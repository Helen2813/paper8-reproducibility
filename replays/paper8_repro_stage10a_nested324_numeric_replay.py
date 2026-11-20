#!/usr/bin/env python3
"""Paper 8 / Stage 10A: independent replay of frozen GSE51981 nested-CV benchmark.

READ-ONLY original project. Reuses the archived 148 x 324 single-CEL MAS5 matrix;
no raw CEL, ComBat, marker selection, or downloads. All outputs go to Downloads.
This is an internal/sample-level retrospective benchmark, NOT external validation.

Invocation (no flags): python3 ~/Downloads/paper8_repro_stage10a_nested324_numeric_replay.py
"""
from __future__ import print_function
import collections
import csv
import datetime
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import traceback
import zipfile

SEED = 20261114
C_GRID = [0.01, 0.1, 1.0, 10.0]
TASKS = [
    ("ENDO_vs_ALL", ("ENDO", "OTHER", "NUPP")),
    ("ENDO_vs_OTHER", ("ENDO", "OTHER")),
    ("ENDO_vs_NUPP", ("ENDO", "NUPP")),
]
ROOT = Path.home() / "Desktop" / "paper8" / "endo-research"
SOURCE = ROOT / "results" / "p8_v422_locked_GSE51981_MAS5_nested_cv_RESUMABLE.zip"
BASE = Path.home() / "Downloads" / "Paper8_repro_stage10A_nested324"
CURRENT_SCRIPT = Path(__file__).resolve()
EXPECTED_ARCHIVE_STATUS = "PASS_V422_LOCKED_INTERNAL_ANALYSIS_WITH_HOLDS"
EXPECTED_TASK_COUNTS = {
    "ENDO_vs_ALL": (148, 77, 71),
    "ENDO_vs_OTHER": (114, 77, 37),
    "ENDO_vs_NUPP": (111, 77, 34),
    "ENDO_vs_ALL_2009_ONLY": (77, 46, 31),
    "ENDO_vs_OTHER_2009_ONLY": (77, 46, 31),
}
# Published internal-benchmark checkpoints, not training/tuning targets.
PUBLISHED_CORE_AUC = {
    "ENDO_vs_ALL": 0.7556246570331078,
    "ENDO_vs_OTHER": 0.7328887328887328,
    "ENDO_vs_NUPP": 0.8934300993124522,
    "ENDO_vs_OTHER_2009_ONLY": 0.7237026647966339,
}
TOL_SCORES = 1e-8
TOL_METRICS = 1e-8

class Hold(Exception):
    pass

def log(*parts):
    print(*parts, flush=True)

def sha_bytes(blob):
    return hashlib.sha256(blob).hexdigest()

def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for b in iter(lambda: stream.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()

def load_table(blob):
    with io.StringIO(blob.decode("utf-8-sig"), newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise Hold("INVALID_TABLE_HEADER")
        return list(reader), list(reader.fieldnames)

def write_csv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def create_output_dir():
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = BASE / ("Run_" + stamp)
    path.mkdir(parents=True, exist_ok=False)
    return path

def load_inputs():
    if not SOURCE.is_file() or SOURCE.is_symlink():
        raise Hold("MISSING_ARCHIVED_V422_ZIP: " + str(SOURCE))
    required = (
        "00_manifest.json", "05_148_x_324_single_CEL_MAS5_intensity.tsv",
        "06_oof_predictions_NOT_EXTERNAL.tsv", "07_outer_fold_manifest.tsv",
        "08_nested_cv_summary.json",
    )
    with zipfile.ZipFile(SOURCE) as zipfile_in:
        unsafe = [n for n in zipfile_in.namelist()
                  if n.startswith("/") or ".." in Path(n).parts or "\\" in n]
        if unsafe: raise Hold("UNSAFE_ARCHIVE_MEMBER")
        if zipfile_in.testzip() is not None: raise Hold("ZIP_CRC_FAILURE")
        for member in required:
            if member not in zipfile_in.namelist():
                raise Hold("MISSING_ARCHIVED_MEMBER_" + member)
        files = {f: zipfile_in.read(f) for f in required}
    manifest = json.loads(files["00_manifest.json"])
    if manifest.get("status") != EXPECTED_ARCHIVE_STATUS:
        raise Hold("UNEXPECTED_ARCHIVED_STATUS_" + str(manifest.get("status")))
    for member in required:
        if not member.endswith(".tsv"): continue
        expected = manifest.get("output_sha256", {}).get(member)
        if not expected or expected != sha_bytes(files[member]):
            raise Hold("ARCHIVED_SHA_MISMATCH_" + member)
    matrix_rows, matrix_header = load_table(files[required[1]])
    pred_rows, pred_header = load_table(files[required[2]])
    folds, fold_header = load_table(files[required[3]])
    summary = json.loads(files[required[4]])
    if matrix_header[0] != "gsm" or len(matrix_header) != 325:
        raise Hold("FROZEN_MATRIX_COLUMN_SCHEMA")
    if len(set(matrix_header[1:])) != 324:
        raise Hold("DUPLICATE_PROBE_IDS")
    if len(matrix_rows) != 148:
        raise Hold("FROZEN_MATRIX_NOT_148")
    features = {}
    for row in matrix_rows:
        gsm = row["gsm"].strip().upper()
        if not gsm or gsm in features: raise Hold("DUPLICATE_OR_EMPTY_GSM")
        values = [float(row[probe]) for probe in matrix_header[1:]]
        if not all(math.isfinite(v) and v > 0 for v in values):
            raise Hold("NONPOSITIVE_OR_MISSING_MAS5_VALUES")
        features[gsm] = values
    if len(pred_rows) < 148 or len(folds) < 5:
        raise Hold("INCOMPLETE_ARCHIVED_OOF_OR_FOLDS")
    groups, dates, phases, donors = {}, {}, {}, {}
    for row in pred_rows:
        gsm = row["gsm"].upper().strip()
        if gsm not in features: raise Hold("OOF_GSM_NOT_IN_MATRIX")
        for field, target in (("group",groups),("date",dates),("phase",phases),
                              ("donor_id_verified",donors)):
            value = row[field]
            if gsm in target and target[gsm] != value:
                raise Hold("INCONSISTENT_" + field + "_FOR_" + gsm)
            target[gsm] = value
    if len(groups) != 148 or len(dates) != 148:
        raise Hold("OOF_NOT_ALL_148_GSM")
    if collections.Counter(groups.values()) != collections.Counter(
            {"ENDO":77,"OTHER":37,"NUPP":34}):
        raise Hold("THREE_GROUP_COUNTS_MISMATCH")
    unique_donors = {v.strip() for v in donors.values() if v.strip()}
    if unique_donors and len(unique_donors) != 148:
        # A donor-grouped approach needs verified crosswalk and grouping, not a partial ID field.
        raise Hold("DONOR_ID_SCHEMA_REQUIRES_SEPARATE_GROUPED_REPLAY")
    if len(set(dates.values())) != 17:
        raise Hold("SCANNER_DATE_COUNT_MISMATCH")
    indexed_pred = {}
    for row in pred_rows:
        key = (row["task"], row["mode"], row["gsm"].upper())
        if key in indexed_pred: raise Hold("DUPLICATED_ARCHIVED_PREDICTION")
        indexed_pred[key] = row
    indexed_fold = {}
    for row in folds:
        key = (row["task"], row["mode"], int(row["fold"]))
        if key in indexed_fold: raise Hold("DUPLICATED_ARCHIVED_FOLD")
        indexed_fold[key] = row
    results = {(row["task"],row["mode"]):row for row in summary["results"]}
    return {
        "source_sha256": sha_file(SOURCE), "manifest": manifest, "features": features,
        "groups": groups, "dates": dates, "phases": phases,
        "predictions": indexed_pred, "folds": indexed_fold, "results": results,
        "source_members": {k:sha_bytes(v) for k,v in files.items()},
        "summary": summary,
    }

def build_folds(y, groupvec, strategy, n_splits=5, offset=0):
    import numpy as np
    from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
    seed = SEED + offset
    y = np.asarray(y)
    if strategy == "sample":
        cv = StratifiedKFold(n_splits=n_splits,shuffle=True,random_state=seed)
        return list(cv.split(np.zeros((len(y),1)),y))
    cv = StratifiedGroupKFold(n_splits=n_splits,shuffle=True,random_state=seed)
    return list(cv.split(np.zeros((len(y),1)),y,np.asarray(groupvec)))

def validate_folds(folds, y, groupvec=None, n_expected=5):
    if len(folds) != n_expected: return False
    alltests = []
    for train, test in folds:
        if len(train)==0 or len(test)==0: return False
        if len(set(y[train])) != 2 or len(set(y[test])) != 2:return False
        if groupvec is not None and set(groupvec[train]) & set(groupvec[test]):
            return False
        alltests.extend(test.tolist())
    return sorted(alltests) == list(range(len(y)))

def evaluate_task(task, mode, allowed, data, gsms_override=None):
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score, balanced_accuracy_score, confusion_matrix, brier_score_loss
    groups, dates, phases = data["groups"], data["dates"], data["phases"]
    if gsms_override is None:
        gsms = sorted(g for g,v in groups.items() if v in allowed)
    else:
        gsms = sorted(g for g in gsms_override if groups[g] in allowed)
    y = np.asarray([int(groups[g]=="ENDO") for g in gsms],dtype=int)
    x = np.log2(np.asarray([data["features"][g] for g in gsms],dtype=float))
    if x.shape != (len(gsms),324) or not np.isfinite(x).all():raise Hold("INVALID_TASK_FEATURES_"+task)
    if mode == "core":
        strategy, groupvec = "sample", None
    else:
        strategy = "scanner_date"
        groupvec = np.asarray([dates[g] for g in gsms])
    try:
        outer = build_folds(y, groupvec,strategy)
    except ValueError as e:
        return {"status":"HOLD_OUTER_FOLD_FAILURE", "reason":str(e)}, [], []
    if not validate_folds(outer,y,groupvec):
        return {"status":"HOLD_OUTER_FOLD_CLASS_OR_GROUP_POSITIVITY"}, [], []
    pred=np.full(len(y),np.nan)
    prior=np.full(len(y),np.nan)
    phasepred=np.full(len(y),np.nan)
    negpred=np.full(len(y),np.nan)
    folds=[]
    phase_all=[phases[g] or "UNKNOWN" for g in gsms]
    def estimator(c):
        return make_pipeline(StandardScaler(),LogisticRegression(
            C=c,class_weight="balanced",solver="liblinear",max_iter=2000,
            random_state=SEED))
    for k,(tr,te) in enumerate(outer,1):
        yt=y[tr]
        gx=groupvec[tr] if groupvec is not None else None
        try:
            inner = build_folds(yt,gx,strategy,4,k)
        except ValueError as e:
            return {"status":"HOLD_INNER_FOLD_FAILURE","reason":str(e)},[],[]
        if not validate_folds(inner,yt,gx,n_expected=4):
            return {"status":"HOLD_INNER_CLASS_OR_GROUP_POSITIVITY"},[],[]
        cand=[]
        for c in C_GRID:
            score=np.full(len(tr),np.nan)
            for a,b in inner:
                model=estimator(c)
                model.fit(x[tr][a],yt[a])
                score[b]=model.predict_proba(x[tr][b])[:,1]
            cand.append((float(roc_auc_score(yt,score)),c))
        selected_c=sorted(cand,key=lambda v:(-v[0],v[1]))[0][1]
        model=estimator(selected_c)
        model.fit(x[tr],yt)
        pred[te]=model.predict_proba(x[te])[:,1]
        prior[te]=float(np.mean(yt))
        nphase=collections.Counter(phase_all[j] for j in tr)
        ncase=collections.Counter(phase_all[j] for j in tr if y[j]==1)
        phasepred[te]=np.asarray([
            (ncase[phase_all[j]]+1)/(nphase[phase_all[j]]+2)
            if nphase[phase_all[j]] else (sum(yt)+1)/(len(yt)+2) for j in te])
        rng=np.random.default_rng(SEED+k)
        yp=rng.permutation(yt)
        neg=estimator(selected_c)
        neg.fit(x[tr],yp)
        negpred[te]=neg.predict_proba(x[te])[:,1]
        folds.append({"task":task,"mode":mode,"fold":k,"n_train":len(tr),
                      "n_test":len(te),"train_endo":int(sum(yt)),
                      "test_endo":int(sum(y[te])),"best_C":selected_c,
                      "heldout_dates":";".join(sorted({dates[gsms[j]] for j in te}))})
    if not all(np.isfinite(p).all() for p in (pred,prior,phasepred,negpred)):
        raise Hold("UNFILLED_PREDICTIONS_"+task)
    def measures(p):
        call=(p>=.5).astype(int)
        tn,fp,fn,tp=confusion_matrix(y,call,labels=[0,1]).ravel()
        return {"roc_auc":float(roc_auc_score(y,p)),
                "BA":float(balanced_accuracy_score(y,call)),
                "sensitivity":float(tp/(tp+fn)),
                "specificity":float(tn/(tn+fp)),
                "precision":float(tp/(tp+fp)) if tp+fp else 0.0,
                "Brier":float(brier_score_loss(y,p)),
                "TP":int(tp),"FP":int(fp),"TN":int(tn),"FN":int(fn)}
    metrics={"status":"PASS_OOF_EXPLORATORY","task":task,"mode":mode,
             "n":len(gsms),"endo":int(sum(y)),"nonendo":int(len(y)-sum(y)),
             "point":measures(pred),"class_prior_baseline":measures(prior),
             "phase_only_baseline":measures(phasepred),
             "permuted_training_negative_control":measures(negpred)}
    rows=[{"task":task,"mode":mode,"gsm":g,
           "OOF_probability":float(pred[i]),"prior_probability":float(prior[i]),
           "phase_only_probability":float(phasepred[i]),
           "permuted_train_probability":float(negpred[i])}
          for i,g in enumerate(gsms)]
    return metrics,rows,folds

def check_pairwise_archive(data, new_results, preds, folds):
    mismatches=[]
    pred_comparison=[]
    fold_comparison=[]
    metrics_comparison=[]
    max_score_error=0.0
    all_score_keys=("OOF_probability","prior_probability",
                    "phase_only_probability","permuted_train_probability")
    expected_keys=set(data["results"])
    if set(new_results)!=expected_keys:
        mismatches.append("RESULT_SET_MISMATCH")
    for key,result in new_results.items():
        stored=data["results"].get(key)
        if stored is None:
            mismatches.append("MISSING_ARCHIVED_RESULT_"+str(key));continue
        entry={"task":key[0],"mode":key[1],"new_status":result["status"],
               "archived_status":stored["status"]}
        if result["status"]!=stored["status"]:
            mismatches.append("STATUS_"+str(key))
        if "point" in result:
            if "point" not in stored:
                mismatches.append("MISSING_ARCHIVED_POINT_"+str(key))
            else:
                for group in ("point","class_prior_baseline","phase_only_baseline",
                              "permuted_training_negative_control"):
                    for name,value in result[group].items():
                        ref=stored[group][name]
                        diff=abs(float(value)-float(ref))
                        if diff > TOL_METRICS:
                            mismatches.append("METRIC_"+str(key)+"_"+group+"_"+name)
                entry["auc"]=result["point"]["roc_auc"]
                entry["ba"]=result["point"]["BA"]
        metrics_comparison.append(entry)
    actual_pred={}
    for row in preds:
        key=(row["task"],row["mode"],row["gsm"])
        if key in actual_pred: raise Hold("DUPLICATE_NEW_OOF_ROW")
        actual_pred[key]=row
    if set(actual_pred)!=set(data["predictions"]):
        mismatches.append("OOF_KEYS_MISMATCH")
    for key in sorted(set(actual_pred)&set(data["predictions"])):
        cur,ref=actual_pred[key],data["predictions"][key]
        record={"task":key[0],"mode":key[1],"gsm":key[2]}
        for field in all_score_keys:
            diff=abs(float(cur[field])-float(ref[field]))
            max_score_error=max(max_score_error,diff)
            record[field+"_abs_diff"]=format(diff,'.17g')
            if diff>TOL_SCORES:
                mismatches.append("OOF_"+str(key)+"_"+field)
        pred_comparison.append(record)
    actual_folds={}
    for row in folds:
        key=(row["task"],row["mode"],int(row["fold"]))
        if key in actual_folds: raise Hold("DUPLICATE_NEW_FOLD_ROW")
        actual_folds[key]=row
    if set(actual_folds)!=set(data["folds"]):
        mismatches.append("FOLD_KEYS_MISMATCH")
    for key in sorted(set(actual_folds)&set(data["folds"])):
        cur,ref=actual_folds[key],data["folds"][key]
        different=[]
        for field in ("n_train","n_test","train_endo","test_endo"):
            if int(cur[field])!=int(ref[field]):different.append(field)
        if abs(float(cur["best_C"])-float(ref["best_C"]))>1e-12:
            different.append("best_C")
        if cur["heldout_dates"]!=ref["heldout_dates"]:
            different.append("heldout_dates")
        if different:
            mismatches.append("FOLD_"+str(key)+"_"+";".join(different))
        fold_comparison.append({"task":key[0],"mode":key[1],"fold":key[2],
                                "match":not bool(different),
                                "differences":";".join(different),
                                "frozen_C":ref["best_C"],"replayed_C":cur["best_C"]})
    for task,value in PUBLISHED_CORE_AUC.items():
        key=(task,"core")
        if key not in new_results or "point" not in new_results[key]:
            mismatches.append("MISSING_PUBLISHED_CHECKPOINT_"+task)
        elif abs(new_results[key]["point"]["roc_auc"]-value)>TOL_METRICS:
            mismatches.append("PUBLISHED_AUC_MISMATCH_"+task)
    # Never silently collapse multiple mismatches in report.
    return mismatches,pred_comparison,fold_comparison,metrics_comparison,max_score_error

def process(out):
    import numpy as np
    import sklearn
    data=load_inputs()
    log("[ARCHIVE] CRC + archived table SHA PASS",SOURCE.name)
    log("[INPUT] 148 GSM x 324 positive MAS5 intensities; three diagnosis groups; 17 scanner dates")
    np.seterr(all="raise")
    new_results={};allpred=[];allfold=[]
    for task,allowed in TASKS:
        for mode in ("core","dateblocked"):
            result,pred,fold=evaluate_task(task,mode,allowed,data)
            new_results[(task,mode)]=result
            allpred.extend(pred);allfold.extend(fold)
            log("[REPLAY]",task,mode,result["status"],
                "AUC="+str(round(result["point"]["roc_auc"],9)) if "point" in result else "")
    # Predeclared 2009 subset from original v4.22. Repeat only when the
    # original fixed class-positivity gate is met (>=10 each and >=40 total).
    year2009={g for g in data["groups"] if data["dates"][g].startswith("2009-")}
    for task,allowed in TASKS:
        subtask=task+"_2009_ONLY"
        eligible=[g for g in year2009 if data["groups"][g] in allowed]
        positives=sum(data["groups"][g]=="ENDO" for g in eligible)
        negatives=len(eligible)-positives
        if min(positives,negatives)<10 or len(eligible)<40:
            result={"status":"HOLD_2009_PREDECLARED_MIN_10_EACH_40_TOTAL",
                    "task":subtask,"mode":"core","n":len(eligible),
                    "endo":positives,"nonendo":negatives}
            pred,fold=[],[]
        else:
            result,pred,fold=evaluate_task(subtask,"core",allowed,data,eligible)
        new_results[(subtask,"core")]=result
        allpred.extend(pred);allfold.extend(fold)
        log("[REPLAY]",subtask,"core",result["status"],
            "AUC="+str(round(result["point"]["roc_auc"],9)) if "point" in result else "")
    errors,comp_pred,comp_fold,comp_metric,maxdiff=check_pairwise_archive(
        data,new_results,allpred,allfold)
    write_csv(out/"01_oof_gsm_score_differences.csv",comp_pred,
              ["task","mode","gsm"]+[k+"_abs_diff" for k in
               ("OOF_probability","prior_probability","phase_only_probability",
                "permuted_train_probability")])
    write_csv(out/"02_outer_fold_comparison.csv",comp_fold,
              ["task","mode","fold","match","differences","frozen_C","replayed_C"])
    write_csv(out/"03_task_metric_comparison.csv",comp_metric,
              ["task","mode","new_status","archived_status","auc","ba"])
    report={"stage":"P8_STAGE10A_324_PROBE_NUMERIC_REPLAY",
            "status":"PASS_V422_INTERNAL_BENCHMARK_NUMERIC_REPLAY" if not errors else "HOLD_V422_REPLAY_DIFFERENCES",
            "scope":"Independent model replay from previously frozen 148x324 single-CEL MAS5 matrix; no raw CEL or MAS5 rerun",
            "n_matrix_rows":148,"n_probes":324,"seed":SEED,"grid":C_GRID,
            "reference_zip_sha256":data["source_sha256"],
            "reference_members_sha256":data["source_members"],
            "worker_python":sys.version.split()[0],"numpy_version":np.__version__,
            "sklearn_version":sklearn.__version__,"max_oof_score_diff":maxdiff,
            "n_compared_oof_gsm_rows":len(comp_pred),"n_compared_folds":len(comp_fold),
            "n_compared_result_cells":len(comp_metric),"n_mismatches":len(errors),
            "mismatches":errors[:200],"result_cells":comp_metric,
            "limitations":["Exploratory sample-level internal CV; independence of donors is not established",
                           "193-gene/324-probe panel may have been informed by prior analyses",
                           "This does not reproduce the published Fogg model or external diagnostic accuracy",
                           "This step does not rerun MAS5 from raw CEL or the v4.23 leave-one-scanner-date-out experiment"]}
    (out/"stage10a_nested_replay.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    log("[SCORE_MAX_ABS]",maxdiff)
    log("[FOLDS]",len(comp_fold),"[GSM OOF rows]",len(comp_pred))
    log("[STATUS]",report["status"])
    return report

def run_worker(out):
    report={"status":"HOLD_STAGE10A","stage":"P8_STAGE10A_324_PROBE_NUMERIC_REPLAY"}
    try:
        report=process(out)
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
        report["traceback_tail"]=traceback.format_exc()[-3000:]
        (out/"stage10a_nested_replay.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")
        log("[HOLD]",report["error"])
    (out/"README_SCOPE.txt").write_text(
        "Paper 8 Stage 10A independent nested-CV replay from v4.22 archived MAS5 matrix.\n"
        "Read-only original project. Code has not performed raw-CEL MAS5 re-extraction.\n"
        "This is NOT independent clinical/patient validation.\n"
        "Interpret status HOLD as a mismatch or unavailable archived prerequisite, not a negative clinical result.\n",
        encoding="utf-8")
    outfile=out/"Paper8_stage10A_NESTED324_REPLAY_PRIVATE.zip"
    with zipfile.ZipFile(outfile,"w",zipfile.ZIP_DEFLATED,compresslevel=6) as zipout:
        for f in sorted(out.iterdir()):
            if f.is_file() and f!=outfile:zipout.write(f,f.name)
    with zipfile.ZipFile(outfile) as z:
        if z.testzip() is not None:raise Hold("OUTPUT_ZIP_CRC_FAILURE")
    log("[PRIVATE_ZIP]",outfile)
    return 0 if report["status"]=="PASS_V422_INTERNAL_BENCHMARK_NUMERIC_REPLAY" else 2

def selftest():
    import numpy as np
    yy=np.asarray([0,1]*30)
    folds=build_folds(yy,None,"sample")
    assert validate_folds(folds,yy)
    assert len(C_GRID)==4 and SEED==20261114
    assert sum((77,37,34))==148
    assert len(PUBLISHED_CORE_AUC)==4
    log("[SELFTEST] PASS_FIXED_NESTED_CV_SCHEMATIC")

def main():
    if len(sys.argv)>1 and sys.argv[1]=="--worker":
        return run_worker(Path(sys.argv[2]))
    log("PAPER 8 | STAGE 10A | INDEPENDENT INTERNAL 324-PROBE NESTED CV REPLAY")
    log("[POLICY] ORIGINAL PROJECT READ ONLY | NO MAS5/CEL | NO COMBAT | NO TUNING")
    out=create_output_dir()
    worker=ROOT/".envs"/"p8_v422_py"/"bin"/"python3"
    if not worker.is_file():
        worker=ROOT/".envs"/"p8_fogg_classifier"/"bin"/"python3"
    log("[WORKER]",worker)
    if not worker.is_file():
        (out/"stage10a_nested_replay.json").write_text(json.dumps(
            {"status":"HOLD_NO_LOCKED_PYTHON_ENV"},indent=2)+"\n")
        return 2
    cmd=[str(worker),str(CURRENT_SCRIPT),"--worker",str(out)]
    process_run=subprocess.run(cmd,capture_output=True,text=True)
    (out/"worker_stdout.log").write_text(process_run.stdout,encoding="utf8")
    (out/"worker_stderr.log").write_text(process_run.stderr,encoding="utf8")
    private_kit = out / "Paper8_stage10A_NESTED324_REPLAY_PRIVATE.zip"
    if private_kit.is_file():
        with zipfile.ZipFile(private_kit,"a",zipfile.ZIP_DEFLATED,compresslevel=6) as kit:
            kit.write(out/"worker_stdout.log","worker_stdout.log")
            kit.write(out/"worker_stderr.log","worker_stderr.log")
        with zipfile.ZipFile(private_kit) as kit:
            if kit.testzip() is not None:
                raise Hold("OUTPUT_ZIP_CRC_FAILURE_AFTER_LOGS")
    print(process_run.stdout,end="")
    if process_run.stderr:
        log("[WORKER_STDERR_TAIL]",process_run.stderr[-2000:])
    log("[EXIT]",process_run.returncode)
    log("[OUTPUT]",out)
    return process_run.returncode

if __name__=="__main__":
    sys.exit(main())
