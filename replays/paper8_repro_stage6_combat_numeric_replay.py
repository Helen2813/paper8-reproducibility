#!/usr/bin/env python3
"""Paper 8 / Stage 6: numeric ComBat replay against locked archival output.

Normal invocation (macOS):
  python3 "$HOME/Downloads/paper8_repro_stage6_combat_numeric_replay.py"

INPUT POLICY: read only from ~/Desktop/paper8/endo-research.
OUTPUT POLICY: write only into a new, unique ~/Downloads/Paper8_repro_stage6_combat/Run_*.
NO installs, downloads, RMA, classifier fitting, hyperparameter search, or research-tree edits.

An independent R worker reconstructs the frozen full 44,341 x 216 ComBat arm
(sva::ComBat, ref.batch=1, model.matrix(~condition)) and compares every value
against the already archived post-ComBat output. This is not a training-only
transformation and is NOT an independent diagnostic validation.

Only a checked result report (not a copied research matrix) is saved. A successful
archived SHA-256 check is NOT itself a numerical replay.
"""
from __future__ import print_function
import csv
import datetime
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

ROOT = Path.home() / "Desktop" / "paper8" / "endo-research"
DOWNLOADS = Path.home() / "Downloads"
BASE_OUT = DOWNLOADS / "Paper8_repro_stage6_combat"
PREFLIGHT = DOWNLOADS / "Paper8_repro_stage5_combat_preflight" / "stage5_combat_preflight.json"
STAGE5_STATUS = "PASS_ARCHIVED_COMBAT_INPUTS_HASHED_NOT_RERUN"
FROZEN = "20261108T015502Z"
MATRIX = "data/p8_cross_series_independence/fogg_original_combat/" + FROZEN
SOURCE_SCRIPT = "data/p8_cross_series_independence/fogg_code/df538583652d3359d71959437837ea3e0d1195e6/data_prep/batch_correction.r"
FILES = {
    "precombat": (MATRIX + "/all_phases_all_genes_precombat.tsv", "2570400220ea4107e38771232f112899cee1f93b07686e5a07ff8d799d303314"),
    "postcombat": (MATRIX + "/all_phases_all_genes_counts.tsv", "45337e9db40226814cfe314c57b4771c8d1993a733c0e145805b7f75cee09b7c"),
    "metadata": (MATRIX + "/all_phases_coldata.tsv", "a2696216f28f6f80a8e7db716f280c45e4d2e4264730723fb3775b94751ebc4a"),
    "author_R": (SOURCE_SCRIPT, "20366a374148d732602894f80a53e126224c278dc9492dfbe560d4097efab6be"),
}
EXPECTED_INPUT_SHA = FILES["precombat"][1]
EXPECTED_OUTPUT_SHA = FILES["postcombat"][1]
R_BIN = ROOT / ".envs/p8_fogg_r44/bin/Rscript"
R_SCRIPT = r'''#!/usr/bin/env Rscript
# Paper 8 / Stage 6 / frozen full-matrix ComBat replay.
# Inputs read-only. Metrics and environment files written only to supplied Downloads directory.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop("STAGE6_EXPECTED_FOUR_ARGUMENTS")
pre_file <- args[1]; post_file <- args[2]; meta_file <- args[3]; out_dir <- args[4]
if (!requireNamespace("sva", quietly=TRUE)) stop("HOLD_SVA_MISSING")
if (as.character(utils::packageVersion("sva")) != "3.54.0") stop("HOLD_SVA_VERSION_NOT_LOCKED_3_54_0")
message("[R_START] P8 Stage6 numeric replay: full 44341 x 216, no models")
message("[R_PACKAGE] sva ", as.character(utils::packageVersion("sva")))

# Deliberately do not use row names: HGNC symbols can occur more than once.
meta <- utils::read.delim(meta_file, sep="\t", header=TRUE,
                          check.names=FALSE, stringsAsFactors=FALSE,
                          quote="", comment.char="")
if (nrow(meta) != 216L || !identical(names(meta), c("sample_name","condition","phase","batch"))) {
    stop("HOLD_METADATA_SCHEMA")
}
if (anyNA(meta) || anyDuplicated(meta$sample_name) != 0L) stop("HOLD_METADATA_MISSING_OR_DUPLICATE")
pre <- utils::read.delim(pre_file, sep="\t", header=TRUE,
                         check.names=FALSE, stringsAsFactors=FALSE,
                         quote="", comment.char="")
if (!identical(dim(pre), c(44341L,217L))) stop("HOLD_PRECOMBAT_SHAPE")
if (names(pre)[1] != "hgnc_symbol") stop("HOLD_PRECOMBAT_SYMBOL_COLUMN")
if (!identical(as.character(names(pre)[-1]), as.character(meta$sample_name))) stop("HOLD_PRECOMBAT_GSM_ORDER")
if (anyNA(pre[[1]])) stop("HOLD_PRECOMBAT_GENE_IDS_MISSING")
# Avoid character coersion across the entire 216-column data frame.
x <- as.matrix(pre[,-1,drop=FALSE])
if (!is.numeric(x) || any(!is.finite(x))) stop("HOLD_PRECOMBAT_NONNUMERIC_OR_NONFINITE")
colnames(x) <- as.character(meta$sample_name)
gene_ids <- as.character(pre[[1]])
rm(pre); invisible(gc(verbose=FALSE))

# Match the previously locked original full-cohort design.
# Both training and test samples enter ComBat; TRUE test conditions are supplied.
# This is not leakage-free prospective preprocessing.
if (!all(as.character(meta$batch) %in% c("1","2","3","4"))) {
    message("[INFO] batch levels ", paste(unique(as.character(meta$batch)),collapse=","))
}
if (!("1" %in% unique(as.character(meta$batch)))) stop("HOLD_REFERENCE_BATCH_1_MISSING")
mod <- stats::model.matrix(~condition, data=meta)
if (nrow(mod) != 216L || any(!is.finite(mod))) stop("HOLD_DESIGN_INVALID")
if (qr(mod)$rank != ncol(mod)) stop("HOLD_DESIGN_RANK_DEFICIENT")
message("[R_DESIGN] ",paste(colnames(mod), collapse=","),
        " | batch_counts=",paste(as.character(table(meta$batch)),collapse=","))
message("[R_COMBAT_START] all genes; outcome-informed; ref.batch=1")
start <- proc.time()[[3]]
corrected <- sva::ComBat(dat=x, batch=meta$batch, mod=mod,
                         ref.batch=1, par.prior=TRUE, prior.plots=FALSE)
spent <- proc.time()[[3]] - start
if (!identical(dim(corrected), c(44341L,216L)) || any(!is.finite(corrected))) {
    stop("HOLD_REPLAY_NONFINITE_OR_SHAPE")
}
rm(x); invisible(gc(verbose=FALSE))
message("[R_COMBAT_DONE] elapsed_sec=",format(spent,digits=8))

# The archived numeric values are compared cell by cell. Text-file SHA256
# equality is NOT required because R numeric serialization can differ.
reference <- utils::read.delim(post_file, sep="\t", header=TRUE,
                               check.names=FALSE, stringsAsFactors=FALSE,
                               quote="", comment.char="")
if (!identical(dim(reference),c(44341L,217L)) || names(reference)[1]!="symbol") {
    stop("HOLD_ARCHIVED_POSTCOMBAT_SCHEMA")
}
if (!identical(names(reference)[-1], as.character(meta$sample_name))) {
    stop("HOLD_POSTCOMBAT_GSM_ORDER")
}
if (!identical(as.character(reference[[1]]),gene_ids)) stop("HOLD_GENE_ROW_ORDER")

max_abs <- 0.0; sum_abs <- 0.0; sum_sq <- 0.0
above_1e6 <- 0; above_1e8 <- 0; above_1e10 <- 0
n_numeric <- as.double(nrow(corrected)) * ncol(corrected)
worst_gene <- NA_character_; worst_gsm <- NA_character_; worst_gene_row <- NA_integer_
for (j in seq_len(ncol(corrected))) {
    archived <- reference[[j+1L]]
    if (!is.numeric(archived) || any(!is.finite(archived))) stop("HOLD_REFERENCE_NONNUMERIC_OR_NONFINITE")
    err <- abs(corrected[,j] - archived)
    local_max <- max(err)
    if (local_max > max_abs) {
        idx <- which.max(err)
        max_abs <- local_max
        worst_gene <- gene_ids[[idx]]
        worst_gsm <- colnames(corrected)[[j]]
        worst_gene_row <- idx
    }
    sum_abs <- sum_abs + sum(err)
    sum_sq <- sum_sq + sum(err*err)
    above_1e6 <- above_1e6 + sum(err > 1e-6)
    above_1e8 <- above_1e8 + sum(err > 1e-8)
    above_1e10 <- above_1e10 + sum(err > 1e-10)
    if ((j %% 54L)==0L) message("[R_COMPARE] ",j,"/216 samples")
}
rm(reference,corrected); invisible(gc(verbose=FALSE))
metrics <- data.frame(
    key=c("status","n_genes","n_samples","n_values","sva_version",
          "R_version","design_columns","reference_batch","computation_seconds",
          "max_abs_error","mean_abs_error","rmse", "n_abs_gt_1e-6", "n_abs_gt_1e-8",
          "n_abs_gt_1e-10","worst_gene","worst_gsm","worst_gene_row"),
    value=as.character(c(if (above_1e6 == 0L) "PASS_NUMERIC_COMBAT_REPLAY_1E6" else "HOLD_NUMERIC_DIFFERENCE",
          "44341","216",format(n_numeric,scientific=FALSE),
          as.character(utils::packageVersion("sva")),R.version.string,
          paste(colnames(mod),collapse=";"),"1",format(spent,digits=16),
          format(max_abs,digits=16,scientific=TRUE),format(sum_abs/n_numeric,digits=16,scientific=TRUE),
          format(sqrt(sum_sq/n_numeric),digits=16,scientific=TRUE),
          format(above_1e6,scientific=FALSE),format(above_1e8,scientific=FALSE),
          format(above_1e10,scientific=FALSE),worst_gene,worst_gsm,worst_gene_row)),
    stringsAsFactors=FALSE
)
utils::write.table(metrics, file.path(out_dir,"numeric_metrics.tsv"),
                   sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
message("[R_RESULT] ",metrics$value[[1]]," | max_abs=",format(max_abs,digits=16),
        " | mismatch_1e-6=",above_1e6)
'''


def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(4*1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def make_out():
    BASE_OUT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = BASE_OUT / ("Run_" + stamp)
    counter = 1
    while out.exists():
        out = BASE_OUT / ("Run_" + stamp + "_" + str(counter))
        counter += 1
    out.mkdir()
    (out / "temp").mkdir()
    return out


def read_metrics(path):
    with path.open("r",encoding="utf-8",newline="") as f:
        return {r["key"]:r["value"] for r in csv.DictReader(f,delimiter="\t")}


def bounded_log_tail(path, limit=2000):
    try:
        with path.open("rb") as f:
            f.seek(0,os.SEEK_END)
            f.seek(max(0,f.tell()-limit))
            return f.read().decode("utf-8",errors="replace")
    except OSError:
        return ""


def validate_stage5():
    if not PREFLIGHT.is_file():
        raise RuntimeError("MISSING_STAGE5_REPORT: " + str(PREFLIGHT))
    data = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    if data.get("report_status") != STAGE5_STATUS:
        raise RuntimeError("HOLD_STAGE5_NOT_PASS")
    candidate = next((x for x in data.get("candidates",[]) if x.get("run")==FROZEN),None)
    if not candidate or candidate.get("report_status")!="PASS_ORIGINAL_FOGG_COMBAT_216":
        raise RuntimeError("HOLD_FROZEN_RUN_REPORT_MISMATCH")
    matched = {Path(item.get("relative_path", "")).name: item for item in candidate.get("files",[])}
    for key,name in [("precombat","all_phases_all_genes_precombat.tsv"),
                     ("postcombat","all_phases_all_genes_counts.tsv"),
                     ("metadata","all_phases_coldata.tsv")]:
        item=matched.get(name)
        if not item or item.get("hash_matches") is not True or item.get("expected_sha256")!=FILES[key][1]:
            raise RuntimeError("HOLD_STAGE5_INPUT_MISMATCH_"+key)
    if (data.get("r_environment",{}).get("author_batch_script",{}).get("sha256")
            != FILES["author_R"][1]):
        raise RuntimeError("HOLD_STAGE5_ORIGINAL_R_SHA_MISMATCH")
    return data


def selftest():
    assert len(FILES)==4
    assert "sva::ComBat" in R_SCRIPT
    assert "model.matrix(~condition" in R_SCRIPT
    assert "44341" in R_SCRIPT and "216" in R_SCRIPT
    assert R_SCRIPT.count("write.table(") == 1
    assert all(len(h)==64 for _,h in FILES.values())
    assert "ref.batch=1" in R_SCRIPT and "prior.plots=FALSE" in R_SCRIPT
    print("[SELFTEST] PASS_STAGE6_TEMPLATE_AND_LOCKS")


def main():
    if sys.argv[1:]==["--selftest"]:
        selftest(); return 0
    if len(sys.argv)!=1:
        raise SystemExit("Run without flags; or use --selftest")
    print("PAPER 8 | STAGE 6 | INDEPENDENT FULL-GENE COMBAT NUMERIC REPLAY",flush=True)
    print("[POLICY] READ ONLY PROJECT | NEW OUTPUTS ONLY IN DOWNLOADS | NO FITS",flush=True)
    out=make_out()
    report={"stage":"P8_STAGE6_FULL_NUMERIC_COMBAT_REPLAY","status":"HOLD_UNSTARTED",
            "created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "policy":"PROJECT_READ_ONLY; NO_DOWNLOADS_OR_INSTALLS; OUTPUTS_DOWNLOADS_ONLY",
            "source":"original author code verified but independent frozen-method R worker executed",
            "original_author_code_executed":False,
            "entire_GEO_RMA_recomputed":False,"ComBat_numerical_replay_performed":False,
            "classifier_fitted":False,"test_outcomes_supplied_to_ComBat":True,
            "n_samples":216,"n_genes":44341,"inputs":{},"errors":[]}
    try:
        if not ROOT.is_dir(): raise RuntimeError("HOLD_PROJECT_ROOT_MISSING")
        validate_stage5()
        print("[GATE] STAGE5_ARCHIVED_SHA PASS",flush=True)
        if not R_BIN.is_file(): raise RuntimeError("HOLD_RSCRIPT_MISSING")
        # Guard against accidentally beginning a full-matrix R job on an exhausted disk.
        free = shutil.disk_usage(out).free
        if free < 2*1024**3:
            raise RuntimeError("HOLD_FREE_DISK_BELOW_2_GiB")
        for key,(rel,wanted) in FILES.items():
            path=ROOT/rel
            if not path.is_file(): raise RuntimeError("HOLD_FILE_MISSING_"+key)
            print("[HASH]",key,"...",flush=True)
            observed=sha256_file(path)
            report["inputs"][key]={"sha256":observed,"expected_sha256":wanted,
                                    "size_bytes":path.stat().st_size,"match":observed==wanted,
                                    "source_relative_path":rel}
            if observed!=wanted: raise RuntimeError("HOLD_HASH_MISMATCH_"+key)
            print("[HASH]",key,"MATCH",flush=True)
        worker=out/"stage6_combat_replay.R"
        worker.write_text(R_SCRIPT,encoding="utf-8")
        report["worker_sha256"]=sha256_file(worker)
        env=os.environ.copy()
        env["TMPDIR"]=str(out/"temp")
        env["TMP"]=str(out/"temp")
        env["TEMP"]=str(out/"temp")
        cmd=[str(R_BIN),"--vanilla",str(worker),
             str(ROOT/FILES["precombat"][0]),str(ROOT/FILES["postcombat"][0]),
             str(ROOT/FILES["metadata"][0]),str(out)]
        report["rscript_binary_sha256"]=sha256_file(R_BIN)
        report["r_command_executable"]=".envs/p8_fogg_r44/bin/Rscript"
        print("[R_START] full-gene outcome-informed ComBat, comparison across all 9,577,656 values",flush=True)
        # Stream R progress to the Terminal and simultaneously retain the full log.
        # Do not store raw expression values in logs.
        with (out/"R_run.log").open("w",encoding="utf-8") as log:
            process=subprocess.Popen(cmd,cwd=str(out),env=env,stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,text=True,bufsize=1,
                                     errors="replace")
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                log.flush()
                if line.startswith("[R_") or line.startswith("[INFO]"):
                    print(line.rstrip(),flush=True)
            rc=process.wait()
        report["R_exit_code"]=rc
        report["ComBat_numerical_replay_performed"]=(rc==0)
        print("[R_EXIT]",rc,flush=True)
        if rc!=0:
            print("[R_ERROR_TAIL]",bounded_log_tail(out/"R_run.log"),flush=True)
            raise RuntimeError("HOLD_R_WORKER_EXIT_"+str(rc))
        metricpath=out/"numeric_metrics.tsv"
        if not metricpath.is_file(): raise RuntimeError("HOLD_NUMERIC_METRICS_NOT_CREATED")
        metrics=read_metrics(metricpath)
        report["metrics"]=metrics
        if metrics.get("n_genes")!="44341" or metrics.get("n_samples")!="216" or metrics.get("n_abs_gt_1e-6") is None:
            raise RuntimeError("HOLD_NUMERIC_METRICS_SCHEMA")
        report["status"] = ("PASS_NUMERIC_COMBAT_REPLAY_1E6" if
            metrics.get("status")=="PASS_NUMERIC_COMBAT_REPLAY_1E6" and
            int(metrics["n_abs_gt_1e-6"])==0 else "HOLD_NUMERIC_DIFFERENCE")
        print("[MAX_ABS]",metrics.get("max_abs_error"),flush=True)
        print("[MISMATCHES_GT_1E6]",metrics.get("n_abs_gt_1e-6"),flush=True)
    except Exception as e:
        report["errors"].append(str(e))
        report["status"]="HOLD_STAGE6_ERROR"
        print("[HOLD]",str(e),flush=True)
    finally:
        report["completed_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
        (out/"stage6_replay.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        (out/"summary.txt").write_text(
            "PAPER 8 / STAGE 6 / COMBAT NUMERIC REPLAY\n"
            "Status: "+report["status"]+"\n"
            "Max absolute error: "+str(report.get("metrics",{}).get("max_abs_error","NOT_COMPUTED"))+"\n"
            "Cells differing > 1e-6: "+str(report.get("metrics",{}).get("n_abs_gt_1e-6","NOT_COMPUTED"))+"\n"
            "Note: full 216-array joint ComBat uses true test diagnoses; not a leakage-free validation.\n",
            encoding="utf-8")
        print("[OUTPUT]",out,flush=True)
        print("[STATUS]",report["status"],flush=True)
    return 0 if report["status"]=="PASS_NUMERIC_COMBAT_REPLAY_1E6" else 1


if __name__=="__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("[INTERRUPTED] no original project data changed",file=sys.stderr)
        sys.exit(130)
