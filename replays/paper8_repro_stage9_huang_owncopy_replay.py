#!/usr/bin/env python3
"""Paper 8 Stage 9A: independently replay v4.30 Huang-inspired RF own-copy omission.

Uses the four pinned historical RESULT ZIPs v4.21, v4.27, v4.28, v4.29 and
compares against the v4.30 historical ZIP. It does NOT need original CEL, perform
MAS5 processing or re-run feature selection. Source archives are READ ONLY.
New reports are created exclusively in ~/Downloads/Paper8_repro_stage9_huang/.

This does NOT claim to reproduce Huang et al.'s original fitted model. The 26
own-copy-out scores come from 26 distinct leave-one-out fitted classifiers; the
pooled AUC is NOT an external AUC for one fixed diagnostic model.

Run: python3 ~/Downloads/paper8_repro_stage9_huang_owncopy_replay.py
Test: python3 ~/Downloads/paper8_repro_stage9_huang_owncopy_replay.py --selftest
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
import traceback
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path.home() / "Desktop" / "paper8" / "endo-research"
DOWNLOADS = Path.home() / "Downloads"
WORK = DOWNLOADS / "Paper8_repro_stage9_huang"
EXPECTED_SOURCE_V429 = "939ecccfda4e65e1f053b0a1e46475de4edff7d09bb9b11e90cceb792b5feb3e"
SOURCE_V430_SHA = "d64c4d0982d01715327ebc0523c7c19c7e4fe41adeaf5d97bf7938c33497a1f7"
GENES = ("ACKR1", "LMNB1", "MFAP4", "NMU", "SEMA3C")
SEED_RF = 1234567
SEED_PLACEBO = 20261115
N_ESTIMATORS = 1000
N_PLACEBO = 199
BUNDLES = {
    "v421": "p8_v421_exact_CEL_overlap_and_scanner_dates_RESUMABLE.zip",
    "v427": "p8_v427_published_markers_date_stratified_RESUMABLE.zip",
    "v428": "p8_v428_huang_5gene_duplicate_sensitivity_RESUMABLE.zip",
    "v429": "p8_v429_huang_training_ablation_RESUMABLE.zip",
    "v430": "p8_v430_targeted_reanalysis_20261115T021900Z.zip",
}

class Hold(Exception):
    pass

def digest(data):
    return hashlib.sha256(data).hexdigest()

def file_digest(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def log(*items):
    print(*items, flush=True)

def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")

def table_from_bytes(blob, label, delimiter="\t"):
    text = blob.decode("utf-8-sig")
    rd = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    if not rd.fieldnames or len(rd.fieldnames) != len(set(rd.fieldnames)):
        raise Hold("Invalid or duplicate columns: " + label)
    result = list(rd)
    if not result:
        raise Hold("Empty table: " + label)
    return result

def write_tsv(path, rows):
    if not rows:
        raise Hold("Attempted to write empty results: " + path.name)
    with path.open("w", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)

class Bundle:
    """Read-only access to pinned archives, using ZIP entries without extraction."""
    def __init__(self, path):
        self.path = path
        if not path.is_file() or path.is_symlink():
            raise Hold("Missing or symbolic-link historical archive: " + path.name)
        self.zip = zipfile.ZipFile(path)
        if self.zip.testzip() is not None:
            raise Hold("Archive CRC failed: " + path.name)
        names = self.zip.namelist()
        if not names or len(set(names)) != len(names):
            raise Hold("Archive names invalid: " + path.name)
        for n in names:
            if n.startswith("/") or ".." in Path(n).parts or "\\" in n:
                raise Hold("Unsafe ZIP member " + n)
        self.names = set(names)

    def data(self, name):
        if name not in self.names:
            raise Hold("Missing " + self.path.name + "!" + name)
        if self.zip.getinfo(name).file_size > 4_000_000:
            raise Hold("Unexpected large archive member: " + name)
        return self.zip.read(name)

    def json(self, name):
        return json.loads(self.data(name).decode("utf-8-sig"))

    def table(self, name):
        return table_from_bytes(self.data(name), self.path.name + "!" + name)

    def close(self):
        self.zip.close()

def find_pairs(rows, train, test):
    pairs = []
    for row in rows:
        gsm = sorted(set(v.strip() for v in row.values()
                         if isinstance(v, str) and re.fullmatch(r"GSM\d+", v.strip())))
        tr = [g for g in gsm if g in train and train[g]["series"] == "GSE6364"]
        te = [g for g in gsm if g in test]
        if len(tr) != 1 or len(te) != 1:
            raise Hold("Ambiguous exact CEL pair")
        pairs.append((tr[0], te[0]))
    if len(pairs) != 26 or len({a for a, _ in pairs}) != 26 or len({b for _, b in pairs}) != 26:
        raise Hold("Exact CEL duplicate pair roster is not one-to-one 26")
    return sorted(pairs, key=lambda r: r[1])

def genes_from_vector(bundle, folder, gsm, probes):
    stem = folder + "/" + gsm
    data = bundle.data(stem + ".tsv")
    meta = bundle.json(stem + ".json")
    if meta.get("output_sha") != digest(data):
        raise Hold("Pinned MAS5 CEL vector hash mismatch: " + gsm)
    rows = table_from_bytes(data, stem)
    vals = {r["probe_id"]: float(r["mas5_signal"]) for r in rows}
    if len(vals) != len(rows) or not probes.issubset(vals):
        raise Hold("Missing or repeated probe in " + gsm)
    if any(not math.isfinite(vals[p]) or vals[p] <= 0 for p in probes):
        raise Hold("Nonfinite or nonpositive MAS5 signal in " + gsm)
    return vals

def initialize(bundles):
    import numpy as np
    b21, b27, b28, b29, b30 = (bundles[k] for k in ("v421", "v427", "v428", "v429", "v430"))
    source = b29.data("07_reproduction_source.py")
    if digest(source) != EXPECTED_SOURCE_V429:
        raise Hold("Wrong historical v4.29 source SHA-256")
    if b28.json("00_manifest.json").get("status") != "PASS_V428_NEW_HUANG_FIVE_GENE_RECONSTRUCTION_SENSITIVITY":
        raise Hold("v4.28 manifest did not PASS")
    if b29.json("00_manifest.json").get("status") != "PASS_V429_HUANG_SAME_TEST_TRAINING_ABLATION":
        raise Hold("v4.29 manifest did not PASS")
    if b30.json("00_Huang_single_copy_status.json").get("status") != "PASS_26_INDIVIDUAL_COPY_OMISSIONS":
        raise Hold("v4.30 archived Huang gate did not PASS")
    if b27.json("01_status_manifest.json").get("status") != "PASS_V427_MARKER_DATE_STRATIFIED_EXPLORATORY":
        raise Hold("v4.27 manifest did not PASS")
    train_rows = b28.table("05_train_59_GEO_title_labels.tsv")
    test_rows = b28.table("03_external_predictions_148.tsv")
    annotation = b28.table("02_fixed_probes_by_gene.tsv")
    if len(train_rows) != 59 or len(test_rows) != 148 or len(annotation) != 5:
        raise Hold("Wrong Huang training/test/annotation sizes")
    train = {r["gsm"]: r for r in train_rows}
    test = {r["gsm"]: r for r in test_rows}
    if len(train) != 59 or len(test) != 148 or set(train) & set(test):
        raise Hold("Duplicate or overlapping GSM roster")
    if Counter(r["series"] for r in train_rows) != Counter({"GSE6364": 37, "GSE25628": 22}):
        raise Hold("Wrong source training roster")
    if Counter(r["class"] for r in train_rows) != Counter({"ENDO": 37, "CONTROL": 22}):
        raise Hold("Wrong training label totals")
    if Counter(r["group"] for r in test_rows) != Counter({"ENDO": 77, "OTHER": 37, "NUPP": 34}):
        raise Hold("Wrong test group totals")
    mapping = {r["gene"]: r["shared_probe_ids"].split(";") for r in annotation}
    if set(mapping) != set(GENES) or any(not p for p in mapping.values()):
        raise Hold("Wrong frozen five-gene panel")
    probes = set(p for lst in mapping.values() for p in lst)
    if len(probes) != 6:
        raise Hold("Expected exactly six shared GPL570/GPL571 probe IDs")
    pairs = find_pairs(b21.table("02_GSE51981_GSE6364_EXACT_CEL_PAIRS.tsv"), train, test)
    tr_copies = {a for a, b in pairs}
    te_copies = {b for a, b in pairs}
    if Counter(train[a]["class"] for a in tr_copies) != Counter({"ENDO": 19, "CONTROL": 7}):
        raise Hold("Wrong 26 train-copy diagnosis roster")
    if Counter(test[b]["group"] for b in te_copies) != Counter({"ENDO": 19, "NUPP": 7}):
        raise Hold("Wrong 26 test-copy diagnosis roster")
    vectors = {}
    for g in train:
        vectors[g] = genes_from_vector(b28, "source_CEL_MAS5_cache", g, probes)
    for g in test:
        vectors[g] = genes_from_vector(b27, "target_probe_MAS5_per_CEL", g, probes)
    if len(vectors) != 207:
        raise Hold("Wrong MAS5 vector count")
    for a, b in pairs:
        if any(not math.isclose(vectors[a][p], vectors[b][p], abs_tol=1e-9, rel_tol=1e-9) for p in probes):
            raise Hold("Identical CELs have differing MAS5 signals: " + a + " " + b)
    order_train, order_test = sorted(train), sorted(test)
    def mx(ids):
        return np.asarray([[float(np.median(np.log2([vectors[g][p] for p in mapping[gene]])))
                            for gene in GENES] for g in ids], dtype=float)
    Xtr = mx(order_train)
    Xt = mx(order_test)
    ytr = np.array([int(train[g]["class"] == "ENDO") for g in order_train], dtype=int)
    yt = np.array([int(test[g]["group"] == "ENDO") for g in order_test], dtype=int)
    if Xtr.shape != (59, 5) or Xt.shape != (148, 5):
        raise Hold("Wrong frozen expression matrix shape")
    log("[INPUT] 59 train; 148 test; 26 copy pairs; 5 genes; 6 probes; pinned MAS5 cache SHA-256 PASS")
    return train, test, order_train, order_test, Xtr, Xt, ytr, yt, pairs

def make_rf():
    from sklearn.ensemble import RandomForestClassifier
    return RandomForestClassifier(n_estimators=N_ESTIMATORS, max_features="sqrt",
                                  random_state=SEED_RF, n_jobs=1, bootstrap=True,
                                  class_weight=None)

def auc(y, scores):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, scores))

def number_equal(a, b, tol=1e-11):
    return math.isclose(float(a), float(b), rel_tol=0, abs_tol=tol)

def compare_reference(b30, individuals, placebo, summary):
    """Compare individual scores plus all 199 archived placebo records; no BA-only gates."""
    ref_ind = b30.table("01_26_own_copy_out_scores.tsv")
    ref_pl = b30.table("02_199_same_class_other_drop_placebos.tsv")
    ref_sm = b30.table("03_own_copy_ablation_summary.tsv")
    if len(ref_ind) != 26 or len(ref_pl) != 199 or len(ref_sm) != 3:
        raise Hold("v4.30 reference row counts differ")
    new_ind = {r["test_gsm"]: r for r in individuals}
    ref_ind = {r["test_gsm"]: r for r in ref_ind}
    if set(new_ind) != set(ref_ind):
        raise Hold("26 test GSM mismatch to archive")
    metrics = {"individual_scores_checked": 0, "placebo_auc_values_checked": 0,
               "summary_auc_values_checked": 0, "max_abs_score_diff": 0.,
               "max_abs_placebo_auc_diff": 0., "max_abs_summary_metric_diff": 0.}
    for g in sorted(new_ind):
        a, b = new_ind[g], ref_ind[g]
        for key in ("train_copy_gsm", "group"):
            if str(a[key]) != str(b[key]):
                raise Hold("Mismatch archived duplicate identity: " + g + " / " + key)
        for key in ("base_score", "own_copy_omitted_score", "placebo_other_drop_mean_score"):
            difference = abs(float(a[key]) - float(b[key]))
            metrics["max_abs_score_diff"] = max(metrics["max_abs_score_diff"], difference)
            metrics["individual_scores_checked"] += 1
    for r, s in zip(placebo, ref_pl):
        if str(r["replicate"]) != str(s["replicate"]):
            raise Hold("Placebo replicate ordering changed")
        for k in ("duplicate_auc", "endo_nupp_hybrid_auc", "all_148_hybrid_auc"):
            difference = abs(float(r[k]) - float(s[k]))
            metrics["max_abs_placebo_auc_diff"] = max(metrics["max_abs_placebo_auc_diff"], difference)
            metrics["placebo_auc_values_checked"] += 1
    old_summary = {r["subset"]: r for r in ref_sm}
    if set(old_summary) != {r["subset"] for r in summary}:
        raise Hold("Historical three-subset summary mismatch")
    for r in summary:
        old = old_summary[r["subset"]]
        for k in ("base_auc_all_59_train", "own_copy_removed_auc", "same_class_other_drop_placebo_mean",
                  "placebo_q025", "placebo_q975", "empirical_tail_add_one_descriptive"):
            difference = abs(float(r[k]) - float(old[k]))
            metrics["max_abs_summary_metric_diff"] = max(metrics["max_abs_summary_metric_diff"], difference)
            metrics["summary_auc_values_checked"] += 1
    metrics["exact_scores_to_1e-10"] = metrics["max_abs_score_diff"] <= 1e-10
    metrics["exact_placebos_to_1e-10"] = metrics["max_abs_placebo_auc_diff"] <= 1e-10
    metrics["exact_summary_to_1e-10"] = metrics["max_abs_summary_metric_diff"] <= 1e-10
    return metrics

def compute(bundles, out):
    import numpy as np
    import sklearn
    train, test, ntr, nte, Xtr, Xt, ytr, yt, pairs = initialize(bundles)
    ti = {g: i for i, g in enumerate(ntr)}
    ei = {g: i for i, g in enumerate(nte)}
    src = sorted(g for g in ntr if train[g]["series"] == "GSE6364")
    rf = make_rf()
    rf.fit(Xtr, ytr)
    base = rf.predict_proba(Xt)[:, list(rf.classes_).index(1)]
    for g, score in zip(nte, base):
        if not number_equal(score, test[g]["score"], 2e-12):
            raise Hold("Full-training RF baseline differs from v4.28: " + g)
    log("[BASELINE] PASS_148_FROZEN_HUANG_RECONSTRUCTION_SCORES")
    excluded = {}
    for i, g in enumerate(src, 1):
        mask = [j for j, name in enumerate(ntr) if name != g]
        model = make_rf()
        model.fit(Xtr[mask], ytr[mask])
        excluded[g] = model.predict_proba(Xt)[:, list(model.classes_).index(1)]
        if i % 10 == 0 or i == len(src):
            log("[SINGLE_OMISSION_MODEL]", str(i) + "/37")
    idx = np.array([ei[b] for a, b in pairs], dtype=int)
    actual_y = yt[idx]
    if Counter(actual_y.tolist()) != Counter({1: 19, 0: 7}):
        raise Hold("Duplicate test class distribution changed")
    own_scores = np.array([excluded[a][ei[b]] for a, b in pairs], dtype=float)
    all_own = base.copy()
    for (a, b), value in zip(pairs, own_scores):
        all_own[ei[b]] = value
    rng = np.random.default_rng(SEED_PLACEBO)
    placebo, placebo_scores = [], np.empty((N_PLACEBO, len(pairs)))
    d_nupp = np.array([test[g]["group"] != "OTHER" for g in nte], bool)
    for rep in range(N_PLACEBO):
        for j, (a, b) in enumerate(pairs):
            candidates = sorted(g for g in src if g != a and train[g]["class"] == train[a]["class"])
            omitted = candidates[int(rng.integers(len(candidates)))]
            placebo_scores[rep, j] = excluded[omitted][ei[b]]
        combined = base.copy()
        combined[idx] = placebo_scores[rep]
        placebo.append({"replicate": rep + 1, "duplicate_auc": auc(actual_y, placebo_scores[rep]),
                        "endo_nupp_hybrid_auc": auc(yt[d_nupp], combined[d_nupp]),
                        "all_148_hybrid_auc": auc(yt, combined),
                        "n_same_class_other_omissions": 26, "n_target_own_copies_kept": 26})
    subsets = {"DUP26": np.array([g in {b for a, b in pairs} for g in nte], bool),
               "ENDO_NUPP111": d_nupp, "ALL148": np.ones(len(nte), bool)}
    column = {"DUP26": "duplicate_auc", "ENDO_NUPP111": "endo_nupp_hybrid_auc",
              "ALL148": "all_148_hybrid_auc"}
    summary = []
    for label, mask in subsets.items():
        values = np.array([r[column[label]] for r in placebo], float)
        before, after = auc(yt[mask], base[mask]), auc(yt[mask], all_own[mask])
        summary.append({"subset": label, "n": int(mask.sum()), "cases": int(yt[mask].sum()),
                        "base_auc_all_59_train": before, "own_copy_removed_auc": after,
                        "delta_auc": after - before,
                        "same_class_other_drop_placebo_mean": float(values.mean()),
                        "placebo_q025": float(np.quantile(values, .025)),
                        "placebo_q975": float(np.quantile(values, .975)),
                        "count_placebo_auc_le_own": int(np.sum(values <= after + 1e-12)),
                        "empirical_tail_add_one_descriptive": float((1 + np.sum(values <= after + 1e-12)) / (N_PLACEBO + 1))})
    individuals = []
    for j, (a, b) in enumerate(pairs):
        i = ei[b]
        individuals.append({"train_copy_gsm": a, "test_gsm": b, "group": test[b]["group"],
                            "test_date": test[b].get("date", ""), "test_phase": test[b].get("phase", ""),
                            "base_score": float(base[i]), "own_copy_omitted_score": float(own_scores[j]),
                            "signed_correct_class_margin_change": (1 if yt[i] else -1) * float(own_scores[j] - base[i]),
                            "own_removed_correct_at_0_5": int((own_scores[j] >= .5) == bool(yt[i])),
                            "placebo_other_drop_mean_score": float(placebo_scores[:, j].mean())})
    write_tsv(out / "01_26_own_copy_out_scores.tsv", individuals)
    write_tsv(out / "02_199_same_class_other_drop_placebos.tsv", placebo)
    write_tsv(out / "03_own_copy_ablation_summary.tsv", summary)
    metrics = compare_reference(bundles["v430"], individuals, placebo, summary)
    metrics["sklearn_version"] = sklearn.__version__
    metrics["numpy_version"] = np.__version__
    metrics["fit_count"] = 38  # one baseline + 37 leave-one-out models
    metrics["n_assembled_placebo_score_sets"] = 199
    metrics["n_own_copy_target_models"] = 26
    metrics["primary_auc_full"] = summary[0]["base_auc_all_59_train"]
    metrics["primary_auc_own_removed"] = summary[0]["own_copy_removed_auc"]
    metrics["archive_comparison_pass"] = all((metrics["exact_scores_to_1e-10"],
                                              metrics["exact_placebos_to_1e-10"],
                                              metrics["exact_summary_to_1e-10"]))
    metrics["interpretation"] = ("Post-hoc Huang-inspired RF reconstruction from pinned single-CEL MAS5 vectors; "
         "26 pooled predictions come from 26 separately fitted omission models. NOT reproduction of Huang et al's original model; "
         "NOT independent clinical validation; no test-outcome-informed training or parameter selection.")
    return metrics

def selftest():
    assert len(GENES) == 5 and len(set(GENES)) == 5
    assert N_ESTIMATORS == 1000 and N_PLACEBO == 199 and SEED_RF == 1234567
    assert digest(b"abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    toy = b"gsm\tscore\nGSM1\t0.45\nGSM2\t0.83\n"
    assert len(table_from_bytes(toy, "toy")) == 2
    train = {"GSM1": {"series": "GSE6364"}}
    assert isinstance(number_equal(0.5, 0.500000000001, 1e-10), bool)
    log("[SELFTEST] PASS_HASH_TABLES_LOCKED_RF_CONFIGURATION")

def run_worker(out):
    log("PAPER 8 | STAGE 9A | HUANG OWN-COPY INDEPENDENT REPLAY")
    log("[POLICY] READ ONLY HISTORICAL ZIPs | LOCKED FIVE GENES | NO CEL OR MAS5 RUNS | OUTPUTS DOWNLOADS")
    report = {"stage": "P8_STAGE9A_HUANG_OWN_COPY_REPLAY", "status": "HOLD_UNSTARTED",
              "project_inputs_read_only": True, "runs_raw_CEL": False, "runs_MAS5": False,
              "published_Huang_model_reproduced": False, "v430_source_sha256_known": SOURCE_V430_SHA,
              "expected_v429_source_sha256": EXPECTED_SOURCE_V429,
              "archive_inputs": {}, "validation": {}}
    bundles = {}
    try:
        import numpy, sklearn
        log("[ENV] Python", sys.version.split()[0], "NumPy", numpy.__version__, "sklearn", sklearn.__version__)
        for key, name in BUNDLES.items():
            p = ROOT / "results" / name
            bundles[key] = Bundle(p)
            report["archive_inputs"][key] = {"name": name, "sha256": file_digest(p)}
            log("[ARCHIVE]", key, "CRC PASS", "size=", p.stat().st_size)
        metrics = compute(bundles, out)
        report["validation"] = metrics
        report["status"] = ("PASS_INDEPENDENT_V430_HUANG_OWN_COPY_REPLAY" if metrics["archive_comparison_pass"]
                            else "HOLD_V430_SCORE_OR_PLACEBO_DISAGREEMENT")
        log("[RESULT] full AUC", round(metrics["primary_auc_full"], 6),
            "own-copy-out pooled AUC", round(metrics["primary_auc_own_removed"], 6))
        log("[SCORE_COMPARE]", metrics["individual_scores_checked"], "values max_abs=", metrics["max_abs_score_diff"])
        log("[PLACEBO_COMPARE]", metrics["placebo_auc_values_checked"], "values max_abs=", metrics["max_abs_placebo_auc_diff"])
    except Exception as e:
        report["status"] = "HOLD_STAGE9A_PRECONDITION_OR_NUMERICAL_ERROR"
        report["error"] = str(e)
        report["traceback_tail"] = traceback.format_exc()[-2500:]
        log("[HOLD]", str(e))
    finally:
        for b in bundles.values():
            b.close()
        save_json(out / "stage9a_replay.json", report)
        (out / "README_STATUS.txt").write_text(
            "Paper 8 Stage 9A: " + report["status"] + "\n"
            "This replays the frozen five-gene Huang-inspired model, not the original Huang et al. model.\n"
            "26 target-specific scores pooled across separate fitted forests are not external validation of a single classifier.\n",
            encoding="utf-8")
        z = out / "Paper8_stage9A_HUANG_REPLAY_PRIVATE.zip"
        with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as archive:
            for p in sorted(out.iterdir()):
                if p.is_file() and p != z:
                    archive.write(p, arcname=p.name)
        with zipfile.ZipFile(z) as archive:
            if archive.testzip() is not None:
                raise RuntimeError("Output ZIP CRC failed")
        log("[STATUS]", report["status"])
        log("[PRIVATE_ZIP]", z)
    return report["status"].startswith("PASS_")

def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        selftest()
        return
    if len(sys.argv) == 3 and sys.argv[1] == "--worker":
        out = Path(sys.argv[2]).resolve()
        if WORK.resolve() not in out.parents:
            raise SystemExit("HOLD_OUTPUT_OUTSIDE_DOWNLOADS")
        selftest()
        sys.exit(0 if run_worker(out) else 2)
    if len(sys.argv) != 1:
        raise SystemExit("Run without arguments or with --selftest")
    if not ROOT.is_dir():
        raise SystemExit("HOLD_PROJECT_ROOT_MISSING " + str(ROOT))
    py = ROOT / ".envs/p8_fogg_classifier/bin/python3"
    if not py.is_file():
        raise SystemExit("HOLD_LOCKED_CLASSIFIER_PYTHON_MISSING " + str(py))
    out = WORK / ("Run_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    out.mkdir(parents=True, exist_ok=False)
    print("[WORKER]", py, flush=True)
    outcome = subprocess.run([str(py), str(Path(__file__).resolve()), "--worker", str(out)], check=False)
    if outcome.returncode:
        raise SystemExit(outcome.returncode)

if __name__ == "__main__":
    main()
