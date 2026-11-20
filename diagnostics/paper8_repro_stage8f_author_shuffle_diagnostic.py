#!/usr/bin/env python3
"""Paper 8 / Stage 8F -- one pre-specified training-order diagnostic.

Usage on macOS:
    python3 "$HOME/Downloads/paper8_repro_stage8f_author_shuffle_diagnostic.py"

No new ComBat; no permutations generated; no parameter search; no changes to the
original research directory. Fits exactly the twelve *previously tested*
ComBat null models once, with a fixed training-row order derived from the
verified v3.50g Fogg reconstruction code:

    rng = np.random.RandomState(123)
    fit = fit.sample(frac=1, random_state=rng)

Since pandas DataFrame.sample(frac=1) uses RandomState.choice with replacement
False, RandomState.permutation(n) yields the same indices. This tests a single
historically motivated hypothesis -- it does NOT verify that the unavailable
v3.71 Python driver used that order. All results, including negative outcomes,
are reported without selection or tuning.

Inputs: frozen v3.71 ComBat null matrices and the user's prior Stage 8D/E
private diagnostics, all read-only. Outputs: Downloads only.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
import warnings
import zipfile

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
DL = Path.home() / 'Downloads'
V371 = ROOT / 'results/p8_v371_full_gene_combat_null_RESUMABLE'
P8PY = ROOT / '.envs/p8_null/bin/python3'
PREVIOUS_D = DL / 'Paper8_repro_stage8D_numeric_pilot'
PREVIOUS_E = DL / 'Paper8_repro_stage8E_score_diagnostics'
OUTROOT = DL / 'Paper8_repro_stage8F_author_order'
AUTHOR_SCRIPT = ROOT / 'paper8_v350g_core_matrisome_refit.py'
AUTHOR_HASH = '6d5284e7a8be20eb890c70fae2986e10942a34957734d913764882f4cec4c86e'
FAMILIES = (('FULL179', 'WITHIN_BATCH'), ('FULL179', 'WITHIN_BATCH_PHASE'), ('DEDUP153', 'WITHIN_BATCH'))
CONFIGS = {'ARCHIVED': (0.6552218193440044, 0.0675608348364084),
           'PRIOR_FIXED_SENSITIVITY': (3.0, 0.8)}


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def csv_rows(path):
    with open(path, 'r', newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def dump(path, val):
    Path(path).write_text(json.dumps(val, indent=2, ensure_ascii=False,
                                    allow_nan=False) + '\n', encoding='utf-8')


def latest_run(base, marker):
    for p in sorted(base.glob('Run_*'), reverse=True):
        if p.is_dir() and (p / marker).is_file():
            yield p


def select_locked_runs():
    for e in latest_run(PREVIOUS_E, 'stage8e_summary.json'):
        e_info = json.loads((e/'stage8e_summary.json').read_text())
        if e_info.get('status') != 'PASS_DIAGNOSTIC_ONLY_ARCHIVED_NULL_SCORES_UNRESOLVED':
            continue
        d = PREVIOUS_D / e_info.get('selected_stage8d', '')
        if not d.is_dir() or not (d/'stage8d_pilot.json').is_file() or not (d/'worker_report.json').is_file():
            continue
        d_info = json.loads((d/'stage8d_pilot.json').read_text())
        w = json.loads((d/'worker_report.json').read_text())
        if d_info.get('status') != 'HOLD_RECONSTRUCTED_NUMERIC_PILOT_DISAGREEMENT':
            continue
        if (w.get('n_fits'), w.get('n_mismatches'), w.get('n_RMA_null_fits'), w.get('n_ComBat_null_fits')) != (42,12,18,12):
            continue
        if len(d_info.get('R_replays', [])) != 3 or len(w.get('pilot_rows', [])) != 30:
            continue
        scores = e/'per_gsm_score_differences_PRIVATE.csv'
        if not scores.is_file():
            continue
        return d, d_info, e, e_info, scores
    raise RuntimeError('HOLD_NO_MATCHED_STAGE8D_STAGE8E_DIAGNOSTIC_RUNS')


def validate_sources(d, d_info, e, e_info, scores):
    if not ROOT.is_dir() or not V371.is_dir() or not P8PY.is_file():
        raise RuntimeError('HOLD_MISSING_PROJECT_OR_LOCKED_PYTHON')
    if not AUTHOR_SCRIPT.is_file() or digest(AUTHOR_SCRIPT) != AUTHOR_HASH:
        raise RuntimeError('HOLD_V350G_HISTORIC_AUTHOR_SCRIPT_HASH_MISMATCH')
    historic = AUTHOR_SCRIPT.read_text(encoding='utf-8')
    for exact in ('fit.sample(frac=1, random_state=rng)',
                  'rng = np.random.RandomState(123)'):
        if exact not in historic:
            raise RuntimeError('HOLD_V350G_SHUFFLE_HYPOTHESIS_NOT_VERIFIED')
    g = json.loads((V371/'gates.json').read_text())
    if g.get('ARCHIVED_37_SCORE_IDENTITY', {}).get('status') != 'HOLD_UNRESOLVED_PRIOR_DIFFERENCE':
        raise RuntimeError('HOLD_HISTORIC_SCORE_PROVENANCE_GATE_CHANGED')
    for r in d_info['R_replays']:
        f = Path(r['output'])
        if not f.is_file() or digest(f) != r['sha256']:
            raise RuntimeError('HOLD_STAGE8D_COMBAT_OUTPUT_HASH_CHANGED')
        if not f.resolve().is_relative_to(PREVIOUS_D.resolve()):
            raise RuntimeError('HOLD_STAGE8D_COMBAT_PATH_OUTSIDE_PREVIOUS_RUN')
    src_scores = csv_rows(scores)
    keys = ['layout','family','preprocessing','configuration','gsm']
    if len(src_scores) != 444 or len({tuple(x[k] for k in keys) for x in src_scores}) != 444:
        raise RuntimeError('HOLD_STAGE8E_NOT_444_UNIQUE_SAMPLE_SCORES')
    hist = V371/'protocol_lock.json'
    if not hist.is_file() or digest(hist) != 'b0b40894a9181e2e65b9afdd98f6fe658a8c0c77e362bb9f9d0d82c2f7286ac1':
        raise RuntimeError('HOLD_HISTORIC_PROTOCOL_HASH_MISMATCH')
    for n in ('verified_216_metadata.csv','removed_train26.txt','core270_genes.txt'):
        if not (V371/n).is_file():
            raise RuntimeError('HOLD_REQUIRED_INPUT_MISSING_'+n)
    return src_scores


def read_mat(path):
    import numpy as np
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        ids = header[1:]
        if len(ids) not in (216,190) or len(ids) != len(set(ids)):
            raise RuntimeError('HOLD_MATRIX_SAMPLE_ROSTER_'+path.name)
        genes, vals = [], []
        for row in reader:
            if len(row) != len(header):
                raise RuntimeError('HOLD_MATRIX_ROW_WIDTH_'+path.name)
            genes.append(row[0])
            vals.append([float(x) for x in row[1:]])
    if len(genes) != 270 or len(set(genes)) != 270:
        raise RuntimeError('HOLD_MATRIX_GENE_ROSTER_'+path.name)
    x = np.asarray(vals, dtype=np.float64).T
    if not np.isfinite(x).all():
        raise RuntimeError('HOLD_NONFINITE_'+path.name)
    return ids, genes, x


def numeric_worker(d, e, scores_file, out):
    import numpy as np
    import sklearn
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.linear_model import LogisticRegression
    from sklearn.exceptions import ConvergenceWarning

    if (sys.version.split()[0], np.__version__, sklearn.__version__) != ('3.12.15','2.5.4','1.9.1'):
        raise RuntimeError('HOLD_NUMERIC_ENV_DIFFERS_FROM_V371')
    source = csv_rows(scores_file)
    by_source = {}
    for x in source:
        key = (x['layout'], x['family'], x['preprocessing'], x['configuration'])
        by_source.setdefault(key, {})[x['gsm']] = x
    if len(by_source) != 12 or any(len(v) != 37 for v in by_source.values()):
        raise RuntimeError('HOLD_PREVIOUS_CASES_MISSING')

    meta = csv_rows(V371/'verified_216_metadata.csv')
    if len(meta) != 216 or len(set(x['gsm'] for x in meta)) != 216:
        raise RuntimeError('HOLD_FROZEN_ROSTER_INVALID')
    ids = {r['gsm']: r for r in meta}
    train = sorted(g for g in ids if ids[g]['split'] == 'TRAIN179')
    test = sorted(g for g in ids if ids[g]['split'] == 'TEST37')
    dropped = set((V371/'removed_train26.txt').read_text().splitlines())
    panel = (V371/'core270_genes.txt').read_text().splitlines()
    if len(train) != 179 or len(test) != 37 or len(dropped) != 26 or not dropped <= set(train) or len(panel) != 270:
        raise RuntimeError('HOLD_SPLIT_OR_PANEL_COUNTS')
    true_y = np.asarray([int(ids[g]['truth']) for g in test], dtype=int)

    rows = []
    details = []
    for layout, family in FAMILIES:
        use_train = train if layout == 'FULL179' else sorted(set(train) - dropped)
        assignment = csv_rows(V371/('assign_'+layout+'_'+family+'.csv'))
        rep1 = {r['gsm']: int(r['permuted_label']) for r in assignment if r['replicate'] == '1'}
        if len(rep1) != len(use_train) or set(rep1) != set(use_train):
            raise RuntimeError('HOLD_NULL_ASSIGNMENT_REP1_ROSTER_'+layout+'_'+family)
        y = np.asarray([rep1[g] for g in use_train],dtype=int)
        if len(set(y.tolist())) != 2:
            raise RuntimeError('HOLD_NULL_LABELS_SINGLE_CLASS')
        # ONE pre-specified alternative: the 123-seeded author shuffle in v3.50g.
        # The score from the previously completed sorted-order fit is read from
        # Stage 8E, never recomputed or optimized in this stage.
        order = np.random.RandomState(123).permutation(len(use_train))
        for mode in ('NO_OUTCOME', 'OUTCOME'):
            src = (V371/('baseline_'+layout+'_NO_OUTCOME.tsv.gz') if mode == 'NO_OUTCOME'
                   else d/'R_combat_null'/(layout+'__'+family)/'rep_001.tsv.gz')
            if not src.is_file():
                raise RuntimeError('HOLD_MISSING_PREVIOUS_COMBAT_MATRIX_'+src.name)
            samples, genes, data = read_mat(src)
            if genes != panel:
                raise RuntimeError('HOLD_FROZEN_GENE_ORDER_'+src.name)
            colidx = {g:i for i,g in enumerate(samples)}
            if not set(use_train).issubset(colidx) or not set(test).issubset(colidx):
                raise RuntimeError('HOLD_MATRIX_TRAIN_TEST_GSM_MISSING_'+src.name)
            X_train = data[[colidx[g] for g in use_train]]
            X_test = data[[colidx[g] for g in test]]
            for config, (C, l1) in CONFIGS.items():
                key = (layout, family, mode, config)
                s = by_source[key]
                archived = np.asarray([float(s[g]['archived_score']) for g in test])
                sorted_score = np.asarray([float(s[g]['reconstructed_score']) for g in test])
                # Independent archive-v371 checksum of the score contents, not
                # just relying on the Stage 8E copy.
                with (V371/('null_'+layout+'_'+family+'_'+mode+'.csv')).open('r',encoding='utf-8-sig',newline='') as h:
                    cand = [r for r in csv.DictReader(h) if r['replicate']=='1' and r['config']==config]
                if len(cand) != 1:
                    raise RuntimeError('HOLD_ARCHIVED_SCORE_SOURCE_ROW_COUNT')
                archive_direct = np.asarray(json.loads(cand[0]['score_by_sorted_GSM']),dtype=float)
                if len(archive_direct)!=37 or not np.allclose(archive_direct,archived,atol=1e-12,rtol=0):
                    raise RuntimeError('HOLD_STAGE8E_ARCHIVED_SCORE_COPY_MISMATCH')
                estimator = make_pipeline(StandardScaler(), LogisticRegression(
                    penalty='elasticnet', solver='saga', C=C, l1_ratio=l1,
                    random_state=123, max_iter=5000, tol=1e-4, n_jobs=1))
                with warnings.catch_warnings(record=True) as observed:
                    warnings.simplefilter('always',ConvergenceWarning)
                    estimator.fit(X_train[order],y[order])
                shuffled_score = estimator.decision_function(X_test)
                ref_max=float(np.max(np.abs(archived-sorted_score)))
                new_max=float(np.max(np.abs(archived-shuffled_score)))
                ref_mean=float(np.mean(np.abs(archived-sorted_score)))
                new_mean=float(np.mean(np.abs(archived-shuffled_score)))
                sorted_pred = sorted_score>=0
                shuffled_pred = shuffled_score>=0
                arch_pred = archived>=0
                cm=lambda p:(int(np.sum((p==1)&(true_y==1))),int(np.sum((p==0)&(true_y==1))),
                             int(np.sum((p==0)&(true_y==0))),int(np.sum((p==1)&(true_y==0))))
                ba=lambda p:float((np.mean(p[true_y==1])+np.mean(~p[true_y==0]))/2)
                counts=cm(shuffled_pred)
                archived_counts=tuple(int(cand[0][k]) for k in ('TP','FN','TN','FP'))
                if cm(arch_pred) != archived_counts:
                    raise RuntimeError('HOLD_ARCHIVE_SCORES_NOT_REPRODUCING_CONFUSION')
                warnings_seen=[str(w.message) for w in observed]
                row={
                    'layout':layout,'family':family,'mode':mode,'configuration':config,
                    'sorted_order_mean_abs_gap':ref_mean,
                    'author_shuffle_mean_abs_gap':new_mean,
                    'sorted_order_max_abs_gap':ref_max,
                    'author_shuffle_max_abs_gap':new_max,
                    'mean_abs_gap_decreased':bool(new_mean<ref_mean),
                    'max_abs_gap_decreased':bool(new_max<ref_max),
                    'sorted_order_class_disagreements':int(np.sum(sorted_pred!=arch_pred)),
                    'author_shuffle_class_disagreements':int(np.sum(shuffled_pred!=arch_pred)),
                    'sorted_order_BA':ba(sorted_pred),'author_shuffle_BA':ba(shuffled_pred),
                    'archived_BA':ba(arch_pred),
                    'archived_confusion':list(archived_counts),
                    'shuffled_confusion':list(counts),
                    'shuffled_n_iter':int(estimator.named_steps['logisticregression'].n_iter_[0]),
                    'shuffled_convergence_warnings':int(sum('converg' in x.lower() for x in warnings_seen)),
                    'sha_provenance_status':'Unresolved for v3.71 assignments',
                    'not_original_v371_driver':True,
                }
                rows.append(row)
                print('[ORDER_TEST]',layout,family,mode,config,
                      'mean_gap sorted=',round(ref_mean,7),'v350g_shuffle=',round(new_mean,7),
                      'class_changes=',row['author_shuffle_class_disagreements'],
                      'n_iter=',row['shuffled_n_iter'],flush=True)
                for i, g in enumerate(test):
                    details.append({'layout':layout,'family':family,'mode':mode,'config':config,
                                    'gsm':g,'truth':int(true_y[i]),
                                    'archived_score':format(float(archived[i]),'.17g'),
                                    'sorted_score_stage8E':format(float(sorted_score[i]),'.17g'),
                                    'v350g_shuffled_score':format(float(shuffled_score[i]),'.17g'),
                                    'archived_pred':int(arch_pred[i]),
                                    'sorted_pred':int(sorted_pred[i]),
                                    'shuffled_pred':int(shuffled_pred[i])})
    if len(rows)!=12 or len(details)!=444:
        raise RuntimeError('HOLD_DIAGNOSTIC_EXPECTED_TWELVE_MODELS_444_SCORES')
    with (out/'all_444_scores_PRIVATE.csv').open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=list(details[0]));w.writeheader();w.writerows(details)
    with (out/'twelve_order_comparisons.csv').open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    final={
        'stage':'P8_STAGE8F_FIXED_V350G_SEEDED_SHUFFLE_DIAGNOSTIC',
        'status':'PASS_DIAGNOSTIC_ONLY_HISTORIC_V371_DRIVER_UNVERIFIED',
        'source':'Verified archived v3.50g fixed-seed row shuffle vs independent Stage8E sorted-row reconstruction',
        'source_v350g_sha256':AUTHOR_HASH,
        'historical_v371_driver_found':False,
        'v371_assignment_byte_hash_pins_verified':False,
        'new_ComBat_rerun':False,
        'new_permutations_generated':False,
        'n_new_classifier_fits':12,
        'n_archived_score_vectors_compared':12,
        'n_cases_per_vector':37,
        'n_models_mean_gap_decreased':sum(r['mean_abs_gap_decreased'] for r in rows),
        'n_models_max_gap_decreased':sum(r['max_abs_gap_decreased'] for r in rows),
        'n_models_shuffled_score_max_abs_below_1e6':sum(r['author_shuffle_max_abs_gap']<=1e-6 for r in rows),
        'n_models_shuffled_classification_identical_to_archive':sum(r['author_shuffle_class_disagreements']==0 for r in rows),
        'models':rows,
        'interpretation':(
            'A single historical training-row-order hypothesis was checked, not selected as a replacement model. '
            'A numerical improvement or deterioration does not establish the original v3.71 implementation; '
            'continuous-score provenance remains unresolved unless other locked evidence identifies the historical pipeline.'),
    }
    dump(out/'stage8f_order_diagnostic.json',final)
    return final


def main():
    run_id=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out=OUTROOT/('Run_'+run_id)
    out.mkdir(parents=True,exist_ok=False)
    overall={'stage':'P8_STAGE8F','status':'HOLD_NOT_STARTED','errors':[]}
    print('PAPER 8 | STAGE 8F | VERIFIED V350G SHUFFLE HYPOTHESIS',flush=True)
    print('[POLICY] PROJECT READ ONLY | NO COMBAT | 12 FIXED FITS | NO TUNING',flush=True)
    try:
        d, d_info, e, e_info, scores=select_locked_runs()
        validate_sources(d,d_info,e,e_info,scores)
        print('[SOURCE] archived v3.50g seeded training-row shuffle verified',flush=True)
        print('[REFERENCE]',d.name,'+',e.name,flush=True)
        env=dict(os.environ)
        env.update(PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
        proc=subprocess.run([str(P8PY),str(Path(__file__).resolve()),'--numeric-worker',str(d),str(e),str(scores),str(out)],
                            cwd=str(out),env=env,capture_output=True,text=True,timeout=500)
        (out/'model_worker_log.txt').write_text(proc.stdout+'\n--- STDERR ---\n'+proc.stderr,encoding='utf-8')
        print(proc.stdout[-9000:],flush=True)
        overall.update(source_stage8D=d.name,source_stage8E=e.name,worker_exit=proc.returncode)
        if proc.returncode:
            raise RuntimeError('HOLD_NUMERIC_WORKER_EXIT_'+str(proc.returncode)+' '+proc.stderr[-550:])
        details=json.loads((out/'stage8f_order_diagnostic.json').read_text())
        overall.update(status=details['status'],n_models=details['n_new_classifier_fits'],
                       n_model_mean_gaps_decreased=details['n_models_mean_gap_decreased'],
                       n_model_shuffled_classification_matches=details['n_models_shuffled_classification_identical_to_archive'])
    except Exception as exc:
        overall['status']='HOLD_STAGE8F_DIAGNOSTIC_FAILURE'
        overall['errors'].append(str(exc))
        traceback.print_exc()
    finally:
        dump(out/'stage8f_summary.json',overall)
        archive=out/'Paper8_stage8F_AUTHOR_SHUFFLE_DIAGNOSTICS_PRIVATE.zip'
        with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for file in ('stage8f_summary.json','stage8f_order_diagnostic.json',
                         'twelve_order_comparisons.csv','all_444_scores_PRIVATE.csv','model_worker_log.txt'):
                f=out/file
                if f.is_file():z.write(f,arcname=file)
        print('[STATUS]',overall['status'],flush=True)
        print('[PRIVATE_ZIP]',archive,flush=True)
        print('[NOTE] Only one historically motivated order was tested. Original v3.71 driver still unverified.',flush=True)
    return 0 if overall['status'].startswith('PASS_') else 2


def selftest():
    import numpy as np
    import pandas as pd
    for n in (153,179):
        reference = pd.DataFrame({'idx':np.arange(n)}).sample(frac=1,random_state=np.random.RandomState(123))['idx'].to_numpy()
        actual = np.random.RandomState(123).permutation(n)
        assert np.array_equal(reference,actual)
    print('[SELFTEST] PASS_PANDAS_RANDOMSTATE123_SAMPLE_EQUALS_NUMPY_PERMUTATION')


if __name__=='__main__':
    if len(sys.argv)==1:
        raise SystemExit(main())
    elif len(sys.argv)==2 and sys.argv[1]=='--selftest':
        selftest()
    elif len(sys.argv)==6 and sys.argv[1]=='--numeric-worker':
        raise SystemExit(0 if numeric_worker(Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]),Path(sys.argv[5])) else 2)
    else:
        raise SystemExit('USAGE: python3 script.py [--selftest]')
