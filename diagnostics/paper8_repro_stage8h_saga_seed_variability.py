#!/usr/bin/env python3
"""Paper 8 / Stage 8H -- prespecified SAGA random-state sensitivity.

Run on the original Mac:
    python3 "$HOME/Downloads/paper8_repro_stage8h_saga_seed_variability.py"

Goal: measure variability of scores in the ALREADY RECONSTRUCTED 12 v3.71
ComBat-null models for replicate 1 when changing ONLY LogisticRegression's
random_state. No fresh ComBat, no new permutations, no hyperparameter tuning.

All historical research files are READ ONLY. Results are written exclusively
under ~/Downloads/Paper8_repro_stage8H_seed_variability/Run_<UTC timestamp>/.

FIXED SEEDS: 123, 11, 23, 37, 53, 71, in this exact order. Seed 123 is the
original reconstructed reference. The other five are sensitivity diagnostics,
not candidates for replacing it. All results are reported, not ranked/chosen.

The historical v3.71 Python driver and the three assignment CSV SHA-256 pins
remain UNVERIFIED. Even a large observed seed spread cannot prove which
historical procedure produced the archived scores.

This implementation requires the exactly pinned Stage 8F script that the user
has already run, so its validated input reconstruction can be safely reused.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import traceback
import warnings
import zipfile

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
DL = Path.home() / 'Downloads'
PREVIOUS = DL / 'paper8_repro_stage8f_author_shuffle_diagnostic.py'
PREVIOUS_SHA = 'd5bc04245cc523017eefdb9d8527cdc642e97884aaa96d92de2a034188a7c4f2'
P8PY = ROOT / '.envs/p8_null/bin/python3'
V371 = ROOT / 'results/p8_v371_full_gene_combat_null_RESUMABLE'
OUTROOT = DL / 'Paper8_repro_stage8H_seed_variability'
SEEDS = (123, 11, 23, 37, 53, 71)
FAMILIES = (('FULL179','WITHIN_BATCH'),('FULL179','WITHIN_BATCH_PHASE'),('DEDUP153','WITHIN_BATCH'))
MODES = ('NO_OUTCOME','OUTCOME')
CONFIGS = {'ARCHIVED':(0.6552218193440044,0.0675608348364084),
           'PRIOR_FIXED_SENSITIVITY':(3.0,0.8)}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                                    allow_nan=False) + '\n', encoding='utf-8')


def write_rows(path, records):
    if not records:
        raise RuntimeError('HOLD_EMPTY_ROWS_'+str(path))
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)


def load_pinned_helpers():
    if not PREVIOUS.is_file() or sha(PREVIOUS) != PREVIOUS_SHA:
        raise RuntimeError('HOLD_STAGE8F_HELPER_MISSING_OR_HASH_CHANGED; expected '+str(PREVIOUS))
    return runpy.run_path(str(PREVIOUS))


def confusion(y, scores):
    import numpy as np
    y = np.asarray(y, dtype=np.int64)
    yp = np.asarray(scores, dtype=np.float64) >= 0
    return (int(np.sum((y==1)&yp)),int(np.sum((y==1)&(~yp))),
            int(np.sum((y==0)&(~yp))),int(np.sum((y==0)&yp)))


def balanced_accuracy(cm):
    tp,fn,tn,fp = cm
    if tp+fn == 0 or tn+fp == 0:
        raise RuntimeError('HOLD_SINGLE_CLASS_TEST')
    return .5 * (tp/(tp+fn) + tn/(tn+fp))


def worker(out):
    import numpy as np
    import sklearn
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.linear_model import LogisticRegression
    from sklearn.exceptions import ConvergenceWarning

    if (sys.version.split()[0], np.__version__, sklearn.__version__) != ('3.12.15','2.5.4','1.9.1'):
        raise RuntimeError('HOLD_NOT_FROZEN_P8_NULL_ENV')
    lib = load_pinned_helpers()
    d, d_info, e, e_info, scores_file = lib['select_locked_runs']()
    lib['validate_sources'](d, d_info, e, e_info, scores_file)
    read_mat=lib['read_mat']
    csv_rows=lib['csv_rows']

    prev = csv_rows(scores_file)
    by_source={}
    for r in prev:
        k=(r['layout'],r['family'],r['preprocessing'],r['configuration'])
        by_source.setdefault(k,{})[r['gsm']]=r
    if len(by_source)!=12 or any(len(v)!=37 for v in by_source.values()):
        raise RuntimeError('HOLD_STAGE8E_ARCHIVE_SCORES_NOT_444')
    roster=csv_rows(V371/'verified_216_metadata.csv')
    by_id={r['gsm']:r for r in roster}
    train=sorted(g for g,v in by_id.items() if v['split']=='TRAIN179')
    test=sorted(g for g,v in by_id.items() if v['split']=='TEST37')
    dropped=set((V371/'removed_train26.txt').read_text().splitlines())
    panel=(V371/'core270_genes.txt').read_text().splitlines()
    if len(roster)!=216 or len(by_id)!=216 or len(train)!=179 or len(test)!=37 or len(dropped)!=26 or not dropped<=set(train) or len(panel)!=270:
        raise RuntimeError('HOLD_FROZEN_ROSTERS_OR_PANEL')
    y_true=np.asarray([int(by_id[g]['truth']) for g in test],dtype=np.int64)
    if set(y_true.tolist())!={0,1}:
        raise RuntimeError('HOLD_TEST_LABELS')

    all_seed_rows=[]
    per_sample=[]
    summaries=[]
    for layout,family in FAMILIES:
        selected=train if layout=='FULL179' else sorted(set(train)-dropped)
        assign=csv_rows(V371/('assign_'+layout+'_'+family+'.csv'))
        rep1={r['gsm']:int(r['permuted_label']) for r in assign if r['replicate']=='1'}
        if len(rep1)!=len(selected) or set(rep1)!=set(selected):
            raise RuntimeError('HOLD_ASSIGNMENT_REP1_'+layout+'_'+family)
        yy=np.asarray([rep1[g] for g in selected],dtype=np.int64)
        if set(yy.tolist())!={0,1}:
            raise RuntimeError('HOLD_ONLY_ONE_TRAINING_CLASS')
        for mode in MODES:
            matrix=(V371/('baseline_'+layout+'_NO_OUTCOME.tsv.gz') if mode=='NO_OUTCOME'
                    else d/'R_combat_null'/(layout+'__'+family)/'rep_001.tsv.gz')
            if not matrix.is_file():
                raise RuntimeError('HOLD_MISSING_EXISTING_MATRIX_'+str(matrix))
            gsm,genes,data=read_mat(matrix)
            if genes!=panel:
                raise RuntimeError('HOLD_MATRIX_GENE_PANEL_NOT_LOCKED')
            ii={g:i for i,g in enumerate(gsm)}
            if len(ii)!=len(gsm) or not set(selected+test)<=set(ii):
                raise RuntimeError('HOLD_MATRIX_SAMPLE_ROSTER')
            XX=data[[ii[g] for g in selected]]
            ZZ=data[[ii[g] for g in test]]
            if XX.shape!=(len(selected),270) or ZZ.shape!=(37,270):
                raise RuntimeError('HOLD_SHAPE')
            archive_rows=csv_rows(V371/('null_'+layout+'_'+family+'_'+mode+'.csv'))
            ar={r['config']:r for r in archive_rows if r['replicate']=='1'}
            if set(ar)!=set(CONFIGS):
                raise RuntimeError('HOLD_ARCHIVED_REP1_NOT_TWO_CONFIGS')
            for config,(C,l1) in CONFIGS.items():
                historic=np.asarray(json.loads(ar[config]['score_by_sorted_GSM']),dtype=np.float64)
                if historic.shape!=(37,) or not np.isfinite(historic).all():
                    raise RuntimeError('HOLD_ARCHIVED_SCORE_FORMAT')
                prior=by_source[(layout,family,mode,config)]
                if not np.allclose(historic,[float(prior[g]['archived_score']) for g in test],rtol=0,atol=1e-12):
                    raise RuntimeError('HOLD_STAGED_REFERENCE_DIFFERENCE')
                y_com=confusion(y_true,historic)
                if y_com!=tuple(int(ar[config][col]) for col in ('TP','FN','TN','FP')):
                    raise RuntimeError('HOLD_ARCHIVED_SCORE_CONFUSION')

                vecs=[]
                seed_rows=[]
                for seed in SEEDS:
                    estimator=make_pipeline(StandardScaler(), LogisticRegression(
                        penalty='elasticnet',solver='saga',C=C,l1_ratio=l1,
                        random_state=seed,max_iter=5000,tol=1e-4,n_jobs=1))
                    with warnings.catch_warnings(record=True) as warning_list:
                        warnings.simplefilter('always',ConvergenceWarning)
                        estimator.fit(XX,yy)
                    scores=estimator.decision_function(ZZ)
                    if not np.isfinite(scores).all() or scores.shape!=(37,):
                        raise RuntimeError('HOLD_BAD_NUMERIC_SCORES')
                    fitted=estimator.named_steps['logisticregression']
                    ni=int(fitted.n_iter_[0])
                    nw=sum(isinstance(w.message,ConvergenceWarning) for w in warning_list)
                    if ni>=5000 or nw:
                        raise RuntimeError('HOLD_NO_CONVERGENCE_'+layout+'_'+family+'_'+mode+'_'+config+'_'+str(seed))
                    cm=confusion(y_true,scores)
                    row=dict(layout=layout,family=family,preprocessing=mode,
                             configuration=config,seed=seed,n_iter=ni,n_convergence_warnings=nw,
                             TP=cm[0],FN=cm[1],TN=cm[2],FP=cm[3],BA=balanced_accuracy(cm),
                             max_abs_vs_archived=float(np.max(np.abs(scores-historic))),
                             mean_abs_vs_archived=float(np.mean(np.abs(scores-historic))))
                    all_seed_rows.append(row)
                    seed_rows.append(row)
                    vecs.append(scores)
                primary=vecs[0]
                old=np.asarray([float(prior[g]['reconstructed_score']) for g in test])
                discrepancy=float(np.max(np.abs(primary-old)))
                if discrepancy>1e-7:
                    raise RuntimeError('HOLD_SEED123_DOES_NOT_REPRODUCE_STAGE8E_'+layout+'_'+family+'_'+mode+'_'+config+':'+str(discrepancy))
                matrix_seed=np.stack(vecs,axis=0)
                lo=matrix_seed.min(axis=0)
                hi=matrix_seed.max(axis=0)
                span=hi-lo
                sigma=matrix_seed.std(axis=0,ddof=1)
                ref_gap=np.abs(primary-historic)
                class_ambiguity=np.any((matrix_seed>=0)!=(matrix_seed[0]>=0),axis=0)
                archive_inside=(historic>=lo)&(historic<=hi)
                for k,g in enumerate(test):
                    per_sample.append(dict(layout=layout,family=family,preprocessing=mode,
                        configuration=config,gsm=g,truth=int(y_true[k]),
                        historic_score=float(historic[k]),seed123_score=float(primary[k]),
                        seed_range_min=float(lo[k]),seed_range_max=float(hi[k]),
                        seed_spread=float(span[k]),seed_sd=float(sigma[k]),
                        historic_vs_seed123_abs=float(ref_gap[k]),
                        historic_in_six_seed_range=bool(archive_inside[k]),
                        class_changed_across_sensitivity_seeds=bool(class_ambiguity[k])))
                summary=dict(layout=layout,family=family,preprocessing=mode,configuration=config,
                    n_seeds=len(SEEDS),n_test=37,seed123_BA=seed_rows[0]['BA'],
                    archived_BA=float(ar[config]['BA']),
                    archived_vs_seed123_max_abs=float(ref_gap.max()),
                    archived_vs_seed123_mean_abs=float(ref_gap.mean()),
                    seed_spread_max=float(span.max()),
                    seed_spread_mean=float(span.mean()),
                    seed_spread_median=float(np.median(span)),
                    seed_sd_median=float(np.median(sigma)),
                    historic_within_seed_range_count=int(archive_inside.sum()),
                    sample_classes_varying_with_seeds=int(class_ambiguity.sum()),
                    n_seed_BA_equal_archive=sum(abs(x['BA']-float(ar[config]['BA']))<1e-12 for x in seed_rows),
                    seed_BA_min=min(x['BA'] for x in seed_rows),
                    seed_BA_max=max(x['BA'] for x in seed_rows),
                    max_iter_used=max(x['n_iter'] for x in seed_rows),
                    exact_v371_driver_recovered=False,assignment_csv_byte_sha_verified=False,
                    diagnostic_only_not_selection=True)
                summaries.append(summary)
                print('[SEED_VARIABILITY]',layout,family,mode,config,
                    'archival_gap=',round(summary['archived_vs_seed123_max_abs'],6),
                    'max_seed_spread=',round(summary['seed_spread_max'],6),
                    'seed_class_changes=',summary['sample_classes_varying_with_seeds'],
                    'archive_inside_range=',str(summary['historic_within_seed_range_count'])+'/37',flush=True)
    if len(summaries)!=12 or len(all_seed_rows)!=72 or len(per_sample)!=444:
        raise RuntimeError('HOLD_UNEXPECTED_OUTPUT_SHAPES')
    write_rows(out/'twelve_model_seed_summary.csv',summaries)
    write_rows(out/'seventy_two_seed_fits.csv',all_seed_rows)
    write_rows(out/'four_hundred_forty_four_gsm_seed_variability_PRIVATE.csv',per_sample)
    result=dict(stage='P8_STAGE8H_FIXED_SEED_SAGA_VARIABILITY',
        status='PASS_NUMERIC_SEED_DIAGNOSTIC_NOT_HISTORICAL_REPRODUCTION',
        random_states=list(SEEDS),
        baseline_random_state=123,
        original_training_row_order='GSM sorted, exactly as Stage8E',
        only_change_across_seed_fits='LogisticRegression random_state',
        params=dict(solver='saga',tol=1e-4,max_iter=5000,n_jobs=1,
                    configurations={k:list(v) for k,v in CONFIGS.items()}),
        based_on_stage8D=d.name,based_on_stage8E=e.name,
        new_combat_runs=0,generated_new_assignments=False,
        classifier_fits=72,
        n_archival_score_vectors=12,
        n_samples_per_vector=37,
        n_models_where_archive_within_seed_range_all_37=sum(s['historic_within_seed_range_count']==37 for s in summaries),
        n_models_with_seed_induced_classification_variability=sum(s['sample_classes_varying_with_seeds']>0 for s in summaries),
        n_models_seed123_BA_equals_archive=sum(abs(s['seed123_BA']-s['archived_BA'])<1e-12 for s in summaries),
        models=summaries,
        limitation='This measures score variability across six preselected random_state values with fixed data, labels, model settings and row order. It cannot identify the historic v3.71 implementation or establish assignment identity. Do not select the seed closest to historic scores.')
    save_json(out/'stage8h_seed_variability.json',result)
    return result


def main():
    now=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out=OUTROOT/('Run_'+now)
    out.mkdir(parents=True,exist_ok=False)
    status='HOLD_STAGE8H_NOT_STARTED'
    print('PAPER 8 | STAGE 8H | PREDECLARED SAGA RANDOM-STATE SENSITIVITY',flush=True)
    print('[POLICY] NO COMBAT | NO NEW PERMUTATIONS | FIXED 12 MODELS x 6 SEEDS',flush=True)
    try:
        helpers=load_pinned_helpers()
        d,d_info,e,e_info,scores=helpers['select_locked_runs']()
        helpers['validate_sources'](d,d_info,e,e_info,scores)
        if not P8PY.is_file():
            raise RuntimeError('HOLD_LOCKED_CLASSIFIER_PYTHON_NOT_FOUND')
        print('[INPUTS]',d.name,'+',e.name,flush=True)
        env=dict(os.environ)
        env.update(PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
        proc=subprocess.run([str(P8PY),str(Path(__file__).resolve()),'--numeric-worker',str(out)],
                            cwd=str(out),env=env,text=True,capture_output=True,timeout=1200)
        (out/'numeric_worker.log').write_text(proc.stdout+'\n--- STDERR ---\n'+proc.stderr,encoding='utf-8')
        print(proc.stdout[-12000:],flush=True)
        if proc.returncode:
            raise RuntimeError('HOLD_STAGE8H_WORKER_EXIT_'+str(proc.returncode)+':'+proc.stderr[-1200:])
        result=json.loads((out/'stage8h_seed_variability.json').read_text(encoding='utf-8'))
        status=result['status']
    except Exception as ex:
        status='HOLD_STAGE8H_DIAGNOSTIC'
        print('[ERROR]',str(ex),flush=True)
        traceback.print_exc()
    finally:
        save_json(out/'stage8h_status.json',dict(status=status,stage='8H',
            fixed_seeds=list(SEEDS),initial_hypothesis='seed-induced SAGA numerical score spread',
            caution='HOLD for historical v3.71 score reproduction remains; diagnostic only'))
        archive=out/'Paper8_stage8H_SAGA_SEED_DIAGNOSTIC_PRIVATE.zip'
        with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for name in ('stage8h_status.json','stage8h_seed_variability.json',
                         'twelve_model_seed_summary.csv','seventy_two_seed_fits.csv',
                         'four_hundred_forty_four_gsm_seed_variability_PRIVATE.csv','numeric_worker.log'):
                p=out/name
                if p.is_file():
                    z.write(p,arcname=name)
        print('[STATUS]',status,flush=True)
        print('[PRIVATE_ZIP]',archive,flush=True)
        print('[LIMIT] Seed spread does not establish historic driver or assignment hash identity',flush=True)
    return 0 if status.startswith('PASS_') else 2


def selftest():
    cm=confusion([1,1,0,0],__import__('numpy').array([0.4,-0.3,-0.1,-0.2]))
    assert cm==(1,1,2,0), cm
    assert abs(balanced_accuracy(cm)-0.75)<1e-12
    assert len(SEEDS)==6 and len(set(SEEDS))==6 and SEEDS[0]==123
    print('[SELFTEST] PASS_PREDECLARED_SEEDS_AND_CONFUSION')


if __name__=='__main__':
    if len(sys.argv)==1:
        raise SystemExit(main())
    elif len(sys.argv)==2 and sys.argv[1]=='--selftest':
        selftest()
    elif len(sys.argv)==3 and sys.argv[1]=='--numeric-worker':
        raise SystemExit(0 if worker(Path(sys.argv[2])) else 2)
    else:
        raise SystemExit('USAGE: python3 paper8_repro_stage8h_saga_seed_variability.py [--selftest]')
