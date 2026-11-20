#!/usr/bin/env python3
"""Paper 8 Stage 8E: investigate Stage 8D continuous-score discrepancies.

Usage (macOS):
  python3 "$HOME/Downloads/paper8_repro_stage8e_score_diagnostics.py"

The original research project is strictly READ ONLY. Repeats only 12 already
completed Stage 8D null CLASSIFIER fits (NO new ComBat, no new permutations,
no hyperparameter search); writes reports to ~/Downloads only.

An exact historical v3.71 source replay is NOT claimed. This is a diagnostic
replay of the INDEPENDENT reconstructed pilot, with original historic scores
used solely for comparison after fitting.
"""
from __future__ import annotations
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback
import warnings
import zipfile

ROOT = Path.home()/'Desktop/paper8/endo-research'
RUNS = Path.home()/'Downloads/Paper8_repro_stage8D_numeric_pilot'
OUTROOT = Path.home()/'Downloads/Paper8_repro_stage8E_score_diagnostics'
P8PY = ROOT/'.envs/p8_null/bin/python3'
V371 = ROOT/'results/p8_v371_full_gene_combat_null_RESUMABLE'
CONFIGS = {'ARCHIVED':(0.6552218193440044,0.0675608348364084),
           'PRIOR_FIXED_SENSITIVITY':(3.0,0.8)}
FAMILIES = [('FULL179','WITHIN_BATCH'),('FULL179','WITHIN_BATCH_PHASE'),('DEDUP153','WITHIN_BATCH')]

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def read_csv(p):
    with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def dump(p,obj):
    Path(p).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')

def find_run():
    candidates=sorted((p for p in RUNS.glob('Run_*') if (p/'worker_report.json').is_file() and
                       (p/'stage8d_pilot.json').is_file()),reverse=True)
    if not candidates:raise RuntimeError('HOLD_NO_COMPLETE_STAGE8D_V3_RUN_FOUND')
    for p in candidates:
        s=json.loads((p/'stage8d_pilot.json').read_text())
        w=json.loads((p/'worker_report.json').read_text())
        if len(w.get('pilot_rows',[]))==30 and len(w.get('baseline_rows',[]))==12 and len(s.get('R_replays',[]))==3:
            return p,s,w
    raise RuntimeError('HOLD_STAGE8D_RUN_WITH_30_NULL_AND_12_BASELINES_NOT_FOUND')

def read_mat(p):
    import numpy as np
    op=gzip.open if p.name.endswith('.gz') else open
    with op(p,'rt',encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f,delimiter='\t')
        header=next(reader)
        if len(header) not in (217,191):raise RuntimeError('HOLD_WRONG_MATRIX_WIDTH_'+str(p.name))
        gsm=header[1:]
        if len(set(gsm))!=len(gsm):raise RuntimeError('HOLD_DUPLICATE_GSM_'+p.name)
        genes=[];arr=[]
        for row in reader:
            if len(row)!=len(header):raise RuntimeError('HOLD_INVALID_MATRIX_ROW_'+p.name)
            genes.append(row[0]);arr.append([float(x) for x in row[1:]])
    if len(genes)!=270 or len(set(genes))!=270:raise RuntimeError('HOLD_MATRIX_270_GENES_'+p.name)
    x=np.array(arr,dtype='float64').T
    if not np.isfinite(x).all():raise RuntimeError('HOLD_NONFINITE_MATRIX_'+p.name)
    return gsm,genes,x

def worker(indir,outdir):
    import numpy as np
    import sklearn
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.exceptions import ConvergenceWarning
    expected_env=json.loads((V371/'environment_lock.json').read_text())
    if (sys.version.split()[0],np.__version__,sklearn.__version__)!=('3.12.15','2.5.4','1.9.1'):
        raise RuntimeError('HOLD_ENVIRONMENT_NOT_PINNED_'+sys.version.split()[0]+'_'+np.__version__+'_'+sklearn.__version__)
    report0=json.loads((indir/'worker_report.json').read_text())
    roster=read_csv(V371/'verified_216_metadata.csv')
    if len(roster)!=216:raise RuntimeError('HOLD_BAD_ROSTER_216')
    expected_genes=(V371/'core270_genes.txt').read_text().splitlines()
    if len(expected_genes)!=270 or len(set(expected_genes))!=270:raise RuntimeError('HOLD_PANEL_LOCK')
    truth={r['gsm']:int(r['truth']) for r in roster}
    all_train=sorted(r['gsm'] for r in roster if r['split']=='TRAIN179')
    test=sorted(r['gsm'] for r in roster if r['split']=='TEST37')
    removed=set((V371/'removed_train26.txt').read_text().splitlines())
    if len(all_train)!=179 or len(test)!=37 or len(removed)!=26:raise RuntimeError('HOLD_ROSTER_COUNTS')
    meta_old={(r['layout'],r['family'],r['preprocessing'],r['configuration'],int(r['replicate'])):r
              for r in report0['pilot_rows'] if r['type'] in ('COMBAT_NO_OUTCOME_NULL','COMBAT_OUTCOME_NULL')}
    if len(meta_old)!=12:raise RuntimeError('HOLD_OLD_REPORT_COMBAT_COUNT_NOT_12')
    old_baseline=list(report0['baseline_rows'])
    old_rma=[r for r in report0['pilot_rows'] if r['type']=='RMA_NULL']
    if len(old_baseline)!=12 or len(old_rma)!=18:raise RuntimeError('HOLD_PREVIOUS_COUNTS')
    prior=list(old_baseline)+old_rma
    if any(not r['agrees'] for r in prior):raise RuntimeError('HOLD_PREVIOUS_NON_COMBAT_MISMATCH')
    if not (V371/'protocol_lock.json').is_file():raise RuntimeError('HOLD_PROTOCOL_MISSING')
    protocol=json.loads((V371/'protocol_lock.json').read_text())
    if protocol.get('fixed_hyperparameters',{})!={k:list(v) for k,v in CONFIGS.items()}:
        raise RuntimeError('HOLD_FROZEN_PARAMETER_DIFFERS')
    # Load the 3 x 199 assignment tables, use ONLY replicate 1: no new labels generated.
    ass={}
    archived={}
    for layout,family in FAMILIES:
        src=V371/('assign_'+layout+'_'+family+'.csv')
        entries=read_csv(src)
        rep1={r['gsm']:int(r['permuted_label']) for r in entries if int(r['replicate'])==1}
        if len(rep1)!=(179 if layout=='FULL179' else 153) or len(set(rep1))!=len(rep1):
            raise RuntimeError('HOLD_ASSIGNED_TRAINING_ROSTER_'+layout+'_'+family)
        expected_train=set(all_train) if layout=='FULL179' else set(all_train)-removed
        if set(rep1)!=expected_train:raise RuntimeError('HOLD_ASSIGNED_SAMPLE_IDS_'+layout+'_'+family)
        ass[(layout,family)]=rep1
        for mode in ('OUTCOME','NO_OUTCOME'):
            name='null_'+layout+'_'+family+'_'+mode+'.csv'
            rows=read_csv(V371/name)
            ref={r['config']:r for r in rows if int(r['replicate'])==1}
            if len(ref)!=2:raise RuntimeError('HOLD_REF_REPLICATE1_'+name)
            archived[(layout,family,mode)]=ref
    mats={}
    for layout in ('FULL179','DEDUP153'):
        p=V371/('baseline_'+layout+'_NO_OUTCOME.tsv.gz')
        gsm,genes,x=read_mat(p)
        if genes!=expected_genes:raise RuntimeError('HOLD_PANEL_ORDER_'+p.name)
        mats[(layout,'NO_OUTCOME')]=(gsm,x)
    for layout,family in FAMILIES:
        p=indir/'R_combat_null'/(layout+'__'+family)/'rep_001.tsv.gz'
        if not p.is_file():raise RuntimeError('HOLD_FRESH_R_COMBAT_MISSING_'+str(p))
        gsm,genes,x=read_mat(p)
        if genes!=expected_genes:raise RuntimeError('HOLD_PANEL_ORDER_'+p.name)
        mats[(layout,family,'OUTCOME')]=(gsm,x)
    details=[];models=[]
    for layout,family in FAMILIES:
        selected=all_train if layout=='FULL179' else sorted(set(all_train)-removed)
        labels=np.array([ass[(layout,family)][g] for g in selected],dtype=int)
        ytrue=np.array([truth[g] for g in test],dtype=int)
        for mode in ('NO_OUTCOME','OUTCOME'):
            gsm,x=(mats[(layout,'NO_OUTCOME')] if mode=='NO_OUTCOME' else mats[(layout,family,'OUTCOME')])
            index={g:i for i,g in enumerate(gsm)}
            if not set(selected).issubset(index) or not set(test).issubset(index):
                raise RuntimeError('HOLD_GSM_NOT_IN_MATRIX_'+layout+'_'+mode)
            X=x[[index[g] for g in selected]]
            Z=x[[index[g] for g in test]]
            for config,(c,l1) in CONFIGS.items():
                old=meta_old[(layout,family,mode,config,1)]
                ref=archived[(layout,family,mode)][config]
                refscores=np.array(json.loads(ref['score_by_sorted_GSM']),dtype='float64')
                if len(refscores)!=37 or not np.isfinite(refscores).all():
                    raise RuntimeError('HOLD_ARCHIVED_SCORE_NOT_37')
                estimator=make_pipeline(StandardScaler(),LogisticRegression(
                    penalty='elasticnet',solver='saga',C=c,l1_ratio=l1,
                    random_state=123,max_iter=5000,tol=1e-4,n_jobs=1))
                with warnings.catch_warnings(record=True) as obs:
                    warnings.simplefilter('always',ConvergenceWarning)
                    estimator.fit(X,labels)
                score=estimator.decision_function(Z)
                absdiff=np.abs(score-refscores)
                threshold_new=(score>=0)
                threshold_ref=(refscores>=0)
                flips=int(np.sum(threshold_new!=threshold_ref))
                yt=ytrue; yp=score>=0
                tp=int(np.sum((yt==1)&yp));fn=int(np.sum((yt==1)&~yp))
                tn=int(np.sum((yt==0)&~yp));fp=int(np.sum((yt==0)&yp))
                ba=.5*(tp/(tp+fn)+tn/(tn+fp))
                err=float(np.max(absdiff))
                if abs(err-float(old['max_score_abs_error']))>1e-7:
                    raise RuntimeError('HOLD_STAGE8E_NOT_REPRODUCING_STAGE8D_SCORE_DELTA_'+layout+'_'+family+'_'+mode+'_'+config)
                if (tp,fn,tn,fp)!=(int(old['TP']),int(old['FN']),int(old['TN']),int(old['FP'])):
                    raise RuntimeError('HOLD_STAGE8E_NOT_REPRODUCING_STAGE8D_CONFUSION')
                warn=[str(w.message) for w in obs]
                n_iter=int(estimator.named_steps['logisticregression'].n_iter_[0])
                row={
                    'layout':layout,'family':family,'preprocessing':mode,'configuration':config,
                    'max_abs_score_difference':err,'mean_abs_score_difference':float(np.mean(absdiff)),
                    'median_abs_score_difference':float(np.median(absdiff)),
                    'n_scores_gt_1e6':int(np.sum(absdiff>1e-6)),
                    'n_prediction_disagreements':flips,
                    'n_convergence_warnings':sum('converge' in w.lower() for w in warn),
                    'n_iter':n_iter,'max_iter':5000,'new_BA':ba,'archived_BA':float(ref['BA']),
                    'new_confusion':[tp,fn,tn,fp],
                    'archived_confusion':[int(ref[k]) for k in ('TP','FN','TN','FP')],
                    'previous_report_max_abs_error':float(old['max_score_abs_error']),
                    'source':'Reconstructed independent fixed-code refit; original v3.71 driver missing',
                }
                models.append(row)
                print('[DIAGNOSTIC]',layout,family,mode,config,
                      'max_abs=',round(err,7),'changed_classes=',flips,'n_iter=',n_iter,flush=True)
                for i,gsm_id in enumerate(test):
                    details.append({'layout':layout,'family':family,'preprocessing':mode,'configuration':config,
                                    'gsm':gsm_id,'truth':int(yt[i]),'archived_score':format(float(refscores[i]),'.17g'),
                                    'reconstructed_score':format(float(score[i]),'.17g'),
                                    'signed_difference':format(float(score[i]-refscores[i]),'.17g'),
                                    'absolute_difference':format(float(absdiff[i]),'.17g'),
                                    'archived_prediction':int(refscores[i]>=0),
                                    'reconstructed_prediction':int(score[i]>=0),
                                    'class_disagreement':int(threshold_new[i]!=threshold_ref[i])})
    if len(models)!=12 or len(details)!=444:raise RuntimeError('HOLD_DIAGNOSTIC_OUTPUT_COUNTS')
    with open(outdir/'per_gsm_score_differences_PRIVATE.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(details[0]));w.writeheader();w.writerows(details)
    info={'stage':'STAGE8E','status':'PASS_DIAGNOSTIC_ONLY_ARCHIVED_NULL_SCORES_UNRESOLVED',
          'not_historic_v371_source_reexecution':True,'not_a_new_null_experiment':True,
          'not_full_3582_refit_replay':True,'no_tuning_or_parameter_search':True,
          'model_parameters':{'max_iter':5000,'random_state':123,'tol':1e-4,'solver':'saga'},
          'archived_gate':json.loads((V371/'gates.json').read_text()).get('ARCHIVED_37_SCORE_IDENTITY'),
          'n_baseline_BA_reproduced':12,'n_RMA_null_BA_reproduced':18,
          'n_combat_models_diagnosed':len(models),
          'n_combat_models_with_score_mismatch':sum(m['max_abs_score_difference']>1e-6 for m in models),
          'n_combat_models_with_BA_mismatch':sum(abs(m['new_BA']-m['archived_BA'])>1e-10 for m in models),
          'models':models,
          'important_limitation':'The historic v3.71 executable and assignment-file hash pins are unresolved. Do not treat this as exact historic model reproducibility.'}
    dump(outdir/'score_diagnostics_PRIVATE.json',info)
    return info

def main():
    timestamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    outdir=OUTROOT/('Run_'+timestamp)
    outdir.mkdir(parents=True,exist_ok=False)
    rep={'stage':'STAGE8E','status':'HOLD_NOT_STARTED','notes':[]}
    try:
        indir,stage8d,worker8d=find_run()
        print('PAPER 8 | STAGE 8E | NULL SCORE DISCREPANCY DIAGNOSTICS',flush=True)
        print('[POLICY] READ ONLY RESEARCH PROJECT; NO COMBAT / NO TUNING; NEW OUTPUTS DOWNLOADS',flush=True)
        print('[SELECTED_STAGE8D_RUN]',indir.name,flush=True)
        rep['selected_stage8d']=indir.name
        if stage8d.get('status')!='HOLD_RECONSTRUCTED_NUMERIC_PILOT_DISAGREEMENT':
            raise RuntimeError('HOLD_NOT_EXPECTED_STAGE8D_STATUS')
        if worker8d.get('n_fits')!=42 or worker8d.get('n_mismatches')!=12:
            raise RuntimeError('HOLD_UNEXPECTED_PREVIOUS_FIT_COUNTS')
        archived=json.loads((V371/'gates.json').read_text())
        old=archived.get('ARCHIVED_37_SCORE_IDENTITY',{})
        print('[HISTORICAL_SCORE_GATE]',old.get('status'),old.get('detail'),flush=True)
        if old.get('status')!='HOLD_UNRESOLVED_PRIOR_DIFFERENCE':
            raise RuntimeError('HOLD_HISTORIC_GATE_CHANGED')
        for row in stage8d['R_replays']:
            fn=Path(row['output'])
            if not fn.is_file() or sha(fn)!=row['sha256']:
                raise RuntimeError('HOLD_FRESH_R_COMBAT_OUTPUT_CHANGED')
        py=P8PY
        if not py.is_file():raise RuntimeError('HOLD_LOCKED_PYTHON_ENV_MISSING')
        env=dict(os.environ)
        env.update(PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        proc=subprocess.run([str(py),str(Path(__file__).resolve()),'--worker',str(indir),str(outdir)],
                            cwd=str(outdir),capture_output=True,text=True,env=env,timeout=400)
        (outdir/'worker_stdout.txt').write_text(proc.stdout+'\n--- STDERR ---\n'+proc.stderr,encoding='utf-8')
        print(proc.stdout[-5500:],flush=True)
        rep['worker_exit']=proc.returncode
        if proc.returncode:
            raise RuntimeError('HOLD_WORKER_EXIT_'+str(proc.returncode)+' '+proc.stderr[-500:])
        diag=json.loads((outdir/'score_diagnostics_PRIVATE.json').read_text())
        rep.update(status=diag['status'],n_models=12,
                   n_score_mismatches=diag['n_combat_models_with_score_mismatch'],
                   n_BA_mismatches=diag['n_combat_models_with_BA_mismatch'],
                   max_score_gap=max(m['max_abs_score_difference'] for m in diag['models']))
    except Exception as e:
        rep['status']='HOLD_STAGE8E_DIAGNOSTIC_FAILED';rep['error']=str(e)
        traceback.print_exc()
    finally:
        dump(outdir/'stage8e_summary.json',rep)
        bundle=outdir/'Paper8_stage8E_NULL_SCORE_DIAGNOSTICS_PRIVATE.zip'
        with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
            for name in ('stage8e_summary.json','score_diagnostics_PRIVATE.json','per_gsm_score_differences_PRIVATE.csv','worker_stdout.txt'):
                path=outdir/name
                if path.is_file():z.write(path,arcname=name)
            # Include the small original Stage 8D reports, but NOT whole CEL matrices.
            if 'indir' in locals():
                for name in ('stage8d_pilot.json','worker_report.json','pilot_comparisons.csv','assignment_semantic_gate.json'):
                    path=indir/name
                    if path.is_file():z.write(path,arcname='prior_stage8D/'+name)
        print('[STATUS]',rep['status'],flush=True)
        print('[PRIVATE_KIT]',bundle,flush=True)
        print('[NOTE] Do not publish this diagnostic archive before removing local-path metadata.',flush=True)
    return 0 if rep['status'].startswith('PASS_') else 2

if __name__=='__main__':
    if len(sys.argv)==2 and sys.argv[1]=='--selftest':
        from io import StringIO
        kit=Path('/mnt/data/Paper8_stage8C_NULL_EVIDENCE_KIT_PRIVATE.zip')
        if kit.is_file():
            with zipfile.ZipFile(kit) as z:
                env=json.loads(z.read('archived_v371/environment_lock.json'))
                assert env['numpy']=='2.5.4' and env['sklearn']=='1.9.1'
                gate=json.loads(z.read('archived_v371/gates.json'))
                assert gate['ARCHIVED_37_SCORE_IDENTITY']['status']=='HOLD_UNRESOLVED_PRIOR_DIFFERENCE'
                n=0
                for layout,family in FAMILIES:
                    for mode in ('NO_OUTCOME','OUTCOME'):
                        raw=z.read('archived_v371/null_'+layout+'_'+family+'_'+mode+'.csv').decode()
                        rows=list(csv.DictReader(StringIO(raw)))
                        assert len(rows)==398
                        rep1=[r for r in rows if r['replicate']=='1']
                        assert len(rep1)==2
                        assert all(len(json.loads(r['score_by_sorted_GSM']))==37 for r in rep1)
                        n+=len(rep1)
                assert n==12
        print('[SELFTEST] PASS_12_FROZEN_NULL_SCORES_AND_PRIOR_HOLD_STATUS')
        sys.exit(0)
    if len(sys.argv)==4 and sys.argv[1]=='--worker':
        try:
            result=worker(Path(sys.argv[2]),Path(sys.argv[3]))
            sys.exit(0)
        except Exception as e:
            traceback.print_exc()
            sys.exit(3)
    sys.exit(main())
