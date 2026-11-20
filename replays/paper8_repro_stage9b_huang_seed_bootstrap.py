#!/usr/bin/env python3
"""Paper 8 Stage 9B: independent 5-seed Huang-inspired RF sensitivity replay
and diagnostic 5000-resample class-stratified conditional bootstrap.

Source: pinned v4.21/27/28/29/30 result ZIPs; no raw CEL, no MAS5, no new
feature selection. Reuses the SHA256-pinned Stage9A *data loading* functions;
forest fits and bootstrap are freshly executed here. Compares against v4.31
five-seed values reported in the submitted manuscript's accompanying evidence.

Important: No historical bootstrap RNG/source is available here, so our 5000
bootstrap interval is an independently seeded conditional sensitivity check,
NOT a claim of exact historical v4.31 bootstrap-sequence reproduction.

Run: python3 ~/Downloads/paper8_repro_stage9b_huang_seed_bootstrap.py
Selftest: python3 ~/Downloads/paper8_repro_stage9b_huang_seed_bootstrap.py --selftest
Outputs only ~/Downloads/Paper8_repro_stage9B_huang/Run_<timestamp>/
"""
from __future__ import annotations
import csv
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
import traceback
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
DL = Path.home() / 'Downloads'
OUT_ROOT = DL / 'Paper8_repro_stage9B_huang'
SOURCE9A = DL / 'paper8_repro_stage9_huang_owncopy_replay.py'
SOURCE9A_SHA = '29ec8dff8472b9e67990381204220d1c6f712c0c3b26d46600a9e0778654c392'
SEEDS = [1234567,112358,314159,271828,161803]
# Source: evidence/v431_five_seed_AUC.tsv in Paper8_Overleaf_Main_Fixed.zip
# Exact numbers are explicit historic reference values; these are never fitted.
REF = {1234567: 0.32330827067669177,
       112358: 0.30451127819548873,
       314159: 0.30827067669172936,
       271828: 0.29323308270676696,
       161803: 0.32330827067669177}
REF_BOOT_CI = (0.10526315789473684, 0.5714285714285714)
NBOOT = 5000
BOOT_SEED = 20261115  # new explicitly fixed diagnostic seed, not a historical claim
TOL = 1e-10

class Hold(Exception): pass

def log(*args): print(*args,flush=True)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for c in iter(lambda:f.read(1<<20),b''):h.update(c)
    return h.hexdigest()

def save_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')

def save_table(path,rows):
    if not rows:raise Hold('Empty output table: '+str(path))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader();w.writerows(rows)

def load_stage9a():
    if not SOURCE9A.is_file() or SOURCE9A.is_symlink():
        raise Hold('Missing pinned Stage9A loader in Downloads: '+str(SOURCE9A))
    if sha(SOURCE9A) != SOURCE9A_SHA:
        raise Hold('Stage9A SHA mismatch: refuse unreviewed loader change')
    spec=importlib.util.spec_from_file_location('p8_stage9a_pinned',str(SOURCE9A))
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def auc(y,s):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y,s))

def worker(out):
    import numpy as np
    import sklearn
    from sklearn.ensemble import RandomForestClassifier
    log('PAPER 8 | STAGE 9B | FIVE-SEED HUANG + CONDITIONAL BOOTSTRAP')
    log('[POLICY] READ ONLY project | 5 frozen RF seeds | no MAS5 / CEL | no tuning')
    report={'stage':'P8_STAGE9B','status':'HOLD_UNSTARTED',
            'input_read_only':True,'historical_bootstrap_rng_recovered':False,
            'not_original_published_Huang_model':True,'n_boot':NBOOT,
            'diagnostic_boot_seed':BOOT_SEED,'seeds':SEEDS,'ref_auc':REF,
            'ref_boot_ci':REF_BOOT_CI,'library_versions':{'numpy':np.__version__,'sklearn':sklearn.__version__}}
    bundles={}
    try:
        a=load_stage9a()
        if a.GENES != ('ACKR1','LMNB1','MFAP4','NMU','SEMA3C') or a.N_ESTIMATORS != 1000:
            raise Hold('Source Stage9A study specification changed')
        for k,n in a.BUNDLES.items():
            bundles[k]=a.Bundle(ROOT/'results'/n)
            log('[SOURCE]',k,'CRC PASS')
        train,test,ntr,nte,Xtr,Xt,ytr,yt,pairs=a.initialize(bundles)
        ntri={g:i for i,g in enumerate(ntr)}
        ntei={g:i for i,g in enumerate(nte)}
        src=sorted(g for g in ntr if train[g]['series']=='GSE6364')
        if len(src)!=37 or len(pairs)!=26:raise Hold('Source roster changed')
        pair_tr=[a for a,b in pairs]
        pair_te=[b for a,b in pairs]
        ids=np.array([ntei[b] for b in pair_te],dtype=int)
        labels=yt[ids]
        if int(labels.sum())!=19 or int(len(labels)-labels.sum())!=7:
            raise Hold('26 duplicated target class counts differ')
        rows=[];ind_rows=[];base_seed_scores=None
        for seed in SEEDS:
            def rf():
                return RandomForestClassifier(n_estimators=1000,max_features='sqrt',random_state=seed,
                                             n_jobs=1,bootstrap=True,class_weight=None)
            mod=rf();mod.fit(Xtr,ytr)
            baseline=mod.predict_proba(Xt)[:, list(mod.classes_).index(1)]
            if seed==SEEDS[0]:
                for g,s in zip(nte,baseline):
                    if abs(float(s)-float(test[g]['score']))>2e-12:
                        raise Hold('Base seed Stage9A historic baseline score does not match '+g)
            omitted={}
            for g in src:
                keep=[j for j,t in enumerate(ntr) if t!=g]
                if len(keep)!=58:raise Hold('Wrong leave-one-out training size')
                m=rf();m.fit(Xtr[keep],ytr[keep]);omitted[g]=m.predict_proba(Xt)[:,list(m.classes_).index(1)]
            own=np.asarray([omitted[g][ntei[b]] for g,b in pairs],float)
            other=np.asarray([np.mean([omitted[g][ntei[b]] for g in src
                   if g!=a_train and train[g]['class']==train[a_train]['class']]) for a_train,b in pairs],float)
            baseauc=auc(labels,baseline[ids]);ownauc=auc(labels,own)
            otherauc=auc(labels,other)
            target=REF[seed]
            diff=abs(ownauc-target)
            rows.append({'seed':seed,'n_models':38,'n_train_each_omission':58,'n_test_duplicate':26,
                         'base_dup26_auc':baseauc,'own_removed_dup26_auc':ownauc,
                         'other_same_class_deleted_mean_probability_auc':otherauc,
                         'reference_v431_own_auc':target,'own_auc_difference':diff,
                         'own_auc_matches_1e10':int(diff<=TOL),
                         'base_auc_equals_one':int(abs(baseauc-1.0)<=TOL)})
            for i,(a_train,b) in enumerate(pairs):
                ind_rows.append({'seed':seed,'train_copy_gsm':a_train,'test_gsm':b,
                     'truth':int(labels[i]),'baseline_score':float(baseline[ids[i]]),
                     'own_copy_omitted_score':float(own[i]),
                     'same_class_other_removed_mean_score':float(other[i])})
            if seed==SEEDS[0]:
                base_seed_scores=own.copy()
                reference={r['test_gsm']:float(r['own_copy_omitted_score'])
                           for r in bundles['v430'].table('01_26_own_copy_out_scores.tsv')}
                maxgap=max(abs(float(own[i])-reference[b]) for i,b in enumerate(pair_te))
                if maxgap>TOL:raise Hold('Base seed individual scores disagree with v430')
            log('[SEED]',seed,'base=',round(baseauc,9),'own=',round(ownauc,9),
                'v431_gap=',round(diff,12),'other_class_mean=',round(otherauc,9))
        save_table(out/'01_independent_five_seed_auc.tsv',rows)
        save_table(out/'02_individual_scores_five_seed.tsv',ind_rows)
        # Class-stratified conditional percentile bootstrap on ORIGINAL base-seed scores.
        # Resamples fixed scores, not forests or verified independent participants.
        idxpos=np.flatnonzero(labels==1);idxneg=np.flatnonzero(labels==0)
        rng=np.random.default_rng(BOOT_SEED)
        boot_scores=[]
        for rep in range(NBOOT):
            i=np.concatenate((rng.choice(idxneg,len(idxneg),replace=True),
                              rng.choice(idxpos,len(idxpos),replace=True)))
            boot_scores.append({'replicate':rep+1,'auc':auc(labels[i],base_seed_scores[i])})
        save_table(out/'03_conditional_bootstrap_5000.tsv',boot_scores)
        samples=np.asarray([z['auc'] for z in boot_scores],float)
        ci=(float(np.quantile(samples,.025)),float(np.quantile(samples,.975)))
        ci_diffs=[abs(ci[i]-REF_BOOT_CI[i]) for i in range(2)]
        report['five_seed_results']=rows
        report['bootstrap_conditional']={'n_samples':26,'n_case':19,'n_control':7,
            'resamples':NBOOT,'bootstrap_unit':'CEL sample (not independently verified patient)',
            'seed':BOOT_SEED,'interval_low':ci[0],'interval_high':ci[1],
            'historical_v431_interval':list(REF_BOOT_CI),'difference_from_historical':ci_diffs,
            'caution':'Resamples frozen predictions only; does not refit 1000-tree RF or account for repeated patients'}
        report['n_forest_fits']=len(SEEDS)*38
        auc_match=all(r['own_auc_matches_1e10'] and r['base_auc_equals_one'] for r in rows)
        report['status']='PASS_FIVE_SEED_NUMERIC_REPLAY_BOOTSTRAP_SENSITIVITY_NOT_HISTORIC_CI' if auc_match else 'HOLD_FIVE_SEED_AUC_DIFFERENCE'
        log('[BOOTSTRAP] independent conditional 95% =',ci,'historic reference=',REF_BOOT_CI)
        log('[STATUS]',report['status'])
    except Exception as e:
        report['status']='HOLD_STAGE9B_ERROR'
        report['error']=str(e)
        report['traceback_tail']=traceback.format_exc()[-2400:]
        log('[HOLD]',repr(e))
    finally:
        for b in bundles.values():b.close()
        save_json(out/'stage9b_huang.json',report)
        (out/'README_STATUS.txt').write_text('Paper 8 Stage9B: '+report['status']+'\n'
          'Independent 5-seed v4.31 forest AUC replay + independently seeded conditional 5000 bootstrap.\n'
          'Not Huang et al original fitted model; 26 scores are from separately fit omission models.\n'
          'This bootstrap does not claim to match the original historical bootstrap random sequence.\n',encoding='utf-8')
        z=out/'Paper8_stage9B_HUANG_FIVE_SEED_BOOTSTRAP_PRIVATE.zip'
        with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
            for p in sorted(out.iterdir()):
                if p.is_file() and p!=z:f.write(p,arcname=p.name)
        with zipfile.ZipFile(z) as f:
            if f.testzip():raise RuntimeError('Output zip CRC failed')
        log('[PRIVATE_ZIP]',str(z))
    return report['status'].startswith('PASS_')

def selftest():
    assert len(SEEDS)==len(REF)==5 and len(set(SEEDS))==5
    assert NBOOT==5000 and REF[1234567]==0.32330827067669177
    assert sha(Path(__file__)) and len(sha(Path(__file__)))==64
    log('[SELFTEST] PASS_HUANG_SEEDS_AND_SOURCE_HASH_GATES')

def main():
    if len(sys.argv)==2 and sys.argv[1]=='--selftest':selftest();return
    if len(sys.argv)==3 and sys.argv[1]=='--worker':
        out=Path(sys.argv[2]).resolve()
        if OUT_ROOT.resolve() not in out.parents:raise SystemExit('HOLD_OUTPUT_NOT_DOWNLOADS')
        selftest();sys.exit(0 if worker(out) else 2)
    if len(sys.argv)!=1:raise SystemExit('No arguments or --selftest only')
    if not ROOT.is_dir():raise SystemExit('HOLD_PROJECT_ROOT_NOT_FOUND')
    python=ROOT/'.envs/p8_fogg_classifier/bin/python3'
    if not python.is_file():raise SystemExit('HOLD_FROZEN_PYTHON_MISSING')
    out=OUT_ROOT/('Run_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    out.mkdir(parents=True,exist_ok=False)
    log('[PYTHON_WORKER]',str(python))
    res=subprocess.run([str(python),str(Path(__file__).resolve()),'--worker',str(out)],check=False)
    if res.returncode:raise SystemExit(res.returncode)
if __name__=='__main__':main()
