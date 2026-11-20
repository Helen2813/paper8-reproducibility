#!/usr/bin/env python3
"""Paper 8 Stage 10B -- independent 2009 date-LODO and within-date AUC replay.

Runs from previously frozen 148 x 324 MAS5 intensities. Independently refits
five outer date-held-out models using training-only inner-date model selection.
Also independently recomputes within-date pairwise AUC from archived v4.22 OOF
scores (no new fit for that portion). Never refits on held-out outcomes, does
not rerun MAS5, select markers, or write to the research project.

Outputs solely to ~/Downloads/Paper8_repro_stage10B_date/Run_.../.
Requires the historical v4.22 ZIP; the original v4.23 ZIP is needed for exact
77-score comparison. If v4.23 is absent, records HOLD instead of guessing.
No CLI options required: python3 ~/Downloads/paper8_repro_stage10b_date_lodo_withindate.py
"""
from __future__ import print_function
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import traceback
import zipfile
from collections import Counter

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
RESULTS = ROOT / 'results'
SOURCE_V422 = RESULTS / 'p8_v422_locked_GSE51981_MAS5_nested_cv_RESUMABLE.zip'
SOURCE_V423 = RESULTS / 'p8_v423_2009_LODO_and_paper_membership_RESUMABLE.zip'
OUTBASE = Path.home() / 'Downloads' / 'Paper8_repro_stage10B_date'
CURRENT = Path(__file__).resolve()
SEED = 20261114
C_GRID = (0.01, 0.1, 1.0, 10.0)
DATE_COUNTS = {
    '2009-01-08': (3, 5),
    '2009-01-14': (10, 11),
    '2009-01-15': (17, 5),
    '2009-01-22': (8, 8),
    '2009-02-06': (8, 2),
}
EXPECTED_DATE_AUC = {
    '2009-01-08': 0.5333333333333333,
    '2009-01-14': 0.7181818181818183,
    '2009-01-15': 0.776470588235294,
    '2009-01-22': 0.734375,
    '2009-02-06': 0.625,
}
EXPECTED_DATE_C = {
    '2009-01-08': 10.0,
    '2009-01-14': 0.01,
    '2009-01-15': 1.0,
    '2009-01-22': 0.01,
    '2009-02-06': 10.0,
}
WITHIN_DATE_CHECKPOINTS = {
    'ENDO_vs_ALL': (400, 0.685),
    'ENDO_vs_OTHER': (338, 0.6449704142011834),
    'ENDO_vs_NUPP': (62, 0.6774193548387096),
    'ENDO_vs_OTHER_2009_ONLY': (290, 0.7),
}
SCORE_TOL = 1e-8
METRIC_TOL = 1e-8

class Hold(RuntimeError):
    pass

def log(*parts):
    print(*parts, flush=True)

def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()

def sha_file(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            h.update(b)
    return h.hexdigest()

def read_tsv(blob):
    with io.StringIO(blob.decode('utf-8-sig'),newline='') as s:
        reader = csv.DictReader(s,delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise Hold('MALFORMED_TSV_HEADER')
        return list(reader),list(reader.fieldnames)

def write_csv(path, rows, columns):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def archive_members(path, names, optional=False):
    if not path.is_file() or path.is_symlink():
        if optional:return None
        raise Hold('MISSING_ARCHIVE_'+str(path.name))
    with zipfile.ZipFile(path) as z:
        listed=z.namelist()
        if len(listed)!=len(set(listed)) or any(
            n.startswith('/') or '\\' in n or '..' in Path(n).parts for n in listed
        ): raise Hold('UNSAFE_OR_DUPLICATE_ARCHIVE_MEMBERS_'+path.name)
        if z.testzip() is not None:raise Hold('ARCHIVE_CRC_FAILED_'+path.name)
        if not all(n in listed for n in names):
            raise Hold('ARCHIVE_REQUIRED_MEMBERS_MISSING_'+path.name+' '+repr(sorted(set(names)-set(listed))))
        return {n:z.read(n) for n in names}

def get_v422():
    wanted=('00_manifest.json','01_date_x_three_groups.tsv',
            '05_148_x_324_single_CEL_MAS5_intensity.tsv',
            '06_oof_predictions_NOT_EXTERNAL.tsv')
    blobs=archive_members(SOURCE_V422,wanted)
    manifest=json.loads(blobs['00_manifest.json'])
    if manifest.get('status')!='PASS_V422_LOCKED_INTERNAL_ANALYSIS_WITH_HOLDS':
        raise Hold('V422_LOCKED_STATUS_MISMATCH')
    for k in wanted[1:]:
        if manifest.get('output_sha256',{}).get(k)!=sha_bytes(blobs[k]):
            raise Hold('V422_MEMBER_SHA256_MISMATCH_'+k)
    rows,cols=read_tsv(blobs[wanted[2]])
    preds,predcols=read_tsv(blobs[wanted[3]])
    daterows,_=read_tsv(blobs[wanted[1]])
    if len(rows)!=148 or len(cols)!=325 or cols[0]!='gsm' or len(set(cols[1:]))!=324:
        raise Hold('V422_148X324_MATRIX_SCHEMA')
    features={}
    for r in rows:
        g=r['gsm'].strip().upper()
        if not g or g in features:raise Hold('V422_DUPLICATE_GSM_MATRIX')
        features[g]=[float(r[k]) for k in cols[1:]]
    if any(not all(math.isfinite(v) and v>0 for v in row) for row in features.values()):
        raise Hold('V422_NONPOSITIVE_MAS5')
    # For every task, v4.22 carries one held-out prediction per GSM in mode=core.
    indexed={}
    for r in preds:
        key=(r['task'],r['mode'],r['gsm'].upper())
        if key in indexed:raise Hold('V422_DUPLICATED_OOF_ROW')
        indexed[key]=r
    if len([1 for k in indexed if k[0]=='ENDO_vs_ALL' and k[1]=='core'])!=148:
        raise Hold('V422_NO_FULL_148_CORE_ROSTER')
    bydate={r['date']:(int(r['ENDO']),int(r['OTHER'])) for r in daterows
            if r['date'].startswith('2009-')}
    if bydate!=DATE_COUNTS:
        raise Hold('FROZEN_2009_DATE_DISTRIBUTION_DIFFERENT_'+repr(bydate))
    log('[ARCHIVE] v422 CRC and internal SHA PASS')
    return {'features':features,'predictions':indexed,'archive_sha':sha_file(SOURCE_V422),
            'member_sha':{k:sha_bytes(v) for k,v in blobs.items()}}

def get_2009_dataset(data):
    import numpy as np
    roster={}
    for (task,mode,g),r in data['predictions'].items():
        if task!='ENDO_vs_OTHER_2009_ONLY' or mode!='core':continue
        if g in roster:raise Hold('DUPLICATE_2009_GSM')
        group=r['group'].upper()
        if group not in ('ENDO','OTHER') or int(r['outcome']) != int(group=='ENDO'):
            raise Hold('2009_LABEL_NOT_CONSISTENT_WITH_GROUP')
        if r['date'] not in DATE_COUNTS: raise Hold('2009_DATE_NOT_FROZEN')
        if g not in data['features']:raise Hold('2009_GSM_NOT_IN_FEATURES')
        if r.get('donor_id_verified','').strip():
            raise Hold('DONOR_ID_PRESENT_REQUIRES_GROUPED_REPLAY')
        roster[g]=r
    if len(roster)!=77:raise Hold('2009_TEST_ROSTER_NOT_77')
    gsms=sorted(roster)
    y=np.asarray([int(roster[g]['outcome']) for g in gsms],dtype=np.int8)
    dates=np.asarray([roster[g]['date'] for g in gsms])
    phases=np.asarray([roster[g]['phase'] or 'UNKNOWN' for g in gsms])
    x=np.log2(np.asarray([data['features'][g] for g in gsms],dtype=np.float64))
    if x.shape!=(77,324) or not np.isfinite(x).all() or int(sum(y))!=46:
        raise Hold('2009_MATRIX_OR_LABEL_SHAPE')
    for d,(a,b) in DATE_COUNTS.items():
        sel=(dates==d)
        if (int(sum(y[sel])),int(sum(1-y[sel])))!=(a,b):
            raise Hold('2009_CLASS_DATE_COUNTS_CHANGED_'+d)
    return gsms,x,y,dates,phases

def date_lodo(x,y,dates,phases):
    import numpy as np
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    def classifier(c):
        return make_pipeline(StandardScaler(),LogisticRegression(
            C=float(c),class_weight='balanced',solver='liblinear',
            max_iter=2000,random_state=SEED))
    obs=np.full(len(y),np.nan)
    prior=np.full(len(y),np.nan)
    phasepred=np.full(len(y),np.nan)
    rows=[]
    for heldout in sorted(DATE_COUNTS):
        train=np.flatnonzero(dates!=heldout)
        test=np.flatnonzero(dates==heldout)
        otherdates=sorted(set(dates[train]))
        if len(otherdates)!=4 or len(set(y[test]))!=2:
            raise Hold('OUTER_DATE_POSITIVITY_'+heldout)
        candidates=[]
        for c in C_GRID:
            date_aucs=[]
            for inner_heldout in otherdates:
                innertrain=train[dates[train]!=inner_heldout]
                innertest=train[dates[train]==inner_heldout]
                if len(set(y[innertrain]))!=2 or len(set(y[innertest]))!=2:
                    raise Hold('INNER_DATE_POSITIVITY_'+inner_heldout)
                fit=classifier(c)
                fit.fit(x[innertrain],y[innertrain])
                p=fit.predict_proba(x[innertest])[:,1]
                date_aucs.append(float(roc_auc_score(y[innertest],p)))
            candidates.append((float(np.mean(date_aucs)),float(c)))
        best=sorted(candidates,key=lambda z:(-z[0],z[1]))[0]
        model=classifier(best[1]);model.fit(x[train],y[train])
        obs[test]=model.predict_proba(x[test])[:,1]
        prior[test]=float(np.mean(y[train]))
        for pos in test:
            eligible=train[phases[train]==phases[pos]]
            phasepred[pos]=(sum(y[eligible])+1)/(len(eligible)+2) if len(eligible) else (sum(y[train])+1)/(len(train)+2)
        auc=float(roc_auc_score(y[test],obs[test]))
        rows.append({'date':heldout,'selected_C':best[1],'heldout_auc':auc,
                     'n_test':len(test),'n_train':len(train),'test_endo':int(sum(y[test])),
                     'test_other':int(sum(1-y[test])),'inner_date_mean_auc':best[0]})
    if not (np.isfinite(obs).all() and np.isfinite(prior).all() and np.isfinite(phasepred).all()):
        raise Hold('NAN_OOF_DATE_MODEL')
    auc=float(roc_auc_score(y,obs))
    macro=float(np.mean([r['heldout_auc'] for r in rows]))
    return obs,prior,phasepred,rows,auc,macro

def within_date(data):
    # Rank concordance uses FIXED molecular OOF scores; no new model is fitted.
    output=[];errors=[]
    for task,(nexpected,vexpected) in WITHIN_DATE_CHECKPOINTS.items():
        rows=[r for (t,m,g),r in data['predictions'].items() if t==task and m=='core']
        if not rows:raise Hold('MISSING_V422_TASK_'+task)
        case=[r for r in rows if r['group']=='ENDO']
        ctrl=[r for r in rows if r['group']!='ENDO']
        pairs=0;points=0.0;ties=0
        for a in case:
            for b in ctrl:
                if a['date']!=b['date']:continue
                x=float(a['OOF_probability']);y=float(b['OOF_probability'])
                if not math.isfinite(x) or not math.isfinite(y):raise Hold('BAD_V422_OOF_SCORE')
                pairs+=1
                if x>y:points+=1.0
                elif x==y:points+=0.5;ties+=1
        val=points/pairs if pairs else float('nan')
        if pairs!=nexpected or not math.isfinite(val) or abs(val-vexpected)>METRIC_TOL:
            errors.append('WITHIN_DATE_PUBLISHED_CHECKPOINT_'+task)
        output.append({'task':task,'pairs':pairs,'within_date_auc':val,'ties':ties,
                       'published_pairs':nexpected,'published_within_date_auc':vexpected,
                       'abs_error':abs(val-vexpected)})
        log('[WITHIN_DATE]',task,'pairs=',pairs,'rank_auc=',round(val,9),
            'match=',pairs==nexpected and abs(val-vexpected)<=METRIC_TOL)
    return output,errors

def compare_v423(gsms, probs, prior, phase, heldout, pooled, macro):
    wanted=('00_manifest.json','01_OOF_2009_by_GSM.tsv',
            '02_heldout_date_summary.tsv','03_2009_LODO_summary.json')
    raw=archive_members(SOURCE_V423,wanted,optional=True)
    if raw is None:
        return {'status':'HOLD_ARCHIVED_V423_ZIP_NOT_FOUND','n_scores_compared':0},[],[]
    m=json.loads(raw['00_manifest.json'])
    if m.get('status')!='PASS_V423_POSTHOC_DATE_LODO_EXPLORATORY':
        raise Hold('ARCHIVED_V423_STATUS_NOT_EXPECTED_'+str(m.get('status')))
    for k in wanted[1:]:
        if m.get('output_sha256',{}).get(k)!=sha_bytes(raw[k]):
            raise Hold('V423_OUTPUT_HASH_MISMATCH_'+k)
    rows,_=read_tsv(raw[wanted[1]])
    dates,_=read_tsv(raw[wanted[2]])
    record=json.loads(raw[wanted[3]])
    if len(rows)!=77 or len(dates)!=5:raise Hold('V423_ARCHIVED_ROW_COUNTS')
    pred_by_gsm={r['gsm'].upper():r for r in rows}
    if set(pred_by_gsm)!=set(gsms):raise Hold('V423_ARCHIVED_OOF_ROSTER')
    comp=[];err=[];mx=0.0
    for i,g in enumerate(gsms):
        r=pred_by_gsm[g]
        ds=[abs(float(r['lodo_score'])-float(probs[i])),
            abs(float(r['prior_score'])-float(prior[i])),
            abs(float(r['phase_only_score'])-float(phase[i]))]
        mx=max(mx,*ds)
        if any(d>SCORE_TOL for d in ds):err.append('V423_INDIVIDUAL_SCORE_'+g)
        comp.append({'gsm':g,'score_abs_error':ds[0],'prior_abs_error':ds[1],
                     'phase_abs_error':ds[2]})
    archive_dates={r['scanner_date']:r for r in dates}
    report_dates=[]
    for r in heldout:
        d=r['date'];a=archive_dates.get(d)
        if a is None:err.append('ARCHIVED_MISSING_DATE_'+d);continue
        delta=abs(float(a['auc'])-r['heldout_auc'])
        cdelta=abs(float(a['selected_C'])-r['selected_C'])
        if delta>METRIC_TOL or cdelta>METRIC_TOL:err.append('ARCHIVED_DATE_METRIC_'+d)
        report_dates.append({'date':d,'auc_abs_error':delta,'selected_C_matches':cdelta<=1e-12})
    archived_obs=record.get('observed',{})
    if abs(float(archived_obs.get('macro_date_auc',float('nan')))-macro)>METRIC_TOL:
        err.append('MACRO_DATE_AUC_V423_DIFFERENT')
    if abs(float(archived_obs.get('pooled',{}).get('auc',float('nan')))-pooled)>METRIC_TOL:
        err.append('POOLED_AUC_V423_DIFFERENT')
    return {'status':'PASS_77_NUMERIC_OOF_V423' if not err else 'HOLD_V423_NUMERIC_DIFFERENCES',
            'n_scores_compared':77,'max_abs_score_difference':mx,
            'n_mismatches':len(err),'mismatches':err[:30],
            'v423_archive_sha256':sha_file(SOURCE_V423),
            'v423_members_sha256':{k:sha_bytes(v) for k,v in raw.items()}},comp,report_dates

def process(out):
    import numpy as np
    import sklearn
    from sklearn.metrics import roc_auc_score
    v422=get_v422()
    gsms,x,y,dates,phases=get_2009_dataset(v422)
    log('[INPUT] 77 GSM, 46 ENDO, 31 OTHER, five dates, 324 frozen probes')
    scores,prior,phase,heldout,pooled,macro=date_lodo(x,y,dates,phases)
    log('[2009_DATE_LODO] pooled AUC=',round(pooled,9),'macro date AUC=',round(macro,9))
    discrepancies=[]
    for r in heldout:
        d=r['date'];expected=EXPECTED_DATE_AUC[d]
        if abs(r['heldout_auc']-expected)>METRIC_TOL or r['selected_C']!=EXPECTED_DATE_C[d]:
            discrepancies.append('PUBLISHED_2009_DATE_CHECKPOINT_'+d)
        log('[HELDOUT]',d,'ENDO=',r['test_endo'],'OTHER=',r['test_other'],
            'AUC=',round(r['heldout_auc'],9),'C=',r['selected_C'])
    within,wn_err=within_date(v422);discrepancies.extend(wn_err)
    historical,individual,date_compare=compare_v423(gsms,scores,prior,phase,heldout,pooled,macro)
    if historical['status']!='PASS_77_NUMERIC_OOF_V423':
        discrepancies.append(historical['status'])
    write_csv(out/'01_2009_heldout_date_fits.csv',heldout,
              ['date','n_train','n_test','test_endo','test_other','selected_C','inner_date_mean_auc','heldout_auc'])
    write_csv(out/'02_2009_GSM_replayed_scores.csv',[
        {'gsm':g,'date':d,'outcome':int(yy),'lodo_score':format(float(s),'.17g'),
         'prior_score':format(float(p),'.17g'),'phase_only_score':format(float(ph),'.17g')}
        for g,d,yy,s,p,ph in zip(gsms,dates,y,scores,prior,phase)],
        ['gsm','date','outcome','lodo_score','prior_score','phase_only_score'])
    write_csv(out/'03_within_date_pairwise_AUC.csv',within,
              ['task','pairs','within_date_auc','ties','published_pairs',
               'published_within_date_auc','abs_error'])
    if individual:
        write_csv(out/'04_V423_OOF_difference_by_GSM.csv',individual,
                  ['gsm','score_abs_error','prior_abs_error','phase_abs_error'])
    if date_compare:
        write_csv(out/'05_V423_fold_comparison.csv',date_compare,
                  ['date','auc_abs_error','selected_C_matches'])
    stat='PASS_V423_DATE_LODO_AND_WITHIN_DATE_REPLAY' if not discrepancies else 'HOLD_DATE_REPLAY_OR_MISSING_REFERENCE'
    summary={'stage':'P8_STAGE10B_2009_LODO_AND_WITHIN_DATE_RANK',
             'status':stat,'scope':'Frozen 148x324 MAS5 data and v4.22 OOF molecular scores; independent date-held-out refit for 2009; no 199-permutation replay',
             'date_2009_n':len(y),'date_2009_endo':int(sum(y)),
             'date_2009_other':int(sum(1-y)),'n_dates':5,
             'date_lodo_pooled_auc':pooled,'date_lodo_macro_auc':macro,
             'date_heldout_fits':heldout,'v423_exact_comparison':historical,
             'within_date_fixed_molecular_scores':within,'mismatches':discrepancies,
             'source_v422_archive_sha256':v422['archive_sha'],
             'source_v422_member_sha256':v422['member_sha'],
             'python':sys.version.split()[0],'numpy':np.__version__,
             'sklearn':sklearn.__version__,
             'limitations':['Exploratory sample-level, not patient-level, cross-validation',
                            'Panel might have been informed by the same GSE51981 outcomes prior to lock',
                            'Within-date statistic uses fixed OOF scores; not a within-date refitted model',
                            'Stage 10B does not rerun 199 within-date label permutations',
                            'Stage 10B does not reproduce date-only metadata OOF predictor from v4.24',
                            'No raw CEL re-extraction or MAS5 rerun']}
    (out/'stage10b_replay.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    log('[V423_ARCHIVED_OOF]',historical['status'])
    log('[STATUS]',stat)
    return summary

def selftest():
    assert len(DATE_COUNTS)==5 and sum(a+b for a,b in DATE_COUNTS.values())==77
    assert len(WITHIN_DATE_CHECKPOINTS)==4 and len(C_GRID)==4
    assert sum(a for a,b in DATE_COUNTS.values())==46
    assert sum(b for a,b in DATE_COUNTS.values())==31
    log('[SELFTEST] PASS_FIXED_2009_ROSTER_AND_PUBLISHED_CHECKPOINTS')

def create_run():
    OUTBASE.mkdir(parents=True,exist_ok=True)
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    p=OUTBASE/('Run_'+stamp)
    p.mkdir(parents=True,exist_ok=False)
    return p

def worker(out):
    try:
        report=process(out)
    except Exception as e:
        report={'status':'HOLD_STAGE10B_ERROR','error':type(e).__name__+': '+str(e),
                'traceback_tail':traceback.format_exc()[-3000:]}
        (out/'stage10b_replay.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
        log('[HOLD]',report['error'])
    (out/'README_LIMITATIONS.txt').write_text(
        'Paper 8 Stage 10B: retrospective 2009 held-out-date independent numerical replay.\n'
        'All raw data and code in Desktop/paper8 are read-only.\n'
        'This is not patient-level or external validation; no prospective claims.\n'
        'No 199-permutation reruns; no date-only model reruns; no CEL/MAS5 reruns.\n'
        'Only PASS_V423_DATE_LODO_AND_WITHIN_DATE_REPLAY proves score-level reproduction.\n',encoding='utf8')
    zpath=out/'Paper8_stage10B_DATE_REPLAY_PRIVATE.zip'
    with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(out.iterdir()):
            if p.is_file() and p!=zpath:z.write(p,p.name)
    with zipfile.ZipFile(zpath) as z:
        if z.testzip() is not None:raise Hold('OUTPUT_ZIP_CRC_FAILED')
    log('[PRIVATE_ZIP]',zpath)
    return 0 if report['status']=='PASS_V423_DATE_LODO_AND_WITHIN_DATE_REPLAY' else 2

def main():
    if len(sys.argv)==3 and sys.argv[1]=='--worker':
        return worker(Path(sys.argv[2]))
    log('PAPER 8 | STAGE 10B | INDEPENDENT 2009 DATE LODO + WITHIN-DATE SCORE REPLAY')
    log('[POLICY] PROJECT READ ONLY | NO CEL/MAS5 | NO RESELECTION | NO PERMUTATIONS | OUTPUTS DOWNLOADS')
    selftest();out=create_run()
    python=ROOT/'.envs'/'p8_v422_py'/'bin'/'python3'
    if not python.is_file():python=ROOT/'.envs'/'p8_fogg_classifier'/'bin'/'python3'
    if not python.is_file():
        (out/'stage10b_replay.json').write_text(json.dumps({'status':'HOLD_PYTHON_ENV_MISSING'},indent=2)+'\n')
        log('[STATUS] HOLD_PYTHON_ENV_MISSING');return 2
    log('[WORKER]',python)
    proc=subprocess.run([str(python),str(CURRENT),'--worker',str(out)],
                        capture_output=True,text=True)
    (out/'worker_stdout.log').write_text(proc.stdout,encoding='utf8')
    (out/'worker_stderr.log').write_text(proc.stderr,encoding='utf8')
    zpath=out/'Paper8_stage10B_DATE_REPLAY_PRIVATE.zip'
    if zpath.is_file():
        with zipfile.ZipFile(zpath,'a',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            z.write(out/'worker_stdout.log','worker_stdout.log')
            z.write(out/'worker_stderr.log','worker_stderr.log')
        with zipfile.ZipFile(zpath) as z:
            if z.testzip() is not None:raise Hold('FINAL_ZIP_CRC_FAILURE')
    print(proc.stdout,end='')
    if proc.stderr:log('[WORKER_STDERR_TAIL]',proc.stderr[-1200:])
    log('[EXIT]',proc.returncode)
    log('[OUTPUT]',out)
    return proc.returncode

if __name__=='__main__':
    sys.exit(main())
