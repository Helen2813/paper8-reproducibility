#!/usr/bin/env python3
"""Paper 8 Stage 10C v3: independent numerical replay of v4.24 date-only controls.

Read-only research project. Writes only to Downloads. No CEL/MAS5, ComBat, molecular
training, external metadata, model tuning or new hypothesis testing. Rebuilds
metadata-only out-of-fold predictions from the locked v4.22 records and verifies
against all archived v4.24 metadata model predictions and published metrics.

Usage on Mac:
    python3 "$HOME/Downloads/paper8_repro_stage10c_date_only_numeric_replay_v3.py"

The original v4.24 script is NOT executed. This is a separate implementation of
its published frozen recipe, with explicit gates to avoid claiming independent
patient or external validation.
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
from collections import Counter, defaultdict

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
RESULTS = ROOT / 'results'
V422 = RESULTS / 'p8_v422_locked_GSE51981_MAS5_nested_cv_RESUMABLE.zip'
V424 = RESULTS / 'p8_v424_scanner_date_negative_controls_20261115T003336Z.zip'
OUTBASE = Path.home() / 'Downloads' / 'Paper8_repro_stage10C_date_only'
THIS = Path(__file__).resolve()
SEED = 20261114
SMOOTH = 5.0
TASKS = {
    'ENDO_vs_ALL': (148, 77, 71),
    'ENDO_vs_OTHER': (114, 77, 37),
    'ENDO_vs_NUPP': (111, 77, 34),
    'ENDO_vs_OTHER_2009_ONLY': (77, 46, 31),
}
MODELS = ('year','date','year_and_phase','date_and_phase',
          'molecular_archived','phase_archived')
NEW_MODELS = MODELS[:4]
TOL = 1e-10

class Hold(RuntimeError):
    pass

def say(*parts):
    print(*parts, flush=True)

def shabytes(raw):
    return hashlib.sha256(raw).hexdigest()

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def open_archive(path, needed, manifest_key, hash_scope="all"):
    """Read an archive; respect the manifest's documented hash scope."""
    if not path.is_file() or path.is_symlink():
        raise Hold('MISSING_PINNED_ZIP_'+path.name)
    with zipfile.ZipFile(path) as z:
        names=z.namelist()
        if len(names)!=len(set(names)) or any(
                n.startswith('/') or '\\' in n or '..' in Path(n).parts
                for n in names):
            raise Hold('UNSAFE_ZIP_'+path.name)
        if z.testzip() is not None:
            raise Hold('ZIP_CRC_ERROR_'+path.name)
        missing=set(needed)-set(names)
        if missing:raise Hold('MISSING_ARCHIVE_MEMBER_'+path.name+'_'+','.join(sorted(missing)))
        raw={n:z.read(n) for n in needed}
    manifest=json.loads(raw['00_manifest.json'])
    expected=manifest.get(manifest_key,{})
    for n in needed:
        if n=='00_manifest.json':continue
        if hash_scope=='tsv' and not n.endswith('.tsv'):
            continue  # v4.22 records hashes only for *.tsv, never for JSON
        wanted=expected.get(n)
        if not isinstance(wanted,str) or len(wanted)!=64 or wanted!=shabytes(raw[n]):
            raise Hold('BAD_INTERNAL_SHA_'+path.name+'_'+n)
    return raw,manifest

def _verify_v422_cross_manifest(v424_manifest, v422_sha, raw_v422):
    """Use original v4.24 source provenance to verify v4.22's unpinned JSON.

    The v4.22 generator uses `OUT.glob('*.tsv')` when building output_sha256;
    therefore its 08_nested_cv_summary.json is NOT listed. The v4.24
    generator recorded both `source_zip_sha256` and `verified_archive_sha256`
    for every input it read, including that JSON. This does not accept an
    unverified JSON or relax TSV hash checking.
    """
    if v424_manifest.get('status') != 'PASS_V424_POSTHOC_METADATA_NEGATIVE_CONTROLS':
        raise Hold('V424_PROVENANCE_STATUS_NOT_PASS')
    exp_zip = v424_manifest.get('source_zip_sha256')
    if not isinstance(exp_zip, str) or len(exp_zip)!=64 or exp_zip != v422_sha:
        raise Hold('V422_WHOLE_ZIP_SHA_MISMATCH_V424_PROVENANCE')
    per_member = v424_manifest.get('verified_archive_sha256')
    if not isinstance(per_member, dict):
        raise Hold('V424_VERIFIED_SOURCE_HASHES_MISSING')
    for n,content in raw_v422.items():
        wanted=per_member.get(n)
        if not isinstance(wanted,str) or len(wanted)!=64 or wanted!=shabytes(content):
            raise Hold('V422_SOURCE_MEMBER_SHA_MISMATCH_V424_PROVENANCE_'+n)
    return True


def verify_v422_provenance_from_v424(raw_v422):
    if not V424.is_file() or V424.is_symlink():
        raise Hold('MISSING_V424_PROVENANCE_ARCHIVE')
    with zipfile.ZipFile(V424) as z:
        names=z.namelist()
        if len(names)!=len(set(names)) or '00_manifest.json' not in names or any(
                n.startswith('/') or '\\' in n or '..' in Path(n).parts for n in names):
            raise Hold('UNSAFE_V424_PROVENANCE_ARCHIVE')
        if z.testzip() is not None:
            raise Hold('V424_PROVENANCE_CRC_ERROR')
        v424_manifest=json.loads(z.read('00_manifest.json'))
    _verify_v422_cross_manifest(v424_manifest, sha(V422), raw_v422)
    say('[V422_PROVENANCE] v4.22 TSV SHA PASS; whole ZIP + unpinned JSON SHA verified via original v4.24 manifest')
    return True


def tsv_bytes(raw):
    with io.StringIO(raw.decode('utf-8-sig'),newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        if not r.fieldnames or len(set(r.fieldnames))!=len(r.fieldnames):
            raise Hold('BAD_TSV_HEADER')
        return list(r)

def csv_bytes(raw):
    with io.StringIO(raw.decode('utf-8-sig'),newline='') as f:
        r=csv.DictReader(f)
        if not r.fieldnames or len(set(r.fieldnames))!=len(r.fieldnames):
            raise Hold('BAD_CSV_HEADER')
        return list(r)

def year_or_date(row,kind):
    date=row['date']
    phase=row['phase']
    if kind=='year':return (date[:4],)
    if kind=='date':return (date,)
    if kind=='year_and_phase':return (date[:4],phase)
    if kind=='date_and_phase':return (date,phase)
    raise Hold('UNEXPECTED_MODEL_'+kind)

def load_v422():
    wanted=('00_manifest.json','01_date_x_three_groups.tsv',
            '06_oof_predictions_NOT_EXTERNAL.tsv','07_outer_fold_manifest.tsv',
            '08_nested_cv_summary.json')
    data,manifest=open_archive(V422,wanted,'output_sha256',hash_scope='tsv')
    verify_v422_provenance_from_v424(data)
    if manifest.get('version')!='4.22' or not str(manifest.get('status','')).startswith('PASS_V422'):
        raise Hold('V422_ARCHIVE_NOT_LOCKED')
    summary=json.loads(data['08_nested_cv_summary.json'])
    if int(summary.get('frozen_seed',-1))!=SEED:
        raise Hold('V422_SEED_CHANGED')
    items=tsv_bytes(data['06_oof_predictions_NOT_EXTERNAL.tsv'])
    bytask=defaultdict(list)
    seen=set()
    for r in items:
        task=r['task']
        if r.get('mode')!='core' or task not in set(TASKS)|{'ENDO_vs_ALL_2009_ONLY'}:
            raise Hold('V422_UNEXPECTED_TASK_OR_MODE')
        k=(task,r['gsm'].upper())
        if k in seen:raise Hold('V422_DUPLICATE_GSM_TASK')
        seen.add(k)
        if r.get('donor_id_verified','').strip():
            raise Hold('V422_DONOR_GATE_CHANGED')
        if (r['outcome']=='1')!=(r['group']=='ENDO'):
            raise Hold('V422_GROUP_LABEL_MISMATCH')
        try: dt.date.fromisoformat(r['date'])
        except ValueError:raise Hold('V422_BAD_DATE')
        for field in ('OOF_probability','phase_only_probability','prior_probability'):
            v=float(r[field])
            if not math.isfinite(v) or not (0<=v<=1):raise Hold('V422_BAD_PROBABILITY')
        bytask[task].append(r)
    a={r['gsm']:r for r in bytask['ENDO_vs_ALL_2009_ONLY']}
    b={r['gsm']:r for r in bytask['ENDO_vs_OTHER_2009_ONLY']}
    if set(a)!=set(b) or any(a[g]['OOF_probability']!=b[g]['OOF_probability'] for g in a):
        raise Hold('V422_DUPLICATE_2009_COHORT_NOT_IDENTICAL')
    for task,(n,case,control) in TASKS.items():
        rr=bytask[task]
        if len(rr)!=n or sum(r['outcome']=='1' for r in rr)!=case or n-case!=control:
            raise Hold('V422_CLASS_COUNTS_'+task)
        bytask[task]=sorted(rr,key=lambda r:r['gsm'])
    # The archive also contains ENDO_vs_ALL_2009_ONLY, which is identical
    # to ENDO_vs_OTHER_2009_ONLY. Check it above, but do NOT analyse it as
    # a fifth independent task (v4.24 archived only four tasks).
    if set(bytask) != set(TASKS) | {'ENDO_vs_ALL_2009_ONLY'}:
        raise Hold('V422_UNEXPECTED_TASK_SET_AFTER_ALIAS_CHECK')
    bytask = {task: bytask[task] for task in TASKS}
    if sum(len(v) for v in bytask.values()) != 450:
        raise Hold('V422_EXPECTED_FOUR_TASKS_450_OOF_ROWS')
    say('[V422] SHA/CRC PASS | 4 tasks | 2009 duplicate alias verified and excluded | archived molecular OOF available')
    return bytask,sha(V422)

def predict_oof(rows):
    """Independent implementation of frozen training-fold-only smoothed prevalence."""
    import numpy as np
    from sklearn.model_selection import StratifiedKFold
    y=np.asarray([int(r['outcome']) for r in rows],dtype=int)
    n=len(rows)
    folds=StratifiedKFold(n_splits=5,shuffle=True,random_state=SEED)
    predictions={model:np.full(n,np.nan) for model in NEW_MODELS}
    assigned=set()
    for tr,te in folds.split(np.zeros((n,1)),y):
        prevalence=float(np.mean(y[tr]))
        cats={k:Counter(year_or_date(rows[j],k) for j in tr) for k in NEW_MODELS}
        cases={k:Counter(year_or_date(rows[j],k) for j in tr if y[j]==1) for k in NEW_MODELS}
        for j in te:
            if j in assigned:raise Hold('DUPLICATE_OOF_ASSIGNMENT')
            assigned.add(int(j))
            if abs(float(rows[j]['prior_probability'])-prevalence)>1e-12:
                raise Hold('V422_PRIOR_NOT_SAME_FOLDS_'+rows[j]['gsm'])
            for k in NEW_MODELS:
                cat=year_or_date(rows[j],k)
                predictions[k][j]=(cases[k][cat]+SMOOTH*prevalence)/(cats[k][cat]+SMOOTH)
    if len(assigned)!=n:raise Hold('INCOMPLETE_OOF_ASSIGNMENTS')
    if any(not np.isfinite(v).all() for v in predictions.values()):
        raise Hold('NONFINITE_DATE_OOF')
    predictions['molecular_archived']=np.asarray([float(r['OOF_probability']) for r in rows])
    predictions['phase_archived']=np.asarray([float(r['phase_only_probability']) for r in rows])
    return y,predictions

def pairwise_auc(y,rows,prob):
    from collections import defaultdict
    groups=defaultdict(lambda: [[],[]])
    for i,r in enumerate(rows):
        groups[r['date']][int(y[i])].append(float(prob[i]))
    tot=0;wins=0.0;dates=0
    for d,(controls,cases) in groups.items():
        if not controls or not cases:continue
        dates+=1
        for a in cases:
            for b in controls:
                tot+=1
                if a>b:wins+=1
                elif a==b:wins+=.5
    return {'pairs':tot,'mixed_dates':dates,'auc':wins/tot if tot else None}

def run_bootstrap(y, molecular, date):
    import numpy as np
    from sklearn.metrics import roc_auc_score
    rng=np.random.default_rng(SEED+424)
    n=len(y);diffs=[]
    for i in range(2500):
        ix=rng.integers(0,n,n)
        if len(set(y[ix]))!=2:continue
        diffs.append(float(roc_auc_score(y[ix],molecular[ix])-roc_auc_score(y[ix],date[ix])))
    if len(diffs)<2375:raise Hold('BOOTSTRAP_NOT_ENOUGH_VALID_DRAWS')
    return {'point':float(roc_auc_score(y,molecular)-roc_auc_score(y,date)),
            'ci95':[float(x) for x in np.quantile(diffs,[.025,.975])],
            'reps':len(diffs),'unit':'SAMPLE_NOT_PATIENT'}

def calculate(bytask):
    from sklearn.metrics import roc_auc_score,brier_score_loss,balanced_accuracy_score
    import numpy as np
    metric=[]; preds=[]; withins=[]; boot={}
    for task,rows in bytask.items():
        y,probs=predict_oof(rows)
        for m in MODELS:
            p=probs[m]
            metric.append({'task':task,'model':m,'n':len(y),'cases':int(sum(y)),
                           'controls':int(len(y)-sum(y)),
                           'auc':float(roc_auc_score(y,p)),
                           'brier':float(brier_score_loss(y,p)),
                           'balanced_accuracy':float(balanced_accuracy_score(y,(p>=0.5).astype(int)))})
        for i,r in enumerate(rows):
            preds.append({'task':task,'gsm':r['gsm'],'group':r['group'],
                          'scanner_date':r['date'],'phase':r['phase'],
                          'outcome':int(r['outcome']),
                          **{m:float(probs[m][i]) for m in MODELS}})
        within=pairwise_auc(y,rows,probs['molecular_archived'])
        withins.append({'task':task,'pooled_oof_auc':float(roc_auc_score(y,probs['molecular_archived'])),
                        'date_conditional_pair_auc':within['auc'],
                        'mixed_dates':within['mixed_dates'],
                        'within_date_case_control_pairs':within['pairs']})
        boot[task]=run_bootstrap(y,probs['molecular_archived'],probs['date'])
        say('[DATE_ONLY]',task,'date_auc=',round(next(x['auc'] for x in metric if x['task']==task and x['model']=='date'),9),
            'molecular_auc=',round(next(x['auc'] for x in metric if x['task']==task and x['model']=='molecular_archived'),9),
            'within_date=',round(within['auc'],9))
    if len(metric)!=24 or len(preds)!=450 or len(withins)!=4 or len(boot)!=4:
        raise Hold('RESULT_COUNT_GATE_metrics_%d_predictions_%d_within_%d_bootstrap_%d' %
                   (len(metric),len(preds),len(withins),len(boot)))
    return metric,preds,withins,boot

def compare_archived(metric,preds,within,boot):
    names=('00_manifest.json','01_metrics.csv','02_oof_metadata_control_predictions.csv',
           '03_date_conditional_AUC.csv','05_paired_bootstrap_exploratory.json')
    raw,manifest=open_archive(V424,names,'outputs_sha256')
    if manifest.get('status')!='PASS_V424_POSTHOC_METADATA_NEGATIVE_CONTROLS':
        raise Hold('ARCHIVED_V424_NOT_PASS')
    if int(manifest.get('baseline_rule',{}).get('seed',-1))!=SEED or float(manifest.get('baseline_rule',{}).get('shrinkage',-1))!=SMOOTH:
        raise Hold('V424_PARAMS_NOT_LOCKED')
    m0=csv_bytes(raw['01_metrics.csv']);p0=csv_bytes(raw['02_oof_metadata_control_predictions.csv'])
    w0=csv_bytes(raw['03_date_conditional_AUC.csv'])
    boot0=json.loads(raw['05_paired_bootstrap_exploratory.json'])
    errors=[];maxerr=0.0;count=0
    def match_numeric(left,right,label):
        nonlocal maxerr,count
        try: a=float(left);b=float(right)
        except (TypeError,ValueError):
            errors.append('NONNUMERIC_'+label);return
        if not math.isfinite(a) or not math.isfinite(b):errors.append('NONFINITE_'+label);return
        err=abs(a-b);maxerr=max(maxerr,err);count+=1
        if err>TOL:errors.append(label+'__diff_'+str(err))
    def index(rows,keys):
        result={}
        for r in rows:
            k=tuple(str(r[x]) for x in keys)
            if k in result:raise Hold('DUPLICATE_ARCHIVE_KEY_'+repr(k))
            result[k]=r
        return result
    for label,a,b,keys,nums,texts in [
      ('metrics',metric,m0,('task','model'),('n','cases','controls','auc','brier','balanced_accuracy'),()),
      ('scores',preds,p0,('task','gsm'),('outcome',)+MODELS,('group','scanner_date','phase')),
      ('within',within,w0,('task',),('pooled_oof_auc','date_conditional_pair_auc','mixed_dates','within_date_case_control_pairs'),()),
    ]:
        ia=index(a,keys); ib=index(b,keys)
        if ia.keys()!=ib.keys():
            errors.append('ARCHIVED_KEYS_DIFFER_'+label)
        for k in set(ia)&set(ib):
            for field in nums:
                match_numeric(ia[k][field],ib[k][field],label+'_'+str(k)+'_'+field)
            for field in texts:
                if str(ia[k][field])!=str(ib[k][field]):
                    errors.append(label+'_'+str(k)+'_TEXT_'+field)
    if set(boot)!=set(boot0):errors.append('BOOTSTRAP_TASK_SET_CHANGED')
    for task in set(boot)&set(boot0):
        match_numeric(boot[task]['point'],boot0[task]['point'],'boot_'+task+'_point')
        match_numeric(boot[task]['ci95'][0],boot0[task]['ci95'][0],'boot_'+task+'_lower')
        match_numeric(boot[task]['ci95'][1],boot0[task]['ci95'][1],'boot_'+task+'_upper')
        match_numeric(boot[task]['reps'],boot0[task]['reps'],'boot_'+task+'_reps')
    return {'status':'PASS_NUMERIC_V424' if not errors else 'HOLD_V424_NUMERIC_DIFFERENCES',
            'n_numeric_values_compared':count,'max_abs_numeric_error':maxerr,
            'n_mismatches':len(errors),'mismatch_preview':errors[:30],
            'archive_sha256':sha(V424)}

def savecsv(path,rows,fields):
    with path.open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        w.writeheader();w.writerows(rows)

def worker(path):
    report={'stage':'P8_STAGE10C_DATE_ONLY_NEGATIVE_CONTROL','status':'HOLD_NOT_RUN',
            'limitations':['Sample-level internal folds, not patient-independent or external validation',
                           'Panel and cohort were investigated previously; post hoc analysis',
                           'No raw CEL/MAS5 replay; molecular scores reused from v4.22',
                           'Stage 10C does not compare missing v4.23 per-GSM scores']}
    try:
        bytask,source_sha=load_v422()
        metric,preds,within,boot=calculate(bytask)
        comp=compare_archived(metric,preds,within,boot)
        savecsv(path/'01_replayed_v424_metrics.csv',metric,
                ['task','model','n','cases','controls','auc','brier','balanced_accuracy'])
        savecsv(path/'02_replayed_v424_predictions.csv',preds,
                ['task','gsm','group','scanner_date','phase','outcome']+list(MODELS))
        savecsv(path/'03_replayed_within_date.csv',within,
                ['task','pooled_oof_auc','date_conditional_pair_auc','mixed_dates','within_date_case_control_pairs'])
        (path/'04_replayed_sample_bootstrap.json').write_text(json.dumps(boot,indent=2)+'\n',encoding='utf8')
        report.update({'v422_v424_cross_manifest_provenance':'PASS',
                       'status':'PASS_INDEPENDENT_V424_METADATA_CONTROL_REPLAY' if comp['status']=='PASS_NUMERIC_V424' else 'HOLD_DATE_ONLY_NUMERIC_DIFFERENCES',
                       'archive_comparison':comp,'source_v422_sha256':source_sha,
                       'recomputed_n_predictions':len(preds),'recomputed_n_metric_rows':len(metric),
                       'key_results':{
                          'ENDO_vs_NUPP_date_auc':next(r['auc'] for r in metric if r['task']=='ENDO_vs_NUPP' and r['model']=='date'),
                          'ENDO_vs_NUPP_molecular_auc':next(r['auc'] for r in metric if r['task']=='ENDO_vs_NUPP' and r['model']=='molecular_archived')},
                       'note':'Independent metadata-only OOF reimplementation; not new feature selection or validation'})
    except Exception as e:
        report.update({'status':'HOLD_STAGE10C_ERROR','error':type(e).__name__+': '+str(e),
                       'traceback_tail':traceback.format_exc()[-2500:]})
    (path/'stage10c_replay.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    say('[STATUS]',report['status'])
    if 'archive_comparison' in report:
        say('[COMPARE]',report['archive_comparison'])
    if 'error' in report:say('[ERROR]',report['error'])
    zipout=path/'Paper8_stage10C_DATE_ONLY_REPLAY_PRIVATE.zip'
    with zipfile.ZipFile(zipout,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(path.iterdir()):
            if p.is_file() and p!=zipout:z.write(p,p.name)
    with zipfile.ZipFile(zipout) as z:
        if z.testzip() is not None:raise Hold('OUTPUT_CRC_ERROR')
    say('[PRIVATE_ZIP]',zipout)
    return 0 if report['status']=='PASS_INDEPENDENT_V424_METADATA_CONTROL_REPLAY' else 2

def selftest():
    assert sum(v[0] for v in TASKS.values())==450
    assert set(TASKS)=={'ENDO_vs_ALL','ENDO_vs_OTHER','ENDO_vs_NUPP','ENDO_vs_OTHER_2009_ONLY'}
    assert 'ENDO_vs_ALL_2009_ONLY' not in TASKS  # duplicate 2009 cohort must not be counted twice
    assert 6*len(TASKS)==24
    assert set(MODELS)==set(NEW_MODELS)|{'molecular_archived','phase_archived'}
    rr={'date':'2009-01-08','phase':'mse'}
    assert year_or_date(rr,'date_and_phase')==('2009-01-08','mse')
    assert year_or_date(rr,'year')==('2009',)
    # Independent check that category boundaries don't mix years or phases.
    assert year_or_date(rr,'year_and_phase')!=year_or_date({'date':'2008-01-08','phase':'mse'},'year_and_phase')
    # Artificial provenance fixture: v4.22 JSON has no own hash, but v4.24
    # binds every source member and the entire source ZIP. Reject tampering.
    contents={'00_manifest.json':b'{}', '01_date_x_three_groups.tsv':b'date\n2009',
              '08_nested_cv_summary.json':b'{"frozen_seed":20261114}'}
    fake_zip_hash='a'*64
    m={'status':'PASS_V424_POSTHOC_METADATA_NEGATIVE_CONTROLS',
       'source_zip_sha256':fake_zip_hash,
       'verified_archive_sha256':{n:shabytes(raw) for n,raw in contents.items()}}
    assert _verify_v422_cross_manifest(m,fake_zip_hash,contents)
    corrupted=dict(contents)
    corrupted['08_nested_cv_summary.json']=b'{}'
    try:
        _verify_v422_cross_manifest(m,fake_zip_hash,corrupted)
        raise AssertionError('JSON_TAMPERING_NOT_DETECTED')
    except Hold as e:
        assert str(e).startswith('V422_SOURCE_MEMBER_SHA_MISMATCH')
    try:
        _verify_v422_cross_manifest(m,'0'*64,contents)
        raise AssertionError('WHOLE_ZIP_TAMPERING_NOT_DETECTED')
    except Hold as e:
        assert 'WHOLE_ZIP_SHA_MISMATCH' in str(e)
    say('[SELFTEST] PASS_DATE_ONLY_FIXED_MODEL_AND_COHORT_RULES + CROSS_MANIFEST_TAMPER_GATES')

def main():
    say('PAPER 8 | STAGE 10C | INDEPENDENT DATE-ONLY NEGATIVE CONTROL REPLAY')
    say('[POLICY] PROJECT READ ONLY | NO RAW CEL/MAS5 OR COMBAT | NO PANEL SELECTION')
    selftest()
    OUTBASE.mkdir(parents=True,exist_ok=True)
    path=OUTBASE/('Run_'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    path.mkdir(parents=True,exist_ok=False)
    workerpy=ROOT/'.envs'/'p8_v422_py'/'bin'/'python3'
    if not workerpy.is_file():
        workerpy=ROOT/'.envs'/'p8_fogg_classifier'/'bin'/'python3'
    if not workerpy.is_file():
        say('[STATUS] HOLD_PYTHON_ENV_MISSING');return 2
    say('[WORKER]',workerpy)
    proc=subprocess.run([str(workerpy),str(THIS),'--worker',str(path)],capture_output=True,text=True)
    (path/'worker_stdout.log').write_text(proc.stdout,encoding='utf8')
    (path/'worker_stderr.log').write_text(proc.stderr,encoding='utf8')
    private_zip=path/'Paper8_stage10C_DATE_ONLY_REPLAY_PRIVATE.zip'
    if private_zip.is_file():
        with zipfile.ZipFile(private_zip,'a',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            z.write(path/'worker_stdout.log','worker_stdout.log')
            z.write(path/'worker_stderr.log','worker_stderr.log')
        with zipfile.ZipFile(private_zip) as z:
            if z.testzip() is not None:
                say('[STATUS] HOLD_FINAL_PRIVATE_ZIP_CRC');return 2
    print(proc.stdout,end='')
    if proc.stderr:say('[WORKER_STDERR_TAIL]',proc.stderr[-1600:])
    say('[EXIT]',proc.returncode)
    say('[OUTPUT]',path)
    return proc.returncode

if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--worker':
        sys.exit(worker(Path(sys.argv[2])))
    sys.exit(main())
