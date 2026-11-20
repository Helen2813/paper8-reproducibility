#!/usr/bin/env python3
"""Paper 8 / Reproducibility Stage 7.

One bounded, frozen core-matrisome model replay from archived post-ComBat
classifier-input matrices. This does NOT reproduce raw CEL -> RMA, and does
NOT estimate independent diagnostic performance. It checks the archived
AUTHOR_ARCHIVED_PARAMETERS arm against all 37 fixed-score references.

Project is READ ONLY. Outputs are exclusively under ~/Downloads/.
Usage: python3 ~/Downloads/paper8_repro_stage7_classifier.py
       python3 ~/Downloads/paper8_repro_stage7_classifier.py --selftest
"""
import base64
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
from datetime import datetime, timezone
import zlib

STAGE = 'P8_STAGE7_FROZEN_CORE_CLASSIFIER_REPLAY'
PROJECT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
DOWNLOADS = Path.home() / 'Downloads'
OUTPUT_ROOT = DOWNLOADS / 'Paper8_repro_stage7_classifier'
R_SCRIPT = 'paper8_v352d_frozen_score.py'
EXPECTED_SOURCES = {
    'paper8_v350f_classifier_inputs.py': '80b5a01d0b47ff5d0112551909110035c7131539e3676d78c4f9039e7dc50080',
    'paper8_v350g_core_matrisome_refit.py': '6d5284e7a8be20eb890c70fae2986e10942a34957734d913764882f4cec4c86e',
    'paper8_v352d_frozen_score.py': '7c479fa3f66be94fd6bbcced37500edc526d86f6003cddb7f44bf9e1f738037e',
}
EXPECTED_ROSTER_SHA = '4b786c25bfd81ed78ac3fee6d141de623e64d1cc3678b931049f0fbc0e99e97b'
EXPECTED_S3_SHA = '865eea2ef379ded51ce574096693e07d1ac99b56f48d01df61f4eaf3d7a4708d'
S2_COMPRESSED_B64 = 'eNqNmL1uk0EQRXsk3uQrdnZ3/koKhCigIOlRCIZYSpTINkh5e764CSYTfNqro5ndvb7jsX/u75ZvV4frm+Ww+3W4WfYPt9vD8nBztd+8ffPh4pO0DJmLLG25/PLu42fxXB5297fbH5vd1WH7+y9KEWWIckL1dkptrna3j1/3m+vd5nC/e3zmBHIdcv+8x932e0UpogxRjqhAVBJqNETJWaqrmY5lrOD/H/YITgrqKVh8QI6YnWKvnc9p22D1EmHW0CVM4Omss7aDtZ0MU9bUWDXqhDEnjDnhjWHCMOaDMx98smrMB2eJcGdYsCskwoLlIWgeorN6zIdgeQg2mYLlIWgeguUhmBPJ8pDCqrE8JPMhWR5S4cslS0SyRCRLRCIfvCEfvMFEeOus7WDYZKdTejpjbZ1hQdsyL4R5IcKwDk8ng4LMDVGGGW3rFAz2zMyN/uRGO4/JKfbq6Xpn9QbDJm2rFDTW2Gm9OAVf+Rj0RG1Hg22HsHodnW4wN8ZkmLKmxjBnGPNhMB9mQ9Umc2GyTMzBmjIXpjKMuTCZCzNY00TVlLmgdDZppyDLgzInlOVB2WRS5oSyPCjLg7FvCWN5MJYHYy4Y/Y4w+h1hLBPGnDDmhLFMOHPCmRPOnHD4d5I73Xac7Z5Od0/2y9iD7Z7Bdk/2u9iDjYlgAzvYmAi6wAQb2cE+nkkXmBRWj47sZMMimRfJHjnZBEj2xIlmcbTzE0CbrEY8/Qt/+f7ict2gnmWp5V7Lo5ZnLWst2yq3l7LXctRyVnJvUsvldXobtTxrWWvZatlrOeoD1teRVsulaV3qW0p9S6lvKfUtpb6l1LeUqOUs5V7fstde9l7Lo5ZnLWstWy2//GiuO3xb97j5/+AdKUdUEGrdp+a56Xfkxin3cg4dKXR+c9gzUM8kPdelhlDwNdalhlQbsNpE1ZS8xroekVrQAUcOOHIg2pmefwBDop6K'
S3_COMPRESSED_B64 = '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'

BASE = PROJECT / 'data' / 'p8_cross_series_independence' / 'fogg_classifier_inputs_v350f' / '20261108T015702Z'
FROZEN = PROJECT / 'results' / 'p8_v352d_frozen_score_20261108T033500Z'
PREFERRED_PYTHON = PROJECT / '.envs' / 'p8_fogg_classifier' / 'bin' / 'python3'
MATRICES = {
    'train_counts': BASE / 'fit' / 'all_phases_core_matrisome_counts.tsv',
    'test_counts': BASE / 'test' / 'all_phases_core_matrisome_counts.tsv',
    'train_meta': BASE / 'fit' / 'all_phases_coldata.tsv',
    'test_meta': BASE / 'test' / 'all_phases_coldata.tsv',
    'frozen_summary': FROZEN / 'frozen_model_summary.csv',
    'frozen_coefficients': FROZEN / 'frozen_model_coefficients.csv',
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for piece in iter(lambda:f.read(4 * 1024 * 1024), b''):
            h.update(piece)
    return h.hexdigest()


def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def frozen_reference_bytes():
    s2 = zlib.decompress(base64.b64decode(S2_COMPRESSED_B64))
    s3 = zlib.decompress(base64.b64decode(S3_COMPRESSED_B64))
    assert hashlib.sha256(s2).hexdigest() == EXPECTED_ROSTER_SHA, 'LOCKED_S2_HASH_MISMATCH'
    assert hashlib.sha256(s3).hexdigest() == EXPECTED_S3_SHA, 'LOCKED_S3_HASH_MISMATCH'
    return s2, s3


def validate_refs(s2, s3):
    roster = list(csv.DictReader(io.StringIO(s2.decode('utf-8-sig'))))
    reference = list(csv.DictReader(io.StringIO(s3.decode('utf-8-sig'))))
    if len(roster) != 216 or len(reference) != 37:
        raise RuntimeError('LOCKED_REFERENCE_ROW_COUNT_MISMATCH')
    if len({r['gsm'] for r in roster}) != 216 or len({r['gsm'] for r in reference}) != 37:
        raise RuntimeError('LOCKED_REFERENCE_DUPLICATED_GSM')
    train = [r for r in roster if r['split'] == 'TRAIN179']
    test = [r for r in roster if r['split'] == 'TEST37']
    if len(train) != 179 or len(test) != 37:
        # Defensive compatibility with a spelling change; actual frozen roster
        # is inspected, not guessed from external data.
        values = sorted({r['split'] for r in roster})
        raise RuntimeError('LOCKED_SPLIT_COUNTS_UNEXPECTED_' + ','.join(values))
    if {r['gsm'] for r in test} != {r['gsm'] for r in reference}:
        raise RuntimeError('LOCKED_TEST_GSM_SET_DISAGREES')
    for r in reference:
        if r.get('group') != 'AUTHOR_ARCHIVED_PARAMETERS':
            raise RuntimeError('REFERENCE_ARM_NOT_ARCHIVED_PARAMETERS')
        if r['truth'] not in ('0','1') or r['baseline_score'] == '':
            raise RuntimeError('REFERENCE_SCHEMA_INVALID')
    return roster, reference


def selftest():
    s2, s3 = frozen_reference_bytes()
    roster, reference = validate_refs(s2, s3)
    print('[SELFTEST] PASS_FROZEN_S2_S3_SHA_ROSTER', 'train=',sum(r['split']=='TRAIN179' for r in roster),
          'test=',len(reference))


def make_report(run_dir):
    return {'stage':STAGE, 'status':'HOLD_NOT_STARTED', 'created_utc':datetime.now(timezone.utc).isoformat(),
            'policy':'PROJECT_READ_ONLY_OUTPUT_DOWNLOADS_ONLY; NO_TEST_TUNING; NO_NEW_SCIENCE',
            'project':str(PROJECT), 'output_dir':str(run_dir),
            'not_raw_cel_replay':True, 'not_leakage_free_validation':True,
            'outcome_informed_combat_preprocessing_already_frozen':True,
            'fixed_test_GSM_count':37, 'frozen_train_GSM_count':179,
            'scoring_reference':'S3_v430_AUTHOR_ARCHIVED_PARAMETERS',
            'reference_shas':{'S2':EXPECTED_ROSTER_SHA,'S3':EXPECTED_S3_SHA},
            'sources':{}, 'inputs':{}, 'warnings':[]}


def worker(out_dir, inputs):
    import warnings
    import numpy as np
    import pandas as pd
    import scipy
    import sklearn
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.exceptions import ConvergenceWarning

    with (out_dir / 'reference_S2.csv').open(newline='',encoding='utf-8-sig') as f:
        roster=list(csv.DictReader(f))
    with (out_dir / 'reference_S3.csv').open(newline='',encoding='utf-8-sig') as f:
        s3=list(csv.DictReader(f))
    roster_by_id={r['gsm']:r for r in roster}
    expected_by_id={r['gsm']:r for r in s3}

    def load_expression_and_metadata(kind):
        count_path=inputs[kind+'_counts']; meta_path=inputs[kind+'_meta']
        mat=pd.read_csv(count_path,sep='\t',low_memory=False)
        metadata=pd.read_csv(meta_path,sep='\t',low_memory=False)
        if len(mat.columns) != (180 if kind == 'train' else 38):
            raise RuntimeError('HOLD_'+kind.upper()+'_MATRIX_COLUMNS_'+str(len(mat.columns)))
        if len(mat) != 270:
            raise RuntimeError('HOLD_'+kind.upper()+'_FEATURE_ROWS_'+str(len(mat)))
        if str(mat.columns[0]).lower() not in ('symbol','hgnc_symbol','gene','gene_symbol'):
            raise RuntimeError('HOLD_'+kind.upper()+'_FIRST_COLUMN_'+str(mat.columns[0]))
        if 'sample_name' not in metadata.columns or 'condition' not in metadata.columns:
            raise RuntimeError('HOLD_'+kind.upper()+'_METADATA_SCHEMA_'+str(list(metadata.columns)))
        sample_names=list(map(str,mat.columns[1:]))
        if sample_names != metadata['sample_name'].astype(str).tolist():
            raise RuntimeError('HOLD_'+kind.upper()+'_METADATA_GSM_ORDER')
        frozen_ids={r['gsm'] for r in roster if r['split']==('TRAIN179' if kind=='train' else 'TEST37')}
        if set(sample_names) != frozen_ids or len(sample_names) != len(frozen_ids):
            raise RuntimeError('HOLD_'+kind.upper()+'_S2_GSM_SET')
        gene_names=mat.iloc[:,0].astype(str).tolist()
        if len(set(gene_names))!=len(gene_names):
            raise RuntimeError('HOLD_'+kind.upper()+'_FEATURE_NAME_DUPLICATES')
        try:
            X=mat.iloc[:,1:].to_numpy(dtype=np.float64).T
        except (TypeError,ValueError) as e:
            raise RuntimeError('HOLD_'+kind.upper()+'_NONNUMERIC_MATRIX') from e
        if X.shape!=(179 if kind=='train' else 37,270) or not np.isfinite(X).all():
            raise RuntimeError('HOLD_'+kind.upper()+'_MATRIX_SHAPE_OR_FINITE')
        truth=np.asarray([int(roster_by_id[g]['truth']) for g in sample_names],dtype=np.int64)
        if not set(truth.tolist()) <= {0,1}:
            raise RuntimeError('HOLD_'+kind.upper()+'_TRUTH_NOT_BINARY')
        # Check expected condition labels, without relying on guessed patient IDs.
        for gsm,cond in zip(sample_names,metadata['condition'].astype(str)):
            v=cond.strip().lower()
            if v in ('healthy','control','normal','nupp') and int(roster_by_id[gsm]['truth']) != 0:
                raise RuntimeError('HOLD_'+kind.upper()+'_CONDITION_TRUTH_MISMATCH')
            if v in ('endometriosis','endo','case') and int(roster_by_id[gsm]['truth']) != 1:
                raise RuntimeError('HOLD_'+kind.upper()+'_CONDITION_TRUTH_MISMATCH')
        return X,truth,sample_names,gene_names

    train,y_train,train_ids,train_genes=load_expression_and_metadata('train')
    test,y_test,test_ids,test_genes=load_expression_and_metadata('test')
    if train_genes != test_genes:
        raise RuntimeError('HOLD_TRAIN_TEST_FEATURE_ORDER')
    summary=pd.read_csv(inputs['frozen_summary'],low_memory=False)
    if len(summary)==0:
        raise RuntimeError('HOLD_FROZEN_MODEL_SUMMARY_EMPTY')
    cdict={str(c).strip().lower():c for c in summary.columns}
    config_key=next((cdict[k] for k in ('arm','configuration','model_arm') if k in cdict),None)
    if config_key is None:
        raise RuntimeError('HOLD_FROZEN_ARM_COLUMN_MISSING columns='+','.join(map(str,summary.columns)))
    selected=summary.loc[summary[config_key].astype(str).str.strip()=='AUTHOR_ARCHIVED_PARAMETERS']
    if len(selected)!=1:
        raise RuntimeError('HOLD_AUTHOR_ARCHIVED_ARM_COUNT_'+str(len(selected))+' columns='+','.join(map(str,summary.columns)))
    pars=selected.iloc[0]
    c_col=next((cdict[k] for k in ('c','c_value','cval') if k in cdict),None)
    l_col=next((cdict[k] for k in ('l1_ratio','l1ratio') if k in cdict),None)
    if c_col is None or l_col is None:
        raise RuntimeError('HOLD_FROZEN_C_OR_L1_RATIO_NOT_EXPLICIT columns='+','.join(map(str,summary.columns)))
    C=float(pars[c_col]); l1=float(pars[l_col])
    if not (0<C<1e9 and 0<=l1<=1):
        raise RuntimeError('HOLD_INVALID_FROZEN_HYPERPARAMETERS')
    mi_col=next((cdict[k] for k in ('max_iter','max_iterations') if k in cdict),None)
    max_iter=int(float(pars[mi_col])) if mi_col is not None and str(pars[mi_col]).strip().lower() not in ('','nan','none') else 100
    if max_iter not in (100,5000):
        raise RuntimeError('HOLD_UNEXPECTED_MAX_ITER_'+str(max_iter))
    source=(PROJECT / R_SCRIPT).read_text(encoding='utf-8')
    if 'StandardScaler' not in source or 'LogisticRegression' not in source or 'random_state=123' not in source:
        raise RuntimeError('HOLD_FROZEN_SOURCE_MODEL_SPEC_NOT_IDENTIFIED')
    print('[FROZEN_AUTHOR_CONFIG] C=',C,' l1_ratio=',l1,' max_iter=',max_iter,sep='',flush=True)
    print('[INPUT_MATRIX] fit=',train.shape,' test=',test.shape,' genes_equal=',train_genes==test_genes,flush=True)

    model=make_pipeline(StandardScaler(),LogisticRegression(penalty='elasticnet', solver='saga',
            C=C,l1_ratio=l1,random_state=123,max_iter=max_iter,n_jobs=1,tol=1e-4))
    with warnings.catch_warnings(record=True) as wrn:
        warnings.simplefilter('always',ConvergenceWarning)
        model.fit(train,y_train)
    fit_warnings=[str(w.message) for w in wrn]
    scores=model.decision_function(test)
    predictions=(scores>=0.0).astype(int)
    reference=np.asarray([float(expected_by_id[g]['baseline_score']) for g in test_ids],dtype=np.float64)
    ref_truth=np.asarray([int(expected_by_id[g]['truth']) for g in test_ids],dtype=int)
    if not np.array_equal(y_test,ref_truth):
        raise RuntimeError('HOLD_S3_TEST_TRUTH_DISAGREES')
    reference_predictions=(reference>=0.0).astype(int)
    err=np.abs(scores-reference)
    tp=int(np.sum((y_test==1)&(predictions==1)))
    tn=int(np.sum((y_test==0)&(predictions==0)))
    fp=int(np.sum((y_test==0)&(predictions==1)))
    fn=int(np.sum((y_test==1)&(predictions==0)))
    ba=0.5*(tp/(tp+fn)+tn/(tn+fp))
    match_all=np.array_equal(predictions,reference_predictions)
    max_abs=float(err.max()); n_wrong=int(np.sum(err>1e-6))
    expected_tp,expected_fn,expected_tn,expected_fp=19,2,16,0
    pass_all=(match_all and max_abs<=1e-6 and (tp,fn,tn,fp)==
                 (expected_tp,expected_fn,expected_tn,expected_fp))
    status='PASS_FROZEN_AUTHOR_MODEL_37_SCORES_1E6' if pass_all else 'HOLD_FROZEN_AUTHOR_MODEL_REPLAY_MISMATCH'
    with (out_dir/'37_predictions_comparison.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['gsm','truth','reference_score','fresh_score','absolute_error',
                         'reference_prediction','fresh_prediction','prediction_match'])
        for gsm,truth,rs,sc,e,rp,p in zip(test_ids,y_test,reference,scores,err,reference_predictions,predictions):
            writer.writerow([gsm,int(truth),format(float(rs),'.17g'),format(float(sc),'.17g'),
                             format(float(e),'.17g'),int(rp),int(p),int(rp)==int(p)])
    results={'status':status,'model':'AUTHOR_ARCHIVED_PARAMETERS',
             'n_train':len(y_train),'n_test':len(y_test),'n_features':len(train_genes),
             'fixed_parameters':{'C':C,'l1_ratio':l1,'max_iter':max_iter,
                'random_state':123,'solver':'saga','penalty':'elasticnet','tol':1e-4},
             'fit_warnings':fit_warnings,'n_fit_warnings':len(fit_warnings),
             'confusion':{'TP':tp,'FN':fn,'TN':tn,'FP':fp},
             'balanced_accuracy':ba,'reference_confusion':
                 {'TP':expected_tp,'FN':expected_fn,'TN':expected_tn,'FP':expected_fp},
             'n_class_disagreements':int(np.sum(predictions!=reference_predictions)),
             'max_score_absolute_error':max_abs,'n_score_disagreements_gt_1e6':n_wrong,
             'scoring_reference':'published v4.30 supplementary S3 per-GSM baseline_score',
             'packages':{'python':sys.version.split()[0], 'numpy':np.__version__,
                         'pandas':pd.__version__, 'scipy':scipy.__version__, 'sklearn':sklearn.__version__},
             'notes':['Archived post-ComBat input matrices are used; the ComBat step already consumed 37 true test diagnoses.',
                      'Fixed test labels are only used to assess predictions, not to select features or parameters.',
                      'No claim of independent clinical validation or reproduction from raw CEL.']}
    write_json(out_dir/'worker_metrics.json',results)
    print('[CONFUSION]',tp,fn,tn,fp,'BA=',format(ba,'.10f'),flush=True)
    print('[SCORE_MATCH] max_abs=',format(max_abs,'.12g'),'n_gt_1e-6=',n_wrong,
          'classification_disagreements=',results['n_class_disagreements'],flush=True)
    print('[WORKER_STATUS]',status,flush=True)
    return 0 if pass_all else 2


def run():
    s2,s3=frozen_reference_bytes()
    validate_refs(s2,s3)
    run_dir=OUTPUT_ROOT/('Run_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    run_dir.mkdir(parents=True,exist_ok=False)
    report=make_report(run_dir)
    (run_dir/'reference_S2.csv').write_bytes(s2)
    (run_dir/'reference_S3.csv').write_bytes(s3)
    try:
        print('PAPER 8 | STAGE 7 | FROZEN CORE CLASSIFIER NUMERIC REPLAY')
        print('[POLICY] PROJECT READ ONLY | RUN OUTPUTS DOWNLOADS ONLY | NO TEST TUNING')
        if not PROJECT.is_dir():
            raise RuntimeError('HOLD_PROJECT_ROOT_MISSING')
        for name,expected in EXPECTED_SOURCES.items():
            src=PROJECT/name
            if not src.is_file():
                raise RuntimeError('HOLD_SOURCE_MISSING_'+name)
            actual=sha256(src)
            report['sources'][name]={'sha256':actual,'match':actual==expected}
            if actual!=expected:
                raise RuntimeError('HOLD_SOURCE_CHANGED_'+name)
            print('[SOURCE]',name,'MATCH',flush=True)
        for name,path in MATRICES.items():
            if not path.is_file():
                raise RuntimeError('HOLD_INPUT_MISSING_'+name)
            record={'relative_path':str(path.relative_to(PROJECT)),
                    'bytes':path.stat().st_size,'sha256':sha256(path)}
            report['inputs'][name]=record
            print('[INPUT]',name,record['bytes'],'bytes',flush=True)
        if not PREFERRED_PYTHON.is_file():
            raise RuntimeError('HOLD_FROZEN_ENVIRONMENT_MISSING')
        print('[RUN] Independent frozen-fit model in locked classifier environment',flush=True)
        argv=[str(PREFERRED_PYTHON),str(Path(__file__).resolve()),'--worker',str(run_dir)]
        env=os.environ.copy()
        env.update({'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
                    'VECLIB_MAXIMUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1',
                    'PYTHONDONTWRITEBYTECODE':'1'})
        completed=subprocess.run(argv,cwd=str(run_dir),env=env,stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,text=True,timeout=900)
        (run_dir/'worker_stdout.log').write_text(completed.stdout,encoding='utf-8')
        (run_dir/'worker_stderr.log').write_text(completed.stderr,encoding='utf-8')
        for ln in completed.stdout.splitlines():
            print(ln,flush=True)
        if completed.stderr:
            print('[WORKER_STDERR_TAIL]',completed.stderr[-1400:],flush=True)
        report['worker_exit_code']=completed.returncode
        if (run_dir/'worker_metrics.json').is_file():
            metrics=json.loads((run_dir/'worker_metrics.json').read_text(encoding='utf-8'))
            report['model_metrics']=metrics
            report['status']=metrics['status']
        else:
            report['status']='HOLD_WORKER_NO_METRICS'
            report['worker_error_tail']=completed.stderr[-2000:]
        if completed.returncode and report['status'].startswith('PASS_'):
            report['status']='HOLD_WORKER_EXIT_INCONSISTENT'
    except Exception as e:
        report['status']='HOLD_STAGE7_GATE'
        report['error']=str(e)
        (run_dir/'stage7_failure_traceback.txt').write_text(traceback.format_exc(),encoding='utf-8')
    finally:
        write_json(run_dir/'stage7_replay.json',report)
        (run_dir/'README_STAGE7_SCOPE.txt').write_text(
            'Paper 8 Stage 7, numerical replay of archived AUTHOR_ARCHIVED_PARAMETERS core classifier.\n'
            'Project is read-only. Outputs and embedded S2/S3 are saved only in Downloads.\n'
            'Does not use raw CEL, re-run RMA or ComBat, claim leakage-free validation, or perform test tuning.\n'
            'If score reproduction fails, HOLD; no parameter search is carried out.\n',encoding='utf-8')
        print('[OUTPUT]',run_dir,flush=True)
        print('[STATUS]',report['status'],flush=True)
    return 0 if report['status']=='PASS_FROZEN_AUTHOR_MODEL_37_SCORES_1E6' else 2


if __name__=='__main__':
    if len(sys.argv)==2 and sys.argv[1]=='--selftest':
        selftest()
    elif len(sys.argv)==3 and sys.argv[1]=='--worker':
        try:
            sys.exit(worker(Path(sys.argv[2]),MATRICES))
        except Exception:
            traceback.print_exc()
            sys.exit(2)
    elif len(sys.argv)==1:
        sys.exit(run())
    else:
        print('Usage: python3 ~/Downloads/paper8_repro_stage7_classifier.py [--selftest]',file=sys.stderr)
        sys.exit(2)
