#!/usr/bin/env python3
"""Portable Paper 8 Stage7 model replay from exact pinned inputs.
This preserves the archived Stage 7 statistical model; no test tuning or hyperparameter search.
Project-specific source SHA and input provenance are checked by portable_cli.py.
"""
import csv, json, sys
from pathlib import Path

def write_json(path,payload):
    Path(path).write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")

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

    expected_versions={'python':'3.12.16','numpy':'2.5.4','pandas':'3.1.1',
                       'scipy':'1.18.1','sklearn':'1.9.1'}
    actual_versions={'python':sys.version.split()[0],'numpy':np.__version__,
                     'pandas':pd.__version__,'scipy':scipy.__version__,
                     'sklearn':sklearn.__version__}
    pinned_environment=(actual_versions==expected_versions)
    print('[ENV_LOCK]', 'MATCH' if pinned_environment else 'HOLD_VERSION_MISMATCH',
          actual_versions,flush=True)

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
    # Portable source check: the frozen input hashes and explicit model configuration
    # are verified before this worker runs. No historical project source tree is required.
    if abs(C-0.6552218193440044)>1e-12 or abs(l1-0.0675608348364084)>1e-12 or max_iter!=100:
        raise RuntimeError('HOLD_FROZEN_MODEL_PARAMETERS_DIFFER')
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
    numeric_match=(match_all and max_abs<=1e-6 and (tp,fn,tn,fp)==
                   (expected_tp,expected_fn,expected_tn,expected_fp))
    pass_all=numeric_match and pinned_environment
    status=('PASS_FROZEN_AUTHOR_MODEL_37_SCORES_1E6' if pass_all else
            'HOLD_UNPINNED_ENVIRONMENT' if numeric_match else
            'HOLD_FROZEN_AUTHOR_MODEL_REPLAY_MISMATCH')
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
             'pinned_environment':pinned_environment,'expected_versions':expected_versions,
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

if __name__ == '__main__':
    import json, traceback
    if len(sys.argv)!=3:
        print('Expected output directory and JSON file of flat input paths',file=sys.stderr)
        sys.exit(2)
    try:
        output=Path(sys.argv[1]); paths=json.loads(Path(sys.argv[2]).read_text(encoding='utf-8'))
        keys=('train_counts','test_counts','train_meta','test_meta','frozen_summary','frozen_coefficients')
        inputs={k:Path(paths[k]) for k in keys}
        sys.exit(worker(output,inputs))
    except Exception:
        traceback.print_exc()
        sys.exit(2)
