#!/usr/bin/env python3
"""Paper 8 portable numerical replay core. Python >=3.9, stdlib for orchestration.
The two compute modes require pinned R and Python scientific environments.

Reruns are from preserved derived matrices. They are not a full raw-GEO pipeline,
not a leakage-free medical validation and not a replay of historical null drivers.
"""
import argparse, csv, datetime, hashlib, json, os, shutil, subprocess, sys, traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CONF=ROOT/'portable'/'input_manifest.json'

def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(4*1024*1024),b''):h.update(part)
    return h.hexdigest()

def load_conf():
    d=json.loads(CONF.read_text(encoding='utf-8'))
    assert len(d['resources'])==9
    return d

def under(path,parent):
    path=Path(path).resolve();parent=Path(parent).resolve()
    return path==parent or parent in path.parents

def safe_out(path, roots):
    p=Path(path).expanduser().resolve()
    if p == Path('/') or under(p, ROOT):
        raise RuntimeError('HOLD_UNSAFE_OUTPUT_INSIDE_RELEASE')
    for r in roots:
        if r is not None and under(p,r):
            raise RuntimeError('HOLD_UNSAFE_OUTPUT_INSIDE_SOURCE_OR_INPUTS')
    return p

def ensure_bundle():
    # No clinical tables are bundled in the public code-only release.
    conf=load_conf()
    refs=conf.get('external_references',{})
    if set(refs)!={'s2','s3'}:
        raise RuntimeError('HOLD_EXTERNAL_REFERENCE_MANIFEST_SCHEMA')
    for key,ref in refs.items():
        if ref.get('portable_filename') not in ('reference_S2.csv','reference_S3.csv'):
            raise RuntimeError('HOLD_EXTERNAL_REFERENCE_FILENAME_'+key)
        if len(str(ref.get('sha256','')))!=64:
            raise RuntimeError('HOLD_EXTERNAL_REFERENCE_HASH_FORMAT_'+key)
    return conf

def check_external_references(conf,inputs):
    # Fail closed: never substitute synthetic, absent or unverified reference tables.
    found={}
    for key,spec in conf['external_references'].items():
        p=inputs/spec['portable_filename']
        if p.is_symlink() or not p.is_file():
            raise RuntimeError('HOLD_PRIVATE_REFERENCE_MISSING_'+key)
        if sha256(p)!=spec['sha256']:
            raise RuntimeError('HOLD_PRIVATE_REFERENCE_HASH_'+key)
        print('[REFERENCE]',key,'HASH_MATCH',flush=True)
        found[key]=p
    return found

def check_inputs(conf,inputs,keys):
    found={}
    for k in keys:
        spec=conf['resources'][k]
        p=inputs/spec['portable_filename']
        if p.is_symlink() or not p.is_file():
            raise RuntimeError('HOLD_PORTABLE_INPUT_MISSING_OR_SYMLINK_'+k)
        actual=sha256(p)
        if actual!=spec['sha256']:
            raise RuntimeError('HOLD_PORTABLE_INPUT_HASH_'+k)
        found[k]={'name':spec['portable_filename'],'sha256':actual,'bytes':p.stat().st_size}
        print('[INPUT]',k,'HASH_MATCH',flush=True)
    return found

def script_run(cmd,cwd,out):
    env=os.environ.copy()
    env.update({'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'})
    cp=subprocess.run([str(s) for s in cmd],cwd=str(cwd),env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (out/'stdout.log').write_text(cp.stdout,encoding='utf-8')
    (out/'stderr.log').write_text(cp.stderr,encoding='utf-8')
    print(cp.stdout[-6000:],flush=True)
    if cp.stderr:print('[STDERR_TAIL]',cp.stderr[-1300:],flush=True)
    return cp.returncode

def export_inputs(a):
    conf=ensure_bundle()
    source=Path(a.source_root).expanduser().resolve()
    if not source.is_dir():raise RuntimeError('HOLD_SOURCE_ROOT_MISSING')
    output=safe_out(a.inputs_dir,[source])
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('HOLD_EXPORT_DIR_NOT_EMPTY: choose a new directory (no overwrites)')
    if shutil.disk_usage(output.parent if output.parent.is_dir() else Path.home()).free < 1_000_000_000:
        raise RuntimeError('HOLD_LOW_FREE_SPACE')
    # Verify ALL 9 source files before copying any bytes.
    matches=[]
    for k,spec in conf['resources'].items():
        p=source/spec['source_relative_path']
        if p.is_symlink() or not p.is_file() or sha256(p)!=spec['sha256']:
            raise RuntimeError('HOLD_ORIGINAL_INPUT_PROVENANCE_'+k)
        matches.append((k,p,spec))
        print('[SOURCE]',k,'HASH_MATCH',flush=True)
    output.mkdir(parents=True,exist_ok=True)
    copied=[]
    for k,p,spec in matches:
        dest=output/spec['portable_filename']
        with p.open('rb') as fi,dest.with_suffix(dest.suffix+'.partial').open('xb') as fo:
            shutil.copyfileobj(fi,fo,1024*1024)
        temp=dest.with_suffix(dest.suffix+'.partial')
        if sha256(temp)!=spec['sha256']:
            raise RuntimeError('HOLD_COPY_HASH_MISMATCH_'+k)
        temp.replace(dest)
        copied.append(k)
        print('[EXPORTED]',k,flush=True)
    (output/'INPUTS_PRIVATE_MANIFEST.json').write_text(json.dumps({'status':'PASS_PORTABLE_COPY_SOURCE_HASHED','resources':copied,'note':'Private expression matrices; do not publish without provenance and permissions review.'},indent=2)+'\n')
    print('[STATUS] PASS_PORTABLE_COPY_SOURCE_HASHED')
    print('[NOTE] Two external published supplementary references must be separately copied and SHA-verified before classifier replay.')
    print('[INPUTS_DIR]',output)

def fresh_out(base,label):
    base=Path(base).expanduser().resolve()
    base.mkdir(parents=True,exist_ok=True)
    d=base/(label+'_'+now())
    i=1
    while d.exists():d=base/(label+'_'+now()+'_'+str(i));i+=1
    d.mkdir()
    return d

def compute(a,mode):
    conf=ensure_bundle()
    inputs=Path(a.inputs_dir).expanduser().resolve()
    if not inputs.is_dir():raise RuntimeError('HOLD_INPUTS_DIRECTORY_MISSING')
    if under(inputs,ROOT):raise RuntimeError('HOLD_PRIVATE_INPUTS_INSIDE_PUBLIC_REPOSITORY')
    results=safe_out(a.out_dir,[inputs])
    keys=list(conf['resources']) if mode=='all' else (['combat_pre','combat_post','combat_meta'] if mode=='combat' else ['train_counts','test_counts','train_meta','test_meta','frozen_summary','frozen_coefficients'])
    check_inputs(conf,inputs,keys)
    references=check_external_references(conf,inputs) if mode in ('classifier','all') else {}
    if mode in ('combat','all'):
        if not a.rscript:raise RuntimeError('HOLD_RSCRIPT_ARGUMENT_REQUIRED')
        if not Path(a.rscript).expanduser().is_file():raise RuntimeError('HOLD_RSCRIPT_NOT_FOUND')
    if mode in ('classifier','all'):
        if not a.python:raise RuntimeError('HOLD_PYTHON_ARGUMENT_REQUIRED')
        if not Path(a.python).expanduser().is_file():raise RuntimeError('HOLD_PYTHON_NOT_FOUND')
    out=fresh_out(results,'p8_'+mode)
    report={'method':'Paper8_portable_code_v6','mode':mode,'status':'HOLD_NOT_FINISHED','runs':[],
            'input_shas_verified':True,'source_root_not_required':True,
            'not_leakage_free':True,'not_raw_CEL_replay':True,'not_full_3582_null_replay':True}
    try:
        if mode in ('combat','all'):
            cr=out/'combat';cr.mkdir()
            argv=[a.rscript,ROOT/'portable/combat_replay.R',inputs/'combat_pre.tsv',inputs/'combat_post.tsv',inputs/'combat_meta.tsv',cr]
            code=script_run(argv,cr,cr)
            metrics={}
            m=cr/'numeric_metrics.tsv'
            if m.exists():
                with m.open(newline='',encoding='utf-8') as f:metrics={r['key']:r['value'] for r in csv.DictReader(f,delimiter='\t')}
            ok=code==0 and metrics.get('status')=='PASS_NUMERIC_COMBAT_REPLAY_1E6' and int(metrics.get('n_abs_gt_1e-6','-1'))==0
            report['runs'].append({'name':'combat','status':'PASS' if ok else 'HOLD','exit_code':code,'metrics':metrics})
            if not ok:raise RuntimeError('HOLD_COMBAT_REPLAY')
        if mode in ('classifier','all'):
            cr=out/'classifier';cr.mkdir()
            for k in ('s2','s3'):
                shutil.copyfile(references[k],cr/('reference_S'+('2' if k=='s2' else '3')+'.csv'))
            keys=['train_counts','test_counts','train_meta','test_meta','frozen_summary','frozen_coefficients']
            paths={k:str(inputs/conf['resources'][k]['portable_filename']) for k in keys}
            (cr/'input_paths.json').write_text(json.dumps(paths,indent=2)+'\n',encoding='utf-8')
            code=script_run([a.python,ROOT/'portable/classifier_replay.py',cr,cr/'input_paths.json'],cr,cr)
            m=cr/'worker_metrics.json';metrics=json.loads(m.read_text()) if m.is_file() else {}
            ok=code==0 and metrics.get('status')=='PASS_FROZEN_AUTHOR_MODEL_37_SCORES_1E6'
            report['runs'].append({'name':'classifier','status':'PASS' if ok else 'HOLD','exit_code':code,'metrics':metrics})
            if not ok:raise RuntimeError('HOLD_CLASSIFIER_REPLAY')
        report['status']='PASS_PORTABLE_REPLAY_PINNED_INPUTS'
    except Exception as e:
        report['status']='HOLD_PORTABLE_REPLAY'
        report['reason']=str(e)
    (out/'portable_replay_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('[STATUS]',report['status'])
    print('[REPORT]',out/'portable_replay_report.json')
    return 0 if report['status'].startswith('PASS') else 2

def verify():
    ensure_bundle()
    calls=[[sys.executable,'-B',ROOT/'scripts/public_code_preflight.py'],
           [sys.executable,'-B','-m','unittest','discover','-s','tests','-v']]
    for c in calls:
        cp=subprocess.run([str(x) for x in c],cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if cp.returncode:
            print(cp.stdout,cp.stderr,file=sys.stderr)
            raise RuntimeError('HOLD_PUBLIC_CODE_SELFTEST')
    print('[STATUS] PASS_PUBLIC_CODE_STATIC_VERIFY')
    print('[SCOPE] SHA-256 and synthetic tests only; scientific input matrices and patient-level references are NOT bundled')
    return 0

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(dest='mode',required=True)
    sub.add_parser('verify')
    s=sub.add_parser('export-inputs',help='One-time local private copy after verifying nine locked SHA256s')
    s.add_argument('--source-root',required=True);s.add_argument('--inputs-dir',required=True)
    for mode in ('combat','classifier','all'):
        s=sub.add_parser(mode,help='Numeric replay of archived frozen inputs')
        s.add_argument('--inputs-dir',required=True);s.add_argument('--out-dir',required=True)
        s.add_argument('--rscript',default=None);s.add_argument('--python',default=None)
    args=ap.parse_args()
    try:
        if args.mode=='verify':return verify()
        if args.mode=='export-inputs':export_inputs(args);return 0
        return compute(args,args.mode)
    except Exception as e:
        print('[STATUS] HOLD_PORTABLE_STAGE13',type(e).__name__,str(e),file=sys.stderr)
        return 2
if __name__=='__main__':
    sys.exit(main())
