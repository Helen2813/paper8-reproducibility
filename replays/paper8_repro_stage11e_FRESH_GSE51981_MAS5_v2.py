#!/usr/bin/env python3
"""Paper 8 Stage 11E -- fresh-namespace raw-CEL -> MAS5 replay.

Recomputes independent affy::mas5(sc=500) for *all* 148 GSE51981 CELs,
extracts the exact frozen 324-probe panel, and compares 47,952 values with
v4.22, the archived sample-level exploratory benchmark input matrix.

No RMA, ComBat, models, downloads, parameter tuning, or project writes.
No CLI arguments required; run from ~/Downloads on macOS (Python 3.9+).

STRICT NO-REUSE VERSION. Every invocation creates a unique, initially empty
checkpoint directory INSIDE its uniquely named Run_... folder. It executes R
separately for all 148 CELs, without consulting checkpoints from other runs.
Interrupted runs retain partial evidence for diagnosis but are NOT resumed:
rerunning this script starts a fresh, separate run. A PASS requires exactly
148 successful new R executions and zero reused checkpoints. All original
research inputs remain read-only. Never publish reports without reviewing
local-path metadata.
"""
from __future__ import print_function
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import traceback
import zipfile

ROOT = Path.home() / 'Desktop/paper8/endo-research'
BASE = Path.home() / 'Downloads/Paper8_repro_stage11E_FRESH_GSE51981_MAS5'
RAW = ROOT / 'data/external/p8_v334_gse51981_cel/GSE51981_RAW.tar'
V422 = ROOT / 'results/p8_v422_locked_GSE51981_MAS5_nested_cv_RESUMABLE.zip'
V410A = ROOT / 'results/p8_v410a_GSE7846_GPL570_source_only_RESUMABLE.zip'
R_BIN = ROOT / '.envs/p8_fogg_r44/bin/Rscript'
ARCHIVE_MATRIX = '05_148_x_324_single_CEL_MAS5_intensity.tsv'
PANEL_BASENAME = '01_frozen_324_probe_ids.txt'
N_CELS, N_PROBES = 148, 324
MAX_CEL_BYTES = 90_000_000
MAX_RAW_MEMBER_BYTES = 90_000_000
MIN_FREE_BYTES = 2 * 1024**3
TOL = 1e-6
GSM_RE = re.compile(r'GSM[0-9]+', re.IGNORECASE)
EXPECTED_RSCRIPT_SHA256 = 'cfcc3adb69274ad4dfc8fff5152a6e476ac5c1b949540a257d011036e363fcbd'

# Matches frozen v4.22's one-array ReadAffy + MAS5(sc=500) protocol.
# This worker was written independently; the original v4.22 driver is NOT executed.
R_WORKER = r'''#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(hgu133plus2cdf))
suppressPackageStartupMessages(library(Biobase))
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=3L)stop('HOLD_EXPECTED_CEL_PANEL_OUTPUT')
panel <- readLines(a[[2L]],warn=FALSE)
if(length(panel)!=324L || anyDuplicated(panel))stop('HOLD_PANEL_324_INVALID')
x <- affy::ReadAffy(filenames=a[[1L]])
if(ncol(Biobase::exprs(x))!=1L)stop('HOLD_NOT_ONE_ARRAY')
y <- affy::mas5(x,sc=500)
v <- Biobase::exprs(y)[,1L]
names(v) <- rownames(Biobase::exprs(y))
if(any(!panel %in% names(v)))stop('HOLD_FROZEN_PROBES_UNAVAILABLE')
v <- v[panel]
if(length(v)!=324L || anyNA(v) || any(!is.finite(v)) || any(v<=0))stop('HOLD_INVALID_MAS5_SIGNALS')
writeLines(c('probe_id\tmas5_signal',paste0(panel,'\t',sprintf('%.17g',as.numeric(v)))),a[[3L]],useBytes=TRUE)
'''

class Hold(RuntimeError):
    pass

def say(*items):
    print(*items, flush=True)

def hash_bytes(data):
    return hashlib.sha256(data).hexdigest()

def hash_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4*1024*1024),b''):
            h.update(b)
    return h.hexdigest()

def as_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)+'\n',encoding='utf-8')

def parse_matrix(data):
    with io.StringIO(data.decode('utf-8-sig'),newline='') as text:
        rows=csv.reader(text,delimiter='\t')
        header=next(rows,None)
        if not header or len(header)!=N_PROBES+1 or header[0]!='gsm':
            raise Hold('V422_MAS5_MATRIX_HEADER')
        panel=[x.lower() for x in header[1:]]
        if len(set(panel))!=N_PROBES or any(not p for p in panel):
            raise Hold('V422_PANEL_SCHEMA')
        values={}
        for row in rows:
            if len(row)!=N_PROBES+1:raise Hold('V422_ROW_WIDTH')
            gsm=row[0].strip().upper()
            if not re.fullmatch('GSM[0-9]+',gsm) or gsm in values:
                raise Hold('V422_INVALID_OR_DUPLICATE_GSM')
            nums=[float(v) for v in row[1:]]
            if any(not math.isfinite(n) or n<=0 for n in nums):
                raise Hold('V422_MAS5_NONPOSITIVE_OR_NAN')
            values[gsm]=nums
        if len(values)!=N_CELS:raise Hold('V422_NOT_148_SAMPLES')
    return panel,values

def parse_single_vector(data, panel):
    with io.StringIO(data.decode('utf-8-sig'),newline='') as text:
        r=list(csv.reader(text,delimiter='\t'))
    if len(r)!=N_PROBES+1 or r[0]!=['probe_id','mas5_signal']:
        raise Hold('R_MAS5_SINGLE_VECTOR_FORMAT')
    probe_ids=[x[0].lower() for x in r[1:]]
    if probe_ids!=panel:raise Hold('R_MAS5_FROZEN_PROBE_ORDER_MISMATCH')
    vals=[float(x[1]) for x in r[1:]]
    if any(not math.isfinite(v) or v<=0 for v in vals):
        raise Hold('R_MAS5_NONPOSITIVE_OR_NAN')
    return vals

def max_diff(ref, obs, panel):
    if len(ref)!=N_PROBES or len(obs)!=N_PROBES:
        raise Hold('VECTOR_WIDTH')
    diffs=[abs(a-b) for a,b in zip(ref,obs)]
    maxval=max(diffs)
    worst=panel[diffs.index(maxval)]
    return {'max_abs':maxval,'mean_abs':sum(diffs)/N_PROBES,
            'n_abs_gt_1e_6':sum(d>TOL for d in diffs),
            'n_abs_gt_1e_10':sum(d>1e-10 for d in diffs),
            'worst_probe':worst}

def load_sources():
    for p in (RAW,V422,V410A,R_BIN):
        if not p.is_file() or p.is_symlink():raise Hold('MISSING_OR_SYMLINK_SOURCE_'+str(p))
    if hash_file(R_BIN)!=EXPECTED_RSCRIPT_SHA256:
        raise Hold('RSCRIPT_BINARY_SHA256_CHANGED')
    hashes={'raw_tar':hash_file(RAW), 'v422_zip':hash_file(V422),
            'v410a_zip':hash_file(V410A), 'Rscript':hash_file(R_BIN),
            'independent_worker':hash_bytes(R_WORKER.encode('utf-8'))}
    with zipfile.ZipFile(V422) as z:
        if z.testzip() is not None:raise Hold('V422_ARCHIVE_CRC_FAILED')
        if ARCHIVE_MATRIX not in z.namelist() or '00_manifest.json' not in z.namelist():
            raise Hold('V422_MISSING_MATRIX_OR_MANIFEST')
        matrix=z.read(ARCHIVE_MATRIX)
        manifest=json.loads(z.read('00_manifest.json').decode('utf-8'))
    if manifest.get('status')!='PASS_V422_LOCKED_INTERNAL_ANALYSIS_WITH_HOLDS':
        raise Hold('V422_UNEXPECTED_STATUS')
    matrix_hash=hash_bytes(matrix)
    if manifest.get('output_sha256',{}).get(ARCHIVE_MATRIX)!=matrix_hash:
        raise Hold('V422_FROZEN_MATRIX_SHA_MISMATCH')
    source_pins=manifest.get('source_files',{})
    if source_pins.get('GSE51981_raw_tar_sha256')!=hashes['raw_tar']:
        raise Hold('V422_SOURCE_GEO_TAR_HASH_MISMATCH')
    panel,reference=parse_matrix(matrix)
    with zipfile.ZipFile(V410A) as z:
        if z.testzip() is not None:raise Hold('V410A_ARCHIVE_CRC_FAILED')
        candidates=[n for n in z.namelist() if Path(n).name==PANEL_BASENAME]
        if len(candidates)!=1:raise Hold('V410A_FROZEN_PANEL_NOT_UNIQUE')
        panel_blob=z.read(candidates[0])
    panel_hash=hash_bytes(panel_blob)
    if source_pins.get('frozen_panel_sha256')!=panel_hash:
        raise Hold('V422_V410A_FROZEN_PANEL_HASH_MISMATCH')
    panel_raw=[x.strip() for x in panel_blob.decode('utf-8-sig').splitlines() if x.strip()]
    if len(panel_raw)!=N_PROBES or len(set(panel_raw))!=N_PROBES:
        raise Hold('V410A_PANEL_COUNT')
    if [x.lower() for x in panel_raw]!=panel:
        raise Hold('V410A_V422_PROBE_ORDER_DIFFER')
    hashes['v422_MAS5_matrix']=matrix_hash
    hashes['v410a_panel']=panel_hash
    say('[PIN] v422 MAS5 matrix SHA PASS; V410a frozen panel SHA PASS; original GEO TAR SHA PASS')
    return panel,panel_blob,reference,hashes

def enumerate_gsm_members():
    result={}
    with tarfile.open(RAW,mode='r:*') as tar:
        for member in tar:
            base=Path(member.name).name
            if not(base.lower().endswith('.cel') or base.lower().endswith('.cel.gz')):
                continue
            if '\\' in member.name or '\x00' in member.name or '..' in Path(member.name).parts:
                raise Hold('UNSAFE_TAR_MEMBER_PATH')
            matches=set(x.upper() for x in GSM_RE.findall(base))
            if not member.isfile() or member.size>MAX_RAW_MEMBER_BYTES or len(matches)!=1:
                raise Hold('BAD_CEL_MEMBER')
            gsm=next(iter(matches))
            if gsm in result:raise Hold('DUPLICATED_GSM_IN_GEO_TAR')
            result[gsm]=member.name
    if len(result)!=N_CELS:raise Hold('GEO_TAR_NOT_148_CELS_'+str(len(result)))
    return result

def r_preflight():
    expr="""pkgs<-c('affy','Biobase','hgu133plus2cdf');for(p in pkgs) {if(!requireNamespace(p,quietly=TRUE))stop(paste('MISSING',p));cat(p,'=',as.character(packageVersion(p)),'\\n',sep='')}"""
    cmd=[str(R_BIN),'-e',expr]
    run=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
    if run.returncode:
        raise Hold('R_PACKAGE_PROBE_FAILED_'+(run.stderr or run.stdout)[-400:])
    values={}
    for line in run.stdout.splitlines():
        if '=' in line:
            k,v=line.split('=',1)
            values[k.strip()]=v.strip()
    expected={'affy':'1.84.0','Biobase':'2.66.0','hgu133plus2cdf':'2.18.0'}
    if any(values.get(k)!=v for k,v in expected.items()):
        raise Hold('R_PACKAGE_VERSION_MISMATCH_'+str(values))
    say('[R_PACKAGES]',', '.join(k+'='+v for k,v in values.items()))
    return values

def uncompress_member(tar, member_name, destination):
    member=tar.getmember(member_name)
    stream=tar.extractfile(member)
    if stream is None:raise Hold('GEO_MEMBER_CANNOT_BE_READ')
    digest=hashlib.sha256(); count=0
    with stream:
        source=gzip.GzipFile(fileobj=stream,mode='rb') if member_name.lower().endswith('.gz') else stream
        try:
            with destination.open('xb') as output:
                while True:
                    block=source.read(1024*1024)
                    if not block:break
                    count+=len(block)
                    if count>MAX_CEL_BYTES:raise Hold('CEL_UNCOMPRESSED_SIZE_LIMIT')
                    digest.update(block);output.write(block)
        finally:
            if source is not stream:source.close()
    if count<100_000:raise Hold('CEL_TRUNCATED')
    return digest.hexdigest(),count

def create_new_empty_run_cache(run_dir):
    """Allocate a fresh output-only namespace, never based on source hashes.

    run_dir itself was created by mkdtemp, so even same-second invocations
    cannot share any CEL checkpoint. An existing directory is a hard error.
    """
    cache = Path(run_dir) / 'FRESH_R_MAS5_THIS_RUN'
    cache.mkdir(mode=0o700, parents=False, exist_ok=False)
    if any(cache.iterdir()):
        raise Hold('NEW_RUN_CACHE_NOT_EMPTY')
    return cache

def write_csv(path,rows,cols):
    with path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)

def selftest():
    global N_PROBES
    old=N_PROBES
    try:
        N_PROBES=2
        matrix=b'gsm\tP1\tP2\nGSM1\t1.5\t2\nGSM2\t3\t4\n'
        # Miniature parser uses actual 148-gate, check this gate rejects incomplete archive.
        try:parse_matrix(matrix)
        except Hold as exc:
            assert 'NOT_148' in str(exc)
        else:raise AssertionError('EXPECTED_148_GATE')
        assert parse_single_vector(b'probe_id\tmas5_signal\nP1\t1.5\nP2\t2\n',['p1','p2'])==[1.5,2]
        assert max_diff([1,2],[1,2],['p1','p2'])['max_abs']==0
        assert max_diff([1,2],[1,2.01],['p1','p2'])['n_abs_gt_1e_6']==1
        with tempfile.TemporaryDirectory(prefix='p8_fresh_check_') as temp:
            check_run=Path(temp)
            fresh=create_new_empty_run_cache(check_run)
            assert fresh.parent==check_run
            assert not list(fresh.iterdir())
            try:
                create_new_empty_run_cache(check_run)
            except FileExistsError:
                pass
            else:
                raise AssertionError('MUST_REJECT_REUSED_RUN_CACHE')
    finally:
        N_PROBES=old
    return True

def run():
    say('PAPER 8 | STAGE 11E v2 | STRICT NO-CACHE-REUSE 148-CEL MAS5 REPLAY')
    say('[POLICY] SOURCE PROJECT READ ONLY | OFFLINE | NO RMA/COMBAT/CLASSIFIER | OUTPUTS DOWNLOADS')
    assert selftest()
    say('[SELFTEST] PASS_PANEL_PROBE_AND_NUMERIC_TAMPER_GATES')
    BASE.mkdir(parents=True,exist_ok=True)
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    # mkdtemp creates a unique output directory even for multiple starts in the same second.
    run_dir=Path(tempfile.mkdtemp(prefix='Run_'+stamp+'_',dir=str(BASE)))
    report={'stage':'P8_STAGE11E_FRESH_ALL148_MAS5','status':'HOLD_START','n_expected':N_CELS,
            'n_probes':N_PROBES,'n_comparisons_expected':N_CELS*N_PROBES,
            'threshold':TOL,'project_read_only':True,
            'scope':'Single-array MAS5 sc=500 for exploratory v4.22 benchmark (not original Fogg model)',
            'RMA_recomputed':False,'ComBat_recomputed':False,'classifier_fitted':False,
            'all_148_new_MAS5_computed':False,
            'n_new_R_MAS5_executions_this_run':0,
            'n_reused_provenance_verified_checkpoints':0,
            'cache_reuse_allowed':False,
            'limitations':[
                'The frozen 324-probe panel reflects prior project analyses and is not independent feature discovery.',
                'These MAS5 intensities are inputs to internal sample-level benchmark, not independent patient validation.',
                'ReadAffy + affy::mas5(sc=500) is repeated separately for each CEL; no full-series pooling.']}
    results=[]
    n_new_R_executed=0
    try:
        free=shutil.disk_usage(str(BASE)).free
        say('[DISK] Downloads free',round(free/1024**3,2),'GiB')
        if free<MIN_FREE_BYTES:raise Hold('DOWNLOADS_DISK_FREE_BELOW_2GiB')
        panel,panel_blob,reference,pins=load_sources()
        report['source_sha256']=pins
        members=enumerate_gsm_members()
        if set(members)!=set(reference):
            raise Hold('SOURCE_TAR_V422_GSM_ROSTER_MISMATCH')
        say('[GEO_TAR]',len(members),'CEL; exact match with 148 archived MAS5 GSM')
        versions=r_preflight();report['R_package_versions']=versions
        panel_file=run_dir/'frozen_324_probes.txt'
        panel_file.write_bytes(panel_blob)
        worker=run_dir/'independent_single_CEL_MAS5.R'
        worker.write_text(R_WORKER,encoding='utf-8')
        # Mandatory one-use-per-run cache. NEVER point at any earlier Stage 11D
        # or Stage 11E checkpoint folder; do not quietly resume older values.
        cache=create_new_empty_run_cache(run_dir)
        proof=cache/'STAGE11E_FRESH_SOURCE_ORIGIN.json'
        proof_record={
            'stage':'STAGE11E_v2_FRESH_ALL148',
            'source_script_sha256':hash_file(Path(__file__).resolve()),
            'raw_tar_sha256':pins['raw_tar'],
            'v422_matrix_sha256':pins['v422_MAS5_matrix'],
            'v410a_panel_sha256':pins['v410a_panel'],
            'worker_sha256':pins['independent_worker'],
            'Rscript_sha256':pins['Rscript'],
            'created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
            'origin':'NEW_EMPTY_FOLDER_UNDER_UNIQUE_RUN_DIRECTORY',
            'checkpoint_reuse_permitted':False,
        }
        as_json(proof,proof_record)
        report['checkpoint_folder']=str(cache)
        report['checkpoint_origin_proof_sha256']=hash_file(proof)
        report['checkpoint_origin']=proof_record
        say('[CACHE_ORIGIN] NEW UNIQUE EMPTY FOLDER; OLD CHECKPOINTS NEVER READ')
        say('[CHECKPOINTS]',cache)
        with tarfile.open(RAW,mode='r:*') as tar:
            for number,gsm in enumerate(sorted(reference),1):
                resultfile=cache/(gsm+'.tsv')
                metafile=cache/(gsm+'.json')
                if resultfile.exists() or metafile.exists():
                    raise Hold('UNEXPECTED_OLD_CHECKPOINT_IN_FRESH_RUN_'+gsm)
                with tempfile.TemporaryDirectory(prefix='one_CEL_',dir=str(run_dir)) as temp:
                    t=Path(temp)
                    src=t/(gsm+'.CEL')
                    cel_sha,size=uncompress_member(tar,members[gsm],src)
                    out=t/'vector.tsv'
                    localtmp=t/'r_temp'
                    localtmp.mkdir()
                    env=os.environ.copy()
                    env['TMPDIR']=str(localtmp)
                    args=[str(R_BIN),str(worker),str(src),str(panel_file),str(out)]
                    # A new R process ALWAYS computes this GSM in THIS invocation.
                    p=subprocess.run(args,cwd=str(t),env=env,capture_output=True,
                                     text=True,timeout=480)
                    if p.returncode!=0:
                        (run_dir/('R_error_'+gsm+'.txt')).write_text(
                            p.stdout+'\n'+p.stderr,encoding='utf-8')
                        raise Hold('MAS5_R_WORKER_FAILED_'+gsm+'_'+str(p.returncode))
                    if not out.is_file():
                        raise Hold('MAS5_WORKER_OUTPUT_MISSING_'+gsm)
                    # Count only successful new R processes producing a file.
                    n_new_R_executed+=1
                    report['n_new_R_MAS5_executions_this_run']=n_new_R_executed
                    values=parse_single_vector(out.read_bytes(),panel)
                    stat=max_diff(reference[gsm],values,panel)
                    if stat['n_abs_gt_1e_6']:
                        (run_dir/('HOLD_vector_'+gsm+'.tsv')).write_bytes(out.read_bytes())
                        raise Hold('MAS5_NUMERIC_DISAGREEMENT_'+gsm+'_MAX_'+str(stat['max_abs']))
                    shutil.copyfile(out,resultfile)
                    as_json(metafile,{
                        'gsm':gsm,'source_tar_sha256':pins['raw_tar'],
                        'panel_sha256':pins['v410a_panel'],
                        'worker_sha256':pins['independent_worker'],
                        'reference_sha256':pins['v422_MAS5_matrix'],
                        'vector_sha256':hash_file(resultfile),
                        'cel_sha256':cel_sha,'cel_uncompressed_bytes':size,
                        'method':'fresh_single_CEL_affy_mas5_sc500',
                        'computed_in_this_run':True,
                    })
                values=parse_single_vector(resultfile.read_bytes(),panel)
                stat=max_diff(reference[gsm],values,panel)
                if stat['n_abs_gt_1e_6']:
                    raise Hold('FRESH_R_MAS5_NUMERIC_COMPARISON_DIFFER_'+gsm)
                results.append({
                    'gsm':gsm,'status':'MATCH','reused_checkpoint':0,
                    'cel_sha256':cel_sha,'vector_sha256':hash_file(resultfile),**stat
                })
                if number%10==0 or number==N_CELS:
                    say('[MAS5_COMPARE]',number,'/',N_CELS,'GSM; max_abs=',
                        format(max(r['max_abs'] for r in results),'.12g'))
        cols=['gsm','status','reused_checkpoint','cel_sha256','vector_sha256',
              'max_abs','mean_abs','n_abs_gt_1e_6','n_abs_gt_1e_10','worst_probe']
        write_csv(run_dir/'148_CEL_comparison.csv',results,cols)
        # Write an independently generated wide matrix, identical layout/order to v4.22.
        with (run_dir/'148x324_MAS5_reconstructed.tsv').open('w',newline='',encoding='utf-8') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n')
            w.writerow(['gsm']+panel)
            for gsm in sorted(reference):
                v=parse_single_vector((cache/(gsm+'.tsv')).read_bytes(),panel)
                w.writerow([gsm]+[format(a,'.17g') for a in v])
        report['comparison']={'n_CEL':N_CELS,'n_probes':N_PROBES,
           'n_values':N_CELS*N_PROBES,'max_abs':max(r['max_abs'] for r in results),
           'n_gt_1e_6':sum(r['n_abs_gt_1e_6'] for r in results),
           'n_gt_1e_10':sum(r['n_abs_gt_1e_10'] for r in results),
           'reused_cached_CEL':sum(r['reused_checkpoint'] for r in results)}
        report['reconstructed_wide_matrix_sha256']=hash_file(run_dir/'148x324_MAS5_reconstructed.tsv')
        n_reused=sum(r['reused_checkpoint'] for r in results)
        report['n_new_R_MAS5_executions_this_run']=n_new_R_executed
        report['n_reused_provenance_verified_checkpoints']=n_reused
        report['all_148_new_MAS5_computed']=(n_new_R_executed==N_CELS and n_reused==0)
        if not report['all_148_new_MAS5_computed']:
            raise Hold('FRESH_R_EXECUTION_COUNT_NOT_148_OR_REUSED_CHECKPOINTS')
        report['status']='PASS_FRESH_148_RAW_CEL_MAS5_324_NUMERIC_REPLAY_1E6'
        say('[R_EXECUTIONS_THIS_RUN]',n_new_R_executed,'[VERIFIED_CHECKPOINTS]',n_reused)
        say('[RESULT]',report['comparison'])
    except BaseException as exc:
        report['status']='HOLD_STAGE11E_'+type(exc).__name__
        report['error']=str(exc)
        (run_dir/'traceback.log').write_text(traceback.format_exc(),encoding='utf-8')
        say('[ERROR]',str(exc))
        if isinstance(exc,KeyboardInterrupt):say('[NOTE] Keyboard interrupt, checkpointed GSM results preserved')
    finally:
        report['n_new_R_MAS5_executions_this_run']=n_new_R_executed
        report['n_GSM_completed_in_this_attempt']=len(results)
        as_json(run_dir/'stage11e_replay.json',report)
        if results:
            cols=['gsm','status','reused_checkpoint','cel_sha256','vector_sha256',
                  'max_abs','mean_abs','n_abs_gt_1e_6','n_abs_gt_1e_10','worst_probe']
            write_csv(run_dir/'completed_CEL_comparisons.csv',results,cols)
        # Private evidence bundle excludes raw CEL and cached per-CEL vector files.
        if report['status'].startswith('PASS_'):
            bundle=run_dir/'Paper8_stage11E_FRESH_MAS5_REPLAY_PRIVATE.zip'
            with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
                for name in ['stage11e_replay.json','148_CEL_comparison.csv',
                             '148x324_MAS5_reconstructed.tsv','independent_single_CEL_MAS5.R','frozen_324_probes.txt']:
                    z.write(run_dir/name,arcname=name)
            with zipfile.ZipFile(bundle) as z:
                if z.testzip() is not None:raise RuntimeError('OUTPUT_ZIP_CRC_ERROR')
            say('[PRIVATE_ZIP]',bundle)
        say('[STATUS]',report['status'])
        say('[OUTPUT]',run_dir)
    return 0 if report['status'].startswith('PASS_') else 2

if __name__=='__main__':
    sys.exit(run())
