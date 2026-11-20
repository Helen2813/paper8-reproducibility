#!/usr/bin/env python3
"""Paper 8 Stage 11A: independent, bounded, per-series RMA numeric replay.

Reads the previously archived public GEO GSE29981 RAW tar, independently runs
Bioconductor affy::rma in the locked local R environment, and compares the whole
probe-by-GSM matrix to the archived pre-ComBat RMA table.

Original research project is read-only. All scratch/output remains in Downloads.
No ComBat, model fitting, tuning, downloading, or library installation.

This is not yet a full four-series (232-sample) RAW-to-results reproduction.
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
import time
import traceback

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
OUTBASE = Path.home() / 'Downloads' / 'Paper8_repro_stage11A_GSE29981_RMA'
RAW_TAR = ROOT / 'data/p8_cross_series_independence/fogg_original_raw/GSE29981_RAW.tar'
REFERENCE = ROOT / 'data/p8_cross_series_independence/fogg_rma_original_commit/counts/GSE29981_counts_unmapped.tsv'
AUTHOR_R = ROOT / 'data/p8_cross_series_independence/fogg_code/df538583652d3359d71959437837ea3e0d1195e6/data_prep/rma_counts.r'
R_BIN = ROOT / '.envs/p8_fogg_r44/bin/Rscript'
HISTORICAL_PY = ROOT / 'paper8_v350b_fogg_original_rma.py'
PIN_HISTORICAL_PY_SHA256 = 'faf08941b9c098255d36e828de56e2f52930bbb57e3ce5332b9cff509e5e36fe'
SERIES = 'GSE29981'
GSM_RE = re.compile(r'GSM\d+', re.IGNORECASE)
MAX_CEL_UNCOMPRESSED = 80_000_000
MAX_ARCHIVE_MEMBER_COMPRESSED = 80_000_000
MAX_CELS = 200
NUMERIC_TOL = 1e-6
EXPECTED_PROBES = 54675


def say(*args):
    print(*args, flush=True)


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def extract_gsm(name):
    matches = set(x.upper() for x in GSM_RE.findall(name))
    if len(matches) != 1:
        raise RuntimeError('HOLD_GSM_PARSING_' + str(name))
    return next(iter(matches))


def validate_archive_name(name):
    # Tar members are not extracted by the tar library; only their bytes are read.
    # Names remain untrusted and must not determine filesystem paths.
    if '\\' in name or '\x00' in name or '..' in Path(name).parts:
        raise RuntimeError('HOLD_UNSAFE_TAR_NAME')
    base = Path(name).name
    if not (base.lower().endswith('.cel') or base.lower().endswith('.cel.gz')):
        return None
    return extract_gsm(base)


def select_members(raw_tar):
    found = {}
    with tarfile.open(str(raw_tar), mode='r:*') as tar:
        for m in tar:
            gsm = validate_archive_name(m.name)
            if gsm is None:
                continue
            if not m.isfile() or m.size > MAX_ARCHIVE_MEMBER_COMPRESSED:
                raise RuntimeError('HOLD_UNSAFE_OR_OVERSIZE_CEL_MEMBER_' + gsm)
            if gsm in found:
                raise RuntimeError('HOLD_DUPLICATE_GSM_IN_RAW_TAR_' + gsm)
            found[gsm] = {'tar_member': m.name, 'compressed_bytes': m.size}
    if not (5 <= len(found) <= MAX_CELS):
        raise RuntimeError('HOLD_CEL_COUNT_OUTSIDE_BOUNDS_' + str(len(found)))
    return found


def gsm_header(header):
    if len(header) < 2:
        raise RuntimeError('HOLD_EMPTY_RMA_HEADER')
    gsms = [extract_gsm(x) for x in header[1:]]
    if len(set(gsms)) != len(gsms):
        raise RuntimeError('HOLD_DUPLICATE_GSM_HEADER')
    return gsms


def read_tsv_header(path):
    with Path(path).open('r', newline='', encoding='utf-8') as f:
        row = next(csv.reader(f, delimiter='\t'), None)
        if row is None:
            raise RuntimeError('HOLD_EMPTY_REFERENCE_TSV')
        return row


def uncompress_cels(raw_tar, selection, stage_dir):
    stage_dir.mkdir(parents=True, exist_ok=False)
    records = []
    with tarfile.open(str(raw_tar), mode='r:*') as tar:
        for gsm in sorted(selection):
            info = selection[gsm]
            member = tar.getmember(info['tar_member'])
            if not member.isfile():
                raise RuntimeError('HOLD_NONFILE_ARCHIVE_MEMBER_' + gsm)
            dest = stage_dir / (gsm + '.CEL')
            archived_stream = tar.extractfile(member)
            if archived_stream is None:
                raise RuntimeError('HOLD_TAR_MEMBER_NOT_READABLE_' + gsm)
            with archived_stream as source:
                if member.name.lower().endswith('.gz'):
                    stream = gzip.GzipFile(fileobj=source, mode='rb')
                else:
                    stream = source
                digest = hashlib.sha256()
                count = 0
                with dest.open('xb') as output:
                    while True:
                        chunk = stream.read(2 * 1024 * 1024)
                        if not chunk:
                            break
                        count += len(chunk)
                        if count > MAX_CEL_UNCOMPRESSED:
                            raise RuntimeError('HOLD_UNCOMPRESSED_CEL_TOO_LARGE_' + gsm)
                        digest.update(chunk)
                        output.write(chunk)
                if stream is not source:
                    stream.close()
            if count < 100_000:
                raise RuntimeError('HOLD_TRUNCATED_CEL_' + gsm)
            records.append({'gsm': gsm, 'source_tar_member': member.name,
                            'cel_uncompressed_sha256': digest.hexdigest(),
                            'uncompressed_bytes': count})
    return records


R_WORKER = r'''#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop('HOLD_EXPECTED_THREE_INPUTS')
list_file <- args[[1]]
out_file <- args[[2]]
version_file <- args[[3]]
required <- c('affy','Biobase','hgu133plus2cdf')
for (p in required) if (!requireNamespace(p,quietly=TRUE)) stop(paste('HOLD_MISSING_R_PACKAGE',p))
versions <- vapply(required,function(p) as.character(utils::packageVersion(p)),character(1))
if (versions[['affy']] != '1.84.0' || versions[['hgu133plus2cdf']] != '2.18.0' ||
    versions[['Biobase']] != '2.66.0') stop('HOLD_PACKAGE_VERSIONS_NOT_FROZEN')
paths <- readLines(list_file,warn=FALSE)
if (length(paths) < 5 || anyDuplicated(paths) != 0L || any(!file.exists(paths))) stop('HOLD_CEL_LIST')
message('[R_INPUT] ',length(paths),' CEL; full original GEO series')
message('[R_PACKAGES] ',paste(paste(names(versions),versions),collapse='; '))
raw <- affy::ReadAffy(filenames=paths)
message('[R_READAFFY] ',length(Biobase::sampleNames(raw)),' arrays')
processed <- affy::rma(raw)
x <- Biobase::exprs(processed)
if (nrow(x) != 54675L || ncol(x) != length(paths) || any(!is.finite(x))) stop('HOLD_RMA_MATRIX_SHAPE_OR_NUMERICS')
gsms <- regmatches(basename(colnames(x)),regexpr('GSM[0-9]+',basename(colnames(x)),ignore.case=TRUE))
if (length(gsms) != ncol(x) || any(!nzchar(gsms)) || anyDuplicated(gsms) != 0L) stop('HOLD_RMA_GSM_NAMES')
colnames(x) <- toupper(gsms)
probe_ids <- rownames(x)
if (anyDuplicated(probe_ids) != 0L) stop('HOLD_DUPLICATE_PROBE_IDS')
message('[R_RMA_COMPLETE] ',nrow(x),' probes x ',ncol(x),' CEL')
# This table is only a sandbox intermediate for cell-by-cell numeric comparison.
utils::write.table(data.frame(probe_id=probe_ids,x,check.names=FALSE),
                   file=out_file,sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
writeLines(c(paste('R',R.version.string),paste(names(versions),versions,sep='=')),version_file)
message('[R_OUTPUT] ',out_file)
'''


def compare_matrices(archived, rerun, tolerance=NUMERIC_TOL):
    """Compare probe rows by ID and sample columns by GSM, streaming both TSVs.

    A changed probe-row order is a HOLD rather than an implicit reordering.
    Sample column order is explicitly aligned by public GSM identifiers.
    """
    stats = {'status': 'HOLD_NOT_COMPARED', 'n_rows': 0, 'n_samples': 0,
             'values_compared': 0, 'max_abs_error': 0., 'sum_abs_error': 0.,
             'n_abs_gt_1e_6': 0, 'n_abs_gt_1e_8': 0,
             'probe_mismatch_count': 0, 'worst_probe_id': None,
             'worst_gsm': None, 'first_mismatches': []}
    with Path(archived).open('r', newline='', encoding='utf-8') as a, Path(rerun).open('r',newline='',encoding='utf-8') as b:
        aa = csv.reader(a, delimiter='\t')
        bb = csv.reader(b, delimiter='\t')
        a_head, b_head = next(aa), next(bb)
        a_gsms, b_gsms = gsm_header(a_head), gsm_header(b_head)
        if set(a_gsms) != set(b_gsms):
            missing = sorted(set(a_gsms).symmetric_difference(set(b_gsms)))
            raise RuntimeError('HOLD_RMA_GSM_ROSTER_DIFFERENCE_' + str(missing[:8]))
        col_b = {gsm: 1 + b_gsms.index(gsm) for gsm in a_gsms}
        cols = [col_b[x] for x in a_gsms]
        stats['n_samples'] = len(a_gsms)
        for ridx, (row_a, row_b) in enumerate(__import__('itertools').zip_longest(aa, bb), start=1):
            if row_a is None or row_b is None:
                raise RuntimeError('HOLD_RMA_ROW_COUNTS_DIFFER')
            if len(row_a) != len(a_head) or len(row_b) != len(b_head):
                raise RuntimeError('HOLD_RMA_ROW_WIDTH_DIFFER_' + str(ridx))
            if row_a[0] != row_b[0]:
                raise RuntimeError('HOLD_RMA_PROBE_ORDER_DIFFER_' + str(ridx) + '_' + row_a[0] + '_' + row_b[0])
            for j,(gsm,idx) in enumerate(zip(a_gsms,cols),start=1):
                av,bv = float(row_a[j]),float(row_b[idx])
                if not math.isfinite(av) or not math.isfinite(bv):
                    raise RuntimeError('HOLD_NONFINITE_RMA_VALUE')
                diff=abs(av-bv)
                stats['values_compared'] += 1
                stats['sum_abs_error'] += diff
                if diff > stats['max_abs_error']:
                    stats['max_abs_error'] = diff
                    stats['worst_probe_id'] = row_a[0]
                    stats['worst_gsm'] = gsm
                if diff>1e-8:stats['n_abs_gt_1e_8']+=1
                if diff>1e-6:
                    stats['n_abs_gt_1e_6']+=1
                    if len(stats['first_mismatches'])<12:
                        stats['first_mismatches'].append({'probe':row_a[0], 'gsm':gsm,
                                                          'archived':av, 'rerun':bv,'abs_diff':diff})
            stats['n_rows']+=1
            if (ridx % 10000)==0:say('[COMPARE_RMA]',ridx,'rows')
    if stats['n_rows'] != EXPECTED_PROBES:
        raise RuntimeError('HOLD_UNEXPECTED_PROBE_COUNT_'+str(stats['n_rows']))
    stats['mean_abs_error']=stats.pop('sum_abs_error')/stats['values_compared']
    stats['status']='PASS_INDEPENDENT_RMA_GSE29981_NUMERIC_REPLAY_1E6' if stats['n_abs_gt_1e_6']==0 else 'HOLD_RMA_NUMERIC_DIFFERENCE'
    return stats


def selftest():
    assert extract_gsm('a_GSM742055_something.CEL.gz')=='GSM742055'
    assert validate_archive_name('nested/GSM742055.CEL.gz')=='GSM742055'
    assert validate_archive_name('README.txt') is None
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as p:
        a=Path(p)/'a.tsv'; b=Path(p)/'b.tsv'
        a.write_text('probe_id\tGSM1\tGSM2\np1\t1\t2\np2\t3\t4\n',encoding='utf-8')
        b.write_text('probe_id\tGSM2\tGSM1\np1\t2\t1\np2\t4\t3\n',encoding='utf-8')
        # Miniature test uses the same comparator with an adjusted expected row count.
        global EXPECTED_PROBES
        prior=EXPECTED_PROBES;EXPECTED_PROBES=2
        try:
            assert compare_matrices(a,b)['status'].startswith('PASS_')
            b.write_text('probe_id\tGSM2\tGSM1\np1\t2\t1.05\np2\t4\t3\n',encoding='utf-8')
            assert compare_matrices(a,b)['status']=='HOLD_RMA_NUMERIC_DIFFERENCE'
        finally:EXPECTED_PROBES=prior
    return True


def main():
    say('PAPER 8 | STAGE 11A | INDEPENDENT GSE29981 RAW CEL -> RMA REPLAY')
    say('[POLICY] ORIGINAL PROJECT READ ONLY | NO DOWNLOADING | NO COMBAT/FITS | OUTPUT DOWNLOADS')
    assert selftest()
    say('[SELFTEST] PASS_MATRIX_COMPARISON_AND_GSM_ALIGNMENT')
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run=OUTBASE/('Run_'+stamp+'_'+str(os.getpid()))
    run.mkdir(parents=True,exist_ok=False)
    report={'stage':'P8_STAGE11A_RMA_GSE29981','scope':'PER_SERIES_GSE29981_ONLY',
            'status':'HOLD_NOT_STARTED','project_read_only':True,'new_downloads':False,
            'models_fitted':False,'ComBat_executed':False,
            'independent_frozen_method_not_original_author_script':True,
            'parameters':{'method':'affy::ReadAffy; affy::rma; Biobase::exprs',
                          'output_tolerance':NUMERIC_TOL, 'expected_probes':EXPECTED_PROBES},
            'raw_input_checksums':{},'warnings':[]}
    start=time.monotonic()
    try:
        # Avoid unpacking a GEO cohort when Downloads has little free space.
        free_bytes=shutil.disk_usage(str(run)).free
        report['downloads_free_bytes_before']=free_bytes
        if free_bytes < 2_000_000_000:
            raise RuntimeError('HOLD_LOW_DOWNLOADS_FREE_SPACE')
        for name,p in [('GEO_RAW_tar',RAW_TAR),('archived_RMA_matrix',REFERENCE),
                       ('original_author_rma_script',AUTHOR_R),('frozen_v350b_python',HISTORICAL_PY),
                       ('locked_Rscript',R_BIN)]:
            if not p.is_file():raise RuntimeError('HOLD_REQUIRED_INPUT_MISSING_'+name)
            report['raw_input_checksums'][name]=file_sha256(p)
            say('[INPUT_HASH]',name,report['raw_input_checksums'][name][:16])
        if report['raw_input_checksums']['frozen_v350b_python'] != PIN_HISTORICAL_PY_SHA256:
            raise RuntimeError('HOLD_HISTORICAL_V350B_SCRIPT_SHA_CHANGED')
        source_code=AUTHOR_R.read_text(encoding='utf-8',errors='replace')
        if not re.search(r'\brma\s*\(',source_code,re.IGNORECASE):
            raise RuntimeError('HOLD_AUTHOR_RMA_METHOD_NOT_VERIFIED')
        report['original_R_source_only_verified_not_executed']=True
        ref_head=read_tsv_header(REFERENCE)
        ref_gsms=gsm_header(ref_head)
        say('[REFERENCE]',len(ref_gsms),'GSM; archive header',ref_head[0])
        selection=select_members(RAW_TAR)
        report['raw_tar_CEL_count']=len(selection)
        say('[GEO_TAR]',len(selection),'CEL members')
        if set(selection)!=set(ref_gsms):
            differences=sorted(set(selection).symmetric_difference(set(ref_gsms)))
            raise RuntimeError('HOLD_SOURCE_REFERENCE_GSM_DIFFER_'+str(differences[:15]))
        say('[ROSTER]',len(ref_gsms),'GSM exact match')
        cel_dir=run/'temporary_uncompressed_CEL'
        recs=uncompress_cels(RAW_TAR,selection,cel_dir)
        (run/'raw_CEL_manifest.json').write_text(json.dumps(recs,indent=2)+'\n',encoding='utf-8')
        say('[EXTRACT]',len(recs),'CEL into sandbox only')
        order=[str(cel_dir/(gsm+'.CEL')) for gsm in sorted(selection)]
        list_file=run/'cel_list.txt'
        list_file.write_text('\n'.join(order)+'\n',encoding='utf-8')
        worker=run/'independent_rma_worker.R'
        worker.write_text(R_WORKER,encoding='utf-8')
        target=run/'GSE29981_independent_RMA.tsv'
        version_file=run/'R_versions.txt'
        cmd=[str(R_BIN),str(worker),str(list_file),str(target),str(version_file)]
        say('[R_START] full GSE29981 series; independent per-series affy::rma')
        with (run/'R_stdout.log').open('w',encoding='utf-8') as out, (run/'R_stderr.log').open('w',encoding='utf-8') as err:
            r=subprocess.run(cmd,stdout=out,stderr=err,cwd=str(run),check=False)
        report['R_exit_code']=r.returncode
        if r.returncode!=0:
            message=(run/'R_stderr.log').read_text(encoding='utf-8',errors='replace')
            say('[R_ERROR_TAIL]',message[-1500:])
            raise RuntimeError('HOLD_R_EXECUTION_EXIT_'+str(r.returncode))
        if not target.is_file():raise RuntimeError('HOLD_NO_NEW_RMA_OUTPUT')
        report['numeric_comparison']=compare_matrices(REFERENCE,target)
        report['status']=report['numeric_comparison']['status']
        metrics=report['numeric_comparison']
        say('[RMA_COMPARE]',metrics['n_rows'],'probes x',metrics['n_samples'],'GSM')
        say('[RMA_MAX_ABS]',metrics['max_abs_error'])
        say('[RMA_GT_1E6]',metrics['n_abs_gt_1e_6'])
        if metrics['status'].startswith('PASS_'):
            say('[RESULT] Independent RMA numerically matches archived GSE29981')
        else:
            say('[RESULT] HOLD_RMA_NUMERIC_DIFFERENCE; see first_mismatches')
    except Exception as exc:
        report['status']='HOLD_STAGE11A_'+type(exc).__name__
        report['error']=str(exc)
        (run/'traceback.log').write_text(traceback.format_exc(),encoding='utf-8')
        say('[ERROR]',str(exc))
    finally:
        # This directory is created by this run only; never delete project data.
        try:
            temp=run/'temporary_uncompressed_CEL'
            if temp.is_dir():
                shutil.rmtree(temp)
                report['scratch_CEL_cleaned_up']=True
        except Exception as exc:
            report['warnings'].append('scratch_cleanup_failed: '+str(exc))
        report['elapsed_seconds']=round(time.monotonic()-start,3)
        report['limitations']=[
            'Only GSE29981 among the four GEO series is replayed, not all 232 arrays.',
            'Input tar is the archived public GEO RAW source; no fresh GEO download is performed.',
            'The independently coded standard affy::rma pipeline is used; original author R script is hash-recorded but not executed.',
            'No classifier/ComBat or clinical validation is executed.',
            'Separate source-only single-CEL MAS5 is not part of this Stage 11A.'
        ]
        (run/'stage11a_replay.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        (run/'summary.txt').write_text('PAPER 8 Stage 11A\nStatus: '+report['status']+'\nScope: independent GSE29981 per-series RMA\n',encoding='utf-8')
        say('[STATUS]',report['status'])
        say('[OUTPUT]',run)
    return 0 if report['status'].startswith('PASS_') else 2

if __name__=='__main__':
    sys.exit(main())
