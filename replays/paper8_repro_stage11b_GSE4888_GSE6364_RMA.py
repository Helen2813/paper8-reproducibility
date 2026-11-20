#!/usr/bin/env python3
"""PAPER 8 Stage 11B: raw-CEL -> RMA independent numeric replay.

Runs GSE4888 and GSE6364 one series at a time against the archived
unmapped RMA matrices, keeping all temporary output in Downloads.

No network, installs, batch correction, tuning, model training, or
writes to the original paper8/endo-research repository.

GSE51981 (148 CEL) is deliberately excluded from this bounded step;
run it later with a separate memory/storage preflight.

The script embeds the tested Stage 11A CEL/RMA comparator instead of
requiring another Downloads script, making the script self-contained.
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
import zipfile

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
OUTBASE = Path.home() / 'Downloads' / 'Paper8_repro_stage11B_two_series_RMA'
AUTHOR_R = ROOT / 'data/p8_cross_series_independence/fogg_code/df538583652d3359d71959437837ea3e0d1195e6/data_prep/rma_counts.r'
R_BIN = ROOT / '.envs/p8_fogg_r44/bin/Rscript'
HISTORICAL_PY = ROOT / 'paper8_v350b_fogg_original_rma.py'
PIN_HISTORICAL_PY_SHA256 = 'faf08941b9c098255d36e828de56e2f52930bbb57e3ce5332b9cff509e5e36fe'
PIN_AUTHOR_R_SHA256 = '78b360703c1c1936'  # first 16 hexadecimal digits, verified by Stage 11A
REF_BASE = ROOT / 'data/p8_cross_series_independence/fogg_rma_original_commit/counts'
HISTORICAL_REPORT_ZIP = ROOT / 'results/p8_v350b_original_rma_20261108T014606Z.zip'
SERIES_CONFIG = (
    ('GSE4888', [ROOT / 'data/p8_cross_series_independence/fogg_original_raw/GSE4888_RAW.tar']),
    ('GSE6364', [ROOT / 'data/p8_cross_series_independence/raw_archives/GSE6364_RAW.tar',
                 ROOT / 'data/external/p8_v347f_gse6364_cel/GSE6364_RAW.tar']),
)
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
    stats['status']='PASS_INDEPENDENT_RMA_SERIES_NUMERIC_REPLAY_1E6' if stats['n_abs_gt_1e_6']==0 else 'HOLD_RMA_NUMERIC_DIFFERENCE'
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



def load_frozen_rma_manifest(path):
    if not path.is_file():
        raise RuntimeError('HOLD_HISTORICAL_RMA_MANIFEST_ZIP_MISSING')
    with zipfile.ZipFile(str(path)) as archive:
        bad = archive.testzip()
        if bad:
            raise RuntimeError('HOLD_RMA_MANIFEST_ZIP_BAD_MEMBER_' + bad)
        names = [n for n in archive.namelist() if Path(n).name == 'rma_report.json']
        if len(names) != 1:
            raise RuntimeError('HOLD_RMA_MANIFEST_NOT_UNIQUE_' + str(len(names)))
        rep = json.loads(archive.read(names[0]).decode('utf-8'))
    if rep.get('status') != 'PASS_ORIGINAL_FOUR_SERIES_RMA':
        raise RuntimeError('HOLD_HISTORICAL_RMA_MANIFEST_STATUS_' + str(rep.get('status')))
    if not isinstance(rep.get('cohorts'), dict):
        raise RuntimeError('HOLD_HISTORICAL_RMA_MANIFEST_COHORTS_SCHEMA')
    return rep


def choose_verified_tar(candidates, record):
    candidates = [p for p in candidates if p.is_file()]
    if not candidates:
        raise RuntimeError('HOLD_RAW_GEO_TAR_NOT_FOUND')
    hashes = [{'relative_path':str(x.relative_to(ROOT)), 'sha256':file_sha256(x),
               'size_bytes':x.stat().st_size} for x in candidates]
    record['tar_candidates'] = hashes
    if len({h['sha256'] for h in hashes}) > 1:
        raise RuntimeError('HOLD_MULTIPLE_RAW_TAR_COPIES_DISAGREE')
    return candidates[0]


def rma_series(series, paths, historical, output):
    output.mkdir(parents=True, exist_ok=False)
    record = {'series':series, 'status':'HOLD_NOT_STARTED',
              'input_checksums':{}, 'original_project_modified':False,
              'rma_numeric_replay':False, 'warnings':[]}
    try:
        free = shutil.disk_usage(str(output)).free
        record['downloads_free_bytes_start'] = free
        if free < 3_000_000_000:
            raise RuntimeError('HOLD_LESS_THAN_3GB_FREE_SPACE_IN_DOWNLOADS')
        raw_tar = choose_verified_tar(paths, record)
        ref = REF_BASE / (series + '_counts_unmapped.tsv')
        if not ref.is_file():
            raise RuntimeError('HOLD_ARCHIVED_RMA_REFERENCE_NOT_FOUND')
        ref_sha = file_sha256(ref)
        record['input_checksums']['archived_RMA_matrix'] = ref_sha
        archived_series = historical['cohorts'].get(series)
        if not isinstance(archived_series, dict):
            raise RuntimeError('HOLD_SERIES_ABSENT_FROM_HISTORICAL_RMA_REPORT')
        record['archived_manifest_cohort'] = archived_series
        expected_ref = archived_series.get('output_sha256')
        if not expected_ref or expected_ref != ref_sha:
            raise RuntimeError('HOLD_FROZEN_RMA_REFERENCE_SHA_MISMATCH')
        say('[SOURCE_HASH]',series,'raw=',record['tar_candidates'][0]['sha256'][:16],
            'archived_RMA=',ref_sha[:16],'MATCH')
        ref_head=read_tsv_header(ref)
        ref_gsms=gsm_header(ref_head)
        selection=select_members(raw_tar)
        record['n_archived_GSM']=len(ref_gsms)
        record['n_raw_CEL']=len(selection)
        say('[GEO_TAR]',series,'CEL=',len(selection),'RMA_GSM=',len(ref_gsms))
        if set(selection) != set(ref_gsms):
            missing=sorted(set(selection).symmetric_difference(set(ref_gsms)))
            raise RuntimeError('HOLD_RAW_VS_REFERENCE_GSM_DIFFER_'+str(missing[:10]))
        cel_dir=output/'scratch_uncompressed_CEL'
        cels=uncompress_cels(raw_tar,selection,cel_dir)
        (output/'uncompressed_CEL_SHA256.json').write_text(json.dumps(cels,indent=2)+'\n',encoding='utf-8')
        paths_in_order=[str(cel_dir/(g+'.CEL')) for g in sorted(selection)]
        listing=output/'cel_list.txt'
        listing.write_text('\n'.join(paths_in_order)+'\n',encoding='utf-8')
        worker=output/'independent_per_series_rma.R'
        worker.write_text(R_WORKER,encoding='utf-8')
        record['r_worker_sha256']=file_sha256(worker)
        out_mat=output/(series+'_independent_RMA.tsv')
        cmd=[str(R_BIN),str(worker),str(listing),str(out_mat),str(output/'R_versions.txt')]
        say('[R_START]',series,'affy::rma with',len(paths_in_order),'CEL')
        with (output/'R_stdout.log').open('w',encoding='utf-8') as stdout, (output/'R_stderr.log').open('w',encoding='utf-8') as stderr:
            run=subprocess.run(cmd,stdout=stdout,stderr=stderr,cwd=str(output),check=False)
        record['R_exit_code']=run.returncode
        if run.returncode != 0:
            tail=(output/'R_stderr.log').read_text(encoding='utf-8',errors='replace')[-1400:]
            say('[R_STDERR_TAIL]',tail)
            raise RuntimeError('HOLD_R_EXIT_'+str(run.returncode))
        if not out_mat.is_file():
            raise RuntimeError('HOLD_RMA_REPLAY_OUTPUT_NOT_CREATED')
        match=compare_matrices(ref,out_mat)
        record['numeric_comparison']=match
        record['independent_rma_output_sha256']=file_sha256(out_mat)
        record['status']=match['status'] if match['status'] != 'PASS_INDEPENDENT_RMA_SERIES_NUMERIC_REPLAY_1E6' else 'PASS_INDEPENDENT_RMA_'+series+'_NUMERIC_REPLAY_1E6'
        record['rma_numeric_replay']=True
        say('[RMA_COMPARE]',series,match['n_rows'],'probes x',match['n_samples'],'GSM',
            'max_abs=',match['max_abs_error'],'n_gt_1e-6=',match['n_abs_gt_1e_6'])
    except Exception as exc:
        record['status']='HOLD_'+series+'_'+type(exc).__name__
        record['error']=str(exc)
        (output/'traceback.log').write_text(traceback.format_exc(),encoding='utf-8')
        say('[SERIES_ERROR]',series,str(exc))
    finally:
        scratch=output/'scratch_uncompressed_CEL'
        try:
            if scratch.exists():
                shutil.rmtree(scratch)
                record['scratch_CEL_removed']=True
        except Exception as exc:
            record['warnings'].append('SCRATCH_CEL_CLEANUP_FAILURE:'+str(exc))
        (output/'series_replay.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        say('[SERIES_STATUS]',series,record['status'])
    return record


def main():
    say('PAPER 8 | STAGE 11B | GSE4888 + GSE6364 INDEPENDENT RAW CEL -> RMA REPLAY')
    say('[POLICY] ORIGINAL PROJECT READ ONLY | GEO SOURCES OFFLINE | NO MODEL FITS/COMBAT')
    assert selftest()
    say('[SELFTEST] PASS_CEL_AND_STREAMING_NUMERIC_COMPARATOR')
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run=OUTBASE/('Run_'+stamp+'_'+str(os.getpid()))
    run.mkdir(parents=True,exist_ok=False)
    report={'stage':'P8_STAGE11B_TWO_SERIES_RAW_RMA','status':'HOLD_NOT_STARTED',
            'scope':['GSE4888','GSE6364'], 'preserved_project':True,
            'no_downloads':True,'no_combat':True,'no_classifier_fits':True,
            'n_probes_expected':EXPECTED_PROBES,'tolerance':NUMERIC_TOL,
            'steps':[],'limitations':['This stage repeats only GSE4888 and GSE6364.',
            'GSE51981 is deliberately excluded for a separate large-series preflight.',
            'Comparisons are to frozen RMA matrices, not independent clinical validation.',
            'Historic author R code verified but independently written affy::rma worker executed.']}
    try:
        for label,path in [('frozen_v350b_Python',HISTORICAL_PY),
                           ('author_rma_counts_r',AUTHOR_R),('Rscript',R_BIN),
                           ('rma_report_zip',HISTORICAL_REPORT_ZIP)]:
            if not path.is_file():raise RuntimeError('HOLD_REQUIRED_SOURCE_MISSING_'+label)
            report[label+'_sha256']=file_sha256(path)
            say('[SOURCE_SHA]',label,report[label+'_sha256'][:16])
        if report['frozen_v350b_Python_sha256'] != PIN_HISTORICAL_PY_SHA256:
            raise RuntimeError('HOLD_ORIGINAL_V350B_SCRIPT_HASH_CHANGED')
        if not report['author_rma_counts_r_sha256'].startswith(PIN_AUTHOR_R_SHA256):
            raise RuntimeError('HOLD_ORIGINAL_AUTHOR_R_CODE_HASH_CHANGED')
        if not re.search(r'\brma\s*\(',AUTHOR_R.read_text(encoding='utf-8',errors='replace'),re.I):
            raise RuntimeError('HOLD_RMA_FUNCTION_NOT_IN_ORIGINAL_SOURCE')
        historical=load_frozen_rma_manifest(HISTORICAL_REPORT_ZIP)
        report['historical_report_status']=historical['status']
        report['historical_report_zip_sha256']=file_sha256(HISTORICAL_REPORT_ZIP)
        for series,tars in SERIES_CONFIG:
            item=rma_series(series,tars,historical,run/series)
            report['steps'].append({'series':series, 'status':item['status'],
                                    'numeric_comparison':item.get('numeric_comparison'),
                                    'GSM':item.get('n_raw_CEL')})
        passed=[x for x in report['steps'] if x['status'].startswith('PASS_')]
        report['status']='PASS_INDEPENDENT_RMA_TWO_SERIES_1E6' if len(passed)==2 else 'HOLD_PARTIAL_RMA_REPLAY_'+str(len(passed))+'_OF_2'
    except Exception as exc:
        report['status']='HOLD_STAGE11B_'+type(exc).__name__
        report['error']=str(exc)
        (run/'traceback.log').write_text(traceback.format_exc(),encoding='utf-8')
        say('[ERROR]',str(exc))
    finally:
        (run/'stage11b_replay.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        (run/'summary.txt').write_text('Paper 8 Stage 11B\nStatus: '+report['status']+'\n'+'\n'.join(
           x['series']+': '+x['status'] for x in report['steps'])+'\n',encoding='utf-8')
        say('[STATUS]',report['status'])
        say('[OUTPUT]',run)
    return 0 if report['status'].startswith('PASS_') else 2

if __name__ == '__main__':
    sys.exit(main())
