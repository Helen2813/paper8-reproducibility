#!/usr/bin/env python3
"""Paper 8 / Stage 8D v3: bounded independent null-refit pilot.

Run from macOS Terminal:
    python3 "$HOME/Downloads/paper8_repro_stage8d_numeric_pilot_v3.py"

Project files are READ ONLY; results are written ONLY under ~/Downloads.
No downloads, installs, test tuning, random-label regeneration, or source modifications.

Replays 3 RMA-only replicates per family (18 model fits), and 1 ComBat replicate
per family under the outcome-informed and no-outcome arms (12 model fits).
The original v3.71 Python driver was not recovered: this is a FIXED-HYPOTHESIS
independent reconstruction, not a claim of exact source-code re-execution.
All numerical mismatches are retained as HOLD, never tuned away.

IMPORTANT: DEDUP153 was independently permuted WITHIN BATCH after removing
26 training arrays. It is NOT the FULL179 permutation restricted to 153 GSMs.
DEDUP153 files have unresolved v3.71 file SHA pins, so an in-range and
stratum-conserving assignment proves protocol compatibility, NOT historical
permutation identity; any pilot match is explicitly conditional.
"""
from __future__ import annotations

import collections
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback
import warnings
import zipfile

STAGE = 'P8_STAGE8D_V3_LAYOUT_SPECIFIC_ASSIGNMENT_GATE'
ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
DOWNLOADS = Path.home() / 'Downloads'
OUTROOT = DOWNLOADS / 'Paper8_repro_stage8D_numeric_pilot'
V371 = ROOT / 'results/p8_v371_full_gene_combat_null_RESUMABLE'
V366 = ROOT / 'results/p8_v366_batch_split_recovery_20261109T002942Z'
V367 = ROOT / 'results/p8_v367_retrained_nulls_20261109T005048Z'
V369 = ROOT / 'results/p8_v369_core270_corrected_lineage_20261109T015634Z'
V352 = ROOT / 'results/p8_v352g_factorial_20261108T042414Z'
ORIG = ROOT / 'data/p8_cross_series_independence/fogg_original_combat/20261108T015502Z'
R_BIN = ROOT / '.envs/p8_fogg_r44/bin/Rscript'
PY_BIN = ROOT / '.envs/p8_null/bin/python3'

EXPECTED_PROTOCOL_SHA = 'b0b40894a9181e2e65b9afdd98f6fe658a8c0c77e362bb9f9d0d82c2f7286ac1'
EXPECTED_R_HELPER_SHA = 'c352cb9980daa3904a9c8a79f6f0027b8ccc1ca0a58223b08fce66b622f95705'
EXPECTED_RMA_CORE_SHA = 'e3b04bbe5f4f9fc163e8c69e697b3c292d3c984d15fe0f02b2b55e581aaa551c'
EXPECTED_V369_ARCHIVED_SHA = 'c96f2ec13105fb3d015f844ed874fe458060627bd8db1d28b2e437a429ce1a72'
EXPECTED_V371_SHA = '22fcef9c8fb1124f293cae634756eacbe7f35940b26b314ed333443abf35109c'
CONFIGS = {
    'ARCHIVED': (0.6552218193440044, 0.0675608348364084),
    'PRIOR_FIXED_SENSITIVITY': (3.0, 0.8),
}
# Convergence setting is an independent fixed reconstruction, NOT verified
# against the unavailable original null driver; no search over max_iter.
RECONSTRUCTED_ITERATIONS = 5000
RMA_REPLICATES = (1, 2, 3)
COMBAT_REPLICATE = 1
FAMILIES = (
    ('FULL179', 'WITHIN_BATCH'),
    ('FULL179', 'WITHIN_BATCH_PHASE'),
    ('DEDUP153', 'WITHIN_BATCH'),
)
EXPECTED_ASSIGNMENT_FILES = {
    ('FULL179','WITHIN_BATCH'): 'assign_FULL179_WITHIN_BATCH.csv',
    ('FULL179','WITHIN_BATCH_PHASE'): 'assign_FULL179_WITHIN_BATCH_PHASE.csv',
    ('DEDUP153','WITHIN_BATCH'): 'assign_DEDUP153_WITHIN_BATCH.csv',
}


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def save_json(path,obj):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def read_dicts(path, delim=','):
    with Path(path).open('r',newline='',encoding='utf-8-sig') as f:
        return list(csv.DictReader(f,delimiter=delim))


def check_sha(path, expected, name, manifest):
    path=Path(path)
    if not path.is_file() or path.is_symlink():
        raise RuntimeError('HOLD_MISSING_OR_SYMLINK_'+name+' '+str(path))
    real=digest(path)
    manifest[name]={'path':str(path),'sha256':real,'expected_sha256':expected,'match':real==expected,'bytes':path.stat().st_size}
    if real!=expected:
        raise RuntimeError('HOLD_SHA256_MISMATCH_'+name)
    print('[HASH]',name,'MATCH',flush=True)
    return real


def canonical_research_file(raw):
    marker='/endo-research/'
    if marker not in str(raw):
        raise RuntimeError('HOLD_PINNED_PATH_OUTSIDE_STUDY_ROOT')
    r=Path(str(raw).split(marker,1)[1])
    if r.is_absolute() or '..' in r.parts:
        raise RuntimeError('HOLD_UNSAFE_PINNED_PATH')
    dest=ROOT/r
    if not dest.resolve().is_relative_to(ROOT.resolve()):
        raise RuntimeError('HOLD_PINNED_SYMLINK_ESCAPE')
    return dest


def scan_downloads_for_pinned_sources(report):
    """Optional local provenance: never execute historically located scripts."""
    hits={}
    scans=0
    for name,pin in [('v369',EXPECTED_V369_ARCHIVED_SHA),('v371',EXPECTED_V371_SHA)]:
        hits[name]={'found':False,'expected_sha256':pin}
    for p in sorted(DOWNLOADS.glob('*.py')):
        if scans>=600: break
        if not p.is_file() or p.is_symlink() or p.stat().st_size>1_000_000: continue
        scans+=1
        if not re.search(r'(?:p8|paper8|v3(?:69|70|71))',p.name,re.I): continue
        h=digest(p)
        for name,pin in [('v369',EXPECTED_V369_ARCHIVED_SHA),('v371',EXPECTED_V371_SHA)]:
            if h==pin: hits[name]={'found':True,'path':str(p),'expected_sha256':pin}
    report['exact_historical_drivers']=hits
    report['searched_downloads_scripts_bounded']=scans
    print('[EXACT_DRIVERS_FOUND]',{k:v['found'] for k,v in hits.items()},flush=True)


def validate_sources(run_dir, report):
    print('PAPER 8 | STAGE 8D | INDEPENDENT BOUNDED NULL NUMERIC PILOT',flush=True)
    print('[POLICY] ORIGINAL PROJECT READ ONLY | OUTPUTS DOWNLOADS ONLY | NO TUNING',flush=True)
    if not ROOT.is_dir(): raise RuntimeError('HOLD_PROJECT_ROOT_NOT_FOUND')
    prot_path=V371/'protocol_lock.json'
    check_sha(prot_path,EXPECTED_PROTOCOL_SHA,'protocol_v371',report['hashes'])
    proto=json.loads(prot_path.read_text(encoding='utf-8'))
    if proto.get('v371_script_sha256')!=EXPECTED_V371_SHA:
        raise RuntimeError('HOLD_ORIGINAL_DRIVER_PIN_DIFFERS')
    if len(proto.get('permutation_sha256',{}))!=3:
        raise RuntimeError('HOLD_PERMUTATION_PIN_COUNT')
    check_sha(V371/'fullgene_combat_helper.R',EXPECTED_R_HELPER_SHA,'original_R_helper',report['hashes'])
    check_sha(V352/'FULL216__NO_COMBAT.tsv',EXPECTED_RMA_CORE_SHA,'RMA_only_CORE270',report['hashes'])
    pinned=proto['input_hashes']
    for raw,exp in pinned.items():
        key='pinned_'+Path(raw).name
        # duplicate basenames may occur; preserve context and avoid overwrites
        rel=str(canonical_research_file(raw).relative_to(ROOT))
        key='pinned_'+rel.replace('/','__')
        check_sha(canonical_research_file(raw),exp,key,report['hashes'])
    check_sha(V366/'12_locked_training_null_assignments_NOT_FITTED.csv',
             'b3a741d39b577a200106eb71a5c753caa34fcf0be4c4e2041f7d2f3571bb567c',
             'v366_full179_null_registry',report['hashes'])
    rost=read_dicts(V366/'08c_verified_216_source_batch_roster.csv')
    if len(rost)!=216 or len({r['gsm'] for r in rost})!=216 or collections.Counter(r['split'] for r in rost)!={'TRAIN179':179,'TEST37':37}:
        raise RuntimeError('HOLD_V366_ROSTER_SCHEMA_OR_COUNTS')
    train_ids={r['gsm'] for r in rost if r['split']=='TRAIN179'}
    test_ids={r['gsm'] for r in rost if r['split']=='TEST37'}
    dropped=read_dicts(V367/'06_dedup_crosswalk_train_removed.csv')
    removed={r['training_removed_gsm'] for r in dropped}
    if len(dropped)!=26 or len(removed)!=26 or not removed<=train_ids:
        raise RuntimeError('HOLD_DEDUP_TRAIN_CROSSWALK')
    # Byte hashes of the v3.71 assignment CSVs have not yet been validated.
    # Do NOT treat their recorded protocol pins as file-content hashes without
    # independently checking every row against the pinned v3.66 source registry.
    # A byte mismatch is recorded prominently; a LABEL mismatch is fatal.
    legacy=collections.defaultdict(dict)
    frozen_source=V366/'12_locked_training_null_assignments_NOT_FITTED.csv'
    roster_by_gsm={r['gsm']:r for r in rost}
    phase_rows=read_dicts(V366/'08_verified_179_training_roster.csv')
    phase_by_gsm={r['gsm']:r['phase'] for r in phase_rows}
    if len(phase_by_gsm)!=179 or set(phase_by_gsm)!=train_ids:
        raise RuntimeError('HOLD_V366_PHASE_ROSTER_INVALID')
    reg_count=0
    for r in read_dicts(frozen_source):
        if r['training_layout']!='FULL179_ONLY' or r['family'] not in ('WITHIN_BATCH','WITHIN_BATCH_PHASE'):
            raise RuntimeError('HOLD_LOCKED_REGISTRY_UNEXPECTED_LAYOUT_OR_FAMILY')
        fam=r['family']; rep=int(r['replicate']); gsm=r['gsm']
        if rep not in range(1,200) or gsm not in train_ids:
            raise RuntimeError('HOLD_LOCKED_REGISTRY_BAD_REP_OR_GSM')
        yy=int(r['permuted_train_label'])
        if yy not in (0,1) or int(r['original_truth'])!=int(roster_by_gsm[gsm]['truth']) or (str(r['batch'])!=str(roster_by_gsm[gsm]['batch'])) or r['phase']!=phase_by_gsm[gsm]:
            raise RuntimeError('HOLD_LOCKED_REGISTRY_ROSTER_OR_TRUTH_MISMATCH')
        k=(fam,rep)
        if gsm in legacy[k]:raise RuntimeError('HOLD_LOCKED_REGISTRY_DUPLICATE')
        legacy[k][gsm]=yy
        reg_count+=1
    if reg_count!=71242 or len(legacy)!=398 or any(set(v)!=train_ids for v in legacy.values()):
        raise RuntimeError('HOLD_LOCKED_REGISTRY_COVERAGE_NOT_71242')
    # Preserve group-level case counts for the frozen FULL179 permutations.
    # The DEDUP153 layout is a post-hoc restriction; its group counts need not
    # remain identical after the 26 removed arrays are excluded.
    for fam in ('WITHIN_BATCH','WITHIN_BATCH_PHASE'):
        grouped=collections.defaultdict(list)
        for g in sorted(train_ids):
            meta=roster_by_gsm[g]
            stratum=(meta['batch'],phase_by_gsm[g]) if fam=='WITHIN_BATCH_PHASE' else (meta['batch'],)
            grouped[stratum].append(g)
        for rep in range(1,200):
            perm=legacy[(fam,rep)]
            for ids in grouped.values():
                if sum(perm[g] for g in ids)!=sum(int(roster_by_gsm[g]['truth']) for g in ids):
                    raise RuntimeError('HOLD_V366_NULL_STRATUM_BALANCE_'+fam+'_REP'+str(rep))

    ass_data={}
    evidence=[]
    blockers=[]
    for layout,family in FAMILIES:
        path=V371/EXPECTED_ASSIGNMENT_FILES[(layout,family)]
        if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(ROOT.resolve()):
            raise RuntimeError('HOLD_MISSING_UNSAFE_ASSIGNMENT_'+layout+'_'+family)
        expected=proto['permutation_sha256'][layout+'/'+family]
        actual=digest(path)
        use_train=train_ids if layout=='FULL179' else train_ids-removed
        byrep=collections.defaultdict(dict)
        rows=read_dicts(path)
        need={'replicate','gsm','permuted_label'}
        if not rows or not need<=set(rows[0]):
            raise RuntimeError('HOLD_ASSIGNMENT_REQUIRED_COLUMNS_'+layout+'_'+family)
        failures=[]
        for r in rows:
            try:
                rep=int(r['replicate']);gsm=r['gsm'];y=int(r['permuted_label'])
            except (KeyError,TypeError,ValueError):
                failures.append({'failure':'unparseable_row','data':str(r)[:150]})
                if len(failures)>=10:break
                continue
            if rep not in range(1,200) or gsm not in use_train or y not in (0,1) or gsm in byrep[rep]:
                failures.append({'failure':'invalid_or_duplicate','replicate':rep,'gsm':gsm})
                if len(failures)>=10:break
            else:
                byrep[rep][gsm]=y
        coverage=(len(rows)==199*len(use_train) and len(byrep)==199
                  and all(set(byrep[r])==use_train for r in range(1,200)))
        n_distinct=(len({tuple(byrep[r][g] for g in sorted(use_train)) for r in range(1,200)})
                    if coverage else 0)
        # Strict compositional test on the ACTUAL training population for each
        # layout. DEDUP153 must preserve the true labels within its 153-sample
        # batch strata; it is NOT expected to match a restricted FULL179 draw.
        strata=collections.defaultdict(list)
        for gsm in sorted(use_train):
            row=roster_by_gsm[gsm]
            key=(row['batch'],phase_by_gsm[gsm]) if family=='WITHIN_BATCH_PHASE' else (row['batch'],)
            strata[key].append(gsm)
        stratum_reference={str(key):{'n':len(gs),'cases':sum(int(roster_by_gsm[g]['truth']) for g in gs)}
                           for key,gs in sorted(strata.items())}
        strata_preserved=0
        if coverage:
            for rep in range(1,200):
                if all(sum(byrep[rep][g] for g in gs)==sum(int(roster_by_gsm[g]['truth']) for g in gs)
                       for gs in strata.values()):
                    strata_preserved+=1
                elif len(failures)<10:
                    failures.append({'failure':'stratum_case_counts_changed','replicate':rep})
        # FULL179 has an independently pinned source roster v3.66 and must
        # match it exactly, element-by-element. DEDUP153 has NO such source
        # in v3.66, so an equality demand would be a false provenance test.
        reference_matches=None
        n_mismatch_labels=None
        projected_matches=None
        if coverage:
            projected_matches=0
            for rep in range(1,200):
                projected_bad=sum(legacy[(family,rep)][g]!=byrep[rep][g] for g in use_train)
                if projected_bad==0:projected_matches+=1
                if layout=='FULL179' and projected_bad and len(failures)<10:
                    for gsm in sorted(use_train):
                        if legacy[(family,rep)][gsm]!=byrep[rep][gsm]:
                            failures.append({'failure':'frozen_v366_mismatch','replicate':rep,'gsm':gsm,
                                             'reference':legacy[(family,rep)][gsm],'current':byrep[rep][gsm]})
                            if len(failures)>=10:break
            if layout=='FULL179':
                reference_matches=projected_matches
                n_mismatch_labels=sum(legacy[(family,rep)][g]!=byrep[rep][g]
                                      for rep in range(1,200) for g in use_train)
        complete=(coverage and n_distinct==199 and strata_preserved==199)
        if layout=='FULL179':
            if complete and reference_matches==199:
                status=('PASS_ORIGINAL_V371_FILE_HASH_AND_V366_LABELS' if actual==expected
                        else 'PASS_V366_SEMANTICS_V371_FILE_HASH_UNRESOLVED')
            else:
                status='HOLD_FULL179_V366_SEMANTIC_OR_STRATUM_MISMATCH'
        elif complete:
            status=('PASS_DEDUP153_PINNED_FILE_HASH_AND_STRATA' if actual==expected
                    else 'PROTOCOL_CONSISTENT_DEDUP153_HISTORICAL_SEED_UNVERIFIED')
        else:
            status='HOLD_DEDUP153_PERMUTATION_ROSTER_OR_STRATUM_MISMATCH'
        rec={'family':layout+'/'+family,'status':status,'file':str(path),'size_bytes':path.stat().st_size,
             'expected_protocol_sha256':expected,'actual_file_sha256':actual,'byte_hash_match':actual==expected,
             'observed_rows':len(rows),'expected_rows':199*len(use_train),'n_replicates':len(byrep),
             'n_distinct_assignments':n_distinct,'n_replicates_preserving_stratum_counts':strata_preserved,
             'stratum_reference':stratum_reference,
             'n_replicates_matching_v366':reference_matches,
             'n_label_mismatches_vs_v366':n_mismatch_labels,
             'n_replicates_equal_to_full179_restriction':projected_matches,
             'coverage_ok':coverage,'first_issues':failures[:10],
             'DEDUP153_is_independently_permuted_not_restricted_FULL179':layout=='DEDUP153',
             'v371_assignment_seed_exactness_proved':actual==expected}
        evidence.append(rec)
        report['hashes']['v371_assign_'+layout+'__'+family]={
            'path':str(path),'sha256':actual,'expected_sha256':expected,'match':actual==expected,
            'protocol_status':status,'n_matching_v366':reference_matches,
            'strata_preserved':strata_preserved}
        print('[ASSIGNMENT]',layout,family,
              'byte_sha=', 'MATCH' if actual==expected else 'MISMATCH',
              'same_as_v366=',(str(reference_matches)+'/199') if reference_matches is not None else 'NOT_APPLICABLE',
              'strata_preserved=',str(strata_preserved)+'/199',
              'distinct=',n_distinct,'status=',status,flush=True)
        if status.startswith('HOLD_'):blockers.append(layout+'/'+family)
        if actual!=expected:
            report.setdefault('warnings',[]).append('UNRESOLVED_V371_FILE_SHA_'+layout+'_'+family)
        if layout=='DEDUP153' and actual!=expected:
            report.setdefault('warnings',[]).append(
                'DEDUP153 historical assignment identity is UNVERIFIED; pilot matches, if any, are conditional.')
        ass_data[(layout,family)]=byrep
    report['assignment_semantic_gate']={
        'source_sha256':digest(frozen_source),'source_rows':reg_count,
        'n_families':3,'entries':evidence,'blockers':blockers,
        'status':'HOLD' if blockers else 'PROTOCOL_COMPATIBLE_WITH_UNRESOLVED_V371_SOURCE_HASHES',
        'interpretation':(
            'FULL179 assignments match the exact v366 register for 199/199 per family. '
            'DEDUP153 is a separately permuted within-batch 153-sample design and must '
            'preserve batch-specific case counts, rather than match FULL179 projected '
            'onto the deduplicated sample set. Three v371 byte-hash pins remain unresolved; '
            'DEDUP153 original seed/replicate-index identity is not established.'),
        'not_proof_of_original_v371_driver':True}
    save_json(run_dir/'assignment_semantic_gate.json',report['assignment_semantic_gate'])
    if blockers:
        diagnostic_zip=run_dir/'Paper8_stage8D_ASSIGNMENT_DIAGNOSTICS_PRIVATE.zip'
        with zipfile.ZipFile(diagnostic_zip,'w',compression=zipfile.ZIP_DEFLATED) as z:
            z.write(run_dir/'assignment_semantic_gate.json',arcname='assignment_semantic_gate.json')
            for layout,family in FAMILIES:
                file=V371/EXPECTED_ASSIGNMENT_FILES[(layout,family)]
                z.write(file,arcname='assignments/'+file.name)
        report['private_assignment_diagnostics_zip']=str(diagnostic_zip)
        raise RuntimeError('HOLD_LAYOUT_SPECIFIC_ASSIGNMENT_CHECK: '+','.join(blockers))
    if report.get('warnings'):
        print('[WARNING] Historical v371 file SHA pins unresolved. DEDUP153 tested against '
              'frozen 153-sample batch-stratum constraints, NOT v366 projected labels.',flush=True)
    report['assignment_equivalence_v366']={a+'/'+b:{'matching_replicates_vs_v366':199,'of':199}
                                             for (a,b) in FAMILIES if a=='FULL179'}
    report['assignment_protocol_consistency']={
        a+'/'+b:{'n_replicates':len(reps),'n_training':len(reps[1]),
                 'source_verified':bool(next(x for x in evidence if x['family']==a+'/'+b)['byte_hash_match']),
                 'status':next(x for x in evidence if x['family']==a+'/'+b)['status']}
        for (a,b),reps in ass_data.items()}
    report['n_test']=len(test_ids)
    report['n_removed_training_exact_cel']=len(removed)
    # Mandatory presence only. Do not hash source files without historic known hashes.
    for p in (V371/'verified_216_metadata.csv', V371/'removed_train26.txt',V371/'core270_genes.txt',
              V371/'baseline_FULL179_NO_OUTCOME.tsv.gz',V371/'baseline_DEDUP153_NO_OUTCOME.tsv.gz',
              V371/'baseline_FULL179_OUTCOME.tsv.gz',V371/'baseline_DEDUP153_OUTCOME.tsv.gz',
              V369/'08_retrained_RMA_null_metrics.csv',V371/'null_summaries.csv'):
        if not p.is_file() or p.is_symlink():raise RuntimeError('HOLD_REQUIRED_HISTORIC_ARTIFACT_MISSING_'+p.name)
    meta=read_dicts(V371/'verified_216_metadata.csv')
    if len(meta)!=216 or {r['gsm'] for r in meta}!={r['gsm'] for r in rost}:
        raise RuntimeError('HOLD_V371_META_ROSTER_MISMATCH')
    orig_by={r['gsm']:r for r in rost}
    if any(str(r['truth'])!=str(orig_by[r['gsm']]['truth']) or str(r['batch'])!=str(orig_by[r['gsm']]['batch']) or str(r['split'])!=str(orig_by[r['gsm']]['split']) for r in meta):
        raise RuntimeError('HOLD_V371_META_TRUTH_BATCH_OR_SPLIT_MISMATCH')
    removed_lines={x.strip() for x in (V371/'removed_train26.txt').read_text().splitlines() if x.strip()}
    if removed_lines!=removed:raise RuntimeError('HOLD_V371_REMOVAL_LIST_MISMATCH')
    allgenes=[x.strip() for x in (V371/'core270_genes.txt').read_text().splitlines() if x.strip()]
    if len(allgenes)!=270 or len(set(allgenes))!=270:raise RuntimeError('HOLD_CORE270_PANEL_LOCK')
    scan_downloads_for_pinned_sources(report)
    report['full_genome_R_combat_helper_verified']=True
    report['pilot_spec']={'RMA_replicates_per_family':list(RMA_REPLICATES),'ComBat_replicate_per_family':COMBAT_REPLICATE,
        'family_set':[a+'/'+b for a,b in FAMILIES],'configurations':CONFIGS,'standardization':'StandardScaler train-only',
        'classifier':'scikit-learn SAGA elastic-net','max_iter_reconstructed_fixed_not_source_verified':RECONSTRUCTED_ITERATIONS,
        'random_state':123,'tol':1e-4,'n_jobs':1,'no_test_tuning':True}
    return proto


def load_rows_map(path,kind):
    rows=read_dicts(path)
    if not rows:raise RuntimeError('HOLD_ARCHIVED_METRICS_EMPTY_'+kind)
    return rows


def prepare_frozen_refs(out_dir):
    paths={
        'v369_RMA_null_metrics':V369/'08_retrained_RMA_null_metrics.csv',
        'v369_RMA_baselines':V369/'08a_RMA_verified_baselines.json',
        'v371_baseline_models':V371/'baseline_combat_models.csv',
    }
    for key,path in paths.items():
        if not path.is_file():raise RuntimeError('HOLD_ARCHIVED_SOURCE_MISSING_'+key)
    # This small set provides an immutable audit trail for the downstream worker.
    saved={}
    for key,path in paths.items():
        dest=out_dir/path.name
        shutil.copyfile(path,dest)
        saved[key]={'path':str(path),'sha256':digest(path),'copied_to':str(dest)}
    return saved


def R_replay_one(family,layout,run_dir,report):
    slug=layout+'__'+family
    out=run_dir/'R_combat_null'/slug
    out.mkdir(parents=True,exist_ok=False)
    expected=out/'rep_001.tsv.gz'
    cmd=[str(R_BIN),str(run_dir/'fullgene_combat_helper.R'),
        'null',layout,'OUTCOME',str(V371/'verified_216_metadata.csv'),
        str(V371/'core270_genes.txt'),str(ORIG/'all_phases_all_genes_precombat.tsv'),
        str(ORIG/'all_phases_all_genes_counts.tsv'),str(V352/'FULL216__OUTCOME_INFORMED_COMBAT.tsv'),
        str(V371/EXPECTED_ASSIGNMENT_FILES[(layout,family)]),str(out),str(COMBAT_REPLICATE),str(COMBAT_REPLICATE)]
    print('[R_COMBAT_REPLAY]',slug,'rep=',COMBAT_REPLICATE,flush=True)
    proc=subprocess.run(cmd,cwd=str(run_dir),capture_output=True,text=True,timeout=300)
    (run_dir/('R_'+slug+'.log')).write_text(proc.stdout+'\n--- STDERR ---\n'+proc.stderr,encoding='utf-8')
    row={'family':slug,'returncode':proc.returncode,'output':str(expected),'log':'R_'+slug+'.log',
         'fresh_r_combat_executed':True,'true_test_outcome_used_by_combat':True}
    report['R_replays'].append(row)
    if proc.returncode!=0 or not expected.is_file() or expected.stat().st_size==0:
        raise RuntimeError('HOLD_R_COMBAT_REPLAY_FAILED_'+slug+' '+proc.stderr[-400:])
    row['sha256']=digest(expected)
    return expected


def numeric_worker(run_dir):
    # Worker is intentionally invoked in the PREVIOUSLY LOCKED v3.71 numpy/sklearn environment.
    import numpy as np
    import sklearn
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    worker_report={'status':'HOLD_NOT_STARTED','python':sys.version.split()[0],
                   'numpy':np.__version__,'sklearn':sklearn.__version__,
                   'model_source':'independent reconstructed fixed specification; original v371 driver unavailable',
                   'expected_environment':'Python 3.12.15, NumPy 2.5.4, scikit-learn 1.9.1',
                   'max_iter':RECONSTRUCTED_ITERATIONS,'baseline_rows':[],'pilot_rows':[],'warnings':[]}
    save_json(run_dir/'worker_report.json',worker_report)
    if (sys.version.split()[0],np.__version__,sklearn.__version__)!=('3.12.15','2.5.4','1.9.1'):
        raise RuntimeError('HOLD_V371_REPLAY_ENVIRONMENT_NOT_IDENTICAL')
    def matrix_from_tsv(path,kind):
        opener=gzip.open if str(path).endswith('.gz') else open
        with opener(path,'rt',encoding='utf-8-sig',newline='') as handle:
            reader=csv.reader(handle,delimiter='\t')
            header=next(reader)
            samples=header[1:]
            if len(samples) not in (216,190) or len(set(samples))!=len(samples):
                raise RuntimeError('HOLD_'+kind+'_SAMPLE_COLUMNS')
            names=[]; vectors=[]
            for r in reader:
                if len(r)!=len(header):raise RuntimeError('HOLD_'+kind+'_TSV_ROW_WIDTH')
                names.append(r[0]);vectors.append([float(x) for x in r[1:]])
        if len(vectors)!=270 or len(set(names))!=270:
            raise RuntimeError('HOLD_'+kind+'_CORE270_ROWS')
        array=np.asarray(vectors,dtype=np.float64).T
        if not np.isfinite(array).all():raise RuntimeError('HOLD_'+kind+'_NONFINITE')
        return samples,names,array
    meta=read_dicts(V371/'verified_216_metadata.csv')
    by_gsm={x['gsm']:x for x in meta}
    actual_train=sorted(x['gsm'] for x in meta if x['split']=='TRAIN179')
    removed={x.strip() for x in (V371/'removed_train26.txt').read_text().splitlines() if x.strip()}
    actual_test=sorted(x['gsm'] for x in meta if x['split']=='TEST37')
    if len(actual_train)!=179 or len(removed)!=26 or len(actual_test)!=37:
        raise RuntimeError('HOLD_META_TRAIN_TEST_COUNTS')
    truth={g:int(r['truth']) for g,r in by_gsm.items()}
    raw_assign={}
    for layout,family in FAMILIES:
        rows=read_dicts(V371/EXPECTED_ASSIGNMENT_FILES[(layout,family)])
        assigns=collections.defaultdict(dict)
        for row in rows:
            rep=int(row['replicate'])
            if rep in RMA_REPLICATES or rep==COMBAT_REPLICATE:
                assigns[rep][row['gsm']]=int(row['permuted_label'])
        raw_assign[(layout,family)]=assigns
    rma_rows=load_rows_map(V369/'08_retrained_RMA_null_metrics.csv','RMA_NULL')
    rma_ref={(x['layout'],x['family'],int(x['replicate']),x['configuration']):x for x in rma_rows}
    if len(rma_ref)!=1194:raise RuntimeError('HOLD_RMA_REF_1194_UNIQUE')
    prior_baseline=json.loads((V369/'08a_RMA_verified_baselines.json').read_text())
    rma_baseline={(x['layout'],x['configuration']):x for x in prior_baseline}
    if len(rma_baseline)!=4:raise RuntimeError('HOLD_RMA_BASELINE_REF_FOUR')
    v371_ref={}
    for layout,family in FAMILIES:
        for prep in ('NO_OUTCOME','OUTCOME'):
            path=V371/('null_'+layout+'_'+family+'_'+prep+'.csv')
            rows=read_dicts(path)
            if len(rows)!=398:
                raise RuntimeError('HOLD_V371_REF_NOT_398_ROWS_'+path.name)
            for x in rows:
                key=(layout,family,prep,int(x['replicate']),x['config'])
                if key in v371_ref:raise RuntimeError('HOLD_V371_REF_DUPLICATE_KEY')
                v371_ref[key]=x
    if len(v371_ref)!=2388:raise RuntimeError('HOLD_V371_REF_NOT_2388_ROWS')
    baseline_combat=read_dicts(V371/'baseline_combat_models.csv')
    baseline_lookup={}
    if baseline_combat:
        for r in baseline_combat:
            config=r.get('configuration') or r.get('config')
            if config and {'layout','preprocessing','TP','FN','TN','FP','BA'}<=set(r):
                key=(r['layout'],r['preprocessing'],config)
                baseline_lookup[key]=r
            else:
                worker_report['warnings'].append('BASELINE_COMBAT_SCHEMA_UNRECOGNIZED_'+','.join(baseline_combat[0]))
                break
    expected_genes=[x.strip() for x in (V371/'core270_genes.txt').read_text().splitlines() if x.strip()]
    all_matrix={'RMA_ONLY': V352/'FULL216__NO_COMBAT.tsv',
               'NO_OUTCOME_FULL179': V371/'baseline_FULL179_NO_OUTCOME.tsv.gz',
               'NO_OUTCOME_DEDUP153': V371/'baseline_DEDUP153_NO_OUTCOME.tsv.gz',
               'OUTCOME_FULL179': V371/'baseline_FULL179_OUTCOME.tsv.gz',
               'OUTCOME_DEDUP153': V371/'baseline_DEDUP153_OUTCOME.tsv.gz'}
    pool={}
    for name,path in all_matrix.items():
        cols,genes,X=matrix_from_tsv(path,name)
        if genes!=expected_genes:
            raise RuntimeError('HOLD_'+name+'_FROZEN_PANEL_GENE_ORDER')
        pool[name]=(cols,X)
        print('[MATRIX]',name,X.shape,flush=True)
    def fit(X,y,Z,C,l1):
        est=make_pipeline(StandardScaler(),LogisticRegression(penalty='elasticnet',solver='saga',
                         C=C,l1_ratio=l1,random_state=123,max_iter=RECONSTRUCTED_ITERATIONS,
                         tol=1e-4,n_jobs=1))
        with warnings.catch_warnings(record=True) as obs:
            warnings.simplefilter('always',ConvergenceWarning)
            est.fit(X,y)
        return est.decision_function(Z),int(est.named_steps['logisticregression'].n_iter_[0]),[str(x.message) for x in obs]
    def compute_matrix_case(Xcols,X,layout,assign):
        train=set(actual_train) if layout=='FULL179' else set(actual_train)-removed
        sample_map={g:j for j,g in enumerate(Xcols)}
        if not train<=sample_map.keys() or not set(actual_test)<=sample_map.keys():
            raise RuntimeError('HOLD_MATRIX_MISSING_REQUIRED_GSM')
        train_order=sorted(train)
        tr=X[[sample_map[g] for g in train_order],:]
        te=X[[sample_map[g] for g in actual_test],:]
        y=np.asarray([assign.get(g,truth[g]) if assign is not None else truth[g] for g in train_order],dtype=int)
        yt=np.asarray([truth[g] for g in actual_test],dtype=int)
        if not set(y.tolist())<={0,1} or not set(yt.tolist())<={0,1}:
            raise RuntimeError('HOLD_TRAIN_TEST_LABELS')
        return tr,y,te,yt
    def metric(score,yt):
        yp=(score>=0).astype(int)
        tp=int(np.sum((yt==1)&(yp==1)));fn=int(np.sum((yt==1)&(yp==0)))
        tn=int(np.sum((yt==0)&(yp==0)));fp=int(np.sum((yt==0)&(yp==1)))
        return {'TP':tp,'FN':fn,'TN':tn,'FP':fp,'BA':float(0.5*(tp/(tp+fn)+tn/(tn+fp)))}
    def comparison(tag,ident,layout,family,prep,config,rep,reference,score,n_iter,model_warnings,yt):
        x=metric(score,yt)
        v={name:int(reference[name]) for name in ('TP','FN','TN','FP')}
        expect_ba=float(reference['BA'])
        valid=all(x[k]==v[k] for k in v) and abs(x['BA']-expect_ba)<1e-10
        score_err=None
        if 'score_by_sorted_GSM' in reference:
            ref=np.asarray(json.loads(reference['score_by_sorted_GSM']),dtype=np.float64)
            if ref.size!=37 or not np.isfinite(ref).all():raise RuntimeError('HOLD_SCORE_REFERENCE_NOT_37')
            score_err=float(np.max(np.abs(ref-score)))
            valid=valid and score_err<1e-6
        row={'type':tag,'identifier':ident,'layout':layout,'family':family,'preprocessing':prep,
             'configuration':config,'replicate':rep,'TP':x['TP'],'FN':x['FN'],'TN':x['TN'],'FP':x['FP'],
             'BA':x['BA'],'archived_TP':v['TP'],'archived_FN':v['FN'],'archived_TN':v['TN'],'archived_FP':v['FP'],
             'archived_BA':expect_ba,'n_iter':n_iter,'convergence_warnings':len(model_warnings),
             'model_warnings':model_warnings[:3],'max_score_abs_error':score_err,'agrees':bool(valid)}
        worker_report['baseline_rows' if tag=='BASELINE' else 'pilot_rows'].append(row)
        print('[MODEL]',tag,ident,config,'BA',round(x['BA'],6),'archived',round(expect_ba,6),
              'SCORE_MAX_DIFF',score_err,'MATCH' if valid else 'HOLD',flush=True)
    # Baselines: detect any mismatch before running fresh full-gene ComBat.
    for layout in ('FULL179','DEDUP153'):
        for config,(C,l1) in CONFIGS.items():
            cols,X=pool['RMA_ONLY']; tr,y,te,yt=compute_matrix_case(cols,X,layout,None)
            score,n_iter,ww=fit(tr,y,te,C,l1)
            comparison('BASELINE','RMA_'+layout,layout,'-', 'RMA_ONLY',config,0,
                       rma_baseline[(layout,config)],score,n_iter,ww,yt)
            for prep in ('NO_OUTCOME','OUTCOME'):
                cols,X=pool[prep+'_'+layout];tr,y,te,yt=compute_matrix_case(cols,X,layout,None)
                score,n_iter,ww=fit(tr,y,te,C,l1)
                if (layout,prep,config) in baseline_lookup:
                    comparison('BASELINE','COMBAT_'+prep+'_'+layout,layout,'-',prep,config,0,
                               baseline_lookup[(layout,prep,config)],score,n_iter,ww,yt)
                else:
                    worker_report['warnings'].append('COMBAT_BASELINE_REF_SCHEMA_UNAVAILABLE_'+layout+'_'+prep+'_'+config)
    for layout,family in FAMILIES:
        for rep in RMA_REPLICATES:
            cols,X=pool['RMA_ONLY'];tr,y,te,yt=compute_matrix_case(cols,X,layout,raw_assign[(layout,family)][rep])
            for config,(C,l1) in CONFIGS.items():
                score,n_iter,ww=fit(tr,y,te,C,l1)
                comparison('RMA_NULL','RMA_'+layout+'_'+family+'_REP'+str(rep),layout,family,'RMA_ONLY',config,
                          rep,rma_ref[(layout,family,rep,config)],score,n_iter,ww,yt)
    # One no-outcome label permutation per family using frozen label-blind full-gene ComBat outputs.
    for layout,family in FAMILIES:
        rep=COMBAT_REPLICATE
        cols,X=pool['NO_OUTCOME_'+layout]
        tr,y,te,yt=compute_matrix_case(cols,X,layout,raw_assign[(layout,family)][rep])
        for config,(C,l1) in CONFIGS.items():
            score,n_iter,ww=fit(tr,y,te,C,l1)
            comparison('COMBAT_NO_OUTCOME_NULL','NO_OUTCOME_'+layout+'_'+family,layout,family,'NO_OUTCOME',config,
                       rep,v371_ref[(layout,family,'NO_OUTCOME',rep,config)],score,n_iter,ww,yt)
    # Freshly run ComBat with permuted TRAIN labels and true TEST labels, on 44,341 rows.
    for layout,family in FAMILIES:
        rep=COMBAT_REPLICATE
        p=run_dir/'R_combat_null'/(layout+'__'+family)/'rep_001.tsv.gz'
        if not p.is_file():
            raise RuntimeError('HOLD_FRESH_COMBAT_OUTPUT_MISSING_'+layout+'_'+family)
        cols,genes,X=matrix_from_tsv(p,'fresh_'+layout+'_'+family)
        if genes!=expected_genes:raise RuntimeError('HOLD_FRESH_COMBAT_PANEL_ORDER')
        tr,y,te,yt=compute_matrix_case(cols,X,layout,raw_assign[(layout,family)][rep])
        for config,(C,l1) in CONFIGS.items():
            score,n_iter,ww=fit(tr,y,te,C,l1)
            comparison('COMBAT_OUTCOME_NULL','OUTCOME_'+layout+'_'+family,layout,family,'OUTCOME',config,
                       rep,v371_ref[(layout,family,'OUTCOME',rep,config)],score,n_iter,ww,yt)
    with (run_dir/'pilot_comparisons.csv').open('w',newline='',encoding='utf-8') as f:
        rows=worker_report['baseline_rows']+worker_report['pilot_rows']
        keys=['type','identifier','layout','family','preprocessing','configuration','replicate',
              'TP','FN','TN','FP','BA','archived_TP','archived_FN','archived_TN','archived_FP','archived_BA',
              'n_iter','convergence_warnings','max_score_abs_error','agrees']
        w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(rows)
    n_fitted=len(worker_report['baseline_rows'])+len(worker_report['pilot_rows'])
    n_mismatches=sum(not x['agrees'] for x in worker_report['baseline_rows']+worker_report['pilot_rows'])
    worker_report['n_fits']=n_fitted
    worker_report['n_compared']=n_fitted
    worker_report['n_mismatches']=n_mismatches
    worker_report['n_RMA_null_fits']=sum(x['type']=='RMA_NULL' for x in worker_report['pilot_rows'])
    worker_report['n_ComBat_null_fits']=sum(x['type'].startswith('COMBAT_') for x in worker_report['pilot_rows'])
    worker_report['status']='PASS_RECONSTRUCTED_NUMERIC_PILOT_DRIVER_UNVERIFIED' if n_mismatches==0 else 'HOLD_RECONSTRUCTED_NUMERIC_PILOT_DISAGREEMENT'
    save_json(run_dir/'worker_report.json',worker_report)
    print('[WORKER_STATUS]',worker_report['status'],'fits',n_fitted,'disagreements',n_mismatches,flush=True)
    return 0 if n_mismatches==0 else 3


def run():
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run_dir=OUTROOT/('Run_'+stamp)
    run_dir.mkdir(parents=True,exist_ok=False)
    report={'stage':STAGE,'status':'HOLD_NOT_STARTED','created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
            'policy':'PROJECT_READ_ONLY; OUTPUTS_ONLY_DOWNLOADS; NO_DOWNLOADS_INSTALLS_TUNING; NO_MODEL_INPUT_MUTATIONS',
            'not_original_v371_python_driver_execution':True,'not_3582_full_refit_replay':True,
            'DEDUP153_historical_assignment_identity_unverified':True,
            'not_external_validation':True,'source_status':'original_v371_driver_missing','hashes':{},
            'R_replays':[],'messages':[],'errors':[],'warnings':[],'output_dir':str(run_dir)}
    try:
        proto=validate_sources(run_dir,report)
        prepare_frozen_refs(run_dir)
        if not R_BIN.is_file() or not PY_BIN.is_file():
            raise RuntimeError('HOLD_R_OR_PYTHON_ENVIRONMENT_MISSING')
        shutil.copyfile(V371/'fullgene_combat_helper.R',run_dir/'fullgene_combat_helper.R')
        (run_dir/'REPRODUCTION_SCOPE.txt').write_text(
            'Independent fixed-method pilot, not original v3.71 Python driver.\n'
            'Three RMA-only replicates per family; one full-gene outcome-informed ComBat\n'
            'per family; uses true TEST diagnoses in outcome-informed ComBat, so\n'
            'these results are NOT independent diagnostic validation.\n'
            'Model max_iter=5000 is an explicit reconstruction assumption until\n'
            'the original Python driver is recovered; the run never tunes it.\n'
            'DEDUP153 uses a separate within-batch permutation after removal of 26 arrays;\n'
            'it must not be compared with a restricted FULL179 permutation.\n'
            'v3.71 assignment file hashes and DEDUP153 original permutation seed\n'
            'remain unverified; pilot outcomes are conditional on observed assignments.\n'
            'Original project is read-only. Do not publish this private run as-is.\n',encoding='utf-8')
        # Replay only three ComBat permutations, not 199*3.
        for layout,family in FAMILIES:
            R_replay_one(family,layout,run_dir,report)
        env=dict(os.environ)
        env.update(PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        cmd=[str(PY_BIN),str(Path(__file__).resolve()),'--numeric-worker',str(run_dir)]
        print('[MODEL_WORKER] pinned Python3.12.15 NumPy2.5.4 sklearn1.9.1',flush=True)
        # Stream worker logs live; retain full log under Downloads.
        with (run_dir/'python_worker.log').open('w',encoding='utf-8') as handle:
            proc=subprocess.run(cmd,cwd=str(run_dir),env=env,stdout=handle,stderr=subprocess.STDOUT,
                                timeout=600)
        report['python_worker_exit']=proc.returncode
        txt=(run_dir/'python_worker.log').read_text(encoding='utf-8',errors='replace')
        print(txt[-9000:],flush=True)
        report['python_worker_log']='python_worker.log'
        wp=run_dir/'worker_report.json'
        if wp.is_file():report['worker']=json.loads(wp.read_text())
        report['status']=report.get('worker',{}).get('status','HOLD_PYTHON_WORKER_FAILED') if proc.returncode in (0,3) else 'HOLD_PYTHON_WORKER_FAILED'
        if report['status']=='PASS_RECONSTRUCTED_NUMERIC_PILOT_DRIVER_UNVERIFIED' and report.get('warnings'):
            report['status']='PASS_CONDITIONAL_NUMERIC_PILOT_V371_ASSIGNMENT_PROVENANCE_UNVERIFIED'
            report['warnings'].append('Numerical pilot matched on supplied v371 assignment files. Exact original driver and v371 file hashes are NOT verified; DEDUP153 seed/replicate-index identity remains unverified.')
    except Exception as exc:
        report['errors'].append(type(exc).__name__+': '+str(exc))
        report['status']='HOLD_STAGE8D_PREFLIGHT_OR_COMPUTATION'
        traceback.print_exc()
    finally:
        save_json(run_dir/'stage8d_pilot.json',report)
        print('[REPORT]',run_dir/'stage8d_pilot.json',flush=True)
        print('[STATUS]',report['status'],flush=True)
        print('[SCOPE] 18 RMA-only + 12 ComBat pilot fits planned; NOT full 3582.',flush=True)
    return 0 if report['status'] in ('PASS_RECONSTRUCTED_NUMERIC_PILOT_DRIVER_UNVERIFIED','PASS_CONDITIONAL_NUMERIC_PILOT_V371_ASSIGNMENT_PROVENANCE_UNVERIFIED') else 2


def selftest():
    from zipfile import ZipFile
    p=Path('/mnt/data/Paper8_stage8C_NULL_EVIDENCE_KIT_PRIVATE.zip')
    if p.is_file():
        with ZipFile(p) as z:
            protocol=json.loads(z.read('known_protocol/protocol_lock.json'))
            assert len(protocol['permutation_sha256'])==3
            s=z.read('selected_provenance/results/p8_v366_batch_split_recovery_20261109T002942Z/12_locked_training_null_assignments_NOT_FITTED.csv')
            r=list(csv.DictReader(io.StringIO(s.decode('utf-8-sig'))))
            assert len(r)==71242
            assert collections.Counter((x['training_layout'],x['family']) for x in r)=={
                ('FULL179_ONLY','WITHIN_BATCH'):35621,('FULL179_ONLY','WITHIN_BATCH_PHASE'):35621}
            assert len([x for x in protocol['assignments'] if x['n']==199])==3
    assert RECONSTRUCTED_ITERATIONS==5000
    assert len(RMA_REPLICATES)==3 and COMBAT_REPLICATE==1 and len(FAMILIES)==3
    print('[SELFTEST] PASS_STAGE8D_PINNED_PROTOCOL_NULL_FAMILY_SCHEMA_AND_BOUNDED_COUNTS')


if __name__=='__main__':
    if len(sys.argv)==2 and sys.argv[1]=='--selftest':
        selftest()
    elif len(sys.argv)==3 and sys.argv[1]=='--numeric-worker':
        try:
            sys.exit(numeric_worker(Path(sys.argv[2])))
        except Exception:
            traceback.print_exc()
            sys.exit(2)
    elif len(sys.argv)==1:
        sys.exit(run())
    else:
        raise SystemExit('USAGE: python3 "$HOME/Downloads/paper8_repro_stage8d_numeric_pilot_v3.py"')
