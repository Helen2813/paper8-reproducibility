#!/usr/bin/env python3
"""Paper 8 Stage 5: read-only preflight for archived full-cohort ComBat replay.

Writes reports only to ~/Downloads/Paper8_repro_stage5_combat_preflight.
Does not run ComBat, RMA, classifiers, installs or downloads.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path.home() / 'Desktop/paper8/endo-research'
OUT = Path.home() / 'Downloads/Paper8_repro_stage5_combat_preflight'
REPORT = {
    'stage': 'P8_STAGE5_COMBAT_PREFLIGHT',
    'timestamp_utc': datetime.now(timezone.utc).isoformat(),
    'policy': 'READ_ONLY_INPUTS; NO_COMPUTATION_REPLAY; REPORT_ONLY_IN_DOWNLOADS',
    'report_status': 'HOLD',
    'candidates': [],
    'r_environment': {},
    'issues': [],
}

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def file_record(path, do_hash=False, expected=None):
    r = {'relative_path': str(path.relative_to(ROOT)), 'exists': path.is_file()}
    if r['exists']:
        r['size_bytes'] = path.stat().st_size
        if do_hash:
            r['sha256'] = sha256(path)
            if expected:
                r['expected_sha256'] = expected
                r['hash_matches'] = r['sha256'] == expected
    return r

def header_info(path):
    with path.open('r', encoding='utf-8', newline='') as f:
        line = f.readline().rstrip('\r\n')
        # Strict bounded header only: never load the matrix into memory.
        cols = line.split('\t')
    return {'column_count': len(cols), 'first_column': cols[0] if cols else '',
            'first_samples': cols[1:4], 'last_samples': cols[-3:]}

def main():
    print('PAPER 8 | STAGE 5 | COMBAT ARCHIVE VS OUTPUT PREFLIGHT')
    print('[POLICY] No RMA, ComBat, fitting, installs or project writes')
    if not ROOT.is_dir():
        raise SystemExit('[HOLD] Paper8 project root not found: ' + str(ROOT))
    relbase = Path('data/p8_cross_series_independence/fogg_original_combat')
    base = ROOT / relbase
    for name in ('20261108T015339Z', '20261108T015502Z'):
        d = base / name
        entry = {'run': name, 'directory_exists': d.is_dir()}
        report_candidates = [ROOT / 'results' / ('p8_v350e_original_combat_' + name) / 'original_combat_report.json',
                             d / 'original_combat_report.json']
        report_path = next((p for p in report_candidates if p.is_file()), None)
        entry['report_path'] = str(report_path.relative_to(ROOT)) if report_path else None
        report = None
        if report_path:
            try:
                report = json.loads(report_path.read_text(encoding='utf-8'))
                entry['report_status'] = report.get('status', 'UNKNOWN')
                entry['report_keys'] = sorted(report.keys())
            except Exception as exc:
                entry['report_error'] = type(exc).__name__ + ': ' + str(exc)[:150]
        expected = {'all_phases_all_genes_precombat.tsv': 'precombat_sha256',
                    'all_phases_all_genes_counts.tsv': 'combat_sha256',
                    'all_phases_coldata.tsv': 'metadata_sha256'}
        entry['files'] = []
        for fname, key in expected.items():
            path = d / fname
            wanted = report.get(key) if isinstance(report, dict) else None
            rec = file_record(path, do_hash=path.is_file(), expected=wanted)
            if path.is_file() and fname.endswith('.tsv'):
                try:
                    rec['header'] = header_info(path)
                except Exception as exc:
                    rec['header_error'] = type(exc).__name__
            entry['files'].append(rec)
            if not path.is_file():
                REPORT['issues'].append(name + ': missing ' + fname)
            if rec.get('hash_matches') is False:
                REPORT['issues'].append(name + ': SHA256 mismatch ' + fname)
        entry['zip_candidates'] = sorted(p.name for p in (ROOT/'results').glob('*original_combat*zip'))
        REPORT['candidates'].append(entry)
        print('[RUN]', name, 'dir=', d.is_dir(), 'report=', entry['report_status'] if report else 'MISSING')
        for f in entry['files']:
            print('  ', Path(f['relative_path']).name, 'exists=', f['exists'],
                  'hash=', f.get('hash_matches', 'NO_EXPECTED_HASH'))
    rscript = ROOT / '.envs/p8_fogg_r44/bin/Rscript'
    original = ROOT / 'data/p8_cross_series_independence/fogg_code/df538583652d3359d71959437837ea3e0d1195e6/data_prep/batch_correction.r'
    REPORT['r_environment']['Rscript'] = file_record(rscript)
    REPORT['r_environment']['author_batch_script'] = file_record(original, do_hash=original.is_file())
    if not rscript.is_file(): REPORT['issues'].append('Rscript absent')
    if not original.is_file(): REPORT['issues'].append('Original author ComBat code absent')
    strong = next((x for x in REPORT['candidates'] if x['run']=='20261108T015502Z'), None)
    valid = (strong and strong.get('report_status')=='PASS_ORIGINAL_FOGG_COMBAT_216'
             and all(x['exists'] and x.get('hash_matches', False) for x in strong['files'])
             and rscript.is_file() and original.is_file())
    REPORT['report_status'] = 'PASS_ARCHIVED_COMBAT_INPUTS_HASHED_NOT_RERUN' if valid else 'HOLD_REPLAY_PREFLIGHT'
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / 'stage5_combat_preflight.json'
    dest.write_text(json.dumps(REPORT, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    (OUT / 'summary.txt').write_text('Status: ' + REPORT['report_status'] + '\n' +
                                     '\n'.join(REPORT['issues']) + '\n', encoding='utf-8')
    print('[STATUS]', REPORT['report_status'])
    print('[OUTPUT]', dest)
    print('[NOTE] Archived data/hash consistency, not a numerical ComBat replay.')

if __name__ == '__main__':
    main()
