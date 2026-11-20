#!/usr/bin/env python3
"""Paper 8 reproducibility stage 2: read-only, bounded local inventory.

Does not modify research files, download datasets, launch old analysis, or hash large files.
Writes only an inventory report under ~/Downloads/Paper8_repro_stage2_inventory/.
"""
from __future__ import annotations
import csv
import datetime as dt
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

HOME = Path.home()
PROJECTS = [HOME / 'Desktop/paper8/endo-research', HOME / 'Desktop/paper8']
DOWNLOADS = HOME / 'Downloads'
OUTPUT = DOWNLOADS / 'Paper8_repro_stage2_inventory'
MAX_DEPTH = 7
MAX_FILES = 120000
MAX_EXAMPLES = 18

PATTERNS = {
    'raw_cel': re.compile(r'\.cel(?:\.gz)?$', re.I),
    'raw_archives': re.compile(r'(?:\.tar(?:\.gz)?|\.tgz|\.zip)$', re.I),
    'scripts': re.compile(r'\.(?:py|r|rmd|sh)$', re.I),
    'frozen_configs': re.compile(r'(?:prespec|freeze|frozen|manifest|seed|config|protocol|lock|params|sessioninfo|requirements|environment|renv\.lock)', re.I),
    'matrix_data': re.compile(r'\.(?:rds|rdata|h5|h5ad|hdf5|parquet|feather|npy|npz)$', re.I),
    'tables': re.compile(r'\.(?:csv|tsv)(?:\.gz)?$', re.I),
}
STAGES = {
    'RMA_and_raw': re.compile(r'(?i)\b(?:rma|affy|mas5|cel|gse51981|gse6364|fogg)\b|(?:_rma|rma_|\.cel)', re.I),
    'ComBat': re.compile(r'(?i)combat|sva|batch|no_outcome|outcome', re.I),
    'null_refits': re.compile(r'(?i)null|refit|permut|label|seed|3582|1194|v3[.]7', re.I),
    'Huang': re.compile(r'(?i)huang|random.?forest|rf_|dedup|leave.one', re.I),
    'benchmark_324': re.compile(r'(?i)324|nested|date|phase|gse51981|v4[.]30|v430', re.I),
}
SKIP = {'.git', 'node_modules', '__pycache__', '.cache', '.pytest_cache', '.mypy_cache', '.next', 'build', 'dist', '.tox', 'site-packages', 'lib', 'include', 'bin', 'share', '.venv', 'venv', '.envs'}


def safe_stat(p):
    try:
        st = p.stat()
        return {'bytes': st.st_size, 'modified_utc': dt.datetime.fromtimestamp(st.st_mtime, dt.timezone.utc).isoformat()}
    except OSError:
        return {'bytes': None, 'modified_utc': None}


def scan(root, max_depth=MAX_DEPTH):
    counts = {x: 0 for x in PATTERNS}
    samples = {x: [] for x in PATTERNS}
    stage_samples = {x: [] for x in STAGES}
    stage_counts = {x: 0 for x in STAGES}
    total = 0
    skipped = 0
    errors = 0
    if not root.is_dir():
        return {'root': str(root), 'exists': False}
    for base, dirs, files in os.walk(root, followlinks=False):
        basep = Path(base)
        depth = len(basep.relative_to(root).parts)
        dirs[:] = sorted(d for d in dirs if d not in SKIP and depth < max_depth and not (basep / d).is_symlink())
        for name in sorted(files):
            f = basep / name
            if f.is_symlink():
                continue
            total += 1
            if total > MAX_FILES:
                skipped += 1
                continue
            rel = str(f.relative_to(root))
            low = rel.lower()
            try:
                st = f.stat()
            except OSError:
                errors += 1
                continue
            details = {'path': str(f), 'size_bytes': st.st_size}
            for kind, regex in PATTERNS.items():
                if regex.search(name):
                    counts[kind] += 1
                    if len(samples[kind]) < MAX_EXAMPLES:
                        samples[kind].append(details)
            for stage, regex in STAGES.items():
                if regex.search(low):
                    stage_counts[stage] += 1
                    if len(stage_samples[stage]) < MAX_EXAMPLES:
                        stage_samples[stage].append(details)
    return {'root': str(root), 'exists': True, 'files_seen': total, 'max_files': MAX_FILES,
            'depth_limit': max_depth, 'skipped_by_limit': skipped, 'stat_errors': errors,
            'counts': counts, 'examples': samples, 'stage_counts': stage_counts, 'stage_examples': stage_samples}


def command_version(args):
    try:
        result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, timeout=12, check=False)
        return {'command': args, 'returncode': result.returncode, 'output': result.stdout.strip()[:650]}
    except (OSError, subprocess.TimeoutExpired) as e:
        return {'command': args, 'error': type(e).__name__ + ': ' + str(e)}


def find_envs(root):
    found = []
    for path in [root / '.envs', root / '.venv', root / 'venv', root.parent / '.envs']:
        if path.exists():
            found.append({'path': str(path), 'Rscript': [str(p) for p in path.glob('*/bin/Rscript')][:20],
                          'python': [str(p) for p in path.glob('*/bin/python*')][:20]})
    return found


def main():
    roots = [p for p in PROJECTS if p.is_dir()]
    if roots and roots[0] == PROJECTS[0]:
        # Avoid scanning the same files twice; inspect nested project once.
        roots = [roots[0]]
    if not roots:
        print('PROJECT_NOT_FOUND: Expected ~/Desktop/paper8/endo-research or ~/Desktop/paper8')
        print('Nothing was changed. Check the project directory and rerun.')
        sys.exit(2)
    print('PAPER 8 | STAGE 2 | READ-ONLY LOCAL INVENTORY')
    print('[POLICY] No raw reads, hashing, downloads, training, or writes to research directories')
    records = [scan(p) for p in roots]
    tools = {'python': command_version([sys.executable, '--version'])}
    for name, args in [('Rscript', ['Rscript', '--version']), ('git', ['git', '--version'])]:
        tools[name] = command_version(args) if shutil.which(args[0]) else {'available_on_PATH': False}
    envs = find_envs(PROJECTS[0])
    # Probe R packages only if an existing Rscript executable is available, without installation.
    r_candidates = [Path(p) for e in envs for p in e.get('Rscript', [])]
    if shutil.which('Rscript'):
        r_candidates.insert(0, Path(shutil.which('Rscript')))
    pkg_expr = 'for(p in c("affy","hgu133plus2cdf","sva","Biobase","glmnet","randomForest","caret")){cat(p,":",if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else"MISSING","\\n",sep="")}'
    r_packages = []
    for script in dict.fromkeys(map(str, r_candidates)):
        r_packages.append({'executable': script, 'probe': command_version([script, '-e', pkg_expr])})
    report = {'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'platform': {'system': platform.system(), 'release': platform.release(), 'arch': platform.machine(), 'python': sys.version.split()[0]},
              'policy': 'READ ONLY RESEARCH DATA; limited directory traversal; no hidden dataset contents; no model execution',
              'project_roots': records, 'environments': envs, 'tools': tools, 'r_packages': r_packages,
              'validation_level': 'Inventory only, not proof of any computational result'}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    outfile = OUTPUT / 'inventory.json'
    outfile.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    lines = ['PAPER 8 REPRODUCIBILITY STAGE 2 INVENTORY', 'Status: INVENTORY_ONLY (no analysis reruns)', '']
    for rec in records:
        lines += [f"Root: {rec['root']}", f"Files inspected: {rec['files_seen']}",
                  f"Bounded traversal: depth <= {MAX_DEPTH}, cap {MAX_FILES}"]
        for typ, cnt in rec['counts'].items():
            lines.append(f'{typ}: {cnt}')
        lines += ['', 'STAGE CANDIDATES:']
        for stage, cnt in rec['stage_counts'].items():
            lines.append(f'{stage}: {cnt} candidate paths')
        lines.append('')
    lines += ['Environment roots:'] + [f"  {e['path']}" for e in envs]
    lines += ['R package probes:'] + [f"  {p['executable']}: {p['probe'].get('output', p['probe'].get('error', ''))}" for p in r_packages]
    lines += ['', 'Detailed sample paths: inventory.json',
              'WARNING: File names and full local paths may be personal. Review before sharing publicly.']
    (OUTPUT / 'summary.txt').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print('\n'.join(lines))
    print('[REPORT]', outfile)
    print('[STATUS] PASS_INVENTORY_ONLY')

if __name__ == '__main__':
    main()
