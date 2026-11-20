#!/usr/bin/env python3
"""Paper 8 code-only release integrity gate. Not a rights/privacy clearance."""
import hashlib, json, re, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SKIP={'.git','__pycache__','.pytest_cache','.DS_Store'}
FORBIDDEN_DIRS={'.envs','.mamba','.conda','.venv','venv','data','results','sources',
                'private_inputs','private_results','checkpoints','outputs','experiments'}
FORBIDDEN_SUFFIX={'.csv','.tsv','.xlsx','.xls','.zip','.tar','.gz','.cel','.rds','.h5',
                  '.h5ad','.npy','.npz','.pkl','.pickle','.sqlite','.db','.pdf','.log','.parquet'}
SECRET_RULES=[
  ('HOME_PATH',re.compile(r'/Users/[A-Za-z0-9._-]+/')),
  ('GH_TOKEN',re.compile('ghp_'+'[A-Za-z0-9]{30,}')),
  ('OPENAI_TOKEN',re.compile('sk-'+'[A-Za-z0-9_-]{20,}')),
  ('PRIVATE_KEY',re.compile('-----BEGIN '+'PRIVATE KEY-----')),
]

def all_files():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and
                  not any(x in SKIP for x in p.relative_to(ROOT).parts))

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def audit():
    errors=[]
    manifest=json.loads((ROOT/'SHA256SUMS.json').read_text(encoding='utf-8'))
    expected=set(manifest)
    files=all_files()
    actual={p.relative_to(ROOT).as_posix() for p in files}-{'SHA256SUMS.json'}
    for key in sorted(expected-actual):errors.append('MISSING '+key)
    for key in sorted(actual-expected):errors.append('EXTRA '+key)
    for p in files:
        rel=p.relative_to(ROOT)
        key=rel.as_posix()
        if p.is_symlink():errors.append('SYMLINK '+key)
        if any(x in FORBIDDEN_DIRS for x in rel.parts):errors.append('PRIVATE_DIR '+key)
        if p.suffix.lower() in FORBIDDEN_SUFFIX:errors.append('DATA_OR_ARCHIVE_FILE '+key)
        if p.stat().st_size > 1000000:errors.append('TOO_LARGE '+key)
        if key in expected and digest(p)!=manifest[key]:errors.append('HASH '+key)
        if key=='SHA256SUMS.json':continue
        try:content=p.read_text(encoding='utf-8')
        except UnicodeError:errors.append('NON_UTF8 '+key);continue
        if key=='scripts/public_code_preflight.py':continue
        for name,regex in SECRET_RULES:
            if regex.search(content):errors.append(name+' '+key)
    if (ROOT/'LICENSE_PENDING.md').exists():
        errors.append('LICENSE_PENDING_STILL_PRESENT')
    lic = ROOT/'LICENSE'
    if not lic.is_file() or not lic.read_text(encoding='utf-8').startswith('MIT License\n\nCopyright (c) 2026 Elena Krikun'):
        errors.append('MIT_LICENSE_MISSING_OR_UNEXPECTED')
    return errors

def main():
    errors=audit()
    for err in errors[:30]: print('[HOLD]',err)
    if errors:
        print('[STATUS] HOLD_PUBLIC_CODE_V7_INTEGRITY',len(errors),'issues')
        return 2
    print('[STATUS] PASS_PUBLIC_CODE_V7_STATIC_INTEGRITY')
    print('[FILES]',len(all_files()),'code-and-doc files verified')
    print('[PUBLICATION_GATE] HOLD_MANUAL_THIRD_PARTY_CODE_PROVENANCE_AND_RELEASE_REVIEW')
    print('[LIMIT] Static scan is not a privacy, security or legal review')
    return 0
if __name__=='__main__':sys.exit(main())
