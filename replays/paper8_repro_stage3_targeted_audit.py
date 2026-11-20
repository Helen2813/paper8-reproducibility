#!/usr/bin/env python3
"""Paper 8 Stage 3. Read-only targeted reproducibility audit; writes only to Downloads.
No computations, raw data reading, installs, network or changes to research project.
Python >=3.9, standard library only.
"""
import ast
import datetime as dt
import hashlib
import json
import os
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path.home() / 'Desktop' / 'paper8' / 'endo-research'
OUT = Path.home() / 'Downloads' / 'Paper8_repro_stage3_targeted'
STAGES = {
 'RMA': ['paper8_v350b_fogg_original_rma.py','paper8_v349z_reproduction_environment.py'],
 'ComBat': ['paper8_v350e_original_combat.py','paper8_v350j_factorial_combat_dedup.py','paper8_v352b_combat_execution_gate.py'],
 'label_and_null': ['paper8_v352c_single_label_intervention.py','paper8_v352e_joint_permutations.py','paper8_v354d_full_panel_nolabel_convergence.py'],
 'Huang': ['paper8_v348j_leave_one_out.py','paper8_v350i_deduplication_sensitivity.py'],
 'locks': ['paper8_v349f_model_spec_lock.py','paper8_v352a_lock_protocol.py','paper8_v352d_frozen_score.py'],
}
RESULT_GLOBS = [
 'p8_v350b_original_rma*.zip','p8_v350g_core_refit*.zip','p8_v352c_single_label*.zip',
 'p8_v352e_joint_permutations*.zip','p8_v428_huang*.zip','p8_v429_huang*.zip',
 'p8_v370*.zip','p8_v371*.zip','p8_v4*.zip',
]
SOURCE_FILES = [
 'data/external/p8_v334_gse51981_cel/GSE51981_RAW.tar',
 'data/external/p8_v347f_gse6364_cel/GSE6364_RAW.tar',
 'data/external/p8_v401_GEO_sources/GSE25628_RAW.tar',
 'data/p8_cross_series_independence/fogg_original_raw/GSE29981_RAW.tar',
 'data/p8_cross_series_independence/fogg_original_raw/GSE4888_RAW.tar',
 'data/p8_cross_series_independence/fogg_rma_original_commit/counts/GSE51981_counts_unmapped.tsv',
 'data/p8_cross_series_independence/fogg_code/df538583652d3359d71959437837ea3e0d1195e6/data_prep/rma_counts.r',
 'data/p8_cross_series_independence/fogg_code/df538583652d3359d71959437837ea3e0d1195e6/data_prep/batch_correction.r',
]
RE_NAME = re.compile(r'(seed|random_state|n_perm|n_reps|repetitions|alpha|lambda|combat|rma|glmnet|randomforest|sklearn|cv|fold|manifest|sha256|frozen|input|result|output|run_id)', re.I)
RE_VERSION = re.compile(r'(v3(?:50|52|54|65|69|70|71)|v4(?:01|07|22|28|29|30)|huang|324|nested|date|null)',re.I)

def file_meta(p):
    try:
        s=p.stat()
        return {'path':str(p),'exists':True,'bytes':s.st_size}
    except (OSError, ValueError):
        return {'path':str(p),'exists':False}

def safe_text(p,limit=180000):
    try:
        if p.stat().st_size>limit: return None
        return p.read_text(encoding='utf-8',errors='replace')
    except (OSError,UnicodeError): return None

def analyze_script(p):
    rec=file_meta(p)
    if not rec['exists']:return rec
    txt=safe_text(p)
    if txt is None: rec['inspection']='SKIPPED_LARGE_OR_UNREADABLE';return rec
    rec['sha256']=hashlib.sha256(txt.encode('utf-8',errors='replace')).hexdigest()
    imports=[]
    if p.suffix=='.py':
        try:
            tree=ast.parse(txt)
            for node in ast.walk(tree):
                if isinstance(node,ast.Import):imports.extend(n.name for n in node.names)
                elif isinstance(node,ast.ImportFrom):imports.append(node.module or '')
            rec['syntax']='PASS'
        except SyntaxError as e:rec['syntax']=f'SYNTAX_ERROR_LINE_{e.lineno}'
    rec['imports']=sorted(set(imports))[:120]
    lines=txt.splitlines()
    rec['relevant_lines']=[{'number':i,'text':line[:320]} for i,line in enumerate(lines,1)
                           if RE_NAME.search(line) and not line.lstrip().startswith('#')][:45]
    return rec

def inspect_zip(p):
    rec=file_meta(p)
    if not rec['exists']:return rec
    try:
        with zipfile.ZipFile(p) as z:
            names=z.namelist()
            rec['entries']=len(names)
            rec['selected_entries']=[{'name':x.filename,'bytes':x.file_size} for x in z.infolist()
                     if RE_VERSION.search(x.filename) or x.filename.lower().endswith(('.json','.r','.py','.md'))][:90]
            rec['bad_member_paths']=[n for n in names if n.startswith('/') or '..' in Path(n).parts][:10]
            # Directory inspection only; do not extract or test raw archives.
    except (OSError,zipfile.BadZipFile) as e:rec['zip_error']=str(e)
    return rec

def main():
    if not ROOT.is_dir():
        sys.exit('HOLD: research root not found: '+str(ROOT))
    scripts={key:[analyze_script(ROOT/x) for x in names] for key,names in STAGES.items()}
    # Only inspect project-root script filenames; no broad scan of private data.
    discovered=sorted((x for x in ROOT.iterdir() if x.is_file() and x.suffix in ('.py','.R','.r') and RE_VERSION.search(x.name)),key=lambda x:x.name)
    extra=[analyze_script(x) for x in discovered if x.name not in {n for ns in STAGES.values() for n in ns}][:80]
    result_dir=ROOT/'results'
    found=[]
    if result_dir.is_dir():
        seen=set()
        for pat in RESULT_GLOBS:
            for p in sorted(result_dir.glob(pat)):
                if p.is_file() and str(p) not in seen:
                    seen.add(str(p));found.append(p)
    # Limit archive inspection to the most specific, compact results, avoid large payload reads.
    zip_recs=[inspect_zip(x) for x in found[:85]]
    sources=[file_meta(ROOT/x) for x in SOURCE_FILES]
    environments=[]
    for rel in ['.envs/p8_fogg_r44/bin/Rscript','.envs/p8_null/bin/python3',
                '.envs/p8_fogg_classifier/bin/python3','.envs/p8_v422_py/bin/python3']:
        environments.append(file_meta(ROOT/rel))
    report={
      'status':'TARGETED_INVENTORY_ONLY',
      'created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
      'policy':'NO raw contents, analysis, installs, downloads or writes in research project',
      'root':str(ROOT),'stages':scripts,'additional_root_scripts':extra,
      'source_and_original_code_paths':sources,'result_archives':zip_recs,
      'result_archive_count_discovered':len(found),'environment_executables':environments,
      'limitations':['Script source text is inspected, not executed',
                     'ZIP central directories are listed, their contents not executed or extracted',
                     'Archived scripts and outputs not proven to be final inputs without provenance checks',
                     'Full data hashes and reruns not attempted']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'targeted_inventory.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    summary=['PAPER 8 | STAGE 3 | TARGETED READ-ONLY PROVENANCE AUDIT',
             '[STATUS] TARGETED_INVENTORY_ONLY','[ROOT] '+str(ROOT),
             '[POLICY] No raw reads, training, downloads, installs, or research-directory writes','']
    for stage,items in scripts.items():
        summary.append(f'[{stage}] '+', '.join(f'{Path(x["path"]).name}:{"FOUND" if x["exists"] else "MISSING"}' for x in items))
    summary += ['','[ADDITIONAL PROJECT-ROOT SCRIPTS] '+str(len(extra)),
               '[RESULT ARCHIVES INSPECTED] '+str(len(zip_recs)),
               '[RESULT ARCHIVES DISCOVERED] '+str(len(found)),
               '[INPUT FILES FOUND] '+str(sum(x['exists'] for x in sources))+'/'+str(len(sources)),
               '[ENVIRONMENT EXECUTABLES FOUND] '+str(sum(x['exists'] for x in environments))+'/'+str(len(environments)),
               '[NEXT] Share targeted_inventory.json privately; do not publish unreviewed absolute paths.',
               '[REPORT] '+str(OUT/'targeted_inventory.json')]
    (OUT/'summary.txt').write_text('\n'.join(summary)+'\n',encoding='utf8')
    print('\n'.join(summary))

if __name__=='__main__':main()
