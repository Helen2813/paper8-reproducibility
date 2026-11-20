#!/usr/bin/env python3
"""Paper 8 Stage 4: non-destructive provenance verification, no analysis runs.

Uses previously generated Stage 3 targeted_inventory.json. Writes reports only
under ~/Downloads/Paper8_repro_stage4_gate. Requires Python >=3.9, stdlib only.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import zipfile

VERSION = 'stage4-v1.0'
MAX_JSON_MEMBER = 2_000_000
MAX_MEMBER_SCANS = 1000
OUTPUT_NAME = 'Paper8_repro_stage4_gate'
PRIORITIES = (
    ('RMA', 'p8_v350b_original_rma_'),
    ('ComBat', 'p8_v350e_original_combat_'),
    ('Single-label ComBat', 'p8_v352c_single_label_20261108T032638Z'),
    ('999 ComBat permutations (NOT 3582 refits)', 'p8_v352e_joint_permutations_999'),
    ('Null 199-refit blocks', 'p8_v371_full_gene_combat_null_'),
    ('Huang deletion', 'p8_v428_huang_5gene_duplicate_sensitivity_'),
    ('Huang ablation', 'p8_v429_huang_training_ablation_'),
    ('Nested 324-probe', 'p8_v422_locked_GSE51981_MAS5_nested_cv_'),
    ('Date comparison', 'p8_v424_scanner_date_negative_controls_'),
    ('Last targeted evidence', 'p8_v430_targeted_reanalysis_'),
)
SCRIPT_TOKENS = {
    'RMA':'paper8_v350b_fogg_original_rma.py',
    'Probe mapping':'paper8_v350c_original_probe_mapping.py',
    'Metadata reconciliation':'paper8_v350d_mapping_metadata_reconciliation.py',
    'ComBat':'paper8_v350e_original_combat.py',
    'Classifier matrices':'paper8_v350f_classifier_inputs.py',
    'Core classifier':'paper8_v350g_core_matrisome_refit.py',
    'Frozen score':'paper8_v352d_frozen_score.py',
    'Single-label intervention':'paper8_v352c_single_label_intervention.py',
    'Joint label permutations':'paper8_v352e_joint_permutations.py',
    'Factorial comparison':'paper8_v352g_factorial.py',
    'Huang reconstruction':'paper8_v348j_leave_one_out.py',
}

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as fp:
        for chunk in iter(lambda:fp.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def safe_inside(root,p):
    try: return p.resolve().is_relative_to(root.resolve())
    except AttributeError:
        try: p.resolve().relative_to(root.resolve()); return True
        except ValueError: return False

def json_member(z, name):
    inf=z.getinfo(name)
    if inf.file_size > MAX_JSON_MEMBER: return {'_status':'TOO_LARGE'}
    try: return json.loads(z.read(name).decode('utf-8-sig'))
    except (ValueError,UnicodeError,RuntimeError) as e: return {'_status':'INVALID_JSON','_error':str(e)[:100]}

def summarize(obj):
    if not isinstance(obj,dict): return {'type':type(obj).__name__}
    result={}
    for k,v in obj.items():
        if k.lower() in ('status','task','version','analysis','method','n_permutations','n_models','n_refits',
            'n_samples','n_cel','n_fits','n_replicates','seed','random_seed','n_success','n_failed',
            'n_test','n_train','n_complete','n_checked','n_predictions','n_models_run',
            'conclusion','overall_status','replay_status'):
            if isinstance(v,(str,int,float,bool)):
                result[k]=v if not isinstance(v,str) else v[:180]
    return result

def archive_inspect(p):
    result={'name':p.name,'exists':p.is_file()}
    if not p.is_file():return result
    try:
        with zipfile.ZipFile(p) as z:
            infos=z.infolist()
            result['entries']=len(infos)
            result['unsafe_members']=[i.filename for i in infos if i.filename.startswith('/') or '..' in Path(i.filename).parts][:5]
            manifests=[]
            for inf in infos[:MAX_MEMBER_SCANS]:
                nm=inf.filename
                base=Path(nm).name.lower()
                if base.endswith('.json') and any(t in base for t in ('manifest','report','gate','summary','status','lock','protocol')):
                    doc=json_member(z,nm)
                    manifests.append({'file':nm,'summary':summarize(doc),'top_keys':list(doc)[:16] if isinstance(doc,dict) else []})
            result['manifests']=manifests[:35]
            result['unread_manifests']=max(0,len(manifests)-35)
    except (OSError,zipfile.BadZipFile,RuntimeError) as e:
        result['error']=str(e)[:220]
    return result

def package_probe(path):
    try:
        res=subprocess.run([str(path),'-c',
            'import importlib.util,json,sys; m=["numpy","pandas","scipy","sklearn","skopt","yaml"]; print(json.dumps({"python":sys.version.split()[0],"modules":{x:importlib.util.find_spec(x) is not None for x in m}}))'],
            capture_output=True,text=True,timeout=15)
        if res.returncode:return {'status':'PROBE_FAILED','stderr':res.stderr[-400:]}
        return {'status':'OK',**json.loads(res.stdout)}
    except Exception as e:return {'status':'PROBE_ERROR','message':str(e)[:200]}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory',default=str(Path.home()/'Downloads/Paper8_repro_stage3_targeted/targeted_inventory.json'))
    parser.add_argument('--output',default=str(Path.home()/'Downloads'/OUTPUT_NAME))
    parser.add_argument('--selftest',action='store_true')
    args=parser.parse_args()
    if args.selftest:
        assert summarize({'status':'PASS','x':[1]})=={'status':'PASS'}
        assert not safe_inside(Path('/tmp/xyzroot'),Path('/tmp/other'))
        print('[SELFTEST] PASS');return 0
    inv_path=Path(args.inventory).expanduser()
    if not inv_path.is_file():
        print('[HOLD] INVENTORY_NOT_FOUND:',inv_path)
        return 2
    inv=json.loads(inv_path.read_text(encoding='utf-8'))
    if inv.get('status')!='TARGETED_INVENTORY_ONLY':
        print('[HOLD] UNEXPECTED_STAGE3_STATUS');return 2
    root=Path(inv['root'])
    if not root.is_dir():
        print('[HOLD] PROJECT_ROOT_UNAVAILABLE:',root);return 2
    output=Path(args.output).expanduser()
    if safe_inside(root,output):
        print('[HOLD] OUTPUT_INSIDE_RESEARCH_DIRECTORY');return 2
    print('PAPER 8 | STAGE 4 | SOURCE-AND-ARCHIVE PROVENANCE GATE')
    print('[POLICY] NO model runs, no RMA/ComBat, no data changes, no downloads, no installs')
    print('[INPUT]',inv_path)
    report={'version':VERSION,'status':'INSPECTION_ONLY','created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
            'source_inventory':str(inv_path),'root':str(root),'scripts':[], 'archive_summaries':[],
            'environment_probes':[],'warnings':[]}
    indexed={}
    for section in inv['stages'].values():
        for item in section:indexed[Path(item['path']).name]=item
    for item in inv['additional_root_scripts']:indexed[Path(item['path']).name]=item
    for stage,filename in SCRIPT_TOKENS.items():
        entry=indexed.get(filename)
        p=root/filename
        record={'stage':stage,'name':filename,'exists':p.is_file(),
                'expected_sha256':entry.get('sha256') if entry else None,'syntax':entry.get('syntax') if entry else 'NOT_IN_STAGE3'}
        if p.is_file():
            record['actual_sha256']=digest(p)
            record['sha256_status']='MATCH' if record['expected_sha256'] and record['actual_sha256']==record['expected_sha256'] else 'NO_EXPECTED_HASH_OR_MISMATCH'
            if record['sha256_status']!='MATCH':report['warnings'].append('SCRIPT_HASH_NOT_MATCHED '+filename)
        else:report['warnings'].append('SCRIPT_MISSING '+filename)
        report['scripts'].append(record)
        print('[SCRIPT]',stage,record.get('sha256_status','MISSING'))
    archives=[Path(x['path']) for x in inv['result_archives']]
    for stage,prefix in PRIORITIES:
        selected=[p for p in archives if p.name.startswith(prefix)]
        selected.sort()
        if not selected:
            report['warnings'].append('ARCHIVE_NOT_FOUND '+stage)
            report['archive_summaries'].append({'stage':stage,'status':'MISSING'})
            print('[ARCHIVE]',stage,'MISSING')
            continue
        for p in selected:
            if not safe_inside(root,p):
                report['warnings'].append('ARCHIVE_OUTSIDE_PROJECT '+p.name);continue
            row=archive_inspect(p);row['stage']=stage
            report['archive_summaries'].append(row)
            print('[ARCHIVE]',stage,p.name,'manifests=',len(row.get('manifests',[])),'error=',row.get('error','none'))
    for env in inv['environment_executables']:
        p=Path(env['path'])
        if p.name.startswith('python'):
            x={'path':str(p),'exists':p.is_file(),'probe':package_probe(p) if p.is_file() else None}
            report['environment_probes'].append(x)
            print('[PYTHON]',p.parent.parent.name,x['probe'].get('status') if x['probe'] else 'MISSING')
    # Explicit publication distinctions: do not equate the 999 label-permutation runs with 3582 model fits.
    report['interpretation']={
        'stage3_does_not_prove_models_reran':True,
        '999_joint_combat_permutations_is_distinct_from_3582_refits':True,
        '3582_comprises_18_conditions_times_199_refits_from_archived_publication':True,
        'full_raw_cel_replay_not_attempted':True,
        'no_experimental_results_claimed':True,
    }
    output.mkdir(parents=True,exist_ok=True)
    (output/'stage4_provenance.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    lines=['PAPER 8 - STAGE 4 PROVENANCE', 'STATUS: INSPECTION_ONLY',
           'Scripts SHA matching: '+str(sum(x.get('sha256_status')=='MATCH' for x in report['scripts']))+'/'+str(len(report['scripts'])),
           'Matched archive records: '+str(len(report['archive_summaries'])),
           'Warnings: '+str(len(report['warnings'])),
           'NOT A RECOMPUTATION: 999 ComBat permutations differ from 3582 classifier refits.',
           'No original research data or scripts were modified.',
           'NEXT: inspect validated chains and choose ONE bounded computation in a sandbox.']
    if report['warnings']:lines.extend(['WARNINGS:']+report['warnings'])
    (output/'summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('[OUTPUT]',output/'stage4_provenance.json')
    print('[STATUS] STAGE4_INSPECTION_ONLY')
    return 0

if __name__=='__main__':sys.exit(main())
