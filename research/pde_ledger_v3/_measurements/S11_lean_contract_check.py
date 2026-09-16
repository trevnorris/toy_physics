#!/usr/bin/env python3
"""H4 bounded S11 checks: one Lean process, durable logs, mathematical controls.

No CAS production runs. Reuse successful builds only with identical commands,
transitive local sources, dependency pins and output-module hashes. Mutations
always run afresh. Resource/import/tactic-environment failures never count.
"""
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

BASE = Path(__file__).resolve().parents[1]
LEAN = BASE/'lean'
SCRATCH = LEAN/'s11/_scratch/verification'
REPORT = BASE/'_measurements/S11_lean_contract_checks.json'
ORDER = ('Action', 'Spectrum', 'Threshold')
PREAMBLE = '''import S11Homogeneous.Threshold
open S10Pilot S11Homogeneous
noncomputable section
abbrev k3 : Vec 3 := ![0, 0, 1]
abbrev e1 : Vec 3 := ![1, 0, 0]
abbrev e3 : Vec 3 := ![0, 0, 1]
theorem k3_ne : k3 ≠ 0 := by
  intro h
  have h2 := congrFun h 2
  norm_num [k3, Matrix.cons_val_two] at h2
'''
COUNT = '''theorem contract_control :
    Module.finrank ℝ (S11Homogeneous.modeSpace 1 1 BULK OMEGA k3) = REL := by
  rw [kernel_census_three (by norm_num) k3_ne]
  norm_num [coneValue, normSq, dot, k3, Fin.sum_univ_succ]
'''
MODE = '''theorem contract_control : REL S11Homogeneous.ModalStationary 1 1 BULK OMEGA WAVE AMP := by
  norm_num [S11Homogeneous.modal_stationary_iff, S11Homogeneous.modalOperator,
    normSq, dot, k3, e1, e3, Fin.sum_univ_succ, funext_iff, Fin.forall_fin_succ]
'''
THRESHOLD = '''theorem contract_control : normalWaveSq 1 BULK 2 k3 REL 0 := by
  norm_num [normalWaveSq, normSq, dot, k3, Fin.sum_univ_succ]
'''
EXAMPLES = [
    ('transverse_count', COUNT.replace('BULK','4').replace('OMEGA','1'), '2', '3'),
    ('longitudinal_count', COUNT.replace('BULK','4').replace('OMEGA','2'), '1', '2'),
    ('coincidence_count', COUNT.replace('BULK','1').replace('OMEGA','1'), '3', '2'),
    ('transverse_independent_of_compression', MODE.replace('BULK','9').replace('OMEGA','1').replace('WAVE','k3').replace('AMP','e1'), '', '¬'),
    ('longitudinal_root_is_not_shear_root', MODE.replace('BULK','4').replace('OMEGA','1').replace('WAVE','k3').replace('AMP','e3'), '¬', ''),
    ('off_roots_excluded', MODE.replace('BULK','4').replace('OMEGA','3').replace('WAVE','k3').replace('AMP','e1'), '¬', ''),
    ('zero_compression_static_mode', MODE.replace('BULK','0').replace('OMEGA','0').replace('WAVE','k3').replace('AMP','e3'), '', '¬'),
    ('zero_wavevector_static_mode', MODE.replace('BULK','4').replace('OMEGA','0').replace('WAVE','0').replace('AMP','e1'), '', '¬'),
    ('zero_wavevector_no_dynamic_mode', MODE.replace('BULK','4').replace('OMEGA','1').replace('WAVE','0').replace('AMP','e1'), '¬', ''),
    ('below_grazing', THRESHOLD.replace('BULK','1'), '<', '>'),
    ('at_grazing', THRESHOLD.replace('BULK','4'), '=', '≠'),
    ('above_grazing', THRESHOLD.replace('BULK','9'), '>', '<'),
    ('grazing_requires_nonzero_wavevector', '''theorem contract_control :
    normalWaveSq 1 1 2 (0 : Vec 3) = 0 ∧ (1 : ℝ) REL 1 * 2 ^ 2 := by
  norm_num [normalWaveSq, normSq, dot]
''', '≠', '='),
]
MUTATIONS = [
    ('compression_action_sign', 'Action',
     'mu / 2 * stiffness J -\n    B / 2', 'mu / 2 * stiffness J +\n    B / 2',
     'lagrangian_split'),
    ('phase_average_half', 'Action',
     'phaseAverage rho mu B omega k a = (1 / 2 : ℝ) * modalAction',
     'phaseAverage rho mu B omega k a = (1 : ℝ) * modalAction', 'phaseAverage_eq'),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes():
    pending = [LEAN/'s11/S11Homogeneous.lean']; seen = set()
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        for module in re.findall(r'^import ((?:S9|S10|S11)\w+(?:\.\w+)*)$',path.read_text(),re.M):
            directory = 's11' if module.startswith('S11') else 's10' if module.startswith('S10') else 's9'
            child = LEAN/directory/(module.replace('.','/')+'.lean')
            assert child.exists(),child
            pending.append(child)
    seen |= {LEAN/name for name in ['lakefile.toml','lake-manifest.json','lean-toolchain']}
    return {str(p.relative_to(BASE)):sha(p) for p in sorted(seen)}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--reuse-build',action='store_true')
    args = parser.parse_args()
    SCRATCH.mkdir(parents=True,exist_ok=True)
    before = source_hashes()
    prior = json.loads(REPORT.read_text()) if REPORT.exists() else {}
    if prior:
        archive = SCRATCH/('prior_'+sha(REPORT)+'.json')
        if not archive.exists():
            archive.write_bytes(REPORT.read_bytes())
    reusable = {r['name']:r for r in prior.get('checks',[]) if r['outcome']=='PASS'} if (
        args.reuse_build and before==prior.get('source_sha256_before')) else {}
    report = {'status':'RUNNING','started_utc':datetime.now(timezone.utc).isoformat(),
              'instrument_sha256':sha(Path(__file__)),'source_sha256_before':before,'checks':[],
              'resources':{'workers':1,'lean_threads':1,'lean_allocator_limit_mib':4096,
                           'process_timeout_seconds':300,'note':'Allocator limit is not an OS RSS cap.'}}

    def save():
        temporary = REPORT.with_suffix('.new')
        temporary.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
        temporary.replace(REPORT)

    def run(name,path,expected,required=None,output=None):
        command=['nice','-n','10','lake','env','lean','-j1','-M4096','-DwarningAsError=true']
        if output:
            output.parent.mkdir(parents=True,exist_ok=True)
            command += ['-o',str(output)]
        rel=str(path.relative_to(LEAN));command.append(rel)
        old=reusable.get(name,{})
        if output and old.get('command')==command and output.exists() and old.get('olean_sha256')==sha(output):
            record={**old,'reused_from_instrument_sha256':prior['instrument_sha256']}
            report['checks'].append(record);save();return record
        log=SCRATCH/(name+'.log');start=time.monotonic()
        with log.open('w') as stream:
            result=subprocess.run(command,cwd=LEAN,stdout=stream,stderr=subprocess.STDOUT,timeout=300,
                                  env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
        text=log.read_text();source=path.read_text()
        decls=[(source[:m.start()].count('\n')+1,m[1]) for m in re.finditer(r'^\s*(?:theorem|def)\s+(\w+)',source,re.M)]
        errors=list(re.finditer(re.escape(rel)+r':(\d+):\d+: error(?:\([^)]*\))?:',text));diagnostics=[]
        for i,error in enumerate(errors):
            line=int(error[1]);decl=next((n for pos,n in reversed(decls) if pos<=line),'')
            message=text[error.start():errors[i+1].start() if i+1<len(errors) else len(text)]
            diagnostics.append({'declaration':decl,'line':line,'message':message})
        invalid=re.search(r'unknown (?:module|namespace|identifier)|unexpected token|maximum (?:recursion|heartbeats)|'
                          r'excessive memory|out of memory|PANIC|No such file|failed to create',text,re.I)
        # Concrete false claims must reduce to False. Source mutations must expose
        # an unsolved mathematical identity in the named declaration, inspected at review.
        intended=[d for d in diagnostics if d['declaration']==required and
                  ('⊢ False' in d['message'] if required=='contract_control' else 'unsolved goals' in d['message'])]
        valid=(result.returncode==0 and not re.search(r'\b(?:error|warning)(?:\([^)]*\))?:',text)) if expected=='PASS' else (
            result.returncode==1 and bool(intended) and not invalid)
        record={'name':name,'expected':expected,'outcome':expected if valid else 'UNEXPECTED',
                'command':command,'exit_status':result.returncode,'wall_seconds':round(time.monotonic()-start,3),
                'source_sha256':sha(path),'required_failure_in':required,'diagnostics':diagnostics,
                'output':text,'output_sha256':sha(log)}
        if output and result.returncode==0:record['olean_sha256']=sha(output)
        if expected=='REJECTED' or name.endswith('_positive'):record['source']=source
        report['checks'].append(record);save()
        if not valid:raise RuntimeError(f'{name}: unexpected outcome; inspect {log}')
        return record

    save()
    try:
        for module in ORDER:
            run('build_'+module,LEAN/'s11/S11Homogeneous'/f'{module}.lean','PASS',
                output=LEAN/'.lake/build/lib/lean/S11Homogeneous'/f'{module}.olean')
        audit=run('axiom_audit',LEAN/'s11/S11Homogeneous.lean','PASS',output=LEAN/'.lake/build/lib/lean/S11Homogeneous.olean')
        axioms=re.findall(r'depends on axioms: \[([^\]]*)\]',audit['output'])
        assert len(axioms)==29
        assert all(set(a.replace(' ','').split(','))<={'propext','Classical.choice','Quot.sound',''} for a in axioms)
        report['axiom_audit_count']=len(axioms)
        for path in (LEAN/'s11').rglob('*.lean'):
            if '_scratch' not in path.parts:
                assert not re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)',path.read_text(),re.M),path
        for name,module,old,new,required in MUTATIONS:
            canonical=LEAN/'s11/S11Homogeneous'/f'{module}.lean';source=canonical.read_text()
            assert source.count(old)==1,(name,source.count(old))
            path=SCRATCH/(name+'.lean');path.write_text(source.replace(old,new))
            record=run(name,path,'REJECTED',required)
            record.update(canonical_source=str(canonical.relative_to(BASE)),replacement={'old':old,'new':new});save()
        for name,statement,correct,wrong in EXAMPLES:
            for suffix,value,expected in [('positive',correct,'PASS'),('mutant',wrong,'REJECTED')]:
                path=SCRATCH/(name+'_'+suffix+'.lean');path.write_text(PREAMBLE+statement.replace('REL',value))
                run(name+'_'+suffix,path,expected,'contract_control')
        after=source_hashes();assert before==after,'source/dependency changed during checks'
        report.update(status='PASS',source_sha256_after=after)
    except Exception as error:
        report.update(status='ERROR',error=repr(error))
        raise
    finally:
        report['finished_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':
    main()
