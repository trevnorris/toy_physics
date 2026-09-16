#!/usr/bin/env python3
"""K4 bounded D3 bulk variation checks: one Lean process, durable logs, mathematical controls.

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
SCRATCH = LEAN/'s11/_scratch/d3_bulk_verification'
REPORT = BASE/'_measurements/S11_lean_d3_bulk_contract_checks.json'
ORDER = ('Action','Calculus','Variation','Bulk','Census')
DEPENDENCIES = [('s10', 'S10Pilot', ('Action',)),
 ('s10', 'S10Pilot', ('PlaneWave',)),
 ('s10', 'S10Pilot', ('Analytic',)),
 ('s10', 'S10Pilot', ('Variation',)),
 ('s11', 'S11D3Invariants', ('Quadratic',)),
 ('s11', 'S11D3Invariants', ('Rotation',)),
 ('s11', 'S11D3Invariants', ('Constraints',)),
 ('s11', 'S11D3Invariants', ('Classification',)),
 ('s10', 'S10Controls', ('Action',)),
 ('s10', 'S10Controls', ('PlaneWave',)),
 ('s10', 'S10Controls', ('Analytic',)),
 ('s10', 'S10Controls', ('Variation',)),
 ('s10', 'S10Pilot', ('PhaseAverage',)),
 ('s10', 'S10Controls', ('PhaseAverage',)),
 ('s11', 'S11Homogeneous', ('Action',))]
PREAMBLE = """import S11D3Bulk.Census
open S10Pilot S11D3Bulk
noncomputable section
"""
EXAMPLES = [('zero_null',
  'theorem contract_control : REL VariationallyNull 0 := by\n  norm_num [variationallyNull_iff, Matrix.cons_val_two]\n',
  '',
  '¬'),
 ('nonzero_null',
  'theorem contract_control : REL VariationallyNull ![1,-1,0] := by\n'
  '  norm_num [variationallyNull_iff, Matrix.cons_val_two]\n',
  '',
  '¬'),
 ('frobenius_not_null',
  'theorem contract_control : REL VariationallyNull ![0,0,1] := by\n'
  '  norm_num [variationallyNull_iff, Matrix.cons_val_two]\n',
  '¬',
  ''),
 ('trace_equivalence',
  'theorem contract_control : REL BulkEquivalent ![1,0,0] ![0,1,0] := by\n'
  '  norm_num [bulkEquivalent_iff, Matrix.cons_val_two]\n',
  '',
  '¬'),
 ('different_c',
  'theorem contract_control : REL BulkEquivalent ![1,0,0] ![1,0,1] := by\n'
  '  norm_num [bulkEquivalent_iff, Matrix.cons_val_two]\n',
  '¬',
  ''),
 ('response_count',
  'theorem contract_control : Module.finrank ℝ (LinearMap.range responseMap) = REL := by\n'
  '  norm_num [bulk_response_dimension]\n',
  '2',
  '3'),
 ('null_count',
  'theorem contract_control : Module.finrank ℝ (LinearMap.ker responseMap) = REL := by\n'
  '  norm_num [null_dimension]\n',
  '1',
  '0'),
 ('zero_wavevector',
  'theorem contract_control : S11D3Bulk.modalOperator ![1,2,3] 0 ![1,0,0] 0 = REL := by\n'
  '  norm_num [modal_zero_wavevector]\n',
  '0',
  '1')]
MUTATIONS = [('density_sign',
  'Action',
  '  -(1 / 2 : ℝ) * (v 0 * divergence J ^ 2',
  '  (1 / 2 : ℝ) * (v 0 * divergence J ^ 2',
  'density_identity'),
 ('density_half',
  'Action',
  '  -(1 / 2 : ℝ) * (v 0 * divergence J ^ 2',
  '  -(1 : ℝ) * (v 0 * divergence J ^ 2',
  'density_identity'),
 ('bulk_coefficient_map',
  'Bulk',
  '(-(v 0 + v 1) * dot k a) • k',
  '(-(v 0 - v 1) * dot k a) • k',
  'homogeneous_operator'),
 ('boundary_current_sign',
  'Bulk',
  'u x i * fieldJet u x j.succ j - u x j * fieldJet u x j.succ i',
  'u x i * fieldJet u x j.succ j + u x j * fieldJet u x j.succ i',
  'boundary_identity')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes():
    pending = [LEAN/'s11/S11D3Bulk.lean']; seen = set()
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
    seen.add(BASE/'_measurements/S11_lean_d3_generate.py')
    return {str(p.relative_to(BASE)):sha(p) for p in sorted(seen)}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--reuse-build',action='store_true')
    args = parser.parse_args()
    SCRATCH.mkdir(parents=True,exist_ok=True)
    subprocess.run(["python3",str(BASE/"_measurements/S11_lean_d3_generate.py"),"--check"],check=True,stdout=subprocess.DEVNULL)
    before = source_hashes()
    dependency_oleans = {}
    for rel in before:
        p=Path(rel)
        if p.suffix=='.lean' and 'S11D3Bulk' not in str(p):
            module=Path(*p.parts[2:]).with_suffix('.olean')
            obj=LEAN/'.lake/build/lib/lean'/module
            assert obj.exists(),obj
            dependency_oleans[str(obj.relative_to(BASE))]=sha(obj)
    prior = json.loads(REPORT.read_text()) if REPORT.exists() else {}
    if prior:
        archive = SCRATCH/('prior_'+sha(REPORT)+'.json')
        if not archive.exists():
            archive.write_bytes(REPORT.read_bytes())
    reusable = {r['name']:r for r in prior.get('checks',[]) if r['outcome']=='PASS'} if (
        args.reuse_build and dependency_oleans==prior.get('dependency_olean_sha256')) else {}
    report = {'status':'RUNNING','started_utc':datetime.now(timezone.utc).isoformat(),
              'instrument_sha256':sha(Path(__file__)),'source_sha256_before':before,'dependency_olean_sha256':dependency_oleans,'checks':[],
              'resources':{'workers':1,'lean_threads':1,'lean_allocator_limit_mib':4096,
                           'process_timeout_seconds':600,'note':'Allocator limit is not an OS RSS cap.'}}

    def save():
        temporary = REPORT.with_suffix('.new')
        temporary.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
        temporary.replace(REPORT)

    def build_inputs(path):
        # A changed downstream proof need not invalidate already checked imports.
        # Require the full local transitive source closure, pins and generator
        # for each reused object, and verify objects in dependency order below.
        pending=[path];seen=set()
        while pending:
            item=pending.pop()
            if item in seen:continue
            seen.add(item)
            for module in re.findall(r'^import ((?:S9|S10|S11)\w+(?:\.\w+)*)$',item.read_text(),re.M):
                directory='s11' if module.startswith('S11') else 's10' if module.startswith('S10') else 's9'
                pending.append(LEAN/directory/(module.replace('.','/')+'.lean'))
        seen|={LEAN/n for n in ['lakefile.toml','lake-manifest.json','lean-toolchain']}
        seen.add(BASE/'_measurements/S11_lean_d3_generate.py')
        return {str(p.relative_to(BASE)):sha(p) for p in sorted(seen)}

    def run(name,path,expected,required=None,output=None):
        command=['nice','-n','10','lake','env','lean','-j1','-M4096','-DwarningAsError=true']
        if output:
            output.parent.mkdir(parents=True,exist_ok=True)
            command += ['-o',str(output)]
        rel=str(path.relative_to(LEAN));command.append(rel)
        old=reusable.get(name,{})
        inputs=build_inputs(path) if output else {}
        input_objects={}
        if output:
            for source_rel in inputs:
                source_path=BASE/source_rel
                if source_path.suffix=='.lean' and source_path!=path:
                    parts=Path(source_rel).parts
                    obj=LEAN/'.lake/build/lib/lean'/Path(*parts[2:]).with_suffix('.olean')
                    assert obj.exists(),obj
                    input_objects[str(obj.relative_to(BASE))]=sha(obj)
        same_inputs=all(prior.get('source_sha256_before',{}).get(rel)==digest for rel,digest in inputs.items())
        if output and same_inputs and old.get('input_olean_sha256')==input_objects and old.get('command')==command and output.exists() and old.get('olean_sha256')==sha(output):
            record={**old,'reused_from_instrument_sha256':prior['instrument_sha256'],
                    'reuse_source_dependency_sha256':inputs}
            report['checks'].append(record);save();return record
        log=SCRATCH/(name+'.log');start=time.monotonic()
        with log.open('w') as stream:
            result=subprocess.run(command,cwd=LEAN,stdout=stream,stderr=subprocess.STDOUT,timeout=600,
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
        if output and result.returncode==0:
            record['olean_sha256']=sha(output)
            record['input_olean_sha256']=input_objects
        if expected=='REJECTED' or name.endswith('_positive'):record['source']=source
        report['checks'].append(record);save()
        if not valid:raise RuntimeError(f'{name}: unexpected outcome; inspect {log}')
        return record

    save()
    try:
        # Bind unchanged local imports to their source bytes sequentially.
        # Completed contracts and historical reports are not modified.
        for directory,namespace,modules in DEPENDENCIES:
            for module in modules:
                run('dependency_'+namespace+'_'+module, LEAN/directory/namespace/f'{module}.lean',
                    'PASS', output=LEAN/'.lake/build/lib/lean'/namespace/f'{module}.olean')
        report['dependency_olean_sha256_before_build']=dependency_oleans
        dependency_oleans={rel:sha(BASE/rel) for rel in dependency_oleans}
        report['dependency_olean_sha256']=dependency_oleans
        save()
        for module in ORDER:
            run('build_'+module,LEAN/'s11/S11D3Bulk'/f'{module}.lean','PASS',
                output=LEAN/'.lake/build/lib/lean/S11D3Bulk'/f'{module}.olean')
        audit=run('axiom_audit',LEAN/'s11/S11D3Bulk.lean','PASS',output=LEAN/'.lake/build/lib/lean/S11D3Bulk.olean')
        roots=re.findall(r'^#print axioms (\S+)',(LEAN/'s11/S11D3Bulk.lean').read_text(),re.M)
        axiom_rows=re.findall(r"^'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)$",audit['output'],re.M)
        assert [name for name,_ in axiom_rows]==roots
        assert all(set(a.replace(' ','').split(','))<={'propext','Classical.choice','Quot.sound',''} for _,a in axiom_rows)
        report['axiom_audit_count']=len(axiom_rows)
        report['axiom_declarations']={name:[x.strip() for x in a.split(',') if x.strip()] for name,a in axiom_rows}
        for path in (LEAN/'s11/S11D3Bulk').rglob('*.lean'):
            if '_scratch' not in path.parts:
                assert not re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)',path.read_text(),re.M),path
        for name,module,old,new,required in MUTATIONS:
            canonical=LEAN/'s11/S11D3Bulk'/f'{module}.lean';source=canonical.read_text()
            assert source.count(old)==1,(name,source.count(old))
            path=SCRATCH/(name+'.lean');path.write_text(source.replace(old,new))
            record=run(name,path,'REJECTED',required)
            record.update(canonical_source=str(canonical.relative_to(BASE)),replacement={'old':old,'new':new});save()
        for name,statement,correct,wrong in EXAMPLES:
            for suffix,value,expected in [('positive',correct,'PASS'),('mutant',wrong,'REJECTED')]:
                path=SCRATCH/(name+'_'+suffix+'.lean');path.write_text(PREAMBLE+statement.replace('REL',value))
                run(name+'_'+suffix,path,expected,'contract_control')
        for name,statement in [('admissible_compact_variation_positive', 'theorem contract_control : ∃ u h : Point 3 → Vec 3, SmoothField u ∧ TestField h ∧ deriv (S11D3Bulk.relativeAction ![-2,3,-5] u h) 0 ≠ 0 := exists_nonzero_firstVariation _ (by norm_num)'), ('longitudinal_positive', 'theorem contract_control : S11D3Bulk.modalOperator ![1,2,3] ![1,0,0] ![1,0,0] = ![-6,0,0] := by\n  ext i\n  fin_cases i <;> norm_num [modal_longitudinal, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two]'), ('transverse_positive', 'theorem contract_control : S11D3Bulk.modalOperator ![1,2,3] ![1,0,0] ![0,1,0] = ![0,-3,0] := by\n  ext i\n  fin_cases i <;> norm_num [S11D3Bulk.modalOperator, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two]')]:
            path=SCRATCH/(name+'.lean');path.write_text(PREAMBLE+statement+'\n')
            run(name,path,'PASS')
        after=source_hashes();assert before==after,'source/dependency changed during checks'
        assert dependency_oleans=={rel:sha(BASE/rel) for rel in dependency_oleans}, 'dependency object changed during checks'
        report.update(status='PASS',source_sha256_after=after)
    except Exception as error:
        report.update(status='ERROR',error=repr(error))
        raise
    finally:
        report['finished_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':
    main()
