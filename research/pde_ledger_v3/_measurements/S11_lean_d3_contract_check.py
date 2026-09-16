#!/usr/bin/env python3
"""J4 bounded D3 invariants checks: one Lean process, durable logs, mathematical controls.

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
SCRATCH = LEAN/'s11/_scratch/d3_verification'
REPORT = BASE/'_measurements/S11_lean_d3_contract_checks.json'
ORDER = ('Quadratic','Rotation','Constraints','Classification','Census','Controls')
DEPENDENCIES = []
PREAMBLE = """import S11D3Invariants.Controls
open S11D3Invariants
noncomputable section
"""
EXAMPLES = [
 ('so_dimension', "theorem contract_control : Module.finrank ℝ soSpace = REL := by\n  norm_num [so_dimension]\n", '3','2'),
 ('o_dimension', "theorem contract_control : Module.finrank ℝ oSpace = REL := by\n  norm_num [o_dimension]\n", '3','4'),
 ('odd_dimension', "theorem contract_control : Module.finrank ℝ oddSpace = REL := by\n  norm_num [odd_dimension]\n", '0','1'),
 ('trace_required', "theorem contract_control : REL ∃ b c : ℝ, invariantForm ![1,0,0] = invariantForm ![0,b,c] := by\n  norm_num [trace_not_omittable]\n", '¬',''),
 ('trace_of_square_required', "theorem contract_control : REL ∃ a c : ℝ, invariantForm ![0,1,0] = invariantForm ![a,0,c] := by\n  norm_num [traceOfSquare_not_omittable]\n", '¬',''),
 ('frobenius_required', "theorem contract_control : REL ∃ a b : ℝ, invariantForm ![0,0,1] = invariantForm ![a,b,0] := by\n  norm_num [frobenius_not_omittable]\n", '¬',''),
 ('no_odd_extra', "theorem contract_control : REL ∃ Q : Quad, SOInvariant Q ∧ ReflectionOdd Q ∧ Q ≠ 0 := by\n  norm_num [nonzero_odd_impossible]\n", '¬',''),
 ('full_group_sensitivity', "theorem contract_control : REL SOInvariant (monomial 0 0) := by\n  norm_num [single_entry_not_invariant]\n", '¬',''),
]
MUTATIONS = [
 ('trace_normalization','Rotation','2 • monomial 0 4','3 • monomial 0 4','traceSquare_apply'),
 ('trace_of_square_sign','Rotation','2 • monomial 1 3','(-2 : ℝ) • monomial 1 3','traceOfSquare_apply'),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes():
    pending = [LEAN/'s11/S11D3Invariants.lean']; seen = set()
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
        if p.suffix=='.lean' and 'S11D3Invariants' not in str(p):
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
                           'process_timeout_seconds':300,'note':'Allocator limit is not an OS RSS cap.'}}

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
        same_inputs=all(prior.get('source_sha256_before',{}).get(rel)==digest for rel,digest in inputs.items())
        if output and same_inputs and old.get('command')==command and output.exists() and old.get('olean_sha256')==sha(output):
            record={**old,'reused_from_instrument_sha256':prior['instrument_sha256'],
                    'reuse_source_dependency_sha256':inputs}
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
        # No completed local proof modules are dependencies of this increment.
        # Reuse is still guarded by all source, pin, command and object hashes.
        for directory,namespace,modules in DEPENDENCIES:
            for module in modules:
                run('dependency_'+namespace+'_'+module, LEAN/directory/namespace/f'{module}.lean',
                    'PASS', output=LEAN/'.lake/build/lib/lean'/namespace/f'{module}.olean')
        report['dependency_olean_sha256_before_build']=dependency_oleans
        dependency_oleans={rel:sha(BASE/rel) for rel in dependency_oleans}
        report['dependency_olean_sha256']=dependency_oleans
        save()
        for module in ORDER:
            run('build_'+module,LEAN/'s11/S11D3Invariants'/f'{module}.lean','PASS',
                output=LEAN/'.lake/build/lib/lean/S11D3Invariants'/f'{module}.olean')
        audit=run('axiom_audit',LEAN/'s11/S11D3Invariants.lean','PASS',output=LEAN/'.lake/build/lib/lean/S11D3Invariants.olean')
        roots=re.findall(r'^#print axioms (\S+)',(LEAN/'s11/S11D3Invariants.lean').read_text(),re.M)
        axiom_rows=re.findall(r"^'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)$",audit['output'],re.M)
        assert [name for name,_ in axiom_rows]==roots
        assert all(set(a.replace(' ','').split(','))<={'propext','Classical.choice','Quot.sound',''} for _,a in axiom_rows)
        report['axiom_audit_count']=len(axiom_rows)
        report['axiom_declarations']={name:[x.strip() for x in a.split(',') if x.strip()] for name,a in axiom_rows}
        for path in (LEAN/'s11/S11D3Invariants').rglob('*.lean'):
            if '_scratch' not in path.parts:
                assert not re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)',path.read_text(),re.M),path
        for name,module,old,new,required in MUTATIONS:
            canonical=LEAN/'s11/S11D3Invariants'/f'{module}.lean';source=canonical.read_text()
            assert source.count(old)==1,(name,source.count(old))
            path=SCRATCH/(name+'.lean');path.write_text(source.replace(old,new))
            record=run(name,path,'REJECTED',required)
            record.update(canonical_source=str(canonical.relative_to(BASE)),replacement={'old':old,'new':new});save()
        for name,statement,correct,wrong in EXAMPLES:
            for suffix,value,expected in [('positive',correct,'PASS'),('mutant',wrong,'REJECTED')]:
                path=SCRATCH/(name+'_'+suffix+'.lean');path.write_text(PREAMBLE+statement.replace('REL',value))
                run(name+'_'+suffix,path,expected,'contract_control')
        for name,statement in [
            ('nonempty_positive','theorem contract_control : ∃ Q : Quad, SOInvariant Q ∧ Q ≠ 0 := nonzero_invariant_exists'),
            ('zero_positive','theorem contract_control : SOInvariant (0 : Quad) ∧ ReflectionOdd (0 : Quad) := zero_invariant'),
            ('negative_coefficients_positive','theorem contract_control : OInvariant (invariantForm ![-2,3,-5]) := invariantForm_O _'),
            ('unique_coefficients_positive','theorem contract_control : ∃! v : Fin 3 → ℝ, invariantForm ![-2,3,-5] = invariantForm v := SO_unique _ (invariantForm_SO _)')]:
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
