#!/usr/bin/env python3
"""D5.1–D5.4 quadratic invariant checks: one Lean process, durable logs, mathematical controls.

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
import signal
import subprocess
import time

BASE = Path(__file__).resolve().parents[1]
LEAN = BASE/'lean'
SCRATCH = LEAN/'s11/_scratch/d5_verification'
REPORT = BASE/'_measurements/S11_lean_d5_contract_checks.json'
ORDER = ('CoordinateAlgebra', 'BilinearExpansion', 'Quadratic', 'Rotation', 'Forms', 'ConstraintPolynomial', 'ConstraintBlock0', 'ConstraintBlock1', 'ConstraintBlock2', 'ConstraintBlock3', 'ConstraintBlock4', 'ConstraintBlock5', 'ConstraintBlock6', 'ConstraintBlock7', 'ConstraintBlock8', 'ConstraintBlock9', 'ConstraintBlock10', 'ConstraintBlock11', 'ConstraintBlock12', 'ConstraintBlock13', 'ConstraintBlock14', 'ConstraintBlock15', 'ConstraintBlock16', 'ReconstructionBlock0', 'ReconstructionBlock1', 'ReconstructionBlock2', 'ReconstructionBlock3', 'ReconstructionBlock4', 'ReconstructionBlock5', 'ReconstructionBlock6', 'ReconstructionBlock7', 'ReconstructionBlock8', 'ReconstructionBlock9', 'ReconstructionBlock10', 'ReconstructionBlock11', 'ReconstructionBlock12', 'ReconstructionBlock13', 'ReconstructionBlock14', 'ReconstructionBlock15', 'ReconstructionBlock16', 'Constraints', 'Classification', 'Census', 'Controls')
DEPENDENCIES = []
PREAMBLE = 'import S11D5Invariants.Controls\nopen S11D5Invariants\nnoncomputable section\n'
EXAMPLES = [('SO_count', 'theorem contract_control : Module.finrank ℝ soSpace = REL := by\n  norm_num only [so_dimension]\n', '3', '4'), ('O_count', 'theorem contract_control : Module.finrank ℝ oSpace = REL := by\n  norm_num only [o_dimension]\n', '3', '2'), ('odd_count', 'theorem contract_control : Module.finrank ℝ oddSpace = REL := by\n  norm_num only [odd_dimension]\n', '0', '1'), ('odd_existence', 'theorem contract_control : REL (∃ Q : Quad, SOInvariant Q ∧ ReflectionOdd Q ∧ Q ≠ 0) := by\n  norm_num only [nonzero_odd_impossible]\n', '¬', ''), ('omit_trace', 'theorem contract_control : REL (∃ b c : ℝ, invariantForm ![1,0,0] = invariantForm ![0,b,c]) := by\n  norm_num only [trace_not_omittable]\n', '¬', ''), ('omit_pair', 'theorem contract_control : REL (∃ a c : ℝ, invariantForm ![0,1,0] = invariantForm ![a,0,c]) := by\n  norm_num only [traceOfSquare_not_omittable]\n', '¬', ''), ('omit_frobenius', 'theorem contract_control : REL (∃ a b : ℝ, invariantForm ![0,0,1] = invariantForm ![a,b,0]) := by\n  norm_num only [frobenius_not_omittable]\n', '¬', ''), ('wrong_span_member', 'theorem contract_control : REL (SOInvariant (monomial 0 0)) := by\n  norm_num only [single_entry_not_invariant]\n', '¬', ''), ('trace_factor', 'theorem contract_control : traceSquare twoDiagonal = REL := by\n  norm_num only [trace_normalization]\n', '4', '2'), ('pair_factor', 'theorem contract_control : traceOfSquare offSymmetric = REL := by\n  norm_num only [pair_normalization]\n', '2', '1'), ('frobenius_factor', 'theorem contract_control : frobeniusSquare offDiagonal = REL := by\n  norm_num only [frobenius_normalization]\n', '1', '2'), ('reflection_sign', 'theorem contract_control : reflection.det = REL := by\n  norm_num only [reflection_det]\n', '-1', '1'), ('SO_O_equality', 'theorem contract_control : soSpace REL oSpace := by\n  norm_num only [so_eq_o, ne_self_iff_false]\n', '=', '≠')]
MUTATIONS = []
EXTRA_POSITIVES = [('zero_form_positive', 'theorem contract_control : SOInvariant (0 : Quad) ∧ ReflectionOdd (0 : Quad) := zero_invariant'), ('nonzero_form_positive', 'theorem contract_control : ∃ Q : Quad, SOInvariant Q ∧ Q ≠ 0 := nonzero_invariant_exists'), ('arbitrary_coefficients_positive', 'theorem contract_control (v : Fin 3 → ℝ) : OInvariant (invariantForm v) := invariantForm_O v'), ('unique_coefficients_positive', 'theorem contract_control (Q : Quad) (h : SOInvariant Q) : ∃! v : Fin 3 → ℝ, Q = invariantForm v := SO_unique Q h')]

def run_owned_process(command,stream,timeout_seconds=600):
    """Bound a launcher and all children to one process group and timeout."""
    process=subprocess.Popen(command,cwd=LEAN,stdout=stream,stderr=subprocess.STDOUT,
        start_new_session=True,
        env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
    try:
        process.wait(timeout=timeout_seconds)
    except subprocess.TimeoutExpired as error:
        try:
            os.killpg(process.pid,signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL)
            process.wait()
        # Stop a child even if the launcher exited first on SIGTERM.
        try:
            os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:
            pass
        error.owned_process_group=process.pid
        raise
    return subprocess.CompletedProcess(command,process.returncode)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes():
    pending = [LEAN/'s11/S11D5Invariants.lean']; seen = set()
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
    seen |= {BASE/'_measurements'/name for name in ['S11_lean_d5_source_check.py','S11_lean_d5_source_checks.json','S11_lean_d5_preserved_inputs.json']}
    seen |= {BASE/rel for rel in historical_record()['read_only_inputs']}
    seen.add(BASE/'_measurements/S11_lean_d5_generate.py')
    return {str(p.relative_to(BASE)):sha(p) for p in sorted(seen)}


def local_build_order():
    visited=set(); ordered=[]
    def visit(module):
        if module in visited:return
        visited.add(module)
        directory='s11' if module.startswith('S11') else 's10' if module.startswith('S10') else 's9'
        path=LEAN/directory/(module.replace('.','/')+'.lean')
        for child in re.findall(r'^import ((?:S9|S10|S11)\w+(?:\.\w+)*)$',path.read_text(),re.M):visit(child)
        if not module.startswith('S11D5Invariants'):
            namespace,leaf=module.split('.',1)
            ordered.append((directory,namespace,[leaf.replace('.','/')]))
    visit('S11D5Invariants')
    return ordered

def external_inputs():
    manifest=json.loads((LEAN/'lake-manifest.json').read_text())
    packages=[]
    for p in manifest['packages']:
        directory=LEAN/'.lake/packages'/p['name']
        if not directory.exists():
            directory=LEAN/'.lake/packages'/p['name'].strip('«»')
        assert directory.exists(),directory
        head=subprocess.check_output(['git','-C',str(directory),'rev-parse','HEAD'],text=True).strip()
        assert head==p['rev'],p['name']
        assert not subprocess.check_output(['git','-C',str(directory),'status','--porcelain','--untracked-files=no'],text=True),p['name']
        packages.append({'name':p['name'],'commit':head,'tracked_source_clean':True})
    objects={}
    sources={}
    for rel in source_hashes():
        path=BASE/rel
        if path.suffix!='.lean':continue
        for module in re.findall(r'^import (Mathlib(?:\.\w+)+)$',path.read_text(),re.M):
            stem=module.replace('.','/')
            source=LEAN/'.lake/packages/mathlib'/(stem+'.lean')
            obj=LEAN/'.lake/packages/mathlib/.lake/build/lib/lean'/(stem+'.olean')
            assert source.exists() and obj.exists(),module
            sources[str(source.relative_to(BASE))]=sha(source)
            objects[str(obj.relative_to(BASE))]=sha(obj)
    return {'package_pins':packages,'direct_import_source_sha256':sources,'direct_import_olean_sha256':objects}


def historical_record():
    return json.loads((BASE/'_measurements/S11_lean_d5_preserved_inputs.json').read_text())

def validate_preserved():
    record=historical_record()
    for field in ('preserved_objects','read_only_inputs','historical_sources'):
        for rel,expected in record[field].items():assert sha(BASE/rel)==expected,rel
    return record

def validate_native():
    path=BASE/'_measurements/S11_lean_d5_source_checks.json'
    report=json.loads(path.read_text())
    assert report['status']=='PASS'
    assert report['instrument_sha256']==sha(BASE/'_measurements/S11_lean_d5_source_check.py')
    assert all(group['same_rref'] for group in report['native_spans'].values())
    assert [report['native_spans'][k]['dimension'] for k in ['SO','O','odd']]==[3,3,0]
    assert report['odd_operator_identity'] and report['PD_zero']
    assert report['controls']['same_count_wrong_span_rejected'] and report['controls']['wrong_orientation_span_rejected']
    for rel,expected in report['source_sha256'].items():assert sha(BASE/rel)==expected,rel
    return {'report_sha256':sha(path),'instrument_sha256':report['instrument_sha256']}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--reuse-build',action='store_true')
    args = parser.parse_args()
    SCRATCH.mkdir(parents=True,exist_ok=True)
    historical = validate_preserved()
    native = validate_native()
    before = source_hashes()
    external = external_inputs()
    dependency_oleans = {}
    for rel in before:
        p=Path(rel)
        if p.suffix=='.lean' and 'S11D5Invariants' not in str(p):
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
              'historical_preservation':historical,'native_identification':native,'instrument_sha256':sha(Path(__file__)),'source_sha256_before':before,'dependency_olean_sha256':dependency_oleans,'checks':[],
              'external_dependencies':external,
              'resources':{'workers':1,'lean_threads':1,'lean_allocator_limit_mib':4096,
                           'process_timeout_seconds':600,'note':'Allocator limit is not an OS RSS cap.'}}

    def save():
        temporary = REPORT.with_suffix('.new')
        temporary.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
        temporary.replace(REPORT)

    def build_inputs(path):
        # A changed downstream proof need not invalidate already checked imports.
        # Require the full local transitive source closure and dependency pins
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
        seen.add(BASE/'_measurements/S11_lean_d5_generate.py')

        return {str(p.relative_to(BASE)):sha(p) for p in sorted(seen)}

    def run(name,path,expected,required=None,output=None):
        command=['nice','-n','10','lake','env','lean','-j1','-M4096','-DwarningAsError=true']
        if output:
            output.parent.mkdir(parents=True,exist_ok=True)
            command += ['-o',str(output)]
        rel=str(path.relative_to(LEAN));command.append(rel)
        old=reusable.get(name,{})
        inputs=build_inputs(path) if output else {}
        if output: inputs.update(external['direct_import_source_sha256'])
        input_objects=dict(external['direct_import_olean_sha256']) if output else {}
        if output:
            for source_rel in inputs:
                source_path=BASE/source_rel
                if source_path.suffix=='.lean' and source_path!=path and '.lake' not in source_path.parts:
                    parts=Path(source_rel).parts
                    obj=LEAN/'.lake/build/lib/lean'/Path(*parts[2:]).with_suffix('.olean')
                    assert obj.exists(),obj
                    input_objects[str(obj.relative_to(BASE))]=sha(obj)
        same_inputs=old.get('source_dependency_sha256')==inputs and prior.get('external_dependencies')==external
        if output and same_inputs and old.get('input_olean_sha256')==input_objects and old.get('command')==command and output.exists() and old.get('olean_sha256')==sha(output):
            record={**old,'reused_from_instrument_sha256':prior['instrument_sha256'],
                    'reuse_source_dependency_sha256':inputs}
            report['checks'].append(record);save();return record
        log=SCRATCH/(name+'.log');start=time.monotonic()
        with log.open('w') as stream:
            try:
                result=run_owned_process(command,stream)
            except subprocess.TimeoutExpired as error:
                report['timeout_cleanup']={'name':name,'process_group':error.owned_process_group,
                    'whole_group_stopped':True,'accepted_as_mutation':False}
                save()
                raise
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
            result.returncode==1 and len(diagnostics)==1 and len(intended)==1 and intended[0]['message'].count('⊢ False')==1 and not invalid and not re.search(r'\bwarning(?:\([^)]*\))?:',text))
        record={'name':name,'expected':expected,'outcome':expected if valid else 'UNEXPECTED',
                'command':command,'exit_status':result.returncode,'wall_seconds':round(time.monotonic()-start,3),
                'source_sha256':sha(path),'required_failure_in':required,'diagnostics':diagnostics,
                'output':text,'output_sha256':sha(log)}
        if output and result.returncode==0:
            record['olean_sha256']=sha(output)
            record['input_olean_sha256']=input_objects
            record['source_dependency_sha256']=inputs
        if expected=='REJECTED' or name.endswith('_positive'):record['source']=source
        report['checks'].append(record);save()
        if not valid:raise RuntimeError(f'{name}: unexpected outcome; inspect {log}')
        return record

    save()
    try:
        # Bind unchanged local imports to their source bytes sequentially.
        # Completed contracts and historical reports are not modified.
        for directory,namespace,modules in local_build_order():
            for module in modules:
                run('dependency_'+namespace+'_'+module, LEAN/directory/namespace/f'{module}.lean',
                    'PASS', output=LEAN/'.lake/build/lib/lean'/namespace/f'{module}.olean')
        report['dependency_olean_sha256_before_build']=dependency_oleans
        dependency_oleans={rel:sha(BASE/rel) for rel in dependency_oleans}
        report['dependency_olean_sha256']=dependency_oleans
        save()
        for module in ORDER:
            run('build_'+module,LEAN/'s11/S11D5Invariants'/f'{module}.lean','PASS',
                output=LEAN/'.lake/build/lib/lean/S11D5Invariants'/f'{module}.olean')
        audit=run('axiom_audit',LEAN/'s11/S11D5Invariants.lean','PASS',output=LEAN/'.lake/build/lib/lean/S11D5Invariants.olean')
        roots=re.findall(r'^#print axioms (\S+)',(LEAN/'s11/S11D5Invariants.lean').read_text(),re.M)
        axiom_rows=re.findall(r"^'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)$",audit['output'],re.M)
        assert [name for name,_ in axiom_rows]==roots
        assert all({x for x in ''.join(a.split()).split(',')}<={'propext','Classical.choice','Quot.sound',''} for _,a in axiom_rows)
        report['axiom_audit_count']=len(axiom_rows)
        report['axiom_declarations']={name:[x for x in ''.join(a.split()).split(',') if x] for name,a in axiom_rows}
        for path in (LEAN/'s11/S11D5Invariants').rglob('*.lean'):
            if '_scratch' not in path.parts:
                assert not re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)',path.read_text(),re.M),path
        for name,module,old,new,required in MUTATIONS:
            canonical=LEAN/'s11/S11D5Invariants'/f'{module}.lean';source=canonical.read_text()
            assert source.count(old)==1,(name,source.count(old))
            path=SCRATCH/(name+'.lean');path.write_text(source.replace(old,new))
            record=run(name,path,'REJECTED',required)
            record.update(canonical_source=str(canonical.relative_to(BASE)),replacement={'old':old,'new':new});save()
        for name,statement,correct,wrong in EXAMPLES:
            for suffix,value,expected in [('positive',correct,'PASS'),('mutant',wrong,'REJECTED')]:
                path=SCRATCH/(name+'_'+suffix+'.lean');path.write_text(PREAMBLE+statement.replace('REL',value))
                run(name+'_'+suffix,path,expected,'contract_control')
        for name,statement in EXTRA_POSITIVES:
            path=SCRATCH/(name+'.lean');path.write_text(PREAMBLE+statement+'\n')
            run(name,path,'PASS')
        assert historical==validate_preserved(),'historical object/source drift'
        assert native==validate_native(),'native identification drift'
        after=source_hashes();assert before==after,'source/dependency changed during checks'
        assert dependency_oleans=={rel:sha(BASE/rel) for rel in dependency_oleans}, 'dependency object changed during checks'
        assert external==external_inputs(),'external dependency changed during checks'
        report.update(status='PASS',source_sha256_after=after)
    except Exception as error:
        report.update(status='ERROR',error=repr(error))
        raise
    finally:
        report['finished_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':
    main()
