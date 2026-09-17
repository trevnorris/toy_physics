#!/usr/bin/env python3
"""D4C.1–D4C.4 full D4 bulk checks: one Lean process, durable logs, mathematical controls.

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
SCRATCH = LEAN/'s11/_scratch/d4_bulk_verification'
REPORT = BASE/'_measurements/S11_lean_d4_bulk_contract_checks.json'
ORDER = ('Even', 'Action', 'Variation', 'Bulk', 'Census', 'Controls')
DEPENDENCIES = []
PREAMBLE = 'import S11D4Bulk.Controls\nopen S11D4Bulk\nnoncomputable section\n'
EXAMPLES = [('response_dimension', 'theorem contract_control : Module.finrank ℝ (LinearMap.range responseMap) = REL := by\n  norm_num only [bulk_response_dimension]\n', '2', '3'), ('null_dimension', 'theorem contract_control : Module.finrank ℝ (LinearMap.ker responseMap) = REL := by\n  norm_num only [null_dimension]\n', '2', '1'), ('null_sum_sign', 'theorem contract_control : REL VariationallyNull ![1,1,0,0] := by\n  norm_num [variationallyNull_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]\n', '¬', ''), ('omit_laplacian', 'theorem contract_control : REL VariationallyNull ![0,0,1,0] := by\n  norm_num [variationallyNull_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]\n', '¬', ''), ('retain_odd_null', 'theorem contract_control : REL VariationallyNull ![0,0,0,7] := by\n  norm_num [variationallyNull_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]\n', '', '¬'), ('even_density_not_zero', 'theorem contract_control : lagrangian ![1,-1,0,0] evenJet = REL := by\n  norm_num only [even_null_density_nonzero]\n', '-1', '0'), ('even_momentum_sign', 'theorem contract_control : momentum ![1,0,0,0] evenJet 1 0 = REL := by\n  norm_num only [even_momentum_normalization]\n', '-2', '2'), ('even_momentum_factor', 'theorem contract_control : momentum ![1,0,0,0] evenJet 1 0 = REL := by\n  norm_num only [even_momentum_normalization]\n', '-2', '-4'), ('odd_density_factor', 'theorem contract_control : lagrangian ![0,0,0,2] S11D4Odd.witnessJet = REL := by\n  norm_num only [odd_density_normalization]\n', '-1', '-2'), ('odd_momentum_sign', 'theorem contract_control : momentum ![0,0,0,2] S11D4Odd.witnessJet 1 1 = REL := by\n  norm_num only [odd_momentum_normalization]\n', '-1', '1'), ('odd_momentum_factor', 'theorem contract_control : momentum ![0,0,0,2] S11D4Odd.witnessJet 1 1 = REL := by\n  norm_num only [odd_momentum_normalization]\n', '-1', '-2'), ('even_current_factor', 'theorem contract_control : (∑ i : Fin 4, S10Pilot.coordDeriv i.succ (fun y => boundaryCurrent affineWitness y i) 0) = REL := by\n  norm_num only [even_current_normalization]\n', '2', '1'), ('odd_current_factor', 'theorem contract_control : S11D4Odd.currentJet ![0,2,0,0] S11D4Odd.witnessJet 0 = REL := by\n  norm_num only [odd_current_normalization]\n', '1', '2'), ('transverse_coefficient', 'theorem contract_control : eulerLagrange ![2,3,7,11] (S10Pilot.planeWave 0 ![1,0,0,0] ![0,1,0,0]) 0 1 = REL := by\n  norm_num only [transverse_response]\n', '-7', '-5'), ('longitudinal_coefficient', 'theorem contract_control : eulerLagrange ![2,3,7,11] (S10Pilot.planeWave 0 ![1,0,0,0] ![1,0,0,0]) 0 0 = REL := by\n  norm_num only [longitudinal_response]\n', '-12', '-7'), ('bulk_equivalence', 'theorem contract_control : REL BulkEquivalent ![2,3,7,11] ![8,-3,7,-9] := by\n  norm_num [bulkEquivalent_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]\n', '', '¬')]
MUTATIONS = []
EXTRA_POSITIVES = [('mixed_null_positive', 'theorem contract_control (t beta : ℝ) : VariationallyNull ![t,-t,0,beta] := mixed_null_positive t beta'), ('zero_wavevector_positive', 'theorem contract_control (v a : S10Pilot.Vec 4) : modalOperator v 0 a = 0 := modal_zero_wavevector v a'), ('negative_coefficients_positive', 'theorem contract_control : BulkEquivalent ![-2,1,-7,3] ![5,-6,-7,-9] := by\n  norm_num [bulkEquivalent_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]'), ('nonzero_variation_positive', 'theorem contract_control : ∃ u h : S10Pilot.Point 4 → S10Pilot.Vec 4, S10Pilot.SmoothField u ∧ S10Pilot.TestField h ∧ deriv (relativeAction ![1,1,0,0] u h) 0 ≠ 0 :=\n  exists_nonzero_firstVariation _ (by norm_num [Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three])')]

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
    pending = [LEAN/'s11/S11D4Bulk.lean']; seen = set()
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
    seen |= {BASE/'_measurements'/name for name in ['S11_lean_d4_bulk_source_check.py','S11_lean_d4_bulk_source_checks.json','S11_lean_d4_bulk_preserved_inputs.json']}
    seen |= {BASE/rel for rel in historical_record()['read_only_inputs']}
    return {str(p.relative_to(BASE)):sha(p) for p in sorted(seen)}


def local_build_order():
    visited=set(); ordered=[]
    def visit(module):
        if module in visited:return
        visited.add(module)
        directory='s11' if module.startswith('S11') else 's10' if module.startswith('S10') else 's9'
        path=LEAN/directory/(module.replace('.','/')+'.lean')
        for child in re.findall(r'^import ((?:S9|S10|S11)\w+(?:\.\w+)*)$',path.read_text(),re.M):visit(child)
        if not module.startswith('S11D4Bulk'):
            namespace,leaf=module.split('.',1)
            ordered.append((directory,namespace,[leaf.replace('.','/')]))
    visit('S11D4Bulk')
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
    return json.loads((BASE/'_measurements/S11_lean_d4_bulk_preserved_inputs.json').read_text())

def validate_preserved():
    record=historical_record()
    for field in ('preserved_objects','read_only_inputs','historical_sources'):
        for rel,expected in record[field].items():assert sha(BASE/rel)==expected,rel
    return record

def validate_native():
    path=BASE/'_measurements/S11_lean_d4_bulk_source_checks.json'
    report=json.loads(path.read_text())
    assert report['status']=='PASS'
    assert report['instrument_sha256']==sha(BASE/'_measurements/S11_lean_d4_bulk_source_check.py')
    assert all(report['identities'].values()) and all(report['target_controls'].values())
    assert all(not c['V5_every_actual_basis_element'] and not c['actual_L_native_EL_equals_negative_Lean_EL'] for c in report['native_mutations'].values())
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
        if p.suffix=='.lean' and 'S11D4Bulk' not in str(p):
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
        seen.add(BASE/'_measurements/S11_lean_d4_generate.py')

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
            run('build_'+module,LEAN/'s11/S11D4Bulk'/f'{module}.lean','PASS',
                output=LEAN/'.lake/build/lib/lean/S11D4Bulk'/f'{module}.olean')
        audit=run('axiom_audit',LEAN/'s11/S11D4Bulk.lean','PASS',output=LEAN/'.lake/build/lib/lean/S11D4Bulk.olean')
        roots=re.findall(r'^#print axioms (\S+)',(LEAN/'s11/S11D4Bulk.lean').read_text(),re.M)
        axiom_rows=re.findall(r"^'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)$",audit['output'],re.M)
        assert [name for name,_ in axiom_rows]==roots
        assert all({x.strip() for x in a.split(',')}<={'propext','Classical.choice','Quot.sound',''} for _,a in axiom_rows)
        report['axiom_audit_count']=len(axiom_rows)
        report['axiom_declarations']={name:[x.strip() for x in a.split(',') if x.strip()] for name,a in axiom_rows}
        for path in (LEAN/'s11/S11D4Bulk').rglob('*.lean'):
            if '_scratch' not in path.parts:
                assert not re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)',path.read_text(),re.M),path
        for name,module,old,new,required in MUTATIONS:
            canonical=LEAN/'s11/S11D4Bulk'/f'{module}.lean';source=canonical.read_text()
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
