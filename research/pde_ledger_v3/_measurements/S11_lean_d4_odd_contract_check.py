#!/usr/bin/env python3
"""D4B.4 bounded constant-coefficient D4 odd variation checks: one Lean process, durable logs, mathematical controls.

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
SCRATCH = LEAN/'s11/_scratch/d4_odd_verification'
REPORT = BASE/'_measurements/S11_lean_d4_odd_contract_checks.json'
ORDER = ('Action','Calculus','Variation','Boundary')
DEPENDENCIES = [('s10','S10Pilot',['Action','PlaneWave','Analytic']),
 ('s11','S11D4Invariants',['Quadratic','Rotation','Orientation','Forms','ConstraintPolynomial','ConstraintBlock0','ConstraintBlock1','ConstraintBlock2','ConstraintBlock3','Constraints','Classification'])]
PREAMBLE = 'import S11D4Odd.Boundary\nopen S10Pilot S11D4Odd\nnoncomputable section\nopen scoped ContDiff\n'
EXAMPLES = [('nonzero_density', 'theorem contract_control : orientationJet witnessJet = REL := by\n  norm_num [nonzero_density]\n', '1', '0'), ('momentum_sign', 'theorem contract_control : S11D4Odd.momentum 2 witnessJet 1 1 = REL := by\n  norm_num [nonzero_momentum]\n', '-1', '1'), ('current_factor', 'theorem contract_control : currentJet ![0,2,0,0] witnessJet 0 = REL := by\n  norm_num [current_normalization]\n', '1', '2'), ('degree_two_contraction', 'theorem contract_control : (∑ i : Fin 4, ∑ j : Fin 4, witnessJet i.succ j * dualCurl witnessJet i j) = REL := by\n  norm_num [dualCurl_contraction, nonzero_density]\n', '2', '1'), ('first_variation_zero', 'theorem contract_control {u h : Point 4 → Vec 4} (hu : SmoothField u) (hh : TestField h) : deriv (S11D4Odd.relativeAction 3 u h) 0 = REL := by\n  norm_num [firstVariation_zero hu hh]\n', '0', '1'), ('local_EL_zero', 'theorem contract_control {u : Point 4 → Vec 4} (hu : SmoothField u) (x : Point 4) : S11D4Odd.eulerLagrange 3 u x 0 = REL := by\n  norm_num [eulerLagrange_zero hu]\n', '0', '1'), ('negative_beta_density', 'theorem contract_control : S11D4Odd.lagrangian (-2) witnessJet = REL := by\n  norm_num [witness_lagrangian]\n', '1', '-1'), ('time_momentum_zero', 'theorem contract_control (J : Jet 4) : S11D4Odd.momentum 3 J 0 1 = REL := by\n  norm_num [S11D4Odd.momentum_eq]\n', '0', '1')]
MUTATIONS = [('density_normalization', 'Action', '-(beta / 2) * orientationJet J', '-beta * orientationJet J', 'density_identity'), ('complementary_orientation_sign', 'Action', '![-antisym J 2 3, 0, antisym J 0 3, -antisym J 0 2]', '![antisym J 2 3, 0, antisym J 0 3, -antisym J 0 2]', 'momentum_eq'), ('mixed_derivative_index', 'Boundary', '(∑ i : Fin 4, coordDeriv i.succ (fun y => dualCurl (fieldJet u y) i j) x) = 0 := by', '(∑ i : Fin 4, coordDeriv (1 : Fin 5) (fun y => dualCurl (fieldJet u y) i j) x) = 0 := by', 'dualCurl_divergence_zero'), ('boundary_current_factor', 'Boundary', '(1 / 2 : ℝ) * ∑ j : Fin 4, v j * dualCurl J i j', '(1 : ℝ) * ∑ j : Fin 4, v j * dualCurl J i j', 'boundary_identity')]
EXTRA_POSITIVES = [('zero_field_admissible_positive', 'theorem contract_control : TestField (fun _ : Point 4 => (0 : Vec 4)) := by\n  constructor\n  · intro i\n    fun_prop\n  · intro i\n    simp [HasCompactSupport]\n'), ('negative_beta_stationary_positive', 'theorem contract_control {u : Point 4 → Vec 4} (hu : SmoothField u) : S11D4Odd.ActionStationary (-7) u := every_background_stationary hu _\n'), ('zero_beta_stationary_positive', 'theorem contract_control {u : Point 4 → Vec 4} (hu : SmoothField u) : S11D4Odd.ActionStationary 0 u := every_background_stationary hu _\n'), ('smooth_affine_background_positive', 'theorem contract_control : SmoothField (fun x : Point 4 => (![0,x 1,0,x 3] : Vec 4)) := by\n  intro i\n  fin_cases i <;> simp [Matrix.cons_val_two, Matrix.cons_val_three] <;> fun_prop\n')]


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
    pending = [LEAN/'s11/S11D4Odd.lean']; seen = set()
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
    seen.add(BASE/'_measurements/S11_lean_d4_generate.py')
    return {str(p.relative_to(BASE)):sha(p) for p in sorted(seen)}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--reuse-build',action='store_true')
    args = parser.parse_args()
    SCRATCH.mkdir(parents=True,exist_ok=True)
    subprocess.run(["python3",str(BASE/"_measurements/S11_lean_d4_generate.py"),"--check"],check=True,stdout=subprocess.DEVNULL)
    before = source_hashes()
    dependency_oleans = {}
    for rel in before:
        p=Path(rel)
        if p.suffix=='.lean' and 'S11D4Odd' not in str(p):
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
            run('build_'+module,LEAN/'s11/S11D4Odd'/f'{module}.lean','PASS',
                output=LEAN/'.lake/build/lib/lean/S11D4Odd'/f'{module}.olean')
        audit=run('axiom_audit',LEAN/'s11/S11D4Odd.lean','PASS',output=LEAN/'.lake/build/lib/lean/S11D4Odd.olean')
        roots=re.findall(r'^#print axioms (\S+)',(LEAN/'s11/S11D4Odd.lean').read_text(),re.M)
        axiom_rows=re.findall(r"^'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)$",audit['output'],re.M)
        assert [name for name,_ in axiom_rows]==roots
        assert all(set(a.replace(' ','').split(','))<={'propext','Classical.choice','Quot.sound',''} for _,a in axiom_rows)
        report['axiom_audit_count']=len(axiom_rows)
        report['axiom_declarations']={name:[x.strip() for x in a.split(',') if x.strip()] for name,a in axiom_rows}
        for path in (LEAN/'s11/S11D4Odd').rglob('*.lean'):
            if '_scratch' not in path.parts:
                assert not re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)',path.read_text(),re.M),path
        for name,module,old,new,required in MUTATIONS:
            canonical=LEAN/'s11/S11D4Odd'/f'{module}.lean';source=canonical.read_text()
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
