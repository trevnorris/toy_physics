#!/usr/bin/env python3
"""D5B.1–D5B.4: isolated, sequential formal verification with guarded D5 dependency reuse.

The seed is a completed fresh portable D5 replay, validated against all live
inputs, external dependencies, logs and objects. No historical local object
is trusted or overwritten. Only dependencies are seeded; every new statement
and paired control is checked in this contract. Resource failures never count.
"""
import argparse
from datetime import datetime,timezone
import hashlib,json,os,re,shutil,subprocess,sys,time
from pathlib import Path
BASE=Path(__file__).resolve().parents[1];LEAN=BASE/'lean'
sys.path.insert(0,str(LEAN))
import verify as portable
SCRATCH=LEAN/'s11/_scratch/d5_bulk_verification'
OBJECTS=SCRATCH/'objects'
REPORT=BASE/'_measurements/S11_lean_d5_bulk_contract_checks.json'
PREAMBLE='import S11D5Bulk.Controls\nopen S11D5Bulk\nnoncomputable section\n'
EXAMPLES = [('response_dimension', 'theorem contract_control : Module.finrank ℝ (LinearMap.range responseMap) = REL := by\n  norm_num only [bulk_response_dimension]\n', '2', '3'), ('null_dimension', 'theorem contract_control : Module.finrank ℝ (LinearMap.ker responseMap) = REL := by\n  norm_num only [null_dimension]\n', '1', '2'), ('null_sum_sign', 'theorem contract_control : REL VariationallyNull ![1,1,0] := by\n  norm_num [variationallyNull_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]\n', '¬', ''), ('omit_laplacian', 'theorem contract_control : REL VariationallyNull ![0,0,1] := by\n  norm_num [variationallyNull_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]\n', '¬', ''), ('even_density_not_zero', 'theorem contract_control : lagrangian ![1,-1,0] evenJet = REL := by\n  norm_num only [even_null_density_nonzero]\n', '-1', '0'), ('even_momentum_sign', 'theorem contract_control : momentum ![1,0,0] evenJet 1 0 = REL := by\n  norm_num only [even_momentum_normalization]\n', '-2', '2'), ('even_momentum_factor', 'theorem contract_control : momentum ![1,0,0] evenJet 1 0 = REL := by\n  norm_num only [even_momentum_normalization]\n', '-2', '-4'), ('even_current_factor', 'theorem contract_control : (∑ i : Fin 5, S10Pilot.coordDeriv i.succ (fun y => boundaryCurrent affineWitness y i) 0) = REL := by\n  norm_num only [even_current_normalization]\n', '2', '1'), ('transverse_coefficient', 'theorem contract_control : eulerLagrange ![2,3,7] (S10Pilot.planeWave 0 ![1,0,0,0,0] ![0,1,0,0,0]) 0 1 = REL := by\n  norm_num only [transverse_response]\n', '-7', '-5'), ('longitudinal_coefficient', 'theorem contract_control : eulerLagrange ![2,3,7] (S10Pilot.planeWave 0 ![1,0,0,0,0] ![1,0,0,0,0]) 0 0 = REL := by\n  norm_num only [longitudinal_response]\n', '-12', '-7'), ('bulk_equivalence', 'theorem contract_control : REL BulkEquivalent ![2,3,7] ![8,-3,7] := by\n  norm_num [bulkEquivalent_iff, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three]\n', '', '¬'), ('time_momentum', 'theorem contract_control : momentum ![2,3,7] evenJet 0 4 = REL := by\n  norm_num only [momentum_time]\n', '0', '1'), ('fifth_transverse', 'theorem contract_control : eulerLagrange ![2,3,7] (S10Pilot.planeWave 0 ![0,0,0,0,1] ![1,0,0,0,0]) 0 0 = REL := by\n  norm_num only [fifth_transverse_response]\n', '-7', '0'), ('fifth_longitudinal', 'theorem contract_control : eulerLagrange ![2,3,7] (S10Pilot.planeWave 0 ![0,0,0,0,1] ![0,0,0,0,1]) 0 4 = REL := by\n  norm_num only [fifth_longitudinal_response]\n', '-12', '0')]
EXTRA_POSITIVES = [('null_positive', 'theorem contract_control (t : ℝ) : VariationallyNull ![t,-t,0] := null_positive t'), ('zero_wavevector_positive', 'theorem contract_control (v : Coeff) (a : S10Pilot.Vec 5) : modalOperator v 0 a = 0 := modal_zero_wavevector v a'), ('negative_coefficients_positive', 'theorem contract_control : BulkEquivalent ![-2,1,-7] ![5,-6,-7] := by\n  norm_num [bulkEquivalent_iff, Matrix.cons_val, Matrix.cons_val_two]'), ('nonzero_variation_positive', 'theorem contract_control : ∃ u h : S10Pilot.Point 5 → S10Pilot.Vec 5, S10Pilot.SmoothField u ∧ S10Pilot.TestField h ∧ deriv (relativeAction ![1,1,0] u h) 0 ≠ 0 :=\n  exists_nonzero_firstVariation _ (by norm_num [Matrix.cons_val, Matrix.cons_val_two])')]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def source_path(module):return BASE/portable.module_path(module)
def obj_path(module):return OBJECTS/(module.replace('.','/')+'.olean')
def closure(module):return portable.build_order([module])
def historical_record():return load(BASE/'_measurements/S11_lean_d5_bulk_preserved_inputs.json')
def validate_preserved():
    record=historical_record()
    for field in ['historical_sources','preserved_objects','read_only_inputs']:
        for rel,h in record[field].items():assert sha(BASE/rel)==h,rel
    assert sha(LEAN/'INSTALL_VALIDATION.json')==record['historical_install_validation_sha256']
    return record

def validate_native():
    path=BASE/'_measurements/S11_lean_d5_bulk_source_checks.json';report=load(path)
    assert report['status']=='PASS'
    assert report['instrument_sha256']==sha(BASE/'_measurements/S11_lean_d5_bulk_source_check.py')
    assert all(report['identities'].values()) and all(report['target_controls'].values())
    assert all(not r['V5_every_actual_basis_element'] and not r['actual_L_native_EL_equals_negative_Lean_EL'] for r in report['native_mutations'].values())
    for rel,h in report['source_sha256'].items():assert sha(BASE/rel)==h,rel
    return {'report_sha256':sha(path),'instrument_sha256':report['instrument_sha256']}

def external_inputs(env,metadata,order):
    direct={};sources={};objects={}
    manifest=load(LEAN/'lake-manifest.json')
    for module in order:
        for child in portable.IMPORT.findall(source_path(module).read_text()):
            if portable.LOCAL.fullmatch(child):continue
            stem=child.replace('.','/')
            obj=next((Path(p)/(stem+'.olean') for p in env['LEAN_PATH'].split(os.pathsep) if (Path(p)/(stem+'.olean')).is_file()),None)
            assert obj is not None,child
            source=next((LEAN/'.lake/packages'/p['name'].strip('«»')/(stem+'.lean') for p in manifest['packages'] if (LEAN/'.lake/packages'/p['name'].strip('«»')/(stem+'.lean')).is_file()),None)
            assert source is not None,child
            sources[str(source.relative_to(BASE))]=sha(source)
            objects[str(obj.relative_to(BASE))]=sha(obj)
            direct[child]={'source':str(source),'source_sha256':sha(source),'path':str(obj),'olean_sha256':sha(obj)}
    return {**metadata,'direct_imports':direct,'direct_import_source_sha256':sources,'direct_import_olean_sha256':objects}

def validate_seed(seed_path,external):
    seed=load(seed_path);run=seed_path.parent
    assert seed['status']=='PASS' and seed['plan']['contracts']==['d5'] and seed['plan']['mode']=='full'
    assert len(seed['objects'])==45 and seed['plan']['control_executions']==30
    assert sum(c.get('axiom_audits',0) for c in seed['checks'])==51
    assert seed['dependencies']['packages']==external['packages']
    assert seed['dependencies']['lean_version']==external['lean_version']
    for rel,h in seed['input_sha256'].items():
        assert sha(BASE/rel)==h,rel
        # The source-check output is newly produced only in the portable snapshot.
        if not rel.endswith('_source_checks.json'):assert sha(run/'snapshot'/rel)==h,rel
    for module,item in seed['dependencies']['direct_imports'].items():
        assert sha(item['path'])==item['olean_sha256'],module
        assert module in external['direct_imports'],module
        assert item['olean_sha256']==external['direct_imports'][module]['olean_sha256']
        assert item['source_sha256']==external['direct_imports'][module]['source_sha256']
    for record in seed['checks']:
        assert record['passed'],record['name']
        assert sha(run/record['log'])==record['log_sha256']
    for rel,h in seed['objects'].items():assert sha(run/'objects'/rel)==h,rel
    controls=[c for c in seed['checks'] if c['name'].startswith('d5_')]
    assert len(controls)==30 and sum(c['exit_status']==1 for c in controls)==13
    return seed

def source_hashes(order,seed_path):
    paths={source_path(m) for m in order}
    paths|={LEAN/n for n in ['lakefile.toml','lake-manifest.json','lean-toolchain','verify.py']}
    paths|={BASE/'_measurements'/n for n in ['S11_lean_d5_generate.py','S11_lean_d5_bulk_source_check.py','S11_lean_d5_bulk_source_checks.json','S11_lean_d5_bulk_preserved_inputs.json']}
    paths|={BASE/p for p in historical_record()['read_only_inputs']}
    paths.add(seed_path)
    return {str(p.relative_to(BASE)):sha(p) for p in sorted(paths)}

def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--seed-report',type=Path,required=True)
    parser.add_argument('--reuse-build',action='store_true')
    args=parser.parse_args();seed_path=args.seed_report.resolve()
    assert seed_path.is_relative_to(BASE/'_scratch/lean_portable'),seed_path
    SCRATCH.mkdir(parents=True,exist_ok=True);OBJECTS.mkdir(exist_ok=True)
    order=closure('S11D5Bulk')
    env,executable,metadata=portable.external_environment()
    external=external_inputs(env,metadata,order)
    seed=validate_seed(seed_path,external)
    preserved=validate_preserved();native=validate_native()
    before=source_hashes(order,seed_path)
    # Exclude all shared ledger objects; only these isolated objects and pinned external caches.
    env['LEAN_PATH']=str(OBJECTS)+os.pathsep+env['LEAN_PATH']
    env.update(LEAN_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    env.pop('PYTHONOPTIMIZE',None)
    old=load(REPORT) if REPORT.exists() else {}
    if old:
        archive=SCRATCH/('prior_'+sha(REPORT)+'.json')
        if not archive.exists():archive.write_bytes(REPORT.read_bytes())
    reusable={r['name']:r for r in old.get('checks',[]) if r['outcome']=='PASS'} if args.reuse_build else {}
    result={'status':'RUNNING','started_utc':datetime.now(timezone.utc).isoformat(),
            'instrument_sha256':sha(Path(__file__)),'source_sha256_before':before,
            'external_dependencies':external,'historical_preservation':preserved,
            'native_identification':native,'dependency_seed':{'report':str(seed_path.relative_to(BASE)),'sha256':sha(seed_path),'policy':'42 needed classification objects copied from the completed portable D5 execution; no proof skipped without its source/dependency/log/object evidence.'},
            'checks':[],'resources':{'workers':1,'lean_threads':1,'lean_allocator_limit_mib':4096,'process_timeout_seconds':600},
            'limits':['Paired controls only; no canonical source-replacement mutants.','External pinned cache baseline; no full Mathlib rebuild.','All local outputs isolated; historical shared objects remain untouched.']}
    def save():
        p=REPORT.with_suffix('.new');p.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n');p.replace(REPORT)
    def build_inputs(module):
        inputs={str(source_path(m).relative_to(BASE)):sha(source_path(m)) for m in closure(module)}
        for rel in ['lean/lakefile.toml','lean/lake-manifest.json','lean/lean-toolchain','_measurements/S11_lean_d5_generate.py','lean/verify.py']:
            inputs[rel]=sha(BASE/rel)
        inputs.update(external['direct_import_source_sha256'])
        objects={str(obj_path(m).relative_to(BASE)):sha(obj_path(m)) for m in closure(module) if m!=module}
        objects.update(external['direct_import_olean_sha256'])
        return inputs,objects
    def execute(name,path,expected='PASS',module=None):
        source=path.read_text();out=obj_path(module) if module else None
        command=[executable,'-j1','-M4096','-DwarningAsError=true']
        if out:
            out.parent.mkdir(parents=True,exist_ok=True);command+=['-o',str(out)]
        command.append(str(path.relative_to(LEAN)))
        inputs,objects=build_inputs(module) if module else ({},{})
        prior=reusable.get(name,{})
        same=(prior.get('command')==command and prior.get('source_dependency_sha256')==inputs and prior.get('input_olean_sha256')==objects and old.get('external_dependencies')==external)
        if out and same and out.exists() and sha(out)==prior.get('olean_sha256'):
            assert sha(SCRATCH/prior['log'])==prior['output_sha256']
            record={**prior,'reused_from_instrument_sha256':old['instrument_sha256']}
            result['checks'].append(record);save();return record
        log=SCRATCH/(name+'.log');start=time.monotonic()
        code=portable.run_process(command,LEAN,log,600,env)
        text=log.read_text()
        valid=(code==0 and not re.search(r'\b(?:error|warning)(?:\([^)]*\))?:',text)) if expected=='PASS' else portable.adjudicate(code,text,source,['contract_control'])
        record={'name':name,'expected':expected,'outcome':expected if valid else 'UNEXPECTED','command':command,'exit_status':code,
                'wall_seconds':round(time.monotonic()-start,3),'source_sha256':sha(path),'log':log.name,'output':text,'output_sha256':sha(log)}
        if module and valid:
            assert out.is_file(),out
            record.update(olean_sha256=sha(out),output_olean=str(out.relative_to(BASE)),source_dependency_sha256=inputs,input_olean_sha256=objects)
            record['axiom_audits']=portable.audit_axioms(text,source)
        if not module:
            record.update(source=source,required_failure_in='contract_control')
        result['checks'].append(record);save()
        if not valid:raise RuntimeError(f'{name}: unexpected outcome; inspect {log}')
        return record
    save()
    try:
        for module in order:
            if module.startswith('S11D5Invariants'):
                # Copy only after its transitive inputs are already present/checked in order.
                rel=module.replace('.','/')+'.olean';src=seed_path.parent/'objects'/rel;dst=obj_path(module)
                prior=next(c for c in seed['checks'] if c['name']=='build_'+module)
                inputs,objects=build_inputs(module)
                assert sha(src)==prior['olean_sha256']==seed['objects'][rel]
                dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
                record={'name':'seed_'+module,'expected':'PASS','outcome':'PASS','kind':'validated_dependency_seed',
                        'source_sha256':sha(source_path(module)),'source_dependency_sha256':inputs,'input_olean_sha256':objects,
                        'output_olean':str(dst.relative_to(BASE)),'olean_sha256':sha(dst),'seed_command':prior['command'],
                        'seed_log':str((seed_path.parent/prior['log']).relative_to(BASE)),'seed_log_sha256':prior['log_sha256']}
                result['checks'].append(record);save()
            else:
                execute('build_'+module,source_path(module),module=module)
        result['axiom_audit_count']=sum(c.get('axiom_audits',0) for c in result['checks'])
        assert result['axiom_audit_count']==54
        for path in (LEAN/'s11/S11D5Bulk').glob('*.lean'):
            assert not re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)',path.read_text(),re.M),path
        for name,statement,good,bad in EXAMPLES:
            for suffix,value,expected in [('positive',good,'PASS'),('mutant',bad,'REJECTED')]:
                path=SCRATCH/(name+'_'+suffix+'.lean');path.write_text(PREAMBLE+statement.replace('REL',value))
                execute(name+'_'+suffix,path,expected)
        for name,statement in EXTRA_POSITIVES:
            path=SCRATCH/(name+'.lean');path.write_text(PREAMBLE+statement+'\n');execute(name,path)
        assert native==validate_native() and preserved==validate_preserved()
        assert before==source_hashes(order,seed_path),'source changed during execution'
        fresh_env,_,fresh_meta=portable.external_environment()
        assert external==external_inputs(fresh_env,fresh_meta,order),'external dependency changed'
        validate_seed(seed_path,external)
        result['objects']={str(obj_path(m).relative_to(BASE)):sha(obj_path(m)) for m in order}
        for c in result['checks']:
            if 'olean_sha256' in c:assert sha(BASE/c['output_olean'])==c['olean_sha256']
        result.update(status='PASS',source_sha256_after=before)
    except BaseException as error:
        result.update(status='ERROR',error=repr(error),instrument_failure_is_mutation=False);raise
    finally:
        result['finished_utc']=datetime.now(timezone.utc).isoformat();save()
if __name__=='__main__':main()
