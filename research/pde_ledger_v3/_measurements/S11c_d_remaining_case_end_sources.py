#!/usr/bin/env python3
"""Fresh case end symbols and native current-generator input joins."""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import textwrap
import time
import sympy as sp
import S11c_d_remaining_case_bindings as binding
import S11c_d_uniform_source as uniform

f,engine=binding.f,binding.engine
PLAN=f.M/'S11c_d_remaining_case_end_sources_plan.md'
MATRICES=f.M/'S11c_d_remaining_case_matrices_checkpoint.json'
UNIFORM=f.M/'S11c_d_uniform_source_checkpoint.json'
PRODUCER=f.STORE/'s11c-thickness-coordinate-20260914/d_full'
BASELINE=binding.BASELINE


def copy_input(path,target,sha,inputs):
    f.require(f.digest(path)==sha,('source input',str(path)))
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
    f.require(f.digest(target)==sha,'byte-identical source copy');inputs[str(path)]=sha
    return f.unpickle(target)


def load(base):
    cp=json.loads(MATRICES.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']=='ACCEPTED_FOUR_CASE_INTERIOR_MATRICES','accepted all-case interiors')
    f.require(f.digest(origin/'checks.json')==cp['checksSha256'],'completed matrix checkpoint')
    for n,h in cp['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(origin/'source'/n)==h,('unchanged matrix source',n))
    inputs={str(origin/'checks.json'):cp['checksSha256']};cases={}
    for label in cp['cases']:
        cases[label]={}
        for kind in ('reduced-action','actions','assembly'):
            name='accepted-cases/'+label+'/'+kind+'.pickle'
            cases[label][kind]=copy_input(origin/name,base/name,cp['artifacts'][name]['sha256'],inputs)
    uc=json.loads(UNIFORM.read_text());uo=Path(uc['runDirectory'])
    f.require(uc['status']=='PUBLISHED_ANNEX_VERIFIED','accepted baseline uniform source')
    f.require(f.digest(f.ROOT/uc['publication']['path'])==uc['publication']['sha256'],'baseline source annex join')
    baseline=copy_input(uo/'uniform-source.pickle',base/'accepted-uniform-source.pickle',uc['checks']['artifacts']['uniform-source.pickle']['sha256'],inputs)
    producer=json.loads((PRODUCER/'manifest.json').read_text())
    f.require(producer['exit_code']==0 and producer['source_hashes_before']==producer['source_hashes_after'],'completed stable four-case symbol producer')
    frozen=PRODUCER/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    f.require(f.digest(frozen)==producer['source_hashes_after']['scripts/S11c_d_mixing_scattering_sympy_audit.py'],'original native producer engine')
    joins=binding.factors.cases.definition_joins(frozen,{'ConstantEndPencil','ReducedPencil','UniformSlabCurrent','ChannelInput'})
    for n,h in producer['source_hashes_after'].items():
        if n.endswith('_exports.py') or n.startswith('directives/') or n=='scripts/ledger_fold.py':
            f.require(f.digest(f.ROOT/n)==h,('original physical producer source',n))
    inputs.update({str(PRODUCER/'manifest.json'):f.digest(PRODUCER/'manifest.json'),str(frozen):f.digest(frozen)})
    cached={}
    for label in cases:
        cached[label]={}
        for end in ('REFERENCE','LEFT','RIGHT'):
            name='symbols/'+end+'_'+label.replace('__','_')+'.pickle'
            cached[label][end]=copy_input(PRODUCER/name,base/'accepted-symbols'/Path(name).name,producer['artifacts'][name]['sha256'],inputs)
    current_sources={};extra=[]
    for end in ('LEFT','RIGHT'):
        checkpoint=f.M/('S11c_d_end_current_source_'+end.lower()+'_thickness_repair_checkpoint.json')
        c=json.loads(checkpoint.read_text());directory=Path(c['runDirectory']);path=directory/'objects.pickle'
        f.require(c['retainedNonzeroScalars']==c['nonzeroCancellationIdentities']==0,'accepted baseline current source')
        current_sources[end]=copy_input(path,base/('accepted-'+end.lower()+'-current-source.pickle'),c['artifacts']['objects.pickle']['sha256'],inputs)
        original=directory/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
        f.require(uniform.class_ast(original.read_text(),'UniformSlabCurrent')==uniform.class_ast(engine.HERE.read_text(),'UniformSlabCurrent'),'unchanged current input constructor')
        inputs[str(original)]=f.digest(original);extra.append(checkpoint)
    specification=json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text())
    f.require(specification==json.loads((origin/'inputs.json').read_text())['input'],'actual approved input')
    pins=dict(cp['sourceFiles'])
    for p in (Path(__file__).resolve(),PLAN,MATRICES,UNIFORM,Path(uniform.__file__),*extra):pins[str(p.relative_to(f.ROOT))]=f.digest(p)
    for n,h in pins.items():
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target)
        f.require(f.digest(target)==h,'frozen source identity')
    manifest={'sourceFiles':pins,'inputPackets':inputs,'nativeDefinitionJoins':joins,'input':specification,
        'scope':'Fresh constant-background end pencils and current-generator operands for all cases; no new mode/current matrix, scattering solve or physical pole result.'}
    manifest['copiedInputs']={str(p.relative_to(base)):f.digest(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    f.save(base/'inputs.json',manifest)
    return cases,cached,baseline,current_sources,manifest


def prefix_function():
    source=textwrap.dedent(inspect.getsource(engine.UniformSlabCurrent.construct))
    original=ast.parse(source).body[0]
    stop=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='virtual' for t in n.targets))
    node=copy.deepcopy(original);node.decorator_list=[];node.name='current_generator_inputs';node.body=node.body[:stop]
    names={'SOURCE_ENERGY':'source','UNIFORM_SOURCE_ENERGY':'uniform_source','SOURCE_PARAMETER_ALIGNMENT':'tuple(parameter_map.items())',
        'FIELD_MAP':'field_map','HARMONIC_ENERGY':'harmonic_energy','TANGENTIAL_ENERGY_REDUCTION':'energy',
        'ZERO_TRANSFER_MASS_ROW':'mass_control','MATERIAL_CONSTRAINT_COEFFICIENTS':'constraint_coefficients',
        'MATERIAL_CONSTRAINT_TRUNCATION_REMAINDER':'constraint_residual'}
    node.body.append(ast.Return(ast.Dict(keys=[ast.Constant(k) for k in names],values=[ast.parse(v,mode='eval').body for v in names.values()])))
    f.require(ast.dump(ast.Module(body=node.body[:-1],type_ignores=[]))==ast.dump(ast.Module(body=original.body[:stop],type_ignores=[])), 'whole native current prefix AST identity')
    module=ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[]));namespace=dict(vars(engine));exec(compile(module,'<native-current-generator-prefix>','exec'),namespace)
    return namespace[node.name],{'nativeFunctionSha256':hashlib.sha256(source.encode()).hexdigest(),
        'nativePrefixAstSha256':hashlib.sha256(ast.dump(ast.Module(body=original.body[:stop],type_ignores=[])).encode()).hexdigest(),
        'nativePrefixStatements':stop,'stopBefore':'virtual variation; no current matrix is reconstructed by this prefix'}


def expanded_difference(a,b):
    if isinstance(a,sp.MatrixBase):return (a-b).applyfunc(sp.expand)
    if isinstance(a,(tuple,list)):return tuple(expanded_difference(x,y) for x,y in zip(a,b))
    return sp.expand(a-b)


def scalars(value):
    if isinstance(value,sp.MatrixBase):yield from value
    elif isinstance(value,dict):
        for v in value.values():yield from scalars(v)
    elif isinstance(value,(tuple,list)):
        for v in value:yield from scalars(v)
    else:yield value


def current_inputs(target,label,end,r,inp,strong,energy,units,function):
    state=copy.copy(r);state.end_values={key:inp.limits[value] for key,value in r.end_values.items()}
    ends=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);ends.r=state;ends.kn=sp.Symbol('s11cdSpectralNormalMomentum',real=True)
    current=engine.UniformSlabCurrent(state,{'value':energy},ends,strong[3,:])
    endpoint={'REFERENCE':None,'LEFT':-sp.oo,'RIGHT':sp.oo}[end]
    anchoring=label.split('__')[0]
    data=function(current,anchoring,endpoint)
    data.update(strong=strong,fieldUnits=tuple(current.field_units),end=end,case=label,
        profileBindings=inp.limits,sourceEnergy=energy,anchoring=anchoring,
        retainedConstraintResidual=data['MATERIAL_CONSTRAINT_TRUNCATION_REMAINDER'].applyfunc(current.retained))
    f.atomic_pickle(target/(end.lower()+'-current-inputs.pickle'),data)
    f.require(all(v==0 for v in data['retainedConstraintResidual']),'actual retained homogeneous material constraint')
    f.require(tuple(current.field_units)==tuple(units),'current/source field units')
    f.require(not data['TANGENTIAL_ENERGY_REDUCTION'].has(sp.Integral),'complete current energy reduction')
    return data


def family_comparison(target,name,current,previous):
    pairs={k:(current[k],previous[k]) for k in ('strong','weak','curl','energy','constraint','mass')}
    f.atomic_pickle(target/(name+'-pairs.pickle'),pairs)
    differences={k:expanded_difference(a,b) for k,(a,b) in pairs.items()}
    result={'differences':differences,'nonzeroScalars':{k:sum(v!=0 for v in scalars(x)) for k,x in differences.items()},
        'fieldUnitsIdentical':binding.same(current['fieldUnits'],previous['fieldUnits']),
        'profileBindingsIdentical':binding.same(current['profileBindings'],previous['profileBindings'])}
    result['sameGeneratorInputs']=not any(result['nonzeroScalars'].values()) and result['fieldUnitsIdentical'] and result['profileBindingsIdentical']
    f.atomic_pickle(target/(name+'-comparison.pickle'),result)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));start=time.monotonic()
    def timeout(*_):raise TimeoutError('remaining end source budget; preserve completed symbols and current inputs')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    def progress(stage):
        with (base/'progress.jsonl').open('a') as s:s.write(json.dumps({'stage':stage,'wallSeconds':time.monotonic()-start})+'\n')
    cases,cached,baseline,accepted_currents,manifest=load(base);function,prefix=prefix_function();f.save(base/'prefix-ast.json',prefix)
    outputs={};current_packets={};families=[];inventory={};focused={}
    for label,packets in cases.items():
        target=base/'cases'/label;target.mkdir(parents=True,exist_ok=True);progress(label+'_started')
        r,dimensions,pencil=binding.factors.context(packets['reduced-action'],packets['actions'],packets['assembly'])
        inp=engine.ChannelInput(r,manifest['input'])
        if label==BASELINE:
            result=baseline
            shutil.copyfile(base/'accepted-uniform-source.pickle',target/'uniform-source.pickle')
            f.require(f.digest(target/'uniform-source.pickle')==f.digest(base/'accepted-uniform-source.pickle'),'unchanged accepted baseline source')
        else:
            result=uniform.construct(target,r,pencil,inp,cached[label],lambda name:progress(label+'_'+name))
            f.atomic_pickle(target/'uniform-source.pickle',result)
        outputs[label]=result;current_packets[label]={}
        field_units=tuple(dimensions.known[v] for v in pencil.fields)
        for end in ('REFERENCE','LEFT','RIGHT'):
            record=result['records'][end];strong=record['strong'];energy=cached[label][end][5]
            f.require(binding.same(energy,cached[BASELINE][end][5]),'actual original energy source identity')
            data=current_inputs(target,label,end,r,inp,strong,energy,field_units,function);current_packets[label][end]=data
            if label==BASELINE and end in accepted_currents:
                old=accepted_currents[end]['conservative'];names=('TANGENTIAL_ENERGY_REDUCTION','ZERO_TRANSFER_MASS_ROW','MATERIAL_CONSTRAINT_COEFFICIENTS')
                pairs={n:(data[n],old[n]) for n in names};f.atomic_pickle(target/(end.lower()+'-accepted-current-pairs.pickle'),pairs)
                proof={n:expanded_difference(a,b) for n,(a,b) in pairs.items()};f.atomic_pickle(target/(end.lower()+'-accepted-current-proofs.pickle'),proof)
                f.require(all(v==0 for v in scalars(proof)),'native prefix reproduces accepted actual current-generator operands')
                # Same-unit coefficient mutation of the actual reduced energy.
                response=data['TANGENTIAL_ENERGY_REDUCTION'];f.require(response!=0,'actual energy coefficient mutation responds')
                focused[end]={'zeroScalars':sum(1 for _ in scalars(proof)),'energyMutationNonzero':True}
                f.save(base/'focused-current-inputs.json',focused)
            signature={'strong':strong,'weak':record['weak'],'curl':result['curl'],
                'energy':data['TANGENTIAL_ENERGY_REDUCTION'],'constraint':data['MATERIAL_CONSTRAINT_COEFFICIENTS'],
                'mass':data['ZERO_TRANSFER_MASS_ROW'],'fieldUnits':field_units,'profileBindings':inp.limits}
            matches=[];comparisons=[]
            for i,(address,previous) in enumerate(families):
                comparison=family_comparison(target,end.lower()+'-family-'+str(i),signature,previous)
                comparisons.append({'family':i,'origin':address,'nonzeroScalars':comparison['nonzeroScalars'],'sameGeneratorInputs':comparison['sameGeneratorInputs']})
                if comparison['sameGeneratorInputs']:matches.append(i);break
            if matches:family=matches[0]
            else:family=len(families);families.append(((label,end),signature))
            inventory[label+'__'+end]={'family':family,'origin':families[family][0],'comparisons':comparisons,
                'strongShape':strong.shape,'weakShape':record['weak'].shape,
                'couplingNonzero':{k:sum(v!=0 for v in x) for k,x in record['coupling'].items()},
                'scope':'Exact symbolic generator-input family; no mode/current/boundary reuse has yet been performed.'}
            f.save(base/'end-inventory.json',inventory);progress(label+'_'+end+'_current_inputs_saved')
        f.atomic_pickle(target/'case-end-sources.pickle',{'source':result,'currentInputs':current_packets[label],
            'case':label,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],
            'dimensionState':dict(vars(dimensions))})
    f.require(set(focused)=={'LEFT','RIGHT'},'both accepted current-generator input checks')
    f.atomic_pickle(base/'remaining-case-end-sources.pickle',{'sources':outputs,'currentInputs':current_packets,
        'families':families,'inventory':inventory,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets']})
    for n,h in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==h,'unchanged current/frozen source')
    for n,h in manifest['inputPackets'].items():f.require(f.digest(Path(n))==h,'unchanged original input')
    for n,h in manifest['copiedInputs'].items():f.require(f.digest(base/n)==h,'unchanged copied inputs')
    artifacts={str(p.relative_to(base)):binding.factors.artifact(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    checks={**manifest,'runDirectory':str(base),'status':'COMPLETED_CASE_END_AND_CURRENT_INPUT_SOURCES','endInventory':inventory,
        'distinctGeneratorFamilies':len(families),'focusedCurrentInputs':focused,'prefixAst':prefix,'artifacts':artifacts,
        'newQuadratureNodes':0,'newModeSolves':0,'newCurrentMatrices':0,'newScatteringSolves':0,
        'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
