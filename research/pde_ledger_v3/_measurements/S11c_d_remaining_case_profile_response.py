#!/usr/bin/env python3
"""Complete the three missing FORM responses from accepted case bindings."""
import argparse,ast,contextlib,copy,gc,json,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_remaining_case_profile_bindings as binding
import S11c_d_remaining_case_response as response
import S11c_d_remaining_case_response_finish as aggregate
f,m,b=response.f,response.m,response.b
matrices=response.matrices
profile=binding.profile
BASELINE=binding.BASELINE
PLAN=f.M/'S11c_d_remaining_case_profile_response_plan.md'
BCP=f.M/'S11c_d_remaining_case_profile_bindings_checkpoint.json'
RCP=f.M/'S11c_d_remaining_case_response_checkpoint.json'
ECP=f.M/'S11c_d_remaining_case_boundary_checkpoint.json'
FCP=f.M/'S11c_d_remaining_case_profile_response_focused.json'


def load(base,resume=None):
    br,bc,_=binding.accepted(BCP,'ACCEPTED_CASE_PROFILE_BINDINGS')
    rr,rc,_=binding.accepted(RCP,'PUBLISHED_ANNEX_VERIFIED')
    er,ec,_=binding.accepted(ECP,'ACCEPTED_FOUR_CASE_BOUNDARY_MAPS')
    pins={}
    for checks in (bc,rc,ec):
        for name,value in checks['sourceFiles'].items():
            f.require(name not in pins or pins[name]==value,'same actual consumed source versions');pins[name]=value
    for path in (Path(__file__),PLAN,BCP,RCP,ECP,Path(aggregate.__file__),Path(response.__file__),Path(profile.__file__)):
        pins[str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':{},'copiedInputs':{},'input':bc['input'],'baselineInput':bc['baselineInput'],'settings':bc['settings'],
        'scope':'One actual FORM comparison per case, with full finite/independent-grade response and unchanged own-case end coordinates. Positive regulator, approximate boundaries and omitted parent pure-second-order terms remain.'}
    for root,checks in ((br,bc),(rr,rc),(er,ec)):
        manifest['inputPackets'][str(root/'checks.json')]=f.digest(root/'checks.json')
        for name,value in checks['inputPackets'].items():
            f.require(name not in manifest['inputPackets'] or manifest['inputPackets'][name]==value,'same original shared operand');manifest['inputPackets'][name]=value
    labels=tuple(bc['cases'])
    if resume:
        accepted=json.loads(FCP.read_text());prior=json.loads((resume/'checks.json').read_text())
        f.require(accepted['status']=='ACCEPTED_CASE_FORM_RESPONSE_INPUTS' and f.digest(resume/'checks.json')==accepted['checksSha256'],'accepted focused input/matrix-tail evidence')
        f.require(prior['sourceFiles']==pins and prior['mode']=='focused','exact focused sources')
        for name,value in prior['artifacts'].items():m.retain(resume/name,base/name,manifest,value['sha256'])
        manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json')
        manifest['inputPackets'][str(FCP)]=f.digest(FCP)
        manifest['completedFocusedReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'artifacts':len(prior['artifacts'])}
    else:
        for name,value in bc['artifacts'].items():m.retain(br/name,base/'profile-inputs'/name,manifest,value['sha256'])
        m.retain(br/'checks.json',base/'profile-inputs/checks.json',manifest)
        for label in labels:
            for kind in ('reduced-action','actions','assembly','factorization'):
                name='accepted-cases/'+label+'/'+kind+'.pickle';m.retain(br/name,base/'interiors'/name,manifest,bc['artifacts'][name]['sha256'])
            if label==BASELINE:
                original=f.unpickle(base/'profile-inputs/accepted-bindings'/label/'case-binding.pickle')
                case={'binding':f.unpickle(base/'profile-inputs/baseline-binding-view.pickle'),'grades':original['grades'],'reusedRows':[{'row':i,'fromCase':BASELINE,'fromRow':i} for i in range(len(original['binding']['bound']['rows']))],'newRows':[]}
                dest=base/'interiors/accepted-bindings'/label/'case-binding.pickle';dest.parent.mkdir(parents=True);f.atomic_pickle(dest,case)
            else:
                name='cases/'+label+'/case-binding.pickle';m.retain(br/name,base/'interiors/accepted-bindings'/label/'case-binding.pickle',manifest,bc['artifacts'][name]['sha256'])
            for name,value in rc['artifacts'].items():
                if name.startswith('cases/'+label+'/') and name.endswith('.pickle'):m.retain(rr/name,base/'unchanged-response'/name,manifest,value['sha256'])
            name='boundary-cases/'+label+'/case-boundary.pickle'
            m.retain(rr/name,base/'unchanged-response'/name,manifest,rc['artifacts'][name]['sha256'])
        for name,value in ec['artifacts'].items():
            if name.startswith('boundary-cases/') or name=='accepted-modes/reference/modal.pickle':m.retain(er/name,base/name,manifest,value['sha256'])
        system=br/'accepted-profile/ablation/finite-system.pickle'
        m.retain(system,base/'interiors/accepted-finite-system.pickle',manifest)
        for name in ('finite-system.pickle','finite-solution.pickle'):
            m.retain(br/'accepted-profile/ablation'/name,base/'cases'/BASELINE/'finite'/name,manifest)
        m.retain(br/'accepted-profile/full.out',base/'cases'/BASELINE/'continuum/full.out',manifest)
        pc=json.loads(binding.PCP.read_text());m.retain(Path(pc['runDirectory'])/'checks.json',base/'cases'/BASELINE/'continuum/checks.json',manifest,pc['checksSha256'])
    for name,value in pins.items():
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest);f.require(f.digest(dest)==value,'frozen FORM response source')
    manifest['settings']=binding.restore_settings(manifest['settings'],f.unpickle(base/'interiors/accepted-finite-system.pickle')['settings'])
    f.save(base/'inputs.json',manifest)
    return manifest,labels


def matrix_tail():
    main=response.function(ast.parse(Path(matrices.__file__).read_text()),'main')
    loop=next(n for n in main.body if isinstance(n,ast.For) and ast.unparse(n.target)=='(label, case)')
    first=next(i for i,n in enumerate(loop.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='(r, adapter, packets)')
    body=copy.deepcopy(loop.body[first:])
    f.require(ast.dump(ast.Module(body=body,type_ignores=[]))==ast.dump(ast.Module(body=loop.body[first:],type_ignores=[])),'whole unchanged native row integration and independent assembly tail')
    names=('base','label','case','system','all_rows','summaries','new_nodes','baseline','target')
    node=ast.FunctionDef(name='complete_case_matrices',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body+[ast.Return(ast.Tuple([ast.Subscript(ast.Name('summaries',ast.Load()),ast.Name('label',ast.Load()),ast.Load()),ast.Name('new_nodes',ast.Load())],ast.Load()))],decorator_list=[])
    return response.compile_function(node,vars(matrices)),{'wholeMatrixTailAstJoin':True,'nativeAccumulatorUnchanged':True,'nativeAssemblerUnchanged':True}


def profile_emitter(label):
    tree=ast.parse(Path(profile.__file__).read_text());main=response.function(tree,'main');emitter=copy.deepcopy(response.function(tree,'emit_result'))
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='summary')
    body=copy.deepcopy(main.body[first:stop]);original=copy.deepcopy(body);prefix='s11cd'+label.replace('__','_')+'ProfileForm'
    class Rename(ast.NodeTransformer):
        def __init__(self,old,new):self.old,self.new,self.count=old,new,0
        def visit_Constant(self,n):
            if n.value==self.old:self.count+=1;return ast.Constant(self.new)
            return n
    rename=Rename('s11cdProfileForm',prefix);body=[rename.visit(n) for n in body];f.require(rename.count==1,'one case-owned FORM export namespace')
    undo=Rename(prefix,'s11cdProfileForm');restored=[undo.visit(copy.deepcopy(n)) for n in body]
    f.require(ast.dump(ast.Module(body=restored,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[])),'whole original FORM emission/replay tail')
    names=('base','result','r','pins','operands','before')
    node=ast.FunctionDef(name='emit_case_form',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body+[ast.Return(ast.Dict(keys=[ast.Constant(k) for k in ('tags','keys','metadataPaths')],values=[ast.Call(ast.Name('len',ast.Load()),[ast.Name('entries',ast.Load())],[]),ast.Name('keys',ast.Load()),ast.Name('paths',ast.Load())]))],decorator_list=[])
    scope=dict(vars(profile),PREFIX='PROFILE_FORM_'+label.replace('__','_'));scope['emit_result']=response.compile_function(emitter,scope)
    return response.compile_function(node,scope)


def preflight(base,manifest,labels):
    cases={label:f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle') for label in labels}
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');old_profile=f.unpickle(base/'profile-inputs/accepted-profile/profile-form.pickle')
    prepared=matrices.preflight(base/'interiors',cases,basis)
    joins={}
    for label,case in cases.items():
        old=f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-system.pickle');ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle')
        old_ends=f.unpickle(base/'unchanged-response/boundary-cases'/label/'case-boundary.pickle')
        f.require(m.same(ends,old_ends) and m.same(old['channels'],ends['finite']),'full own-case unchanged finite and continuum end/current/phase inputs')
        f.require(np.array_equal(old['nodes'],basis['nodes']) and m.same(old['derivativeMatrices'],basis['derivativeMatrices']) and old['settings']==basis['settings'],'complete actual finite trial and quadrature settings')
        f.require(m.same(ends['fieldUnits'],case['binding']['fieldUnits']) and m.same(ends['rowUnits'],case['binding']['equationUnits']),'actual field/equation maps before response')
        if label!=BASELINE:
            endpoints=case['endpoints'];f.require(all(a==c for a,c in endpoints.values()),'all actual constant-end profile limits')
        else:f.require(m.same(basis['channels'],ends['finite']),'accepted baseline FORM end maps')
        r,adapter,_=matrices.context(base/'interiors',label,case,manifest['input'])
        f.require(binding.native.same(adapter.input.profiles,old_profile['alteredShapes']),'actual profile before moment reuse')
        f.require(profile.shape_input(manifest['baselineInput'])==manifest['input'],'only authorized profile FORM input change')
        wrong=copy.deepcopy(ends['finite']);wrong['LEFT']['incomingBoundaryData']*=-1
        f.require(not m.same(wrong,old['channels']),'actual changed incoming boundary input rejected')
        joins[label]={'sameCompleteEndPackets':True,'sameFiniteChannels':True,'sameFiniteBasis':True,'sameFieldUnits':True,'sameMomentProfiles':True,'incomingMutationRejected':True}
        if label!=BASELINE:profile_emitter(label)
    _,matrix_join=matrix_tail();_,_,finite_join=response.finite_tail()
    result={'rowReuse':prepared,'endJoins':joins,'matrixTail':matrix_join,'finiteTail':finite_join,'emissionTailJoins':3,'baselineFormSolvesRepeated':0,'newQuadratureNodes':0,'newModes':0,'newSolves':0}
    f.save(base/'preflight.json',result);return result


def compare_case(base,label,manifest):
    directory=base/'cases'/label;oldbase=base/'unchanged-response/cases'/label
    finite=f.unpickle(directory/'finite/finite-solution.pickle');new=f.unpickle(directory/'finite/observable.pickle');old=f.unpickle(oldbase/'finite/observable.pickle')
    current=f.unpickle(directory/'continuum/continuum-response.pickle');original=f.unpickle(oldbase/'continuum/continuum-response.pickle')
    fresh=current['response'];prior=original['response']
    f.require(m.same(new['positions'],old['positions']) and m.same(new['incomingPhase'],old['incomingPhase']) and m.same(new['outgoingPhase'],old['outgoingPhase']),'same actual finite common-origin coordinates')
    flux_keys=('amplitudeHomotopy','outgoingCurrentHomotopy','incomingCurrentHomotopy','openOutgoingFluxHomotopy','openOutgoingFractionCoefficients')
    f.require(all(current['flux'][key].keys()==original['flux'][key].keys() for key in flux_keys),'complete actual current homotopy support')
    for key in ('incomingOriginPhase','outgoingOriginPhase','incomingCurrentOrigin','outgoingCurrentOrigin','incomingCurrentRootInverse','fieldIncomingMap','labels'):
        f.require(m.same(fresh[key],prior[key]),'actual identical continuum incident/current/channel coordinates')
    f.require(m.same(new['originCurrent'],old['originCurrent']),'actual identical finite current coordinates')
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle')
    grid=f.polynomial_basis(new['positions'],basis['settings']['sourceBound'],len(basis['nodes']))
    def fields(value):
        maps=b.J.multiply(value['response']['incomingOriginPhase'],value['response']['incomingCurrentRootInverse']);size=len(basis['nodes'])
        parts=[b.J.multiply({g:grid@array[i*size:(i+1)*size] for g,array in value['solve']['coefficients'].items()},maps) for i in range(5)]
        return {g:np.stack([part[g] for part in parts]) for g in b.G}
    new_fields,old_fields=fields(current),fields(original)
    comparisons={'finiteScattering':new['originScattering']-old['originScattering'],'finiteCurrent':new['totalCurrentRatio']-old['totalCurrentRatio'],'finiteFields':new['originFields']-old['originFields'],
        'continuum':{key:b.subtract(fresh[key],prior[key]) for key in ('fieldOriginScattering','fluxOriginScattering')},
        'closedMatching':{end:b.subtract(fresh['outgoingFieldAmplitudes'][end],prior['outgoingFieldAmplitudes'][end]) for end in ('LEFT','RIGHT')},
        'continuumCurrent':{key:{g:current['flux'][key][g]-original['flux'][key][g] for g in current['flux'][key]} for key in flux_keys},
        'continuumOriginFields':{'altered':new_fields,'baseline':old_fields,'difference':b.subtract(new_fields,old_fields)},'fieldPositions':new['positions']}
    f.atomic_pickle(directory/'form-comparisons.pickle',comparisons)
    coefficients=f.unpickle(base/'interiors/cases'/label/'interior-matrices.pickle');old_profile=f.unpickle(base/'profile-inputs/accepted-profile/profile-form.pickle')
    computed={'finite':finite,'finiteOriginScattering':new['originScattering'],'baselineFiniteOriginScattering':old['originScattering'],'solve':current['solve'],'response':fresh,'baselineResponse':prior,'comparisons':comparisons,
        'settings':current['settings'],'independentFiniteDifference':new['independentDifference'],'recombination':coefficients['comparisons']['approved_total']}
    result={'computed':computed,'moments':old_profile['moments'],'baselineShapes':old_profile['baselineShapes'],'alteredShapes':old_profile['alteredShapes'],'fieldUnits':current['fieldUnits'],'currentUnit':current['currentUnit'],'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']+' Closed amplitudes stay boundary anchored. Depth-integrated normal current is not bulk-depth escape. Separate bump moment is not bump scattering.'}
    f.atomic_pickle(directory/'continuum/profile-form.pickle',result)
    summary={'finiteAmplitudeChange':b.norm(comparisons['finiteScattering']),'finiteCurrentChange':b.norm(comparisons['finiteCurrent']),'commonGridFieldChange':b.norm(comparisons['finiteFields']),
        'baselineFieldMaximum':b.norm(old['originFields']),'baselineAmplitudeMaximum':b.norm(old['originScattering']),
        'continuumCoefficientChanges':{k:{str(g):b.norm(v) for g,v in series.items()} for k,series in comparisons['continuum'].items()},
        'continuumCurrentChanges':{key:{str(g):b.norm(v) for g,v in series.items()} for key,series in comparisons['continuumCurrent'].items()},
        'continuumCommonGridFieldChanges':{str(g):b.norm(v) for g,v in comparisons['continuumOriginFields']['difference'].items()},
        'sameActualEndCoordinates':True,'closedMatchingChanges':{end:{str(g):b.norm(v) for g,v in series.items()} for end,series in comparisons['closedMatching'].items()}}
    f.save(directory/'form-checks.json',summary);return summary


def emit_case(base,label,manifest):
    directory=base/'cases'/label/'continuum';result=f.unpickle(directory/'profile-form.pickle');case=f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle')
    r,adapter,_=matrices.context(base/'interiors',label,case,manifest['input']);current=f.unpickle(directory/'continuum-response.pickle')
    f.engine.PHYSICAL_METADATA.dimensions.__dict__.update(current['dimensionState'])
    emitted=profile_emitter(label)(directory,result,r,manifest['sourceFiles'],manifest['inputPackets'],f.digest(directory/'profile-form.pickle'))
    f.save(directory/'emission-checks.json',emitted);return emitted


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--focused',action='store_true');ap.add_argument('--resume-focused',type=Path);args=ap.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    f.require(args.focused!=(args.resume_focused is not None),'focused preparation or completed focused reuse')
    manifest,labels=load(base,args.resume_focused)
    focused=json.loads((base/'preflight.json').read_text()) if args.resume_focused else preflight(base,manifest,labels)
    matrix_counts={};solves={};comparisons={};emissions={};new_nodes=0;combined=None
    if not args.focused:
        basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');all_rows={BASELINE:f.unpickle(base/'profile-inputs/baseline-profile-row-matrices.pickle')};complete,_=matrix_tail()
        for label in labels:
            if label==BASELINE:continue
            target=base/'interiors/cases'/label;target.mkdir(parents=True,exist_ok=True);case=f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle')
            _,new_nodes=complete(base/'interiors',label,case,basis,all_rows,matrix_counts,new_nodes,None,target);gc.collect()
        # All native matrix packets precede any solve; all complete responses
        # precede comparison and emission. Output failures cannot require reruns.
        for label in labels:
            if label==BASELINE:continue
            solves[label]=response.construct_case(base,label,manifest);f.save(base/'response-inventory.json',solves);gc.collect()
        for label in labels:
            if label!=BASELINE:comparisons[label]=compare_case(base,label,manifest);f.save(base/'comparison-inventory.json',comparisons);gc.collect()
        for label in labels:
            if label!=BASELINE:emissions[label]=emit_case(base,label,manifest);gc.collect()
        with (base/'original-combined.out').open('xb') as stream:
            for label in labels:stream.write((base/'cases'/label/'continuum/full.out').read_bytes())
        combined=aggregate.aggregate(base,labels);f.save(base/'aggregation-checks.json',combined)
        f.atomic_pickle(base/'remaining-case-profile-form.pickle',{'cases':solves,'comparisons':comparisons,'baseline':str(base/'profile-inputs/accepted-profile/profile-form.pickle'),'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen source pre/post')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original input pre/post')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'unchanged completed operand copy')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,'mode':'focused' if args.focused else 'construct','status':'PASSED_CASE_FORM_RESPONSE_INPUTS' if args.focused else 'COMPLETED_FOUR_CASE_PROFILE_FORM','preflight':focused,'matrixCases':matrix_counts,'responseCases':solves,'comparisons':comparisons,'emissions':emissions,'aggregate':combined,'newNodes':new_nodes,'newFiniteResponses':len(solves),'newContinuumResponses':len(solves),'newModeConstructions':0,'newCurrentClosures':0,'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
