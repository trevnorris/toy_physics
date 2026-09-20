#!/usr/bin/env python3
"""Complete case responses from accepted interior matrices and actual end maps."""
import argparse,ast,copy,contextlib,json,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import scipy.linalg as la
import sympy as sp
import S11c_d_remaining_case_boundary as maps
import S11c_d_remaining_case_matrices as matrices
import S11c_d_continuum_response as response
import S11c_d_finite_scattering_domain as domain
f,m,b=maps.f,maps.modes,maps.boundary
PLAN=f.M/'S11c_d_remaining_case_response_plan.md'
BCP=f.M/'S11c_d_remaining_case_boundary_checkpoint.json'
ICP=f.M/'S11c_d_remaining_case_matrices_checkpoint.json'
RCP=f.M/'S11c_d_continuum_response_checkpoint.json'
DCP=f.M/'S11c_d_finite_scattering_domain_checkpoint.json'
BASELINE=m.BASELINE


def unsplit_arrays(packet):
    # The native direct binding has already evaluated the physical contrasts.
    # Its sole collection slot is storage for that full unsplit matrix, not
    # the independent zero-contrast coefficient of the operator.
    keys=('local','nonlocal','total')
    f.require(all(set(packet[k])=={(0,0,0)} for k in keys),'one actual fully bound unsplit slot')
    result={k:packet[k][(0,0,0)] for k in keys}
    f.require(all(isinstance(v,np.ndarray) and v.shape==(645,645) and np.isfinite(v).all() for v in result.values()),'complete actual unsplit matrices')
    f.require(np.array_equal(result['local']+result['nonlocal'],result['total']),'literal full unsplit decomposition')
    return {**packet,**result}


def function(tree,name):return next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)


def compile_function(node,namespace):
    scope=dict(namespace);exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),__file__,'exec'),scope);return scope[node.name]


def finite_tail():
    original=function(ast.parse(Path(f.__file__).read_text()),'construct')
    start=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='original' for t in n.targets))
    body=copy.deepcopy(original.body[start:]);names=('base','matrix','channels','size','derivative','x','local_matrix','setting','pins','operands','differentiation','groups')
    node=ast.FunctionDef(name='finish_finite',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    f.require(ast.dump(ast.Module(body=node.body,type_ignores=[]))==ast.dump(ast.Module(body=original.body[start:],type_ignores=[])),'literal complete finite solve/output tail')
    # The prefix up to the first packet write is available for baseline replay without solving.
    stop=next(i for i,n in enumerate(body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='atomic_pickle')
    prepare=copy.deepcopy(node);prepare.name='finite_system';prepare.body=prepare.body[:stop]+[ast.Return(ast.Tuple([ast.Name(n,ast.Load()) for n in ('matrix','rhs','original')],ast.Load()))]
    return compile_function(node,vars(f)),compile_function(prepare,vars(f)),{'wholeFiniteSolveTailJoin':True,'firstTailStatement':start,'boundaryPrefixStatements':stop,'changedNumericalStatements':0}


def emitter(label):
    tree=ast.parse(Path(response.__file__).read_text());emit=copy.deepcopy(function(tree,'emit_result'));main=function(tree,'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='summary' for t in n.targets))
    body=copy.deepcopy(main.body[first:stop]);original=copy.deepcopy(body);changes=[]
    key_prefix='s11cd'+label.replace('__','_')+'ContinuumResponse'
    class Keys(ast.NodeTransformer):
        def visit_Constant(self,n):
            if n.value=='s11cdContinuumResponse':changes.append(True);return ast.Constant(key_prefix)
            return n
    body=[Keys().visit(n) for n in body];f.require(len(changes)==1,'one case-owned export key namespace')
    class Undo(ast.NodeTransformer):
        def visit_Constant(self,n):return ast.Constant('s11cdContinuumResponse') if n.value==key_prefix else n
    restored=[Undo().visit(copy.deepcopy(n)) for n in body]
    f.require(ast.dump(ast.Module(body=restored,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[])),'whole original emission and validation tail')
    node=ast.FunctionDef(name='emit_case',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','r','pins','operands','before')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body+[ast.Return(ast.Dict(keys=[ast.Constant(k) for k in ('tags','keys','metadataPaths')],values=[ast.Call(ast.Name('len',ast.Load()),[ast.Name('entries',ast.Load())],[]),ast.Name('keys',ast.Load()),ast.Name('metadata_paths',ast.Load())]))],decorator_list=[])
    namespace=dict(vars(response),PREFIX='CONTINUUM_RESPONSE_'+label.replace('__','_'));namespace['emit_result']=compile_function(emit,namespace)
    return compile_function(node,namespace)


def checked_checkpoint(path,status):
    cp=json.loads(path.read_text());root=Path(cp['runDirectory'])
    f.require(cp['status']==status and f.digest(root/'checks.json')==cp['checksSha256'],'accepted final checkpoint')
    for n,v in cp['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(root/'source'/n)==v,'accepted source/frozen joins')
    for n,v in cp['artifacts'].items():f.require(f.digest(root/n)==v['sha256'],'accepted artifact')
    return cp,root


def load(base,resume=None):
    bc,br=checked_checkpoint(BCP,'ACCEPTED_FOUR_CASE_BOUNDARY_MAPS');ic,ir=checked_checkpoint(ICP,'ACCEPTED_FOUR_CASE_INTERIOR_MATRICES')
    pins={**ic['sourceFiles'],**bc['sourceFiles']}
    f.require(all(pins[n]==v for n,v in ic['sourceFiles'].items()),'shared source versions')
    for path in (Path(__file__),PLAN,BCP,ICP,RCP,DCP,Path(response.__file__),Path(domain.__file__)):
        pins[str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':{str(br/'checks.json'):bc['checksSha256'],str(ir/'checks.json'):ic['checksSha256']},'copiedInputs':{},'settings':json.loads((ir/'inputs.json').read_text())['settings'],'input':json.loads((br/'inputs.json').read_text())['input'],'scope':'Four actual case finite and independent-grade responses; finite regulator and approximate modal boundaries.'}
    if resume:
        old=json.loads((resume/'checks.json').read_text());f.require(old['mode']=='focused' and old['status']=='PASSED_BASELINE_RESPONSE_WIRING','accepted focused baseline replay')
        f.require(old['sourceFiles']==pins,'same focused constructor and sources')
        for n,v in old['artifacts'].items():m.retain(resume/n,base/n,manifest,v['sha256'])
        manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json');manifest['completedFocusedReuse']=str(resume)
    else:
        for n,v in ic['artifacts'].items():m.retain(ir/n,base/'interiors'/n,manifest,v['sha256'])
        for n,v in bc['artifacts'].items():
            if n.startswith('boundary-cases/') or n=='accepted-modes/reference/modal.pickle':m.retain(br/n,base/n,manifest,v['sha256'])
        # The accepted baseline response is copied; no baseline solve or emission repeats.
        dc=json.loads(DCP.read_text());f.require(dc['status']=='VALIDATED_FINITE_BOUNDARY_REGULATOR_RESPONSE','accepted baseline finite response')
        droot=Path(dc['runDirectory'])/'regulator/complete';record=next(v for v in dc['savedOperandVerification'] if v['name']=='regulator')
        for name in ('finite-system.pickle','finite-solution.pickle'):m.retain(droot/name,base/'cases'/BASELINE/'finite'/name,manifest,record['packetHashes'][name]['sha256'])
        rc,rr=checked_checkpoint(RCP,'PUBLISHED_ANNEX_VERIFIED')
        for n,v in rc['artifacts'].items():m.retain(rr/n,base/'cases'/BASELINE/'continuum'/n,manifest,v['sha256'])
    for n,v in pins.items():
        dst=base/'source'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dst);f.require(f.digest(dst)==v,'frozen response source')
    f.save(base/'inputs.json',manifest)
    return manifest,tuple(ic['cases'])


def finite_view(system,solution,independent):
    a=system['matrix'];x=solution['coefficients'];rhs=system['rhs'];size=len(system['nodes'])
    residual=a@x-rhs;f.require(np.array_equal(residual,solution['equationResidual']),'literal finite equation residual')
    f.require(solution['rank']==len(a) and b.norm(residual/solution['rowScale'][:,None])<1e-9,'full rank and scaled equations')
    other=np.linalg.solve(a,rhs) if independent else None
    difference=other-x if independent else None
    if independent:f.require(b.norm(difference)/(1+b.norm(x))<1e-8,'independent finite solve')
    boundary_residuals={};trace_residuals={};offset=0
    values=np.stack([system['derivativeMatrices'][0]@x[i*size:(i+1)*size] for i in range(5)])
    derivatives=np.stack([system['derivativeMatrices'][1]@x[i*size:(i+1)*size] for i in range(5)])
    for end,index in [('LEFT',0),('RIGHT',size-1)]:
        c=system['channels'][end];count=len(c['incoming']);forcing=np.zeros_like(values[:,index]);forcing[:,offset:offset+count]=c['incomingBoundaryData']
        boundary_residuals[end]=derivatives[:,index]-c['traceMap']@values[:,index]-forcing
        trace=c['right']@solution['modalAmplitudes'][end];trace[:,offset:offset+count]+=c['incomingValues'];trace_residuals[end]=trace-values[:,index];offset+=count
    f.require(b.norm(boundary_residuals)<1e-8 and b.norm(trace_residuals)<1e-8,'all actual finite boundary and evanescent amplitudes')
    inc,out=domain.phases(system,solution);s=out[:,None]*solution['boundaryAnchoredFluxBasisScattering']*inc[None,:];inverse=np.diag(1/out);current=inverse.conj().T@solution['outgoingChannelCurrent']@inverse
    ratio=np.diag(s.conj().T@current@s).real/solution['incomingFlux'];phase_residual=ratio-solution['outgoingFluxRatio'];f.require(b.norm(phase_residual)<1e-12,'full origin-current quotient')
    positions=np.linspace(-48.,48.,257);basis=f.polynomial_basis(positions,system['settings']['sourceBound'],size);fields=np.stack([basis@x[i*size:(i+1)*size] for i in range(5)])*inc[None,None,:]
    return {'originScattering':s,'originCurrent':current,'incomingPhase':inc,'outgoingPhase':out,'originCurrentResidual':phase_residual,'totalCurrentRatio':ratio,'positions':positions,'originFields':fields,'boundaryResiduals':boundary_residuals,'traceResiduals':trace_residuals,'independentCoefficients':other,'independentDifference':difference,'independentStatus':'new direct solve' if independent else 'accepted baseline independent-solve evidence retained in checkpoint','scope':'Each density family uses its own actual finite channel coordinates. Evanescent amplitudes remain boundary anchored.'}


def baseline_check(base,manifest):
    target=base/'cases'/BASELINE;system=f.unpickle(target/'finite/finite-system.pickle');solution=f.unpickle(target/'finite/finite-solution.pickle');basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');ends=f.unpickle(base/'boundary-cases'/BASELINE/'case-boundary.pickle');coefficients=f.unpickle(base/'interiors/cases'/BASELINE/'interior-matrices.pickle');direct=f.unpickle(base/'interiors/cases'/BASELINE/'direct-unsplit.pickle')
    direct=unsplit_arrays(direct)
    f.require(m.same(basis,system) and m.same(system['channels'],ends['finite']),'full accepted finite basis and actual boundary input')
    finish,prepare,join=finite_tail();a,rhs,original=prepare(target/'finite',direct['total'].copy(),ends['finite'],len(basis['nodes']),basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],solution['polynomialDerivativeResiduals'],solution['groups'])
    ar=(a-system['matrix'])/solution['rowScale'][:,None];rr=rhs-system['rhs'];actual=(a@solution['coefficients']-rhs)/solution['rowScale'][:,None]
    f.require(b.norm(ar)<1e-10 and b.norm(rr)==0 and b.norm(actual)<1e-9,'accepted full finite assembly and solution reuse')
    matrix,forcing=response.systems(coefficients,ends,basis);old=f.unpickle(target/'continuum/coefficient-systems.pickle');old_solution=f.unpickle(target/'continuum/coefficient-solutions.pickle')
    differences={g:matrix[g]-old['matrices'][g] for g in b.G};forcing_difference=b.subtract(forcing,old['rhs']);equation=b.subtract(b.J.multiply(matrix,old_solution['coefficients']),forcing)
    f.require(b.norm(differences)/(1+b.norm(matrix))<1e-10 and b.norm(forcing_difference)==0 and b.norm({g:v/old_solution['rowScale'][:,None] for g,v in equation.items()})<1e-9,'accepted four-grade systems and solution reuse')
    view=finite_view(system,solution,False);f.atomic_pickle(target/'finite/observable.pickle',view)
    # The actual boundary forcing changes when its incident sign is reversed.
    sign_mutation=-rhs-rhs;f.require(b.norm(sign_mutation)>1e-10,'actual one-sided incident sign mutation')
    f.atomic_pickle(base/'baseline-response-comparisons.pickle',{'finiteScaledMatrix':ar,'finiteRhs':rr,'finiteReusedSolutionEquation':actual,'continuumMatrix':differences,'continuumRhs':forcing_difference,'continuumReusedSolutionEquation':equation,'incidentSignMutation':sign_mutation})
    return {**join,'finiteMatrixResidual':b.norm(ar),'finiteReusedSolutionResidual':b.norm(actual),'continuumMatrixResidual':b.norm(differences),'incidentSignMutation':b.norm(sign_mutation),'baselineSolvesRepeated':0}


def construct_case(base,label,manifest):
    target=base/'cases'/label;finite_dir=target/'finite';cont_dir=target/'continuum';finite_dir.mkdir(parents=True);cont_dir.mkdir()
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');source=base/'interiors/cases'/label;coefficients=f.unpickle(source/'interior-matrices.pickle');direct=f.unpickle(source/'direct-unsplit.pickle');rows=f.unpickle(source/'row-matrices.pickle');binding=f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle');ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle');reference=f.unpickle(base/'accepted-modes/reference/modal.pickle')[0]
    direct=unsplit_arrays(direct)
    f.require(coefficients['settings']==rows['settings']==basis['settings']==binding['binding']['settings'],'complete finite setting identity')
    f.require(coefficients['fieldUnits']==ends['fieldUnits'] and coefficients['equationUnits']==ends['rowUnits'],'actual field and equation units')
    native=f.unpickle(source/'direct-native-cells.pickle')
    f.require(set(rows['rows'])==set(range(len(binding['binding']['bound']['rows']))) and len(binding['grades']['termJoins'])==native['terms'],'complete accepted case row and native cell census')
    f.require(b.norm(direct['total']-native['total'])/(1+b.norm(native['total']))<1e-10,'accepted direct unsplit/native cell matrix identity')
    finish,_,_=finite_tail();original_solution=f.unpickle(base/'cases'/BASELINE/'finite/finite-solution.pickle')
    # The tail only retains this field. It contains the complete accepted case
    # row packet, including its genuine new groups and every exact row reuse;
    # it is not a purported new integration or the legacy pilot group schema.
    result,_=finish(finite_dir,direct['total'].copy(),ends['finite'],len(basis['nodes']),basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],original_solution['polynomialDerivativeResiduals'],[rows])
    system=f.unpickle(finite_dir/'finite-system.pickle');view=finite_view(system,result,True);f.atomic_pickle(finite_dir/'observable.pickle',view)
    matrix,rhs=response.systems(coefficients,ends,basis);f.atomic_pickle(cont_dir/'coefficient-systems.pickle',{'matrices':matrix,'rhs':rhs})
    solved=response.solve(matrix,rhs);f.atomic_pickle(cont_dir/'coefficient-solutions.pickle',solved)
    channels=response.channels(solved,ends,basis,reference);f.atomic_pickle(cont_dir/'channel-response.pickle',channels);f.require(b.norm(channels['residuals'])<1e-8,'complete current/phase/normalization response')
    ratio=sp.Rational(manifest['input']['parameters']['W_0'])/sp.Rational(manifest['input']['parameters']['L_W']);flux=response.open_flux(channels,ratio);remainders={}
    for eta,sigma in ((0.01,0.001),(0.005,0.0005),(0.0025,0.00025)):
        a=b.evaluate(matrix,eta,sigma);rhs_at=b.evaluate(rhs,eta,sigma);actual=la.solve(a,rhs_at);predicted=b.evaluate(solved['coefficients'],eta,sigma);difference=actual-predicted
        remainders[eta,sigma]={'direct':actual,'retained':predicted,'difference':difference,'maximumReferenceFrame':b.norm(difference),'directEquationResidual':a@actual-rhs_at}
    f.atomic_pickle(cont_dir/'formal-remainders.pickle',remainders)
    packet={'solve':solved,'response':channels,'flux':flux,'remainders':remainders,'size':len(basis['nodes']),'fieldUnits':ends['fieldUnits'],'rowUnits':ends['rowUnits'],'currentUnit':ends['currentUnit'],'ratio':ratio,'settings':basis['settings'],'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'dimensionState':coefficients['dimensionState'],'scope':'Actual case finite-domain independent-grade response; positive regulator, approximate modal boundaries and omitted parent pure-second-order terms.'}
    f.atomic_pickle(cont_dir/'continuum-response.pickle',packet)
    summary={'case':label,'rows':len(rows['rows']),'terms':len(binding['grades']['termJoins']),'sources':len(binding['binding']['jets']),'unknowns':len(system['matrix']),'incidentColumns':system['rhs'].shape[1],'finiteRank':result['rank'],'finiteCondition':result['balancedCondition'],'finiteScaledResidual':b.norm(result['scaledEquationResidual']),'finiteIndependentDifference':b.norm(view['independentDifference']),'finiteBoundaryResidual':b.norm(view['boundaryResiduals']),'finiteCurrentRatios':view['totalCurrentRatio'].tolist(),'continuumRank':solved['rank'],'continuumCondition':solved['condition'],'continuumScaledResidual':b.norm(solved['scaledResidual']),'continuumIndependentDifference':b.norm(solved['independentDifference']),'continuumMapResidual':b.norm(channels['residuals']),'mixedForcingMutation':b.norm(solved['mixedForcingMutation']),'formalRemainderMaxima':{str(k):v['maximumReferenceFrame'] for k,v in remainders.items()},'newQuadratureNodes':0}
    f.save(target/'response-checks.json',summary)
    return summary


def emit_saved_case(base,label,manifest):
    directory=base/'cases'/label/'continuum';result=f.unpickle(directory/'continuum-response.pickle');binding=f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle')
    r,adapter,packets=matrices.context(base/'interiors',label,binding,manifest['input'])
    f.engine.PHYSICAL_METADATA.dimensions.__dict__.update(result['dimensionState'])
    emitted=emitter(label)(directory,result,r,manifest['sourceFiles'],manifest['inputPackets'],f.digest(directory/'continuum-response.pickle'))
    f.save(directory/'emission-checks.json',emitted);return emitted


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--focused',action='store_true');ap.add_argument('--resume-focused',type=Path);args=ap.parse_args();start=time.monotonic()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    manifest,labels=load(base,args.resume_focused);join=json.loads((args.resume_focused/'checks.json').read_text())['baselineReplay'] if args.resume_focused else baseline_check(base,manifest)
    counts={};emissions={}
    if not args.focused:
        for label in labels:
            if label==BASELINE:continue
            with (base/'progress.jsonl').open('a') as stream:stream.write(json.dumps({'case':label,'phase':'solve','wallSeconds':time.monotonic()-start})+'\n')
            counts[label]=construct_case(base,label,manifest);f.save(base/'case-inventory.json',counts)
        # Finish all numerical packets before output; an output failure never requires another solve.
        for label in labels:
            if label!=BASELINE:emissions[label]=emit_saved_case(base,label,manifest)
        same_basis={}
        for density in ('RHO4_CONSTANT','RHOBR_CONSTANT'):
            lab='LAB_HELD__'+density;mat='MATERIAL_ADVECTED__'+density
            le=f.unpickle(base/'boundary-cases'/lab/'case-boundary.pickle');me=f.unpickle(base/'boundary-cases'/mat/'case-boundary.pickle');f.require(m.same(le['ends'],me['ends']) and m.same(le['finite'],me['finite']),'same actual channel coordinates before anchoring comparison')
            lv=f.unpickle(base/'cases'/lab/'finite/observable.pickle');mv=f.unpickle(base/'cases'/mat/'finite/observable.pickle');lc=f.unpickle(base/'cases'/lab/'continuum/channel-response.pickle');mc=f.unpickle(base/'cases'/mat/'continuum/channel-response.pickle')
            same_basis[density]={'finiteAmplitudeDifference':mv['originScattering']-lv['originScattering'],'finiteCurrentDifference':mv['totalCurrentRatio']-lv['totalCurrentRatio'],'finiteFieldDifference':mv['originFields']-lv['originFields'],'continuumFluxAmplitudeDifferences':b.subtract(mc['fluxOriginScattering'],lc['fluxOriginScattering']),'scope':'Anchoring comparison at fixed density and identical end coordinates; no cross-density S-matrix subtraction.'}
        f.atomic_pickle(base/'anchoring-comparisons.pickle',same_basis)
        with (base/'full.out').open('xb') as stream:
            for label in labels:stream.write((base/'cases'/label/'continuum/full.out').read_bytes())
        tags=set();keys=set()
        for line in response.grades.decoded_lines(base/'full.out'):
            tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in tags,'four-case output tag uniqueness');tags.add(tag)
            if tag.endswith('_WRITE_KEYS') and not tag.startswith('PY_S11CD_METADATA_'):
                values={str(k):str(v) for k,v in response.grades._restore(body)}
                f.require(len(values)==len(set(values.values())) and not keys&set(values.values()),'distinct actual case export keys');keys.update(values.values())
        manifest['combinedOutput']={'tags':len(tags),'keys':len(keys),'caseParts':{label:f.digest(base/'cases'/label/'continuum/full.out') for label in labels},'scope':'Lossless concatenation of four individually validated continuum transcripts; each original emission index remains local to its own case part.'}
        f.atomic_pickle(base/'remaining-case-response.pickle',{'cases':counts,'baseline':str(base/'cases'/BASELINE),'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'source pre/post identity')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original pre/post identity')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'unchanged copied operand')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts}
    checks={**manifest,'mode':'focused' if args.focused else 'construct','status':'PASSED_BASELINE_RESPONSE_WIRING' if args.focused else 'COMPLETED_FOUR_CASE_RESPONSES','baselineReplay':join,'cases':counts,'emissions':emissions,'artifacts':artifacts,'newFiniteCases':len(counts),'newContinuumCases':len(counts),'newQuadratureNodes':0,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
