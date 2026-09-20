#!/usr/bin/env python3
"""Solve missing derivative controls from saved complete boundary systems."""
import argparse, ast, builtins, copy, dis, gc, json, resource, shutil, signal, time
from pathlib import Path
import numpy as np
import S11c_d_remaining_case_first_jet_response_inputs as inputs
import S11c_d_remaining_case_profile_response as profile
import S11c_d_remaining_case_flux as flux_helper

f,m,b=inputs.f,inputs.m,inputs.b
native=inputs.response
current=flux_helper.native
BASELINE=inputs.BASELINE
FCP=f.M/'S11c_d_remaining_case_first_jet_response_focused.json'
CCP=f.M/'S11c_d_remaining_case_flux_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_first_jet_response_production_plan.md'
SCOPE=('Four actual finite-contrast and independent-grade first-w-derivative controls in their own case coordinates. '
       'Four finite and three continuum solves are new; the baseline continuum is reused. '
       'Positive regulator, approximate modal boundaries and omitted parent pure-second-order terms remain. '
       'One-sided closed-operator sensitivity, not a consistent new profile or isolated advection channel. '
       'Full transcript validation/publication follows these saved numerical results.')


def finite_solver():
    source=native.function(ast.parse(Path(f.__file__).read_text()),'construct')
    start=next(i for i,n in enumerate(source.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
               and ast.unparse(n.value.func)=='atomic_pickle' and 'finite-system.pickle' in ast.unparse(n))
    body=copy.deepcopy(source.body[start:])
    names=('base','matrix','rhs','original','channels','size','derivative','x','local_matrix','setting','pins','operands','differentiation','groups')
    node=ast.FunctionDef(name='solve_prepared_finite',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    f.require(ast.dump(ast.Module(body=node.body,type_ignores=[]))==ast.dump(ast.Module(body=source.body[start:],type_ignores=[])), 'whole unchanged native finite save/solve/current tail')
    fn=native.compile_function(node,vars(f))
    missing={v.argval for v in dis.get_instructions(fn) if v.opname=='LOAD_GLOBAL' and v.argval not in fn.__globals__ and not hasattr(builtins,v.argval)}
    f.require(not missing,'complete prepared finite native namespace')
    return fn,{'wholeNativeFiniteSolveTail':True,'firstStatement':start,'changedNumericalStatements':0,'nativeFileSha256':f.digest(Path(f.__file__))}


def continuum_solver():
    source=native.function(ast.parse(Path(native.__file__).read_text()),'construct_case')
    start=next(i for i,n in enumerate(source.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='solved')
    names=('base','label','manifest','target','cont_dir','basis','ends','reference','coefficients','matrix','rhs','system','result','view','rows','binding')
    node=ast.FunctionDef(name='solve_prepared_continuum',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=copy.deepcopy(source.body[start:]),decorator_list=[])
    f.require(ast.dump(ast.Module(body=node.body,type_ignores=[]))==ast.dump(ast.Module(body=source.body[start:],type_ignores=[])), 'whole unchanged continuum solve/maps/flux/formal-check tail')
    fn=native.compile_function(node,vars(native))
    missing={v.argval for v in dis.get_instructions(fn) if v.opname=='LOAD_GLOBAL' and v.argval not in fn.__globals__ and not hasattr(builtins,v.argval)}
    f.require(not missing,'complete prepared continuum native namespace')
    return fn,{'wholeNativeContinuumSolveTail':True,'firstStatement':start,'changedNumericalStatements':0,'nativeFileSha256':f.digest(Path(native.__file__))}


def comparison_prefix():
    source=native.function(ast.parse(Path(profile.__file__).read_text()),'compare_case')
    stop=next(i for i,n in enumerate(source.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='coefficients')
    original=copy.deepcopy(source.body[:stop]);body=copy.deepcopy(original)
    class Rename(ast.NodeTransformer):
        def __init__(self,old,new):self.old,self.new,self.count=old,new,0
        def visit_Constant(self,n):
            if n.value==self.old:self.count+=1;return ast.Constant(self.new)
            return n
    rename=Rename('form-comparisons.pickle','first-jet-comparisons.pickle');body=[rename.visit(n) for n in body]
    undo=Rename('first-jet-comparisons.pickle','form-comparisons.pickle');back=[undo.visit(copy.deepcopy(n)) for n in body]
    f.require(rename.count==1 and ast.dump(ast.Module(body=back,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[])), 'exact generic own-case comparison body with one packet address change')
    names=('base','label','manifest')
    node=ast.FunctionDef(name='compare_prepared_case',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),
        body=body+[ast.Return(ast.Tuple([ast.Name(n,ast.Load()) for n in ('comparisons','new','old','current','original')],ast.Load()))],decorator_list=[])
    return native.compile_function(node,vars(profile)),{'wholeOwnCaseComparisonPrefix':True,'packetAddressChanges':1,'changedNumericalStatements':0,'nativeFileSha256':f.digest(Path(profile.__file__))}


def load(base,focus):
    root,checks,cp=inputs.provenance.accepted(FCP,'ACCEPTED_FIRST_JET_RESPONSE_INPUTS')
    f.require(root==focus,'exact accepted prepared-system directory')
    cr,cc,_=inputs.provenance.accepted(CCP,'PUBLISHED_ANNEX_VERIFIED')
    pins=dict(checks['sourceFiles'])
    for name,value in cc['sourceFiles'].items():f.require(name not in pins or pins[name]==value,'same consumed current source');pins[name]=value
    for path in (Path(__file__),PLAN,FCP,CCP,Path(profile.__file__),Path(flux_helper.__file__),Path(current.__file__)):
        pins[str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':dict(checks['inputPackets']),'copiedInputs':{},'input':checks['input'],'settings':checks['settings'],'scope':SCOPE,
              'completedFocusedReuse':{'directory':str(root),'checksSha256':cp['checksSha256'],'artifacts':len(checks['artifacts'])}}
    for rr,ch in ((root,checks),(cr,cc)):
        manifest['inputPackets'][str(rr/'checks.json')]=f.digest(rr/'checks.json')
        for n,v in ch['inputPackets'].items():f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n]==v,'same original consumed input');manifest['inputPackets'][n]=v
    for n,v in checks['artifacts'].items():m.retain(root/n,base/n,manifest,v['sha256'])
    labels=tuple(checks['cases'])
    for label in labels:
        # Complete prepared systems are copied, not rebuilt in the numerical run.
        m.retain(root/'preparation'/label/'continuum-prepared.pickle',base/'cases'/label/'continuum/coefficient-systems.pickle',manifest)
        for name in ('channel-selectors','bulk-kinematics','continuum-currents','open-channel-currents','finite-end-currents'):
            path='cases/'+label+'/continuum/'+name+'.pickle';m.retain(cr/path,base/'unchanged-currents'/path,manifest,cc['artifacts'][path]['sha256'])
        source='case-response/'+label+'/continuum-response.pickle'
        f.require(f.digest(cr/source)==f.digest(base/'unchanged-response/cases'/label/'continuum/continuum-response.pickle')==cc['artifacts'][source]['sha256'],'actual original current response input')
        manifest['inputPackets'][str(cr/source)]=cc['artifacts'][source]['sha256']
        address='boundary-cases/'+label+'/case-boundary.pickle'
        f.require(f.digest(cr/address)==f.digest(base/address)==cc['artifacts'][address]['sha256'],'actual complete current end input')
        manifest['inputPackets'][str(cr/address)]=cc['artifacts'][address]['sha256']
    for source,target in (('first-jet-solutions','coefficient-solutions'),('first-jet-channels','channel-response'),('first-jet-response','continuum-response')):
        m.retain(root/'accepted-first-jet'/(source+'.pickle'),base/'cases'/BASELINE/'continuum'/(target+'.pickle'),manifest)
    m.retain(root/'accepted-first-jet/full.out',base/'cases'/BASELINE/'continuum/full.out',manifest)
    jc=json.loads(inputs.JCP.read_text());m.retain(Path(jc['runDirectory'])/'checks.json',base/'cases'/BASELINE/'continuum/checks.json',manifest,jc['checksSha256'])
    for n,v in pins.items():
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dest);f.require(f.digest(dest)==v,'frozen numerical response source')
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');manifest['settings']=inputs.provenance.restore_settings(manifest['settings'],basis['settings'])
    f.save(base/'inputs.json',manifest);return manifest,labels


def wiring(base,manifest,labels):
    finite,fp=finite_solver();_,cp=continuum_solver();_,compare=comparison_prefix()
    pre=json.loads((base/'preflight.json').read_text())
    for name,sha in pre['continuumFunctions'].items():f.require(inputs.matrices.binding.native_body(getattr(native.response,name))==sha,'unchanged accepted continuum numerical functions')
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');routes={}
    class Captured(Exception):pass
    for label in labels:
        prepared=f.unpickle(base/'preparation'/label/'finite-prepared.pickle');ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle')
        direct=native.unsplit_arrays(f.unpickle(base/'interiors/cases'/label/'direct-unsplit.pickle'))
        rows=f.unpickle(base/'interiors/cases'/label/'row-matrices.pickle');old=f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-solution.pickle')
        held={}
        def capture(path,payload):held.update(path=path,payload=payload);raise Captured()
        scope=dict(finite.__globals__,atomic_pickle=capture)
        probe=__import__('types').FunctionType(finite.__code__,scope,finite.__name__)
        try:probe(base/'cases'/label/'finite',prepared['matrix'],prepared['rhs'],prepared['unreplacedOperator'],ends['finite'],129,basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],old['polynomialDerivativeResiduals'],[rows])
        except Captured:pass
        f.require(held['path']==base/'cases'/label/'finite/finite-system.pickle','native first write address before solve')
        value=held['payload']
        for name in ('matrix','rhs','unreplacedOperator'):f.require(m.same(value[name],prepared[name]),'whole prepared finite array in native first write')
        f.require(m.same(value['channels'],ends['finite']) and m.same(value['settings'],basis['settings']),'actual native first-write end and settings')
        changed=value['rhs'].copy();changed[0,0]+=1
        f.require(not m.same(changed,prepared['rhs']),'actual prepared forcing mutation rejected')
        routes[label]={'firstWriteCapturedBeforeSolve':True,'preparedMatrixShape':list(value['matrix'].shape),'preparedForcingShape':list(value['rhs'].shape),'changedForcingRejected':True}
        del prepared,ends,direct,rows,old,held,value,probe;gc.collect()
    result={'finiteTail':fp,'continuumTail':cp,'comparisonPrefix':compare,'routes':routes,'newSolves':0,'newQuadratureNodes':0}
    f.save(base/'solver-wiring.json',result);return result


def construct(base,label,manifest):
    target=base/'cases'/label;fd=target/'finite';fd.mkdir();cd=target/'continuum'
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle')
    direct=native.unsplit_arrays(f.unpickle(base/'interiors/cases'/label/'direct-unsplit.pickle'));rows=f.unpickle(base/'interiors/cases'/label/'row-matrices.pickle')
    selected=f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle');co=f.unpickle(base/'interiors/cases'/label/'interior-matrices.pickle')
    saved=f.unpickle(base/'preparation'/label/'finite-prepared.pickle');old=f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-solution.pickle')
    finish,_=finite_solver()
    result,_=finish(fd,saved['matrix'],saved['rhs'],saved['unreplacedOperator'],ends['finite'],129,basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],old['polynomialDerivativeResiduals'],[rows])
    system=f.unpickle(fd/'finite-system.pickle');view=native.finite_view(system,result,True);f.atomic_pickle(fd/'observable.pickle',view)
    if label==BASELINE:
        continuum=f.unpickle(cd/'continuum-response.pickle')
        summary={'case':label,'rows':len(rows['rows']),'terms':len(selected['grades']['termJoins']),'sources':len(selected['binding']['jets']),
                 'unknowns':645,'incidentColumns':4,'finiteRank':result['rank'],'finiteCondition':result['balancedCondition'],'finiteScaledResidual':b.norm(result['scaledEquationResidual']),
                 'finiteIndependentDifference':b.norm(view['independentDifference']),'finiteBoundaryResidual':b.norm(view['boundaryResiduals']),'finiteCurrentRatios':view['totalCurrentRatio'].tolist(),
                 'continuumRank':continuum['solve']['rank'],'continuumScaledResidual':b.norm(continuum['solve']['scaledResidual']),
                 'newFiniteControl':True,'newContinuumControl':False,'reusedHistoricalContinuum':True,'formalRemainders':'Historical baseline first-derivative control did not compute direct formal-point solves; none is invented here.'}
        f.save(target/'response-checks.json',summary)
    else:
        p=f.unpickle(cd/'coefficient-systems.pickle');ref=f.unpickle(base/'accepted-modes/reference/modal.pickle')[0];complete,_=continuum_solver()
        summary=complete(base,label,manifest,target,cd,basis,ends,ref,co,p['matrices'],p['rhs'],system,result,view,rows,selected)
        summary.update(newFiniteControl=True,newContinuumControl=True,reusedHistoricalContinuum=False)
        f.save(target/'response-control-disposition.json',summary)
    return summary


def finite_open_current(system,solution):
    result={};offset=0
    for end in ('LEFT','RIGHT'):
        channel=system['channels'][end];metric=channel['current'];open_rows=[(i,v) for i,v in enumerate(channel['outgoing']) if v['kind']=='open'];outgoing=[v['MATRIX_COLUMN'] for i,v in open_rows];incoming=[v['MATRIX_COLUMN'] for v in channel['incoming']]
        f.require(len(set(outgoing+incoming))==len(outgoing)+len(incoming)==metric.shape[0] and set(outgoing+incoming)==set(range(len(metric))),'complete finite open-current matrix column census')
        out=np.zeros((len(metric),4),complex);inc=np.zeros_like(out);out[outgoing]=np.vstack([solution['modalAmplitudes'][end][i] for i,v in open_rows])
        inc[np.ix_(incoming,range(offset,offset+len(incoming)))]=np.eye(len(incoming));offset+=len(incoming)
        amplitude=out+inc;full=amplitude.conj().T@metric@amplitude;outward=out.conj().T@metric@out;inward=inc.conj().T@metric@inc;cross=out.conj().T@metric@inc+inc.conj().T@metric@out
        scalar=np.array([[sum(amplitude[i,a].conjugate()*metric[i,j]*amplitude[j,z] for i in range(len(metric)) for j in range(len(metric))) for z in range(4)] for a in range(4)])
        residuals={'direct':scalar-full,'decomposition':full-outward-inward-cross}
        f.require(b.norm(residuals)<1e-8,'full finite open-current contraction including incoming/outgoing cross terms')
        result[end]={'amplitude':amplitude,'metric':metric,'full':full,'outgoing':outward,'incoming':inward,'interference':cross,'residuals':residuals,'coordinate':'actual boundary anchored open-current basis','closedMatchingAmplitudes':{i:solution['modalAmplitudes'][end][i] for i,v in enumerate(channel['outgoing']) if v['kind']=='evanescent'},'scope':'This native finite metric covers open directions only. Closed amplitudes are retained, not assigned zero current. Full closed/cross-mode current bookkeeping is computed separately for the retained continuum response.'}
    f.require(offset==4,'complete finite incident columns');return result


def compare(base,label,manifest):
    prefix,_=comparison_prefix();differences,new,old,value,original=prefix(base,label,manifest)
    ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle');directory=base/'cases'/label
    keys=('openBoundaryScattering','openOriginScattering','fieldOriginScattering','fluxOriginScattering','incomingCurrentOrigin','outgoingCurrentOrigin')
    channel={key:b.subtract(value['response'][key],original['response'][key]) for key in keys}
    operator=f.unpickle(base/'interiors/cases'/label/'interior-matrices.pickle')
    eta=float(operator['gradeOrigin'][operator['generators'][1]]);sigma=float(operator['gradeOrigin'][operator['generators'][2]])
    f.require(abs(float(value['ratio'])-float(manifest['input']['parameters']['W_0'])/float(manifest['input']['parameters']['L_W']))<1e-15,'actual accepted homotopy ratio')
    del operator
    def quotient(packet):
        s=b.evaluate(packet['response']['openOriginScattering'],eta,sigma);ji=b.evaluate(packet['response']['incomingCurrentOrigin'],eta,sigma);jo=b.evaluate(packet['response']['outgoingCurrentOrigin'],eta,sigma)
        direct=np.array([sum(s[i,a].conjugate()*jo[i,j]*s[j,a] for i in range(4) for j in range(4))/ji[a,a] for a in range(4)])
        ratio=np.diag(s.conj().T@jo@s)/np.diag(ji);f.require(b.norm(ratio-direct)<1e-12,'independent actual retained-polynomial physical current quotient')
        return {'amplitude':s,'incomingCurrent':ji,'outgoingCurrent':jo,'ratio':ratio,'independent':direct,'residual':ratio-direct}
    retained={'reversed':quotient(value),'baseline':quotient(original)};retained['difference']=retained['reversed']['ratio']-retained['baseline']['ratio']
    selection=f.unpickle(base/'unchanged-currents/cases'/label/'continuum/channel-selectors.pickle');prior=f.unpickle(base/'unchanged-currents/cases'/label/'continuum/continuum-currents.pickle')
    reference=f.unpickle(base/'accepted-modes/reference/modal.pickle')[0]
    f.require(m.same(selection,current.selectors(value,ends,reference)),'complete source-derived selected directions and classifiers')
    for key in ('incomingOriginPhase','outgoingOriginPhase','fieldIncomingMap','fieldOpenOutgoingMap'):
        f.require(m.same(value['response'][key],original['response'][key]),'exact unchanged current metric coordinates')
    # The old actual metric was derived from the identical end and phase/field
    # operands. Only the reversed response amplitudes need new contractions.
    metrics=prior['metrics'];opened,mutation=current.construct_open(value,metrics,selection)
    closed={}
    for end in ('LEFT','RIGHT'):
        ix=[i for row in selection['census'][end] if row['kind']=='evanescent' and row['direction']=='outgoing' for i in row['columns']]
        closed[end]=current.amplitude_bookkeeping({g:a[ix] for g,a in value['response']['outgoingFieldAmplitudes'][end].items()})
    end_current=current.end_currents(value,ends)
    bookkeeping={'selection':selection,'metrics':metrics,'metricResiduals':prior['metricResiduals'],'open':opened,'closedAmplitudes':closed,'endCurrents':end_current,'bulkDomains':prior['bulkDomains'],'currentMutation':mutation,
                 'currentUnit':value['currentUnit'],'ratio':value['ratio'],'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'settings':value['settings'],'dimensionState':value['dimensionState'],
                 'scope':'Actual reversed amplitudes and inherited full slab/bulk normal currents, including closed matching and cross terms. Depth-integrated normal current is not depth escape.'}
    f.atomic_pickle(directory/'continuum/continuum-currents.pickle',bookkeeping)
    finite={}
    for name,root in (('reversed',directory),('baseline',base/'unchanged-response/cases'/label)):
        finite[name]=finite_open_current(f.unpickle(root/'finite/finite-system.pickle'),f.unpickle(root/'finite/finite-solution.pickle'))
    finite['differences']={end:{k:finite['reversed'][end][k]-finite['baseline'][end][k] for k in ('full','outgoing','incoming','interference')} for end in ('LEFT','RIGHT')}
    f.atomic_pickle(directory/'finite/open-end-currents.pickle',finite)
    sensitivity={'case':label,'evaluatedGradeOrigin':{'eta':eta,'sigma':sigma},'channelCoefficients':channel,'evaluatedChannels':{key:b.evaluate(series,eta,sigma) for key,series in channel.items()},'retainedPolynomialCurrent':retained,
                 'finiteAndCommonGrid':differences,'fieldCoefficientDifference':b.subtract(value['solve']['coefficients'],original['solve']['coefficients']),
                 'openCurrentDifference':{name:{part:current.subtract(opened[name]['currentParts'][part],prior['open'][name]['currentParts'][part]) for part in ('slab','bulk','total')} for name in opened},
                 'endCurrentDifference':{end:{part:{kind:current.subtract(end_current[end]['parts'][part][kind],prior['endCurrents'][end]['parts'][part][kind]) for kind in ('full','outgoing','incoming','interference')} for part in end_current[end]['parts']} for end in ('LEFT','RIGHT')},
                 'fieldUnits':value['fieldUnits'],'currentUnit':value['currentUnit'],'settings':value['settings'],'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':SCOPE,
                 'materialComparison':{'status':'accepted baseline comparison retained in accepted-first-jet/first-jet-comparisons.pickle' if label==BASELINE else 'not yet constructed for this case','zeroSubstitute':False}}
    f.atomic_pickle(directory/'first-jet-sensitivity.pickle',sensitivity)
    summary={'finiteAmplitudeChange':b.norm(differences['finiteScattering']),'finiteCurrentChange':b.norm(differences['finiteCurrent']),'finiteFieldChange':b.norm(differences['finiteFields']),
             'continuumAmplitudeChange':b.norm(sensitivity['evaluatedChannels']['fluxOriginScattering']),'retainedPolynomialCurrentChange':b.norm(retained['difference']),
             'openCurrentResidual':max(b.norm(v['residuals']) for v in opened.values()),'endCurrentResidual':max(b.norm(v['residuals']) for v in end_current.values()),'currentMutation':b.norm(mutation),
             'finiteOpenEndCurrentResidual':max(b.norm(finite[name][end]['residuals']) for name in ('reversed','baseline') for end in ('LEFT','RIGHT')),
             'amplitudeResolutionGoal':1e-4,'currentResolutionGoal':1e-6,'sameOwnCaseCoordinates':True,'scope':SCOPE}
    f.save(directory/'sensitivity-checks.json',summary);return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--resume-focused',type=Path,required=True);p.add_argument('--wiring-only',action='store_true');args=p.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    manifest,labels=load(base,args.resume_focused);wire=wiring(base,manifest,labels);solves={};comparisons={}
    if not args.wiring_only:
        for label in labels:
            solves[label]=construct(base,label,manifest);f.save(base/'response-inventory.json',solves);gc.collect()
        # Every new solution is durable before any comparison or output work.
        for label in labels:
            comparisons[label]=compare(base,label,manifest);f.save(base/'comparison-inventory.json',comparisons);gc.collect()
        f.atomic_pickle(base/'remaining-case-first-jet-response.pickle',{'cases':solves,'comparisons':comparisons,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':SCOPE})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'source/current/frozen pre/post')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original input pre/post')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'copied prepared/accepted operand pre/post')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    result={**manifest,'status':'PASSED_PREPARED_FIRST_JET_SOLVER_WIRING' if args.wiring_only else 'COMPLETED_FIRST_JET_NUMERICAL_RESPONSES','solverWiring':wire,'cases':solves,'comparisons':comparisons,
            'newFiniteControlSolves':len(solves),'newContinuumControlSolves':sum(v['newContinuumControl'] for v in solves.values()),'reusedBaselineContinuum':True,'newQuadratureNodes':0,'newInteriorAssemblies':0,'newModes':0,'newCurrentClosures':0,
            'newCurrentContractions':len(comparisons),'physicalTranscriptEmitted':False,'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
