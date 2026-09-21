#!/usr/bin/env python3
"""Saved-system material-coordinate controls, with unchanged native solvers."""
import argparse,ast,copy,gc,json,resource,shutil,signal,time,types
from pathlib import Path
import numpy as np
import S11c_d_remaining_case_coordinate_response_inputs as inputs
import S11c_d_remaining_case_first_jet_response as solver
import S11c_d_remaining_case_flux as flux

f,m,b=inputs.f,inputs.m,inputs.b
native=inputs.response;current=flux.native;matrices=inputs.matrices
BASELINE=inputs.BASELINE
FCP=f.M/'S11c_d_remaining_case_coordinate_response_focused.json'
CCP=f.M/'S11c_d_remaining_case_flux_checkpoint.json'
JCP=f.M/'S11c_d_remaining_case_first_jet_response_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_coordinate_response_production_plan.md'
SCOPE=('Own-case material-coordinate controls constructed in common Eulerian coordinates before solving. '
       'Four new finite controls and three new continuum controls; historical baseline continuum reused unchanged. '
       'Positive regulator, approximate boundaries, unresolved tiny signals and omitted parent pure-second-order scope. '
       'Finite open-only current and separate closed matching amplitudes; full seven-direction continuum currents. '
       'No post-hoc S conjugation, unlike-density subtraction, isolated advection or c2 closure claim.')


def load(base,focus,resume=None):
    fr,fc,fp=matrices.accepted(FCP,'ACCEPTED_CASE_MATERIAL_RESPONSE_INPUTS')
    f.require(fr==focus,'exact accepted saved response systems')
    cr,cc,cp=matrices.accepted(CCP,'PUBLISHED_ANNEX_VERIFIED')
    jr,jc,jp=matrices.accepted(JCP,'ACCEPTED_FIRST_JET_NUMERICAL_RESPONSES')
    pins=dict(fc['sourceFiles']);operands=dict(fc['inputPackets'])
    for root,ch in ((fr,fc),(cr,cc),(jr,jc)):
        for n,sha in ch['sourceFiles'].items():
            f.require(n not in pins or pins[n]==sha,'unchanged shared scientific source');pins[n]=sha
        for n,sha in ch['inputPackets'].items():
            f.require(n not in operands or operands[n]==sha,'unchanged shared input');operands[n]=sha
        operands[str(root/'checks.json')]=f.digest(root/'checks.json')
    for path in (Path(__file__),PLAN,FCP,CCP,JCP,Path(inputs.__file__),Path(solver.__file__),Path(current.__file__)):
        pins[str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'copiedInputs':{},
        'settings':fc['settings'],'input':fc['input'],'scope':SCOPE,'materialUnitSupplement':fc['materialUnitSupplement'],
        'completedFocusedReuse':{'directory':str(fr),'checksSha256':fp['checksSha256'],'artifacts':len(fc['artifacts'])}}
    labels=tuple(fc['cases'])
    if resume:
        old=json.loads((resume/'checks.json').read_text());run=resume.parent
        inv=json.loads((run/'coordinate_response_wiring.invocation.json').read_text());guard=json.loads((run/'resource-guard/outcome.json').read_text())
        f.require(inv['exitCode']==guard['exitCode']==guard['childOutcome']['exitCode']==0 and inv['stderrBytes']==guard['stderrBytes']==0 and guard['limitsVerified'] and guard['childOutcome']['guardReason'] is None,'actual clean solver-routing preflight')
        f.require((run/'coordinate_response_wiring.stdout').read_bytes()==(resume/'checks.json').read_bytes(),'preflight final checks/stdout identity')
        f.require(old['status']=='PASSED_MATERIAL_RESPONSE_SOLVER_WIRING' and old['sourceFiles']==pins,'exact numerical helper and completed preflight')
        for n,v in old['artifacts'].items():m.retain(resume/n,base/n,manifest,v['sha256'])
        manifest['inputPackets'].update(old['inputPackets']);manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json')
        manifest['completedSolverPreflightReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'artifacts':len(old['artifacts'])}
    else:
        for n,v in fc['artifacts'].items():m.retain(fr/n,base/n,manifest,v['sha256'])
        m.retain(fr/'checks.json',base/'accepted-focus-checks.json',manifest,fp['checksSha256'])
        vr=Path(fp['validator']['runDirectory'])
        m.retain(vr/'checks.json',base/'accepted-focus-validation.json',manifest,fp['validator']['checksSha256'])
        for label in labels:
            m.retain(fr/'preparation'/label/'continuum-prepared.pickle',base/'cases'/label/'continuum/coefficient-systems.pickle',manifest)
            # Original Eulerian currents and their complete input coordinates are operands.
            for name in ('channel-selectors','bulk-kinematics','continuum-currents'):
                address='cases/'+label+'/continuum/'+name+'.pickle'
                m.retain(cr/address,base/'unchanged-currents'/address,manifest,cc['artifacts'][address]['sha256'])
            for address,local in (('case-response/'+label+'/continuum-response.pickle','unchanged-response/cases/'+label+'/continuum/continuum-response.pickle'),
                                  ('boundary-cases/'+label+'/case-boundary.pickle','boundary-inputs/boundary-cases/'+label+'/case-boundary.pickle')):
                sha=cc['artifacts'][address]['sha256'];f.require(f.digest(cr/address)==f.digest(base/local)==sha,'actual Eulerian current response/end inputs')
                manifest['inputPackets'][str(cr/address)]=sha
            prior='cases/'+label+'/finite/open-end-currents.pickle'
            m.retain(jr/prior,base/'unchanged-finite-currents'/(label+'.pickle'),manifest,jc['artifacts'][prior]['sha256'])
            for name in ('finite-system','finite-solution'):
                address='unchanged-response/cases/'+label+'/finite/'+name+'.pickle'
                f.require(f.digest(jr/address)==f.digest(base/address)==jc['artifacts'][address]['sha256'],'exact original finite current operands before reuse')
                manifest['inputPackets'][str(jr/address)]=jc['artifacts'][address]['sha256']
        for source,target in (('material-coefficient-solutions','coefficient-solutions'),('material-channel-response','channel-response'),('coordinate-response','continuum-response')):
            m.retain(fr/'accepted-coordinate'/(source+'.pickle'),base/'cases'/BASELINE/'continuum'/(target+'.pickle'),manifest)
        for name in ('full.out','checks.json'):
            m.retain(fr/'accepted-coordinate'/name,base/'cases'/BASELINE/'continuum'/name,manifest)
    for n,sha in pins.items():
        path=base/'source'/n;path.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,path)
        f.require(f.digest(path)==sha,'frozen numerical helper and physical sources')
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle')
    f.require(json.loads(json.dumps(basis['settings']))==manifest['settings'],'native settings and exact JSON view')
    manifest['settings']=basis['settings'];f.save(base/'inputs.json',manifest);matrices.hash_check(base,manifest)
    return manifest,labels


def wiring(base,manifest,labels):
    finite,fj=solver.finite_solver();_,cj=solver.continuum_solver()
    tree=ast.parse(Path(__file__).read_text())
    functions={v.name:v for v in tree.body if isinstance(v,ast.FunctionDef)}
    call=lambda fn,name:next(v for v in ast.walk(functions[fn]) if isinstance(v,ast.Call) and isinstance(v.func,ast.Name) and v.func.id==name)
    checked=copy.deepcopy(call('wiring','probe'));actual=call('construct','finish')
    class Names(ast.NodeTransformer):
        def visit_Name(self,node):return ast.copy_location(ast.Name({'probe':'finish','p':'saved'}.get(node.id,node.id),node.ctx),node)
    checked=Names().visit(checked)
    # The production fd local is exactly this saved first-write directory.
    checked.args[0]=ast.Name('fd',ast.Load())
    f.require(ast.dump(checked)==ast.dump(actual),'actual production finite call matches captured first-write operands')
    saved=json.loads((base/'native-response-wiring.json').read_text())
    f.require(fj==saved['finiteSolver'],'unchanged accepted finite solve suffix')
    for name,sha in saved['nativeContinuum'].items():f.require(matrices.body(getattr(native.response,name))==sha,'unchanged native continuum function')
    joins={'finite':fj,'continuum':cj,'actualFiniteCallMatchesCapture':True,'finiteCurrent':matrices.body(solver.finite_open_current),
        'currentFunctions':{n:matrices.body(getattr(current,n)) for n in ('selectors','open_metrics','construct_open','end_currents','amplitude_bookkeeping','direct_quadratic')}}
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');routes={}
    class Captured(Exception):pass
    for label in labels:
        p=f.unpickle(base/'preparation'/label/'finite-prepared.pickle');ends=f.unpickle(base/'boundary-inputs/material-cases'/label/'case-material-boundary.pickle')
        direct=native.unsplit_arrays(f.unpickle(base/'interiors/cases'/label/'direct-unsplit.pickle'));rows=f.unpickle(base/'interiors/cases'/label/'row-matrices.pickle')
        old=f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-solution.pickle');held={}
        def capture(path,value):held.update(path=path,value=value);raise Captured()
        probe=types.FunctionType(finite.__code__,dict(finite.__globals__,atomic_pickle=capture),finite.__name__)
        try:probe(base/'cases'/label/'finite',p['matrix'],p['rhs'],p['unreplacedOperator'],ends['finite'],129,basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],old['polynomialDerivativeResiduals'],[rows])
        except Captured:pass
        f.require(held['path']==base/'cases'/label/'finite/finite-system.pickle','actual solver first write before solve')
        for name in ('matrix','rhs','unreplacedOperator'):f.require(m.same(held['value'][name],p[name]),'actual saved numerical first-write arrays')
        f.require(m.same(held['value']['channels'],ends['finite']),'actual material/common-Eulerian boundary before solve')
        changed=p['rhs'].copy();changed[0,0]+=1;f.require(not m.same(changed,held['value']['rhs']),'actual forcing mutation rejects')
        for end in ('LEFT','RIGHT'):
            e=ends['finite'][end];columns=[v['MATRIX_COLUMN'] for v in e['incoming']]+[v['MATRIX_COLUMN'] for v in e['outgoing'] if v['kind']=='open']
            f.require(sorted(columns)==list(range(4)) and e['current'].shape==(4,4),'actual finite open-current coverage')
            f.require(len([v for v in e['outgoing'] if v['kind']=='evanescent'])==3,'actual separate closed matching amplitudes')
        co=f.unpickle(base/'interiors/cases'/label/'interior-matrices.pickle')
        for gen in co['generators'][1:]:f.require(gen in co['gradeOrigin'],'actual independent physical contrast origin')
        if label==BASELINE:
            original=f.unpickle(base/'accepted-coordinate/coordinate-response.pickle');view=f.unpickle(base/'cases'/label/'continuum/continuum-response.pickle')
            f.require(m.same(original,view),'full historical material continuum payload unchanged')
            for key,name in (('solve','coefficient-solutions'),('response','channel-response')):
                f.require(m.same(view[key],f.unpickle(base/'cases'/label/'continuum'/(name+'.pickle'))),'actual original baseline continuum subpacket')
            del original,view
        routes[label]={'nativeFirstWrite':True,'changedForcingRejected':True,'finiteOpenColumns':4,'separateClosedDirections':3,
            'newFiniteSolve':True,'newContinuumSolve':label!=BASELINE,'unitSupplementMandatory':True}
        del p,ends,direct,rows,old,held,probe,co;gc.collect()
    answer={'nativeJoins':joins,'cases':routes,'newSolves':0,'newContractions':0}
    f.save(base/'solver-wiring.json',answer);return answer


def prohibit_completed_work():
    def forbidden(*a,**kw):raise RuntimeError('completed scientific construction prohibited in material response solve')
    h=inputs.boundary_inputs
    for cls in (h.engine.NumericalReducedAction,h.engine.ModalCurrentSubspaces,h.engine.TwoEndedMatchingChannels):cls.__init__=forbidden
    h.engine.NumericalReducedAction.bind=forbidden
    for cls in ('ClosedAcousticEnergy','ClosedCurrentPairing','AdjointCurrentMap'):setattr(h.engine,cls,forbidden)
    for name in ('__init__','maps','image','basis'):setattr(h.c.Chart,name,forbidden)
    for module,names in ((h.c,('bind','material_ends')),(h.h,('restore','restore_chart')),
        (h.c.source,('coordinate_change',)),(h.c.MaterialMomentum,('prepare_basis',)),
        (matrices.interior,('assemble',)),(matrices.matrices,('direct_cells',)),
        (native.response,('systems',)),(f,('boundary_map','source_jets'))):
        for name in names:setattr(module,name,forbidden)
    np.polynomial.legendre.leggauss=forbidden;native.sp.lambdify=forbidden


def construct(base,label,manifest):
    target=base/'cases'/label;fd=target/'finite';fd.mkdir();cd=target/'continuum'
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle');ends=f.unpickle(base/'boundary-inputs/material-cases'/label/'case-material-boundary.pickle')
    co=f.unpickle(base/'interiors/cases'/label/'interior-matrices.pickle');direct=native.unsplit_arrays(f.unpickle(base/'interiors/cases'/label/'direct-unsplit.pickle'))
    rows=f.unpickle(base/'interiors/cases'/label/'row-matrices.pickle')
    original=f.unpickle(base/'boundary-inputs/bindings/sources/accepted-bindings'/label/'case-binding.pickle')
    saved=f.unpickle(base/'preparation'/label/'finite-prepared.pickle');old=f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-solution.pickle')
    finish,_=solver.finite_solver()
    result,_=finish(fd,saved['matrix'],saved['rhs'],saved['unreplacedOperator'],ends['finite'],129,basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],old['polynomialDerivativeResiduals'],[rows])
    system=f.unpickle(fd/'finite-system.pickle');view=native.finite_view(system,result,True);f.atomic_pickle(fd/'observable.pickle',view)
    if label==BASELINE:
        value=f.unpickle(cd/'continuum-response.pickle')
        answer={'case':label,'rows':len(rows['rows']),'terms':len(original['grades']['termJoins']),'sources':len(original['binding']['jets']),
            'unknowns':645,'incidentColumns':4,'finiteRank':result['rank'],'finiteCondition':result['balancedCondition'],
            'finiteScaledResidual':b.norm(result['scaledEquationResidual']),'finiteIndependentDifference':b.norm(view['independentDifference']),
            'finiteBoundaryResidual':b.norm(view['boundaryResiduals']),'finiteCurrentRatios':view['totalCurrentRatio'].tolist(),
            'continuumRank':value['solve']['rank'],'continuumScaledResidual':b.norm(value['solve']['scaledResidual']),
            'formalRemainders':'Historical continuum has no direct formal-point solves; none invented or repeated.'}
        f.save(target/'response-checks.json',answer)
    else:
        p=f.unpickle(cd/'coefficient-systems.pickle');ref=f.unpickle(base/'accepted-modes/reference/modal.pickle')[0]
        complete,_=solver.continuum_solver()
        answer=complete(base,label,manifest,target,cd,basis,ends,ref,co,p['matrices'],p['rhs'],system,result,view,rows,original)
    answer.update(newFiniteControl=True,newContinuumControl=label!=BASELINE,reusedHistoricalContinuum=label==BASELINE)
    f.save(target/'response-disposition.json',answer);return answer


def comparison(base,label,manifest):
    target=base/'cases'/label;prior=base/'unchanged-response/cases'/label
    material=f.unpickle(target/'continuum/continuum-response.pickle');eulerian=f.unpickle(prior/'continuum/continuum-response.pickle')
    new=f.unpickle(target/'finite/observable.pickle');old=f.unpickle(prior/'finite/observable.pickle')
    ends=f.unpickle(base/'boundary-inputs/material-cases'/label/'case-material-boundary.pickle')
    physical=f.unpickle(base/'boundary-inputs/boundary-cases'/label/'case-boundary.pickle')
    co=f.unpickle(base/'interiors/cases'/label/'interior-matrices.pickle');basis=f.unpickle(base/'interiors/accepted-finite-system.pickle')
    eta,sigma=(float(co['gradeOrigin'][v]) for v in co['generators'][1:])
    for key in ('fieldUnits','rowUnits','currentUnit','settings'):f.require(m.same(material[key],eulerian[key]),'actual own-case response units/settings')
    f.require(m.same(new['positions'],old['positions']),'full same common-Eulerian field grid')
    # Coordinate identity is witnessed by actual accepted map constructions,
    # physical column addresses and raw common-basis/phase residuals, not S.
    coordinate={}
    for end in ('LEFT','RIGHT'):
        owner=ends['materialFamilyRoutes'][end];folder=base/'boundary-inputs/material-families'/owner
        sig=f.unpickle(folder/'input.pickle');f.require(m.same(sig['finite'],physical['finite'][end]) and m.same(sig['continuum'],physical['ends'][end]),'actual own physical family operands')
        finite=f.unpickle(folder/'finite/finite-material-boundary.pickle')
        f.require(m.same(finite['commonEulerian'],ends['finite'][end]),'material/common route before replacement')
        addresses={}
        for direction in ('incoming','outgoing'):
            left,right=ends['finite'][end][direction],physical['finite'][end][direction]
            fields=('RECORD_INDEX','BASIS_COLUMN','kind','k')
            f.require(len(left)==len(right),'same complete coordinate count')
            addresses[direction]=[(tuple(v[n] for n in fields),tuple(w[n] for n in fields)) for v,w in zip(left,right)]
            for v,w in zip(left,right):
                for n in fields[:-1]:f.require(m.same(v[n],w[n]),'same actual physical channel address')
        momentum=np.array([v['k']-w['k'] for direction in ('incoming','outgoing') for v,w in zip(ends['finite'][end][direction],physical['finite'][end][direction])])
        coordinate[end]={'addresses':addresses,'normalMomentum':momentum,'rightBasis':ends['finite'][end]['right']-physical['finite'][end]['right'],
            'incomingValues':ends['finite'][end]['incomingValues']-physical['finite'][end]['incomingValues'],
            'traceMap':ends['finite'][end]['traceMap']-physical['finite'][end]['traceMap'],
            'phaseIncoming':new['incomingPhase']-old['incomingPhase'],'phaseOutgoing':new['outgoingPhase']-old['outgoingPhase'],
            'finiteCurrent':ends['finite'][end]['current']-physical['finite'][end]['current']}
        f.require(m.same(material['response']['labels'][end],eulerian['response']['labels'][end]),'full continuum record and basis addresses')
    phasekeys=('incomingOriginPhase','outgoingOriginPhase','fieldIncomingMap','fieldOpenOutgoingMap')
    coordinate['continuum']={k:(b.subtract(material['response'][k],eulerian['response'][k]) if isinstance(material['response'][k],dict) else material['response'][k]-eulerian['response'][k]) for k in phasekeys}
    f.atomic_pickle(target/'common-coordinate-joins.pickle',coordinate)
    coordinate_residual=max(b.norm({k:v for k,v in coordinate[e].items() if k!='addresses'}) for e in ('LEFT','RIGHT'))
    f.require(coordinate_residual<1e-10 and b.norm(coordinate['continuum'])<1e-10,'computed common-coordinate basis/phase/current joins before comparison')
    grid=f.polynomial_basis(new['positions'],basis['settings']['sourceBound'],len(basis['nodes']))
    def fields(value):
        maps=b.J.multiply(value['response']['incomingOriginPhase'],value['response']['incomingCurrentRootInverse'])
        parts=[b.J.multiply({g:grid@a[i*129:(i+1)*129] for g,a in value['solve']['coefficients'].items()},maps) for i in range(5)]
        return {g:np.stack([p[g] for p in parts]) for g in b.G}
    newfields,oldfields=fields(material),fields(eulerian)
    keys=('openBoundaryScattering','openOriginScattering','fieldOriginScattering','fluxOriginScattering','incomingCurrentOrigin','outgoingCurrentOrigin')
    channels={k:b.subtract(material['response'][k],eulerian['response'][k]) for k in keys}
    def quotient(value):
        a=b.evaluate(value['response']['openOriginScattering'],eta,sigma);ji=b.evaluate(value['response']['incomingCurrentOrigin'],eta,sigma);jo=b.evaluate(value['response']['outgoingCurrentOrigin'],eta,sigma)
        q=np.diag(a.conj().T@jo@a)/np.diag(ji)
        independent=np.array([sum(a[i,j].conjugate()*jo[i,k]*a[k,j] for i in range(4) for k in range(4))/ji[j,j] for j in range(4)])
        f.require(b.norm(q-independent)<1e-12,'independent physical retained-polynomial current quotient')
        return {'amplitude':a,'incomingCurrent':ji,'outgoingCurrent':jo,'ratio':q,'independent':independent,'residual':q-independent}
    quotients={'material':quotient(material),'eulerian':quotient(eulerian)}
    quotients['difference']=quotients['material']['ratio']-quotients['eulerian']['ratio']
    differences={'finiteScattering':new['originScattering']-old['originScattering'],'finiteCurrent':new['totalCurrentRatio']-old['totalCurrentRatio'],
        'finiteFields':new['originFields']-old['originFields'],'continuumChannels':channels,
        'evaluatedChannels':{k:b.evaluate(v,eta,sigma) for k,v in channels.items()},'continuumFields':{'material':newfields,'eulerian':oldfields,'difference':b.subtract(newfields,oldfields)},
        'continuumFieldCoefficients':b.subtract(material['solve']['coefficients'],eulerian['solve']['coefficients']),
        'closedMatching':{e:b.subtract(material['response']['outgoingFieldAmplitudes'][e],eulerian['response']['outgoingFieldAmplitudes'][e]) for e in ('LEFT','RIGHT')},
        'retainedPolynomialCurrents':quotients,'gradeOrigin':co['gradeOrigin'],'generators':co['generators'],'fieldPositions':new['positions'],
        'fieldUnits':material['fieldUnits'],'currentUnit':material['currentUnit'],'scope':SCOPE}
    f.atomic_pickle(target/'coordinate-comparisons.pickle',differences)
    # New material-response amplitude bookkeeping uses its actual common maps.
    reference=f.unpickle(base/'accepted-modes/reference/modal.pickle')[0]
    selection=current.selectors(material,ends,reference)
    previous=f.unpickle(base/'unchanged-currents/cases'/label/'continuum/continuum-currents.pickle')
    f.require(m.same(selection,previous['selection']),'full actual selector and candidate dispositions')
    metrics,mres=current.open_metrics(material,ends,selection)
    f.atomic_pickle(target/'continuum/current-metric-inputs.pickle',{'selection':selection,'metrics':metrics,'residuals':mres})
    opened,mutation=current.construct_open(material,metrics,selection)
    f.atomic_pickle(target/'continuum/open-channel-currents.pickle',{'metrics':metrics,'records':opened,'residuals':mres,'mutation':mutation})
    closed={}
    for end in ('LEFT','RIGHT'):
        ix=[i for row in selection['census'][end] if row['direction']=='outgoing' and row['kind']=='evanescent' for i in row['columns']]
        closed[end]=current.amplitude_bookkeeping({g:a[ix] for g,a in material['response']['outgoingFieldAmplitudes'][end].items()})
    end_current=current.end_currents(material,ends);f.atomic_pickle(target/'continuum/finite-end-currents.pickle',end_current)
    domains=f.unpickle(base/'unchanged-currents/cases'/label/'continuum/bulk-kinematics.pickle')
    f.require(m.same(domains,previous['bulkDomains']),'exact accepted own-case physical depth domains in common Eulerian coordinates')
    bookkeeping={'selection':selection,'metrics':metrics,'metricResiduals':mres,'open':opened,'closedAmplitudes':closed,'endCurrents':end_current,
        'bulkDomains':domains,'currentMutation':mutation,'currentUnit':material['currentUnit'],'ratio':material['ratio'],
        'settings':material['settings'],'dimensionState':material['dimensionState'],'materialUnitSupplement':manifest['materialUnitSupplement'],
        'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':SCOPE+' Depth-integrated normal current is not depth escape.'}
    f.atomic_pickle(target/'continuum/continuum-currents.pickle',bookkeeping)
    fin=solver.finite_open_current(f.unpickle(target/'finite/finite-system.pickle'),f.unpickle(target/'finite/finite-solution.pickle'))
    old_fin=f.unpickle(base/'unchanged-finite-currents'/(label+'.pickle'))['baseline']
    fdiff={e:{k:fin[e][k]-old_fin[e][k] for k in ('full','outgoing','incoming','interference')} for e in ('LEFT','RIGHT')}
    f.atomic_pickle(target/'finite/open-end-currents.pickle',{'material':fin,'eulerian':old_fin,'differences':fdiff})
    current_diff={'open':{name:{part:current.subtract(opened[name]['currentParts'][part],previous['open'][name]['currentParts'][part]) for part in ('slab','bulk','total')} for name in opened},
        'ends':{end:{part:{kind:current.subtract(end_current[end]['parts'][part][kind],previous['endCurrents'][end]['parts'][part][kind]) for kind in ('full','outgoing','incoming','interference')} for part in end_current[end]['parts']} for end in ('LEFT','RIGHT')}}
    f.atomic_pickle(target/'coordinate-current-comparisons.pickle',current_diff)
    summary={'commonCoordinateResidual':coordinate_residual,'finiteAmplitudeDifference':b.norm(differences['finiteScattering']),
        'finiteCurrentDifference':b.norm(differences['finiteCurrent']),'finiteFieldDifferenceDescriptive':b.norm(differences['finiteFields']),
        'continuumAmplitudeDifference':b.norm(differences['evaluatedChannels']['fluxOriginScattering']),
        'retainedPolynomialCurrentDifference':b.norm(quotients['difference']),
        'openCurrentResidual':max(b.norm(v['residuals']) for v in opened.values()),'endCurrentResidual':max(b.norm(v['residuals']) for v in end_current.values()),
        'finiteOpenCurrentResidual':max(b.norm(v['residuals']) for v in fin.values()),'currentMutation':b.norm(mutation),
        'amplitudeGoal':1e-4,'currentGoal':1e-6,'scope':SCOPE}
    f.save(target/'coordinate-checks.json',summary);return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--resume-focused',type=Path,required=True)
    p.add_argument('--preflight',action='store_true');p.add_argument('--resume-preflight',type=Path);args=p.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    manifest,labels=load(base,args.resume_focused,args.resume_preflight)
    if args.preflight:
        def forbidden(*a,**kw):raise RuntimeError('numerical solve prohibited in material solver-routing preflight')
        for module in (np.linalg,native.la):
            for name in ('solve','inv','pinv','lstsq','svd','eig','eigh','eigvals','eigvalsh','lu_factor','lu_solve','solve_sylvester'):
                if hasattr(module,name):setattr(module,name,forbidden)
    wire=json.loads((base/'solver-wiring.json').read_text()) if args.resume_preflight else wiring(base,manifest,labels)
    prohibit_completed_work();solves={};comparisons={}
    if not args.preflight:
        for label in labels:
            solves[label]=construct(base,label,manifest);f.save(base/'response-inventory.json',solves);gc.collect()
        for label in labels:
            comparisons[label]=comparison(base,label,manifest);f.save(base/'comparison-inventory.json',comparisons);gc.collect()
        f.atomic_pickle(base/'remaining-case-material-response.pickle',{'cases':solves,'comparisons':comparisons,'sourceFiles':manifest['sourceFiles'],
            'inputPackets':manifest['inputPackets'],'materialUnitSupplement':manifest['materialUnitSupplement'],'scope':SCOPE})
    matrices.hash_check(base,manifest)
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    result={**manifest,'status':'PASSED_MATERIAL_RESPONSE_SOLVER_WIRING' if args.preflight else 'COMPLETED_CASE_MATERIAL_RESPONSES',
        'solverWiring':wire,'cases':solves,'comparisons':comparisons,'newFiniteControlSolves':len(solves),
        'newContinuumControlSolves':sum(v['newContinuumControl'] for v in solves.values()),'reusedHistoricalContinuum':True,
        'newMaterialAmplitudeBookkeeping':len(comparisons),'newQuadratureNodes':0,'newInteriorAssemblies':0,'newModes':0,'newCurrentClosures':0,
        'physicalTranscriptEmitted':False,'artifacts':artifacts,'wallSeconds':time.monotonic()-start,'maximumRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
