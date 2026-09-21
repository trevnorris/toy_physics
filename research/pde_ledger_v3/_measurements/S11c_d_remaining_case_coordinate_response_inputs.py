#!/usr/bin/env python3
"""Prepare actual material-coordinate response systems from accepted packets."""
import argparse
import copy
import gc
import json
from pathlib import Path
import resource
import shutil
import signal
import time
import types

import numpy as np
import S11c_d_remaining_case_coordinate_boundary as boundary_inputs
import S11c_d_remaining_case_response as response
import S11c_d_remaining_case_first_jet_response as prepared_solver

f, m, b = boundary_inputs.f, boundary_inputs.modes, response.b
matrices = boundary_inputs.matrices
BASELINE = matrices.BASELINE
ICP = f.M/'S11c_d_remaining_case_coordinate_matrices_checkpoint.json'
BCP = f.M/'S11c_d_remaining_case_coordinate_boundary_checkpoint.json'
RCP = f.M/'S11c_d_remaining_case_response_checkpoint.json'
CCP = f.M/'S11c_d_coordinate_response_checkpoint.json'
PLAN = f.M/'S11c_d_remaining_case_coordinate_response_plan.md'
WIRING = f.M/'S11c_d_remaining_case_coordinate_response_wiring.json'
SCOPE = ('Saved material operators and actual common-Eulerian finite/continuum end/current/phase routes before solving. '
         'Four finite controls and three continuum controls remain new work; historical baseline continuum is reused. '
         'No quadrature, source binding, interior assembly, boundary/current construction or response solve in this focus.')


def load(base):
    accepted = [matrices.accepted(path, status) for path, status in (
        (ICP,'ACCEPTED_CASE_MATERIAL_INTERIOR_MATRICES'),
        (BCP,'ACCEPTED_CASE_MATERIAL_BOUNDARY_MAPS'),
        (RCP,'PUBLISHED_ANNEX_VERIFIED'), (CCP,'PUBLISHED_ANNEX_VERIFIED'))]
    (ir,ic,ip),(br,bc,bp),(rr,rc,rp),(cr,cc,cp) = accepted
    manifest = {'runDirectory':str(base),'sourceFiles':{},'inputPackets':{},'copiedInputs':{},
                'input':ic['input'],'settings':ic['settings'],'scope':SCOPE}
    for root, checks, checkpoint in accepted:
        for n,sha in checks['sourceFiles'].items():
            f.require(n not in manifest['sourceFiles'] or manifest['sourceFiles'][n]==sha,'unchanged shared source')
            f.require(f.digest(f.ROOT/n)==f.digest(root/'source'/n)==sha,'actual accepted current/frozen source')
            manifest['sourceFiles'][n]=sha
        for n,sha in checks['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n]==sha,'identical inherited input')
            manifest['inputPackets'][n]=sha
        manifest['inputPackets'][str(root/'checks.json')]=checkpoint['checksSha256']
    labels=tuple(bc['result']['cases'])
    def retain(root, checks, name, destination):
        f.require(name in checks['artifacts'],('actual accepted artifact',name))
        m.retain(root/name,base/destination,manifest,checks['artifacts'][name]['sha256'])
    # Keep all raw family/current/chart/source proof packets and failed-run lineage.
    for name in bc['artifacts']: retain(br,bc,name,'boundary-inputs/'+name)
    m.retain(br/'checks.json',base/'accepted-boundary-checks.json',manifest,bp['checksSha256'])
    m.retain(ir/'checks.json',base/'accepted-interior-checks.json',manifest,ip['checksSha256'])
    supplement=bp['metadataSupplement'];vr=Path(bp['validator']['runDirectory'])
    f.require(supplement['mandatoryForDownstreamUnitConsumers'] and Path(supplement['runDirectory'])==vr,
              'mandatory accepted material unit supplement')
    f.require(f.digest(vr/'checks.json')==bp['validator']['checksSha256'],'final saved boundary validation')
    m.retain(vr/'checks.json',base/'accepted-boundary-validation.json',manifest,bp['validator']['checksSha256'])
    for n,v in supplement['artifacts'].items():
        m.retain(vr/n,base/'boundary-unit-supplement'/n,manifest,v['sha256'])
    for n,sha in bp['validator']['sourceFiles'].items():
        f.require(n not in manifest['sourceFiles'] or manifest['sourceFiles'][n]==sha,'unit supplement source join')
        f.require(f.digest(f.ROOT/n)==sha,'unchanged supplemental unit helper')
        manifest['sourceFiles'][n]=sha
    manifest['materialUnitSupplement']={'checkpoint':str(BCP),'validatorChecksSha256':bp['validator']['checksSha256'],
        'runDirectory':str(vr),'artifacts':supplement['artifacts'],'mandatory':True}
    retain(ir,ic,'eulerian/accepted-finite-system.pickle','interiors/accepted-finite-system.pickle')
    retain(ir,ic,'preflight.json','accepted-matrix-preflight.json')
    retain(ir,ic,'transported-trial-reuse.pickle','transported-trial-reuse.pickle')
    retain(ir,ic,'eulerian/cases/'+BASELINE+'/interior-matrices.pickle','baseline/eulerian-interior.pickle')
    retain(ir,ic,'bindings/baseline-material-row-matrices.pickle','baseline/material-row-matrices.pickle')
    for n in ('material-binding.pickle','material-coefficient-matrices.pickle','material-coefficient-systems.pickle',
              'material-coefficient-solutions.pickle','material-channel-response.pickle','coordinate-response.pickle',
              'material-boundary.pickle','coordinate-comparisons.pickle','full.out'):
        retain(cr,cc,n,'accepted-coordinate/'+n)
    m.retain(cr/'checks.json',base/'accepted-coordinate/checks.json',manifest,cp['checksSha256'])
    for label in labels:
        for part,names in (('finite',('finite-system','finite-solution','observable')),
                           ('continuum',('coefficient-systems','coefficient-solutions','channel-response','continuum-response'))):
            for name in names:
                address='cases/'+label+'/'+part+'/'+name+'.pickle'
                retain(rr,rc,address,'unchanged-response/'+address)
        if label!=BASELINE:
            for n in ('interior-matrices','direct-unsplit','direct-native-cells','row-matrices','comparisons'):
                address='cases/'+label+'/'+n+'.pickle';retain(ir,ic,address,'interiors/'+address)
            retain(ir,ic,'bindings/material-cases/'+label+'/case-material-binding.pickle','material-bindings/'+label+'.pickle')
        for n in ('integration-plan','coefficient-routes','native-cell-routes'):
            retain(ir,ic,'preparation/'+label+'/'+n+'.pickle','accepted-preparation/'+label+'/'+n+'.pickle')
        retain(ir,ic,'eulerian/cases/'+label+'/direct-unsplit.pickle','eulerian/'+label+'/direct-unsplit.pickle')
    retain(rr,rc,'accepted-modes/reference/modal.pickle','accepted-modes/reference/modal.pickle')
    for path in (Path(__file__),PLAN,WIRING,ICP,BCP,RCP,CCP,Path(response.__file__),Path(prepared_solver.__file__),f.ACCEPTANCE):
        manifest['sourceFiles'][str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    for n,sha in manifest['sourceFiles'].items():
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target)
        f.require(f.digest(target)==sha,'frozen response preparation source')
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle')
    f.require(json.loads(json.dumps(basis['settings']))==manifest['settings'],'exact native settings and serialized view')
    manifest['settings']=basis['settings']
    manifest['acceptedMaterialBoundaryArtifacts']=len(bc['artifacts'])
    f.save(base/'inputs.json',manifest);matrices.hash_check(base,manifest)
    return manifest,labels


def baseline_views(base,manifest):
    """Retain every saved coefficient grade; evaluate the full polynomial once.

    This is a storage/arithmetic view, not a new native-cell assembly or finite
    solve. The baseline producer saved all source grades, not just the response
    rectangle. Compare the full recombination to the accepted direct operator.
    """
    old=f.unpickle(base/'baseline/eulerian-interior.pickle')
    keys=('size','generators','gradeOrigin','blockUnits','fieldUnits','equationUnits','dimensionState','settings')
    co={k:old[k] for k in keys};del old;gc.collect()
    packet=f.unpickle(base/'accepted-coordinate/material-coefficient-matrices.pickle')
    co['matrices']=packet['matrices'];del packet
    original=f.unpickle(base/'boundary-inputs/bindings/sources/accepted-bindings'/BASELINE/'case-binding.pickle')
    material=f.unpickle(base/'accepted-coordinate/material-binding.pickle')
    support=set(material['support']);f.require(all(set(co['matrices'][n])==support for n in ('local','nonlocal','total')),
                                                  'all original material grades, including outside the response rectangle')
    f.require(support=={g for v in original['grades']['records'].values() for g in v['record']['COEFFICIENTS']},
              'complete original source-grade support before baseline finite recombination')
    full={n:matrices.interior.recombine(co['matrices'][n],co['gradeOrigin'],co['generators']) for n in ('local','nonlocal','total')}
    saved_total=full['total'];full['total']=full['local']+full['nonlocal']
    old_direct=response.unsplit_arrays(f.unpickle(base/'eulerian'/BASELINE/'direct-unsplit.pickle'))
    differences={n:matrices.interior.differences(full[n],old_direct[n],129) for n in full}
    addition=saved_total-full['total']
    target=base/'interiors/cases'/BASELINE;target.mkdir(parents=True)
    f.atomic_pickle(target/'interior-matrices.pickle',dict(co,scope='All accepted baseline material coefficient arrays; no reassembly.'))
    f.atomic_pickle(target/'direct-unsplit.pickle',{n:{(0,0,0):v} for n,v in full.items()})
    f.atomic_pickle(target/'baseline-recombination.pickle',{'differences':differences,'additionOrderDifference':addition,
        'gradeOrigin':co['gradeOrigin'],'support':tuple(sorted(support)),
        'scope':'Full saved coefficient-polynomial evaluation. This is not a new direct native-cell contraction.'})
    rows=f.unpickle(base/'baseline/material-row-matrices.pickle')
    f.atomic_pickle(target/'row-matrices.pickle',{'rows':rows,'settings':manifest['settings'],'wholeBaselineRowsReused':True})
    f.require(max(v['maximumScaledReferenceFrame'] for v in differences.values())<1e-9 and
              b.norm(addition)/(1+b.norm(saved_total))<1e-12,'complete baseline material arithmetic view and own Eulerian direct operator')
    result=f.unpickle(base/'accepted-coordinate/coordinate-response.pickle')
    f.require(m.same(result['matrices'],co['matrices']) and m.same(result['binding'],material) and m.same(result['rows'],rows),
              'actual original full material response operands')
    f.require(m.same(result['generators'],co['generators']) and m.same(result['settings'],co['settings']),
              'actual baseline coordinate grades and settings')
    f.save(base/'baseline-storage-views.json',{'newAssemblies':0,'newNumericalRows':0,'finiteBaselineSolveExists':False,
        'completeSourceGrades':[list(g) for g in sorted(support)],'arithmeticOnly':True,
        'operatorMaximum':max(v['maximumScaledReferenceFrame'] for v in differences.values()),
        'parents':{str(p.relative_to(base)):f.digest(p) for p in (base/'accepted-coordinate/material-coefficient-matrices.pickle',
        base/'baseline/eulerian-interior.pickle',base/'accepted-coordinate/coordinate-response.pickle',base/'baseline/material-row-matrices.pickle')}})
    del result,material,original,co,full,old_direct,differences,rows;gc.collect()


def disable_constructors():
    # Keep the captured pure boundary-replacement systems function available;
    # disable its module entry point along with all scientific constructors.
    boundary_inputs.inputs.disable_constructors()
    def forbidden(*a,**kw): raise RuntimeError('scientific construction prohibited in saved coordinate-response input focus')
    for module in (np.linalg,response.la):
        for name in ('solve','inv','pinv','lstsq','svd','eig','eigh','eigvals','eigvalsh','lu_factor','lu_solve','solve_sylvester'):
            if hasattr(module,name):setattr(module,name,forbidden)
    np.polynomial.legendre.leggauss=forbidden;response.sp.lambdify=forbidden
    for cls in ('UniformSlabCurrent','ClosedAcousticEnergy','ClosedCurrentPairing','ModalCurrentSubspaces','AdjointCurrentMap'):
        setattr(f.engine,cls,forbidden)


def prepare(base,manifest,labels,systems,finite_prefix,finite_solver,wiring):
    baseline_views(base,manifest)
    basis=f.unpickle(base/'interiors/accepted-finite-system.pickle')
    state=f.unpickle(base/'boundary-inputs/bindings/chart-state.pickle');g=state['values']
    f.require(len(basis['nodes'])==129 and basis['nodes'][0]==-64 and basis['nodes'][-1]==64,'actual complete approved trial interval')
    trial=f.unpickle(base/'transported-trial-reuse.pickle')
    for n in ('nodes','derivativeMatrices','settings'):f.require(m.same(trial[n],basis[n]),'accepted actual transported trial array')
    f.require(m.same(trial['chartState'],state) and max(trial['acceptedLocalWitness']['derivativeResiduals'].values())<1e-9,
              'actual complete chart and saved native transported derivative witness')
    f.save(base/'trial-input-joins.json',{'fullSavedTrialAndChartIdentity':True,'newBasisConstruction':False,
        'inheritedDerivativeMaximum':max(trial['acceptedLocalWitness']['derivativeResiduals'].values()),
        'evidenceSha256':f.digest(base/'transported-trial-reuse.pickle')})
    del trial;gc.collect()
    inventory={};row_count=term_count=source_count=0
    for label in labels:
        target=base/'preparation'/label;target.mkdir(parents=True)
        original=f.unpickle(base/'boundary-inputs/bindings/sources/accepted-bindings'/label/'case-binding.pickle')
        if label==BASELINE: material=f.unpickle(base/'accepted-coordinate/material-binding.pickle')
        else:
            bound=f.unpickle(base/'material-bindings'/(label+'.pickle'))
            f.require(m.same(bound['originalCase'],original),'exact own-case physical source and material binding')
            material=bound['material'];del bound
        co=f.unpickle(base/'interiors/cases'/label/'interior-matrices.pickle')
        plan=f.unpickle(base/'accepted-preparation'/label/'integration-plan.pickle')
        for n in ('gradeOrigin','fieldUnits','equationUnits','settings'):f.require(m.same(co[n],plan[n]),'actual accepted typed coefficient preparation')
        f.require(m.same(plan['nodes'],basis['nodes']),'full own-case prepared node identity')
        if label!=BASELINE:f.require(np.array_equal(co['materialCoordinates'],basis['nodes']/float(g['d'])),
                                    'material coefficients were evaluated at the actual transported coordinates')
        direct=response.unsplit_arrays(f.unpickle(base/'interiors/cases'/label/'direct-unsplit.pickle'))
        rows=f.unpickle(base/'interiors/cases'/label/'row-matrices.pickle')
        ends=f.unpickle(base/'boundary-inputs/material-cases'/label/'case-material-boundary.pickle')
        physical=f.unpickle(base/'boundary-inputs/boundary-cases'/label/'case-boundary.pickle')
        old_system=f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-system.pickle')
        old_solution=f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-solution.pickle')
        f.require(m.same(old_system['channels'],physical['finite']),'actual unchanged response owns these original finite coordinates')
        for n in ('nodes','derivativeMatrices','settings'): f.require(m.same(old_system[n],basis[n]),'complete own-case common trial basis/settings')
        f.require(material['settings']==co['settings']==rows['settings']==basis['settings'],'all full native material settings')
        for n,cn in (('fieldUnits','fieldUnits'),('rowUnits','equationUnits')):
            f.require(m.same(ends[n],physical[n]) and m.same(ends[n],co[cn]) and m.same(ends[n],original['binding'][cn]),'actual field/equation unit route')
        f.require(m.same(ends['currentUnit'],physical['currentUnit']) and m.same(material['geometry'],g['g']), 'actual chart/current units')
        f.require(set(rows['rows'])==set(range(len(material['rows']))) and len(material['rows'])==len(original['binding']['bound']['rows']), 'complete material row addresses')
        row_routes=[]
        for row,old in zip(material['rows'],original['binding']['bound']['rows']):
            for n in ('index','original','unit'):f.require(m.same(row[n],old[n]),'physical source/address/row unit')
            mapped=tuple((v,g['d']*a+g['kappa'],g['d']*z+g['kappa']) for v,a,z in old['limits'])
            source=(old['sourceLimit'][0],old['sourceLimit'][1]/g['d'],old['sourceLimit'][2]/g['d'])
            f.require(row['limits']==mapped and row['sourceLimit']==source,'ordered material momentum and source endpoints')
            row_routes.append({'row':row['index'],'originalLimits':old['limits'],'materialLimits':row['limits'],
                'originalSourceLimit':old['sourceLimit'],'materialSourceLimit':row['sourceLimit'],'unit':row['unit']})
        for term in original['grades']['termJoins']:
            f.require(m.same(term['originalIntegral'],material['rows'][term['integralIndex']]['original']),'all original physical cell-term addresses')
        f.atomic_pickle(target/'source-row-joins.pickle',row_routes)
        boundary_routes={}
        for end in ('LEFT','RIGHT'):
            owner=ends['materialFamilyRoutes'][end];folder=base/'boundary-inputs/material-families'/owner
            signature=f.unpickle(folder/'input.pickle')
            f.require(m.same(signature['chart'],state) and m.same(signature['finite'],physical['finite'][end]) and
                      m.same(signature['continuum'],physical['ends'][end]),'full own-case original source at the accepted material family')
            finite=f.unpickle(folder/'finite/finite-material-boundary.pickle')
            continuum=f.unpickle(folder/'continuum'/(end.lower()+'-material-boundary.pickle'))
            f.require(m.same(ends['finite'][end],finite['commonEulerian']) and m.same(ends['ends'][end],continuum['commonEulerian']) and
                      m.same(ends['finitePhases'][end],finite['phases']), 'actual material-to-common route before boundary replacement')
            f.require(ends['finite'][end]['current'].shape==(4,4) and len(finite['closedMatchingDirections'])==3 and
                      int(ends['ends'][end]['offsets'][-1])==7,'actual open/closed current domains')
            boundary_routes[end]={'owner':owner,'finite':f.digest(folder/'finite/finite-material-boundary.pickle'),
                'continuum':f.digest(folder/'continuum'/(end.lower()+'-material-boundary.pickle')),
                'unitSupplementRequired':owner=='LAB_HELD__RHOBR_CONSTANT__RIGHT'}
        if any(v['unitSupplementRequired'] for v in boundary_routes.values()):
            units=f.unpickle(base/'boundary-unit-supplement/unit-overlay/unit-source-pairs.pickle')
            view=f.unpickle(base/'boundary-unit-supplement/unit-overlay/dimension-state-view.pickle')
            for atom,unit in units['joinedKnown'].items():f.require(tuple(view['known'][atom])==tuple(unit)==(-1,0,0),'actual accepted supplemental current units')
        f.atomic_pickle(target/'boundary-input-joins.pickle',boundary_routes)
        a,rhs,retained=finite_prefix(target,direct['total'].copy(),ends['finite'],129,basis['derivativeMatrices'],
            basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],old_solution['polynomialDerivativeResiduals'],[rows])
        f.atomic_pickle(target/'finite-prepared.pickle',{'matrix':a,'rhs':rhs,'unreplacedOperator':retained})
        ca,crhs=systems(co,ends,basis)
        f.atomic_pickle(target/'continuum-prepared.pickle',{'matrices':ca,'rhs':crhs})
        f.require(a.shape==(645,645) and rhs.shape==(645,4) and np.array_equal(retained,direct['total']) and np.isfinite(a).all() and np.isfinite(rhs).all(),'full finite system before solve')
        f.require(set(ca)==set(crhs)==set(b.G) and all(v.shape==(645,645) and np.isfinite(v).all() for v in ca.values()) and
                  all(v.shape==(645,4) and np.isfinite(v).all() for v in crhs.values()),'full independent-grade systems and forcing')
        # Actual one-sided sign and omitted material covector reach the native RHS.
        wrong=copy.deepcopy(ends['finite']);wrong['LEFT']['incomingBoundaryData'] *= -1
        _,wrong_rhs,_=finite_prefix(target,direct['total'].copy(),wrong,129,basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],{},{},old_solution['polynomialDerivativeResiduals'],[rows])
        left=base/'boundary-inputs/material-families'/ends['materialFamilyRoutes']['LEFT']
        lf=f.unpickle(left/'finite/finite-material-boundary.pickle');X,U,T=f.unpickle(left/'saved-chart-map-reuse.pickle')['maps']
        omitted=T@lf['material']['incomingDerivative']/float(g['d'])-lf['commonEulerian']['traceMap']@lf['commonEulerian']['incomingValues']
        wrong_covector=copy.deepcopy(ends['finite']);wrong_covector['LEFT']['incomingBoundaryData']=omitted
        _,cov_rhs,_=finite_prefix(target,direct['total'].copy(),wrong_covector,129,basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],{},{},old_solution['polynomialDerivativeResiduals'],[rows])
        controls={'incidentSign':wrong_rhs-rhs,'omittedIncomingCovector':cov_rhs-rhs,'omittedForcing':omitted,
                  'correctForcing':ends['finite']['LEFT']['incomingBoundaryData']}
        f.atomic_pickle(target/'forcing-controls.pickle',controls)
        f.require(b.norm(controls['incidentSign'])>1e-10 and b.norm(controls['omittedIncomingCovector'])>1e-10,'actual forcing sign and material covector omissions respond')
        # Stop the unchanged finite solver at its first durable write.
        class Captured(Exception):pass
        held={}
        def capture(path,payload):held.update(path=path,payload=payload);raise Captured()
        fn=types.FunctionType(finite_solver.__code__,dict(finite_solver.__globals__,atomic_pickle=capture),finite_solver.__name__)
        try:fn(target/'future-finite',a,rhs,retained,ends['finite'],129,basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],manifest['sourceFiles'],manifest['inputPackets'],old_solution['polynomialDerivativeResiduals'],[rows])
        except Captured:pass
        f.require(held['path']==target/'future-finite/finite-system.pickle','unchanged native first write before solve')
        for n,val in (('matrix',a),('rhs',rhs),('unreplacedOperator',retained),('channels',ends['finite'])):f.require(m.same(held['payload'][n],val),'actual first-write array and own boundary input')
        changed=rhs.copy();changed[0,0]+=1;f.require(not m.same(changed,held['payload']['rhs']),'changed native forcing rejects')
        f.save(target/'first-write-checks.json',{'nativeTail':wiring['finiteSolver'],'capturedBeforeSolve':True,'changedForcingRejected':True})
        baseline_replay=None
        if label==BASELINE:
            old=f.unpickle(base/'accepted-coordinate/material-coefficient-systems.pickle')
            solved=f.unpickle(base/'accepted-coordinate/material-coefficient-solutions.pickle')
            delta=b.subtract(ca,old['matrices']);fd=b.subtract(crhs,old['rhs'])
            residual=b.subtract(b.J.multiply(ca,solved['coefficients']),crhs)
            scaled={grade:v/solved['rowScale'][:,None] for grade,v in residual.items()}
            f.atomic_pickle(target/'baseline-continuum-replay.pickle',{'matrixDifferences':delta,'forcingDifferences':fd,'savedSolutionScaledResidual':scaled})
            f.require(b.norm(delta)==b.norm(fd)==0 and b.norm(scaled)<1e-9,'exact original baseline material equations and saved solution reuse')
            baseline_replay={'matrixMaximum':b.norm(delta),'forcingMaximum':b.norm(fd),'scaledResidual':b.norm(scaled)}
            del old,solved,delta,fd,residual,scaled
        row_count+=len(material['rows']);term_count+=len(original['grades']['termJoins']);source_count+=len(material['jets'])
        inventory[label]={'rows':len(material['rows']),'terms':len(original['grades']['termJoins']),'sources':len(material['jets']),
            'unknowns':645,'incidentColumns':4,'nativeFirstWriteCaptured':True,'newFiniteSolveRequired':True,
            'newContinuumSolveRequired':label!=BASELINE,'baselineContinuumReplay':baseline_replay,
            'incidentSignMutation':b.norm(controls['incidentSign']),'covectorOmission':b.norm(controls['omittedIncomingCovector']),
            'materialFamilies':ends['materialFamilyRoutes'],'fullCurrentAndPhaseInputJoins':True}
        f.save(base/'case-inventory.json',inventory)
        del original,material,co,plan,direct,rows,ends,physical,old_system,old_solution,a,rhs,retained,ca,crhs,held,fn,wrong,wrong_rhs,wrong_covector,cov_rhs,lf,controls;gc.collect()
    f.require((row_count,term_count,source_count)==(300,647,120),'full original four-case source census')
    result={'cases':inventory,'rowAddresses':row_count,'nativeTerms':term_count,'sourceAmplitudes':source_count,
            'remainingFiniteControlSolves':4,'remainingContinuumControlSolves':3,'baselineContinuumReused':True,
            'newQuadratureNodes':0,'newInteriorAssemblies':0,'newBoundaryConstructions':0,'newCurrentContractions':0,'newSolves':0,
            'nativeWiring':wiring,'scope':SCOPE}
    f.save(base/'preflight.json',result);return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    _,finite_prefix,prefix_join=response.finite_tail();finite,solver_join=prepared_solver.finite_solver()
    systems=response.response.systems
    wiring={'finitePrefix':prefix_join,'finiteSolver':solver_join,
            'nativeContinuum':{n:matrices.body(getattr(response.response,n)) for n in ('systems','solve','channels','open_flux')},
            'recombine':matrices.body(matrices.interior.recombine),'newScientificSolves':0}
    f.save(base/'native-response-wiring.json',wiring)
    disable_constructors();manifest,labels=load(base)
    result=prepare(base,manifest,labels,systems,finite_prefix,finite,wiring)
    matrices.hash_check(base,manifest)
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
               if p.is_file() and 'source' not in p.relative_to(base).parts and p.name not in ('inputs.json','checks.json')}
    # Include copied nested checks/input packets; only this run's final manifests are excluded.
    for n in manifest['copiedInputs']:
        p=base/n;artifacts[n]={'sha256':f.digest(p),'bytes':p.stat().st_size}
    checks={**manifest,'status':'COMPLETED_CASE_MATERIAL_RESPONSE_INPUTS','preflight':result,'cases':result['cases'],
            'artifacts':artifacts,'wallSeconds':time.monotonic()-start,'maximumRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
