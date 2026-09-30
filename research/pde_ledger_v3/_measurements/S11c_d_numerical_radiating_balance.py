#!/usr/bin/env python3
"""Central finite balance after inspected saved integration/refinement checks.

No completed source, endpoint, current, path or reference-integral replay.
"""
import argparse
import gc
import json
from pathlib import Path
import resource
import sys
import time
import traceback
import warnings
import S11c_d_numerical_radiating_end_maps_v2 as storage
import S11c_d_numerical_radiating_integration_v3 as I
import S11c_d_numerical_radiating_finite as F
from S11c_d_numerical_radiating_blob_store import BlobStore
from S11c_d_numerical_radiating_saved_operations import exact_structure
from S11c_d_numerical_radiating_uniform_refinement_v2 import restored_blob
require=I.require


def matrix_group_setting(setting,variables):
    result=dict(setting)
    if len(variables)==1:result['outerOrder']=setting['singleLegOuterOrder']
    return result


def run(manifest,J):
    prior=manifest['priorStartup'];root=Path(prior['path'])
    for rel,v in prior['files'].items():require(storage.digest(root/rel)==v['sha256'],'saved pilot file '+rel)
    operations=[json.loads(x) for x in (root/'operation-index.jsonl').read_text().splitlines()]
    require(len(operations)==36 and len({x['name'] for x in operations})==36,'complete saved pilot journal')
    old={x['name']:x for x in operations};source=BlobStore(root/'operations.sqlite');provenance={'database':str(root/'operations.sqlite'),'sha256':prior['files']['operations.sqlite']['sha256']}
    def restore(name):return restored_blob(J,source,old[name]['return'],name,dict(provenance,operation=old[name]))
    try:
        frequency=restore('restore/frequency');native=restore('restore/native')['bound'];endmaps=restore('restore/completed-end-maps')
        bindings={a:restore(f'bind/contrast-{a}') for a in (0.,.25,.5,1.)}
        roots=restore('integration/radical-inventory');settings,coverage=restore('discretization/coverage')
        endpoint_proof=restore('integration/endpoint-order-census');branch_joins=restore('integration/source-branch-joins')
        endpoint_values=restore('integration/two-sided-endpoint-values')
        for name in ('transformed-base','direct-part-0','direct-part-1','direct-part-2','transformed-refined','wrong-middle-sheet','omitted-transformed-jacobians'):
            restore('integration/middle/row-70/'+name)
    finally:source.close()
    for name in ('middle-integration-comparison.json','integration-mutations.json','two-sided-endpoint-values.json'):
        (J.base/('prior-'+name)).write_bytes((root/name).read_bytes())
        require(storage.digest(J.base/('prior-'+name))==prior['files'][name]['sha256'],'completed integration summary copy')
    def complex_saved(x):return complex(x['real'],x['imag'])
    middle=json.loads((root/'middle-integration-comparison.json').read_text());mutation=json.loads((root/'integration-mutations.json').read_text())
    require(all(abs(complex_saved(d))<=t for d,t in zip(middle['difference'],middle['tolerance'])) and middle['referenceError']<=max(middle['tolerance']),'saved selected middle integral meets original tolerance')
    require(all(max(abs(complex_saved(v)) for v in values)>1e-10 for values in mutation.values()),'saved actual sheet/Jacobian controls respond')
    require(all(v['maximumScaledResidual']<1e-11 for v in branch_joins),'saved source-root joins')
    require(all(min(v['orders'].values())>=-1 for v in endpoint_proof['orders']),'saved integrable endpoint powers')
    require(all(v['residual']==0 for v in endpoint_proof['denominatorJoins']),'saved endpoint denominator identities')
    require(len(endmaps['summary'])==8 and all(v['status']=='NUMERICAL_TRANSVERSE_CURRENT_FACE_TRACE_SUPPORTED' for v in endmaps['summary']),'accepted saved central end maps')
    uniform=manifest['uniformRefinement'];ur=Path(uniform['path'])
    for rel,v in uniform['files'].items():require(storage.digest(ur/rel)==v['sha256'],'saved uniform refinement file '+rel)
    uc=json.loads((ur/'checks.json').read_text());require(uc['status']=='UNIFORM_REFINEMENT_SUPPORTS_ORIGINAL_ACTION_TOLERANCE' and all(v['passesOriginalTolerance'] for v in uc['trials']),'saved fixed-tolerance uniform support')
    us=BlobStore(ur/'operations.sqlite')
    try:refinement=restored_blob(J,us,uniform['result'],'completed-uniform-refinement',{'database':str(ur/'operations.sqlite'),'sha256':uniform['files']['operations.sqlite']['sha256']})
    finally:us.close()
    for trial in refinement['results']:
        require(trial['passesOriginalTolerance'],'saved trial support')
        for c in trial['comparisons']:
            require(np.isfinite(c['assembled']).all() and np.all(abs(c['difference'])<=c['tolerance']) and np.array_equal(c['assembled']-c['reference'],c['difference']),'actual saved uniform difference/tolerance arrays')
    require(refinement['sourceNodeMovementInToleranceUnits']<=1,'saved source-node refinement stability')
    require({(v['setting']['outerOrder'],v['setting']['sourceNodes']) for v in refinement['results']}=={(32,512),(48,512),(64,512),(64,768)},'actual tested quadrature settings')
    H=I.native_helpers(manifest['nativeEngine'],manifest['nativeFinite'],np,sp);F.initialize(np,sp,I);r=I.context(native,frequency);zero=bindings[0.]
    require(endmaps['fieldUnits']==frequency['fieldUnits']==zero['fieldUnits'],'actual end/source/operator units')
    for a,binding in bindings.items():require(all(exact_structure(binding[key],zero[key]) for key in ('rows','jets','sources')),'saved actual contrast-independent kernel sharing')
    for name,setting in settings.items():
        setting['singleLegOuterOrder']=48 if name=='base' else 64
    J.report('numerical-settings',{'settings':settings,'coverage':coverage,'singleLegRule':'Order48 base,64 refinements from actual fixed-tolerance Gaussian convergence; nested rules unchanged.'})
    J.report('consumed-integration-checks',{'middle':middle,'mutations':mutation,'uniformSummary':uc,'completeReturnsReplayed':False,'sourceRootJoins':len(branch_joins),'endpointOrders':len(endpoint_proof['orders']),'physicalLossInterpretation':'RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED','independentMethodInterpretation':'AWAIT_REQUIRED_REPLACEMENT_METHOD_REVIEW'})
    Radiating=I.radiating_class(H)
    del native,frequency,refinement,endpoint_proof,endpoint_values;gc.collect()
    for a in bindings:J.refs[f'bind/contrast-{a}']=J.refs[f'restore/bind/contrast-{a}']
    cases={};summaries=[];solves=0
    for setting_name,setting in settings.items():
        N=setting['size'];L=setting['sourceBound'];x=np.sort(L*np.cos(np.arange(N)*np.pi/(N-1)))
        worker=Radiating(zero,r,roots,setting);worker.prepare_basis(zero['jets'],x,L,N,setting['sourceNodes'])
        source_arrays={i:J.blob(setting_name+f'/source-basis-{i}.pickle',value) for i,value in worker.amplitudes.items() if worker.amplitude_aliases[i]==i}
        source_arrays={i:source_arrays[worker.amplitude_aliases[i]] for i in worker.amplitudes}
        J.blob(setting_name+'/source-basis-operands.pickle',{'jets':zero['jets'],'nodes':worker.source_nodes,'weights':worker.source_weights,'amplitudeReceipts':source_arrays,'setting':setting,'positions':x})
        row_matrices={}
        groups=sorted({tuple(l[0] for l in row['limits']) for row in zero['rows']},key=lambda v:(len(v),tuple(map(str,v))))
        for index,variables in enumerate(groups):
            prefix=setting_name+f'/group-{index}'
            group_setting=matrix_group_setting(setting,variables);worker.setting=group_setting
            data=J.call(prefix,{'rows':[v for v in zero['rows'] if tuple(l[0] for l in v['limits'])==variables],'jets':zero['jets'],'variables':variables,'setting':group_setting,'positions':x},lambda variables=variables,prefix=prefix,group_setting=group_setting:F.matrix_group(worker,variables,group_setting,x,J,prefix))
            row_matrices.update({int(i):storage.decode(J.store.get(ref)) for i,ref in data['matrixReceipts'].items()});J.report(prefix.replace('/','-'),{'rowIndices':data['rowIndices'],'nodes':data['nodes'],'batches':data['batches'],'seconds':data['seconds'],'maximumActionResidual':F.norm(data['actionResidual']),'massResidual':data['massResidual']})
        require(set(row_matrices)==set(range(80)),'all actual eighty nonlocal matrices')
        worker.compiler.compiled.cache_clear();del worker,data;gc.collect()
        for a in ((0.,.25,.5,1.) if setting_name=='base' else (0.,1.)):
            prefix=setting_name+f'/contrast-{a}';binding=bindings[a]
            system=J.call(prefix+'/assemble',{'binding':J.refs[f'bind/contrast-{a}'],'rowMatrixGroups':[J.refs[setting_name+f'/group-{i}'] for i in range(len(groups))],'positions':x,'setting':setting},lambda:F.matrices(binding,row_matrices,x,setting,r,H))
            selected={e:values[a] for e,values in endmaps['ends'].items()}
            solves+=1;J.report('finite-solve-count',{'central':solves,'completed':len(summaries),'pilotMaximum':30})
            solved=J.call(prefix+'/solve',{'system':J.refs[prefix+'/assemble'],'endMaps':selected},lambda:F.solve(system,selected,J,prefix))
            cases[setting_name,a]={k:v for k,v in solved.items() if k in ('deficit','normalizedCrossCurrents','incident','outwardTransverseByEnd','reflected','transmitted','solveChecks','setting')}
            summary={'setting':setting_name,'contrast':a,'incident':solved['incident'],'reflected':solved['reflected'],'transmitted':solved['transmitted'],'deficit':solved['deficit'],'crossCurrents':solved['normalizedCrossCurrents'],'rank':solved['solveChecks']['rank'],'condition':solved['solveChecks']['condition'],'maximumScaledResidual':solved['solveChecks']['maximumScaledResidual'],'physicalLossInterpretation':'RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED'}
            summaries.append(summary);J.report('completed-cases',summaries);J.report('finite-solve-count',{'central':solves,'pilotMaximum':30})
            require(a!=0 or np.all(abs(solved['deficit'])<=1e-6),'UNIFORM_BACKGROUND_CONTROL_FAILED',summary)
            del solved,system;gc.collect()
        del row_matrices;gc.collect()
    report=J.call('finite-controls',cases,lambda:F.control_report(cases));J.report('central-transverse-balance',report)
    artifact=J.blob('completed-central-pilot.pickle',{'report':report,'cases':cases,'frequency':3,'fieldUnits':zero['fieldUnits'],'rowUnits':zero['rowUnits'],'currentUnit':endmaps['currentUnit']});J.report('result-artifacts',{'centralPilot':artifact})
    bad=any(s in ('NEGATIVE_DEFICIT_UNRESOLVED','PARTITION_UNRESOLVED') for row in report['cases'] for s in row['statuses'])
    return {'status':'CENTRAL_FINITE_PILOT_UNRESOLVED_STOP_NEIGHBORS' if bad else 'CENTRAL_FINITE_PILOT_COMPLETED_INSPECTION_PENDING','frequency':3,'case':'LAB_HELD/RHO4_CONSTANT','finiteSolves':solves,'summaries':summaries,'report':report,'neighborsRun':False,'analogLightCalibration':'OPEN','independentMethodInterpretation':'AWAIT_REQUIRED_REPLACEMENT_METHOD_REVIEW'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--run-directory',type=Path,required=True);a=p.parse_args()
    manifest=json.loads(a.manifest.read_text());gate=json.loads(a.gate.read_text())
    require(gate['status']=='READY_FOR_AUTHORIZED_CENTRAL_NUMERICAL_PILOT','actual launch readiness')
    require(gate['workerSha256']==storage.digest(__file__) and gate['manifestSha256']==storage.digest(a.manifest),'worker and manifest pins')
    for path,h in gate['sourcePins'].items():require(storage.digest(path)==h,'source/helper pin '+path)
    require(all(gate[k] is None for k in ('wallDeadlineSeconds','nativeDeadlineSeconds','inactivityDeadlineSeconds')),'standing unlimited runtime')
    require(gate['independentMethodClearance'] is False and gate['proceedAuthority']=='EXPLICIT_USER_NO_FURTHER_GROK_MINOR_CORRECTIONS_AND_RUN','literal review/authority status')
    base=a.run_directory.resolve();require(str(base)==manifest['resultDirectory'],'declared result directory');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));started=time.monotonic();code=0;J=None
    try:
        global np,sp
        import numpy as np
        import sympy as sp
        storage.np=np;storage.sp=sp;warnings.simplefilter('error',RuntimeWarning)
        from scipy.linalg import LinAlgWarning
        warnings.simplefilter('error',LinAlgWarning)
        J=storage.Journal(base);result=run(manifest,J)
    except BaseException as error:
        code=1;result={'status':'STOPPED_UNRESOLVED','exceptionType':type(error).__name__,'message':str(error),'traceback':traceback.format_exc(),'incompleteOperation':J.active if J else 'initialization'}
        count=base/'finite-solve-count.json';result['finiteSolves']=json.loads(count.read_text())['central'] if count.exists() else 0
        if J and getattr(error,'evidence',None) is not None:
            ref=J.blob('final-failed-check-evidence.pickle',error.evidence);result['failedCheckEvidence']=ref
        storage.save(base/'failure.json',result)
    pins=dict(gate['sourcePins']);pins.update({v['path']:v['sha256'] for v in manifest['packets'].values()});pins[manifest['endMaps']['database']]=manifest['endMaps']['sha256']
    pins.update({str(Path(manifest['priorStartup']['path'])/p):v['sha256'] for p,v in manifest['priorStartup']['files'].items()})
    pins.update({str(Path(manifest['uniformRefinement']['path'])/p):v['sha256'] for p,v in manifest['uniformRefinement']['files'].items()})
    post={path:{'expected':h,'actual':storage.digest(path)} for path,h in pins.items()};intact=all(v['expected']==v['actual'] for v in post.values());storage.save(base/'posthashes.json',post)
    if not intact:code=1;result['status']='INTEGRITY_FAILURE'
    if J:J.store.integrity_check();J.store.close()
    inventory={str(path.relative_to(base)):{'sha256':storage.digest(path),'bytes':path.stat().st_size} for path in sorted(base.rglob('*')) if path.is_file()};storage.save(base/'artifact-index.json',inventory)
    result.update(wallSeconds=time.monotonic()-started,completeOperations=J.count if J else 0,restoredOperations=J.restored if J else 0,posthashesIntact=intact,peakRssKiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,independentMethodClearance=False,physicalLossInterpretation='RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED',automaticRetry=False)
    text=json.dumps(storage.readable(result),indent=2,allow_nan=False)+'\n';(base/'checks.json').write_text(text);sys.stdout.write(text);sys.stdout.flush();return code

if __name__=='__main__':sys.exit(main())
