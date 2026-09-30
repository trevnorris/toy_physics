#!/usr/bin/env python3
"""Central real-frequency finite pilot; source restoration only under guard.

No source producers, root/path replay, symbolic epsilon extraction or time cap.
Stops on a substantive integration/background/solve failure; no retries.
"""
import argparse
import pickle
import gc
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback
import warnings
import S11c_d_numerical_radiating_end_maps_v2 as storage
import S11c_d_numerical_radiating_integration_v3 as I
import S11c_d_numerical_radiating_finite as F
from S11c_d_numerical_radiating_blob_store import BlobStore

require=I.require

from S11c_d_numerical_radiating_saved_operations import SavedOperationsJournal as SavedStartupJournal, exact_structure


def load_packet(item,base):
    path=Path(item['path']);require(storage.digest(path)==item['sha256'],'accepted source packet hash')
    target=base/'inputs'/path.name;target.parent.mkdir(exist_ok=True);shutil.copyfile(path,target)
    require(storage.digest(target)==item['sha256'],'opaque source copy hash')
    return storage.decode(target.read_bytes())


def saved_end_maps(manifest,J):
    item=manifest['endMaps'];require(storage.digest(item['database'])==item['sha256'],'accepted end-map database hash')
    store=BlobStore(item['database'])
    try:raw=store.get(item['result'])
    finally:store.close()
    receipt=J.store.put('inputs/completed-end-maps.pickle',raw);maps=storage.decode(raw)
    require(len(maps['summary'])==8 and all(v['status']=='NUMERICAL_TRANSVERSE_CURRENT_FACE_TRACE_SUPPORTED' for v in maps['summary']),'saved bounded central end support')
    evidence=[]
    for end,cases in maps['ends'].items():
        require(set(cases)=={0.,.25,.5,1.},'all four saved contrast end maps')
        for contrast,case in cases.items():
            current=case['currents'];checks=case['checks']
            evidence.append({'end':end,'contrast':contrast,'fullCurrentForms':current['forms'],'hermitianResiduals':current['hermitianResiduals'],'lossSideMixed':current['lossSideMixed'],'traceResidual':case['traceMap']['residual'],'selectedControls':{str(i):{k:v for k,v in c.items() if k in ('faceProjections','bulkProjection','interfacePowerProjection','nativeVelocityOmissions','currentOffDiagonalOmission','orientationFlipMovement','fluxNormalizationResidual')} for i,c in checks.items()}})
    J.report('consumed-end-map-evidence',{'receipt':receipt,'savedValuesOnly':True,'evidence':evidence})
    return maps


def run(manifest,J):
    packets={name:J.call('restore/'+name,item,lambda item=item:load_packet(item,J.base)) for name,item in manifest['packets'].items()}
    endmaps=J.call('restore/completed-end-maps',manifest['endMaps'],lambda:saved_end_maps(manifest,J))
    native=packets['native']['bound'];frequency=packets['frequency'];grades=packets['grades']
    H=I.native_helpers(manifest['nativeEngine'],manifest['nativeFinite'],np,sp);F.initialize(np,sp,I)
    r=I.context(native,frequency)
    J.report('reused-numerical-helper-asts',H.proofs)
    require(endmaps['fieldUnits']==frequency['fieldUnits'],'actual source/end field-unit identity')
    bindings={}
    for a in (0.,.25,.5,1.):
        bindings[a]=J.call(f'bind/contrast-{a}',{'frequencyPacket':manifest['packets']['frequency'],'nativePacket':manifest['packets']['native'],'gradesPacket':manifest['packets']['grades'],'frequency':3,'contrast':a},lambda a=a:I.bind_sources(frequency,native,grades,a,r,H))
    del packets,native,grades;gc.collect()
    zero=bindings[0.]
    for a,binding in bindings.items():
        require(all(exact_structure(binding[key],zero[key]) for key in ('rows','jets','sources')),'actual row/source operands independent of contrast')
    J.report('kernel-contrast-reuse',{'contrasts':list(bindings),'actualRowsJetsAndSourceCharactersEqual':True,'localAndCellCoefficientsBoundSeparately':True})
    roots=J.call('integration/radical-inventory',frequency['sources']['rowCensus'],lambda:I.radical_inventory(frequency['sources'],J))
    settings,coverage=J.call('discretization/coverage',endmaps['ends'],lambda:F.settings(endmaps['ends']));J.report('numerical-settings',{'settings':settings,'coverage':coverage})
    physical=json.loads(Path(manifest['physicalInput']['path']).read_text())
    proof=J.call('integration/endpoint-order-census',{'boundRows':zero['rows'],'roots':roots,'parameters':physical['parameters']},lambda:I.endpoint_orders(zero,roots,physical['parameters'],J))
    Radiating=I.radiating_class(H);worker=Radiating(zero,r,roots,settings['base']);positions=np.asarray([-2.,0.,3.])
    J.call('integration/source-branch-joins',{'rows':zero['rows'],'roots':roots,'positions':positions},lambda:I.source_branch_joins(worker,zero,roots,positions,J))
    J.call('integration/transformed-measure',{'setting':settings['base'],'roots':roots},lambda:I.check_rule(worker,roots,settings['base']))
    boundary=J.call('integration/Gaussian-source',{'jets':zero['jets'],'sourceBound':64.,'sourceNodes':512,'width':8.,'momentum':0.},lambda:worker.prepare_gaussian(zero['jets'],64.,512))
    J.blob('integration/Gaussian-source-arrays.pickle',{'nodes':worker.source_nodes,'weights':worker.source_weights,'amplitudes':worker.amplitudes,'boundaryValues':boundary})
    endpoint=J.call('integration/two-sided-endpoint-values',{'rows':zero['rows'],'roots':roots,'positions':positions,'setting':settings['base']},lambda:I.endpoint_values(worker,zero,roots,positions,J));J.report('two-sided-endpoint-values',endpoint)
    middle=I.middle_check(worker,zero,roots,positions,J)
    uniform=I.uniform_gaussian(worker,zero,frequency,roots,positions,J)
    J.blob('integration/completed-checks.pickle',{'endpointProof':proof,'endpointValues':endpoint,'middle':middle,'uniform':uniform})
    J.report('integration-status',{'status':'CENTRAL_INTEGRATION_CHECKS_SUPPORTED','finiteSolves':0,'physicalLossInterpretation':'RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED'})
    worker.compiler.compiled.cache_clear();del worker,proof,uniform,middle,frequency;gc.collect()
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
            data=J.call(prefix,{'rows':[v for v in zero['rows'] if tuple(l[0] for l in v['limits'])==variables],'jets':zero['jets'],'variables':variables,'setting':setting,'positions':x},lambda variables=variables,prefix=prefix:F.matrix_group(worker,variables,setting,x,J,prefix))
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
    return {'status':'CENTRAL_FINITE_PILOT_UNRESOLVED_STOP_NEIGHBORS' if bad else 'CENTRAL_FINITE_PILOT_COMPLETED_INSPECTION_PENDING','frequency':3,'case':'LAB_HELD/RHO4_CONSTANT','finiteSolves':solves,'summaries':summaries,'report':report,'neighborsRun':False,'analogLightCalibration':'OPEN'}


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
        J=SavedStartupJournal(base,manifest['priorStartup']);result=run(manifest,J)
    except BaseException as error:
        code=1;result={'status':'STOPPED_UNRESOLVED','exceptionType':type(error).__name__,'message':str(error),'traceback':traceback.format_exc(),'incompleteOperation':J.active if J else 'initialization'}
        count=base/'finite-solve-count.json';result['finiteSolves']=json.loads(count.read_text())['central'] if count.exists() else 0
        if J and getattr(error,'evidence',None) is not None:
            ref=J.blob('final-failed-check-evidence.pickle',error.evidence);result['failedCheckEvidence']=ref
        storage.save(base/'failure.json',result)
    pins=dict(gate['sourcePins']);pins.update({v['path']:v['sha256'] for v in manifest['packets'].values()});pins[manifest['endMaps']['database']]=manifest['endMaps']['sha256']
    pins.update({str(Path(manifest['priorStartup']['path'])/p):v['sha256'] for p,v in manifest['priorStartup']['files'].items()})
    post={path:{'expected':h,'actual':storage.digest(path)} for path,h in pins.items()};intact=all(v['expected']==v['actual'] for v in post.values());storage.save(base/'posthashes.json',post)
    if not intact:code=1;result['status']='INTEGRITY_FAILURE'
    if J:J.store.integrity_check();J.store.close()
    inventory={str(path.relative_to(base)):{'sha256':storage.digest(path),'bytes':path.stat().st_size} for path in sorted(base.rglob('*')) if path.is_file()};storage.save(base/'artifact-index.json',inventory)
    result.update(wallSeconds=time.monotonic()-started,completeOperations=J.count if J else 0,restoredOperations=J.restored if J else 0,posthashesIntact=intact,peakRssKiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,independentMethodClearance=False,physicalLossInterpretation='RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED',automaticRetry=False)
    text=json.dumps(storage.readable(result),indent=2,allow_nan=False)+'\n';(base/'checks.json').write_text(text);sys.stdout.write(text);sys.stdout.flush();return code

if __name__=='__main__':sys.exit(main())
