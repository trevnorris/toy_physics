#!/usr/bin/env python3
"""Focused uniform-action refinement; saved reference, unchanged equations.

Scientific imports/restoration only inside shared guard and supervisor.
No path/root/current/middle-integral replay, finite solve or time limit.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback
import warnings
import S11c_d_numerical_radiating_end_maps_v2 as storage
import S11c_d_numerical_radiating_integration_v3 as I
from S11c_d_numerical_radiating_blob_store import BlobStore
from S11c_d_numerical_radiating_saved_operations import exact_structure
require=I.require


def restored_blob(J,source,ref,name,prior):
    raw=source.get(ref)
    inp=J.blob('restore/'+name+'/input.pickle',{'source':prior,'receipt':ref})
    out=J.store.put('restore/'+name+'/return.pickle',raw)
    value=storage.decode(raw)
    record={'name':'restore/'+name,'input':inp,'return':out,'status':'RESTORED_PRIOR_COMPLETE_RETURN','functionCalled':False,'sourceReceipt':ref,'sourceIdentity':'EXACT_HASH_VERIFIED_RETURN_BYTES','prior':prior}
    with (J.base/'operation-index.jsonl').open('a') as f:
        f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
    J.refs['restore/'+name]=out;J.count+=1;J.restored+=1
    return value


def gaussian_setup(worker,jets,bound,order,width,momentum):
    boundary=worker.prepare_gaussian(jets,bound,order,width,momentum)
    return {'nodes':worker.source_nodes,'weights':worker.source_weights,'amplitudes':worker.amplitudes,'size':worker.size,'boundaryValues':boundary,'scope':'Previously unsaved ephemeral source quadrature setup; unchanged helper, not a producer or completed integral.'}


def assemble(binding,row_values,positions,column,momentum,r):
    assembled=np.zeros((5,len(positions)),complex)
    for n,m in binding['local'].items():
        d=I.gaussian_derivative(positions,n,8.,momentum)
        for i in range(5):assembled[i]+=np.asarray(sp.lambdify(r.z,m[i,column],'numpy',cse=True)(positions),complex)*d
    for term in binding['termJoins']:
        if term['column']!=column:continue
        c=binding['actual']['cell',term['row'],column,term['term']]
        if c!=0:assembled[term['row']]+=np.asarray(sp.lambdify(r.z,c,'numpy',cse=True)(positions),complex)*row_values[term['integralIndex']]
    return assembled


def run(manifest,J):
    prior=manifest['priorStartup'];root=Path(prior['path'])
    for rel,v in prior['files'].items():require(storage.digest(root/rel)==v['sha256'],'preserved complete input '+rel)
    operations=[json.loads(x) for x in (root/'operation-index.jsonl').read_text().splitlines()]
    require(len(operations)==36 and all(v['status'] in ('COMPLETE','RESTORED_PRIOR_COMPLETE_RETURN') for v in operations),'all complete prior operation receipts')
    old={v['name']:v for v in operations};source=BlobStore(root/'operations.sqlite')
    provenance={'database':str(root/'operations.sqlite'),'sha256':prior['files']['operations.sqlite']['sha256']}
    def restore(name):return restored_blob(J,source,old[name]['return'],name,dict(provenance,operation=old[name]))
    try:
        frequency=restore('restore/frequency');native=restore('restore/native')['bound'];binding=restore('bind/contrast-0.0')
        roots=restore('integration/radical-inventory');settings,coverage=restore('discretization/coverage')
        refs={};inputs={};baselines={}
        for column in (3,4):
            name='uniform-Gaussian/0.6000000000000001/reference-column-'+str(column)
            refs[column]=restore(name)
            inputs[column]=restored_blob(J,source,old[name]['input'],name+'-operands',provenance)
            item=manifest['baselineComparisons'][str(column)]
            baselines[column]=restored_blob(J,source,item,'baseline-comparison-'+str(column),provenance)
        for filename in manifest['evidenceFiles']:
            raw=(root/filename).read_bytes();(J.base/('prior-'+filename)).write_bytes(raw)
            require(storage.digest(J.base/('prior-'+filename))==prior['files'][filename]['sha256'],'saved readable evidence copy')
    finally:source.close()
    H=I.native_helpers(manifest['nativeEngine'],manifest['nativeFinite'],np,sp);r=I.context(native,frequency)
    require(binding['contrast']==0 and binding['fieldUnits']==frequency['fieldUnits'],'actual zero-contrast/source unit joins')
    base=settings['base'];positions=inputs[4]['positions'];momentum=inputs[4]['momentum'];width=inputs[4]['width']
    require(width==8 and abs(momentum-.6)<1e-15 and exact_structure(positions,np.asarray([-2.,0.,3.])),'selected saved Gaussian operands')
    matrix=frequency['ends']['REFERENCE']['pencil']['livePencil'].subs(frequency['ends']['REFERENCE']['pencil']['frequency'],3)
    scale=next(iter(roots.values()))['scale']
    for column in (3,4):
        x=inputs[column];b=baselines[column]
        require(x['column']==column and x['K']==base['momentumBound'] and x['width']==width and x['momentum']==momentum and exact_structure(x['positions'],positions),'saved reference argument joins')
        require(x['matrix']==matrix and x['physicalScale']==scale,'actual saved reference pencil and physical scale')
        require(exact_structure(b['reference'],refs[column]['value']) and b['column']==column and b['momentum']==momentum,'saved reference/baseline identity')
    J.report('source-context-joins',{'frequency':3,'contrast':0,'momentum':momentum,'width':width,'positions':positions,'fieldUnits':binding['fieldUnits'],'rowUnits':binding['rowUnits'],'baseSetting':base,'referenceColumns':[3,4],'referenceRecomputed':False,'middleIntegralsRecomputed':False,'sourceSetup':'Only unsaved ephemeral source arrays recreated from saved jets at the declared quadrature orders.'})
    active=[v for v in binding['termJoins'] if binding['actual']['cell',v['row'],v['column'],v['term']]!=0]
    used=sorted({v['integralIndex'] for v in active});require(used==[14,27,29,35,36],'same five active uniform rows')
    expressions=[binding['actual']['cell',v['row'],v['column'],v['term']] for v in active]
    expressions += [f['coefficient'] for i in used for f in binding['rows'][i]['factors']]
    expressions += [x for matrix in binding['local'].values() for x in matrix]
    require(all(not x.has(r.regulator) for x in expressions),'unchanged regulator-free uniform operator')
    Radiating=I.radiating_class(H);worker=Radiating(binding,r,roots,base);setup_cache={};results=[]
    for trial in manifest['trials']:
        setting=dict(base,outerOrder=trial['outerOrder'],sourceNodes=trial['sourceNodes']);worker.setting=setting;prefix=trial['label'];source_order=trial['sourceNodes']
        if source_order not in setup_cache:
            setup_cache[source_order]=J.call('source-setup/'+str(source_order),{'jets':binding['jets'],'sourceBound':setting['sourceBound'],'sourceNodes':source_order,'width':width,'momentum':momentum},lambda:gaussian_setup(worker,binding['jets'],setting['sourceBound'],source_order,width,momentum))
        setup=setup_cache[source_order];worker.source_nodes=setup['nodes'];worker.source_weights=setup['weights'];worker.amplitudes=setup['amplitudes'];worker.size=setup['size']
        row_values={}
        for i in used:
            row=binding['rows'][i]
            item=J.call(prefix+'/row-'+str(i),{'row':row,'setting':setting,'positions':positions,'sourceSetup':J.refs['source-setup/'+str(source_order)]},lambda row=row:worker.row_integral(row,positions))
            row_values[i]=item['value']
        comparisons=[]
        for column in (3,4):
            assembled=J.call(prefix+'/assemble-column-'+str(column),{'binding':J.refs['restore/bind/contrast-0.0'],'rowReturns':{i:J.refs[prefix+'/row-'+str(i)] for i in used},'positions':positions,'column':column,'momentum':momentum},lambda column=column:assemble(binding,row_values,positions,column,momentum,r))
            reference=refs[column]['value'];tolerance=1e-8+1e-6*abs(reference);difference=assembled-reference
            require(np.isfinite(assembled).all() and np.isfinite(reference).all(),'finite uniform-action comparison')
            passed=bool(np.max(abs(reference))>1e-10 and np.all(abs(difference)<=tolerance) and refs[column]['error']<=np.max(tolerance))
            record={'column':column,'momentum':momentum,'setting':setting,'assembled':assembled,'reference':reference,'difference':difference,'tolerance':tolerance,'referenceError':refs[column]['error'],'passedOriginalTolerance':passed,'maximumAbsoluteDifference':float(np.max(abs(difference))),'maximumToleranceRatio':float(np.max(abs(difference)/tolerance))}
            J.blob(prefix+'/comparison-column-'+str(column)+'.pickle',record);J.report(prefix+'-column-'+str(column),record);comparisons.append(record)
        result={'label':prefix,'setting':setting,'comparisons':comparisons,'passesOriginalTolerance':all(v['passedOriginalTolerance'] for v in comparisons)};results.append(result)
        J.report('refinement-progress',[{'label':v['label'],'passesOriginalTolerance':v['passesOriginalTolerance'],'maximumToleranceRatio':max(c['maximumToleranceRatio'] for c in v['comparisons'])} for v in results])
    last=results[-1];previous=results[-2]
    source_movement=max(float(np.max(abs(a['assembled']-b['assembled'])/a['tolerance'])) for a,b in zip(last['comparisons'],previous['comparisons']))
    adjacent_pass=any(a['passesOriginalTolerance'] and b['passesOriginalTolerance'] for a,b in zip(results[:-2],results[1:-1]))
    supported=last['passesOriginalTolerance'] and previous['passesOriginalTolerance'] and adjacent_pass and source_movement<=1
    artifact=J.blob('completed-uniform-refinement.pickle',{'results':results,'baseline':baselines,'references':refs,'sourceNodeMovementInToleranceUnits':source_movement,'fieldUnits':binding['fieldUnits'],'rowUnits':binding['rowUnits']});J.report('result-artifacts',{'uniformRefinement':artifact})
    summary={'status':'UNIFORM_REFINEMENT_SUPPORTS_ORIGINAL_ACTION_TOLERANCE' if supported else 'UNIFORM_REFINEMENT_ORIGINAL_TOLERANCE_UNRESOLVED','frequency':3,'contrast':0,'momentum':momentum,'trials':[{'label':v['label'],'outerOrder':v['setting']['outerOrder'],'sourceNodes':v['setting']['sourceNodes'],'passesOriginalTolerance':v['passesOriginalTolerance'],'maximumToleranceRatio':max(c['maximumToleranceRatio'] for c in v['comparisons'])} for v in results],'sourceNodeMovementInToleranceUnits':source_movement,'finiteSolves':0,'operatorAccuracyIsFluxErrorBound':False,'requestedCoarseFluxResolution':1e-5,'coarseResolutionStatus':'Requires actual finite uniform/scaling/refinement/domain/regulator/sign controls and independent method review; operator-action accuracy alone is not a flux bound.','physicalLossInterpretation':'RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED'}
    J.report('uniform-refinement-summary',summary);return summary


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
    post={path:{'expected':h,'actual':storage.digest(path)} for path,h in pins.items()};intact=all(v['expected']==v['actual'] for v in post.values());storage.save(base/'posthashes.json',post)
    if not intact:code=1;result['status']='INTEGRITY_FAILURE'
    if J:J.store.integrity_check();J.store.close()
    inventory={str(path.relative_to(base)):{'sha256':storage.digest(path),'bytes':path.stat().st_size} for path in sorted(base.rglob('*')) if path.is_file()};storage.save(base/'artifact-index.json',inventory)
    result.update(wallSeconds=time.monotonic()-started,completeOperations=J.count if J else 0,restoredOperations=J.restored if J else 0,posthashesIntact=intact,peakRssKiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,independentMethodClearance=False,physicalLossInterpretation='RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED',automaticRetry=False)
    text=json.dumps(storage.readable(result),indent=2,allow_nan=False)+'\n';(base/'checks.json').write_text(text);sys.stdout.write(text);sys.stdout.flush();return code

if __name__=='__main__':sys.exit(main())
