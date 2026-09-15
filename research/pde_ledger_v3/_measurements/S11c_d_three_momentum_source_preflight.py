#!/usr/bin/env python3
"""Compare computed finite prefixes with saved native operands; estimate cost."""
import argparse,itertools,json,time,resource
import numpy as np
import S11c_d_three_momentum_source_check as s

class PrefixComplete(Exception):pass

def prefix(worker,ti,variables,setting,bound,width):
    saved=[]
    def stop(state):
        if state['batchCount']==64:saved.append(state);raise PrefixComplete()
    started=time.monotonic()
    try:worker.group(ti,variables,setting,bound['pairs'],width,bound['positions'],stop)
    except PrefixComplete:pass
    if len(saved)!=1:raise ValueError('prefix not completed')
    return saved[0],time.monotonic()-started

def main():
    p=argparse.ArgumentParser();p.add_argument('--directory',type=s.Path,required=True);a=p.parse_args()
    base=a.directory.resolve();base.relative_to(s.STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    r,bound,rows,variables,finest,provenance,pins=s.load(base)
    worker=s.engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
    native=s.engine.BoundedSourceFourierQuadrature.CachedMomentum(bound['rows'],bound['sources'],r)
    accepted,old=s.momentum.source.accepted(s.M/'S11c_d_three_momentum_checkpoint.json')
    records=[];partial_sources=[];timings=[]
    for ti in range(2):
        setting=dict(finest[ti]['setting'],innerOrders=(24,24));width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
        path=old/'partials'/f'{ti}-4-000064.pickle'
        record=next(v for v in accepted['partialArtifacts'] if v['path']==str(path.relative_to(old)))
        if s.digest(path)!=record['sha256']:raise ValueError('saved prefix hash')
        original=s.unpickle(path);partial_sources.append(dict(record,absolutePath=str(path)))
        if original['variables']!=variables or original['rowIndices']!=[v['index'] for v in rows] or original['setting']!=finest[ti]['setting'] or accepted['provenance']['BOUND_PACKET_SHA256']!=provenance['BOUND_PACKET_SHA256']:raise ValueError('original prefix joins')
        current,seconds=prefix(worker,ti,variables,setting,bound,width);timings.append({'test':ti,'kind':'baselinePrefix','seconds':seconds,'nodes':current['nodeCount']})
        node_difference=weight_difference=0.
        for new,prior in zip(itertools.islice(worker.batches(variables,setting,bound['pairs'],width),64),itertools.islice(native.batches(variables,dict(setting,panelOrder=24),bound['pairs'],width),64),strict=True):
            node_difference=max(node_difference,float(np.max(abs(new[0]-prior[0]))));weight_difference=max(weight_difference,float(np.max(abs(new[1]-prior[1]))))
        group=next(g for g in finest[ti]['groups'] if g['variables']==variables)
        action=s.full_action(r,bound,finest[ti],group)
        records.append({'test':ti,'currentPrefix':current,'savedPrefix':original,'prefixValueResidual':current['values']-original['values'],
            'prefixMutationResidual':current['measureMutationValues']-original['measureMutationValues'],
            'prefixMassResidual':current['quadratureMass']-original['quadratureMass'],
            'nodesResidualMax':node_difference,'weightsResidualMax':weight_difference,'retainedBaselineAction':action})
    order_controls=[]
    for orders in ((32,24),(24,32)):
        altered=dict(setting,innerOrders=orders)
        changed=next(worker.batches(variables,altered,bound['pairs'],width))[0]
        baseline=next(worker.batches(variables,setting,bound['pairs'],width))[0]
        order_controls.append({'orders':orders,'baselineNodes':baseline,'changedNodes':changed,'residual':changed-baseline})
    highest=dict(setting,outerOrder=144,innerOrders=(24,24),sourceNodes=256,profileNodes=256)
    for ti in range(2):
        current,seconds=prefix(worker,ti,variables,highest,bound,width)
        timings.append({'test':ti,'kind':'finestPrefix','seconds':seconds,'nodes':current['nodeCount']})
        s.atomic_pickle(base/f'finest-prefix-{ti}.pickle',current)
    caches=[]
    for (points,order),(nodes,weights) in worker.fixed_rules.items():
        original=s.engine.BoundedSourceFourierQuadrature.rule(points,order)
        caches.append({'points':points,'order':order,'profileRule':(points,order) in worker.profile_rule_keys,'nodesEqual':np.array_equal(nodes,original[0]),'weightsEqual':np.array_equal(weights,original[1]),'readOnly':not nodes.flags.writeable and not weights.flags.writeable})
    emission_records=[]
    for ti in range(2):
        group=next(g for g in finest[ti]['groups'] if g['variables']==variables)
        emission_records.append({'test':ti,'index':0,'method':'retained','setting':dict(finest[ti]['setting'],innerOrders=(24,24)),
            'group':group,**s.full_action(r,bound,finest[ti],group)})
    cache_arrays=[]
    for (points,order),(nodes,weights) in worker.fixed_rules.items():
        nn,ww=s.engine.BoundedSourceFourierQuadrature.rule(points,order)
        cache_arrays.append({'points':points,'order':order,'nodesResidual':nodes-nn,'weightsResidual':weights-ww,
            'readOnly':not nodes.flags.writeable and not weights.flags.writeable,'bytes':nodes.nbytes+weights.nbytes,
            'profileRule':(points,order) in worker.profile_rule_keys})
    smoke={'records':emission_records,'cacheChecks':cache_arrays,
        'phaseWorkspaceEstimateBytes':worker.peak_phase_workspace_estimate,'batchCacheEstimateBytes':worker.peak_batch_cache_estimate,
        'sourceRuleCacheBytes':sum(v['bytes'] for v in cache_arrays),'workspaceBudgetBytes':worker.workspace_bytes,
        'evaluatedProfileIntegrals':tuple(sorted(worker.profile_integrals,key=s.sp.default_sort_key)),
        'sourceFrequencyCensus':worker.source_frequency_census}
    s.atomic_pickle(base/'emission-operands.pickle',smoke)
    entries,keys,metadata_paths=s.emit_and_replay(base,smoke,r,bound,rows,finest,provenance)
    s.atomic_pickle(base/'focused.pickle',{'records':records,'orderControls':order_controls,'caches':caches,'timings':timings,'partialSources':partial_sources})
    checks={'runDirectory':str(base),'sourceFiles':pins,'provenance':provenance,'nativeEngineAstJoin':True,'unchangedEngineSourceJoin':True,
        'rows':len(rows),'tests':2,'prefixNodesPerTest':16384,'prefixSources':partial_sources,'timings':timings,
        'maxPrefixValueResidual':max(float(np.max(abs(v['prefixValueResidual']))) for v in records),
        'maxPrefixMutationResidual':max(float(np.max(abs(v['prefixMutationResidual']))) for v in records),
        'maxPrefixMassResidual':max(abs(v['prefixMassResidual']) for v in records),
        'maxNodeResidual':max(v['nodesResidualMax'] for v in records),'maxWeightResidual':max(v['weightsResidualMax'] for v in records),
        'maxRetainedActionResidual':max(float(np.max(abs(v['retainedBaselineAction']['differenceFromReference']))) for v in records),
        'minPrefixMaxMutation':min(float(np.max(abs(v['currentPrefix']['measureMutationValues']-v['currentPrefix']['values']))) for v in records),
        'orderControlMaxima':[float(np.max(abs(v['residual']))) for v in order_controls],
        'cacheChecks':caches,'focusedPacketSha256':s.digest(base/'focused.pickle'),'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'wallSeconds':time.monotonic()-started,
        'emissionTags':len(entries),'emissionWriteKeys':len(keys),'emissionMetadataPaths':metadata_paths,'emissionTranscriptSha256':s.digest(base/'full.out'),
        'scope':'Computed prefixes and retained full-baseline joins only; no new full three-momentum grid or convergence result.'}
    s.save(base/'checks.json',checks)
    if any(checks[k]!=0 for k in ('maxPrefixValueResidual','maxPrefixMutationResidual','maxPrefixMassResidual','maxNodeResidual','maxWeightResidual','maxRetainedActionResidual')):raise ValueError('prefix/baseline mismatch')
    if checks['minPrefixMaxMutation']<=0 or min(checks['orderControlMaxima'])<=0 or any(not v[k] for v in caches for k in ('nodesEqual','weightsEqual','readOnly')):raise ValueError('cache or sensitivity check')
    if pins!={n:s.digest(s.ROOT/n) for n in pins}:raise ValueError('focused sources changed')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
