#!/usr/bin/env python3
"""Refine two ordered momenta while retaining accepted other action layouts."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
import sympy as sp
import S11c_d_momentum_action_check as momentum
from S11c_d_momentum_action_check import ROOT,STORE,engine,digest,save,atomic_pickle,unpickle,number,progress
from S11c_d_output_codec import decoded_lines,restore_emission_index
from ledger_fold import _restore

M=ROOT/'_measurements';CHECKPOINT=M/'S11c_d_single_momentum_checkpoint.json'
PLAN=M/'S11c_d_pair_momentum_plan.md';PREFIX='PAIR_MOMENTUM_LAB_HELD_RHO4_CONSTANT'


def load(base):
    accepted,previous=momentum.source.accepted(CHECKPOINT)
    for name,sha in accepted['sourceFiles'].items():
        if digest(previous/'source'/name)!=sha or (name!=str(engine.HERE.relative_to(ROOT)) and digest(ROOT/name)!=sha):
            raise ValueError(('accepted momentum source changed',name))
    for key in ('recordArtifacts',):
        for record in accepted[key]:
            if digest(previous/record['path'])!=record['sha256']:raise ValueError('accepted momentum saved operand changed')
    old=ast.parse((previous/'source'/engine.HERE.relative_to(ROOT)).read_text());current=ast.parse(engine.HERE.read_text())
    cls=next(n for n in current.body if getattr(n,'name',None)=='BoundedSourceFourierQuadrature')
    additions=[n for n in cls.body if getattr(n,'name',None)=='CachedMomentum']
    if len(additions)!=1:raise ValueError('cached momentum helper count')
    helper=additions[0]
    profile=next(n for n in helper.body if getattr(n,'name',None)=='profile_value')
    native=next(n for n in cls.body if getattr(n,'name',None)=='FiniteMomentum')
    old_profile=next(n for n in native.body if getattr(n,'name',None)=='profile_value')
    calls=[n for n in ast.walk(profile) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='profile_rule']
    if len(calls)!=1:raise ValueError('profile cache call census')
    calls[0].func=ast.Attribute(value=ast.Name(id='BoundedSourceFourierQuadrature',ctx=ast.Load()),attr='rule',ctx=ast.Load())
    if ast.dump(profile)!=ast.dump(old_profile):raise ValueError('profile evaluator changed beyond exact rule cache')
    cls.body.remove(additions[0])
    if ast.dump(old)!=ast.dump(current):raise ValueError('engine changed beyond cached momentum helper')
    bound_packet=unpickle(previous/'accepted-bound-momentum.pickle');numerical=unpickle(previous/'accepted-momentum-action.pickle')
    single=unpickle(previous/'single-momentum.pickle')
    if numerical['boundPacketSha256']!=digest(previous/'accepted-bound-momentum.pickle'):raise ValueError('bound/numerical packet join')
    cp,source_base=momentum.source.accepted(momentum.source.SOURCE)
    r,dimensions=momentum.source.native.source.restore_context(unpickle(source_base/'reduced-action.pickle'))
    dimensions.__dict__.update(bound_packet['dimensionState'])
    bound=bound_packet['result'];values=numerical['result']
    if single['boundPacketSha256']!=digest(previous/'accepted-bound-momentum.pickle') or accepted['provenance']['BOUND_PACKET_SHA256']!=single['boundPacketSha256']:raise ValueError('single/momentum packet provenance join')
    rows=[row for row in bound['rows'] if len(row['limits'])==2];variables=tuple(l[0] for l in rows[0]['limits'])
    if len(rows)!=30 or any(row['limits']!=rows[0]['limits'] for row in rows):raise ValueError('two-momentum native ordered-layout census')
    source_ids=sorted({f['sourceIndex'] for row in rows for f in row['factors']})
    original={ti:next(item for item in numerical['result']['records'] if item['test']==ti and item['grid']==4) for ti in range(2)}
    finest={}
    for ti,item in original.items():
        refined=next(v for v in single['result']['records'] if v['test']==ti and v['method']=='gauss' and v['index']==5)
        if refined['setting']['outerOrder']!=216 or refined['setting']['sourceNodes']!=256:raise ValueError('accepted single rule join')
        reference=dict(item,action=refined['action'],contributions=[{'index':v['index'],'terms':v['values']} for v in refined['terms']])
        reference['groups']=[refined['group'] if len(g['variables'])==1 else g for g in item['groups']]
        reference['layoutSettings']={len(g['variables']):(refined['setting'] if len(g['variables'])==1 else item['setting']) for g in item['groups']}
        group=next(g for g in item['groups'] if g['variables']==variables)
        if group['rowIndices']!=[row['index'] for row in rows]:raise ValueError('reference paired row ordering')
        finest[ti]=reference
    provenance={'SINGLE_CHECKPOINT_SHA256':digest(CHECKPOINT),'BOUND_PACKET_SHA256':digest(previous/'accepted-bound-momentum.pickle'),
        'NUMERICAL_PACKET_SHA256':digest(previous/'accepted-momentum-action.pickle'),'APPROVED_INPUT_SHA256':accepted['provenance']['APPROVED_INPUT_SHA256'],'SINGLE_PACKET_SHA256':digest(previous/'single-momentum.pickle')}
    paths=tuple(dict.fromkeys((*(ROOT/n for n in accepted['sourceFiles']),CHECKPOINT,PLAN,Path(__file__).resolve(),M/'S11c_d_momentum_refinement_diagnostic.json',M/'S11c_d_pair_momentum_preflight.py')))
    pins={str(p.relative_to(ROOT)):digest(p) for p in paths}
    for path in paths:
        target=base/'source'/path.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
    for old_name,new_name in (('accepted-bound-momentum.pickle','accepted-bound-momentum.pickle'),('accepted-momentum-action.pickle','accepted-momentum-action.pickle'),('single-momentum.pickle','accepted-single-momentum.pickle')):
        shutil.copyfile(previous/old_name,base/new_name)
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'nativeEngineAstJoin':True,'profileEvaluatorAstJoin':True,
        'sourceRunDirectory':str(previous),'rowIndices':[row['index'] for row in rows],'sourceIndices':source_ids})
    return r,bound,rows,variables,finest,provenance,pins


def full_action(r,bound,reference,group):
    values={row:group['values'][j] for j,row in enumerate(group['rowIndices'])}
    old_terms={tuple(record['index']):record['terms'] for record in reference['contributions']}
    compiler=engine.BoundedActionQuadrature({});action=reference['local'].copy();records=[]
    for cell in (c for c in bound['cells'] if c['test']==reference['test']):
        for pi,z in enumerate(bound['positions']):
            index=(pi,cell['column'],cell['row']);fresh=[];changed=[]
            for ni,(row,coef) in enumerate(cell['terms']):
                if row in values:
                    entry=complex(compiler.evaluate(coef,{r.z:z,r.regulator:reference['setting']['regulator']}))*values[row][pi]
                    changed.append(True)
                else:entry=old_terms[index][ni];changed.append(False)
                fresh.append(entry)
            action[index]+=sum(fresh)
            records.append({'index':index,'rows':[row for row,_ in cell['terms']],'changed':changed,
                'reference':old_terms[index],'values':fresh,'differences':np.asarray(fresh)-np.asarray(old_terms[index])})
    return {'action':action,'terms':records,'differenceFromReference':action-reference['action']}


def construct(base,r,bound,rows,variables,finest):
    contractor=engine.BoundedSourceFourierQuadrature.CachedMomentum(bound['rows'],bound['sources'],r)
    records=[];inventory=[];partial_inventory=[];adaptive_inventory=[]
    choices=((64,16,128,128),(96,16,128,128),(144,16,128,128),(216,16,128,128),(216,24,128,128),(216,32,128,128),(216,32,192,128),(216,32,256,128),(216,32,256,192),(216,32,256,256))
    for ti in range(2):
        previous=None
        for index,(order,panel_order,source_order,profile_order) in enumerate(choices):
            setting=dict(finest[ti]['setting'],outerOrder=order,panelOrder=panel_order,sourceNodes=source_order,profileNodes=profile_order)
            width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
            def partial(state):
                if state['batchCount']%64==0:
                    path=base/'partials'/f'{ti}-{index}-{state["batchCount"]:06}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,state)
                    partial_inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)});save(base/'partial-inventory.json',partial_inventory)
            group=contractor.group(ti,variables,setting,bound['pairs'],width,bound['positions'],partial)

            item={'test':ti,'index':index,'method':'gauss','setting':setting,'group':group,
                **full_action(r,bound,finest[ti],group)}
            reference_group=next(g for g in finest[ti]['groups'] if g['variables']==variables)
            item['referenceIntegralDifference']=group['values']-reference_group['values']
            if previous is not None:
                item['precedingIntegralDifference']=group['values']-previous['group']['values']
                item['precedingActionDifference']=item['action']-previous['action']
            path=base/'records'/f'{ti}-{index}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,item)
            inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)});save(base/'record-inventory.json',inventory)
            if abs(group['volumeResidual'])>1e-11:raise ValueError('pair momentum mass')
            if np.any((np.abs(group['values'])>1e-9)&(np.abs(group['measureMutationResidual'])<1e-12)):
                raise ValueError('resolved paired value lost measure sensitivity')
            if index==0 and (np.max(abs(item['referenceIntegralDifference']))>1e-10 or np.max(abs(item['differenceFromReference']))>1e-10):
                raise ValueError('cached-rule reference layout/action join')
            for record in item['terms']:
                if any(v!=0 for changed,v in zip(record['changed'],record['differences']) if not changed):raise ValueError('held native term changed')
            records.append(item);previous=item;progress(base,'pair_momentum_saved',test=ti,index=index,outerOrder=order,sourceOrder=source_order)
        point_count=0
        def save_outer(k,values):
            nonlocal point_count
            path=base/'outer-points'/f'{ti}-{point_count:06}.pickle';path.parent.mkdir(exist_ok=True)
            atomic_pickle(path,{'test':ti,'variables':variables,'outerValue':k,'innerValues':values,'setting':setting,'positions':bound['positions']})
            adaptive_inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)});save(base/'outer-point-inventory.json',adaptive_inventory);point_count+=1
        adaptive=contractor.adaptive_outer(ti,variables,setting,bound['pairs'],width,bound['positions'],checkpoint=save_outer)

        item={'test':ti,'index':len(choices),'method':'adaptive','setting':setting,'group':adaptive,
            **full_action(r,bound,finest[ti],adaptive),
            'precedingIntegralDifference':adaptive['values']-previous['group']['values']}
        item['precedingActionDifference']=item['action']-previous['action']
        path=base/'records'/f'{ti}-adaptive.pickle';atomic_pickle(path,item)
        inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)});save(base/'record-inventory.json',inventory)
        if not adaptive['success'] or not np.all(np.isfinite(adaptive['values'])) or not np.isfinite(adaptive['unitFrameErrorEstimate']):
            raise ValueError('adaptive outer quadrature unfinished')
        records.append(item);progress(base,'adaptive_momentum_saved',test=ti,evaluations=int(adaptive['evaluations']))
    cache_checks=[]
    for (points,order),(nodes,weights) in contractor.fixed_rules.items():
        native_nodes,native_weights=engine.BoundedSourceFourierQuadrature.rule(points,order)
        cache_checks.append({'points':points,'order':order,'nodesResidual':nodes-native_nodes,'weightsResidual':weights-native_weights,
            'readOnly':not nodes.flags.writeable and not weights.flags.writeable,'bytes':nodes.nbytes+weights.nbytes,'profileRule':(points,order) in contractor.profile_rule_keys})
    if any(np.any(v['nodesResidual']) or np.any(v['weightsResidual']) or not v['readOnly'] for v in cache_checks):raise ValueError('exact source rule cache join')
    expected_sources={(ti,f['sourceIndex']) for ti in range(2) for row in rows for f in row['factors']}
    if set(contractor.source_frequency_census)!=expected_sources:raise ValueError('single source evaluation census')
    profiles=set().union(*(f['coefficient'].atoms(sp.Integral) for row in rows for f in row['factors']))
    if contractor.profile_integrals!=profiles or not profiles<=set(bound['profileUnits']):raise ValueError('paired nested profile census')
    result={'records':records,'recordArtifacts':inventory,'cacheChecks':cache_checks,'sourceFrequencyCensus':contractor.source_frequency_census,
        'evaluatedProfileIntegrals':tuple(sorted(profiles,key=sp.default_sort_key)),
        'partialArtifacts':partial_inventory,'outerPointArtifacts':adaptive_inventory,'sourceRuleCacheBytes':sum(v['bytes'] for v in cache_checks),'phaseWorkspaceEstimateBytes':contractor.peak_phase_workspace_estimate,
        'batchCacheEstimateBytes':contractor.peak_batch_cache_estimate,'workspaceBudgetBytes':contractor.workspace_bytes}
    atomic_pickle(base/'pair-momentum.pickle',{'result':result,'boundPacketSha256':digest(base/'accepted-bound-momentum.pickle'),
        'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    return result


def emit(result,r,bound,rows,finest,provenance):
    dimensions=engine.PHYSICAL_METADATA.dimensions;zero=dimensions.zero
    metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=r
    def numeric(name,value,unit,literal=False):
        body=number(value);engine.emit(PREFIX+'_'+name,body if literal else metadata.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    engine.physical(PREFIX+'_PROVENANCE',provenance)
    numeric('CHANGED_ROW_INDICES',[row['index'] for row in rows],lambda p:zero,True)
    numeric('TEST_POSITIONS',bound['positions'],lambda p:dimensions.measure(r.z),True)
    numeric('TEST_WIDTHS_MOMENTA',bound['tests'],lambda p:tuple((1 if p[1]==0 else -1)*v for v in dimensions.measure(r.z)),True)
    for row in rows:engine.fingerprinted(PREFIX+'_ORIGINAL_INTEGRAL_'+str(row['index']),row['original'])
    for ti,reference in finest.items():
        numeric('REFERENCE_ACTION_'+str(ti),reference['action'],lambda p:bound['equationUnits'][p[2]])
        engine.physical(PREFIX+'_HELD_RULE_KIND_'+str(ti),reference['setting']['kind'])
        numeric('REFERENCE_RULE_ORDERS_'+str(ti),[[n,*[setting[k] for k in ('panelOrder','outerOrder','sourceNodes','profileNodes')]] for n,setting in sorted(reference['layoutSettings'].items())],lambda p:zero,True)
        for gi,g in enumerate(reference['groups']):
            numeric('REFERENCE_ROWS_'+str(ti)+'_'+str(gi),g['rowIndices'],lambda p:zero,True)
            numeric('REFERENCE_INTEGRALS_'+str(ti)+'_'+str(gi),g['values'],lambda p:bound['rows'][g['rowIndices'][p[0]]]['unit'])
    for item in result['records']:
        suffix=str(item['test'])+'_'+str(item['index']);group=item['group'];setting=item['setting']
        engine.physical(PREFIX+'_METHOD_'+suffix,item['method'])
        numeric('CHANGED_RULE_ORDERS_'+suffix,(setting['outerOrder'],setting['panelOrder'],setting['sourceNodes'],setting['profileNodes']),lambda p:zero,True)
        numeric('SOURCE_BOUNDS_'+suffix,(-setting['sourceBound'],setting['sourceBound']),lambda p:dimensions.measure(r.zp),True)
        numeric('MOMENTUM_BOUNDS_'+suffix,(-setting['momentumBound'],setting['momentumBound']),lambda p:bound['momentumUnit'],True)
        numeric('PROFILE_BOUNDS_'+suffix,(-setting['profileBound'],setting['profileBound']),lambda p:dimensions.measure(r.xi),True)
        numeric('REGULATOR_'+suffix,setting['regulator'],lambda p:dimensions.measure(r.regulator),True)
        integral_unit=lambda p:bound['rows'][group['rowIndices'][p[0]]]['unit']
        numeric('INTEGRAL_VALUES_'+suffix,group['values'],integral_unit)
        for name in ('referenceIntegralDifference','precedingIntegralDifference'):
            if name in item:numeric(name.upper()+'_'+suffix,item[name],integral_unit,True)
        for name in ('action','differenceFromReference','precedingActionDifference'):
            if name in item:numeric(name.upper()+'_'+suffix,item[name],lambda p:bound['equationUnits'][p[2]],name!='action')
        if item['method']=='gauss':
            numeric('MEASURE_MUTATION_OPERAND_'+suffix,group['measureMutationValues'],integral_unit)
            numeric('MEASURE_MUTATION_RESIDUAL_'+suffix,group['measureMutationResidual'],integral_unit,True)
            numeric('BOX_MASS_OPERANDS_'+suffix,(group['quadratureMass'],group['boxVolume'],group['volumeResidual']),lambda p:tuple(2*v for v in bound['momentumUnit']),True)
            numeric('NODE_COUNT_'+suffix,group['nodeCount'],lambda p:zero,True)
        else:
            numeric('ADAPTIVE_ERROR_TOLERANCE_'+suffix,(group['unitFrameErrorEstimate'],group['unitFrameTolerance']),lambda p:zero,True)
            numeric('ADAPTIVE_EVALUATIONS_STATUS_'+suffix,(group['evaluations'],group['status'],int(group['success'])),lambda p:zero,True)
            numeric('ADAPTIVE_INTERVALS_'+suffix,group['intervals'],lambda p:bound['momentumUnit'],True)
            numeric('ADAPTIVE_INTERVAL_VALUES_'+suffix,group['intervalValues'],lambda p:bound['rows'][group['rowIndices'][p[1]]]['unit'])
            numeric('ADAPTIVE_INTERVAL_UNIT_FRAME_ERRORS_'+suffix,group['unitFrameIntervalErrorEstimates'],lambda p:zero,True)
        for ni,term in enumerate(item['terms']):
            key=suffix+'_'+str(ni);unit=lambda p:bound['equationUnits'][term['index'][2]]
            numeric('CELL_'+key,term['index'],lambda p:zero,True)
            numeric('TERM_ROWS_'+key,term['rows'],lambda p:zero,True)
            numeric('CHANGED_TERM_MASK_'+key,list(map(int,term['changed'])),lambda p:zero,True)
            numeric('REFERENCE_TERMS_'+key,term['reference'],unit)
            numeric('TERMS_'+key,term['values'],unit)
            numeric('TERM_DIFFERENCES_'+key,term['differences'],unit,True)
    for index,cache in enumerate(result['cacheChecks']):
        suffix=str(index)
        cache_unit=dimensions.measure(r.xi if cache['profileRule'] else r.zp)
        numeric('CACHED_RULE_KIND_'+suffix,int(cache['profileRule']),lambda p:zero,True)
        numeric('CACHED_PANELS_'+suffix,cache['points'],lambda p:cache_unit,True)
        numeric('CACHED_SOURCE_ORDER_'+suffix,cache['order'],lambda p:zero,True)
        for name in ('nodesResidual','weightsResidual'):numeric(name.upper()+'_'+suffix,cache[name],lambda p:cache_unit,True)
        numeric('SOURCE_RULE_READONLY_'+suffix,int(cache['readOnly']),lambda p:zero,True)
    numeric('WORKSPACE_CACHE_BYTES',(result['phaseWorkspaceEstimateBytes'],result['batchCacheEstimateBytes'],result['sourceRuleCacheBytes'],result['workspaceBudgetBytes']),lambda p:zero,True)
    for i,profile in enumerate(result['evaluatedProfileIntegrals']):
        engine.emit(PREFIX+'_EVALUATED_PROFILE_'+str(i),engine.carrier_fingerprint(profile))
        engine.emit('METADATA_'+PREFIX+'_EVALUATED_PROFILE_'+str(i),metadata.numeric_metadata(profile,lambda p:bound['profileUnits'][profile]))
    for (ti,si),value in sorted(result['sourceFrequencyCensus'].items()):
        numeric('SOURCE_FREQUENCIES_'+str(ti)+'_'+str(si),value,lambda p:zero if p[0]==2 else bound['momentumUnit'],True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    r,bound,rows,variables,finest,provenance,pins=load(base)
    before={n:digest(base/n) for n in ('accepted-bound-momentum.pickle','accepted-momentum-action.pickle','accepted-single-momentum.pickle')}
    result=construct(base,r,bound,rows,variables,finest);before['pair-momentum.pickle']=digest(base/'pair-momentum.pickle')
    with (base/'full.out').open('x') as out,contextlib.redirect_stdout(out):
        emit(result,r,bound,rows,finest,provenance)
        keys={tag:'s11cd'+''.join(w.title() for w in tag.removeprefix('PY_S11CD_').split('_')) for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')}
        momentum.source.emit_manifest(PREFIX+'_WRITE_KEYS',keys)
        index=engine.emission_index(engine.EMISSION_LINES);zero_units={p:(0,0,0) for p,_ in engine.leaves(engine.cas(index))}
        momentum.source.emit_manifest(PREFIX+'_EMISSION_LINES',index,zero_dimensions=zero_units)
    entries={}
    for line in decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ')
        if tag in entries:raise ValueError('duplicate pair momentum tag')
        entries[tag]=_restore(body)
    seen=set();old=engine.emit
    def compare(name,value):
        tag='PY_S11CD_'+name
        if tag in seen or entries.get(tag)!=engine.cas(value):raise ValueError(('pair momentum replay',tag))
        seen.add(tag)
    engine.emit=compare
    try:
        emit(result,r,bound,rows,finest,provenance);momentum.source.emit_manifest(PREFIX+'_WRITE_KEYS',keys)
        momentum.source.emit_manifest(PREFIX+'_EMISSION_LINES',index,zero_dimensions=zero_units)
    finally:engine.emit=old
    if seen!=set(entries) or len(keys)!=len(set(keys.values())) or set(keys.values())&set(engine.IMPORT_KEYS):raise ValueError('pair momentum key/replay census')
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    metadata_paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        for record in body:
            if len(record)==2 and isinstance(record[0],sp.Tuple):fields={str(k):v for k,v in record[1]};count=1
            else:fields={str(k):v for k,v in record};count=len(fields['PATHS'])
            unit=fields['DIMENSION_L_T_M']
            if len(unit)!=3 or any(v.free_symbols for v in unit) or any(k not in fields for k in ('MULTIGRADE','EPSILON_LAMBDA_SUPPORT')):raise ValueError('pair momentum metadata')
            metadata_paths+=count
    norms=[]
    for item in result['records']:
        row={'test':item['test'],'index':item['index'],'method':item['method'],'outerOrder':item['setting']['outerOrder'],'sourceOrder':item['setting']['sourceNodes'],'panelOrder':item['setting']['panelOrder'],'profileOrder':item['setting']['profileNodes']}
        for name in ('referenceIntegralDifference','precedingIntegralDifference','differenceFromReference','precedingActionDifference'):
            if name in item:row[name]=float(np.max(np.abs(item[name])))
        if item['method']=='adaptive':row.update(errorEstimate=float(item['group']['unitFrameErrorEstimate']),evaluations=int(item['group']['evaluations']),success=bool(item['group']['success']))
        norms.append(row)
    summary={'runDirectory':str(base),'sourceFiles':pins,'provenance':provenance,'nativeEngineAstJoin':True,'profileEvaluatorAstJoin':True,'rowCount':len(rows),'profileCount':len(result['evaluatedProfileIntegrals']),
        'sourceCount':len(result['sourceFrequencyCensus']),'records':len(result['records']),'norms':norms,'ruleCaches':len(result['cacheChecks']),
        'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':metadata_paths,'recordArtifacts':result['recordArtifacts'],'partialArtifacts':result['partialArtifacts'],'outerPointArtifacts':result['outerPointArtifacts'],
        'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'scope':'Two-momentum outer/panel/source/profile refinements and adaptive outer comparison, with accepted single-momentum and original finest three-momentum terms explicitly held. No full-action tail, regulator limit, scattering or pole result.'}
    save(base/'checks.json',summary)
    for record in (*result['recordArtifacts'],*result['partialArtifacts'],*result['outerPointArtifacts']):
        if digest(base/record['path'])!=record['sha256']:raise ValueError('pair momentum saved record changed')
    if pins!={n:digest(ROOT/n) for n in pins} or before!={n:digest(base/n) for n in before} or engine.PHYSICAL_METADATA.dimensions.constraints:
        raise ValueError('pair momentum final source/packet/dimension guard')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
