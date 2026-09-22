#!/usr/bin/env python3
"""Bounded saved structural source-coefficient review; no coefficient arithmetic."""
import argparse
import json
from pathlib import Path
import resource
import signal
import time

import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_remaining_case_frequency_end_continuation_inputs as io

M,REPO=io.M,io.REPO
ROOT=REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921/source-coefficients'
require,same=io.require,io.same


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic();reader=io.Reader()
    inspection=reader.json(ROOT/'completion-inspection.json');checks=reader.json(ROOT/'complete/checks.json',inspection['checksSha256'])
    require(inspection['actualExits']==[0,0,0] and inspection['emptyStrictStderr'] and inspection['checksStdoutIdentity'] and inspection['zeroCapOOMSwap'],'clean producer final inspection')
    for name,r in inspection['evidence'].items():reader.retain(r['path'],r['sha256'])
    guard=reader.json(ROOT/'resource-guard/outcome.json');child=reader.json(ROOT/'resource-guard/child-outcome.json');inv=reader.json(ROOT/'frequency_source_coefficients.invocation.json')
    require(guard['exitCode']==child['exitCode']==inv['exitCode']==0 and child['guardReason'] is None,'actual producer exits')
    for p in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_source_coefficients_review_plan.md'):reader.retain(p)
    cache={}
    def packet(record):
        p=record.get('logical',record.get('path'));r=reader.retain(p,record['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(p,record['sha256'])
        return cache[r['canonical']]
    def artifact(name):
        r=checks['artifacts'][name]
        return reader.json(ROOT/'complete'/name,r['sha256']) if name.endswith('.json') else packet(r)
    for name,r in checks['artifacts'].items():require(reader.retain(ROOT/'complete'/name,r['sha256'])['bytes']==r['bytes'],'producer artifact')
    inputs=artifact('inputs.json')
    for p,r in inputs['consumedRoutes'].items():require(reader.retain(p,r['sha256'])==r,'all consumed original routes')
    source=artifact('native-caller-and-new-projection.json');require(not source['nativeSourceJetsCalled'],'explicit new projection identity')
    calls=[]
    for n in range(checks['counts']['coefficientOperationCalls']):
        prefix='coefficient-operations/'+str(n);inp=artifact(prefix+'/input.pickle');value=artifact(prefix+'/value.pickle');receipt=artifact(prefix+'/completed.json')
        require(receipt['input']==checks['artifacts'][prefix+'/input.pickle'] and receipt['value']==checks['artifacts'][prefix+'/value.pickle'],'immediate actual operation receipts')
        require(inp['call']['function'] in ('Add','Mul') and inp['call']['kwargs']=={} and receipt['function']==inp['call']['function'],'only new declared coefficient operations')
        if receipt['kind']=='REUSED_NEW_COMPLETE_CALL':
            owner=receipt['owner'];require(0<=owner<n and same(calls[owner]['input']['call'],inp['call']) and same(calls[owner]['value'],value),'saved reuse owner full arguments and actual value')
        else:require(receipt['kind']=='NEW_COEFFICIENT_OPERATION' and receipt['owner']==n,'actual new operation owner')
        calls.append({'input':inp,'value':value,'receipt':receipt})
    families=[]
    for family in range(checks['counts']['newFullProjections']):
        prefix='families/'+str(family);inp=artifact(prefix+'/input.pickle');full=inp['call'];value=artifact(prefix+'/value.pickle');receipt=artifact(prefix+'/completed.json')
        request=inp['consumer']['fullInputRoutes'];old=packet(request['nativeSource']['packet']);bound=packet(request['boundAmplitude']['packet']);si=request['sourceIndex']
        source_record=old['bound']['sources'][0,si];probe=old['jets'][si]['probe'];coordinate=probe.args[0]
        require(same(full,{'expression':bound['actual']['source',si],'probe':probe,'coordinate':coordinate,
                           'amplitudeUnit':source_record['amplitudeUnit'],'integralUnit':source_record['integralUnit']}),'complete saved physical input')
        require(receipt['input']==checks['artifacts'][prefix+'/input.pickle'] and receipt['value']==checks['artifacts'][prefix+'/value.pickle'],'full projection receipt')
        records=artifact(prefix+'/projection-tree.pickle');nodes={};used=set()
        for record in records:
            node=record['node'];actual=full['expression']
            for index in record['argsAddress']:actual=actual.args[index]
            require(same(actual,node),'exact input subtree address')
            result=record['result'];kind=record['kind'];operations=record['operations']
            if kind=='actual-saved-probe-node':
                if same(node,probe):order=0
                else:
                    require(isinstance(node,sp.Derivative) and node.expr==probe and all(x==coordinate for x,n in node.variable_count),'saved actual derivative node')
                    order=int(sum(n for x,n in node.variable_count))
                require(same(result,{order:sp.S.One}),'structural probe coefficient identity')
            elif kind=='unchanged-independent-subtree':
                require(not node.has(AppliedUndef,sp.Derivative) and same(result,{None:node}),'verbatim independent source subtree')
            else:
                children=[nodes[v] for v in node.args]
                if kind=='linear-add':
                    require(isinstance(node,sp.Add),'native saved sum shape');orders=set().union(*(v.keys() for v in children))
                    args={k:tuple(v[k] for v in children if k in v) for k in orders};function='Add'
                else:
                    require(kind=='linear-product' and isinstance(node,sp.Mul),'native saved product shape')
                    dep=[i for i,v in enumerate(children) if any(k is not None for k in v)]
                    require(len(dep)==1,'exact homogeneous linear product');index=dep[0]
                    args={k:tuple(coefficient if i==index else children[i][None] for i in range(len(children))) for k,coefficient in children[index].items()};function='Mul'
                require(set(result)==set(args)==set(operations),'complete projected derivative orders')
                for order,operands in args.items():
                    op=operations[order]
                    if len(operands)==1:
                        require(op=={'kind':'single-saved-operand'} and same(result[order],operands[0]),'no single-operand arithmetic')
                    else:
                        n=next(i for i,c in enumerate(calls) if c['receipt']['input']==op['input']);call=calls[n];used.add(n)
                        context={k:full[k] for k in ('probe','coordinate','amplitudeUnit','integralUnit')};context['jetOrder']=order
                        require(call['receipt']==op and same(call['input']['call'],{'function':function,'args':operands,'kwargs':{},'unitContext':context}) and same(result[order],call['value']),'complete before/after coefficient wiring without arithmetic replay')
                        require(call['input']['consumer']=={'family':family,'argsAddress':record['argsAddress'],'jetOrder':order},'own actual coefficient consumer')
            nodes[node]=result
        root=nodes[full['expression']]
        require(None not in root and set(root)==set(range(len(value['coefficients']))) and all(same(root[n],a) for n,a in enumerate(value['coefficients'])),'complete homogeneous source coefficient return')
        require(same(value['originalBoundAmplitude'],full['expression']) and same(value['probe'],probe) and same(value['amplitudeUnit'],full['amplitudeUnit']) and same(value['integralUnit'],full['integralUnit']),'full returned source/unit identity')
        require(value['column']==int(probe.func.__name__.removeprefix('s11cdPencilProbe')) and 'residual' not in value,'actual field identity; no fabricated native residual')
        ci=artifact(prefix+'/comparison-input.pickle');cv=artifact(prefix+'/comparison-value.pickle');summary=artifact(prefix+'/comparison.json')
        require(same(ci['source'],full) and same(ci['coefficients'],value['coefficients']) and ci['tolerance']==2e-12,'actual numeric comparison full inputs')
        require(len(cv)==summary['actualComparisonCount']==len(ci['points'])*len(ci['jetVectors'])==12,'bounded comparison census')
        require(summary['maximumScaledDifference']==max(v['scaledDifference'] for v in cv)<2e-12 and summary['coefficientMutationResponse']>1e-5 and summary['jetOrderOrSignMutationResponse']>1e-12,'saved direct comparisons and actual mutation responses')
        for row in cv:
            require(row['point'] in ci['points'] and row['jetVector'] in ci['jetVectors'] and len(row['coefficientValues'])==len(value['coefficients']),'actual coordinate and jet comparison routes')
            require(row['coefficientMutationValues'][1:]==row['coefficientValues'][1:] and row['coefficientMutationValues'][0]!=row['coefficientValues'][0],'actual changed numeric coefficient and unaffected fields')
            require(row['jetOrderMutationValues']!=row['coefficientValues'],'actual saved derivative-order/sign mutation')
        families.append({'input':full,'value':value,'used':used})
        io.save(base,'validated-family-'+str(family)+'.json',{'input':receipt['input'],'value':receipt['value'],'projectionNodes':len(records),'actualOperationSlots':sorted(used),'comparison':summary,'newScientificCalls':0})
    require(set().union(*(v['used'] for v in families))==set(range(len(calls))),'all actual arithmetic slots consumed')
    counts={'sourceUses':0,'savedNativeUses':0,'newFullProjections':0,'sharedNewUses':0}
    for name in checks['artifacts']:
        if not name.endswith('/source-coefficient-routes.json'):continue
        routes=artifact(name)
        for route in routes:
            counts['sourceUses']+=1;request=route['consumer']['fullInputRoutes'];si=request['sourceIndex'];old=packet(request['nativeSource']['packet']);bound=packet(request['boundAmplitude']['packet']);s=old['bound']['sources'][0,si];probe=old['jets'][si]['probe']
            full={'expression':bound['actual']['source',si],'probe':probe,'coordinate':probe.args[0],'amplitudeUnit':s['amplitudeUnit'],'integralUnit':s['integralUnit']}
            target=route['result'];value=packet(target['packet'])
            for k in target['keys']:value=value[k]
            require(same((value['originalBoundAmplitude'],value['probe'],value['amplitudeUnit'],value['integralUnit']),
                         (full['expression'],probe,full['amplitudeUnit'],full['integralUnit'])),'complete own source coefficient consumer')
            if route['kind']=='SAVED_NATIVE_COMPLETE_JET':
                counts['savedNativeUses']+=1;require(any(same(target,{'packet':v['packet'],'keys':v['keys']}) for v in request['matches']),'literal full accepted return')
            else:
                matches=[v for v in families if same(v['input'],full) and same(v['value'],value)];require(len(matches)==1,'full new coefficient input and return')
                counts['newFullProjections' if route['kind']=='NEW_STRUCTURAL_PROJECTION' else 'sharedNewUses']+=1
        io.save(base,name.replace('source-coefficient-routes','validated-routes'),{'sourceUses':len(routes),'allFullInputsJoined':True})
    require(all(counts[k]==checks['counts'][k] for k in counts),'independent consumer counts')
    reader.postcheck();io.save(base,'validated-paths.json',reader.routes)
    result={'status':'PASSED_BOUNDED_SAVED_SOURCE_COEFFICIENT_REVIEW','producerChecksSha256':inspection['checksSha256'],
        'counts':counts,'actualOperationReceipts':len(calls),'maximumSavedScaledComparison':checks['maximumScaledDirectComparison'],
        'newScientificCalls':0,'allConsumedHashesUnchanged':True,'consumedLogicalPaths':len(reader.routes),'wallSeconds':time.monotonic()-start,
        'scope':'Saved structural source coefficients at1-0.01i. No symbolic/numeric coefficient, differentiation, polynomial, quadrature, map or response recomputation.'}
    io.save(base,'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
