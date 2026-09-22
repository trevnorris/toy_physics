#!/usr/bin/env python3
"""Bounded saved scalar review and source-jet input catalogue; no science calls."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time

import S11c_d_remaining_case_frequency_end_continuation_inputs as io

M,REPO=io.M,io.REPO
ROOT=REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921/scalar-bindings'
SHA='87a03a6ec9d2f89e9a0cb9e146d0d3f5d616bd3975e9723878c36c8ff8d17b46'
require,same=io.require,io.same


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    base=parser.parse_args().run_directory.resolve();base.relative_to(REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    reader=io.Reader();checks=reader.json(ROOT/'complete/checks.json',SHA)
    def metadata(name,value):io.save(base,name,value)
    def artifact(name):
        p=ROOT/'complete'/name;record=checks['artifacts'][name]
        return reader.json(p,record['sha256']) if name.endswith('.json') else reader.packet(p,record['sha256'])
    for name,v in checks['artifacts'].items():
        require(reader.retain(ROOT/'complete'/name,v['sha256'])['bytes']==v['bytes'],'all producer artifact bytes')
    guard=reader.json(ROOT/'resource-guard/outcome.json');child=reader.json(ROOT/'resource-guard/child-outcome.json');inv=reader.json(ROOT/'frequency_scalar_bindings.invocation.json')
    require(guard['exitCode']==child['exitCode']==inv['exitCode']==0 and child['guardReason'] is None and guard['limitsVerified'],'actual clean producer exits')
    limits=reader.json(ROOT/'resource-guard/effective-limits.json');validation=reader.json(ROOT/'resource-guard/limit-validation.json')
    require(validation['verified'] and limits['memory.max']=='2147483648' and limits['memory.swap.max']=='0' and limits['pids.max']=='32'
            and limits['nice']==15 and len(limits['affinity'])==1 and set(limits['threads'].values())=={'1'},'mandatory containment')
    for n in ('frequency_scalar_bindings.stderr','guard.stderr','resource-guard/stderr'):
        require(reader.retain(ROOT/n)['bytes']==0,'empty strict stderr')
    require(reader.retain(ROOT/'frequency_scalar_bindings.stdout')['sha256']==SHA,'checks/stdout identity')
    telemetry=ROOT/'resource-guard/resource-samples.jsonl';reader.retain(telemetry)
    for s in map(json.loads,telemetry.read_text().splitlines()):
        events=dict(line.split() for line in s['memory.events'].splitlines())
        require(int(s['memory.swap.current'])==0 and int(s['memory.peak'])<=2147483648 and
                all(int(events[k])==0 for k in ('max','oom','oom_kill')),'zero cap/OOM/swap')
    inputs=artifact('inputs.json')
    for p,v in inputs['consumedRoutes'].items():require(reader.retain(p,v['sha256'])==v,'consumed original route')
    for p in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_scalar_review_plan.md'):
        reader.retain(p)
    cp=reader.json(inputs['acceptedAnalyticCheckpoint']['logical'],inputs['acceptedAnalyticCheckpoint']['sha256']);origin=Path(cp['runDirectory'])
    mcp=reader.json(inputs['baselineFrequencyMatrix']['logical'],inputs['baselineFrequencyMatrix']['sha256']);mr=Path(mcp['runDirectory'])
    def accepted(name):return reader.packet(origin/name,cp['artifacts'][name]['sha256'])
    def route(name):return reader.retain(origin/name,cp['artifacts'][name]['sha256'])
    baseline_name='complex/frequency-binding.pickle';baseline_route=reader.retain(mr/baseline_name,mcp['artifacts'][baseline_name]['sha256'])
    bound=reader.packet(mr/baseline_name,baseline_route['sha256']);frequency=bound['frequency']
    sourcefile=M/'S11c_d_finite_scattering.py';native_sha=mcp['sourceFiles']['_measurements/S11c_d_finite_scattering.py']
    reader.retain(sourcefile,native_sha);reader.retain(mr/'source/_measurements/S11c_d_finite_scattering.py',native_sha)
    source_tree=ast.parse(sourcefile.read_text());source_jets=next(n for n in source_tree.body if isinstance(n,ast.FunctionDef) and n.name=='source_jets')
    metadata('native-source-jet-caller.json',{'file':reader.retain(sourcefile),'body':ast.unparse(source_jets),'called':False})
    native=artifact('native-scalar-binding-join.json');require(native['wholeNativeBindReverseAST'],'whole native scalar caller receipt')
    baseline_analytic=accepted('analytic-cases/LAB_HELD__RHO4_CONSTANT/frequency-analytic.pickle')
    baseline_by_key={v['key']:v for v in bound['sourceJoins']}
    new_inputs={};new_values={}
    for number in range(checks['newScalarCalls']):
        prefix=f'operations/{number}';arg=artifact(prefix+'/input.pickle');value=artifact(prefix+'/value.pickle');receipt=artifact(prefix+'/completed.json')
        require(receipt['value']==checks['artifacts'][prefix+'/value.pickle'],'immediate saved scalar return')
        require(arg['call']['kwargs']=={'simultaneous':True},'literal native substitution keyword')
        new_inputs[checks['artifacts'][prefix+'/value.pickle']['path']]=arg
        new_values[checks['artifacts'][prefix+'/value.pickle']['path']]=value
        metadata(f'validated-operation-{number}.json',{'input':checks['artifacts'][prefix+'/input.pickle'],'value':receipt['value'],
                 'ownerCase':arg['consumer']['case'],'ownerKey':arg['consumer']['key'],'savedValueType':str(type(value)),
                 'evidence':'Actual before/after packets under unchanged native call; substitution not repeated.'})
    def requested(record,variables):
        return {'expression':record['analytic'],'mapping':{variables['frequency']:frequency,**variables['origin']},'kwargs':{'simultaneous':True},
                'liveSource':record['originalLive'],'unit':record['unit'],'limits':record['limits']}
    jet_atlas=[]
    for index,jet in bound['jets'].items():
        jet_atlas.append({'jet':jet,'owner':{'packet':baseline_route,'keys':['jets',index],'kind':'accepted-baseline-complex-source-jet'}})
    native_rows={}
    for label in cp['cases']:
        name='accepted/frequency-cases/'+label+'/original-row-inputs.pickle';old=accepted(name);native_rows[label]=(old,route(name))
        for index,jet in old['jets'].items():
            jet_atlas.append({'jet':jet,'owner':{'packet':route(name),'keys':['jets',index],'kind':'accepted-own-seed-source-jet'}})
    summaries={};jet_requests=[];total_new_uses=0
    for label in cp['cases']:
        chart=baseline_analytic if label=='LAB_HELD__RHO4_CONSTANT' else accepted('analytic-cases/'+label+'/frequency-analytic.pickle')
        values=artifact('cases/'+label+'/scalar-bindings.pickle');routes=artifact('cases/'+label+'/binding-routes.pickle')
        by_key={v['key']:v for v in values['joins']}
        require(set(by_key)==set(chart['records']) and len(routes)==len(by_key),'all own scalar addresses')
        counts={'records':len(by_key),'newUses':0,'reusedUses':0}
        for r in routes:
            key=r['consumer']['key'];record=chart['records'][key];join=by_key[key];owner=r['owner']
            request=requested(record,chart['variables'])
            require(r['consumer']['case']==label and same(r['consumer']['address'],record['address']) and
                    same(join['analytic'],record['analytic']) and same(join['unit'],record['unit']) and same(join['limits'],record['limits']), 'full actual scalar consumer')
            require(same(values['actual'][tuple(record['address'])],join['bound']) and same(values['mapping'][join['original']],join['bound']), 'native result maps/joins')
            if owner['kind']=='accepted-baseline-frequency-binding':
                old=baseline_by_key[owner['key']];old_record=baseline_analytic['records'][owner['key']]
                require(same(request,requested(old_record,baseline_analytic['variables'])) and same(join['bound'],old['bound']), 'complete actual baseline scalar call/return')
            else:
                p=owner['value']['path'];require(p in new_inputs and same(request,new_inputs[p]['call']) and same(join['bound'],new_values[p]), 'complete actual new scalar call/return')
                if r['kind']=='NEW_COMPLETE_CALL':
                    require(new_inputs[p]['consumer']['case']==label and new_inputs[p]['consumer']['key']==key,'actual first new owner')
            counts['newUses' if r['kind']=='NEW_COMPLETE_CALL' else 'reusedUses']+=1
        total_new_uses+=counts['newUses']
        end_routes=artifact('cases/'+label+'/end-map-routes.json')
        require(set(end_routes)=={'LEFT','RIGHT'},'both actual end routes')
        for v in end_routes.values():
            p=v['map']['file'];require(reader.retain(p['logical'],p['sha256'])==p,'accepted end-map bytes')
        old,old_route=native_rows[label]
        desired=[]
        for si,prior_jet in old['jets'].items():
            source=old['bound']['sources'][0,si];expr=values['actual']['source',si]
            candidates=[item for item in jet_atlas if same((expr,source['amplitudeUnit'],source['integralUnit'],prior_jet['probe']),
                        (item['jet']['originalBoundAmplitude'],item['jet']['amplitudeUnit'],item['jet']['integralUnit'],item['jet']['probe']))]
            if candidates:require(all(same(candidates[0]['jet'],v['jet']) for v in candidates),'complete matching saved jet returns agree')
            request={'case':label,'sourceIndex':si,'boundAmplitude':{'packet':checks['artifacts']['cases/'+label+'/scalar-bindings.pickle'],
                     'keys':['actual',['source',si]],'keyEncoding':'source key is the exact tuple (source,index)'},
                     'nativeSource':{'packet':old_route,'keys':['bound','sources',[0,si]],'keyEncoding':'source key is the exact tuple (0,index)'},
                     'ownProbeAndCoordinate':{'packet':old_route,'keys':['jets',si,'probe'],'coordinateKey':'args[0]'},
                     'matches':[v['owner'] for v in candidates], 'status':'SAVED_COMPLETE_JET' if candidates else 'NO_SAVED_COMPLETE_JET_IN_INSPECTED_ATLAS',
                     'sourceBody':reader.retain(sourcefile),'newExtractionAuthorized':False}
            desired.append(request);jet_requests.append(request)
        counts.update(sourceJets=len(desired),savedJetUses=sum(bool(v['matches']) for v in desired),unmatchedJetUses=sum(not v['matches'] for v in desired))
        summaries[label]=counts
        metadata('cases/'+label+'/validated-scalars-and-end-routes.json',{'summary':counts,'endRoutes':end_routes})
        metadata('cases/'+label+'/source-jet-input-routes.json',desired)
    require(total_new_uses==checks['newScalarCalls']==63 and sum(v['records'] for v in summaries.values())==1467,'actual new/saved census')
    controls=artifact('scalar-route-mutation-operands.pickle');responses=artifact('scalar-route-controls.json')
    for name in ('coefficient','unit','frequency'):
        changed=controls[name];require(responses[name] and not same(changed,controls['original']),'saved changed routing input')
        modified={'coefficient':'expression','unit':'unit','frequency':'mapping'}[name]
        require(all(same(v,changed[k]) for k,v in controls['original'].items() if k!=modified),'unchanged routing control fields')
    metadata('validated-routing-controls.json',responses)
    metadata('source-jet-input-catalogue.json',{'availableJetPackets':len(jet_atlas),'cases':summaries,'requests':jet_requests,
             'scope':'Full saved source-jet inputs only. An unmatched whole return does not establish that individual derivative/Poly/extraction operations never ran.'})
    reader.postcheck();metadata('validated-paths.json',reader.routes)
    result={'status':'PASSED_BOUNDED_SCALAR_BINDINGS_REVIEW_AND_JET_INPUTS','producerChecksSha256':SHA,'cases':summaries,
            'newScalarCallsReviewed':63,'physicalScalarRecords':1467,'endRoutes':8,'newScientificCalls':0,'allConsumedHashesUnchanged':True,
            'consumedLogicalPaths':len(reader.routes),'consumedPhysicalFiles':len(reader.physical),'wallSeconds':time.monotonic()-start,
            'scope':'Saved scalar bindings and existing source-jet routing at1-0.01i. No jet extraction, quadrature, matrix, solve, pole or full response acceptance.'}
    metadata('checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
