#!/usr/bin/env python3
"""Saved four-row numerical receipt review; no coefficient or row recomputation."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_remainder_rows as producer

saved,p,sp,np,io=producer.saved,producer.p,producer.sp,producer.np,producer.io
M,F,require,same=producer.M,producer.F,producer.require,producer.same
ROOT=F/'remainder-rows';CHECKS_SHA='7b1166e22be054dac47ce59804bf376c77d5963d4f7b7b6a97bfabf327041618'
HELPER_SHA='79704fbe36e5d2acbc62f218b655b6958e7a3d91b152113991d51107a9d93a84'
NAME='S11c_d_remaining_case_frequency_remainder_rows_review'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal=saved.Reader(),producer.prep.inputs.MetadataJournal(base)
    c=reader.json(ROOT/'complete/checks.json',CHECKS_SHA);arts=c['artifacts'];require(c['status']=='COMPLETED_BOUNDED_FOUR_REMAINDER_ROWS','actual complete row batch')
    g=reader.json(ROOT/'resource-guard/outcome.json');ch=reader.json(ROOT/'resource-guard/child-outcome.json');iv=reader.json(ROOT/'frequency_remainder_rows.invocation.json');limits=reader.json(ROOT/'resource-guard/effective-limits.json')
    require(g['exitCode']==ch['exitCode']==iv['exitCode']==0 and ch['guardReason'] is None and g['limitsVerified'],'actual final producer outcomes')
    require(limits['memory.max']=='2147483648' and limits['memory.swap.max']=='0' and limits['pids.max']=='32' and limits['nice']==15 and len(limits['affinity'])==1 and set(limits['threads'].values())=={'1'},'mandatory containment')
    for n in ('frequency_remainder_rows.stderr','guard.stderr','resource-guard/stderr'):require(reader.retain(ROOT/n)['bytes']==0,'empty producer strict stderr')
    require(reader.retain(ROOT/'frequency_remainder_rows.stdout')['sha256']==CHECKS_SHA,'actual checks/stdout identity')
    samples=ROOT/'resource-guard/resource-samples.jsonl';reader.retain(samples)
    for line in samples.read_text().splitlines():
        s=json.loads(line);e=dict(x.split() for x in s['memory.events'].splitlines());require(int(s['memory.swap.current'])==0 and int(s['memory.peak'])<=2147483648 and all(int(e[k])==0 for k in ('max','oom','oom_kill')),'actual resource events')
    reader.retain(producer.__file__,HELPER_SHA)
    def packet(rec):return reader.packet(rec.get('path',rec.get('logical')),rec['sha256'])
    def meta(name):return reader.json(arts[name]['path'],arts[name]['sha256'])
    def value(name):return packet(arts[name])
    def array(v,shape,dtype=complex):require(isinstance(v,np.ndarray) and v.dtype==np.dtype(dtype) and v.shape==shape and np.isfinite(v).all(),'full finite saved array/type/shape')
    inputs=meta('inputs.json')
    for logical,rec in inputs['consumedRoutes'].items():reader.retain(logical,rec['sha256'])
    cp=reader.json(inputs['acceptedPreparation']['logical'],inputs['acceptedPreparation']['sha256']);require(cp['status']=='ACCEPTED_CASE_FREQUENCY_REMAINDER_NUMERICAL_PREPARATION','accepted full source/profile inputs')
    owned={v['rowIndex']:packet(v['ownInput']) for v in cp['sourceActions']};owner=owned[51];unit=owner['row']['factors'][0]['unit'];context=owner['context'];k,q=(lim[0] for lim in owner['row']['limits']);z,reg=context['z'],context['regulator']
    rows={};positions=None
    for ri in (51,52,53,54):
        own=value('rows/'+str(ri)+'/input.pickle');source=next(v for v in cp['sourceActions'] if v['rowIndex']==ri)
        require(same(own['own'],owned[ri]) and own['sourceInput']==source['ownInput'] and own['sourceAction']==source['actionValue'] and own['sharedCoefficientOwner']==cp['sourceActions'][0]['ownInput'],'full own source action/row input')
        require(same(owned[ri]['coefficient'],owner['coefficient']) and same(owned[ri]['row']['factors'][0]['unit'],unit) and same(owned[ri]['context'],context) and same(owned[ri]['fieldUnits'],owner['fieldUnits']) and same(owned[ri]['equationUnits'],owner['equationUnits']) and same(owned[ri]['row']['limits'],owner['row']['limits']) and same(owned[ri]['settings'],owner['settings']),'full coefficient unit/context equality with own independent source')
        positions=own['positions'] if positions is None else positions;require(same(positions,own['positions']),'same actual trial grid');array(packet(source['actionValue']),(1024,129));rows[ri]=own
    partition=value('literal-factor-partition.pickle');require(same(partition['expression'],owner['coefficient']) and same(partition['coefficientUnit'],unit) and partition['ownSourceInputs']==[v['ownInput'] for v in cp['sourceActions']],'actual literal factor partition input')
    flattened=[v for group in partition['groups'].values() for v in group];require(len(flattened)==len(owner['coefficient'].args) and set(flattened)==set(owner['coefficient'].args),'all literal factors retained without reconstruction')
    profile=packet(cp['profileInput']);integral=profile['integral'];pcat=reader.json(cp['profileFineCatalogue']['path'],cp['profileFineCatalogue']['sha256'])
    profile_blocks={};profile_values=[]
    for route in pcat['values']:
        key=route['packet']['path']
        if key not in profile_blocks:profile_blocks[key]=packet(route['packet'])
        v=profile_blocks[key]
        for key in route['keys']:v=v[key]
        profile_values.append(v)
    profile_values=np.asarray(profile_values)
    require(same(profile['originalCoefficient'],owner['coefficient']) and same(profile['context'],context),'full finite profile and retained Abel coefficient')
    fine_rule=value('grids/2048/rule-value.pickle');fine=fine_rule['nodes'];array(fine,(2049,),float);array(fine_rule['weights'],(2049,),float)
    # Scientific functions remain disabled while their typed pickle classes survive.
    def forbidden(*args,**kwargs):raise RuntimeError('saved row review disables scientific computation')
    producer.main=producer.partition=producer.prep.evaluator=producer.prep.main=forbidden
    p.source_matrix=p.independent_source=io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden;journal.write=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    catalogue=meta('operation-catalogue.json');require(c['acceptedFactorSlots']==0,'actual zero accepted factor-slot receipt count')
    factors=0;factor_slots={}
    for name,rec in arts.items():
        if name.startswith('factors/') and name.endswith('/input.pickle'):
            folder=name.rsplit('/',1)[0];arg=packet(rec);receipt=meta(folder+'/completed.json');group=arg['group'];env=arg['environment']
            require(arg['expression'] in partition['groups'][group] and same(arg['coefficientUnit'],unit) and arg['partition']==arts['literal-factor-partition.pickle'] and receipt=={'input':rec,'value':arts[folder+'/value.pickle']},'whole numerical factor input/units/receipt')
            shape=() if group=='constant' else (len(env[q]),) if group=='input' else (len(env[k]),129)
            if group=='output':require(same(env[z],positions[None,:]),'actual full output position argument')
            array(packet(receipt['value']),shape)
            fi=partition['groups'][group].index(arg['expression']);coordinates=[None] if group=='constant' else list(env[q] if group=='input' else env[k][:,0])
            for index,x in enumerate(coordinates):factor_slots[group,fi,None if x is None else float(x)]={'packet':receipt['value'],'keys':[] if group=='constant' else [index]}
            factors+=1;journal.json(folder+'.json',{'input':rec,'value':receipt['value'],'fullTypedInputJoined':True})
    phase_owners={};fourier_count=0;fourier_value_count=0;fourier_slots={}
    for entry in catalogue['fourierBatches']:
        arg=packet(entry['input']);folder=str(Path(entry['input']['path']).relative_to(ROOT/'complete').parent);receipt=meta(folder+'/completed.json');parg=value(folder+'/phase-input.pickle');phase_receipt=meta(folder+'/phase-completed.json')
        source=next(v for v in cp['sourceActions'] if v['sourceIndex']==entry['sourceIndex']);own=owned[source['rowIndex']];length=len(arg['frequencies'])
        require(arg['sourceAction']==source['actionValue'] and arg['sourceNodes']==cp['sourceNodes'] and arg['ownSourceInput']==source['ownInput'] and same(arg['integralUnit'],own['jet']['integralUnit']) and parg['sourceNodes']==cp['sourceNodes'] and same(parg['frequencies'],arg['frequencies']) and parg['sign']==-1,'full own Fourier and phase arguments')
        require(receipt=={'input':entry['input'],'phase':phase_receipt['value'],'value':entry['value']} and phase_receipt['input']==arts[folder+'/phase-input.pickle'],'full Fourier and phase receipt')
        if phase_receipt['disposition']=='NEW_PHASE':
            require(phase_receipt['ownerInput']==phase_receipt['input'],'actual new phase owner');phase_owners[phase_receipt['input']['path']]=(parg,phase_receipt['value'])
        else:
            require(phase_receipt['disposition']=='COMPLETED_NEW_PHASE' and phase_receipt['ownerInput']['path'] in phase_owners,'actual preceding phase owner')
            previous,returned=phase_owners[phase_receipt['ownerInput']['path']];require(same(previous,parg) and returned==phase_receipt['value'],'exact saved full phase call reuse')
        array(packet(phase_receipt['value']),(length,1024));array(packet(entry['value']),(length,129));require(entry['pointCount']==length,'actual Fourier count')
        for index,x in enumerate(arg['frequencies']):
            key=(entry['sourceIndex'],float(x));require(key not in fourier_slots,'no repeated full own Fourier point');fourier_slots[key]={'packet':entry['value'],'keys':[index]}
        journal.json(folder+'.json',{'input':entry['input'],'phase':phase_receipt,'value':entry['value'],'sourceIndex':entry['sourceIndex'],'count':length,'newScienceInReview':0});fourier_count+=1;fourier_value_count+=length
    covered=np.zeros((2049,2049),bool);mixed_count=0
    fine_mixed=value('grids/2048/mixed-value.pickle');coarse_mixed=value('grids/1024/mixed-value.pickle')
    array(fine_mixed,(2049,2049));array(coarse_mixed,(1025,1025))
    for entry in catalogue['mixedBatches']:
        arg=packet(entry['input']);folder=str(Path(entry['input']['path']).relative_to(ROOT/'complete').parent);receipt=meta(folder+'/completed.json');ri,ci=arg['masterRows'],arg['masterColumns'];shape=(len(ri),len(ci))
        require(same(arg['k'],fine[ri]) and same(arg['q'],fine[ci]) and same(arg['coefficientUnit'],unit) and arg['partition']==arts['literal-factor-partition.pickle'] and same(arg['literalFactors'],partition['groups']['mixed']) and arg['regulator']==owner['settings']['regulator'],'full actual mixed coordinates/coefficient/units/Abel arguments')
        require(np.array_equal(arg['profileIndices'],ri[:,None]-ci[None,:]+2048) and not covered[np.ix_(ri,ci)].any(),'literal metadata indices and no repeated completed mixed pair')
        pv=profile_values[arg['profileIndices']];factor_refs=[]
        for index,expr in enumerate(partition['groups']['mixed']):
            prefix=folder+'/factors/'+str(index);a=value(prefix+'/input.pickle');r=meta(prefix+'/completed.json')
            require(a['batchInput']==entry['input'] and same(a['expression'],expr) and same(a['coefficientUnit'],unit) and same(a['environment'][k],arg['k'][:,None]) and same(a['environment'][q],arg['q'][None,:]) and a['environment'][reg]==arg['regulator'] and same(a['environment'][integral],pv),'exact saved mixed factor environment and profile returns')
            require(r=={'input':arts[prefix+'/input.pickle'],'value':arts[prefix+'/value.pickle']},'actual mixed factor receipt');array(packet(r['value']),shape);factor_refs.append(r['value'])
        product=value(folder+'/product-input.pickle');require(product=={'factorValues':factor_refs,'batchInput':entry['input']} and receipt=={'input':entry['input'],'productInput':arts[folder+'/product-input.pickle'],'value':entry['value']},'full mixed product input and receipt')
        returned=packet(entry['value']);array(returned,shape);require(same(fine_mixed[np.ix_(ri,ci)],returned),'fine assembled array elements are actual saved batch returns')
        er=np.flatnonzero(ri%2==0);ec=np.flatnonzero(ci%2==0)
        if len(er) and len(ec):require(same(coarse_mixed[np.ix_(ri[er]//2,ci[ec]//2)],returned[np.ix_(er,ec)]),'coarse assembled array elements are actual saved batch returns')
        covered[np.ix_(ri,ci)]=True;mixed_count+=1;journal.json(folder+'.json',{'input':entry['input'],'value':entry['value'],'shape':shape,'typedValuesReadWithoutRecomputation':True})
    require(covered.all(),'every requested fine mixed pair has actual saved return')
    fine_rows=[]
    for panels in (1024,2048):
        label=str(panels);folder='grids/'+label;summary=meta(folder+'/summary.json');rule=packet(summary['rule']);array(rule['nodes'],(panels+1,),float);array(rule['weights'],(panels+1,),float)
        mi=value(folder+'/mixed-input.pickle');array(value(folder+'/mixed-value.pickle'),(panels+1,panels+1))
        require(mi['partition']==arts['literal-factor-partition.pickle'] and meta(folder+'/mixed-completed.json')=={'input':arts[folder+'/mixed-input.pickle'],'value':arts[folder+'/mixed-value.pickle']},'actual completed mixed array assembly receipt')
        kernel=value(folder+'/kernel-input.pickle');array(value(folder+'/kernel-value.pickle'),(panels+1,panels+1));require(kernel['mixed']==arts[folder+'/mixed-value.pickle'] and kernel['rule']==summary['rule'] and kernel['bothMeasuresOnce'] and meta(folder+'/kernel-completed.json')=={'input':arts[folder+'/kernel-input.pickle'],'value':arts[folder+'/kernel-value.pickle']},'whole recorded kernel input/measure receipt')
        array(packet(kernel['inputProduct']),(panels+1,));array(packet(kernel['constantProduct']),())
        for group in ('constant','input','output'):
            product_folder='products/'+('1024' if group=='constant' else label)+'/'+group;arg=value(product_folder+'/input.pickle');receipt=meta(product_folder+'/completed.json')
            wanted=[None] if group=='constant' else list(map(float,rule['nodes']))
            expected=[[factor_slots[group,fi,x] for x in wanted] for fi in range(len(partition['groups'][group]))]
            require(arg['factors']==expected and arg['group']==group and same(arg['coefficientUnit'],unit) and arg['partition']==arts['literal-factor-partition.pickle'] and receipt=={'input':arts[product_folder+'/input.pickle'],'value':arts[product_folder+'/value.pickle']},'actual saved full product factor arguments/receipts')
            if group in ('constant','input'):require(kernel[group+'Product']==receipt['value'],'kernel uses actual saved factor product')
        for ri in (51,52,53,54):
            prefix='rows/'+str(ri)+'/'+label;inner=value(prefix+'/inner-input.pickle');array(value(prefix+'/inner-value.pickle'),(panels+1,129))
            source=next(v for v in cp['sourceActions'] if v['rowIndex']==ri)
            require(inner['kernel']==arts[folder+'/kernel-value.pickle'] and inner['ownSourceInput']==source['ownInput'] and len(inner['sourceFourier'])==panels+1 and meta(prefix+'/inner-completed.json')=={'input':arts[prefix+'/inner-input.pickle'],'value':arts[prefix+'/inner-value.pickle']},'actual whole own inner contraction receipt')
            require(inner['sourceFourier']==[fourier_slots[source['sourceIndex'],float(x)] for x in rule['nodes']],'complete own source Fourier consumer sequence')
            row_input=value(prefix+'/row-input.pickle');require(row_input['outputProduct']==arts['products/'+label+'/output/value.pickle'],'actual saved output product consumer');array(packet(row_input['outputProduct']),(panels+1,129));array(value(prefix+'/row-value.pickle'),(129,129))
            require(row_input['inner']==arts[prefix+'/inner-value.pickle'] and row_input['operation']=='ordinary transpose @' and row_input['ownRowInput']==arts['rows/'+str(ri)+'/input.pickle'] and meta(prefix+'/row-completed.json')=={'input':arts[prefix+'/row-input.pickle'],'value':arts[prefix+'/row-value.pickle']},'actual full own row input/value receipt')
            row_summary=meta(prefix+'/summary.json');require(row_summary in c['rows'],'original whole row summary');journal.json(prefix+'.json',{'source':source['ownInput'],'summary':row_summary,'fullTypedInputAndFiniteResultJoined':True})
            if panels==2048:fine_rows.append(row_summary)
    for comparison in c['comparisons']:
        ri=comparison['rowIndex'];ci=packet(comparison['input']);array(packet(comparison['difference']),(129,129))
        require(ci=={'coarse':arts['rows/'+str(ri)+'/1024/row-value.pickle'],'fine':arts['rows/'+str(ri)+'/2048/row-value.pickle']} and comparison==meta('rows/'+str(ri)+'/comparison.json') and comparison['targetMet'] and comparison['absoluteDifference']<=max(1e-4,.01*comparison['fineNorm']),'actual saved row refinement result without difference recomputation')
    for index in range(3):
        prefix='coefficient-checks/'+str(index);arg=value(prefix+'/input.pickle');r=meta(prefix+'/completed.json');direct=packet(r['direct']);composed=packet(r['composed'])
        require(same(arg['expression'],owner['coefficient']) and arg['partition']==arts['literal-factor-partition.pickle'] and np.isfinite(direct) and np.isfinite(composed) and r['scaledDifference']<2e-12,'actual saved whole coefficient comparison')
    contraction=meta('selected-contraction/completed.json');array(packet(contraction['value']),(3,3));require(contraction['scaledDifference']<2e-12,'saved selected direct contraction comparison')
    control=meta('abel-control/completed.json');arg=packet(control['input']);require(arg['wholeInput']==arts['coefficient-checks/0/input.pickle'] and arg['changedCarrier'] in partition['groups']['mixed'] and control['response']>1e-14 and np.isfinite(packet(control['value'])),'actual omitted-Abel control response')
    require(factors==c['newFactorBatches']==8 and fourier_count==c['newFourierBatches']==132 and fourier_value_count==c['newFourierValues']==8196 and len(phase_owners)==c['newFourierPhaseBatches']==33 and mixed_count==c['newMixedBatches']==82 and int(covered.sum())==c['mixedPointPairs']==4198401,'actual receipt counts')
    journal.json('validated-result.json',{'fineRows':fine_rows,'comparisons':c['comparisons'],'coefficientCalls':factors,'fourierCalls':fourier_count,'phaseCalls':len(phase_owners),'mixedCalls':mixed_count,'newScientificCalls':0,'scope':c['scope']})
    for rec in arts.values():reader.retain(rec['path'],rec['sha256'])
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        reader.retain(path);dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    reader.postcheck();journal.json('validated-paths.json',reader.routes)
    result={'status':'PASSED_BOUNDED_SAVED_FOUR_REMAINDER_ROWS_REVIEW','rows':4,'fullRowReturns':8,'newScientificCalls':0,'consumedLogicalPaths':len(reader.routes),'allConsumedHashesUnchanged':True,'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started}
    journal.json('checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
