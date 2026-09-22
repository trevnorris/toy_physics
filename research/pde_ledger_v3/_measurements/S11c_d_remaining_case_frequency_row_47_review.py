#!/usr/bin/env python3
"""Read actual saved paired-row arguments and results; no numerical replay."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time

import S11c_d_remaining_case_frequency_row_47 as producer
p,pilot,np,sp,io=producer.p,producer.pilot,producer.np,producer.sp,producer.io
M,F,require,same=producer.M,producer.F,producer.require,producer.same
ROOT=F/'row-47'
CHECKS_SHA='f8e31e662d0f812047bb23fc50da9594853b967e5d5914aa2847e21259afb156'
HELPER_SHA='fc4a90c73245caaa5c55d78d2dd059c9e73a234a8fddbcbeab25f1eaee399b74'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=pilot.digest
    reader=pilot.Reader();packets={}
    def packet(rec):
        path=rec.get('path',rec.get('logical'))
        if path not in packets:packets[path]=reader.packet(path,rec['sha256'])
        return packets[path]
    def address(route):
        value=packet(route['packet'])
        for key in route['keys']:value=value[key]
        return value
    def save(name,value):io.save(base,name,value)
    def finite(v,shape):require((isinstance(v,np.ndarray) or (shape==() and isinstance(v,np.complexfloating))) and v.shape==shape and v.dtype==np.dtype(complex) and np.isfinite(v).all(),'finite actual saved array/scalar')
    inspection=reader.json(ROOT/'completion-inspection.json');require(inspection['actualExits']==[0,0,0] and inspection['emptyStrictStderr'] and inspection['checksStdoutIdentity'] and inspection['zeroCapOOMSwap'],'clean final producer')
    for rec in inspection['evidence'].values():reader.retain(rec['path'],rec['sha256'])
    outcomes=[reader.json(ROOT/n) for n in ('resource-guard/outcome.json','resource-guard/child-outcome.json','frequency_row_47.invocation.json')]
    require(all(v['exitCode']==0 for v in outcomes) and outcomes[1]['guardReason'] is None,'actual final producer exits')
    checks=reader.json(ROOT/'complete/checks.json',CHECKS_SHA);art=checks['artifacts']
    for rec in art.values():require(reader.retain(rec['path'],rec['sha256'])['bytes']==rec['bytes'],'entire actual new artifacts')
    def saved(name):
        rec=art[name];return reader.json(rec['path'],rec['sha256']) if name.endswith('.json') else packet(rec)
    for logical,rec in saved('inputs.json')['consumedRoutes'].items():require(reader.retain(logical,rec['sha256'])==rec,'all consumed native/current/frozen/logical/canonical inputs')
    reader.retain(Path(producer.__file__),HELPER_SHA)
    for path in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_row_47_review_plan.md',producer.PLAN,M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    callers=saved('native-and-new-numerical-callers.json');oldcaller=reader.json(callers['acceptedCallers']['path'],callers['acceptedCallers']['sha256'])
    for item in oldcaller['native'].values():
        r=reader.retain(item['file']['logical'],item['file']['sha256']);tree=ast.parse(Path(r['canonical']).read_text())
        for name,body in item['bodies'].items():require(ast.dump(next(n for n in tree.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'whole native caller body')
    module,join=pilot.recovery.adapted(Path(p.__file__).read_text());compile(module,'<unexecuted-native-adapter>','exec')
    require(join==callers['evaluatorJoin']==oldcaller['evaluatorJoin'],'actual entire evaluator adapter join')
    def forbidden(*args,**kwargs):raise RuntimeError('saved paired row review forbids science')
    producer.main=pilot.main=pilot.partition=pilot.trapezoid=pilot.Journal.write=p.evaluate=p.source_matrix=p.independent_source=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    for name in ('exp','tanh','sin','cos','linspace'):setattr(np,name,forbidden)
    cp=reader.json(producer.CP,producer.CP_SHA);oldchecks=reader.json(Path(cp['runDirectory'])/'checks.json',cp['checksSha256']);oldart=oldchecks['artifacts']
    def old(name):
        rec=oldart[name];return reader.json(rec['path'],rec['sha256']) if name.endswith('.json') else packet(rec)
    selected=next(v['route'] for v in reader.json(cp['nextSavedInputs']['path'],cp['nextSavedInputs']['sha256']) if v['route']['rowIndex']==47)
    raw=packet(selected['sourceInputs']['packet']);scalars=packet(selected['scalarInputs']);system=packet(selected['basis'])
    row=saved('row-input.pickle');source=saved('source-basis-route.json');pair=saved('literal-factor-input-pair.pickle');part=pair['requested'];oldpart=old('literal-factor-partition.pickle')
    physical=reader.json(selected['ownPhysicalRoutes']['packet']['logical'],selected['ownPhysicalRoutes']['packet']['sha256']);common=packet(physical['context']);context=common['contextPair'][0];jet=address(selected['jets']['23'])
    require(same(row['row'],raw['bound']['rows'][47]) and same(row['coefficient'],scalars['actual']['factor',47,0]) and same(row['source'],raw['bound']['sources'][0,23]) and same(row['jet'],jet),'actual whole requested row/source/scalar')
    require(same(row['context'],context) and same(*common['contextPair']) and same(*common['basisPair']) and same(row['positions'],system['nodes']) and
        same(row['physicalRoutes'],physical) and same(row['frequency'],scalars['frequency']),'full own context/basis/units/frequency')
    for name in ('fieldUnits','equationUnits','settings'):require(same(row[name],raw[name]),'field-equation/settings inputs')
    for name in ('abel','pairs','profileUnits'):require(same(row[name],raw['bound'][name]),'profile/measure/ordered-limit inputs')
    require(same(jet['originalBoundAmplitude'],scalars['actual']['source',23]) and same(jet['amplitudeUnit'],row['source']['amplitudeUnit']) and same(jet['integralUnit'],row['source']['integralUnit']),'actual source coefficients and units')
    action=address(source['chosen']);nodes=address(source['nodeRoute']);weights=address(source['weightRoute']);finite(action,(1024,129))
    key=(jet['probe'],tuple(jet['coefficients']),jet['amplitudeUnit'],jet['integralUnit'])
    for candidate in source['matches']:
        j=address(candidate['input']);require(same(key,(j['probe'],tuple(j['coefficients']),j['amplitudeUnit'],j['integralUnit'])) and same(address(candidate['result']),action) and
            same(address(candidate['nodeRoute']),nodes) and same(address(candidate['weightRoute']),weights),'full accepted source-basis input/value/measure match')
    original=packet(cp['rowInput']);require(same(row['context'],original['context']) and same(row['row']['limits'],original['row']['limits']) and same(row['positions'],original['positions']) and same(pair['accepted'],oldpart),'entire accepted common context and ordered limits')
    groups=part['groups'];k,q=(lim[0] for lim in row['row']['limits']);z=context['z'];integral=part['integral']
    require(same(part['coefficient'],row['coefficient']) and same(part['integral'],oldpart['profile']) and same(part['profileUnit'],row['profileUnits'][integral]) and same(part['coefficientUnit'],row['row']['factors'][0]['unit']),'same actual profile and own coefficient unit')
    factors=[v for group in groups.values() for v in group];require(len(factors)==len(row['coefficient'].args) and set(factors)==set(row['coefficient'].args),'complete literal factors without rewriting')
    save('validated-row-source-profile.json',{'rowInput':art['row-input.pickle'],'sourceAction':source['chosen'],'sameActualAcceptedProfile':True,'newSourceProfileRuleCalls':0})
    catalogue=saved('operation-catalogue.json');slots={};new_batches=set();accepted_slots=0
    for index,entry in enumerate(catalogue['factorSlots']):
        group,fi,coordinate=entry['group'],entry['factorIndex'],entry['coordinate'];arg=packet(entry['input']);value=address(entry['value'])
        require(same(arg['expression'],groups[group][fi]) and same(arg['coefficientUnit'],part['coefficientUnit']) and arg['group']==group,'full actual factor expression and physical unit')
        env=arg['environment'];indices=entry['value']['keys']
        if group=='constant':require(coordinate is None and env=={} and indices==[],'actual constant argument')
        else:
            idx=indices[0];require(len(indices)==1 and float(env[q][idx] if group=='input' else env[k][idx,0])==coordinate,'actual full factor coordinate slot')
            if group=='output':require(same(env[z],row['positions'][None,:]),'actual full output position argument')
        key=(group,fi,coordinate);require(key not in slots,'unique actual factor-call slot');slots[key]=entry['value']
        if entry['disposition']=='NEW':
            path=Path(entry['input']['path']);folder=str(path.parent.relative_to(ROOT/'complete'));receipt=saved(folder+'/completed.json')
            require(receipt=={'input':entry['input'],'value':entry['value']['packet']} and arg['ownInput']==art['row-input.pickle'],'immediate actual new factor receipt');new_batches.add(folder)
        else:
            require(entry['disposition']=='ACCEPTED_VALUE','explicit accepted factor route');path=Path(entry['input']['path']);folder=str(path.parent.relative_to(Path(cp['runDirectory'])))
            require(old(folder+'/completed.json')=={'input':entry['input'],'value':entry['value']['packet']},'actual accepted factor receipt');accepted_slots+=1
        require(np.isfinite(value).all(),'finite actual scalar/array factor slot')
        if index%256==255 or index==len(catalogue['factorSlots'])-1:save('validated-factors/'+str(index)+'.json',{'throughInclusive':index,'newBatches':len(new_batches),'acceptedSlots':accepted_slots})
    require(len(new_batches)==checks['newFactorBatches'],'actual factor batch count')
    grids={};products={}
    for label in ('512','1024'):
        rule=old('momentum-grids/'+label+'/rule/value.pickle');grid=rule['nodes'];grids[label]=rule
        for group in ('constant','input','output'):
            prefix='products/'+label+'/'+group;arg=saved(prefix+'/input.pickle');receipt=saved(prefix+'/completed.json');value=packet(receipt['value'])
            wanted=[None] if group=='constant' else [float(v) for v in grid]
            require(arg['factors']==[[slots[group,fi,v] for v in wanted] for fi in range(len(groups[group]))] and
                same(arg['coefficientUnit'],part['coefficientUnit']) and receipt['input']==art[prefix+'/input.pickle'],'all exact factor product operands')
            finite(value,() if group=='constant' else ((len(grid),129) if group=='output' else (len(grid),)))
            if receipt['disposition']=='NEW_PRODUCT':require(receipt['value']==art[prefix+'/value.pickle'],'actual new product return')
            elif receipt['disposition']=='COMPLETED_NEW_PRODUCT':require(group=='constant' and receipt['value']==products['512',group],'actual preceding new constant product')
            else:
                oldprefix='factor-products/'+('512' if group=='constant' else label)+'/'+group
                require(receipt['disposition']=='ACCEPTED_PRODUCT' and receipt['value']==oldart[oldprefix+'/value.pickle'] and arg['factors']==old(oldprefix+'/input.pickle')['factors'],'actual saved full product')
            products[label,group]=receipt['value']
    prior=reader.json(cp['upstreamCheckpoint']['logical'],cp['upstreamCheckpoint']['sha256']);prior_checks=reader.json(Path(prior['runDirectory'])/'checks.json',prior['checksSha256']);prior_art=prior_checks['artifacts']
    catref=prior['savedReview']['pointCatalogue'];cat=reader.json(catref['path'],catref['sha256']);oldpoints={v['momentumHex']:v for v in cat if v['rowIndex']==30}
    fourier={};new=old_count=0
    for index,entry in enumerate(catalogue['fourier']):
        arg=packet(entry['input']);v=packet(entry['value']);h=entry['frequencyHex'];finite(v,(129,))
        require(h not in fourier and type(arg['frequency']) is complex and arg['frequency'].imag==0 and arg['frequency'].real.hex()==h and
            same(arg['sourceAction'],source['chosen']) and same(arg['sourceNodes'],source['nodeRoute']) and same(arg['unit'],jet['integralUnit']) and
            arg['method']=='literal minus-phase vector times accepted weighted source action','full actual Fourier arguments')
        if entry['disposition']=='NEW':
            prefix='fourier/'+str(new);require(saved(prefix+'/completed.json')=={'input':entry['input'],'value':entry['value']},'actual new Fourier receipt');new+=1
        else:
            require(entry['disposition']=='ACCEPTED' and h in oldpoints,'actual accepted Fourier owner');owner=oldpoints[h];prefix='rows/30/points/'+str(owner['pointIndex'])
            r=prior_art[prefix+'/fourier-completed.json'];require(entry['input']==prior_art[prefix+'/fourier-input.pickle'] and entry['value']==owner['operands']['sourceFourier'] and
                reader.json(r['path'],r['sha256'])=={'input':entry['input'],'value':entry['value']},'actual original full Fourier input/return/receipt');old_count+=1
        fourier[h]=entry['value']
        if index%128==127 or index==len(catalogue['fourier'])-1:save('validated-fourier/'+str(index)+'.json',{'throughInclusive':index,'new':new,'accepted':old_count})
    require(new==checks['newFourierCalls'] and old_count==checks['acceptedFourierCalls'],'actual Fourier counts')
    for expected in checks['rows']:
        label=str(expected['panels']);prefix='grids/'+label;grid=grids[label]['nodes'];n=len(grid);ki=saved(prefix+'/kernel-input.pickle');kv=saved(prefix+'/kernel-value.pickle');finite(kv,(n,n))
        oldki=old('momentum-grids/'+label+'/kernel-input.pickle')
        require(ki['rule']==oldart['momentum-grids/'+label+'/rule/value.pickle'] and ki['profileRoutes']==oldki['profileValues'] and ki['inputProduct']==products[label,'input'] and
            ki['constantProduct']==products[label,'constant'] and ki['ownInput']==art['row-input.pickle'] and ki['bothMeasuresOnce'] and
            saved(prefix+'/kernel-completed.json')=={'input':art[prefix+'/kernel-input.pickle'],'value':art[prefix+'/kernel-value.pickle']},'whole same saved profile/rule and own new kernel operands')
        ii=saved(prefix+'/inner-input.pickle');iv=saved(prefix+'/inner-value.pickle');finite(iv,(n,129))
        require(ii=={'kernel':art[prefix+'/kernel-value.pickle'],'sourceValues':[fourier[float(v).hex()] for v in grid],'sourceAction':source['chosen']} and
            saved(prefix+'/inner-completed.json')=={'input':art[prefix+'/inner-input.pickle'],'value':art[prefix+'/inner-value.pickle']},'actual full inner contraction operands and value')
        ri=saved(prefix+'/row-input.pickle');rv=saved(prefix+'/row-value.pickle');finite(rv,(129,129))
        require(ri=={'output':products[label,'output'],'inner':art[prefix+'/inner-value.pickle'],'operation':'ordinary transpose then matmul','rowInput':art['row-input.pickle']} and
            saved(prefix+'/row-completed.json')=={'input':art[prefix+'/row-input.pickle'],'value':art[prefix+'/row-value.pickle']} and saved(prefix+'/summary.json')==expected,'whole new final row and receipt')
        save('validated-grids/'+label+'.json',expected)
    control=saved('selected-contraction/completed.json');finite(saved('selected-contraction/value.pickle'),(3,3));require(control['scaledDifference']<2e-12,'actual saved selected contraction check')
    comparison=saved('comparison/value.json');finite(saved('comparison/difference.pickle'),(129,129));require(comparison==checks['comparison'],'saved full row comparison')
    save('validated-comparison.json',{'comparison':comparison,'selectedContraction':control,'noNumericalRecomputation':True})
    reader.postcheck();save('validated-paths.json',reader.routes)
    result={'status':'PASSED_BOUNDED_SAVED_ROW_47_REVIEW','producerChecksSha256':CHECKS_SHA,'case':producer.CASE,'rowIndex':47,'sourceIndex':23,
        'rows':checks['rows'],'comparison':comparison,'newFactorBatches':len(new_batches),'acceptedFactorSlots':accepted_slots,'newFourierReceipts':new,
        'acceptedFourierReturns':old_count,'newScientificCalls':0,'allConsumedHashesUnchanged':True,'consumedLogicalPaths':len(reader.routes),'wallSeconds':time.monotonic()-started,
        'scope':'One own finite2D row47 saved input/value/native/source/unit/receipt review, no numerical or old science replay. Fixed-row refinement only; no scattering accuracy or pole/domain acceptance.'}
    save('checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
