#!/usr/bin/env python3
"""Bounded saved-packet review of the single finite two-momentum row pilot."""
import argparse
import ast
import gc
import json
from pathlib import Path
import resource
import signal
import time

import S11c_d_remaining_case_frequency_row_2d_pilot as producer
p, io, np, sp = producer.p, producer.io, producer.np, producer.sp
M, F, require, same = producer.M, producer.F, producer.require, producer.same
ROOT = F/'row-2d-pilot'
CHECKS_SHA = '5c7b7e222e5b4c5863050dfef9a6621858c52516d6a48cd38a627c91fc0b31b5'
HELPER_SHA = 'bd0328d418ed191e882d2817f2a10b926bab0ea74707de365e2eeeccd9585b47'


class Reader(producer.Reader):
    def json(self, path, expected=None):
        route = self.retain(path, expected)
        with Path(route['canonical']).open('r') as stream:
            value = json.load(stream); producer.release(stream)
        return value


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    io.digest=producer.digest;reader=Reader();large={}
    def packet(record):return reader.packet(record.get('path',record.get('logical')),record['sha256'])
    def address(route):
        r=route['packet'];key=r.get('path',r.get('logical'))
        if key not in large:large[key]=packet(r)
        value=large[key]
        for index in route['keys']:value=value[index]
        return value
    def save(name,value):io.save(base,name,value)
    inspection=reader.json(ROOT/'completion-inspection.json')
    require(inspection['actualExits']==[0,0,0] and inspection['emptyStrictStderr'] and inspection['checksStdoutIdentity'] and inspection['zeroCapOOMSwap'],'actual final clean producer')
    for rec in inspection['evidence'].values():reader.retain(rec['path'],rec['sha256'])
    outcomes=[reader.json(ROOT/n) for n in ('resource-guard/outcome.json','resource-guard/child-outcome.json','frequency_row_2d_pilot.invocation.json')]
    require(all(v['exitCode']==0 for v in outcomes) and outcomes[1]['guardReason'] is None,'actual final guard/supervisor/child exits')
    checks=reader.json(ROOT/'complete/checks.json',CHECKS_SHA);artifacts=checks['artifacts']
    for rec in artifacts.values():require(reader.retain(rec['path'],rec['sha256'])['bytes']==rec['bytes'],'actual complete new packets')
    def saved(name):
        rec=artifacts[name]
        return reader.json(rec['path'],rec['sha256']) if name.endswith('.json') else packet(rec)
    inputs=saved('inputs.json')
    for logical,rec in inputs['consumedRoutes'].items():require(reader.retain(logical,rec['sha256'])==rec,'exact consumed logical/canonical/source input')
    reader.retain(Path(producer.__file__),HELPER_SHA)
    for path in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_row_2d_review_plan.md',producer.PLAN,
                 M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    caller=saved('native-and-new-numerical-callers.json')
    for item in caller['native'].values():
        r=reader.retain(item['file']['logical'],item['file']['sha256']);tree=ast.parse(Path(r['canonical']).read_text())
        for name,body in item['bodies'].items():require(ast.dump(next(n for n in tree.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'entire actual native caller source')
    module,join=producer.recovery.adapted(Path(p.__file__).read_text());compile(module,'<unexecuted-old-adapter>','exec')
    require(join==caller['evaluatorJoin'] and ast.unparse(ast.Module(body=[module.body[0]],type_ignores=[]))==caller['evaluator'],'only saved accepted numerical evaluator source')
    def forbidden(*args,**kwargs):raise RuntimeError('scientific call forbidden in saved2D row review')
    producer.main=producer.partition=producer.trapezoid=producer.Journal.write=forbidden
    p.evaluate=p.source_matrix=p.independent_source=p.integrate.quad_vec=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    for name in ('exp','tanh','sin','cos','linspace'):setattr(np,name,forbidden)
    cp=reader.json(producer.CP,producer.CP_SHA);next_view=reader.json(cp['nextSavedInput']['path'],cp['nextSavedInput']['sha256']);selected=next_view['route']
    raw=packet(selected['sourceInputs']['packet']);scalars=packet(selected['scalarInputs']);system=packet(selected['basis'])
    physical=reader.json(selected['ownPhysicalRoutes']['packet']['logical'],selected['ownPhysicalRoutes']['packet']['sha256']);common=packet(physical['context']);context=common['contextPair'][0]
    row=saved('row-input.pickle');si=row['row']['factors'][0]['sourceIndex'];jet=address(selected['jets'][str(si)])
    require(si==20 and same(row['row'],raw['bound']['rows'][46]) and same(row['coefficient'],scalars['actual']['factor',46,0]) and
            same(row['source'],raw['bound']['sources'][0,si]) and same(row['jet'],jet),'full actual row/source/scalar inputs')
    require(same(row['context'],context) and same(*common['contextPair']) and same(*common['basisPair']) and
            same(row['physicalRoutes'],physical) and same(row['positions'],system['nodes']) and same(row['frequency'],scalars['frequency']),'full own context/basis/frequency/address inputs')
    for name in ('fieldUnits','equationUnits','settings'):require(same(row[name],raw[name]),'own field-equation/settings units')
    for name in ('abel','pairs','profileUnits'):require(same(row[name],raw['bound'][name]),'actual profile/measure/ordered limits')
    require(same(jet['originalBoundAmplitude'],scalars['actual']['source',si]) and same(jet['integralUnit'],row['source']['integralUnit']) and
            same(jet['amplitudeUnit'],row['source']['amplitudeUnit']),'full actual source amplitude and units')
    sr=saved('source-basis-route.json');action=address(sr['chosen']);nodes=address(sr['nodeRoute']);weights=address(sr['weightRoute'])
    key=(jet['probe'],tuple(jet['coefficients']),jet['amplitudeUnit'],jet['integralUnit'])
    for candidate in sr['matches']:
        j=address(candidate['input']);require(same(key,(j['probe'],tuple(j['coefficients']),j['amplitudeUnit'],j['integralUnit'])) and
            same(address(candidate['result']),action) and same(address(candidate['nodeRoute']),nodes) and same(address(candidate['weightRoute']),weights),'exact full accepted source-action inputs/returns')
    require(action.shape==(1024,129) and action.dtype==np.dtype(complex) and np.isfinite(action).all() and same(row['sourceAction'],sr['chosen']),'complete actual weighted source array')
    old_row=packet(sr['actualPriorFourierOwner']);require(same(old_row['sourceAction'],row['sourceAction']) and same(old_row['context'],context) and
        same(key,(old_row['jet']['probe'],tuple(old_row['jet']['coefficients']),old_row['jet']['amplitudeUnit'],old_row['jet']['integralUnit'])),'accepted computational Fourier owner stays distinct from consumer')
    part=saved('literal-factor-partition.pickle');groups=part['groups'];k,q,z,xi=part['variables'];integral=part['profile'];difference=part['differenceSubtree']
    require(same(part['coefficient'],row['coefficient']) and same(part['profileUnit'],row['profileUnits'][integral]) and
            same(part['coefficientUnit'],row['row']['factors'][0]['unit']) and same(row['source']['frequency'],q),'entire coefficient/profile unit/source momentum joins')
    actual=[]
    for group,values in groups.items():
        actual.extend(values)
        for expression in values:
            allowed={'constant':set(),'output':{k,z},'input':{q},'profile':{k,q}}[group]
            require(expression.free_symbols<=allowed and (isinstance(expression,sp.Integral)==(group=='profile')),'literal factor group covers actual free variables')
    require(len(actual)==len(row['coefficient'].args) and set(actual)==set(row['coefficient'].args) and groups['profile']==[integral],'all original literal factors retained exactly once')
    require(isinstance(difference,sp.Add) and k in difference.args and len(difference.args)==2,'saved literal difference subtree')
    neg=next(v for v in difference.args if v!=k)
    require(isinstance(neg,sp.Mul) and set(neg.args)=={sp.S.NegativeOne,q},'actual output-minus-input sign')
    def covered(node):
        if node==difference:return
        require(node not in (k,q),'all raw profile momentum occurrences covered')
        for child in node.args:covered(child)
    covered(integral.function)
    save('validated-input-source-and-partition.json',{'rowInput':artifacts['row-input.pickle'],'sourceAction':sr['chosen'],'completeLiteralFactors':len(actual),'newSourceCalls':0,'wholeOwnSourceUnitsAndNativeCallerJoined':True})
    def finite(value,shape):
        # NumPy multiplication of zero-dimensional arrays may return its scalar
        # class. Preserve that actual saved class rather than coercing it.
        require((isinstance(value,np.ndarray) or (shape==() and isinstance(value,np.complexfloating))) and
            value.shape==shape and value.dtype==np.dtype(complex) and np.isfinite(value).all(),'actual complete finite saved array/scalar')
    rules={}
    for name,count in [('profile-rules/2048',2048),('profile-rules/4096',4096),('momentum-grids/512/rule',512),('momentum-grids/1024/rule',1024)]:
        arg=saved(name+'/input.pickle');v=saved(name+'/value.pickle');done=saved(name+'/completed.json');guard=saved(name+'/measure-check.json')
        require(done=={'input':artifacts[name+'/input.pickle'],'value':artifacts[name+'/value.pickle']} and arg['panels']==count and
            v['nodes'].shape==v['weights'].shape==(count+1,) and np.isfinite(v['nodes']).all() and np.isfinite(v['weights']).all() and
            np.all(v['weights']>0) and v['nodes'][0]==arg['lower'] and v['nodes'][-1]==arg['upper'] and guard['mass']==guard['expectedMass'],'saved rule, measure guard and literal endpoints')
        rules[artifacts[name+'/value.pickle']['path']]=v
    catalogue=saved('operation-catalogue.json');profile_calls={};profile_values={}
    for index,rec in enumerate(catalogue['profileBatches']):
        arg=packet(rec['input']);value=packet(rec['value']);folder='profiles/'+str(index);done=saved(folder+'/completed.json');rv=rules[arg['xiRule']['path']]
        require(done=={'input':rec['input'],'integrand':artifacts[folder+'/integrand-value.pickle'],'contractionInput':artifacts[folder+'/contraction-input.pickle'],'value':rec['value']} and
                same(arg['integral'],integral) and same(arg['differenceSubtree'],difference) and same(arg['profileUnit'],part['profileUnit']) and
                len(arg['deltas'])==rec['count']<=64,'full actual profile arguments and immediate receipt')
        integrand_value=saved(folder+'/integrand-value.pickle');finite(integrand_value,(len(arg['deltas']),len(rv['nodes'])));finite(value,(len(arg['deltas']),))
        contraction=saved(folder+'/contraction-input.pickle');require(contraction=={'integrand':done['integrand'],'rule':arg['xiRule'],'axis':1,'measureIncludedOnce':True},'actual weighted profile consumer operands')
        for j,d in enumerate(arg['deltas']):
            key=(rec['panels'],float(d).hex());require(key not in profile_calls,'each actual profile input invoked once')
            profile_calls[key]={'packet':rec['value'],'keys':[j]};profile_values[rec['value']['path'],j]=value[j]
        save('validated-profiles/'+str(index)+'.json',{'input':rec['input'],'value':rec['value'],'count':len(value),'panels':rec['panels'],'actualIntegrandAndWeightedReturnRead':True})
        del integrand_value
    require(len(profile_calls)==checks['uniqueProfileCalls']==2056 and len(catalogue['profileBatches'])==34,'actual profile call counts')
    profile_compare=saved('profile-comparison/value.json');pcin=saved('profile-comparison/input.json')
    require(profile_compare['absoluteSpread']==checks['profileComparisonAbsolute']<profile_compare['tolerance']==1e-9,'saved actual finite-profile comparison result')
    for panel,label in ((2048,'coarse'),(4096,'fine')):
        for d,r in zip(saved('pilot-scope.json')['profileComparisonDeltas'],pcin[label]):require(r==profile_calls[panel,float(d).hex()],'complete profile comparison operand route')
    fourier_values={};new_count=old_count=0
    oldchecks=reader.json(Path(cp['runDirectory'])/'checks.json',cp['checksSha256']);oldarts=oldchecks['artifacts']
    oldcatref=cp['savedReview']['pointCatalogue'];oldcat=reader.json(oldcatref['path'],oldcatref['sha256'])
    oldpoints={v['momentumHex']:v for v in oldcat if v['rowIndex']==29}
    for index,route in enumerate(catalogue['sourceFourier']):
        arg=packet(route['input']);value=packet(route['value']);h=route['frequencyHex']
        require(h not in fourier_values and type(arg['frequency']) is complex and arg['frequency'].imag==0 and arg['frequency'].real.hex()==h and
            same(arg['sourceAction'],sr['chosen']) and same(arg['sourceNodes'],sr['nodeRoute']) and same(arg['unit'],jet['integralUnit']) and
            arg['method']=='literal minus-phase vector times accepted weighted source action','full typed actual source Fourier input')
        finite(value,(129,))
        if route['disposition']=='NEW':
            prefix='source-fourier/'+str(new_count);require(route['input']==artifacts[prefix+'/input.pickle'] and route['value']==artifacts[prefix+'/value.pickle'] and
                saved(prefix+'/completed.json')=={'input':route['input'],'value':route['value']},'actual new Fourier operation receipt');new_count+=1
        else:
            require(route['disposition']=='ACCEPTED_SAVED' and h in oldpoints,'recorded accepted Fourier route')
            owner=oldpoints[h];prefix='rows/29/points/'+str(owner['pointIndex']);receipt_ref=oldarts[prefix+'/fourier-completed.json']
            require(route['input']==oldarts[prefix+'/fourier-input.pickle'] and route['value']==owner['operands']['sourceFourier'] and
                reader.json(receipt_ref['path'],receipt_ref['sha256'])=={'input':route['input'],'value':route['value']},'accepted actual complete input/receipt/return match');old_count+=1
        fourier_values[h]=route['value']
        if index%64==63 or index==len(catalogue['sourceFourier'])-1:save('validated-fourier/'+str(index)+'.json',{'throughInclusive':index,'newCount':new_count,'savedCount':old_count,'fullInputsValuesReceiptsRead':True})
    require(new_count==checks['newSourceFourierCalls']==842 and old_count==checks['acceptedSavedFourierCalls']==183,'actual new/saved Fourier counts')
    del oldchecks,oldarts,oldcat,oldpoints;gc.collect()
    factor_values={};factor_requests={}
    for name in artifacts:
        if not name.startswith('coefficient-factors/') or not name.endswith('/input.pickle'):continue
        folder=name.rsplit('/',1)[0];arg=saved(name);v=saved(folder+'/value.pickle');done=saved(folder+'/completed.json');group=arg['group'];fi=arg['factorIndex']
        require(same(arg['expression'],groups[group][fi]) and same(arg['coefficientUnit'],part['coefficientUnit']) and arg['partition']==artifacts['literal-factor-partition.pickle'] and
            done=={'input':artifacts[name],'value':artifacts[folder+'/value.pickle']},'actual numerical factor arguments/value/receipt')
        env=arg['environment']
        if group=='constant':require(env=={} and v.shape==(),'actual constant factor');points=[None]
        elif group=='input':require(set(env)=={q},'only actual input momentum environment');points=list(env[q]);finite(v,(len(points),))
        else:
            require(set(env)=={k,z} and same(env[z],row['positions'][None,:]) and env[k].shape[1:]==(1,),'complete actual output factor position environment')
            points=list(env[k][:,0]);finite(v,(len(points),129))
        for j,point in enumerate(points):
            h=None if point is None else float(point).hex();key=(group,fi,h);require(key not in factor_requests,'actual factor input evaluated once across nested grids')
            route={'packet':done['value'],'keys':[] if group=='constant' else [j]};factor_requests[key]=route
        factor_values[done['value']['path']]=v
    for product in catalogue['factorProducts']:
        if product.get('disposition')=='EXACT_COMPLETED_PRODUCT':
            require(product['group']=='constant' and product['label']=='1024' and product['value']==artifacts['factor-products/512/constant/value.pickle'],'constant product uses exact completed return');continue
        arg=packet(product['input']);value=packet(product['value']);group=product['group'];label=product['label'];prefix='factor-products/'+label+'/'+group
        require(saved(prefix+'/completed.json')=={'input':product['input'],'value':product['value']} and arg['group']==group and
            arg['partition']==artifacts['literal-factor-partition.pickle'],'actual full factor product consumer')
        grid=rules[artifacts['momentum-grids/'+label+'/rule/value.pickle']['path']]['nodes']
        hs=[None] if group=='constant' else [float(v).hex() for v in grid]
        require(len(arg['factors'])==len(groups[group]) and all(routes==[factor_requests[group,fi,h] for h in hs] for fi,routes in enumerate(arg['factors'])),'complete fine/coarse factor routes without reevaluation')
        finite(value,() if group=='constant' else ((len(grid),129) if group=='output' else (len(grid),)))
    save('validated-factor-calls.json',{'uniqueFullCoordinateInputs':len(factor_requests),'productReceipts':len(catalogue['factorProducts']),'coarseValuesAndConstantProductReused':True})
    summaries=[]
    for expected in checks['rows']:
        label=str(expected['panels']);folder='momentum-grids/'+label;gridref=artifacts[folder+'/rule/value.pickle'];grid=rules[gridref['path']]['nodes'];n=len(grid)
        ki=saved(folder+'/kernel-input.pickle');kv=saved(folder+'/kernel-value.pickle');kd=saved(folder+'/kernel-completed.json');finite(kv,(n,n))
        require(kd=={'input':artifacts[folder+'/kernel-input.pickle'],'value':artifacts[folder+'/kernel-value.pickle']} and ki['grid']==gridref and
            ki['differenceIndex']=='outputIndex-inputIndex+panels' and ki['bothMomentumMeasuresIncludedOnce'] and ki['rowInput']==artifacts['row-input.pickle'] and
            len(ki['deltas'])==len(ki['profileValues'])==2*expected['panels']+1,'actual complete kernel arguments and returned matrix')
        for d,route in zip(ki['deltas'],ki['profileValues']):require(route==profile_calls[4096,float(d).hex()],'actual saved profile values cover kernel difference addresses')
        sv=saved(folder+'/source-value-routes.pickle');require(sv['grid']==gridref and sv['returns']==[fourier_values[float(qv).hex()] for qv in grid],'all full source Fourier column routes')
        ii=saved(folder+'/inner-input.pickle');iv=saved(folder+'/inner-value.pickle');idone=saved(folder+'/inner-completed.json');finite(iv,(n,129))
        require(ii=={'kernel':artifacts[folder+'/kernel-value.pickle'],'sourceValues':artifacts[folder+'/source-value-routes.pickle'],'operation':'kernel @ sourceValues'} and
            idone=={'input':artifacts[folder+'/inner-input.pickle'],'value':artifacts[folder+'/inner-value.pickle']},'actual first full contraction before/after packet')
        ri=saved(folder+'/row-input.pickle');rv=saved(folder+'/row-value.pickle');rd=saved(folder+'/row-completed.json');finite(rv,(129,129))
        require(ri=={'outputFactors':artifacts['factor-products/'+label+'/output/value.pickle'],'inner':artifacts[folder+'/inner-value.pickle'],
            'operation':'outputFactors.T @ inner','transposeConjugates':False} and rd=={'input':artifacts[folder+'/row-input.pickle'],'value':artifacts[folder+'/row-value.pickle']},'actual final complete contraction with ordinary transpose')
        summary=saved(folder+'/summary.json');require(summary==expected and summary['rowReturn']==rd['value'],'actual complete row summary')
        save('validated-grids/'+label+'.json',summary);summaries.append(summary)
        del kv,iv,rv
    for index in range(3):
        prefix='coefficient-checks/'+str(index);arg=saved(prefix+'/input.pickle');comparison=saved(prefix+'/comparison-input.pickle');receipt=saved(prefix+'/completed.json')
        direct=saved(prefix+'/direct-value.pickle');grouped=saved(prefix+'/grouped-value.pickle');finite(direct,(129,));finite(grouped,(129,))
        require(same(arg['expression'],row['coefficient']) and same(arg['environment'][z],row['positions']) and
            receipt=={'input':artifacts[prefix+'/input.pickle'],'direct':artifacts[prefix+'/direct-value.pickle'],'comparisonInput':artifacts[prefix+'/comparison-input.pickle'],
            'grouped':artifacts[prefix+'/grouped-value.pickle'],'scaledDifference':receipt['scaledDifference']} and receipt['scaledDifference']<2e-12,'saved full coefficient evaluation and actual comparison receipt')
    control=saved('contraction-check/completed.json');direct=saved('contraction-check/direct-value.pickle');mutations=saved('contraction-check/mutation-values.pickle');finite(direct,(3,3))
    for v in mutations.values():finite(v,(3,3));require(not same(v,direct),'actual changed contraction operand responds')
    require(control['scaledDifference']<2e-12 and all(v>1e-14 for v in control['responses'].values()),'saved selected literal-sum and measure/index control results')
    comparison=saved('row-comparison/value.json');delta=saved('row-comparison/difference.pickle');finite(delta,(129,129))
    require(comparison==checks['comparison'] and saved('row-comparison/input.json')=={'coarse':summaries[0]['rowReturn'],'fine':summaries[1]['rowReturn']},'actual saved final two-grid comparison')
    save('validated-comparisons-and-controls.json',{'profile':profile_compare,'selectedContraction':control,'rowComparison':comparison,'noComparisonArithmeticRepeated':True})
    # Only a small metadata view of the five remaining LAB_RHOBR2D rows, already
    # present in consumed source/scalar packets; no new profile or row evaluation.
    inventory=reader.json(p.READY/'completed-input-artifact-inventory.json',p.INVENTORY_SHA);next_rows=[]
    for ri in (47,51,52,53,54):
        name='cases/'+producer.CASE+'/row-'+str(ri)+'.json';route=reader.json(p.READY/'complete'/name,inventory[name]['sha256']);r=raw['bound']['rows'][ri]
        entries=[]
        for fi,factor in enumerate(r['factors']):
            expression=scalars['actual']['factor',ri,fi]
            entries.append({'sourceIndex':factor['sourceIndex'],'sourceFrequency':str(raw['bound']['sources'][0,factor['sourceIndex']]['frequency']),
                'coefficient':str(expression),'integrals':[{'expression':str(v),'unit':repr(raw['bound']['profileUnits'].get(v))} for v in sorted(expression.atoms(sp.Integral),key=str)]})
        next_rows.append({'route':route,'coefficientOperands':entries})
    save('next-case-2d-inputs.json',next_rows)
    reader.postcheck();save('validated-paths.json',reader.routes)
    result={'status':'PASSED_BOUNDED_SAVED_2D_ROW_REVIEW','producerChecksSha256':CHECKS_SHA,'case':producer.CASE,'rowIndex':46,'sourceIndex':20,
        'rows':summaries,'comparison':comparison,'profileComparison':profile_compare,'newFourierReceipts':new_count,'acceptedFourierReturns':old_count,
        'profileCalls':len(profile_calls),'profileBatches':len(catalogue['profileBatches']),'newScientificCalls':0,'allConsumedHashesUnchanged':True,
        'consumedLogicalPaths':len(reader.routes),'wallSeconds':time.monotonic()-started,
        'scope':'Saved full actual input/intermediate/value/receipt/source/unit review of one finite2D row. No profile/evaluator/source/Fourier/factor/product/contraction/rule/comparison or old science replay. Fixed-row stability only, not scattering accuracy or pole/domain certification.'}
    save('checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
