#!/usr/bin/env python3
"""Finish the paired finite2D row using accepted profile/rule/factor returns."""
import argparse
import ast
import gc
import json
from pathlib import Path
import resource
import signal
import time

import S11c_d_remaining_case_frequency_row_2d_pilot as pilot
p,np,sp,io=pilot.p,pilot.np,pilot.sp,pilot.io
M,F,require,same=pilot.M,pilot.F,pilot.require,pilot.same
CP=M/'S11c_d_remaining_case_frequency_row_2d_checkpoint.json'
CP_SHA='a4e1a9cc50f126a40730dac2d5f8886ad527f8d5b667be5ee3be83ccd8706229'
PLAN=M/'S11c_d_remaining_case_frequency_row_47_plan.md'
CASE=pilot.CASE


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    io.digest=pilot.digest;reader,journal=pilot.Reader(),pilot.Journal(base);packets={}
    def packet(rec):
        path=rec.get('path',rec.get('logical'))
        if path not in packets:packets[path]=reader.packet(path,rec['sha256'])
        return packets[path]
    def address(route):
        v=packet(route['packet'])
        for k in route['keys']:v=v[k]
        return v
    cp=reader.json(CP,CP_SHA);require(cp['status']=='ACCEPTED_BOUNDED_CASE_FREQUENCY_2D_ROW_PILOT','accepted paired2D pilot')
    oldchecks=reader.json(Path(cp['runDirectory'])/'checks.json',cp['checksSha256']);oldarts=oldchecks['artifacts']
    def old(name):
        rec=oldarts[name]
        return reader.json(rec['path'],rec['sha256']) if name.endswith('.json') else packet(rec)
    next_rows=reader.json(cp['nextSavedInputs']['path'],cp['nextSavedInputs']['sha256'])
    selected=next(v['route'] for v in next_rows if v['route']['rowIndex']==47)
    require(selected['layoutDimension']==2 and selected['firstOwner']==[CASE,47] and not selected['savedCompleteMatches'],'actual unmatched whole row47 input')
    own=packet(selected['sourceInputs']['packet']);scalars=packet(selected['scalarInputs']);system=packet(selected['basis'])
    physical=reader.json(selected['ownPhysicalRoutes']['packet']['logical'],selected['ownPhysicalRoutes']['packet']['sha256']);common=packet(physical['context']);context=common['contextPair'][0]
    original=packet(cp['rowInput']);row=own['bound']['rows'][47];require(row['index']==47 and len(row['factors'])==1,'whole single factor paired row')
    factor=row['factors'][0];si=factor['sourceIndex'];source=own['bound']['sources'][0,si];jet=address(selected['jets'][str(si)])
    coefficient=scalars['actual']['factor',47,0];positions=system['nodes'];settings=own['settings'];k,q=(v[0] for v in row['limits']);z,xi=context['z'],context['xi']
    require(si==23 and same(source['frequency'],q) and same(jet['originalBoundAmplitude'],scalars['actual']['source',si]) and
        same(jet['integralUnit'],source['integralUnit']) and same(jet['amplitudeUnit'],source['amplitudeUnit']),'own full source and units')
    require(same(*common['contextPair']) and same(*common['basisPair']) and same(context,original['context']) and
        same(settings,original['settings']) and same(settings,system['settings']) and same(positions,original['positions']) and
        same(own['fieldUnits'],original['fieldUnits']) and same(own['equationUnits'],original['equationUnits']) and
        same(physical,original['physicalRoutes']) and same(scalars['frequency'],original['frequency']) and
        same(row['limits'],original['row']['limits']),'whole same physical caller and ordered limits while preserving own source')
    inventory=reader.json(p.READY/'completed-input-artifact-inventory.json',p.INVENTORY_SHA)
    def metadata(name):return reader.json(p.READY/'complete'/name,inventory[name]['sha256'])
    candidates=[v for v in next(v for v in metadata('cases/'+CASE+'/source-basis-inputs.json') if v['sourceIndex']==si)['savedCoefficientBasisCandidates'] if v['owner'][0]==CASE]
    require(candidates,'actual whole source-basis candidates')
    key=(jet['probe'],tuple(jet['coefficients']),jet['amplitudeUnit'],jet['integralUnit']);action_route=candidates[0]['result'];node_route=candidates[0]['nodeRoute'];weight_route=candidates[0]['weightRoute']
    action=address(action_route);nodes=address(node_route);weights=address(weight_route)
    for candidate in candidates:
        j=address(candidate['input']);require(same(key,(j['probe'],tuple(j['coefficients']),j['amplitudeUnit'],j['integralUnit'])) and
            same(address(candidate['result']),action) and same(address(candidate['nodeRoute']),nodes) and same(address(candidate['weightRoute']),weights),'full actual source-basis input/value matches')
    require(action.shape==(1024,129) and action.dtype==np.dtype(complex) and np.isfinite(action).all() and
        same(nodes,address(old('source-basis-route.json')['nodeRoute'])) and same(weights,address(old('source-basis-route.json')['weightRoute'])),'exact whole source rule and saved finite array')
    input_record=journal.write('row-input.pickle',{'row':row,'coefficient':coefficient,'source':source,'jet':jet,'context':context,'settings':settings,
        'positions':positions,'physicalRoutes':physical,'scalarInputs':selected['scalarInputs'],'frequency':scalars['frequency'],'sourceAction':action_route,
        'fieldUnits':own['fieldUnits'],'equationUnits':own['equationUnits'],'profileUnits':own['bound']['profileUnits'],'abel':own['bound']['abel'],'pairs':own['bound']['pairs']})
    journal.json('source-basis-route.json',{'matches':candidates,'chosen':action_route,'nodeRoute':node_route,'weightRoute':weight_route,'newSourceActions':0})
    # New literal-factor metadata only; all matching scientific values stay saved.
    reader.retain(Path(pilot.__file__),'bd0328d418ed191e882d2817f2a10b926bab0ea74707de365e2eeeccd9585b47')
    groups,integral,difference=pilot.partition(coefficient,k,q,z,xi);oldpart=old('literal-factor-partition.pickle')
    require(same(integral,oldpart['profile']) and same(difference,oldpart['differenceSubtree']) and
        same(own['bound']['profileUnits'][integral],oldpart['profileUnit']),'entire same accepted finite profile call and physical unit')
    part_ref=journal.write('literal-factor-input-pair.pickle',{'requested':{'coefficient':coefficient,'groups':groups,'integral':integral,'differenceSubtree':difference,
        'coefficientUnit':factor['unit'],'profileUnit':own['bound']['profileUnits'][integral]},'accepted':oldpart,'acceptedInput':cp['rowInput'],'ownInput':input_record})
    caller=old('native-and-new-numerical-callers.json')
    for item in caller['native'].values():
        r=reader.retain(item['file']['logical'],item['file']['sha256']);tree=ast.parse(Path(r['canonical']).read_text())
        for name,body in item['bodies'].items():require(ast.dump(next(n for n in tree.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'whole original native row/source caller')
    for logical,rec in old('inputs.json')['consumedRoutes'].items():
        if logical.endswith(('.py','.md')):reader.retain(logical,rec['sha256'])
    module,join=pilot.recovery.adapted(Path(p.__file__).read_text());namespace=dict(vars(p));numerical=ast.Module(body=[module.body[0]],type_ignores=[])
    exec(compile(numerical,'<accepted-numerical-evaluator>','exec'),namespace);evaluate=namespace['evaluate']
    require(join==caller['evaluatorJoin'] and ast.unparse(numerical)==caller['evaluator'],'entire unchanged principal numerical evaluator')
    for path in (Path(__file__).resolve(),PLAN,M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    for path in (Path(__file__).resolve(),PLAN):
        target=base/'source'/path.name;target.parent.mkdir(exist_ok=True)
        with target.open('xb') as out:out.write(path.read_bytes())
        reader.retain(target,pilot.digest(path))
    journal.json('native-and-new-numerical-callers.json',{'acceptedCallers':oldarts['native-and-new-numerical-callers.json'],'evaluatorJoin':join,
        'actualSourceCandidatePairs':candidates,'newHelper':reader.retain(Path(__file__).resolve()),'profileAndRulesAreSavedReturns':True})
    scope={'case':CASE,'rowIndex':47,'sourceIndex':si,'frequency':{'real':1.,'imag':-.01},'panels':[512,1024],
        'momentumBounds':[-4.,4.],'profileBounds':[-14.,14.],'chartAndSheet':'Exact same accepted principal coefficient branches; own source23 and end-map context retained.',
        'newProfileRulesOrValues':0,'newSourceActions':0,'measuredPriorFineGridSeconds':1.5236331930063898,
        'nextCallCostGate':'3x actual/prior grid cost +60 inside remaining900s','comparisonIndication':{'relative':.01,'absolute':1e-4},
        'scope':'Second paired finite2D129x129 row at1-.01i. Saved grids/profile and full matching factor values, own source/Fourier. No native rule/derivative/compiler/basis/end/science replay. Fixed-row refinement only, not scattering accuracy or pole/domain claim.'}
    journal.json('batch-scope.json',scope)
    def forbidden(*args,**kwargs):raise RuntimeError('completed native science disabled in paired2D row')
    pilot.main=pilot.trapezoid=p.source_matrix=p.independent_source=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    # Full actual saved factor input/value slots, including fine nested-grid slots.
    available={}
    for name,rec in oldarts.items():
        if name.startswith('coefficient-factors/') and name.endswith('/input.pickle'):
            folder=name.rsplit('/',1)[0];arg=old(name);value=old(folder+'/value.pickle');receipt=old(folder+'/completed.json')
            require(receipt=={'input':rec,'value':oldarts[folder+'/value.pickle']},'accepted actual factor receipt')
            group=arg['group'];env=arg['environment']
            if group=='output':require(same(env[z],positions[None,:]),'full saved output-factor position argument')
            coords=[None] if group=='constant' else list(env[q] if group=='input' else env[k][:,0])
            for i,coordinate in enumerate(coords):
                address_key=(group,None if coordinate is None else float(coordinate))
                available.setdefault(address_key,[]).append({'expression':arg['expression'],'unit':arg['coefficientUnit'],
                    'value':value if group=='constant' else value[i],'route':{'packet':receipt['value'],'keys':[] if group=='constant' else [i]},'input':rec})
    factor_cache={};factor_routes=[];product_cache={};new_factor_batches=0
    def products(grid,label):
        nonlocal new_factor_batches
        result={}
        for group in ('constant','input','output'):
            all_routes=[];all_values=[]
            for fi,expression in enumerate(groups[group]):
                wanted=[None] if group=='constant' else [float(v) for v in grid];missing=[]
                for coordinate in wanted:
                    key=(group,fi,coordinate)
                    if key in factor_cache:continue
                    matches=[v for v in available.get((group,coordinate),[]) if same(v['expression'],expression) and same(v['unit'],factor['unit'])]
                    if matches:
                        require(all(same(matches[0]['value'],v['value']) for v in matches),'all exact full factor candidate returns agree')
                        v=matches[0];factor_cache[key]=(v['value'],v['route'])
                        factor_routes.append({'group':group,'factorIndex':fi,'coordinate':coordinate,'disposition':'ACCEPTED_VALUE','input':v['input'],'value':v['route']})
                    else:missing.append(coordinate)
                if missing:
                    coords=None if group=='constant' else np.asarray(missing)
                    env={} if group=='constant' else ({q:coords} if group=='input' else {k:coords[:,None],z:positions[None,:]})
                    folder='factor-calls/'+str(new_factor_batches);new_factor_batches+=1
                    ar=journal.write(folder+'/input.pickle',{'expression':expression,'environment':env,'coefficientUnit':factor['unit'],'group':group,'factorIndex':fi,'ownInput':input_record})
                    with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
                        value=np.asarray(evaluate(expression,env),complex)
                        value=value.reshape(()) if group=='constant' else np.broadcast_to(value,(len(missing),129) if group=='output' else (len(missing),))
                    vr=journal.write(folder+'/value.pickle',value);journal.json(folder+'/completed.json',{'input':ar,'value':vr});require(np.isfinite(value).all(),'finite genuinely missing factor values')
                    for i,coordinate in enumerate(missing):
                        route={'packet':vr,'keys':[] if group=='constant' else [i]};factor_cache[group,fi,coordinate]=(value if group=='constant' else value[i],route)
                        factor_routes.append({'group':group,'factorIndex':fi,'coordinate':coordinate,'disposition':'NEW','input':ar,'value':route})
                all_routes.append([factor_cache[group,fi,v][1] for v in wanted]);all_values.append([factor_cache[group,fi,v][0] for v in wanted])
            prefix='products/'+label+'/'+group
            arg=journal.write(prefix+'/input.pickle',{'group':group,'factors':all_routes,'ownPartition':part_ref,'coefficientUnit':factor['unit']})
            old_product='factor-products/'+('512' if group=='constant' else label)+'/'+group
            old_arg=old(old_product+'/input.pickle')
            reusable=same(factor['unit'],oldpart['coefficientUnit']) and same(groups[group],oldpart['groups'][group]) and all_routes==old_arg['factors']
            if reusable:vr=oldarts[old_product+'/value.pickle'];value=packet(vr);disposition='ACCEPTED_PRODUCT'
            elif group=='constant' and group in product_cache:
                value,vr,previous_routes=product_cache[group];require(previous_routes==all_routes,'exact new constant product reuse');disposition='COMPLETED_NEW_PRODUCT'
            else:
                value=np.ones(() if group=='constant' else ((len(grid),129) if group=='output' else (len(grid),)),complex)
                for values in all_values:value=value*(values[0] if group=='constant' else np.asarray(values))
                vr=journal.write(prefix+'/value.pickle',value);disposition='NEW_PRODUCT'
                if group=='constant':product_cache[group]=(value,vr,all_routes)
            journal.json(prefix+'/completed.json',{'input':arg,'value':vr,'disposition':disposition});result[group]=(value,vr)
        return result
    # The other accepted1D row supplies exactly matching source23 Fourier calls.
    prior_cp=reader.json(cp['upstreamCheckpoint']['logical'],cp['upstreamCheckpoint']['sha256'])
    prior_checks=reader.json(Path(prior_cp['runDirectory'])/'checks.json',prior_cp['checksSha256']);prior_art=prior_checks['artifacts']
    prior_owner=next(v for v in prior_cp['rows'] if v['rowIndex']==30);prior_row=packet(prior_owner['rowInput'])
    require(same(prior_row['sourceAction'],action_route) and same(key,(prior_row['jet']['probe'],tuple(prior_row['jet']['coefficients']),prior_row['jet']['amplitudeUnit'],prior_row['jet']['integralUnit'])) and same(prior_row['context'],context),'whole accepted source30/ownsource23 Fourier caller join')
    catref=prior_cp['savedReview']['pointCatalogue'];cat=reader.json(catref['path'],catref['sha256']);oldpoints={v['momentumHex']:v for v in cat if v['rowIndex']==30}
    fourier_cache={};fourier_routes=[];new_fourier=0
    def fouriers(grid):
        nonlocal new_fourier
        for qv in map(float,grid):
            h=qv.hex()
            if h in fourier_cache:continue
            request={'frequency':complex(qv),'sourceAction':action_route,'sourceNodes':node_route,'unit':jet['integralUnit'],'method':'literal minus-phase vector times accepted weighted source action'}
            if h in oldpoints:
                owner=oldpoints[h];prefix='rows/30/points/'+str(owner['pointIndex']);ar=prior_art[prefix+'/fourier-input.pickle'];vr=owner['operands']['sourceFourier'];receipt=prior_art[prefix+'/fourier-completed.json']
                require(same(packet(ar),request) and reader.json(receipt['path'],receipt['sha256'])=={'input':ar,'value':vr},'actual complete accepted Fourier call');value=packet(vr);disposition='ACCEPTED'
            else:
                prefix='fourier/'+str(new_fourier);new_fourier+=1;ar=journal.write(prefix+'/input.pickle',request)
                value=np.exp(-1j*qv*nodes)@action;vr=journal.write(prefix+'/value.pickle',value);journal.json(prefix+'/completed.json',{'input':ar,'value':vr});disposition='NEW'
            require(value.shape==(129,) and value.dtype==np.dtype(complex) and np.isfinite(value).all(),'full actual source Fourier return')
            fourier_cache[h]=(value,vr);fourier_routes.append({'frequencyHex':h,'input':ar,'value':vr,'disposition':disposition})
        return np.asarray([fourier_cache[float(v).hex()][0] for v in grid]),[fourier_cache[float(v).hex()][1] for v in grid]
    costs=[];results=[];summaries=[]
    for panels in (512,1024):
        label=str(panels);folder='grids/'+label;oldfolder='momentum-grids/'+label
        reserve=3*max([1.5236331930063898]+costs)+60;remaining=900-(time.monotonic()-started)
        journal.json(folder+'/cost-decision.json',{'completedGridSeconds':costs,'requiredReserveSeconds':reserve,'remainingSeconds':remaining});require(reserve<remaining,'measured bounded next paired row grid fits')
        tick=time.monotonic();rule=old(oldfolder+'/rule/value.pickle');grid,mass=rule['nodes'],rule['weights'];old_kernel=old(oldfolder+'/kernel-input.pickle')
        profile=np.asarray([address(v) for v in old_kernel['profileValues']]);source_values,source_refs=fouriers(grid);parts=products(grid,label)
        a,ar=parts['output'];b,br=parts['input'];c,cr=parts['constant']
        arg=journal.write(folder+'/kernel-input.pickle',{'rule':oldarts[oldfolder+'/rule/value.pickle'],'profileRoutes':old_kernel['profileValues'],
            'profileInput':oldarts[oldfolder+'/kernel-input.pickle'],'inputProduct':br,'constantProduct':cr,'ownInput':input_record,'bothMeasuresOnce':True})
        index=np.arange(panels+1)[:,None]-np.arange(panels+1)[None,:]+panels
        kernel=profile[index]*mass[:,None]*mass[None,:]*b[None,:]*c
        kr=journal.write(folder+'/kernel-value.pickle',kernel);journal.json(folder+'/kernel-completed.json',{'input':arg,'value':kr})
        arg=journal.write(folder+'/inner-input.pickle',{'kernel':kr,'sourceValues':source_refs,'sourceAction':action_route})
        inner=kernel@source_values;ir=journal.write(folder+'/inner-value.pickle',inner);journal.json(folder+'/inner-completed.json',{'input':arg,'value':ir})
        arg=journal.write(folder+'/row-input.pickle',{'output':ar,'inner':ir,'operation':'ordinary transpose then matmul','rowInput':input_record})
        value=a.T@inner;vr=journal.write(folder+'/row-value.pickle',value);journal.json(folder+'/row-completed.json',{'input':arg,'value':vr})
        require(value.shape==(129,129) and np.isfinite(value).all(),'full finite paired row')
        if panels==512:
            ri,cj=(0,64,128),(0,11,128);inp=journal.write('selected-contraction/input.pickle',{'output':ar,'kernel':kr,'sourceValues':source_refs,'row':vr,'positions':ri,'columns':cj})
            direct=np.asarray([[np.sum(a[:,r,None]*kernel*source_values[None,:,s]) for s in cj] for r in ri]);dr=journal.write('selected-contraction/value.pickle',direct)
            spread=float(np.max(abs(direct-value[np.ix_(ri,cj)]))/(1+np.max(abs(direct))))
            journal.json('selected-contraction/completed.json',{'input':inp,'value':dr,'scaledDifference':spread});require(spread<2e-12,'new own source contraction check')
        cost=time.monotonic()-tick;costs.append(cost);results.append(value)
        summary={'panels':panels,'rowInput':input_record,'rowReturn':vr,'wallSeconds':cost,'rowMaxNorm':float(np.max(abs(value)))}
        journal.json(folder+'/summary.json',summary);summaries.append(summary);del kernel,inner,index;gc.collect()
    journal.json('comparison/input.json',{'first':summaries[0]['rowReturn'],'second':summaries[1]['rowReturn']})
    delta=results[1]-results[0];norm=float(np.max(abs(results[1])));absolute=float(np.max(abs(delta)))
    journal.write('comparison/difference.pickle',delta)
    comparison={'rowMaxNorm':norm,'absoluteDifference':absolute,'relativeDifference':absolute/norm if norm else None,'pilotTargetMet':absolute<=max(1e-4,.01*norm)}
    journal.json('comparison/value.json',comparison);journal.json('operation-catalogue.json',{'factorSlots':factor_routes,'fourier':fourier_routes})
    reader.postcheck()
    for rec in journal.artifacts.values():require(pilot.digest(rec['path'])==rec['sha256'],'stable new input/value/receipt bytes')
    journal.json('inputs.json',{'accepted2DPilot':reader.retain(CP,CP_SHA),'sourceAction':action_route,'consumedRoutes':reader.routes})
    checks={'status':'COMPLETED_PAIRED_FINITE_2D_ROW','case':CASE,'rowIndex':47,'sourceIndex':si,'rows':summaries,'comparison':comparison,
        'newSourceActionCalls':0,'newProfileOrRuleCalls':0,'newFactorBatches':new_factor_batches,'newFourierCalls':new_fourier,
        'acceptedFourierCalls':sum(v['disposition']=='ACCEPTED' for v in fourier_routes),'allConsumedHashesUnchanged':True,
        'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,'scope':scope['scope']}
    journal.json('checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
