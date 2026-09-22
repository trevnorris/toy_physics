#!/usr/bin/env python3
"""Two missing full 1D rows, consuming saved source actions and scalar returns."""
import argparse
import ast
import gc
import inspect
import json
from pathlib import Path
import resource
import signal
import time

import S11c_d_remaining_case_frequency_row_pilot_recover as recovery
p=recovery.prior
np,sp,io=p.np,p.sp,p.io
M,F,require,same=p.M,p.F,p.require,p.same
CP=M/'S11c_d_remaining_case_frequency_row_pilot_checkpoint.json'
CP_SHA='0c2448600fbd33b82b5e9c1a0d6b16ddbf5dbb423decff540ceec7d221098e15'
CATALOGUE_SHA='2332f982586c546d50db118825748921537eecc175f35cee5e4bd6fb15cb96f4'
PLAN=M/'S11c_d_remaining_case_frequency_rows_1d_plan.md'
ROWS=(29,30)
CASE='LAB_HELD__RHOBR_CONSTANT'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    reader,journal=io.Reader(),p.storage.Journal(base);packets={}
    def packet(record):
        route=reader.retain(record.get('path',record.get('logical')),record['sha256'])
        if route['canonical'] not in packets:packets[route['canonical']]=reader.packet(route['canonical'],record['sha256'])
        return packets[route['canonical']]
    def address(route):
        v=packet(route['packet'])
        for k in route['keys']:v=v[k]
        return v
    cp=reader.json(CP,CP_SHA);require(cp['status']=='ACCEPTED_BOUNDED_CASE_FREQUENCY_ROW_PILOT','accepted fixed-setting pilot')
    pc=reader.json(Path(cp['runDirectory'])/'checks.json',cp['checksSha256']);pilot=packet(cp['rowInput'])
    catalogue=reader.json(F/'row-pilot-review/complete/validated-point-catalogue.json',CATALOGUE_SHA)
    old_points={v['momentumHex']:v for v in catalogue};require(len(old_points)==6336,'actual reviewed completed point owners')
    reader.retain(p.READY/'complete/checks.json',p.READY_SHA)
    inventory=reader.json(p.READY/'completed-input-artifact-inventory.json',p.INVENTORY_SHA)
    def metadata(name):return reader.json(p.READY/'complete'/name,inventory[name]['sha256'])
    prior_inputs=metadata('inputs.json')
    for logical,record in prior_inputs['consumedRoutes'].items():
        if logical.endswith(('.py','.md')):reader.retain(logical,record['sha256'])
    callers=metadata('native-callers.json')
    for entry in callers.values():
        src=reader.retain(entry['file']['logical'],entry['file']['sha256']);tree=ast.parse(Path(src['canonical']).read_text())
        for name,body in entry['bodies'].items():
            node=next(n for n in tree.body if getattr(n,'name',None)==name)
            require(ast.dump(node)==ast.dump(ast.parse(body).body[0]),'whole native row/source caller source')
    source_meta=metadata('saved-prepared-basis/'+CASE+'.json');prepared=packet(source_meta['packet']);saved=prepared['original']
    selected_rows=[metadata('cases/'+CASE+'/row-'+str(ri)+'.json') for ri in ROWS]
    basis_meta=metadata('cases/'+CASE+'/source-basis-inputs.json')
    raw=packet(selected_rows[0]['sourceInputs']['packet']);scalars=packet(selected_rows[0]['scalarInputs']);system=packet(selected_rows[0]['basis'])
    physical=reader.json(selected_rows[0]['ownPhysicalRoutes']['packet']['logical'],selected_rows[0]['ownPhysicalRoutes']['packet']['sha256'])
    common=packet(physical['context']);context=common['contextPair'][0]
    require(same(*common['contextPair']) and same(*common['basisPair']),'actual complete own context and units')
    for r in physical.values():
        if isinstance(r,dict) and 'logical' in r:reader.retain(r['logical'],r['sha256'])
    settings=raw['settings'];positions=system['nodes'];nodes=saved['source_nodes'];weights=saved['source_weights']
    require(same(settings,system['settings']) and same(settings,context['settings']) and same(json.loads(json.dumps(settings)),source_meta['sourceSettings']),'full original source-rule settings join')
    require(complex(scalars['frequency'])==1-.01j and len(positions)==saved['size']==129 and nodes.shape==weights.shape==(1024,),'actual fixed frequency and full trial/source sizes')
    # Compile just the previously reviewed numerical evaluator. No old main or
    # source constructor is executed. Its principal-complex rule stays literal.
    reader.retain(Path(p.__file__),recovery.HELPER_SHA)
    reader.retain(Path(recovery.__file__),'7d9d446165f6e3301dabd5dc225f1585930f5879a1feb65fa58e96ba24f82a40')
    module,join=recovery.adapted(Path(p.__file__).read_text())
    numerical=ast.Module(body=[module.body[0]],type_ignores=[]);namespace=dict(vars(p));exec(compile(numerical,'<accepted-row-numerical-evaluator>','exec'),namespace);evaluate=namespace['evaluate']
    branch_record=pc['artifacts']['recovery-native-power-join.json'];branch=reader.json(branch_record['path'],branch_record['sha256'])
    for key in ('current','frozen'):reader.retain(branch[key]['logical'],branch[key]['sha256'])
    require(branch['literalPower']=='lambda x, p: np.asarray(x, dtype=complex) ** p','accepted actual numerical power convention')
    for path in (Path(__file__).resolve(),PLAN,Path(io.__file__),Path(p.storage.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md',Path(inspect.getsourcefile(p.integrate.quad_vec))):reader.retain(path)
    journal.json('native-and-numerical-callers.json',{'native':callers,'sourceBasisCaller':source_meta['nativeCaller'],'evaluatorJoin':join,'evaluator':ast.unparse(numerical),'principalPower':branch,'newHelper':reader.retain(Path(__file__).resolve())})
    scope={'case':CASE,'rows':list(ROWS),'frequency':{'real':1.,'imag':-.01},'chartAndSheet':'Same accepted expression branches, finite real Fourier momentum interval [-4,4], and own end pilot at1-0.01i; no root reselection or end call.',
        'sourceBound':64,'momentumBound':4,'profileBound':14,'regulator':.1,'sourceOrder':256,'profileOrder':512,
        'rules':{'29':['gk21'],'30':['gk21','gk15']},'epsabs':[1e-8,5e-9],'epsrel':[1e-6,5e-7],'norm':'max','intervalLimit':256,'maximumDistinctPointsPerRow':8192,
        'comparisonTolerance':1e-6,'measuredPriorGK21Seconds':15.934622110013152,'measuredPriorGK15Seconds':12.482678062006016,
        'batchCostDecision':'Three new full integration calls, estimated at3x the measured prior cost plus60s reserve, within900s; recheck actual remaining cost before each next call.',
        'newSourceActions':0,'scope':'Two missing full rows in the analog toy model. Fixed-setting numerical integrals only; no scattering accuracy, tiny-effect, pole/domain or numerical family reuse acceptance.'}
    journal.json('batch-scope.json',scope)
    # Preserve native pickle decoding, including Basic.xreplace. Disable actual
    # old constructors, symbolic compilation and differentiation explicitly.
    def forbidden(*args,**kwargs):raise RuntimeError('completed native science disabled in new 1D rows')
    p.source_matrix=p.independent_source=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    scalar_cache={};fourier_cache={};action_owners=[];summaries=[];costs=[]
    for selected in selected_rows:
        ri=selected['rowIndex'];root='rows/'+str(ri)
        require(selected['status']=='UNMATCHED_FULL_ROW_INPUT' and not selected['savedCompleteMatches'] and selected['layoutDimension']==1 and selected['firstOwner']==[CASE,ri],'actual previously missing row only')
        require(same(selected['sourceInputs'],selected_rows[0]['sourceInputs']) and same(selected['scalarInputs'],selected_rows[0]['scalarInputs']) and same(selected['basis'],selected_rows[0]['basis']) and same(selected['ownPhysicalRoutes'],selected_rows[0]['ownPhysicalRoutes']),'full same own caller packets')
        row=raw['bound']['rows'][ri];require(row['index']==ri and len(row['factors'])==1,'single complete factor in actual row')
        factor=row['factors'][0];si=factor['sourceIndex'];jet=address(selected['jets'][str(si)]);source=raw['bound']['sources'][0,si]
        coefficient=scalars['actual']['factor',ri,0];variable=row['limits'][0][0]
        require(same(jet['originalBoundAmplitude'],scalars['actual']['source',si]),'full actual bound source amplitude joins accepted jet return')
        require(same(jet['probe'].args,(context['zp'],)) and same(jet['amplitudeUnit'],source['amplitudeUnit']) and same(jet['integralUnit'],source['integralUnit']),'actual probe and full amplitude/measure units')
        require(not coefficient.has(sp.Integral) and coefficient.free_symbols<={variable,context['z'],context['regulator']} and source['frequency'].free_symbols<={variable},'only complete supported 1D row, no omitted profile integral')
        item=next(v for v in basis_meta if v['sourceIndex']==si);candidates=[v for v in item['savedCoefficientBasisCandidates'] if v['owner'][0]==CASE]
        require(candidates and same(item['coefficientRoute'],selected['jets'][str(si)]),'saved full source-basis candidates')
        key=(jet['probe'],tuple(jet['coefficients']),jet['amplitudeUnit'],jet['integralUnit'])
        require(not same(key,(pilot['jet']['probe'],tuple(pilot['jet']['coefficients']),pilot['jet']['amplitudeUnit'],pilot['jet']['integralUnit'])),'no completed first-pilot Fourier action is replayed')
        arrays=[]
        for candidate in candidates:
            oldjet=address(candidate['input']);require(same(key,(oldjet['probe'],tuple(oldjet['coefficients']),oldjet['amplitudeUnit'],oldjet['integralUnit'])),'full saved source coefficients and physical units')
            require(same(address(candidate['nodeRoute']),nodes) and same(address(candidate['weightRoute']),weights),'exact saved source nodes and weights')
            arrays.append(address(candidate['result']))
        matrix=arrays[0];require(matrix.shape==(1024,129) and matrix.dtype==np.dtype(complex) and np.isfinite(matrix).all() and all(same(matrix,a) for a in arrays),'all full saved basis candidate values agree')
        source_route=candidates[0]['result']
        action_matches=[i for i,a in enumerate(action_owners) if same(a['key'],key) and same(a['matrix'],matrix)]
        action_owner=action_matches[0] if action_matches else len(action_owners)
        if not action_matches:action_owners.append({'key':key,'matrix':matrix,'route':source_route})
        action_route=action_owners[action_owner]['route']
        inp={'row':row,'coefficient':coefficient,'source':source,'jet':jet,'context':context,'settings':settings,'physicalRoutes':physical,'scalarInputs':selected['scalarInputs'],'positions':positions,'fieldUnits':raw['fieldUnits'],'equationUnits':raw['equationUnits'],'sourceAction':source_route,'frequency':scalars['frequency'],'abel':raw['bound']['abel'],'pairs':raw['bound']['pairs'],'profileUnits':raw['bound']['profileUnits']}
        row_record=journal.write(root+'/input.pickle',inp)
        journal.json(root+'/source-basis-route.json',{'requestedJet':selected['jets'][str(si)],'matches':candidates,'chosen':source_route,'nodeRoute':candidates[0]['nodeRoute'],'weightRoute':candidates[0]['weightRoute'],'fullRuleAndCallerJoined':True,'newSourceActionCalls':0})
        # Reuse an old scalar return only on the same complete expression,
        # environment and physical context; computational owner stays separate.
        old_common=same(context,pilot['context']) and same(positions,pilot['positions']) and same(raw['fieldUnits'],pilot['fieldUnits']) and same(raw['equationUnits'],pilot['equationUnits']) and same(settings,pilot['settings']) and same(variable,pilot['row']['limits'][0][0])
        pilot_matches={'coefficient':old_common and same(coefficient,pilot['coefficient']) and same(factor['unit'],pilot['row']['factors'][0]['unit']),
                       'frequency':old_common and same(source['frequency'],pilot['source']['frequency'])}
        journal.write(root+'/saved-scalar-input-pairs.pickle',{'requested':inp,'acceptedRowInput':cp['rowInput'],'accepted':pilot,'matches':pilot_matches})
        units={'coefficient':factor['unit'],'frequency':('source-frequency',variable)}
        expressions={'coefficient':coefficient,'frequency':source['frequency']}
        point_cache={};counts={'newPoints':0,'savedPoints':0,'newCoefficients':0,'newFrequencies':0,'savedScalars':0,'newFourier':0,'savedFourier':0};results=[];rule_routes=[]
        def point(momentum):
            h=float(momentum).hex()
            if h in point_cache:
                counts['savedPoints']+=1;rec=point_cache[h]
                return reader.packet(rec['path'],rec['sha256'])
            require(counts['newPoints']<8192,'declared per-row point budget')
            prefix=root+'/points/'+str(counts['newPoints'])
            point_input=journal.write(prefix+'/input.pickle',{'momentum':float(momentum),'variable':variable,'positions':selected['basis'],'rowInput':row_record,'sourceAction':source_route})
            scalar_values={};scalar_routes={}
            for kind in ('coefficient','frequency'):
                expression=expressions[kind];cache_key=(kind,expression,tuple(units[kind]),h)
                if cache_key in scalar_cache:
                    rec=scalar_cache[cache_key];value=reader.packet(rec['path'],rec['sha256']);counts['savedScalars']+=1
                elif pilot_matches[kind] and h in old_points:
                    owner=old_points[h];old_input=packet(owner['input'])
                    require(old_input['momentum'].hex()==h and same(old_input['variable'],variable) and old_input['rowInput']==cp['rowInput'],'actual old full scalar caller point')
                    rec=pc['artifacts']['integrand/'+str(owner['index'])+'/'+kind+'-value.pickle'];value=reader.packet(rec['path'],rec['sha256']);counts['savedScalars']+=1
                    scalar_cache[cache_key]=rec
                else:
                    environment={variable:float(momentum)}
                    if kind=='coefficient':environment.update({context['z']:positions,context['regulator']:settings['regulator']})
                    argument=journal.write(prefix+'/'+kind+'-input.pickle',{'expression':expression,'environment':environment,'unit':units[kind],'rowInput':row_record})
                    with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
                        raw_value=evaluate(expression,environment)
                        value=np.broadcast_to(np.asarray(raw_value,complex),positions.shape) if kind=='coefficient' else complex(raw_value)
                    rec=journal.write(prefix+'/'+kind+'-value.pickle',value);journal.json(prefix+'/'+kind+'-completed.json',{'input':argument,'value':rec})
                    scalar_cache[cache_key]=rec;counts['newCoefficients' if kind=='coefficient' else 'newFrequencies']+=1
                scalar_values[kind]=value;scalar_routes[kind]=rec
            c,frequency=scalar_values['coefficient'],scalar_values['frequency']
            require(c.shape==(129,) and np.isfinite(c).all() and abs(frequency.imag)<1e-14,'finite coefficient and real source Fourier frequency')
            fourier_key=(action_owner,float(frequency.real).hex())
            if fourier_key in fourier_cache:
                fr=fourier_cache[fourier_key];source_value=reader.packet(fr['path'],fr['sha256']);counts['savedFourier']+=1
            else:
                argument=journal.write(prefix+'/fourier-input.pickle',{'frequency':frequency,'sourceAction':action_route,'sourceNodes':candidates[0]['nodeRoute'],'unit':jet['integralUnit'],'method':'literal minus-phase vector times accepted weighted source action'})
                with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):source_value=np.exp(-1j*frequency.real*nodes)@matrix
                fr=journal.write(prefix+'/fourier-value.pickle',source_value);journal.json(prefix+'/fourier-completed.json',{'input':argument,'value':fr});fourier_cache[fourier_key]=fr;counts['newFourier']+=1
            journal.json(prefix+'/operand-routes.json',{'coefficient':scalar_routes['coefficient'],'frequency':scalar_routes['frequency'],'sourceFourier':fr})
            value=c[:,None]*source_value[None,:]
            output=journal.write(prefix+'/value.pickle',value);journal.json(prefix+'/completed.json',{'input':point_input,'value':output})
            require(value.shape==(129,129) and np.isfinite(value).all(),'finite complete full row integrand')
            point_cache[h]=output;counts['newPoints']+=1
            return value
        specifications=(('gk21',1e-8,1e-6),) if ri==29 else (('gk21',1e-8,1e-6),('gk15',5e-9,5e-7))
        for rule,ea,er in specifications:
            reserve=3*max([15.934622110013152]+costs)+60;remaining=900-(time.monotonic()-started)
            prefix=root+'/quadrature/'+rule
            journal.json(prefix+'/cost-decision.json',{'priorCompletedCosts':costs,'requiredReserveSeconds':reserve,'remainingSeconds':remaining});require(reserve<remaining,'measured next small row call fits remaining guarded budget')
            argument=journal.write(prefix+'/input.pickle',{'rowInput':row_record,'sourceAction':source_route,'interval':(-4.,4.),'rule':rule,'epsabs':ea,'epsrel':er,'workers':1,'norm':'max','limit':256,'cacheSize':32*1024**2})
            tick=time.monotonic();value,error,info=p.integrate.quad_vec(point,-4.,4.,epsabs=ea,epsrel=er,norm='max',quadrature=rule,workers=1,limit=256,cache_size=32*1024**2,full_output=True);cost=time.monotonic()-tick;costs.append(cost)
            output=journal.write(prefix+'/value.pickle',{'row':value,'errorEstimate':error,'success':info.success,'status':info.status,'message':info.message,'evaluations':info.neval,'intervals':info.intervals,'intervalValues':info.integrals,'intervalErrors':info.errors,'wallSeconds':cost})
            journal.json(prefix+'/completed.json',{'input':argument,'value':output,'counts':dict(counts)})
            require(info.success and np.isfinite(value).all(),'complete converged own row')
            results.append(value);rule_routes.append({'rule':rule,'value':output,'errorEstimate':error,'evaluations':info.neval,'wallSeconds':cost})
        comparison=None
        if len(results)==2:
            journal.json(root+'/comparison-input.json',{'first':rule_routes[0]['value'],'second':rule_routes[1]['value']})
            delta=results[1]-results[0];norm=float(np.max(abs(results[1])));absolute=float(np.max(abs(delta)));scaled=absolute/(1+norm)
            journal.write(root+'/comparison-value.pickle',{'difference':delta,'rowMaxNorm':norm,'absoluteDifference':absolute,'scaledDifference':scaled})
            comparison={'absoluteDifference':absolute,'scaledDifference':scaled,'rowMaxNorm':norm};require(scaled<1e-6,'selected second-row integration rule comparison')
        summary={'rowIndex':ri,'sourceIndex':si,'rules':rule_routes,'counts':counts,'comparison':comparison,'rowMaxNorm':float(np.max(abs(results[-1]))),'rowReturn':{'packet':rule_routes[-1]['value'],'keys':['row']},'rowInput':row_record,'sourceAction':source_route}
        journal.json(root+'/summary.json',summary);summaries.append(summary)
        del point_cache,results;gc.collect()
    reader.postcheck();journal.json('inputs.json',{'acceptedPilot':reader.retain(CP,CP_SHA),'savedRowInspection':reader.retain(p.READY/'complete/checks.json',p.READY_SHA),'consumedRoutes':reader.routes})
    checks={'status':'COMPLETED_BOUNDED_MISSING_1D_FREQUENCY_ROWS','case':CASE,'rows':summaries,'newSourceActions':0,'newNativeBasisCalls':0,'newEndMapCalls':0,'allConsumedHashesUnchanged':True,'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,'scope':scope['scope']}
    journal.json('checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
