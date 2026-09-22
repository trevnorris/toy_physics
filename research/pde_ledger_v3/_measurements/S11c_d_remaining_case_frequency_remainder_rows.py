#!/usr/bin/env python3
"""Four actual remainder/Abel rows by a bounded shared numerical contraction."""
import argparse
import ast
import gc
import json
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_remainder_prepare as prep

saved,p,sp,np,io=prep.saved,prep.p,prep.sp,prep.np,prep.io
M,F,require,same=prep.M,prep.F,prep.require,prep.same
NAME='S11c_d_remaining_case_frequency_remainder_rows'
CP=M/'S11c_d_remaining_case_frequency_remainder_preparation_checkpoint.json'
CP_SHA='3df3fb8b3bf9796eb32f22b8d86ab68f3708811bc24119e2757b93866a3a3a25'


def partition(expression,k,q,z,regulator):
    require(expression.func is sp.Mul,'actual literal coefficient factors')
    groups={'constant':[],'input':[],'output':[],'mixed':[]}
    for term in expression.args:
        free=term.free_symbols
        group='constant' if not free else 'output' if free<={k,z} else 'input' if free<={q} else 'mixed'
        require(group!='mixed' or free<={k,q,regulator},'covered actual mixed coefficient variables')
        groups[group].append(term)
    require(groups['mixed'] and groups['output'],'actual output phase and mixed coefficient retained')
    return groups


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal,packets=saved.Reader(),prep.inputs.MetadataJournal(base),{}
    def packet(rec):
        path=rec.get('path',rec.get('logical'))
        if path not in packets:packets[path]=reader.packet(path,rec['sha256'])
        return packets[path]
    def meta(rec):return reader.json(rec.get('path',rec.get('logical')),rec['sha256'])
    def address(route):
        value=packet(route['packet'])
        for key in route['keys']:value=value[key]
        return value
    cp=reader.json(CP,CP_SHA);require(cp['status']=='ACCEPTED_CASE_FREQUENCY_REMAINDER_NUMERICAL_PREPARATION','accepted own source/profile preparation')
    pc=reader.json(Path(cp['runDirectory'])/'checks.json',cp['checksSha256']);pa=pc['artifacts']
    sources=[packet(v['ownInput']) for v in cp['sourceActions']];owner=sources[0];context=owner['context'];settings=owner['settings'];coefficient=owner['coefficient']
    k,q=(v[0] for v in owner['row']['limits']);z,reg=context['z'],context['regulator'];unit=owner['row']['factors'][0]['unit']
    system=packet(owner['sourceInputRoute']['basis']);positions=system['nodes'];nodes=address(cp['sourceNodes']);actions=[packet(v['actionValue']) for v in cp['sourceActions']]
    require(len(positions)==129 and nodes.shape==(1024,) and same(settings,system['settings']) and settings['regulator']==.1 and settings['momentumBound']==4,'actual full finite trial/source setting')
    for own,source in zip(sources,cp['sourceActions']):
        require(same(own['coefficient'],coefficient) and same(own['row']['factors'][0]['unit'],unit) and same(own['row']['limits'],owner['row']['limits']) and
            same(own['context'],context) and same(own['settings'],settings) and same(own['fieldUnits'],owner['fieldUnits']) and same(own['equationUnits'],owner['equationUnits']) and
            same(own['physicalRoutes'],owner['physicalRoutes']) and same(own['source']['frequency'],q),'full coefficient/unit/context matches with separate actual source consumers')
        require(source['sourceIndex']==own['row']['factors'][0]['sourceIndex'] and own['size']==129,'own source index and full action')
    profile=packet(cp['profileInput']);profile_cat=meta(cp['profileFineCatalogue']);difference_grid=packet(cp['differenceGrid'])
    profile_values=np.asarray([address(v) for v in profile_cat['values']]);integral=profile['integral']
    require(same(profile['originalCoefficient'],coefficient) and same(profile['context'],context) and profile_values.shape==difference_grid.shape==(4097,), 'whole saved finite profile and coefficient join')
    groups=partition(coefficient,k,q,z,reg)
    partitions=journal.write('literal-factor-partition.pickle',{'expression':coefficient,'groups':groups,'profileInput':cp['profileInput'],'coefficientUnit':unit,'ownSourceInputs':[v['ownInput'] for v in cp['sourceActions']]})
    # Actual prior coordinate-only factor calls; no match is inferred from family counts.
    old_cps=[reader.json(M/'S11c_d_remaining_case_frequency_row_2d_checkpoint.json','a4e1a9cc50f126a40730dac2d5f8886ad527f8d5b667be5ee3be83ccd8706229'),
             reader.json(M/'S11c_d_remaining_case_frequency_row_47_checkpoint.json','2d2c98959f64f2edf48c10decd7f25051f8ac660bcd96ff4855966987eb7ad18')]
    available={}
    for oldcp in old_cps:
        oldchecks=reader.json(Path(oldcp['runDirectory'])/'checks.json',oldcp['checksSha256']);arts=oldchecks['artifacts'];oldown=packet(oldcp['rowInput'])
        require(same(oldown['context'],context) and same(oldown['positions'],positions) and same(oldown['fieldUnits'],owner['fieldUnits']) and same(oldown['equationUnits'],owner['equationUnits']),'actual saved factor physical context')
        for name,rec in arts.items():
            if name.endswith('/input.pickle') and name.startswith(('factor-calls/','coefficient-factors/')):
                arg=packet(rec);group=arg['group'];folder=name.rsplit('/',1)[0];value=packet(arts[folder+'/value.pickle']);receipt=meta(arts[folder+'/completed.json'])
                require(receipt=={'input':rec,'value':arts[folder+'/value.pickle']},'actual prior factor receipt')
                coordinates=[None] if group=='constant' else list(arg['environment'][q] if group=='input' else arg['environment'][k][:,0])
                if group=='output':require(same(arg['environment'][z],positions[None,:]),'actual full output-factor positions')
                for index,x in enumerate(coordinates):
                    available.setdefault((group,None if x is None else float(x)),[]).append({'expression':arg['expression'],'unit':arg['coefficientUnit'],
                        'value':value if group=='constant' else value[index],'input':rec,'route':{'packet':arts[folder+'/value.pickle'],'keys':[] if group=='constant' else [index]}})
    reader.retain(Path(prep.__file__),'685b35715c212ea8cc5a050db276bd849de8c992c3285647d583fac5c201ee9e')
    evaluate,evaluation_join=prep.evaluator();caller=meta(pa['native-and-new-numerical-callers.json'])
    require(evaluation_join==caller['evaluator'],'entire accepted numerical evaluator')
    for path,rec in cp['sourceFiles'].items():reader.retain(path,rec['sha256'])
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md'),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    journal.json('native-and-numerical-callers.json',{'acceptedPreparation':reader.retain(CP,CP_SHA),'nativeAndNumericalSource':pa['native-and-new-numerical-callers.json'],
        'evaluator':evaluation_join,'fullSharedCoefficientAndUnitJoin':True,'newProfileSourceBasisOrCompilerCalls':0,'actualOwnSources':[v['ownInput'] for v in cp['sourceActions']]})
    scope={'case':p.CASE,'rows':[51,52,53,54],'sources':[5,7,9,11],'frequency':{'real':1.,'imag':-.01},'panels':[1024,2048],
        'momentumBounds':[-4.,4.],'regulator':.1,'profileBounds':[-14.,14.],'profileRule':'accepted finite split768 per leg; zero new calls',
        'chart':'actual principal complex coefficient tree at1-0.01i with retained positive Abel regulator; no root reselection',
        'method':'Both-measure matrix contraction, whole shared coefficient/unit/context; own four source actions; ordinary output transpose.',
        'comparisonTarget':{'absolute':1e-4,'relative':.01},'priorMeasuredPairedGridSeconds':1.5236331930063898,
        'nextGridBudget':'3x measured coarse full-batch cost +60 inside remaining900s','mixedBatchRows':64,
        'reasonForSelectedResolution':'Positive Abel factor has momentum difference scale0.01; use one fixed1024/2048 pair, no automatic extension.',
        'scope':'Four fixed finite-frequency rows only; no scattering error bound, tiny effect, physical pole, outgoing domain or regulator limit claim.'}
    journal.json('batch-scope.json',scope)
    for source,own in zip(cp['sourceActions'],sources):journal.write('rows/'+str(source['rowIndex'])+'/input.pickle',{'sourceInput':source['ownInput'],'sourceAction':source['actionValue'],'own':own,'positions':positions,'sharedCoefficientOwner':cp['sourceActions'][0]['ownInput']})
    def forbidden(*args,**kwargs):raise RuntimeError('completed native or preparation science disabled in row batch')
    prep.main=p.source_matrix=p.independent_source=prep.special.roots_legendre=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    factor_cache={};factor_routes=[];factor_calls=0
    def products(grid,label):
        nonlocal factor_calls
        result={}
        for group in ('constant','input','output'):
            wanted=[None] if group=='constant' else list(map(float,grid));parts=[];routes=[]
            for fi,expression in enumerate(groups[group]):
                missing=[]
                for x in wanted:
                    key=(group,fi,x)
                    if key in factor_cache:continue
                    matches=[v for v in available.get((group,x),[]) if same(v['expression'],expression) and same(v['unit'],unit)]
                    if matches:
                        require(all(same(matches[0]['value'],v['value']) for v in matches),'all actual old factor returns agree')
                        v=matches[0];factor_cache[key]=(v['value'],v['route']);factor_routes.append({'group':group,'factorIndex':fi,'coordinate':x,'disposition':'ACCEPTED','input':v['input'],'value':v['route']})
                    else:missing.append(x)
                if missing:
                    env={} if group=='constant' else {q:np.asarray(missing)} if group=='input' else {k:np.asarray(missing)[:,None],z:positions[None,:]}
                    prefix='factors/'+str(factor_calls);factor_calls+=1
                    ar=journal.write(prefix+'/input.pickle',{'expression':expression,'environment':env,'coefficientUnit':unit,'partition':partitions,'group':group})
                    value=np.asarray(evaluate(expression,env),complex);value=value.reshape(()) if group=='constant' else np.broadcast_to(value,(len(missing),129) if group=='output' else (len(missing),))
                    vr=journal.write(prefix+'/value.pickle',value);journal.json(prefix+'/completed.json',{'input':ar,'value':vr});require(np.isfinite(value).all(),'finite new factor array')
                    for i,x in enumerate(missing):
                        route={'packet':vr,'keys':[] if group=='constant' else [i]};factor_cache[group,fi,x]=(value if group=='constant' else value[i],route)
                        factor_routes.append({'group':group,'factorIndex':fi,'coordinate':x,'disposition':'NEW','input':ar,'value':route})
                parts.append([factor_cache[group,fi,x][0] for x in wanted]);routes.append([factor_cache[group,fi,x][1] for x in wanted])
            if group=='constant' and group in constant_product:
                result[group]=constant_product[group];continue
            ar=journal.write('products/'+label+'/'+group+'/input.pickle',{'factors':routes,'group':group,'coefficientUnit':unit,'partition':partitions})
            value=np.ones(() if group=='constant' else (len(grid),129) if group=='output' else (len(grid),),complex)
            for item in parts:value=value*(item[0] if group=='constant' else np.asarray(item))
            vr=journal.write('products/'+label+'/'+group+'/value.pickle',value);journal.json('products/'+label+'/'+group+'/completed.json',{'input':ar,'value':vr});result[group]=(value,vr)
            if group=='constant':constant_product[group]=result[group]
        return result
    constant_product={};fourier_cache=[{} for _ in sources];fourier_routes=[];fourier_batches=0;phase_cache={}
    def fouriers(grid,label):
        nonlocal fourier_batches
        outputs=[]
        for number,(source,action) in enumerate(zip(cp['sourceActions'],actions)):
            cache=fourier_cache[number];missing=[float(qv) for qv in grid if float(qv) not in cache]
            for start in range(0,len(missing),64):
                points=np.asarray(missing[start:start+64]);prefix='fourier/'+str(fourier_batches);fourier_batches+=1
                arg=journal.write(prefix+'/input.pickle',{'frequencies':points,'sourceAction':source['actionValue'],'sourceNodes':cp['sourceNodes'],
                    'integralUnit':sources[number]['jet']['integralUnit'],'ownSourceInput':source['ownInput'],'method':'new batched minus-phase array @ new numerical recurrence source action'})
                phase_key=tuple(map(float,points))
                phase_input=journal.write(prefix+'/phase-input.pickle',{'frequencies':points,'sourceNodes':cp['sourceNodes'],'sign':-1,'method':'literal exp minus-i frequency times source position'})
                if phase_key in phase_cache:
                    phase,pr,phase_owner=phase_cache[phase_key];phase_disposition='COMPLETED_NEW_PHASE'
                else:
                    phase=np.exp(-1j*points[:,None]*nodes[None,:]);pr=journal.write(prefix+'/phase-value.pickle',phase)
                    phase_owner=phase_input;phase_cache[phase_key]=(phase,pr,phase_owner);phase_disposition='NEW_PHASE'
                journal.json(prefix+'/phase-completed.json',{'input':phase_input,'value':pr,'disposition':phase_disposition,'ownerInput':phase_owner})
                value=phase@action;vr=journal.write(prefix+'/value.pickle',value);journal.json(prefix+'/completed.json',{'input':arg,'phase':pr,'value':vr});require(np.isfinite(value).all(),'full new Fourier batch')
                for index,qv in enumerate(points):cache[float(qv)]=(value[index],{'packet':vr,'keys':[index]})
                fourier_routes.append({'sourceIndex':source['sourceIndex'],'input':arg,'value':vr,'pointCount':len(points)})
            outputs.append((np.asarray([cache[float(v)][0] for v in grid]),[cache[float(v)][1] for v in grid]))
        return outputs
    mixed_master=np.empty((2049,2049),complex);mixed_done=np.zeros((2049,2049),bool);mixed_routes=[];mixed_calls=0
    def mixed_values(grid,label,indices):
        nonlocal mixed_calls
        for start in range(0,len(indices),64):
            subset=indices[start:start+64]
            for is_even in (True,False):
                rows=subset[(subset%2==0)==is_even]
                if not len(rows):continue
                needed=indices[~mixed_done[rows[0],indices]]
                if not len(needed):continue
                require(not mixed_done[np.ix_(rows,needed)].any(),'only missing full coefficient point pairs')
                kv=-4+rows/256.;qv=-4+needed/256.;profile_indices=rows[:,None]-needed[None,:]+2048
                pv=profile_values[profile_indices];prefix='mixed/'+str(mixed_calls);mixed_calls+=1
                arg=journal.write(prefix+'/input.pickle',{'k':kv,'q':qv,'profileIndices':profile_indices,'profileCatalogue':cp['profileFineCatalogue'],
                    'regulator':settings['regulator'],'coefficientUnit':unit,'partition':partitions,'masterRows':rows,'masterColumns':needed,'literalFactors':groups['mixed']})
                env={k:kv[:,None],q:qv[None,:],reg:settings['regulator'],integral:pv};values=[];refs=[]
                for fi,expression in enumerate(groups['mixed']):
                    ar=journal.write(prefix+'/factors/'+str(fi)+'/input.pickle',{'expression':expression,'batchInput':arg,'environment':env,'coefficientUnit':unit})
                    value=np.broadcast_to(np.asarray(evaluate(expression,env),complex),(len(rows),len(needed)))
                    vr=journal.write(prefix+'/factors/'+str(fi)+'/value.pickle',value);journal.json(prefix+'/factors/'+str(fi)+'/completed.json',{'input':ar,'value':vr});values.append(value);refs.append(vr)
                ar=journal.write(prefix+'/product-input.pickle',{'factorValues':refs,'batchInput':arg})
                value=np.ones((len(rows),len(needed)),complex)
                for item in values:value=value*item
                vr=journal.write(prefix+'/product-value.pickle',value);journal.json(prefix+'/completed.json',{'input':arg,'productInput':ar,'value':vr});require(np.isfinite(value).all(),'finite actual mixed coefficient including Abel')
                mixed_master[np.ix_(rows,needed)]=value;mixed_done[np.ix_(rows,needed)]=True;mixed_routes.append({'input':arg,'value':vr})
        require(mixed_done[np.ix_(indices,indices)].all(),'all actual requested mixed pairs completed')
        return mixed_master[np.ix_(indices,indices)]
    costs=[];row_values={v['rowIndex']:[] for v in cp['sourceActions']};all_summaries=[]
    old46=old_cps[0];oc=reader.json(Path(old46['runDirectory'])/'checks.json',old46['checksSha256']);old_rule=oc['artifacts']['momentum-grids/1024/rule/value.pickle']
    for panels in (1024,2048):
        label=str(panels);prefix='grids/'+label;remaining=900-(time.monotonic()-started);reserve=3*max(costs+[1.5236331930063898])+60
        journal.json(prefix+'/cost-decision.json',{'remainingSeconds':remaining,'requiredSeconds':reserve,'completedGridSeconds':costs});require(reserve<remaining,'bounded measured next batch fits')
        tick=time.monotonic()
        if panels==1024:
            rule=packet(old_rule);grid,mass=rule['nodes'],rule['weights'];rule_ref=old_rule
        else:
            arg=journal.write(prefix+'/rule-input.pickle',{'method':'new explicit uniform composite trapezoid','bounds':(-4.,4.),'panels':panels})
            grid=np.linspace(-4,4,panels+1);mass=np.full(panels+1,8/panels);mass[[0,-1]]*=.5
            rule_ref=journal.write(prefix+'/rule-value.pickle',{'nodes':grid,'weights':mass});journal.json(prefix+'/rule-completed.json',{'input':arg,'value':rule_ref})
        indices=np.arange(0,2049,2 if panels==1024 else 1)
        require(np.array_equal(grid,-4+indices/256.),'literal saved/fine coordinate addresses')
        parts=products(grid,label);ft=fouriers(grid,label);mixed=mixed_values(grid,label,indices)
        a,ar=parts['output'];b,br=parts['input'];c,cr=parts['constant']
        mi=journal.write(prefix+'/mixed-input.pickle',{'partition':partitions,'masterIndices':indices,'completedBatches':mixed_routes})
        mr=journal.write(prefix+'/mixed-value.pickle',mixed);journal.json(prefix+'/mixed-completed.json',{'input':mi,'value':mr})
        ki=journal.write(prefix+'/kernel-input.pickle',{'mixed':mr,'inputProduct':br,'constantProduct':cr,'rule':rule_ref,'bothMeasuresOnce':True})
        kernel=mixed*b[None,:]*c*mass[:,None]*mass[None,:];kr=journal.write(prefix+'/kernel-value.pickle',kernel);journal.json(prefix+'/kernel-completed.json',{'input':ki,'value':kr})
        for number,source in enumerate(cp['sourceActions']):
            ri=source['rowIndex'];fourier,fourier_refs=ft[number];folder='rows/'+str(ri)+'/'+label
            ii=journal.write(folder+'/inner-input.pickle',{'kernel':kr,'sourceFourier':fourier_refs,'ownSourceInput':source['ownInput']})
            inner=kernel@fourier;ir=journal.write(folder+'/inner-value.pickle',inner);journal.json(folder+'/inner-completed.json',{'input':ii,'value':ir})
            vi=journal.write(folder+'/row-input.pickle',{'outputProduct':ar,'inner':ir,'operation':'ordinary transpose @','ownRowInput':journal.artifacts['rows/'+str(ri)+'/input.pickle']})
            value=a.T@inner;vr=journal.write(folder+'/row-value.pickle',value);journal.json(folder+'/row-completed.json',{'input':vi,'value':vr});require(value.shape==(129,129) and np.isfinite(value).all(),'full own finite129x129 row')
            if panels==1024 and number==0:
                rows,columns=[0,64,128],[0,11,128]
                ci=journal.write('selected-contraction/input.pickle',{'output':ar,'kernel':kr,'sourceFourier':fourier_refs,'row':vr,'positions':rows,'columns':columns})
                direct=np.asarray([[np.sum(a[:,r,None]*kernel*fourier[None,:,s]) for s in columns] for r in rows]);dv=journal.write('selected-contraction/value.pickle',direct)
                spread=float(np.max(abs(direct-value[np.ix_(rows,columns)]))/(1+np.max(abs(direct))))
                journal.json('selected-contraction/completed.json',{'input':ci,'value':dv,'scaledDifference':spread});require(spread<2e-12,'selected literal double-sum contraction')
            summary={'rowIndex':ri,'sourceIndex':source['sourceIndex'],'panels':panels,'rowInput':journal.artifacts['rows/'+str(ri)+'/input.pickle'],'value':vr,'norm':float(np.max(abs(value)))}
            journal.json(folder+'/summary.json',summary);all_summaries.append(summary);row_values[ri].append(value)
        if panels==1024:
            controls=[]
            for gi,gj,zi in ((512,560,32),(600,480,64),(440,500,96)):
                kv,qv=grid[gi],grid[gj];pv=profile_values[indices[gi]-indices[gj]+2048];env={k:kv,q:qv,z:positions[zi],reg:settings['regulator'],integral:pv}
                index=len(controls);ci=journal.write('coefficient-checks/'+str(index)+'/input.pickle',{'expression':coefficient,'environment':env,'partition':partitions,'factorProducts':[ar,br,cr,mr],'indices':(gi,gj,zi)})
                direct=complex(evaluate(coefficient,env));dv=journal.write('coefficient-checks/'+str(index)+'/direct-value.pickle',direct)
                composed=c*a[gi,zi]*b[gj]*mixed[gi,gj];cv=journal.write('coefficient-checks/'+str(index)+'/composed-value.pickle',composed)
                error=float(abs(direct-composed)/(1+abs(direct)));journal.json('coefficient-checks/'+str(index)+'/completed.json',{'input':ci,'direct':dv,'composed':cv,'scaledDifference':error});require(error<2e-12,'whole coefficient agrees with literal factor contraction')
                controls.append({'input':ci,'value':dv,'scaledDifference':error})
            # A new actual omitted-Abel control uses the same saved finite profile.
            carrier=next(v for v in groups['mixed'] if v.has(sp.Integral))
            ci=journal.write('abel-control/input.pickle',{'wholeInput':controls[0]['input'],'changedCarrier':carrier,'changedValue':profile_values[indices[512]-indices[560]+2048]})
            inp=packet(controls[0]['input']);env=dict(inp['environment']);env[carrier]=profile_values[indices[512]-indices[560]+2048]
            wrong=complex(evaluate(coefficient,env));wv=journal.write('abel-control/value.pickle',wrong);response=float(abs(wrong-packet(controls[0]['value'])))
            journal.json('abel-control/completed.json',{'input':ci,'value':wv,'response':response});require(response>1e-14,'retained Abel term responds to omission')
        cost=time.monotonic()-tick;costs.append(cost);journal.json(prefix+'/summary.json',{'panels':panels,'wallSeconds':cost,'rowCount':4,'rule':rule_ref})
        del kernel,mixed;gc.collect()
    comparisons=[]
    for ri,values in row_values.items():
        ci=journal.write('rows/'+str(ri)+'/comparison-input.pickle',{'coarse':journal.artifacts['rows/'+str(ri)+'/1024/row-value.pickle'],'fine':journal.artifacts['rows/'+str(ri)+'/2048/row-value.pickle']})
        delta=values[1]-values[0];dv=journal.write('rows/'+str(ri)+'/difference.pickle',delta);norm=float(np.max(abs(values[1])));absolute=float(np.max(abs(delta)))
        comparison={'rowIndex':ri,'input':ci,'difference':dv,'fineNorm':norm,'absoluteDifference':absolute,'relativeDifference':absolute/norm if norm else None,'targetMet':absolute<=max(1e-4,.01*norm)}
        journal.json('rows/'+str(ri)+'/comparison.json',comparison);comparisons.append(comparison)
    journal.json('operation-catalogue.json',{'factors':factor_routes,'fourierBatches':fourier_routes,'mixedBatches':mixed_routes})
    reader.postcheck()
    for rec in journal.artifacts.values():require(saved.digest(rec['path'])==rec['sha256'],'new inputs/results/receipts unchanged')
    journal.json('inputs.json',{'acceptedPreparation':reader.retain(CP,CP_SHA),'consumedRoutes':reader.routes})
    checks={'status':'COMPLETED_BOUNDED_FOUR_REMAINDER_ROWS','rows':all_summaries,'comparisons':comparisons,'newFactorBatches':factor_calls,
        'acceptedFactorSlots':sum(v['disposition']=='ACCEPTED' for v in factor_routes),'newFourierBatches':fourier_batches,'newFourierValues':sum(v['pointCount'] for v in fourier_routes),
        'newMixedBatches':mixed_calls,'mixedPointPairs':int(mixed_done.sum()),'newFourierPhaseBatches':len(phase_cache),'newSourceProfileOrNativeCalls':0,'gridSeconds':costs,'allConsumedHashesUnchanged':True,
        'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,'scope':scope['scope']}
    journal.json('checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
