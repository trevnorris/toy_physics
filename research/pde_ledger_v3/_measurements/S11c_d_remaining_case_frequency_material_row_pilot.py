#!/usr/bin/env python3
"""One material row on two saved finite rules, with new batched numerical calls."""
import argparse
import ast
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_remainder_prepare as prior

saved,p,io,sp,np=prior.saved,prior.p,prior.io,prior.sp,prior.np
M,F,require,same=prior.M,prior.F,prior.require,prior.same
NAME='S11c_d_remaining_case_frequency_material_row_pilot'
CASE='MATERIAL_ADVECTED__RHO4_CONSTANT'
ROW=2
SOURCE_CP_SHA='866bf623bae5bc054bb2c603bd5de10bd7c839f8ce7e0afaacb675e3fe32a809'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal,cache=saved.Reader(),prior.inputs.MetadataJournal(base),{}
    def packet(rec):
        r=reader.retain(rec.get('logical',rec.get('path')),rec['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(r['logical'],r['sha256'])
        return cache[r['canonical']]
    def meta(rec):return reader.json(rec.get('logical',rec.get('path')),rec['sha256'])
    def address(route):
        v=packet(route['packet'])
        for key in route['keys']:v=v[tuple(key) if isinstance(key,list) else key]
        return v
    cp=reader.json(M/'S11c_d_remaining_case_frequency_material_sources_checkpoint.json',SOURCE_CP_SHA)
    require(cp['status']=='ACCEPTED_BOUNDED_MATERIAL_NUMERICAL_SOURCES','accepted own new material source actions')
    ready=reader.json(F/'material-inputs-recovery-01/complete/checks.json','e86d67717cac30c693bf2055875283b9efbaadd24597fc74c3f7ff3a3ad34f23')
    view_ref=next(v['input'] for v in ready['rows'] if v['case']==CASE and v['rowIndex']==ROW);view=meta(view_ref);selected=view['route']
    require(selected['firstOwner']==[CASE,ROW] and not selected['savedCompleteMatches'] and selected['layoutDimension']==1 and selected['status']=='UNMATCHED_FULL_ROW_INPUT','actual missing full material row')
    raw=packet(selected['sourceInputs']['packet']);scalars=packet(selected['scalarInputs']);system=packet(selected['basis']);physical=view['physicalRoutes'];common=packet(physical['context']);context=common['contextPair'][0]
    require(same(*common['contextPair']) and same(*common['basisPair']),'actual whole context and field unit basis')
    for r in physical.values():
        if isinstance(r,dict) and 'logical' in r:reader.retain(r['logical'],r['sha256'])
    row=raw['bound']['rows'][ROW];require(row['index']==ROW and len(row['factors'])==1 and len(row['limits'])==1,'entire one-factor one-dimensional row')
    factor=row['factors'][0];si=factor['sourceIndex'];source=raw['bound']['sources'][0,si];jet=address(selected['jets'][str(si)]);coefficient=scalars['actual']['factor',ROW,0]
    variable=row['limits'][0][0];positions=system['nodes'];settings=raw['settings']
    require(si==2 and same(source['frequency'],variable),'literal source frequency is the same integration variable; no compiler or new frequency evaluation')
    require(same(jet['originalBoundAmplitude'],scalars['actual']['source',si]) and same(jet['amplitudeUnit'],source['amplitudeUnit']) and same(jet['integralUnit'],source['integralUnit']),'full own source jet and amplitude/measure units')
    require(same(settings,context['settings']) and same(settings,system['settings']) and settings['momentumBound']==4 and settings['sourceBound']==64 and positions.shape==(129,) and complex(scalars['frequency'])==1-.01j,'actual fixed frequency and full native typed settings')
    require(not coefficient.has(sp.Integral) and coefficient.free_symbols<={variable,context['z'],context['regulator']},'whole supported coefficient without omitted profile or parameter')
    sa=next(v for v in cp['sourceActions'] if v['case']==CASE and v['sourceIndex']==si);own=packet(sa['ownInput']);ai=packet(sa['actionInput']);action=packet(sa['actionValue']);nodes=address(ai['nodes']);weights=address(ai['weights'])
    require(same(own['source'],source) and same(own['jet'],jet) and same(own['context'],context) and same(own['settings'],settings) and same(own['fieldUnits'],raw['fieldUnits']) and same(own['equationUnits'],raw['equationUnits']), 'full own accepted material source input')
    require(action.shape==(1024,129) and action.dtype==np.dtype(complex) and np.isfinite(action).all() and nodes.shape==weights.shape==(1024,) and ai['bound']==64 and ai['size']==129,'actual complete accepted source action and rule')
    sources_checks=meta(cp['artifacts']['manifest']);sr=sources_checks['artifacts']['cases/'+CASE+'/sources/'+str(si)+'/action-completed.json'];require(meta(sr)=={'input':sa['actionInput'],'value':sa['actionValue']},'actual saved source action receipt')
    inventory=reader.json(p.READY/'completed-input-artifact-inventory.json',p.INVENTORY_SHA);reader.retain(p.READY/'complete/checks.json',p.READY_SHA)
    native=meta(inventory['native-callers.json'])
    for entry in native.values():
        r=reader.retain(entry['file']['logical'],entry['file']['sha256']);tree=ast.parse(Path(r['canonical']).read_text())
        for name,body in entry['bodies'].items():require(ast.dump(next(n for n in tree.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'whole native source_basis/matrix_group/frequency caller')
    old_methods=reader.json(M/'S11c_d_remaining_case_frequency_remainder_preparation_checkpoint.json','3df3fb8b3bf9796eb32f22b8d86ab68f3708811bc24119e2757b93866a3a3a25');om=meta(old_methods['manifestReferences']['artifacts']['file']);nj=meta(om['artifacts']['native-and-new-numerical-callers.json'])
    reader.retain(prior.__file__,'685b35715c212ea8cc5a050db276bd849de8c992c3285647d583fac5c201ee9e');reader.retain(p.__file__,'c579b8881e890f0257fb7cebb5b50ca30dce356b3b2141f622f0b35deb7ad2ec');reader.retain(saved.recovery.__file__,'7d9d446165f6e3301dabd5dc225f1585930f5879a1feb65fa58e96ba24f82a40')
    evaluate,ej=prior.evaluator();require(same(ej,nj['evaluator']),'whole accepted principal-complex numerical evaluator')
    oldcp=reader.json(M/'S11c_d_remaining_case_frequency_row_2d_checkpoint.json','a4e1a9cc50f126a40730dac2d5f8886ad527f8d5b667be5ee3be83ccd8706229');oldchecks=reader.json(Path(oldcp['runDirectory'])/'checks.json',oldcp['checksSha256']);oa=oldchecks['artifacts']
    oldcaller=meta(oa['native-and-new-numerical-callers.json'])
    for key in ('engineCurrent','engineFrozen'):reader.retain(oldcaller[key]['logical'],oldcaller[key]['sha256'])
    rules={}
    for panels in (512,1024):
        prefix='momentum-grids/'+str(panels)+'/rule/';arg=packet(oa[prefix+'input.pickle']);value=packet(oa[prefix+'value.pickle']);receipt=meta(oa[prefix+'completed.json'])
        require(arg=={'lower':-4.,'upper':4.,'panels':panels,'method':'uniform composite trapezoid'} and receipt=={'input':oa[prefix+'input.pickle'],'value':oa[prefix+'value.pickle']},'full existing numerical rule input and actual return')
        require(value['nodes'].shape==value['weights'].shape==(panels+1,) and value['nodes'][0]==-4 and value['nodes'][-1]==4,'complete actual finite rule')
        rules[panels]=(value,oa[prefix+'value.pickle'])
    require(same(rules[512][0]['nodes'],rules[1024][0]['nodes'][::2]),'exact actual coarse nodes occur at fine even indices')
    full=journal.write('row-input.pickle',{'row':row,'coefficient':coefficient,'source':source,'jet':jet,'context':context,'settings':settings,'positions':positions,'physicalRoutes':physical,'scalarInputs':selected['scalarInputs'],
        'fieldUnits':raw['fieldUnits'],'equationUnits':raw['equationUnits'],'abel':raw['bound']['abel'],'pairs':raw['bound']['pairs'],'profileUnits':raw['bound']['profileUnits'],'frequency':scalars['frequency'],'sourceAction':sa['actionValue'],'sourceInput':sa['ownInput'],'rowView':view_ref})
    journal.json('saved-source-and-rule-routes.json',{'sourceAction':sa,'sourceCompletion':sr,'sourceNodes':ai['nodes'],'sourceWeights':ai['weights'],'rules':{str(n):r for n,(_,r) in rules.items()},'newSourceAndRuleCalls':0})
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md'),Path(saved.__file__),Path(prior.inputs.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    journal.json('native-and-new-numerical-callers.json',{'native':native,'evaluatorJoin':ej,'previousEvaluator':om['artifacts']['native-and-new-numerical-callers.json'],'newHelper':reader.retain(__file__),
        'nativeContraction':'(c * weights[:, None]).T @ s','sourceFourier':'np.exp(-1j * unique[:, None] * source_nodes[None, :]) @ amplitudes','transposeConjugates':False,'newNativeCompilerBasisRuleCalls':0})
    scope={'case':CASE,'rowIndex':ROW,'sourceIndex':si,'frequency':{'real':1.,'imag':-.01},'chartAndSheet':'Accepted principal complex powers in the full own saved expression; real finite momentum interval[-4,4], unchanged own frequency/end state.',
        'panels':[512,1024],'rule':'Exact saved uniform composite trapezoid arrays; an explicit alternative to native16/4/4 rules.','sourceNodes':1024,'sourceBound':64,'trialSize':129,'nativeSettings':settings,
        'rowTarget':{'absolute':1e-4,'relative':.01},'methodCheckTarget':2e-12,'cost':'Measure first grid; next grid only if3x maximum measured cost(initial2s planning allowance)+60s fits900s.',
        'scope':'One material1D row pilot. Fixed-row comparison only, not scattering accuracy or tiny-effect/pole/domain evidence. No automatic extension; own645x645/fourincident material response and other case work remain.'}
    journal.json('pilot-scope.json',scope)
    def forbidden(*args,**kwargs):raise RuntimeError('completed native source/profile/rule/mode/map science disabled')
    prior.main=saved.main=p.main=p.source_matrix=p.independent_source=saved.trapezoid=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    grid_results=[];costs=[];coarse=None;counts={'newCoefficientBatches':0,'newFourierBatches':0,'newCoefficientPoints':0,'newFourierPoints':0,'savedCoarsePointUses':0,'newRowContractions':0,'newSourceProfileRuleEndCalls':0}
    for panels in (512,1024):
        folder='grids/'+str(panels);rule,rr=rules[panels];mom=rule['nodes'];w=rule['weights'];remaining=900-(time.monotonic()-started);reserve=3*max([2.]+costs)+60
        journal.json(folder+'/cost-decision.json',{'remainingSeconds':remaining,'requiredReserveSeconds':reserve,'previousGridSeconds':costs});require(remaining>reserve,'next small grid fits measured budget')
        tick=time.monotonic();indices=np.arange(len(mom)) if coarse is None else np.arange(1,len(mom),2);new_mom=mom[indices]
        env={variable:new_mom[:,None],context['z']:positions[None,:],context['regulator']:settings['regulator']}
        ci=journal.write(folder+'/coefficient-input.pickle',{'expression':coefficient,'environment':env,'unit':factor['unit'],'rowInput':full,'momentumRule':rr,'indices':indices})
        with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):new_c=np.broadcast_to(np.asarray(evaluate(coefficient,env),complex),(len(indices),129))
        cv=journal.write(folder+'/coefficient-value.pickle',new_c);journal.json(folder+'/coefficient-completed.json',{'input':ci,'value':cv});counts['newCoefficientBatches']+=1;counts['newCoefficientPoints']+=len(indices)
        require(new_c.shape==(len(indices),129) and np.isfinite(new_c).all(),'complete finite new coefficient values')
        fi=journal.write(folder+'/fourier-input.pickle',{'frequencyExpression':source['frequency'],'variable':variable,'frequencies':new_mom,'sourceNodes':ai['nodes'],'sourceAction':sa['actionValue'],'sourceInput':sa['ownInput'],'unit':jet['integralUnit'],'rowInput':full,'minusPhase':True})
        with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):phase=np.exp(-1j*new_mom[:,None]*nodes[None,:])
        phase_ref=journal.write(folder+'/phase-value.pickle',phase);fourier_contraction=journal.write(folder+'/fourier-contraction-input.pickle',{'phase':phase_ref,'sourceAction':sa['actionValue'],'operation':'phase @ weightedSourceAction'})
        new_s=phase@action
        fv=journal.write(folder+'/fourier-value.pickle',new_s);journal.json(folder+'/fourier-completed.json',{'input':fi,'phase':phase_ref,'contractionInput':fourier_contraction,'value':fv});counts['newFourierBatches']+=1;counts['newFourierPoints']+=len(indices)
        require(new_s.shape==(len(indices),129) and np.isfinite(new_s).all(),'complete finite new source Fourier values')
        if coarse is None:c,s,cr,fr=new_c,new_s,cv,fv
        else:
            merge=journal.write(folder+'/exact-coarse-reuse-input.pickle',{'coarseRule':rules[512][1],'fineRule':rr,'coarseCoefficient':coarse['coefficient'],'coarseFourier':coarse['fourier'],'newCoefficient':cv,'newFourier':fv,'newIndices':indices,'savedIndices':np.arange(0,len(mom),2),'rowInput':full})
            c=np.empty((len(mom),129),complex);s=np.empty_like(c);c[::2]=coarse['c'];s[::2]=coarse['s'];c[indices]=new_c;s[indices]=new_s;counts['savedCoarsePointUses']+=len(coarse['c'])
            cr=journal.write(folder+'/full-coefficient.pickle',c);fr=journal.write(folder+'/full-fourier.pickle',s)
            journal.json(folder+'/exact-coarse-reuse-completed.json',{'input':merge,'coefficient':cr,'fourier':fr,'oldComputationalCallsRepeated':0,'savedPoints':len(coarse['c'])})
        wi=journal.write(folder+'/weighted-coefficient-input.pickle',{'coefficient':cr,'rule':rr,'operation':'c * weights[:, None]','measureIncludedOnce':True,'rowInput':full});weighted=c*w[:,None]
        wv=journal.write(folder+'/weighted-coefficient-value.pickle',weighted);journal.json(folder+'/weighted-coefficient-completed.json',{'input':wi,'value':wv})
        ri=journal.write(folder+'/row-input.pickle',{'weightedCoefficient':wv,'sourceFourier':fr,'operation':'weightedCoefficient.T @ sourceFourier','transposeConjugates':False,'rowInput':full});value=weighted.T@s
        rv=journal.write(folder+'/row-value.pickle',value);journal.json(folder+'/row-completed.json',{'input':ri,'value':rv});counts['newRowContractions']+=1
        require(value.shape==(129,129) and np.isfinite(value).all(),'full finite own129by129 material row')
        if coarse is None:
            selections=[(i,j) for i in (0,64,128) for j in (0,64,128)];fourier_selections=[(1,0),(len(indices)//2+5,64),(len(indices)-2,128)]
            check_input=journal.write('method-checks/input.pickle',{'row':rv,'weightedCoefficient':wv,'sourceFourier':fr,'phase':phase_ref,'sourceAction':sa['actionValue'],'rule':rr,'rowSelections':selections,'fourierSelections':fourier_selections,'method':'Selected literal scalar sums from the actual saved contraction operands','target':2e-12})
            direct=np.asarray([sum(weighted[t,i]*s[t,j] for t in range(len(mom))) for i,j in selections]);direct_fourier=np.asarray([sum(phase[i,t]*action[t,j] for t in range(len(nodes))) for i,j in fourier_selections])
            bad_measure=np.asarray([sum(weighted[t,i]*1.001*s[t,j] for t in range(len(mom))) for i,j in selections]);bad_index=np.asarray([sum(weighted[t,i]*s[t,(j+1)%129] for t in range(len(mom))) for i,j in selections])
            check_value=journal.write('method-checks/value.pickle',{'literalRows':direct,'literalFourier':direct_fourier,'measureMutation':bad_measure,'sourceColumnMutation':bad_index})
            selected_value=np.asarray([value[i,j] for i,j in selections]);selected_fourier=np.asarray([new_s[i,j] for i,j in fourier_selections])
            checks={'rowScaledSpread':float(np.max(abs(direct-selected_value))/(1+np.max(abs(selected_value)))),'fourierScaledSpread':float(np.max(abs(direct_fourier-selected_fourier))/(1+np.max(abs(selected_fourier)))),
                'measureResponse':float(np.max(abs(bad_measure-direct))),'sourceColumnResponse':float(np.max(abs(bad_index-direct))),'input':check_input,'value':check_value}
            journal.json('method-checks/completed.json',checks);require(checks['rowScaledSpread']<2e-12 and checks['fourierScaledSpread']<2e-12 and checks['measureResponse']>0 and checks['sourceColumnResponse']>0,'selected literal contractions and actual measure/index controls')
            coarse={'c':c,'s':s,'coefficient':cr,'fourier':fr}
        restored=reader.packet(rv['path'],rv['sha256']);require(same(restored,value),'exact saved full row readback, no recomputation')
        cost=time.monotonic()-tick;costs.append(cost);summary={'panels':panels,'rule':rr,'rowInput':ri,'rowValue':rv,'rowMaxNorm':float(np.max(abs(value))),'wallSeconds':cost,'newPoints':len(indices),'savedPoints':0 if panels==512 else len(coarse['c'])}
        journal.json(folder+'/summary.json',summary);grid_results.append((value,summary))
    comparison_input=journal.write('comparison-input.pickle',{'coarse':grid_results[0][1]['rowValue'],'fine':grid_results[1][1]['rowValue'],'absoluteTarget':1e-4,'relativeTarget':.01,'rowInput':full})
    diff=grid_results[1][0]-grid_results[0][0];norm=grid_results[1][1]['rowMaxNorm'];absolute=float(np.max(abs(diff)));target=max(1e-4,.01*norm)
    comparison_value=journal.write('comparison-value.pickle',{'difference':diff,'absoluteSpread':absolute,'fineMaxNorm':norm,'relativeSpread':absolute/norm if norm else None,'target':target,'withinPilotTarget':absolute<=target})
    comparison={'input':comparison_input,'value':comparison_value,'absoluteSpread':absolute,'fineMaxNorm':norm,'relativeSpread':absolute/norm if norm else None,'target':target,'withinPilotTarget':absolute<=target}
    journal.json('row-comparison.json',comparison)
    result={'case':CASE,'rowIndex':ROW,'sourceIndex':si,'rowInput':full,'rowReturn':grid_results[-1][1]['rowValue'],'sourceAction':sa['actionValue'],'grids':[v for _,v in grid_results],'comparison':comparison,'methodChecks':checks,'counts':counts,'exactSavedReadback':True}
    journal.json('material-row-pilot-result.json',result);journal.json('inputs.json',{'consumedRoutes':dict(reader.routes),'acceptedSources':reader.retain(M/'S11c_d_remaining_case_frequency_material_sources_checkpoint.json',SOURCE_CP_SHA),'savedRowView':view_ref})
    reader.postcheck();final={'status':'COMPLETED_BOUNDED_MATERIAL_1D_ROW_PILOT','result':result,'consumedRoutes':dict(reader.routes),'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,'scope':scope['scope']}
    journal.json('checks.json',final);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
