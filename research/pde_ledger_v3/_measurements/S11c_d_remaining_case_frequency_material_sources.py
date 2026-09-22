#!/usr/bin/env python3
"""Missing material numerical source actions with exact full-call reuse."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_remainder_prepare as prior

saved,p,io,sp,np=prior.saved,prior.p,prior.io,prior.sp,prior.np
M,F,require,same=prior.M,prior.F,prior.require,prior.same
NAME='S11c_d_remaining_case_frequency_material_sources'
READY=F/'material-inputs-recovery-01/complete'
READY_SHA='e86d67717cac30c693bf2055875283b9efbaadd24597fc74c3f7ff3a3ad34f23'


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
    ready=reader.json(READY/'checks.json',READY_SHA);require(ready['newScientificCalls']==ready['newPhysicalPickleWrites']==0,'completed saved material input handoff')
    views=[meta(v['input']) for v in ready['sources']]
    cp=reader.json(M/'S11c_d_remaining_case_frequency_remainder_preparation_checkpoint.json','3df3fb8b3bf9796eb32f22b8d86ab68f3708811bc24119e2757b93866a3a3a25')
    pc=meta(cp['manifestReferences']['artifacts']['file']);nc=meta(pc['artifacts']['native-and-new-numerical-callers.json']);native=meta(nc['savedCaller'])
    for entry in native['native'].values():
        r=reader.retain(entry['file']['logical'],entry['file']['sha256']);parsed=ast.parse(Path(r['canonical']).read_text())
        for name,body in entry['bodies'].items():require(ast.dump(next(n for n in parsed.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'whole native source/basis caller')
    for rec in (native['complexCaller'],native['nativeHeavisidePrinter']['file'],native['nativeHeavisidePrinter']['namespaceFile']):reader.retain(rec['logical'],rec['sha256'])
    reader.retain(prior.__file__,'685b35715c212ea8cc5a050db276bd849de8c992c3285647d583fac5c201ee9e')
    reader.retain(p.__file__,'c579b8881e890f0257fb7cebb5b50ca30dce356b3b2141f622f0b35deb7ad2ec')
    reader.retain(saved.recovery.__file__,'7d9d446165f6e3301dabd5dc225f1585930f5879a1feb65fa58e96ba24f82a40')
    evaluate,ej=prior.evaluator();require(same(ej,nc['evaluator']),'whole previously reviewed principal complex numerical evaluator')
    parsed=ast.parse(Path(p.__file__).read_text());functions=[next(n for n in parsed.body if getattr(n,'name',None)==name) for name in ('source_matrix','independent_source')]
    for node,key in zip(functions,('newSourceRecurrence','selectedIndependentSource')):require(ast.dump(node)==ast.dump(ast.parse(nc[key]).body[0]),'whole previously reviewed numerical source method')
    ns={'np':np,'require':require};exec(compile(ast.Module(body=functions,type_ignores=[]),'<two numerical source functions only>','exec'),ns)
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md'),Path(saved.__file__),Path(prior.inputs.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    scope={'frequency':{'real':1.,'imag':-.01},'sourceBound':64.,'sourceNodes':1024,'columns':129,'maximumDerivativeOrder':2,
        'method':'Previously reviewed new numerical expression-tree evaluation and simultaneous weighted Chebyshev recurrence; native symbolic preparation is not reconstructed.',
        'sourceRequests':len(views),'savedInputFamilies':ready['sourceInputFamilies'],'comparison':'One independent closed-trigonometric comparison for each actual encountered order0,1,2; target2e-10 scaled.',
        'stopping':'One finite source-preparation batch only; later calls need3x measured maximum action cost plus60s inside900s. No automatic expansion on a failed comparison.',
        'precisionContext':{'resolvedRelativeObservable':.01,'absoluteAmplitude':1e-4,'absoluteCurrent':1e-6},
        'scope':'Source preparation only. No row/profile/Fourier/rule/end/matrix/response calculation or accuracy/domain/pole claim. Own source, units, grades and all ordered physical limits remain referenced.'}
    journal.json('source-preparation-scope.json',scope)
    journal.json('native-and-numerical-method-joins.json',{'acceptedMethods':pc['artifacts']['native-and-new-numerical-callers.json'],'nativeCaller':nc['savedCaller'],
        'wholeEvaluatorAndRecurrenceBodiesEqual':True,'oldNativePrepareBasisRan':native['oldPrepareBasisRan'],'newNativeBasisCalls':0})
    def forbidden(*args,**kwargs):raise RuntimeError('completed native/source/row science disabled in new material preparation')
    prior.main=p.main=p.source_matrix=p.independent_source=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    # Compare genuinely saved coefficient calls, including unit and context inputs.
    old_coefficients=[]
    for entry in cp['sourceActions']:
        old=packet(entry['ownInput'])
        for order,vr in enumerate(entry['coefficientValues']):
            folder='sources/'+str(entry['sourceIndex'])+'/coefficients/'+str(order);ar=pc['artifacts'][folder+'/input.pickle'];arg=packet(ar)['requested'];receipt=meta(pc['artifacts'][folder+'/completed.json'])
            require(receipt['input']==ar and receipt['value']==vr,'actual saved coefficient input/value/receipt')
            key={'expression':arg['expression'],'variable':arg['variable'],'nodes':address(arg['nodes']),'units':arg['units'],'context':old['context']}
            old_coefficients.append({'key':key,'input':ar,'value':vr,'receipt':pc['artifacts'][folder+'/completed.json']})
    first=reader.json(M/'S11c_d_remaining_case_frequency_row_pilot_checkpoint.json','0c2448600fbd33b82b5e9c1a0d6b16ddbf5dbb423decff540ceec7d221098e15')
    fc=reader.json(Path(first['runDirectory'])/'checks.json',first['checksSha256']);fa=fc['artifacts'];old=packet(fa['source-action/input.pickle'])
    for name,ar in fa.items():
        if name.startswith('source-action/coefficients/') and name.endswith('/input.pickle'):
            arg=packet(ar);folder=name.rsplit('/',1)[0];vr=fa[folder+'/value.pickle'];receipt=meta(fa[folder+'/completed.json'])
            rr=reader.retain(receipt['value']['path'],receipt['value']['sha256']);rv=reader.retain(vr['path'],vr['sha256']);require(all(rr[k]==rv[k] for k in ('canonical','sha256','bytes')),'actual linked source22 coefficient return')
            old_coefficients.append({'key':{'expression':arg['expression'],'variable':arg['variable'],'nodes':old['nodes'],
                'units':(old['jet']['amplitudeUnit'],old['jet']['integralUnit']),'context':old['context']},'input':ar,'value':vr,'receipt':fa[folder+'/completed.json']})
    coefficients_cache=[];actions=[];results=[];checked_orders=set();costs=[];counts={'newCoefficientArrays':0,'savedCoefficientUses':0,'newActions':0,'savedActionUses':0,'independentComparisons':0}
    for view in views:
        case,si=view['case'],view['sourceIndex'];folder='cases/'+case+'/sources/'+str(si);source=address(view['nativeSource']);jet=address(view['sourceJet']);raw=packet(view['nativeSource']['packet']);physical=view['physicalRoutes'];common=packet(physical['context']);context=common['contextPair'][0]
        require(same(*common['contextPair']) and same(*common['basisPair']) and same(context['settings'],view['settings']) and same(raw['settings'],view['settings']),'full own context/unit basis/settings')
        for rec in physical.values():
            if isinstance(rec,dict) and 'logical' in rec:reader.retain(rec['logical'],rec['sha256'])
        nodes,weights=address(view['nodeRoute']),address(view['weightRoute']);units=(jet['amplitudeUnit'],jet['integralUnit']);order=len(jet['coefficients'])-1
        require(view['order']==order<=2 and view['size']==129 and nodes.shape==weights.shape==(1024,) and same(jet['amplitudeUnit'],source['amplitudeUnit']) and
                same(jet['integralUnit'],source['integralUnit']) and jet['probe'].args==(context['zp'],),'complete actual source and rule/unit inputs')
        require(not view['originalPreparedArrayCandidates']['savedCoefficientBasisCandidates'] and not any(v['wholeSourceCallMatch'] for v in view['newSavedActionComparisons']),'genuinely missing complete action in inspected saved inputs')
        key={'jet':{k:jet[k] for k in ('probe','column','coefficients','amplitudeUnit','integralUnit','originalBoundAmplitude')},'context':context,'nodes':nodes,'weights':weights,'bound':view['settings']['sourceBound'],'size':view['size']}
        matches=[v for v in actions if same(v['key'],key)]
        full=journal.write(folder+'/input.pickle',{'source':source,'jet':jet,'context':context,'fieldUnits':raw['fieldUnits'],'equationUnits':raw['equationUnits'],
            'profileUnits':raw['bound']['profileUnits'],'abel':raw['bound']['abel'],'pairs':raw['bound']['pairs'],'settings':view['settings'],'physicalRoutes':physical,
            'sourceView':next(v['input'] for v in ready['sources'] if v['case']==case and v['sourceIndex']==si),'nodeRoute':view['nodeRoute'],'weightRoute':view['weightRoute'],'bound':key['bound'],'size':key['size'],
            'actualCompletedActionMatches':[v['receipt'] for v in matches]})
        if matches:
            entry=matches[0];matrix=entry['matrix'];vr=entry['value'];croutes=entry['coefficients'];counts['savedActionUses']+=1;action_input=entry['actionInput'];comparison=None
            journal.json(folder+'/accepted-new-action-reuse.json',{'ownInput':full,'computationalInput':entry['input'],'actionInput':action_input,'value':vr,'receipt':entry['receipt'],'fullTypedInputMatch':True})
        else:
            values=[];croutes=[]
            for n,expression in enumerate(jet['coefficients']):
                require(expression.free_symbols<={context['zp']},'entire numerical source coefficient has only its own source variable')
                ckey={'expression':expression,'variable':context['zp'],'nodes':nodes,'units':units,'context':context}
                cm=[v for v in coefficients_cache if same(v['key'],ckey)];om=[v for v in old_coefficients if same(v['key'],ckey)]
                ar=journal.write(folder+'/coefficients/'+str(n)+'/input.pickle',{'expression':expression,'variable':context['zp'],'nodes':view['nodeRoute'],'units':units,'ownInput':full,'order':n,
                    'actualNewMatches':[v['receipt'] for v in cm],'actualOldMatches':[v['receipt'] for v in om]})
                if cm:value=cm[0]['array'];cr=cm[0]['value'];disposition='REUSED_COMPLETED_NEW';counts['savedCoefficientUses']+=1
                elif om:value=packet(om[0]['value']);cr=om[0]['value'];disposition='REUSED_ACCEPTED';counts['savedCoefficientUses']+=1
                else:
                    with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):value=np.broadcast_to(np.asarray(evaluate(expression,{context['zp']:nodes}),complex),nodes.shape)
                    cr=journal.write(folder+'/coefficients/'+str(n)+'/value.pickle',value);disposition='NEW';counts['newCoefficientArrays']+=1
                receipt_name=folder+'/coefficients/'+str(n)+'/completed.json';journal.json(receipt_name,{'input':ar,'value':cr,'disposition':disposition})
                require(value.shape==(1024,) and value.dtype==np.dtype(complex) and np.isfinite(value).all(),'full finite coefficient array')
                if not cm:coefficients_cache.append({'key':ckey,'array':value,'value':cr,'receipt':journal.artifacts[receipt_name]})
                values.append(value);croutes.append(cr)
            remaining=900-(time.monotonic()-started);reserve=3*max([.1]+costs)+60
            journal.json(folder+'/cost-decision.json',{'remainingSeconds':remaining,'previousActionSeconds':costs,'requiredReserveSeconds':reserve});require(remaining>reserve,'bounded next source action fits measured budget')
            action_input=journal.write(folder+'/action-input.pickle',{'coefficients':croutes,'nodes':view['nodeRoute'],'weights':view['weightRoute'],'bound':key['bound'],'size':key['size'],'units':units,'ownInput':full,'method':'unchanged reviewed simultaneous numerical Chebyshev recurrence'})
            tick=time.monotonic();matrix=ns['source_matrix'](values,nodes,weights,key['bound'],key['size']);costs.append(time.monotonic()-tick)
            vr=journal.write(folder+'/action-value.pickle',matrix);receipt_name=folder+'/action-completed.json';journal.json(receipt_name,{'input':action_input,'value':vr});counts['newActions']+=1
            require(matrix.shape==(1024,129) and matrix.dtype==np.dtype(complex) and np.isfinite(matrix).all(),'full finite material numerical source action')
            comparison=None
            if order not in checked_orders:
                ci=journal.write(folder+'/comparison-input.pickle',{'coefficients':croutes,'actionInput':action_input,'actionValue':vr,'method':'independent closed trigonometric Chebyshev derivatives','scaledTarget':2e-10})
                other=ns['independent_source'](values,nodes,weights,key['bound'],key['size']);cv=journal.write(folder+'/comparison-value.pickle',other)
                spread=float(np.max(abs(matrix-other))/(1+np.max(abs(matrix))));comparison={'input':ci,'value':cv,'scaledDifference':spread,'target':2e-10};journal.json(folder+'/comparison.json',comparison)
                require(spread<2e-10,'selected actual source recurrence/trigonometric comparison');checked_orders.add(order);counts['independentComparisons']+=1
            actions.append({'key':key,'matrix':matrix,'input':full,'actionInput':action_input,'value':vr,'coefficients':croutes,'receipt':journal.artifacts[receipt_name]})
        # Bounded exact saved readback, no numerical source operation repeated.
        restored=reader.packet(vr['path'],vr['sha256']);require(same(restored,matrix),'exact saved full source action readback')
        summary={'case':case,'sourceIndex':si,'order':order,'ownInput':full,'actionInput':action_input,'actionValue':vr,'coefficientValues':croutes,'comparison':comparison,
            'newAction':not bool(matches),'fullSavedReadback':True,'scope':scope['scope']}
        journal.json(folder+'/summary.json',summary);results.append(summary)
    journal.json('material-source-results.json',results)
    journal.json('inputs.json',{'consumedRoutes':dict(reader.routes),'savedMaterialInputChecks':reader.retain(READY/'checks.json',READY_SHA),'sourceScope':scope})
    reader.postcheck()
    checks={'status':'COMPLETED_BOUNDED_MATERIAL_NUMERICAL_SOURCES','counts':counts,'sources':results,'checkedOrders':sorted(checked_orders),'nativeBasisCompilerDerivativeCalls':0,'newRowProfileRuleEndCalls':0,
        'consumedPaths':len(reader.routes),'consumedRoutes':dict(reader.routes),'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,
        'boundedSavedReadback':True,'finalGuardRequired':True,'scope':scope['scope']}
    journal.json('checks.json',checks);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
