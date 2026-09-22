#!/usr/bin/env python3
"""One response sensitivity to an already saved coarser row, no quadrature."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time
from types import SimpleNamespace
import scipy.linalg as la
import S11c_d_remaining_case_frequency_lab_operator as native

saved,io,np,sp=native.saved,native.io,native.np,native.sp
M,F,require,same=native.M,native.F,native.require,native.same
NAME='S11c_d_remaining_case_frequency_lab_row_sensitivity'
CP_SHA='992b97501ed1ef96d45a66d67d6d451e8b97d338d71216872283bc176e127d30'
ROW_SHA='a4e1a9cc50f126a40730dac2d5f8886ad527f8d5b667be5ee3be83ccd8706229'
CASE='LAB_HELD__RHOBR_CONSTANT'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(io.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal=saved.Reader(),native.inputs.metadata.MetadataJournal(base);cache={}
    def packet(rec):
        r=reader.retain(rec.get('logical',rec.get('path')),rec['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(r['logical'],r['sha256'])
        return cache[r['canonical']]
    def meta(rec):return reader.json(rec.get('logical',rec.get('path')),rec['sha256'])
    cp=reader.json(M/'S11c_d_remaining_case_frequency_lab_operator_checkpoint.json',CP_SHA)
    require(cp['status']=='ACCEPTED_BOUNDED_LAB_CASE_FREQUENCY_RESPONSE','accepted complete own response')
    checks=meta(cp['artifacts']['manifest']);arts=checks['artifacts'];case=cp['cases'][CASE]
    for module in (native,native.inputs,native.inputs.metadata,saved):
        path=str(Path(module.__file__).resolve());rec=checks['consumedRoutes'][path]
        reader.retain(path,rec['sha256'])
    review=cp['independentSavedReview']['checks'];reader.retain(review['path'],review['sha256'])
    def old(name):return packet(arts[name])
    def oldmeta(name):return meta(arts[name])
    original_input=packet(case['input']);fine_interior=packet(case['interior']);fine_system=packet(case['system']);fine_response=packet(case['response'])
    original_solve=old('cases/'+CASE+'/solve/input.pickle');system=packet(original_solve['referenceSystem']);seed=packet(original_solve['fixedScales']);ends=original_solve['unitFrame']
    reader.retain(native.__file__,'06aa307c7b43c97db370cbebbb7a179042d41a374203ca16ee76bc31bf5869db')
    caller=oldmeta('native-and-numerical-callers.json')
    for entry in caller['native'].values():
        for name in ('current','frozen'):reader.retain(entry[name]['logical'],entry[name]['sha256'])
    for rec in caller['linearAlgebraSources'].values():reader.retain(rec['logical'],rec['sha256'])
    normrec=caller['nativeNorm']['file'];reader.retain(normrec['logical'],normrec['sha256'])
    norm_ns={'np':np};exec(compile(ast.parse(caller['nativeNorm']['body']),'<native norm only>','exec'),norm_ns)
    source=Path(caller['native']['S11c_d_frequency_matrix.py']['current']['canonical']).read_text();module,join=native.solve_adapter(source)
    require(same(join,caller['solveAdapter']),'whole previously reviewed native solve adapter')
    rowcp=reader.json(M/'S11c_d_remaining_case_frequency_row_2d_checkpoint.json',ROW_SHA);ra=packet(rowcp['rowInput'])
    rowchecks=reader.json(Path(rowcp['runDirectory'])/'checks.json',rowcp['checksSha256']);rowarts=rowchecks['artifacts']
    coarse_ref=rowarts['momentum-grids/512/row-value.pickle'];fine_ref=rowcp['rowReturn'];coarse=packet(coarse_ref);fine=packet(fine_ref)
    for panels in (512,1024):
        folder='momentum-grids/'+str(panels);arg=packet(rowarts[folder+'/row-input.pickle']);receipt=meta(rowarts[folder+'/row-completed.json'])
        require(receipt=={'input':rowarts[folder+'/row-input.pickle'],'value':rowarts[folder+'/row-value.pickle']} and
                arg['operation']=='outputFactors.T @ inner' and arg['transposeConjugates'] is False,'actual completed finite row return and literal contraction')
    consumed=next(v for v in original_input['rows'] if v['rowIndex']==46)
    require(consumed['input']['sha256']==rowcp['rowInput']['sha256'] and consumed['value']['packet']['sha256']==fine_ref['sha256'] and consumed['value']['keys']==[], 'exact accepted own row46 input and fine consumer return')
    require(rowcp['case']==CASE and rowcp['rowIndex']==46 and rowcp['sourceIndex']==20 and same(ra['frequency'],original_input['frequency']) and
            same(ra['settings'],system['settings']) and same(ra['positions'],system['nodes']),'own row source/frequency/settings/trial nodes')
    require(tuple(map(tuple,ra['fieldUnits']))==tuple(map(tuple,ends['fieldUnits'])) and tuple(map(tuple,ra['equationUnits']))==tuple(map(tuple,ends['rowUnits'])),'own row/equation/end physical units')
    require(coarse.shape==fine.shape==(129,129) and coarse.dtype==fine.dtype==np.dtype(complex) and np.isfinite(coarse).all() and np.isfinite(fine).all(),'complete actual coarse and fine row arrays')
    maps={};map_routes=meta(case['endMaps'])
    for side,v in map_routes.items():
        current=packet(v['actualMap']['file'])
        for key in v['actualMap']['keys']:current=current[key]
        maps[side]=dict(current,openIndices=v['openIndices'])
        for name in ('selection','physicalInput','existingWholeCallerJoin'):
            rec=v[name];reader.retain(rec['logical'],rec['sha256'])
    scope={'case':CASE,'frequency':{'real':1.,'imag':-.01},'chart':'unchanged accepted principal complex coefficient branch and continuously tracked own outgoing/incoming end maps',
        'seed':'same accepted real-seed current coordinates and fixed diagonal scales','frequencyCount':1,'newResponseCountMaximum':1,'unknowns':645,'incidents':4,
        'variation':'Only row46/source20 uses its existing512-panel finite momentum return instead of1024. All other rows, profiles, positive regulator0.1, source action, maps and physical parameters stay at their exact accepted settings.',
        'why':'Row46 has the largest observed absolute paired-grid spread among the new2D LAB rows; propagate it to actual amplitudes before deciding further refinement.',
        'targets':{'resolvedRelativeAmplitude':.01,'absoluteAmplitude':1e-4,'absoluteCurrent':1e-6},
        'priorMeasuredNativeSolveSeconds':.562564978026785,'stopping':'One saved-row substitution and at most one native solve; no rule/evaluator/map calls, no automatic extension.',
        'interpretation':'Observed sensitivity to this one row resolution only, not a rigorous error bound or complete observable/domain/regulator stability; complex amplitudes are not Hermitian flux probabilities.'}
    journal.json('sensitivity-scope.json',scope)
    inp=journal.write('sensitivity-input.pickle',{'acceptedResponseCheckpoint':reader.retain(M/'S11c_d_remaining_case_frequency_lab_operator_checkpoint.json',CP_SHA),
        'ownInput':case['input'],'fineSystem':case['system'],'fineResponse':case['response'],'fineInterior':case['interior'],'rowInput':rowcp['rowInput'],
        'coarseRow':coarse_ref,'fineRow':fine_ref,'rowCheckpoint':reader.retain(M/'S11c_d_remaining_case_frequency_row_2d_checkpoint.json',ROW_SHA),
        'maps':case['endMaps'],'nativeWholeSolve':join,'physicalContext':original_input['physicalRoutes'],'fixedScales':original_solve['fixedScales']})
    for path in (Path(__file__),M/(NAME+'_plan.md'),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    for path in (Path(__file__),M/(NAME+'_plan.md')):
        dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as stream:stream.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    journal.json('native-caller-and-reuse.json',{'acceptedCaller':arts['native-and-numerical-callers.json'],'wholeSolveReverseAST':True,
        'oldSourceEvaluatorRowProfileFourierMapCalls':0,'coefficientValues':'Exact prior full scalar returns; no new coefficient evaluator.',
        'accumulation':'Reuse the last unchanged accumulated prefix, then apply native ordered additions with saved unaffected products and genuinely new row46 products.'})
    def forbidden(*args,**kwargs):raise RuntimeError('completed source/row/end/evaluator science disabled in response sensitivity')
    native.main=native.prep.main=native.prep.evaluator=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    nonlocal_matrix=fine_interior['nonlocal'].copy();changed=[];new_products=0;reused_products=0;new_additions=0;product_cache=[]
    for i in range(5):
        for j in range(5):
            prefix='cases/'+CASE+'/blocks/nonlocal/'+str(i)+'-'+str(j);arg=old(prefix+'/input.pickle');items=arg['items']
            indexes=[n for n,item in enumerate(items) if item['term']['integralIndex']==46]
            if not indexes:continue
            require(all(same(items[n]['operand'],fine) and items[n]['operandAddress']['packet']['sha256']==fine_ref['sha256'] for n in indexes),'every actual row46 use is the accepted own fine operand')
            original_block=oldmeta(prefix+'/completed.json');require(original_block['disposition']=='NEW_ORDERED_NATIVE_BLOCK','actual own saved product/accumulation history available')
            require(same(arg['ownInput'],case['input']) and same(arg['fieldUnits'],ends['fieldUnits']) and
                    same(arg['rowUnits'],ends['rowUnits']) and same(arg['settings'],system['settings']),'whole own block caller and units')
            folder='blocks/'+str(i)+'-'+str(j);br=journal.write(folder+'/input.pickle',{'originalBlockInput':arts[prefix+'/input.pickle'],'row46ItemPositions':indexes,
                'replacementRow':coarse_ref,'wholeSensitivityInput':inp,'originalBlockReturn':original_block['value']})
            first=indexes[0]
            if first:
                acc_ref=arts[prefix+'/terms/'+str(first-1)+'/accumulated.pickle'];acc=packet(acc_ref).copy()
            else:
                zero=journal.write(folder+'/initial-input.pickle',{'shape':(129,129),'dtype':'complex128','nativeInitialAccumulator':True})
                acc=np.zeros((129,129),complex);acc_ref=journal.write(folder+'/initial-value.pickle',acc);journal.json(folder+'/initial-completed.json',{'input':zero,'value':acc_ref})
            for n in range(first,len(items)):
                item=items[n];oldtp=prefix+'/terms/'+str(n);oldarg=old(oldtp+'/input.pickle');oldreceipt=oldmeta(oldtp+'/completed.json')
                require(oldreceipt['input']==arts[oldtp+'/input.pickle'] and oldarg['item']==n and same(oldarg['operand'],item['operandAddress']),'actual ordered saved product input/receipt')
                if n in indexes:
                    coef=packet(oldarg['coefficient']);require(coef.shape==(129,) and coef.dtype==np.dtype(complex),'full accepted native coefficient return')
                    product_input={'coefficient':coef,'operand':coarse_ref,'coefficientUnit':item['unit'],
                        'fieldUnit':arg['fieldUnits'][j],'rowUnit':arg['rowUnits'][i],'physicalContext':original_input['physicalRoutes']}
                    matches=[entry for entry in product_cache if same(entry['input'],product_input)]
                    nr=journal.write(folder+'/terms/'+str(n)+'/product-input.pickle',{'coefficient':oldarg['coefficient'],'oldFullInput':arts[oldtp+'/input.pickle'],
                        'operand':coarse_ref,'originalRowInput':rowcp['rowInput'],'wholeBlockInput':br,'operation':'coefficient[:,None]*operand',
                        'fullTypedCall':product_input,'actualMatches':[entry['receipt'] for entry in matches]})
                    if matches:
                        contribution,pr=matches[0]['value'],matches[0]['reference'];reused_products+=1
                    else:
                        contribution=coef[:,None]*coarse;pr=journal.write(folder+'/terms/'+str(n)+'/product-value.pickle',contribution);new_products+=1
                    receipt_name=folder+'/terms/'+str(n)+'/product-completed.json'
                    journal.json(receipt_name,{'input':nr,'value':pr,'newCall':not bool(matches)})
                    receipt=journal.artifacts[receipt_name]
                    if not matches:product_cache.append({'input':product_input,'value':contribution,'reference':pr,'receipt':receipt})
                else:
                    pr=oldreceipt['product'];contribution=packet(pr)
                ar=journal.write(folder+'/terms/'+str(n)+'/sum-input.pickle',{'left':acc_ref,'right':pr,'nativeOrder':n,'wholeBlockInput':br})
                acc+=contribution;acc_ref=journal.write(folder+'/terms/'+str(n)+'/accumulated.pickle',acc)
                journal.json(folder+'/terms/'+str(n)+'/sum-completed.json',{'input':ar,'value':acc_ref});new_additions+=1
            require(np.isfinite(acc).all(),'actual finite new coarse-row block')
            nonlocal_matrix[i*129:(i+1)*129,j*129:(j+1)*129]=acc
            entry={'block':[i,j],'input':br,'value':acc_ref,'changedItemPositions':indexes};journal.json(folder+'/completed.json',entry);changed.append(entry)
    journal.json('changed-block-catalogue.json',{'blocks':changed,'newProducts':new_products,'reusedProducts':reused_products,'newAdditions':new_additions})
    require(changed,'one row resolution actually enters the native own operator')
    ni=journal.write('interior-input.pickle',{'fineInterior':case['interior'],'changedNonlocalBlocks':changed,'unchangedLocal':{'packet':case['interior'],'keys':['local']}})
    total=fine_interior['local']+nonlocal_matrix;ir=journal.write('interior-value.pickle',{'nonlocal':nonlocal_matrix,'total':total,'unchangedLocalReference':case['interior']})
    journal.json('interior-completed.json',{'input':ni,'value':ir});require(np.isfinite(total).all(),'full actual variant operator')
    exact_same=same(total,fine_interior['total'])
    journal.json('whole-response-call-comparison.json',{'sameFullInterior':exact_same,'sameFullMapsAndScales':True,'ownInput':inp})
    native_calls=0;solve_seconds=0.
    if exact_same:
        solved,result=fine_system,fine_response;result_ref=case['response'];system_ref=case['system']
        journal.json('complete-response-reuse.json',{'input':inp,'system':system_ref,'response':result_ref,'newSolveCalls':0})
    else:
        remain=900-(time.monotonic()-started);decision={'remainingSeconds':remain,'priorNativeSolveSeconds':.562564978026785,'reserveSeconds':3*.562564978026785+60}
        journal.json('solve-cost-decision.json',decision);require(remain>decision['reserveSeconds'],'one measured native solve fits declared budget')
        solve_base=base/'solve';solve_base.mkdir();sr=journal.write('solve/input.pickle',{'interior':ir,'wholeSensitivityInput':inp,'samePhysicalInput':case['input'],'sameEndMaps':case['endMaps'],'sameFixedScales':original_solve['fixedScales']})
        def record_call(site,name,fn,args,kwargs):
            nonlocal native_calls
            folder='solve/operations/'+str(native_calls);native_calls+=1
            ar=journal.write(folder+'/input.pickle',{'site':site,'function':name,'args':args,'kwargs':kwargs,'wholeInput':sr})
            value=fn(*args,**kwargs);vr=journal.write(folder+'/value.pickle',value);journal.json(folder+'/completed.json',{'input':ar,'value':vr});return value
        def observe(site,names,values):journal.write('solve/locals/'+site+'.pickle',{k:values[k] for k in names if k in values})
        def atomic(path,value):return journal.write(str(Path(path).relative_to(base)),value)
        ns={'np':np,'la':la,'f':SimpleNamespace(require=require,atomic_pickle=atomic),'end':SimpleNamespace(norm=norm_ns['norm']),'record_call':record_call,'observe':observe}
        exec(compile(module,'<unchanged native sensitivity solve>','exec'),ns)
        tick=time.monotonic();solved,result=ns['system_and_solve'](solve_base,{'system':system,'ends':ends},{'frequency':original_input['frequency']},{'total':total},maps,(seed['fixedRowScale'],seed['fixedColumnScale']))
        solve_seconds=time.monotonic()-tick;result_ref=journal.artifacts['solve/frequency-solution.pickle'];system_ref=journal.artifacts['solve/frequency-system.pickle']
    ci=journal.write('amplitude-comparison-input.pickle',{'fine':case['response'],'coarseVariant':result_ref,'unitFrame':original_solve['unitFrame'],
        'observable':'openOriginScattering','relativeTarget':.01,'absoluteTarget':1e-4,'variedRow':46})
    a=fine_response['openOriginScattering'];b=result['openOriginScattering'];difference=b-a;absolute=np.abs(difference);scale=np.abs(a);target=np.maximum(1e-4,.01*scale)
    value={'fine':a,'coarseVariant':b,'difference':difference,'absoluteDifference':absolute,'fineMagnitude':scale,'perEntryTarget':target,'perEntryWithinTarget':absolute<=target}
    cv=journal.write('amplitude-comparison-value.pickle',value);journal.json('amplitude-comparison-completed.json',{'input':ci,'value':cv})
    require(a.shape==b.shape==(4,4) and np.isfinite(b).all(),'complete actual four-incident open amplitudes')
    summary={'rowIndex':46,'sourceIndex':20,'newRowOrMapOrCoefficientCalls':0,'changedBlocks':len(changed),'newProducts':new_products,'reusedProducts':reused_products,'newAdditions':new_additions,
        'nativeLinearAlgebraCalls':native_calls,'nativeSolveSeconds':solve_seconds,'rank':result['rank'],'condition':result['fixedFrameCondition'],
        'maximumScaledEquationResidual':norm_ns['norm'](result['scaledEquationResidual']),'maximumBoundaryResidual':{k:norm_ns['norm'](v) for k,v in result['boundaryResiduals'].items()},
        'maximumAbsoluteAmplitudeSpread':float(np.max(absolute)),'maximumFineAmplitude':float(np.max(scale)),'allEntriesWithinDeclaredTarget':bool(np.all(absolute<=target)),
        'amplitudeEntries':[{'outgoing':i,'incident':j,'fineReal':float(a[i,j].real),'fineImag':float(a[i,j].imag),'coarseReal':float(b[i,j].real),'coarseImag':float(b[i,j].imag),
                            'absoluteSpread':float(absolute[i,j]),'relativeSpreadIfAboveAbsoluteResolution':float(absolute[i,j]/scale[i,j]) if scale[i,j]>1e-4 else None,
                            'target':float(target[i,j]),'withinTarget':bool(absolute[i,j]<=target[i,j])} for i in range(4) for j in range(4)],
        'scope':scope['interpretation']}
    journal.json('response-sensitivity-summary.json',summary)
    # Bounded saved-result readback under this guard; no arithmetic or solve repeated.
    saved_value=reader.packet(cv['path'],cv['sha256']);saved_summary=reader.json(base/'response-sensitivity-summary.json')
    require(same(saved_value,value) and saved_summary==summary,'exact saved comparison and response summary readback')
    journal.json('saved-result-readback.json',{'comparison':cv,'summary':reader.retain(base/'response-sensitivity-summary.json'),'exactSavedReadback':True,
        'nativeSourceAndSolvePreviouslyReviewed':cp['independentSavedReview']['checks'],'scienceRecomputed':0})
    journal.json('inputs.json',{'consumedRoutes':dict(reader.routes),'scope':scope,'response':result_ref,'system':system_ref})
    reader.postcheck()
    checks={'status':'COMPLETED_BOUNDED_LAB_ROW46_RESPONSE_SENSITIVITY','acceptedResponseCheckpointSha256':CP_SHA,'summary':summary,
        'response':result_ref,'system':system_ref,'consumedPaths':len(reader.routes),'consumedRoutes':dict(reader.routes),'artifacts':dict(journal.artifacts),
        'wallSeconds':time.monotonic()-started,'savedResultReadback':True,'acceptanceRequiresFinalGuard':True,'scope':scope['interpretation']}
    journal.json('checks.json',checks);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
