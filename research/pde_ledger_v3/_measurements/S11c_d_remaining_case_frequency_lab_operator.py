#!/usr/bin/env python3
"""One own LAB frequency response from accepted full rows and saved end maps."""
import argparse
import ast
import copy
import json
import inspect
from pathlib import Path
import resource
import signal
import time
from types import SimpleNamespace
import scipy.linalg as la
import S11c_d_remaining_case_frequency_lab_operator_inputs as inputs
import S11c_d_remaining_case_frequency_remainder_prepare as prep

saved, io = inputs.saved, inputs.io
np, sp, M, F, same, require = inputs.np, inputs.sp, inputs.M, inputs.F, inputs.same, inputs.require
NAME = 'S11c_d_remaining_case_frequency_lab_operator'
READY = F/'lab-operator-inputs-finish'
INVENTORY_SHA = 'c7d3e333bc5adb8afb7eb7c2997e236a8eb2b63d37240cd8c58b6a6ff623b821'
READY_SHA = '73956d81f6985fa027899949be54cd57480ff3cd2f6b4a6fc8c57a64ba708f75'
BASELINE = 'LAB_HELD__RHO4_CONSTANT'
CASE = 'LAB_HELD__RHOBR_CONSTANT'


def solve_adapter(text):
    """Keep the whole native solve; add durable operation/local observations."""
    original = next(n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == 'system_and_solve')
    changed = copy.deepcopy(original)
    calls = {}
    class Calls(ast.NodeTransformer):
        def visit_Call(self, node):
            node = self.generic_visit(node)
            name = ast.unparse(node.func)
            if name in ('np.linalg.svd', 'la.lu_factor', 'la.lu_solve'):
                site = str(len(calls)); calls[site] = name
                return ast.copy_location(ast.Call(func=ast.Name(id='record_call', ctx=ast.Load()),
                    args=[ast.Constant(site), ast.Constant(name), node.func, ast.Tuple(elts=node.args,ctx=ast.Load()),
                          ast.Dict(keys=[ast.Constant(k.arg) for k in node.keywords],values=[k.value for k in node.keywords])],keywords=[]),node)
            return node
    changed = Calls().visit(changed)
    observed = []; body = []
    for index, node in enumerate(changed.body):
        body.append(node)
        if isinstance(node, (ast.Assign, ast.AugAssign, ast.For, ast.If)):
            # Entire native loops are operations; snapshots precede the next guard.
            names = sorted({n.id for n in ast.walk(node) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)} -
                           {'i','j','d','v','block','label','index'})
            for target in ast.walk(node):
                if isinstance(target,ast.Subscript) and isinstance(target.ctx,ast.Store):
                    root=target.value
                    while isinstance(root,ast.Subscript):root=root.value
                    if isinstance(root,ast.Name):names.append(root.id)
                if isinstance(target,ast.Call) and isinstance(target.func,ast.Attribute) and target.func.attr=='append' and isinstance(target.func.value,ast.Name):
                    names.append(target.func.value.id)
            names = sorted(set(n for n in names if n not in ('current','reference','derivative')))
            if names:
                observed.append({'site':index,'names':names})
                body.extend(ast.parse('observe('+repr(str(index))+', '+repr(names)+', locals())').body)
    changed.body = body
    class Reverse(ast.NodeTransformer):
        def visit_Expr(self,node):
            if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id == 'observe': return None
            return self.generic_visit(node)
        def visit_Call(self,node):
            if isinstance(node.func,ast.Name) and node.func.id == 'record_call':
                return self.visit(ast.copy_location(ast.Call(func=node.args[2],args=node.args[3].elts,
                    keywords=[ast.keyword(arg=k.value,value=v) for k,v in zip(node.args[4].keys,node.args[4].values)]),node))
            return self.generic_visit(node)
    restored = Reverse().visit(copy.deepcopy(changed))
    require(ast.dump(restored)==ast.dump(original),'whole native system/solve reverse AST')
    return ast.fix_missing_locations(ast.Module(body=[changed],type_ignores=[])), {
        'original':ast.unparse(original),'adapted':ast.unparse(changed),'wholeReverseAST':True,
        'numericalCalls':calls,'localObservations':observed,'nativeArithmeticAndGuardsUnchanged':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(io.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    io.digest=saved.digest
    reader,journal= saved.Reader(), inputs.metadata.MetadataJournal(base)
    cache={}
    def packet(rec):
        path=rec.get('logical',rec.get('path'));r=reader.retain(path,rec['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(path,rec['sha256'])
        return cache[r['canonical']]
    def meta(rec):return reader.json(rec.get('logical',rec.get('path')),rec['sha256'])
    def address(route):
        v=packet(route['packet'])
        for k in route['keys']:v=v[tuple(k) if isinstance(k,list) else k]
        return v
    inv=reader.json(READY/'completed-input-artifact-inventory.json',INVENTORY_SHA)
    ready=reader.json(READY/'complete/checks.json',READY_SHA)
    require(ready['newScientificCalls']==0,'saved metadata handoff')
    def prior(name):return meta(inv[name])
    callers=prior('native-operator-callers.json')
    for entry in callers.values():
        for kind in ('current','frozen'):reader.retain(entry[kind]['logical'],entry[kind]['sha256'])
        parsed=ast.parse(Path(entry['current']['canonical']).read_text())
        for name,body in entry['bodies'].items():
            require(ast.dump(next(n for n in parsed.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'whole saved native operator caller')
    end_source=reader.retain(M/'S11c_d_frequency_end.py')
    norm_node=next(n for n in ast.parse(Path(end_source['canonical']).read_text()).body if getattr(n,'name',None)=='norm')
    norm_namespace={'np':np};exec(compile(ast.Module(body=[norm_node],type_ignores=[]),'<native norm only>','exec'),norm_namespace)
    native_norm=norm_namespace['norm']
    operation_sources={name:reader.retain(inspect.getsourcefile(fn)) for name,fn in [('np.linalg.svd',np.linalg.svd),('la.lu_factor',la.lu_factor),('la.lu_solve',la.lu_solve)]}
    original_text=Path(callers['S11c_d_frequency_matrix.py']['current']['canonical']).read_text()
    solve_module,solve_join=solve_adapter(original_text)
    for path,h in ((M/'S11c_d_remaining_case_frequency_row_pilot.py','c579b8881e890f0257fb7cebb5b50ca30dce356b3b2141f622f0b35deb7ad2ec'),
                   (M/'S11c_d_remaining_case_frequency_row_pilot_recover.py','7d9d446165f6e3301dabd5dc225f1585930f5879a1feb65fa58e96ba24f82a40'),
                   (Path(prep.__file__),'685b35715c212ea8cc5a050db276bd849de8c992c3285647d583fac5c201ee9e')):reader.retain(path,h)
    evaluate,evaluator_join=prep.evaluator()
    for path in (Path(__file__),M/(NAME+'_plan.md'),Path(inputs.__file__),Path(inputs.metadata.__file__),
                 Path(saved.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        reader.retain(path)
    for path in (Path(__file__),M/(NAME+'_plan.md')):
        dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    journal.json('native-and-numerical-callers.json',{'native':callers,'nativeNorm':{'file':end_source,'body':ast.unparse(norm_node)},'linearAlgebraSources':operation_sources,'solveAdapter':solve_join,'numericalEvaluator':evaluator_join,
        'assembly':'Ordered per-block native values(coefficient)[:,None]*operand accumulation; full matching blocks route saved native results.',
        'scope':'New numerical expression-tree evaluation replaces no saved symbolic compiler object. Complete own grades remain referenced.'})
    bcp=reader.json(M/'S11c_d_frequency_matrix_checkpoint.json','0c350c5f9c554c2e585920745157be826daf51775e05676a0554d8670a750663')
    def accepted(cp,name):
        rec=reader.retain(Path(cp['runDirectory'])/name,cp['artifacts'][name]['sha256']);return packet(rec),rec
    bb,bbref=accepted(bcp,'complex/frequency-binding.pickle')
    brow,browref=accepted(bcp,'complex/frequency-rows.pickle')
    bi,biref=accepted(bcp,'complex/frequency-interior.pickle')
    bs,bsref=accepted(bcp,'complex/frequency-system.pickle')
    bx,bxref=accepted(bcp,'complex/frequency-solution.pickle')
    bm,bmref=accepted(bcp,'complex/frequency-end-maps.pickle')
    seed,seedref=accepted(bcp,'seed/frequency-system.pickle')
    bgrade_path=io.REPO/'_scratch/s11c/s11c-continuum-grades-20260919/production/complete/continuum-grades.pickle'
    bgrade_ref=reader.retain(bgrade_path,'e39412750fe58734ba4db2bd84d40908d18706ae055f6ddd6e3f53f560aa9b71');bgrade=packet(bgrade_ref)
    system_path=io.REPO/'_scratch/s11c/s11c-finite-domain-20260919/production/complete/regulator/complete/finite-system.pickle'
    sysref=reader.retain(system_path,'2a043895b5e8a203cc4dd0b517ba3feec61300ab812a0404f19195cc39c0319a');system=packet(sysref)
    boundary_cp=reader.json(M/'S11c_d_remaining_case_boundary_checkpoint.json','51f7a9458895d15fa4347c124e8127598b9967c55f0ebf246096d4933c81b009')
    matrices_cp=reader.json(M/'S11c_d_remaining_case_matrices_checkpoint.json','198fd1f1c9659fe47888368e6a598c477275811a70de4eb425a493eab2e41562')
    previous=prior('saved-baseline-whole-operator-returns.json')
    require(all(same(previous[n]['packet'],reader.retain(Path(bcp['runDirectory'])/n,bcp['artifacts'][n]['sha256'])) for n in previous),'actual baseline complete returns')
    baseline_scalar={tuple(v['address']):v for v in bb['sourceJoins']}
    cases={}
    for label in (BASELINE,CASE):
        prefix='cases/'+label
        physical=prior(prefix+'/physical-operator-inputs.json');scalar=packet(physical['scalar']);routes=physical['routes']
        raw=packet(routes['originalNativeRowsSourcesProfilesSettings']);terms=packet(routes['physicalRowsTermsCharacters'])['terms']
        source=packet(routes['source']);analytic=packet(routes['analytic']);context=packet(routes['context'])
        graded,graded_ref=accepted(matrices_cp,'accepted-bindings/'+label+'/case-binding.pickle')
        ends,ends_ref=accepted(boundary_cp,'boundary-cases/'+label+'/case-boundary.pickle')
        require(same(*context['contextPair']) and same(*context['basisPair']),'own full physical context pairs')
        ctx=context['contextPair'][0]
        require(same(raw['settings'],system['settings']) and same(ctx['settings'],system['settings']) and
                same(scalar['fieldUnits'],raw['fieldUnits']) and same(scalar['equationUnits'],raw['equationUnits']) and
                tuple(map(tuple,ends['fieldUnits']))==tuple(map(tuple,scalar['fieldUnits'])) and
                tuple(map(tuple,ends['rowUnits']))==tuple(map(tuple,scalar['equationUnits'])),'full own unit and setting frame')
        require(complex(scalar['frequency'])==1-.01j and len(system['nodes'])==129,'one actual fixed frequency and full trial size')
        require(same(terms,graded['grades']['termJoins']),'whole own independent-grade cell term source')
        native_cells={(cell['row'],cell['column'],ti):(ri,coefficient) for cell in raw['bound']['cells'] if cell['test']==0 for ti,(ri,coefficient) in enumerate(cell['terms'])}
        require(set(native_cells)=={(t['row'],t['column'],t['term']) for t in terms},'all own native cell occurrences')
        for term in terms:
            ri,coef=native_cells[term['row'],term['column'],term['term']]
            require(ri==term['integralIndex'],'native cell/term integral addresses')
        rows={};row_routes=prior(prefix+'/complete-row-routes.json')
        for v in row_routes:
            rows[v['rowIndex']]=address(v['value']);packet(v['input']);meta(v['savedFullRowInputJoin'])
        require(set(rows)==set(range(len(raw['bound']['rows']))),'every own physical row retained')
        joins={tuple(v['address']):v for v in scalar['joins']}
        records={tuple(v['address']):v['record'] for v in graded['grades']['records'].values()}
        for adr,v in joins.items():
            if adr[0] in ('local','cell'):
                require(same(scalar['actual'][adr],v['bound']) and same(records[adr]['ORIGINAL'],v['original']),'actual raw coefficient -> accepted frequency binding')
        maps={};map_metadata={};end_routes=prior(prefix+'/saved-forcing-observation-end-inputs.json')
        for side,v in end_routes.items():
            rt=v['route'];current=packet(rt['map']['file'])
            for key in rt['map']['keys']:current=current[key]
            full=packet(rt['physicalSources'][0]);selection=packet(rt['selection']);meta(rt['wholeInputJoin'])
            require(full['address']==(label,side) and tuple(map(tuple,full['signature']['fieldUnits']))==tuple(map(tuple,ends['fieldUnits'])) and
                    tuple(map(tuple,full['signature']['equationUnits']))==tuple(map(tuple,ends['rowUnits'])),'actual case end source and unit frame')
            open_indices=[i for i,item in enumerate(selection['channel']['outgoing']) if item['kind']=='open']
            if 'openIndices' in current:require(current['openIndices']==open_indices,'saved complete outgoing open addresses')
            else:
                ecp=reader.json(M/'S11c_d_remaining_case_frequency_end_pilot_checkpoint.json','df9ebba15ab7c4e02bd438a70675fac1bfb13924914384dea727e560b541cf2d')
                mi,miref=accepted(ecp,'maps/fine-input.pickle')
                kinds=[c['seed']['kind'] for c in mi['clusters'] if c['seed']['direction']=='outgoing' for _ in range(c['seed']['R'].shape[1])]
                require(kinds==[x['kind'] for x in selection['channel']['outgoing']] and [i for i,k in enumerate(kinds) if k=='open']==open_indices,
                        'actual map stacking order and selected channel order')
                require(all(c['state']['frequency']==1-.01j for c in mi['clusters']),'own saved fine endpoint')
            require(current['trace'].shape==(5,5) and current['forcingAtCommonOrigin'].shape==(5,2) and len(open_indices)==2,'complete physical end dimensions')
            maps[side]=dict(current,openIndices=open_indices)
            map_metadata[side]={'actualMap':rt['map'],'openIndices':open_indices,'selection':rt['selection'],'physicalInput':rt['physicalSources'][0],
                                'existingWholeCallerJoin':rt['wholeInputJoin'],'arrayOperations':0}
        journal.json(prefix+'/end-map-consumer-routes.json',map_metadata)
        context_ref=journal.write(prefix+'/operator-input.pickle',{'case':label,'frequency':scalar['frequency'],'scalar':physical['scalar'],
            'physicalRoutes':routes,'independentGrades':graded_ref,'boundaryUnitFrame':ends_ref,'trialSystem':sysref,'rows':row_routes,'endMaps':map_metadata,
            'fixedScales':seedref,'nativeCallers':inv['native-operator-callers.json'],'settings':system['settings']})
        cases[label]={'scalar':scalar,'joins':joins,'terms':terms,'raw':raw,'records':records,'context':ctx,'rows':rows,'ends':ends,'maps':maps,
                      'input':context_ref,'rowRoutes':row_routes,'physical':physical,'mapMetadata':map_metadata}
    # All input/metadata and source joins above precede any new numerical operation.
    scope={'cases':[BASELINE,CASE],'frequency':{'real':1.,'imag':-.01},'chart':'accepted principal complex coefficient tree and saved continuously tracked own end maps',
        'sourceSeed':'actual real-seed current coordinates and saved fixed diagonal scales','frequencyCount':1,'unknowns':645,'incidentColumns':4,
        'precision':{'resolvedRelative':.01,'absoluteAmplitude':1e-4,'absoluteCurrent':1e-6},
        'stopping':'One own regular finite response; native full-rank/equation/SVD guards. No frequency sweep or automatic extension.',
        'pilotCostContext':'Baseline complete finite matrix job already completed; this one645 SVD/LU response has no quadrature or new end work. Actual cost measured below.',
        'scope':'Fixed finite-regulator approximate-boundary amplitudes; no Hermitian complex-frequency flux probabilities, continuum expansion completion or pole claim.'}
    journal.json('operator-scope.json',scope)
    def forbidden(*args,**kwargs):raise RuntimeError('completed symbolic/quadrature/end science disabled in own finite response')
    prep.main=io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    # Full consumed scalar/term/row inputs and units determine block reuse, not source labels/counts.
    def block_items(case,kind,i,j):
        scalar=case['scalar']['actual'];seq=[]
        if kind=='local':
            for adr in scalar:
                if adr[0]=='local' and adr[2:]==(i,j):
                    v=case['joins'][adr];seq.append({'address':adr,'original':v['original'],'coefficient':scalar[adr],'unit':v['unit'],
                        'limits':v['limits'],'operand':system['derivativeMatrices'][adr[1]],'operandAddress':{'packet':sysref,'keys':['derivativeMatrices',adr[1]]}})
        else:
            for term in case['terms']:
                if (term['row'],term['column'])!=(i,j):continue
                adr=('cell',i,j,term['term']);v=case['joins'][adr];ri=term['integralIndex']
                seq.append({'address':adr,'original':v['original'],'coefficient':scalar[adr],'unit':v['unit'],'limits':v['limits'],
                    'term':term,'operand':case['rows'][ri],'operandAddress':case['rowRoutes'][ri]['value']})
        return seq
    oldcase={'scalar':{'actual':{a:v['bound'] for a,v in baseline_scalar.items()}},'joins':baseline_scalar,
             'terms':bgrade['termJoins'],'rows':brow['rows'],
             'rowRoutes':[{'value':{'packet':browref,'keys':['rows',i]}} for i in range(80)]}
    def full_match(a,b):
        return len(a)==len(b) and all(same({k:v for k,v in x.items() if k!='operandAddress'},
            {k:v for k,v in y.items() if k!='operandAddress'}) for x,y in zip(a,b))
    numerical_cache=[];block_routes=[];case_blocks={};counts={'newBlocks':0,'savedBlocks':0,'newCoefficientCalls':0,'newWholeSolves':0}
    for label,case in cases.items():
        matched={};sequences={}
        for kind in ('local','nonlocal'):
            for i in range(5):
                for j in range(5):
                    key=(kind,i,j);a=block_items(case,kind,i,j);b=block_items(oldcase,kind,i,j)
                    requested=journal.write('cases/'+label+'/blocks/'+kind+f'/{i}-{j}/input.pickle',{
                        'ownInput':case['input'],'items':a,'baselineItems':b,'baselineReturn':biref,'block':key,
                        'fieldUnits':case['ends']['fieldUnits'],'rowUnits':case['ends']['rowUnits'],'settings':system['settings']})
                    unit_match=tuple(map(tuple,case['ends']['fieldUnits']))==tuple(map(tuple,bb['fieldUnits'])) and tuple(map(tuple,case['ends']['rowUnits']))==tuple(map(tuple,bb['rowUnits']))
                    matched[key]=full_match(a,b) and unit_match and same(system['settings'],bb['settings']);sequences[key]=(a,requested)
        journal.json('cases/'+label+'/whole-baseline-operator-match.json',{'all50Blocks':all(matched.values()),'blockMatches':[[*k,v] for k,v in matched.items()],
            'frequencyMatches':same(case['scalar']['frequency'],bb['frequency']),'sameTrialInput':sysref,'originalIndependentGradeInputs':case['input']})
        whole=all(matched.values()) and same(case['scalar']['frequency'],bb['frequency'])
        if label==BASELINE:
            require(whole,'complete baseline LAB operator caller before any baseline result reuse')
            require(all(all(same(case['maps'][side][key],bm[side][key]) for key in ('trace','forcingAtCommonOrigin','observationAtCommonOrigin','directIncomingSubtraction','openIndices')) for side in bm),'whole saved baseline forcing/observation/end caller')
            require(same(seed['fixedRowScale'],bs['fixedRowScale']) and same(seed['fixedColumnScale'],bs['fixedColumnScale']),'complete baseline solve fixed scales')
            journal.json('cases/'+label+'/completed-response-reuse.json',{'input':case['input'],'interior':biref,'system':bsref,'solution':bxref,'newScientificCalls':0})
            case_blocks[label]=(bi,bs,bx);continue
        require(not whole,'genuinely missing full RHOBR operator')
        matrices={k:np.empty((645,645),complex) for k in ('local','nonlocal')}
        for key,(items,arg) in sequences.items():
            kind,i,j=key;prefix='cases/'+label+'/blocks/'+kind+f'/{i}-{j}';sl=(slice(i*129,(i+1)*129),slice(j*129,(j+1)*129))
            if matched[key]:
                value=bi[kind][sl];vr={'packet':biref,'keys':[kind],'rowSlice':[i*129,(i+1)*129],'columnSlice':[j*129,(j+1)*129]}
                disposition='ACCEPTED_COMPLETE_BLOCK';counts['savedBlocks']+=1
            else:
                value=np.zeros((129,129),complex);terms_values=[]
                for ti,item in enumerate(items):
                    expr=item['coefficient'];unit=item['unit'];ctx=case['context'];env={ctx['z']:system['nodes'],ctx['regulator']:system['settings']['regulator']}
                    request={'expression':expr,'environment':env,'unit':unit,'fieldUnits':case['ends']['fieldUnits'],'rowUnits':case['ends']['rowUnits']}
                    matches=[v for v in numerical_cache if same(request,v['input'])]
                    if matches:coef,cr=matches[0]['value'],matches[0]['reference']
                    else:
                        cp='coefficient-values/'+str(counts['newCoefficientCalls']);counts['newCoefficientCalls']+=1
                        ar=journal.write(cp+'/input.pickle',{'call':request,'ownInput':case['input'],'callerBlock':arg,'address':item['address']})
                        coef=np.broadcast_to(np.asarray(evaluate(expr,env),complex),system['nodes'].shape)
                        cr=journal.write(cp+'/value.pickle',coef);journal.json(cp+'/completed.json',{'input':ar,'value':cr})
                        require(not (expr.free_symbols-set(env)) and np.isfinite(coef).all(),'full finite numerical coefficient')
                        numerical_cache.append({'input':request,'value':coef,'reference':cr})
                    ar=journal.write(prefix+'/terms/'+str(ti)+'/input.pickle',{'blockInput':arg,'item':ti,'coefficient':cr,'operand':item['operandAddress'],'nativeExpression':'coefficient[:,None]*operand'})
                    contribution=coef[:,None]*item['operand'];pr=journal.write(prefix+'/terms/'+str(ti)+'/value.pickle',contribution)
                    value+=contribution;acc=journal.write(prefix+'/terms/'+str(ti)+'/accumulated.pickle',value)
                    journal.json(prefix+'/terms/'+str(ti)+'/completed.json',{'input':ar,'product':pr,'accumulated':acc});terms_values.append((coef,item['operand']))
                vr=journal.write(prefix+'/value.pickle',value);counts['newBlocks']+=1;disposition='NEW_ORDERED_NATIVE_BLOCK'
                # Three genuinely new direct scalar contractions detect index/broadcast omissions.
                comparison=[]
                for rr,cc in ((0,0),(64,17),(128,128)):
                    expected=sum((complex(c[rr])*complex(d[rr,cc]) for c,d in terms_values),0j)
                    comparison.append({'row':rr,'column':cc,'direct':expected,'assembled':complex(value[rr,cc]),'absolute':abs(expected-value[rr,cc])})
                pr=journal.write(prefix+'/selected-direct-comparison.pickle',comparison)
                require(all(v['absolute']<1e-11*(1+abs(v['direct'])) for v in comparison),'selected literal scalar block contractions')
            matrices[kind][sl]=value
            entry={'block':list(key),'input':arg,'value':vr,'disposition':disposition};journal.json(prefix+'/completed.json',entry);block_routes.append(entry)
        ar=journal.write('cases/'+label+'/interior-input.pickle',{'ownInput':case['input'],'blocks':block_routes,'nativeOperation':'local + nonlocal'})
        matrices['total']=matrices['local']+matrices['nonlocal'];ir=journal.write('cases/'+label+'/frequency-interior.pickle',matrices)
        journal.json('cases/'+label+'/interior-completed.json',{'input':ar,'value':ir})
        require(np.isfinite(matrices['total']).all(),'full finite own interior matrix')
        cost={'elapsedInputAssemblySeconds':time.monotonic()-started,'remainingSeconds':900-(time.monotonic()-started),'solveCount':1,'minimumReserveSeconds':60}
        journal.json('solve-cost-decision.json',cost);require(cost['remainingSeconds']>60,'bounded one complete native solve budget')
        solve_base=base/'cases'/label/'solve';solve_base.mkdir()
        sr=journal.write('cases/'+label+'/solve/input.pickle',{'ownInput':case['input'],'interior':ir,'endMaps':case['mapMetadata'],'unitFrame':case['ends'],
            'referenceSystem':sysref,'fixedScales':seedref,'nativeBody':solve_join['original']})
        op_count=0
        def record_call(site,name,fn,args,kwargs):
            nonlocal op_count
            prefix='cases/'+label+'/solve/operations/'+str(op_count);op_count+=1
            ar=journal.write(prefix+'/input.pickle',{'site':site,'function':name,'args':args,'kwargs':kwargs,'wholeInput':sr})
            value=fn(*args,**kwargs);vr=journal.write(prefix+'/value.pickle',value)
            journal.json(prefix+'/completed.json',{'input':ar,'value':vr});return value
        def observe(site,names,values):
            journal.write('cases/'+label+'/solve/locals/'+site+'.pickle',{k:values[k] for k in names if k in values})
        def atomic(path,value):return journal.write(str(Path(path).relative_to(base)),value)
        namespace={'np':np,'la':la,'f':SimpleNamespace(require=require,atomic_pickle=atomic),
                   'end':SimpleNamespace(norm=native_norm),'record_call':record_call,'observe':observe}
        exec(compile(solve_module,'<whole-native-frequency-solve>','exec'),namespace)
        tick=time.monotonic();solved,result=namespace['system_and_solve'](solve_base,{'system':system,'ends':case['ends']},
                {'frequency':case['scalar']['frequency']},matrices,case['maps'],(seed['fixedRowScale'],seed['fixedColumnScale']))
        counts['newWholeSolves']+=1
        journal.json('cases/'+label+'/response-summary.json',{'nativeSolveSeconds':time.monotonic()-tick,'rank':result['rank'],
            'fixedFrameCondition':result['fixedFrameCondition'],'maximumScaledEquationResidual':float(np.max(np.abs(result['scaledEquationResidual']))),
            'maximumBoundaryResidual':{k:float(np.max(np.abs(v))) for k,v in result['boundaryResiduals'].items()},
            'openAmplitudeShape':list(result['openOriginScattering'].shape),'scope':scope['scope']})
        case_blocks[label]=(matrices,solved,result)
    delta=case_blocks[CASE][2]['openOriginScattering']-bx['openOriginScattering']
    journal.write('LAB-response-contrast.pickle',{'own':case_blocks[CASE][2]['openOriginScattering'],'baseline':bx['openOriginScattering'],'difference':delta,
        'scope':'Actual finite complex amplitudes in own real-seed current frames; no probability or error-bound interpretation.'})
    journal.json('inputs.json',{'consumedRoutes':dict(reader.routes),'scope':scope,'completedMetadataHandoff':reader.retain(READY/'complete/checks.json',READY_SHA)})
    reader.postcheck()
    checks={'status':'COMPLETED_BOUNDED_LAB_FREQUENCY_RESPONSE','counts':counts,'cases':[BASELINE,CASE],'unknowns':645,'incidentColumns':4,
        'frequency':{'real':1.,'imag':-.01},'rows':{k:len(v['rows']) for k,v in cases.items()},'consumedPaths':len(reader.routes),
        'artifacts':dict(journal.artifacts),'consumedRoutes':dict(reader.routes),'wallSeconds':time.monotonic()-started,'scope':scope['scope'],
        'acceptancePendingSavedReview':True,'completedRowOrEndCalls':0}
    journal.json('checks.json',checks);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
