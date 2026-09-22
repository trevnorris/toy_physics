#!/usr/bin/env python3
"""Missing native rational determinant with exact saved row and return routing."""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import resource
import signal
import textwrap
import time
from types import SimpleNamespace
import S11c_d_remaining_case_frequency_end_rational as previous

h,f,sp,engine=previous.h,previous.f,previous.sp,previous.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_rational_checkpoint.json'
CP_SHA='d9f275cc6c266698d421ca91d0b2e895ff654d506b5243eb48d86bc1b2f4f2d7'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_determinant_plan.md'
OWNER=previous.OWNER
same=previous.same
equal=previous.equal
source_body=previous.source_body


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(f.digest(CP)==CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_RATIONAL_ENTRIES','accepted rational entry checkpoint')
    for root,stage,sha in ((origin.parent,'frequency_end_rational',cp['checksSha256']),(vr,'validate',cp['validation']['checksSha256'])):
        h.source.receipts.inspect_guard(root,stage)
        checks=origin/'checks.json' if stage!='validate' else vr/'checks.json'
        f.require(f.digest(checks)==sha and checks.read_bytes()==(root/(stage+'.stdout')).read_bytes(),'clean accepted checks/stdout')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],
              'acceptedRationalEntries':{'checkpoint':str(CP),'checkpointSha256':CP_SHA,'runDirectory':str(origin),
                  'checksSha256':cp['checksSha256'],'validationDirectory':str(vr),'validatorChecksSha256':cp['validation']['checksSha256']}}
    ref=lambda path,name,sha:h.source.reference(base,manifest,path,name,sha)
    for name,item in cp['artifacts'].items():ref(origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        p=origin/name;f.require(p.is_symlink() and str(p.readlink())==item['original'] and str(p.resolve())==item['resolvedOriginal'] and p.stat().st_size==item['bytes'] and f.digest(p)==item['sha256'],'original rational reference identity')
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'accepted current/frozen source')
        ref(origin/'source'/name,'source/'+name,value)
    for path,name in ((CP,'accepted-rational-checkpoint.json'),(origin/'checks.json','accepted-rational-checks.json'),
                      (origin/'inputs.json','accepted-rational-manifest.json'),(vr/'checks.json','accepted-rational-validation.json')):ref(path,name,f.digest(path))
    for name,item in cp['validation']['artifacts'].items():ref(vr/name,'accepted-rational-validation/'+name,item['sha256'])
    diagnostic=f.STORE/'s11c-frequency-source-20260919/threshold-repair/elimination.pickle'
    f.require(manifest['inputPackets'][str(diagnostic)]=='e27be9461983fcb9a3bfc9e531d9479c9df85ec6f7f68a1fe2502ea39155ff27','accepted original determinant diagnostic hash')
    ref(diagnostic,'original-determinant-diagnostic-elimination.pickle',manifest['inputPackets'][str(diagnostic)])
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));value=f.digest(path);f.require(name not in manifest['sourceFiles'],'fresh determinant constructor source')
        manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return cp,manifest


def native_adapter():
    """Keep the whole determinant body; route completed rows/tails and calls."""
    original=ast.parse(textwrap.dedent(inspect.getsource(engine.FullPencilModes.rational_determinant))).body[0]
    sites=[]
    class Calls(ast.NodeTransformer):
        def visit_Call(self,node):
            old=copy.deepcopy(node);name=ast.unparse(node.func);node=self.generic_visit(node)
            if name not in ('sp.cancel','sp.denom','sp.lcm','sp.ImmutableMatrix','sp.prod','sp.fraction','polynomial_matrix.det'):return node
            site=len(sites);sites.append({'site':site,'call':ast.unparse(old),'nativeAST':ast.dump(old)})
            operation=name.split('.')[-1];args=[ast.Constant(site),ast.Constant(operation)]
            if name=='polynomial_matrix.det':args.append(node.func.value)
            args.extend(node.args)
            return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='router',ctx=ast.Load()),attr='call',ctx=ast.Load()),args=args,keywords=node.keywords),node)
    routed=Calls().visit(copy.deepcopy(original))
    loop=next(n for n in routed.body if isinstance(n,ast.For));observed=[]
    for stmt in loop.body:
        observed.append(stmt)
        if isinstance(stmt,ast.Assign):
            for target in stmt.targets:
                if isinstance(target,ast.Name):observed.extend(ast.parse(f"router.observe('{target.id}', {target.id})").body)
        elif isinstance(stmt,ast.Expr) and isinstance(stmt.value,ast.Call) and ast.unparse(stmt.value.func)=='rows.append':
            observed.extend(ast.parse("router.observe('clearedRow', rows[-1])").body)
    observed.extend(ast.parse('router.row_complete(i,rows[-1],denominators[-1])').body)
    branch=ast.parse('if router.restore_row(matrix,i,rows,denominators):\n pass\nelse:\n pass').body[0];branch.orelse=observed;loop.body=[branch]
    for idx,stmt in enumerate(routed.body):
        if isinstance(stmt,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='polynomial_matrix' for t in stmt.targets):matrix_index=idx;break
    tail=[]
    for stmt in routed.body[matrix_index+1:]:
        tail.append(stmt)
        if isinstance(stmt,ast.Assign):
            for target in stmt.targets:
                if isinstance(target,ast.Name):tail.extend(ast.parse(f"router.observe('{target.id}', {target.id})").body)
        if isinstance(stmt,ast.Return):stmt.value=ast.Call(func=ast.Attribute(value=ast.Name(id='router',ctx=ast.Load()),attr='finish',ctx=ast.Load()),args=[stmt.value],keywords=[])
    branch_tail=ast.parse('if router.restore_tail(polynomial_matrix,denominators):\n return router.finish(router.saved_tail)\nelse:\n pass').body[0];branch_tail.orelse=tail
    routed.body[matrix_index+1:]=ast.parse("router.observe('polynomial_matrix', polynomial_matrix)").body+[branch_tail]
    class Reverse(ast.NodeTransformer):
        def visit_If(self,node):
            if isinstance(node.test,ast.Call) and ast.unparse(node.test.func) in ('router.restore_row','router.restore_tail'):
                return [self.visit(n) for n in node.orelse if not (isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func) in ('router.observe','router.row_complete'))]
            return self.generic_visit(node)
        def visit_Expr(self,node):
            if isinstance(node.value,ast.Call) and ast.unparse(node.value.func)=='router.observe':return None
            return self.generic_visit(node)
        def visit_Call(self,node):
            if ast.unparse(node.func)=='router.call':return ast.parse(sites[node.args[0].value]['call'],mode='eval').body
            if ast.unparse(node.func)=='router.finish':return self.visit(node.args[0])
            return self.generic_visit(node)
    restored=Reverse().visit(copy.deepcopy(routed));f.require(ast.dump(restored)==ast.dump(original),'whole native rational_determinant reverse AST')
    routed.name='native_determinant';routed.decorator_list=[];routed.args.args.append(ast.arg(arg='router'))
    module=ast.fix_missing_locations(ast.Module(body=[routed],type_ignores=[]));env={}
    exec(compile(module,'<native rational determinant with exact saved operand routing>','exec'),env)
    return env['native_determinant'],{'wholeNativeAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'entireOriginalBodyReverseAST':True,'nativeOperationSites':sites,'rowReuseGuard':'full literal row/unit/coordinate/caller input',
        'tailReuseGuard':'full literal cleared matrix/denominators/unit/coordinate/caller input','adapterSource':ast.unparse(module)}


def joined_sources(base,manifest):
    sources={'rationalDeterminant':engine.FullPencilModes.rational_determinant,'endSources':h.source.q.end_sources,
        'endTables':h.inputs.chart.end_tables,'acceptedThresholdAdapter':previous.previous.accepted_threshold.end_adapter,
        'acceptedRealAxisAnalysis':previous.previous.accepted_threshold.real_axis_analysis,
        'rationalEntryAdapter':previous.native_adapter,'rationalEntryConstructor':previous.construct}
    result={name:source_body(value) for name,value in sources.items()}
    old=json.loads((base/'native-end-rational-joins.json').read_text())
    for name,item in result.items():
        relative=str(Path(item['path']).relative_to(f.ROOT));f.require(manifest['sourceFiles'][relative]==item['sha256'],'accepted whole native source')
    for name in ('rationalDeterminant','endSources','endTables','acceptedThresholdAdapter','acceptedRealAxisAnalysis'):same(result[name],old[name])
    fn,proof=native_adapter();result['determinantAdapter']=proof
    for name in ('cancel','denom','lcm','ImmutableMatrix','prod','fraction'):result['sympy_'+name]=source_body(getattr(sp,name))
    result['matrixDet']=source_body(sp.MatrixBase.det);result['sympyVersion']=sp.__version__
    f.save(base/'native-end-determinant-joins.json',result)
    return fn,result


def forbidden(*args,**kwargs):raise RuntimeError('only genuinely missing native determinant operations are enabled')


def prohibit():
    functions={name:getattr(sp,name) for name in ('cancel','denom','lcm','ImmutableMatrix','prod','fraction')}
    functions['det']=sp.MatrixBase.det
    previous.previous.prohibit()
    for name,value in functions.items():
        if name!='det':setattr(sp,name,value)
    for name in ('load','construct','native_adapter','main'):setattr(previous,name,forbidden)
    for name in ('__init__','operation','start','observe'):setattr(previous.OperationRouter,name,forbidden)
    return functions


class Router:
    def __init__(self,base,own,baselines,functions):
        self.base,self.own,self.baselines,self.functions=base,own,baselines,functions
        self.units=own['fullBindingInput']['units']['PENCIL_PLUS'];self.coordinates=(own['frequencyCoordinate'],own['radicalCoordinate'])
        self.rows=[];self.calls=[];self.atlas=[];self.row=None;self.current_folder=None;self.saved_tail=None
        self.unit_context={'physicalInputUnits':self.units,'coefficientCoordinateUnits':((0,0,0),(0,0,0)),
            'coordinateMap':own['coordinateMap'],'coordinates':self.coordinates,
            'nativeMatrixConvention':'dimensionless coefficients between explicit physical equation/field unit reference bases'}

    def row_key(self,matrix,i,units,coordinates):
        return {'entries':tuple(matrix[i,j] for j in range(matrix.cols)),'columns':matrix.cols,
            'physicalEntryUnits':units[i],'coefficientCoordinates':coordinates,'coefficientCoordinateUnits':((0,0,0),(0,0,0))}

    def restore_row(self,matrix,i,rows,denominators):
        self.row=i;folder=self.base/'determinant-rows'/str(i);folder.mkdir(parents=True);self.current_folder=folder
        requested=self.row_key(matrix,i,self.units,self.coordinates);candidates=[]
        for label,b in self.baselines.items():
            t=b['threshold'];u=b['matrixUnits'];coords=(t['frequencyCoordinate'],t['radicalCoordinate'])
            for j in range(t['coefficientMatrix'].rows):
                key=self.row_key(t['coefficientMatrix'],j,u,coords)
                if equal(key,requested):candidates.append({'input':key,'value':([t['clearedMatrix'][j,k] for k in range(t['clearedMatrix'].cols)],t['rowDenominators'][j]),
                    'owner':(label,j),'path':b['thresholdPath'],'sha256':b['thresholdSha256']})
        f.atomic_pickle(folder/'input.pickle',{'requested':requested,'matrix':matrix,'row':i,'candidates':candidates,
            'wholeOwnInputReference':{'path':str(self.base/'remaining-case-end-threshold-inputs.pickle'),'sha256':f.digest(self.base/'remaining-case-end-threshold-inputs.pickle')}})
        self.row_start=len(self.calls)
        if not candidates:return False
        value=candidates[0]['value']
        for other in candidates[1:]:same(value,other['value'])
        f.atomic_pickle(folder/'value.pickle',value)
        receipt={'row':i,'reused':True,'owner':candidates[0]['owner'],'inputSha256':f.digest(folder/'input.pickle'),
            'valueSha256':f.digest(folder/'value.pickle'),'operationIndices':[]}
        f.save(folder/'completed.json',receipt);self.rows.append(receipt)
        rows.append(value[0]);denominators.append(value[1]);return True

    def row_complete(self,i,values,denominator):
        f.require(i==self.row,'actual row completion owner')
        folder=self.current_folder;f.atomic_pickle(folder/'value.pickle',(values,denominator))
        receipt={'row':i,'reused':False,'owner':(OWNER,i),'inputSha256':f.digest(folder/'input.pickle'),
            'valueSha256':f.digest(folder/'value.pickle'),'operationIndices':list(range(self.row_start,len(self.calls)))}
        f.save(folder/'completed.json',receipt);self.rows.append(receipt)

    def restore_tail(self,matrix,denominators):
        self.row=None;self.current_folder=self.base/'determinant-tail';self.current_folder.mkdir()
        requested={'polynomialMatrix':matrix,'rowDenominators':denominators,'inputUnits':self.units,'coefficientCoordinates':self.coordinates}
        candidates=[];matrix_only=[]
        for label,b in self.baselines.items():
            t=b['threshold'];key={'polynomialMatrix':t['clearedMatrix'],'rowDenominators':list(t['rowDenominators']),
                'inputUnits':b['matrixUnits'],'coefficientCoordinates':(t['frequencyCoordinate'],t['radicalCoordinate'])}
            full={'input':key,'value':((t['numerator'],t['denominator']),t['clearedMatrix'],t['rowDenominators']),
                'owner':label,'path':b['thresholdPath'],'sha256':b['thresholdSha256']}
            if equal(requested,key):candidates.append(full)
            if equal(matrix,t['clearedMatrix']) and equal(self.units,b['matrixUnits']) and equal(self.coordinates,key['coefficientCoordinates']):matrix_only.append(full)
        f.atomic_pickle(self.current_folder/'input.pickle',{'requested':requested,'candidates':candidates,'matrixOnlyCandidates':matrix_only})
        self.matrix_only=matrix_only;self.tail_start=len(self.calls);self.tail_reused=bool(candidates)
        if not candidates:return False
        self.saved_tail=candidates[0]['value']
        for other in candidates[1:]:same(other['value'],self.saved_tail)
        self.tail_owner=candidates[0]['owner'];return True

    def observe(self,name,value):
        folder=self.current_folder if name in ('entries','denominator','clearedRow') else self.base/'determinant-values'
        folder.mkdir(parents=True,exist_ok=True);f.atomic_pickle(folder/(name+'.pickle'),value)

    def call(self,site,name,*args,**kwargs):
        index=len(self.calls);folder=self.base/'determinant-operations'/str(index);folder.mkdir(parents=True)
        key={'function':name,'args':args,'kwargs':kwargs,'unitContext':self.unit_context}
        matches=[v for v in self.atlas if equal(v['input'],key)]
        matrix_candidates=[]
        if name=='ImmutableMatrix':
            for label,b in self.baselines.items():
                t=b['threshold'];rows=[[t['clearedMatrix'][i,j] for j in range(t['clearedMatrix'].cols)] for i in range(t['clearedMatrix'].rows)]
                if equal(args,(rows,)) and not kwargs and equal(self.units,b['matrixUnits']) and equal(self.coordinates,(t['frequencyCoordinate'],t['radicalCoordinate'])):
                    matrix_candidates.append({'input':key,'value':t['clearedMatrix'],'owner':{'kind':'accepted-cleared-matrix','end':label,
                        'path':b['thresholdPath'],'sha256':b['thresholdSha256']}})
            matches+=matrix_candidates
        f.atomic_pickle(folder/'input.pickle',{'requested':key,'matches':matches,'site':site,'row':self.row,
            'nativeCaller':'FullPencilModes.rational_determinant','sourceOwner':OWNER})
        if matches:
            value=matches[0]['value'];owner=matches[0]['owner']
            for other in matches[1:]:same(value,other['value'])
        else:
            if name=='det':
                f.require(type(args[0]).det is self.functions['det'],'actual native immutable matrix method dispatch')
                # The historical complete tail is reusable, but no standalone
                # determinant value was serialized by those old producers.
                # Never recover it by algebraic reconstruction or repeat it.
                if self.matrix_only:
                    f.atomic_pickle(folder/'unavailable-completed-determinant-value.pickle',{'requested':key,'completedMatrixInputs':self.matrix_only})
                    raise RuntimeError('completed standalone determinant input has no saved standalone return; preserve unfinished tail')
            value=self.functions[name](*args,**kwargs)
            owner={'kind':'new-native-determinant-operation','call':index,'inputPath':str(folder/'input.pickle'),'valuePath':str(folder/'value.pickle')}
        f.atomic_pickle(folder/'value.pickle',value)
        receipt={'index':index,'function':name,'site':site,'row':self.row,'reused':bool(matches),'owner':owner,
            'inputSha256':f.digest(folder/'input.pickle'),'valueSha256':f.digest(folder/'value.pickle')}
        f.save(folder/'completed.json',receipt);self.calls.append(receipt);self.atlas.append({'input':key,'value':value,'owner':owner,'receipt':receipt})
        return value

    def finish(self,value):
        f.atomic_pickle(self.base/'native-determinant-return.pickle',value)
        f.save(self.base/'determinant-tail/completed.json',{'reused':self.tail_reused,'owner':self.tail_owner if self.tail_reused else OWNER,
            'inputSha256':f.digest(self.base/'determinant-tail/input.pickle'),'valueSha256':f.digest(self.base/'native-determinant-return.pickle'),
            'operationIndices':list(range(self.tail_start,len(self.calls)))})
        return value


def construct(base,cp,native,functions):
    own=f.unpickle(base/'remaining-case-end-threshold-inputs.pickle');rational=f.unpickle(base/'remaining-case-end-rational.pickle')
    same(rational['sourceInput'],own);same(rational['table']['source'],own['source']);same(own['owner'],OWNER)
    baseline_packet=f.unpickle(base/'saved-native-threshold-input-values.pickle');baselines=baseline_packet['baselines']
    same(baseline_packet['own'],own['fullBindingInput'])
    pairs=f.unpickle(base/'threshold-native-call-input-pairs.pickle');units=own['fullBindingInput']['units']['PENCIL_PLUS']
    requested={'coefficientMatrix':own['coefficientMatrix'],'inputUnits':units,'frequencyCoordinate':own['frequencyCoordinate'],'radicalCoordinate':own['radicalCoordinate']}
    same(pairs['requestedDeterminant'],requested);f.require(not any(v['exactInput'] for v in pairs['determinantCandidates']),'accepted missing full determinant call')
    # The original accepted diagnostic saved the final fraction/cleared rows,
    # not a standalone determinant. Read its actual full record as well.
    dp=base/'original-determinant-diagnostic-elimination.pickle';diagnostic=f.unpickle(dp);reference=baselines['REFERENCE']['threshold']
    for key in ('coefficientMatrix','clearedMatrix','rowDenominators','numerator','denominator','wave','elimination'):same(diagnostic[key],reference[key])
    same(diagnostic['coordinate'],reference['frequencyCoordinate']);same(diagnostic['radicalCoordinate'],reference['radicalCoordinate'])
    f.require('determinant' not in diagnostic and all('determinant' not in b['threshold'] for b in baselines.values()),'historical packets retain full fraction return, not an independent determinant slot')
    f.atomic_pickle(base/'saved-native-determinant-returns.pickle',{'baselines':baselines,'diagnostic':diagnostic,'diagnosticPath':str(dp),'diagnosticSha256':f.digest(dp),
        'standaloneDeterminantValueSaved':False,'wholeOwnInput':requested})
    f.atomic_pickle(base/'native-determinant-input.pickle',{'requested':requested,'fullOwnInputPath':str(base/'remaining-case-end-threshold-inputs.pickle'),
        'fullOwnInputSha256':f.digest(base/'remaining-case-end-threshold-inputs.pickle'),'acceptedRationalTablePath':str(base/'new-rational-end/right-rational-end.pickle'),
        'acceptedRationalTableSha256':f.digest(base/'new-rational-end/right-rational-end.pickle'),'coordinateMap':own['coordinateMap'],'owner':OWNER})
    router=Router(base,own,baselines,functions);value=native(own['coefficientMatrix'],router)
    (numerator,denominator),cleared,row_denominators=value
    output={'coefficientMatrix':own['coefficientMatrix'],'wave':own['wave'],'numerator':numerator,'denominator':denominator,
        'clearedMatrix':cleared,'rowDenominators':row_denominators,'frequencyCoordinate':own['frequencyCoordinate'],
        'radicalCoordinate':own['radicalCoordinate'],'coordinateMap':own['coordinateMap'],'owner':OWNER}
    f.atomic_pickle(base/'new-end-determinant.pickle',output)
    f.require(cleared.shape==(5,5) and len(row_denominators)==5 and denominator!=0,'actual full native determinant return')
    # Route future calls using only actual completed determinant outputs.
    # No resultant, Poly, threshold analysis or root operation executes here.
    elimination_input={'numerator':numerator,'wave':own['wave'],'radicalCoordinate':own['radicalCoordinate'],'coefficientUnit':(0,0,0)}
    em=[];am=[];bm=[]
    for label,b in baselines.items():
        t=b['threshold'];old={'numerator':t['numerator'],'wave':t['wave'],'radicalCoordinate':t['radicalCoordinate'],'coefficientUnit':(0,0,0)}
        if equal(old,elimination_input):em.append({'input':old,'value':t['elimination'],'owner':label,'path':b['thresholdPath'],'sha256':b['thresholdSha256']})
        branch_input={'wave':t['wave'],'radicalCoordinate':t['radicalCoordinate'],'frequencyCoordinate':t['frequencyCoordinate'],'unit':(0,0,0)}
        if equal(branch_input,pairs['requestedBranches']):bm.append({'input':branch_input,'value':t['bulkBranchFrequencies'],'owner':label,'path':b['thresholdPath'],'sha256':b['thresholdSha256']})
    f.atomic_pickle(base/'determinant-elimination-call-pairs.pickle',{'requested':elimination_input,'matches':em})
    f.atomic_pickle(base/'determinant-branch-call-pairs.pickle',{'requested':pairs['requestedBranches'],'matches':bm})
    for x in em[1:]:same(x['value'],em[0]['value'])
    for x in bm[1:]:same(x['value'],bm[0]['value'])
    pending={'elimination':[] if em else [elimination_input],'analysis':[],'branches':[] if bm else [pairs['requestedBranches']]}
    if em:
        key={'elimination':em[0]['value'],'coordinate':own['frequencyCoordinate'],'analysisSource':previous.previous.source_body_from_join(base,'acceptedRealAxisAnalysis'),
            'callerSource':previous.previous.source_body_from_join(base,'acceptedThresholdAdapter'),'polynomialConstructor':'sp.Poly(elimination, coordinate)'}
        for label,b in baselines.items():
            t=b['threshold'];old=dict(key,elimination=t['elimination'],coordinate=t['frequencyCoordinate'])
            if equal(key,old):am.append({'input':old,'value':t['realAxisAnalysis'],'owner':label,'path':b['thresholdPath'],'sha256':b['thresholdSha256']})
        f.atomic_pickle(base/'determinant-analysis-call-pairs.pickle',{'requested':key,'matches':am})
        for x in am[1:]:same(x['value'],am[0]['value'])
        if not am:pending['analysis'].append(key)
    f.atomic_pickle(base/'pending-end-threshold-operations.pickle',pending)
    routes={'eliminationCandidateOwners':[v['owner'] for v in em],'analysisCandidateOwners':[v['owner'] for v in am],
        'branchCandidateOwners':[v['owner'] for v in bm],'analysisDeferredUntilElimination':not em,'newEliminationAnalysisRootCalls':0}
    f.save(base/'determinant-threshold-call-routing.json',routes)
    aliases={}
    for label,case in own['aliases'].items():
        aliases[label]={}
        for end,prior in case.items():
            shared=tuple(prior['owner'])==OWNER
            item={'address':(label,end),'owner':prior['owner'],'fullInputRoute':prior,
                'determinantPath':str(base/'new-end-determinant.pickle') if shared else prior['thresholdPath'],
                'mode':'shared-new-native-determinant' if shared else 'accepted-baseline-threshold-return',
                'endDomainContinuationOrNumericalReuseAccepted':False}
            folder=base/'determinant-cases'/label/end.lower();folder.mkdir(parents=True);f.save(folder/'route.json',item);aliases[label][end]=item
        f.save(base/'determinant-cases'/label/'case-summary.json',aliases[label])
    changed=sp.MutableDenseMatrix(own['coefficientMatrix']);changed[0,0]+=1;changed=sp.ImmutableMatrix(changed)
    wrong_unit=tuple(tuple((u[0]+1,u[1],u[2]) if (i,j)==(0,0) else u for j,u in enumerate(row)) for i,row in enumerate(units))
    mutations={'actualMatrix':own['coefficientMatrix'],'changedMatrix':changed,'actualUnits':units,'changedUnits':wrong_unit,
        'actualOwner':OWNER,'changedOwner':(OWNER[0],'LEFT'),'actualNumerator':numerator,'changedNumerator':numerator+1}
    f.atomic_pickle(base/'determinant-mutation-operands.pickle',mutations)
    controls={'matrixCoefficient':not equal(mutations['actualMatrix'],changed),'physicalUnit':not equal(units,wrong_unit),
        'owner':OWNER!=mutations['changedOwner'],'returnedNumerator':not equal(numerator,mutations['changedNumerator'])}
    f.save(base/'determinant-mutation-controls.json',controls);f.require(all(controls.values()),'actual determinant route/value mutations')
    counts={name:{'new':sum(c['function']==name and not c['reused'] for c in router.calls),'reused':sum(c['function']==name and c['reused'] for c in router.calls)} for name in functions}
    f.save(base/'determinant-operation-inventory.json',{'rows':router.rows,'calls':router.calls,'counts':counts,'tailReused':router.tail_reused})
    f.atomic_pickle(base/'remaining-case-end-determinant.pickle',{'owner':OWNER,'sourceInput':own,'determinant':output,'rows':router.rows,
        'operations':router.calls,'aliases':aliases,'thresholdRoutes':routes,'pending':pending,'thresholdsAndEndDomainsUnaccepted':True})
    return {'physicalEnds':8,'sharedNewUses':2,'baselineAliases':6,'newEndOwners':1,'cases':aliases,
        'newNativeRowCalls':sum(not r['reused'] for r in router.rows),'reusedNativeRowCalls':sum(r['reused'] for r in router.rows),
        'tailReused':router.tail_reused,'operationCounts':counts,'thresholdRoutes':routes,'actualMutationControls':len(controls),
        'newEliminationAnalysisBindingDerivativeModeCurrentCalls':0,'thresholdDomainContinuationOrNumericalReuseAccepted':False}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    previous.previous.previous.previous.previous.protect(base)
    cp,manifest=load(base);native,joins=joined_sources(base,manifest);functions=prohibit();result=construct(base,cp,native,functions)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen source posthash')
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'reference postidentity')
    for key,item in joins.items():
        if key.startswith('sympy_') or key=='matrixDet':f.require(f.digest(Path(item['path']))==item['sha256'],'native operation source posthash')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_DETERMINANT','nativeJoins':joins,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
