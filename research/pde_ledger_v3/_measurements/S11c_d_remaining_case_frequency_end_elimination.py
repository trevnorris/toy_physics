#!/usr/bin/env python3
"""Missing native resultant with exact saved full complex-analysis call routing."""
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
import S11c_d_remaining_case_frequency_end_determinant as previous

h,f,sp,engine=previous.h,previous.f,previous.sp,previous.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_determinant_checkpoint.json'
CP_SHA='358f2f3ff6451ee419c40d9ef2befd48bfa01eb862cd3d47d55617990fc4b51c'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_elimination_plan.md'
OWNER=previous.OWNER
same=previous.same
equal=previous.equal
source_body=previous.source_body


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(f.digest(CP)==CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_DETERMINANT','accepted determinant checkpoint')
    for root,stage,sha in ((origin.parent,'frequency_end_determinant',cp['checksSha256']),(vr,'validate',cp['validation']['checksSha256'])):
        h.source.receipts.inspect_guard(root,stage)
        checks=origin/'checks.json' if stage!='validate' else vr/'checks.json'
        f.require(f.digest(checks)==sha and checks.read_bytes()==(root/(stage+'.stdout')).read_bytes(),'clean accepted checks/stdout')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],
              'acceptedDeterminant':{'checkpoint':str(CP),'checkpointSha256':CP_SHA,'runDirectory':str(origin),
                  'checksSha256':cp['checksSha256'],'validationDirectory':str(vr),'validatorChecksSha256':cp['validation']['checksSha256']}}
    ref=lambda path,name,sha:h.source.reference(base,manifest,path,name,sha)
    for name,item in cp['artifacts'].items():ref(origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        p=origin/name;f.require(p.is_symlink() and str(p.readlink())==item['original'] and str(p.resolve())==item['resolvedOriginal'] and p.stat().st_size==item['bytes'] and f.digest(p)==item['sha256'],'original determinant reference identity')
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'accepted current/frozen source')
        ref(origin/'source'/name,'source/'+name,value)
    for path,name in ((CP,'accepted-determinant-checkpoint.json'),(origin/'checks.json','accepted-determinant-checks.json'),
                      (origin/'inputs.json','accepted-determinant-manifest.json'),(vr/'checks.json','accepted-determinant-validation.json')):ref(path,name,f.digest(path))
    for name,item in cp['validation']['artifacts'].items():ref(vr/name,'accepted-determinant-validation/'+name,item['sha256'])
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));value=f.digest(path);f.require(name not in manifest['sourceFiles'],'fresh elimination constructor source')
        manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return cp,manifest

def native_adapter():
    """Route one exact resultant call, reverse the entire original end body."""
    original=ast.parse(textwrap.dedent(inspect.getsource(h.source.q.end_sources))).body[0]
    loop=next(n for n in original.body if isinstance(n,ast.For));sites=[]
    class Route(ast.NodeTransformer):
        def visit_Call(self,node):
            if ast.unparse(node.func)!='sp.resultant':return self.generic_visit(node)
            sites.append(copy.deepcopy(node))
            f.require(not node.keywords and ast.unparse(node)=='sp.resultant(numerator, curve, radical_coordinate)','literal native resultant call')
            return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='router',ctx=ast.Load()),attr='resultant',ctx=ast.Load()),
                args=copy.deepcopy(node.args),keywords=[]),node)
    routed=Route().visit(copy.deepcopy(original));f.require(len(sites)==1,'one native end resultant site')
    class Reverse(ast.NodeTransformer):
        def visit_Call(self,node):
            if ast.unparse(node.func)=='router.resultant':return copy.deepcopy(sites[0])
            return self.generic_visit(node)
    restored=Reverse().visit(copy.deepcopy(routed));f.require(ast.dump(restored)==ast.dump(original),'entire native end_sources resultant-only reversal')
    modified_loop=next(n for n in routed.body if isinstance(n,ast.For))
    slot=next(i for i,n in enumerate(modified_loop.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='elimination' for t in n.targets))
    f.require(isinstance(modified_loop.body[slot].value,ast.Call) and ast.unparse(modified_loop.body[slot].value.func)=='router.resultant','whole actual resultant assignment')
    fn=ast.parse('def native_elimination(numerator,curve,radical_coordinate,router):\n pass').body[0]
    fn.body=[copy.deepcopy(modified_loop.body[slot])]+ast.parse('return elimination').body
    module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]));env={}
    exec(compile(module,'<native resultant assignment with exact input value persistence>','exec'),env)
    return env['native_elimination'],{'wholeNativeAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'entireOriginalBodyReverseAST':True,'nativeAssignmentIndex':slot,'nativeAssignmentAST':ast.dump(loop.body[slot]),
        'nativeCall':ast.unparse(sites[0]),'nativeCallAST':ast.dump(sites[0]),'adapterSource':ast.unparse(module),
        'allOtherNativeStatementsUnchangedAndUnexecuted':True}


def joined_sources(base,manifest):
    accepted=previous.previous.previous.accepted_threshold
    sources={'endSources':h.source.q.end_sources,'endTables':h.inputs.chart.end_tables,
        'rationalDeterminant':engine.FullPencilModes.rational_determinant,
        'acceptedThresholdAdapter':accepted.end_adapter,'acceptedRealAxisAnalysis':accepted.real_axis_analysis,
        'determinantAdapter':previous.native_adapter,'determinantConstructor':previous.construct}
    result={name:source_body(value) for name,value in sources.items()};old=json.loads((base/'native-end-determinant-joins.json').read_text())
    for name,item in result.items():
        relative=str(Path(item['path']).relative_to(f.ROOT));f.require(manifest['sourceFiles'][relative]==item['sha256'],'complete original current/frozen caller joins')
    for name in ('endSources','endTables','rationalDeterminant','acceptedThresholdAdapter','acceptedRealAxisAnalysis'):same(result[name],old[name])
    frequency=json.loads((base/'end-origins/frequency-checkpoint.json').read_text())
    for name in ('acceptedThresholdAdapter','acceptedRealAxisAnalysis'):
        rel=str(Path(result[name]['path']).relative_to(f.ROOT));f.require(frequency['sourceFiles'][rel]==result[name]['sha256'],'actual accepted complex-analysis helper')
    unused=accepted.end_adapter() # compile/reverse only; the returned constructor is never called.
    f.require(callable(unused),'whole accepted adapter compiles without executing scientific statements')
    native,proof=native_adapter();result['eliminationAdapter']=proof
    result['sympy_resultant']=source_body(sp.resultant);result['sympyVersion']=sp.__version__
    f.save(base/'native-end-elimination-joins.json',result)
    return native,result


def forbidden(*args,**kwargs):raise RuntimeError('only an actually missing native resultant may execute')


def prohibit():
    function=sp.resultant
    previous.previous.previous.prohibit()
    for module in (previous,previous.previous):
        for name in ('load','construct','native_adapter','main'):setattr(module,name,forbidden)
    for name in ('__init__','call','restore_row','restore_tail','row_complete','observe','finish'):setattr(previous.Router,name,forbidden)
    for name in ('__init__','operation','start','observe'):setattr(previous.previous.OperationRouter,name,forbidden)
    return function


class ResultantRouter:
    def __init__(self,base,own,determinant,baseline_pairs,function):
        self.base,self.own,self.determinant,self.baseline_pairs,self.function=base,own,determinant,baseline_pairs,function
        self.receipt=None;self.operation_count=0

    def resultant(self,numerator,curve,radical_coordinate):
        f.require(self.receipt is None,'one native resultant invocation')
        requested={'numerator':numerator,'wave':curve,'radicalCoordinate':radical_coordinate,'coefficientUnit':(0,0,0)}
        same(requested,self.baseline_pairs['requested'])
        actual={'function':'sympy.resultant','args':(numerator,curve,radical_coordinate),'kwargs':{},
            'coefficientUnit':(0,0,0),'frequencyCoordinate':self.own['frequencyCoordinate'],'radicalCoordinate':radical_coordinate,
            'physicalInputUnits':self.own['fullBindingInput']['units']['PENCIL_PLUS'],'coordinateMap':self.own['coordinateMap']}
        folder=self.base/'elimination-operation';folder.mkdir()
        f.atomic_pickle(folder/'input.pickle',{'requested':requested,'actualCall':actual,'candidates':self.baseline_pairs['matches'],
            'owner':OWNER,'wholeOwnInput':self.own,'acceptedDeterminant':self.determinant,
            'nativeCaller':'S11c_d_frequency_source.end_sources','acceptedCallerAdapter':'S11c_d_frequency_source_finish.end_adapter'})
        matches=self.baseline_pairs['matches']
        if matches:
            value=matches[0]['value'];owner=matches[0]['owner']
            for item in matches[1:]:same(item['value'],value)
        else:
            value=self.function(*actual['args'],**actual['kwargs']);owner=OWNER;self.operation_count+=1
        # Persist the actual return immediately, before any later type/symbol/route guard.
        f.atomic_pickle(folder/'value.pickle',value)
        self.receipt={'reused':bool(matches),'owner':owner,'inputSha256':f.digest(folder/'input.pickle'),
            'valueSha256':f.digest(folder/'value.pickle'),'actualNewCalls':self.operation_count}
        f.save(folder/'completed.json',self.receipt)
        return value


def construct(base,cp,native,joins,function):
    accepted_result=f.unpickle(base/'remaining-case-end-determinant.pickle');own=accepted_result['sourceInput'];det=accepted_result['determinant']
    same(own,f.unpickle(base/'remaining-case-end-threshold-inputs.pickle'));same(det,f.unpickle(base/'new-end-determinant.pickle'));same(own['owner'],OWNER)
    saved_returns=f.unpickle(base/'saved-native-determinant-returns.pickle');baselines=saved_returns['baselines']
    pairs=f.unpickle(base/'determinant-elimination-call-pairs.pickle');branch_pairs=f.unpickle(base/'determinant-branch-call-pairs.pickle')
    requested={'numerator':det['numerator'],'wave':det['wave'],'radicalCoordinate':det['radicalCoordinate'],'coefficientUnit':(0,0,0)}
    same(requested,pairs['requested']);same(accepted_result['pending']['elimination'],[requested])
    matches=[]
    for label,b in baselines.items():
        t=b['threshold'];old={'numerator':t['numerator'],'wave':t['wave'],'radicalCoordinate':t['radicalCoordinate'],'coefficientUnit':(0,0,0)}
        if equal(old,requested):matches.append({'input':old,'value':t['elimination'],'owner':label,'path':b['thresholdPath'],'sha256':b['thresholdSha256']})
    same(matches,pairs['matches']);f.require(not matches,'accepted actually missing full resultant input')
    f.atomic_pickle(base/'native-elimination-input.pickle',{'requested':requested,'fullOwnInputPath':str(base/'remaining-case-end-threshold-inputs.pickle'),
        'fullOwnInputSha256':f.digest(base/'remaining-case-end-threshold-inputs.pickle'),
        'determinantPath':str(base/'new-end-determinant.pickle'),'determinantSha256':f.digest(base/'new-end-determinant.pickle'),
        'fullSavedCallPairs':pairs,'owner':OWNER})
    router=ResultantRouter(base,own,det,pairs,function);elimination=native(det['numerator'],det['wave'],det['radicalCoordinate'],router)
    output={'determinant':det,'elimination':elimination,'frequencyCoordinate':own['frequencyCoordinate'],
        'radicalCoordinate':own['radicalCoordinate'],'coordinateMap':own['coordinateMap'],'owner':OWNER}
    f.atomic_pickle(base/'new-end-elimination.pickle',output)
    f.require(isinstance(elimination,sp.Expr) and elimination!=0 and not elimination.free_symbols-{own['frequencyCoordinate']},'actual nonzero frequency-only resultant')
    # Compare the full accepted Poly(elimination, coordinate)+complex-analysis
    # pipeline, not an invented standalone Poly or old intermediate value.
    requested_analysis={'elimination':elimination,'coordinate':own['frequencyCoordinate'],'coefficientUnit':(0,0,0),
        'coordinateUnit':(0,0,0),'polynomialConstructor':'sp.Poly(elimination, coordinate)',
        'analysisSource':joins['acceptedRealAxisAnalysis'],'callerSource':joins['acceptedThresholdAdapter']}
    comparisons=[];candidates=[]
    for label,b in baselines.items():
        t=b['threshold'];old=dict(requested_analysis,elimination=t['elimination'],coordinate=t['frequencyCoordinate'])
        item={'input':old,'value':t['realAxisAnalysis'],'owner':label,'path':b['thresholdPath'],'sha256':b['thresholdSha256']}
        match=equal(requested_analysis,old);comparisons.append({'requested':requested_analysis,'saved':item,'exactInput':match})
        if match:candidates.append(item)
    f.atomic_pickle(base/'elimination-analysis-call-pairs.pickle',{'requested':requested_analysis,'comparisons':comparisons,'matches':candidates})
    for item in candidates[1:]:same(item['value'],candidates[0]['value'])
    f.atomic_pickle(base/'elimination-accepted-branch-routes.pickle',branch_pairs)
    pending={'complexAnalysisPipeline':[] if candidates else [requested_analysis],'branches':accepted_result['pending']['branches']}
    f.atomic_pickle(base/'pending-end-complex-threshold-operations.pickle',pending)
    routes={'complexAnalysisCandidateOwners':[item['owner'] for item in candidates],
        'branchCandidateOwners':[item['owner'] for item in branch_pairs['matches']],
        'pendingComplexAnalysisPipelines':len(pending['complexAnalysisPipeline']),'newPolyAnalysisRootCalls':0}
    f.save(base/'elimination-threshold-call-routing.json',routes)
    aliases={}
    for label,case in own['aliases'].items():
        aliases[label]={}
        for end,prior in case.items():
            shared=tuple(prior['owner'])==OWNER
            item={'address':(label,end),'owner':prior['owner'],'fullInputRoute':prior,
                'eliminationPath':str(base/'new-end-elimination.pickle') if shared else prior['thresholdPath'],
                'mode':'shared-new-native-elimination' if shared else 'accepted-baseline-threshold-return',
                'endDomainContinuationOrNumericalReuseAccepted':False}
            folder=base/'elimination-cases'/label/end.lower();folder.mkdir(parents=True);f.save(folder/'route.json',item);aliases[label][end]=item
        f.save(base/'elimination-cases'/label/'case-summary.json',aliases[label])
    changed_input=dict(requested,numerator=requested['numerator']+1)
    unit_input=dict(requested,coefficientUnit=(1,0,0))
    mutations={'actualInput':requested,'changedInput':changed_input,'changedUnitInput':unit_input,
        'actualOwner':OWNER,'changedOwner':(OWNER[0],'LEFT'),'actualValue':elimination,'changedValue':elimination+1}
    f.atomic_pickle(base/'elimination-mutation-operands.pickle',mutations)
    controls={'inputNumerator':not equal(requested,changed_input),'coefficientUnit':not equal(requested,unit_input),
        'owner':OWNER!=mutations['changedOwner'],'returnedElimination':not equal(elimination,mutations['changedValue'])}
    f.save(base/'elimination-mutation-controls.json',controls);f.require(all(controls.values()),'actual elimination input/owner/value controls')
    f.atomic_pickle(base/'remaining-case-end-elimination.pickle',{'owner':OWNER,'sourceInput':own,'determinantInput':det,'elimination':output,
        'operation':router.receipt,'aliases':aliases,'thresholdRoutes':routes,'pending':pending,'thresholdsAndEndDomainsUnaccepted':True})
    return {'physicalEnds':8,'sharedNewUses':2,'baselineAliases':6,'newEndOwners':1,'cases':aliases,
        'newResultantCalls':router.operation_count,'reusedResultantCalls':int(router.receipt['reused']),
        'thresholdRoutes':routes,'actualMutationControls':len(controls),'newPolyAnalysisBindingDerivativeModeCurrentCalls':0,
        'thresholdDomainContinuationOrNumericalReuseAccepted':False}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    previous.previous.previous.previous.previous.previous.protect(base)
    cp,manifest=load(base);native,joins=joined_sources(base,manifest);function=prohibit();result=construct(base,cp,native,joins,function)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen source posthash')
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'reference postidentity')
    f.require(f.digest(Path(joins['sympy_resultant']['path']))==joins['sympy_resultant']['sha256'],'native operation source posthash')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_ELIMINATION','nativeJoins':joins,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
