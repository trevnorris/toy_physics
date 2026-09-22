#!/usr/bin/env python3
"""Missing accepted complex threshold analysis with durable native operations."""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import operator
import resource
import signal
import textwrap
import time
from types import SimpleNamespace
import S11c_d_remaining_case_frequency_end_elimination as previous

h,f,sp,engine=previous.h,previous.f,previous.sp,previous.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_elimination_checkpoint.json'
CP_SHA='64bdd5aa0d56e5572815e8d314b059aab89d4547b3e5c2c658004496d53618fe'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_complex_threshold_plan.md'
OWNER=previous.OWNER
accepted=previous.previous.previous.previous.accepted_threshold
same=previous.same
equal=previous.equal
source_body=previous.source_body


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(f.digest(CP)==CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_ELIMINATION','accepted elimination checkpoint')
    for root,stage,sha in ((origin.parent,'frequency_end_elimination',cp['checksSha256']),(vr,'validate',cp['validation']['checksSha256'])):
        h.source.receipts.inspect_guard(root,stage)
        checks=origin/'checks.json' if stage!='validate' else vr/'checks.json'
        f.require(f.digest(checks)==sha and checks.read_bytes()==(root/(stage+'.stdout')).read_bytes(),'clean accepted checks/stdout')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],
              'acceptedElimination':{'checkpoint':str(CP),'checkpointSha256':CP_SHA,'runDirectory':str(origin),
                  'checksSha256':cp['checksSha256'],'validationDirectory':str(vr),'validatorChecksSha256':cp['validation']['checksSha256']}}
    ref=lambda path,name,sha:h.source.reference(base,manifest,path,name,sha)
    for name,item in cp['artifacts'].items():ref(origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        p=origin/name;f.require(p.is_symlink() and str(p.readlink())==item['original'] and str(p.resolve())==item['resolvedOriginal'] and p.stat().st_size==item['bytes'] and f.digest(p)==item['sha256'],'original elimination reference identity')
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'accepted current/frozen source')
        ref(origin/'source'/name,'source/'+name,value)
    for path,name in ((CP,'accepted-elimination-checkpoint.json'),(origin/'checks.json','accepted-elimination-checks.json'),
                      (origin/'inputs.json','accepted-elimination-manifest.json'),(vr/'checks.json','accepted-elimination-validation.json')):ref(path,name,f.digest(path))
    for name,item in cp['validation']['artifacts'].items():ref(vr/name,'accepted-elimination-validation/'+name,item['sha256'])
    focused_cp=f.M/'S11c_d_frequency_source_threshold_acceptance.json'
    focused=json.loads(focused_cp.read_text());focused_root=Path(focused['focusedRun'])
    f.require(focused['status']=='ACCEPTED_COMPLEX_THRESHOLD_REPAIR' and focused['exitCode']==0 and focused['stderrBytes']==0 and focused['checksStdoutIdentical'],'original accepted focused threshold evidence')
    ref(focused_cp,'complex-original-focused/acceptance.json',f.digest(focused_cp))
    focused_checks=focused_root/'checks.json';f.require(f.digest(focused_checks)==focused['checksSha256'],'original focused checks hash')
    ref(focused_checks,'complex-original-focused/checks.json',focused['checksSha256'])
    for name in ('actual-analysis.pickle','actual-imaginary-mutation.pickle','focused-result.pickle'):
        ref(focused_root/name,'complex-original-focused/'+name,focused['checks']['artifacts'][name])
    for name in ('_measurements/S11c_d_frequency_source_finish.py',):
        value=focused['checks']['sourceFiles'][name];f.require(f.digest(f.ROOT/name)==manifest['sourceFiles'][name]==value,'original focused whole native caller source')
    focused_helper=f.M/'S11c_d_frequency_source_threshold_focused.py'
    f.require(f.digest(focused_helper)==focused['focusedScriptSha256'],'actual original focused caller source')
    ref(focused_helper,'complex-original-focused/focused.py',focused['focusedScriptSha256'])
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));value=f.digest(path);f.require(name not in manifest['sourceFiles'],'fresh complex threshold constructor source')
        manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return cp,manifest

def instrument(original,prefix):
    """Reverse every operation route and observation to the whole original AST."""
    sites={};observations=[]
    def tag(node,kind):
        name=prefix+str(len(sites));sites[name]={'kind':kind,'source':ast.unparse(node),'ast':ast.dump(node)};return name
    class Route(ast.NodeTransformer):
        def visit_Call(self,node):
            old=copy.deepcopy(node);node=self.generic_visit(node)
            if not isinstance(old.func,ast.Attribute) or ast.unparse(old.func.value)=='f':return node
            name=tag(old,'call');global_call=ast.unparse(old.func.value)=='sp'
            receiver=ast.Constant(None) if global_call else node.func.value
            f.require(all(k.arg is not None for k in node.keywords),'explicit native keyword arguments')
            return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='router',ctx=ast.Load()),attr='call',ctx=ast.Load()),
                args=[ast.Constant(name),ast.Constant(('sp.' if global_call else 'method.')+old.func.attr),receiver,
                    ast.Tuple(elts=node.args,ctx=ast.Load()),ast.Dict(keys=[ast.Constant(k.arg) for k in node.keywords],values=[k.value for k in node.keywords])],keywords=[]),node)
        def visit_BinOp(self,node):
            old=copy.deepcopy(node);node=self.generic_visit(node);name=tag(old,'arithmetic')
            return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='router',ctx=ast.Load()),attr='arithmetic',ctx=ast.Load()),
                args=[ast.Constant(name),ast.Constant(type(old.op).__name__),node.left,node.right],keywords=[]),node)
    routed=Route().visit(copy.deepcopy(original))
    def observe(body):
        out=[]
        for statement in body:
            if isinstance(statement,(ast.If,ast.For)):
                statement.body=observe(statement.body);statement.orelse=observe(statement.orelse)
            out.append(statement)
            if isinstance(statement,ast.Assign):
                names=[n.id for target in statement.targets for n in ast.walk(target) if isinstance(n,ast.Name)]
                for name in names:
                    site=prefix+'observation'+str(len(observations));observations.append({'site':site,'name':name})
                    out.extend(ast.parse(f"router.observe('{site}', '{name}', {name})").body)
            if isinstance(statement,ast.Expr) and isinstance(statement.value,ast.Call) and ast.unparse(statement.value.func)=='f.require':
                site=prefix+'guard'+str(len(observations));observations.append({'site':site,'guard':ast.unparse(statement)})
                out.extend(ast.parse(f"router.guard('{site}')").body)
        return out
    routed.body=observe(routed.body)
    class Reverse(ast.NodeTransformer):
        def visit_Expr(self,node):
            if isinstance(node.value,ast.Call) and ast.unparse(node.value.func) in ('router.observe','router.guard'):return None
            return self.generic_visit(node)
        def visit_Call(self,node):
            if ast.unparse(node.func) in ('router.call','router.arithmetic'):
                return ast.parse(sites[node.args[0].value]['source'],mode='eval').body
            return self.generic_visit(node)
    f.require(ast.dump(Reverse().visit(copy.deepcopy(routed)))==ast.dump(original),'whole native operation/persistence reversal')
    routed.args.args.append(ast.arg(arg='router'))
    module=ast.fix_missing_locations(ast.Module(body=[routed],type_ignores=[]));env={'sp':sp,'f':f}
    exec(compile(module,'<whole saved complex-analysis operation adapter>','exec'),env)
    return env[routed.name],{'entireOriginalBodyReverseAST':True,'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'sites':sites,'observations':observations,'adapterSource':ast.unparse(module)}


def native_adapter():
    original=ast.parse(textwrap.dedent(inspect.getsource(accepted.real_axis_analysis))).body[0]
    analysis,proof=instrument(original,'analysis-')
    # The original full end body is retained; extract only its initial Poly and
    # nonzero guard plus its threshold metadata assignment. All earlier science
    # is supplied as accepted values, and the obsolete real-only block is unused.
    whole=ast.parse(textwrap.dedent(inspect.getsource(h.source.q.end_sources))).body[0]
    loop=next(n for n in whole.body if isinstance(n,ast.For))
    index=next(i for i,n in enumerate(loop.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='poly' for t in n.targets))
    initial=ast.parse('def initial_poly(elimination,coordinate):\n pass').body[0]
    initial.body=copy.deepcopy(loop.body[index:index+2])+ast.parse('return poly').body
    poly,poly_proof=instrument(initial,'initial-')
    threshold=next(n for n in loop.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='threshold' for t in n.targets))
    caller=ast.parse('def threshold_packet(zero,cleared,row_denominators,numerator,denominator,curve,elimination,normalized,squarefree,factorization,reconstruction,intervals,branches,coordinate,radical_coordinate,inp,data,analysis):\n pass').body[0]
    caller.body=[copy.deepcopy(threshold)]+ast.parse("threshold['realAxisAnalysis']=analysis\nreturn threshold").body
    packet,packet_proof=instrument(caller,'packet-')
    restored=copy.deepcopy(whole);restored_loop=next(n for n in restored.body if isinstance(n,ast.For))
    restored_loop.body[index:index+2]=copy.deepcopy(initial.body[:-1])
    f.require(ast.dump(restored)==ast.dump(whole),'full original end body with exact initial Poly/guard slots')
    f.require(ast.unparse(initial.body[0])=='poly = sp.Poly(elimination, coordinate)' and ast.dump(caller.body[0])==ast.dump(threshold),'exact original Poly and threshold packet assignment')
    return (poly,analysis,packet),{'wholeEndAST':hashlib.sha256(ast.dump(whole).encode()).hexdigest(),
        'allEarlierEndStatementsUnchangedAndUnexecuted':True,'initial':poly_proof,'analysis':proof,'packet':packet_proof,
        'onlyAcceptedComplexAnalysisUsed':True,'thresholdAddition':"threshold['realAxisAnalysis']=analysis"}


def joined_sources(base,manifest):
    result={name:source_body(value) for name,value in {
        'endSources':h.source.q.end_sources,'endTables':h.inputs.chart.end_tables,
        'rationalDeterminant':engine.FullPencilModes.rational_determinant,
        'acceptedThresholdAdapter':accepted.end_adapter,'acceptedRealAxisAnalysis':accepted.real_axis_analysis,
        'eliminationAdapter':previous.native_adapter,'eliminationConstructor':previous.construct}.items()}
    old=json.loads((base/'native-end-elimination-joins.json').read_text())
    for name,item in result.items():
        rel=str(Path(item['path']).relative_to(f.ROOT));f.require(manifest['sourceFiles'][rel]==item['sha256'],'full current/frozen helper source')
    for name in ('endSources','endTables','rationalDeterminant','acceptedThresholdAdapter','acceptedRealAxisAnalysis'):same(result[name],old[name])
    frequency=json.loads((base/'end-origins/frequency-checkpoint.json').read_text())
    for name in ('acceptedThresholdAdapter','acceptedRealAxisAnalysis'):
        rel=str(Path(result[name]['path']).relative_to(f.ROOT));f.require(frequency['sourceFiles'][rel]==result[name]['sha256'],'actual accepted complex analysis implementation')
    f.require(callable(accepted.end_adapter()),'accepted whole end adapter compile/reversal only')
    functions,proof=native_adapter();result['complexAdapter']=proof
    names=('Poly','re','im','expand','gcd','div','gcdex','Rational','factor_list','prod')
    dispatch={name:getattr(sp,name) for name in names}
    result['sympyFunctions']={name:source_body(value) for name,value in dispatch.items()};result['sympyVersion']=sp.__version__
    f.save(base/'native-end-complex-threshold-joins.json',result)
    return functions,result,dispatch


def forbidden(*args,**kwargs):raise RuntimeError('only the actually missing accepted complex-analysis pipeline may execute')


def prohibit(dispatch):
    previous.prohibit()
    # SymPy algorithms use these public helpers internally. Restore only the
    # required exact originals; native dispatch is instrumented, prior producers
    # and every unrelated scientific entry point remain disabled.
    for name,value in dispatch.items():setattr(sp,name,value)
    for module in (previous,previous.previous,accepted):
        for name in ('load','construct','main','native_adapter','end_adapter','real_axis_analysis','reuse'):
            if hasattr(module,name):setattr(module,name,forbidden)


class AnalysisRouter:
    def __init__(self,base,own,dispatch,baselines,joins):
        self.base,self.own,self.dispatch,self.joins=base,own,dispatch,joins
        self.context={'coefficientUnit':(0,0,0),'coordinateUnit':(0,0,0),'coordinate':own['frequencyCoordinate']}
        self.atlas=[];self.missing=[];self.receipts=[];self.observation_counts={};self.native_sources={}
        for label,item in baselines.items():
            t=item['threshold'];a=t['realAxisAnalysis'];w=t['frequencyCoordinate']
            context=dict(self.context,coordinate=w);origin={'owner':label,'path':item['thresholdPath'],'sha256':item['thresholdSha256']}
            def add(name,receiver,args,kwargs,value,field):
                self.atlas.append({'input':self.key(name,receiver,args,kwargs,context),'value':value,'origin':dict(origin,field=field)})
            add('sp.Poly',None,(t['normalizedElimination'],w),{'domain':sp.QQ_I},a['normalized'],'realAxisAnalysis.normalized')
            add('method.as_expr',a['normalized'],(),{},t['normalizedElimination'],'normalizedElimination')
            add('method.sqf_part',a['normalized'],(),{},a['complexSquarefree'],'realAxisAnalysis.complexSquarefree')
            add('method.as_expr',a['complexSquarefree'],(),{},t['squarefreeElimination'],'squarefreeElimination')
            add('sp.factor_list',None,(t['normalizedElimination'],w),{'extension':sp.I},a['factorization'],'realAxisAnalysis.factorization')
            # These actual constructor inputs are saved, but their intermediate
            # Poly returns are not. An exact collision must stop, never replay.
            for expression,kwargs,field in ((t['elimination'],{},'initialPoly'),(a['realPart'],{'domain':sp.QQ},'realPoly'),(a['imaginaryPart'],{'domain':sp.QQ},'imaginaryPoly')):
                self.missing.append({'input':self.key('sp.Poly',None,(expression,w),kwargs,context),'origin':dict(origin,field=field)})
        # Inspect the other actually completed original diagnostic outputs too.
        # Their raw normalized/squarefree Poly objects are saved. Do not invent
        # the unsaved pure-real regression result or mutation's initial Poly.
        for name in ('actual-analysis.pickle','actual-imaginary-mutation.pickle'):
            path=base/'complex-original-focused'/name;a=f.unpickle(path)
            origin={'owner':'original-focused','path':str(path),'sha256':f.digest(path),'field':'complexSquarefree'}
            self.atlas.append({'input':self.key('method.sqf_part',a['normalized'],(),{},self.context),'value':a['complexSquarefree'],'origin':origin})
        f.atomic_pickle(base/'completed-complex-operation-atlas.pickle',{'available':self.atlas,'unavailableIntermediateReturns':self.missing,
            'scope':'Only actual saved typed full-call inputs/results; no old Poly or other intermediate constructed from expression metadata.'})

    @staticmethod
    def key(name,receiver,args,kwargs,context):
        return {'function':name,'receiver':receiver,'args':args,'kwargs':kwargs,'coefficientContext':context}

    def operation(self,site,name,receiver,args,kwargs,function):
        key=self.key(name,receiver,args,kwargs,self.context)
        matches=[item for item in self.atlas if equal(item['input'],key)]
        unavailable=[item for item in self.missing if equal(item['input'],key)]
        folder=self.base/'complex-operations'/str(len(self.receipts));folder.mkdir(parents=True)
        f.atomic_pickle(folder/'input.pickle',{'site':site,'input':key,'matches':matches,'unavailableMatches':unavailable,
            'owner':OWNER,'wholePipelineInputPath':str(self.base/'native-complex-threshold-input.pickle'),
            'wholePipelineInputSha256':f.digest(self.base/'native-complex-threshold-input.pickle')})
        if unavailable and not matches:
            f.save(folder/'unavailable-completed-value.json',{'site':site,'function':name,'actualSavedOrigins':[x['origin'] for x in unavailable]})
            raise RuntimeError('completed exact native call has no saved intermediate return; preserve input and inspect without replay')
        if matches:
            value=matches[0]['value'];origin=matches[0]['origin']
            for item in matches[1:]:same(value,item['value'])
        else:
            value=function();origin={'owner':OWNER,'path':str(folder/'value.pickle'),'call':len(self.receipts)}
        f.atomic_pickle(folder/'value.pickle',value)
        receipt={'site':site,'function':name,'reused':bool(matches),'owner':origin,
            'inputSha256':f.digest(folder/'input.pickle'),'valueSha256':f.digest(folder/'value.pickle')}
        f.save(folder/'completed.json',receipt);self.receipts.append(receipt)
        self.atlas.append({'input':key,'value':value,'origin':dict(origin,receiptPath=str(folder/'completed.json'),receiptSha256=f.digest(folder/'completed.json'))})
        return value

    def call(self,site,name,receiver,args,kwargs):
        if name=='sp.prod':
            # Consume the native generator once, in native order. Routed power
            # operations persist each factor before this product's input saves.
            args=(tuple(args[0]),)+args[1:]
        if name.startswith('sp.'):
            function=self.dispatch[name[3:]]
        else:
            function=getattr(receiver,name[7:])
            source=source_body(function);source_key=type(receiver).__module__+'.'+type(receiver).__qualname__+'.'+name[7:]
            if source_key not in self.native_sources:
                self.native_sources[source_key]=source;f.save(self.base/'complex-method-source'/('method-'+str(len(self.native_sources))+'.json'),{'method':source_key,**source})
            else:same(source,self.native_sources[source_key])
        return self.operation(site,name,receiver,args,kwargs,lambda:function(*args,**kwargs))

    def arithmetic(self,site,name,left,right):
        function={'Add':operator.add,'Sub':operator.sub,'Mult':operator.mul,'Div':operator.truediv,'Pow':operator.pow}[name]
        return self.operation(site,'operator.'+name,None,(left,right),{},lambda:function(left,right))

    def observe(self,site,name,value):
        count=self.observation_counts.get(site,0);self.observation_counts[site]=count+1
        folder=self.base/'complex-locals'/site/str(count);folder.mkdir(parents=True)
        f.atomic_pickle(folder/'value.pickle',value)
        f.save(folder/'completed.json',{'site':site,'name':name,'valueSha256':f.digest(folder/'value.pickle')})

    def guard(self,site):
        count=self.observation_counts.get(site,0);self.observation_counts[site]=count+1
        f.save(self.base/'complex-guards'/(site+'-'+str(count)+'.json'),{'passed':True,'site':site,'completedOperationCount':len(self.receipts)})


def construct(base,cp,functions,joins,dispatch):
    previous_result=f.unpickle(base/'remaining-case-end-elimination.pickle');own=previous_result['sourceInput']
    elimination=previous_result['elimination'];same(elimination,f.unpickle(base/'new-end-elimination.pickle'))
    same(own,f.unpickle(base/'remaining-case-end-threshold-inputs.pickle'));same(own['owner'],OWNER)
    pairs=f.unpickle(base/'elimination-analysis-call-pairs.pickle');same(previous_result['pending']['complexAnalysisPipeline'],[pairs['requested']])
    f.require(not pairs['matches'],'actual missing whole Poly-plus-complex-analysis pipeline')
    baselines=f.unpickle(base/'saved-native-determinant-returns.pickle')['baselines']
    branch_pairs=f.unpickle(base/'elimination-accepted-branch-routes.pickle')
    f.atomic_pickle(base/'native-complex-threshold-input.pickle',{'requested':pairs['requested'],'fullSavedPairs':pairs,
        'owner':OWNER,'elimination':elimination,'wholeOwnInput':own,'acceptedBranchCalls':branch_pairs})
    same(pairs['requested']['elimination'],elimination['elimination']);same(pairs['requested']['coordinate'],own['frequencyCoordinate'])
    same(pairs['requested']['analysisSource'],joins['acceptedRealAxisAnalysis']);same(pairs['requested']['callerSource'],joins['acceptedThresholdAdapter'])
    router=AnalysisRouter(base,own,dispatch,baselines,joins);initial,analyze,packet=functions
    poly=initial(elimination['elimination'],own['frequencyCoordinate'],router)
    analysis=analyze(poly,own['frequencyCoordinate'],router)
    f.atomic_pickle(base/'new-complex-analysis.pickle',analysis)
    matches=branch_pairs['matches'];f.require(matches and not previous_result['pending']['branches'],'all actual branch inputs have saved full returns')
    branches=matches[0]['value']
    for item in matches:same(item['input'],branch_pairs['requested']);same(item['value'],branches)
    f.atomic_pickle(base/'complex-accepted-branch-return.pickle',{'requested':branch_pairs['requested'],'matches':matches,'value':branches})
    det=elimination['determinant'];cm=own['coordinateMap']
    # Restore only metadata caller state read by the unchanged threshold packet.
    inp=SimpleNamespace(frame=cm['unitFrame']);data={'ends':{'fieldUnits':cm['fieldReferenceUnits'],'rowUnits':cm['equationReferenceUnits']}}
    threshold=packet(own['coefficientMatrix'],det['clearedMatrix'],det['rowDenominators'],det['numerator'],det['denominator'],own['wave'],
        elimination['elimination'],analysis['normalized'],analysis['complexSquarefree'],analysis['factorization'],analysis['factorReconstructionResidual'],
        analysis['intervals'],branches,own['frequencyCoordinate'],own['radicalCoordinate'],inp,data,analysis,router)
    f.atomic_pickle(base/'new-end-threshold/right-threshold-candidates.pickle',threshold)
    same(threshold['coordinateMap'],cm)
    # Native real-axis guards already ran inside the unchanged instrumented body.
    # Full-result persistence and exact metadata/source routes add no root proof.
    aliases={}
    for label,case in own['aliases'].items():
        aliases[label]={}
        for end,prior in case.items():
            shared=tuple(prior['owner'])==OWNER
            item={'address':(label,end),'owner':prior['owner'],'fullInputRoute':prior,
                'thresholdPath':str(base/'new-end-threshold/right-threshold-candidates.pickle') if shared else prior['thresholdPath'],
                'mode':'shared-new-complex-threshold-candidates' if shared else 'accepted-baseline-threshold-return',
                'physicalThresholdOrEndDomainOrNumericalReuseAccepted':False}
            folder=base/'complex-threshold-cases'/label/end.lower();folder.mkdir(parents=True);f.save(folder/'route.json',item);aliases[label][end]=item
        f.save(base/'complex-threshold-cases'/label/'case-summary.json',aliases[label])
    mutations={'actualInput':pairs['requested'],'changedInput':dict(pairs['requested'],elimination=pairs['requested']['elimination']+1),
        'changedUnit':dict(pairs['requested'],coefficientUnit=(1,0,0)),'actualOwner':OWNER,'changedOwner':(OWNER[0],'LEFT'),
        'actualAnalysis':analysis,'changedRealPart':dict(analysis,realPart=analysis['realPart']+1),
        'changedBezoutResidual':dict(analysis,bezoutResidual=analysis['bezoutResidual']+1)}
    f.atomic_pickle(base/'complex-threshold-mutation-operands.pickle',mutations)
    controls={'eliminationCoefficient':not equal(mutations['actualInput'],mutations['changedInput']),
        'coefficientUnit':not equal(mutations['actualInput'],mutations['changedUnit']),'owner':OWNER!=mutations['changedOwner'],
        'realPart':not equal(analysis,mutations['changedRealPart']),'bezoutResidual':not equal(analysis,mutations['changedBezoutResidual'])}
    f.save(base/'complex-threshold-mutation-controls.json',controls);f.require(all(controls.values()),'actual saved complex threshold input/owner/proof changes')
    counts={}
    for receipt in router.receipts:
        item=counts.setdefault(receipt['function'],{'new':0,'reused':0});item['reused' if receipt['reused'] else 'new']+=1
    f.save(base/'complex-operation-inventory.json',{'calls':router.receipts,'counts':counts,'methodSources':router.native_sources})
    f.atomic_pickle(base/'remaining-case-end-complex-threshold.pickle',{'owner':OWNER,'sourceInput':own,'eliminationInput':elimination,
        'analysis':analysis,'threshold':threshold,'operations':router.receipts,'aliases':aliases,
        'physicalThresholdsComplexEndDomainsAndContinuationUnaccepted':True})
    return {'physicalEnds':8,'sharedNewUses':2,'baselineAliases':6,'newEndOwners':1,'cases':aliases,'operationCounts':counts,
        'actualOperationCalls':len(router.receipts),'actualMutationControls':len(controls),'realCommonRootIntervals':len(analysis['intervals']),
        'newComplexAnalysisPipelines':1,'newBranchCalls':0,'physicalThresholdOrEndDomainContinuationOrNumericalReuseAccepted':False},router


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    previous.previous.previous.previous.previous.previous.previous.protect(base)
    cp,manifest=load(base);functions,joins,dispatch=joined_sources(base,manifest);prohibit(dispatch)
    result,router=construct(base,cp,functions,joins,dispatch)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen source posthash')
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'reference postidentity')
    for item in list(joins['sympyFunctions'].values())+list(router.native_sources.values()):f.require(f.digest(Path(item['path']))==item['sha256'],'actual native operation source posthash')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_COMPLEX_THRESHOLD','nativeJoins':joins,'artifacts':artifacts,
        'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
