#!/usr/bin/env python3
"""Resume after a metadata-directory fault, preserving the completed initial Poly."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import sys
import textwrap

sys.dont_write_bytecode=True
import S11c_d_remaining_case_frequency_end_complex_threshold as original

f=original.f
PREVIOUS=f.STORE/'s11c-remaining-case-frequency-20260921/end-complex-threshold'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_complex_threshold_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_end_complex_threshold_directory_repair.json'
STATE={}


def metadata_writer(base):
    """Create only fresh local JSON parents; preserve the original serializer."""
    saved=f.save
    def write(path,value):
        path=Path(path);path.relative_to(base);path.resolve().relative_to(base)
        f.require(not path.exists() and not path.is_symlink(),'fresh metadata file, never an existing reference')
        for parent in path.parents:
            if parent==base:break
            f.require(not parent.is_symlink(),'no metadata writes through reference parents')
        path.parent.mkdir(parents=True,exist_ok=True)
        return saved(path,value)
    f.save=write


def load(base):
    metadata_writer(base)
    repair=json.loads(REPAIR.read_text());old=PREVIOUS/'complete';inventory_path=Path(repair['completedInventory']['path'])
    f.require(f.digest(Path(original.__file__))==repair['originalHelperSha256'],'immutable failed complex helper')
    f.require(f.digest(old/'inputs.json')==repair['originalInputsSha256'] and f.digest(inventory_path)==repair['completedInventory']['sha256'],'completed original manifest/inventory hashes')
    inventory=json.loads(inventory_path.read_text());prior=json.loads((old/'inputs.json').read_text());manifest=copy.deepcopy(prior)
    manifest['runDirectory']=str(base);manifest['referencedInputs']={}
    ref=lambda p,n,s:original.h.source.reference(base,manifest,p,n,s)
    for name,item in inventory.items():
        path=old/name
        f.require(f.digest(path)==item['sha256'] and path.stat().st_size==item['bytes'] and str(path.resolve())==item['resolved']
            and (str(path.readlink()) if path.is_symlink() else None)==item['rawLink'],('immutable completed complex file',name))
        ref(path,'original-complex-inputs.json' if name=='inputs.json' else name,item['sha256'])
    for name,item in prior['referencedInputs'].items():
        saved=inventory[name]
        f.require(saved['rawLink']==item['original'] and saved['resolved']==item['resolvedOriginal'] and saved['sha256']==item['sha256'] and saved['bytes']==item['bytes'],'full original reference route')
    for name,item in repair['previousLogs'].items():ref(PREVIOUS/name,'original-complex-logs/'+name,item['sha256'])
    for key,name in (('completedInventory','completed-complex-file-inventory.json'),('outcomeInspection','failed-complex-threshold-outcome.json')):
        item=repair[key];ref(Path(item['path']),name,item['sha256'])
    for name,value in prior['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(old/'source'/name)==f.digest(base/'source'/name)==value,'original current/frozen source identity')
    for path in (Path(__file__).resolve(),PLAN,REPAIR):
        name=str(path.relative_to(f.ROOT));f.require(name not in manifest['sourceFiles'],'fresh recovery source pin')
        value=f.digest(path);manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
        f.require(not target.exists() and not target.is_symlink(),'fresh recovery frozen source');target.write_bytes(path.read_bytes())
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'fresh original input prehash')
    cp=json.loads((base/'accepted-elimination-checkpoint.json').read_text())
    f.require(f.digest(base/'accepted-elimination-checkpoint.json')==original.CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_ELIMINATION','same accepted resultant input')
    tree=ast.parse(Path(original.__file__).read_text())
    for name,sha in repair['wholeOriginalBodies'].items():
        node=next(n for n in tree.body if getattr(n,'name',None)==name)
        f.require(hashlib.sha256(ast.dump(node).encode()).hexdigest()==sha,('whole immutable original helper body',name))
    STATE.update(repair=repair,inventory=inventory,manifest=manifest)
    manifest['complexDirectoryRecovery']={'previousDirectory':str(PREVIOUS),'originalHelperSha256':repair['originalHelperSha256'],
        'originalInputsSha256':repair['originalInputsSha256'],'completedFilesReferenced':len(inventory),
        'previousOutcome':repair['previousOutcome'],'oldLoaderCalled':False,'completedInitialPolyCalls':1,
        'originalAnalysisCalls':0,'originalFinalChecksOrAcceptance':False,'metadataWriterOnly':'create fresh checked JSON parents; original bytes/serializer unchanged'}
    f.save(base/'completed-complex-reference-reuse.json',manifest['complexDirectoryRecovery']);f.save(base/'inputs.json',manifest)
    return cp,manifest


def packet(base,name):
    item=STATE['inventory'][name];path=base/name
    f.require(path.is_symlink() and path.readlink()==PREVIOUS/'complete'/name and str(path.resolve())==item['resolved']
        and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],('actual completed packet address/hash',name))
    return f.unpickle(path)


def restore_prefix(base,functions,joins,dispatch):
    saved=packet(base,'native-complex-threshold-input.pickle');own=saved['wholeOwnInput'];elimination=saved['elimination']
    pairs=saved['fullSavedPairs'];branch_pairs=saved['acceptedBranchCalls'];previous_result=packet(base,'remaining-case-end-elimination.pickle')
    op=packet(base,'complex-operations/0/input.pickle')
    poly=packet(base,'complex-operations/0/value.pickle');local=packet(base,'complex-locals/initial-observation0/0/value.pickle')
    original.same(poly,local)
    receipt=json.loads((base/'complex-operations/0/completed.json').read_text())
    f.require(receipt['inputSha256']==f.digest(base/'complex-operations/0/input.pickle') and receipt['valueSha256']==f.digest(base/'complex-operations/0/value.pickle')
        and receipt['site']=='initial-0' and receipt['function']=='sp.Poly' and not receipt['reused'],'completed exact first Poly receipt')
    original.same(op['input'],{'function':'sp.Poly','receiver':None,'args':(elimination['elimination'],own['frequencyCoordinate']),'kwargs':{},
        'coefficientContext':{'coefficientUnit':(0,0,0),'coordinateUnit':(0,0,0),'coordinate':own['frequencyCoordinate']}})
    original.same(saved['owner'],original.OWNER);original.same(op['owner'],original.OWNER)
    f.require(op['wholePipelineInputPath']==str(PREVIOUS/'complete/native-complex-threshold-input.pickle')
        and op['wholePipelineInputSha256']==f.digest(base/'native-complex-threshold-input.pickle'),'full actual initial caller input')
    original.same(pairs['requested']['analysisSource'],joins['acceptedRealAxisAnalysis']);original.same(pairs['requested']['callerSource'],joins['acceptedThresholdAdapter'])
    original.same(own,packet(base,'remaining-case-end-threshold-inputs.pickle'))
    original.same(elimination,packet(base,'new-end-elimination.pickle'))
    atlas=packet(base,'completed-complex-operation-atlas.pickle')
    router=object.__new__(original.AnalysisRouter)
    router.base,router.own,router.dispatch,router.joins=base,own,dispatch,joins
    router.context=op['input']['coefficientContext'];router.atlas=list(atlas['available']);router.missing=list(atlas['unavailableIntermediateReturns'])
    router.observation_counts={'initial-observation0':1,'initial-guard1':1};router.native_sources={}
    origin={'owner':op['owner'],'path':str(PREVIOUS/'complete/complex-operations/0/value.pickle'),'call':0}
    f.require(json.loads(json.dumps(origin))==receipt['owner'],'actual native JSON owner boundary only')
    native_receipt=dict(receipt,owner=origin);router.receipts=[native_receipt]
    f.require(json.loads(json.dumps(native_receipt))==receipt,'exact native receipt serializer, physical operands unchanged')
    f.atomic_pickle(base/'completed-initial-poly-owner-view.pickle',{'rawTypedOwnerFromActualInput':origin,
        'originalJsonReceipt':receipt,'nativeTypedReceipt':native_receipt,'nativeJsonDecodedReceipt':json.loads(json.dumps(native_receipt))})
    atlas_origin=dict(origin,receiptPath=str(PREVIOUS/'complete/complex-operations/0/completed.json'),receiptSha256=f.digest(base/'complex-operations/0/completed.json'))
    router.atlas.append({'input':op['input'],'value':poly,'origin':atlas_origin})
    f.save(base/'completed-initial-poly-and-guard-reuse.json',{'originalHelperSha256':STATE['repair']['originalHelperSha256'],
        'originalOperationReceipt':receipt,'fullInputValueLocalSha256':{name:f.digest(base/name) for name in ('native-complex-threshold-input.pickle','completed-complex-operation-atlas.pickle','complex-operations/0/input.pickle','complex-operations/0/value.pickle','complex-operations/0/completed.json','complex-locals/initial-observation0/0/value.pickle')},
        'traceSha256':STATE['repair']['previousLogs']['frequency_end_complex_threshold.stderr']['sha256'],
        'guardCompletionEvidence':'Literal native nonzero f.require returned before original AnalysisRouter.guard reached failed JSON parent write. No independent old guard receipt existed.',
        'initialPolyConstructorRepeated':False,'initialNonzeroComparisonRepeated':False,'oldAtlasAndPrefixRepeated':False,
        'restoredAtlasUsesActualSavedInputsValues':True,'restoredLocals':STATE['constructJoin']['restoredLocals']})
    return own,elimination,pairs,branch_pairs,router,functions[1],functions[2],poly,previous_result


def construct_adapter():
    before=ast.parse(textwrap.dedent(inspect.getsource(original.construct))).body[0]
    stop=next(i for i,n in enumerate(before.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='analysis' for t in n.targets))
    f.require(ast.unparse(before.body[stop].value)=='analyze(poly, own[\'frequencyCoordinate\'], router)','resume only at previously unstarted analysis')
    prefix=copy.deepcopy(before.body[:stop]);restored=ast.parse('own,elimination,pairs,branch_pairs,router,analyze,packet,poly,previous_result=restore_prefix(base,functions,joins,dispatch)').body[0]
    changed=copy.deepcopy(before);changed.body[:stop]=[restored]
    reverse=copy.deepcopy(changed);reverse.body[:1]=prefix
    f.require(ast.dump(reverse)==ast.dump(before),'whole original construct reversal: completed initial prefix only')
    defined={n.id for s in prefix for n in ast.walk(s) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
    used={n.id for s in before.body[stop:] for n in ast.walk(s) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Load)}
    later={n.id for s in before.body[stop:] for n in ast.walk(s) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
    required=(defined & used)-later;assigned={n.id for n in ast.walk(restored.targets[0]) if isinstance(n,ast.Name)}
    f.require(required==assigned,'all actual unfinished suffix locals restored')
    module=ast.fix_missing_locations(ast.Module(body=[changed],type_ignores=[]));env=dict(vars(original),restore_prefix=restore_prefix)
    exec(compile(module,'<saved complex initial Poly prefix continuation>','exec'),env)
    return env['construct'],{'entireOriginalConstructReverseAST':True,'originalConstructAST':hashlib.sha256(ast.dump(before).encode()).hexdigest(),
        'completedPrefixStatements':stop,'restoredLocals':sorted(required),'allLaterStatementsAndNativeGuardsUnchanged':True,'adapterSource':ast.unparse(module)}


def joined_sources(base,manifest):
    # Restore completed source observations; only compile/reverse adapters for
    # the explicit new continuation. No original joined_sources or constructor.
    old=json.loads((base/'native-end-complex-threshold-joins.json').read_text())
    functions,proof=original.native_adapter();original.same(proof,old['complexAdapter'])
    dispatch={name:getattr(original.sp,name) for name in old['sympyFunctions']}
    for name,value in dispatch.items():original.same(original.source_body(value),old['sympyFunctions'][name])
    fn,join=construct_adapter();STATE['construct']=fn;STATE['constructJoin']=join
    result=copy.deepcopy(old);result['directoryRecovery']={'originalHelperSha256':STATE['repair']['originalHelperSha256'],
        'wholeOriginalMainAST':STATE['repair']['wholeOriginalBodies']['main'],'wholeOriginalAnalysisRouterAST':STATE['repair']['wholeOriginalBodies']['AnalysisRouter'],
        'constructPrefixJoin':join,'writerSource':original.source_body(metadata_writer),'previousOutcome':STATE['repair']['previousOutcome'],
        'originalNativeJoinsSha256':f.digest(base/'native-end-complex-threshold-joins.json'),'noScientificOperationChanged':True}
    f.save(base/'recovery-native-end-complex-threshold-joins.json',result)
    return functions,result,dispatch


def construct(base,cp,functions,joins,dispatch):
    result,router=STATE['construct'](base,cp,functions,joins,dispatch)
    result.update(completedInitialPolyCallsReused=1,newInitialPolyCalls=0,operationCountsIncludeOriginalCompletedPoly=True,
        originalFailureProvenance=STATE['repair']['previousOutcome'],directoryRecovery=joins['directoryRecovery'])
    return result,router


def main():
    node=ast.parse(textwrap.dedent(inspect.getsource(original.main))).body[0]
    f.require(hashlib.sha256(ast.dump(node).encode()).hexdigest()==json.loads(REPAIR.read_text())['wholeOriginalBodies']['main'],'whole unchanged original main/final guards')
    env=dict(vars(original),load=load,joined_sources=joined_sources,construct=construct)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),'<whole original main with explicit continuation slots>','exec'),env)
    env['main']()

if __name__=='__main__':main()
