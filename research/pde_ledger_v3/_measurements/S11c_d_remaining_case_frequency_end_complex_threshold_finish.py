#!/usr/bin/env python3
"""Finish saved complex threshold output, aliases and controls without science replay."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import sys
import textwrap
from types import SimpleNamespace

sys.dont_write_bytecode=True
import S11c_d_remaining_case_frequency_end_complex_threshold_recover as previous
original=previous.original
f=original.f
PREVIOUS=f.STORE/'s11c-remaining-case-frequency-20260921/end-complex-threshold-recovery-01'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_complex_threshold_finish_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_end_complex_threshold_packet_repair.json'
STATE={}


def fresh_writers(base):
    saved_json,saved_pickle=f.save,f.atomic_pickle
    def checked(path):
        path=Path(path);path.relative_to(base);path.resolve().relative_to(base)
        f.require(not path.exists() and not path.is_symlink(),'fresh output, never existing data or reference')
        for parent in path.parents:
            if parent==base:break
            f.require(not parent.is_symlink(),'no reference-parent writes')
        path.parent.mkdir(parents=True,exist_ok=True)
        return path
    f.save=lambda path,value:saved_json(checked(path),value)
    f.atomic_pickle=lambda path,value:saved_pickle(checked(path),value)


def load(base):
    fresh_writers(base)
    repair=json.loads(REPAIR.read_text());old=PREVIOUS/'complete';inventory_path=Path(repair['completedInventory']['path'])
    f.require(f.digest(Path(original.__file__))==repair['originalHelperSha256'],'immutable original complex helper')
    f.require(f.digest(Path(previous.__file__))==repair['previousHelperSha256'],'immutable first continuation helper')
    f.require(f.digest(old/'inputs.json')==repair['originalInputsSha256'] and f.digest(inventory_path)==repair['completedInventory']['sha256'],'completed original manifest/inventory hashes')
    inventory=json.loads(inventory_path.read_text());prior=json.loads((old/'inputs.json').read_text());manifest=copy.deepcopy(prior)
    manifest['runDirectory']=str(base);manifest['referencedInputs']={}
    ref=lambda p,n,s:original.h.source.reference(base,manifest,p,n,s)
    for name,item in inventory.items():
        path=old/name
        f.require(f.digest(path)==item['sha256'] and path.stat().st_size==item['bytes'] and str(path.resolve())==item['resolved']
            and (str(path.readlink()) if path.is_symlink() else None)==item['rawLink'],('immutable completed complex file',name))
        ref(path,'previous-complex-recovery-inputs.json' if name=='inputs.json' else name,item['sha256'])
    # Fulfil the failed destination with the exact already saved final packet
    # bytes. No reconstruction or serialization of that scientific value.
    ref(old/'complex-locals/packet-observation1/0/value.pickle','new-end-threshold/right-threshold-candidates.pickle',repair['previousOutcome']['actualThresholdPacketSha256'])
    for name,item in prior['referencedInputs'].items():
        saved=inventory[name]
        f.require(saved['rawLink']==item['original'] and saved['resolved']==item['resolvedOriginal'] and saved['sha256']==item['sha256'] and saved['bytes']==item['bytes'],'full original reference route')
    for name,item in repair['previousLogs'].items():ref(PREVIOUS/name,'previous-complex-recovery-logs/'+name,item['sha256'])
    for key,name in (('completedInventory','completed-complex-recovery-file-inventory.json'),('outcomeInspection','failed-complex-packet-outcome.json')):
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
    previous_tree=ast.parse(Path(previous.__file__).read_text())
    for name,sha in repair['wholePreviousBodies'].items():
        node=next(n for n in previous_tree.body if getattr(n,'name',None)==name)
        f.require(hashlib.sha256(ast.dump(node).encode()).hexdigest()==sha,('whole first continuation helper body',name))
    STATE.update(repair=repair,inventory=inventory,manifest=manifest)
    manifest['complexPacketFinish']={'previousDirectory':str(PREVIOUS),'originalHelperSha256':repair['originalHelperSha256'],
        'previousHelperSha256':repair['previousHelperSha256'],'originalInputsSha256':repair['originalInputsSha256'],
        'completedFilesReferenced':len(inventory),'previousOutcome':repair['previousOutcome'],'oldLoaderCalled':False,
        'completedNativeOperationReceipts':58,'newAnalysisOrNativeProofCalls':0,'originalFinalChecksOrAcceptance':False,
        'thresholdOutputUsesExactSavedBytes':True,'originalScientificValuesReconstructed':False}
    f.save(base/'completed-complex-packet-reference-reuse.json',manifest['complexPacketFinish']);f.save(base/'inputs.json',manifest)
    return cp,manifest


def packet(base,name):
    item=STATE['inventory'][name];path=base/name
    f.require(path.is_symlink() and path.readlink()==PREVIOUS/'complete'/name and str(path.resolve())==item['resolved']
        and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],('actual completed packet address/hash',name))
    return f.unpickle(path)



def restore_completed(base,functions,joins,dispatch):
    saved=packet(base,'native-complex-threshold-input.pickle');own=saved['wholeOwnInput'];elimination=saved['elimination'];pairs=saved['fullSavedPairs']
    analysis=packet(base,'new-complex-analysis.pickle');threshold=packet(base,'complex-locals/packet-observation1/0/value.pickle');cm=own['coordinateMap']
    destination=base/'new-end-threshold/right-threshold-candidates.pickle';target=PREVIOUS/'complete/complex-locals/packet-observation1/0/value.pickle'
    f.require(destination.is_symlink() and destination.readlink()==target and f.digest(destination)==STATE['repair']['previousOutcome']['actualThresholdPacketSha256']
        and destination.read_bytes()==target.read_bytes(),'exact saved full native threshold output bytes')
    original.same(threshold['realAxisAnalysis'],analysis)
    original.same(analysis,packet(base,'complex-locals/analysis-observation21/0/value.pickle'))
    branch=packet(base,'complex-accepted-branch-return.pickle');original.same(threshold['bulkBranchFrequencies'],branch['value'])
    original.same(own,packet(base,'remaining-case-end-threshold-inputs.pickle'));original.same(elimination,packet(base,'new-end-elimination.pickle'))
    original.same(pairs['requested']['analysisSource'],joins['acceptedRealAxisAnalysis']);original.same(pairs['requested']['callerSource'],joins['acceptedThresholdAdapter'])
    receipt_names=sorted([n for n in STATE['inventory'] if n.startswith('complex-operations/') and n.endswith('/completed.json')],key=lambda n:int(Path(n).parent.name))
    f.require(len(receipt_names)==58,'all completed native operation receipts')
    receipts=[];owner_views=[]
    for name in receipt_names:
        folder=Path(name).parent;input_name=str(folder/'input.pickle');value_name=str(folder/'value.pickle')
        actual=packet(base,input_name);value=packet(base,value_name);receipt=json.loads((base/name).read_text())
        f.require(receipt['inputSha256']==f.digest(base/input_name) and receipt['valueSha256']==f.digest(base/value_name),'actual completed before/after hashes')
        if receipt['reused']:
            f.require(actual['matches'],'captured live matching call exists');owner=actual['matches'][0]['origin']
        else:
            owner={'owner':actual['owner'],'path':str(PREVIOUS/'complete'/folder/'value.pickle'),'call':int(folder.name)}
            if folder.name=='0':
                # Call zero predates this recovery and its raw typed owner view
                # is already saved under the original-source/receipt join.
                first=packet(base,'completed-initial-poly-owner-view.pickle');owner=first['nativeTypedReceipt']['owner']
        native=dict(receipt,owner=owner)
        owner_views.append({'call':int(folder.name),'rawTypedOwner':owner,'originalJsonReceipt':receipt,
            'nativeJsonReceipt':json.loads(json.dumps(native)),'inputPath':str(base/input_name),'valuePath':str(base/value_name),
            'actualValueType':type(value).__module__+'.'+type(value).__qualname__})
        receipts.append(native)
    f.atomic_pickle(base/'completed-complex-operation-owner-views.pickle',owner_views)
    for view in owner_views:f.require(view['nativeJsonReceipt']==view['originalJsonReceipt'],'native receipt JSON owner boundary, no physical conversion')
    methods={}
    for name in sorted(n for n in STATE['inventory'] if n.startswith('complex-method-source/') and n.endswith('.json')):
        item=json.loads((base/name).read_text());key=item.pop('method')
        f.require(key not in methods and f.digest(Path(item['path']))==item['sha256'],'completed native method-source identity')
        methods[key]=item
    router=SimpleNamespace(receipts=receipts,native_sources=methods)
    guards={name:json.loads((base/name).read_text()) for name in STATE['inventory'] if name.startswith('complex-guards/') and name.endswith('.json')}
    f.require(guards=={'complex-guards/'+k:v for k,v in STATE['repair']['previousOutcome']['completedNativeGuards'].items()},'all saved native guard receipts')
    # Read representation views of completed values only. No Poly/as_expr,
    # normalization, root operation, proof arithmetic or physical conversion.
    overview={'normalizedElimination':str(threshold['normalizedElimination']),
        'realCommonFactor':str(analysis['realCommonFactor']),'realSquarefree':str(analysis['realSquarefree']),
        'realPart':str(analysis['realPart']),'imaginaryPart':str(analysis['imaginaryPart']),
        'decompositionResidual':str(analysis['decompositionResidual']),'divisionResiduals':[str(v) for v in analysis['divisionResiduals']],
        'bezoutResidual':str(analysis['bezoutResidual']),'factorReconstructionResidual':str(analysis['factorReconstructionResidual']),
        'intervals':[[str(a),str(b),str(n)] for (a,b),n in analysis['intervals']],
        'scope':analysis['scope'],'representationViewOnly':True,'nativeAnalysisAndGuardProofsReusedNotRepeated':True}
    f.save(base/'completed-complex-analysis-view.json',overview)
    f.save(base/'completed-complex-analysis-and-packet-reuse.json',{'originalHelperSha256':STATE['repair']['originalHelperSha256'],
        'previousHelperSha256':STATE['repair']['previousHelperSha256'],'analysisSha256':f.digest(base/'new-complex-analysis.pickle'),
        'thresholdSavedViewSha256':f.digest(destination),'thresholdOutputReference':str(target),'thresholdBytesPreserved':True,
        'allCompletedNativeOperationReceipts':58,'completedNativeGuards':guards,'originalGuardReceiptsReusedWithoutRecheckingProofs':True,
        'originalAnalysisPacketAndBranchCallsRepeated':False,'newScientificConstructorCalls':0,
        'traceSha256':STATE['repair']['previousLogs']['frequency_end_complex_threshold_recover.stderr']['sha256'],
        'failedOutputWriteCompletedBySavedReference':True,'originalOutputWriteWasNotCompleted':True,
        'restoredLocals':STATE['constructJoin']['restoredLocals']})
    return threshold,cm,own,analysis,router,pairs,elimination


def construct_adapter():
    before=ast.parse(textwrap.dedent(inspect.getsource(original.construct))).body[0]
    boundary=next(i for i,n in enumerate(before.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
        and ast.unparse(n.value.func)=='f.atomic_pickle' and ast.unparse(n.value.args[0])=="base / 'new-end-threshold/right-threshold-candidates.pickle'")
    stop=boundary+1;prefix=copy.deepcopy(before.body[:stop])
    f.require(ast.unparse(before.body[boundary].value.args[1])=='threshold','exact failed output value')
    replacement=ast.parse('threshold,cm,own,analysis,router,pairs,elimination=restore_completed(base,functions,joins,dispatch)').body[0]
    changed=copy.deepcopy(before);changed.body[:stop]=[replacement]
    reverse=copy.deepcopy(changed);reverse.body[:1]=prefix
    f.require(ast.dump(reverse)==ast.dump(before),'whole original construct: saved complete prefix plus failed output reference only')
    defined={n.id for s in prefix for n in ast.walk(s) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
    used={n.id for s in before.body[stop:] for n in ast.walk(s) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Load)}
    later={n.id for s in before.body[stop:] for n in ast.walk(s) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
    required=(defined & used)-later;assigned={n.id for n in ast.walk(replacement.targets[0]) if isinstance(n,ast.Name)}
    f.require(required==assigned,'all unfinished alias/control/hash suffix locals restored')
    module=ast.fix_missing_locations(ast.Module(body=[changed],type_ignores=[]));env=dict(vars(original),restore_completed=restore_completed)
    exec(compile(module,'<completed complex threshold packet output continuation>','exec'),env)
    return env['construct'],{'entireOriginalConstructReverseAST':True,'originalConstructAST':hashlib.sha256(ast.dump(before).encode()).hexdigest(),
        'completedScientificPrefixAndFailedOutputStatements':stop,'restoredLocals':sorted(required),
        'remainingNativeCoordinateMapAliasControlAndFinalStatementsUnchanged':True,'adapterSource':ast.unparse(module)}


def joined_sources(base,manifest):
    old=json.loads((base/'native-end-complex-threshold-joins.json').read_text())
    functions,proof=original.native_adapter();original.same(proof,old['complexAdapter'])
    unused,prior_proof=previous.construct_adapter()
    prior_join=json.loads((base/'recovery-native-end-complex-threshold-joins.json').read_text())
    original.same(prior_proof,prior_join['directoryRecovery']['constructPrefixJoin'])
    dispatch={name:getattr(original.sp,name) for name in old['sympyFunctions']}
    for name,value in dispatch.items():original.same(original.source_body(value),old['sympyFunctions'][name])
    fn,join=construct_adapter();STATE['construct']=fn;STATE['constructJoin']=join
    result=copy.deepcopy(old);result['packetFinish']={'originalHelperSha256':STATE['repair']['originalHelperSha256'],
        'previousHelperSha256':STATE['repair']['previousHelperSha256'],'wholeOriginalMainAST':STATE['repair']['wholeOriginalBodies']['main'],
        'wholePreviousMainAST':STATE['repair']['wholePreviousBodies']['main'],'constructPrefixJoin':join,
        'writerSource':original.source_body(fresh_writers),'previousOutcome':STATE['repair']['previousOutcome'],
        'originalNativeJoinsSha256':f.digest(base/'native-end-complex-threshold-joins.json'),
        'previousContinuationJoinsSha256':f.digest(base/'recovery-native-end-complex-threshold-joins.json'),
        'allEarlierNativeAndContinuationAdaptersCompileOnly':True,'noScientificOperationChangedOrRepeated':True}
    f.save(base/'finish-native-end-complex-threshold-joins.json',result)
    return functions,result,dispatch


def forbidden(*args,**kwargs):raise RuntimeError('saved threshold packet finish cannot repeat any native analysis or proof operation')


def prohibit(dispatch):
    original.prohibit({})
    for name in dispatch:
        if name!='Poly':setattr(original.sp,name,forbidden)
    for name in ('call','operation','arithmetic','observe','guard'):setattr(original.AnalysisRouter,name,forbidden)


def construct(base,cp,functions,joins,dispatch):
    result,router=STATE['construct'](base,cp,functions,joins,dispatch)
    result.update(completedInitialPolyCallsReused=1,newInitialPolyCalls=0,completedNativeOperationReceiptsReused=58,
        newNativeAnalysisOrProofCalls=0,operationCountsIncludePreviousCompletedScience=True,
        originalFailureProvenance=STATE['repair']['previousOutcome'],packetFinish=joins['packetFinish'])
    return result,router


def main():
    node=ast.parse(textwrap.dedent(inspect.getsource(original.main))).body[0]
    f.require(hashlib.sha256(ast.dump(node).encode()).hexdigest()==json.loads(REPAIR.read_text())['wholeOriginalBodies']['main'],'whole unchanged original main/final guards')
    env=dict(vars(original),load=load,joined_sources=joined_sources,construct=construct,prohibit=prohibit)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),'<whole original main with saved packet finish slots>','exec'),env)
    env['main']()

if __name__=='__main__':main()
