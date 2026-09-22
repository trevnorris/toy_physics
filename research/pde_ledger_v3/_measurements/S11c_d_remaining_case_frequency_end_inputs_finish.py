#!/usr/bin/env python3
"""Finish the saved first-end summary and the remaining end input catalogue."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import shutil
import sys
import textwrap

sys.dont_write_bytecode=True
import S11c_d_remaining_case_frequency_end_inputs_recover as prior

h,f,source=prior.h,prior.f,prior.source
PREVIOUS=prior.PREVIOUS.parent/'end-inputs-recovery-01'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_inputs_finish_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_end_inputs_locals_repair.json'
STATE={}
same,array_same=prior.same,prior.array_same


def packet(base,name):
    item=STATE['repair']['completedFiles'][name];path=base/name
    f.require(path.is_symlink() and path.readlink()==PREVIOUS/'complete'/name
              and str(path.resolve())==item['resolved'] and path.stat().st_size==item['bytes']
              and f.digest(path)==item['sha256'],('completed actual input packet',name))
    return f.unpickle(path)


def evidence(base,name):
    item=STATE['repair']['completedFiles'][name];path=base/name
    f.require(path.is_symlink() and path.readlink()==PREVIOUS/'complete'/name
              and str(path.resolve())==item['resolved'] and f.digest(path)==item['sha256'],'completed JSON evidence identity')
    return json.loads(path.read_text())


def route_reuse(base,receipt):
    result=[]
    for route in receipt['inputRoutes']:
        relative=Path(route['current']).relative_to(PREVIOUS/'complete');path=base/relative
        item=STATE['old']['referencedInputs'][str(relative)]
        f.require(path.is_symlink() and path.readlink()==Path(route['current'])
                  and str(path.resolve())==route['resolved']==item['resolvedOriginal']
                  and f.digest(path)==route['sha256']==item['sha256']
                  and path.stat().st_size==route['bytes']==item['bytes'],'complete saved-prefix actual input route')
        for old in (route['original'],route['acceptedSource'],route['current']):
            p=Path(old);f.require(str(p.resolve())==route['resolved'] and f.digest(p)==route['sha256'],
                                 'all actual original prefix source addresses retained')
        result.append({'current':str(path),'savedRoute':route,'fullHashAddressJoin':True})
    return result


def restore_case_prefix(base,manifest,baseline_context,label,context,state,frame,analytic,case,uniform,response):
    f.require(label==h.BASELINE,'only the completed first case prefix')
    saved=packet(base,'end-input-cases/'+label+'/context-pairs.pickle')
    previous=evidence(base,'completed-end-case-prefix-reuse.json')
    requested={'context':context,'baselineContext':baseline_context,'seedInput':state,'frame':frame,
        'analyticVariables':analytic['variables'],
        'ownSourcePath':str(prior.PREVIOUS/'complete/end-accepted/sources/cases'/label/'case-end-sources.pickle'),
        'analyticPacketPath':str(prior.PREVIOUS/'complete/analytic-cases'/label/'frequency-analytic.pickle')}
    same(saved,requested)
    f.require(previous['fullTypedRequestedSavedIdentity'] and not previous['originalGuardsRepeated']
              and previous['case']==label,'actual completed case-prefix receipt')
    routes=route_reuse(base,previous)
    f.save(base/'finish-completed-case-prefix-reuse.json',{'previousReceiptSha256':f.digest(base/'completed-end-case-prefix-reuse.json'),
        'restoredInputRoutes':routes,'fullTypedRequestedSavedIdentity':True,
        'originalCaseGuardsRepeated':False,'previousPrefixChecksRepeated':False,'completedCaseClaimed':False})


def restore_end_prefix(base,label,end,groups):
    f.require(label==h.BASELINE and end=='LEFT' and not groups,'only completed first end proof/control prefix')
    folder='end-input-cases/'+label+'/left/'
    raw=packet(base,folder+'full-inputs.pickle');selection=packet(base,folder+'whole-cluster-inputs.pickle')
    pair=packet(base,folder+'baseline-source-input-pair.pickle')
    controls=evidence(base,folder+'mutation-controls.json')
    previous=evidence(base,'completed-end-selection-prefix-reuse.json')
    routes=route_reuse(base,previous)
    proof_write=evidence(base,'baseline-proof-pair-write-reuse.json')
    f.require(proof_write['fullTypedRequestedSavedIdentity'] and proof_write['rawBytesPreserved']
              and proof_write['sha256']==f.digest(base/folder/'baseline-accepted-proof-pair.pickle'), 'actual prior proof-pair write/input join')
    f.require(set(controls)=={'coefficient','unit','frequency','address','selectedVector'}
              and all(v is True for v in controls.values()),'completed original mutation receipts')
    f.require(pair['sameFullSourceInput'] is True and previous['restoredSourceFamily']==0
              and previous['savedCandidate'] is True,'actual saved first family result')
    # Restore all eight locals consumed by the unchanged summary/hash/signature
    # suffix. In particular, selected is the captured complete list, not rebuilt
    # from vectors or a numerical selection algorithm.
    signature=raw['signature'];modes={item['info']['INDEX']:item for item in selection['allCandidateModes']}
    by_index=selection['groups'];selected=selection['selectedItems'];candidate=pair['sameFullSourceInput']
    groups.append(((label,end),signature))
    f.save(base/'finish-completed-left-end-reuse.json',{'case':label,'end':end,
        'capturedInputs':{name:f.digest(base/folder/name) for name in ('full-inputs.pickle','whole-cluster-inputs.pickle',
            'baseline-source-input-pair.pickle','baseline-accepted-proof-pair.pickle','mutation-operands.pickle','mutation-controls.json')},
        'inputRoutes':routes,'selectedField':('whole-cluster-inputs.pickle','selectedItems'),
        'restoredLocals':STATE['adapterJoin']['requiredRestoredLocals'],
        'previousHelperSha256':STATE['repair']['previousHelperSha256'],
        'traceSha256':STATE['repair']['originalLogs']['frequency_end_inputs_recover.stderr']['sha256'],
        'completionEvidence':'Exact original guard/control sequence precedes failed summary len(selected); no old end/case summary existed.',
        'originalProofOrMutationGuardsRepeated':False,'newEndSummaryPending':True})
    return base/folder,signature,modes,by_index,0,candidate,selected,controls


def adapter():
    original=ast.parse(textwrap.dedent(inspect.getsource(prior.prepare_adapter)))
    boundary=next(x for x in ast.walk(original) if isinstance(x,ast.Assign) and isinstance(x.targets[0],ast.Name) and x.targets[0].id=='boundary')
    replacement=ast.parse("boundary=next(i for i,x in enumerate(inner.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='summary[end]')").body[0]
    changed=copy.deepcopy(original);cb=next(x for x in ast.walk(changed) if isinstance(x,ast.Assign) and isinstance(x.targets[0],ast.Name) and x.targets[0].id=='boundary')
    cb.value=copy.deepcopy(replacement.value)
    old_literal='target,signature,background,channel,modes,by_index,group,candidate,selection=restore_end_prefix(base,label,end,groups)'
    new_literal='target,signature,modes,by_index,group,candidate,selected,controls=restore_end_prefix(base,label,end,groups)'
    strings=[x for x in ast.walk(changed) if isinstance(x,ast.Constant) and x.value==old_literal]
    f.require(len(strings)==1,'one exact restored-local assignment boundary');strings[0].value=new_literal
    reverse=copy.deepcopy(changed)
    rb=next(x for x in ast.walk(reverse) if isinstance(x,ast.Assign) and isinstance(x.targets[0],ast.Name) and x.targets[0].id=='boundary')
    rb.value=copy.deepcopy(boundary.value)
    next(x for x in ast.walk(reverse) if isinstance(x,ast.Constant) and x.value==new_literal).value=old_literal
    f.require(ast.dump(reverse)==ast.dump(original),'whole prior adapter reverse identity: prefix extent and consumed local list only')
    ns=dict(vars(prior),restore_case_prefix=restore_case_prefix,restore_end_prefix=restore_end_prefix)
    exec(compile(ast.fix_missing_locations(changed),'<completed-end-summary-adapter>','exec'),ns)
    prepare,join=ns['prepare_adapter']()
    whole=ast.parse(textwrap.dedent(inspect.getsource(h.prepare)))
    outer=next(x for x in whole.body[0].body if isinstance(x,ast.For) and ast.unparse(x.iter)=='labels')
    inner=next(x for x in outer.body if isinstance(x,ast.For) and ast.unparse(x.target)=='end')
    index=next(i for i,x in enumerate(inner.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='summary[end]')
    local_names={x.id for x in ast.walk(whole) if isinstance(x,ast.Name) and isinstance(x.ctx,ast.Store)}
    used={x.id for stmt in inner.body[index:] for x in ast.walk(stmt) if isinstance(x,ast.Name) and isinstance(x.ctx,ast.Load)}
    required=(used & local_names)-{'summary','end','groups','signatures','label'}
    assigned={x.id for x in ast.walk(ast.parse(new_literal).body[0].targets[0]) if isinstance(x,ast.Name)}
    f.require(required==assigned,'every original summary/signature suffix local is explicitly restored')
    result={'wholePriorAdapterReverseAST':True,'priorAdapterAST':prior.body(prior.prepare_adapter),
        'wholeOriginalPrepareJoin':join,'requiredRestoredLocals':sorted(required),'assignedRestoredLocals':sorted(assigned),
        'originalSummarySuffixAST':hashlib.sha256(ast.dump(ast.Module(body=inner.body[index:],type_ignores=[])).encode()).hexdigest(),
        'originalSummaryAndRemainingEndGuardsUnchanged':True}
    return prepare,result


def protect(base):
    prior.original_protect(base);save_json,save_pickle=f.save,f.atomic_pickle
    def checked(path):
        path=Path(path);path.relative_to(base)
        f.require(not path.is_symlink(),'no reference writes')
        f.require(not path.exists() or path==base/'end-input-case-inventory.json','fresh output or native growing case inventory')
        return path
    def json_writer(path,value):
        path=Path(path)
        if path==base/'native-end-input-joins.json':
            f.require(value==json.loads(path.read_text())==json.loads((base/'recovery-native-end-input-joins.json').read_text()),
                      'original and previous native whole-source joins unchanged')
            return save_json(checked(base/'finish-native-end-input-joins.json'),value)
        return save_json(checked(path),value)
    def pickle_writer(path,value):return save_pickle(checked(path),value)
    f.save,f.atomic_pickle=json_writer,pickle_writer


def load(base):
    repair=json.loads(REPAIR.read_text());oldbase=PREVIOUS/'complete';old=json.loads((oldbase/'inputs.json').read_text())
    STATE.update(repair=repair,old=old,base=base)
    for path,key in ((Path(__file__),'newHelperSha256'),(PLAN,'newPlanSha256'),(Path(prior.__file__),'previousHelperSha256'),
                     (Path(h.__file__),'originalHelperSha256'),(prior.REPAIR,'previousRepairSha256')):
        f.require(f.digest(path)==repair[key],('reviewed immutable continuation source',str(path)))
    f.require(prior.body(h.main)==repair['adapterJoin']['wholeOriginalPrepareJoin']['originalMainAST'],'whole original main unchanged')
    f.require(f.digest(oldbase/'inputs.json')==repair['originalInputsSha256'],'completed prior loader unchanged')
    actual={str(p.relative_to(oldbase)) for p in oldbase.rglob('*') if p.is_file()
            and 'source' not in p.relative_to(oldbase).parts and p!=oldbase/'inputs.json'}
    f.require(actual==set(repair['completedFiles']),'all previous completed files preserved')
    manifest=dict(old,runDirectory=str(base),sourceFiles=dict(old['sourceFiles']),inputPackets=dict(old['inputPackets']),referencedInputs={})
    for name,item in repair['completedFiles'].items():
        path=oldbase/name
        f.require(f.digest(path)==item['sha256'] and path.stat().st_size==item['bytes'] and str(path.resolve())==item['resolved']
                  and (str(path.readlink()) if path.is_symlink() else None)==item['rawLink'],('previous file bytes and address',name))
        source.reference(base,manifest,path,name,item['sha256'])
    source.reference(base,manifest,oldbase/'inputs.json','previous-end-recovery-inputs.json',repair['originalInputsSha256'])
    for name,item in repair['originalLogs'].items():source.reference(base,manifest,PREVIOUS/name,'previous-end-recovery-logs/'+name,item['sha256'])
    for path in (Path(__file__).resolve(),PLAN,REPAIR):
        name=str(path.relative_to(f.ROOT));f.require(name not in manifest['sourceFiles'],'separate recovery source pin');manifest['sourceFiles'][name]=f.digest(path)
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==value,'current source prehash')
        if name in old['sourceFiles']:f.require(f.digest(oldbase/'source'/name)==value,'previous frozen source prehash')
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,target)
        f.require(f.digest(target)==value,'new frozen source hash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'all original/reference input prehashes')
    trace=(base/'previous-end-recovery-logs/frequency_end_inputs_recover.stderr').read_text()
    f.require("UnboundLocalError: local variable 'selected' referenced before assignment" in trace
              and 'line 95, in prepare' in trace,'actual failed summary after completed left-end guards')
    guard=json.loads((base/'previous-end-recovery-logs/resource-guard/outcome.json').read_text())
    inv=json.loads((base/'previous-end-recovery-logs/frequency_end_inputs_recover.invocation.json').read_text())
    f.require(guard['exitCode']==guard['childOutcome']['exitCode']==inv['exitCode']==1
              and guard['childOutcome']['guardReason'] is None,'actual prior guarded runtime failure')
    array_receipt=evidence(base,'array-reader-mutation-controls.json')
    f.require(set(array_receipt)=={'value','dtype','shape','container','originalReaderArrayBoundary','exactCopyAccepted'}
              and all(v is True for v in array_receipt.values()),'completed exact array reader controls')
    f.save(base/'finish-completed-array-controls-reuse.json',{'controls':array_receipt,
        'operandSha256':f.digest(base/'array-reader-mutation-operands.pickle'),
        'controlsSha256':f.digest(base/'array-reader-mutation-controls.json'),
        'actualPreviousSourceTraceJoin':True,'controlsOrOriginalFailureInstrumentRepeated':False})
    manifest['summaryLocalsRecovery']={'previousRun':str(PREVIOUS),'previousFailure':repair['failure'],
        'previousInputsSha256':repair['originalInputsSha256'],'completedFilesReused':len(repair['completedFiles']),
        'adapterJoin':STATE['adapterJoin'],'completedLeftProofAndControlsReused':True,
        'originalOrNewScientificConstructionRepeated':False}
    f.save(base/'inputs.json',manifest)
    f.save(base/'summary-locals-reference-reuse.json',manifest['referencedInputs'])
    f.save(base/'summary-locals-continuation-join.json',manifest['summaryLocalsRecovery'])
    return manifest,tuple(repair['cases'])


if __name__=='__main__':
    prepare,join=adapter();STATE['adapterJoin']=join
    f.require(join==json.loads(REPAIR.read_text())['adapterJoin'],'reviewed whole source/body and complete suffix local joins')
    h.inputs.protect_references=protect;h.load=load;h.prepare=prepare;h.same=same
    # This is the original producer prohibition; the six completed array controls
    # are reused above and must not be called through prior.prohibit again.
    h.prohibit=prior.original_prohibit
    h.main()
