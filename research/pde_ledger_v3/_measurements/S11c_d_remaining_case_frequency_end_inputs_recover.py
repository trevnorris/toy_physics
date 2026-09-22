#!/usr/bin/env python3
"""Finish saved end inputs with the accepted exact array reader and saved prefix."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import shutil
import sys
import textwrap

sys.dont_write_bytecode = True
import S11c_d_remaining_case_frequency_end_inputs as h
import S11c_d_remaining_case_modes as modes_reader

f,source=h.f,h.source
PREVIOUS=Path('/var/projects/toy_physics/_scratch/s11c/s11c-remaining-case-frequency-20260921/end-inputs')
PLAN=f.M/'S11c_d_remaining_case_frequency_end_inputs_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_end_inputs_array_repair.json'
STATE={}
original_protect=h.inputs.protect_references
original_prohibit=h.prohibit
array_same=modes_reader.same


def body(value):return hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(value)))).encode()).hexdigest()


def same(a,b):f.require(array_same(a,b),'full typed saved end input identity including exact array dtype/shape/value')


def saved(base,name):
    item=STATE['repair']['completedFiles'][name];path=base/name
    f.require(path.is_symlink() and str(path.readlink())==str(PREVIOUS/'complete'/name)
              and str(path.resolve())==item['resolved'] and path.stat().st_size==item['bytes']
              and f.digest(path)==item['sha256'],('completed packet input/address identity',name))
    return f.unpickle(path)


def input_route(base,name):
    item=STATE['old']['referencedInputs'][name];old=PREVIOUS/'complete'/name;path=base/name
    f.require(old.is_symlink() and str(old.readlink())==item['original']
              and str(old.resolve())==item['resolvedOriginal'] and f.digest(old)==item['sha256'], 'original accepted reference route')
    f.require(path.is_symlink() and path.readlink()==old and path.resolve()==old.resolve()
              and f.digest(path)==item['sha256'] and path.stat().st_size==item['bytes'], 'recovery actual input joins original source')
    return {'original':str(old),'current':str(path),'acceptedSource':item['original'],
            'resolved':str(path.resolve()),'sha256':item['sha256'],'bytes':item['bytes']}


def restore_case_prefix(base,manifest,baseline_context,label,context,state,frame,analytic,case,uniform,response):
    f.require(label==h.BASELINE,'only the completed first case prefix')
    name='end-input-cases/'+label+'/context-pairs.pickle';packet=saved(base,name)
    own='end-accepted/sources/cases/'+label+'/case-end-sources.pickle'
    analytic_name='analytic-cases/'+label+'/frequency-analytic.pickle'
    requested={'context':context,'baselineContext':baseline_context,'seedInput':state,'frame':frame,
               'analyticVariables':analytic['variables'],'ownSourcePath':str(PREVIOUS/'complete'/own),
               'analyticPacketPath':str(PREVIOUS/'complete'/analytic_name)}
    same(packet,requested)
    names=['accepted/contexts/'+label+'/native-input.pickle','accepted/contexts/'+label+'/seed.pickle',
           'accepted/accepted-contexts/'+label+'/frame.pickle',analytic_name,own,
           'end-accepted/sources/remaining-case-end-sources.pickle',
           'end-accepted/uniform/cases/'+label+'/continuum/uniform-response.pickle']
    routes=[input_route(base,n) for n in names]
    # No original context/physical/unit proof is executed here. Its complete
    # input packet and exact accepted producer routes restore that saved prefix.
    f.save(base/'completed-end-case-prefix-reuse.json',{'case':label,'packet':name,'sha256':f.digest(base/name),
        'inputRoutes':routes,'originalHelperSha256':STATE['repair']['originalHelperSha256'],
        'failedTraceSha256':STATE['repair']['originalLogs']['frequency_end_inputs.stderr']['sha256'],
        'originalCasePrefixAST':STATE['prepareJoin']['casePrefixAST'],
        'fullTypedRequestedSavedIdentity':True,'originalGuardsRepeated':False,'completedCaseClaimed':False})


def restore_end_prefix(base,label,end,groups):
    f.require(label==h.BASELINE and end=='LEFT' and not groups,'only completed first physical end prefix')
    folder='end-input-cases/'+label+'/left/';raw=saved(base,folder+'full-inputs.pickle')
    selection=saved(base,folder+'whole-cluster-inputs.pickle')
    pair=saved(base,folder+'baseline-source-input-pair.pickle')
    # These are captured locals, never a new mode, current, source or binding.
    signature=raw['signature'];background=raw['uniformBackground'];channel=selection['channel']
    modes={item['info']['INDEX']:item for item in selection['allCandidateModes']}
    by_index=selection['groups'];candidate=pair['sameFullSourceInput']
    same(pair['actualAddress'],(label,end));same(raw['address'],(label,end))
    same(pair['actual'],signature);same(pair['baseline'],signature)
    same(selection['channel'],background['response']['channels'][end])
    f.require(type(candidate) is bool and candidate,'saved original first self-candidate result')
    names=['end-accepted/uniform/remaining-case-currents.pickle',
           'end-accepted/uniform/cases/'+label+'/left/modal.pickle',
           'end-accepted/uniform/cases/'+label+'/left/mode-inputs.pickle',
           'end-accepted/uniform/input-routes.json']
    routes=[input_route(base,n) for n in names]
    f.require(raw['fullCurrentPacketPath']==routes[0]['original'] and raw['fullModalPacketPath']==routes[1]['original'],
              'saved actual complete current/modal input addresses')
    groups.append(((label,end),signature))
    f.save(base/'completed-end-selection-prefix-reuse.json',{'case':label,'end':end,
        'packets':{n:f.digest(base/folder/n) for n in ('full-inputs.pickle','whole-cluster-inputs.pickle','baseline-source-input-pair.pickle')},
        'inputRoutes':routes,'restoredSourceFamily':0,'savedCandidate':candidate,
        'originalEndPrefixAST':STATE['prepareJoin']['endPrefixAST'],
        'originalSelectionAndSourceGuardsRepeated':False,'completedEndClaimed':False,
        'unfinished':'Entire baseline table/background/channel guard group and all later guards; no sub-call completion is inferred.'})
    return base/folder,signature,background,channel,modes,by_index,0,candidate,selection


def prepare_adapter():
    original=ast.parse(textwrap.dedent(inspect.getsource(h.prepare)));changed=copy.deepcopy(original)
    outer=next(x for x in changed.body[0].body if isinstance(x,ast.For) and ast.unparse(x.iter)=='labels')
    # The first directory already contains saved inputs, but is a real directory.
    mkdir=next(x for x in ast.walk(outer) if isinstance(x,ast.Call) and ast.unparse(x.func)=='directory.mkdir')
    mkdir.keywords.append(ast.keyword(arg='exist_ok',value=ast.Constant(True)))
    start=next(i for i,x in enumerate(outer.body) if isinstance(x,ast.Expr) and isinstance(x.value,ast.Call)
               and ast.unparse(x.value.func)=='f.atomic_pickle')
    stop=next(i for i,x in enumerate(outer.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='summary')
    case_prefix=copy.deepcopy(outer.body[start:stop])
    restored=ast.parse('restore_case_prefix(base,manifest,baseline_context,label,context,state,frame,analytic,case,uniform,response)').body
    outer.body[start:stop]=[ast.If(test=ast.parse('label==BASELINE',mode='eval').body,body=restored,orelse=copy.deepcopy(case_prefix))]
    inner=next(x for x in outer.body if isinstance(x,ast.For) and ast.unparse(x.target)=='end')
    boundary=next(i for i,x in enumerate(inner.body) if isinstance(x,ast.If) and ast.unparse(x.test)=='label == BASELINE')
    end_prefix=copy.deepcopy(inner.body[:boundary])
    restored=ast.parse('target,signature,background,channel,modes,by_index,group,candidate,selection=restore_end_prefix(base,label,end,groups)').body
    inner.body[:boundary]=[ast.If(test=ast.parse("label==BASELINE and end=='LEFT'",mode='eval').body,
                                    body=restored,orelse=copy.deepcopy(end_prefix))]
    class Reader(ast.NodeTransformer):
        count=0
        def visit_Call(self,node):
            self.generic_visit(node)
            if ast.unparse(node.func)=='native.same':node.func=ast.Name(id='array_same',ctx=ast.Load());self.count+=1
            return node
    reader=Reader();changed=reader.visit(changed)
    reverse=copy.deepcopy(changed)
    class Reverse(ast.NodeTransformer):
        def visit_Call(self,node):
            self.generic_visit(node)
            if ast.unparse(node.func)=='array_same':node.func=ast.Attribute(value=ast.Name(id='native',ctx=ast.Load()),attr='same',ctx=ast.Load())
            return node
    reverse=Reverse().visit(reverse)
    ro=next(x for x in reverse.body[0].body if isinstance(x,ast.For) and ast.unparse(x.iter)=='labels')
    ri=next(x for x in ro.body if isinstance(x,ast.For) and ast.unparse(x.target)=='end')
    ri.body[:1]=copy.deepcopy(end_prefix)
    ro.body[start:start+1]=copy.deepcopy(case_prefix)
    rm=next(x for x in ast.walk(ro) if isinstance(x,ast.Call) and ast.unparse(x.func)=='directory.mkdir')
    rm.keywords.pop()
    f.require(ast.dump(reverse)==ast.dump(original),'entire original prepare reverse identity')
    ns=dict(vars(h),array_same=array_same,same=same,restore_case_prefix=restore_case_prefix,restore_end_prefix=restore_end_prefix)
    exec(compile(ast.fix_missing_locations(changed),'<saved-end-prefix-array-reader>','exec'),ns)
    return ns['prepare'],{'wholePrepareReverseAST':True,'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'casePrefixAST':hashlib.sha256(ast.dump(ast.Module(body=case_prefix,type_ignores=[])).encode()).hexdigest(),
        'endPrefixAST':hashlib.sha256(ast.dump(ast.Module(body=end_prefix,type_ignores=[])).encode()).hexdigest(),
        'arrayReaderCalls':reader.count,'existingDirectoryKeyword':1,'arrayComparatorAST':body(array_same),
        'originalMainAST':body(h.main),'originalSameAST':body(h.same),'allRemainingGuardsAndMainUnchanged':True}


def protect(base):
    original_protect(base);save_json,save_pickle=f.save,f.atomic_pickle
    baseline=base/'end-input-cases'/h.BASELINE/'left/baseline-accepted-proof-pair.pickle'
    def checked(path):
        path=Path(path);path.relative_to(base)
        f.require(not path.is_symlink(),'no reference writer')
        # The original main has one growing fresh case inventory by design.
        f.require(not path.exists() or path==base/'end-input-case-inventory.json','fresh output or native growing inventory')
        return path
    def json_writer(path,value):
        path=Path(path)
        if path==base/'native-end-input-joins.json':
            f.require(value==json.loads(path.read_text()),'whole original native input join unchanged')
            save_json(checked(base/'recovery-native-end-input-joins.json'),value)
            return
        return save_json(checked(path),value)
    def pickle_writer(path,value):
        path=Path(path)
        if path==baseline:
            old=saved(base,str(path.relative_to(base)));same(old,value)
            save_json(checked(base/'baseline-proof-pair-write-reuse.json'),{'path':str(path),'sha256':f.digest(path),
                'fullTypedRequestedSavedIdentity':True,'rawBytesPreserved':True,'originalWriteRepeated':False,
                'proofGroupStatus':'All comparisons following the saved pair are unfinished here.'})
            return
        return save_pickle(checked(path),value)
    f.save,f.atomic_pickle=json_writer,pickle_writer


def array_controls(base):
    packet=saved(base,'end-input-cases/'+h.BASELINE+'/left/whole-cluster-inputs.pickle')
    array=packet['selectedItems'][0]['vector'];np=modes_reader.np
    f.require(isinstance(array,np.ndarray) and array.size>1 and array.dtype.kind in 'fc','actual saved selected mode array')
    value=array.copy();value.flat[0]+=1
    dtype=np.asarray(array,dtype=np.complex64 if array.dtype!=np.dtype('complex64') else np.complex128)
    shape=array.reshape(array.shape+(1,))
    original={'items':[('selected-vector',array)]}
    changed={'value':{'items':[('selected-vector',value)]},'dtype':{'items':[('selected-vector',dtype)]},
             'shape':{'items':[('selected-vector',shape)]},'container':{'items':(('selected-vector',array),)}}
    f.atomic_pickle(base/'array-reader-mutation-operands.pickle',{'original':original,'changed':changed,
        'sourcePacket':str(base/'end-input-cases'/h.BASELINE/'left/whole-cluster-inputs.pickle'),
        'sourceField':('selectedItems',0,'vector')})
    controls={name:not array_same(original,v) for name,v in changed.items()}
    try:h.native.same(original,copy.deepcopy(original))
    except ValueError as error:original_failure='truth value of an array' in str(error)
    else:original_failure=False
    controls['originalReaderArrayBoundary']=original_failure
    controls['exactCopyAccepted']=array_same(original,copy.deepcopy(original))
    f.save(base/'array-reader-mutation-controls.json',controls)
    f.require(all(controls.values()),'actual array reader dtype/shape/value/container controls')


def load(base):
    repair=json.loads(REPAIR.read_text());oldbase=PREVIOUS/'complete';old=json.loads((oldbase/'inputs.json').read_text())
    STATE.update(base=base,repair=repair,old=old)
    f.require(f.digest(Path(__file__))==repair['newHelperSha256'] and f.digest(PLAN)==repair['newPlanSha256'],'reviewed recovery helper/plan')
    f.require(f.digest(Path(h.__file__))==repair['originalHelperSha256'] and body(h.main)==repair['prepareJoin']['originalMainAST'],
              'original helper/main immutable')
    f.require(f.digest(Path(modes_reader.__file__))==repair['arrayReaderSha256'] and body(array_same)==repair['prepareJoin']['arrayComparatorAST'],
              'existing accepted whole array comparator')
    f.require(f.digest(oldbase/'inputs.json')==repair['originalInputsSha256'],'completed original loader')
    actual={str(p.relative_to(oldbase)) for p in oldbase.rglob('*') if p.is_file()
            and 'source' not in p.relative_to(oldbase).parts and p!=oldbase/'inputs.json'}
    f.require(actual==set(repair['completedFiles']),'all completed failed catalogue files')
    manifest=dict(old,runDirectory=str(base),sourceFiles=dict(old['sourceFiles']),inputPackets=dict(old['inputPackets']),referencedInputs={})
    for name,item in repair['completedFiles'].items():
        path=oldbase/name
        f.require(f.digest(path)==item['sha256'] and path.stat().st_size==item['bytes'] and str(path.resolve())==item['resolved']
                  and (str(path.readlink()) if path.is_symlink() else None)==item['rawLink'],('completed original file identity',name))
        source.reference(base,manifest,path,name,item['sha256'])
    source.reference(base,manifest,oldbase/'inputs.json','original-end-inputs.json',repair['originalInputsSha256'])
    for name,item in repair['originalLogs'].items():
        source.reference(base,manifest,PREVIOUS/name,'original-end-input-logs/'+name,item['sha256'])
    for path in (Path(__file__).resolve(),PLAN,REPAIR):
        name=str(path.relative_to(f.ROOT));f.require(name not in manifest['sourceFiles'],'distinct recovery source namespace')
        manifest['sourceFiles'][name]=f.digest(path)
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==value,'current source prehash')
        if name in old['sourceFiles']:f.require(f.digest(oldbase/'source'/name)==value,'original frozen source prehash')
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,target)
        f.require(f.digest(target)==value,'recovery frozen source')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'original/reference input prehash')
    guard=json.loads((base/'original-end-input-logs/resource-guard/outcome.json').read_text())
    inv=json.loads((base/'original-end-input-logs/frequency_end_inputs.invocation.json').read_text())
    trace=(base/'original-end-input-logs/frequency_end_inputs.stderr').read_text()
    f.require(guard['exitCode']==guard['childOutcome']['exitCode']==inv['exitCode']==1
              and guard['childOutcome']['guardReason'] is None and 'line 204, in prepare' in trace
              and 'truth value of an array' in trace,'actual failed callsite after saved first-end prefix')
    manifest['arrayReaderRecovery']={'originalRun':str(PREVIOUS),'originalFailure':repair['failure'],
        'originalInputsSha256':repair['originalInputsSha256'],'completedFilesReused':len(repair['completedFiles']),
        'prepareJoin':STATE['prepareJoin'],'comparatorSource':str(Path(modes_reader.__file__)),
        'comparatorSha256':repair['arrayReaderSha256'],'scienceOrCompletedCaseRepeated':False}
    f.save(base/'inputs.json',manifest)
    f.save(base/'completed-end-input-reference-reuse.json',manifest['referencedInputs'])
    f.save(base/'end-input-array-reader-join.json',manifest['arrayReaderRecovery'])
    return manifest,tuple(repair['cases'])


def prohibit():
    original_prohibit()
    # Original native source inspection has now finished; controls run only
    # after every scientific producer has been disabled.
    array_controls(STATE['base'])


if __name__=='__main__':
    prepare,join=prepare_adapter();STATE['prepareJoin']=join
    repair=json.loads(REPAIR.read_text());f.require(join==repair['prepareJoin'],'reviewed whole prepare/main/source join')
    h.inputs.protect_references=protect;h.load=load;h.prepare=prepare;h.same=same;h.prohibit=prohibit
    h.main()
