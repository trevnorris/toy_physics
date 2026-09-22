#!/usr/bin/env python3
"""Continue saved end-call inputs after a cached-method source-file reader fault."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import sys
import textwrap

sys.dont_write_bytecode=True
import S11c_d_remaining_case_frequency_end_calls as original

f=original.f
PREVIOUS=f.STORE/'s11c-remaining-case-frequency-20260921/end-call-inputs'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_calls_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_end_calls_reader_repair.json'
NATIVE_BODY=original.body
NATIVE_LOAD=original.load


def adapter():
    text=textwrap.dedent(inspect.getsource(NATIVE_BODY));before=ast.parse(text).body[0]
    node=copy.deepcopy(before);changed=[]
    class Reader(ast.NodeTransformer):
        def visit_Call(self,n):
            if ast.unparse(n.func)=='inspect.getsourcefile':
                f.require(len(n.args)==1 and ast.unparse(n.args[0])=='value' and not n.keywords,
                          'exact original source-file reader operand')
                changed.append(True)
                n.args[0]=ast.Call(ast.Attribute(ast.Name('inspect',ast.Load()),'unwrap',ast.Load()),[n.args[0]],[])
            return self.generic_visit(n)
    node=Reader().visit(node);f.require(len(changed)==1,'one decorated-object source-file reader')
    restored=copy.deepcopy(node)
    class Undo(ast.NodeTransformer):
        def visit_Call(self,n):
            if ast.unparse(n.func)=='inspect.getsourcefile':
                f.require(ast.unparse(n.args[0])=='inspect.unwrap(value)','exact inspected wrapper')
                n.args=n.args[0].args
            return self.generic_visit(n)
    restored=Undo().visit(restored)
    f.require(ast.dump(restored)==ast.dump(before),'whole original body reverse AST')
    namespace=dict(vars(original));module=ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[]))
    exec(compile(module,'<saved-end-call-source-file-reader>','exec'),namespace)
    return namespace['body'],{'wholeOriginalBodyAST':hashlib.sha256(ast.dump(before).encode()).hexdigest(),
                             'onlyChange':'inspect.getsourcefile(inspect.unwrap(value))',
                             'wholeBodyReverseAST':True,'sourceTextIncludingDecoratorsUnchanged':True}


def load(base):
    repair=json.loads(REPAIR.read_text());old=PREVIOUS/'complete'
    f.require(f.digest(Path(original.__file__))==repair['originalHelperSha256'], 'original helper remains immutable')
    f.require(f.digest(old/'inputs.json')==repair['originalInputsSha256'],'completed input manifest')
    inventory_path=Path(repair['completedInventory']['path'])
    f.require(f.digest(inventory_path)==repair['completedInventory']['sha256'],'completed reference inventory')
    inventory=json.loads(inventory_path.read_text())
    prior=json.loads((old/'inputs.json').read_text());manifest=copy.deepcopy(prior)
    manifest['runDirectory']=str(base);manifest['referencedInputs']={}
    # The original loader completed; preserve its manifest and all exact
    # raw/resolved/source/hash routes without redoing its preparation.
    for name,item in inventory.items():
        path=old/name
        f.require(f.digest(path)==item['sha256'] and path.stat().st_size==item['bytes']
                  and str(path.resolve())==item['resolved']
                  and (str(path.readlink()) if path.is_symlink() else None)==item['rawLink'],
                  ('immutable completed file',name))
        destination='original-end-call-inputs.json' if name=='inputs.json' else name
        original.h.source.reference(base,manifest,path,destination,item['sha256'])
    for name,item in prior['referencedInputs'].items():
        saved=inventory[name]
        f.require(saved['rawLink']==item['original'] and saved['resolved']==item['resolvedOriginal']
                  and saved['sha256']==item['sha256'] and saved['bytes']==item['bytes'],
                  'full completed original reference route')
    for name,item in repair['previousLogs'].items():
        original.h.source.reference(base,manifest,PREVIOUS/name,str(Path('original-end-call-logs')/name),item['sha256'])
    for path,name in ((inventory_path,'completed-end-call-file-inventory.json'),
                      (Path(repair['outcomeInspection']['path']),'failed-end-call-outcome.json')):
        expected=repair['completedInventory']['sha256'] if path==inventory_path else repair['outcomeInspection']['sha256']
        original.h.source.reference(base,manifest,path,name,expected)
    for name,value in prior['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(old/'source'/name)==f.digest(base/'source'/name)==value,
                  'original current/frozen source identity')
    for path in (Path(__file__).resolve(),PLAN,REPAIR):
        name=str(path.relative_to(f.ROOT));f.require(name not in manifest['sourceFiles'],'distinct recovery source pin')
        value=f.digest(path);manifest['sourceFiles'][name]=value
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
        f.require(not target.exists() and not target.is_symlink(),'fresh frozen recovery source')
        target.write_bytes(path.read_bytes());f.require(f.digest(target)==value,'frozen recovery source hash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'fresh original input prehash')
    accepted=json.loads((base/'accepted-end-input-checkpoint.json').read_text())
    f.require(accepted['status']=='ACCEPTED_CASE_FREQUENCY_END_INPUTS'
              and f.digest(base/'accepted-end-input-checkpoint.json')==prior['acceptedEndInputs']['checkpointSha256'],
              'same accepted input checkpoint')
    fixed,join=adapter();original.body=fixed
    # The cached constructor is inspected only. Its actual wrapper and
    # underlying callable must resolve to the same full decorated source.
    value=original.engine.ClosedCurrentPairing.construct;unwrapped=inspect.unwrap(value)
    f.require(value is not unwrapped and hasattr(value,'__wrapped__'),'actual native cache wrapper route')
    f.require(inspect.getsource(value)==inspect.getsource(unwrapped),'same full original decorated body')
    joined=fixed(value)
    f.require(joined['path']==str(original.engine.HERE) and joined['sha256']==f.digest(original.engine.HERE),
              'actual native engine source file')
    join.update(originalHelperSha256=repair['originalHelperSha256'],
                originalMainAST=repair['originalMainAST'],originalPrepareAST=repair['originalPrepareAST'],
                cachedMethodSourceSha256=joined['sha256'],cachedMethodBodyAST=joined['wholeBodyAST'],
                actualWrapperType=type(value).__name__,nativeConstructorCalled=False)
    for name in ('main','prepare'):
        function=next(n for n in ast.parse(Path(original.__file__).read_text()).body if getattr(n,'name',None)==name)
        f.require(hashlib.sha256(ast.dump(function).encode()).hexdigest()==repair['original'+name.title()+'AST'],
                  'whole original main/prepare remains unchanged')
    manifest['sourceFileReaderRecovery']={'previousDirectory':str(PREVIOUS),
        'originalInputsSha256':repair['originalInputsSha256'],'completedFilesReused':len(inventory),
        'originalSourcePins':len(prior['sourceFiles']),'originalReferenceRoutes':len(prior['referencedInputs']),
        'previousOutcome':repair['previousOutcome'],'completedCatalogueCalls':0,
        'oldLoaderCalled':False,'originalNativeAlgorithmsChanged':False,'readerJoin':join}
    f.save(base/'source-file-reader-join.json',join)
    f.save(base/'completed-end-call-reference-reuse.json',manifest['sourceFileReaderRecovery'])
    f.save(base/'inputs.json',manifest)
    return manifest,tuple(accepted['cases'])


if __name__=='__main__':
    original.load=load
    original.main()
