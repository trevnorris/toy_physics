"""Stdlib-only restoration/storage tests; no scientific imports or payloads."""
import ast
import hashlib
import json
import math
from pathlib import Path
import pickle
import shutil
import sys
import tempfile
import time
import os
import types
import zipfile
from S11c_d_numerical_radiating_blob_store import BlobStore

M=Path(__file__).parent
WORKER=M/'S11c_d_numerical_radiating_boundary_continue_v2.py'
OLD=M/'S11c_d_numerical_radiating_boundary.py'
GUARD=M/'S11c_d_numerical_radiating_boundary_continue_v2_guard.py'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class FalseEquality:
    def __eq__(self,other):return False

def main():
    tree=ast.parse(WORKER.read_text());prior_tree=ast.parse(OLD.read_text())
    checks={}
    old={n.name:n for n in prior_tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
    new={n.name:n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
    for name in ('polynomial','Pair','bind_source','trace_path','original_seed','synthetic_check','compare','projector'):
        assert ast.dump(old[name])==ast.dump(new[name]),name
    prior_table={n.name:n for n in old['Table'].body if isinstance(n,ast.FunctionDef)}
    current_table={n.name:n for n in new['Table'].body if isinstance(n,ast.FunctionDef)}
    for name in prior_table:assert ast.dump(prior_table[name])==ast.dump(current_table[name])
    # run differs only in the path argument needed to join the old copied packet.
    old_run=ast.unparse(old['run'])
    new_run=ast.unparse(new['run'])
    old_restore="inputs[key] = J.call('restore/' + key, {'path': str(target), 'sha256': item['sha256']}, lambda p=target: pickle.loads(p.read_bytes()))"
    new_restore="restore_args = {'path': str(Path(manifest['priorComplete']['path']) / 'inputs' / target.name), 'sha256': item['sha256']}\n        inputs[key] = J.call('restore/' + key, restore_args, lambda p=target: pickle.loads(p.read_bytes()))"
    assert old_restore in old_run and new_restore in new_run
    assert old_run.replace(old_restore,'SOURCE_RESTORE')==new_run.replace(new_restore,'SOURCE_RESTORE')
    checks['scientificEquationPathAndRunASTUnchangedExceptRestorePath']=True
    calls=[ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)]
    assert not any(x in calls for x in ('signal.alarm','signal.setitimer'))
    guard=GUARD.read_text();assert "'--property=RuntimeMaxSec=infinity'" in guard and "'--property=Restart=no'" in guard
    assert "reason = 'wall-time limit'" not in guard and 'progressStall' not in guard
    checks['noWorkerAlarmServiceOrInactivityDeadline']=True
    nodes=[new[n] for n in ('require','save_json','same_saved','Journal')]
    fake_np=types.SimpleNamespace(ndarray=type('SyntheticArray',(),{}),number=type('SyntheticNumber',(),{}),isnan=math.isnan)
    space={'np':fake_np,'sp':types.SimpleNamespace(true=object()),'hashlib':hashlib,'Path':Path,'json':json,'pickle':pickle,'zipfile':zipfile,'shutil':shutil,'time':time,'os':os,'BlobStore':BlobStore,'digest':sha,'numeric_json':lambda v:v}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(WORKER),'exec'),space)
    Journal=space['Journal'];same=space['same_saved']
    assert same({'x':(1,b'v')},{'x':(1,b'v')}) and not same({'x':(1,b'v')},{'x':[1,b'v']})
    checks['exactContainerJoinRejectsChangedStructure']=True
    with tempfile.TemporaryDirectory(prefix='s11c-restore-synthetic-') as directory:
        root=Path(directory);prior=root/'prior';prior.mkdir();records=[]
        def blob(z,name,value):
            raw=pickle.dumps(value,protocol=5);z.writestr(name,raw)
            return {'member':name,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
        with zipfile.ZipFile(prior/'operations.zip','w') as z:
            for i in range(5001):
                name='synthetic/'+str(i)
                records.append({'name':name,'status':'COMPLETE','input':blob(z,name+'/input.pickle',{'index':i,'container':(2,b'k'),'opaque':FalseEquality()}),'return':blob(z,name+'/return.pickle',{'value':i+1}),'seconds':0})
            pending={'name':'unfinished','input':blob(z,'unfinished/input.pickle',{'index':5001}),'status':'STARTED'}
            extra=blob(z,'synthetic-context.pickle',{'context':'opaque'})
        (prior/'operation-index.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in records))
        (prior/'active-operation.json').write_text(json.dumps(pending))
        files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in prior.iterdir()}
        manifest={'priorComplete':{'path':str(prior),'files':files,'extraContext':extra}}
        bad=root/'bad';bad.mkdir();j=Journal(bad,manifest)
        def forbidden():raise AssertionError('A completed function was replayed')
        try:j.call('synthetic/0',{'index':999},forbidden)
        except ValueError:pass
        else:raise AssertionError('Argument mismatch accepted')
        assert (bad/'argument-mismatch.json').exists() and j.count==0
        j.finish();checks['mismatchSavedAndStoppedBeforeFunction']=True
        good=root/'good';good.mkdir();j=Journal(good,manifest)
        for i in range(5001):assert j.call('synthetic/'+str(i),{'index':i,'container':(2,b'k'),'opaque':FalseEquality()},forbidden)=={'value':i+1}
        assert j.restored==j.count==5001
        j.blob(extra['member'],{'context':'opaque'})
        assert j.call('unfinished',{'index':5001},lambda:{'fresh':True})=={'fresh':True}
        assert j.call('next',{'index':5002},lambda:7)==7
        assert j.restored==5001 and j.count==5003
        j.finish()
        current=[json.loads(line) for line in (good/'operation-index.jsonl').read_text().splitlines()]
        assert all(r['status']=='RESTORED_PRIOR_COMPLETE_RETURN' and r['functionCalled'] is False for r in current[:5001])
        assert all(r['status']=='COMPLETE' for r in current[5001:])
        assert json.loads((good/'incomplete-input-join.json').read_text())['equal'] is True
        assert all(sha(prior/p)==sha(good/'prior-complete'/p)==v['sha256'] for p,v in files.items())
        checks['all5001ReturnsRestoredWithoutCallingFunctions']=True
        assert all(r['argumentJoin']=='EXACT_PRIOR_SERIALIZED_BYTES' for r in current[:5001])
        checks['exactBytesAcceptFalseObjectEqualityWithoutChangingInputs']=True
        checks['pendingInputJoinedBeforeFirstNewReturn']=True
        checks['originalArchiveCopiedByteIdentically']=True
    print(json.dumps({'status':'PASS','checks':checks,'scientificImportsOrPayloads':False,'workerSha256':sha(WORKER),'guardSha256':sha(GUARD),'testSha256':sha(__file__)},indent=2))

if __name__=='__main__':main()
