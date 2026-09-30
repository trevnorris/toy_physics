"""Saved-return continuation storage; no calculation functions are replayed."""
import json
import os
from pathlib import Path
import shutil
import S11c_d_numerical_radiating_end_maps_v2 as storage
from S11c_d_numerical_radiating_blob_store import BlobStore


def exact_structure(a,b):
    if type(a) is not type(b):return False
    if isinstance(a,dict):return a.keys()==b.keys() and all(exact_structure(a[k],b[k]) for k in a)
    if isinstance(a,(tuple,list)):return len(a)==len(b) and all(exact_structure(x,y) for x,y in zip(a,b))
    if hasattr(a,'dtype') and hasattr(a,'shape'):
        return a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes()
    return bool(a==b)


class SavedOperationsJournal(storage.Journal):
    def __init__(self,base,prior):
        super().__init__(base);source=Path(prior['path']);copies={}
        for rel,item in prior['files'].items():
            p=source/rel;storage.require(storage.digest(p)==item['sha256'],'prior file pin '+rel)
            target=base/'prior-complete'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
            storage.require(storage.digest(target)==item['sha256'],'prior byte copy '+rel);copies[rel]=item
        storage.save(base/'prior-copy-index.json',copies)
        rows=[json.loads(x) for x in (source/'operation-index.jsonl').read_text().splitlines()]
        storage.require(len(rows)==prior['completeOperations'] and all(v['status'] in ('COMPLETE','RESTORED_PRIOR_COMPLETE_RETURN') for v in rows),'all prior complete operations')
        self.saved_rows={v['name']:v for v in rows};self.saved_order=list(self.saved_rows);self.saved_store=BlobStore(base/'prior-complete/operations.sqlite')
        self.pending=json.loads((source/'active-operation.json').read_text());storage.require(self.pending['status']=='STARTED','one unfinished operation input')
        for name,h,n in self.saved_store.connection.execute('SELECT name,sha256,bytes FROM blobs'):
            if name==self.pending['input']['member']:continue
            raw=self.saved_store.get({'member':name,'sha256':h,'bytes':n});self.store.put(name,raw)
        exclude={'checks.json','failure.json','posthashes.json','artifact-index.json','active-operation.json','failed-check-evidence.json'}
        for path in source.glob('*.json'):
            if path.name not in exclude and not path.name.endswith('copy-index.json'):shutil.copyfile(path,base/path.name)
    def join(self,old,actual):
        if old==actual:return 'EXACT_SERIALIZED_BYTES'
        storage.require(exact_structure(storage.decode(old),storage.decode(actual)),'actual saved input structural identity')
        return 'EXACT_STRUCTURAL_IDENTITY'
    def call(self,name,args,fn):
        if name not in self.saved_rows:
            storage.require(self.restored==len(self.saved_rows),'complete prior prefix restored')
            if self.pending is not None:
                storage.require(name==self.pending['name'],'first resumed unfinished operation')
                import pickle
                actual=pickle.dumps(args,protocol=5);join=self.join(self.saved_store.get(self.pending['input']),actual)
                storage.save(self.base/'unfinished-input-identity.json',{'operation':name,'join':join,'priorInput':self.pending['input']});self.pending=None
            return super().call(name,args,fn)
        storage.require(name==self.saved_order[self.restored],'prior complete operation order')
        row=self.saved_rows[name];self.active=name;actual=self.blob('restore-argument-checks/'+name+'.pickle',args)
        storage.save(self.base/'active-operation.json',{'name':name,'actualInput':actual,'status':'RESTORING'})
        join=self.join(self.saved_store.get(row['input']),self.store.get(actual));result=storage.decode(self.saved_store.get(row['return']))
        record={'name':name,'input':row['input'],'actualInput':actual,'return':row['return'],'status':'RESTORED_PRIOR_COMPLETE_RETURN','functionCalled':False,'argumentJoin':join,'prior':row}
        with (self.base/'operation-index.jsonl').open('a') as f:f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
        self.refs[name]=row['return'];self.count+=1;self.restored+=1;self.active=None;storage.save(self.base/'active-operation.json',record)
        return result
