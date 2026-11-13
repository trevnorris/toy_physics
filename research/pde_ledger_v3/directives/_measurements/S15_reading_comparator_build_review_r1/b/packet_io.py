from pathlib import Path
import gzip,hashlib,json,tarfile,sys
B=Path(__file__).resolve().parent
P=B.parent/'packet'
R=Path.cwd()
S=P/'research/pde_ledger_v3/scripts'
DETAIL='S15_reading_comparator_details.out'
archive=S/'out'/DETAIL
sys.path.insert(0,str(S))
def sha(x): return hashlib.sha256(x).hexdigest()
def member(tf,name):
    data=tf.extractfile(name).read()
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)
class Records:
    def __init__(self):
        self.tar=tarfile.open(archive,'r:')
        self.cfg=member(self.tar,'configuration.json')
        self.idx={r['id']:r for r in member(self.tar,'coverage.json.gz')}
        self.f={stream:(R/i['path']).open('rb') for stream,i in self.cfg['inputs'].items()}
        self.n={}
        for r in self.idx.values():
            self.n.setdefault((r['id'].split(':')[0],r['tag']),[]).append(r['id'])
    def raw(self,rid):
        r=self.idx[rid]; stream=rid.split(':')[0]; f=self.f[stream]; span=r['span']
        f.seek(span['offset']); payload=f.read(span['length'])
        assert payload.endswith(b'\n') and sha(payload[:-1])==r['payload'],rid
        return payload[:-1].decode()
    def get(self,e,t):
        rid,=self.n[e,t]
        return self.raw(rid)
    def obj(self,i):
        return member(self.tar,f'objects/{i:03d}.json.gz')
    def ref(self,name,field):
        return self.row(name,field)['left'][0],self.row(name,field)['right'][0]
    def row(self,name,field):
        obj=self.obj(self.cfg['required_objects'].index(name))
        return next(r for r in obj['rows'] if r['field']==field and r['role']=='cross')
