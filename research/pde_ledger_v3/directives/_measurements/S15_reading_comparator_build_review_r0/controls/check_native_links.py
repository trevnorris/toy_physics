"""Check the referenced signed-sum inputs directly in the four native output files."""
from collections import Counter
from contextlib import ExitStack
import gzip
import hashlib
import json
from pathlib import Path
import re
import tarfile

root=Path(__file__).resolve().parents[5]
packet=Path(__file__).resolve().parents[1]/'packet'
archive=packet/'research/pde_ledger_v3/scripts/out/S15_reading_comparator_details.out'
roles=Counter()
by_stream=Counter()
with tarfile.open(archive) as tar, ExitStack() as files:
    def read(name):
        data=tar.extractfile(name).read()
        return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)
    cfg=read('configuration.json')
    coverage={x['id']:x for x in read('coverage.json.gz')}
    sources={k:files.enter_context((root/v['path']).open('rb')) for k,v in cfg['inputs'].items()}
    for rec in coverage.values():
        if 'payload' not in rec:continue
        side=rec['id'].split(':')[0]
        start,length=rec['span']['offset'],rec['span']['length']
        assert start+length<=cfg['inputs'][side]['bytes']
        f=sources[side]
        f.seek(start)
        native=f.read(length)
        assert native[-1:]==b'\n'
        assert hashlib.sha256(native[:-1]).hexdigest()==rec['payload']
        by_stream[side]+=1
    for member in sorted(tar.getnames()):
        if not re.fullmatch(r'objects/\d{3}\.json\.gz',member):continue
        obj=read(member)
        m=re.fullmatch(r'(inner|outer|both_surfaces)/sum_(rest|transport|crossing)/(component_[123])',obj['object'])
        if not m:continue
        location,kind,axis=m.groups()
        surfaces=['inner','outer'] if location=='both_surfaces' else [location]
        pieces=[]
        if kind in ('crossing','rest'):pieces+=['carried','pressure','viscous','other']
        if kind in ('crossing','transport'):pieces+=['transport_centre','transport_radius']
        expected=Counter('/'.join([s,p,axis]) for s in surfaces for p in pieces)
        for row in obj['rows']:
            if row['field']!='original_Q':continue
            engine=row['role'].split('_')[0]
            assert coverage[row['left'][0]]['tag']=='reading/'+obj['object']+'/Q'
            assert all(i.startswith(engine+'_original:') for i in row['right'])
            assert Counter(coverage[i]['tag'] for i in row['right'])==expected
            roles[engine,kind]+=1
print('Payload digests match native bytes before their respective accepted boundaries:',dict(sorted(by_stream.items())))
print('Each referenced sum contains exactly its native component, surface(s), and original pressure/transport terms:',sorted(roles.items()))
