import ast
import gzip
import io
import re
import json
import hashlib
import tarfile
from collections import Counter,defaultdict
from pathlib import Path

LOCAL=Path(__file__).resolve().parent
PACKET=LOCAL.parent/'packet/research/pde_ledger_v3'
REPO=Path('/var/projects/toy_physics')
manifest=json.loads((PACKET/'directives/_measurements/S15_reading_comparator_run_manifest.json').read_text())
files={k:REPO/v['path'] for k,v in manifest['repository_inputs'].items()}

def scan():
    indices={}
    for k,v in manifest['repository_inputs'].items():
        index=defaultdict(list)
        hashed=hashlib.sha256(); stop=0
        with files[k].open('rb') as f:
            current=None
            while stop < v['bytes']:
                line=f.readline()
                if not line:raise RuntimeError('truncated')
                hashed.update(line)
                h, sep, tail=line.partition(b': ')
                if sep and re.fullmatch(rb'[A-Za-z][a-zA-Z0-9_{},+ /-]*', h):
                    current=[stop+len(h)+2,len(tail)]
                    index[h.decode()].append(current)
                else:
                    current[1]+=len(line)
                stop+=len(line)
        assert stop==v['bytes'] and hashed.hexdigest()==v['sha256'],(k,stop,hashed.hexdigest())
        indices[k]=index
    return indices

def text_at(k,ref):
    with files[k].open('rb') as f:
        f.seek(ref[0]);return f.read(ref[1]).decode().rstrip('\n')

def raw(indices,k,tag):
    ref,=indices[k][tag]
    return text_at(k,ref)

def canon(s):
    return re.sub(r'labels_\{([^}]+)\}',lambda m:'labels_'+m[1].replace(', ','+'),s)

def roster():
    labs=['drain1','drain2','remainder']
    import itertools
    groups=['+'.join(g) for i in (1,2,3) for g in itertools.combinations(labs,i)]
    names=set()
    for i in (1,2,3):
        comp='component_'+str(i)
        for surface in ('inner','outer'):
            for piece in ('other','pressure','viscous','carried','transport_centre','transport_radius'):
                key=f'{surface}/{piece}/{comp}'; names.add(key)
                if piece=='other':continue
                prefix='/labels_'
                if piece=='pressure':
                    names.update([key+'/first_order_choice',key+'/whole_minus_expansion']); prefix='/expansion_labels_'
                names.update(key+prefix+g for g in groups)
        for surface in ('inner','outer','both_surfaces'):
            for piece in ('rest','transport','crossing'):names.add(f'{surface}/sum_{piece}/{comp}')
        for q in ['momentum','momentum_rate','integrated_force']+['momentum_labels_'+g for g in groups]+['rate_labels_'+g for g in groups]:
            names.add(f'region/{comp}/{q}')
    return names

def detail(tar,name):
    data=tar.extractfile(name).read()
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)

if __name__=='__main__':
    indices=scan()
    for e in ('py','wl'):
        k=e+'_reading';orig=e+'_original'
        names={canon(tag[8:-2]):tag[8:-2] for tag in indices[k] if tag.startswith('reading/') and tag.endswith('/Q')}
        direct={q for q,v in names.items() if v in indices[orig] and raw(indices,k,'reading/'+v+'/Q')==raw(indices,orig,v)}
        print(k,'bytes',manifest['repository_inputs'][k]['bytes'],'roster/missing/extra/direct/sums',
              len(names),len(roster()-set(names)),len(set(names)-roster()),len(direct),sum('/sum_' in q for q in set(names)-direct))
        for q in ('inner/other/component_1','outer/other/component_1'):
            print(e,q,raw(indices,orig,q))
        for q in ('inner/outward_normal','outer/outward_normal','inner/area_measure','outer/area_measure'):
            if q in indices[orig]:print(e,q,raw(indices,orig,q))
    with tarfile.open(PACKET/'scripts/out/S15_reading_comparator_details.out') as tar:
        names=detail(tar,'configuration.json')['required_objects']
        for name in ['both_surfaces/sum_crossing/component_1','inner/other/component_1','region/component_1/rate_labels_drain2']:
            o=detail(tar,f'objects/{names.index(name):03d}.json.gz')
            print('ARCHIVE',name,[(r['field'],r['status']) for r in o['rows'] if r['field'] in ['original_Q','Q','E_geo','E_resp','velocity_part_1/residual','velocity_part_2/residual','item3/family_integral','item3/object']])
        for name in ['inner/other/component_1','both_surfaces/sum_crossing/component_1']:
            obj=detail(tar,f'objects/{names.index(name):03d}.json.gz')
            for row in obj['rows']:
                if row['field']=='E_resp' and row['role']=='cross':
                    print('ARCHIVE RESPONSE',name,'ids',row['left'],row['right'],'first',row.get('samples',[None])[0])
        aux=detail(tar,'context_checks.json.gz')
        print('domain endpoints',*[q for q in aux['kernels'] if 'positive endpoint' in q['field']])
        print('angular example',next(q for q in aux['kernels'] if 'angular moment' in q['field']))
        print('argument sample',[q for q in aux['argument_slots'] if q.get('status')=='no_wl_declared_slot'][:1])
