"""Independent source/operand, grid, roster, provenance, and branch-coverage audit."""
import ast, collections, gzip, hashlib, itertools, json, re, subprocess, tarfile
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'packet'
R=next(d for d in Path(__file__).resolve().parents if (d/'.git').exists())
L=P/'research/pde_ledger_v3'
A=L/'scripts/out/S15_reading_comparator_details.out'
MAN=L/'directives/_measurements/S15_reading_comparator_run_manifest.json'
OLD=R/'research/pde_ledger_v3/scripts/out/preserved/S15_reading_comparator_details.build-r0.out'
H=lambda x:hashlib.sha256(x).hexdigest()
def member(tf,n):
    v=tf.extractfile(n).read()
    return json.loads(gzip.decompress(v) if n.endswith('.gz') else v)
def main():
    man=json.loads(MAN.read_text());tf=tarfile.open(A);cf=member(tf,'configuration.json')
    print('archive bytes, sha',A.stat().st_size,H(A.read_bytes()))
    assert H(A.read_bytes())==man['detail_archive']['sha256'] and H(OLD.read_bytes())==cf['reused_exact']['archive_sha256']
    assert H((P/'prior/S15_reading_comparator_run_manifest.json').read_bytes())==cf['reused_exact']['manifest_sha256']
    for k,h in man['source_sha256'].items(): assert H((L/'scripts'/k).read_bytes())==h
    for k,h in man['engine_sources'].items(): assert H((P/k).read_bytes())==h
    for k,a in man['accepted_artifacts'].items(): assert H(subprocess.check_output(['git','show',a['git_object']],cwd=R))==a['sha256']
    print('prior archive sha',H(OLD.read_bytes()),'pinned scripts engines artifacts',*[len(man[k]) for k in ('source_sha256','engine_sources','accepted_artifacts')])
    cov=member(tf,'coverage.json.gz');raw={};ids={r['id']:r for r in cov};assert len(ids)==len(cov)
    blob={v['sha256']:v['native_operand'].encode() for v in map(json.loads,gzip.GzipFile(fileobj=tf.extractfile('operands.jsonl.gz')))}
    for key,spec in cf['inputs'].items():
        file=(R/spec['path']); end=spec['bytes']
        with file.open('rb') as f: text=f.read(end); extra=f.read(1)
        assert H(text)==spec['sha256']
        assert bool(extra)==(key=='wl_reading')
        recs=[q for q in cov if q['id'].startswith(key+':')]
        positions=sorted((q['span']['offset']-len(q['tag'].encode())-2,q['span']['offset'],q) for q in recs)
        for i,(start,off,q) in enumerate(positions):
            assert text[start:off]==q['tag'].encode()+b': '
            nxt=positions[i+1][0] if i+1<len(positions) else end
            assert nxt==off+q['span']['length'] and text[nxt-1:nxt]==b'\n'
            if 'payload' in q:
                native=text[off:nxt-1]; assert H(native)==q['payload'] and native==blob[q['payload']]
        print('stream',key,end,spec['sha256'],len(recs),positions[-1][2]['tag'])
    shards=[member(tf,f'objects/{i:03}.json.gz') for i in range(336)]
    prior=tarfile.open(OLD)
    fields=collections.Counter(); sources=collections.Counter(); dense=collections.Counter()
    groups=['+'.join(q) for n in (1,2,3) for q in itertools.combinations(('drain1','drain2','remainder'),n)]
    required=set()
    for i in (1,2,3):
        for surf in ('inner','outer'):
            for kind in ('carried','pressure','viscous','other','transport_centre','transport_radius'):
                k=f'{surf}/{kind}/component_{i}';required.add(k)
                if kind=='pressure':
                    required.update([k+'/first_order_choice',k+'/whole_minus_expansion',*[k+'/expansion_labels_'+q for q in groups]])
                elif kind!='other':required.update(k+'/labels_'+q for q in groups)
        required.update(f'{surf}/sum_{kind}/component_{i}' for surf in ('inner','outer','both_surfaces') for kind in ('rest','transport','crossing'))
        k=f'region/component_{i}/';required.update(k+q for q in ('momentum','momentum_rate','integrated_force'))
        required.update(k+p+g for g in groups for p in ('momentum_labels_','rate_labels_'))
    assert required==set(cf['required_objects'])
    for obj in shards:
        k=obj['object']; old=member(prior,f"objects/{obj['index']:03}.json.gz")
        oldrows={(tuple(r['left']),tuple(r['right'])):r for r in old['rows'] if r.get('residual')=='0'}
        for i,row in enumerate(obj['rows']):
            if row.get('provenance'):
                prv=oldrows[tuple(row['left']),tuple(row['right'])]
                assert prv['status']==row['status']=='exact_identity' and prv['method']==row['method'];sources['prior']+=1
            if 'same_comparison_as' in row:
                earliest=obj['rows'][row['same_comparison_as']]
                assert row['same_comparison_as']<i and earliest['status']==row['status']
                for side in ('left','right'): assert [(x[:2],ids[x]['payload']) for x in row[side]]==[(x[:2],ids[x]['payload']) for x in earliest[side]]
                sources['alias']+=1
            if 'samples' in row:
                for s in row['samples']:
                    x,y=s['values'][-1];assert x-y==s['residual']
                sources['initial samples']+=len(row['samples'])
            if row['role'].endswith('_same_engine'):
                left,=row['left'];right=row['right']
                if '/sum_' in k:
                    surf,kind,i=k.split('/');kind=kind.removeprefix('sum_')
                    sfs=['inner','outer'] if surf=='both_surfaces' else [surf]
                    pcs=({'rest':['carried','pressure','viscous','other'],'transport':['transport_centre','transport_radius'],
                         'crossing':['carried','pressure','viscous','other','transport_centre','transport_radius']})[kind]
                    assert set(ids[r]['tag'] for r in right)=={f'{t}/{p}/{i}' for t in sfs for p in pcs};sources['sum links']+=1
                else: assert len(right)==1 and ids[left]['payload']==ids[right[0]]['payload'];sources['direct links']+=1
            if row['item']==6:
                for rid in row['left']+row['right']:
                    st=blob[ids[rid]['payload']].decode().replace('{','[').replace('}',']')
                    val=ast.literal_eval(st);dense[str(tuple(tuple(x) for x in val))]+=1
            fields[row['field'],row['status']]+=1
    print('roster',len(required),'counts',dict(sources),'dimension native collections',dict(dense))
    aux=member(tf,'context_checks.json.gz')
    print('argument-slot gaps',[(q['argument'],q['expected']) for q in aux['argument_slots'] if q['status']=='no_wl_declared_slot'])
    print('Wolfram source registration',*[f'{i} {s}' for i,s in enumerate((L/'mathematica/S15_separation_departure.wl').read_text().splitlines(),1) if 'Do[register[centreVelocityChange' in s])
    obj=shards[19]
    for name in ('item3/exponent_branch_inventory','item3/Q_closed_operand','item3/object'):
        row=next(q for q in obj['rows'] if q['field']==name)
        if name.endswith('inventory'):
            for rid in row['left']+row['right']: print('branch inventory',rid,blob[ids[rid]['payload']].decode())
        else:
            origin=obj['rows'][row['same_comparison_as']] if 'same_comparison_as' in row else row
            planes=[list(map(F,g)) for g in obj['exponent_loci']]
            point=[F(-1),F(9,4),F(-2),F(9,2)]
            hit=lambda p:tuple(i for i,h in enumerate(planes) if h[4]+sum(a*b for a,b in zip(h[:4],p))==0)
            observed={hit(list(map(F,q['exponents']))) for q in origin['samples']}
            print('sampled simultaneous branch evidence',obj['object'],name,row['left'],row['right'],'individual sampled plane sets',sorted(observed),'specified point',list(map(str,point)),hit(point),'covered',hit(point) in observed)
    selected=member(tf,'refinement/019.json.gz')['fields'][0]['rechecks'][0]
    print('first precision rows',selected)
if __name__=='__main__': main()
