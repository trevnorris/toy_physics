import ast,collections,gzip,hashlib,itertools,json,re,subprocess,tarfile
from packet_io import *
d=Records(); cf=d.cfg; ar=d.tar
mf=json.loads((S.parent/'directives/_measurements/S15_reading_comparator_run_manifest.json').read_text())
for key,path in [('detail_archive',archive),('transcript',S/'out/S15_reading_comparator.out')]:
    m=mf[key]; print(key,path.stat().st_size,sha(path.read_bytes())); assert path.stat().st_size==m['bytes'] and sha(path.read_bytes())==m['sha256']
for k,v in mf['accepted_artifacts'].items():
    h=sha(subprocess.check_output(['git','show',v['git_object']]));print(k,h);assert h==v['sha256']
for name,raw in cf['inputs'].items():
    f=d.f[name]; blob=f.read(raw['bytes']); assert len(blob)==raw['bytes'] and sha(blob)==raw['sha256']
    print('native',name,len(blob),sha(blob),'last record',blob.rstrip().split(b'\n')[-1].decode(),'trailing bytes', (R/raw['path']).stat().st_size-len(blob))
    # Independent physical-line scanner with separate tag grammar and payload endpoints.
    brackets=re.compile(rb'^([a-zA-Z][a-zA-Z0-9_,+{} /-]*): ',re.M)
    cnt=collections.Counter(); matches=list(brackets.finditer(blob))
    for j,mt in enumerate(matches):
        tag=mt[1].decode(); pos=mt.end(); end=matches[j+1].start() if j+1<len(matches) else len(blob)
        r=d.idx[f'{name}:{pos}']; assert (r['tag'],r['occurrence'],r['span']['length'])==(tag,cnt[tag],end-pos)
        if 'payload' in r: assert sha(blob[pos:end-1])==r['payload']
        cnt[tag]+=1
    print('occurrences',len(matches),'multiline',sum(blob.count(b'\n',m.end(),matches[i+1].start() if i+1<len(matches) else len(blob))>1 for i,m in enumerate(matches)),'distinct tags',len(cnt))
    del blob
wanted={r['payload'] for r in d.idx.values() if 'payload' in r}
with gzip.GzipFile(fileobj=ar.extractfile('operands.jsonl.gz')) as f:
    for ln in f:
        ent=json.loads(ln); h=sha(ent['native_operand'].encode());assert h==ent['sha256'];wanted.remove(h)
assert not wanted
pg=P.parent.parent/'r0/packet/research/pde_ledger_v3/scripts/out'/DETAIL
old=json.loads((P/'prior/S15_reading_comparator_run_manifest.json').read_text()); assert sha(pg.read_bytes())==mf['reused_exact']['archive_sha256']==old['detail_archive']['sha256']
print('prior',sha(pg.read_bytes()))
eligible={}; reused=collections.Counter(); samps=collections.Counter(); cases=collections.Counter(); rows=collections.Counter(); inventory=collections.Counter()
oldtar=tarfile.open(pg); a=member(oldtar,'configuration.json'); assert a['inputs']==cf['inputs']
parts=['carried','pressure','viscous','other','transport_centre','transport_radius']
groups=['+'.join(g) for n in [1,2,3] for g in itertools.combinations(('drain1','drain2','remainder'),n)]
roster=set()
for c in ('component_1','component_2','component_3'):
    for side,piece in itertools.product(['inner','outer'],parts):
        stem=f'{side}/{piece}/{c}';roster.add(stem)
        if piece in ('carried','viscous','transport_centre','transport_radius'):roster.update(stem+'/labels_'+g for g in groups)
        if piece=='pressure':roster.update(stem+'/expansion_labels_'+g for g in groups);roster.update(stem+'/'+p for p in ['first_order_choice','whole_minus_expansion'])
    roster.update(f'{side}/sum_{s}/{c}' for side,s in itertools.product(['inner','outer','both_surfaces'],['rest','transport','crossing']))
    roster.update('region/'+c+'/'+g for g in ('momentum','momentum_rate','integrated_force'))
    roster.update('region/'+c+'/'+p+g for g,p in itertools.product(groups,['momentum_labels_','rate_labels_']))
assert set(cf['required_objects'])==roster
for i,name in enumerate(cf['required_objects']):
    prev=member(oldtar,f'objects/{i:03d}.json.gz')
    for r in prev['rows']:
        if r.get('residual')=='0' and r.get('method') in ('identical native expression on same engine','exact algebra on translated operands') and r['status']=='exact_identity':eligible[r['left'][0],tuple(r['right'])]=r['method']
    o=d.obj(i); assert o['object']==name and o['index']==i and o['build_id']==cf['build_id']
    for j,r in enumerate(o['rows']):
        rows[r['role'],r['status']]+=1
        if 'provenance' in r:
            key=(r['left'][0],tuple(r['right']));assert eligible[key]==r['method'] and r['status']=='exact_identity';reused[r['method']]+=1
        ids=r['left']+r['right']
        if r['role'].endswith('same_engine'):
            if '/sum_' in name:
                loc,typ,c=name.split('/'); terms=parts[:4] if typ=='sum_rest' else parts[4:] if typ=='sum_transport' else parts
                expected=[f'{l}/{t}/{c}' for l in (('inner','outer') if loc=='both_surfaces' else (loc,)) for t in terms]
                assert [d.idx[x]['tag'] for x in r['right']]==expected
            else:assert d.raw(ids[0])==d.raw(ids[1])
            inventory[r['role'],len(r['right'])>1]+=1
        if 'same_comparison_as' in r:
            parent=o['rows'][r['same_comparison_as']];assert r['same_comparison_as']<j
            assert [d.idx[x]['payload'] for x in ids]==[d.idx[x]['payload'] for x in parent['left']+parent['right']]
        for s in r.get('samples',[]):
            samps[s['status']]+=1; cases[tuple(s['case'][:2])]+=1
            if 'values' in s:
                x,y=s['values'][-1]; assert x-y==s['residual']
    if i==0:
        for name2 in ('Q','E_geo','E_resp'):
            r=next(item for item in o['rows'] if item['field']==name2 and item['role']=='cross'); origin=o['rows'][r.get('same_comparison_as',o['rows'].index(r))]
            print('native sample',name,name2,'ids',r['left'],r['right'],'first',origin.get('samples',['exact'])[0])
print('roster',len(roster),'direct and signed-reference links',dict(inventory),'eligible',len(eligible),'reused',dict(reused),'samples',dict(samps),'case counts',dict(cases))
for tag in ('family/centre_substitution','profile_family/separation_chart','profile_family/chart_condition'):
    print('context',tag,d.get('py_reading',tag), d.get('wl_reading',{'family/centre_substitution':'family/separation_vector','profile_family/chart_condition':'profile_family/domain'}.get(tag,tag)))
cc=member(ar,'context_checks.json.gz')
for row in cc['kernels']:
    if row['field'].startswith(('positive endpoint','angular moment (3, 0, 0)')): print('kernel',row)
for tag in ('E_geo/dimension_operands','E_geo/dimension_reference_residuals','velocity_part_2/dimension_reference_residuals'):
    for n in ('region/component_1/momentum','both_surfaces/sum_crossing/component_1'):
        ids=d.ref(n,tag);print('dimension',n,tag,[(i,d.raw(i)) for i in ids])
for tag in ('departure_existence_condition','item3/existence_condition'):
    ids=d.ref('both_surfaces/sum_crossing/component_1',tag);print('existence',tag,[(i,d.raw(i)) for i in ids])
print('all scoped occurrences',len(set(d.idx)-{k for k,v in d.idx.items() if v['role']!='scoped'}),'unpaired kernel counts',sum(map(len,cc['unpaired_internal_kernel_occurrences'].values())))
print('refinement sample')
obj=member(ar,'refinement/019.json.gz')
for f in obj['fields']:
    if f['rechecks']:print(obj['object'],f['field'],f['rechecks'][0]);break
