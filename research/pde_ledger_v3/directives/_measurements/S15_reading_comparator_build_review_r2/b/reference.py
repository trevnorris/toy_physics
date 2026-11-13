"""Independent native inventory, exact affine geometry, archive sampling and slot audit."""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import gzip
import json
from pathlib import Path
import re
import tarfile

import sympy

BASE = Path('research/pde_ledger_v3')
PACKET = Path('_scratch/s15/reading_comparator_build_review/r2/packet')
ROOT = PACKET / BASE
ARCH = ROOT / 'scripts/out/S15_reading_comparator_details.out'
COV = ROOT / 'directives/_measurements/S15_reading_comparator_coverage.json'
READ = {'py': ('scripts', 192065031, '10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6'),
        'wl': ('mathematica',333550018,'75622128d04da2d934c2a52468505b4438fff9f17a902d7a6cd35993d769b6c4')}
COORDS = sympy.symbols('a_1 a_2 b_1 b_2')
FIELD = 'item3/object'
PY_EXPR = re.compile(r'(?:Ne|Eq)\(([^,()]+), ([-+]?\d+)\)')
WL_EXPR = re.compile(r'(?:, | && )([^,{}\n\[\]A-Za-z]*?(?:[ab][ab]\[[12]\][+*/ ()\-\d]*)+) == ([-\d]+)')
WL_INV = re.compile(r'([^,&{}]+) == ([-\d]+)')

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def plane(coefficients):
    values = [F(x) for x in coefficients]
    leading = next(v for v in values[:4] if v)
    return tuple(v/leading for v in values)

@lru_cache(None)
def affine(lhs, rhs):
    lhs = re.sub(r'([ab])\1\[([12])\]',r'\1_\2',lhs.strip())
    expr = sympy.Poly(sympy.sympify(lhs) - int(rhs), *COORDS)
    assert expr.total_degree() == 1
    terms = [expr.coeff_monomial(x) for x in COORDS] + [expr.TC()]
    return plane(map(str, terms))

def extract(expr, pattern, count):
    entries = pattern.findall(expr)
    assert len(entries) == count, (len(entries),count)
    return frozenset(affine(a,b) for a,b in entries)

def native(engine):
    folder, size, expected = READ[engine]
    obj, inventory, rec = {}, {}, {}
    h = sha256()
    with (BASE / folder / 'out/S15_separation_departure.out').open('rb') as f:
        while f.tell() < size:
            start = f.tell()
            line = f.readline()
            h.update(line)
            tag,sep,val = line.partition(b': ')
            if not sep: continue
            if tag in {b'family/centre_velocity_changes', b'family/velocity_changes',
                       b'reading/dimension_argument_audit', b'family/dimension_signatures_including_changes'}:
                rec[tag.decode()] = (start+len(tag)+2,val.decode().strip())
            if not tag.startswith(b'reading/') or not tag.endswith((b'/item3/object',b'/item3/exponent_branch_inventory')): continue
            name,field = tag.decode()[8:].split('/item3/')
            name = re.sub(r'\{([^}]*)\}',lambda m:m[1].replace(', ','+'),name)
            value = val.decode().strip()
            if field == 'object':
                if engine=='py':
                    flat = extract(value,PY_EXPR,value.count('Ne(')+value.count('Eq('))
                else: flat = extract(value,WL_EXPR,value.count('=='))
                assert name not in obj
                obj[name] = (flat,start+len(tag)+2)
            else:
                if engine=='py':
                    flat=extract(value,PY_EXPR,value.count('Ne(')+value.count('Eq('))
                else:
                    flat=extract(value.replace(' && ',','),WL_INV,value.count('=='))
                assert name not in inventory
                inventory[name]=flat
        assert f.tell() == size and h.hexdigest()==expected
    assert {k:v[0] for k,v in obj.items()} == inventory
    return obj, rec

def flat_witness(equations, planes):
    # Form a reduced augmented system E x = -e over the rationals.
    rows=[list(p[:4])+[-p[-1]] for p in equations]
    pivots=[]
    for col in range(4):
        k=len(pivots)
        pivot=next((j for j in range(k,len(rows)) if rows[j][col]),None)
        if pivot is None: continue
        rows[k],rows[pivot]=rows[pivot],rows[k]
        a=rows[k][col]
        rows[k]=[x/a for x in rows[k]]
        for j in range(len(rows)):
            if j!=k and rows[j][col]:
                c=rows[j][col]
                rows[j]=[x-c*y for x,y in zip(rows[j],rows[k])]
        pivots.append(col)
    if len(pivots)!=len(equations): return None
    free=[j for j in range(4) if j not in pivots]
    # On this flat the free variables are t,...,t^d. Each unintended equality
    # is a nonzero polynomial of degree <= d. More than d*m distinct t suffice.
    polys=[]
    for coeffs in planes:
        z=list(coeffs)
        for col,row in zip(pivots,rows):
            multiple=z[col]
            z[-1]+=multiple*row[-1]
            for j in range(4):
                z[j]-=multiple*row[j]
        polys.append([z[-1]]+[z[j] for j in free])
    implied=tuple(i for i,p in enumerate(polys) if not any(p))
    for t in range(len(free)*len(planes)+1):
        if all(any(p) is False or sum(c*t**i for i,c in enumerate(p))!=0 for p in polys): break
    else: raise AssertionError('no witness')
    return implied,len(pivots)

@lru_cache(None)
def arrangement(planes):
    realized={}
    # Every nonempty affine flat in 4-space has an independent defining
    # subfamily of at most four members. Dependent consistent subsets therefore
    # contribute no additional flat. Each closure has its own rational witness.
    for dim in range(min(len(planes),4)+1):
        for subset in combinations(planes,dim):
            w=flat_witness(subset,planes)
            if w is not None:
                code,codim=w
                assert realized.setdefault(code,codim)==codim
    return realized

@lru_cache(None)
def hits(planes, coords):
    p=tuple(F(i) for i in coords)+(F(1),)
    return tuple(i for i,row in enumerate(planes) if sum(a*b for a,b in zip(row,p))==0)

def member(t,name):
    data=t.extractfile(name).read()
    return json.loads(gzip.decompress(data) if name.endswith('gz') else data)

def observation_set(t,path,d,counters):
    row=next(r for r in d['rows'] if r['field']==FIELD)
    ids=row['left']+row['right']
    visited=set()
    while 'same_comparison_as' in row:
        idx=row['same_comparison_as']
        assert idx not in visited
        visited.add(idx)
        row=d['rows'][idx]
    if visited: counters['alias']+=1
    observations={}
    def incorporate(record,source):
        key=tuple(record['case'])
        newpoint=tuple(record['exponents'])
        if key in observations:
            assert observations[key][0]==newpoint
        observations[key]=(newpoint,record['status'])
        counters[source+':'+record['status']]+=1
    for s in row.get('samples',[]): incorporate(s,'base')
    for prefix in ('refinement','log_branches'):
        sname=path.replace('objects/',prefix+'/')
        if sname not in t.getnames(): continue
        add=member(t,sname)
        assert add['object']==d['object'] and add['build_id']==d['build_id']
        if prefix=='refinement':
            for r in add['fields']:
                if r['field']!=FIELD: continue
                assert ids==r['left']+r['right']
                for s in r['additional_samples']: incorporate(s,'additional')
                for s in r['rechecks']: incorporate(s,'high_precision')
        else:
            for s in add['rechecks']:
                if s['field']==FIELD: incorporate(s,'log')
    return ids,observations

def slots(py,wl,rows,source):
    line=next(l for l in source.splitlines() if 'Do[register[centreVelocityChange[j]' in l)
    decl=re.fullmatch(r'Do\[register\[centreVelocityChange\[j\],lenDim-timeDim,\{(\w+)\}\],\{j,(\d+)\}\];',line)
    assert decl
    unit=re.search(decl[1]+r'=\{(\d+),(\d+),(\d+)\}',source)
    dims=list(map(int,unit.groups()))
    reg=wl['family/velocity_changes']
    sign=wl['family/dimension_signatures_including_changes'][1]
    audit=py['reading/dimension_argument_audit']
    changes=py['family/centre_velocity_changes'][1]
    matched=[]
    for i in range(1,int(decl[2])+1):
        head=f'centreVelocityChange[{i}]'
        change=f'X_1_{i}_velocity_change'
        assert f'Derivative[1][xx[1, {i}]][t] -> {head}[t]' in reg[1]
        assert f'"{head}" -> {{0, 1, -1}}' in sign
        assert f'X_1_{i}_velocity(t): {change}(t)' in changes
        key=[change,0,'t']
        entry=next(r for r in rows if r['argument']==key)
        native=re.search(re.escape(repr(tuple(key)))+r': \(frozenset\(\{\(([^()]*)\)\}\), frozenset\(\{\(([^()]*)\)\}\)\)',audit[1])
        expected,observed=([int(n) for n in v.split(',')] for v in native.groups())
        assert dims==entry['wl_expected']==expected==observed, (head,dims,expected,observed)
        assert entry['expected']==[expected] and entry['observed']==[observed]
        assert entry['wl_slot']==1 and entry['wl_head']==head
        assert entry['wl_velocity_rules']==f'wl_reading:{reg[0]}' and entry['py_record']==f'py_reading:{audit[0]}'
        src=entry['wl_source']
        assert src['sha256']==sha256(source.encode()).hexdigest()
        assert src['declaration']==line==source.splitlines()[src['line']-1] and src['revision']=='cc62dc9f'
        assert entry['status']=='matches_source_declared_slot'
        matched.append((change,head,1,dims))
    return matched

def compare_sets(pub,planes,realized,points):
    got=tuple(plane(r) for r in pub['planes'])
    assert len(set(got))==len(got) and set(got)==set(planes), 'native plane mismatch'
    def convert(key):
        return tuple(sorted(planes.index(got[i]) for i in key))
    sampled=set(hits(planes,x) for x,status in points.values() if status in ('within_estimate','within_precision_estimate','separated'))
    assert sampled<=set(realized)
    untested=set(realized)-sampled
    for label,ref in [('sampled',sampled),('untested',untested)]:
        listed=[convert(i) for i in pub[label]]
        assert len(listed)==len(ref) and set(listed)==ref, (label,len(listed),len(ref),sorted(set(listed)^ref))
    counts=Counter(str(realized[x]) for x in untested)
    assert dict(counts)==pub['untested_by_codimension']
    return len(sampled),len(untested),counts

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--packet',type=Path,default=PACKET)
    args=parser.parse_args()
    directory=args.packet/BASE
    check=subprocess_check(args.packet)
    print('packet SHA256SUMS entries verified',check,flush=True)
    py, pyrec=native('py'); wl, wlrec=native('wl')
    print('native hashes, bytes, object and inventory counts',[(e,r[1],r[2]) for e,r in READ.items()],len(py),len(wl),flush=True)
    assert set(py)==set(wl)
    for k in py: assert py[k][0]==wl[k][0],k
    coverage=json.loads((directory/'directives/_measurements/S15_reading_comparator_coverage.json').read_bytes())
    index={r['object']:r for r in coverage['objects']}
    assert set(index)=={k for k,v in py.items() if v[0]}
    counts=Counter(); aliases=Counter(); status=Counter(); group=Counter(); codims=Counter()
    reclassified=set(); examples={}
    with tarfile.open(ARCH) as t:
        cfg=member(t,'configuration.json')
        print('native input file hash matches for originals')
        for key,d in cfg['inputs'].items():
            if key.endswith('original'):
                assert Path(d['path']).stat().st_size==d['bytes'] and digest(Path(d['path']))==d['sha256']
                print(key,d['bytes'],d['sha256'])
        assert cfg['build_id']==coverage['build_id'] and digest(ARCH)==coverage['archive_sha256']
        for path in t.getnames():
            if not path.startswith('objects/'): continue
            d=member(t,path); key=d['object']; nativep=py[key][0]
            c=counts if nativep else aliases
            ids, points=observation_set(t,path,d,c)
            assert ids==[f'py_reading:{py[key][1]}', f'wl_reading:{wl[key][1]}']
            if not nativep: continue
            assert all(s in ('within_estimate','within_precision_estimate') for p,s in points.values())
            pub=index[key]; assert pub['archive_index']==d['index'] and ids==pub['operand_ids']
            planes=tuple(sorted(nativep))
            feasible=arrangement(planes)
            m,n,cs=compare_sets(pub,planes,feasible,points)
            codims.update(cs)
            group[(len(planes),len(feasible),m,n)]+=1
            status.update(s for p,s in points.values())
            if len(planes)==8 and any(key.startswith(x) for x in ['region/']): examples.setdefault('compound',(key,pub,planes,feasible,points))
            if len(planes)==14: examples.setdefault('mixed',(key,pub,planes,feasible,points))
        slots_arch=member(t,'context_checks.json.gz')['argument_slots']
        assert len(slots_arch)==len(coverage['argument_slots'])
        for orig,new in zip(slots_arch,coverage['argument_slots']):
            if orig!=new: reclassified.add(tuple(orig['argument']))
            else: assert orig['status']!='no_wl_declared_slot'
        assert len(reclassified)==3
    total={'objects_with_planes':sum(group.values()),'objects_with_gaps':sum(v for (*_,n),v in group.items() if n),'sampled':sum(m*v for (_,_,m,_),v in group.items()),'untested':sum(n*v for (_,_,_,n),v in group.items())}
    assert total==coverage['counts']
    print('object groups (native planes, feasible selections, with samples, without samples) : object counts',sorted(group.items()))
    print('all objects with planes archive observations',sorted(counts.items()))
    print('without planes alias count',aliases['alias'],'with planes alias count',counts['alias'])
    print('effective successful observations by status',sorted(status.items()))
    print('missing selections by exact codimension',sorted(codims.items()),'totals',total)
    source=(directory/'mathematica/S15_separation_departure.wl').read_text()
    responses=slots(pyrec,wlrec,coverage['argument_slots'],source)
    print('native velocity argument evidence',responses)
    assert reclassified=={(x[0],0,'t') for x in responses}
    argcounts=dict(Counter(r['status'] for r in coverage['argument_slots']))
    assert argcounts==coverage['argument_counts']
    print('final item 6 status counts',argcounts)
    return examples,source,pyrec,wlrec,coverage

def subprocess_check(directory):
    count=0
    for line in (directory/'SHA256SUMS').read_text().splitlines():
        checksum,relative=line.split('  ')
        assert digest(directory/relative)==checksum
        count+=1
    return count

if __name__=='__main__': main()
