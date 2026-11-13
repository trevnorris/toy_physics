"""One-sided component substitution on the actual total shell response, using the frozen reviewer packet."""
import gzip
import json
from pathlib import Path
import sys
import tarfile
from collections import Counter
from itertools import product

base=Path(__file__).resolve().parents[1]
repo=base.parents[3]
scripts=base/'packet/research/pde_ledger_v3/scripts'
sys.path.insert(0,str(scripts))
import S15_reading_math as math
import sympy as s

with tarfile.open(scripts/'out/S15_reading_comparator_details.out') as tar:
    def read(name):
        return json.loads(gzip.decompress(tar.extractfile(name).read()))
    c=json.load(tar.extractfile('configuration.json'))
    loc={r['id']: r for r in read('coverage.json.gz')}
    obj=read('objects/000.json.gz')
    row=next(r for r in obj['rows'] if r['field']=='E_resp')
    raw=[]
    for rid in (row['left']+row['right']):
        stream=rid.split(':')[0]
        span=loc[rid]['span']
        with (repo/c['inputs'][stream]['path']).open('rb') as f:
            f.seek(span['offset']);raw.append(f.read(span['length']).decode().strip())
assert raw[1].count('change[oo[1, 1]]')==2
changed=raw[1].replace('change[oo[1, 1]]','change[oo[1, 2]]')
left,right,mutant=math.parse(raw[0],'py'),math.parse(raw[1],'wl'),math.parse(changed,'wl')
print('object',obj['object'],row['field'],'native source occurrence IDs',row['left'],row['right'],flush=True)
print('Two Wolfram occurrences of change[oo[1,1]] replaced with change[oo[1,2]], preserving the other fields and measures.',flush=True)
# Evaluate the change in isolation, using the native integral and derivative of the interpolated input.
key=s.Symbol('corruption')
subs={x:x+key*(s.Function('Pi_o_1_2_change')(*x.args)-x) for x in right.atoms(math.AppliedUndef) if x.func.__name__=='Pi_o_1_1_change'}
assert len(subs)==2
def plane(v):
    if not v.has(key):return s.S.Zero
    if isinstance(v,s.Integral):return s.Integral(plane(v.function),*v.limits)
    if v.is_Add:return s.Add(*(plane(q) for q in v.args))
    return s.diff(v,key).subs(key,0)
diff=plane(right.xreplace(subs))
assert diff!=0
out=math.sample_values(left,mutant,())
print('Frozen published status',row['status'],row['sample_counts'],'one-sided mutation',out['status'],out['sample_counts'],flush=True)
direct=[]
for case in row['samples']:
    form,stratum,point,_=case['case']
    exps=tuple(map(s.Rational,case['exponents']))
    direct.append(float(math.Sample(form,stratum,point,exps,24).value(diff)))
print('Isolated mutated Wolfram minus original Wolfram at quadrature order 24, first sample:',direct[0],flush=True)
print('Min/max of isolated signed changes across all 30 cases:',min(direct),max(direct),flush=True)
print('Nonzero isolated samples:',sum(v!=0 for v in direct),'of',len(direct),flush=True)
print('First-case full baseline then mutated operands (order 24):',row['samples'][0]['values'][-1],out['samples'][0]['values'][-1],flush=True)
