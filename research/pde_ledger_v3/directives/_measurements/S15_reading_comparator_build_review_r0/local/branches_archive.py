"""Compare actual family outputs with a meridional integral from the surface flux.
Audit every published raw roster, dimension table, and first-stage separation tally.
"""
from native import records,PACKET
import ast,gzip,json,sys,tarfile,collections
import mpmath as mp
import sympy as s
sys.path.insert(0,str(PACKET/'scripts'))
import S15_reading_math as m
from S15_reading_refinement import rational
import S15_reading_comparator as c
wlo='reading/inner/transport_radius/component_1/labels_{drain1}/'
pyo='reading/inner/transport_radius/component_1/labels_drain1/'
fields=('item3/Q_closed_operand','item3/object','item3/exponent_branch_inventory')
need={e+'_reading':{pre+f for f in fields} for e,pre in [('py',pyo),('wl',wlo)]}
raw,hashes=records(need)
def val(e,f):return raw[e+'_reading',(pyo if e=='py' else wlo)+f]
print('domain d=5, r=2, dr/dt=3/2, N=(1,0,0); A1=2 and B1=3')
print('branch raw',val('py',fields[2]),val('wl',fields[2]))
# Independent reduction under rotation along Cartesian x_1 makes e_1=mu.
# Q = -2 pi r^2 rdot A B int (r mu-d)*(r^2+d^2-2drmu)^((a+b-1)/2) dmu.
def evaluate(expr,a,b):
    t=s.Symbol('t')
    sub={s.Symbol('d'):5,**dict(zip(m.DIRECTION,(1,0,0))),m.EXP[0]:a,m.EXP[2]:b,
         s.Derivative(s.Function('r_in')(t),t):s.Rational(3,2),s.Function('r_in')(t):2,
         s.Function('A_1')(t):2,s.Function('B_1')(t):3}
    ans=expr.xreplace(sub)
    return mp.mpf(str(s.N(ans,58)))
mp.mp.dps=55
def independent(d,a,b):
    return -2*mp.pi*4*mp.mpf('1.5')*2*3*mp.quad(lambda mu:(2*mu-d)*(4+d*d-4*d*mu)**(mp.mpf(a+b-1)/2),[-1,1])
closed={}
for f in fields[:2]:
    closed[f]=[m.parse(val(e,f)[1],e) for e in ('py','wl')]
wtext=val('wl',fields[0])[1]
token='(Pi*Log[dist - rin[t]])'
assert token in wtext
bad=m.parse(wtext.replace(token,'(-Pi*Log[dist - rin[t]])',1),'wl')
for a,b in [(2,1),(-5,0),(-3,0),(-1,0),(-7,0)]:
    d=mp.mpf(5)
    q=independent(d,a,b); effect=d*mp.diff(lambda z:independent(z,a,b),d)+2*q
    print('sum',a+b)
    for f,ref in zip(fields[:2],(q,effect)):
        py,wl=[evaluate(q,a,b) for q in closed[f]]
        print(f,'analytic 1d integral PY WL, max native deviation',*[mp.nstr(n,19) for n in (ref,py,wl,max(abs(py-ref),abs(wl-ref)))])
    py,wl=[evaluate(q,a,b) for q in closed[fields[0]]]
    changed=evaluate(bad,a,b)
    print('one native log sign in WL closed Q: changed WL, residual PY-WL',mp.nstr(changed,19),mp.nstr(py-changed,19))
print('original and changed raw span',val('wl',fields[0])[0],repr(token),'->','(-Pi*Log[dist - rin[t]])')
archive=PACKET/'scripts/out/S15_reading_comparator_details.out'
with tarfile.open(archive) as tf:
    def load(fname):return json.load(gzip.GzipFile(fileobj=tf.extractfile(fname))) if fname.endswith('.gz') else json.load(tf.extractfile(fname))
    index=load('configuration.json')
    native={r['id']:r for r in load('coverage.json.gz')}
    names=index['required_objects']; print('unique required objects',len(set(names)))
    same=collections.defaultdict(collections.Counter); dimensions=collections.Counter()
    physical=collections.defaultdict(collections.Counter)
    representative=[]
    offsets={}
    for idx,name in enumerate(names):
        rows=load(f'objects/{idx:03d}.json.gz')['rows']
        for r in rows:
            if r['role'].endswith('_same_engine'):
                if len(r['right'])==1:
                    eq=native[r['left'][0]]['payload']==native[r['right'][0]]['payload']
                    same[r['role']]['direct same payload' if eq else 'changed payload']+=1
                else:
                    same[r['role']]['surface sum '+r['status']]+=1
            elif r['role']=='cross' and r['item']==6 and not r['field'].endswith('comparison'):
                dimensions[r['status']]+=1
            elif r['role']=='cross' and r['field'] in ('Q','E_geo','E_resp','item3/object','velocity_part_1/residual','velocity_part_2/residual'):
                physical[r['field']][r['status']]+=1
            if name=='region/component_1/rate_labels_drain2' and r['field']=='Q':
                representative.append((name,r['status'],r['samples'][0]['status'],r['samples'][0]['residual']))
    print('source/reference links',dict(same))
    print('dimension rows',dict(dimensions))
    print('raw object field counts',dict(physical))
    print('archive example unresolved convergence',representative)
    refine,extra,rechecked=collections.Counter(),collections.Counter(),{}
    for member in tf.getnames():
        if member.startswith('log_branches/') and member.endswith('.gz'):
            data=load(member)
            for r in data['rechecks']: rechecked[(data['object'],r['field'],tuple(r['case']))]=r['status']
    for member in tf.getnames():
        if member.startswith('refinement/') and member.endswith('.gz'):
            data=load(member)
            for f in data['fields']:
                for r in f['rechecks']: refine[r['status']]+=1
                for r in f['additional_samples']:
                    key=(data['object'],f['field'],tuple(r['case']))
                    extra[rechecked.get(key,r['status'])]+=1
    print('independent recheck counts',dict(refine),'extra cases after log replacements',dict(extra))
