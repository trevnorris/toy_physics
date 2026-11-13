"""Native object checks with physical fields defined independently in reference.py.
Use the reviewed translators on their own byte streams; the reference integrals
and this evaluator do not use the builders' test polynomials, tolerances, or jets.
"""
import sys
import hashlib
import gzip
import tarfile
import json
from collections import Counter
from functools import lru_cache
import sympy as s
import numpy as np
from native import PACKET, records
sys.path.insert(0,str(PACKET/'scripts'))
import S15_reading_comparator as c
import S15_reading_math as m

def tag(e,object_name,field):
    if e=='wl':
        import re
        object_name=re.sub(r'labels_([^/]+)',lambda q:'labels_{'+q[1].replace('+',', ')+'}',object_name)
    return f'reading/{object_name}/{field}'
def request(obj,fields):
    for engine in ('py','wl'):
        need[engine+'_reading'].update(tag(engine,obj,field) for field in fields)
need={e+'_'+stage:set() for e in ('py','wl') for stage in ('reading','original')}
pieces=('carried','pressure','viscous','other','transport_centre','transport_radius')
for stream in ('py_original','wl_original'):
    need[stream].update(f'{loc}/{piece}/component_1' for piece in pieces for loc in ('inner','outer'))
fields=['Q','E_geo','E_resp','separation_derivative','departure_partition/operand_a','departure_partition/operand_b','departure_partition/residual','item3/Q_closed_operand','item3/object','velocity_part_1/operand_a','velocity_part_1/operand_b','velocity_part_1/residual','velocity_part_2/operand_a','velocity_part_2/operand_b','velocity_part_2/residual','E_geo/dimension_reference_operand','E_geo/dimension_pair_residuals','E_geo/dimension_reference_residuals']
for loc in ('inner','outer'):
    request(loc+'/transport_centre/component_1/labels_drain2',fields)
    request(loc+'/other/component_1',['Q'])
request('inner/viscous/component_1/labels_remainder',['E_resp'])
main='both_surfaces/sum_crossing/component_1'
request(main,['Q','E_resp'])
raw,hs=records(need)
def native(engine,obj,field):return raw[engine+'_reading',tag(engine,obj,field)]
def pair(obj,field):
    return [m.parse(native(engine,obj,field)[1],engine) for engine in ('py','wl')]
# Spherical angular product rule, Gauss in cos(theta), periodic trapezoid in phi.
# Native measures already contain sin(theta), so the Jacobian is divided out here.
mu,weights=np.polynomial.legendre.leggauss(18)
theta=np.arccos(mu[:,None]); phi=np.arange(40)[None,:]*(2*np.pi/40)
coords={'theta':theta,'phi':phi,'t':1.,'a_2':1.,'b_2':2.,'a_1':2.,'b_1':1.,
        'N_1':2/3,'N_2':1/3,'N_3':2/3,'d':5.}
vec=[2,-3,5]
matrix=[[2,3,5],[7,11,13],[17,19,23]]
def headval(name,orders,args):
    if name.startswith('X_2_'):
        j=int(name[-1])-1
        return [1,4,-2][j] if not sum(orders) else vec[j]
    if name in ('r_in','R'):return 2 if name=='r_in' else 3
    if name.startswith('Pi_o_1_'):
        j=int(name[-1])-1
        return [1,4,-2][j]+sum(a*b for a,b in zip(matrix[j],args[:3]))
    spec={'rho_2':lambda r,t:5+r**2, 'u_2':lambda r,t:3*r,
          'rho_2_change':lambda r,t:7-r, 'u_2_change':lambda r,t:r+4,
          'A_2':lambda t:3.,'B_2':lambda t:2.}
    return spec[name](*args)
@lru_cache(None)
def independent_value(v):
    ev=independent_value
    if v==s.pi:return np.pi
    if v.is_Number:return float(v)
    if v.is_Symbol:return coords[str(v)]
    if v.is_Add:return sum(ev(q) for q in v.args)
    if v.is_Mul:
        ans=1.
        for a in v.args:ans=ans*ev(a)
        return ans
    if v.is_Pow:return ev(v.base)**ev(v.exp)
    if v.func in (s.sin,s.cos):return (np.sin if v.func==s.sin else np.cos)(ev(v.args[0]))
    if isinstance(v,s.Integral):
        assert set(map(str,v.variables))=={'theta','phi'}
        return float(np.sum(weights[:,None]*np.broadcast_to(ev(v.function),(18,40))/np.sin(theta))*(2*np.pi/40))
    if isinstance(v,m.AppliedUndef):return headval(v.func.__name__,[0]*len(v.args),tuple(ev(q) for q in v.args))
    h,orders,points=m.original.partial_signature(v)
    return headval(h.__name__,orders,tuple(ev(q) for q in points))
for loc in ('inner','outer'):
    obj=loc+'/transport_centre/component_1/labels_drain2'
    print(obj,'(values divided by pi; PY,WL;', 'offsets are native payload start)')
    for f in fields:
        aa,bb=[native(e,obj,f) for e in ('py','wl')]
        if '/dimension_' in f:
            print(f,(aa[0],bb[0]),aa[1],bb[1])
            continue
        norms=[independent_value(expr)/np.pi for expr in pair(obj,f)]
        print(f,(aa[0],bb[0]),*(format(x,'.11g') for x in norms),'diff/pi',format(norms[0]-norms[1],'.3g'))

def check_localized():
    obj='inner/other/component_1'
    py,wl=pair(obj,'Q')
    text=native('wl',obj,'Q')[1]
    for label,changed in [('original',text),('one-sided sign','-('+text+')'),('one-sided offdiagonal',text.replace('oo[1, 2]','oo[1, 3]'))]:
        result=m.parse(changed,'wl')
        lhs,rhs=[independent_value(v)/np.pi for v in (py,result)]
        print(label,obj,'Q/pi',format(lhs,'.12g'),format(rhs,'.12g'),'residual/pi',format(lhs-rhs,'.12g'),'exact_zero',m.exact_residual(py,result)==0)
    obj='inner/viscous/component_1/labels_remainder'
    py,wl=pair(obj,'E_resp'); key='Derivative[0, 1, 0, 0][change[ww[1]]]'
    txt=native('wl',obj,'E_resp')[1]; assert txt.count(key)==1
    bad=m.parse(txt.replace(key,'Derivative[0, 0, 1, 0][change[ww[1]]]'),'wl')
    ctx=m.Sample('mixed_quadratic','generic',0,(1,2,3,4),24)
    print('native derivative slot collapse y->z only on WL, component 1 response', [native(e,obj,'E_resp')[0] for e in ('py','wl')])
    print('unchanged PY, WL baseline, WL changed',*[float(ctx.value(v)) for v in (py,wl,bad)])
    for kind,e in [('good',wl),('bad',bad)]:print(kind,'difference',float(ctx.value(py)-ctx.value(e)),'exact_zero',m.exact_residual(py,e)==0)
check_localized()
print('signed original terms via exact difference with native reading Q')
for e in ('py','wl'):
    parts=[m.parse(raw[e+'_original',f'{loc}/{piece}/component_1'][1],e) for loc in ('inner','outer') for piece in pieces]
    q=m.parse(native(e,main,'Q')[1],e)
    print(e,'original 12 piece residual ==0',m.exact_residual(q,sum(parts))==0)
    if e=='wl':
        wrong=sum(parts)-2*parts[3] # inner other piece only; rest retains native values
        rs=m.exact_residual(q,wrong)
        print('WL original reference inner other only sign reversed, exact difference zero?',rs==0,
              'equals twice native inner other?',m.exact_residual(rs,2*parts[3])==0)

def outer_alias():
    py,wl=pair(main,'E_resp')
    text=native('wl',main,'E_resp')[1]
    before='change[oo[1, 2]]'; after='change[oo[1, 3]]'
    assert text.count(before)==2 and 'rin[t]' in text[:text.index(before)]
    bad=m.parse(text.replace(before,after,1),'wl')
    with tarfile.open(PACKET/'scripts/out/S15_reading_comparator_details.out') as tar:
        archive=json.load(gzip.GzipFile(fileobj=tar.extractfile('objects/000.json.gz')))
    baseline=next(r for r in archive['rows'] if r['field']=='E_resp')
    result=m.sample_values(py,bad,())
    print('native net E_resp only inner change[oo[1,2]] -> change[oo[1,3]];',
          'native offsets',*[native(e,main,'E_resp')[0] for e in ('py','wl')])
    print('raw cross status',baseline['status'],baseline['sample_counts'],'->',result['status'],result['sample_counts'])
    orig=baseline['samples'][0]; mut=result['samples'][0]
    print('case',mut['case'],'old',orig['values'][-1],orig['residual'],'changed',mut['values'][-1],mut['residual'],'tolerance',mut['tolerance'])
    print('unchanged PY value at all 30 points',all(a['values'][-1][0]==b['values'][-1][0] for a,b in zip(baseline['samples'],result['samples'])))
    print('number of changed rounded WL 24-node values',sum(a['values'][-1][1]!=b['values'][-1][1] for a,b in zip(baseline['samples'],result['samples'])))
    # Difference of just the two WL payloads, preserving the net object's
    # actual algebra and isolating the perturbed residue before roundoff.
    isolated=m.exact_residual(wl,bad)
    ctx=m.Sample('mixed_quadratic','generic',0,(1,2,3,4),24)
    print('isolated original minus mutated native WL at same point',ctx.value(isolated))
outer_alias()
