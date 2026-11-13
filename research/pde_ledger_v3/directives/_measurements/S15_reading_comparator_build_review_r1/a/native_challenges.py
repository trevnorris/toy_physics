"""Ground accepted native samples against independent source definitions and mutate one native input at a time."""
import gzip, hashlib, itertools, json, math, re, sys, tarfile
from pathlib import Path
from unittest.mock import patch
from functools import lru_cache
from collections import Counter, defaultdict
import numpy as np
import sympy as s
import mpmath as mp
from scipy.special import roots_legendre

from archive_audit import P,R,L,A,member
import native_reference as n
sys.path.insert(0,str(L/'scripts'))
import S15_reading_math as m
import S15_reading_refinement as r
import S15_reading_comparator as c
import S15_reading_context_checks as k

sys.setrecursionlimit(14000)
tf=tarfile.open(A)
cfg=member(tf,'configuration.json')
lookup={v['id']:v for v in member(tf,'coverage.json.gz')}
def read(rid):
    v=lookup[rid]; off,size=v['span']['offset'],v['span']['length']
    with (R/cfg['inputs'][rid.split(':')[0]]['path']).open('rb') as f:
        f.seek(off);payload=f.read(size-1)
    assert hashlib.sha256(payload).hexdigest()==v['payload']
    return payload.decode()
def obj(index): return member(tf,f'objects/{index:03}.json.gz')
def findrow(o,field): return next(v for v in o['rows'] if v['role']=='cross' and v['field']==field)
def getobj(name): return obj(cfg['required_objects'].index(name))
def operands(o,field):
    row=findrow(o,field)
    refs=row['left'][0],row['right'][0]
    return refs,[read(q) for q in refs]
def parsed(pair):return [m.parse(v,eng) for v,eng in zip(pair,('py','wl'))]
def ctx(form=m.FORMS[0],stratum='generic',point=0,order=32,exp=None):
    return m.Sample(form,stratum,point,exp if exp is not None else next(m.exponent_samples(point,()))[1],order)
def change_one(value,old,new,count=1):
    assert value.count(old)>=count
    return value.replace(old,new,count)
def archivecases():
    o=getobj('both_surfaces/sum_crossing/component_2')
    refs,v=operands(o,'E_resp')
    print('raw net response',refs,*(len(q) for q in v))
    py,wl=parsed(v)
    # Two different components and signs; with every term and all native bounds present.
    mutations={'wrong inner orientation':change_one(v[1],'Cos[phi]*rin[t]^2*Sin[theta]^2*change[oo[2, 1]]','-Cos[phi]*rin[t]^2*Sin[theta]^2*change[oo[2, 1]]'),
               'off diagonal same row':change_one(v[1],'change[oo[2, 1]]','change[oo[2, 3]]'),
               'transpose off diagonal':change_one(v[1],'change[oo[2, 1]]','change[oo[1, 2]]')}
    mutants={k:m.parse(y,'wl') for k,y in mutations.items()}
    base=findrow(o,'E_resp')
    if 'same_comparison_as' in base: base=o['rows'][base['same_comparison_as']]
    first=base['samples'][0]
    print('archive actual accepted raw sample',first)
    detections=Counter()
    for f,st in itertools.product(m.FORMS,m.STRATA):
        C=ctx(f,st)
        x,y=C.value(py),C.value(wl)
        archived=next(a for a in base['samples'] if a['case']==[f,st,0,'generic'])
        assert [float(x),float(y)]==archived['values'][-1]
        for name,val in mutants.items():
            z=C.value(val)
            if abs(x-z)>archived['tolerance']:detections[name]+=1
            if f==m.FORMS[0] and st==m.STRATA[0]:
                print('one-sided net native ablation',name,float(x),float(z),float(x-z),'baseline empirical allowance',archived['tolerance'])
    print('one-sided detections across all ten form/stratum pairs',dict(detections))
    assert set(detections.values())=={10}

def reference_shells():
    u=m.original.FORMAL; t=u[0]
    x,y,z,time=u
    # Analytic jets of the reference realization at a common instant; retain original bounds and both measures.
    tr=t-m.POINTS[0][-1]
    values={'r_in':n.inner+tr*n.radial_rates[0]+tr**2,'R':n.outer+tr*n.radial_rates[1]+2*tr**2,
            'eta':2+y,'zeta':3+x,'varrho':5+x-2*z,
            'u_1':2*u[0]/3,'u_2':-u[0]/5,'rho_1':2+u[0]**2/7,'rho_2':3+u[0]**2/11,'p_br':u[0]**2}
    for j in range(3):
        values[f'X_1_{j+1}']=n.source[j]+tr*n.v1[j]+(j+1)*tr**2
        values[f'X_2_{j+1}']=n.target[j]+tr*n.v2[j]+(j+4)*tr**2
    w=s.Matrix([[1,2,-1],[3,1,2],[-2,4,1]])*s.Matrix([x,y,z])+s.Matrix([2,3,5])
    for i in range(3):
        values[f'W_{i+1}']=w[i]
        for j in range(3):
            code=3*i+j+1;values[f'Pi_o_{i+1}_{j+1}']=code+sum(code**(d+1)*v for d,v in enumerate((x,y,z)))
    old=m.jet_function
    @lru_cache(None)
    def jet(name,arity,orders,form,stratum,point):
        if name not in values:return old(name,arity,orders,form,stratum,point)
        ds=[(v,k) for v,k in zip(u[:arity],orders) if k]
        value=s.diff(values[name],*ds) if ds else values[name]
        return s.lambdify(u[:arity],value,'numpy')
    # The exact rational values were derived, emitted and frozen by native_reference.py first.
    exact={ 'inner':[[s.Rational(v) for v in q] for q in [('1420487/15400','33333901/207900','9235679/69300'),('6538/825','216433/9900','-119983/1650')]],
            'outer':[[s.Rational(v) for v in q] for q in [('-278508387/107800','-235871043/53900','-182140611/53900'),('2774547/5390','4339809/5390','886122/539')]]}
    exact['both_surfaces']=[[a+b for a,b in zip(u,v)] for u,v in zip(exact['inner'],exact['outer'])]
    maxerr=0.
    with patch.object(m,'jet_function',jet):
        C=ctx()
        for side in exact:
            a,b=exact[side]
            for i,field in itertools.product(range(3),('rest','transport','crossing')):
                name=f'{side}/sum_{field}/component_{i+1}'
                _,raw=operands(getobj(name),'Q'); values_pair=[float(C.value(v))/math.pi for v in parsed(raw)]
                ref={'rest':a[i],'transport':b[i],'crossing':a[i]+b[i]}[field]
                err=max(abs(float(ref)-v) for v in values_pair)
                maxerr=max(maxerr,err)
                if side=='both_surfaces' and field=='crossing':print('native Q/pi',name,ref,values_pair)
                assert err<1.e-9
    print('54 independent signed native sums max abs error in units of pi',maxerr)

def jets():
    for f in m.FORMS:
        for st in m.STRATA:
            checked=0
            for p in range(3):
                origin=m.POINTS[p][-1]
                for name in [f'X_{k}_{j}' for k in (1,2) for j in (1,2,3)]+['r_in','R']:
                    stopped=(name.startswith('X_1_') and st=='X1_velocity_zero' or
                        name.startswith('X_2_') and st in ('X2_velocity_zero','X2_and_radii_velocity_zero') or
                        name in ('R','r_in') and st in ('radii_velocity_zero','X2_and_radii_velocity_zero'))
                    actual=[m.jet_function(name,1,(k,),f,st,p)(float(origin)) for k in (1,2,3)]
                    assert (actual[0]==0)==stopped and all(v!=0 for v in actual[1:])
                    checked+=1
                v=[m.jet_function(f'X_1_{j}_velocity_change',1,(0,),f,st,p)(float(origin)) for j in (1,2,3)]
                assert all(q!=0 for q in v) and len(set(v))==3
            print('geometry cases exact stationary slopes, nonzero orders 2,3 and live velocity responses',f,st,checked)
    for form in m.FORMS:
        vals=[m.jet_function(f'Pi_o_{i}_{j}_change',4,(0,0,0,0),form,'generic',0)(.11,.28,-.57,1/6) for i in (1,2,3) for j in (1,2,3)]
        print('nonsymmetric native tensor response samples',form,len(set(vals)),vals)
        assert len(set(vals))==9
def family():
    rd=defaultdict(lambda:defaultdict(list))
    for q in lookup.values():
        st=q['id'].split(':')[0]
        if st in ('py_reading','wl_reading') and 'payload' in q:rd[st[:2]][q['tag']].append(q['id'])
    for eng,tag in [('py','family/centre_substitution'),('wl','profile_family/separation_chart'),('wl','profile_family/domain')]:
        rid,=rd[eng][tag];print('native binding',rid,read(rid))
    # The frozen candidate locates Git blobs relative to its own checkout. Use the actual review repo as root.
    c.inputs.ROOT=R;c.inputs.LEDGER=R/'research/pde_ledger_v3'
    for eng,tag,old,new in [('wl','profile_family/separation_chart','xx[1, 2][t] ->','xx[2, 2][t] ->'),
                          ('py','family/centre_substitution',' + X_1_2(t),',' + X_2_2(t),'),
                          ('wl','profile_family/domain','Reals','Complexes')]:
        rid,=rd[eng][tag]
        changes=[row['field'] for row in c.family_context(rd['py'],rd['wl'],lambda x:change_one(read(x),old,new) if x==rid else read(x)) if row['status']!='same_values']
        print('one-sided context mutation',eng,tag,old,'=>',new,changes)
        assert changes
    rows=member(tf,'context_checks.json.gz')['kernels']
    for q in rows:
        if q['field'].startswith(('positive endpoint','power primitive,')) or q['field']=='angular moment (0, 2, 0)':print('native kernel',q)
    print('full azimuth measured rows',sum(q.get('py_bounds')==q.get('wl_bounds')=='(kernel_az, 0, 2*pi)' for q in rows))

if __name__=='__main__':
    jets()
    family()
    reference_shells()
    archivecases()
