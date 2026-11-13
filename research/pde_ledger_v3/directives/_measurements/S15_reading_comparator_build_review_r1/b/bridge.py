"""Evaluate native printouts at the analytic data in native_reference.py.

Independent field fixture uses SymPy function derivatives and chart quadrature.
The independent angular reference uses exact rational monomial moments.
"""
from functools import lru_cache
import re
import numpy as np
import sympy as s
from packet_io import *
import native_reference as nr
import S15_reading_math as m
from S15_reading_refinement import PrecisionSample
rec=Records()
REFERENCE=nr.surface
t0=m.POINTS[0][-1]; t=s.symbols('t'); args=m.original.FORMAL
order0=tuple(next(m.exponent_samples(0,()))[1])
source=(3,-1,2);target=(1,3,-2)
def case_values(st):
    v1=(0,0,0) if st=='X1_velocity_zero' else source
    v2=(0,0,0) if st in ('X2_velocity_zero','X2_and_radii_velocity_zero') else target
    bdot=(0,0) if st in ('radii_velocity_zero','X2_and_radii_velocity_zero') else (-2,4)
    return v1,v2,bdot
@lru_cache(None)
def fixture(name,arity,orders,form,stratum,point):
    v1,v2,rates=case_values(stratum);t0=m.POINTS[point][-1]
    coeff=s.Symbol('lambda_c');a=args[:arity]; dt=a[-1]-t0
    base=name.removesuffix('_change'); c2=(1,-2,3);D=(4,3,0)
    if base.endswith('_velocity'):
        v=5+int(base[4]); q=v+(v+1)*dt
    elif m.geometry_head(base):
        if base in ('r_in','R'):
            k=int(base=='R'); q=(1,2)[k]+rates[k]*dt+(11,-13)[k]*dt**2
        else:
            _,ks,js=base.split('_'); k=int(ks);j=int(js)-1
            q=c2[j]+(D[j] if k==1 else 0)+(v1 if k==1 else v2)[j]*dt+(7+k+3*j)*dt**2
    else:
        if base in ('u_1','u_2','rho_1','rho_2'):
            q={'u_1':(1+3*coeff)*a[0],'u_2':(-2+5*coeff)*a[0],'rho_1':(2+7*coeff)*a[0]**2,'rho_2':(3-11*coeff)*a[0]**2}[base]
        elif base=='varrho':q=8+13*coeff+sum(b*x for b,x in zip((2+coeff,3-2*coeff,1+4*coeff),a[:3]))
        elif base=='p_br':q=a[0]**2+a[0]+1
        elif base=='eta':q=9+a[0]+coeff*(8+a[2])
        elif base=='zeta':q=2+a[1]+coeff*(4-a[0])
        elif base.startswith('Pi_o_'):
            i,j=[int(x)-1 for x in base.split('_')[-2:]]
            q=sum(((i+2)**2*(j+1)*(k+1)+coeff*(2+i+3*j)**2*(k+1)**2)*a[k] for k in range(3))+(1+i)*j+coeff*(2+j-i)
        elif base.startswith('W_'):
            i=int(base[-1])-1;q=i+1+(i+5)*coeff+sum(((-1)**(i+j)*(1+3*i+j)+(2+i+2*j)*coeff)*a[j] for j in range(3))
        else:raise ValueError(name)
        q=q.expand().coeff(coeff,int(name.endswith('_change')))
        if arity>1:q+=(2+len(name))*dt + (5+len(name))*dt**2
    q=s.diff(q, *((k,n) for k,n in zip(a,orders) if n)) if any(orders) else q
    return s.lambdify(a,q,modules='numpy')
def exact(num):return float(num)*4*np.pi
raw_jet=m.jet_function
m.jet_function=fixture
@lru_cache(None)
def expr(name,field,e):
    ids=rec.ref(name,field)
    return m.parse(rec.raw(ids[int(e=='wl')]),e)
def val(a,st,order=40):
    return float(m.Sample('mixed_quadratic',st,0,order0,order).value(a))
def close(x,y):
    return abs(x-y)<max(2e-6,5e-12*max(abs(x),abs(y)))
number=0
for st in m.STRATA:
    v1,v2,rates=case_values(st)
    sides=[(1,-1,rates[0]),(2,1,rates[1])]
    local={mode:[REFERENCE(b,n,r,v2,*opts) for b,n,r in sides] for mode,opts in [('geo',(1,0)),('resp',(0,1))]}
    case_names=['both_surfaces'] if st!='generic' else ['inner','outer','both_surfaces']
    worst=0.
    for side in case_names:
        take={'inner':[0],'outer':[1],'both_surfaces':[0,1]}[side]
        for typ in ('rest','transport','crossing'):
            for i in range(3):
                nm=f'{side}/sum_{typ}/component_{i+1}'
                q=sum(local['geo'][k][typ][i][0] for k in take)
                g=sum(local['geo'][k][typ][i][1] for k in take)
                h=sum(local['resp'][k][typ][i][1] for k in take)
                for fld,ref in [('Q',q),('E_geo',g+2*q),('E_resp',h),('separation_derivative',h+g)]:
                    vals=[val(expr(nm,fld,e),st) for e in ('py','wl')]
                    want=exact(ref); res=[x-want for x in vals]
                    worst=max(worst,max(map(abs,res))); number+=2
                    if nm=='both_surfaces/sum_crossing/component_2' and fld in ('E_resp','Q'):
                        print('fixture',st,nm,fld,rec.ref(nm,fld),'normalized by 4pi',ref,'native values',vals,'native-minus-reference',res,flush=True)
                    assert all(close(want,x) for x in vals),(nm,fld,ref,vals,res)
    print(st,'cumulative operand count',number,'max abs error',worst,flush=True)
print('one-sided actual native changes')
nm='both_surfaces/sum_crossing/component_2';field='E_resp'
ids=rec.ref(nm,field); a=expr(nm,field,'py'); w=rec.raw(ids[1]);print('WL form snippet',re.search(r'.{30}change\[oo\[2, 1\]\].{90}',w)[0],flush=True)
for label,raw in [('baseline',w),('sign', '-('+w+')'),('form',w.replace('change[oo[2, 1]]','change[oo[1, 2]]'))]:
    b=m.parse(raw,'wl'); diffs=[]
    for st in m.STRATA:
        x=val(a,st); y=val(b,st);diffs.append(x-y)
        if st=='generic': print('native',label,ids, 'number of changed tensor calls',w.count('change[oo[2, 1]]') if label=='form' else 0,'PY',x,'WL',y,'difference',x-y)
    print(label,'all five differences',diffs);assert all((abs(q)<1e-5 if label=='baseline' else abs(q)>1) for q in diffs)
m.jet_function=raw_jet
# Retrieve as-stored numerical operands at a non-generic stationary branch.
nm='both_surfaces/sum_crossing/component_2';field='E_resp';obj=rec.obj(rec.cfg['required_objects'].index(nm));row=rec.row(nm,field)
if 'same_comparison_as' in row:row=obj['rows'][row['same_comparison_as']]
case=next(x for x in row['samples'] if x['case']==['mixed_cubic','X2_and_radii_velocity_zero',1,'generic'])
print('raw sample',case)
f,st,p,_=case['case'];ctx=m.Sample(f,st,p,tuple(map(s.Rational,case['exponents'])),32)
vals=[float(ctx.value(expr(nm,field,e))) for e in ('py','wl')];print('recalculated',vals,vals[0]-vals[1]);assert all(close(x,y) for x,y in zip(vals,case['values'][-1]))
