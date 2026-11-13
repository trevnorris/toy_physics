"""Evaluate both native integrals at an independent smooth fixture with all drain groups live."""
import sys,io,contextlib,re,itertools,functools,tarfile,collections
from streams import scan,raw,PACKET,detail,text_at
from native_controls import c,m,compare,describe
import sympy as s
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import reference as ref

variables=s.symbols('xx yy zz tt')
xx,yy,zz,tt=variables
x=s.Matrix(variables[:3])
C1=s.Matrix([7,2,1])+tt*s.Matrix([1,1,2])+tt**2*s.Matrix([1,-2,3])
C2=ref.X+tt*ref.C+tt**2*s.Matrix([2,1,3])
d1=sum((a-b)**2 for a,b in zip(x,C1))
d2=sum((a-b)**2 for a,b in zip(x,C2))
W=ref.v.subs(dict(zip(ref.x,x)))-(1+tt)*(x-C1)-(2+tt)*(x-C2)
special={'r_in':(1+tt/2+tt**2,(tt,)), 'R':(2-2*tt/3+tt**2,(tt,)),
         'p_br':(xx**2,(xx,)), 'eta':(2+xx,variables),'zeta':(3+yy,variables),
         'varrho':(ref.rho.subs(dict(zip(ref.x,x)))-(1+d1+tt)-(3+2*d2+tt),variables)}
for i in range(3):
    special[f'X_1_{i+1}']=(C1[i],(tt,))
    special[f'X_2_{i+1}']=(C2[i],(tt,))
    special[f'W_{i+1}']=(W[i],variables)
    for j in range(3):special[f'Pi_o_{i+1}_{j+1}']=(ref.other[i,j].subs(dict(zip(ref.x,x))),variables)
for k in (1,2):
    special[f'u_{k}']=((k+tt)*xx,(xx,tt))
    special[f'rho_{k}']=((1 if k==1 else 3)+k*xx**2+tt,(xx,tt))
for k in (1,2):
    for prefix in ('A','B'):
        special[f'{prefix}_{k}']=((1 if prefix=='A' else 3)+k+tt,(tt,))
@functools.lru_cache(None)
def jet(name,arity,orders):
    base,slots=special[name]
    v=s.diff(base,*[(slot,n) for slot,n in zip(slots,orders) if n]) if any(orders) else base
    return s.lambdify(slots,v,'numpy')

class Native(m.Sample):
    def __init__(self):
        self.order=16
        self.envs=[{}];self.grids={};self.cache={}
        sep=np.array([float(c.subs(tt,0)) for c in C1-C2]); d=np.linalg.norm(sep)
        self.constants={'t':0,'d':d,**dict(zip(map(str,m.DIRECTION),sep/d)),**dict(zip(map(str,m.EXP),[1,-2,3,0]))}
    def compute(self,v,environment):
        sig=None
        if isinstance(v,m.AppliedUndef):sig=(v.func,(0,)*len(v.args),v.args)
        elif isinstance(v,(s.Subs,s.Derivative)):sig=m.original.partial_signature(v)
        if sig:
            f,orders,args=sig
            return jet(f.__name__,len(args),orders)(*(self.value(a,environment) for a in args))
        return super().compute(v,environment)
def compare_references(ix):
    ctx=Native()
    for loc in ('inner','outer','both_surfaces'):
        shell=ref.shells[0] if loc=='inner' else ref.shells[1] if loc=='outer' else None
        for piece in ('rest','transport','crossing'):
            expected=shell[piece] if shell else tuple(map(sum,zip(*(f[piece] for f in ref.shells))))
            for i,truth in enumerate(expected,1):
                tag=f'reading/{loc}/sum_{piece}/component_{i}/Q'
                numbers=[ctx.value(m.parse(raw(ix,k+'_reading',tag),k))/(4*np.pi) for k in ('py','wl')]
                print('native sums/(4*pi)',loc,piece,i,'ref',truth,'py,wl',*numbers,'largest error',max(abs(float(truth)-z) for z in numbers),flush=True)
    for k,l in [('py','drain2'),('wl','{drain2}')]:
        pre=f'reading/outer/transport_centre/component_1/labels_{l}/'
        # Spherical second moment int n_i n_j = 4*pi/3 delta_ij. Geometric
        # derivative is zero for this target-only current. Use nonzero target velocity.
        # At a2=-2, b2=0, R=2, A2=3, B2=5, v2_x=2, Q/(4*pi)=10.
        for field,truth in [('item3/Q_closed_operand',10),('item3/family_integral',20),('item3/object',20)]:
            val=ctx.value(m.parse(raw(ix,k+'_reading',pre+field),k))/(4*np.pi)
            print(k,'family',field,'ref',truth,'observed',val,flush=True)
    for k,l in [('py','drain2'),('wl','{drain2}')]:
        pre=f'reading/outer/transport_centre/component_1/labels_{l}/'
        for drain in [1,2]:
            vals=[ctx.value(m.parse(raw(ix,k+'_reading',pre+f'velocity_part_{drain}/'+r),k))/(4*np.pi) for r in ['operand_a','operand_b','residual']]
            print(k,'actual velocity part',drain,'original,cut,residual / (4*pi)',vals,flush=True)
def controls(ix):
    target1='reading/inner/viscous/component_1/labels_remainder/Q';target2=target1.replace('labels_remainder','labels_{remainder}')
    a=raw(ix,'py_reading',target1);b=raw(ix,'wl_reading',target2)
    for label,wrong in [('baseline slot',b),('z slot to y',b.replace('Derivative[0, 0, 1, 0][ww[1]]','Derivative[0, 1, 0, 0][ww[1]]',1))]:
        describe(label,compare(a,wrong,'Q',1))
    # Target a native logarithmic exponent branch; ordinary exponents share the same operand.
    key='reading/outer/transport_centre/component_1/labels_'
    a=raw(ix,'py_reading',key+'drain1/item3/object');b=raw(ix,'wl_reading',key+'{drain1}/item3/object')
    changed=b.replace('Log[dist - rout[t]]','(Log[dist - rout[t]]+1)',1)
    print('branch mutation, WL first log changed, remaining log count',b.count('Log[dist - rout[t]]')-1)
    loci=m.branch_loci(raw(ix,'py_reading',key+'drain1/item3/exponent_branch_inventory'),'py')
    a,b,bad=[m.parse(v,e) for v,e in [(a,'py'),(b,'wl'),(changed,'wl')]]
    for branch,exp in m.exponent_samples(0,loci):
        ctx=m.Sample(m.FORMS[0],m.STRATA[0],0,exp,24)
        x,y,z=[ctx.value(o) for o in (a,b,bad)]
        print('exp',branch,exp,'PY,WL,WLmut',x,y,z,'old,new residual',x-y,x-z,flush=True)
if __name__=='__main__':
    ix=scan()
    compare_references(ix)
    controls(ix)
