"""Exact family directional derivative at fixed target from the P5/P6 surface tensors."""
from reference_bridge import *
from sympy import Matrix,eye,Rational as R
P=s.Symbol('familyparameter')
# Smooth radial velocity responses proportional to radius keep the total fields polynomial in Cartesian x.
for k in (1,2):
    special[f'u_{k}_change']=((1+2*k+tt)*xx,(xx,tt))
    special[f'rho_{k}_change']=(2+k+tt+k*xx**2,(xx,tt))
for k in range(3):
    special[f'X_1_{k+1}_velocity_change']=(3+k+tt,(tt,))
    special[f'W_{k+1}_change']=((k+1)*xx-(k+2)*yy+3*zz+tt,variables)
    for j in range(3):
        special[f'Pi_o_{k+1}_{j+1}_change']=((k+1)*xx+(j+1)*yy+(2*k+j+1)*zz+tt,variables)
for k,expr in [('varrho_change',xx-2*yy+zz+tt),('eta_change',2+xx+yy),('zeta_change',1-zz+tt)]:special[k]=(expr,variables)
def family_variation(response):
    # Work at epsilon=t=0 on both shells. The positional part moves X1 by D.
    q=lambda h:special[h][0].subs(tt,0)
    c1,c2=(cc.subs(tt,0) for cc in (C1,C2))
    dvec1=x-c1;dvec2=x-c2;sep=c1-c2
    rho=ref.rho.subs(dict(zip(ref.x,x)))
    v=ref.v.subs(dict(zip(ref.x,x)))
    grad=v.jacobian(x)
    geom_v=-sep;geom_rho=-2*dvec1.dot(sep)
    dv=geom_v if not response else sum(((1+2*k)*dv0 for k,dv0 in enumerate((dvec1,dvec2),1)),Matrix([q(f'W_{i}_change') for i in (1,2,3)]))
    drho=geom_rho if not response else q('varrho_change')+sum(2+k+k*dv0.dot(dv0) for k,dv0 in enumerate((dvec1,dvec2),1))
    de, dz=(q('eta_change'),q('zeta_change')) if response else (0,0)
    dgrad=dv.jacobian(x)
    shear=lambda g:g+g.T-R(2,3)*s.trace(g)*eye(3)
    dshear=de*shear(grad)+q('eta')*shear(dgrad)+(dz*s.trace(grad)+q('zeta')*s.trace(dgrad))*eye(3)
    other=Matrix(3,3,lambda i,j:q(f'Pi_o_{i+1}_{j+1}_change')) if response else s.zeros(3)
    dpi=drho*v*v.T+rho*(dv*v.T+v*dv.T)+2*rho*drho*eye(3)-dshear+other
    dG=drho*v+rho*dv
    out=[]
    for sign,b,rate in [(-1,1,R(1,2)),(1,2,-R(2,3))]:
        n=ref.N*sign
        integrand=-dpi*n+dG*(ref.C+rate*ref.N).dot(n)
        prim=integrand.subs(dict(zip(x,ref.x)))
        out.append(ref.integrate(prim,b))
    return list(map(sum,zip(*out)))
if __name__=='__main__':
    ix=scan()
    geom=family_variation(False);resp=family_variation(True)
    Q=list(map(sum,zip(*(r['crossing'] for r in ref.shells))))
    truth={'E_geo':[g+2*q for g,q in zip(geom,Q)],'E_resp':resp,'separation_derivative':[g+r for g,r in zip(geom,resp)]}
    print('exact crossing reference / (4*pi)',truth,flush=True)
    ctx=Native()
    for i in (1,2,3):
        for field in truth:
            tag=f'reading/both_surfaces/sum_crossing/component_{i}/{field}'
            a,b=[raw(ix,k+'_reading',tag) for k in ('py','wl')]
            vals=[ctx.value(m.parse(v,k))/(4*np.pi) for v,k in [(a,'py'),(b,'wl')]]
            print('component',i,field,'truth / (4*pi)',truth[field][i-1],'py,wl',vals,flush=True)
            if i==1 and field=='E_resp':
                mut=b.replace('change[oo[1, 2]]','change[oo[1, 1]]',1)
                changed=ctx.value(m.parse(mut,'wl'))/(4*np.pi)
                print('inner off-diagonal collapse, py, wl_mut / (4*pi)',vals[0],changed,'residual',vals[0]-changed,'analytic residual',R(1,3),flush=True)
