"""Independent of the comparator. Direct tensor contraction and spherical moments.

Weights below give the spherical mean of a Cartesian normal monomial. The
radial ansatz u_k=U_k r and rho_k=B_k r^2 makes every test field smooth on the
stated domain and leaves enough degree to distinguish the Cartesian slots.
Units are divided by 4*pi to leave rational exact surface integrals.
"""
from math import prod
from fractions import Fraction as F
from sympy.polys.rings import ring
from sympy import QQ

P, n0,n1,n2,h = ring('n0,n1,n2,h', QQ)
n=(n0,n1,n2)
I=range(3)

def df(k): return prod(range(k,0,-2))

def mean(p,order):
    return sum((F(int(c.numerator),int(c.denominator))*prod(df(i-1) for i in k[:3])/df(sum(k[:3])+1)
        for k,c in P(p).items() if k[3]==order and all(i%2==0 for i in k[:3])), F(0))

def dot(a,b): return sum(x*y for x,y in zip(a,b))

def surface(b, ns, rate, vc, geom, resp, swapped=False):
    c2=(1,-2,3); distance=(4,3,0)  # d=5 > 2 > 1
    c1=[a+(1+geom*h)*d for a,d in zip(c2,distance)]
    x=[a+b*v for a,v in zip(c2,n)]
    dx=[[a-c for a,c in zip(x,cent)] for cent in (c1,c2)]
    change=resp*h
    U=[1+3*change,-2+5*change]; B=[2+7*change,3-11*change]
    M=[[(-1)**(i+j)*(1+3*i+j) + (2+i+2*j)*change for j in I] for i in I]
    V=[sum(U[k]*dx[k][i] for k in (0,1)) +i+1 +(i+5)*change+dot(M[i],x) for i in I]
    rho=sum(B[k]*dot(dx[k],dx[k]) for k in (0,1))+8+13*change+dot((2+change,3-2*change,1+4*change),x)
    p=rho*rho+rho+1 # fixed pressure function, changing density
    grad=[[M[i][j]+(sum(U) if i==j else 0) for j in I] for i in I]
    div=sum(grad[i][i] for i in I)
    eta=9+x[0]+change*(8+x[2]); zeta=2+x[1]+change*(4-x[0])
    def other(i,j):
        return sum(((i+2)**2*(j+1)*(k+1) + change*(2+i+3*j)**2*(k+1)**2)*x[k] for k in I)+(1+i)*j+change*(2+j-i)
    rows={k:[] for k in ('carried','pressure','viscous','other','centre','radius','rest','transport','crossing')}
    for i in I:
        pieces=[-rho*V[i]*dot(V,n), -p*n[i],
                eta*(sum((grad[i][j]+grad[j][i])*n[j] for j in I)- QQ(2,3)*div*n[i])+zeta*div*n[i],
                -sum((other(j,i) if swapped else other(i,j))*n[j] for j in I),
                rho*V[i]*dot(vc,n),rho*V[i]*rate]
        allpieces=pieces+[sum(pieces[:4]),sum(pieces[4:]),sum(pieces)]
        for key,poly in zip(rows,allpieces):
            rows[key].append(tuple(mean(ns*b*b*poly,o) for o in (0,1)))
    return rows

def run():
    strata={'moving':((1,3,-2),(-2,4)),'stationary target':((0,0,0),(-2,4)),
            'stationary radii':((1,3,-2),(0,0)),'stationary target and radii':((0,0,0),(0,0)),
            'stationary source':((1,3,-2),(-2,4))}
    def fmt(vals):return ' '.join(str(v) for v in vals)
    for case,(v,rates) in strata.items():
        sums=[[surface(b,ns,rate,v,*mode) for b,ns,rate in zip((1,2),(-1,1),rates)] for mode in ((1,0),(0,1))]
        print(case)
        for key in ['rest','transport','crossing']:
            for side,(a,c) in enumerate([(0,1),(1,2),(0,2)]):
                vals=lambda mode,order:[sum(sums[mode][s][key][i][order] for s in range(a,c)) for i in I]
                q=vals(0,0); geo=[2*x+y for x,y in zip(q,vals(0,1))]; resp=vals(1,1)
                print(['in','out','net'][side],key,'Q',fmt(q),'E_geo',fmt(geo),'E_resp',fmt(resp))
    v=(1,3,-2)
    for reversed_slots in (False,True):
        rows=[surface(b,s,r,v,0,1,reversed_slots) for b,s,r in ((1,-1,-2),(2,1,4))]
        print('net other transposed',reversed_slots,'Q',fmt([sum(s['other'][i][0] for s in rows) for i in I]),'E_resp',fmt([sum(s['other'][i][1] for s in rows) for i in I]))
    # For a rate of a radial scalar sampled in a moving chart the convective
    # derivative of |x-X1|² is 2(x-X1).(v2+bdot*n-v1). Select each v separately.
    b=2; D=(4,3,0); direction=(F(2,3),F(-2,3),F(1,3)); dr=[b*x-d for x,d in zip(direction,D)]
    src=(3,-1,2); tgt=(1,3,-2); bd=4
    def material(v1,v2): return 2*dot(dr,[z+bd*t-y for z,t,y in zip(v2,direction,v1)])
    base=material(src,tgt)
    print('rate centre checks Q, no source, no target, parts',fmt([base,material((0,)*3,tgt),material(src,(0,)*3),base-material((0,)*3,tgt),base-material(src,(0,)*3)]))
    # Choice of different nonzero accelerations even when first jets vanish.
    # Taylor jets f''(0)=2*c2 with f'(0) controlled solely by c1.
    print('stationary first and higher jets (source,target,inner,outer)',[(0,2*c) for c in (7,-5,11,-13)])
    # Positive affine bases on signed orientations. phi in [0,2*pi], spherical
    # averages n_y and n_y^2 are 0 and 1/3. The half azimuth has a nonzero n_y.
    import mpmath as mp
    mp.mp.dps=65
    for sign in (-1,1):
        f=lambda r: mp.log(5+sign*r)/r
        exact=lambda r: mp.log(5)*mp.log(r)-mp.polylog(2,-sign*r/5)
        v=mp.quad(f,[1,2]); print('radial log sign',sign,'integral',mp.nstr(v,53),'primitive residual',mp.nstr(v-exact(mp.mpf(2))+exact(mp.mpf(1)),4))
    print('azimuth sin integral full/half',mp.nstr(mp.quad(mp.sin,[0,2*mp.pi]),4),mp.nstr(mp.quad(mp.sin,[0,mp.pi]),4))
    print('dimension base D,L,T surface: rho*vel^2*area',(1,4,-2),'momentum rho*vel*volume',(1,4,-1),'rate',(1,4,-2),'literal zero slots',[],'equal residual triples',[(0,0,0)])
if __name__=='__main__':run()
