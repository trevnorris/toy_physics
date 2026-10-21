#!/usr/bin/env python3
"""Two exact solution probes, not an S14a construction or acceptance test."""
import resource
resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
import sympy as s

t, r, w = s.symbols('t r w', real=True)
T, om = s.symbols('T omega', positive=True)
c = s.Rational(1, 2)

def out(name, x):
    print(name + ' = ' + str(s.simplify(x)))

def div(ur, uw):
    return s.diff(r**2 * ur, r) / r**2 + s.diff(uw, w)

print('A: t>=0, omega>0, 1<=r<=2, faces 0 and h=1+t/4; n=1 in all regions')
print('A: chi=1/2 inside, zero outside; J_chi=0; same velocity on both sides of each face')
h = 1 + t/4
n = s.S.One
ur = (-2+s.cos(om*t))/(r**2*h)
uw = s.S.Zero
gr = gd = 1/(r**2*h)
out('A u_r', ur)
out('A Gamma_return=Gamma_drain in the interior', gr)
out('A material regular residual', s.diff(n,t)+div(n*ur,n*uw))
out('A order regular residual', s.diff(c*n,t)+div(c*n*ur,c*n*uw)-n*(gr-gd))
aminus = uw
aplus = uw - s.diff(h,t)
Cminus = c*n*aminus
Cplus = -c*n*aplus
out('A material face jump residuals', s.Tuple(n*aminus-n*aminus,n*aplus-n*aplus))
out('A order face coefficients', s.Tuple(Cminus,Cplus))
sigma = s.integrate(n,(w,0,h))
F = s.integrate(n*ur,(w,0,h))
K = s.diff(s.integrate(r**2*sigma,(r,1,2)),t)
face = s.integrate(r**2*(-n*aminus+n*aplus),(r,1,2))
core = -(r**2*F).subs(r,1)
phi = (r**2*F).subs(r,2)
J = face+core
src = s.integrate(r**2*(s.integrate(n*(gr-gd),(w,0,h))+Cminus+Cplus),(r,1,2))
ret = s.integrate(r**2*s.integrate(n*gr,(w,0,h)),(r,1,2))
drain = s.integrate(r**2*s.integrate(n*gd,(w,0,h)),(r,1,2))
out('A return,drain totals/(4*pi)', s.Tuple(ret,drain))
out('A storage,face,core,J,outer flux/(4*pi)', s.Tuple(K,face,core,J,phi))
out('A material shell residual/(4*pi)', K+J+phi)
out('A order shell residual/(4*pi)', c*K+c*core+c*phi-src)
meanJ = s.integrate(J,(t,0,T))/T
dcJ = s.limit(meanJ,T,s.oo)
slowJ = s.limit(J,om,0,dir='+')
out('A finite-window mean J/(4*pi)', meanJ)
out('A DC J/(4*pi), omega>0', dcJ)
out('A fixed-t omega->0 J/(4*pi)', slowJ)
out('A DC minus slow J/(4*pi)', dcJ-slowJ)
out('A discarded J/(4*pi) in omega->0', s.limit(J-dcJ,om,0,dir='+'))
out('A DC projected material storage', s.diff(sigma,t))
out('A mean projected density', s.integrate(sigma,(t,0,T))/T)
out('A mean projected density T->infinity', s.limit(s.integrate(sigma,(t,0,T))/T,T,s.oo))

print('B: fixed faces 0,1; n=1+cos(omega*t)/4; chi=1/2 inside, zero outside; J_chi=0')
n = 1+s.cos(om*t)/4
ndot = s.diff(n,t)
ur = (-1/r**2-r*ndot/3)/n
gr = gd = 1/(r**2*n)
zz = s.symbols('zz',real=True)
out('B density range', s.calculus.util.function_range(1+zz/4,zz,s.Interval(-1,1)))
out('B u_r', ur)
out('B material regular residual', s.diff(n,t)+div(n*ur,0))
out('B order regular residual', s.diff(c*n,t)+div(c*n*ur,0)-n*(gr-gd))
fixed_face_factor = uw-s.diff(s.S.Zero,t)
out('B material and order face coefficients', s.Tuple(n*fixed_face_factor-n*fixed_face_factor,n*fixed_face_factor-n*fixed_face_factor,c*n*fixed_face_factor,-c*n*fixed_face_factor))
out('B return,drain totals/(4*pi)', s.Tuple(s.integrate(r**2*n*gr,(r,1,2)),s.integrate(r**2*n*gd,(r,1,2))))
volume = r**2*ur
out('B volume flux/(4*pi)', volume)
out('B log(n) derivative minus n_dot/n', s.diff(s.log(n),t)-ndot/n)
out('B log(n) change over one period', (s.log(n).subs(t,2*s.pi/om)-s.log(n).subs(t,0)))
z = s.symbols('z',real=True)
# Tangent half-angle substitution in the period integral of 1/n.
inverse_mean = s.integrate(2/(s.Rational(5,4)+s.Rational(3,4)*z**2),(z,-s.oo,s.oo))/(2*s.pi)
out('B period mean 1/n', inverse_mean)
out('B DC volume flux/(4*pi), omega>0', -inverse_mean)
out('B fixed-t omega->0 volume flux/(4*pi)', s.limit(volume,om,0,dir='+'))
out('B DC minus slow volume flux/(4*pi)', -inverse_mean-s.limit(volume,om,0,dir='+'))
out('B period mean n', s.integrate(n,(t,0,2*s.pi/om))/(2*s.pi/om))
out('B period mean outer number flux/(4*pi)', s.integrate((r**2*n*ur).subs(r,2),(t,0,2*s.pi/om))/(2*s.pi/om))
