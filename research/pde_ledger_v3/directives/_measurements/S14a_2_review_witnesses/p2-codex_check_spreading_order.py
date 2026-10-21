#!/usr/bin/env python3
"""An expanding, finite-excess slab with constant global rates."""
import resource
resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
import sympy as s

t = s.symbols('t', nonnegative=True)
r, w = s.symbols('r w', real=True)
T, L = s.symbols('T L', positive=True)
z = s.symbols('z', nonnegative=True)
q = s.Integer(2)
e = s.Rational(1,2)
c = s.Rational(1,2)
def out(name,x):
    print(name+' = '+str(s.simplify(x)))

print('domain: t>=0, r>=a=1; L=R^3-1>0; lower face 0, upper face h')
print('n,u are the same in all three regions; chi is zero outside; J_chi=0 in all regions')
h = 1+e*s.exp(-(r**3-1)/(1+t))
ur = -q/r**2
uw = s.S.Zero
gr = gd = s.exp(-(r**3-1))/h
out('n,chi_inside,J_chi,u_r,u_w', s.Tuple(1,c,0,ur,uw))
out('h',h)
out('Gamma_return=Gamma_drain in the interior',gr)
out('material regular residual',s.diff(r**2*ur,r)/r**2)
out('order regular residual',s.diff(r**2*c*ur,r)/r**2-(gr-gd))
out('material face jump coefficients',s.Tuple(0,0))
relative_plus = uw-s.diff(h,t)-ur*s.diff(h,r)
Cplus = -c*relative_plus
out('order lower,upper face source coefficients',s.Tuple(0,Cplus))
# The exact substitution z=r^3-1 gives r^2 dr=dz/3.
excess = s.integrate(e*s.exp(-z/(1+t))/3,(z,0,L))
K = s.diff(excess,t)
ha = h.subs(r,1)
hR = 1+e*s.exp(-L/(1+t))
face = -K+q*(hR-ha)
core = q*ha
phi = -q*hR
J = core+face
face_source = c*(K-q*(hR-ha))
out('finite-shell excess material/(4*pi)',excess)
out('finite-shell storage/(4*pi)',K)
out('face,core,J,outer flux/(4*pi)',s.Tuple(face,core,J,phi))
out('material shell residual/(4*pi)',K+J+phi)
out('order shell residual/(4*pi)',c*K+c*core+c*phi-face_source)
regular_total = s.integrate(s.exp(-z)/3,(z,0,L))
out('regular return,drain shell totals/(4*pi)',s.Tuple(regular_total,regular_total))
out('global regular return,drain totals/(4*pi)',s.Tuple(s.limit(regular_total,L,s.oo),s.limit(regular_total,L,s.oo)))
out('global storage/(4*pi)',s.limit(K,L,s.oo))
out('global net face source/(4*pi)',s.limit(face_source,L,s.oo))
out('global J/(4*pi)',s.limit(J,L,s.oo))
out('time derivatives of global storage,face source,J',s.Tuple(*(s.diff(s.simplify(s.limit(x,L,s.oo)),t) for x in (K,face_source,J))))
meanK = (excess.subs(t,T)-excess.subs(t,0))/T
dcK = s.limit(meanK,T,s.oo)
out('DC finite-shell storage/(4*pi)',dcK)
out('DC then far storage/(4*pi)',s.limit(dcK,L,s.oo))
out('far then DC storage/(4*pi)',s.limit(K,L,s.oo))
primitive = (1+t)*s.exp(-L/(1+t))+L*s.Ei(-L/(1+t))
out('mean-exp primitive derivative residual',s.diff(primitive,t)-s.exp(-L/(1+t)))
meanExp = (primitive.subs(t,T)-primitive.subs(t,0))/T
dcPhi = -q*(1+e*s.limit(meanExp,T,s.oo))
farPhi = s.limit(phi,L,s.oo)
out('DC then far outer flux/(4*pi)',s.limit(dcPhi,L,s.oo))
out('far then DC outer flux/(4*pi)',farPhi)
out('order-of-limits difference in outer flux/(4*pi)',s.limit(dcPhi,L,s.oo)-farPhi)
out('DC finite-shell face channel/(4*pi)',-dcK+q*(1+e-ha))
out('DC finite-shell face source/(4*pi)',c*(dcK-q*(1+e-ha)))
out('far-first face channel/(4*pi)',s.limit(face,L,s.oo))
