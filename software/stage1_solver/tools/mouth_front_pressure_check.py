"""Independent polynomial volume test for a relaxed, pressure-loaded three-ball."""
import sympy as sp


r, mu = sp.symbols('r mu', real=True)
b, c, sigma, receiver = sp.symbols('b c sigma receiver', real=True)
xi = r*mu  # Unit radius; polar axis points toward the companion.
surface = 2*sp.pi
volume = 4*sp.pi/3


def lap(h):
    radial = sp.diff(r**2*h.diff(r),r)
    angular = sp.diff((1-mu**2)*h.diff(mu),mu)
    return sp.simplify((radial+angular)/r**2)


def sideways(h):
    return mu*h.diff(r)+(1-mu**2)*h.diff(mu)/r


def integral(f):
    return sp.factor(surface*sp.integrate(sp.expand(r**2*f),(mu,-1,1),(r,0,1)))


h0 = b*(1-r**4)
h1 = c*xi*(1-r**2)
p0, p1 = -sigma*lap(h0), -sigma*lap(h1)
load = integral(-p1*sideways(h0))
shape = integral(-p0*sideways(h1))
edge = -h0.diff(r).subs(r,1)*integral(xi*p1)
print('P0, dP =',p0,p1)
print('direct, shape, sum, edge per ball volume =', *(sp.factor(term/volume) for term in (load,shape,load+shape,edge)))
assert sp.simplify(load+shape-edge)==0 and shape!=0 and sp.simplify(edge/volume-8*sigma*b*c)==0
# A receiving pressure proportional to the same odd component changes its
# moment without changing the mean pressure or the integrated linear intake.
print('mean receiver perturbation =',integral(receiver*p1)/volume)
net = sp.factor((1-receiver)*edge/volume)
print('net load per volume; receiver gains 0,1,2 =',net,*(net.subs(receiver,k) for k in (0,1,2)))
assert integral(receiver*p1)==0 and net.subs(receiver,1)==0 and net.subs(receiver,0)==-net.subs(receiver,2)
