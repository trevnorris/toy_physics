"""Signed layer count on a smooth, tilted receiving flare at fixed mouth radius."""
import sympy as sp


def check(label, expr, expected):
    result = sp.factor(expr)
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise AssertionError(label)


# r/R=u, cos(theta)=c. Project the continuous face normal into the
# original brane directions, then normalize it only where nonzero.
u, beta, eps = sp.symbols("u beta epsilon", positive=True)
c, x = sp.symbols("c x", real=True)
nz = (u*c-beta)/sp.sqrt(u**2+beta**2-2*u*beta*c)
check("angular Jacobian per dc dphi", nz.diff(c),
      u**2*(u-beta*c)/(u**2+beta**2-2*u*beta*c)**sp.Rational(3, 2))
for uv, bv in [(2, 0), (2, 1), (2, 3), (3, 4)]:
    f = nz.subs({u: uv, beta: bv})
    check(f"degree u={uv} beta={bv}", (f.subs(c, 1)-f.subs(c, -1))/2, int(uv > bv))

loss = 30*x**2*(1-x)**2
lower = sp.symbols("lower", real=True)
remaining = sp.integrate(loss, (x, lower, 1))
for bv in [0, sp.Rational(9, 4), sp.Rational(5, 2), sp.Rational(11, 4), 4]:
    cutoff = min(1, max(0, bv-2))
    print(f"tilt={bv} layer count={remaining.subs(lower, cutoff)}")
check("change of count inside the layer", remaining.diff(lower), -loss.subs(x, lower))

# The physical face is a graph with Hessian epsilon/R times the identity.
# Rotate the instantaneous gradient to the third coordinate at any point.
slope = sp.symbols("slope", real=True)
R, r, density, mass, rate = sp.symbols("R r n0 m D", positive=True)
p = sp.Matrix([0, 0, slope])
metric = sp.eye(3)+p*p.T
w = sp.sqrt(metric.det())
mean = sp.trace(metric.inv()*sp.eye(3)*eps/(R*w))/3
check("volume determinant at zero projection", metric.det().subs(slope, 0), 1)
check("mean curvature magnitude at zero projection", mean.subs(slope, 0), eps/R)
check("mean curvature magnitude divided by epsilon/R",
      mean*R/eps, (3+2*slope**2)/(3*(1+slope**2)**sp.Rational(3, 2)))

# Follow the graph as it moves, with physical surface density n0/W and
# positive radial velocity in material coordinates. Use number throughput D.
metric_factor = sp.sqrt(1 + eps**2*((r/R)**2 + beta**2 - 2*beta*r*c/R))
surf_density = density/metric_factor
v = rate/(4*sp.pi*density*r**2)
flux = 2*sp.pi*sp.integrate(r**2*metric_factor*surf_density*v, (c, -1, 1))
check("number rate across every section", flux, rate)
check("material number divergence in graph coordinates",
      sp.diff(r**2*metric_factor*surf_density*v, r), 0)
radial_slope = eps*(r/R-beta*c)
static_ke = sp.integrate(sp.pi*mass*density*r**2*v**2*
                         sp.integrate(1+radial_slope**2, (c, -1, 1)), (r, 2*R, 3*R))
check("static collar kinetic energy in units m D^2/(8 pi n0 R)",
      static_ke*8*sp.pi*density*R/(mass*rate**2),
      sp.Rational(1, 6)+eps**2*(1+beta**2/18))
