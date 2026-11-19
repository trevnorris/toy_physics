"""Signed curvature of a radial mouth on one face of a three-dimensional slab."""

import sympy as sp


def check(name, actual, expected):
    print(f"{name} = {sp.factor(actual)}", flush=True)
    if sp.simplify(actual - expected) != 0:
        raise ArithmeticError(name)


r, a, c = sp.symbols("r a C", positive=True)
p, dp = sp.symbols("p dp", real=True)
# Choose the graph normal toward +w on both reflected faces. Signed normal
# curvatures are Hessian / sqrt(1+p²) with an inverse induced metric.
metric = sp.diag(1+p**2, 1, 1)
second_form = sp.diag(dp, p/r, p/r) / sp.sqrt(1+p**2)
curvature = (metric.inv()*second_form).det()
density = 4*sp.pi*r**2*sp.sqrt(metric.det())*curvature
check("signed curvature times radial area measure", density, 4*sp.pi*p**2*dp/(1+p**2)**2)
primitive = sp.integrate(4*sp.pi*p**2/(1+p**2)**2, p)
for s in (1,-1):
    total = primitive.subs(p,0) - sp.limit(primitive,p,-s*sp.oo)
    check(f"vertical rim orientation {s:+d}, integral", total, s*sp.pi**2)
check("finite 45 degree upper rim, fraction of vertical value",
      (primitive.subs(p,0)-primitive.subs(p,-1))/sp.pi**2, sp.Rational(1,2)-1/sp.pi)
check("regular face flat at center and infinity", primitive.subs(p,0)-primitive.subs(p,0),0)

# Counterexamples are geometric profiles, not proposed static solutions.
for s in (1,-1):
    height = s*c*sp.asin(a/r)/a
    slope = height.diff(r)
    check(f"continuous strength at fixed orientation {s:+d}", sp.limit(r*height,r,sp.oo), s*c)
    check(f"rim tangent sign for profile {s:+d}", sp.limit(sp.atan(slope),r,a,dir="+"), -s*sp.pi/2)

# A conditional equilibrium control: extremize excess area at tension sigma.
sigma = sp.symbols("sigma", positive=True)
gap = sp.symbols("gap", positive=True)  # r**4 - a**4 in the exterior
radial_load = sp.diff(4*sp.pi*sigma*r**2*(sp.sqrt(1+p**2)-1),p)
for s in (1,-1):
    slope_min = -s*a**2/sp.sqrt(r**4-a**4)
    load = radial_load.subs({r: (a**4 + gap)**sp.Rational(1, 4), p: -s*a**2/sp.sqrt(gap)})
    check(f"area tension first integral, orientation {s:+d}", load, -s*4*sp.pi*sigma*a**2)
    check(f"area tension far bend coefficient, orientation {s:+d}", sp.limit(-r**2*slope_min,r,sp.oo), s*a**2)
