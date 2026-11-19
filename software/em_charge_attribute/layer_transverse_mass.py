"""Low-frequency entrance traction from distributed momentum in one free layer."""

import sympy as sp


def check(label, measured, reference):
    error = sp.factor(measured - reference)
    print(f"{label}: measured={measured}, reference={reference}, residual={error}", flush=True)
    if error != 0:
        raise ValueError(label)


# U(a)=1.  M(q) is mass beyond q, M'=-rho, and B(a)=0 with
# B'=M/mu for U=1+omega^2 B+... . Integration by parts then yields
# the quartic entrance traction from the squared downstream mass.
q = sp.symbols("q", real=True)
M, B, mu = [sp.Function(name)(q) for name in ("M", "B", "mu")]
rho = -sp.diff(M, q)
integrand = (rho*B + sp.diff(M*B, q)).subs(sp.diff(B, q), M/mu)
check("quartic traction integrand after integration by parts", integrand, M**2/mu)

k, c, rho0 = sp.symbols("k c rho0", positive=True)
a, mass, compliance, tip_mass = sp.symbols("a Sigma H tip_mass", real=True)
omega = c*k
traction = omega**2*mass + omega**4*compliance
edge_ratio = traction/(rho0*c**2*k)
apparent = a + sp.series(sp.atan(edge_ratio), k, 0, 5).removeO()/k
check("leading apparent position", apparent.expand().coeff(k,0), a+mass/rho0)
print("first dispersive coefficient =", sp.factor(apparent.expand().coeff(k,2)))
velocity = sp.symbols("v_t")
kinetic_per_area = mass*velocity**2/2
check("momentum per area per uniform transverse velocity",
      sp.diff(kinetic_per_area, velocity)/velocity, mass)
print("apparent position with inertial end load =",
      sp.limit(a+sp.atan((traction+omega**2*tip_mass)/(rho0*c**2*k))/k, k, 0))

# Test two different densities with the same column on the same linear
# stiffness ramp. x=0 is the free outer edge, x=1 is the layer entrance.
x, u, h = sp.symbols("x u h", positive=True)
for name, density_weight in (("linear", x), ("smoothstep", x**2*(3-2*x))):
    density = rho0*density_weight
    downstream = h*sp.integrate(density.subs(x,u), (u,0,x))
    column = downstream.subs(x,1)
    stiffness = rho0*c**2*x
    fourth = h*sp.integrate(sp.cancel(downstream**2/stiffness), (x,0,1))
    word = sp.factor(apparent.subs({mass:column, compliance:fourth}))
    shift = sp.simplify(word-a)
    print(f"{name}: Sigma={column}, apparent penetration through k^2={shift.expand()}")
    check(f"{name} leading penetration", shift.expand().coeff(k,0), column/rho0)
