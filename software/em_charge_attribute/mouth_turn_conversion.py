"""Follow the signed turning through an axisymmetric order-conversion layer."""

import sympy as sp


def check(label, actual, expected):
    result = sp.trigsimp(sp.factor(actual))
    print(f"{label} = {result}", flush=True)
    if sp.simplify(sp.expand_trig(result-expected)) != 0:
        raise AssertionError(label)


# The meridian runs from the flat feed toward the bulk. Its normal is
# continued through a vertical tangent without resetting its w component.
theta, theta_rate, r = sp.symbols("theta theta_rate r", real=True, positive=True)
for s in (1,-1):
    tangent = sp.Matrix([-sp.cos(theta), s*sp.sin(theta)])
    normal = sp.Matrix([s*sp.sin(theta), sp.cos(theta)])
    axial = normal.dot(tangent.diff(theta))*theta_rate
    angular = -normal[0]/r
    density = 4*sp.pi*r**2*axial*angular**2/sp.pi**2
    check(f"turn per unit arclength, orientation {s:+}", density,
          4*s*sp.sin(theta)**2*theta_rate/sp.pi)


def count(angle):
    return (2*angle-sp.sin(2*angle))/sp.pi


for label, angle, expected in [("cap",0,0), ("vertical",sp.pi/2,1), ("fully flared",sp.pi,2)]:
    check(f"upper mouth through {label}", count(angle), expected)
check("inclination susceptibility at vertical", count(theta).diff(theta).subs(theta,sp.pi/2),
      4/sp.pi)

# The order fraction is one before the layer and zero after it. u runs
# from -1 to +1 across the layer. Integrate chi d(turn) directly and add
# the fully ordered approach from the feed. The skew changes where the
# loss occurs, preserving its integral and keeping it nonnegative.
u, b = sp.symbols("u b", real=True)
loss = sp.Rational(3,4)*(1-u**2)*(1+b*u)
chi = 1-sp.integrate(loss,(u,-1,u))
check("order normalization", chi.subs(u,1), 0)
check("vertical collar with arbitrary skew", sp.integrate(loss,(u,-1,1))*count(sp.pi/2), 1)
bend = sp.pi/2+sp.pi*u/6
ordered = sp.simplify(count(sp.pi/3)+sp.integrate(chi*count(bend).diff(u),(u,-1,1)))
via_loss = sp.integrate(loss*count(bend),(u,-1,1))
check("two expressions for ordered count", ordered-via_loss, 0)
check("symmetric loss across inclined layer", ordered.subs(b,0), 1)
for skew in (-sp.S.Half, sp.S.Zero, sp.S.Half):
    print(f"ordered turn at skew {skew} = {float(ordered.subs(b,skew)):.12f}", flush=True)

# On a smooth chart a current cut off by order is generally not conserved.
t, q = sp.symbols("t q",real=True)
primitive = sp.Function("B")(q,t)
order = sp.Function("chi")(q,t)
rho, current = primitive.diff(q), -primitive.diff(t)
check("order-weighted continuity residual", (order*rho).diff(t)+(order*current).diff(q),
      rho*order.diff(t)+current*order.diff(q))

# Static probe derivative at fixed matter, order and geometry. Write the
# time part of the gauged matter action before differentiating. Other
# stored or gradient energies have no explicit temporal probe.
amp, phase_dot, hbar, qstar, a0, probe, shear, order_energy = sp.symbols(
    "amp phase_dot hbar qstar A0 f shear order_energy",real=True)
cov_dt = sp.I*amp*(phase_dot+qstar*probe*a0/hbar)
lagrangian = sp.I*hbar*amp*(cov_dt-sp.conjugate(cov_dt))/2-order_energy-order*shear
check("electric source for probe f", -lagrangian.diff(a0), qstar*amp**2*probe)
charge, f_lip, f_bulk = sp.symbols("Q f_lip f_bulk",real=True)
projection = charge*(f_lip-f_bulk)
check("opposite physical sources, flat probe", projection.subs(f_bulk,f_lip), 0)
