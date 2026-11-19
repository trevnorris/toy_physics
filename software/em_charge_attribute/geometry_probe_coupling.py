"""Geometric source on a breathing, once-covered conversion layer."""
import sympy as sp


def check(label, value, expected):
    result = sp.factor(value)
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise ValueError(label)


# Radial material chart. The two boundaries are 2R and 2R+ell. The
# geometric current follows level sets of u, independently of the intake.
r, t, u = sp.symbols("r t u", real=True)
radius, width, action = [sp.Function(name)(t) for name in ("R", "ell", "I")]
sign = sp.Symbol("s", real=True)
coord = (r - 2*radius)/width
weight = sp.Function("h")(coord)
velocity = -coord.diff(t)/coord.diff(r)
density = sign*weight/(4*sp.pi*r**2*width)
flux = velocity*density


def balance(rho, current):
    return rho.diff(t) + (r**2*current).diff(r)/r**2


check("geometric current local continuity", balance(density, flux), 0)
g = sp.Function("g")(radius, action)
check("electric continuity residual per unit geometric density",
      balance(g*density, g*flux)/density, g.diff(t))

# The test profile is the order drop used on the relaxed path. Integrate the
# weak-probe source at arbitrary instantaneous R and positive width.
potential, coefficient = sp.symbols("phi g0", real=True)
zint = sp.Symbol("Zint", positive=True)
pullback = (4*sp.pi*r**2*density*width).subs(r, 2*radius+width*u)
count = sp.integrate(pullback.subs(sp.Function("h")(u), 30*u**2*(1-u)**2), (u, 0, 1))
interaction = -potential*coefficient*count/sp.sqrt(zint)
check("conditional normalized monopole", -interaction.diff(potential), coefficient*sign/sp.sqrt(zint))
check("monopole derivative along breathing at fixed coupling and medium",
      -sp.diff(interaction, potential, t), 0)

# Source of a common static electric mode on a mirror pair of depth cells.
# The two cell densities can be arbitrary functions of shape and wave action.
up, down = [sp.Function(name)(radius, action) for name in ("n_plus", "n_minus")]
q, hbar, phase = sp.symbols("q hbar phase", real=True)
probe = sp.symbols("probe", real=True)


def temporal_term(n):
    rate = sp.I*(phase+q*probe/hbar)
    return sp.I*hbar*n*(rate-sp.conjugate(rate))/2


def pair_action(first, second):
    return temporal_term(first)+temporal_term(second)


odd_source = -sp.diff(pair_action(up,down)-pair_action(down,up), probe)/2
check("matched native carrier odd monopole", odd_source, 0)
omega = sp.Function("omega")(radius)
wave_energy = action*omega
check("direct static probe derivative of neutral wave", wave_energy.diff(probe), 0)
check("wave radial load at fixed action", -wave_energy.diff(radius), -action*omega.diff(radius))

# Two flux boundaries of an exact current. A smooth polarization of
# compact support cannot create a net monopole by moving charge locally.
lo, hi = sp.symbols("r_in r_out", positive=True)
surface_moment = sp.Function("M")  # M(r) = 4 pi r^2 P_r(r)
polar_charge = sp.integrate(-surface_moment(r).diff(r), (r,lo,hi))
check("total radial polarization charge with both boundaries",
      polar_charge, surface_moment(lo)-surface_moment(hi))
