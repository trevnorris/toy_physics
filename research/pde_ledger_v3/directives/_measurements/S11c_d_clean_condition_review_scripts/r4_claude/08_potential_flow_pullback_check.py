#!/usr/bin/env python3
"""Leg script 08: independent derivation of the rest-frame potential-flow pullback on class P
from the supplied bulk laws (S11c-b spec :95-97 = S11c-a §1b):
   v_bulk = grad4 phi,  dp = -rho_m d_t phi,  d_t^2 phi = c_s0^2 lap4 phi,
   j = rho_4D v_bulk,   d_t rho_4D + div4 j = 0,
with phi = Phi(x1,w) exp(i(k2 x2 - omega t)).  Compares each derived trace with directive B
§2.4's displayed equations (residual printed), and derives the density/current traces that
B §2.4 does not list although the engine codomain objects depend on them (leg script 01).
Prints computed objects only.
"""
import sympy as sp

x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
k2, om, rho_m, c = sp.symbols("k2 omega rho_m c_s0", positive=True)
Phi = sp.Function("Phi")(x1, w)
E = sp.exp(sp.I * (k2 * x2 - om * t))
phi = Phi * E
X4 = (x1, x2, x3, w)
v = [sp.diff(phi, q) for q in X4]
dp = -rho_m * sp.diff(phi, t)
# wave equation solved for d_w^2 Phi
wave = sp.diff(phi, t, 2) - c**2 * sum(sp.diff(phi, q, 2) for q in X4)
Phi_ww = sp.solve(sp.Eq(sp.simplify(wave / E), 0), sp.diff(Phi, w, 2))[0]
Psi = sp.diff(Phi, w)

def strip(e_):
    return sp.simplify(sp.expand(e_ / E))

derived = {
    "delta_p": strip(dp),
    "v1": strip(v[0]), "v2": strip(v[1]), "v3": strip(v[2]), "v4": strip(v[3]),
    "d_w delta_p": strip(sp.diff(dp, w)),
    "d_w v1": strip(sp.diff(v[0], w)), "d_w v2": strip(sp.diff(v[1], w)),
    "d_w v3": strip(sp.diff(v[2], w)),
    "d_w v4": strip(sp.diff(v[3], w)).subs(sp.diff(Phi, w, 2), Phi_ww),
}
B_eq = {
    "delta_p": sp.I * rho_m * om * Phi,
    "v1": sp.diff(Phi, x1), "v2": sp.I * k2 * Phi, "v3": sp.Integer(0), "v4": Psi,
    "d_w delta_p": sp.I * rho_m * om * Psi,
    "d_w v1": sp.diff(Psi, x1), "d_w v2": sp.I * k2 * Psi, "d_w v3": sp.Integer(0),
    "d_w v4": (k2**2 - om**2 / c**2) * Phi - sp.diff(Phi, x1, 2),
}
for key in derived:
    print("PULLBACK", key, "DERIVED", derived[key], "DIRECTIVE_B", B_eq[key],
          "RESIDUAL", sp.simplify(derived[key] - B_eq[key]))

# Density/current traces from the supplied continuity law (bulk background rho_m, no flow):
drho = sp.Function("drho")(x1, w) * E
cont = sp.diff(drho, t) + rho_m * sum(sp.diff(v[i], X4[i]) for i in range(4))
drho_amp = sp.solve(sp.Eq(sp.simplify(cont / E), 0), sp.Function("drho")(x1, w))[0]
drho_amp = sp.simplify(drho_amp.subs(sp.diff(Phi, w, 2), Phi_ww))
print("DERIVED_DENSITY_TRACE drho =", drho_amp, " ; drho - dp/c^2 =",
      sp.simplify(drho_amp - derived["delta_p"] / c**2))
j = [sp.simplify(rho_m * derived[k]) for k in ("v1", "v2", "v3", "v4")]
print("DERIVED_CURRENT_TRACE j =", j)
print("DERIVED_CURRENT_TRACE_COMPONENT_3", j[2])
