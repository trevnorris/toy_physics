#!/usr/bin/env python3
"""
Preliminary SymPy workspace for closing the missing one-body 2PN response sector.

This script does four things:
1) Reproduces the current 2PN test-mass mismatch from the existing harness.
2) Shows that a single quadratic correction to the Bernoulli denominator,
       D(u) = 1 - 4 u + chi u^2,
   changes only the U^2 v^2 slot at 2PN.
3) Solves for the unique value chi = 8 that reproduces the exact
   Schwarzschild-isotropic test-mass target through 2PN when the static sector is
   kept in the current pair-counted form.
4) Extends the 1DOF throat-response closure to second order in SymPy and computes the
   breathing series a(rho), the equilibrium energy series F_eq(rho), and a minimal
   quadratic geometry dressing coefficient needed to realize chi = 8.

Conventions:
    u  = U/c^2
    nu = v^2/c^2

All expansions use weighted PN counting with wt(u)=wt(nu)=1.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def weighted_series(expr: sp.Expr, vars_with_weight1: list[sp.Symbol], order: int) -> sp.Expr:
    """Series-expand with a single bookkeeping parameter lam."""
    lam = sp.symbols("lam", real=True)
    subs = {var: lam * var for var in vars_with_weight1}
    return sp.expand(sp.series(sp.expand(expr.subs(subs)), lam, 0, order + 1).removeO().subs(lam, 1))


def coeff(expr: sp.Expr, *powers: tuple[sp.Symbol, int]) -> sp.Expr:
    out = sp.expand(expr)
    for sym, power in powers:
        out = sp.expand(sp.Poly(out, sym).coeff_monomial(sym ** power))
    return sp.simplify(out)


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

u, nuv, chi = sp.symbols("u nu chi", real=True)
eps = sp.symbols("eps", real=True)
rho, a = sp.symbols("rho a", positive=True)
A1, A2, A3 = sp.symbols("A1 A2 A3")


# ---------------------------------------------------------------------------
# 1) Baseline harness mismatch and exact isotropic target
# ---------------------------------------------------------------------------

banner("1) Baseline 2PN one-body candidate and the isotropic mismatch")

L_baseline = -(1 - u) * sp.sqrt(1 - nuv / (1 - 4 * u)) - sp.Rational(1, 2) * u ** 2 + sp.Rational(1, 4) * u ** 3
L_baseline_2pn = weighted_series(L_baseline, [u, nuv], 3)

L_iso_exact = -sp.sqrt(((1 - u / 2) / (1 + u / 2)) ** 2 - (1 + u / 2) ** 4 * nuv)
L_iso_2pn = weighted_series(L_iso_exact, [u, nuv], 3)

baseline_residual = sp.expand(L_baseline_2pn - L_iso_2pn)

print("Baseline candidate through 2PN:")
print(L_baseline_2pn)
print("\nExact isotropic target through 2PN:")
print(L_iso_2pn)
print("\nResidual baseline - isotropic:")
print(baseline_residual)
print("\nU^2 v^2 coefficient mismatch:")
print(coeff(baseline_residual, (u, 2), (nuv, 1)))


# ---------------------------------------------------------------------------
# 2) Minimal quadratic denominator correction
# ---------------------------------------------------------------------------

banner("2) Minimal denominator correction D(u) = 1 - 4u + chi u^2")

L_corr = -(1 - u) * sp.sqrt(1 - nuv / (1 - 4 * u + chi * u ** 2)) - sp.Rational(1, 2) * u ** 2 + sp.Rational(1, 4) * u ** 3
L_corr_2pn = weighted_series(L_corr, [u, nuv], 3)

print("Response-corrected candidate through 2PN:")
print(L_corr_2pn)
print("\nKey coefficients as functions of chi:")
print("free v^6 coefficient      =", coeff(L_corr_2pn, (nuv, 3)))
print("U v^4 coefficient         =", coeff(L_corr_2pn, (u, 1), (nuv, 2)))
print("U^2 v^2 coefficient       =", coeff(L_corr_2pn, (u, 2), (nuv, 1)))
print("static U^3 coefficient    =", coeff(L_corr_2pn.subs(nuv, 0), (u, 3)))

chi_solution = sp.solve(sp.Eq(coeff(L_corr_2pn, (u, 2), (nuv, 1)), coeff(L_iso_2pn, (u, 2), (nuv, 1))), chi)[0]
print("\nSolve for isotropic match in the U^2 v^2 slot:")
print("chi =", chi_solution)

L_corr_iso_2pn = sp.expand(L_corr_2pn.subs(chi, chi_solution))
print("\nResidual after substituting chi = 8:")
print(sp.simplify(L_corr_iso_2pn - L_iso_2pn))


# ---------------------------------------------------------------------------
# 3) One-DOF throat response closure extended to second order
# ---------------------------------------------------------------------------

banner("3) One-DOF throat response: a(rho), F_eq(rho), and the pure response factor")

# Normalized n=5 closure with the 11:2:5 operating-point partition at rho=a=1.
F = 11 * rho ** 2 / a + 2 / (rho * a ** 2) + 5 * rho ** 5 * a ** 3

a_series = 1 + A1 * eps + A2 * eps ** 2 + A3 * eps ** 3
stationary = sp.expand(sp.series(sp.diff(F, a).subs({rho: 1 + eps, a: a_series}), eps, 0, 4).removeO())
poly = sp.Poly(stationary, eps)
sol = sp.solve(
    [
        sp.Eq(poly.coeff_monomial(eps), 0),
        sp.Eq(poly.coeff_monomial(eps ** 2), 0),
        sp.Eq(poly.coeff_monomial(eps ** 3), 0),
    ],
    [A1, A2, A3],
    dict=True,
)[0]

print("a(rho) = 1 + A1 eps + A2 eps^2 + A3 eps^3 + ...")
for key in [A1, A2, A3]:
    print(f"  {key} = {sp.simplify(sol[key])}  ≈  {sp.N(sol[key], 12)}")

F_eq = sp.expand(sp.series(F.subs({rho: 1 + eps, a: a_series}).subs(sol), eps, 0, 4).removeO())
F_eq_norm = sp.expand(sp.simplify(F_eq / F_eq.subs(eps, 0)))

print("\nEquilibrium energy ratio F_eq/F0:")
print(F_eq_norm)

# The base density-driven rest-mass scaling kappa_rho=1 is carried separately in the 1PN ledger.
# Dividing by the background density ratio isolates the pure internal-response factor.
rho_ratio = 1 + eps
R_pv_eps = sp.expand(sp.series(F_eq_norm / rho_ratio, eps, 0, 4).removeO())
print("\nPure response factor R_PV(eps) = (F_eq/F0)/(rho/rho0):")
print(R_pv_eps)

print("\nChecks at the operating point:")
print("  linear response coefficient in eps  =", coeff(R_pv_eps, (eps, 1)))
print("  expected kappa_PV                   = 3/2")


# ---------------------------------------------------------------------------
# 4) Compose with exact n=5 Bernoulli density and extract a(u)
# ---------------------------------------------------------------------------

banner("4) Compose the 1DOF closure with exact n=5 Bernoulli density")

rho_u = sp.expand(sp.series((1 - 4 * u) ** sp.Rational(1, 4), u, 0, 4).removeO())
a_u = sp.expand(sp.series(a_series.subs(sol).subs(eps, rho_u - 1), u, 0, 4).removeO())
R_pv_u = sp.expand(sp.series(R_pv_eps.subs(eps, rho_u - 1), u, 0, 4).removeO())

print("rho/rho0 from exact Bernoulli:")
print(rho_u)
print("\na(u) from the 1DOF closure + Bernoulli composition:")
print(a_u)
print("\nR_PV(u) from the same composition:")
print(R_pv_u)


# ---------------------------------------------------------------------------
# 5) A minimal geometry dressing that realizes chi = 8
# ---------------------------------------------------------------------------

banner("5) Minimal quadratic geometry dressing for the missing dynamic slot")

# Let delta_a(u) = a(u) - 1. If the effective dynamic denominator takes the form
#     D_eff(u) = (1 - 4u) * [1 + mu * delta_a(u)^2] + O(u^3),
# then mu is uniquely fixed by the isotropic target chi = 8.
delta_a_u = sp.expand(a_u - 1)
mu = sp.simplify(chi_solution / coeff(delta_a_u ** 2, (u, 2)))
D_eff = sp.expand(sp.series((1 - 4 * u) * (1 + mu * delta_a_u ** 2), u, 0, 3).removeO())

print("delta_a(u) = a(u) - 1:")
print(delta_a_u)
print("\nRequired quadratic geometry coefficient mu in")
print("  D_eff(u) = (1 - 4u) [1 + mu (delta_a)^2] + O(u^3)")
print("is")
print("  mu =", mu, "≈", sp.N(mu, 12))
print("\nResulting D_eff(u) through u^2:")
print(D_eff)

# For comparison, the raw squared resonance frequency of the monopole channel scales like
# omega_0^2 ~ c_s^2 / L^2 with L ∝ a.  This already gives a positive u^2 shift, but not the full 8.
D_res_raw = sp.expand(sp.series((1 - 4 * u) / a_u ** 2, u, 0, 3).removeO())
print("\nRaw resonance-style proxy c_s^2/L^2 through u^2:")
print(D_res_raw)
print("u^2 coefficient in c_s^2/L^2:")
print(coeff(D_res_raw, (u, 2)), "≈", sp.N(coeff(D_res_raw, (u, 2)), 12))

# A smooth port-normalization factor can upgrade the raw resonance proxy to the exact
# denominator target without touching the already-frozen static sector.
p1, p2 = sp.symbols("p1 p2", real=True)
P_port = 1 + p1 * u + p2 * u ** 2
D_res_norm = sp.expand(sp.series(D_res_raw * P_port, u, 0, 3).removeO())
port_solution = sp.solve(
    [
        sp.Eq(coeff(D_res_norm, (u, 1)), -4),
        sp.Eq(coeff(D_res_norm, (u, 2)), 8),
    ],
    [p1, p2],
    dict=True,
)[0]
print("\nOne explicit port-normalization route with D_eff = (c_s^2/L^2) P_port(u):")
print("  p1 =", sp.simplify(port_solution[p1]), "≈", sp.N(port_solution[p1], 12))
print("  p2 =", sp.simplify(port_solution[p2]), "≈", sp.N(port_solution[p2], 12))
print("  check ->", sp.expand(D_res_norm.subs(port_solution)))


# ---------------------------------------------------------------------------
# 6) Compact conclusion block
# ---------------------------------------------------------------------------

banner("6) Provisional conclusion from the SymPy prototype")
print("- The current self+static 2PN candidate misses the isotropic target only by +4*u^2*nu.")
print("- Keeping the static sector exactly as in the harness, the unique minimal one-body fix is")
print("      D(u) = 1 - 4u + 8u^2 + O(u^3).")
print("- In the 1DOF throat closure, the breathing series is already nontrivial:")
print("      a(u) = 1 + (57/64)u + (298821/131072)u^2 + ...")
print("- A minimal quadratic geometry dressing of the dynamic denominator, proportional to")
print("      (delta_a)^2, can realize the needed chi = 8 without disturbing the 1PN slot.")
print("- The raw resonance proxy c_s^2/L^2 has the correct sign for the u^2 shift but only gives")
print("      324075/65536 ≈ 4.945, so a pure sound-speed + length shift is not yet enough.")
print("- This points to a concrete next derivation target: compute the quadratic mouth/port")
print("  normalization or DtN reduction term that upgrades the raw dynamic correction to chi = 8.")
