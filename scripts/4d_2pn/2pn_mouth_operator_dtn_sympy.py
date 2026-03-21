#!/usr/bin/env python3
"""
Preliminary SymPy workspace for the low-frequency monopole mouth operator.

Goal:
    Turn the next-step target from the 2PN inner-throat prototype into a more
    structural DtN statement.

What this script does:
  1) Re-derives the general one-body 2PN denominator conditions.
  2) Parameterizes the low-frequency monopole mouth operator
         Z00(omega;u) = Z2(u) omega^2 + Z4(u) omega^4 + O(omega^6)
     and identifies two normalized DtN invariants:
         C_s(u) := [Z4/Z2^3] / [Z4/Z2^3]_0   -> sound-speed-like factor,
         G(u)   := [Z4/Z2^2] / [Z4/Z2^2]_0   -> geometry-like factor.
  3) Shows that for the cylinder / Neumann-bottom unit-test branch,
         Z2 = -L/c_s^2,
         Z4 = -L^3/(3 c_s^4),
     so C_s = c_s^2/c_{s0}^2 and G = L/L0 exactly.
  4) Composes those exact DtN identities with the already-fixed n=5 Bernoulli
     continuation and the known 1PN breathing slope.
  5) Solves the unique minimal quadratic conservative response built from the
     DtN geometry invariant that preserves the frozen 1PN slot and closes the
     missing one-body 2PN U^2 v^2 term.
  6) Shows how the earlier raw resonance proxy / port-normalization fit is just
     a refactorization of the same DtN-invariant closure.

Conventions:
    u  = U/c^2    (dimensionless gravitational potential)
    nu = v^2/c^2  (dimensionless speed square)
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def weighted_series(expr: sp.Expr, vars_with_weight1: list[sp.Symbol], order: int) -> sp.Expr:
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

u, nuv = sp.symbols("u nu", real=True)
d1, d2 = sp.symbols("d1 d2", real=True)

z1, z2, w1, w2 = sp.symbols("z1 z2 w1 w2", real=True)
g1, g2 = sp.symbols("g1 g2", real=True)
alpha, beta, mu = sp.symbols("alpha beta mu", real=True)

ell1, ell2 = sp.symbols("ell1 ell2", real=True)
A1, A2 = sp.symbols("A1 A2", real=True)

eps = sp.symbols("eps", real=True)
rho, a = sp.symbols("rho a", positive=True)


# ---------------------------------------------------------------------------
# 1) General denominator conditions
# ---------------------------------------------------------------------------

banner("1) General one-body denominator conditions at 2PN")

D_general = 1 + d1 * u + d2 * u ** 2
L_general = -(1 - u) * sp.sqrt(1 - nuv / D_general)
L_general_2pn = weighted_series(L_general, [u, nuv], 3)

c1_uv = coeff(L_general_2pn, (u, 1), (nuv, 1))
c2_u2v = coeff(L_general_2pn, (u, 2), (nuv, 1))

print("For D(u) = 1 + d1 u + d2 u^2, the worldline expansion gives:")
print("  coeff[u v^2]   =", c1_uv)
print("  coeff[u^2 v^2] =", c2_u2v)
print("\nTo keep the frozen 1PN self coefficient +3/2, we need d1 = -4.")
print("Then coeff[u^2 v^2] =", sp.simplify(c2_u2v.subs(d1, -4)), ", so the exact isotropic target coeff 2 forces d2 = 8.")


# ---------------------------------------------------------------------------
# 2) Generic low-frequency DtN invariants
# ---------------------------------------------------------------------------

banner("2) Generic low-frequency DtN invariants from Z2(u), Z4(u)")

Z2hat = 1 + z1 * u + z2 * u ** 2
Z4hat = 1 + w1 * u + w2 * u ** 2

Ghat = sp.expand(sp.series(Z4hat / Z2hat ** 2, u, 0, 3).removeO())
Cshat = sp.expand(sp.series(Z4hat / Z2hat ** 3, u, 0, 3).removeO())

print("Take normalized DtN data:")
print("  Z2/Z20 =", Z2hat)
print("  Z4/Z40 =", Z4hat)
print("\nThen the normalized invariants")
print("  G(u)   = [Z4/Z2^2]/[Z4/Z2^2]_0")
print("  C_s(u) = [Z4/Z2^3]/[Z4/Z2^3]_0")
print("expand as:")
print("  G(u)   =", Ghat)
print("  C_s(u) =", Cshat)
print("\nSo the measurable invariant coefficients are")
print("  g1 =", coeff(Ghat, (u, 1)))
print("  g2 =", coeff(Ghat, (u, 2)))
print("  c1 =", coeff(Cshat, (u, 1)))
print("  c2 =", coeff(Cshat, (u, 2)))


# ---------------------------------------------------------------------------
# 3) Exact cylinder / Neumann-bottom branch
# ---------------------------------------------------------------------------

banner("3) Exact cylinder branch: Z2 = -L/c_s^2, Z4 = -L^3/(3 c_s^4)")

Lhat = 1 + ell1 * u + ell2 * u ** 2
Cs_exact = 1 - 4 * u

Z2hat_cyl = sp.expand(sp.series(Lhat / Cs_exact, u, 0, 3).removeO())
Z4hat_cyl = sp.expand(sp.series(Lhat ** 3 / Cs_exact ** 2, u, 0, 3).removeO())

Ghat_cyl = sp.expand(sp.series(Z4hat_cyl / Z2hat_cyl ** 2, u, 0, 3).removeO())
Cshat_cyl = sp.expand(sp.series(Z4hat_cyl / Z2hat_cyl ** 3, u, 0, 3).removeO())

print("Normalized cylinder DtN coefficients:")
print("  Z2/Z20 =", Z2hat_cyl)
print("  Z4/Z40 =", Z4hat_cyl)
print("\nExtracted DtN invariants:")
print("  G(u)   =", Ghat_cyl)
print("  C_s(u) =", Cshat_cyl)
print("\nChecks:")
print("  G(u) - L/L0   =", sp.simplify(Ghat_cyl - Lhat))
print("  C_s(u) - c_s^2/c_s0^2 =", sp.simplify(Cshat_cyl - Cs_exact))


# ---------------------------------------------------------------------------
# 4) Rebuild the throat breathing slope from the 1DOF closure
# ---------------------------------------------------------------------------

banner("4) Rebuild the known breathing slope from the 1DOF throat closure")

# Same normalized 11:2:5 closure used previously.
F = 11 * rho ** 2 / a + 2 / (rho * a ** 2) + 5 * rho ** 5 * a ** 3

stationary = sp.expand(sp.series(sp.diff(F, a).subs({rho: 1 + eps, a: 1 + A1 * eps + A2 * eps ** 2}), eps, 0, 3).removeO())
poly = sp.Poly(stationary, eps)
sol = sp.solve(
    [
        sp.Eq(poly.coeff_monomial(eps), 0),
        sp.Eq(poly.coeff_monomial(eps ** 2), 0),
    ],
    [A1, A2],
    dict=True,
)[0]

A1sol = sp.simplify(sol[A1])
A2sol = sp.simplify(sol[A2])

print("a(rho) = 1 + A1 eps + A2 eps^2 + ...")
print("  A1 =", A1sol, "≈", sp.N(A1sol, 12))
print("  A2 =", A2sol, "≈", sp.N(A2sol, 12))
print("\nSo d ln a / d ln rho at the operating point is", A1sol)
print("and the exact n=5 Bernoulli map rho/rho0 = (1-4u)^(1/4) gives eps = -u + O(u^2).")
print("Therefore the linear geometry coefficient in G(u)=L/L0=a(u) is")
print("  g1 = -A1 =", sp.simplify(-A1sol))


# ---------------------------------------------------------------------------
# 5) Minimal quadratic DtN-invariant closure for the missing 2PN slot
# ---------------------------------------------------------------------------

banner("5) Minimal quadratic DtN-invariant response closure")

G_generic = 1 + g1 * u + g2 * u ** 2
D_invariant = sp.expand(sp.series((1 - 4 * u) * (1 + alpha * (G_generic - 1) + beta * (G_generic - 1) ** 2), u, 0, 3).removeO())

print("Take the minimal conservative local ansatz built from the geometry invariant G(u):")
print("  D_eff(u) = C_s(u) [1 + alpha (G-1) + beta (G-1)^2]")
print("with C_s(u) = 1 - 4u through 2PN.")
print("\nIts expansion is")
print("  D_eff(u) =", D_invariant)
print("\nThe linear coefficient is", coeff(D_invariant, (u, 1)), ", so preserving the frozen 1PN slot forces alpha = 0.")
print("Then the quadratic coefficient is", coeff(D_invariant.subs(alpha, 0), (u, 2)), ", so matching d2 = 8 forces")
print("  beta = 8/g1^2")

# Specialize to the already-fixed 1PN slope g1 = 57/64.
g1_value = sp.simplify(-A1sol)
mu_value = sp.simplify(8 / g1_value ** 2)

print("\nUsing the known 1PN throat slope g1 = 57/64 gives")
print("  mu = beta =", mu_value, "≈", sp.N(mu_value, 12))

# Important structural point: only g1 matters at 2PN.
D_mu = sp.expand(sp.series((1 - 4 * u) * (1 + mu_value * (G_generic - 1) ** 2), u, 0, 3).removeO())
print("\nSubstituting mu = 8/g1^2 leaves")
print("  D_eff(u) =", D_mu)
print("showing that only g1 enters the 2PN slot; g2 does not appear until 3PN.")


# ---------------------------------------------------------------------------
# 6) Specialize to the current throat closure and verify exact isotropic match
# ---------------------------------------------------------------------------

banner("6) Specialize to the current throat closure and verify the isotropic one-body match")

ell1_value = g1_value
ell2_value = sp.simplify(A2sol + sp.Rational(3, 2) * (-A1sol))  # not actually used below; kept only for display if desired.
# Better to build a(u) directly from the exact Bernoulli composition to O(u^2).
rho_u = sp.expand(sp.series((1 - 4 * u) ** sp.Rational(1, 4), u, 0, 3).removeO())
a_u = sp.expand(sp.series((1 + A1sol * eps + A2sol * eps ** 2).subs(eps, rho_u - 1), u, 0, 3).removeO())

G_u = a_u
D_u = sp.expand(sp.series((1 - 4 * u) * (1 + mu_value * (G_u - 1) ** 2), u, 0, 3).removeO())

print("a(u) = G(u) from the current throat closure:")
print("  ", G_u)
print("\nThe DtN-invariant denominator becomes")
print("  D_eff(u) =", D_u)

# Verify the one-body worldline closes the isotropic target.
L_dtn = -(1 - u) * sp.sqrt(1 - nuv / D_u) - sp.Rational(1, 2) * u ** 2 + sp.Rational(1, 4) * u ** 3
L_dtn_2pn = weighted_series(L_dtn, [u, nuv], 3)

L_iso_exact = -sp.sqrt(((1 - u / 2) / (1 + u / 2)) ** 2 - (1 + u / 2) ** 4 * nuv)
L_iso_2pn = weighted_series(L_iso_exact, [u, nuv], 3)

print("\nResponse-corrected one-body candidate through 2PN:")
print("  ", L_dtn_2pn)
print("\nExact isotropic target through 2PN:")
print("  ", L_iso_2pn)
print("\nResidual candidate - isotropic:")
print("  ", sp.simplify(sp.expand(L_dtn_2pn - L_iso_2pn)))


# ---------------------------------------------------------------------------
# 7) Relation to the earlier raw resonance proxy / port-normalization fit
# ---------------------------------------------------------------------------

banner("7) Relation to the earlier raw resonance proxy / port-normalization fit")

Draw = sp.expand(sp.series((1 - 4 * u) / G_u ** 2, u, 0, 3).removeO())
Pport = sp.expand(sp.series(G_u ** 2 * (1 + mu_value * (G_u - 1) ** 2), u, 0, 3).removeO())

print("Raw resonance proxy from the DtN coefficients:")
print("  D_raw(u) = [Z2/Z4]/[Z20/Z40] =", Draw)
print("\nRequired multiplicative port factor:")
print("  P_port(u) =", Pport)
print("\nCheck that D_raw * P_port reproduces the exact target denominator:")
print("  ", sp.expand(sp.series(Draw * Pport, u, 0, 3).removeO()))

print("\nSo the earlier fitted port-normalization coefficients are not arbitrary:")
print("  p1 =", coeff(Pport, (u, 1)), "≈", sp.N(coeff(Pport, (u, 1)), 12))
print("  p2 =", coeff(Pport, (u, 2)), "≈", sp.N(coeff(Pport, (u, 2)), 12))
print("and they factor exactly as")
print("  P_port = G^2 [1 + mu (G-1)^2] + O(u^3).")


# ---------------------------------------------------------------------------
# 8) Provisional conclusion
# ---------------------------------------------------------------------------

banner("8) Provisional conclusion from the DtN reduction prototype")
print("- The low-frequency monopole DtN data naturally split into a sound-speed invariant")
print("      C_s = [Z4/Z2^3]/[Z4/Z2^3]_0")
print("  and a geometry invariant")
print("      G   = [Z4/Z2^2]/[Z4/Z2^2]_0.")
print("- In the cylinder unit-test branch, these are exactly")
print("      C_s = c_s^2/c_{s0}^2 = 1 - 4u,")
print("      G   = L/L0 = a(u).")
print("- The minimal conservative local response built from G that preserves the frozen 1PN slot")
print("  must start quadratically in (G-1).")
print("- Matching the missing one-body 2PN U^2 v^2 coefficient then fixes uniquely")
print("      D_eff(u) = C_s(u) [1 + mu (G(u)-1)^2],")
print("      mu = 8/g1^2 = 32768/3249,   g1 = 57/64.")
print("- Under the current throat closure, this reproduces the exact isotropic test-mass target through 2PN.")
print("- The earlier raw resonance proxy / port-normalization fit is exactly the same closure written as")
print("      D_eff = D_raw P_port,   P_port = G^2 [1 + mu (G-1)^2].")
