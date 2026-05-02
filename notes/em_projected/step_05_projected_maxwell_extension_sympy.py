#!/usr/bin/env python3
"""Extend the projected-Maxwell derivation to include a weighted gauge-fixing profile H(w).

This script does three jobs:
1. writes the exact projected brane equation for a generalized gauge-fixing weight H(w),
2. derives the zero-mode projection-first effective coupling and effective gauge parameter,
3. compares the key cases H=1 and H=Z, including Gaussian matched-kernel examples.

The main reduced formulas are
    mu0_eff^(proj) = mu0 * I_WS / I_WZ
    1/xi_eff^(proj) = I_WH / (xi * I_WZ)
with
    I_WZ = ∫ W Z dw,
    I_WH = ∫ W H dw,
    I_WS = ∫ W S dw.

For H = Z, xi_eff^(proj) = xi exactly for any normalized projection kernel W.
"""
from __future__ import annotations

import sympy as sp


def line(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.simplify(expr)
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


# -----------------------------------------------------------------------------
# Symbols and helpers
# -----------------------------------------------------------------------------

t, x, y, z, w = sp.symbols("t x y z w", real=True)
lam, xi, mu0, R = sp.symbols("lambda xi mu0 R", positive=True, finite=True, nonzero=True)
Zint, IWZ = sp.symbols("Z_int I_WZ", positive=True, nonzero=True)


# -----------------------------------------------------------------------------
# 1) Exact projected law with weighted gauge fixing
# -----------------------------------------------------------------------------
line("1) Exact projected inhomogeneous law with a weighted gauge-fixing profile H(w)")

print("Start from the weighted bulk equation")
print("  ∂_M(Z F^{MN}) + (1/ξ) ∂^N(H B) = μ0 J^N")
print("with B := ∂·A and H(w) an arbitrary gauge-fixing weight.")
print()
print("For a brane component ν ∈ {0,1,2,3}, projection with W(w) gives")
print("  ∂_μ Proj_W[Z F^{μν}] + Boundary[W Z F^{wν}] - Proj_{W'}[Z F^{wν}] + (1/ξ) ∂^ν Proj_W[H B] = μ0 Proj_W[J^ν]")
print()
print("Under boundary decay / compact support in w:")
print("  ∂_μ Proj_W[Z F^{μν}] - Proj_{W'}[Z F^{wν}] + (1/ξ) ∂^ν Proj_W[H B] = μ0 Proj_W[J^ν]")
print()
print("So H(w) changes only the projected gauge-driver term.")
print("The homogeneous laws and the transverse leakage term stay the same.")


# -----------------------------------------------------------------------------
# 2) Zero-mode projection-first effective law
# -----------------------------------------------------------------------------
line("2) Zero-mode projection-first effective law")

I_WZ, I_WH, I_WS = sp.symbols("I_WZ I_WH I_WS", positive=True, nonzero=True)
mu0_eff_proj = sp.simplify(mu0 * I_WS / I_WZ)
xi_eff_proj = sp.simplify(xi * I_WZ / I_WH)
inv_xi_eff_proj = sp.simplify(I_WH / (xi * I_WZ))
assert_zero("effective gauge inverse is reciprocal", inv_xi_eff_proj * xi_eff_proj - 1)

print("Assume the zero-mode / far-field brane ansatz")
print("  A_μ(x,w) = a_μ(x),   A_w = 0,   ∂_w A_μ = 0,   F^{wν} = 0,   J^ν(x,w) = j^ν(x) S(w).")
print("Then B = ∂·a is independent of w, and the projected equation becomes")
print("  I_WZ ∂_μ f^{μν} + (I_WH/ξ) ∂^ν(∂·a) = μ0 I_WS j^ν")
print("with")
print("  I_WZ := ∫ W(w) Z(w) dw")
print("  I_WH := ∫ W(w) H(w) dw")
print("  I_WS := ∫ W(w) S(w) dw")
print()
print("After dividing by I_WZ:")
print("  ∂_μ f^{μν} + (1/xi_eff_proj) ∂^ν(∂·a) = mu0_eff_proj j^ν")
print()
print("where")
print("  mu0_eff_proj  =", sp.sstr(mu0_eff_proj))
print("  1/xi_eff_proj =", sp.sstr(inv_xi_eff_proj))
print("  xi_eff_proj   =", sp.sstr(xi_eff_proj))
print()
print("This shows that projection-first produces BOTH an effective coupling and an effective gauge parameter.")


# -----------------------------------------------------------------------------
# 3) General H=1 versus H=Z comparison
# -----------------------------------------------------------------------------
line("3) General comparison of H=1 versus H=Z")

# normalized W: ∫W = 1
xi_eff_from_IWH = xi * I_WZ / I_WH
xi_eff_H1 = sp.simplify(xi_eff_from_IWH.subs(I_WH, 1))
xi_eff_HZ = sp.simplify(xi_eff_from_IWH.subs(I_WH, I_WZ))

print("Use the normalization ∫W(w) dw = 1.")
print()
print("Case H = 1 (unweighted gauge fixing):")
print("  I_WH = 1")
print("  xi_eff_proj(H=1) =", sp.sstr(xi_eff_H1))
print("So the projected gauge parameter depends on the observer kernel through I_WZ = ∫WZ.")
print()
print("Case H = Z (localized gauge fixing):")
print("  I_WH = I_WZ")
print("  xi_eff_proj(H=Z) =", sp.sstr(xi_eff_HZ))
print("So H=Z preserves the same gauge parameter ξ in the zero-mode projected equation for ANY normalized W.")
print()
print("Important: H affects the gauge-driver sector only. The effective coupling μ0_eff_proj depends on the source profile S(w), not on H(w).")

W_generic = sp.Function("W")(w)
Z_generic = sp.Function("Z")(w)
H_generic = sp.Function("H")(w)
B_generic = sp.Function("B")(t, x, y, z)
I_WZ_generic = sp.Integral(W_generic * Z_generic, (w, -sp.oo, sp.oo))
I_WH_generic = sp.Integral(W_generic * H_generic, (w, -sp.oo, sp.oo))
assert_zero(
    "H=Z integral identification before cancellation",
    I_WH_generic.subs(H_generic, Z_generic) - I_WZ_generic,
)

gauge_driver_projected = sp.Integral(W_generic * H_generic * B_generic, (w, -sp.oo, sp.oo)) / xi
gauge_driver_reduced = B_generic * I_WZ_generic / xi
assert_zero(
    "zero-mode H=Z gauge-driver projection",
    gauge_driver_projected.subs(H_generic, Z_generic) - gauge_driver_reduced,
)
assert_zero(
    "H=Z effective gauge via symbolic substitution",
    sp.simplify((xi * I_WZ_generic / I_WH_generic).subs(H_generic, Z_generic)) - xi,
)


# -----------------------------------------------------------------------------
# 4) Source-profile comparison at the zero-mode level
# -----------------------------------------------------------------------------
line("4) Source-profile comparison: when projection-first matches reduction-first")

# Algebraic matching channel: S = Z/Z_int
mu0_eff_match_source = sp.simplify(mu0 * (I_WZ / Zint) / I_WZ)
assert_zero("matched source coupling", mu0_eff_match_source - mu0 / Zint)

print("If the source profile is the normalized localization profile")
print("  S(w) = Z(w) / Z_int,")
print("then I_WS = I_WZ / Z_int and")
print("  mu0_eff_proj =", sp.sstr(mu0_eff_match_source))
print()
print("So projection-first reproduces the reduction-first Maxwell coupling exactly:")
print("  mu0_eff_proj = mu0 / Z_int")
print("for ANY normalized projection kernel W.")
print()
print("Together with H = Z, the zero-mode projected equation becomes")
print("  ∂_μ f^{μν} + (1/ξ) ∂^ν(∂·a) = (μ0 / Z_int) j^ν")
print("which is exactly the same reduced gauge-fixed Maxwell form as the localized-gauge reduction.")
print()
print("So there is an exact matching channel:")
print("  source profile S = Z/Z_int   and   gauge-fixing weight H = Z.")


# -----------------------------------------------------------------------------
# 5) Explicit Gaussian example with a matched observer kernel
# -----------------------------------------------------------------------------
line("5) Explicit Gaussian example with matched observer kernel W = Z/Z_int")

Z = sp.exp(-w**2 / lam**2)
Z_int_gauss = sp.simplify(sp.integrate(Z, (w, -sp.oo, sp.oo)))
Z2_int_gauss = sp.simplify(sp.integrate(Z**2, (w, -sp.oo, sp.oo)))
W_match = sp.simplify(Z / Z_int_gauss)
I_WZ_match = sp.simplify(sp.integrate(W_match * Z, (w, -sp.oo, sp.oo)))
I_WH_H1_match = sp.simplify(sp.integrate(W_match, (w, -sp.oo, sp.oo)))
I_WH_HZ_match = sp.simplify(sp.integrate(W_match * Z, (w, -sp.oo, sp.oo)))

# delta-localized source sampled by the projection kernel
I_WS_delta_match = sp.simplify(W_match.subs(w, 0))
mu0_eff_delta_match = sp.simplify(mu0 * I_WS_delta_match / I_WZ_match)

# matched distributed source S = Z/Z_int
I_WS_source_match = sp.simplify(sp.integrate(W_match * (Z / Z_int_gauss), (w, -sp.oo, sp.oo)))
mu0_eff_source_match = sp.simplify(mu0 * I_WS_source_match / I_WZ_match)

xi_eff_H1_match = sp.simplify(xi * I_WZ_match / I_WH_H1_match)
xi_eff_HZ_match = sp.simplify(xi * I_WZ_match / I_WH_HZ_match)
assert_zero("Gaussian Z_int", Z_int_gauss - sp.sqrt(sp.pi) * lam)
assert_zero("Gaussian Z2_int", Z2_int_gauss - sp.sqrt(2 * sp.pi) * lam / 2)
assert_zero("matched Gaussian I_WZ", I_WZ_match - sp.sqrt(2) / 2)
assert_zero("Gaussian H=Z gauge parameter", xi_eff_HZ_match - xi)
assert_zero("Gaussian H=1 gauge parameter", xi_eff_H1_match - xi / sp.sqrt(2))
assert_zero("Gaussian matched source coupling", mu0_eff_source_match - mu0 / Z_int_gauss)
assert_zero("Gaussian delta-source mismatch ratio", mu0_eff_delta_match / (mu0 / Z_int_gauss) - sp.sqrt(2))

print("Take Z(w) = exp(-w^2/λ^2) and the matched observer kernel W = Z/Z_int.")
print("  Z_int  =", sp.sstr(Z_int_gauss))
print("  Z2_int =", sp.sstr(Z2_int_gauss))
print("  I_WZ   =", sp.sstr(I_WZ_match), "= Z2_int / Z_int")
print()
print("Gauge sector:")
print("  xi_eff_proj(H=1) =", sp.sstr(xi_eff_H1_match))
print("  xi_eff_proj(H=Z) =", sp.sstr(xi_eff_HZ_match))
print()
print("Source sector, case A: delta-localized source S = δ(w)")
print("  I_WS = W(0) =", sp.sstr(I_WS_delta_match))
print("  mu0_eff_proj =", sp.sstr(mu0_eff_delta_match))
print("  ratio to reduction-first mu0/Z_int =", sp.sstr(sp.simplify(mu0_eff_delta_match / (mu0 / Z_int_gauss))))
print()
print("Source sector, case B: matched distributed source S = Z/Z_int")
print("  I_WS =", sp.sstr(I_WS_source_match))
print("  mu0_eff_proj =", sp.sstr(mu0_eff_source_match))
print("  ratio to reduction-first mu0/Z_int =", sp.sstr(sp.simplify(mu0_eff_source_match / (mu0 / Z_int_gauss))))


# -----------------------------------------------------------------------------
# 6) Comparison with reduction-first gauge fixing
# -----------------------------------------------------------------------------
line("6) Comparison with reduction-first gauge fixing")

xi4_general = sp.simplify(xi * Zint / sp.Symbol("H_int", positive=True))
xi4_unweighted_reg = sp.simplify(xi * Z_int_gauss / (2 * R))
xi4_unweighted_limit = sp.limit(xi4_unweighted_reg, R, sp.oo)
assert_zero("reduction-first H=1 regulator limit", xi4_unweighted_limit)

print("Localized action reduction with a finite gauge-fixing weight H(w) gives")
print("  xi4 = xi * Z_int / H_int")
print("with H_int = ∫H(w)dw.")
print()
print("So for H = Z:")
print("  xi4 = xi")
print()
print("For H = 1 on a finite regulator interval [-R,R]:")
print("  xi4(R) =", sp.sstr(xi4_unweighted_reg))
print("  lim_{R→∞} xi4(R) =", sp.sstr(xi4_unweighted_limit))
print()
print("Interpretation:")
print("  - reduction-first with H=1 is unsafe as a noncompact zero-mode gauge-fixed action unless Lorenz gauge is imposed before reduction,")
print("  - projection-first with H=1 stays finite because W is normalized, but it produces an observer-kernel-dependent xi_eff_proj,")
print("  - H=Z aligns the projected and reduced gauge sectors cleanly.")


# -----------------------------------------------------------------------------
# 7) Final summary
# -----------------------------------------------------------------------------
line("7) Final summary")
print("Exact weighted projected law:")
print("  ∂_μ Proj_W[Z F^{μν}] - Proj_{W'}[Z F^{wν}] + (1/ξ) ∂^ν Proj_W[H B] = μ0 Proj_W[J^ν]")
print()
print("Zero-mode projected effective parameters:")
print("  mu0_eff_proj = mu0 (∫WS)/(∫WZ)")
print("  xi_eff_proj  = xi (∫WZ)/(∫WH)")
print()
print("Key structural results:")
print("  1. H = Z gives xi_eff_proj = xi for any normalized W.")
print("  2. S = Z/Z_int gives mu0_eff_proj = mu0/Z_int for any normalized W.")
print("  3. With BOTH H = Z and S = Z/Z_int, projection-first reproduces the same zero-mode gauge-fixed Maxwell law as reduction-first.")
print("  4. With a delta-localized source instead, the coupling mismatch remains, even if H = Z.")
print("STATUS: PASS")
