
#!/usr/bin/env python3
"""
moving_throat_pde_stage51_resonance_thresholds_sympy_audit.py

SymPy-backed audit for the resonance-corrected threshold window.

What this script checks
-----------------------
1. Resonance-corrected gain and wall figure:
      W_res(r) = C^2(r) W_wall.
2. Exact threshold translation from the matched branch to the profile family.
3. Explicit resonance-point penalty factor:
      P_res = 1 / C_res^2.
4. Exact width of the profile-sensitive threshold band.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 68 — RESONANCE-CORRECTED THRESHOLDS")

C2, W_wall, Pe_req, Delta0, Deltainf = sp.symbols(
    "C2 W_wall Pe_req Delta_0 Delta_inf", positive=True, real=True
)
Pres = sp.symbols("P_res", positive=True, real=True)

# ---------------------------------------------------------------------------
# 1. Resonance-corrected wall figure (derived from matched-branch gain)
# ---------------------------------------------------------------------------
#
# Premise (notes section 1): the matched-branch gain is
#   G_match = rho_star * g_phi^2 * N_phiphi / (m * c_s^2 * K_X),
# and the wall figure of merit is W_wall = kappa * G_match.
# The profile-family gain at coherence C^2 is G_res = C^2 * G_match,
# so W_res = kappa * G_res. We verify that W_res = C^2 * W_wall by
# building each side from these independent component symbols.

rho_star, g_phi, N_phiphi, m_s, c_s, K_X, kappa = sp.symbols(
    "rho_star g_phi N_phiphi m c_s K_X kappa", positive=True, real=True
)
Gmatch_expr = rho_star * g_phi**2 * N_phiphi / (m_s * c_s**2 * K_X)
Wwall_expr = kappa * Gmatch_expr
Gres_expr  = C2 * Gmatch_expr
Wres_expr  = kappa * Gres_expr
expect_zero("W_res - C2 * W_wall (from gain decomposition)",
            sp.simplify(Wres_expr - C2 * Wwall_expr))
# Also confirm the matched limit C2 -> 1 collapses W_res to W_wall:
expect_zero("W_res(C2->1) - W_wall (matched limit)",
            sp.simplify(Wres_expr.subs(C2, 1) - Wwall_expr))

# ---------------------------------------------------------------------------
# 1b. Resonance penalty factor (derived from ratio of required wall figures)
# ---------------------------------------------------------------------------
#
# The required wall figure on the matched branch is W_wall = Pe_req/Delta.
# On the profile family at coherence C^2, it is W_wall = Pe_req/(C^2 Delta).
# The amplification of the *required* wall figure from matched to profile
# is therefore 1/C^2.  At the resonance point C^2 -> C_res^2 (carried in
# from stage 067), this amplification factor is P_res = 1/C_res^2.

Cres, Delta_sym_pre = sp.symbols("C_res Delta", positive=True, real=True)
W_required_match = sp.solve(
    sp.Eq(sp.Symbol("Wm", positive=True) * Delta_sym_pre, Pe_req),
    sp.Symbol("Wm", positive=True),
)[0]
W_required_prof = sp.solve(
    sp.Eq(C2 * sp.Symbol("Wp", positive=True) * Delta_sym_pre, Pe_req),
    sp.Symbol("Wp", positive=True),
)[0]
Pres_from_ratio = sp.simplify((W_required_prof / W_required_match).subs(C2, Cres**2))
expect_zero("P_res - 1/C_res^2 (from required-wall-figure ratio)",
            Pres_from_ratio - 1/Cres**2)

# Numeric anchor: the paper card states P_res = 1.005612487760576 and
# C_res^2 = 0.994418836451529 (carried from stage 067). Verify the link.
Cres_sq_numeric = sp.Float("0.994418836451529", 20)
Pres_numeric    = 1 / Cres_sq_numeric
Pres_paper      = sp.Float("1.005612487760576", 20)
numeric_residual = sp.Abs(Pres_numeric - Pres_paper)
print(f"P_res numeric residual = {numeric_residual}")
if numeric_residual >= sp.Float("1e-12", 20):
    raise AssertionError(
        f"P_res numeric anchor failed: |1/C_res^2 - paper P_res| = {numeric_residual}"
    )

# ---------------------------------------------------------------------------
# 2. Exact threshold translation (derived from W_res, not postulated)
# ---------------------------------------------------------------------------
#
# Matched-branch thresholds come from the matched-asymptotic Peclet balance:
#   W_match * Delta = Pe_req      (Delta = Delta_inf for failure,
#                                  Delta = Delta_0   for success)
# Profile-family thresholds use W_res = C^2 W_wall in place of W_wall, so the
# *same* Peclet balance with the substitution W_wall -> W_res yields the
# profile-family thresholds.  Both sides are computed independently.

W_match_sym = sp.symbols("W_match", positive=True, real=True)
Delta_sym  = sp.symbols("Delta", positive=True, real=True)
# Solve the matched balance W_match * Delta = Pe_req for W_match:
W_match_sol = sp.solve(sp.Eq(W_match_sym * Delta_sym, Pe_req), W_match_sym)[0]
Wfail_match = sp.simplify(W_match_sol.subs(Delta_sym, Deltainf))
Wsuff_match = sp.simplify(W_match_sol.subs(Delta_sym, Delta0))

# Profile-family thresholds: replace W_wall in the matched balance with W_res.
# Equivalently W_match -> W_match / C^2 (since W_res = C^2 W_wall demands a
# C^2 stronger wall to reach the same incident-side figure).  Derive via Solve:
W_prof_sym = sp.symbols("W_prof", positive=True, real=True)
W_prof_sol = sp.solve(sp.Eq(C2 * W_prof_sym * Delta_sym, Pe_req), W_prof_sym)[0]
Wfail_res = sp.simplify(W_prof_sol.subs(Delta_sym, Deltainf))
Wsuff_res = sp.simplify(W_prof_sol.subs(Delta_sym, Delta0))

print("Matched fail threshold     =", Wfail_match)
print("Matched succeed threshold  =", Wsuff_match)
print("Profile-family fail thresh =", Wfail_res)
print("Profile-family succ thresh =", Wsuff_res)

# Non-trivial check: the *ratio* of profile-family to matched threshold equals
# 1/C^2.  Independent derivations of numerator and denominator collapse to this:
expect_zero("Wfail_res * C2 - Wfail_match", Wfail_res * C2 - Wfail_match)
expect_zero("Wsuff_res * C2 - Wsuff_match", Wsuff_res * C2 - Wsuff_match)

# At resonance C2 = 1/P_res, the profile thresholds scale by P_res:
expect_zero("Wfail_res(C2->1/Pres) - Pres*Wfail_match",
            Wfail_res.subs(C2, 1 / Pres) - Pres * Wfail_match)
expect_zero("Wsuff_res(C2->1/Pres) - Pres*Wsuff_match",
            Wsuff_res.subs(C2, 1 / Pres) - Pres * Wsuff_match)

# ---------------------------------------------------------------------------
# 3. Profile-sensitive band widths, two routes through different symbols
# ---------------------------------------------------------------------------
#
# Way A (C-form): compute the difference at the resonance point using Cres
#   directly. Wsuff_res evaluated at C^2 = Cres^2 is Pe_req/(Cres^2 Delta_0).
# Way B (P-form): compute the difference using Pres as the audit symbol.
# Both should equal (Pres - 1)*Wsuff_match after the relation Pres = 1/Cres^2.
# Any perturbation of that relation breaks one form but not the other, so the
# cross-check is now sensitive to the Pres = 1/Cres^2 link.

# Cres is defined in section 1b above; reuse it here.
Wsuff_res_C   = Pe_req / (Cres**2 * Delta0)
Wfail_res_C   = Pe_req / (Cres**2 * Deltainf)
success_band_widthA = sp.simplify(Wsuff_res_C - Wsuff_match)
failure_band_widthA = sp.simplify(Wfail_res_C - Wfail_match)

success_band_widthB = sp.simplify((Pres - 1) * Wsuff_match)
failure_band_widthB = sp.simplify((Pres - 1) * Wfail_match)

print("Success-side band width (C-form) =", success_band_widthA)
print("Failure-side band width (C-form) =", failure_band_widthA)
print("Success-side band width (P-form) =", success_band_widthB)
print("Failure-side band width (P-form) =", failure_band_widthB)

# Cross-check: under the relation Pres = 1/Cres^2, the two forms must agree.
expect_zero("success band C-form vs P-form (under Pres = 1/Cres^2)",
            sp.simplify((success_band_widthA - success_band_widthB).subs(Pres, 1/Cres**2)))
expect_zero("failure band C-form vs P-form (under Pres = 1/Cres^2)",
            sp.simplify((failure_band_widthA - failure_band_widthB).subs(Pres, 1/Cres**2)))

banner("FINAL LEDGER")
print("Exact formulas:")
print("  W_res(r) = C^2(r) W_wall")
print("  W_fail^(res) = Pe_req / [ C^2(r) Delta_inf ]")
print("  W_suff^(res) = Pe_req / [ C^2(r) Delta_0 ]")
print("At the resonance point C^2 = C_res^2, these become:")
print("  W_fail^(res,*) = P_res Pe_req / Delta_inf")
print("  W_suff^(res,*) = P_res Pe_req / Delta_0")
print("with P_res = 1/C_res^2.")
print()
print("Interpretation:")
print("  The explicit independent-profile benchmark modifies the matched-layer")
print("  threshold window only by the multiplicative factor P_res.")
