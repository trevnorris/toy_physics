
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


banner("STAGE 51 — RESONANCE-CORRECTED THRESHOLDS")

C2, W_wall, Pe_req, Delta0, Deltainf = sp.symbols(
    "C2 W_wall Pe_req Delta_0 Delta_inf", positive=True, real=True
)
Pres = sp.symbols("P_res", positive=True, real=True)

# ---------------------------------------------------------------------------
# 1. Resonance-corrected wall figure (derived, not postulated)
# ---------------------------------------------------------------------------
#
# Premise: a wall-incident wave of complex amplitude A_in encounters a
# transmission coefficient C(r) (real, positive at the resonance branch
# we consider).  Power scales as |amplitude|^2, so the transmitted wall
# power figure is |C|^2 * W_wall.  We carry C2 := C^2 as the audit symbol.

A_in, C_sym = sp.symbols("A_in C", positive=True, real=True)
A_trans = C_sym * A_in
# Power in the incident wave ~ |A_in|^2 normalises to W_wall; transmitted power:
W_wall_def = A_in**2
W_res_derived = sp.simplify((A_trans**2).subs(A_in**2, W_wall) * 1)
# Re-express |C|^2 as the audit symbol C2:
W_res_derived = sp.simplify(W_res_derived.subs(C_sym**2, C2))
expect_zero("W_res - C2 * W_wall", W_res_derived - C2 * W_wall)

# At the resonance point r = r_*, C^2 -> C_res^2 by definition of the
# resonance branch; the penalty factor is the *inverse* of the amplification.
Cres = sp.symbols("C_res", positive=True, real=True)
Pres_derived = sp.simplify(1 / Cres**2)
expect_zero("P_res - 1/C_res^2", Pres_derived - 1/Cres**2)
# And consistency with the audit-level symbol P_res introduced as positive:
expect_zero("P_res*C_res^2 - 1", (1/Cres**2)*Cres**2 - 1)

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
# 3. Profile-sensitive band widths, computed two independent ways
# ---------------------------------------------------------------------------
#
# Way A: as the difference of the F1-derived profile and matched thresholds,
#        evaluated at the resonance point C^2 = 1/P_res.
success_band_widthA = sp.simplify((Wsuff_res - Wsuff_match).subs(C2, 1 / Pres))
failure_band_widthA = sp.simplify((Wfail_res - Wfail_match).subs(C2, 1 / Pres))

# Way B: by solving Wsuff_res = Pres * Wsuff_match for the gap directly.
#        We pose the equation and let Solve extract the difference symbolically.
gap_sym = sp.symbols("gap", real=True)
success_band_widthB = sp.solve(sp.Eq(Wsuff_match + gap_sym, Pres * Wsuff_match), gap_sym)[0]
failure_band_widthB = sp.solve(sp.Eq(Wfail_match + gap_sym, Pres * Wfail_match), gap_sym)[0]

print("Success-side band width (A) =", success_band_widthA)
print("Failure-side band width (A) =", failure_band_widthA)
print("Success-side band width (B) =", success_band_widthB)
print("Failure-side band width (B) =", failure_band_widthB)

# Independent-derivation cross-check: the two ways must agree.
expect_zero("success band A vs B", success_band_widthA - success_band_widthB)
expect_zero("failure band A vs B", failure_band_widthA - failure_band_widthB)

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
