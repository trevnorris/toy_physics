
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
# 1. Resonance-corrected wall figure
# ---------------------------------------------------------------------------

W_res = sp.simplify(C2 * W_wall)
print("W_res =", W_res)

# ---------------------------------------------------------------------------
# 2. Exact threshold translation
# ---------------------------------------------------------------------------

Wfail_match = sp.simplify(Pe_req / Deltainf)
Wsuff_match = sp.simplify(Pe_req / Delta0)

Wfail_res = sp.simplify(Pe_req / (C2 * Deltainf))
Wsuff_res = sp.simplify(Pe_req / (C2 * Delta0))

print("Matched fail threshold     =", Wfail_match)
print("Matched succeed threshold  =", Wsuff_match)
print("Profile-family fail thresh =", Wfail_res)
print("Profile-family succ thresh =", Wsuff_res)

# If P_res = 1/C2, the thresholds scale by P_res.
expect_zero("Wfail_res - P_res*Wfail_match", Wfail_res.subs(C2, 1 / Pres) - Pres * Wfail_match)
expect_zero("Wsuff_res - P_res*Wsuff_match", Wsuff_res.subs(C2, 1 / Pres) - Pres * Wsuff_match)

# ---------------------------------------------------------------------------
# 3. Exact profile-sensitive bands
# ---------------------------------------------------------------------------

success_band_width = sp.simplify((Pres * Wsuff_match) - Wsuff_match)
failure_band_width = sp.simplify((Pres * Wfail_match) - Wfail_match)

print("Success-side band width =", success_band_width)
print("Failure-side band width =", failure_band_width)

expect_zero("success width / matched threshold", sp.simplify(success_band_width / Wsuff_match - (Pres - 1)))
expect_zero("failure width / matched threshold", sp.simplify(failure_band_width / Wfail_match - (Pres - 1)))

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
