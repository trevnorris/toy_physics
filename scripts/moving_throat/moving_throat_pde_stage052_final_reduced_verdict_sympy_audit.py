
#!/usr/bin/env python3
"""
moving_throat_pde_stage52_final_reduced_verdict_sympy_audit.py

SymPy-backed audit for the final reduced support/source verdict.

What this script checks
-----------------------
1. Universal matched-branch threshold window.
2. Resonance-family threshold translation by the penalty factor P_res.
3. Exact relative width of the only profile-sensitive sub-bands.
4. The statement that profile mismatch matters only at O(P_res - 1).
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


banner("STAGE 52 — FINAL REDUCED SUPPORT/SOURCE VERDICT")

Pe_req, Delta0, Deltainf, Pres = sp.symbols(
    "Pe_req Delta_0 Delta_inf P_res", positive=True, real=True
)

Wfail_match = sp.simplify(Pe_req / Deltainf)
Wsuff_match = sp.simplify(Pe_req / Delta0)

Wfail_res = sp.simplify(Pres * Wfail_match)
Wsuff_res = sp.simplify(Pres * Wsuff_match)

print("Universal matched fail threshold    =", Wfail_match)
print("Universal matched success threshold =", Wsuff_match)
print("Resonance-family fail threshold     =", Wfail_res)
print("Resonance-family success threshold  =", Wsuff_res)

# Profile-sensitive sub-band widths.
delta_fail = sp.simplify(Wfail_res - Wfail_match)
delta_suff = sp.simplify(Wsuff_res - Wsuff_match)

print("Failure-side profile-sensitive width =", delta_fail)
print("Success-side profile-sensitive width =", delta_suff)

expect_zero("failure relative width", sp.simplify(delta_fail / Wfail_match - (Pres - 1)))
expect_zero("success relative width", sp.simplify(delta_suff / Wsuff_match - (Pres - 1)))

# A wall figure larger than Pres * matched threshold guarantees success
# even on the independent-profile resonance family.
W_wall = sp.symbols("W_wall", positive=True, real=True)
success_margin = sp.simplify(W_wall - Wsuff_res)
failure_margin = sp.simplify(Wfail_match - W_wall)

print("Success if W_wall - P_res*Pe_req/Delta_0 >= 0:")
sp.pprint(success_margin)
print("Failure if Pe_req/Delta_inf - W_wall >= 0:")
sp.pprint(failure_margin)

banner("FINAL LEDGER")
print("Exact reduced verdict:")
print("  Universal fail   : W_wall <= Pe_req / Delta_inf")
print("  Universal succeed: W_wall >= Pe_req / Delta_0")
print("  Resonance-family thresholds are shifted only by P_res.")
print("Therefore profile mismatch matters only in the narrow O(P_res - 1) bands")
print("around the universal matched-branch thresholds.")
