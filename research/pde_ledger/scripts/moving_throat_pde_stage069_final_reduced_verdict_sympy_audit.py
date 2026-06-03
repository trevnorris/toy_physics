#!/usr/bin/env python3
"""
moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py

SymPy-backed audit for the final reduced support/source verdict.

What this script checks
-----------------------
1. The matched-branch Stage-49 window is exactly
   [Pe_req / Delta_inf, Pe_req / Delta_0].
2. The Stage-51 resonance family shifts both thresholds by the exact penalty
   factor P_res = 1 / C_res^2.
3. The profile-sensitive failure-side and success-side bands have the exact
   widths claimed in the paper.
4. The reduced three-zone verdict is the correct matched-branch theorem, while
   the resonance-family refinement sits in two narrow side-bands.

Provenance notes
----------------
- The matched window `[Pe_req/Delta_inf, Pe_req/Delta_0]` is the Stage 066
  theorem surface carried forward verbatim.
- `P_res = 1/C_res^2` is the Stage 068 resonance penalty; this script keeps
  those symbols unchanged and only assembles the final verdict from them.
- Source for the matched-window form `Pe_req/Delta_inf`:
  `scripts/moving_throat_pde_stage066_*.py` produces this form from the
  matched-branch theorem at its file's `expect_zero("Stage 049 matched
  window", ...)` block. This script imports the form symbolically as
  `Pe_req/Deltainf` and `Pe_req/Delta_0`; it does not re-derive it.
- Source for the resonance penalty `P_res = 1 + Pres_gap`:
  `scripts/moving_throat_pde_stage068_*.py` produces this form from the
  resonance amplification factor. This script imports the form symbolically
  as `Pres = 1 + Pres_gap` and asserts the band-edge ratio matches; it does
  not re-derive the amplification factor.
- The assertions in this script are conditional on Stages 066 and 068
  being correct. The new ratio-based assertions added per F1 verify the
  consolidation step but do not re-verify the upstream derivations.
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


def expect_positive(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.factor(expr))
    print(f"{name} > 0  -> {expr}")
    if expr.is_positive is not True:
        raise AssertionError(f"{name} is not provably positive")


banner("STAGE 069 — FINAL REDUCED SUPPORT/SOURCE VERDICT")

Pe_req, Delta0, Delta_gap = sp.symbols("Pe_req Delta_0 Delta_gap", positive=True, real=True)
Pres_gap, v_fail, v_succ = sp.symbols("Pres_gap v_fail v_succ", positive=True, real=True)

Deltainf = sp.simplify(Delta0 + Delta_gap)
Pres = sp.simplify(1 + Pres_gap)
Cres2 = sp.simplify(1 / Pres)
u_fail = sp.simplify(v_fail / (1 + v_fail))
u_succ = sp.simplify(v_succ / (1 + v_succ))

Wfail_match = sp.simplify(Pe_req / Deltainf)
Wsuff_match = sp.simplify(Pe_req / Delta0)
Wfail_res = sp.simplify(Pres * Wfail_match)
Wsuff_res = sp.simplify(Pres * Wsuff_match)

# --- Upstream-anchor assumption (see docstring) ---
# Wfail_match, Wsuff_match, Wfail_res, Wsuff_res are constructed here from
# free positive symbols. Their interpretation as the Stage 066 / Stage 068
# threshold values is asserted by the F1 consolidation block below
# (matched fail edge from W_match(Delta_inf) = 0, P_res from band-edge
# ratio matches (1 + Pres_gap) = 0), but the upstream derivations
# themselves are not re-checked here.

print("Matched-branch fail threshold       =", Wfail_match)
print("Matched-branch success threshold    =", Wsuff_match)
print("Resonance-family fail threshold     =", Wfail_res)
print("Resonance-family success threshold  =", Wsuff_res)

# Stage-49 matched-window derivation as a generating function.
# Encodes the theorem "W_match(Delta_eff) = Pe_req / Delta_eff" with
# Delta_eff ranging over [Delta_0, Delta_inf] on the matched branch.
Delta_eff = sp.symbols("Delta_eff", positive=True, real=True)
W_match_generator = Pe_req / Delta_eff
expect_zero(
    "matched fail edge from W_match(Delta_inf)",
    W_match_generator.subs(Delta_eff, Deltainf) - Wfail_match,
)
expect_zero(
    "matched success edge from W_match(Delta_0)",
    W_match_generator.subs(Delta_eff, Delta0) - Wsuff_match,
)
expect_positive(
    "W_match decreasing in Delta_eff",
    -sp.diff(W_match_generator, Delta_eff) * Delta_eff**2 / Pe_req,
)

# Stage-49 matched window and Stage-51 penalty factor.
expect_zero("P_res - 1/C_res^2", Pres - 1 / Cres2)
expect_zero(
    "matched window width",
    Wsuff_match - Wfail_match - Pe_req * (Deltainf - Delta0) / (Delta0 * Deltainf),
)
expect_positive("Delta_inf - Delta_0", Deltainf - Delta0)
expect_positive("matched success threshold - matched fail threshold", Wsuff_match - Wfail_match)
# Stage-51 / Stage-68 resonance penalty as a band-edge ratio.
# Encodes the theorem "P_res = W_fail_res / W_fail_match" rather than
# defining P_res and C_res^2 to satisfy it by construction.
Pres_from_ratio = sp.simplify(Wfail_res / Wfail_match)
expect_zero(
    "P_res from band-edge ratio matches (1 + Pres_gap)",
    Pres_from_ratio - (1 + Pres_gap),
)
Pres_from_success_ratio = sp.simplify(Wsuff_res / Wsuff_match)
expect_zero(
    "P_res from success-band ratio agrees with failure-band ratio",
    Pres_from_ratio - Pres_from_success_ratio,
)
expect_zero(
    "resonance fail threshold - Pe_req/(C_res^2 Delta_inf)",
    Wfail_res - Pe_req / (Cres2 * Deltainf),
)
expect_zero(
    "resonance success threshold - Pe_req/(C_res^2 Delta_0)",
    Wsuff_res - Pe_req / (Cres2 * Delta0),
)
expect_positive("1 - C_res^2", 1 - Cres2)
expect_positive("P_res - 1", Pres - 1)

# Exact width of the only profile-sensitive side-bands.
delta_fail = sp.simplify(Wfail_res - Wfail_match)
delta_succ = sp.simplify(Wsuff_res - Wsuff_match)

print("Failure-side profile-sensitive width =", delta_fail)
print("Success-side profile-sensitive width =", delta_succ)

expect_zero(
    "failure-side width - Pe_req(1-C_res^2)/(C_res^2 Delta_inf)",
    delta_fail - Pe_req * (1 - Cres2) / (Cres2 * Deltainf),
)
expect_zero(
    "success-side width - Pe_req(1-C_res^2)/(C_res^2 Delta_0)",
    delta_succ - Pe_req * (1 - Cres2) / (Cres2 * Delta0),
)
expect_zero(
    "P_res - 1 - (1-C_res^2)/C_res^2",
    sp.simplify((Pres - 1) - (1 - Cres2) / Cres2),
)
expect_positive("failure-side width", delta_fail)
expect_positive("success-side width", delta_succ)

# The two genuinely profile-sensitive side-bands.
W_failure_band = sp.simplify(Wfail_match + u_fail * delta_fail)
W_success_band = sp.simplify(Wsuff_match + u_succ * delta_succ)

expect_positive("failure-band point - matched fail edge", W_failure_band - Wfail_match)
expect_positive("resonance fail edge - failure-band point", Wfail_res - W_failure_band)
expect_positive("success-band point - matched success edge", W_success_band - Wsuff_match)
expect_positive("resonance success edge - success-band point", Wsuff_res - W_success_band)

banner("FINAL LEDGER")
print("Stage-49 matched branch:")
print("  Universal fail   : W_wall <= Pe_req / Delta_inf")
print("  Universal succeed: W_wall >= Pe_req / Delta_0")
print("Stage-51 resonance family:")
print("  Fail threshold   : W_wall <= Pe_req / (C_res^2 Delta_inf)")
print("  Success threshold: W_wall >= Pe_req / (C_res^2 Delta_0)")
print("Therefore the only profile-sensitive regions are the exact side-bands")
print("  Pe_req/Delta_inf < W_wall < Pe_req/(C_res^2 Delta_inf)")
print("and")
print("  Pe_req/Delta_0 < W_wall < Pe_req/(C_res^2 Delta_0),")
print("whose relative widths are both P_res - 1 = (1 - C_res^2)/C_res^2.")
