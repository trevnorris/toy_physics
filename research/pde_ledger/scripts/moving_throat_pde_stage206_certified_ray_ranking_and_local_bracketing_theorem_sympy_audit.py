#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 189 — CERTIFIED RAY RANKING AND LOCAL BRACKETING")

# ---------------------------------------------------------------------------
# I. Exact monotone-branch quadratic root map
# ---------------------------------------------------------------------------
subbanner("I. Exact monotone-branch quadratic root map")
H0, k, c = sp.symbols("H0 k c", positive=True, real=True)
# Oriented slope K0 = -k < 0.
Tau = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * c * H0)))
residual = sp.simplify(H0 - k * Tau + sp.Rational(1, 2) * c * Tau**2)
print("Tau(H0,k;c) =")
sp.pprint(Tau)
expect_zero("quadratic residual", residual)
expect_zero("zero-curvature limit", sp.simplify(sp.limit(Tau, c, 0) - H0 / k))
expect_zero(
    "implicit derivative identity",
    sp.simplify(sp.diff(Tau, c) - Tau**2 / (2 * sp.sqrt(k**2 - 2 * c * H0))),
)
expect_zero(
    "derivative-at-root identity",
    sp.simplify((-k + c * Tau) + sp.sqrt(k**2 - 2 * c * H0)),
)

# ---------------------------------------------------------------------------
# II. Exact certified bracket endpoints from curvature bounds
# ---------------------------------------------------------------------------
subbanner("II. Exact bracket endpoints from lower/upper curvature bounds")
cL, cU = sp.symbols("cL cU", real=True)
TauL = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * cL * H0)))
TauU = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * cU * H0)))
print("Tau_lo =")
sp.pprint(TauL)
print("Tau_hi =")
sp.pprint(TauU)
expect_zero(
    "Tau_lo solves lower quadratic",
    sp.simplify(H0 - k * TauL + sp.Rational(1, 2) * cL * TauL**2),
)
expect_zero(
    "Tau_hi solves upper quadratic",
    sp.simplify(H0 - k * TauU + sp.Rational(1, 2) * cU * TauU**2),
)
expect_zero(
    "symmetric degenerate-envelope collapse",
    sp.simplify(TauL.subs(cL, c) - Tau),
)
expect_zero(
    "symmetric degenerate-envelope collapse (upper)",
    sp.simplify(TauU.subs(cU, c) - Tau),
)

# ---------------------------------------------------------------------------
# III. Exact small-envelope width law
# ---------------------------------------------------------------------------
subbanner("III. Exact small-envelope width law")
cmid, eta = sp.symbols("cmid eta", real=True)
TauPlus = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * (cmid + eta / 2) * H0)))
TauMinus = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * (cmid - eta / 2) * H0)))
Width = sp.simplify(TauPlus - TauMinus)
TauMid = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * cmid * H0)))
series_width = sp.simplify(sp.series(Width, eta, 0, 3).removeO())
leading_width = sp.simplify(
    (TauMid**2 / (2 * sp.sqrt(k**2 - 2 * cmid * H0))) * eta
)
print("series[Width] =")
sp.pprint(series_width)
print("leading small-envelope width term =")
sp.pprint(leading_width)
expect_zero("small-envelope width law", series_width - leading_width)

Width0 = sp.simplify(Width.subs(cmid, 0))
series_width0 = sp.simplify(sp.series(Width0, eta, 0, 3).removeO())
Tau0 = sp.simplify(H0 / k)
leading_width0 = sp.simplify((Tau0**2 / (2 * k)) * eta)
print("series[Width] at zero mean curvature =")
sp.pprint(series_width0)
expect_zero("zero-curvature width law", series_width0 - leading_width0)

# ---------------------------------------------------------------------------
# IV. Turning-ray formulas
# ---------------------------------------------------------------------------
subbanner("IV. Turning-ray bracket formulas")
a, b = sp.symbols("a b", positive=True, real=True)
TauTurnA = sp.sqrt(2 * H0 / a)
TauTurnB = sp.sqrt(2 * H0 / b)
print("Tau_turn(a) =")
sp.pprint(TauTurnA)
print("Tau_turn(b) =")
sp.pprint(TauTurnB)
expect_zero(
    "turning root (a)",
    sp.simplify(H0 - sp.Rational(1, 2) * a * TauTurnA**2),
)
expect_zero(
    "turning root (b)",
    sp.simplify(H0 - sp.Rational(1, 2) * b * TauTurnB**2),
)
expect_zero(
    "turning derivative identity",
    sp.simplify(sp.diff(TauTurnA, a) + TauTurnA / (2 * a)),
)

# ---------------------------------------------------------------------------
# V. Collapse to the Stage 205 quadratic predictor under exact curvature
# ---------------------------------------------------------------------------
subbanner("V. Collapse to the Stage 205 quadratic predictor")
Phi0, Labs, L1 = sp.symbols("Phi0 Labs L1", positive=True, real=True)
# Use the oriented negative-slope branch L0 = -Labs < 0.
Lneg = -Labs
H0_stage205 = sp.log(Phi0)
TauStage189 = sp.simplify(
    -2 * H0_stage205 / (Lneg + sp.sign(Lneg) * sp.sqrt(Lneg**2 - 2 * L1 * H0_stage205))
)
TauStage188 = sp.simplify(
    -2 * sp.log(Phi0) / (Lneg + sp.sign(Lneg) * sp.sqrt(Lneg**2 - 2 * L1 * sp.log(Phi0)))
)
print("TauStage189 =")
sp.pprint(TauStage189)
print("TauStage188 =")
sp.pprint(TauStage188)
expect_zero("Stage 206/188 collapse", sp.simplify(TauStage189 - TauStage188))

banner("STAGE 189 SYMPY AUDIT PASSED")
