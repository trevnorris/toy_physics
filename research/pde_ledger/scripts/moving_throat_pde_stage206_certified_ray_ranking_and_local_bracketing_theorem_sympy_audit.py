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
# V. Collapse to the Stage 239 quadratic log predictor under exact curvature
# ---------------------------------------------------------------------------
subbanner("V. Collapse to the Stage 239 quadratic log predictor")
Phi0, L1 = sp.symbols("Phi0 L1", positive=True, real=True)
L0 = sp.symbols("L0", real=True)
TauBracketLog = sp.simplify(Tau.subs({H0: sp.log(Phi0), k: -L0, c: L1}))
TauLog2 = -2 * sp.log(Phi0) / (
    L0 + sp.sign(L0) * sp.sqrt(L0**2 - 2 * L1 * sp.log(Phi0))
)
print("T_bracket(H0=log(Phi0), K0=L0; c=L1) =")
sp.pprint(TauBracketLog)
print("Stage-239 tau_log2 =")
sp.pprint(TauLog2)
expect_zero(
    "Stage 206/239 log-predictor collapse",
    sp.simplify(TauBracketLog - sp.refine(TauLog2, sp.Q.negative(L0))),
)

# ---------------------------------------------------------------------------
# VI. Pairwise ray ordering and local search-sieve admissibility
# ---------------------------------------------------------------------------
subbanner("VI. Pairwise ray ordering and local search-sieve admissibility")

slack_a_hi, slack_sep, slack_b_lo = sp.symbols(
    "slack_a_hi slack_sep slack_b_lo", nonnegative=True, real=True
)
strict_sep = sp.symbols("strict_sep", positive=True, real=True)
tau_star_a = sp.symbols("tau_star_a", real=True)
tau_hi_a = tau_star_a + slack_a_hi
tau_lo_b = tau_hi_a + strict_sep
tau_star_b = tau_lo_b + slack_b_lo
ordering_gap = sp.simplify(tau_star_b - tau_star_a)
ordering_proven = sp.ask(sp.Q.positive(ordering_gap))
negation_unsat = sp.ask(sp.Q.nonpositive(ordering_gap))
print(f"pairwise ordering gap = {ordering_gap}")
print(f"pairwise ordering theorem = {ordering_proven}")
print(f"pairwise ordering negation satisfiable = {negation_unsat}")
if ordering_proven is not True or negation_unsat is not False:
    raise AssertionError("pairwise ordering implication was not discharged")

relaxed_counterexample = {
    "tau_lo_a": sp.Integer(0),
    "tau_star_a": sp.Integer(1),
    "tau_hi_a": sp.Integer(2),
    "tau_lo_b": sp.Integer(0),
    "tau_star_b": sp.Rational(1, 2),
    "tau_hi_b": sp.Integer(2),
}
relaxed_hypotheses_hold = (
    relaxed_counterexample["tau_lo_a"]
    <= relaxed_counterexample["tau_star_a"]
    <= relaxed_counterexample["tau_hi_a"]
    and relaxed_counterexample["tau_lo_b"]
    <= relaxed_counterexample["tau_star_b"]
    <= relaxed_counterexample["tau_hi_b"]
)
relaxed_conclusion = (
    relaxed_counterexample["tau_star_a"] < relaxed_counterexample["tau_star_b"]
)
print(f"pairwise theorem without separation counterexample hypotheses = {relaxed_hypotheses_hold}")
print(f"pairwise theorem without separation counterexample conclusion = {relaxed_conclusion}")
if not relaxed_hypotheses_hold or relaxed_conclusion:
    raise AssertionError("dropping separation did not produce the required counterexample")


def local_sieve_admissible(H, K, c_lo, c_hi, T_valid, tau_hi, tau_hi_tp):
    monotone = sp.And(
        H > 0,
        K < 0,
        c_lo <= c_hi,
        K**2 - 2 * c_lo * H >= 0,
        K**2 - 2 * c_hi * H >= 0,
        tau_hi <= T_valid,
    )
    turning = sp.And(
        H > 0,
        sp.Eq(K, 0),
        c_lo <= c_hi,
        c_hi < 0,
        tau_hi_tp <= T_valid,
    )
    return sp.simplify_logic(sp.Or(monotone, turning))


monotone_good_tau = sp.simplify(TauU.subs({H0: 1, k: 3, cU: 1}))
turning_good_tau = sp.sqrt(sp.Rational(2) * 2 / 1)
monotone_good = local_sieve_admissible(1, -3, 0, 1, 1, monotone_good_tau, 0)
turning_good = local_sieve_admissible(2, 0, -3, -1, 3, 0, turning_good_tau)
monotone_bad = local_sieve_admissible(1, -3, 0, 1, 1, 2, 0)
print(f"local sieve monotone admissible case = {monotone_good}")
print(f"local sieve turning admissible case = {turning_good}")
print(f"local sieve single-clause violation case = {monotone_bad}")
if monotone_good is not sp.S.true or turning_good is not sp.S.true:
    raise AssertionError("local sieve rejected an admissible bracket")
if monotone_bad is not sp.S.false:
    raise AssertionError("local sieve accepted a bracket with tau_hi > T")

banner("STAGE 189 SYMPY AUDIT PASSED")
