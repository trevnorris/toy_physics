#!/usr/bin/env python3
"""
5pn_stage52_final_reduced_verdict.py

Stage 52 audit: final reduced verdict for the support/source program.
"""

from __future__ import annotations

import math
import mpmath as mp
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
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 52 — FINAL REDUCED VERDICT FOR THE SUPPORT/SOURCE PROGRAM")

Pe_req, Delta0, Deltainf = sp.symbols("Pe_req Delta0 Deltainf", positive=True, real=True)
P_res = sp.Symbol("P_res", positive=True, real=True)

subbanner("52.1 — Universal matched-branch theorem window")
W_fail = sp.simplify(Pe_req / Deltainf)
W_suff = sp.simplify(Pe_req / Delta0)
print("Universal matched-branch theorem:")
print("  W_wall <= W_fail  -> fail")
print("  W_wall >= W_suff  -> succeed")
print("with")
print("  W_fail =")
sp.pprint(W_fail)
print("  W_suff =")
sp.pprint(W_suff)

subbanner("52.2 — Explicit profile-family refinement")
W_fail_res = sp.simplify(P_res * Pe_req / Deltainf)
W_suff_res = sp.simplify(P_res * Pe_req / Delta0)
print("Independent-profile resonance family thresholds:")
print("  W_wall <= P_res Pe_req / Delta_inf  -> fail on the family")
print("  W_wall >= P_res Pe_req / Delta_0    -> succeed on the family")
print("with")
print("  W_fail^(res) =")
sp.pprint(W_fail_res)
print("  W_suff^(res) =")
sp.pprint(W_suff_res)

subbanner("52.3 — Exact narrow profile-sensitive bands")
band_fail = sp.simplify((P_res - 1) * Pe_req / Deltainf)
band_suff = sp.simplify((P_res - 1) * Pe_req / Delta0)
print("Only the narrow sub-bands of width")
print("  (P_res - 1) Pe_req / Delta_inf")
print("and")
print("  (P_res - 1) Pe_req / Delta_0")
print("can be profile-sensitive at reduced-theorem level.")
print("band_fail =")
sp.pprint(band_fail)
print("band_suff =")
sp.pprint(band_suff)

subbanner("52.4 — Numerical resonance penalty from Stage 50")
mp.mp.dps = 60
C_res2 = mp.mpf("0.9944188364515293487")
P_res_num = 1 / C_res2
print("C_res^2 =", C_res2)
print("P_res   =", P_res_num)
print("P_res - 1 =", P_res_num - 1)

banner("STAGE 52 FINAL LEDGER")
print("Stage 52 closes the reduced support/source audit.")
print("The universal theorem remains the matched-branch window [Pe_req/Delta_inf, Pe_req/Delta_0].")
print("The explicit sech–Gaussian benchmark only perturbs those thresholds by the tiny factor")
print("  P_res = 1.005612487760576...")
print("so the genuinely profile-sensitive region is only about 0.56% wide.")
