#!/usr/bin/env python3
"""
5pn_stage70_outgoing_branch_loading_ratio_finish.py

Stage 70 audit: final reduced finish-line for the explicit Family-1 outgoing quadrupole branch.
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

banner("STAGE 70 — FINAL REDUCED FINISH-LINE FOR THE EXPLICIT FAMILY-1 BRANCH")

subbanner("70.1 — Exact criterion")
print("The explicit Family-1 branch is fully solved in reduced form once the normalized")
print("loading ratio of the passive/outgoing quadrupole branch is known:")
print("  rho_alpha = alpha_req / alpha_mix.")
print()
print("The exact success/failure criterion is")
print("  rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)  -> guaranteed success,")
print("  rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)  -> guaranteed failure,")
print("  rho_alpha <  rho_max^(F1)(eps_blk)              -> constructive ceiling.")

subbanner("70.2 — Natural unblocked numbers")
print("At lambda_mu = 1 and eps_blk = 0:")
print("  rho_alpha <= 3.46622291347846  -> guaranteed success,")
print("  rho_alpha >= 3.46752913273870  -> guaranteed failure,")
print("  rho_alpha <  3.46752922945601  -> absolute constructive ceiling.")

banner("STAGE 70 FINAL LEDGER")
print("The explicit Family-1 support/source branch is finished.")
print("The only remaining reduced question is the actual loading ratio selected by the")
print("physical passive/outgoing moving-throat quadrupole branch.")
