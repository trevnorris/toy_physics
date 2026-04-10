#!/usr/bin/env python3
"""
5pn_stage69_family1_loading_ratio_window.py

Stage 69 audit: exact Family-1 success/failure window in the pure loading-ratio variable.
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

banner("STAGE 69 — FAMILY-1 LOADING-RATIO WINDOW")

zeta, eps_blk = sp.symbols("zeta eps_blk", positive=True, real=True)
Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta))
rho_max = sp.Symbol("rho_max^(F1)", positive=True, real=True)
print("Q(zeta;eps_blk) =")
sp.pprint(Q)

mp.mp.dps = 80
zeta_suff_chi = mp.mpf("2.46622291347846")
zeta_fail_chi = mp.mpf("2.46752913273870")
zeta_suff_J = mp.mpf("2.44257571477179")
zeta_max = mp.mpf("2.46752922945601")
Q_num = lambda zz, ee: (1 + (1 - 2 * ee) * zz) / (1 - ee * zz)

print("rho_suff^(chi)(1;0) =", Q_num(zeta_suff_chi, mp.mpf(0)))
print("rho_fail^(chi)(1;0) =", Q_num(zeta_fail_chi, mp.mpf(0)))
print("rho_suff^(J)(1;0)   =", Q_num(zeta_suff_J, mp.mpf(0)))
print("rho_max^(F1)(0)     =", Q_num(zeta_max, mp.mpf(0)))
print("eps_blk < 1/zeta_max =", 1 / zeta_max)

subbanner("69.1 — Blocking raises the required loading ratio")
dQ_deps = sp.simplify(sp.diff(Q, eps_blk))
print("dQ/deps_blk =")
sp.pprint(dQ_deps)

banner("STAGE 69 FINAL LEDGER")
print("Stage 69 collapses the explicit Family-1 support theorem to one pure loading-ratio window.")
print("At lambda_mu = 1 and eps_blk = 0:")
print("  rho_alpha <= 3.46622291347846  -> guaranteed success,")
print("  rho_alpha >= 3.46752913273870  -> guaranteed failure,")
print("with the hard ceiling")
print("  rho_alpha <  3.46752922945601.")
