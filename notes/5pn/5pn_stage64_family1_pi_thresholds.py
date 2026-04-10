#!/usr/bin/env python3
"""
5pn_stage64_family1_pi_thresholds.py

Stage 64 audit: final explicit Family-1 quadrupole-demand window in Pi_tr.
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

banner("STAGE 64 — FAMILY-1 Pi_tr THRESHOLDS")

zeta, eps_blk, C_mix = sp.symbols("zeta eps_blk C_mix", positive=True, real=True)

Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta))
print("Q(zeta;eps_blk) =")
sp.pprint(Q)
expect_zero("Q(0)-1", Q.subs(zeta, 0) - 1)
expect_zero("Q(1)-2", sp.simplify(Q.subs(zeta, 1) - 2))
dQ = sp.simplify(sp.diff(Q, zeta))
print("dQ/dzeta =")
sp.pprint(dQ)

mp.mp.dps = 80
zeta_suff_chi = mp.mpf("2.46622291347846")
zeta_fail_chi = mp.mpf("2.46752913273870")
zeta_suff_J = mp.mpf("2.44257571477179")
zeta_max = mp.mpf("2.46752922945601")
eps0 = mp.mpf(0)

Q_num = lambda zz, ee: (1 + (1 - 2 * ee) * zz) / (1 - ee * zz)

print("Pi_suff^(chi)/C_mix at eps_blk=0 =", Q_num(zeta_suff_chi, eps0))
print("Pi_fail^(chi)/C_mix at eps_blk=0 =", Q_num(zeta_fail_chi, eps0))
print("Pi_suff^(J)  /C_mix at eps_blk=0 =", Q_num(zeta_suff_J, eps0))
print("Pi_max^(F1)  /C_mix at eps_blk=0 =", Q_num(zeta_max, eps0))
print("blocking bound eps_blk < 1/zeta_max =", 1 / zeta_max)

banner("STAGE 64 FINAL LEDGER")
print("Stage 64 rewrites the explicit Family-1 theorem directly in the selected-product variable.")
print("At eps_blk = 0 and lambda_mu = 1:")
print("  Pi_tr / C_mix <= 3.46622291347846  -> guaranteed success,")
print("  Pi_tr / C_mix >= 3.46752913273870  -> guaranteed failure,")
print("with the hard explicit ceiling")
print("  Pi_tr / C_mix <  3.46752922945601.")
