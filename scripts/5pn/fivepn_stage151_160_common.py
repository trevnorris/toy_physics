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
    if isinstance(expr, sp.MatrixBase):
        simp = expr.applyfunc(lambda e: sp.simplify(sp.expand(e)))
        print(f"{name} =")
        sp.pprint(simp)
        if any(e != 0 for e in simp):
            raise AssertionError(f"{name} is not zero")
    else:
        simp = sp.simplify(sp.expand(expr))
        print(f"{name} = {simp}")
        if simp != 0:
            raise AssertionError(f"{name} is not zero")

def grouped_triplet(x20: sp.Expr, x21: sp.Expr, x22: sp.Expr):
    xbar = sp.simplify((x20 + 2*x21 + 2*x22)/5)
    a_x = sp.simplify((2*x20 - x21 - x22)/10)
    b_x = sp.simplify((x21 - x22)/2)
    return xbar, a_x, b_x

def grouped_metric():
    return sp.diag(1, 2, 2)

def root1p(r):
    return sp.sqrt(1 + r**2)

# Renormalized co-evolving Family-1 constants used in the notes.
rF1 = sp.Float('1.77799353547498')
gF1 = sp.Float('0.758035078944663')
Sigma0_can = sp.Float('4.651033550168876')
T_hat_can = sp.Float('1.4467083664567624')
S_can = sp.Float('0.6703621156734617')
Pi_can = sp.Float('3.8715643774790087')
