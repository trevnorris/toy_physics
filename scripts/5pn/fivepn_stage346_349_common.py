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
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")

def dlog(expr, var, dlog_var):
    return sp.simplify((var * sp.diff(expr, var) / expr) * dlog_var)

def linearized_log(expr, subst_map):
    t = sp.symbols('t', real=True)
    out = sp.log(expr)
    out = out.subs({k: k*sp.exp(v*t) for k, v in subst_map.items()})
    return sp.simplify(sp.diff(out, t).subs(t, 0))
