
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

def lower_branch_g(r):
    return sp.simplify(r - sp.sqrt(1 + r**2) / 2)

def lower_branch_gprime(r):
    return sp.simplify(sp.diff(lower_branch_g(r), r))

def root1p(r):
    return sp.sqrt(1 + r**2)

rF1 = sp.Float('1.77799353547498')
gF1 = sp.Float('0.758035078944663')
Sigma0_can = sp.Float('4.651033550168876')
T_hat_can = sp.Float('1.4467083664567624')
S_can = sp.Float('0.6703621156734617')
Pi_can = sp.Float('3.8715643774790087')
