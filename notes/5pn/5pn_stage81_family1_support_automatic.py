
#!/usr/bin/env python3
"""
5pn_stage81_family1_support_automatic.py

Stage 81 audit: the explicit Family-1 support test is automatic on the actual
isotropic branch.
"""

from __future__ import annotations
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

banner("STAGE 81 — THE EXPLICIT FAMILY-1 SUPPORT TEST IS AUTOMATIC")

eps_blk, zeta_max = sp.symbols("eps_blk zeta_max", positive=True, real=True)

subbanner("81.1 — Exact blocked demand on the actual isotropic branch")
zeta_req = sp.simplify((sp.Rational(4, 3) - 1) / (1 - eps_blk * (2 - sp.Rational(4, 3))))
zeta_req_expected = sp.simplify(1 / (3 - 2 * eps_blk))
expect_zero("zeta_req^(act) - 1/(3 - 2 eps_blk)", zeta_req - zeta_req_expected)
dz = sp.diff(zeta_req_expected, eps_blk)
print("d zeta_req^(act) / d eps_blk =", sp.simplify(dz))

subbanner("81.2 — Automatic support theorem for any branch with zeta_max > 1")
bound = sp.simplify(1 / (3 - 2 / zeta_max))
difference = sp.simplify(zeta_max - bound)
print("zeta_max - worst-case bound =", sp.factor(difference))
# Positive for zeta_max > 1, zeta_max > 2/3.
print("Hence any explicit branch with zeta_max > 1 already passes the support test throughout")
print("the admissible blocked regime 0 <= eps_blk < 1/zeta_max.")

subbanner("81.3 — Family-1 specialization")
mp.mp.dps = 50
zeta_max_F1 = mp.mpf("2.46752922945601")
bound_F1 = 1 / (3 - 2 / zeta_max_F1)
print("zeta_max^(F1) =", zeta_max_F1)
print("worst-case blocked demand bound =", bound_F1)
print("margin =", zeta_max_F1 - bound_F1)

banner("STAGE 81 FINAL LEDGER")
print("On the actual isotropic branch, the explicit Family-1 support/source side is automatic:")
print("  zeta_req^(act)(eps_blk) = 1 / (3 - 2 eps_blk),")
print("and any branch with zeta_max > 1 passes throughout the admissible blocked regime.")
print("Since zeta_max^(F1) ≈ 2.46752922945601 > 1, Family-1 is safely inside the support window.")
