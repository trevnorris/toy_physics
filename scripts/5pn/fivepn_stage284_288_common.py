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


def g_minus(r):
    r = sp.sympify(r)
    return sp.simplify(r - sp.sqrt(1 + r**2) / 2)


def g_minus_prime(r):
    r = sp.sympify(r)
    return sp.simplify(sp.diff(g_minus(r), r))


def R_of_g(g, r):
    g = sp.sympify(g)
    r = sp.sympify(r)
    return sp.simplify((g - r) ** 2 / (1 + r**2))


def family1_reference_values():
    return {
        "Sigma0_can": sp.Float("4.651033550168876"),
        "Tmhat_can": sp.Float("1.4467083664567624"),
        "g_star": sp.Float("0.758035078944663"),
        "R_can": sp.Rational(1, 4),
        "S_can": sp.Float("0.6703621156734617"),
        "Pi_can": sp.Float("3.8715643774790087"),
        "r_F1": sp.Float("1.77799353547498"),
    }


def tangent_log_channels():
    # Stage-147 channel names.
    dln_Zq, dln_rho_w, dln_csw, dln_cs, dln_Tm, dln_vw0, dln_a, dln_LW = sp.symbols(
        "dln_Zq dln_rho_w dln_csw dln_cs dln_Tm dln_vw0 dln_a dln_LW", real=True
    )
    return {
        "dln_Zq": dln_Zq,
        "dln_rho_w": dln_rho_w,
        "dln_csw": dln_csw,
        "dln_cs": dln_cs,
        "dln_Tm": dln_Tm,
        "dln_vw0": dln_vw0,
        "dln_a": dln_a,
        "dln_LW": dln_LW,
    }
