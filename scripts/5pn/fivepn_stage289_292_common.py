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


def grouped_trace_anomaly(x20, x21, x22):
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    return {"bar": xbar, "a": ax, "b": bx}


def grouped_bilinear(x20, x21, x22, y20, y21, y22):
    gx = grouped_trace_anomaly(x20, x21, x22)
    gy = grouped_trace_anomaly(y20, y21, y22)
    return sp.simplify(4 * gx["a"] * gy["a"] + sp.Rational(4, 5) * gx["b"] * gy["b"])
