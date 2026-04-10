
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

def expect_zero(name: str, expr):
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

def grouped_trace_anomaly(x20, x21, x22):
    xbar = sp.simplify((x20 + 2*x21 + 2*x22)/5)
    ax = sp.simplify((2*x20 - x21 - x22)/10)
    bx = sp.simplify((x21 - x22)/2)
    return {"bar": xbar, "a": ax, "b": bx}

def grouped_inverse(xbar, ax, bx):
    x20 = sp.simplify(xbar + 4*ax)
    x21 = sp.simplify(xbar - ax + bx)
    x22 = sp.simplify(xbar - ax - bx)
    return {"20": x20, "21": x21, "22": x22}

def response_moments_from_D(D0, D2, D4):
    u2 = sp.simplify(-D2 / D0)
    u4 = sp.simplify((D2**2 - D0*D4) / D0**2)
    return {"u2": u2, "u4": u4}

def prefactor_moments(D0, D2, D4, N0, N2, N4):
    P0 = sp.simplify(N0 / D0)
    P2 = sp.simplify((D0*N2 - 2*D2*N0) / D0**2)
    P4 = sp.simplify((D0**2*N4 - 2*D0*(D2*N2 + D4*N0) + 3*D2**2*N0) / D0**3)
    return {"P0": P0, "P2": P2, "P4": P4}

def pole_defect_from_D(D0, D2, D4):
    moms = response_moments_from_D(D0, D2, D4)
    return sp.simplify(moms["u4"] - 4*moms["u2"]**2)

def orbit_symbols(prefix: str = ""):
    chi0s = sp.symbols(f"{prefix}chi0_star", positive=True, real=True)
    Fstar = sp.symbols(f"{prefix}F_star", real=True)
    return {"chi0s": chi0s, "Fstar": Fstar}

def q_from_mismatch(chi0s, Fstar, mT, mK, mMu):
    qtr = sp.simplify((1 + chi0s) * sp.log(mT))
    qeta = sp.simplify(-sp.log(mK))
    qnt = sp.simplify(sp.log(mMu) - sp.log(mK) - Fstar*sp.log(mT))
    return {"qtr": qtr, "qnt": qnt, "qeta": qeta}

def mismatch_from_q(chi0s, Fstar, qtr, qnt, qeta):
    mT = sp.simplify(sp.exp(qtr/(1 + chi0s)))
    mK = sp.simplify(sp.exp(-qeta))
    mMu = sp.simplify(sp.exp(qnt - qeta + Fstar*qtr/(1 + chi0s)))
    return {"mT": mT, "mK": mK, "mMu": mMu}

def mismatch_from_invariants(chi0s, Fstar, Rtr, Rnt, Reta):
    mT = sp.simplify(Rtr ** (sp.Integer(1)/(1 + chi0s)))
    mK = sp.simplify(1 / Reta)
    mMu = sp.simplify(Rnt * mK * mT**Fstar)
    return {"mT": mT, "mK": mK, "mMu": mMu}

def q_from_invariants(Rtr, Rnt, Reta):
    return {"qtr": sp.simplify(sp.log(Rtr)), "qnt": sp.simplify(sp.log(Rnt)), "qeta": sp.simplify(sp.log(Reta))}

def grouped_distance(a2, b2, a4, b4):
    return sp.simplify(a2**2 + b2**2 + a4**2 + b4**2)
