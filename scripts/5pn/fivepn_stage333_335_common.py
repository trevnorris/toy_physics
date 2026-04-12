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

def c_mix(Lambda, eps):
    return sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2)

def s_req_from_products(Pi_tr, C_mix):
    return sp.simplify(Pi_tr / C_mix)

def zeta_req_from_products(Pi_tr, C_mix, eps):
    return sp.simplify((Pi_tr - C_mix) / (C_mix - eps * (2*C_mix - Pi_tr)))

def zeta_req_from_ratio(rho_alpha, eps):
    return sp.simplify((rho_alpha - 1) / (1 - eps * (2 - rho_alpha)))

def precursor_coeffs_from_ratio(rho_alpha):
    c0 = sp.simplify(1 / rho_alpha)
    c1 = sp.simplify((rho_alpha - 1) / rho_alpha)
    return {"c0": c0, "c1": c1}

def ratio_from_precursor(c0, c1):
    return {
        "rho_alpha_from_c0": sp.simplify(1 / c0),
        "rho_alpha_from_c1": sp.simplify(1 / (1 - c1)),
        "zeta_req": sp.simplify(c1 / c0),
    }
