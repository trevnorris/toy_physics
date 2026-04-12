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

def G_tr(xi, delta, R):
    return sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))

def F_tr(xi, delta, R):
    num = (9 * delta + (9 + 2 * R**2) * xi)**2 * (9 * delta + (9 + 2 * R) * xi)**2
    den = 81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2)**2
    return sp.simplify(num / den)

def Pi_tr(xi, delta, R):
    return sp.simplify(F_tr(xi, delta, R) * G_tr(xi, delta, R))

def c_mix(Lambda, eps):
    return sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2)

def rho_alpha_from_products(Pi_tr_value, C_mix_value):
    return sp.simplify(Pi_tr_value / C_mix_value)

def zeta_req_from_rho(rho_alpha, eps):
    return sp.simplify((rho_alpha - 1) / (1 - eps * (2 - rho_alpha)))

def xi_threshold(k, Mmix, delta, R):
    return sp.simplify(
        (
            k * Mmix * (9 + 2 * R**2)
            - 9 * delta
            + sp.sqrt((k * Mmix * (9 + 2 * R**2) - 9 * delta)**2 + 324 * k * Mmix * delta)
        ) / 18
    )

def ZW_threshold(k, xi, delta, R, eps_eta, eps, chi0):
    return sp.simplify(
        sp.pi**2 * (1 - eps_eta) * (1 - eps) * G_tr(xi, delta, R)
        / (8 * k * (1 + chi0)**2)
    )

def Lambda_threshold(k, xi, delta, R, eps):
    return sp.simplify(sp.pi**2 * Pi_tr(xi, delta, R) / (8 * k * (1 - eps)))
