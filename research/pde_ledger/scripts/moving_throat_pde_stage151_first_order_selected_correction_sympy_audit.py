#!/usr/bin/env python3
"""
moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py

Audit for the first-order source correction selected by the full mouth profile.

The SymPy engine is an exact multi-point cross-check: it verifies M1-M7 exactly
at sampled rational Pi_star values, while keeping r1, r2, A_T, B_T, and gprime
symbolic.  The full all-Pi_star symbolic proof is carried by the Mathematica
engine.
"""

# === DO NOT "fix" this into a fully-symbolic SymPy proof. =======================
# SymPy CANNOT evaluate  integral_0^1 e^(-Pi_star*x) * {cos,cosh}(...) * x^n dx
# with a SYMBOLIC Pi_star -- it hangs indefinitely (confirmed 2026-05-28; the
# attempt was killed at 35 min and again at 19 min). This script is therefore an
# intentional EXACT, MULTI-POINT cross-check at concrete rational Pi_star values
# (symbolic in r1, r2, A_T, B_T, gprime). The full all-Pi_star symbolic proof is
# carried by the Mathematica engine:
#   mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl
# If you think you can make the SymPy side fully symbolic in Pi_star: you cannot
# (tried definite-integrate, indefinite+bounds, and trig->exp rewrite -- all hang).
# ================================================================================

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name, expr):
    res = sp.simplify(sp.cancel(sp.together(expr)))
    print(f"{name} = {res}")
    assert res == 0, f"{name} nonzero: {res}"


banner("FIRST-ORDER SELF-CONSISTENT SOURCE CORRECTION")

x = sp.symbols("x", real=True)
eps = sp.symbols("epsilon", real=True)
r1, r2 = sp.symbols("r1 r2", real=True)
gprime = sp.symbols("gprime", real=True, nonzero=True)
AT, BT = sp.symbols("A_T B_T", real=True)
k = sp.pi / 2

_sympy_integrate = sp.integrate


def _poly_exp_moment(rate, degree):
    if rate == 0:
        return sp.Rational(1, degree + 1)
    if degree == 0:
        return (sp.exp(rate) - 1) / rate
    previous = _poly_exp_moment(rate, degree - 1)
    return sp.exp(rate) / rate - sp.Integer(degree) * previous / rate


def _expand_linear_exponentials(expr):
    expr = expr.rewrite(sp.exp)
    if expr.is_Atom:
        return expr
    if expr.is_Add:
        return sp.Add(*(_expand_linear_exponentials(arg) for arg in expr.args), evaluate=False)
    if expr.is_Mul:
        terms = [sp.Integer(1)]
        for factor in expr.args:
            expanded = _expand_linear_exponentials(factor)
            factor_terms = expanded.args if expanded.is_Add else (expanded,)
            terms = [term * factor_term for term in terms for factor_term in factor_terms]
        return sp.Add(*terms, evaluate=False)
    if expr.is_Pow:
        return expr.func(
            _expand_linear_exponentials(expr.base),
            _expand_linear_exponentials(expr.exp),
            evaluate=False,
        )
    return expr.func(*(_expand_linear_exponentials(arg) for arg in expr.args), evaluate=False)


def _integrate_exp_poly_term(term):
    coeff, dependent = term.as_independent(x, as_Add=False)
    dependent = sp.powsimp(dependent, force=True)
    rate = sp.Integer(0)

    for exp_factor in list(dependent.atoms(sp.exp)):
        exponent = exp_factor.args[0]
        if not exponent.has(x):
            continue
        slope = sp.diff(exponent, x)
        constant = sp.simplify(exponent - slope * x)
        if slope.has(x) or constant.has(x):
            raise ValueError("nonlinear exponential argument")
        rate += slope
        coeff *= sp.exp(constant)
        dependent = dependent / exp_factor

    dependent = sp.cancel(sp.powsimp(dependent, force=True))
    poly = sp.Poly(dependent, x)
    return sum(
        coeff * poly_coeff * _poly_exp_moment(rate, degree)
        for (degree,), poly_coeff in poly.terms()
    )


def _exact_unit_integrate(expr, *args, **kwargs):
    if len(args) == 1 and args[0] == (x, 0, 1) and not kwargs:
        rewritten = _expand_linear_exponentials(expr)
        terms = rewritten.args if rewritten.is_Add else (rewritten,)
        return sp.cancel(sp.together(sum(_integrate_exp_poly_term(term) for term in terms)))
    return _sympy_integrate(expr, *args, **kwargs)


sp.integrate = _exact_unit_integrate

PI_SAMPLES = [
    sp.Rational(1, 2),
    sp.Integer(1),
    sp.Rational(3, 2),
    sp.Integer(2),
    sp.Rational(5, 3),
]


for Pi_star in PI_SAMPLES:
    print(f"\n[Pi={Pi_star}]")

    w0 = sp.exp(-Pi_star * x)
    pert = r1 * x + r2 * x**2

    num = sp.series(w0 * sp.exp(-eps * pert), eps, 0, 2).removeO()
    Z0 = sp.integrate(num.coeff(eps, 0), (x, 0, 1))
    Z1 = sp.integrate(num.coeff(eps, 1), (x, 0, 1))
    ser = sp.series(num / (Z0 + eps * Z1), eps, 0, 2).removeO()
    Sigma_star = sp.simplify(ser.coeff(eps, 0))
    delta_Sigma = sp.simplify(ser.coeff(eps, 1))

    c_kernel = sp.cos(k * x)
    K_kernel = sp.cosh(k * (1 - x)) / sp.cosh(k)

    mean = lambda f: sp.integrate(Sigma_star * f, (x, 0, 1))
    Rbar = mean(pert)
    cbar = mean(c_kernel)
    Kbar = mean(K_kernel)
    CovcR = mean(c_kernel * pert) - cbar * Rbar
    CovKR = mean(K_kernel * pert) - Kbar * Rbar

    expect_zero(
        f"[Pi={Pi_star}] M1 delta_Sigma + Sigma_star*(R - <R>)",
        delta_Sigma + Sigma_star * (pert - Rbar),
    )
    expect_zero(
        f"[Pi={Pi_star}] M2 int Sigma_star dx - 1",
        sp.integrate(Sigma_star, (x, 0, 1)) - 1,
    )
    expect_zero(
        f"[Pi={Pi_star}] M3 int delta_Sigma dx",
        sp.integrate(delta_Sigma, (x, 0, 1)),
    )

    dg = sp.integrate(c_kernel * delta_Sigma, (x, 0, 1))
    dS = sp.integrate(K_kernel * delta_Sigma, (x, 0, 1))
    expect_zero(
        f"[Pi={Pi_star}] M4 dg + CovcR",
        dg + CovcR,
    )
    expect_zero(
        f"[Pi={Pi_star}] M5 dS + CovKR",
        dS + CovKR,
    )

    deltaPi = -dg / gprime
    deltaT = AT * dg + BT * dS
    expect_zero(
        f"[Pi={Pi_star}] M6 deltaPi - CovcR/gprime",
        deltaPi - CovcR / gprime,
    )
    expect_zero(
        f"[Pi={Pi_star}] M7 deltaT + A_T*CovcR + B_T*CovKR",
        deltaT + AT * CovcR + BT * CovKR,
    )

print("\nTheorem:")
print("  The selected first-order source correction is determined by")
print("  Cov_*(c,R_*) and Cov_*(K_q,R_*) at every exact Pi_star sample above.")
