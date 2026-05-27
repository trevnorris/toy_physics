#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from sympy.functions.special.bessel import jn, yn

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")

z = sp.symbols("z")
h2 = sp.expand_func(jn(2, z) + sp.I * yn(2, z))
Lam = sp.simplify(z * sp.diff(h2, z) / h2)
Y = sp.simplify(sp.Integer(-3) / Lam)

print("h_2^(1)(z) =")
sp.pprint(h2)

Lam_series = sp.expand(sp.series(Lam, z, 0, 8).removeO())
Y_series = sp.expand(sp.series(Y, z, 0, 8).removeO())

print("\nLambda_2^out(z) =")
sp.pprint(Lam_series)

print("\nY_2^out(z) =")
sp.pprint(Y_series)

expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
expect_zero("Y z^2 coefficient", Y_series.coeff(z, 2) - sp.Rational(1, 9))
expect_zero("Y z^4 coefficient", Y_series.coeff(z, 4) - sp.Rational(4, 81))
expect_zero("Y imag z^5 coefficient", Y_series.coeff(z, 5)/sp.I - sp.Rational(1, 27))
expect_zero("Y z^6 coefficient", Y_series.coeff(z, 6) + sp.Rational(11, 729))
expect_zero("Y imag z^7 coefficient", Y_series.coeff(z, 7)/sp.I + sp.Rational(1, 243))

a, omega, c_s = sp.symbols("a omega c_s", positive=True, real=True)
zsub = a * omega / c_s
Yw = sp.expand(sp.series(Y.subs(z, zsub), omega, 0, 6).removeO())
print("\nY_2^out(omega) through O(omega^5) =")
sp.pprint(Yw)

expect_zero("omega^2 coefficient", Yw.coeff(omega, 2) - a**2/(9*c_s**2))
expect_zero("omega^4 coefficient", Yw.coeff(omega, 4) - 4*a**4/(81*c_s**4))
expect_zero("imag omega^5 coefficient", Yw.coeff(omega, 5)/sp.I - a**5/(27*c_s**5))

print("\nRESULT:")
print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
