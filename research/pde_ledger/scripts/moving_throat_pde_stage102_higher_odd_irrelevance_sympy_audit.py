#!/usr/bin/env python3
"""Stage 102 SymPy audit: higher odd terms are irrelevant at omega^5.

The 2.5PN observable corresponds to the imaginary part of the omega^5
coefficient of the canonical compact-source quadrupole response. We verify:

  D1 — the omega^5 coefficient is independent of tauQ (the omega^7 inverse-
       denominator perturbation cannot back-propagate to omega^5).
  D2 — the omega^7 coefficient's tauQ-derivative equals 1/4 (the leading
       imaginary place where tauQ surfaces).
  D3 — the canonical odd coefficient at omega^5 equals chiQ * 9/(32 Omega^5).

These three asserts gate the result that would otherwise be a printed-only
inspection. They mirror the Mathematica audit's three expectZero anchors.
"""
from __future__ import annotations
import sympy as sp


def expect_zero(name: str, expr, tol=None) -> None:
    val = sp.simplify(expr)
    if tol is None:
        if val != 0:
            raise AssertionError(f"{name}: residual {val} not zero")
    else:
        try:
            num = float(sp.N(sp.Abs(val), 30))
        except (TypeError, ValueError):
            raise AssertionError(f"{name}: residual {val} non-numeric") from None
        if num > float(tol):
            raise AssertionError(f"{name}: residual {val} > tol {tol}")
    print(f"PASS: {name} (residual = {val})")


omega = sp.symbols('omega', real=True)
Omega, chiQ, tauQ = sp.symbols('Omega chiQ tauQ', positive=True, real=True)

sigma_can = sp.Rational(9, 8) / Omega**5
Y = sp.Rational(3,4) + sp.Rational(1,4) / (1 - omega**2/Omega**2 - sp.I*chiQ*sigma_can*omega**5 - sp.I*tauQ*omega**7)
Yser5 = sp.expand(sp.series(Y, omega, 0, 6).removeO())
Yser8 = sp.expand(sp.series(Y, omega, 0, 8).removeO())

print('series through O(omega^5) =', Yser5)
print('series through O(omega^7) =', Yser8)

tau5 = sp.simplify(sp.diff(sp.im(Yser5.coeff(omega, 5)), tauQ))
tau7 = sp.simplify(sp.diff(sp.im(Yser8.coeff(omega, 7)), tauQ))
print('tauQ coefficient in omega^5 term =', tau5)
print('tauQ coefficient in omega^7 term =', tau7)

expect_zero("D1: tauQ irrelevance at omega^5", tau5)
expect_zero("D2: tauQ coefficient at omega^7 - 1/4", tau7 - sp.Rational(1, 4))
expect_zero("D3: canonical odd coefficient at omega^5",
            sp.im(Yser5.coeff(omega, 5)) - chiQ * sp.Rational(9, 32) / Omega**5)

print("\nSTAGE 102 AUDIT PASSED")
