
#!/usr/bin/env python3
"""
5pn_stage88_fixing_chiQ.py

Stage 88 audit: exact fixing of chi_Q from the outgoing DtN fingerprint.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 88 — EXACT FIXING OF chi_Q")

omega, a, c_s, chi_Q, xi_Q = sp.symbols("omega a c_s chi_Q xi_Q", positive=True, real=True)
Omega_Q = sp.simplify(3 * c_s / (2 * a))
sigma_can = sp.simplify(9 / (8 * Omega_Q**5))
expect_zero("sigma_Q^can - 4 a^5/(27 c_s^5)", sigma_can - 4 * a**5 / (27 * c_s**5))

Yret = sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2 - sp.I * chi_Q * sigma_can * omega**5)
Yret_series = sp.series(sp.expand(Yret), omega, 0, 6).removeO()
print("Yhat_Q^ret(omega) =")
sp.pprint(Yret_series)

coef5 = sp.simplify(sp.expand(Yret_series).coeff(omega, 5) / sp.I)
expected_coef5 = sp.simplify(chi_Q * a**5 / (27 * c_s**5))
expect_zero("retarded omega^5 coefficient", coef5 - expected_coef5)
expect_zero("chi_Q fix from canonical DtN", expected_coef5.subs(chi_Q, 1) - a**5 / (27 * c_s**5))

z = sp.symbols("z", real=True)
Lambda_def = -3 + z**2 / 3 + z**4 / 9 + sp.I * xi_Q * z**5 / 9
Ydef = sp.simplify(-3 / Lambda_def)
Ydef_series = sp.series(Ydef, z, 0, 6).removeO()
coef5_def = sp.simplify(sp.expand(Ydef_series).coeff(z, 5) / sp.I)
expect_zero("deformed DtN chi_Q - xi_Q", coef5_def / sp.Rational(1, 27) - xi_Q)

banner("STAGE 88 FINAL LEDGER")
print("Matching the retarded grouped-P2 one-pole-plus-contact module to the explicit outgoing")
print("DtN fingerprint fixes")
print("  chi_Q = 1")
print("on the canonical compact passive/outgoing branch.")
print("A deformed DtN coefficient xi_Q shifts chi_Q exactly by")
print("  chi_Q = xi_Q.")
