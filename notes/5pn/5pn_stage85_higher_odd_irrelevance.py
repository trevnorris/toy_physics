
#!/usr/bin/env python3
"""
5pn_stage85_higher_odd_irrelevance.py

Stage 85 audit: at 2.5PN the only live retarded obstruction is the leading
omega^5 outgoing normalization.
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

banner("STAGE 85 — HIGHER ODD IRRELEVANCE AT 2.5PN")

omega, Omega_Q, chi_Q, sigma_can, tau_Q = sp.symbols(
    "omega Omega_Q chi_Q sigma_can tau_Q", real=True
)
Yret = sp.Rational(3, 4) + sp.Rational(1, 4) / (
    1 - omega**2 / Omega_Q**2 - sp.I * chi_Q * sigma_can * omega**5 - sp.I * tau_Q * omega**7
)
series_Y = sp.series(sp.expand(Yret), omega, 0, 6).removeO()
print("Yhat_Q^ret through O(omega^5) =")
sp.pprint(series_Y)

coef_w5 = sp.simplify(sp.expand(series_Y).coeff(omega, 5) / sp.I)
expected_coef_w5 = sp.simplify(chi_Q * sigma_can / 4)
expect_zero("omega^5 coefficient", coef_w5 - expected_coef_w5)
expect_zero("tau_Q contribution through O(omega^5)", sp.diff(coef_w5, tau_Q))

banner("STAGE 85 FINAL LEDGER")
print("Any extra retarded structure that first enters at O(omega^7) or higher is irrelevant")
print("to the 2.5PN theorem. The only live retarded obstruction at 2.5PN is the leading")
print("omega^5 outgoing-normalization factor chi_Q.")
