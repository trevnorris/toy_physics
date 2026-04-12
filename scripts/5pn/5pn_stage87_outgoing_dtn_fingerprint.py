
#!/usr/bin/env python3
"""
5pn_stage87_outgoing_dtn_fingerprint.py

Stage 87 audit: exact outgoing l=2 Dirichlet-to-Neumann fingerprint.
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

banner("STAGE 87 — EXACT OUTGOING l=2 DtN FINGERPRINT")

z = sp.symbols("z")
j2 = ((3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2)
y2 = (-(3 / z**3 - 1 / z) * sp.cos(z) - 3 * sp.sin(z) / z**2)
h2 = sp.simplify(j2 + sp.I * y2)
Lambda_out = sp.simplify(z * sp.diff(sp.log(h2), z))
Lambda_series = sp.series(Lambda_out, z, 0, 8).removeO()
print("Lambda_2^out(z) =")
sp.pprint(Lambda_series)

Lambda_target = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9 - 2 * z**6 / 27 - sp.I * z**7 / 27
expect_zero("Lambda_2^out - target", sp.expand(Lambda_series - Lambda_target))

Yhat_out = sp.simplify(-3 / Lambda_out)
Yhat_series = sp.series(Yhat_out, z, 0, 8).removeO()
print("Yhat_2^out(z) =")
sp.pprint(Yhat_series)

Yhat_target = 1 + z**2 / 9 + 4 * z**4 / 81 + sp.I * z**5 / 27 - 11 * z**6 / 729 - sp.I * z**7 / 243
expect_zero("Yhat_2^out - target", sp.expand(Yhat_series - Yhat_target))

banner("STAGE 87 FINAL LEDGER")
print("The exact outgoing spherical l=2 DtN model gives")
print("  Lambda_2^out = -3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + O(z^8),")
print("and therefore the normalized outgoing quadrupole branch is")
print("  Yhat_2^out = 1 + z^2/9 + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243 + O(z^8).")
