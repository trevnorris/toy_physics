#!/usr/bin/env python3
"""
Stage 61 SymPy audit.

Checks:
1. Insert the explicit Stage-60 Theta values into the Stage-58 threshold window.
2. Compute the exact Pe_req success/failure windows for both the natural quadratic datum
   and the conservative Jensen floor.
3. Verify the ordering Pe_suff < Pe_fail in both cases.
"""
from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


banner("STAGE 61 — FAMILY-1 BRANCH VERDICT")

lambda_mu, Pe_req = sp.symbols("lambda_mu Pe_req", positive=True, real=True)
Theta_chi = sp.Float("4.06863235008162") * lambda_mu**2
Theta_J = sp.Float("0.927552032539308") * lambda_mu**2
Theta_fail = sp.Float("3.62605617972939e-4") * Pe_req
Theta_suff = sp.Float("4.21495341569977e-2") * Pe_req

Pe_suff_chi = sp.simplify(Theta_chi / sp.Float("4.21495341569977e-2") / lambda_mu**2)
Pe_fail_chi = sp.simplify(Theta_chi / sp.Float("3.62605617972939e-4") / lambda_mu**2)
Pe_suff_J = sp.simplify(Theta_J / sp.Float("4.21495341569977e-2") / lambda_mu**2)
Pe_fail_J = sp.simplify(Theta_J / sp.Float("3.62605617972939e-4") / lambda_mu**2)

print("Pe_suff^(chi) / lambda_mu^2 =", sp.N(Pe_suff_chi, 20))
print("Pe_fail^(chi) / lambda_mu^2 =", sp.N(Pe_fail_chi, 20))
print("Pe_suff^(J)   / lambda_mu^2 =", sp.N(Pe_suff_J, 20))
print("Pe_fail^(J)   / lambda_mu^2 =", sp.N(Pe_fail_J, 20))

if not (sp.N(Pe_suff_chi) < sp.N(Pe_fail_chi)):
    raise AssertionError("Expected Pe_suff^(chi) < Pe_fail^(chi)")
if not (sp.N(Pe_suff_J) < sp.N(Pe_fail_J)):
    raise AssertionError("Expected Pe_suff^(J) < Pe_fail^(J)")
