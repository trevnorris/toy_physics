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
# Stage-77 family-1 Theta extraction:
#   Theta_w^(chi) / lambda_mu^2 = 4.06863235008161550927... (chi^2 weighted datum)
#   Theta_w^(J)   / lambda_mu^2 = 0.92755203253930797184... (Jensen floor)
# Source: scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt:22-23
Theta_chi = sp.Float("4.06863235008162") * lambda_mu**2
Theta_J = sp.Float("0.927552032539308") * lambda_mu**2
# Stage-75 family-1 threshold window:
#   Theta_fail / Pe_req = 0.00036260561797293886969 (= (37 cosh(111 sqrt(5)/5)
#       + 111 sqrt(5) sinh(111 sqrt(5)/5)/5) / (136900 (-1 + sqrt(5) sinh(111 sqrt(5)/5)/3
#       + cosh(111 sqrt(5)/5))))
#   Upsilon_suff / Pe_req = 4.2149534156997728721 ; Theta_suff = Upsilon_suff / 100
# Source: scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt:30,34
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

# --- Branch verdict (chi vs Jensen-floor) ---------------------------------
# Since Theta_J < Theta_chi and both Pe-thresholds share positive denominators,
# the Jensen-floor success threshold and failure ceiling both lie below the
# chi-datum's corresponding thresholds. The two windows nest with overlap.
if not (sp.N(Pe_suff_J) < sp.N(Pe_suff_chi)):
    raise AssertionError(
        "Expected Jensen-floor success threshold to lie below chi-datum's "
        "(requires Theta_J < Theta_chi)"
    )
if not (sp.N(Pe_fail_J) < sp.N(Pe_fail_chi)):
    raise AssertionError(
        "Expected Jensen-floor failure ceiling to lie below chi-datum's "
        "(requires Theta_J < Theta_chi)"
    )
if not (sp.N(Pe_suff_chi) < sp.N(Pe_fail_J)):
    raise AssertionError(
        "Expected chi-datum success threshold below Jensen-floor failure "
        "ceiling (window overlap)"
    )
print("Pe_suff_J < Pe_suff_chi  :", sp.N(Pe_suff_J) < sp.N(Pe_suff_chi))
print("Pe_fail_J < Pe_fail_chi  :", sp.N(Pe_fail_J) < sp.N(Pe_fail_chi))
print("Pe_suff_chi < Pe_fail_J  :", sp.N(Pe_suff_chi) < sp.N(Pe_fail_J))
