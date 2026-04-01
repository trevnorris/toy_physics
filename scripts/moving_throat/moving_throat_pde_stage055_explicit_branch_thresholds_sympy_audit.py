#!/usr/bin/env python3
"""
moving_throat_pde_stage55_explicit_branch_thresholds_sympy_audit.py

SymPy audit for Stage 55:
explicit branch placement map and threshold surfaces.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 55 — EXPLICIT BRANCH THRESHOLD SURFACES")

chi_s, Lambda_ell, Upsilon_w, Pe_req = sp.symbols("chi_s Lambda_ell Upsilon_w Pe_req", positive=True, real=True)
rho_w, c_sw, V0, hbar = sp.symbols("rho_w c_sw V0 hbar", positive=True, real=True)

kappa = sp.simplify(4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2)
eta = Lambda_ell
alpha = sp.sqrt(kappa)

Delta0 = sp.simplify(
    eta * (sp.cosh(alpha) - 1) /
    (alpha**2 * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))
)
DeltaInf = sp.simplify(
    (sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) /
    (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))
)

Upsilon_fail = sp.simplify(Pe_req / (Lambda_ell**2 * DeltaInf))
Upsilon_suff = sp.simplify(Pe_req / (Lambda_ell**2 * Delta0))

print("kappa =", kappa)
print("eta =", eta)
print("Delta_0 =", Delta0)
print("Delta_inf =", DeltaInf)
print("Upsilon_fail =", Upsilon_fail)
print("Upsilon_suff =", Upsilon_suff)

V0_fail_sq = sp.simplify(hbar**2 * c_sw**2 * Upsilon_fail / (4 * rho_w**2))
V0_suff_sq = sp.simplify(hbar**2 * c_sw**2 * Upsilon_suff / (4 * rho_w**2))
print("V0_fail^2 =", V0_fail_sq)
print("V0_suff^2 =", V0_suff_sq)

# Shell-gradient dominated asymptotics: kappa ~ (4/5) Lambda_ell^2.
c = sp.Rational(2) / sp.sqrt(5)
Upsilon_fail_shell = sp.simplify(2 * Pe_req / (sp.sqrt(5) * Lambda_ell))
Upsilon_suff_shell = sp.simplify(sp.Rational(4, 5) * (1 + sp.Rational(2) / sp.sqrt(5)) * Pe_req)

# Check against leading asymptotics derived directly from Delta0, DeltaInf with alpha ~ c Lambda_ell.
Delta0_shell = sp.simplify(Lambda_ell / ((c * Lambda_ell)**2 * (c * Lambda_ell + Lambda_ell)))
DeltaInf_shell = sp.simplify((1 + Lambda_ell / (c * Lambda_ell)) / (c * Lambda_ell + Lambda_ell))
expect_zero("shell fail asymptotic", sp.simplify(Pe_req / (Lambda_ell**2 * DeltaInf_shell) - Upsilon_fail_shell))
expect_zero("shell suff asymptotic", sp.simplify(Pe_req / (Lambda_ell**2 * Delta0_shell) - Upsilon_suff_shell))

# Compression dominated asymptotics: kappa ~ 4 chi_s^2.
Delta0_comp = sp.simplify(Lambda_ell / ((2 * chi_s)**2 * (2 * chi_s + Lambda_ell)))
DeltaInf_comp = sp.simplify((1 + Lambda_ell / (2 * chi_s)) / (2 * chi_s + Lambda_ell))
Upsilon_fail_comp = sp.simplify(Pe_req / (Lambda_ell**2 * DeltaInf_comp))
Upsilon_suff_comp = sp.simplify(Pe_req / (Lambda_ell**2 * Delta0_comp))

print("Upsilon_fail_shell =", Upsilon_fail_shell)
print("Upsilon_suff_shell =", Upsilon_suff_shell)
print("Upsilon_fail_comp =", Upsilon_fail_comp)
print("Upsilon_suff_comp =", Upsilon_suff_comp)
expect_zero("compression fail asymptotic", Upsilon_fail_comp - 2 * Pe_req * chi_s / Lambda_ell**2)
expect_zero("compression suff asymptotic",
            Upsilon_suff_comp - 4 * Pe_req * chi_s**2 * (Lambda_ell + 2 * chi_s) / Lambda_ell**3)

banner("STAGE 55 THEOREM LEDGER")
print("Explicit branch fail surface:")
print("  Upsilon_fail = Pe_req / [Lambda_ell^2 Delta_inf(4 chi_s^2 + 4/5 Lambda_ell^2, Lambda_ell)]")
print("Explicit branch success surface:")
print("  Upsilon_suff = Pe_req / [Lambda_ell^2 Delta_0(4 chi_s^2 + 4/5 Lambda_ell^2, Lambda_ell)]")
print("Shell-gradient dominated limit:")
print("  Upsilon_fail ~ 2 Pe_req / (sqrt(5) Lambda_ell)")
print("  Upsilon_suff ~ (4/5)(1 + 2/sqrt(5)) Pe_req")
print("Compression dominated limit:")
print("  Upsilon_fail ~ 2 Pe_req chi_s / Lambda_ell^2")
print("  Upsilon_suff ~ 4 Pe_req chi_s^2 (Lambda_ell + 2 chi_s) / Lambda_ell^3")
