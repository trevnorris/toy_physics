#!/usr/bin/env python3
"""
moving_throat_pde_stage57_family1_healing_lock_sympy_audit.py

SymPy audit for Stage 074:
- use the exact GNLS healing/compliance width ell = hbar/(2 m c_s),
- verify chi_s = Lambda_ell/2,
- verify kappa = (9/5) Lambda_ell^2 on the explicit Family-1 branch,
- evaluate the reference-branch values.
"""

from __future__ import annotations
import sympy as sp

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

banner("STAGE 074 — HEALING-LENGTH LOCK AND SUPPORT SCALE")

Lambda_ell = sp.symbols("Lambda_ell", positive=True)
hbar, m_psi, c_s, ell, L = sp.symbols("hbar m_psi c_s ell L", positive=True)

# Physical definition of the dimensionless support scale chi_s = m c_s L / hbar.
chi_def = m_psi * c_s * L / hbar

# Apply the GNLS healing/compliance width: ell = hbar / (2 m c_s),
# equivalently c_s = hbar / (2 m_psi ell).
chi_in_ell = sp.simplify(chi_def.subs(c_s, hbar / (2 * m_psi * ell)))
print("chi (after healing-length substitution) =", chi_in_ell)

# Re-express L/ell as the dimensionless ratio Lambda_ell.
chi_lock = sp.simplify(chi_in_ell.subs(L, Lambda_ell * ell))

# Family-1 branch coefficient: kappa = 4 chi_s^2 + (4/5) Lambda_ell^2.
# Coefficients 4 and 4/5 come from the Family-1 Euler-Lagrange branch
# (carried forward from the earlier Family-1 stages); this stage only
# verifies that, with chi_s locked to Lambda_ell/2, kappa reduces to
# (9/5) Lambda_ell^2.
kappa_lock = sp.simplify(4 * chi_lock**2 + sp.Rational(4, 5) * Lambda_ell**2)

print("chi_s (locked) =", chi_lock)
print("kappa(Lambda_ell) =", kappa_lock)

expect_zero("chi_s - Lambda_ell/2", chi_lock - Lambda_ell / 2)
expect_zero("kappa - (9/5) Lambda_ell^2", kappa_lock - sp.Rational(9, 5) * Lambda_ell**2)

Lambda_ref = sp.Integer(37)
chi_ref = sp.simplify(chi_lock.subs(Lambda_ell, Lambda_ref))
kappa_ref = sp.simplify(kappa_lock.subs(Lambda_ell, Lambda_ref))
alpha_ref = sp.simplify(sp.sqrt(kappa_ref))

print("\nReference branch:")
print("Lambda_ell =", Lambda_ref)
print("chi_s      =", chi_ref)
print("kappa      =", kappa_ref)
print("alpha      =", alpha_ref)
print("alpha (numeric) =", sp.N(alpha_ref, 20))

expect_zero("chi_ref - 37/2", chi_ref - sp.Rational(37, 2))
expect_zero("kappa_ref - 12321/5", kappa_ref - sp.Rational(12321, 5))
expect_zero("alpha_ref - 111/sqrt(5)", alpha_ref - sp.Rational(111) / sp.sqrt(5))

print("\nFinal ledger:")
print("  chi_s = 37/2")
print("  kappa = 12321/5")
