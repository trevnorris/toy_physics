#!/usr/bin/env python3
"""
5pn_stage96_outlet_model_status.py

Stage 96 audit: status update across the explicit outlet-model classes.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 96 — OUTLET-MODEL STATUS UPDATE")

rho_R, sigma_W, gamma_W = sp.symbols("rho_R sigma_W gamma_W", real=True)

chi_R = sp.simplify(3 / (3 - rho_R))
chi_mix = sp.simplify(3 * (1 - 9 * sigma_W * gamma_W) / (3 + sigma_W))
chi_hyb = sp.simplify((1 - 9 * sigma_W * gamma_W) / (1 - sigma_W))

print("chi_Q^R   =", sp.factor(chi_R))
print("chi_Q^mix =", sp.factor(chi_mix))
print("chi_Q^hyb =", sp.factor(chi_hyb))

# Pure Robin preserves the canonical outgoing normalization only trivially.
sol_robin = sp.solve(sp.Eq(chi_R, 1), rho_R)
print("Robin-preserving rho_R branches =", sol_robin)
expect_zero("Robin canonical branch is trivial", sol_robin[0])

# Standalone mixed branch: Stage 94 showed canonical-even preservation forces sigma_W = 0.
chi_mix_series = sp.series(chi_mix, sigma_W, 0, 2).removeO()
print("small-sigma_W mixed expansion =", chi_mix_series)
expect_zero(
    "mixed linear shift formula",
    chi_mix_series - (1 - sigma_W * (sp.Rational(1, 3) + 9 * gamma_W)),
)
chi_mix_at_zero = sp.simplify(chi_mix.subs(sigma_W, 0))
expect_zero("mixed branch at sigma_W=0 gives chi_Q=1", chi_mix_at_zero - 1)

# Hybrid branch: exact nontrivial canonical-even surface exists, and gamma_W=1/9 preserves odd normalization.
gamma_pres = sp.solve(sp.Eq(chi_hyb, 1), gamma_W)[0]
print("hybrid-preserving gamma_W =", gamma_pres)
expect_zero("hybrid gamma_W preservation value", gamma_pres - sp.Rational(1, 9))

banner("STAGE 96 FINAL LEDGER")
print("Explicit outlet audit summary:")
print("  1. Pure Robin loading generically shifts chi_Q and is harmless only when rho_R=0.")
print("  2. A standalone isotropic mixed pole cannot sit on the canonical even branch unless it vanishes.")
print("  3. The hybrid Robin–mixed outlet admits one nontrivial compensated canonical-even branch,")
print("     and its odd normalization is preserved exactly when gamma_W = 1/9.")
