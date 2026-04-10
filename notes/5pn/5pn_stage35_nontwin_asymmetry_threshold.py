#!/usr/bin/env python3
"""
5pn_stage35_nontwin_asymmetry_threshold.py

SymPy audit for Moving-Throat PDE Stage 35.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 35 — EXACT NON-TWIN ASYMMETRY REQUIREMENT")

# Symbols
Pi_tr, C_mix = sp.symbols("Pi_tr C_mix", positive=True, real=True)
eps = sp.symbols("eps", positive=True, real=True)
Omega0, K_phi_eff, K_W_eff = sp.symbols("Omega_0 K_phi_eff K_W_eff", positive=True, real=True)

# Required coherent support ratio
S_req = sp.simplify(Pi_tr / C_mix)
zeta_req = sp.simplify((S_req - 1) / (1 + eps * (S_req - 2)))

print("S_req    =", S_req)
print("zeta_req =", zeta_req)

dzeta = sp.simplify(sp.diff(zeta_req, Pi_tr))
dzeta_expected = sp.simplify(C_mix * (1 - eps) / (C_mix - eps * (2 * C_mix - Pi_tr)) ** 2)
expect_zero("dzeta_req/dPi_tr - expected", dzeta - dzeta_expected)

expect_zero("zeta_req(Pi_tr=C_mix)", sp.simplify(zeta_req.subs(Pi_tr, C_mix)))
expect_zero("zeta_req(Pi_tr=2C_mix) - 1", sp.simplify(zeta_req.subs(Pi_tr, 2 * C_mix) - 1))

Delta_zeta = sp.simplify(zeta_req - 1)
Delta_expected = sp.simplify((1 - eps) * (Pi_tr - 2 * C_mix) / (C_mix - eps * (2 * C_mix - Pi_tr)))
expect_zero("Delta_zeta - expected", Delta_zeta - Delta_expected)

print("\nDelta_zeta =", Delta_zeta)

# Exact physical lowest-lane ratio
zeta_phys = sp.simplify((K_W_eff / K_phi_eff) * Omega0**2)
print("\nzeta_0^(phys) =", zeta_phys)

Omega0_req_sq = sp.simplify(zeta_req * K_phi_eff / K_W_eff)
K_phi_req = sp.simplify(K_W_eff * Omega0**2 / zeta_req)

print("Omega_(0,req)^2 =", Omega0_req_sq)
print("K_(phi,0)^(req) =", K_phi_req)

# Pure-overlap and pure-softening diagnostics
Omega_equal_stiff = sp.simplify(sp.sqrt(zeta_req))
Kphi_equal_overlap = sp.simplify(K_W_eff / zeta_req)
softening_fraction = sp.simplify(1 - Kphi_equal_overlap / K_W_eff)
softening_fraction_expected = sp.simplify((1 - eps) * (Pi_tr - 2 * C_mix) / (Pi_tr - C_mix))

print("\nPure-overlap rescue at equal stiffness:")
print("  Omega_0 >= sqrt(zeta_req) =", Omega_equal_stiff)
print("Pure-softening rescue at equal overlap:")
print("  K_phi <= K_W / zeta_req   =", Kphi_equal_overlap)
expect_zero("equal-overlap softening fraction", softening_fraction - softening_fraction_expected)

banner("STAGE 35 THEOREM LEDGER")
print("Exact regime split:")
print("  Pi_tr <= C_mix         : mixed-only enough")
print("  C_mix < Pi_tr <= 2C_mix: symmetric lowest twin enough")
print("  Pi_tr > 2C_mix         : non-twin asymmetry required")
print()
print("Once Pi_tr > 2C_mix, the exact required lowest-lane support ratio is")
print("  zeta_req = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ].")
print()
print("The two equivalent rescue thresholds are")
print("  Omega_0^2 >= zeta_req K_phi^(eff) / K_W^(eff),")
print("  K_phi^(eff) <= K_W^(eff) Omega_0^2 / zeta_req.")
