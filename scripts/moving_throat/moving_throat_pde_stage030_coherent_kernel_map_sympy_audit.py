#!/usr/bin/env python3
"""
Stage 30 SymPy audit.

Checks:
1. Coherent-kernel identities rho_0 = sigma_0 = chi_0.
2. Exact proportionality Z_phi = zeta Z_W and eps_phi = zeta eps_W.
3. Exact total baseline formula M_tr = M_mix * S(zeta;eps).
4. Exact product law R_target * M_tr.
5. Monotonicity of the support-enhancement factor S in zeta.
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


banner("STAGE 30 — COHERENT KERNEL MAP AUDIT")

# Microscopic coherent-kernel symbols.
Keta_eff, KU, KW_eff, Kphi_eff = sp.symbols("Keta_eff KU KW_eff Kphi_eff", positive=True, real=True)
lamW, lamphi, gamma, c_etaU = sp.symbols("lamW lamphi gamma c_etaU", positive=True, real=True)
Tw, L, TU = sp.symbols("T_w L T_U", positive=True, real=True)
G, cs, a, c, muW = sp.symbols("G c_s a c mu_W", positive=True, real=True)
sigma = sp.Rational(88, 9) / sp.pi ** 2

# Exact coherent identities.
rho0 = sp.simplify((gamma * lamW) * c_etaU / (KU * lamW))
sigma0 = sp.simplify(c_etaU * (gamma * lamphi) / (KU * lamphi))
chi0 = sp.simplify(gamma * c_etaU / KU)

banner("1. Coherent interference ratios")
print("rho_0 =", rho0)
print("sigma_0 =", sigma0)
print("chi_0 =", chi0)
expect_zero("rho_0 - chi_0", rho0 - chi0)
expect_zero("sigma_0 - chi_0", sigma0 - chi0)

# Dimensionless coherent ratios.
eps_eta = sp.simplify(c_etaU ** 2 / (KU * Keta_eff))
epsW = sp.simplify(gamma ** 2 * lamW ** 2 * sigma / (KU * KW_eff))
ZW = sp.simplify(lamW ** 2 / (Keta_eff * KW_eff))
zeta = sp.symbols("zeta", positive=True, real=True)
zeta_def = sp.simplify(lamphi ** 2 * KW_eff / (lamW ** 2 * Kphi_eff))
delta0 = sp.simplify(sp.pi ** 2 * Tw / (L ** 2 * Keta_eff))
deltaU = sp.simplify(sp.pi ** 2 * TU / (L ** 2 * KU))
Lambda = sp.simplify(27 * sp.pi ** 2 * G * cs ** 5 * KW_eff / (20 * a ** 5 * c ** 5 * muW))

eps_phi = sp.simplify(gamma ** 2 * lamphi ** 2 * sigma / (KU * Kphi_eff))
Zphi = sp.simplify(lamphi ** 2 / (Keta_eff * Kphi_eff))

banner("2. Support/mixed proportionality")
print("zeta (symbolic) =", zeta)
print("zeta_def (microscopic) =", zeta_def)
print("eps_W =", epsW)
print("eps_phi =", eps_phi)
print("Z_W =", ZW)
print("Z_phi =", Zphi)
expect_zero("eps_phi - zeta_def*eps_W", eps_phi - zeta_def * epsW)
expect_zero("Z_phi - zeta_def*Z_W", Zphi - zeta_def * ZW)

# Split quantities.
eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
delta = sp.simplify((delta0 + eps_eta * deltaU / (1 + deltaU)) / (1 - eps_eta))

banner("3. Tracking factor and split quantities")
print("R_tr =", Rtr)
print("eps = eps_W^(split) =", eps)
print("delta =", delta)

# Baselines and support-enhancement factor.
Mmix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (sp.pi ** 2 * (1 - eps_eta) * (1 - eps)))
Msupp = sp.simplify(8 * zeta * ZW * (1 + chi0) ** 2 / (sp.pi ** 2 * (1 - eps_eta) * (1 - zeta * eps)))
S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
Mtr = sp.simplify(Mmix + Msupp)
Rtarget = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))

banner("4. Exact coherent baseline formulas")
print("M_mix =", Mmix)
print("M_supp =", Msupp)
print("S(zeta;eps) =", S)
print("M_tr =", Mtr)
print("R_target =", Rtarget)
expect_zero("M_tr - M_mix*S", Mtr - Mmix * S)

# Product law.
product_expected = sp.simplify(8 * Lambda * (1 - eps) / sp.pi ** 2 * S)
product_actual = sp.simplify(Rtarget * Mtr)

banner("5. Product law and monotonicity")
print("R_target * M_tr =", product_actual)
expect_zero("product law", product_actual - product_expected)
expect_zero("dR_target/dzeta", sp.diff(Rtarget, zeta))
expect_zero("dS/dzeta - (1-eps)/(1-zeta eps)^2", sp.diff(S, zeta) - (1 - eps) / (1 - zeta * eps) ** 2)
expect_zero("S(zeta=0)-1", S.subs(zeta, 0) - 1)

print("\nAll Stage-30 symbolic checks passed.")
