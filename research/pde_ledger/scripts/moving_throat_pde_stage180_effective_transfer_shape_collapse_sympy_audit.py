#!/usr/bin/env python3
"""
moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.py

SymPy-backed audit for Stage 180.

Checks:
1. Exact multi-port collapse from Xi_1 = 2 sum rho_r tau_r to a single effective
   transfer shape T_eff^2 = N_0 / K.
2. Exact one-port continuum formula
      T^2 = beta_0/K_0
          = Z_W (1+rho)^2 / [Omega_W^2 (1-eps_W)^2].
3. Exact selected-branch reformulation
      T^2 = (27 pi^2 G c_s^5 / (20 a^5 c^5)) * (1-eps_eta)/R_target.
4. Exact weak-axisymmetric slope laws in both direct-port and selected-branch variables.
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

banner("STAGE 163 — EFFECTIVE TRANSFER-SHAPE COLLAPSE")

# ---------------------------------------------------------------------------
# 1. Multi-port collapse to one effective transfer shape
# ---------------------------------------------------------------------------

T1, T2, tau1, tau2, eps, lam = sp.symbols("T1 T2 tau1 tau2 eps lam", positive=True, real=True)

Teff2 = T1**2 * sp.exp(2*eps*lam*tau1) + T2**2 * sp.exp(2*eps*lam*tau2)
Xi_eff = sp.simplify(sp.diff(sp.log(Teff2), eps).subs(eps, 0) / lam)

rho1 = sp.simplify(T1**2 / (T1**2 + T2**2))
rho2 = sp.simplify(T2**2 / (T1**2 + T2**2))
Xi_expected = sp.simplify(2 * (rho1*tau1 + rho2*tau2))

expect_zero("multi-port effective-shape identity", Xi_eff - Xi_expected)

# ---------------------------------------------------------------------------
# 2. One-port continuum transfer shape
# ---------------------------------------------------------------------------

banner("One-port continuum transfer shape")

muW, mueta = sp.symbols("mu_W mu_eta", positive=True, real=True)
Keta, KW = sp.symbols("K_eta_eff K_W_eff", positive=True, real=True)
ZW, rho, epsW = sp.symbols("Z_W rho eps_W", real=True)
OmegaW2 = sp.symbols("Omega_W2", positive=True, real=True)

K0 = Keta / mueta
beta0 = (muW / mueta) * (Keta / KW) * ZW * (1 + rho)**2 / (1 - epsW)**2

T2_direct = sp.simplify(beta0 / K0)
T2_expected = sp.simplify((muW / KW) * ZW * (1 + rho)**2 / (1 - epsW)**2)
expect_zero("T^2 = beta0/K0 -> muW/KW form", T2_direct - T2_expected)

expect_zero(
    "T^2 = ZW(1+rho)^2 / [OmegaW^2 (1-epsW)^2]",
    T2_direct - (ZW * (1 + rho)**2) / (OmegaW2 * (1 - epsW)**2)
        .subs({OmegaW2: KW / muW})
)

# ---------------------------------------------------------------------------
# 3. Selected-branch reformulation
# ---------------------------------------------------------------------------

banner("Selected-branch reformulation")

G, cs, a, c = sp.symbols("G c_s a c", positive=True, real=True)
epseta, Rtarget = sp.symbols("eps_eta R_target", real=True)
Lambda = sp.simplify(27 * sp.pi**2 * G * cs**5 * KW / (20 * a**5 * c**5 * muW))
Rtarget_def = sp.simplify(Lambda * (1 - epseta) * (1 - epsW)**2 / (ZW * (1 + rho)**2))

T2_selected = sp.simplify(27 * sp.pi**2 * G * cs**5 / (20 * a**5 * c**5) * (1 - epseta) / Rtarget)
expect_zero(
    "selected-branch T^2 identity",
    T2_direct.subs({Rtarget: Rtarget_def}) - T2_selected.subs({Rtarget: Rtarget_def})
)

# ---------------------------------------------------------------------------
# 4. Weak-axisymmetric slope laws
# ---------------------------------------------------------------------------

banner("Weak-axisymmetric slope laws")

zetaW, omegaW, rho1s, epsW1 = sp.symbols("zeta_W omega_W rho_1 epsW_1", real=True)
e = sp.symbols("e", real=True)

T2_pert = (
    ZW * sp.exp(e*lam*zetaW)
    * (1 + rho + e*lam*rho1s)**2
    / ((KW/muW) * sp.exp(e*lam*omegaW) * (1 - epsW - e*lam*epsW1)**2)
)
Xi_direct = sp.simplify(sp.diff(sp.log(T2_pert), e).subs(e, 0) / lam)
Xi_direct_expected = sp.simplify(
    zetaW - omegaW + 2*rho1s/(1 + rho) + 2*epsW1/(1 - epsW)
)
expect_zero("direct slope law", Xi_direct - Xi_direct_expected)

eta1, R1 = sp.symbols("eta_1 R_1", real=True)
T2_sel_pert = (
    27 * sp.pi**2 * G * cs**5 / (20 * a**5 * c**5)
    * (1 - epseta - e*lam*eta1)
    / (Rtarget * sp.exp(e*lam*R1))
)
Xi_sel = sp.simplify(sp.diff(sp.log(T2_sel_pert), e).subs(e, 0) / lam)
Xi_sel_expected = sp.simplify(-eta1/(1 - epseta) - R1)
expect_zero("selected-branch slope law", Xi_sel - Xi_sel_expected)

print("\nCarry-forward formulas:")
print("  T_eff^2 = sum_r T_r^2 = N_0/K")
print("  Xi_1    = delta ln(T_eff^2)/(eps lambda_A)")
print("  On the one-port continuum branch,")
print("    T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-eps_W)^2]")
print("        = (27 pi^2 G c_s^5 / (20 a^5 c^5)) * (1-eps_eta)/R_target")
print("  Hence")
print("    Xi_1 = zeta_W - omega_W + 2 rho_1/(1+rho) + 2 epsW_1/(1-eps_W)")
print("        = - eta_1/(1-eps_eta) - R_1")
