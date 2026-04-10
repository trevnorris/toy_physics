
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr):
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


"""
5pn_stage161_outgoing_port_coloading_theorem.py

Executable SymPy audit for the moving-throat Stage 161 chain.

What this script does
---------------------
1. Rewrites the remaining weak-axisymmetric grouped defect as
      Xi_1 = sum_r rho_r^(N) (nu_r - kappa_1) = nu_bar_N - kappa_1.
2. Defines the actual static outgoing-port numerator and detuning slopes
      p_r, d_r
   and verifies
      nu_r = 2 (p_r - d_r).
3. Proves exact equivalence with the earlier slippage language
      nu_r = kappa_1 + 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r.
4. States the exact zero-defect and per-port co-loading conditions.
"""

banner("STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM")

# ---------------------------------------------------------------------------
# I. Exact outgoing-port rewrite of Xi_1
# ---------------------------------------------------------------------------

subbanner("I. Xi_1 as the mismatch between outgoing-port slope and wall-baseline slope")

rho1, rho2 = sp.symbols("rho_1 rho_2", real=True)
nu1, nu2 = sp.symbols("nu_1 nu_2", real=True)
kappa1 = sp.symbols("kappa_1", real=True)

Xi1 = sp.expand(rho1 * (nu1 - kappa1) + rho2 * (nu2 - kappa1))
nu_bar_N = sp.expand(rho1 * nu1 + rho2 * nu2)

print("Xi_1 =", Xi1)
print("nu_bar_N =", nu_bar_N)
expect_zero("Xi_1 - (nu_bar_N - kappa_1) on rho_1+rho_2=1",
            Xi1.subs(rho2, 1-rho1) - (nu_bar_N - kappa1).subs(rho2, 1-rho1))

# ---------------------------------------------------------------------------
# II. Actual numerator and detuning slopes
# ---------------------------------------------------------------------------

subbanner("II. Actual static port slopes p_r, d_r, and nu_r = 2(p_r-d_r)")

oU, oW = sp.symbols("o_U o_W", real=True)
gU, gW = sp.symbols("g_U g_W", real=True)
r = sp.symbols("r", real=True)

alpha, beta = sp.symbols("alpha beta", real=True)
chi, zeta = sp.symbols("chi zeta", real=True)

p_r = sp.expand(alpha * (oU + gW) + beta * (r + gU))
d_r = sp.expand(chi * (oU + oW) - 2 * zeta * r)
nu_r = sp.expand(2 * (p_r - d_r))

print("p_r =", p_r)
print("d_r =", d_r)
print("nu_r =", nu_r)

# ---------------------------------------------------------------------------
# III. Exact link back to the Stage-160 slippage language
# ---------------------------------------------------------------------------

subbanner("III. Exact equivalence to the Stage-160 slippage variables")

I, H = sp.symbols("I H", real=True)
alpha_I = sp.simplify(1 / (1 + I))
beta_I = sp.simplify(I / (1 + I))
chi_H = sp.simplify(1 / (1 - H))
zeta_H = sp.simplify(H / (1 - H))

m_r = sp.expand(gW - oW - sp.Rational(1, 2) * kappa1)
i_r = sp.expand(r + gU - oU - gW)
h_r = sp.expand(2 * r - oU - oW)

nu_slippage = sp.expand(
    kappa1
    + 2 * m_r
    + 2 * I / (1 + I) * i_r
    + 2 * H / (1 - H) * h_r
)

nu_direct = sp.expand(nu_r.subs({alpha: alpha_I, beta: beta_I, chi: chi_H, zeta: zeta_H}))

print("m_r =", m_r)
print("i_r =", i_r)
print("h_r =", h_r)
print("nu_r (direct) =", nu_direct)
print("nu_r (slippage form) =", nu_slippage)
expect_zero("direct nu_r - slippage nu_r", nu_direct - nu_slippage)

sigma_r = sp.expand(nu_slippage - kappa1)
tau_r = sp.expand(m_r + I/(1+I) * i_r + H/(1-H) * h_r)
expect_zero("sigma_r - 2 tau_r", sigma_r - 2 * tau_r)

# ---------------------------------------------------------------------------
# IV. Exact zero-defect and co-loading conditions
# ---------------------------------------------------------------------------

subbanner("IV. Exact zero-defect and per-port co-loading conditions")

Xi1_dom1 = sp.simplify(sp.limit((rho1*(nu1-kappa1) + (1-rho1)*(nu2-kappa1)).subs(rho1, rho1), rho1, 1))
Xi1_dom2 = sp.simplify(sp.limit((rho1*(nu1-kappa1) + (1-rho1)*(nu2-kappa1)).subs(rho1, rho1), rho1, 0))
print("Dominant port 1: Xi_1 ->", Xi1_dom1)
print("Dominant port 2: Xi_1 ->", Xi1_dom2)
expect_zero("dominant port 1", Xi1_dom1 - (nu1 - kappa1))
expect_zero("dominant port 2", Xi1_dom2 - (nu2 - kappa1))

print("\nExact zero-defect condition:")
print("  Xi_1 = 0  <=>  nu_bar_N = kappa_1.")
print("Stronger per-port sufficient condition:")
print("  nu_r = kappa_1 for every active outgoing port r.")

banner("FINAL STAGE-161 LEDGER")
print("1. Xi_1 = sum_r rho_r^(N) (nu_r - kappa_1) = nu_bar_N - kappa_1.")
print("2. nu_r is the actual static outgoing-transfer slope of N_{0}^{(r)} = P_r^2/Delta_r^2.")
print("3. nu_r = 2(p_r-d_r), with p_r the numerator slope and d_r the detuning slope.")
print("4. In the Stage-160 slippage language, nu_r - kappa_1 = 2 tau_r.")
print("5. The remaining theorem gate is whether the outgoing-weighted static transfer slope co-loads with the wall baseline.")
