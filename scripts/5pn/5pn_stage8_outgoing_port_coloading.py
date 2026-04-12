
#!/usr/bin/env python3
"""
5pn_stage8_outgoing_port_coloading.py

Eighth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Rewrites the remaining weak-axisymmetric grouped defect directly in terms of
   the actual outgoing-port numerator and detuning slopes.
2. Proves the exact Stage-161 identity
      Xi_1 = sum_r rho_r^(N) (nu_r - kappa_1) = \bar{nu}_N - kappa_1,
   where nu_r is the logarithmic slope of N_{0}^{(r)} = P_r^2 / Delta_r^2.
3. Derives the direct portwise formulas for the numerator slope p_r and the
   detuning slope d_r, and verifies that they are exactly equivalent to the
   Stage-160 slippage variables (m_r, i_r, h_r).
4. Introduces the wall-normalized transfer shape T_r by N_{0}^{(r)} = K T_r^2
   and proves the exact Stage-162 collapse
      Xi_1 = 2 sum_r rho_r^(N) tau_r.
5. Records the per-port co-loading theorem, dominant-port limit, and the
   recovery of the square-root mixed-leg law as a special case.

Interpretation
--------------
After this stage the remaining grouped weak-axisymmetric defect is not carried by
every microscopic outgoing slippage separately. It is the mismatch between the
outgoing-weighted static transfer slope and the wall-baseline slope, or
equivalently twice the outgoing-weighted transfer-shape slope.
"""

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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
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


# ---------------------------------------------------------------------------
# I. Direct outgoing-port rewrite of Xi_1
# ---------------------------------------------------------------------------

banner("I. DIRECT OUTGOING-PORT REWRITE OF THE REMAINING DEFECT")

rho1, rho2 = sp.symbols("rho_1 rho_2", real=True)
kappa1 = sp.symbols("kappa_1", real=True)
nu1, nu2 = sp.symbols("nu_1 nu_2", real=True)
nu_bar = sp.symbols("nu_bar", real=True)

Xi1 = sp.simplify(rho1 * (nu1 - kappa1) + rho2 * (nu2 - kappa1))
Xi1_bar = sp.simplify((rho1 * nu1 + rho2 * nu2) - kappa1)

print("Xi_1 =", Xi1)
print("nu_bar - kappa_1 =", Xi1_bar)

subbanner("Normalized outgoing weights")
Xi1_norm = sp.simplify(Xi1.subs(rho2, 1 - rho1))
Xi1_bar_norm = sp.simplify(Xi1_bar.subs(rho2, 1 - rho1))
expect_zero("Stage-161 identity on rho_1 + rho_2 = 1", Xi1_norm - Xi1_bar_norm)


# ---------------------------------------------------------------------------
# II. Actual outgoing-port numerator and detuning slopes
# ---------------------------------------------------------------------------

banner("II. ACTUAL PORT SLOPES p_r, d_r, nu_r")

# Portwise weak-axisymmetric slopes.
oU, oW = sp.symbols("o_U o_W", real=True)
gU, gW = sp.symbols("g_U g_W", real=True)
r = sp.symbols("r", real=True)

# Stage-161 convex weights from the static numerator and detuning.
alpha, beta = sp.symbols("alpha beta", real=True)
chi, zeta = sp.symbols("chi zeta", real=True)

p = sp.simplify(alpha * (oU + gW) + beta * (r + gU))
d = sp.simplify(chi * (oU + oW) - 2 * zeta * r)
nu = sp.simplify(2 * (p - d))

print("p_r =", p)
print("d_r =", d)
print("nu_r =", nu)

subbanner("Rewrite in the Stage-160 slippage language")
I, H = sp.symbols("I H", real=True)
alpha_I = sp.simplify(1 / (1 + I))
beta_I = sp.simplify(I / (1 + I))
chi_H = sp.simplify(1 / (1 - H))
zeta_H = sp.simplify(H / (1 - H))

m = sp.simplify(gW - oW - kappa1 / 2)
i = sp.simplify(r + gU - oU - gW)
h = sp.simplify(2 * r - oU - oW)

nu_slippage = sp.simplify(
    kappa1
    + 2 * m
    + 2 * I / (1 + I) * i
    + 2 * H / (1 - H) * h
)
nu_direct = sp.simplify(nu.subs({alpha: alpha_I, beta: beta_I, chi: chi_H, zeta: zeta_H}))

print("m_r =", m)
print("i_r =", i)
print("h_r =", h)
print("nu_r in slippage form =", nu_slippage)
expect_zero("direct port-slope formula == Stage-160 slippage formula", nu_direct - nu_slippage)

tau = sp.simplify(m + I / (1 + I) * i + H / (1 - H) * h)
sigma = sp.simplify(nu_slippage - kappa1)
expect_zero("sigma_r - 2 tau_r", sigma - 2 * tau)


# ---------------------------------------------------------------------------
# III. Transfer-shape factorization
# ---------------------------------------------------------------------------

banner("III. WALL-NORMALIZED TRANSFER SHAPES")

T1, T2 = sp.symbols("T_1 T_2", positive=True, real=True)
tau1, tau2 = sp.symbols("tau_1 tau_2", real=True)

# With N0^(r) = K T_r^2, the outgoing weights are exactly T_r^2 / sum_s T_s^2.
rho1_T = sp.simplify(T1**2 / (T1**2 + T2**2))
rho2_T = sp.simplify(T2**2 / (T1**2 + T2**2))

Xi1_T = sp.simplify(2 * (rho1_T * tau1 + rho2_T * tau2))
print("rho_1^(N) =", rho1_T)
print("rho_2^(N) =", rho2_T)
print("Xi_1 from transfer-shape slopes =", Xi1_T)

# Effective transfer shape T_eff^2 = T1^2 + T2^2.
eps, lamA = sp.symbols("epsilon lambda_A", real=True)
t = sp.symbols("t", real=True)

T1_t = T1 * sp.exp(t * eps * lamA * tau1)
T2_t = T2 * sp.exp(t * eps * lamA * tau2)
Teff2_t = sp.simplify(T1_t**2 + T2_t**2)
Xi1_eff = sp.simplify(sp.diff(sp.log(Teff2_t), t).subs(t, 0) / (eps * lamA))

print("T_eff^2 =", T1**2 + T2**2)
print("delta ln T_eff^2 / (eps lambda_A) =", Xi1_eff)
expect_zero("Stage-162 effective-transfer identity", Xi1_eff - Xi1_T)

subbanner("Per-port co-loading theorem")
print("Exact zero-defect condition: Xi_1 = 0 iff weighted transfer-shape slope vanishes.")
expect_zero("per-port rigidity tau_1=tau_2=0 => Xi_1=0", Xi1_T.subs({tau1: 0, tau2: 0}))

subbanner("Dominant-port limit")
Xi1_dom1 = sp.simplify(sp.limit(Xi1_T, T2, 0, dir="+"))
Xi1_dom2 = sp.simplify(sp.limit(Xi1_T, T1, 0, dir="+"))
print("If port 1 dominates, Xi_1 ->", Xi1_dom1)
print("If port 2 dominates, Xi_1 ->", Xi1_dom2)
expect_zero("dominant port 1", Xi1_dom1 - 2 * tau1)
expect_zero("dominant port 2", Xi1_dom2 - 2 * tau2)


# ---------------------------------------------------------------------------
# IV. Recovery of the square-root mixed-leg law
# ---------------------------------------------------------------------------

banner("IV. RECOVERY OF THE SQUARE-ROOT MIXED-LEG LAW")

u_r, c_r = sp.symbols("u_r c_r", real=True)

# Under normalized upper-leg rigidity and normalized coupling rigidity,
# tau_r collapses to the mixed-leg drift w_r = m_r.
tau_rigid = sp.simplify(tau.subs({i: 0, h: 0}))
expect_zero("rigid-branch tau_r - m_r", tau_rigid - m)

print("Rigid branch (i_r = h_r = 0): tau_r =", tau_rigid)
print("So Xi_1 is carried only by the mixed-leg slope w_r = g_W - o_W - kappa_1/2.")
print("That is exactly the square-root mixed-leg law corridor.")


# ---------------------------------------------------------------------------
# V. Final theorem ledger
# ---------------------------------------------------------------------------

banner("V. FINAL THEOREM LEDGER")
print("1. The remaining grouped weak-axisymmetric defect is")
print("      Xi_1 = sum_r rho_r^(N) (nu_r - kappa_1) = nu_bar - kappa_1.")
print("2. The actual static outgoing-transfer slope is")
print("      nu_r = 2(p_r - d_r).")
print("3. In the Stage-160 slippage language,")
print("      nu_r - kappa_1 = 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r = 2 tau_r.")
print("4. Therefore")
print("      Xi_1 = 2 sum_r rho_r^(N) tau_r.")
print("5. The many-port weighted sum collapses exactly to one effective transfer-shape slope")
print("      Xi_1 = delta ln T_eff^2 / (eps lambda_A).")
print("6. A strong sufficient condition for zero defect is per-port co-loading:")
print("      tau_r = 0 for every active port r.")
print("7. In the dominant-port limit, Xi_1 is just twice that port's transfer-shape slope.")
print("8. Under upper-leg and coupling rigidity, tau_r collapses to the square-root")
print("   mixed-leg law corridor carried only by the raw mixed leg.")
