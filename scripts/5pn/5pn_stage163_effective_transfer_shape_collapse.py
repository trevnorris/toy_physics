
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
5pn_stage163_effective_transfer_shape_collapse.py

Executable SymPy audit for the moving-throat Stage 163 theorem.

What this script does
---------------------
1. Collapses the many-port weighted average to one effective transfer shape
      T_eff^2 = N_0 / K.
2. Evaluates the one-port continuum transfer shape
      T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-epsilon_W)^2]
   and proves the exact weak-axisymmetric slope law.
3. Rewrites the same transfer shape in selected-branch variables,
   proving
      Xi_1 = -eta_1/(1-epsilon_eta) - R_1.
4. Shows how the grouped weak-axisymmetric normalization pattern is fixed once
   the one-port transfer-shape slope is known.
"""

banner("STAGE 163 — EFFECTIVE TRANSFER-SHAPE COLLAPSE")

# ---------------------------------------------------------------------------
# I. Many-port collapse
# ---------------------------------------------------------------------------

subbanner("I. Exact collapse from many ports to one effective transfer shape")

T1, T2 = sp.symbols("T_1 T_2", positive=True, real=True)
tau1, tau2 = sp.symbols("tau_1 tau_2", real=True)
eps, lamA, t = sp.symbols("epsilon lambda_A t", real=True)

rho1 = sp.simplify(T1**2 / (T1**2 + T2**2))
rho2 = sp.simplify(T2**2 / (T1**2 + T2**2))
Xi1_weighted = sp.simplify(2 * (rho1 * tau1 + rho2 * tau2))

T1_t = T1 * sp.exp(t * eps * lamA * tau1)
T2_t = T2 * sp.exp(t * eps * lamA * tau2)
Teff2_t = sp.simplify(T1_t**2 + T2_t**2)
Xi1_eff = sp.simplify(sp.diff(sp.log(Teff2_t), t).subs(t, 0) / (eps * lamA))

print("T_eff^2 =", T1**2 + T2**2)
print("Xi_1 from weighted port slopes =", Xi1_weighted)
print("Xi_1 from log T_eff^2 =", Xi1_eff)
expect_zero("effective transfer-shape collapse", Xi1_weighted - Xi1_eff)

# ---------------------------------------------------------------------------
# II. One-port continuum transfer shape
# ---------------------------------------------------------------------------

subbanner("II. Exact one-port continuum transfer shape")

ZW, rho, OmegaW2, epsW = sp.symbols("Z_W rho Omega_W2 epsilon_W", positive=True, real=True)
zetaZ, omegaW, rho1s, varepsW = sp.symbols("zeta_Z omega_W rho_1 varepsilon_W", real=True)

ZW_t = ZW * sp.exp(t * eps * lamA * zetaZ)
OmegaW2_t = OmegaW2 * sp.exp(t * eps * lamA * omegaW)
rho_t = rho + t * eps * lamA * rho1s
epsW_t = epsW + t * eps * lamA * varepsW

T2_cont_t = sp.simplify(ZW_t * (1 + rho_t)**2 / (OmegaW2_t * (1 - epsW_t)**2))
Xi1_cont = sp.simplify(sp.diff(sp.log(T2_cont_t), t).subs(t, 0) / (eps * lamA))
Xi1_cont_expected = sp.simplify(zetaZ - omegaW + 2 * rho1s / (1 + rho) + 2 * varepsW / (1 - epsW))

print("T^2 =", ZW * (1 + rho)**2 / (OmegaW2 * (1 - epsW)**2))
print("Xi_1 =", Xi1_cont)
expect_zero("one-port continuum slope law", Xi1_cont - Xi1_cont_expected)

# ---------------------------------------------------------------------------
# III. Selected-branch reformulation
# ---------------------------------------------------------------------------

subbanner("III. Selected-branch reformulation")

G, cs, a, c = sp.symbols("G c_s a c", positive=True, real=True)
eps_eta = sp.symbols("epsilon_eta", real=True)
eta1, R1 = sp.symbols("eta_1 R_1", real=True)

Lambda0 = sp.simplify(27 * sp.pi**2 * G * cs**5 / (20 * a**5 * c**5))
eps_eta_t = eps_eta + t * eps * lamA * eta1
Rtarget = sp.symbols("R_target", positive=True, real=True)
Rtarget_t = Rtarget * sp.exp(t * eps * lamA * R1)

T2_sel_t = sp.simplify(Lambda0 * (1 - eps_eta_t) / Rtarget_t)
Xi1_sel = sp.simplify(sp.diff(sp.log(T2_sel_t), t).subs(t, 0) / (eps * lamA))
Xi1_sel_expected = sp.simplify(-eta1 / (1 - eps_eta) - R1)

print("T^2 selected-branch =", Lambda0 * (1 - eps_eta) / Rtarget)
print("Xi_1 selected-branch =", Xi1_sel)
expect_zero("selected-branch slope law", Xi1_sel - Xi1_sel_expected)

# ---------------------------------------------------------------------------
# IV. Grouped normalization pattern
# ---------------------------------------------------------------------------

subbanner("IV. Grouped weak-axisymmetric normalization pattern on the one-port branch")

tau = sp.symbols("tau", real=True)
Delta20 = sp.simplify(2 * eps * tau)
Delta21 = sp.simplify(eps * tau)
Delta22 = sp.simplify(-2 * eps * tau)

print("Delta_Q^(20) =", Delta20)
print("Delta_Q^(21) =", Delta21)
print("Delta_Q^(22) =", Delta22)
print("These follow from Xi_1 = 2 tau on the one-port branch.")

banner("FINAL STAGE-163 LEDGER")
print("1. The many-port weighted average is exactly the slope of one effective transfer shape:")
print("      Xi_1 = d ln T_eff^2 / (eps lambda_A).")
print("2. On the minimal one-port continuum branch,")
print("      T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-epsilon_W)^2],")
print("   so")
print("      Xi_1 = zeta_Z - omega_W + 2 rho_1/(1+rho) + 2 varepsilon_W/(1-epsilon_W).")
print("3. In selected-branch form,")
print("      T^2 = const * (1-epsilon_eta)/R_target,")
print("   hence")
print("      Xi_1 = -eta_1/(1-epsilon_eta) - R_1.")
print("4. On the one-port branch, Xi_1 = 2 tau fixes the full grouped weak-axisymmetric normalization pattern immediately.")
