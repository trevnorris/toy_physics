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


def expect_zero(name: str, expr) -> None:
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


banner("STAGE 022 — RIGID-MOUTH PHYSICAL NORMAL FORM AND THE CARTESIAN ORBIT-LOCK THEOREM")

# ---------------------------------------------------------------------------
# I. Exact rigid-mouth physical logarithmic chart
# ---------------------------------------------------------------------------
subbanner("I. Exact rigid-mouth physical logarithmic chart")

Lam0 = sp.symbols("Lambda_0", positive=True, real=True)
Tsq_ref, eps_eta_ref = sp.symbols("T_sq_ref epsilon_eta_ref", positive=True, real=True)
U, V = sp.symbols("U V", real=True)

Tsq = sp.simplify(Tsq_ref * sp.exp(U))
eps_eta = sp.simplify(eps_eta_ref * sp.exp(V))
Rtarget = sp.simplify(Lam0 * (1 - eps_eta) / Tsq)
Rtarget_ref = sp.simplify(Lam0 * (1 - eps_eta_ref) / Tsq_ref)
Rratio = sp.simplify(Rtarget / Rtarget_ref)

print("T^2 =")
sp.pprint(Tsq)
print("epsilon_eta =")
sp.pprint(eps_eta)
print("R_target/R_target_ref =")
sp.pprint(Rratio)
expect_zero(
    "R_target/R_ref - exp(-U)*(1-epsilon_eta_ref*exp(V))/(1-epsilon_eta_ref)",
    Rratio - sp.exp(-U) * (1 - eps_eta_ref * sp.exp(V)) / (1 - eps_eta_ref),
)

q_phys = sp.Matrix([U, V])
print("rigid-mouth packet q_phys =")
sp.pprint(q_phys)

# ---------------------------------------------------------------------------
# II. Exact physical projectors and commuting finite legs
# ---------------------------------------------------------------------------
subbanner("II. Exact physical projectors")

P_T = sp.Matrix([[1, 0], [0, 0]])
P_eta = sp.Matrix([[0, 0], [0, 1]])

expect_zero("P_T^2 - P_T", P_T * P_T - P_T)
expect_zero("P_eta^2 - P_eta", P_eta * P_eta - P_eta)
expect_zero("P_T P_eta", P_T * P_eta)
expect_zero("P_eta P_T", P_eta * P_T)
expect_zero("P_T + P_eta - I_2", P_T + P_eta - sp.eye(2))

q_T = P_T * q_phys
q_eta_axis = P_eta * q_phys
print("pure transfer-shape axis =")
sp.pprint(q_T)
print("pure dressing axis =")
sp.pprint(q_eta_axis)
expect_zero("q_phys - q_T - q_eta_axis", q_phys - q_T - q_eta_axis)

Rratio_T = sp.simplify(sp.exp(-U))
Rratio_eta = sp.simplify((1 - eps_eta_ref * sp.exp(V)) / (1 - eps_eta_ref))
expect_zero("Rratio - Rratio_T * Rratio_eta", Rratio - Rratio_T * Rratio_eta)

# ---------------------------------------------------------------------------
# III. Exact physical-to-microscopic dependent-plane compiler
# ---------------------------------------------------------------------------
subbanner("III. Exact physical-to-microscopic dependent-plane compiler")

C_phys_dep = sp.Matrix(
    [
        [0, 0],
        [0, -1],
        [1, -1],
    ]
)
L_phys_dep = sp.Matrix(
    [
        [0, -1, 1],
        [0, -1, 0],
    ]
)

y_dep = C_phys_dep * q_phys
print("y_dep(U,V) =")
sp.pprint(y_dep)
expect_zero("L_phys_dep * C_phys_dep - I_2", L_phys_dep * C_phys_dep - sp.eye(2))
expect_zero("y_dep - (0,-V,U-V)", y_dep - sp.Matrix([0, -V, U - V]))

P_T_dep = sp.simplify(C_phys_dep * P_T * L_phys_dep)
P_eta_dep = sp.simplify(C_phys_dep * P_eta * L_phys_dep)
print("P_T^(dep) =")
sp.pprint(P_T_dep)
print("P_eta^(dep) =")
sp.pprint(P_eta_dep)

expect_zero("(P_T_dep)^2 - P_T_dep", P_T_dep * P_T_dep - P_T_dep)
expect_zero("(P_eta_dep)^2 - P_eta_dep", P_eta_dep * P_eta_dep - P_eta_dep)
expect_zero("P_T_dep P_eta_dep", P_T_dep * P_eta_dep)
expect_zero("P_eta_dep P_T_dep", P_eta_dep * P_T_dep)
expect_zero(
    "P_T_dep + P_eta_dep - identity on Delta_T=0 plane",
    P_T_dep + P_eta_dep - sp.Matrix([[0, 0, 0], [0, 1, 0], [0, 0, 1]]),
)

y_T = sp.simplify(C_phys_dep * q_T)
y_eta = sp.simplify(C_phys_dep * q_eta_axis)
print("y_T =")
sp.pprint(y_T)
print("y_eta =")
sp.pprint(y_eta)
expect_zero("y_T - (0,0,U)", y_T - sp.Matrix([0, 0, U]))
expect_zero("y_eta + V*(0,1,1)", y_eta + V * sp.Matrix([0, 1, 1]))
expect_zero("y_dep - y_T - y_eta", y_dep - y_T - y_eta)
expect_zero("packet recovery U - (Delta_mu-Delta_Keta)", L_phys_dep * y_dep - q_phys)

# ---------------------------------------------------------------------------
# IV. Exact correction compilers
# ---------------------------------------------------------------------------
subbanner("IV. Exact correction compilers")

delta_y_static = sp.simplify(-y_T)
delta_y_eta_rest = sp.simplify(V * sp.Matrix([0, 1, 1]))
delta_y_orbit = sp.simplify(-y_dep)

print("Delta y_static =")
sp.pprint(delta_y_static)
print("Delta y_eta_rest =")
sp.pprint(delta_y_eta_rest)
print("Delta y_orbit =")
sp.pprint(delta_y_orbit)

expect_zero("y_dep + Delta y_static - y_eta", y_dep + delta_y_static - y_eta)
expect_zero("Delta y_static + Delta y_eta_rest - Delta y_orbit", delta_y_static + delta_y_eta_rest - delta_y_orbit)
expect_zero("y_dep + Delta y_orbit", y_dep + delta_y_orbit)

# ---------------------------------------------------------------------------
# V. Exact support-blindness of the physical normal form
# ---------------------------------------------------------------------------
subbanner("V. Exact support-blindness")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
ZW, OmW2 = sp.symbols("Z_W Omega_W_sq", positive=True, real=True)
eps = sp.symbols("epsilon", real=True)
zeta, Mmix = sp.symbols("zeta M_mix", real=True)

Tsq_direct = sp.simplify(ZW * (1 + chi0) ** 2 / (OmW2 * (1 - eps) ** 2))
Mtr = sp.simplify(Mmix * (1 + zeta * (1 - eps) / (1 - zeta * eps)))
U_direct = sp.simplify(sp.log(Tsq_direct / Tsq_ref))
V_direct = sp.simplify(sp.log(eps_eta / eps_eta_ref))

print("M_tr =")
sp.pprint(Mtr)
expect_zero("d T^2 / d zeta", sp.diff(Tsq_direct, zeta))
expect_zero("d T^2 / d M_mix", sp.diff(Tsq_direct, Mmix))
expect_zero("d epsilon_eta / d zeta", sp.diff(eps_eta, zeta))
expect_zero("d epsilon_eta / d M_mix", sp.diff(eps_eta, Mmix))
expect_zero("d U / d zeta", sp.diff(U_direct, zeta))
expect_zero("d U / d M_mix", sp.diff(U_direct, Mmix))
expect_zero("d V / d zeta", sp.diff(V_direct, zeta))
expect_zero("d V / d M_mix", sp.diff(V_direct, Mmix))
expect_zero("d y_dep / d zeta", sp.diff(C_phys_dep * sp.Matrix([U_direct, V_direct]), zeta))
expect_zero("d y_dep / d M_mix", sp.diff(C_phys_dep * sp.Matrix([U_direct, V_direct]), Mmix))

# ---------------------------------------------------------------------------
# VI. First-order rigid-mouth form
# ---------------------------------------------------------------------------
subbanner("VI. First-order rigid-mouth form")

lam = sp.symbols("lam", real=True)
U1, V1 = sp.symbols("U_1 V_1", real=True)
Tsq_lin = sp.simplify(Tsq_ref * sp.exp(lam * U1))
eps_eta_lin = sp.simplify(eps_eta_ref * sp.exp(lam * V1))
U_lin = sp.simplify(sp.diff(sp.log(Tsq_lin / Tsq_ref), lam).subs(lam, 0))
V_lin = sp.simplify(sp.diff(sp.log(eps_eta_lin / eps_eta_ref), lam).subs(lam, 0))
y_lin = sp.simplify(C_phys_dep * sp.Matrix([U_lin, V_lin]))

print("first-order U =")
sp.pprint(U_lin)
print("first-order V =")
sp.pprint(V_lin)
print("first-order y_dep =")
sp.pprint(y_lin)
expect_zero("first-order U - U_1", U_lin - U1)
expect_zero("first-order V - V_1", V_lin - V1)
expect_zero("first-order y_dep - (0,-V_1,U_1-V_1)", y_lin - sp.Matrix([0, -V1, U1 - V1]))

print("\nAll Stage-022 audit checks passed.")
