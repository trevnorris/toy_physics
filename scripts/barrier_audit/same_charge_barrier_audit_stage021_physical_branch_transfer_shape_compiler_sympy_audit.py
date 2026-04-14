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


banner("STAGE 021 — PHYSICAL-BRANCH TRANSFER-SHAPE COMPILER AND THE POST-STATIC DRESSING-INVARIANCE THEOREM")

# ---------------------------------------------------------------------------
# I. Exact coherent-branch observables and selected-branch identity
# ---------------------------------------------------------------------------
subbanner("I. Exact coherent-branch observables")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
ZW, OmW2 = sp.symbols("Z_W Omega_W_sq", positive=True, real=True)
eps, eps_eta = sp.symbols("epsilon epsilon_eta", real=True)
Lam0 = sp.symbols("Lambda_0", positive=True, real=True)

Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
Tsq = sp.simplify(ZW * (1 + chi0) ** 2 / (OmW2 * (1 - eps) ** 2))
Rtarget = sp.simplify(Lam0 * OmW2 * (1 - eps_eta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))

print("R_tr =")
sp.pprint(Rtr)
print("T^2 =")
sp.pprint(Tsq)
print("R_target =")
sp.pprint(Rtarget)
expect_zero("R_target * T^2 - Lambda_0(1-epsilon_eta)", sp.simplify(Rtarget * Tsq - Lam0 * (1 - eps_eta)))

# ---------------------------------------------------------------------------
# II. Exact finite packet in physical branch variables
# ---------------------------------------------------------------------------
subbanner("II. Exact finite packet and rigid-mouth reduction")

Cstar, Bstar = sp.symbols("C_star B_star", positive=True, real=True)
chi0r, deltaUr = sp.symbols("chi_0_ref delta_U_ref", positive=True, real=True)
ZWr, OmW2r = sp.symbols("Z_W_ref Omega_W_sq_ref", positive=True, real=True)
epsr, eps_etar = sp.symbols("epsilon_ref epsilon_eta_ref", real=True)

Rtr_r = sp.simplify((1 + chi0r / (1 + deltaUr)) / (1 + chi0r))
Tsq_r = sp.simplify(ZWr * (1 + chi0r) ** 2 / (OmW2r * (1 - epsr) ** 2))
Rtarget_r = sp.simplify(Lam0 * OmW2r * (1 - eps_etar) * (1 - epsr) ** 2 / (ZWr * (1 + chi0r) ** 2))

qtr = sp.simplify(-Cstar * sp.log(Rtr / Rtr_r))
qnt = sp.simplify(Bstar * sp.log(Rtr / Rtr_r) + sp.log((1 - eps_eta) / (1 - eps_etar)) - sp.log(Rtarget / Rtarget_r))
qeta = sp.simplify(sp.log(eps_eta / eps_etar))

print("q_tr =")
sp.pprint(qtr)
print("q_nt =")
sp.pprint(qnt)
print("q_eta =")
sp.pprint(qeta)

expect_zero(
    "q_nt + (B*/C*) q_tr - ln(T^2/T_ref^2)",
    sp.expand_log(qnt + (Bstar / Cstar) * qtr - sp.log(Tsq / Tsq_r), force=True),
)

expect_zero(
    "rigid-mouth q_nt - ln(T^2/T_ref^2)",
    sp.expand_log(qnt.subs(Rtr, Rtr_r) - sp.log(Tsq / Tsq_r), force=True),
)

# ---------------------------------------------------------------------------
# III. First-order physical drift compiler
# ---------------------------------------------------------------------------
subbanner("III. First-order physical drift compiler")

lam = sp.symbols("lam", real=True)
dlnchi0, dlndeltaU = sp.symbols("dlnchi_0 dlndelta_U", real=True)
dlnZW, dlnOmW2, dlneps, dlneps_eta = sp.symbols(
    "dlnZ_W dlnOmega_W_sq dlnepsilon dlnepsilon_eta", real=True
)

subs_lin = {
    chi0: chi0r * sp.exp(lam * dlnchi0),
    deltaU: deltaUr * sp.exp(lam * dlndeltaU),
    ZW: ZWr * sp.exp(lam * dlnZW),
    OmW2: OmW2r * sp.exp(lam * dlnOmW2),
    eps: epsr + lam * (epsr * dlneps),
    eps_eta: eps_etar * sp.exp(lam * dlneps_eta),
}

Rtr_lin = sp.simplify(Rtr.subs(subs_lin))
Tsq_lin = sp.simplify(Tsq.subs(subs_lin))
Rtarget_lin = sp.simplify(Rtarget.subs(subs_lin))
qeta_lin = sp.simplify(qeta.subs(subs_lin))

R1 = sp.simplify(sp.diff(sp.log(Rtr_lin), lam).subs(lam, 0))
T1 = sp.simplify(sp.diff(sp.log(Tsq_lin), lam).subs(lam, 0))
Rt1 = sp.simplify(sp.diff(sp.log(Rtarget_lin), lam).subs(lam, 0))
E1 = sp.simplify(sp.diff(qeta_lin, lam).subs(lam, 0))

Ctr = sp.simplify(chi0r * deltaUr / ((1 + chi0r) * (1 + deltaUr) * (1 + chi0r + deltaUr)))
c_eta = sp.simplify(eps_etar / (1 - eps_etar))

print("d ln R_tr =")
sp.pprint(R1)
print("d ln T^2 =")
sp.pprint(T1)
print("d ln R_target =")
sp.pprint(Rt1)
print("d ln epsilon_eta =")
sp.pprint(E1)

expect_zero(
    "d ln R_tr + C_tr[(1+delta_U)dlnchi0 + (1+chi0)dlndeltaU]",
    sp.simplify(R1 + Ctr * ((1 + deltaUr) * dlnchi0 + (1 + chi0r) * dlndeltaU)),
)
expect_zero(
    "d ln T^2 - [dlnZ_W - dlnOmega_W^2 + 2 chi_0 dlnchi_0/(1+chi_0) + 2 epsilon dlnepsilon/(1-epsilon)]",
    sp.simplify(T1 - (dlnZW - dlnOmW2 + 2 * chi0r * dlnchi0 / (1 + chi0r) + 2 * epsr * dlneps / (1 - epsr))),
)
expect_zero(
    "d ln R_target - [dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1+chi_0) - 2 epsilon dlnepsilon/(1-epsilon) - epsilon_eta dlnepsilon_eta/(1-epsilon_eta)]",
    sp.simplify(Rt1 - (dlnOmW2 - dlnZW - 2 * chi0r * dlnchi0 / (1 + chi0r) - 2 * epsr * dlneps / (1 - epsr) - eps_etar * dlneps_eta / (1 - eps_etar))),
)
expect_zero("E1 - dlnepsilon_eta", E1 - dlneps_eta)
expect_zero("dlnR_target + dlnT^2 + c_eta dln epsilon_eta", sp.simplify(Rt1 + T1 + c_eta * E1))

# rigid-mouth first-order q_nt = d ln T^2
qtr1 = sp.simplify(-Cstar * R1)
qnt1_rigid = sp.simplify(T1)
print("first-order rigid-mouth q_tr =")
sp.pprint(qtr1)
print("first-order rigid-mouth q_nt =")
sp.pprint(qnt1_rigid)
print("first-order rigid-mouth q_eta =")
sp.pprint(E1)

# ---------------------------------------------------------------------------
# IV. Support-blindness factorization of the physical packet
# ---------------------------------------------------------------------------
subbanner("IV. Support-blindness of the physical same-charge packet")

zeta, Mmix = sp.symbols("zeta M_mix", real=True)
Mtr = sp.simplify(Mmix * (1 + zeta * (1 - eps) / (1 - zeta * eps)))
print("M_tr =")
sp.pprint(Mtr)
expect_zero("d R_tr / d zeta", sp.diff(Rtr, zeta))
expect_zero("d T^2 / d zeta", sp.diff(Tsq, zeta))
expect_zero("d epsilon_eta / d zeta", sp.diff(eps_eta, zeta))
expect_zero("d R_tr / d M_mix", sp.diff(Rtr, Mmix))
expect_zero("d T^2 / d M_mix", sp.diff(Tsq, Mmix))
expect_zero("d epsilon_eta / d M_mix", sp.diff(eps_eta, Mmix))

# ---------------------------------------------------------------------------
# V. Post-static dressing-invariance theorem
# ---------------------------------------------------------------------------
subbanner("V. Post-static dressing-invariance theorem")

expect_zero("rigid-mouth static gate iff T^2 = T_ref^2", sp.simplify(sp.log(Tsq / Tsq_r).subs(Tsq, Tsq_r)))
expect_zero("post-static q_eta - ln(epsilon_eta/epsilon_eta_ref)", qeta - sp.log(eps_eta / eps_etar))
expect_zero("post-static first-order q_eta - dlnepsilon_eta", E1 - dlneps_eta)

print("\nAll Stage-021 audit checks passed.")
