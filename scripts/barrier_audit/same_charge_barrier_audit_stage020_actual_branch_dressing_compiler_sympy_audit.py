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


banner("STAGE 020 — ACTUAL-BRANCH DRESSING COMPILER AND THE FINITE STATIC-BLIND CURVE")

# ---------------------------------------------------------------------------
# I. Finite rigid-mouth packet in the direct actual-branch observables
# ---------------------------------------------------------------------------
subbanner("I. Finite rigid-mouth packet and the exact static-blind curve")

eps_ref, qeta = sp.symbols("epsilon_eta_ref q_eta", positive=True, real=True)
Rratio = sp.symbols("R_target_ratio", positive=True, real=True)

eps = sp.symbols("epsilon_eta", positive=True, real=True)
qnt = sp.log((1 - eps) / (1 - eps_ref)) - sp.log(Rratio)
qeta_from_eps = sp.log(eps / eps_ref)
print("q_nt(epsilon_eta, R_target/R_target,ref) =")
sp.pprint(qnt)
print("q_eta(epsilon_eta) =")
sp.pprint(qeta_from_eps)

Rratio_static = sp.simplify((1 - eps) / (1 - eps_ref))
print("static-blind finite curve: R_target/R_target,ref =")
sp.pprint(Rratio_static)
expect_zero("q_nt on static curve", qnt.subs(Rratio, Rratio_static))

eps_of_qeta = sp.simplify(eps_ref * sp.exp(qeta))
Rratio_of_qeta = sp.simplify((1 - eps_ref * sp.exp(qeta)) / (1 - eps_ref))
print("epsilon_eta(q_eta) =")
sp.pprint(eps_of_qeta)
print("R_target/R_target,ref as a function of q_eta =")
sp.pprint(Rratio_of_qeta)
expect_zero(
    "q_nt(q_eta) on finite static curve",
    sp.simplify(
        qnt.subs({eps: eps_of_qeta, Rratio: Rratio_of_qeta})
    ),
)

qeta_from_Rratio = sp.simplify(sp.log((1 - (1 - eps_ref) * Rratio) / eps_ref))
print("q_eta recovered from R_target/R_target,ref once q_nt = 0 =")
sp.pprint(qeta_from_Rratio)
expect_zero(
    "q_eta(Rratio_static) - q_eta",
    sp.simplify(qeta_from_Rratio.subs(Rratio, Rratio_of_qeta) - qeta),
)

# ---------------------------------------------------------------------------
# II. Linearization and the tangent of the finite static-blind curve
# ---------------------------------------------------------------------------
subbanner("II. Linearization and the tangent relation")

c_eta = sp.simplify(eps_ref / (1 - eps_ref))
print("c_eta =")
sp.pprint(c_eta)

# derivative of ln(Rratio(q_eta)) at q_eta = 0
slope = sp.simplify(sp.diff(sp.log(Rratio_of_qeta), qeta).subs(qeta, 0))
print("d/dq_eta ln(R_target/R_target,ref) at q_eta = 0 =")
sp.pprint(slope)
expect_zero("slope + c_eta", sp.simplify(slope + c_eta))

R1, E1 = sp.symbols("R_1 E_1", real=True)
qnt_lin = -R1 - c_eta * E1
print("first-order rigid-mouth packet:")
print("q_nt =")
sp.pprint(qnt_lin)
print("q_eta = E_1")
expect_zero(
    "post-static tangent relation R_1 + c_eta q_eta",
    sp.simplify((R1 + c_eta * E1).subs(R1, -c_eta * E1)),
)
qeta_post_static = sp.simplify(-R1 / c_eta)
print("q_eta recovered from R_1 once q_nt = 0 =")
sp.pprint(qeta_post_static)
expect_zero(
    "q_eta(post-static) - [-(1-eps_ref)/eps_ref R_1]",
    sp.simplify(qeta_post_static + R1 * (1 - eps_ref) / eps_ref),
)

# ---------------------------------------------------------------------------
# III. Exact microscopic compiler for q_eta
# ---------------------------------------------------------------------------
subbanner("III. Exact microscopic compiler for q_eta")

c_etaU, c_etaU_ref = sp.symbols("c_etaU c_etaU_ref", positive=True, real=True)
K_U, K_U_ref = sp.symbols("K_U K_U_ref", positive=True, real=True)
Keta_eff, Keta_eff_ref = sp.symbols("K_eta_eff K_eta_eff_ref", positive=True, real=True)

eps_micro = sp.simplify(c_etaU**2 / (K_U * Keta_eff))
eps_micro_ref = sp.simplify(c_etaU_ref**2 / (K_U_ref * Keta_eff_ref))
qeta_micro = sp.simplify(sp.log(eps_micro / eps_micro_ref))
qeta_micro_split = sp.simplify(
    2 * sp.log(c_etaU / c_etaU_ref)
    - sp.log(K_U / K_U_ref)
    - sp.log(Keta_eff / Keta_eff_ref)
)
print("epsilon_eta =")
sp.pprint(eps_micro)
print("q_eta microscopic finite ratio =")
sp.pprint(qeta_micro)
expect_zero("q_eta finite split formula", qeta_micro - qeta_micro_split)

# First-order microscopic drift extractor
lam = sp.symbols("lam", real=True)
c1, kappaU, kappaEta = sp.symbols("c_1 kappa_U kappa_eta", real=True)
subs_lin = {
    c_etaU: c_etaU_ref * sp.exp(lam * c1),
    K_U: K_U_ref * sp.exp(lam * kappaU),
    Keta_eff: Keta_eff_ref * sp.exp(lam * kappaEta),
}
qeta_micro_lin = sp.simplify(sp.expand(qeta_micro.subs(subs_lin)))
print("q_eta microscopic linearized exactly =")
sp.pprint(qeta_micro_lin)
coeff_lam = sp.simplify(sp.diff(qeta_micro_lin, lam).subs(lam, 0))
print("first-order microscopic drift coefficient =")
sp.pprint(coeff_lam)
expect_zero("microscopic drift law - (2 c1 - kappaU - kappaEta)", coeff_lam - (2 * c1 - kappaU - kappaEta))

# ---------------------------------------------------------------------------
# IV. Support-blindness of q_eta
# ---------------------------------------------------------------------------
subbanner("IV. Support-blindness")

lambda_phi, Kphi_eff = sp.symbols("lambda_phi K_phi_eff", positive=True, real=True)
lambda_W, KW_eff = sp.symbols("lambda_W K_W_eff", positive=True, real=True)
Mmix = sp.symbols("M_mix", positive=True, real=True)
eps_blk = sp.symbols("epsilon", real=True)
zeta_expr = sp.simplify(lambda_phi**2 * KW_eff / (lambda_W**2 * Kphi_eff))
zeta_sym = sp.symbols('zeta', real=True)
Mtr = sp.simplify(Mmix * (1 + zeta_sym * (1 - eps_blk) / (1 - zeta_sym * eps_blk)))
print("zeta =")
sp.pprint(zeta_expr)
print("M_tr =")
sp.pprint(Mtr)
expect_zero("d q_eta / d lambda_phi", sp.diff(qeta_micro, lambda_phi))
expect_zero("d q_eta / d K_phi_eff", sp.diff(qeta_micro, Kphi_eff))
# Direct observable q_eta is also blind to the support variables.
qeta_obs = sp.log(eps / eps_ref)
expect_zero("d q_eta(obs) / d zeta", sp.diff(qeta_obs, zeta_sym))
expect_zero("d q_eta(obs) / d M_mix", sp.diff(qeta_obs, Mmix))

# ---------------------------------------------------------------------------
# V. Post-static orbit-lock theorem on the actual branch
# ---------------------------------------------------------------------------
subbanner("V. Post-static orbit-lock theorem")

# On the static-blind curve, q_eta = 0 iff epsilon_eta = epsilon_ref iff Rratio = 1.
expect_zero("Rratio(q_eta=0) - 1", sp.simplify(Rratio_of_qeta.subs(qeta, 0) - 1))
expect_zero(
    "epsilon_eta(q_eta=0) - epsilon_eta_ref",
    sp.simplify(eps_of_qeta.subs(qeta, 0) - eps_ref),
)

# Linear dependent-plane carrier from Stage 019 for the post-static residue.
y_eta = sp.Matrix([0, -qeta, -qeta])
print("post-static dependent-plane carrier y_eta =")
sp.pprint(y_eta)
carrier_norm_sq = sp.simplify((y_eta.T * y_eta)[0])
print("||y_eta||^2 =")
sp.pprint(carrier_norm_sq)
expect_zero("||y_eta||^2 - 2 q_eta^2", carrier_norm_sq - 2 * qeta**2)

print("\nAll Stage-020 audit checks passed.")
