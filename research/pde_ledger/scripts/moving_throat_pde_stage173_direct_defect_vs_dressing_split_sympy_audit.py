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


banner("STAGE 173 — DIRECT DEFECT VS DRESSING SPLIT, SUPPORT-BLINDNESS, AND THE SCALAR NO-GO FILTERS")

# ---------------------------------------------------------------------------
# I. Support-blindness of the direct transfer shape and corrected composite
# ---------------------------------------------------------------------------
subbanner("I. Exact support-blindness of T^2, R_target, and N_*")

zeta = sp.symbols("zeta", real=True)
ZW, chi0, OmW2, eps, epseta, Lambda0 = sp.symbols(
    "Z_W chi0 OmegaW2 epsilon epsiloneta Lambda_0", positive=True, real=True
)
deltaU = sp.symbols("deltaU", positive=True, real=True)

T2 = sp.simplify(ZW * (1 + chi0) ** 2 / (OmW2 * (1 - eps) ** 2))
Rtarget = sp.simplify(Lambda0 * OmW2 * (1 - epseta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))
Rtr = sp.simplify((1 + chi0 + deltaU) / ((1 + chi0) * (1 + deltaU)))
Bstar = sp.simplify(2 * (1 + chi0 + deltaU) / deltaU)
Nstar = sp.simplify(T2 * Rtr**Bstar)
Mmix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - epseta) * (1 - eps)))
Msupp = sp.simplify(8 * zeta * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - epseta) * (1 - zeta * eps)))
Ssupport = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
Mtr = sp.simplify(Mmix + Msupp)
Rtarget_support = sp.simplify(8 * Lambda0 * OmW2 * (1 - eps) / sp.pi**2 * Ssupport / Mtr)
T2_support = sp.simplify(Lambda0 * (1 - epseta) / Rtarget_support)
Rtr_sb = sp.Function("Rtr_sb")(zeta)
support_blind_subs = {sp.diff(Rtr_sb, zeta): 0}
Nstar_support = sp.simplify(T2_support * Rtr_sb**Bstar)
Ecomp = sp.simplify(Rtarget_support * T2_support / Lambda0)

print("T^2 =")
sp.pprint(T2)
print("R_target =")
sp.pprint(Rtarget)
print("N_* =")
sp.pprint(Nstar)
print("E = R_target T^2 / Lambda_0 =")
sp.pprint(Ecomp)
expect_zero("support-loaded T^2 reconstruction", T2_support - T2)
expect_zero("support-loaded R_target reconstruction", Rtarget_support - Rtarget)
expect_zero("support-loaded N_* reconstruction", Nstar_support.subs(Rtr_sb, Rtr) - Nstar)
expect_zero("E - (1 - epsiloneta)", Ecomp - (1 - epseta))
expect_zero("dln(T^2_loaded)/dzeta", sp.diff(sp.log(T2_support), zeta))
expect_zero("dln(R_target_loaded)/dzeta", sp.diff(sp.log(Rtarget_support), zeta))
expect_zero(
    "dln(N_* loaded)/dzeta",
    sp.diff(sp.log(Nstar_support), zeta).subs(support_blind_subs),
)

bad = sp.symbols("bad", nonzero=True, real=True)
Msupp_spoiled = sp.simplify(Msupp + bad * zeta * Mmix)
Rtarget_spoiled = sp.simplify(8 * Lambda0 * OmW2 * (1 - eps) / sp.pi**2 * Ssupport / (Mmix + Msupp_spoiled))
dlnR_spoiled = sp.simplify(sp.diff(sp.log(Rtarget_spoiled), zeta).subs(bad, 1))
print("spoiled dln(R_target)/dzeta =")
sp.pprint(dlnR_spoiled)
if dlnR_spoiled == 0:
    raise AssertionError("Expected a spoiled support packet to break the Stage-173 blind route.")

# ---------------------------------------------------------------------------
# II. Microscopic slippage ledger from primitive multiplicative drifts
# ---------------------------------------------------------------------------
subbanner("II. Microscopic slippage ledger from primitive drifts")

s = sp.symbols("s", real=True)
lamW, muW, gamma, cetaU, KU, Keta, KW, TU, sigma, L = sp.symbols(
    "lambdaW muW gamma cetaU K_U K_eta K_W T_U sigma L", positive=True, real=True
)
lambda1, mu1, gamma1, c1, kappaU, kappaEta, kappaW, tau1 = sp.symbols(
    "lambda_1 mu_1 gamma_1 c_1 kappa_U kappa_eta kappa_W tau_1", real=True
)

subs_mult = {
    lamW: lamW * sp.exp(s * lambda1),
    muW: muW * sp.exp(s * mu1),
    gamma: gamma * sp.exp(s * gamma1),
    cetaU: cetaU * sp.exp(s * c1),
    KU: KU * sp.exp(s * kappaU),
    Keta: Keta * sp.exp(s * kappaEta),
    KW: KW * sp.exp(s * kappaW),
    TU: TU * sp.exp(s * tau1),
}

chi0_expr = sp.simplify(gamma * cetaU / KU)
deltaU_expr = sp.simplify(sp.pi**2 * TU / (L**2 * KU))
epseta_expr = sp.simplify(cetaU**2 / (KU * Keta))
epsW_expr = sp.simplify(gamma**2 * lamW**2 * sigma / (KU * KW))
Zratio_expr = sp.simplify(lamW**2 * muW / (Keta * KW**2))


def log_drift(expr: sp.Expr) -> sp.Expr:
    expr_s = sp.simplify(expr.subs(subs_mult))
    return sp.simplify(sp.diff(sp.log(expr_s), s).subs(s, 0))

Sigma_chi_expr = log_drift(chi0_expr)
Sigma_delta_expr = log_drift(deltaU_expr)
Sigma_eta_expr = log_drift(epseta_expr)
Sigma_eps_expr = log_drift(epsW_expr)
Sigma_Z_expr = log_drift(Zratio_expr)

print("Sigma_chi =")
sp.pprint(Sigma_chi_expr)
print("Sigma_delta =")
sp.pprint(Sigma_delta_expr)
print("Sigma_eta =")
sp.pprint(Sigma_eta_expr)
print("Sigma_epsilon =")
sp.pprint(Sigma_eps_expr)
print("Sigma_Z =")
sp.pprint(Sigma_Z_expr)

expect_zero("Sigma_chi - (gamma1 + c1 - kappaU)", Sigma_chi_expr - (gamma1 + c1 - kappaU))
expect_zero("Sigma_delta - (tau1 - kappaU)", Sigma_delta_expr - (tau1 - kappaU))
expect_zero("Sigma_eta - (2 c1 - kappaU - kappaEta)", Sigma_eta_expr - (2 * c1 - kappaU - kappaEta))
expect_zero("Sigma_epsilon - (2 gamma1 + 2 lambda1 - kappaU - kappaW)", Sigma_eps_expr - (2 * gamma1 + 2 * lambda1 - kappaU - kappaW))
expect_zero("Sigma_Z - (2 lambda1 + mu1 - kappaEta - 2 kappaW)", Sigma_Z_expr - (2 * lambda1 + mu1 - kappaEta - 2 * kappaW))

# ---------------------------------------------------------------------------
# III. Exact direct-defect decomposition into tracking + nontracking pieces
# ---------------------------------------------------------------------------
subbanner("III. Exact direct-defect decomposition")

Sigma_Z, Sigma_chi, Sigma_eps, Sigma_delta, Sigma_eta = sp.symbols(
    "Sigma_Z Sigma_chi Sigma_epsilon Sigma_delta Sigma_eta", real=True
)
epsW = sp.symbols("epsilon_W", positive=True, real=True)

Xi_direct = sp.simplify(
    Sigma_Z
    + 2 * chi0 / (1 + chi0) * Sigma_chi
    + 2 * epsW / (1 - eps)
    * (
        (11 + 9 * deltaU) / (11 * (1 + deltaU)) * Sigma_eps
        - 2 * deltaU / (11 * (1 + deltaU) ** 2) * Sigma_delta
    )
)
Sigma_tr = sp.simplify((1 + chi0) * Sigma_delta + (1 + deltaU) * Sigma_chi)
Atr = sp.simplify(2 * chi0 / ((1 + chi0) * (1 + deltaU)))
Ctr = sp.simplify(
    chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
)
Sigma_nt = sp.simplify(
    Sigma_Z
    + 2 * epsW / (1 - eps) * (11 + 9 * deltaU) / (11 * (1 + deltaU)) * Sigma_eps
    - (
        2 * chi0 / (1 + deltaU)
        + 4 * epsW * deltaU / (11 * (1 - eps) * (1 + deltaU) ** 2)
    )
    * Sigma_delta
)

print("Xi_1 (direct microscopic law) =")
sp.pprint(Xi_direct)
print("Sigma_tr =")
sp.pprint(Sigma_tr)
print("Sigma_nt =")
sp.pprint(Sigma_nt)
print("A_tr =")
sp.pprint(Atr)
print("C_tr =")
sp.pprint(Ctr)
expect_zero("Xi_direct - (A_tr Sigma_tr + Sigma_nt)", Xi_direct - (Atr * Sigma_tr + Sigma_nt))

# ---------------------------------------------------------------------------
# IV. Exact direct-defect / dressing split and inverse reconstruction
# ---------------------------------------------------------------------------
subbanner("IV. Exact triangular direct-defect / dressing compiler")

Str, Snt, Seta = sp.symbols("Sigma_tr Sigma_nt Sigma_eta_ba", real=True)
Theta = sp.simplify(-Ctr * Str)
Xi = sp.simplify(Atr * Str + Snt)
RplusXi = sp.simplify(-epseta / (1 - epseta) * Seta)
Rcal = sp.simplify(RplusXi - Xi)

compiler = sp.Matrix(
    [
        [-Ctr, 0, 0],
        [Atr, 1, 0],
        [0, 0, -epseta / (1 - epseta)],
    ]
)
print("compiler[(Sigma_tr,Sigma_nt,Sigma_eta) -> (Theta, Xi, R+Xi)] =")
sp.pprint(compiler)
print("det(compiler) =")
sp.pprint(sp.simplify(sp.factor(compiler.det())))
expect_zero("dXi/dSigma_eta", sp.diff(Xi, Seta))
expect_zero("d(R+Xi)/dSigma_tr", sp.diff(RplusXi, Str))
expect_zero("d(R+Xi)/dSigma_nt", sp.diff(RplusXi, Snt))

Sigma_tr_rec = sp.simplify(-((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) / (chi0 * deltaU) * Theta)
Sigma_nt_rec = sp.simplify(Xi + Atr / Ctr * Theta)
Sigma_eta_rec = sp.simplify(-(1 - epseta) / epseta * (Rcal + Xi))

expect_zero("Sigma_tr reconstruction", Sigma_tr_rec - Str)
expect_zero("A_tr/C_tr - 2(1+chi0+deltaU)/deltaU", sp.simplify(Atr / Ctr - 2 * (1 + chi0 + deltaU) / deltaU))
expect_zero("Sigma_nt reconstruction", Sigma_nt_rec - Snt)
expect_zero("Sigma_eta reconstruction", Sigma_eta_rec - Seta)
expect_zero("tracking-rigid Xi - Sigma_nt", Xi.subs(Str, 0) - Snt)
expect_zero("grouped-defect cancellation theorem", Xi.subs(Snt, -Atr * Str))
expect_zero("selected-branch rigidity theorem", RplusXi.subs(Seta, 0))

# ---------------------------------------------------------------------------
# V. Pure grouped real P2 anisotropy has no linear scalar feed-down
# ---------------------------------------------------------------------------
subbanner("V. Weak-axisymmetric grouped signature and the scalar no-go filter")

epsax, x0, x1, y0, y1 = sp.symbols("epsilon_ax x0 x1 y0 y1", real=True)

x20 = x0 + epsax * x1
x21 = x0 + sp.Rational(1, 2) * epsax * x1
x22 = x0 - epsax * x1

y20 = y0 + epsax * y1
y21 = y0 + sp.Rational(1, 2) * epsax * y1
y22 = y0 - epsax * y1

xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
ax = sp.simplify((2 * x20 - x21 - x22) / 10)
bx = sp.simplify((x21 - x22) / 2)

ybar = sp.simplify((y20 + 2 * y21 + 2 * y22) / 5)
ay = sp.simplify((2 * y20 - y21 - y22) / 10)
by = sp.simplify((y21 - y22) / 2)

Ixy = sp.simplify(4 * ax * ay + sp.Rational(4, 5) * bx * by)
Ixx = sp.simplify(4 * ax**2 + sp.Rational(4, 5) * bx**2)

expect_zero("xbar - x0", xbar - x0)
expect_zero("ybar - y0", ybar - y0)
expect_zero("b_x - 3 a_x", bx - 3 * ax)
expect_zero("b_y - 3 a_y", by - 3 * ay)
expect_zero("I[x,y] - 7/10 eps^2 x1 y1", Ixy - sp.Rational(7, 10) * epsax**2 * x1 * y1)
expect_zero("I[x,x] - 7/10 eps^2 x1^2", Ixx - sp.Rational(7, 10) * epsax**2 * x1**2)
expect_zero("linear term in xbar", sp.diff(xbar, epsax).subs(epsax, 0))
expect_zero("linear term in I[x,x] at eps=0", sp.diff(Ixx, epsax).subs(epsax, 0))

banner("STAGE 173 LEDGER")
print("1. The coherent first-order problem splits exactly into:")
print("      direct packet   (Sigma_tr, Sigma_nt)")
print("      dressing packet Sigma_eta")
print("   with an exact triangular compiler to (Theta_1, Xi_1, R_1+Xi_1).")
print("2. The direct transfer shape T^2, the selected-branch demand ratio R_target,")
print("   and the corrected nontracking composite N_* are all support-blind.")
print("3. The direct grouped defect Xi_1 depends only on the direct slippages")
print("      (Sigma_Z, Sigma_chi, Sigma_epsilon, Sigma_delta),")
print("   while the selected-branch residual adds only Sigma_eta.")
print("4. Pure grouped real P2 anisotropy has no linear scalar feed-down:")
print("   the grouped weighted trace has no O(epsilon) term and the first scalar")
print("   invariant is quadratic, I ~ (7/10) epsilon^2.")
print("5. So the remaining first-order gate is already sharply classified as tracking,")
print("   direct transfer-shape, or dressing mismatch.")
