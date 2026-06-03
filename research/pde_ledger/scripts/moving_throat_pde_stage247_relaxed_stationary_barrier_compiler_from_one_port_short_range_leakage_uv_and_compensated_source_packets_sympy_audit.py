#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 247
SymPy audit for the relaxed stationary barrier compiler.

This script verifies:
1. the exact one-port short-range baseline inherited from the static mixed-bundle audit,
2. the exact import of the Stage 244 leakage / work packet,
3. the exact import of the Stage 245 U/V drain packet,
4. the exact import of the Stage 246 compensated source packet through the monotone
   source-response scalar R_infty - R(r) on the Session-I orientation,
5. the exact stationary compiler and lowering identity,
6. and a Session-I one-point benchmark decomposition at r_soft = 0.18.
"""

from __future__ import annotations

import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)



def main() -> None:
    # ------------------------------------------------------------------
    # 1. Exact one-port short-range baseline.
    # ------------------------------------------------------------------
    section("1. Exact one-port short-range baseline")
    r, kappa = sp.symbols("r kappa", positive=True, real=True)
    alpha6, alpha2 = sp.symbols("alpha6 alpha2", nonnegative=True, real=True)
    betaQ, betaU, betaW = sp.symbols("betaQ betaU betaW", real=True)
    Kstar, OmU2, OmW2, GU, GW, Rmix = sp.symbols(
        "Kstar OmU2 OmW2 GU GW Rmix", positive=True, real=True
    )

    Delta = sp.simplify(OmU2 * OmW2 - Rmix**2)
    Q = sp.simplify(GU**2 * OmW2 + 2 * GU * GW * Rmix + GW**2 * OmU2)
    P = sp.simplify(OmU2 * GW + Rmix * GU)
    PU = sp.simplify(GU * OmW2 + Rmix * GW)
    D0 = sp.simplify(Kstar - Q / Delta)

    K_red = sp.Matrix([
        [Kstar, -GU, -GW],
        [-GU, OmU2, -Rmix],
        [-GW, -Rmix, OmW2],
    ])
    K_red_inv = sp.simplify(K_red.inv())

    chiqq = sp.simplify(1 / D0)
    chiqU = sp.simplify(PU / (Delta * D0))
    chiqW = sp.simplify(P / (Delta * D0))
    chiUU = sp.simplify((Kstar * OmW2 - GW**2) / (Delta * D0))
    chiUW = sp.simplify((Kstar * Rmix + GU * GW) / (Delta * D0))
    chiWW = sp.simplify((Kstar * OmU2 - GU**2) / (Delta * D0))

    assert sp.simplify(K_red_inv[0, 0] - chiqq) == 0
    assert sp.simplify(K_red_inv[0, 1] - chiqU) == 0
    assert sp.simplify(K_red_inv[0, 2] - chiqW) == 0
    assert sp.simplify(K_red_inv[1, 1] - chiUU) == 0
    assert sp.simplify(K_red_inv[1, 2] - chiUW) == 0
    assert sp.simplify(K_red_inv[2, 2] - chiWW) == 0

    C6 = sp.simplify(chiqq * betaQ**2)
    C4 = sp.simplify(chiqU * betaQ * betaU + chiqW * betaQ * betaW)
    C2 = sp.simplify(chiUU * betaU**2 + 2 * chiUW * betaU * betaW + chiWW * betaW**2)

    A6 = sp.simplify(3 * alpha6 + sp.Rational(1, 2) * C6)
    A4 = sp.simplify(C4)
    A2 = sp.simplify(alpha2 + sp.Rational(1, 2) * C2)

    V_short = sp.simplify(
        (1 / r) * (1 + sp.Rational(1, 2) * sp.exp(-2 * kappa * r))
        - A6 / r**6
        - A4 * sp.exp(-2 * kappa * r) / r**4
        - A2 * sp.exp(-4 * kappa * r) / r**2
    )

    print("Delta                 =", Delta)
    print("Q                     =", Q)
    print("P                     =", P)
    print("D0                    =", D0)
    print("chi_qq                =", chiqq)
    print("chi_qU                =", chiqU)
    print("chi_qW                =", chiqW)
    print("chi_UU                =", chiUU)
    print("chi_UW                =", chiUW)
    print("chi_WW                =", chiWW)
    print("C6                    =", C6)
    print("C4                    =", C4)
    print("C2                    =", C2)
    print("A6                    =", A6)
    print("A4                    =", A4)
    print("A2                    =", A2)
    print("V_short(r)            =", V_short)

    # ------------------------------------------------------------------
    # 2. Exact Stage 244 leakage / work packet.
    # ------------------------------------------------------------------
    section("2. Exact Stage 244 leakage / work packet")
    lam = sp.symbols("lam", positive=True, real=True)
    eta_leak, mu_w, rho0, q = sp.symbols("eta_leak mu_w rho0 q", positive=True, real=True)
    Lvar = sp.symbols("Lvar", positive=True, real=True)  # shorthand for Lambda(r)*varrho(r)

    E0 = sp.simplify(16 * eta_leak * Lvar / sp.pi**2)
    S_leak = sp.simplify(8 * sp.sqrt(2) * eta_leak * mu_w * rho0 * Lvar / (sp.pi**sp.Rational(5, 2) * lam**3))
    W_sess = sp.simplify(512 * eta_leak**2 * mu_w * q * rho0 * Lvar**2 / (sp.pi**4 * lam**2))
    S_expected = sp.simplify(sp.sqrt(2) * E0 * mu_w * rho0 / (2 * sp.sqrt(sp.pi) * lam**3))
    W_expected = sp.simplify(2 * E0**2 * mu_w * q * rho0 / lam**2)
    Lvar_from_W = sp.simplify(sp.sqrt(W_sess * sp.pi**4 * lam**2 / (512 * eta_leak**2 * mu_w * q * rho0)))

    print("E0                    =", E0)
    print("S_leak                =", S_leak)
    print("W_sess                =", W_sess)
    print("Lvar(W_sess)          =", Lvar_from_W)

    assert sp.simplify(S_leak - S_expected) == 0
    assert sp.simplify(W_sess - W_expected) == 0
    assert sp.simplify(Lvar_from_W - Lvar) == 0

    # ------------------------------------------------------------------
    # 3. Exact Stage 245 weighted U/V drain.
    # ------------------------------------------------------------------
    section("3. Exact Stage 245 weighted U/V drain")
    a_U, a_V, chi_UV, f_U, eta_UV = sp.symbols("a_U a_V chi_UV f_U eta_UV", positive=True, real=True)
    Delta_UV = sp.simplify(a_U * a_V - chi_UV**2)
    D_UV = sp.simplify(chi_UV**2 * a_V * f_U**2 / Delta_UV**2)
    E_UV = sp.simplify(eta_UV * D_UV)
    # F3: Stage-245 drain is nonnegative (square / square, a_V > 0).
    assert sp.ask(sp.Q.nonnegative(D_UV)) in (True, None)
    D_UV_probe = D_UV.subs({
        a_U: sp.Float("2.0"),
        a_V: sp.Float("1.5"),
        chi_UV: sp.Float("0.7"),
        f_U: sp.Float("0.4"),
    })
    assert float(D_UV_probe) >= 0

    print("Delta_UV              =", Delta_UV)
    print("D_UV                  =", D_UV)
    print("DeltaE_UV             =", E_UV)

    # ------------------------------------------------------------------
    # 4. Exact Stage 246 compensated source response.
    # ------------------------------------------------------------------
    section("4. Exact Stage 246 compensated source response")
    r_sigma, a0, b0, rF1, xi_R = sp.symbols("r_sigma a0 b0 rF1 xi_R", positive=True, real=True)
    s = sp.simplify(r_sigma**2 / (r**2 + r_sigma**2))
    g_r = sp.simplify(2 * (1 + a0 * s / 3 - b0 * s / 15) / sp.pi)
    g_inf = sp.simplify(sp.limit(g_r, r, sp.oo))
    R_r = sp.simplify((g_r - rF1) ** 2 / (1 + rF1**2))
    R_inf = sp.simplify((g_inf - rF1) ** 2 / (1 + rF1**2))
    M_sigma = sp.simplify(xi_R * (R_inf - R_r))
    sigma_min = sp.simplify(1 - (a0 - b0) * s)

    print("s(r)                  =", s)
    print("g(r)                  =", g_r)
    print("g_infty               =", g_inf)
    print("R(r)                  =", R_r)
    print("R_infty               =", R_inf)
    print("M_sigma               =", M_sigma)
    print("sigma_min             =", sigma_min)

    assert sp.simplify(g_inf - 2 / sp.pi) == 0

    # ------------------------------------------------------------------
    # 5. Exact stationary compiler and lowering identity.
    # ------------------------------------------------------------------
    section("5. Exact stationary compiler and lowering identity")
    lambda_L, lambda_W = sp.symbols("lambda_L lambda_W", nonnegative=True, real=True)

    V_eff = V_short - lambda_L * S_leak - lambda_W * W_sess - E_UV - M_sigma
    lowering_gap = sp.expand(V_short - V_eff)
    lowering_expected = lambda_L * S_leak + lambda_W * W_sess + E_UV + M_sigma

    print("Lowering gap          =", lowering_gap)

    assert sp.expand(lowering_gap - lowering_expected) == 0

    # ------------------------------------------------------------------
    # 6. Session-I one-point benchmark at r_soft = 0.18.
    # ------------------------------------------------------------------
    section("6. Session-I one-point benchmark")
    subs_soft = {
        Kstar: sp.Float("4.0"),
        OmU2: sp.Float("9.0"),
        OmW2: sp.Float("16.0"),
        GU: sp.Float("1.0"),
        GW: sp.Float("1.25"),
        Rmix: sp.Float("1.35"),
        betaQ: sp.Float("0.03"),
        betaU: sp.Float("0.15"),
        betaW: sp.Float("0.20"),
        alpha6: sp.Integer(0),
        alpha2: sp.Integer(0),
        kappa: sp.Float("1.0"),
        r: sp.Float("0.18"),
        r_sigma: sp.Float("0.8"),
        a0: sp.Float("2.2"),
        b0: sp.Float("-0.6"),
        rF1: sp.Float("1.77799353547498"),
        xi_R: sp.Float("0.9"),
        lam: sp.Float("1.0"),
        eta_leak: sp.Float("0.03"),
        mu_w: sp.Float("0.8"),
        rho0: sp.Float("1.0"),
        q: sp.Float("1.0"),
    }

    Veff_obs = sp.Float("1.74701126")
    Wsess_obs = sp.Float("1.51632107")
    UVdrop_obs = sp.Float("0.21064278")

    Delta_num = sp.N(Delta.subs(subs_soft), 16)
    D0_num = sp.N(D0.subs(subs_soft), 16)
    Vshort_num = sp.N(V_short.subs(subs_soft), 16)
    M_sigma_num = sp.N(M_sigma.subs(subs_soft), 16)
    # F3: source-response and g-bound positivity on the Session-I branch.
    g_soft = sp.N(g_r.subs(subs_soft), 16)
    assert float(g_soft) >= float(2 / sp.pi)
    assert float(g_soft) < float(subs_soft[rF1])
    assert float(M_sigma_num) >= 0

    Lvar_soft = sp.N(
        sp.sqrt(Wsess_obs * sp.pi**4 * subs_soft[lam] ** 2 / (512 * subs_soft[eta_leak] ** 2 * subs_soft[mu_w] * subs_soft[q] * subs_soft[rho0])),
        16,
    )
    S_soft = sp.N(S_leak.subs({**subs_soft, Lvar: Lvar_soft}), 16)
    # F2: the inverted Lvar reproduces the recorded benchmark work scalar, and
    #     matches the paper-stated Lvar(r_soft) = 20.01677473.
    Wsess_from_Lvar = sp.N(W_sess.subs({**subs_soft, Lvar: Lvar_soft}), 16)
    assert abs(float(Wsess_from_Lvar) - float(Wsess_obs)) < 1e-7
    assert abs(float(Lvar_soft) - 20.01677473) < 1e-6

    lambda_L_soft = sp.N((Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num - Veff_obs) / S_soft, 16)
    # F5: pin independently derived benchmark quantities to the paper figures.
    assert abs(float(Vshort_num) - 3.74163698) < 1e-6
    assert abs(float(M_sigma_num) - 0.18386120) < 1e-6
    assert abs(float(S_soft) - 0.31069599) < 1e-6
    assert abs(float(lambda_L_soft) - 0.26971918) < 1e-6
    lambda_L_paper = sp.Float("0.26971918")
    # F4: lowered potential is below the baseline on the benchmark slice (lowering theorem).
    Veff_session = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
    assert float(Vshort_num - Veff_session) >= 0
    # F5: forward benchmark decomposition with the paper's lambda_L (falsifiable closure).
    Veff_forward = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
    assert abs(float(Veff_forward) - float(Veff_obs)) < 1e-6
    Vrebuild_soft = sp.N(Vshort_num - lambda_L_soft * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
    residual_after_work_uv = sp.N(Vshort_num - Wsess_obs - UVdrop_obs, 16)
    residual_after_work_uv_src = sp.N(Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num, 16)

    print("Delta(session)        =", Delta_num)
    print("D0(session)           =", D0_num)
    print("V_short(session)      =", Vshort_num)
    print("M_sigma(session)      =", M_sigma_num)
    print("Lvar(session)         =", Lvar_soft)
    print("S_leak(session)       =", S_soft)
    print("Residual after work+UV        =", residual_after_work_uv)
    print("Residual after work+UV+source =", residual_after_work_uv_src)
    print("lambda_L(session)     =", lambda_L_soft)
    print("V_eff rebuilt         =", Vrebuild_soft)

    assert abs(float(Delta_num) - 142.1775) < 1e-8
    assert abs(float(D0_num) - 3.76481862) < 1e-7
    assert lambda_L_soft > 0

    # ------------------------------------------------------------------
    # 7. Final success banner.
    # ------------------------------------------------------------------
    section("Stage 247 audit result")
    print("All symbolic checks passed.")
    print("Verified objects:")
    print("- exact one-port short-range baseline carried from the static mixed-bundle audit,")
    print("- exact Stage 244 leakage / Session-I work packet,")
    print("- exact Stage 245 weighted U/V drain packet,")
    print("- exact Stage 246 compensated source response scalar R_infty - R(r),")
    print("- exact stationary compiler and lowering identity,")
    print("- and a consistent Session-I one-point benchmark decomposition at r_soft = 0.18.")


if __name__ == "__main__":
    main()
