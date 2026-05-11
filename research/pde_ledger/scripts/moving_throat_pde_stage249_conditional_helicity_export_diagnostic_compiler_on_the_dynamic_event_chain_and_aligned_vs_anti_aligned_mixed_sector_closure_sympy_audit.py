#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 249
SymPy audit for the conditional helicity-export diagnostic attached to the
Stage 248 dynamic event chain.

This script verifies:
1. the exact subscale-helicity transfer equation obtained by subtracting the
   projected and resolved helicity ledgers,
2. the aligned/anti-aligned two-state closure and its reduction to one
   asymmetry scalar alpha_h,
3. the exact peak-ratio and integrated-ratio Möbius compilers,
4. the independence of those ratios from the overall export scale eta_h,
5. and the main Session-II benchmark numbers reported in the barrier write-up.
"""

from __future__ import annotations

import math

import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Exact projected-minus-resolved helicity transfer equation.
    # ------------------------------------------------------------------
    section("1. Exact projected-minus-resolved helicity compiler")
    t = sp.symbols("t", real=True)
    hbar = sp.Function("hbar")
    hres = sp.Function("hres")
    Fbar = sp.Function("Fbar")
    Fres = sp.Function("Fres")
    Cfull = sp.Function("Cfull")
    Cres = sp.Function("Cres")

    eq_proj = sp.diff(hbar(t), t) + Fbar(t) + 2 * Cfull(t)
    eq_res = sp.diff(hres(t), t) + Fres(t) + 2 * Cres(t)

    hsub = sp.simplify(hbar(t) - hres(t))
    Fsub = sp.simplify(Fbar(t) - Fres(t))
    Csub = sp.simplify(Cfull(t) - Cres(t))

    eq_sub = sp.simplify(eq_proj - eq_res)
    eq_sub_expected = sp.simplify(sp.diff(hsub, t) + Fsub + 2 * Csub)

    print("h_sub                 =", hsub)
    print("F_h,sub               =", Fsub)
    print("C_sub                 =", Csub)
    print("projected-resolved eq =", eq_sub)
    print("expected sub eq       =", eq_sub_expected)

    assert sp.simplify(eq_sub - eq_sub_expected) == 0

    # ------------------------------------------------------------------
    # 2. Volume-integrated aligned/anti-aligned closure.
    # ------------------------------------------------------------------
    section("2. Exact aligned/anti-aligned closure")
    sigma = sp.symbols("sigma", real=True)
    Phi0, Phi1 = sp.symbols("Phi0 Phi1", real=True)
    C0, C1 = sp.symbols("C0 C1", real=True)
    Gamma0, Gamma1 = sp.symbols("Gamma0 Gamma1", real=True)
    alpha_h = sp.symbols("alpha_h", real=True)
    eta_h, G0 = sp.symbols("eta_h G0", positive=True, real=True)

    Phi_sigma = Phi0 + sigma * Phi1
    C_sigma = C0 + sigma * C1
    Hdot_sigma = sp.expand(-Phi_sigma - 2 * C_sigma)
    Hdot_expected = sp.expand((-(Phi0 + 2 * C0)) + sigma * (-(Phi1 + 2 * C1)))

    Gamma0_def = sp.simplify(-(Phi0 + 2 * C0))
    Gamma1_def = sp.simplify(-(Phi1 + 2 * C1))
    Hdot_factored = sp.simplify((Gamma0 + sigma * Gamma1).subs(Gamma1, alpha_h * Gamma0))
    Hdot_scale = sp.simplify(Hdot_factored.subs(Gamma0, eta_h * G0))

    print("Phi_sigma             =", Phi_sigma)
    print("C_sigma               =", C_sigma)
    print("Hdot_sigma            =", Hdot_sigma)
    print("Hdot expected         =", Hdot_expected)
    print("Gamma0                =", Gamma0_def)
    print("Gamma1                =", Gamma1_def)
    print("Gamma0(1+sigma alpha) =", Hdot_factored)
    print("eta_h-scaled form     =", Hdot_scale)

    assert sp.simplify(Hdot_sigma - Hdot_expected) == 0
    assert sp.simplify(Hdot_factored - Gamma0 * (1 + sigma * alpha_h)) == 0
    assert sp.simplify(Hdot_scale - eta_h * G0 * (1 + sigma * alpha_h)) == 0

    # ------------------------------------------------------------------
    # 3. Peak-ratio Möbius compiler and inverse.
    # ------------------------------------------------------------------
    section("3. Exact peak-ratio compiler")
    Rpk = sp.symbols("Rpk", positive=True, real=True)
    Gplus = sp.simplify(Gamma0 * (1 + alpha_h))
    Gminus = sp.simplify(Gamma0 * (1 - alpha_h))
    Rpk_formula = sp.simplify(Gplus / Gminus)
    alpha_from_Rpk = sp.simplify(sp.solve(sp.Eq(Rpk, Rpk_formula), alpha_h)[0])

    print("G_plus                =", Gplus)
    print("G_minus               =", Gminus)
    print("R_pk                  =", Rpk_formula)
    print("alpha(R_pk)           =", alpha_from_Rpk)

    assert sp.simplify(Rpk_formula - (1 + alpha_h) / (1 - alpha_h)) == 0
    assert sp.simplify(alpha_from_Rpk - (Rpk - 1) / (Rpk + 1)) == 0

    # ------------------------------------------------------------------
    # 4. Integrated-ratio Möbius compiler and eta_h cancellation.
    # ------------------------------------------------------------------
    section("4. Exact integrated-ratio compiler")
    I0, I1, Rint = sp.symbols("I0 I1 Rint", positive=True, real=True)
    abar = sp.symbols("abar", real=True)

    Hplus_int = sp.simplify(eta_h * (I0 + I1))
    Hminus_int = sp.simplify(eta_h * (I0 - I1))
    Rint_formula = sp.simplify(Hplus_int / Hminus_int)
    Rint_expected = sp.simplify((1 + I1 / I0) / (1 - I1 / I0))
    abar_from_Rint = sp.simplify(sp.solve(sp.Eq(Rint, (1 + abar) / (1 - abar)), abar)[0])

    print("H_plus(int)           =", Hplus_int)
    print("H_minus(int)          =", Hminus_int)
    print("R_int                 =", Rint_formula)
    print("R_int expected        =", Rint_expected)
    print("abar(R_int)           =", abar_from_Rint)

    assert sp.simplify(Rint_formula - Rint_expected) == 0
    assert eta_h not in sp.simplify(Rint_formula).free_symbols
    assert sp.simplify(abar_from_Rint - (Rint - 1) / (Rint + 1)) == 0

    # ------------------------------------------------------------------
    # 5. Session-II benchmark specialization.
    # ------------------------------------------------------------------
    section("5. Session-II benchmark specialization")
    peak_aligned = 281.79830789
    peak_antialigned = 56.96878122
    hint_aligned = 20.58070146
    hint_antialigned = 5.00843357
    ratio_integrated_report = 4.10920923
    Xi_turn = 0.34437471
    lambda_th = 0.42826825
    v_cross = 2.59221845

    ratio_peak = peak_aligned / peak_antialigned
    alpha_peak_num = (ratio_peak - 1.0) / (ratio_peak + 1.0)
    ratio_final = hint_aligned / hint_antialigned
    alpha_int_num = (ratio_integrated_report - 1.0) / (ratio_integrated_report + 1.0)
    alpha_final_num = (ratio_final - 1.0) / (ratio_final + 1.0)

    peak_difference = peak_aligned - peak_antialigned
    final_difference = hint_aligned - hint_antialigned
    peak_sum = peak_aligned + peak_antialigned
    final_sum = hint_aligned + hint_antialigned

    print("peak aligned          =", peak_aligned)
    print("peak anti-aligned     =", peak_antialigned)
    print("R_pk(session)         =", ratio_peak)
    print("alpha_pk(session)     =", alpha_peak_num)
    print("H_sub,final aligned   =", hint_aligned)
    print("H_sub,final anti      =", hint_antialigned)
    print("R_final(session)      =", ratio_final)
    print("R_int(report)         =", ratio_integrated_report)
    print("alpha_int(report)     =", alpha_int_num)
    print("alpha_final(session)  =", alpha_final_num)
    print("peak difference       =", peak_difference)
    print("peak sum              =", peak_sum)
    print("final difference      =", final_difference)
    print("final sum             =", final_sum)
    print("Xi_turn(session)      =", Xi_turn)
    print("lambda_th(session)    =", lambda_th)
    print("v_cross(session)      =", v_cross)

    assert peak_aligned > 0.0 and peak_antialigned > 0.0
    assert hint_aligned > 0.0 and hint_antialigned > 0.0
    assert ratio_peak > 1.0
    assert ratio_final > 1.0
    assert 0.0 < alpha_peak_num < 1.0
    assert 0.0 < alpha_int_num < 1.0
    assert abs(ratio_peak - 4.94653917) < 5e-8
    assert abs(alpha_peak_num - 0.663669919237628) < 5e-13
    assert abs(ratio_final - ratio_integrated_report) < 5e-9
    assert abs(alpha_int_num - 0.6085499908172678) < 5e-13
    assert abs(alpha_final_num - alpha_int_num) < 1e-10
    assert alpha_peak_num > alpha_int_num
    assert Xi_turn > 0.0
    assert lambda_th > 0.0
    assert v_cross > 0.0

    # ------------------------------------------------------------------
    # 6. Final result banner.
    # ------------------------------------------------------------------
    section("Stage 249 audit result")
    print("All symbolic and numerical checks passed.")
    print("Verified objects:")
    print("- exact projected-minus-resolved subscale-helicity transfer law,")
    print("- exact linear aligned/anti-aligned export closure and its reduction to one asymmetry scalar,")
    print("- exact peak-ratio and integrated-ratio Möbius compilers,")
    print("- exact cancellation of the overall export scale eta_h from the preference ratios,")
    print("- and the main Session-II benchmark packet: R_pk, R_int, alpha_pk, and average asymmetry.")


if __name__ == "__main__":
    main()
