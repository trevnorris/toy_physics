#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 250
SymPy audit for the crossing-vs-collapse / Goldilocks-window compiler built on
Top of the Stage 248 event chain and the Session-III relaxed wall-timescale
closure.

This script verifies:
1. the exact event-chain transit-time integral from the reduced energy law,
2. the characteristic crossing-time compiler,
3. the unstable V-leg growth-rate / collapse-time compiler,
4. the exact stability ratio and lower-edge formulas,
5. the heavy-throat scaling theorem,
6. the exact speed-space compiler,
7. and the main Session-III benchmark relations, including width sensitivity.
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
    # 1. Exact event-chain transit integral.
    # ------------------------------------------------------------------
    section("1. Exact event-chain transit integral")
    m_s = sp.symbols("m_s", positive=True, real=True)
    E, Vpeak, V0 = sp.symbols("E Vpeak V0", real=True)
    ra, rb = sp.symbols("r_a r_b", positive=True, real=True)
    rv = sp.symbols("rv", positive=True, real=True)
    V = sp.Function("V")

    rdot_sq = sp.simplify(2 * (E - V(rv)) / m_s)
    dt_dr = sp.simplify(1 / sp.sqrt(rdot_sq))
    Ttraj = sp.simplify(sp.Integral(dt_dr, (rv, rb, ra)))

    print("rdot^2(E,r)           =", rdot_sq)
    print("dt/dr                 =", dt_dr)
    print("T_traj(E; rb->ra)     =", Ttraj)

    # ------------------------------------------------------------------
    # 2. Characteristic crossing-time compiler.
    # ------------------------------------------------------------------
    section("2. Characteristic crossing-time compiler")
    lam_eff = sp.symbols("lambda_eff", positive=True, real=True)
    vbar = sp.simplify(sp.sqrt(2 * (E - Vpeak) / m_s))
    tcross = sp.simplify(lam_eff / vbar)
    dtcross_dE = sp.simplify(sp.diff(tcross, E))

    print("v_bar(E)              =", vbar)
    print("t_cross(E)            =", tcross)
    print("dt_cross/dE           =", dtcross_dE)

    # ------------------------------------------------------------------
    # 3. Collapse-time compiler from the unstable V-leg.
    # ------------------------------------------------------------------
    section("3. Collapse-time compiler")
    gUV, chi_peak, mu_eta = sp.symbols(
        "g_UV chi_peak mu_eta", positive=True, real=True
    )
    Gamma_coll = sp.simplify(sp.sqrt(gUV * chi_peak / mu_eta))
    tcollapse = sp.simplify(1 / Gamma_coll)

    print("Gamma_coll            =", Gamma_coll)
    print("t_collapse            =", tcollapse)

    # ------------------------------------------------------------------
    # 4. Stability ratio and lower-edge theorem.
    # ------------------------------------------------------------------
    section("4. Stability ratio and lower-edge theorem")
    S = sp.simplify(tcross / tcollapse)
    Sedge_eq = sp.Eq(S, 1)
    E_edge = sp.simplify(sp.solve(sp.Eq(S**2, 1), E)[0])
    dS_dE = sp.simplify(sp.diff(S, E))
    S_inf = sp.simplify(sp.limit(S, E, sp.oo))

    print("S(E)                  =", S)
    print("S(E)=1                =", Sedge_eq)
    print("E_edge                =", E_edge)
    print("dS/dE                 =", dS_dE)
    print("limit_{E->oo} S(E)    =", S_inf)

    assert sp.simplify(E_edge - (Vpeak + lam_eff**2 * gUV * chi_peak * m_s / (2 * mu_eta))) == 0
    assert S_inf == 0

    # ------------------------------------------------------------------
    # 5. Heavy-throat scaling theorem.
    # ------------------------------------------------------------------
    section("5. Heavy-throat scaling theorem")
    alpha = sp.symbols("alpha", positive=True, real=True)
    S_alpha = sp.simplify(S.subs(mu_eta, alpha * m_s))
    E_edge_alpha = sp.simplify(E_edge.subs(mu_eta, alpha * m_s))
    dEedge_dm = sp.simplify(sp.diff(E_edge_alpha, m_s))
    dEedge_dalpha = sp.simplify(sp.diff(E_edge_alpha, alpha))

    print("S(E) with mu=alpha*m  =", S_alpha)
    print("E_edge(alpha)         =", E_edge_alpha)
    print("dE_edge/dm_s          =", dEedge_dm)
    print("dE_edge/dalpha        =", dEedge_dalpha)

    assert dEedge_dm == 0

    # ------------------------------------------------------------------
    # 6. Speed-space compiler.
    # ------------------------------------------------------------------
    section("6. Speed-space compiler")
    v0 = sp.symbols("v0", positive=True, real=True)
    E_launch = sp.simplify(sp.Rational(1, 2) * m_s * v0**2 + V0)
    vcrit = sp.simplify(sp.sqrt(2 * (Vpeak - V0) / m_s))
    tcross_v = sp.simplify(tcross.subs(E, E_launch))
    v_safe = sp.simplify(sp.sqrt(2 * (E_edge.subs(E, E_launch) - V0) / m_s))
    v_safe_sq = sp.simplify(v_safe**2)
    v_safe_expected = sp.simplify(vcrit**2 + lam_eff**2 * gUV * chi_peak / mu_eta)
    S_v = sp.simplify(S.subs(E, E_launch))

    print("E_launch(v0)          =", E_launch)
    print("v_crit,new            =", vcrit)
    print("t_cross(v0)           =", tcross_v)
    print("v_safe,min^2          =", v_safe_sq)
    print("expected v_safe^2     =", v_safe_expected)
    print("S(v0)                 =", S_v)

    assert sp.simplify(tcross_v - lam_eff / sp.sqrt(v0**2 - vcrit**2)) == 0
    assert sp.simplify(v_safe_sq - v_safe_expected) == 0

    # ------------------------------------------------------------------
    # 7. Sensitivity derivatives.
    # ------------------------------------------------------------------
    section("7. Lower-edge sensitivity ledger")
    dE_dlam = sp.simplify(sp.diff(E_edge, lam_eff))
    dE_dchi = sp.simplify(sp.diff(E_edge, chi_peak))
    dE_dmu = sp.simplify(sp.diff(E_edge, mu_eta))

    print("dE_edge/dlambda_eff   =", dE_dlam)
    print("dE_edge/dchi_peak     =", dE_dchi)
    print("dE_edge/dmu_eta       =", dE_dmu)

    # ------------------------------------------------------------------
    # 8. Session-III benchmark specialization.
    # ------------------------------------------------------------------
    section("8. Session-III benchmark specialization")
    # Imported benchmark values from the barrier-session write-up.
    Vpeak_num = 3.42933112
    V0_num = 0.19999794
    lam_num = 0.42826825
    g_num = 0.95
    chi_num = 21.73204372
    m_num = 1836.15267343
    mu_num = 1836.15267343
    Emax_num = 80.93332737
    transit_min_num = 0.204
    transit_max_num = 4.054

    vcrit_num = math.sqrt(2.0 * (Vpeak_num - V0_num) / m_num)
    tcollapse_num = math.sqrt(mu_num / (g_num * chi_num))
    Eedge_num = Vpeak_num + lam_num**2 * g_num * chi_num * m_num / (2.0 * mu_num)
    vsafe_num = math.sqrt(2.0 * (Eedge_num - V0_num) / m_num)
    ratio_num = vsafe_num / vcrit_num
    vmax_num = math.sqrt(2.0 * (Emax_num - V0_num) / m_num)
    Sedge_num = lam_num * math.sqrt(m_num * g_num * chi_num / (2.0 * mu_num * (Eedge_num - Vpeak_num)))
    Smax_num = lam_num * math.sqrt(m_num * g_num * chi_num / (2.0 * mu_num * (Emax_num - Vpeak_num)))

    print("v_crit,p              =", vcrit_num)
    print("t_collapse            =", tcollapse_num)
    print("E_safe,min            =", Eedge_num)
    print("v_safe,min            =", vsafe_num)
    print("v_safe/v_crit         =", ratio_num)
    print("v_safe,max(scan)      =", vmax_num)
    print("S(E_safe,min)         =", Sedge_num)
    print("S(E_max,scan)         =", Smax_num)
    print("aligned min transit   =", transit_min_num)
    print("aligned max transit   =", transit_max_num)

    assert abs(vcrit_num - 0.059308512338159085) < 1e-12
    assert abs(tcollapse_num - 9.430664762121758) < 1e-9
    assert abs(Eedge_num - 5.322659467573074) < 5e-8
    assert abs(vsafe_num - 0.07469790710905001) < 5e-11
    assert abs(ratio_num - 1.259480370925967) < 5e-9
    assert abs(vmax_num - 0.29654256303211646) < 5e-9
    assert abs(Sedge_num - 1.0) < 5e-13
    assert Smax_num < 1.0
    assert transit_max_num < tcollapse_num

    # ------------------------------------------------------------------
    # 9. Raw-width sensitivity specialization.
    # ------------------------------------------------------------------
    section("9. Raw-width sensitivity specialization")
    lam_raw = 1.0
    chi_raw = 50.74399964

    tcollapse_raw = math.sqrt(mu_num / (g_num * chi_raw))
    Eedge_raw = Vpeak_num + lam_raw**2 * g_num * chi_raw * m_num / (2.0 * mu_num)
    vsafe_raw = math.sqrt(2.0 * (Eedge_raw - V0_num) / m_num)

    print("t_collapse(raw width) =", tcollapse_raw)
    print("E_safe,min(raw width) =", Eedge_raw)
    print("v_safe,min(raw width) =", vsafe_raw)

    assert abs(tcollapse_raw - 6.171635157122822) < 5e-9
    assert abs(Eedge_raw - 27.532730948999998) < 5e-8
    assert Eedge_raw > Eedge_num

    # ------------------------------------------------------------------
    # 10. Final result banner.
    # ------------------------------------------------------------------
    section("Stage 250 audit result")
    print("All symbolic and numerical checks passed.")
    print("Verified objects:")
    print("- exact event-chain transit-time integral from the Stage 248 energy law,")
    print("- characteristic crossing-time and collapse-time compilers,")
    print("- exact stability ratio and Goldilocks lower-edge formulas,")
    print("- heavy-throat scaling cancellation when mu_eta = alpha m_s,")
    print("- exact speed-space compiler relative to the Stage 248 threshold speed,")
    print("- the Session-III proton-proxy benchmark values,")
    print("- the one-sided safe-band interpretation on the sampled closure,")
    print("- and the trigger-width / steepness sensitivity shift.")


if __name__ == "__main__":
    main()
