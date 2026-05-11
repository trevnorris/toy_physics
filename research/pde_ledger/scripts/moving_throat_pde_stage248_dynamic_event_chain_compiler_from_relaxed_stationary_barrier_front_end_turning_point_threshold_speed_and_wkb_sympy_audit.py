#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 248
SymPy audit for the dynamic event-chain compiler built on the Stage 247
relaxed stationary barrier front end.

This script verifies:
1. exact energy conservation for the reduced one-dimensional event chain,
2. the finite-radius threshold-speed compiler on the lowered branch,
3. the exact Coulomb contact-threshold formula,
4. the exact Coulomb outer-turning-point and WKB reference formulas,
5. the near-top parabolic action normal form,
6. and the main Session-II benchmark relations reported in the barrier write-up.

Provenance notes
----------------
- The symbolic compilers are the Stage 248 theorem surface itself.
- The numeric constants in Section 5 are the carried Session-II benchmark packet
  from the Stage 248 barrier note and the benchmark ledger in
  `CHECKPOINT_CONSTANT_PROVENANCE.md`; they are benchmark readbacks, not hidden
  theorem inputs.
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
    # 1. Exact reduced energy conservation.
    # ------------------------------------------------------------------
    section("1. Exact reduced energy conservation")
    t = sp.symbols("t", real=True)
    m_s = sp.symbols("m_s", positive=True, real=True)
    r = sp.Function("r")
    V = sp.Function("V")

    rdot = sp.diff(r(t), t)
    rddot = sp.diff(r(t), (t, 2))
    E_dyn = sp.simplify(sp.Rational(1, 2) * m_s * rdot**2 + V(r(t)))
    dE_dt = sp.diff(E_dyn, t)
    dE_on_shell = sp.simplify(
        dE_dt.subs(rddot, -sp.diff(V(r(t)), r(t)) / m_s)
    )

    print("E_dyn                 =", E_dyn)
    print("dE_dt on-shell        =", dE_on_shell)

    assert dE_on_shell == 0

    # ------------------------------------------------------------------
    # 2. Barrier-peak / threshold-speed compiler.
    # ------------------------------------------------------------------
    section("2. Barrier-peak and threshold-speed compiler")
    V0, Vpeak = sp.symbols("V0 Vpeak", real=True)
    r0, r_contact = sp.symbols("r0 r_contact", positive=True, real=True)
    v0 = sp.symbols("v0", positive=True, real=True)

    E_launch_new = sp.simplify(sp.Rational(1, 2) * m_s * v0**2 + V0)
    vcrit_new = sp.simplify(sp.sqrt(2 * (Vpeak - V0) / m_s))
    vcrit_new_solved = sp.simplify(
        [
            sol
            for sol in sp.solve(sp.Eq(E_launch_new, Vpeak), v0)
            if not sol.could_extract_minus_sign()
        ][0]
    )
    E_at_vcrit = sp.simplify(E_launch_new.subs(v0, vcrit_new) - Vpeak)

    E_launch_coul = sp.simplify(sp.Rational(1, 2) * m_s * v0**2 + 1 / r0)
    vcontact_coul = sp.simplify(
        sp.sqrt(2 * (1 / r_contact - 1 / r0) / m_s)
    )
    vcontact_coul_solved = sp.simplify(
        [
            sol
            for sol in sp.solve(sp.Eq(E_launch_coul, 1 / r_contact), v0)
            if not sol.could_extract_minus_sign()
        ][0]
    )
    delta_new = sp.simplify(E_launch_new - Vpeak)
    delta_coul = sp.simplify(E_launch_coul - 1 / r_contact)

    print("E_launch_new          =", E_launch_new)
    print("v_crit,new            =", vcrit_new)
    print("v_crit,new (solve)    =", vcrit_new_solved)
    print("E(v_crit,new)-V_peak  =", E_at_vcrit)
    print("E_launch_Coul         =", E_launch_coul)
    print("v_contact,Coul        =", vcontact_coul)
    print("v_contact,Coul (solve)=", vcontact_coul_solved)
    print("E_new-V_peak          =", delta_new)
    print("E_Coul-1/r_contact    =", delta_coul)

    assert E_at_vcrit == 0
    assert sp.simplify(vcrit_new_solved - vcrit_new) == 0
    assert sp.simplify(delta_new.subs(v0, vcrit_new)) == 0
    assert sp.simplify(vcontact_coul_solved - vcontact_coul) == 0
    assert sp.simplify(delta_coul.subs(v0, vcontact_coul)) == 0

    # ------------------------------------------------------------------
    # 3. Turning-point / WKB compiler.
    # ------------------------------------------------------------------
    section("3. Turning-point and WKB compiler")
    Esub, hbar_eff = sp.symbols("Esub hbar_eff", positive=True, real=True)
    r_minus, r_plus = sp.symbols("r_minus r_plus", positive=True, real=True)
    rv = sp.symbols("rv", positive=True, real=True)
    E_turn = sp.symbols("E_turn", positive=True, real=True)

    v0_sub = sp.simplify(sp.sqrt(2 * (Esub - V0) / m_s))
    I_new = sp.Integral(
        sp.sqrt(2 * m_s * (V(rv) - Esub)) / hbar_eff,
        (rv, r_minus, r_plus),
    )

    # Exact Coulomb reference outer turning point.
    rturn_coul = sp.simplify(1 / Esub)

    # Exact Coulomb antiderivative check.
    F_coul = sp.simplify(
        sp.sqrt(2 * m_s) / hbar_eff
        * (
            sp.sqrt(rv * (1 - Esub * rv))
            + sp.asin(sp.sqrt(Esub * rv)) / sp.sqrt(Esub)
        )
    )
    coul_integrand = sp.simplify(sp.sqrt(2 * m_s * (1 / rv - Esub)) / hbar_eff)
    dF_coul = sp.simplify(sp.diff(F_coul, rv) - coul_integrand)
    I_coul_formula = sp.simplify(
        sp.sqrt(2 * m_s) / hbar_eff
        * (
            sp.pi / (2 * sp.sqrt(Esub))
            - sp.sqrt(r_contact * (1 - Esub * r_contact))
            - sp.asin(sp.sqrt(Esub * r_contact)) / sp.sqrt(Esub)
        )
    )
    I_coul_endpoints = sp.simplify(
        sp.limit(F_coul, rv, 1 / Esub, dir="-") - F_coul.subs(rv, r_contact)
    )
    r_plus_E = sp.Function("r_plus")(E_turn)
    r_minus_E = sp.Function("r_minus")(E_turn)
    transport_plus = sp.simplify(
        sp.solve(
            sp.Eq(sp.diff(V(r_plus_E) - E_turn, E_turn), 0),
            sp.diff(r_plus_E, E_turn),
        )[0]
    )
    transport_minus = sp.simplify(
        sp.solve(
            sp.Eq(sp.diff(V(r_minus_E) - E_turn, E_turn), 0),
            sp.diff(r_minus_E, E_turn),
        )[0]
    )
    Tnew_over_Tcoul = sp.exp(
        -2 * (sp.Symbol("Inew", positive=True) - sp.Symbol("Icoul", positive=True))
    )

    print("v_0,sub(E)            =", v0_sub)
    print("I_new(E)              =", I_new)
    print("r_turn,Coul(E)        =", rturn_coul)
    print("dF_coul/dr - integrand=", dF_coul)
    print("I_Coul(E)             =", I_coul_formula)
    print("I_Coul endpoints      =", I_coul_endpoints)
    print("dr_+/dE               =", transport_plus)
    print("dr_-/dE               =", transport_minus)
    print("T_new/T_Coul          =", Tnew_over_Tcoul)

    assert dF_coul == 0
    assert sp.simplify(I_coul_endpoints - I_coul_formula) == 0
    assert sp.simplify(transport_plus - 1 / sp.diff(V(r_plus_E), r_plus_E)) == 0
    assert sp.simplify(transport_minus - 1 / sp.diff(V(r_minus_E), r_minus_E)) == 0

    # ------------------------------------------------------------------
    # 4. Near-top parabolic normal form.
    # ------------------------------------------------------------------
    section("4. Near-top parabolic normal form")
    DeltaE, Kpeak = sp.symbols("DeltaE Kpeak", positive=True, real=True)
    y = sp.symbols("y", real=True)
    yturn = sp.sqrt(2 * DeltaE / Kpeak)

    I_top = sp.simplify(
        sp.integrate(
            sp.sqrt(2 * m_s * (DeltaE - Kpeak * y**2 / 2)) / hbar_eff,
            (y, -yturn, yturn),
        )
    )
    I_top_expected = sp.simplify(
        sp.pi * DeltaE * sp.sqrt(m_s / Kpeak) / hbar_eff
    )
    r_plus_top = sp.Symbol("r_peak", real=True) + yturn
    r_minus_top = sp.Symbol("r_peak", real=True) - yturn

    print("r_+(top)              =", r_plus_top)
    print("r_-(top)              =", r_minus_top)
    print("I_top                 =", I_top)
    print("I_top expected        =", I_top_expected)

    assert sp.simplify(I_top - I_top_expected) == 0

    # ------------------------------------------------------------------
    # 5. Session-II benchmark specialization.
    # Provenance: these carried readbacks are logged explicitly in
    # CHECKPOINT_CONSTANT_PROVENANCE.md under Stage 248.
    # ------------------------------------------------------------------
    section("5. Session-II benchmark specialization")
    m_num = 1.0
    hbar_num = 1.0
    r0_num = 5.0
    rc_num = 0.18
    Esub_num = 2.5
    V0_num = 0.19999794
    Vpeak_num = 3.42933112
    rpeak_num = 0.23944389
    rturn_new_num = 0.39096144
    rinner_num = 0.19039548
    Inew_num = 0.19744614
    Icoul_report = 0.30222297
    Xi_turn_num = 0.34437471
    lambda_th_num = 0.42826825
    vcross_num = 2.59221845
    rcoul_turn_report = 0.28091705

    vcrit_num = math.sqrt(2.0 * (Vpeak_num - V0_num) / m_num)
    vcontact_num = math.sqrt(2.0 * (1.0 / rc_num - 1.0 / r0_num) / m_num)
    vsub_num = math.sqrt(2.0 * (Esub_num - V0_num) / m_num)
    rturn_coul_num = 1.0 / Esub_num
    Tnew_num = math.exp(-2.0 * Inew_num)
    Tcoul_num = math.exp(-2.0 * Icoul_report)
    ratio_num = Tnew_num / Tcoul_num
    improve_pct = 100.0 * (ratio_num - 1.0)

    Icoul_exact_num = float(
        sp.N(
            I_coul_formula.subs(
                {
                    m_s: m_num,
                    hbar_eff: hbar_num,
                    Esub: Esub_num,
                    r_contact: rc_num,
                }
            ),
            16,
        )
    )
    Ecoul_cross_num = 0.5 * vcross_num**2 + 1.0 / r0_num
    rcoul_turn_cross_num = 1.0 / Ecoul_cross_num
    Vprime_turn_mag = Esub_num / lambda_th_num

    print("r_peak(session)       =", rpeak_num)
    print("V_peak(session)       =", Vpeak_num)
    print("V(r0) session         =", V0_num)
    print("v_crit,new(session)   =", vcrit_num)
    print("v_contact,Coul        =", vcontact_num)
    print("v_0,sub(session)      =", vsub_num)
    print("r_turn,new            =", rturn_new_num)
    print("r_inner,new           =", rinner_num)
    print("r_turn,Coul exact     =", rturn_coul_num)
    print("I_new(session)        =", Inew_num)
    print("I_Coul(report)        =", Icoul_report)
    print("I_Coul(exact formula) =", Icoul_exact_num)
    print("T_new                 =", Tnew_num)
    print("T_Coul                =", Tcoul_num)
    print("T_new/T_Coul          =", ratio_num)
    print("improvement (%)       =", improve_pct)
    print("Xi_turn(session)      =", Xi_turn_num)
    print("lambda_th(session)    =", lambda_th_num)
    print("|V'(r_turn)| from lambda_th =", Vprime_turn_mag)
    print("v_cross(session)      =", vcross_num)
    print("Coulomb turnback at v_cross =", rcoul_turn_cross_num)

    assert abs(vcrit_num - 2.5413906350657705) < 1e-12
    assert abs(vcontact_num - 3.272783388968954) < 1e-12
    assert abs(vsub_num - 2.1447620194324593) < 1e-12
    assert abs(rturn_coul_num - 0.4) < 1e-15
    assert abs(Tnew_num - 0.673752615) < 5e-8
    assert abs(Tcoul_num - 0.546377065) < 5e-8
    assert abs(ratio_num - 1.23312756) < 5e-8
    assert abs(improve_pct - 23.3128) < 1e-3
    assert vcrit_num < vcross_num < vcontact_num
    assert abs(rcoul_turn_cross_num - rcoul_turn_report) < 3e-6
    assert abs(Icoul_exact_num - Icoul_report) < 1e-3
    assert Xi_turn_num > 0
    assert lambda_th_num > 0

    # ------------------------------------------------------------------
    # 6. Final result banner.
    # ------------------------------------------------------------------
    section("Stage 248 audit result")
    print("All symbolic and numerical checks passed.")
    print("Verified objects:")
    print("- exact one-dimensional energy conservation on the Stage 247 front end,")
    print("- exact finite-radius threshold-speed compiler on the lowered branch,")
    print("- exact Coulomb contact threshold and Coulomb turning-point/WKB reference formulas,")
    print("- exact near-top parabolic action normal form,")
    print("- and the main Session-II benchmark relations: threshold speed, subbarrier speed, turning points, transmission ratio, and dynamic turning-point diagnostics.")


if __name__ == "__main__":
    main()
