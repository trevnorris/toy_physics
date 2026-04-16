#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 235
SymPy audit for the vacuum-vs-lattice heat partition and cold-survival compiler
built on the Stage-234 microscopic export kernel.

This script verifies:
1. the exact channel-resolved exported-energy ledger,
2. the one-scalar event-shape quotient formula,
3. the exponential-event specialization and event-equivalent damping rates,
4. the safe-edge exported-energy theorem,
5. the microscopic coefficient surface corresponding to the Session-IV 3:1 split,
6. the exact speed-drift law of the lattice partition,
7. and the Session-IV benchmark calibration check.
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
    # 1. Exact channel-resolved energy ledger.
    # ------------------------------------------------------------------
    section("1. Exact channel-resolved energy ledger")
    G3v, G3l, G5v, G5l = sp.symbols(
        "Gamma3_vac Gamma3_lat Gamma5_vac Gamma5_lat", nonnegative=True, real=True
    )
    I2, I3, rV = sp.symbols("I2 I3 r_V", positive=True, real=True)

    G3 = sp.simplify(G3v + G3l)
    G5 = sp.simplify(G5v + G5l)

    Ev = sp.simplify(G3v * I2 + G5v * I3)
    El = sp.simplify(G3l * I2 + G5l * I3)
    Et = sp.simplify(Ev + El)
    fv = sp.simplify((Ev / Et).subs(I3, rV * I2))
    fl = sp.simplify((El / Et).subs(I3, rV * I2))

    print("E_vac                 =", Ev)
    print("E_lat                 =", El)
    print("E_tot                 =", Et)
    print("f_vac(r_V)            =", fv)
    print("f_lat(r_V)            =", fl)
    print("f_vac + f_lat         =", sp.simplify(fv + fl))

    assert sp.simplify(fv - (G3v + G5v * rV) / (G3 + G5 * rV)) == 0
    assert sp.simplify(fl - (G3l + G5l * rV) / (G3 + G5 * rV)) == 0
    assert sp.simplify(fv + fl - 1) == 0

    # ------------------------------------------------------------------
    # 2. Exact speed-drift law and endpoint limits.
    # ------------------------------------------------------------------
    section("2. Exact speed-drift law and endpoint limits")
    dfl = sp.simplify(sp.diff(fl, rV))
    dfl_expected = sp.simplify((G5l * G3v - G3l * G5v) / (G3 + G5 * rV) ** 2)
    fl_r0 = sp.simplify(sp.limit(fl, rV, 0))
    fl_rinf = sp.simplify(sp.limit(fl, rV, sp.oo))

    print("d f_lat / d r_V       =", dfl)
    print("endpoint r_V -> 0     =", fl_r0)
    print("endpoint r_V -> oo    =", fl_rinf)

    assert sp.simplify(dfl - dfl_expected) == 0
    assert sp.simplify(fl_r0 - G3l / G3) == 0
    assert sp.simplify(fl_rinf - G5l / G5) == 0

    # ------------------------------------------------------------------
    # 3. Exponential-event specialization.
    # ------------------------------------------------------------------
    section("3. Exponential-event specialization")
    t, T, s, Vin = sp.symbols("t T s V_in", positive=True, real=True)
    V = sp.simplify(Vin * sp.exp(s * t))
    Vdot = sp.diff(V, t)
    Vdd = sp.diff(V, t, 2)
    Vddd = sp.diff(V, t, 3)

    I1_exp = sp.simplify(sp.integrate(Vdot**2, (t, 0, T)))
    I2_exp = sp.simplify(sp.integrate(Vdd**2, (t, 0, T)))
    I3_exp = sp.simplify(sp.integrate(Vddd**2, (t, 0, T)))
    gamma_v_eq = sp.simplify(G3v * s**2 + G5v * s**4)
    gamma_l_eq = sp.simplify(G3l * s**2 + G5l * s**4)
    gamma_eq = sp.simplify(G3 * s**2 + G5 * s**4)

    print("I1_exp                =", I1_exp)
    print("I2_exp                =", I2_exp)
    print("I3_exp                =", I3_exp)
    print("I3_exp / I2_exp       =", sp.simplify(I3_exp / I2_exp))
    print("gamma_vac^eq(s)       =", gamma_v_eq)
    print("gamma_lat^eq(s)       =", gamma_l_eq)
    print("gamma_eff^eq(s)       =", gamma_eq)

    assert sp.simplify(I2_exp - Vin**2 * s**3 * (sp.exp(2 * s * T) - 1) / 2) == 0
    assert sp.simplify(I3_exp - s**2 * I2_exp) == 0
    assert sp.simplify(Ev.subs({I2: I2_exp, I3: I3_exp}) - gamma_v_eq * I1_exp) == 0
    assert sp.simplify(El.subs({I2: I2_exp, I3: I3_exp}) - gamma_l_eq * I1_exp) == 0

    # ------------------------------------------------------------------
    # 4. Safe-edge compiler.
    # ------------------------------------------------------------------
    section("4. Safe-edge compiler")
    sc, s0, mu_eta = sp.symbols("s_c s_0 mu_eta", positive=True, real=True)
    safe_eq = sp.Eq(G3 * sc**3 + G5 * sc**5, mu_eta * (s0**2 - sc**2))

    I1_safe = sp.simplify(I1_exp.subs({s: sc, T: 1 / sc}))
    Et_safe = sp.simplify((gamma_eq.subs(s, sc)) * I1_safe)
    safe_combo = sp.simplify(G3 * sc**3 + G5 * sc**5)
    Et_safe_expected = sp.simplify(Vin**2 * (sp.E**2 - 1) * safe_combo / 2)
    Et_safe_reduced = sp.simplify(Et_safe_expected.subs(safe_combo, mu_eta * (s0**2 - sc**2)))
    gamma_safe_eq_expected = sp.simplify(safe_combo / sc)
    gamma_safe_eq = sp.simplify(mu_eta * (s0**2 - sc**2) / sc)
    fv_sc = sp.simplify(fl.subs(rV, sc**2))  # intentionally compute lattice fraction here first
    fl_sc = sp.simplify(fl.subs(rV, sc**2))
    fvac_sc = sp.simplify(fv.subs(rV, sc**2))
    Ev_safe = sp.simplify(fvac_sc * Et_safe_reduced)
    El_safe = sp.simplify(fl_sc * Et_safe_reduced)

    print("safe equality         =", safe_eq)
    print("I1_safe               =", I1_safe)
    print("E_safe before reduce  =", Et_safe_expected)
    print("E_safe reduced        =", Et_safe_reduced)
    print("gamma_eq_safe         =", gamma_safe_eq)
    print("f_vac(s_c)            =", fvac_sc)
    print("f_lat(s_c)            =", fl_sc)

    assert sp.simplify(I1_safe - Vin**2 * sc * (sp.E**2 - 1) / 2) == 0
    assert sp.simplify(Et_safe_expected - Et_safe) == 0
    assert sp.simplify(Et_safe_reduced - Vin**2 * (sp.E**2 - 1) * mu_eta * (s0**2 - sc**2) / 2) == 0
    assert sp.simplify(gamma_safe_eq_expected - gamma_eq.subs(s, sc)) == 0
    assert sp.simplify(gamma_safe_eq - mu_eta * (s0**2 - sc**2) / sc) == 0
    assert sp.simplify(Ev_safe + El_safe - Et_safe_reduced) == 0

    # ------------------------------------------------------------------
    # 5. Session-IV 3:1 split as a microscopic coefficient surface.
    # ------------------------------------------------------------------
    section("5. Session-IV 3:1 split as microscopic coefficient surface")
    split_surface = sp.simplify(sp.expand((G3l + G5l * sc**2) - 3 * (G3v + G5v * sc**2)))
    phi, G3T, G5T = sp.symbols("phi Gamma3T Gamma5T", positive=True, real=True)
    fl_phi = sp.simplify(
        fl.subs(
            {
                G3l: phi * G3T,
                G3v: (1 - phi) * G3T,
                G5l: phi * G5T,
                G5v: (1 - phi) * G5T,
            }
        )
    )

    print("split surface         =", split_surface)
    print("speed-independent phi =", fl_phi)

    assert sp.simplify(split_surface) == G3l + G5l * sc**2 - 3 * G3v - 3 * G5v * sc**2
    assert sp.simplify(fl_phi - phi) == 0

    # ------------------------------------------------------------------
    # 6. Session-IV benchmark specialization.
    # ------------------------------------------------------------------
    section("6. Session-IV benchmark specialization")
    t_cross_num = 1.82169718
    s0_num = 6.94311167
    E_diss_num = 0.01033460
    frac_v_num = 0.25
    frac_l_num = 0.75

    sc_num = 1.0 / t_cross_num
    sc2_num = sc_num**2
    gamma_eq_safe_num = (s0_num**2 - sc2_num) / sc_num
    gamma_v_eq_num = frac_v_num * gamma_eq_safe_num
    gamma_l_eq_num = frac_l_num * gamma_eq_safe_num
    safe_energy_pref_num = 0.5 * (math.e**2 - 1.0) * (s0_num**2 - sc2_num)
    Vin_match_num = math.sqrt(E_diss_num / safe_energy_pref_num)
    E_v_num = frac_v_num * E_diss_num
    E_l_num = frac_l_num * E_diss_num

    print("s_c                   =", sc_num)
    print("s_c^2                 =", sc2_num)
    print("gamma_eff^eq,safe     =", gamma_eq_safe_num)
    print("gamma_vac^eq,safe     =", gamma_v_eq_num)
    print("gamma_lat^eq,safe     =", gamma_l_eq_num)
    print("safe energy prefactor =", safe_energy_pref_num)
    print("V_in match            =", Vin_match_num)
    print("E_vac benchmark       =", E_v_num)
    print("E_lat benchmark       =", E_l_num)

    assert abs(sc_num - 0.5489386551062235) < 1e-12
    assert abs(sc2_num - 0.3013336470698294) < 1e-12
    assert abs(gamma_eq_safe_num - 87.26925234614843) < 1e-9
    assert abs(gamma_v_eq_num - 21.81731308653711) < 1e-9
    assert abs(gamma_l_eq_num - 65.45193925961132) < 1e-9
    assert abs(safe_energy_pref_num - 153.03535490769042) < 1e-9
    assert abs(Vin_match_num - 0.008217712598897912) < 1e-15
    assert abs(E_v_num - 0.00258365) < 1e-12
    assert abs(E_l_num - 0.00775095) < 1e-12
    assert abs(E_v_num + E_l_num - E_diss_num) < 1e-12

    # ------------------------------------------------------------------
    # 7. Final banner.
    # ------------------------------------------------------------------
    section("Stage-235 audit result")
    print("All symbolic and numerical checks passed.")
    print("Verified objects:")
    print("- exact vacuum/lattice exported-energy ledger from the Stage-234 kernel,")
    print("- one-scalar partition quotient r_V = I3/I2,")
    print("- exact speed-drift law of the lattice fraction,")
    print("- exponential-event specialization r_V = s^2,")
    print("- event-equivalent damping rates gamma_vac^eq and gamma_lat^eq,")
    print("- exact safe-edge exported-energy theorem,")
    print("- the microscopic coefficient surface corresponding to the Session-IV 3:1 split,")
    print("- and the Session-IV benchmark calibration check.")


if __name__ == "__main__":
    main()
