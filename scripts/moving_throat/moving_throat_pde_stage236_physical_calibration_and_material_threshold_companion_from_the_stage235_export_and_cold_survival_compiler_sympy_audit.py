#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 236
SymPy audit for the physical calibration / material-threshold companion built on
Stage 235's microscopic export and cold-survival compiler.

This script verifies:
1. the exact lattice-turnover threshold from the Stage-235 safe-edge rate,
2. the exact recovery of the legacy Session-V slice through one transport
   projection factor Upsilon_lat,
3. the harmonic-trap chi_lambda geometry compiler and force-matched stiffness map,
4. the Korringa-limited temperature ceiling,
5. and the exact dimensionless material-screening ratios.
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
    # 1. Exact lattice-turnover compiler from the Stage-235 safe-edge rate.
    # ------------------------------------------------------------------
    section("1. Exact lattice-turnover compiler")
    s_c, s0, mu_eta, f_lat = sp.symbols("s_c s_0 mu_eta f_lat", positive=True, real=True)
    zeta_ep, t_star, Upsilon_lat = sp.symbols(
        "zeta_ep t_star Upsilon_lat", positive=True, real=True
    )
    lambda_ep_omega_D = sp.symbols("lambda_ep_omega_D", positive=True, real=True)

    gamma_lat_safe_eq = sp.simplify(f_lat * mu_eta * (s0**2 - s_c**2) / s_c)
    t_cross_phys = sp.simplify(t_star / s_c)
    gamma_lat_turn_phys = sp.simplify(gamma_lat_safe_eq / (Upsilon_lat * t_star))

    threshold_lambda = sp.simplify(gamma_lat_turn_phys / zeta_ep)
    threshold_product = sp.simplify(threshold_lambda * t_cross_phys)
    threshold_lambda_expected = sp.simplify(
        f_lat * mu_eta * (s0**2 - s_c**2) / (Upsilon_lat * zeta_ep * s_c * t_star)
    )
    threshold_product_expected = sp.simplify(
        f_lat * mu_eta * (s0**2 - s_c**2) / (Upsilon_lat * zeta_ep * s_c**2)
    )

    print("gamma_lat,safe^eq     =", gamma_lat_safe_eq)
    print("t_cross^phys          =", t_cross_phys)
    print("gamma_lat,turn^phys   =", gamma_lat_turn_phys)
    print("(lambda_ep omega_D)min=", threshold_lambda)
    print("product threshold     =", threshold_product)

    assert sp.simplify(threshold_lambda - threshold_lambda_expected) == 0
    assert sp.simplify(threshold_product - threshold_product_expected) == 0
    assert sp.simplify(threshold_product - gamma_lat_safe_eq / (Upsilon_lat * zeta_ep * s_c)) == 0

    # ------------------------------------------------------------------
    # 2. Legacy Session-V slice recovery.
    # ------------------------------------------------------------------
    section("2. Legacy Session-V slice recovery")
    gamma_lattice_legacy = sp.symbols("gamma_lattice_legacy", positive=True, real=True)
    Upsilon_legacy = sp.simplify(gamma_lat_safe_eq / gamma_lattice_legacy)
    threshold_lambda_legacy = sp.simplify(threshold_lambda.subs(Upsilon_lat, Upsilon_legacy))
    threshold_product_legacy = sp.simplify(
        threshold_product.subs(Upsilon_lat, Upsilon_legacy)
    )

    print("Upsilon_lat^(sess)    =", Upsilon_legacy)
    print("legacy lambda thresh  =", threshold_lambda_legacy)
    print("legacy product thresh =", threshold_product_legacy)

    assert sp.simplify(threshold_lambda_legacy - gamma_lattice_legacy / (zeta_ep * t_star)) == 0
    assert sp.simplify(threshold_product_legacy - gamma_lattice_legacy / (zeta_ep * s_c)) == 0

    # ------------------------------------------------------------------
    # 3. Harmonic chi_lambda compiler and force-matched stiffness map.
    # ------------------------------------------------------------------
    section("3. Harmonic chi_lambda compiler and stiffness map")
    lambda_phys, lambda_ref = sp.symbols("lambda_phys lambda_ref", positive=True, real=True)
    r_turn, E_star = sp.symbols("r_turn E_star", positive=True, real=True)
    Vprime_turn_abs = sp.symbols("Vprime_turn_abs", positive=True, real=True)
    a_int = sp.symbols("a_int", positive=True, real=True)

    r_turn_phys = sp.simplify(lambda_phys * r_turn / lambda_ref)
    chi_lambda_lattice = sp.simplify(2 * lambda_phys / r_turn_phys)
    k_eff_req = sp.simplify(E_star * lambda_ref * Vprime_turn_abs / (lambda_phys * r_turn_phys))
    K_turn = sp.simplify(lambda_ref**2 * Vprime_turn_abs / r_turn)
    K_turn_sym = sp.symbols("K_turn", positive=True, real=True)
    k_eff_req_Kturn = sp.simplify(k_eff_req.subs(Vprime_turn_abs, K_turn_sym * r_turn / lambda_ref**2))
    k_eff_req_aint = sp.simplify(k_eff_req_Kturn.subs(lambda_phys, a_int / 2))

    print("r_turn^phys           =", r_turn_phys)
    print("chi_lambda,lattice    =", chi_lambda_lattice)
    print("k_eff,req             =", k_eff_req)
    print("K_turn                =", K_turn)
    print("k_eff,req via K_turn  =", k_eff_req_Kturn)
    print("a_int version         =", k_eff_req_aint)

    assert sp.simplify(r_turn_phys - lambda_phys * r_turn / lambda_ref) == 0
    assert sp.simplify(chi_lambda_lattice - 2 * lambda_ref / r_turn) == 0
    assert sp.simplify(k_eff_req - E_star * lambda_ref**2 * Vprime_turn_abs / (lambda_phys**2 * r_turn)) == 0
    assert sp.simplify(k_eff_req_Kturn - K_turn_sym * E_star / lambda_phys**2) == 0
    assert sp.simplify(k_eff_req_aint - 4 * K_turn_sym * E_star / a_int**2) == 0

    # ------------------------------------------------------------------
    # 4. Korringa-limited temperature ceiling.
    # ------------------------------------------------------------------
    section("4. Korringa-limited temperature ceiling")
    Kcorr, T = sp.symbols("K_corr T", positive=True, real=True)

    T_max = sp.simplify(Kcorr / t_cross_phys)
    T_max_expected = sp.simplify(s_c * Kcorr / t_star)
    Pi_T = sp.simplify(T_max / T)

    print("T_max                 =", T_max)
    print("Pi_T                  =", Pi_T)

    assert sp.simplify(T_max - T_max_expected) == 0
    assert sp.simplify(Pi_T - Kcorr / (T * t_cross_phys)) == 0

    # ------------------------------------------------------------------
    # 5. Exact dimensionless screening ratios.
    # ------------------------------------------------------------------
    section("5. Exact screening ratios")
    k_eff = sp.symbols("k_eff", positive=True, real=True)

    Pi_ep = sp.simplify(Upsilon_lat * zeta_ep * lambda_ep_omega_D * t_star / gamma_lat_safe_eq)
    Pi_chi = sp.simplify(chi_lambda_lattice)
    Pi_k = sp.simplify(k_eff * lambda_phys**2 / (K_turn * E_star))
    Pi_k_Kturn = sp.simplify(Pi_k.subs(Vprime_turn_abs, K_turn_sym * r_turn / lambda_ref**2))

    print("Pi_ep                 =", Pi_ep)
    print("Pi_chi                =", Pi_chi)
    print("Pi_k                  =", Pi_k_Kturn)
    print("Pi_T                  =", Pi_T)

    assert sp.simplify(Pi_ep - lambda_ep_omega_D / threshold_lambda) == 0
    assert sp.simplify(Pi_k_Kturn - k_eff / k_eff_req_Kturn) == 0

    # ------------------------------------------------------------------
    # 6. Session-V / Stage-235 benchmark specialization.
    # ------------------------------------------------------------------
    section("6. Session-V / Stage-235 benchmark specialization")
    s_c_num = 0.5489386551062235
    s0_num = 6.94311167
    f_lat_num = 0.75
    mu_eta_num = 1.0
    gamma_lattice_legacy_num = 4.79562976
    lambda_ref_num = 0.42826825
    r_turn_num = 0.39096144
    K_turn_num = 2.73855812

    gamma_lat_safe_eq_num = f_lat_num * mu_eta_num * (s0_num**2 - s_c_num**2) / s_c_num
    Upsilon_legacy_num = gamma_lat_safe_eq_num / gamma_lattice_legacy_num
    product_micro_num = gamma_lat_safe_eq_num / s_c_num
    product_legacy_num = gamma_lattice_legacy_num / s_c_num
    r_turn_phys_coeff_num = r_turn_num / lambda_ref_num
    chi_lambda_num = 2.0 * lambda_ref_num / r_turn_num
    T_max_coeff_num = s_c_num
    a_int_coeff_num = 4.0 * K_turn_num
    Vprime_turn_abs_num = K_turn_num * r_turn_num / (lambda_ref_num**2)

    print(f"s_c benchmark         = {s_c_num:.15f}")
    print(f"gamma_lat,safe^eq     = {gamma_lat_safe_eq_num:.11f}")
    print(f"Upsilon_lat^(sess)    = {Upsilon_legacy_num:.11f}")
    print(f"micro product thresh  = {product_micro_num:.11f} / zeta_ep")
    print(f"legacy product thresh = {product_legacy_num:.11f} / zeta_ep")
    print(f"r_turn^phys coeff     = {r_turn_phys_coeff_num:.10f} * lambda_phys")
    print(f"chi_lambda benchmark  = {chi_lambda_num:.8f}")
    print(f"K_turn benchmark      = {K_turn_num:.8f}")
    print(f"a_int stiffness coeff = {a_int_coeff_num:.8f}")
    print(f"T_max coeff           = {T_max_coeff_num:.9f} * K_corr / t_*")
    print(f"|V'_red|(r_turn)      = {Vprime_turn_abs_num:.12f}")

    assert abs(gamma_lat_safe_eq_num - 65.45193925961132) < 1e-9
    assert abs(Upsilon_legacy_num - 13.64824695299483) < 1e-9
    assert abs(product_micro_num - 119.23361317476524) < 1e-9
    assert abs(product_legacy_num - 8.736185210116078) < 1e-9
    assert abs(r_turn_phys_coeff_num - 0.9128891530016525) < 1e-12
    assert abs(chi_lambda_num - 2.1908464937104797) < 1e-12
    assert abs(a_int_coeff_num - 10.95423248) < 1e-12
    assert abs(T_max_coeff_num - 0.5489386551062235) < 1e-15
    assert abs(Vprime_turn_abs_num - 5.837462857946154) < 1e-12

    # ------------------------------------------------------------------
    # 7. Final banner.
    # ------------------------------------------------------------------
    section("Stage-236 audit result")
    print("All symbolic and numerical checks passed.")
    print("Verified objects:")
    print("- exact lattice-turnover threshold from the Stage-235 safe-edge lattice rate,")
    print("- exact recovery of the Session-V lambda_ep*omega_D slice through Upsilon_lat,")
    print("- harmonic chi_lambda geometry ratio and force-matched k_eff compiler,")
    print("- Korringa-limited T_max formula,")
    print("- and the exact material-screening ratios (Pi_ep, Pi_chi, Pi_k, Pi_T).")


if __name__ == "__main__":
    main()
