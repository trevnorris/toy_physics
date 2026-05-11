#!/usr/bin/env python3
"""Numerical stress harness for the Stage 253 material-threshold companion.

This harness starts from primitive host inputs and verifies the Stage 253
screening stack forward. It deliberately avoids back-solving host values from
target `Pi_*` margins.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


FLAG_ORDER = ("Pi_ep", "Pi_chi", "Pi_k", "Pi_t")


def fmt(value: float) -> str:
    return f"{value:.12g}"


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def near(lhs: float, rhs: float, tol: float) -> bool:
    return abs(lhs - rhs) <= tol * (1.0 + abs(rhs))


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage253_v2":
        raise ValueError("unexpected sample schema")
    return data


def resolve_upsilon(case: dict, gamma_lat_safe_eq: float) -> float:
    mode = case["mode"]
    if mode == "micro":
        return 1.0
    if mode == "legacy":
        gamma_legacy = float(case["gamma_lattice_legacy"])
        return gamma_lat_safe_eq / gamma_legacy
    if mode == "custom":
        return float(case["Upsilon_lat"])
    raise ValueError(f"unexpected mode: {mode}")


def ordered_flags(raw: dict) -> dict[str, bool]:
    return {name: bool(raw[name]) for name in FLAG_ORDER}


def passes_threshold(value: float, tol: float) -> bool:
    return value >= 1.0 - tol


def stage253_block(config: dict) -> None:
    print("\n=== Stage 253: material-threshold stress ===")
    ratio_tol = float(config["tolerances"]["ratio_tol"])
    consistency_tol = float(config["tolerances"].get("consistency_tol", ratio_tol))

    for case in config["cases"]:
        s_c = float(case["s_c"])
        s0 = float(case["s_0"])
        f_lat = float(case["f_lat"])
        mu_eta = float(case["mu_eta"])
        zeta_ep = float(case["zeta_ep"])
        t_star = float(case["t_star"])
        lambda_phys = float(case["lambda_phys"])
        lambda_ref = float(case["lambda_ref"])
        e_star = float(case["E_star"])
        k_turn = float(case["K_turn"])
        k_corr = float(case["K_corr"])

        host = case["host"]
        lambda_ep_omega_d = float(host["lambda_ep_omega_D"])
        r_turn = float(host["r_turn"])
        r_turn_phys = float(host["r_turn_phys"])
        vprime_turn_abs = float(host["Vprime_turn_abs"])
        k_eff = float(host["k_eff"])
        temperature = float(host["temperature"])

        expected = ordered_flags(case["expected_flags"])

        gamma_lat_safe_eq = f_lat * mu_eta * (s0 * s0 - s_c * s_c) / s_c
        upsilon_lat = resolve_upsilon(case, gamma_lat_safe_eq)
        threshold_lambda = gamma_lat_safe_eq / (upsilon_lat * zeta_ep * t_star)
        t_cross_phys = t_star / s_c
        threshold_product = threshold_lambda * t_cross_phys
        t_max = k_corr / t_cross_phys

        r_turn_phys_expected = lambda_phys * r_turn / lambda_ref
        k_turn_from_vprime = lambda_ref * lambda_ref * vprime_turn_abs / r_turn
        k_eff_req_direct = e_star * lambda_ref * vprime_turn_abs / (lambda_phys * r_turn_phys)
        k_eff_req_reduced = k_turn * e_star / (lambda_phys * lambda_phys)

        pi_ep_direct = upsilon_lat * zeta_ep * lambda_ep_omega_d * t_star / gamma_lat_safe_eq
        pi_ep_ratio = lambda_ep_omega_d / threshold_lambda
        pi_chi_phys = 2.0 * lambda_phys / r_turn_phys
        pi_chi_reduced = 2.0 * lambda_ref / r_turn
        pi_k_direct = k_eff * lambda_phys * lambda_phys / (k_turn * e_star)
        pi_k_ratio = k_eff / k_eff_req_reduced
        pi_t_direct = k_corr / (temperature * t_cross_phys)
        pi_t_ratio = t_max / temperature

        actual = {
            "Pi_ep": passes_threshold(pi_ep_direct, ratio_tol),
            "Pi_chi": passes_threshold(pi_chi_phys, ratio_tol),
            "Pi_k": passes_threshold(pi_k_direct, ratio_tol),
            "Pi_t": passes_threshold(pi_t_direct, ratio_tol),
        }

        print(
            f"\n{case['name']} ({case['mode']}): gamma_safe={fmt(gamma_lat_safe_eq)}, "
            f"Upsilon={fmt(upsilon_lat)}, lambda_min={fmt(threshold_lambda)}"
        )
        print(
            "  host inputs: "
            f"lambda_ep*omega_D={fmt(lambda_ep_omega_d)}, r_turn={fmt(r_turn)}, "
            f"r_turn_phys={fmt(r_turn_phys)}, |V'|={fmt(vprime_turn_abs)}, "
            f"k_eff={fmt(k_eff)}, T={fmt(temperature)}"
        )
        print(
            "  computed margins: "
            f"Pi_ep={fmt(pi_ep_direct)}, Pi_chi={fmt(pi_chi_phys)}, "
            f"Pi_k={fmt(pi_k_direct)}, Pi_t={fmt(pi_t_direct)}"
        )

        require(
            f"{case['name']} primitive admissibility",
            s0 > s_c > 0.0 and f_lat > 0.0 and mu_eta > 0.0 and zeta_ep > 0.0 and upsilon_lat > 0.0,
            (
                f"s0={fmt(s0)}, s_c={fmt(s_c)}, f_lat={fmt(f_lat)}, mu_eta={fmt(mu_eta)}, "
                f"zeta_ep={fmt(zeta_ep)}, Upsilon={fmt(upsilon_lat)}"
            ),
        )
        require(
            f"{case['name']} host-input positivity",
            lambda_ep_omega_d > 0.0 and r_turn > 0.0 and r_turn_phys > 0.0 and vprime_turn_abs > 0.0 and k_eff > 0.0 and temperature > 0.0,
            (
                f"lambda_ep*omega_D={fmt(lambda_ep_omega_d)}, r_turn={fmt(r_turn)}, "
                f"r_turn_phys={fmt(r_turn_phys)}, |V'|={fmt(vprime_turn_abs)}, "
                f"k_eff={fmt(k_eff)}, T={fmt(temperature)}"
            ),
        )
        require(
            f"{case['name']} turnover threshold positivity",
            gamma_lat_safe_eq > 0.0 and threshold_lambda > 0.0 and threshold_product > 0.0,
            (
                f"gamma_safe={fmt(gamma_lat_safe_eq)}, lambda_min={fmt(threshold_lambda)}, "
                f"product_min={fmt(threshold_product)}"
            ),
        )
        require(
            f"{case['name']} turning-radius dictionary",
            near(r_turn_phys, r_turn_phys_expected, consistency_tol),
            f"r_turn_phys={fmt(r_turn_phys)}, expected={fmt(r_turn_phys_expected)}",
        )
        require(
            f"{case['name']} K_turn / Vprime consistency",
            near(k_turn_from_vprime, k_turn, consistency_tol),
            f"K_turn_from_Vprime={fmt(k_turn_from_vprime)}, K_turn={fmt(k_turn)}",
        )
        require(
            f"{case['name']} stiffness compiler consistency",
            near(k_eff_req_direct, k_eff_req_reduced, consistency_tol),
            f"k_req_direct={fmt(k_eff_req_direct)}, k_req_reduced={fmt(k_eff_req_reduced)}",
        )
        require(
            f"{case['name']} product-threshold factorization",
            near(threshold_product, gamma_lat_safe_eq / (upsilon_lat * zeta_ep * s_c), consistency_tol),
            f"product={fmt(threshold_product)}, reduced={fmt(gamma_lat_safe_eq / (upsilon_lat * zeta_ep * s_c))}",
        )
        require(
            f"{case['name']} Pi_ep direct vs threshold ratio",
            near(pi_ep_direct, pi_ep_ratio, ratio_tol),
            f"direct={fmt(pi_ep_direct)}, ratio={fmt(pi_ep_ratio)}",
        )
        require(
            f"{case['name']} Pi_chi physical vs reduced ratio",
            near(pi_chi_phys, pi_chi_reduced, ratio_tol),
            f"physical={fmt(pi_chi_phys)}, reduced={fmt(pi_chi_reduced)}",
        )
        require(
            f"{case['name']} Pi_k direct vs reduced ratio",
            near(pi_k_direct, pi_k_ratio, ratio_tol),
            f"direct={fmt(pi_k_direct)}, ratio={fmt(pi_k_ratio)}",
        )
        require(
            f"{case['name']} Pi_t direct vs Tmax ratio",
            near(pi_t_direct, pi_t_ratio, ratio_tol),
            f"direct={fmt(pi_t_direct)}, ratio={fmt(pi_t_ratio)}",
        )
        require(
            f"{case['name']} screening classification",
            actual == expected,
            f"actual={actual}, expected={expected}",
        )
        require(
            f"{case['name']} overall screen verdict",
            all(actual.values()) == all(expected.values()),
            f"actual={actual}, expected={expected}",
        )

        if case["mode"] == "legacy":
            gamma_legacy = float(case["gamma_lattice_legacy"])
            require(
                f"{case['name']} legacy Upsilon recovery",
                near(upsilon_lat, gamma_lat_safe_eq / gamma_legacy, consistency_tol),
                f"Upsilon={fmt(upsilon_lat)}, gamma_safe/gamma_legacy={fmt(gamma_lat_safe_eq / gamma_legacy)}",
            )
            require(
                f"{case['name']} legacy lambda threshold",
                near(threshold_lambda, gamma_legacy / (zeta_ep * t_star), consistency_tol),
                f"threshold={fmt(threshold_lambda)}, legacy={fmt(gamma_legacy / (zeta_ep * t_star))}",
            )
            require(
                f"{case['name']} legacy product threshold",
                near(threshold_product, gamma_legacy / (zeta_ep * s_c), consistency_tol),
                f"product={fmt(threshold_product)}, legacy={fmt(gamma_legacy / (zeta_ep * s_c))}",
            )


def main(argv: list[str]) -> int:
    default_config = Path(__file__).resolve().with_name("stage253_material_threshold_samples.json")
    config_path = Path(argv[1]) if len(argv) > 1 else default_config
    config = load_config(config_path)
    print(f"Loaded config from {config_path}")
    stage253_block(config)
    print("\nAll Stage 253 numerical stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
