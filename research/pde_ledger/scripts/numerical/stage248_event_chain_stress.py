#!/usr/bin/env python3
"""Numerical stress harness for the Stage 248 event-chain compiler."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp


def fmt(value: float | mp.mpf) -> str:
    if isinstance(value, mp.mpc):
        value = mp.re(value)
    return f"{float(value):.12g}"


def realify(value: mp.mpf | mp.mpc) -> mp.mpf:
    if isinstance(value, mp.mpc):
        if abs(mp.im(value)) > mp.mpf("1.0e-20"):
            raise ValueError(f"unexpected complex residual: {value}")
        return mp.re(value)
    return value


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def near(lhs: float | mp.mpf, rhs: float | mp.mpf, tol: float) -> bool:
    if isinstance(lhs, mp.mpc):
        lhs = mp.re(lhs)
    if isinstance(rhs, mp.mpc):
        rhs = mp.re(rhs)
    lhs_f = float(lhs)
    rhs_f = float(rhs)
    return abs(lhs_f - rhs_f) <= tol * (1.0 + abs(rhs_f))


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage248_v2":
        raise ValueError("unexpected sample schema")
    return data


def lowered_barrier_potential(radius: mp.mpf, strength: mp.mpf, scale: mp.mpf) -> mp.mpf:
    return 1 / radius - strength * mp.exp(-radius / scale)


def lowered_barrier_prime(radius: mp.mpf, strength: mp.mpf, scale: mp.mpf) -> mp.mpf:
    return -1 / radius**2 + strength * mp.exp(-radius / scale) / scale


def lowered_barrier_double_prime(radius: mp.mpf, strength: mp.mpf, scale: mp.mpf) -> mp.mpf:
    return 2 / radius**3 - strength * mp.exp(-radius / scale) / scale**2


def xi_profile(radius: mp.mpf, xi_floor: mp.mpf, xi_amp: mp.mpf, xi_scale: mp.mpf) -> mp.mpf:
    return xi_floor + xi_amp * mp.exp(-radius / xi_scale)


def bisect_root(function, left: mp.mpf, right: mp.mpf, iterations: int = 80) -> mp.mpf:
    f_left = function(left)
    f_right = function(right)
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right > 0:
        raise ValueError("interval does not bracket a root")
    for _ in range(iterations):
        midpoint = (left + right) / 2
        f_mid = function(midpoint)
        if f_left * f_mid <= 0:
            right = midpoint
            f_right = f_mid
        else:
            left = midpoint
            f_left = f_mid
    return (left + right) / 2


def scan_roots(function, lower: mp.mpf, upper: mp.mpf, samples: int = 4000) -> list[mp.mpf]:
    roots: list[mp.mpf] = []
    step = (upper - lower) / samples
    left = lower
    f_left = function(left)
    for index in range(1, samples + 1):
        right = lower + step * index
        f_right = function(right)
        root: mp.mpf | None = None
        if f_left == 0:
            root = left
        elif f_left * f_right < 0:
            root = bisect_root(function, left, right)
        if root is not None and lower < root < upper:
            if all(abs(root - existing) > mp.mpf("1.0e-8") for existing in roots):
                roots.append(root)
        left = right
        f_left = f_right
    return sorted(roots)


def coulomb_action_closed(m_s: mp.mpf, hbar_eff: mp.mpf, energy: mp.mpf, r_contact: mp.mpf) -> mp.mpf:
    return mp.sqrt(2 * m_s) / hbar_eff * (
        mp.pi / (2 * mp.sqrt(energy))
        - mp.sqrt(r_contact * (1 - energy * r_contact))
        - mp.asin(mp.sqrt(energy * r_contact)) / mp.sqrt(energy)
    )


def coulomb_action_quad(m_s: mp.mpf, hbar_eff: mp.mpf, energy: mp.mpf, r_contact: mp.mpf) -> mp.mpf:
    turning_point = 1 / energy

    def integrand(radius: mp.mpf) -> mp.mpf:
        radicand = 2 * m_s * (1 / radius - energy)
        if radicand < 0 and abs(radicand) < mp.mpf("1.0e-30"):
            radicand = mp.mpf("0.0")
        return mp.sqrt(radicand) / hbar_eff

    midpoint = (r_contact + turning_point) / 2
    return mp.quad(integrand, [r_contact, midpoint, turning_point])


def lowered_action_quad(
    m_s: mp.mpf,
    hbar_eff: mp.mpf,
    energy: mp.mpf,
    strength: mp.mpf,
    scale: mp.mpf,
    r_minus: mp.mpf,
    r_peak: mp.mpf,
    r_plus: mp.mpf,
) -> mp.mpf:
    def integrand(radius: mp.mpf) -> mp.mpf:
        radicand = 2 * m_s * (lowered_barrier_potential(radius, strength, scale) - energy)
        if radicand < 0 and abs(radicand) < mp.mpf("1.0e-30"):
            radicand = mp.mpf("0.0")
        return mp.sqrt(radicand) / hbar_eff

    return mp.quad(integrand, [r_minus, r_peak, r_plus])


def classify_barrier_case(case: dict, tolerances: dict) -> dict:
    m_s = mp.mpf(case["m_s"])
    hbar_eff = mp.mpf(case["hbar_eff"])
    r0 = mp.mpf(case["r0"])
    r_contact = mp.mpf(case["r_contact"])
    strength = mp.mpf(case["lowering_strength"])
    scale = mp.mpf(case["lowering_scale"])
    margin = mp.mpf(case["energy_margin_fraction"])
    xi_floor = mp.mpf(case["xi_floor"])
    xi_amp = mp.mpf(case["xi_amp"])
    xi_scale = mp.mpf(case["xi_scale"])
    epsilon = mp.mpf("1.0e-5")

    result: dict[str, object] = {
        "name": case["name"],
        "status": "failure",
        "reason": "unclassified",
    }

    peak_roots = scan_roots(
        lambda radius: lowered_barrier_prime(radius, strength, scale),
        r_contact + epsilon,
        r0 - epsilon,
    )
    result["peak_roots"] = peak_roots
    if len(peak_roots) == 0:
        result["reason"] = "no_peak"
        return result
    if len(peak_roots) > 1:
        result["reason"] = "multiple_stationary_points"
        return result

    r_peak = peak_roots[0]
    v_launch = lowered_barrier_potential(r0, strength, scale)
    v_contact = lowered_barrier_potential(r_contact, strength, scale)
    v_peak = lowered_barrier_potential(r_peak, strength, scale)
    k_peak = -lowered_barrier_double_prime(r_peak, strength, scale)

    result["r_peak"] = r_peak
    result["V0"] = v_launch
    result["V_contact"] = v_contact
    result["V_peak"] = v_peak
    result["K_peak"] = k_peak

    if not (v_launch < v_contact < v_peak):
        result["reason"] = "contact_not_below_peak"
        return result
    if not (mp.mpf("0.0") < margin < mp.mpf("1.0")):
        result["reason"] = "invalid_energy_margin"
        return result

    energy = v_contact + margin * (v_peak - v_contact)
    result["E_sub"] = energy
    if energy * r_contact >= 1:
        result["reason"] = "coulomb_inadmissible"
        return result

    turning_roots = scan_roots(
        lambda radius: lowered_barrier_potential(radius, strength, scale) - energy,
        r_contact + epsilon,
        r0 - epsilon,
    )
    result["turning_roots"] = turning_roots
    if len(turning_roots) != 2:
        result["reason"] = "missing_turning_pair"
        return result

    r_minus, r_plus = turning_roots
    if not (r_contact < r_minus < r_peak < r_plus < r0):
        result["reason"] = "turning_pair_ordering"
        return result

    v_crit = mp.sqrt(2 * (v_peak - v_launch) / m_s)
    v_contact_coul = mp.sqrt(2 * (1 / r_contact - 1 / r0) / m_s)
    v_sub = mp.sqrt(2 * (energy - v_launch) / m_s)
    v_cross = v_crit + mp.mpf(case["v_cross_fraction"]) * (v_contact_coul - v_crit)
    i_new = realify(
        lowered_action_quad(m_s, hbar_eff, energy, strength, scale, r_minus, r_peak, r_plus)
    )
    i_coul_closed = coulomb_action_closed(m_s, hbar_eff, energy, r_contact)
    i_coul_quad = coulomb_action_quad(m_s, hbar_eff, energy, r_contact)
    transmission_ratio = realify(mp.exp(-2 * (i_new - i_coul_closed)))
    xi_turn = realify(xi_profile(r_plus, xi_floor, xi_amp, xi_scale))
    xi_contact = realify(xi_profile(r_contact, xi_floor, xi_amp, xi_scale))
    xi_launch = realify(xi_profile(r0, xi_floor, xi_amp, xi_scale))
    slope_exact = realify(abs(lowered_barrier_prime(r_plus, strength, scale)))

    h_step = min(
        (r_plus - r_peak) / 80,
        (r0 - r_plus) / 80,
        (r_plus - r_minus) / 160,
        r_plus / 500,
    )
    derivative_fd = realify(abs(
        (
            -lowered_barrier_potential(r_plus + 2 * h_step, strength, scale)
            + 8 * lowered_barrier_potential(r_plus + h_step, strength, scale)
            - 8 * lowered_barrier_potential(r_plus - h_step, strength, scale)
            + lowered_barrier_potential(r_plus - 2 * h_step, strength, scale)
        )
        / (12 * h_step)
    ))
    lambda_th = energy / slope_exact

    transport_window = min(energy - v_contact, v_peak - energy)
    transport_delta = mp.mpf("0.05") * transport_window

    upper_turnings = scan_roots(
        lambda radius: lowered_barrier_potential(radius, strength, scale) - (energy + transport_delta),
        r_contact + epsilon,
        r0 - epsilon,
    )
    upper2_turnings = scan_roots(
        lambda radius: lowered_barrier_potential(radius, strength, scale) - (energy + 2 * transport_delta),
        r_contact + epsilon,
        r0 - epsilon,
    )
    lower_turnings = scan_roots(
        lambda radius: lowered_barrier_potential(radius, strength, scale) - (energy - transport_delta),
        r_contact + epsilon,
        r0 - epsilon,
    )
    lower2_turnings = scan_roots(
        lambda radius: lowered_barrier_potential(radius, strength, scale) - (energy - 2 * transport_delta),
        r_contact + epsilon,
        r0 - epsilon,
    )
    if (
        len(upper2_turnings) != 2
        or len(upper_turnings) != 2
        or len(lower_turnings) != 2
        or len(lower2_turnings) != 2
    ):
        result["reason"] = "transport_window_breakdown"
        return result
    transport_fd = realify(
        (
            -upper2_turnings[1]
            + 8 * upper_turnings[1]
            - 8 * lower_turnings[1]
            + lower2_turnings[1]
        )
        / (12 * transport_delta)
    )
    transport_exact = realify(1 / lowered_barrier_prime(r_plus, strength, scale))

    result.update(
        {
            "status": "success",
            "reason": "success",
            "r_minus": r_minus,
            "r_plus": r_plus,
            "v_sub": v_sub,
            "v_crit": v_crit,
            "v_cross": v_cross,
            "v_contact_coul": v_contact_coul,
            "I_new": i_new,
            "I_coul_closed": i_coul_closed,
            "I_coul_quad": i_coul_quad,
            "transmission_ratio": transmission_ratio,
            "Xi_turn": xi_turn,
            "Xi_contact": xi_contact,
            "Xi_launch": xi_launch,
            "lambda_th": lambda_th,
            "dV_turn_exact": slope_exact,
            "dV_turn_fd": derivative_fd,
            "transport_fd": transport_fd,
            "transport_exact": transport_exact,
        }
    )
    return result


def stage248_event_chain_block(config: dict) -> None:
    print("\n=== Stage 248: dynamic event-chain stress ===")
    action_tol = float(config["tolerances"]["action_rel_tol"])
    derivative_tol = float(config["tolerances"]["derivative_rel_tol"])
    transport_tol = float(config["tolerances"]["transport_rel_tol"])

    for case in config["barrier_cases"]:
        result = classify_barrier_case(case, config["tolerances"])
        expected_status = case["expected_status"]
        expected_reason = case.get("expected_reason")

        print(f"\n{case['name']}: expected={expected_status}")
        print(
            f"  raw barrier data: strength={fmt(case['lowering_strength'])}, "
            f"scale={fmt(case['lowering_scale'])}, r_contact={fmt(case['r_contact'])}, "
            f"r0={fmt(case['r0'])}"
        )
        print(
            f"  classification: status={result['status']}, reason={result['reason']}"
        )

        if expected_status == "failure":
            require(
                f"{case['name']} expected failure status",
                result["status"] == "failure",
                f"status={result['status']}",
            )
            require(
                f"{case['name']} expected failure reason",
                result["reason"] == expected_reason,
                f"reason={result['reason']}, expected={expected_reason}",
            )
            continue

        require(
            f"{case['name']} expected success status",
            result["status"] == "success",
            f"status={result['status']}, reason={result['reason']}",
        )

        print(
            f"  peak/turning data: r_peak={fmt(result['r_peak'])}, "
            f"r_minus={fmt(result['r_minus'])}, r_plus={fmt(result['r_plus'])}"
        )
        print(
            f"  barrier heights: V0={fmt(result['V0'])}, V_contact={fmt(result['V_contact'])}, "
            f"V_peak={fmt(result['V_peak'])}, Esub={fmt(result['E_sub'])}"
        )
        print(
            f"  actions: I_new={fmt(result['I_new'])}, I_coul_quad={fmt(result['I_coul_quad'])}, "
            f"I_coul_closed={fmt(result['I_coul_closed'])}, T_ratio={fmt(result['transmission_ratio'])}"
        )
        print(
            f"  diagnostics: Xi_turn={fmt(result['Xi_turn'])}, lambda_th={fmt(result['lambda_th'])}, "
            f"|V'(r_+)|={fmt(result['dV_turn_exact'])}"
        )

        require(
            f"{case['name']} dynamic corridor ordering",
            result["v_sub"] < result["v_crit"] < result["v_cross"] < result["v_contact_coul"],
            (
                f"v_sub={fmt(result['v_sub'])}, v_crit={fmt(result['v_crit'])}, "
                f"v_cross={fmt(result['v_cross'])}, v_contact={fmt(result['v_contact_coul'])}"
            ),
        )
        require(
            f"{case['name']} Coulomb action quadrature",
            near(result["I_coul_quad"], result["I_coul_closed"], action_tol),
            f"quad={fmt(result['I_coul_quad'])}, closed={fmt(result['I_coul_closed'])}",
        )
        require(
            f"{case['name']} lowered branch beats Coulomb",
            result["I_new"] < result["I_coul_closed"] and result["transmission_ratio"] > 1,
            (
                f"I_new={fmt(result['I_new'])}, I_coul={fmt(result['I_coul_closed'])}, "
                f"ratio={fmt(result['transmission_ratio'])}"
            ),
        )
        require(
            f"{case['name']} Xi_turn derived from monotone profile",
            result["Xi_launch"] < result["Xi_turn"] < result["Xi_contact"],
            (
                f"Xi_launch={fmt(result['Xi_launch'])}, Xi_turn={fmt(result['Xi_turn'])}, "
                f"Xi_contact={fmt(result['Xi_contact'])}"
            ),
        )
        require(
            f"{case['name']} derivative at r_plus",
            near(result["dV_turn_fd"], result["dV_turn_exact"], derivative_tol),
            f"fd={fmt(result['dV_turn_fd'])}, exact={fmt(result['dV_turn_exact'])}",
        )
        require(
            f"{case['name']} turning transport law",
            near(result["transport_fd"], result["transport_exact"], transport_tol),
            f"fd={fmt(result['transport_fd'])}, exact={fmt(result['transport_exact'])}",
        )


def parabolic_action_closed(m_s: float, hbar_eff: float, delta_e: float, k_peak: float) -> mp.mpf:
    m = mp.mpf(m_s)
    hbar = mp.mpf(hbar_eff)
    delta = mp.mpf(delta_e)
    k_val = mp.mpf(k_peak)
    return mp.pi * delta * mp.sqrt(m / k_val) / hbar


def parabolic_action_quad(m_s: float, hbar_eff: float, delta_e: float, k_peak: float) -> mp.mpf:
    m = mp.mpf(m_s)
    hbar = mp.mpf(hbar_eff)
    delta = mp.mpf(delta_e)
    k_val = mp.mpf(k_peak)
    y_turn = mp.sqrt(2 * delta / k_val)

    def integrand(y_val: mp.mpf) -> mp.mpf:
        radicand = 2 * m * (delta - k_val * y_val * y_val / 2)
        if radicand < 0 and abs(radicand) < mp.mpf("1.0e-30"):
            radicand = mp.mpf("0.0")
        return mp.sqrt(radicand) / hbar

    return mp.quad(integrand, [-y_turn, mp.mpf("0.0"), y_turn])


def near_top_block(config: dict) -> None:
    print("\n=== Stage 248: near-top parabolic action stress ===")
    tol = float(config["tolerances"]["near_top_rel_tol"])

    for case in config["near_top_cases"]:
        m_s = float(case["m_s"])
        hbar_eff = float(case["hbar_eff"])
        delta_e = float(case["DeltaE"])
        k_peak = float(case["Kpeak"])
        closed = parabolic_action_closed(m_s, hbar_eff, delta_e, k_peak)
        quad = parabolic_action_quad(m_s, hbar_eff, delta_e, k_peak)
        y_turn = mp.sqrt(2 * mp.mpf(delta_e) / mp.mpf(k_peak))

        print(
            f"\n{case['name']}: m_s={fmt(m_s)}, hbar_eff={fmt(hbar_eff)}, "
            f"DeltaE={fmt(delta_e)}, Kpeak={fmt(k_peak)}, y_turn={fmt(y_turn)}"
        )
        require(
            f"{case['name']} near-top quadrature",
            near(quad, closed, tol),
            f"quad={fmt(quad)}, closed={fmt(closed)}",
        )


def main(argv: list[str]) -> int:
    mp.mp.dps = 60
    default_config = Path(__file__).resolve().with_name("stage248_event_chain_samples.json")
    config_path = Path(argv[1]) if len(argv) > 1 else default_config
    config = load_config(config_path)
    print(f"Loaded config from {config_path}")
    stage248_event_chain_block(config)
    near_top_block(config)
    print("\nAll Stage 248 numerical stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
