#!/usr/bin/env python3
"""Foundational resonance/outgoing stress checks for Stages 003 and 004."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


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
    if data.get("schema") != "moving_throat_numerical_stage003_004_v1":
        raise ValueError("unexpected sample schema")
    return data


def stage003_block(config: dict) -> None:
    print("\n=== Stage 003: one-mode resonance stress ===")
    tol = config["stage003"]["tolerances"]
    for case in config["stage003"]["cases"]:
        m = float(case["M"])
        om2 = float(case["Omega_eta2"])
        delta = float(case["delta"])
        g = float(case["g"])
        varpi2 = om2 + delta
        disc = math.sqrt((om2 - varpi2) ** 2 + 4.0 * g * g / m)
        minus = 0.5 * (om2 + varpi2 - disc)
        plus = 0.5 * (om2 + varpi2 + disc)
        approx_minus = om2 - g * g / (m * delta)
        approx_plus = varpi2 + g * g / (m * delta)
        wall_shift = minus - om2
        matter_shift = plus - varpi2
        lambda_ratio = g * g / (m * delta * delta)

        print(f"\n{case['name']} ({case['kind']}): M={fmt(m)}, Omega_eta^2={fmt(om2)}, delta={fmt(delta)}, g={fmt(g)}")
        print(f"  resonance parameter g^2/(M delta^2) = {fmt(lambda_ratio)}")

        require(
            f"{case['name']} exact root ordering",
            minus < om2 < varpi2 < plus,
            f"omega_-^2={fmt(minus)}, Omega_eta^2={fmt(om2)}, varpi^2={fmt(varpi2)}, omega_+^2={fmt(plus)}",
        )
        require(
            f"{case['name']} exact sum of roots",
            near(minus + plus, om2 + varpi2, tol["sum_product_tol"]),
            f"sum={fmt(minus + plus)}, target={fmt(om2 + varpi2)}",
        )
        require(
            f"{case['name']} exact product of roots",
            near(minus * plus, om2 * varpi2 - g * g / m, tol["sum_product_tol"]),
            f"product={fmt(minus * plus)}, target={fmt(om2 * varpi2 - g * g / m)}",
        )
        require(
            f"{case['name']} wall/matter shift signs",
            wall_shift < 0.0 and matter_shift > 0.0,
            f"wall_shift={fmt(wall_shift)}, matter_shift={fmt(matter_shift)}",
        )

        wall_rel = abs(wall_shift - (approx_minus - om2)) / abs(approx_minus - om2)
        matter_rel = abs(matter_shift - (approx_plus - varpi2)) / abs(approx_plus - varpi2)
        print(
            f"  exact shifts: wall={fmt(wall_shift)}, matter={fmt(matter_shift)}; "
            f"perturbative: wall={fmt(approx_minus - om2)}, matter={fmt(approx_plus - varpi2)}"
        )
        if case["regime"] == "perturbative_valid":
            tol_rel = float(case["perturbative_rel_tol"])
            require(
                f"{case['name']} wall-like perturbative validity",
                wall_rel <= tol_rel,
                f"relative error={fmt(wall_rel)}, tol={fmt(tol_rel)}",
            )
            require(
                f"{case['name']} matter-like perturbative validity",
                matter_rel <= tol_rel,
                f"relative error={fmt(matter_rel)}, tol={fmt(tol_rel)}",
            )
        else:
            floor = float(case["breakdown_floor"])
            require(
                f"{case['name']} perturbative breakdown is visible",
                wall_rel >= floor and matter_rel >= floor,
                f"wall_rel={fmt(wall_rel)}, matter_rel={fmt(matter_rel)}, floor={fmt(floor)}",
            )


def stage004_block(config: dict) -> None:
    print("\n=== Stage 004: outgoing mixed-sector stress ===")

    def h1(z: float) -> complex:
        j1 = math.sin(z) / (z * z) - math.cos(z) / z
        y1 = -math.cos(z) / (z * z) - math.sin(z) / z
        return complex(j1, y1)

    def h2(z: float) -> complex:
        j2 = ((3.0 / z**3) - 1.0 / z) * math.sin(z) - 3.0 * math.cos(z) / z**2
        y2 = -((3.0 / z**3) - 1.0 / z) * math.cos(z) - 3.0 * math.sin(z) / z**2
        return complex(j2, y2)

    def y2hat(radius: float, c_s: float, omega: float) -> complex:
        z = omega * radius / c_s
        hh2 = h2(z)
        hh2p = h1(z) - 3.0 * hh2 / z
        lam = (omega / c_s) * hh2p / hh2
        y2_exact = 1.0 / lam
        return y2_exact / (-radius / 3.0)

    for case in config["stage004"]["cases"]:
        oa = float(case["Omega_A"])
        ow = float(case["Omega_W"])
        r = float(case["R"])
        ga = float(case["g_A"])
        gw = float(case["g_W"])
        radius = float(case["radius"])
        c_s = float(case["c_s"])
        denom0 = oa * oa * ow * ow - r * r
        n0 = (oa * oa * gw + r * ga) ** 2 / (denom0 * denom0)
        gamma5 = radius**5 / (27.0 * c_s**5)

        print(f"\n{case['name']} ({case['kind']}): Omega_A={fmt(oa)}, Omega_W={fmt(ow)}, R={fmt(r)}")
        print(f"  denominator Omega_A^2 Omega_W^2 - R^2 = {fmt(denom0)}")
        require(
            f"{case['name']} positive mixed denominator",
            denom0 > 0.0,
            f"denominator={fmt(denom0)}",
        )
        for omega in case["omegas"]:
            omega = float(omega)
            aker = oa * oa - omega * omega
            wker = ow * ow - omega * omega
            delta = aker * wker - r * r
            n_exact = (aker * gw + r * ga) ** 2 / (delta * delta)
            sigma_cons = (ga * ga * wker + 2.0 * ga * gw * r + gw * gw * aker) / delta
            pi_out = 1j * gamma5 * omega**5
            sigma_full = (ga * ga * (wker - pi_out) + 2.0 * ga * gw * r + gw * gw * aker) / (aker * (wker - pi_out) - r * r)
            d_corr = -(sigma_full - sigma_cons)
            y2_exact = y2hat(radius, c_s, omega)
            y2_even = 1.0 + radius * radius * omega * omega / (9.0 * c_s * c_s) + 4.0 * radius**4 * omega**4 / (81.0 * c_s**4)
            n_rel = abs(n_exact - n0) / (1.0 + abs(n0))
            dcorr_rel = abs((d_corr.imag / omega**5) - (-gamma5 * n0)) / (1.0 + abs(gamma5 * n0))
            y2_imag_rel = abs((y2_exact.imag / omega**5) - gamma5) / (1.0 + abs(gamma5))
            y2_real_abs = abs(y2_exact.real - y2_even)

            print(
                f"  omega={fmt(omega)}: N(omega)={fmt(n_exact)}, "
                f"Im(delta D)/omega^5={fmt(d_corr.imag / omega**5)}, "
                f"Im(Y2hat)/omega^5={fmt(y2_exact.imag / omega**5)}"
            )
            require(
                f"{case['name']} N0 convergence at omega={fmt(omega)}",
                n_rel <= float(case["N0_rel_tol"]),
                f"relative error={fmt(n_rel)}, tol={fmt(float(case['N0_rel_tol']))}",
            )
            require(
                f"{case['name']} wall odd coefficient at omega={fmt(omega)}",
                dcorr_rel <= float(case["Dcorr_rel_tol"]),
                f"relative error={fmt(dcorr_rel)}, tol={fmt(float(case['Dcorr_rel_tol']))}",
            )
            require(
                f"{case['name']} Y2hat odd coefficient at omega={fmt(omega)}",
                y2_imag_rel <= float(case["Y2_imag_rel_tol"]),
                f"relative error={fmt(y2_imag_rel)}, tol={fmt(float(case['Y2_imag_rel_tol']))}",
            )
            require(
                f"{case['name']} Y2hat even branch at omega={fmt(omega)}",
                y2_real_abs <= float(case["Y2_real_abs_tol"]),
                f"abs error={fmt(y2_real_abs)}, tol={fmt(float(case['Y2_real_abs_tol']))}",
            )


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[3]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/moving_throat/numerical/stage003_004_foundational_samples.json"
    config = load_config(config_path)
    print(f"Loaded config from {config_path}")
    stage003_block(config)
    stage004_block(config)
    print("\nAll stage-003/004 foundational stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
