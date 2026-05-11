#!/usr/bin/env python3
"""Foundational numerical spot checks for Stages 041 and 153."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


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
    if data.get("schema") != "moving_throat_numerical_stage058_170_v1":
        raise ValueError("unexpected sample schema")
    return data


def stage058_block(config: dict) -> None:
    print("\n=== Stage 058: coupled support/source fixed point ===")
    block = config["stage058"]
    tol = block["tolerances"]
    quadrature_points = int(block["quadrature_points"])
    bisection_iterations = int(block["bisection_iterations"])
    x = np.linspace(0.0, 1.0, quadrature_points)

    def w(alpha: float, eta: float) -> float:
        return alpha * math.sinh(alpha) + eta * math.cosh(alpha)

    def kernel(alpha: float, eta: float, xvals: np.ndarray) -> np.ndarray:
        denom = w(alpha, eta)
        return (
            np.cosh(alpha * xvals)
            + (eta / alpha) * np.sinh(alpha * xvals)
            - np.cosh(alpha * (1.0 - xvals))
        ) / denom

    def sigma(pe: float, xvals: np.ndarray) -> np.ndarray:
        return pe * np.exp(pe * xvals) / (math.exp(pe) - 1.0)

    def delta_formula(pe: float, alpha: float, eta: float) -> float:
        denom = pe * pe - alpha * alpha
        ic = (math.exp(pe) * (pe * math.cosh(alpha) - alpha * math.sinh(alpha)) - pe) / denom
        is_ = (math.exp(pe) * (pe * math.sinh(alpha) - alpha * math.cosh(alpha)) + alpha) / denom
        num = (1.0 - math.cosh(alpha)) * ic + (eta / alpha + math.sinh(alpha)) * is_
        return pe / (math.exp(pe) - 1.0) * num / w(alpha, eta)

    def delta_quadrature(pe: float, alpha: float, eta: float) -> float:
        return float(np.trapz(kernel(alpha, eta, x) * sigma(pe, x), x))

    def delta0(alpha: float, eta: float) -> float:
        return eta * (math.cosh(alpha) - 1.0) / (alpha * alpha * w(alpha, eta))

    def delta_inf(alpha: float, eta: float) -> float:
        return (math.cosh(alpha) + (eta / alpha) * math.sinh(alpha) - 1.0) / w(alpha, eta)

    def bisect(alpha: float, eta: float, xi: float) -> tuple[float, float, float]:
        lo = xi * delta0(alpha, eta)
        hi = xi * delta_inf(alpha, eta)
        flo = lo - xi * delta_formula(lo, alpha, eta)
        fhi = hi - xi * delta_formula(hi, alpha, eta)
        require(
            "stage058 bracket sign",
            flo <= tol["bracket_tol"] and fhi >= -tol["bracket_tol"],
            f"flo={fmt(flo)}, fhi={fmt(fhi)}",
        )
        for _ in range(bisection_iterations):
            mid = 0.5 * (lo + hi)
            fm = mid - xi * delta_formula(mid, alpha, eta)
            if fm <= 0.0:
                lo = mid
            else:
                hi = mid
        root = 0.5 * (lo + hi)
        residual = root - xi * delta_formula(root, alpha, eta)
        return root, xi * delta0(alpha, eta), xi * delta_inf(alpha, eta), residual

    for case in block["cases"]:
        alpha = float(case["alpha"])
        eta = float(case["eta"])
        xi = float(case["Xi"])
        print(f"\n{case['name']} ({case['kind']}): alpha={fmt(alpha)}, eta={fmt(eta)}, Xi={fmt(xi)}")
        d0_formula = delta0(alpha, eta)
        d0_quad = float(np.trapz(kernel(alpha, eta, x), x))
        dinf = delta_inf(alpha, eta)
        root, pe_lo, pe_hi, residual = bisect(alpha, eta, xi)
        delta_root_formula = delta_formula(root, alpha, eta)
        delta_root_quad = delta_quadrature(root, alpha, eta)
        require(
            f"{case['name']} Delta0 quadrature check",
            near(d0_quad, d0_formula, tol["quadrature_tol"]),
            f"quadrature={fmt(d0_quad)}, formula={fmt(d0_formula)}",
        )
        require(
            f"{case['name']} Delta(root) quadrature check",
            near(delta_root_quad, delta_root_formula, tol["quadrature_tol"]),
            f"quadrature={fmt(delta_root_quad)}, formula={fmt(delta_root_formula)}",
        )
        require(
            f"{case['name']} fixed point residual",
            abs(residual) <= tol["residual_tol"],
            f"residual={fmt(residual)}",
        )
        require(
            f"{case['name']} Pe bracket",
            pe_lo - tol["bracket_tol"] <= root <= pe_hi + tol["bracket_tol"],
            f"Pe_lo={fmt(pe_lo)}, root={fmt(root)}, Pe_hi={fmt(pe_hi)}",
        )
        print(
            f"  Delta0={fmt(d0_formula)}, Delta_inf={fmt(dinf)}, "
            f"Pe_*={fmt(root)}, Xi Delta(Pe_*)={fmt(xi * delta_root_formula)}"
        )


def stage170_block(config: dict) -> None:
    print("\n=== Stage 170: grouped outlet map ===")
    block = config["stage170"]
    tol = block["tolerances"]
    eps = float(block["epsilon_step"])
    baseline = block["baseline"]
    d0 = float(baseline["D0"])
    n0 = float(baseline["N0"])
    sigma = float(baseline["sigma"])
    p0 = n0 / d0
    d2 = -d0 / 9.0
    d4 = -d0 / 27.0

    def finite_derivatives(d_d0: float, d_d2: float, d_d4: float, d_n0: float) -> dict[str, float]:
        u2_full = -(d2 + eps * d_d2) / (d0 + eps * d_d0)
        u4_full = ((d2 + eps * d_d2) ** 2 - (d0 + eps * d_d0) * (d4 + eps * d_d4)) / (d0 + eps * d_d0) ** 2
        p0_full = (n0 + eps * d_n0) / (d0 + eps * d_d0)
        du2 = (u2_full - 1.0 / 9.0) / eps
        du4 = (u4_full - 4.0 / 81.0) / eps
        dp0 = (p0_full - p0) / eps
        dkappa_direct = -3.0 * (1.0 - sigma) * du2 / sigma
        dgamma_direct = -(1.0 - sigma) * dp0 / (9.0 * sigma * p0)
        return {"du2": du2, "du4": du4, "dP0": dp0, "dkappa": dkappa_direct, "dgamma": dgamma_direct}

    def expected_map(d_d0: float, d_d2: float, d_n0: float) -> dict[str, float]:
        k_a = d_d2 + d_d0 / 9.0
        g_a = d_n0 - p0 * d_d0
        dkappa = 3.0 * (1.0 - sigma) * k_a / (sigma * d0)
        dgamma = -(1.0 - sigma) * g_a / (9.0 * sigma * n0)
        return {"K_A": k_a, "G_A": g_a, "dkappa": dkappa, "dgamma": dgamma}

    pair_results: dict[str, dict[str, float]] = {}
    for case in block["consistent_pairs"]:
        d_d0 = float(case["dD0"])
        d_d2 = float(case["dD2"])
        d_n0 = float(case["dN0"])
        d_d4 = 2.0 * d_d2 / 3.0 + d_d0 / 27.0
        exact = finite_derivatives(d_d0, d_d2, d_d4, d_n0)
        expected = expected_map(d_d0, d_d2, d_n0)
        print(f"\n{case['name']}: dD0={fmt(d_d0)}, dD2={fmt(d_d2)}, dD4={fmt(d_d4)}, dN0={fmt(d_n0)}")
        require(
            f"{case['name']} direct even consistency",
            near(exact["du4"], (8.0 / 9.0) * exact["du2"], tol["pair_tol"]),
            f"du4={fmt(exact['du4'])}, 8/9 du2={fmt((8.0 / 9.0) * exact['du2'])}",
        )
        require(
            f"{case['name']} delta kappa map",
            near(exact["dkappa"], expected["dkappa"], tol["linear_tol"]),
            f"direct={fmt(exact['dkappa'])}, expected={fmt(expected['dkappa'])}",
        )
        require(
            f"{case['name']} delta gamma map",
            near(exact["dgamma"], expected["dgamma"], tol["linear_tol"]),
            f"direct={fmt(exact['dgamma'])}, expected={fmt(expected['dgamma'])}",
        )
        pair_results[case["name"]] = {**exact, **expected}
        print(
            f"  K_A={fmt(expected['K_A'])}, G_A={fmt(expected['G_A'])}, "
            f"delta_kappa={fmt(exact['dkappa'])}, delta_gamma={fmt(exact['dgamma'])}"
        )

    first, second = [pair_results[item["name"]] for item in block["consistent_pairs"]]
    require(
        "consistent pair K_A collapse",
        near(first["K_A"], second["K_A"], tol["pair_tol"]),
        f"K_A(1)={fmt(first['K_A'])}, K_A(2)={fmt(second['K_A'])}",
    )
    require(
        "consistent pair G_A collapse",
        near(first["G_A"], second["G_A"], tol["pair_tol"]),
        f"G_A(1)={fmt(first['G_A'])}, G_A(2)={fmt(second['G_A'])}",
    )
    require(
        "consistent pair delta kappa equality",
        near(first["dkappa"], second["dkappa"], tol["pair_tol"]),
        f"delta_kappa(1)={fmt(first['dkappa'])}, delta_kappa(2)={fmt(second['dkappa'])}",
    )
    require(
        "consistent pair delta gamma equality",
        near(first["dgamma"], second["dgamma"], tol["pair_tol"]),
        f"delta_gamma(1)={fmt(first['dgamma'])}, delta_gamma(2)={fmt(second['dgamma'])}",
    )

    for case in block["inconsistent_cases"]:
        d_d0 = float(case["dD0"])
        d_d2 = float(case["dD2"])
        d_d4 = float(case["dD4"])
        d_n0 = float(case["dN0"])
        exact = finite_derivatives(d_d0, d_d2, d_d4, d_n0)
        mismatch = exact["du4"] - (8.0 / 9.0) * exact["du2"]
        print(f"\n{case['name']}: mismatch={fmt(mismatch)}")
        require(
            f"{case['name']} violates hidden-even consistency",
            abs(mismatch) >= tol["inconsistency_floor"],
            f"du4 - 8/9 du2={fmt(mismatch)}",
        )


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[2]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/numerical/stage058_170_foundational_samples.json"
    config = load_config(config_path)
    print(f"Loaded config from {config_path}")
    stage058_block(config)
    stage170_block(config)
    print("\nAll stage-058/170 foundational stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
