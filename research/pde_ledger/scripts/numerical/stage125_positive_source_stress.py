#!/usr/bin/env python3
"""Positive-source numerical stress harness for Stage 125."""

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


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage125_v1":
        raise ValueError("unexpected sample schema")
    return data


def stage125_block(config: dict) -> None:
    print("=== Stage 125: positive local mouth-source stress ===")
    n = int(config["quadrature_points"])
    iters = int(config["bisection_iterations"])
    tol = config["tolerances"]
    g_minus = float(config["g_minus"])
    g_plus = float(config["g_plus"])
    x = np.linspace(0.0, 1.0, n)
    kernel = np.cos(np.pi * x / 2.0)

    def trapz(values: np.ndarray) -> float:
        return float(np.trapz(values, x))

    def normalized_profile(profile: dict) -> np.ndarray:
        kind = profile["kind"]
        if kind == "uniform":
            raw = np.ones_like(x)
        elif kind == "left_linear":
            raw = 2.0 * (1.0 - x)
        elif kind == "left_beta":
            p = float(profile["p"])
            raw = (p + 1.0) * (1.0 - x) ** p
        elif kind == "right_beta":
            p = float(profile["p"])
            raw = (p + 1.0) * x ** p
        elif kind == "exp_left":
            beta = float(profile["beta"])
            raw = beta * np.exp(-beta * x) / (1.0 - math.exp(-beta))
        else:
            raise ValueError(f"unknown profile kind: {kind}")
        return raw / trapz(raw)

    def g_value(profile: np.ndarray) -> float:
        return trapz(profile * kernel)

    values: dict[str, float] = {}
    for profile in config["profiles"]:
        sigma = normalized_profile(profile)
        norm = trapz(sigma)
        g_sigma = g_value(sigma)
        values[profile["name"]] = g_sigma
        print(f"\n{profile['name']} ({profile['kind']}): g[sigma]={fmt(g_sigma)}")
        require(
            f"{profile['name']} stays nonnegative",
            float(np.min(sigma)) >= -1e-12,
            f"min sigma={fmt(float(np.min(sigma)))}",
        )
        require(
            f"{profile['name']} normalization",
            abs(norm - 1.0) <= float(tol["normalization_tol"]),
            f"integral={fmt(norm)}",
        )
        require(
            f"{profile['name']} moment stays in [0,1]",
            -float(tol["moment_tol"]) <= g_sigma <= 1.0 + float(tol["moment_tol"]),
            f"g[sigma]={fmt(g_sigma)}",
        )

    require("g_- lies in (0,1)", 0.0 < g_minus < 1.0, f"g_-={fmt(g_minus)}")
    require("g_+ lies above 1", g_plus > 1.0, f"g_+={fmt(g_plus)}")
    require(
        "uniform < g_- < left_linear",
        values["uniform"] < g_minus < values["left_linear"],
        f"uniform={fmt(values['uniform'])}, g_-={fmt(g_minus)}, left_linear={fmt(values['left_linear'])}",
    )
    require(
        "right-localized profiles drive g toward 0",
        values["right_beta_20"] < values["right_beta_8"] < values["uniform"],
        f"right_beta_20={fmt(values['right_beta_20'])}, right_beta_8={fmt(values['right_beta_8'])}, uniform={fmt(values['uniform'])}",
    )
    require(
        "left-localized profiles drive g toward 1",
        values["left_beta_20"] > values["left_beta_8"] > values["left_linear"] > g_minus,
        f"left_beta_20={fmt(values['left_beta_20'])}, left_beta_8={fmt(values['left_beta_8'])}, left_linear={fmt(values['left_linear'])}, g_-={fmt(g_minus)}",
    )

    bracket = config["gminus_realization"]["beta_bracket"]
    lo = float(bracket[0])
    hi = float(bracket[1])

    def exp_profile(beta: float) -> np.ndarray:
        raw = beta * np.exp(-beta * x) / (1.0 - math.exp(-beta))
        return raw / trapz(raw)

    def f(beta: float) -> float:
        return g_value(exp_profile(beta)) - g_minus

    flo = f(lo)
    fhi = f(hi)
    require(
        "g_- realization bracket",
        flo <= 0.0 and fhi >= 0.0,
        f"f(lo)={fmt(flo)}, f(hi)={fmt(fhi)}",
    )
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm <= 0.0:
            lo = mid
        else:
            hi = mid
    beta_star = 0.5 * (lo + hi)
    g_star = g_value(exp_profile(beta_star))
    require(
        "exp_left family realizes g_-",
        abs(g_star - g_minus) <= float(tol["root_tol"]),
        f"beta*={fmt(beta_star)}, g(beta*)={fmt(g_star)}, g_-={fmt(g_minus)}",
    )

    print(f"\nSolved exp_left beta* = {fmt(beta_star)} with g[beta*] = {fmt(g_star)}")
    print("\nAll stage-125 positive-source stress checks passed.")


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[2]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/numerical/stage125_positive_source_samples.json"
    config = load_config(config_path)
    print(f"Loaded config from {config_path}")
    stage125_block(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
