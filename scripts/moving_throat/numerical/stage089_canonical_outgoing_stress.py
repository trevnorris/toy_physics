#!/usr/bin/env python3
"""Branch-sensitivity stress harness for Stage 089."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def fmt(value: float) -> str:
    return f"{value:.12g}"


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage089_v1":
        raise ValueError("unexpected sample schema")
    return data


def source_map_model(a_over_r: float) -> float:
    return 1.0 + a_over_r * a_over_r


def nq_value(mhat0: float, chi_q: float) -> float:
    return 1.0 / (mhat0 * mhat0 * chi_q)


def stage089_block(config: dict) -> None:
    print("=== Stage 089: canonical outgoing closure stress ===")
    print(f"source-map model: {config['source_map_model']}")

    canonical_errors: list[tuple[str, float, float]] = []
    for case in config["canonical_samples"]:
        a_over_r = float(case["a_over_r"])
        mhat0 = source_map_model(a_over_r)
        chi_q = 1.0
        nq = nq_value(mhat0, chi_q)
        gamma_eff_ratio = mhat0 * mhat0 * nq
        coeff_rel_error = abs(nq - 1.0)
        canonical_errors.append((case["name"], a_over_r, coeff_rel_error))
        print(f"\n{case['name']} ({case['kind']}): a/r={fmt(a_over_r)}, mhat0={fmt(mhat0)}")
        print(f"  N_Q={fmt(nq)}, coefficient ratio={fmt(nq)}, gamma_eff ratio={fmt(gamma_eff_ratio)}")
        require(
            f"{case['name']} canonical gamma_eff hits target",
            abs(gamma_eff_ratio - 1.0) <= 1e-12,
            f"gamma_eff/target={fmt(gamma_eff_ratio)}",
        )
        require(
            f"{case['name']} canonical branch keeps N_Q <= 1",
            nq <= 1.0,
            f"N_Q={fmt(nq)}",
        )

    for (name_small, a_small, err_small), (name_large, a_large, err_large) in zip(canonical_errors[:-1], canonical_errors[1:], strict=True):
        require(
            f"{name_small} is closer to target than {name_large}",
            err_small < err_large,
            f"|N_Q-1|({fmt(a_small)})={fmt(err_small)}, |N_Q-1|({fmt(a_large)})={fmt(err_large)}",
        )

    for case in config["branch_perturbations"]:
        a_over_r = float(case["a_over_r"])
        mhat0 = source_map_model(a_over_r)
        chi_q = float(case["chi_Q"])
        nq = nq_value(mhat0, chi_q)
        gamma_eff_ratio = mhat0 * mhat0 * nq
        separation = abs(gamma_eff_ratio - 1.0)
        print(f"\n{case['name']} ({case['kind']}): a/r={fmt(a_over_r)}, mhat0={fmt(mhat0)}, chi_Q={fmt(chi_q)}")
        print(f"  N_Q={fmt(nq)}, coefficient ratio={fmt(nq)}, gamma_eff ratio={fmt(gamma_eff_ratio)}")
        require(
            f"{case['name']} branch deformation moves gamma_eff off target",
            separation >= float(case["separation_floor"]),
            f"|gamma_eff/target - 1|={fmt(separation)}, floor={fmt(float(case['separation_floor']))}",
        )
        require(
            f"{case['name']} gamma_eff follows 1/chi_Q",
            abs(gamma_eff_ratio - 1.0 / chi_q) <= 1e-12,
            f"gamma_eff/target={fmt(gamma_eff_ratio)}, 1/chi_Q={fmt(1.0 / chi_q)}",
        )

    print("\nAll stage-089 canonical outgoing stress checks passed.")


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[3]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/moving_throat/numerical/stage089_canonical_outgoing_samples.json"
    config = load_config(config_path)
    print(f"Loaded config from {config_path}")
    stage089_block(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
