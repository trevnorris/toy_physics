#!/usr/bin/env python3
"""Branch-sensitivity stress harness for Stage 089."""

from __future__ import annotations

import cmath
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


def outgoing_probe_z(a_over_r: float) -> float:
    # Probe-only low-frequency sample used to extract the outgoing odd DtN
    # coefficient from the exact spherical branch.
    return min(2.0e-2, max(5.0e-4, 10.0 * a_over_r))


def spherical_hankel_h1_l1(z: float) -> complex:
    return -cmath.exp(1j * z) * (1.0 + 1j / z) / z


def spherical_hankel_h1_l2(z: float) -> complex:
    return 1j * cmath.exp(1j * z) * (1.0 + 3j / z - 3.0 / (z * z)) / z


def outgoing_yhat_exact(z: float) -> complex:
    h1 = spherical_hankel_h1_l1(z)
    h2 = spherical_hankel_h1_l2(z)
    lambda_out = z * (h1 - 3.0 * h2 / z) / h2
    return -3.0 / lambda_out


def canonical_chi_from_outgoing_dtn(z: float) -> float:
    yhat = outgoing_yhat_exact(z)
    even_part = 1.0 + z * z / 9.0 + 4.0 * z**4 / 81.0
    return 27.0 * (yhat - even_part).imag / z**5


def canonical_point_particle_approx(a_over_r: float) -> float:
    return 1.0 - 2.0 * a_over_r * a_over_r


def stage089_block(config: dict) -> None:
    print("=== Stage 089: canonical outgoing closure stress ===")
    print(f"source-map model: {config['source_map_model']}")

    canonical_errors: list[tuple[str, float, float]] = []
    perturbation_pairs: dict[float, dict[str, tuple[float, float, float]]] = {}

    for case in config["canonical_samples"]:
        a_over_r = float(case["a_over_r"])
        mhat0 = source_map_model(a_over_r)
        chi_q = 1.0
        z_probe = outgoing_probe_z(a_over_r)
        chi_exact = canonical_chi_from_outgoing_dtn(z_probe)
        nq = nq_value(mhat0, chi_q)
        pp_approx = canonical_point_particle_approx(a_over_r)
        remainder = nq - pp_approx
        scaled_remainder = remainder / (a_over_r**4)
        coeff_rel_error = abs(nq - 1.0)
        canonical_errors.append((case["name"], a_over_r, coeff_rel_error))

        print(
            f"\n{case['name']} ({case['kind']}): a/r={fmt(a_over_r)}, "
            f"mhat0={fmt(mhat0)}, z_probe={fmt(z_probe)}"
        )
        print(
            f"  N_Q={fmt(nq)}, point-particle approx={fmt(pp_approx)}, "
            f"chi_Q^out(exact)={fmt(chi_exact)}"
        )

        require(
            f"{case['name']} canonical branch keeps N_Q <= 1",
            nq <= 1.0,
            f"N_Q={fmt(nq)}",
        )
        require(
            f"{case['name']} exact outgoing DtN fixes chi_Q = 1",
            abs(chi_exact - 1.0) <= 5.0e-5,
            f"chi_Q^out(exact)={fmt(chi_exact)}",
        )
        require(
            f"{case['name']} outgoing DtN and source-map canonical N_Q agree",
            abs(nq_value(mhat0, chi_exact) - nq) <= 5.0e-5,
            f"N_Q(chi_exact)={fmt(nq_value(mhat0, chi_exact))}, N_Q(chi=1)={fmt(nq)}",
        )
        require(
            f"{case['name']} exact N_Q stays above the O((a/r)^2) point-particle truncation",
            remainder > 0.0,
            f"N_Q - [1 - 2(a/r)^2]={fmt(remainder)}",
        )
        if a_over_r >= 1e-3:
            require(
                f"{case['name']} point-particle remainder is O((a/r)^4)",
                2.9 <= scaled_remainder <= 3.05,
                f"(N_Q - [1 - 2(a/r)^2])/(a/r)^4={fmt(scaled_remainder)}",
            )
        else:
            require(
                f"{case['name']} point-particle remainder stays below floating-point noise",
                remainder <= 1e-14,
                f"N_Q - [1 - 2(a/r)^2]={fmt(remainder)}",
            )

    for (name_small, a_small, err_small), (name_large, a_large, err_large) in zip(
        canonical_errors[:-1], canonical_errors[1:], strict=True
    ):
        require(
            f"{name_small} is closer to target than {name_large}",
            err_small < err_large,
            f"|N_Q-1|({fmt(a_small)})={fmt(err_small)}, |N_Q-1|({fmt(a_large)})={fmt(err_large)}",
        )

    for case in config["branch_perturbations"]:
        a_over_r = float(case["a_over_r"])
        mhat0 = source_map_model(a_over_r)
        chi_q = float(case["chi_Q"])
        z_probe = outgoing_probe_z(a_over_r)
        chi_exact = canonical_chi_from_outgoing_dtn(z_probe)
        nq = nq_value(mhat0, chi_q)
        nq_canonical = nq_value(mhat0, chi_exact)

        print(
            f"\n{case['name']} ({case['kind']}): a/r={fmt(a_over_r)}, "
            f"mhat0={fmt(mhat0)}, chi_Q={fmt(chi_q)}, "
            f"chi_Q^out(exact)={fmt(chi_exact)}"
        )
        print(
            f"  N_Q={fmt(nq)}, canonical N_Q={fmt(nq_canonical)}, "
            f"chi gap={fmt(chi_q - chi_exact)}"
        )

        if chi_q > chi_exact:
            branch_key = "positive"
            require(
                f"{case['name']} positive chi_Q shift lowers N_Q",
                nq < nq_canonical,
                f"N_Q={fmt(nq)}, canonical={fmt(nq_canonical)}",
            )
            require(
                f"{case['name']} positive chi_Q shift stays above the outgoing DtN branch",
                chi_q > chi_exact,
                f"chi_Q={fmt(chi_q)}, chi_Q^out(exact)={fmt(chi_exact)}",
            )
        else:
            branch_key = "negative"
            require(
                f"{case['name']} negative chi_Q shift raises N_Q",
                nq > nq_canonical,
                f"N_Q={fmt(nq)}, canonical={fmt(nq_canonical)}",
            )
            require(
                f"{case['name']} negative chi_Q shift stays below the outgoing DtN branch",
                chi_q < chi_exact,
                f"chi_Q={fmt(chi_q)}, chi_Q^out(exact)={fmt(chi_exact)}",
            )

        perturbation_pairs.setdefault(a_over_r, {})[branch_key] = (chi_q, chi_exact, nq)

    for a_over_r, pair in sorted(perturbation_pairs.items()):
        if set(pair) != {"positive", "negative"}:
            raise AssertionError(f"missing symmetric perturbation pair at a/r={a_over_r}")
        chi_plus, chi_exact_plus, nq_plus = pair["positive"]
        chi_minus, chi_exact_minus, nq_minus = pair["negative"]
        chi_exact = 0.5 * (chi_exact_plus + chi_exact_minus)
        canonical_nq = nq_value(source_map_model(a_over_r), chi_exact)
        slope = (nq_plus - nq_minus) / (chi_plus - chi_minus)
        curvature = nq_plus + nq_minus - 2.0 * canonical_nq

        require(
            f"a/r={fmt(a_over_r)} symmetric branch slope around the outgoing DtN branch is negative",
            slope < 0.0,
            f"slope={fmt(slope)}",
        )
        require(
            f"a/r={fmt(a_over_r)} symmetric branch response is centered on chi_Q^out",
            abs(chi_exact_plus - chi_exact_minus) <= 1.0e-4,
            f"chi_Q^out(+)= {fmt(chi_exact_plus)}, chi_Q^out(-)= {fmt(chi_exact_minus)}",
        )
        require(
            f"a/r={fmt(a_over_r)} reciprocal conservative response is convex",
            curvature > 0.0,
            f"N_Q(chi_+) + N_Q(chi_-) - 2 N_Q(chi_out)={fmt(curvature)}",
        )
        require(
            f"a/r={fmt(a_over_r)} symmetric chi_Q shifts stay on opposite sides of chi_Q^out",
            chi_minus < chi_exact < chi_plus,
            f"chi_-={fmt(chi_minus)}, chi_out={fmt(chi_exact)}, chi_+={fmt(chi_plus)}",
        )
        require(
            f"a/r={fmt(a_over_r)} symmetric chi_Q shifts separate the conservative coefficient",
            nq_minus - nq_plus > 0.0,
            f"N_Q(chi_-)-N_Q(chi_+)={fmt(nq_minus - nq_plus)}",
        )

    print("\nAll stage-089 canonical outgoing stress checks passed.")


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[2]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/numerical/stage089_canonical_outgoing_samples.json"
    config = load_config(config_path)
    print(f"Loaded config from {config_path}")
    stage089_block(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
