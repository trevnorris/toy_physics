#!/usr/bin/env python3
"""Stage 140 numerical capstone stress harness.

Reuses the Stage-138/139 fixed-point grid to certify three facts:
1. Frozen canonical traction does not keep the co-evolving branch compensated.
2. The renormalized canonical point is unique on the analyzed scan window.
3. The capstone status point really is a higher-traction, higher-bias repair.
"""

from __future__ import annotations

import sys
from pathlib import Path

from stage138_139_fixedpoint_stress import (
    FixedPointHarness,
    check_seed_consistency,
    check_stage139,
    fmt,
    load_config,
    near,
    require,
)


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[2]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/numerical/stage138_139_fixedpoint_samples.json"
    config = load_config(config_path)
    h = FixedPointHarness(config)

    print(f"Loaded config from {config_path}")

    frozen_case = next(case for case in config["seed_check_sigmas"] if case["name"] == "frozen_canonical_traction")
    renorm_case = next(case for case in config["seed_check_sigmas"] if case["name"] == "near_renormalized_branch")
    check_seed_consistency(h, frozen_case)
    check_seed_consistency(h, renorm_case)

    print("\n=== stage140_frozen_vs_renormalized ===")
    frozen = h.fixed_point_summary(h.Sigma0_star, h._canonical_seed)
    require(
        "frozen traction loses exact compensation",
        frozen["g"] < h.g_star and frozen["R"] > 0.25,
        f"g={fmt(frozen['g'])}, g_star={fmt(h.g_star)}, R={fmt(frozen['R'])}",
    )
    require(
        "frozen traction is not the renormalized root",
        abs(frozen["g"] - h.g_star) >= h.moment_tol and abs(frozen["R"] - 0.25) >= h.moment_tol,
        f"|g-g_star|={fmt(abs(frozen['g'] - h.g_star))}, |R-1/4|={fmt(abs(frozen['R'] - 0.25))}",
    )

    check_stage139(h, config)

    renorm = h.fixed_point_summary(h.Sigma0_can_expected, h._canonical_seed)
    require(
        "renormalized traction is strictly above the original canonical point",
        h.Sigma0_can_expected > h.Sigma0_star and renorm["Pi"] > h.Pi_star and renorm["T_hat"] > h.T_hat_star,
        (
            f"Sigma0_can={fmt(h.Sigma0_can_expected)}, Sigma0_star={fmt(h.Sigma0_star)}, "
            f"Pi_can={fmt(renorm['Pi'])}, Pi_star={fmt(h.Pi_star)}, "
            f"T_hat_can={fmt(renorm['T_hat'])}, T_hat_star={fmt(h.T_hat_star)}"
        ),
    )
    require(
        "renormalized point matches the carried Stage-139 tuple",
        near(renorm["S"], h.S_can_expected, h.moment_tol)
        and near(renorm["Pi"], h.Pi_can_expected, h.moment_tol)
        and near(renorm["T_hat"], h.T_hat_can_expected, h.moment_tol),
        (
            f"S={fmt(renorm['S'])}, expected={fmt(h.S_can_expected)}; "
            f"Pi={fmt(renorm['Pi'])}, expected={fmt(h.Pi_can_expected)}; "
            f"T_hat={fmt(renorm['T_hat'])}, expected={fmt(h.T_hat_can_expected)}"
        ),
    )

    print("\nStage 140 numerical capstone checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
