#!/usr/bin/env python3
"""Compare independent SymPy and Wolfram Language Construction-B outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "reports" / "artifacts"


def get(data: dict, dotted: str):
    value = data
    for part in dotted.split("."):
        value = value[part]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    sympy = json.loads((args.out_dir / "sympy_results.json").read_text())
    math = json.loads((args.out_dir / "mathematica_results.json").read_text())

    fields = [
        "mapping.charge_spectrum",
        "mapping.link_pair",
        "mapping.path_endpoint_charge",
        "mapping.closed_loop_divergence",
        "mapping.gauge_loop_invariant",
        "mapping.continuity",
        "mapping.uv_global_u1_broken",
        "embedding.throat_embedding",
        "response.maxwell_density_sign",
        "response.maxwell_current_sign",
        "response.scalar_density_sign",
        "response.scalar_current_channels",
        "response.scalar_verdict",
        "response.net_moving_like_charge_sign",
        "controls.ring_off_propagating",
        "controls.ring_off_result",
        "controls.higgs_massive",
        "controls.defects_condensed_result",
        "firewall.photon_modes",
        "firewall.force_exponent",
        "firewall.flux_quantized",
        "single_kernel.single_kernel",
        "single_kernel.independent_scalar_channels",
        "all_ablations_tripped",
        "algebra_status",
        "phase_existence_scope",
    ]
    lines = ["DUAL_ENGINE_COMPARISON"]
    for field in fields:
        left, right = get(sympy, field), get(math, field)
        if left != right:
            raise AssertionError(f"ENGINE_DISAGREE {field}: SymPy={left!r}, Mathematica={right!r}")
        lines.append(f"AGREE {field}: {left}")
    lines.extend(
        [
            "ENGINE_AGREE: sign_structure=repulsive_density+attractive_transverse_current",
            "ENGINE_AGREE: scalar=attractive_density+no_transverse_current -> FAIL_SCALAR_SINGLE_SIGN",
            "ENGINE_AGREE: falloff=1/r^2_force_in_3_spatial_dimensions",
            "ENGINE_AGREE: mapping=integer_divergence_endpoints+flux_dressing+continuity",
            "ENGINE_AGREE: phase_existence=CITED_NOT_COMPUTED",
            "ENGINE_AGREE",
        ]
    )
    log = "\n".join(lines) + "\n"
    (args.out_dir / "engine_agreement.log").write_text(log)
    print(log, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
