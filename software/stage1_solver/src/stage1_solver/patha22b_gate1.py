"""PathA 22b Gate 1 Z_int quadrature.

This module evaluates only the exported Maxwell localization integral from the
frozen M1c background.  It does not solve or update any branch/profile data.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from stage1_solver.dimensional_check import D, DIMENSIONLESS, expect_dim


FORBIDDEN_TARGET_STRINGS = ("54" + "/5", "10.8" + "/P0")
FREEZE_HASH = "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8"
FROZEN_DIR = Path("software/stage1_solver/frozen/m1c") / FREEZE_HASH
PRIMARY_BACKGROUND = FROZEN_DIR / "wp1_background_10x8.json"
SOLVE_MANIFEST = FROZEN_DIR / "wp1_solve_manifest.json"


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _as_float_array(values: object, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def exported_z_quadrature(background_path: Path = PRIMARY_BACKGROUND) -> dict[str, object]:
    """Integrate the exported cell-centered Z_w over the exported w grid."""

    packet = _load_json(background_path)
    grid = packet["grid"]
    derived = packet["derived"]
    config = packet.get("config", {})
    if not isinstance(grid, Mapping) or not isinstance(derived, Mapping):
        raise ValueError("background packet must carry grid and derived mappings")

    w_centers = _as_float_array(grid["w_centers"], "grid.w_centers")
    w_faces = _as_float_array(grid["w_faces"], "grid.w_faces")
    z_w = _as_float_array(derived["Z_w"], "derived.Z_w")
    if w_faces.size != w_centers.size + 1:
        raise ValueError("w_faces must have one more entry than w_centers")
    if z_w.size != w_centers.size:
        raise ValueError("Z_w length must match w_centers length")
    if not np.all(np.diff(w_faces) > 0.0):
        raise ValueError("w_faces must be strictly increasing")

    widths = np.diff(w_faces)
    z_int = float(np.dot(z_w, widths))
    left_edge = float(z_w[0])
    right_edge = float(z_w[-1])
    edge_max = float(max(abs(left_edge), abs(right_edge)))
    peak_abs = float(np.max(np.abs(z_w)))
    edge_cell_sum = float(abs(z_w[0]) * widths[0] + abs(z_w[-1]) * widths[-1])
    edge_cell_fraction = float(edge_cell_sum / abs(z_int)) if z_int else math.inf

    branch = {}
    if isinstance(config, Mapping) and isinstance(config.get("branch"), Mapping):
        branch = dict(config["branch"])
    localization_floor = branch.get("localization_floor")
    localization_amplitude = branch.get("localization_amplitude")
    localization_width = branch.get("localization_width")
    w_min = float(grid["w_min"])
    w_max = float(grid["w_max"])
    domain_width = w_max - w_min
    floor_integral_finite_domain = None
    gaussian_integral_finite_domain = None
    floor_fraction_finite_domain = None
    gaussian_fraction_finite_domain = None
    if isinstance(localization_floor, (int, float)):
        floor = float(localization_floor)
        floor_integral_finite_domain = float(floor * domain_width)
        gaussian_integral_finite_domain = float(np.dot(z_w - floor, widths))
        if z_int:
            floor_fraction_finite_domain = float(floor_integral_finite_domain / z_int)
            gaussian_fraction_finite_domain = float(gaussian_integral_finite_domain / z_int)
    gaussian_tail_outside_domain = None
    if (
        isinstance(localization_amplitude, (int, float))
        and isinstance(localization_width, (int, float))
        and float(localization_width) > 0.0
    ):
        midpoint = 0.5 * (w_min + w_max)
        width = float(localization_width)
        amp = float(localization_amplitude)
        left_arg = (midpoint - w_min) / width
        right_arg = (w_max - midpoint) / width
        gaussian_tail_outside_domain = float(
            amp
            * width
            * math.sqrt(math.pi)
            * (1.0 - 0.5 * (math.erf(left_arg) + math.erf(right_arg)))
        )

    return {
        "background_path": str(background_path),
        "schema": packet.get("schema"),
        "grid": {
            "measure": grid.get("measure"),
            "nw": int(grid["nw"]),
            "w_min": w_min,
            "w_max": w_max,
            "domain_width": domain_width,
            "dw": float(grid["dw"]),
            "w_centers": w_centers.tolist(),
            "w_faces": w_faces.tolist(),
        },
        "z_int_finite_domain": z_int,
        "quadrature_rule": "cell-centered midpoint sum over exported w_faces widths",
        "z_edge_left": left_edge,
        "z_edge_right": right_edge,
        "z_edge_max_abs": edge_max,
        "z_peak_abs": peak_abs,
        "edge_to_peak_abs_ratio": float(edge_max / peak_abs) if peak_abs else math.inf,
        "edge_cell_integral_sum_abs": edge_cell_sum,
        "edge_cell_integral_fraction_abs": edge_cell_fraction,
        "localization_floor": localization_floor,
        "floor_integral_finite_domain": floor_integral_finite_domain,
        "floor_fraction_finite_domain": floor_fraction_finite_domain,
        "gaussian_integral_finite_domain": gaussian_integral_finite_domain,
        "gaussian_fraction_finite_domain": gaussian_fraction_finite_domain,
        "gaussian_tail_outside_exported_domain_if_source_formula_continues": gaussian_tail_outside_domain,
        "ideal_infinite_integral_status": (
            "BLOCKED_NEEDS_DECAYING_OR_ZERO_FLOOR_EXTENSION"
            if isinstance(localization_floor, (int, float)) and float(localization_floor) != 0.0
            else "TAIL_DIAGNOSTIC_ONLY"
        ),
    }


def resolution_diagnostics(manifest_path: Path = SOLVE_MANIFEST) -> dict[str, object]:
    manifest = _load_json(manifest_path)
    rows: list[dict[str, object]] = []
    for raw_path in manifest["background_paths"]:
        path = Path(str(raw_path))
        quad = exported_z_quadrature(path)
        grid = quad["grid"]
        assert isinstance(grid, Mapping)
        rows.append(
            {
                "path": str(path),
                "nw": int(grid["nw"]),
                "dw": float(grid["dw"]),
                "z_int_finite_domain": float(quad["z_int_finite_domain"]),
            }
        )
    rows.sort(key=lambda row: int(row["nw"]))
    values = [float(row["z_int_finite_domain"]) for row in rows]
    finest = rows[-1]
    next_finest = rows[-2] if len(rows) >= 2 else None
    nearest_grid_delta = (
        abs(float(finest["z_int_finite_domain"]) - float(next_finest["z_int_finite_domain"]))
        if next_finest is not None
        else 0.0
    )
    return {
        "manifest_path": str(manifest_path),
        "rows": rows,
        "finest_nw": int(finest["nw"]),
        "nearest_grid_delta": float(nearest_grid_delta),
        "ladder_spread": float(max(values) - min(values)) if values else 0.0,
    }


def measure_determination(primary: Mapping[str, object]) -> dict[str, object]:
    grid = primary["grid"]
    assert isinstance(grid, Mapping)
    return {
        "pde_z_int_measure": "flat_dw_no_sqrt_g_w",
        "pde_provenance": "research/pde/paper/pde.tex:289-295",
        "localized_maxwell_provenance": "research/pde/paper/pde.tex:357-416",
        "effective_coupling_provenance": "research/pde/paper/pde.tex:541-563",
        "export_provenance": "software/stage1_solver/src/stage1_solver/m1c_background_export.py:166-170",
        "export_grid_measure_label": grid.get("measure"),
        "sqrt_g_w_exported": False,
        "gate0_structural_measure": "sqrt_g_w_times_Z_dw",
        "discrepancy_flag": (
            "Gate 0's structural carrier included sqrt_g_w, but pde.tex defines Z_int with flat dw. "
            "The frozen export provides no independent sqrt_g_w profile; in the flat/densitized convention "
            "sqrt_g_w=1 and the numeric quadrature is identical."
        ),
        "sqrt_g_w_variant_status": "NOT_INDEPENDENTLY_COMPUTABLE_FROM_EXPORTED_Z_W",
        "sqrt_g_w_variant_if_flat_or_densitized": float(primary["z_int_finite_domain"]),
    }


def dimensional_checks() -> dict[str, object]:
    z_dim = DIMENSIONLESS
    z_int_dim = z_dim * D["dw"]
    return {
        "checks": [
            expect_dim(
                "pathA_22b_gate1_dimensional",
                "exported Z_w localization weight",
                z_dim,
                DIMENSIONLESS,
                "The code exports Z_w from dimensionless floor/amplitude constants and a dimensionless Gaussian argument.",
            ).as_dict(),
            expect_dim(
                "pathA_22b_gate1_dimensional",
                "Z_int=int Z(w) dw",
                z_int_dim,
                D["dw"],
                "pde.tex defines Z_int with flat dw, so the integral carries the w-coordinate length dimension.",
            ).as_dict(),
        ],
        "unit_symbols": ["M", "L", "T"],
    }


def target_blindness_guard(paths: Iterable[Path]) -> dict[str, object]:
    hits: list[str] = []
    for path in paths:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for forbidden in FORBIDDEN_TARGET_STRINGS:
            if forbidden in text:
                hits.append(f"{path}:{forbidden}")
    return {
        "status": "TARGET_BLIND_PASS" if not hits else "TARGET_BLIND_FAILURE",
        "hits": hits,
    }


def gate1_report_payload(background_path: Path = PRIMARY_BACKGROUND) -> dict[str, object]:
    primary = exported_z_quadrature(background_path)
    resolution = resolution_diagnostics()
    measure = measure_determination(primary)
    dims = dimensional_checks()
    guard = target_blindness_guard(
        [
            Path(__file__),
            Path("software/stage1_solver/tools/pathA_22b_gate1_crosscheck.wl"),
        ]
    )
    grid_error = float(resolution["nearest_grid_delta"])
    residuals = [
        "The exported primary grid is only nw=8; the nearest-grid finite-domain change is recorded as a grid-resolution delta only, not as the total uncertainty.",
        "The exported Z_w is not small at the domain edges, so the finite interval is not a controlled small-tail approximation to an ideal infinite integral.",
        "Because the frozen source has a nonzero localization_floor, an infinite continuation of that same formula would not give a finite Z_int; an external decaying/zero-floor continuation would be needed to bound the omitted tail.",
        "Floor provenance: localization_floor=0.8 is an undocumented solver config constant (coupled_branch.py:188-192; m1c_physical_run.py:121-123), differs across presets (patha_closed_newton.py:61-63), and has no source support; pde.tex:277,289-295 uses a localized Z(w) over (-infinity,+infinity). Next step to unblock: export/derive a genuinely decaying Z(w), or obtain documented physical provenance for the floor plus a physical w-extent.",
        "No photon-cone speed or lambda_gamma normalization is derived or modified in Gate 1.",
    ]
    return {
        "schema": "stage1_pathA_22b_gate1/v1",
        "primary_quadrature": primary,
        "resolution_diagnostics": resolution,
        "measure_determination": measure,
        "dimensional_checks": dims,
        "target_blindness": guard,
        "z_int_value": float(primary["z_int_finite_domain"]),
        "z_int_error_bar_grid": grid_error,
        "z_int_units": "L",
        "scope_statement": (
            "Gate 1 reports Z_int only as a coupling-normalization artifact: mu0_eff=mu0/Z_int "
            "and q_eff=q_star/sqrt(Z_int). Z_int is not a factor in P0*chi_Q*g_mhat^2*lambda_gamma^5/g_G, "
            "does not gate the xi verdict, and does not promote Z_int to lambda_gamma or alter photon-cone speed. "
            "If needed downstream, carry it symbolically as mu0_eff=mu0/Z_int and q_eff=q_star/sqrt(Z_int), never as the numeric finite-box value."
        ),
        "gate1_outcome": "BLOCKED_NEEDS_DECAYING_Z_PROFILE_OR_FLOOR_PROVENANCE",
        "residuals": residuals,
    }


def _fmt(value: object, digits: int = 12) -> str:
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def render_gate1_markdown(payload: Mapping[str, object]) -> str:
    primary = payload["primary_quadrature"]
    resolution = payload["resolution_diagnostics"]
    measure = payload["measure_determination"]
    dims = payload["dimensional_checks"]
    guard = payload["target_blindness"]
    assert isinstance(primary, Mapping)
    assert isinstance(resolution, Mapping)
    assert isinstance(measure, Mapping)
    assert isinstance(dims, Mapping)
    assert isinstance(guard, Mapping)
    grid = primary["grid"]
    assert isinstance(grid, Mapping)

    lines = [
        "## Gate 1",
        "",
        "### Measure determination",
        "",
        f"- PDE definition used: `{measure['pde_z_int_measure']}`. `pde.tex:289-295` defines `Z_int=int Z(w) dw`; it does not include a `sqrt(g_w)` factor.",
        f"- Discrepancy flag: {measure['discrepancy_flag']}",
        f"- Numeric `sqrt(g_w)` variant: `{measure['sqrt_g_w_variant_status']}`. If the export is read in the flat/densitized convention, the variant equals `{_fmt(measure['sqrt_g_w_variant_if_flat_or_densitized'])} L`.",
        f"- Exported grid measure label: `{measure['export_grid_measure_label']}`; w quadrature used the exported `w_faces` widths.",
        "",
        "### Quadrature result",
        "",
        f"- Headline: `Z_int ~= {_fmt(payload['z_int_value'])} L` (finite-box, floor-dominated, domain-dependent; ideal +/-infinity integral divergent).",
        f"- Rule: {primary['quadrature_rule']}. Primary background: `{primary['background_path']}`.",
        f"- Domain/grid: `w in [{_fmt(grid['w_min'])}, {_fmt(grid['w_max'])}]`, `nw={grid['nw']}`, `dw={_fmt(grid['dw'])}`.",
        f"- Resolution diagnostic: nearest-grid change `{_fmt(resolution['nearest_grid_delta'])} L` (grid-resolution delta only -- NOT total uncertainty); full 4/6/8-point ladder spread `{_fmt(resolution['ladder_spread'])} L`.",
        f"- Tail/truncation diagnostic: edge values `Z_left={_fmt(primary['z_edge_left'])}`, `Z_right={_fmt(primary['z_edge_right'])}`, max edge/peak ratio `{_fmt(primary['edge_to_peak_abs_ratio'])}`.",
        f"- Edge-cell contribution diagnostic: the two boundary cells contribute `{_fmt(primary['edge_cell_integral_sum_abs'])} L`, fraction `{_fmt(primary['edge_cell_integral_fraction_abs'])}` of the finite-domain integral.",
        f"- Ideal infinite-domain status: `{primary['ideal_infinite_integral_status']}`. The exported nonzero floor means the finite-domain result is the faithful exported-grid integral, but the omitted ideal tail is not bounded as a small correction by this export.",
    ]
    if primary.get("floor_integral_finite_domain") is not None:
        lines.append(
            f"- Floor decomposition: `localization_floor={_fmt(primary['localization_floor'])}` over domain width `{_fmt(grid['domain_width'])} L` contributes `{_fmt(primary['floor_integral_finite_domain'])} L`, "
            f"which is `{_fmt(100.0 * float(primary['floor_fraction_finite_domain']), 3)}%` of finite-box `Z_int`; the localized Gaussian part contributes `{_fmt(primary['gaussian_integral_finite_domain'])} L` "
            f"(`{_fmt(100.0 * float(primary['gaussian_fraction_finite_domain']), 3)}%`), so the Gaussian is the minority."
        )
    gaussian_tail = primary.get("gaussian_tail_outside_exported_domain_if_source_formula_continues")
    if gaussian_tail is not None:
        lines.append(
            f"- Gaussian-only omitted-tail diagnostic, ignoring the nonzero floor extension: `{_fmt(gaussian_tail)} L`."
        )

    lines.extend(
        [
            "",
            "### Dimensions and scope",
            "",
        ]
    )
    for raw in dims["checks"]:
        assert isinstance(raw, Mapping)
        lines.append(
            f"- `{raw['name']}`: **{raw['status']}** (expected `{raw['expected']}`, actual `{raw['actual']}`). {raw['note']}"
        )
    lines.extend(
        [
            f"- Scope: {payload['scope_statement']}",
            "",
            "### Provenance",
            "",
            f"- Measure: {measure['pde_provenance']}.",
            f"- Localized Maxwell source of `Z(w)`: {measure['localized_maxwell_provenance']}.",
            f"- Coupling reduction: {measure['effective_coupling_provenance']}.",
            f"- Exported `Z_w`: {measure['export_provenance']}.",
            f"- Frozen data: `{primary['background_path']}`.",
            "",
            "### Target-blindness",
            "",
            f"- `{guard['status']}` over the new Gate-1 module and Mathematica cross-check.",
            "",
            "### Residual ledger",
            "",
        ]
    )
    for residual in payload["residuals"]:
        lines.append(f"- {residual}")
    lines.extend(["", f"- Gate 1 outcome: `{payload['gate1_outcome']}`.", ""])
    return "\n".join(lines)


def write_gate1_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    payload = gate1_report_payload()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_22b_gate1.json"
    scratch_md_path = out_dir / "pathA_22b_gate1.md"
    report_path = report_dir / "pathA_22b_minimal_combination_xi.md"
    rendered = render_gate1_markdown(payload)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered + "\n", encoding="utf-8")

    existing = report_path.read_text(encoding="utf-8") if report_path.exists() else "# PathA 22b minimal combination xi\n"
    marker = "\n## Gate 1\n"
    if marker in existing:
        existing = existing.split(marker, 1)[0].rstrip() + "\n"
    elif not existing.endswith("\n"):
        existing += "\n"
    report_path.write_text(existing.rstrip() + "\n\n" + rendered + "\n", encoding="utf-8")
    return json_path, report_path, payload


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="software/stage1_solver/_scratch")
    parser.add_argument("--report-dir", default="software/stage1_solver/reports")
    args = parser.parse_args(list(argv) if argv is not None else None)
    json_path, report_path, payload = write_gate1_report(Path(args.out_dir), Path(args.report_dir))
    print(f"wrote {json_path}")
    print(f"wrote {report_path}")
    print(
        "Gate 1 Z_int="
        f"{payload['z_int_value']:.12g} {payload['z_int_units']} "
        f"(grid-resolution delta only {payload['z_int_error_bar_grid']:.12g} {payload['z_int_units']}; not total uncertainty)"
    )
    print(payload["gate1_outcome"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
