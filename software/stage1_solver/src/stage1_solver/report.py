"""Markdown report generation for Step 1 validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6e}"
    return str(value)


def _table(headers: list[str], rows: list[dict[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(h)) for h in headers) + " |")
    return "\n".join(lines)


def write_report(results: dict[str, Any], path: str) -> Path:
    if "step1" in results and "mms" in results:
        return _write_step2_report(results, path)
    return _write_step1_report(results, path)


def _write_step1_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    config = results["config"]
    linear = results["linear"]
    cubic = results["cubic"]

    lines: list[str] = []
    lines.append("# Step 1 GPE Benchmark Validation")
    lines.append("")
    lines.append(f"Overall gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append("")
    lines.append("## Solver Controls")
    lines.append("")
    lines.append("```yaml")
    lines.append(f"dtype: {config['backend']['dtype']}")
    lines.append(f"device: {config['backend']['device']}")
    lines.append(f"seed: {config['backend']['seed']}")
    lines.append(f"deterministic_algorithms: {config['backend']['deterministic_algorithms']}")
    for key, value in config["newton"].items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Discretization")
    lines.append("")
    lines.append(
        "Cell-centered finite volume on radial shells. Radial divergence is "
        "`(F_{i+1/2}-F_{i-1/2})/V_i` with face flux "
        "`F=4*pi*r_face^2*d_r psi`; the inner face has zero area, enforcing "
        "the regular r=0 flux condition. The reusable tensor-product `(r,w)` "
        "operator extends this with cell volume `shell_volume*dw`, radial flux "
        "`4*pi*r_face^2*dw*d_r psi`, and w-flux `shell_volume*d_w psi`. "
        "Dirichlet and Robin boundaries use the same ghost-cell boundary operator."
    )
    lines.append("")
    lines.append("## Linear Benchmark")
    lines.append("")
    lines.append(linear["description"])
    lines.append(f"Reference: {linear['reference']}.")
    lines.append("")
    lines.append(
        _table(
            [
                "nr",
                "dr",
                "computed_mu",
                "eigenvalue_error",
                "eigenvalue_order",
                "field_l2_error",
                "field_l2_order",
                "discrete_eigen_residual_linf",
                "mass_drift",
                "origin_flux_abs",
                "current_max_abs",
            ],
            linear["rows"],
        )
    )
    lines.append("")
    lines.append("Linear checks:")
    for key, value in linear["pass_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("## Cubic GPE Benchmark")
    lines.append("")
    lines.append(cubic["description"])
    lines.append(f"Reference: {cubic['reference']}.")
    ref = cubic["reference_details"]
    lines.append(
        "Reference details: "
        f"mu={ref['mu']:.12e}, mass={ref['mass']:.12e}, "
        f"nodes={ref['solver_nodes']}, status={ref['solver_status']}."
    )
    lines.append("")
    lines.append(
        _table(
            [
                "nr",
                "dr",
                "computed_mu",
                "mu_error",
                "mu_order",
                "field_l2_error",
                "field_l2_order",
                "newton_iterations",
                "newton_final_residual_linf",
                "mass_drift",
                "origin_flux_abs",
                "current_max_abs",
            ],
            cubic["rows"],
        )
    )
    lines.append("")
    jac = cubic["jacobian_check"]
    lines.append(
        "Jacobian JVP cross-check: "
        f"relative={jac['relative_residual']:.6e}, "
        f"absolute={jac['absolute_residual']:.6e}, "
        f"epsilon={jac['epsilon']:.6e}."
    )
    lines.append("")
    lines.append("Cubic checks:")
    for key, value in cubic["pass_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("## Manifests")
    lines.append("")
    for section in (linear, cubic):
        for row in section["rows"]:
            lines.append(f"- {section['name']} nr={row['nr']}: `{row['manifest']}`")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def _append_checks(lines: list[str], checks: dict[str, bool]) -> None:
    for key, value in checks.items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")


def _write_step2_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    config = results["config"]
    step1 = results["step1"]
    mms = results["mms"]

    lines: list[str] = []
    lines.append("# Step 2 Manufactured-Solution Validation")
    lines.append("")
    lines.append(f"Overall gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append("")
    lines.append("## Config Hash Hygiene")
    lines.append("")
    lines.append(
        "`HarnessConfig.config_hash()` excludes `run_root` and `report_path`; manifests still record the full config."
    )
    lines.append("")
    lines.append("## Step 1 Regression")
    lines.append("")
    lines.append(f"Linear regression: {'PASS' if step1['linear']['passed'] else 'FAIL'}.")
    lines.append(f"Cubic GPE regression: {'PASS' if step1['cubic']['passed'] else 'FAIL'}.")
    cubic_jac = step1["cubic"]["jacobian_check"]
    lines.append(
        "Cubic JVP check: "
        f"relative={cubic_jac['relative_residual']:.6e}, "
        f"absolute={cubic_jac['absolute_residual']:.6e}."
    )
    lines.append("")
    lines.append("## MMS Solver Controls")
    lines.append("")
    lines.append("```yaml")
    lines.append(f"dtype: {config['backend']['dtype']}")
    lines.append(f"device: {config['backend']['device']}")
    lines.append(f"seed: {config['backend']['seed']}")
    for key, value in config["newton"].items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Operator Forms")
    lines.append("")
    lines.append(
        "Matter: `-(hbar^2/2m) r^-2 d_r(r^2 d_r psi) + V psi + (5K/4)|psi|^8 psi - mu psi`."
    )
    lines.append(
        "Tensor Laplacian: cell FV divergence with radial flux `4*pi*r_face^2*dw*d_r u` and w flux `shell_volume*d_w u`."
    )
    lines.append("Current: `(hbar/m) Im(conj(psi) d_r psi)` with `A_r=0`.")
    lines.append(
        "Maxwell: stationary axisymmetric `d_M(Z F^{MN}) + xi^-1 d^N(Z d.A)` with components `A0, Ar, Aw` and `H=Z`."
    )
    lines.append(
        "Wall: flat-`dw` stationary modal `-d_w(T_w d_w eta) + [K_eta + ell(ell+1)T_Omega] eta`."
    )
    lines.append("")

    for key in ("matter", "tensor", "current", "maxwell", "wall"):
        section = mms["sections"][key]
        lines.append(f"## {section['name']}")
        lines.append("")
        lines.append(section["description"])
        lines.append(f"Source: {section['continuum_source']}")
        lines.append(f"Manufactured field: `{section['manufactured_field']}`")
        lines.append(f"Forcing: {section['forcing_derivation']}")
        if "methodology_note" in section:
            lines.append(section["methodology_note"])
        if "placeholder_coefficients" in section:
            placeholders = section["placeholder_coefficients"]
            lines.append(f"Placeholder caveat: {placeholders['label']}")
            lines.append("Placeholder coefficients:")
            for coeff_key, coeff_value in placeholders.items():
                if coeff_key != "label":
                    lines.append(f"- {coeff_key}: `{coeff_value}`")
        lines.append("")
        lines.append(
            _table(
                ["grid", "spacing", "error", "observed_order", "reference_norm"],
                section["rows"],
            )
        )
        lines.append("")
        final_row = section["rows"][-1]
        lines.append(
            "Finest-grid error: "
            f"{final_row['error']:.6e} at `{final_row['grid']}`."
        )
        lines.append("")
        if key == "matter":
            jac = section["jacobian_check"]
            lines.append(
                "Matter JVP check: "
                f"relative={jac['relative_residual']:.6e}, "
                f"absolute={jac['absolute_residual']:.6e}, "
                f"epsilon={jac['epsilon']:.6e}."
            )
            lines.append("")
        lines.append("Checks:")
        _append_checks(lines, section["pass_checks"])
        lines.append("")
    lines.append("## Manifests")
    lines.append("")
    for section in mms["sections"].values():
        if "manifest" in section:
            lines.append(f"- {section['name']}: `{section['manifest']}`")
        for row in section.get("rows", []):
            lines.append(f"- {section['name']} {row['grid']}: `{row['manifest']}`")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
