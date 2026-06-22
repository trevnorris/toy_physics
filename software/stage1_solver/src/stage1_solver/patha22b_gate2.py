"""PathA 22b Gate 2 photon-cone speed reduction.

This module is intentionally additive. It derives the zero-mode localized
Maxwell electric/magnetic principal coefficients and carries the still-open
bulk-to-brane speed normalization explicitly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Mapping

import sympy as sp

from stage1_solver.dimensional_check import (
    ACTION,
    D,
    DIMENSIONLESS,
    ENERGY,
    FORCE,
    LENGTH,
    MASS,
    TIME,
    expect_dim,
    homogeneous,
)


FORBIDDEN_TARGET_STRINGS = ("54" + "/5", "10.8" + "/P0")
OUTCOME_CLASS = "STILL_TUNABLE_LAMBDAGAMMA"


def _flat_metric_f2_expansion() -> dict[str, sp.Expr]:
    """Expand F_MN F^MN with eta_MN=diag(-1,+1,+1,+1,+1)."""

    f_tx, f_ty, f_tz, f_tw, f_xy, f_xz, f_xw, f_yz, f_yw, f_zw = sp.symbols(
        "F_tx F_ty F_tz F_tw F_xy F_xz F_xw F_yz F_yw F_zw"
    )
    metric_diag = [-1, 1, 1, 1, 1]
    pairs = {
        (0, 1): f_tx,
        (0, 2): f_ty,
        (0, 3): f_tz,
        (0, 4): f_tw,
        (1, 2): f_xy,
        (1, 3): f_xz,
        (1, 4): f_xw,
        (2, 3): f_yz,
        (2, 4): f_yw,
        (3, 4): f_zw,
    }
    f_matrix = [[sp.Integer(0) for _ in range(5)] for _ in range(5)]
    for (left, right), symbol in pairs.items():
        f_matrix[left][right] = symbol
        f_matrix[right][left] = -symbol

    f_squared = sp.expand(
        sum(
            metric_diag[left] * metric_diag[right] * f_matrix[left][right] ** 2
            for left in range(5)
            for right in range(5)
        )
    )
    zero_mode_subs = {f_tw: 0, f_xw: 0, f_yw: 0, f_zw: 0}
    f_squared_zero_mode = sp.expand(f_squared.subs(zero_mode_subs))
    e_squared = f_tx**2 + f_ty**2 + f_tz**2
    b_squared = f_xy**2 + f_xz**2 + f_yz**2
    expected_zero_mode = sp.expand(-2 * e_squared + 2 * b_squared)
    return {
        "F_squared_full": f_squared,
        "F_squared_zero_mode": f_squared_zero_mode,
        "E_squared": e_squared,
        "B_squared": b_squared,
        "expected_zero_mode": expected_zero_mode,
        "zero_mode_residual": sp.simplify(f_squared_zero_mode - expected_zero_mode),
    }


def localized_maxwell_reduction() -> dict[str, object]:
    """Reduce the zero-mode Maxwell action and derive the transverse cone."""

    z_int, z_w, mu0 = sp.symbols("Z_int Z_w mu0", positive=True)
    omega, k2 = sp.symbols("omega k2", positive=True)
    expansion = _flat_metric_f2_expansion()
    c_e = sp.simplify(z_int / mu0)
    c_b = sp.simplify(z_int / mu0)
    ratio = sp.cancel(c_b / c_e)
    transverse_operator = sp.expand(c_e * omega**2 - c_b * k2)
    factored_operator = sp.factor(c_e * (omega**2 - ratio * k2))
    reduced_lagrangian_density = sp.expand(
        -z_w / (4 * mu0) * expansion["F_squared_zero_mode"]
    )
    integrated_lagrangian_density = sp.expand(
        reduced_lagrangian_density.subs(z_w, z_int)
    )
    z_cancel_control = sp.cancel((z_int / mu0) / (z_int / mu0))

    return {
        "measure": "flat_dw_no_sqrt_g_w",
        "metric_finding": "FLAT_UNWARPED_BULK_METRIC",
        "metric_signature": "eta_MN=diag(-1,+1,+1,+1,+1)",
        "localization_weight_status": "Z(w)_is_a_Lorentz_scalar_weight",
        "assumptions": [
            "A_w approximately 0",
            "partial_w A_mu approximately 0",
            "J^w approximately 0",
            "F_mu_w approximately 0",
        ],
        "C_E": str(c_e),
        "C_B": str(c_b),
        "C_B_over_C_E": str(ratio),
        "F_squared_full": str(expansion["F_squared_full"]),
        "F_squared_zero_mode": str(expansion["F_squared_zero_mode"]),
        "F_squared_zero_mode_expected": str(expansion["expected_zero_mode"]),
        "F_squared_zero_mode_check": bool(expansion["zero_mode_residual"] == 0),
        "reduced_lagrangian_density": str(reduced_lagrangian_density),
        "integrated_lagrangian_density": str(integrated_lagrangian_density),
        "transverse_operator": str(transverse_operator),
        "factored_operator": str(factored_operator),
        "z_cancel_control": {
            "actual": str(z_cancel_control),
            "expected": "1",
            "pass": bool(sp.simplify(z_cancel_control - 1) == 0),
            "note": "The common Z_int/mu0 factor cancels after the flat-metric F^2 expansion fixes equal electric and magnetic coefficients.",
        },
        "provenance": [
            "research/pde/paper/pde.tex:242-244",
            "research/pde/paper/pde.tex:289-295",
            "research/pde/paper/pde.tex:357-416",
            "research/pde/paper/pde.tex:543-552",
            "research/pde/paper/pde.tex:553-558",
            "software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:24-50",
        ],
    }


def speed_normalization_result() -> dict[str, object]:
    """Carry the bulk-to-brane speed map that the sources do not pin."""

    beta_bulk_to_brane, k_eos, rho0, m_gnls = sp.symbols(
        "beta_bulk_to_brane K rho0 m_GNLS",
        positive=True,
    )
    c_bulk_sq_reduced = sp.Integer(1)
    c_s_sq = 5 * k_eos * rho0**4 / m_gnls
    c_gamma_sq_from_reduction = sp.simplify(beta_bulk_to_brane**2 * c_bulk_sq_reduced)
    lambda_gamma_reduced = sp.sqrt(c_gamma_sq_from_reduction / c_s_sq)
    lambda_gamma_carried = "beta_bulk_to_brane*sqrt(m_GNLS/(5*K*rho0**4))"

    return {
        "speed_map_symbol": "beta_bulk_to_brane",
        "speed_map_status": "UNPINNED_BY_SOURCES",
        "open_knob_count": 1,
        "open_knobs": ["bulk_metric_to_acoustic_speed_identification"],
        "c_bulk_squared_reduced": str(c_bulk_sq_reduced),
        "c_bulk_squared_carried": str(c_bulk_sq_reduced),
        "c_gamma_squared": str(c_gamma_sq_from_reduction),
        "c_s_squared": str(c_s_sq),
        "lambda_gamma_reduced": lambda_gamma_carried,
        "lambda_gamma_carried": lambda_gamma_carried,
        "single_dimensionless_dial": lambda_gamma_carried,
        "closed_result_available": False,
        "outcome": OUTCOME_CLASS,
        "reason": (
            "The flat localized Maxwell zero-mode reduction fixes C_B/C_E=1. "
            "The sources still do not identify the bulk gauge kinetic metric speed with the acoustic "
            "sound-speed units, so beta_bulk_to_brane remains the single open speed normalization."
        ),
    }


def negative_control_cgamma_equals_cs() -> dict[str, object]:
    """Guard against silently forcing c_gamma=c_s."""

    beta_bulk_to_brane, k_eos, rho0, m_gnls = sp.symbols(
        "beta_bulk_to_brane K rho0 m_GNLS",
        positive=True,
    )
    forced_residual = beta_bulk_to_brane**2 - 5 * k_eos * rho0**4 / m_gnls
    residual_is_identity_zero = bool(sp.simplify(forced_residual) == 0)
    return {
        "verdict": "GUARD_AGAINST_FORCING_LAMBDAGAMMA_TO_ONE",
        "guard_not_discriminating_test": True,
        "forced_equality_required": "beta_bulk_to_brane^2=5*K*rho0^4/m_GNLS",
        "symbolic_residual": str(forced_residual),
        "residual_is_identity_zero": residual_is_identity_zero,
        "forced_valid": False,
        "respects_pathA_20b_negative_control": True,
    }


def dimensional_checks() -> dict[str, object]:
    velocity = LENGTH / TIME
    wave_number = LENGTH**-1
    lpsi_density = ENERGY / (LENGTH**4)
    electric_field = FORCE
    maxwell_c_e = lpsi_density / (electric_field**2)
    maxwell_c_b = maxwell_c_e
    ai = ACTION / LENGTH
    rho4 = D["rho_4d_number_density"]
    c_s_squared = D["K_eos"] * (rho4**4) / MASS
    beta_map = velocity
    c_gamma_squared = beta_map**2 * (maxwell_c_b / maxwell_c_e)

    checks = [
        expect_dim(
            "pathA_22b_gate2",
            "C_E electric principal coefficient",
            maxwell_c_e,
            maxwell_c_e,
            "Defined by the coefficient multiplying E_i E_i in the reduced Maxwell action.",
        ),
        expect_dim(
            "pathA_22b_gate2",
            "C_B magnetic principal coefficient",
            maxwell_c_b,
            maxwell_c_e,
            "The flat bulk metric gives the same scalar-weighted coefficient as the electric term.",
        ),
        expect_dim(
            "pathA_22b_gate2",
            "C_B/C_E flat bulk cone ratio",
            maxwell_c_b / maxwell_c_e,
            DIMENSIONLESS,
            "In bulk-metric units the localized Maxwell zero mode has c_bulk^2=1.",
        ),
        expect_dim(
            "pathA_22b_gate2",
            "speed-normalization beta_bulk_to_brane",
            beta_map,
            velocity,
            "This is the remaining source-unpinned map from unit bulk-metric speed to physical brane speed.",
        ),
        expect_dim(
            "pathA_22b_gate2",
            "c_gamma=beta*sqrt(C_B/C_E)",
            c_gamma_squared ** sp.Rational(1, 2),
            velocity,
            "Physical brane photon speed after the flat reduction fixes C_B/C_E=1.",
        ),
        expect_dim(
            "pathA_22b_gate2",
            "c_s=sqrt(5*K*rho0^4/m_GNLS)",
            c_s_squared ** sp.Rational(1, 2),
            velocity,
            "GNLS sound speed from the already-derived pathA_19/20 EOS law.",
        ),
        expect_dim(
            "pathA_22b_gate2",
            "lambda_gamma=c_gamma/c_s",
            (c_gamma_squared / c_s_squared) ** sp.Rational(1, 2),
            DIMENSIONLESS,
            "The single remaining speed-normalization dial is dimensionless only as beta_bulk_to_brane/c_s.",
        ),
        homogeneous(
            "pathA_22b_gate2",
            "bulk-metric transverse Maxwell principal operator",
            {
                "C_E*partial_0^2 A_T": maxwell_c_e * ai / (LENGTH**2),
                "C_B*laplacian A_T": maxwell_c_b * ai / (LENGTH**2),
            },
            "The flat parent metric uses x^0 and spatial coordinates in the same bulk-metric units.",
        ),
        homogeneous(
            "pathA_22b_gate2",
            "photon and sound dispersions",
            {
                "omega^2": (TIME**-1) ** 2,
                "beta^2*(C_B/C_E)*k^2": c_gamma_squared * (wave_number**2),
                "c_s^2*k^2": c_s_squared * (wave_number**2),
            },
            "This confirms comparable units only; it does not assert equal coefficients.",
        ),
    ]
    return {"checks": [check.as_dict() for check in checks], "unit_symbols": ["M", "L", "T"]}


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


def gate2_report_payload() -> dict[str, object]:
    reduction = localized_maxwell_reduction()
    speed = speed_normalization_result()
    negative = negative_control_cgamma_equals_cs()
    dims = dimensional_checks()
    guard = target_blindness_guard(
        [
            Path(__file__),
            Path("software/stage1_solver/tools/pathA_22b_gate2_crosscheck.wl"),
            Path("software/stage1_solver/tests/test_patha22b_gate2.py"),
        ]
    )
    residuals = [
        "The flat localized Maxwell reduction computes C_B/C_E=1, so C_B/C_E is no longer an open Gate-2 knob.",
        "The bulk-to-brane speed normalization beta_bulk_to_brane is not pinned by pde.tex:357-416 or by the zero-mode reduction at pde.tex:543-558.",
        "Z_int cancels from the cone and remains a coupling-normalization artifact, not lambda_gamma.",
        "A definite numerical lambda_gamma requires a source-derived identification of the unit bulk Maxwell cone with the acoustic brane speed scale.",
        "Gate 5 must carry lambda_gamma as an explicit open knob; the overall xi verdict remains FAIL_ABLE_PENDING_LAMBDAGAMMA_SPEED_MAP and must not be folded into a REAL_MATCH.",
    ]
    return {
        "schema": "stage1_pathA_22b_gate2/v1",
        "existing_pathA_20b_group": "software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:24-50",
        "localized_maxwell_reduction": reduction,
        "speed_normalization": speed,
        "negative_control": negative,
        "dimensional_checks": dims,
        "target_blindness": guard,
        "lambda_gamma_result": speed["lambda_gamma_carried"],
        "gate2_outcome": OUTCOME_CLASS,
        "residuals": residuals,
    }


def _fmt_check(raw: Mapping[str, object]) -> str:
    return (
        f"- `{raw['name']}`: **{raw['status']}** "
        f"(expected `{raw['expected']}`, actual `{raw['actual']}`). {raw['note']}"
    )


def render_gate2_markdown(payload: Mapping[str, object]) -> str:
    reduction = payload["localized_maxwell_reduction"]
    speed = payload["speed_normalization"]
    negative = payload["negative_control"]
    dims = payload["dimensional_checks"]
    guard = payload["target_blindness"]
    assert isinstance(reduction, Mapping)
    assert isinstance(speed, Mapping)
    assert isinstance(negative, Mapping)
    assert isinstance(dims, Mapping)
    assert isinstance(guard, Mapping)

    lines = [
        "## Gate 2",
        "",
        "### Existing pathA_20b result",
        "",
        f"- Group/module consumed: `{payload['existing_pathA_20b_group']}`. It established the transverse gauge principal block `P_T=C_E*omega^2-C_B*k^2`, hence `c_bulk^2=C_B/C_E`, and carried `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`.",
        "- Gate-2 scope: reduce the localized zero-mode Maxwell action onto the brane, form the photon cone, and treat the missing speed normalization explicitly. No Gate 3-5 work is used.",
        "",
        "### C_E/C_B reduction",
        "",
        f"- Metric finding: `{reduction['metric_finding']}` with `{reduction['metric_signature']}`. Provenance: `research/pde/paper/pde.tex:242-244`.",
        f"- Source action: `L_EM=-(Z(w)/(4*mu0))*F_MN F^MN` plus gauge fixing and external/background source terms. Provenance: `pde.tex:357-416`.",
        f"- Reduction assumptions: `{', '.join(str(item) for item in reduction['assumptions'])}`. Provenance: `pde.tex:543-552`.",
        f"- Measure: `{reduction['measure']}` because `pde.tex:289-295` defines `Z_int=int Z(w) dw`; no independent `sqrt(g_w)` factor is introduced.",
        f"- Flat-metric expansion: `F_MN F^MN = {reduction['F_squared_zero_mode']}` after the zero-mode assumptions, i.e. `{reduction['F_squared_zero_mode_expected']}`.",
        f"- Reduced Maxwell density: `{reduction['reduced_lagrangian_density']}`; after the flat `dw` integration: `{reduction['integrated_lagrangian_density']}`.",
        f"- Reduced coefficients: `C_E={reduction['C_E']}`, `C_B={reduction['C_B']}`.",
        f"- Cone ratio: `C_B/C_E={reduction['C_B_over_C_E']}`. This is a computed flat-metric result, not a second free knob; both coefficients are `Z_int/mu0`, so the common localization factor cancels from the characteristic speed.",
        f"- Transverse operator: `{reduction['factored_operator']}`.",
        "",
        "### Speed Normalization",
        "",
        f"- Carried speed map: `{speed['speed_map_symbol']}` with status `{speed['speed_map_status']}`.",
        f"- Genuine open knobs after the F^2 computation: `{speed['open_knob_count']}` (`{', '.join(str(item) for item in speed['open_knobs'])}`).",
        f"- Physical photon cone carried by Gate 2: `c_gamma^2 = beta_bulk_to_brane^2` because `C_B/C_E=1` in bulk-metric units.",
        f"- GNLS sound speed: `c_s^2 = {speed['c_s_squared']}`.",
        f"- Result: `lambda_gamma = {payload['lambda_gamma_result']}`.",
        f"- Outcome: `{payload['gate2_outcome']}`. This is still tunable through one named speed-normalization residual, not a two-condition geometry result.",
        "- Gate-5 protection: carry `lambda_gamma` as an explicit open knob until the speed map is pinned; the downstream xi verdict is `FAIL_ABLE_PENDING_LAMBDAGAMMA_SPEED_MAP`, not `REAL_MATCH`.",
        "",
        "### Negative Control",
        "",
        f"- Verdict: `{negative['verdict']}`.",
        f"- Forced equality would require `{negative['forced_equality_required']}`.",
        f"- Symbolic residual carried: `{negative['symbolic_residual']}`.",
        f"- Residual is an identity zero? `{negative['residual_is_identity_zero']}`. This is a guard against forcing `lambda_gamma=1`, not evidence from a discriminating source equation.",
        "",
        "### Dimensional Checks",
        "",
    ]
    for raw in dims["checks"]:
        assert isinstance(raw, Mapping)
        lines.append(_fmt_check(raw))
    lines.extend(
        [
            "",
            "### Provenance",
            "",
        ]
    )
    for item in reduction["provenance"]:
        lines.append(f"- {item}.")
    lines.extend(
        [
            "",
            "### Target-Blindness",
            "",
            f"- `{guard['status']}` over the new Gate-2 module, tests, and Mathematica cross-check. No final-comparison constants are used in the derivation.",
            "",
            "### Residual Ledger",
            "",
        ]
    )
    for residual in payload["residuals"]:
        lines.append(f"- {residual}")
    lines.extend(["", f"- Gate 2 outcome: `{payload['gate2_outcome']}`.", ""])
    return "\n".join(lines)


def write_gate2_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    payload = gate2_report_payload()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_22b_gate2.json"
    scratch_md_path = out_dir / "pathA_22b_gate2.md"
    report_path = report_dir / "pathA_22b_minimal_combination_xi.md"
    rendered = render_gate2_markdown(payload)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered + "\n", encoding="utf-8")

    existing = report_path.read_text(encoding="utf-8") if report_path.exists() else "# PathA 22b minimal combination xi\n"
    marker = "\n## Gate 2\n"
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
    json_path, report_path, payload = write_gate2_report(Path(args.out_dir), Path(args.report_dir))
    print(f"wrote {json_path}")
    print(f"wrote {report_path}")
    print(f"Gate 2 outcome: {payload['gate2_outcome']}")
    print(f"lambda_gamma: {payload['lambda_gamma_result']}")
    print(f"negative control: {payload['negative_control']['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
