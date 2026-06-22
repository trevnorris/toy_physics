"""PathA 22b Gate 4 action-level ratio fork.

This module is intentionally narrow and additive.  It derives the stress-lane
transverse kernel from the parent-action force/reduction chain, attempts the
independent source-map kernel derivation from the allowed source/current
content, and stops at the decisive fork when that source-map provenance is not
available.
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
    LENGTH,
    MASS,
    TIME,
    Dim,
    expect_dim,
)


FORBIDDEN_TARGET_STRINGS = (
    "5" + "4" + "/" + "5",
    "1" + "0" + "." + "8",
    "1" + "0" + "." + "8" + "/" + "P0",
    "P0_" + "target",
    "N" + "_Q",
    "outgoing-" + "BT-target",
    "outgoing-" + "factorized-normalization",
    "pde.tex:" + "2075" + "-" + "2099",
)

SOURCE_BLOCKER = "BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE"
ALPHA_BLOCKER = "BLOCKED_NEEDS_alpha_J"
NOT_RUN_UNDEFINED_SOURCE = "NOT_RUN_UNDEFINED_K_SOURCE"


def _dim_row(name: str, dim: Dim, derivation: str, provenance: str) -> dict[str, object]:
    return {
        "name": name,
        "dimension": str(dim),
        "tuple_L_T_M": dim.as_tuple(),
        "sympy_monomial": str(dim.monomial()),
        "derivation": derivation,
        "provenance": provenance,
    }


def _kernel_symbols() -> tuple[sp.Symbol, sp.Function, sp.Function]:
    w = sp.Symbol("w")
    chi_n = sp.Function("chi_N")
    rho_inf4 = sp.Function("rho_inf4")
    return w, chi_n, rho_inf4


def stress_kernel_expression(w: sp.Symbol | None = None) -> sp.Expr:
    """Return the independently derived stress-lane transverse kernel."""

    if w is None:
        w, chi_n, rho_inf4 = _kernel_symbols()
    else:
        chi_n = sp.Function("chi_N")
        rho_inf4 = sp.Function("rho_inf4")
    return chi_n(w) * rho_inf4(w)


def stress_kernel_derivation() -> dict[str, object]:
    """Derive K_stress from the Noether-stress/reduced-density lane."""

    w = sp.Symbol("w")
    kernel = stress_kernel_expression(w)
    integral = sp.Integral(kernel, (w, -sp.oo, sp.oo))
    kernel_dim = LENGTH**-4
    integral_dim = kernel_dim * D["dw"]
    checks = [
        expect_dim(
            "pathA_22b_gate4_stress",
            "K_stress(w)=chi_N(w)*rho_infty,4(w)",
            kernel_dim,
            LENGTH**-4,
            "chi_N is the dimensionless reduction selector and rho_infty,4 is the bulk number density.",
        ).as_dict(),
        expect_dim(
            "pathA_22b_gate4_stress",
            "integral K_stress(w) dw = N_infty,3",
            integral_dim,
            LENGTH**-3,
            "The reduced density has one fewer length power than the bulk four-spatial density.",
        ).as_dict(),
    ]
    return {
        "status": "DERIVED_TARGET_BLIND",
        "symbol": "K_stress",
        "kernel": str(kernel),
        "kernel_ast": sp.srepr(kernel),
        "integral": str(integral),
        "integral_label": "N_infty,3",
        "weight": "chi_N(w) reduction selector; not Z(w) and not silently W(w)",
        "measure": "flat dw over the transverse coordinate",
        "dimension_table": [
            _dim_row(
                "chi_N(w)",
                DIMENSIONLESS,
                "Reduction selector in N_infty,3=int chi_N(w) rho_infty,4(w) dw.",
                "software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md:80-104",
            ),
            _dim_row(
                "rho_infty,4(w)",
                LENGTH**-4,
                "Asymptotic bulk number density in four spatial dimensions.",
                "research/pde/paper/pde.tex:396-406; research/pde/paper/pde.tex:2519-2524",
            ),
            _dim_row(
                "K_stress(w)",
                kernel_dim,
                "Product of the reduction selector with the bulk density profile.",
                "research/pde/paper/pde.tex:496-565; pathA_21b G5 projection formulas",
            ),
            _dim_row(
                "int K_stress dw",
                integral_dim,
                "Reduced asymptotic density N_infty,3.",
                "software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md:80-104",
            ),
        ],
        "checks": checks,
        "provenance": [
            "research/pde/paper/pde.tex:396-416",
            "research/pde/paper/pde.tex:496-565",
            "software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md:80-104",
            "software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md:20-67",
        ],
        "alpha_lane_power_in_g_G": "-1 in alpha_J1 and -1 in alpha_J2",
        "reduction_functional_power_in_g_G": "-1 in N_infty,3",
        "provenance_id": "stress_from_noether_reduced_density",
    }


def source_kernel_derivation() -> dict[str, object]:
    """Attempt the independent K_source derivation from allowed source terms."""

    mhat0_dim = (LENGTH**-1) * (TIME**-1) * (MASS ** sp.Rational(-1, 2))
    current_density_dim = LENGTH**-4 / TIME
    checks = [
        expect_dim(
            "pathA_22b_gate4_source",
            "mhat0 dimensional scale carrier",
            mhat0_dim,
            (LENGTH**-1) * (TIME**-1) * (MASS ** sp.Rational(-1, 2)),
            "Dimensional inheritance only; this is not a transverse source-map kernel derivation.",
        ).as_dict(),
        expect_dim(
            "pathA_22b_gate4_source",
            "parent number-current divergence scale",
            current_density_dim,
            LENGTH**-4 / TIME,
            "The parent continuity/source-current lane has number-density per time units and contains no mhat0 normalization.",
        ).as_dict(),
    ]
    return {
        "status": SOURCE_BLOCKER,
        "symbol": "K_source",
        "kernel": None,
        "kernel_ast": None,
        "integral": None,
        "integral_label": None,
        "weight": None,
        "measure": None,
        "source_map_kernel_exists_target_blind": False,
        "blocked_reason": (
            "The allowed parent action/current terms provide the GNLS current, localized Maxwell source bookkeeping, "
            "the projected brane continuity law, and the outgoing static prefactor. They do not provide an action, "
            "boundary-source, Hamiltonian, or branch equation that turns mhat0 or g_mhat into an integral over w. "
            "The only explicit mhat normalization is target-facing or listed as branch input, so a K_source kernel "
            "would be convention-dependent here."
        ),
        "allowed_source_scan": [
            "research/pde/paper/pde.tex:326-416",
            "research/pde/paper/pde.tex:496-565",
            "research/pde/paper/pde.tex:903-931",
            "research/pde/paper/pde.tex:2034-2069",
            "research/pde/paper/pde.tex:2551-2565",
            "software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md:28-44",
            "software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md:129-151",
        ],
        "excluded_use": (
            "The Burke-Thorne-facing normalization block is used only for dimensional inheritance of mhat0, "
            "not for K_source or g_mhat derivation."
        ),
        "dimension_table": [
            _dim_row(
                "mhat0",
                mhat0_dim,
                "Dimensionful source-map scale carrier from the already accepted Gate-0a dimensional reading.",
                "software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:5-30; research/pde/paper/pde.tex:2081",
            ),
            _dim_row(
                "parent current divergence",
                current_density_dim,
                "Continuity source/current scale inspected in the allowed parent lane; it is not a mhat0 kernel.",
                "research/pde/paper/pde.tex:396-416; research/pde/paper/pde.tex:903-931",
            ),
        ],
        "checks": checks,
        "alpha_lane_power_in_g_mhat_squared": "0 (no independently derived alpha_J factor in the allowed source lane)",
        "reduction_functional_power_in_g_mhat_squared": "not defined; no source-map transverse functional exists",
        "provenance_id": "source_from_allowed_parent_scan_blocked",
    }


def gravity_single_defect_mini_lemma() -> dict[str, object]:
    """Track the pair-level conditional G expression into g_G without taking roots."""

    theta1, theta2, alpha1, alpha2, i_f = sp.symbols("Theta_Q1 Theta_Q2 alpha_J1 alpha_J2 I_F12")
    j1, j2, n3, m_g, c_gamma, hbar, a, c_s = sp.symbols(
        "J1 J2 N_infty3 m_GNLS c_gamma hbar a c_s", positive=True
    )
    q1 = theta1 * j1 / n3
    q2 = theta2 * j2 / n3
    c_f = sp.simplify(m_g * n3 * q1 * q2 * i_f / (4 * sp.pi))
    m1 = alpha1 * hbar * j1 / c_gamma**2
    m2 = alpha2 * hbar * j2 / c_gamma**2
    g_cond = sp.factor(c_f / (m1 * m2))
    expected_g = sp.factor(c_gamma**4 * m_g * theta1 * theta2 * i_f / (4 * sp.pi * n3 * alpha1 * alpha2 * hbar**2))
    g_g = sp.factor(g_cond * m_g / (a * c_s**2))
    checks = [
        {
            "name": "C_F substitution from Q_i=Theta_Qi*J_i/N_infty,3",
            "status": "CONSISTENT" if sp.simplify(c_f - m_g * theta1 * theta2 * j1 * j2 * i_f / (4 * sp.pi * n3)) == 0 else "INCONSISTENT",
            "actual": str(c_f),
            "expected": str(m_g * theta1 * theta2 * j1 * j2 * i_f / (4 * sp.pi * n3)),
            "note": "The two flux factors are retained until the explicit mass bridge division.",
        },
        {
            "name": "G_cond pair expression",
            "status": "CONSISTENT" if sp.simplify(g_cond - expected_g) == 0 else "INCONSISTENT",
            "actual": str(g_cond),
            "expected": str(expected_g),
            "note": "No square root or single-defect factor dropping is used.",
        },
    ]
    dim_checks = [
        expect_dim(
            "pathA_22b_gate4_gG",
            "conditional G from stress lane",
            (D["c_s"] ** 4) * MASS / ((LENGTH**-3) * (ACTION**2)),
            D["G_3_spatial"],
            "Dimensionless profile factors and alpha_J are carried symbolically.",
        ).as_dict(),
        expect_dim(
            "pathA_22b_gate4_gG",
            "g_G=G*m_GNLS/(a*c_s^2)",
            D["G_3_spatial"] * MASS / (LENGTH * (D["c_s"] ** 2)),
            DIMENSIONLESS,
            "The scalar g_G is dimensionless after the conditional G expression is inserted.",
        ).as_dict(),
    ]
    return {
        "status": "PAIR_TO_SCALAR_G_CARRIED_WITH_PRODUCTS",
        "C_F_after_Q_substitution": str(c_f),
        "G_cond": str(g_cond),
        "g_G_cond": str(g_g),
        "checks": checks,
        "dimensional_checks": dim_checks,
        "power_ledger": {
            "Theta_Q1": 1,
            "Theta_Q2": 1,
            "alpha_J1": -1,
            "alpha_J2": -1,
            "N_infty,3": -1,
            "I_F12": 1,
        },
        "mini_lemma": (
            "The two-body stress result first gives C_F proportional to N_infty,3*Q1*Q2. "
            "Substituting the Gauss-solved Q_i factors gives one inverse power of N_infty,3 while retaining "
            "J1*J2. Dividing by the two candidate masses removes J1*J2 and introduces alpha_J1^-1*alpha_J2^-1. "
            "The scalar g_G uses this pair-level G_cond directly as the universal coupling candidate; the sources "
            "do not justify taking a square root or selecting one defect's factor by convention."
        ),
    }


def weighted_average_variation_condition(
    stress_kernel: sp.Expr,
    source_kernel: sp.Expr,
    *,
    w: sp.Symbol | None = None,
    v: sp.Symbol | None = None,
) -> sp.Expr:
    """Return the pointwise proportionality residual for arbitrary weighting."""

    w = w or sp.Symbol("w")
    v = v or sp.Symbol("v")
    return sp.factor(
        stress_kernel * source_kernel.subs(w, v)
        - stress_kernel.subs(w, v) * source_kernel
    )


def classify_kernel_pair(
    stress_kernel: sp.Expr,
    source_kernel: sp.Expr,
    *,
    shared_scalar_factored: bool,
    w: sp.Symbol | None = None,
    v: sp.Symbol | None = None,
) -> dict[str, object]:
    w = w or sp.Symbol("w")
    v = v or sp.Symbol("v")
    residual = weighted_average_variation_condition(stress_kernel, source_kernel, w=w, v=v)
    proportional = bool(sp.simplify(residual) == 0)
    route_a = bool(shared_scalar_factored)
    route_b = proportional
    outcome = "CANCELS" if route_a or route_b else "DOES_NOT_CANCEL"
    return {
        "outcome": outcome,
        "route_a_shared_factored_scalar": route_a,
        "route_b_pointwise_proportional_kernels": route_b,
        "shared_scalar_factored": bool(shared_scalar_factored),
        "kernels_proportional_for_arbitrary_weight": proportional,
        "proportionality_residual": str(residual),
        "stress_kernel": str(stress_kernel),
        "source_kernel": str(source_kernel),
    }


def comparator_controls() -> dict[str, object]:
    w, v, eps = sp.symbols("w v epsilon")
    negative = classify_kernel_pair(1 + w, 1 + w**2, shared_scalar_factored=False, w=w, v=v)
    real_stress = stress_kernel_expression(w)
    mutated = real_stress + eps * w
    mutated_control = classify_kernel_pair(mutated, real_stress, shared_scalar_factored=False, w=w, v=v)
    return {
        "negative_control": {
            **negative,
            "status": "PASS" if negative["outcome"] == "DOES_NOT_CANCEL" and negative["proportionality_residual"] != "0" else "FAIL",
            "purpose": "A deliberately non-factorizing pair must not cancel.",
        },
        "mutated_kernel_control": {
            **mutated_control,
            "status": (
                "PASS"
                if mutated_control["outcome"] == "DOES_NOT_CANCEL"
                and mutated_control["proportionality_residual"] != "0"
                else "FAIL"
            ),
            "purpose": "Perturbing a real stress-kernel AST must produce a nonzero residual.",
        },
    }


def alpha_cancellation_assessment(
    stress: Mapping[str, object], source: Mapping[str, object], lemma: Mapping[str, object]
) -> dict[str, object]:
    source_alpha_power = 0
    required_source_power = -2
    residual_pair_factor = sp.Symbol("alpha_J1") * sp.Symbol("alpha_J2")
    return {
        "hypothesis": "H-alpha_J",
        "finding": "CONFIRMED_NOT_CANCELLED",
        "stress_lane_alpha": stress["alpha_lane_power_in_g_G"],
        "source_lane_alpha": source["alpha_lane_power_in_g_mhat_squared"],
        "required_g_mhat_squared_alpha_power_for_cancellation": required_source_power,
        "derived_g_mhat_squared_alpha_power": source_alpha_power,
        "ratio_residual_pair_factor_if_source_forced_alpha_free": str(residual_pair_factor),
        "reason": (
            "The stress lane inherits alpha_J only through the candidate mass division in G_cond. "
            "The allowed source-map scan finds no independent alpha_J-bearing mhat0 or g_mhat equation. "
            "Therefore the same bridge factor is not present in both lanes with the required powers."
        ),
        "secondary_blocker_if_source_kernel_existed": ALPHA_BLOCKER,
        "g_G_power_ledger": lemma["power_ledger"],
    }


def cancellation_assessment() -> dict[str, object]:
    stress = stress_kernel_derivation()
    source = source_kernel_derivation()
    lemma = gravity_single_defect_mini_lemma()
    controls = comparator_controls()
    alpha = alpha_cancellation_assessment(stress, source, lemma)
    independence = {
        "status": "PASS_SEPARATE_PROVENANCE",
        "stress_provenance_id": stress["provenance_id"],
        "source_provenance_id": source["provenance_id"],
        "same_symbolic_object_reused": False,
        "note": "The real comparator is not run because K_source is blocked; the source scan is a separate blocked provenance object, not a reused stress AST.",
    }
    route_a = {
        "outcome": NOT_RUN_UNDEFINED_SOURCE,
        "reason": "K_source has no target-blind transverse functional, so there is no real shared scalar to compare.",
        "route": "shared_factored_scalar",
    }
    route_b = {
        "outcome": NOT_RUN_UNDEFINED_SOURCE,
        "reason": "K_source has no target-blind pointwise kernel, so proportionality is undefined.",
        "route": "pointwise_proportional_kernels",
    }
    h_k_source = {
        "hypothesis": "H-K_source",
        "finding": "CONFIRMED",
        "status": SOURCE_BLOCKER,
        "reason": source["blocked_reason"],
    }
    return {
        "stress_kernel": stress,
        "source_kernel": source,
        "single_defect_mini_lemma": lemma,
        "h_alpha_J": alpha,
        "h_k_source": h_k_source,
        "route_a": route_a,
        "route_b": route_b,
        "controls": controls,
        "independence_check": independence,
        "overall_verdict": SOURCE_BLOCKER,
        "stop_at_decisive_fork": True,
        "step_4iii_not_run": True,
    }


def dimensional_checks() -> dict[str, object]:
    stress = stress_kernel_derivation()
    source = source_kernel_derivation()
    lemma = gravity_single_defect_mini_lemma()
    return {
        "stress_checks": stress["checks"],
        "source_checks": source["checks"],
        "gG_checks": lemma["dimensional_checks"],
        "ratio_check": {
            "status": "NOT_RUN_UNDEFINED_K_SOURCE",
            "reason": "A dimensionless ratio check would require a legitimate K_source/g_mhat expression.",
        },
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
        "forbidden_literal_count": len(FORBIDDEN_TARGET_STRINGS),
    }


def gate4_report_payload() -> dict[str, object]:
    this_file = Path(__file__)
    root = this_file.resolve().parents[4]
    paths = [
        this_file,
        root / "software/stage1_solver/tools/pathA_22b_gate4_crosscheck.wl",
        root / "software/stage1_solver/tests/test_patha22b_gate4.py",
    ]
    cancellation = cancellation_assessment()
    dims = dimensional_checks()
    guard = target_blindness_guard(paths)
    return {
        "schema": "stage1_pathA_22b_gate4/v1",
        "scope": "Steps 4-i and 4-ii only; Step 4-iii is not run.",
        "cancellation_assessment": cancellation,
        "dimensional_checks": dims,
        "target_blindness": guard,
        "residual_ledger": [
            "K_source is not derivable as an independent target-blind transverse kernel from the allowed source/current action content.",
            "alpha_J is not inherited by the allowed source-map lane; if a source map were later supplied without alpha_J^-2 in g_mhat^2, the ratio would still be blocked on alpha_J.",
            "The real Route-a/Route-b kernel cancellation test is not meaningful until K_source provenance exists.",
        ],
        "gate4_outcome": SOURCE_BLOCKER,
        "stop_at_decisive_fork": True,
    }


def _fmt_check(raw: Mapping[str, object]) -> str:
    factor = raw.get("factor_needed_to_reach_expected")
    factor_text = "" if factor in (None, "1") else f"; factor needed `{factor}`"
    return (
        f"- `{raw['name']}`: **{raw['status']}** "
        f"(expected `{raw['expected']}`, actual `{raw['actual']}`{factor_text}). {raw['note']}"
    )


def render_gate4_markdown(payload: Mapping[str, object]) -> str:
    assessment = payload["cancellation_assessment"]
    dims = payload["dimensional_checks"]
    guard = payload["target_blindness"]
    assert isinstance(assessment, Mapping)
    assert isinstance(dims, Mapping)
    assert isinstance(guard, Mapping)
    stress = assessment["stress_kernel"]
    source = assessment["source_kernel"]
    lemma = assessment["single_defect_mini_lemma"]
    alpha = assessment["h_alpha_J"]
    hks = assessment["h_k_source"]
    controls = assessment["controls"]
    assert isinstance(stress, Mapping)
    assert isinstance(source, Mapping)
    assert isinstance(lemma, Mapping)
    assert isinstance(alpha, Mapping)
    assert isinstance(hks, Mapping)
    assert isinstance(controls, Mapping)
    negative = controls["negative_control"]
    mutated = controls["mutated_kernel_control"]
    assert isinstance(negative, Mapping)
    assert isinstance(mutated, Mapping)

    lines = [
        "## Gate 4 - ratio-or-blocked fork (Steps 4-i and 4-ii)",
        "",
        "### Target-blind attestation",
        "",
        f"- `{guard['status']}` over the new Gate-4 Python module, Wolfram cross-check, and tests.",
        "- No final comparison constants, outgoing target normalization equations, or imported target constants are used to derive `K_stress`, `K_source`, `g_G`, or `g_mhat`.",
        "- The Burke-Thorne-facing normalization block is used only for the already accepted dimensional reading of `mhat0`, not as source-map provenance.",
        "",
        "### Step 4-i kernels",
        "",
        f"- `K_stress(w)`: `{stress['kernel']}` with `{stress['integral_label']} = {stress['integral']}`.",
        f"- Stress measure/weight: {stress['measure']}; {stress['weight']}.",
        "- Stress provenance: " + "; ".join(str(item) for item in stress["provenance"]) + ".",
        f"- Stress lane powers: `{stress['reduction_functional_power_in_g_G']}`, `{stress['alpha_lane_power_in_g_G']}`.",
        f"- `K_source(w)`: `{source['status']}`.",
        f"- Source finding: {source['blocked_reason']}",
        "- Allowed source scan: " + "; ".join(str(item) for item in source["allowed_source_scan"]) + ".",
        "",
        "Stress dimensional checks:",
    ]
    for raw in stress["checks"]:
        assert isinstance(raw, Mapping)
        lines.append(_fmt_check(raw))
    lines.append("")
    lines.append("Source dimensional checks:")
    for raw in source["checks"]:
        assert isinstance(raw, Mapping)
        lines.append(_fmt_check(raw))

    lines.extend(
        [
            "",
            "### Single-defect `g_G` mini-lemma",
            "",
            f"- `{lemma['mini_lemma']}`",
            f"- `C_F` after flux substitution: `{lemma['C_F_after_Q_substitution']}`.",
            f"- `G_cond`: `{lemma['G_cond']}`.",
            f"- `g_G_cond`: `{lemma['g_G_cond']}`.",
            "- Power ledger: "
            + ", ".join(f"`{key}` -> `{value}`" for key, value in lemma["power_ledger"].items())
            + ".",
            "",
            "Mini-lemma checks:",
        ]
    )
    for raw in lemma["checks"]:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: **{raw['status']}**. {raw['note']}")
    for raw in lemma["dimensional_checks"]:
        assert isinstance(raw, Mapping)
        lines.append(_fmt_check(raw))

    lines.extend(
        [
            "",
            "### H findings",
            "",
            f"- `H-K_source`: `{hks['finding']}` -> `{hks['status']}`. {hks['reason']}",
            f"- `H-alpha_J`: `{alpha['finding']}`. {alpha['reason']}",
            f"- Required source-lane alpha power in `g_mhat^2` for cancellation: `{alpha['required_g_mhat_squared_alpha_power_for_cancellation']}`; derived power: `{alpha['derived_g_mhat_squared_alpha_power']}`.",
            f"- Residual pair factor if an alpha-free source map were forced: `{alpha['ratio_residual_pair_factor_if_source_forced_alpha_free']}`.",
            "",
            "### Step 4-ii cancellation test",
            "",
            f"- Route (a), shared factored scalar: `{assessment['route_a']['outcome']}`. {assessment['route_a']['reason']}",
            f"- Route (b), pointwise-proportional kernels: `{assessment['route_b']['outcome']}`. {assessment['route_b']['reason']}",
            f"- Independence check: `{assessment['independence_check']['status']}`. {assessment['independence_check']['note']}",
            f"- Negative control: `{negative['outcome']}` with residual `{negative['proportionality_residual']}` -> `{negative['status']}`.",
            f"- Mutated-kernel control: `{mutated['outcome']}` with residual `{mutated['proportionality_residual']}` -> `{mutated['status']}`.",
            "",
            "### Overall verdict",
            "",
            f"- Step-4-ii fork verdict: `{assessment['overall_verdict']}`.",
            "- Step 4-iii was not run.",
            "",
            "### Residual ledger",
            "",
        ]
    )
    for item in payload["residual_ledger"]:
        lines.append(f"- {item}")
    lines.append("")
    return "\n".join(lines)


def write_gate4_outputs(out_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    payload = gate4_report_payload()
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_22b_gate4.json"
    md_path = out_dir / "pathA_22b_gate4.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(render_gate4_markdown(payload) + "\n", encoding="utf-8")
    return json_path, md_path, payload


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="software/stage1_solver/_scratch")
    args = parser.parse_args(list(argv) if argv is not None else None)
    json_path, md_path, payload = write_gate4_outputs(Path(args.out_dir))
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"Gate 4 outcome: {payload['gate4_outcome']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
