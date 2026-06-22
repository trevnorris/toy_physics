"""PathA 22b Gate 0 foundation checks.

This module is intentionally additive.  It reconciles the source-map scale
dimension and checks the structural Z-kernel cancellation condition without
running any branch/profile solve.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Mapping

import sympy as sp

from stage1_solver.dimensional_check import (
    D,
    DIMENSIONLESS,
    LENGTH,
    MASS,
    TIME,
    Dim,
    expect_dim,
    homogeneous,
)


FORBIDDEN_TARGET_STRINGS = ("54" + "/5", "10.8" + "/P0")


def _dim_row(name: str, dim: Dim, derivation: str, provenance: str) -> dict[str, object]:
    return {
        "name": name,
        "dimension": str(dim),
        "tuple_L_T_M": dim.as_tuple(),
        "sympy_monomial": str(dim.monomial()),
        "derivation": derivation,
        "provenance": provenance,
    }


def mhat_reconciliation() -> dict[str, object]:
    """Reconcile the dimensionful normalization carrier with the natural map."""

    gamma5 = TIME**5
    bt_rhs = D["G_3_spatial"] / (D["c"] ** 5)
    factorized_rhs = D["G_3_spatial"] * (D["c_s"] ** 5) / ((D["a"] ** 5) * (D["c"] ** 5))
    mhat_required = (bt_rhs / gamma5) ** sp.Rational(1, 2)
    mhat0 = (LENGTH**-1) * (TIME**-1) * (MASS ** sp.Rational(-1, 2))
    natural_map_factor = DIMENSIONLESS
    natural_correction = (D["a"] ** 2) / (LENGTH**2)

    checks = [
        expect_dim(
            "pathA_22b_gate0a",
            "Gamma5 outgoing coefficient",
            gamma5,
            TIME**5,
            "Gamma5=chi_Q*P0*a^5/(const*c_s^5) with (c_s/a)^2-normalized (dimensionless) P0 per pathA_22a and dimensionless chi_Q.",
        ),
        expect_dim(
            "pathA_22b_gate0a",
            "BT odd-coefficient right side G/c^5",
            bt_rhs,
            (LENGTH**-2) * (TIME**3) * (MASS**-1),
            "The invariant odd-coefficient target carries the dimension that mhat^2*Gamma5 must match.",
        ),
        expect_dim(
            "pathA_22b_gate0a",
            "required mhat from mhat^2*Gamma5=G/c^5",
            mhat_required,
            mhat0,
            "Solving the unit equation fixes the source-map scale carrier.",
        ),
        expect_dim(
            "pathA_22b_gate0a",
            "factorized normalization right side",
            factorized_rhs,
            mhat0**2,
            "After substituting Gamma5, the right side has the dimension of mhat0^2.",
        ),
        homogeneous(
            "pathA_22b_gate0a",
            "dimensionful law",
            {
                "mhat0^2*Gamma5": (mhat0**2) * gamma5,
                "G/c^5": bt_rhs,
            },
            "This is the non-tautological unit reconciliation for the odd-coefficient law.",
        ),
        expect_dim(
            "pathA_22b_gate0a",
            "direct dimensionless mhat reading in odd-coefficient law",
            natural_map_factor**2 * gamma5,
            bt_rhs,
            "This is expected to fail if the dimensionless natural-map sentence is read as the whole mhat.",
        ),
        expect_dim(
            "pathA_22b_gate0a",
            "natural source-map correction a^2/r^2",
            natural_correction,
            DIMENSIONLESS,
            "The natural-map sentence is dimensionless and therefore belongs to the profile/source-map factor.",
        ),
    ]
    direct_dimensionless_check = checks[5]
    outcome = (
        "MHAT_DIMENSIONFUL_CONFIRMED"
        if checks[2].status == "CONSISTENT" and direct_dimensionless_check.status == "INCONSISTENT"
        else "MHAT_DIMENSIONLESS_REFRAME_SCALE_UNRESOLVED"
    )
    return {
        "outcome": outcome,
        "checks": [check.as_dict() for check in checks],
        "dimension_table": [
            _dim_row(
                "Gamma5",
                gamma5,
                "Odd outgoing coefficient dimension from the explicit a^5/c_s^5 carrier.",
                "(c_s/a)^2-normalized (dimensionless) P0 per pathA_22a; research/pde/paper/pde.tex:2053-2062",
            ),
            _dim_row(
                "G/c^5",
                bt_rhs,
                "Right side of the invariant odd-coefficient normalization.",
                "research/pde/paper/pde.tex:2075-2079",
            ),
            _dim_row(
                "mhat0",
                mhat0,
                "Scale carrier forced by [mhat0]^2*[Gamma5]=[G/c^5].",
                "research/pde/paper/pde.tex:2075-2083; software/stage1_solver/src/stage1_solver/dimensional_check.py:4280-4288",
            ),
            _dim_row(
                "natural source-map factor",
                natural_map_factor,
                "Dimensionless profile/source-map factor; the correction a^2/r^2 is dimensionless.",
                "research/pde/paper/pde.tex:2095-2099",
            ),
        ],
        "reasoning": (
            "The normalization law itself forces a dimensionful source-map scale. "
            "The later natural-map sentence is dimensionless, so it cannot be the whole object in that law. "
            "The consistent reading is mhat = mhat0*g_mhat, with mhat0 carrying units and g_mhat the "
            "dimensionless profile/source-map factor."
        ),
        "pathA_22a_effect": (
            "No flip to the pathA_22a decomposition: mhat0^2 remains the scale carrier and g_mhat remains "
            "a dimensionless branch/source-map residual until a natural-source-map derivation pins it."
        ),
        "provenance": [
            "research/pde/paper/pde.tex:2053-2062",
            "research/pde/paper/pde.tex:2075-2083",
            "research/pde/paper/pde.tex:2095-2099",
            "software/stage1_solver/src/stage1_solver/dimensional_check.py:4280-4288",
        ],
    }


def weighted_average_variation_condition(
    stress_kernel: sp.Expr,
    source_kernel: sp.Expr,
    *,
    w: sp.Symbol | None = None,
    v: sp.Symbol | None = None,
) -> sp.Expr:
    """Return the arbitrary-Z cancellation condition for two weighted kernels.

    For R = int Z*K_stress / int Z*K_source, the first variation with respect
    to Z is independent of Z only when the two kernels are proportional.  A
    polynomial identity form of that condition is
    K_stress(w)*K_source(v) - K_stress(v)*K_source(w) == 0 for all w,v.
    """

    w = w or sp.Symbol("w")
    v = v or sp.Symbol("v")
    return sp.factor(
        stress_kernel.subs(w, w) * source_kernel.subs(w, v)
        - stress_kernel.subs(w, v) * source_kernel.subs(w, w)
    )


def classify_kernel_pair(
    stress_kernel: sp.Expr,
    source_kernel: sp.Expr,
    *,
    shared_scalar_factored: bool,
    w: sp.Symbol | None = None,
    v: sp.Symbol | None = None,
) -> dict[str, object]:
    """Classify whether a Z-weighted kernel pair permits cancellation."""

    w = w or sp.Symbol("w")
    v = v or sp.Symbol("v")
    condition = weighted_average_variation_condition(stress_kernel, source_kernel, w=w, v=v)
    proportional = sp.simplify(condition) == 0
    route_a_shared_scalar = bool(shared_scalar_factored)
    route_b_proportional = bool(proportional)
    z_independent = route_a_shared_scalar or route_b_proportional
    outcome = "CANCELS" if z_independent else "DOES_NOT_CANCEL"
    return {
        "outcome": outcome,
        "shared_scalar_factored": shared_scalar_factored,
        "route_a_shared_factored_scalar": route_a_shared_scalar,
        "route_b_pointwise_proportional_kernels": route_b_proportional,
        "z_independent_under_available_route": z_independent,
        "kernels_proportional_for_arbitrary_Z": bool(proportional),
        "proportionality_condition": str(condition),
        "stress_kernel": str(stress_kernel),
        "source_kernel": str(source_kernel),
    }


DOES_NOT_CANCEL_SOURCE_LABEL = (
    "DOES_NOT_CANCEL (NOT_ESTABLISHED — sources do not establish either cancellation route; "
    "a later Gate-4 action-level derivation could still find one)"
)


def z_kernel_cancellation_source_assessment() -> dict[str, object]:
    """Assess whether available sources establish Z-kernel cancellation."""

    w, v = sp.symbols("w v")
    factored_control = classify_kernel_pair(
        sp.Symbol("K_stress"),
        sp.Symbol("K_source"),
        shared_scalar_factored=True,
        w=w,
        v=v,
    )
    negative_control = classify_kernel_pair(1 + w, 1 + w**2, shared_scalar_factored=False, w=w, v=v)
    proportional_weighted_control = classify_kernel_pair(2 * (1 + w**2), 1 + w**2, shared_scalar_factored=False, w=w, v=v)

    I_Z, K_stress, K_source = sp.symbols("I_Z K_stress K_source", nonzero=True)
    factored_ratio = sp.cancel((I_Z * K_stress) / (I_Z * K_source))
    source_establishes_shared_scalar = False
    source_establishes_kernel_proportionality = False
    sources_establish_condition = source_establishes_shared_scalar or source_establishes_kernel_proportionality
    outcome = "CANCELS" if sources_establish_condition else "DOES_NOT_CANCEL"
    outcome_label = "CANCELS" if outcome == "CANCELS" else DOES_NOT_CANCEL_SOURCE_LABEL

    checks = [
        {
            "name": "common scalar functional algebra",
            "status": "CONSISTENT" if factored_ratio == K_stress / K_source else "INCONSISTENT",
            "actual": str(factored_ratio),
            "expected": str(K_stress / K_source),
            "note": "If both factors are proven to contain the same factored I_Z, the scalar cancels algebraically.",
        },
        {
            "name": "non-factorizing weighted negative control",
            "status": "CONSISTENT" if negative_control["outcome"] == "DOES_NOT_CANCEL" else "INCONSISTENT",
            "actual": negative_control["outcome"],
            "expected": "DOES_NOT_CANCEL",
            "note": "Distinct kernels in int Z*K_stress / int Z*K_source leave dependence on Z.",
        },
        {
            "name": "weighted proportional control",
            "status": "CONSISTENT" if proportional_weighted_control["outcome"] == "CANCELS" else "INCONSISTENT",
            "actual": proportional_weighted_control["outcome"],
            "expected": "CANCELS",
            "note": "Pointwise proportional weighted kernels make int Z*K_stress / int Z*K_source independent of Z; this is route (b), distinct from shared scalar route (a).",
        },
    ]

    return {
        "outcome": outcome,
        "outcome_label": outcome_label,
        "structural_condition": (
            "Z-independence of g_mhat^2/g_G can be established by either route (a) shared factored scalar, "
            "where both g_G and g_mhat^2 are I_Z times separate w-independent field-content kernels with the "
            "same I_Z=int sqrt(g_w)*Z(w) dw, or route (b) pointwise-proportional kernels, where "
            "K_stress(w)=const*K_source(w) for all w. Route (b) is equivalently the identity "
            "K_stress(w)*K_source(v)=K_stress(v)*K_source(w) for all w,v, making the weighted-average ratio "
            "independent of the Z/W_eff profile."
        ),
        "sources_establish_condition": sources_establish_condition,
        "source_establishes_shared_scalar_route": source_establishes_shared_scalar,
        "source_establishes_pointwise_proportional_route": source_establishes_kernel_proportionality,
        "source_assessment": (
            "Gate 0 is a source-availability assessment; the action-level derivation of the g_G/g_mhat^2 "
            "kernels is deferred to Gate 4. "
            "The parent action places Z(w) in the localized Maxwell kinetic term and Maxwell equation, while "
            "the GNLS matter sector, exact current, source-map statements, and brane projection W(w) are distinct. "
            "The cited sources do not establish route (a), a shared factored I_Z functional for g_G and "
            "g_mhat^2, nor route (b), proportional stress/source kernels."
        ),
        "checks": checks,
        "factored_control": {
            "outcome": "CANCELS" if factored_ratio == K_stress / K_source else "DOES_NOT_CANCEL",
            "ratio_after_cancel": str(factored_ratio),
        },
        "negative_control": negative_control,
        "proportional_weighted_control": proportional_weighted_control,
        "implication_for_gate4": (
            "Before declaring W_eff/full transverse profile irreducible, Gate 4 must test both cancellation "
            "routes: route (a) shared factored scalar and route (b) pointwise-proportional kernels. Unless "
            "Gate 4 proves one of those routes from the action-level kernels, W_eff/full transverse profile "
            "remains on the critical path for the ratio."
        ),
        "provenance": [
            "research/pde/paper/pde.tex:277",
            "research/pde/paper/pde.tex:289-295",
            "research/pde/paper/pde.tex:357-416",
            "research/pde/paper/pde.tex:496-565",
            "software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md:80-104",
            "software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md",
        ],
    }


def dimensional_checks() -> dict[str, object]:
    mhat = (LENGTH**-1) * (TIME**-1) * (MASS ** sp.Rational(-1, 2))
    gamma5 = TIME**5
    rhs_bt = D["G_3_spatial"] / (D["c"] ** 5)
    iz = LENGTH
    z_weight = DIMENSIONLESS
    sqrt_gw_dw = LENGTH
    return {
        "checks": [
            homogeneous(
                "pathA_22b_gate0_dimensional",
                "mhat0^2*Gamma5 vs G/c^5",
                {
                    "mhat0^2*Gamma5": (mhat**2) * gamma5,
                    "G/c^5": rhs_bt,
                },
            ).as_dict(),
            expect_dim(
                "pathA_22b_gate0_dimensional",
                "I_Z=int sqrt(g_w) Z(w) dw scalar carrier",
                sqrt_gw_dw * z_weight,
                iz,
                "The exact length dimension depends on the w-coordinate convention; it cancels only after a proven common-factor derivation.",
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
        "forbidden_literal_count": len(FORBIDDEN_TARGET_STRINGS),
        "hits": hits,
    }


def gate0_report_payload() -> dict[str, object]:
    this_file = Path(__file__)
    wl_file = Path("software/stage1_solver/tools/pathA_22b_gate0_crosscheck.wl")
    mhat = mhat_reconciliation()
    zassessment = z_kernel_cancellation_source_assessment()
    dims = dimensional_checks()
    guard = target_blindness_guard([this_file, wl_file])
    residuals = [
        "The natural dimensionless source-map factor g_mhat is not dynamically pinned by Gate 0.",
        "The parent sources do not establish either route (a) a shared factored I_Z for g_G and g_mhat^2 or route (b) pointwise-proportional kernels.",
        "Gate 4 cannot use a reduced field-content ratio unless a later action-level derivation proves one of the two cancellation routes; otherwise it needs W_eff/full transverse kernels.",
    ]
    return {
        "schema": "stage1_pathA_22b_gate0/v1",
        "mhat_reconciliation": mhat,
        "z_kernel_cancellation_source_assessment": zassessment,
        "dimensional_checks": dims,
        "target_blindness": guard,
        "gate0_outcomes": {
            "0a": mhat["outcome"],
            "0b": zassessment["outcome"],
            "0b_label": zassessment["outcome_label"],
        },
        "residuals": residuals,
    }


def render_gate0_markdown(payload: Mapping[str, object]) -> str:
    mhat = payload["mhat_reconciliation"]
    zassessment = payload["z_kernel_cancellation_source_assessment"]
    dims = payload["dimensional_checks"]
    guard = payload["target_blindness"]
    assert isinstance(mhat, Mapping)
    assert isinstance(zassessment, Mapping)
    assert isinstance(dims, Mapping)
    assert isinstance(guard, Mapping)

    lines = [
        "# PathA 22b minimal combination xi",
        "",
        "## Gate 0",
        "",
        "### 0a mhat reconciliation",
        "",
        f"- Outcome: `{mhat['outcome']}`.",
        f"- Reasoning: {mhat['reasoning']}",
        f"- Effect on pathA_22a: {mhat['pathA_22a_effect']}",
        "- Provenance: " + "; ".join(str(item) for item in mhat["provenance"]) + ".",
        "",
        "| ingredient | dimension | SymPy monomial | derivation | provenance |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in mhat["dimension_table"]:
        assert isinstance(row, Mapping)
        lines.append(
            f"| `{row['name']}` | `{row['dimension']}` | `{row['sympy_monomial']}` | "
            f"{row['derivation']} | {row['provenance']} |"
        )
    lines.extend(["", "0a checks:"])
    for raw in mhat["checks"]:
        assert isinstance(raw, Mapping)
        factor = raw["factor_needed_to_reach_expected"]
        factor_text = "" if factor in (None, "1") else f"; factor needed `{factor}`"
        expected_note = " (expected negative)" if raw["name"] == "direct dimensionless mhat reading in odd-coefficient law" else ""
        lines.append(
            f"- `{raw['name']}`: **{raw['status']}**{expected_note} "
            f"(expected `{raw['expected']}`, actual `{raw['actual']}`{factor_text}). {raw['note']}"
        )

    lines.extend(
        [
            "",
            "### 0b Z-kernel cancellation source assessment",
            "",
            f"- Outcome: `{zassessment['outcome_label']}`.",
            "- Scope: Gate 0 is a source-availability assessment; the action-level derivation of the `g_G`/`g_mhat²` kernels is deferred to Gate 4.",
            f"- Exact structural condition: {zassessment['structural_condition']}",
            f"- Do the sources establish either route? `{zassessment['sources_establish_condition']}`. {zassessment['source_assessment']}",
            f"- Gate 4 implication: {zassessment['implication_for_gate4']}",
            "- Provenance: " + "; ".join(str(item) for item in zassessment["provenance"]) + ".",
            "",
            "0b checks and controls:",
        ]
    )
    for raw in zassessment["checks"]:
        assert isinstance(raw, Mapping)
        lines.append(
            f"- `{raw['name']}`: **{raw['status']}** (expected `{raw['expected']}`, actual `{raw['actual']}`). {raw['note']}"
        )
    negative = zassessment["negative_control"]
    assert isinstance(negative, Mapping)
    lines.append(
        "- Negative control: "
        f"`{negative['stress_kernel']}` vs `{negative['source_kernel']}` gives `{negative['outcome']}`; "
        f"condition residual `{negative['proportionality_condition']}`."
    )
    proportional = zassessment["proportional_weighted_control"]
    assert isinstance(proportional, Mapping)
    lines.append(
        "- Proportional-kernel control: "
        f"`{proportional['stress_kernel']}` vs `{proportional['source_kernel']}` gives `{proportional['outcome']}` by route (b); "
        f"condition residual `{proportional['proportionality_condition']}`."
    )
    factored = zassessment["factored_control"]
    assert isinstance(factored, Mapping)
    lines.append(
        f"- Positive algebra control: common scalar functional gives `{factored['outcome']}` with ratio `{factored['ratio_after_cancel']}`."
    )

    lines.extend(["", "### Dimensional checks", ""])
    for raw in dims["checks"]:
        assert isinstance(raw, Mapping)
        factor = raw["factor_needed_to_reach_expected"]
        factor_text = "" if factor in (None, "1") else f"; factor needed `{factor}`"
        lines.append(
            f"- `{raw['name']}`: **{raw['status']}** "
            f"(expected `{raw['expected']}`, actual `{raw['actual']}`{factor_text}). {raw['note']}"
        )
    lines.extend(
        [
            "",
            "### Target-blindness",
            "",
            f"- `{guard['status']}` over the new Gate-0 module and Mathematica cross-check.",
            "",
            "### Residual ledger",
            "",
        ]
    )
    for item in payload["residuals"]:
        lines.append(f"- {item}")
    lines.append("")
    return "\n".join(lines)


def write_gate0_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    payload = gate0_report_payload()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_22b_gate0.json"
    scratch_md_path = out_dir / "pathA_22b_gate0.md"
    report_path = report_dir / "pathA_22b_minimal_combination_xi.md"
    rendered = render_gate0_markdown(payload)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered + "\n", encoding="utf-8")
    report_path.write_text(rendered + "\n", encoding="utf-8")
    return json_path, report_path, payload


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="software/stage1_solver/_scratch")
    parser.add_argument("--report-dir", default="software/stage1_solver/reports")
    args = parser.parse_args(list(argv) if argv is not None else None)
    json_path, report_path, payload = write_gate0_report(Path(args.out_dir), Path(args.report_dir))
    outcomes = payload["gate0_outcomes"]
    print(f"wrote {json_path}")
    print(f"wrote {report_path}")
    print(f"Gate 0 outcomes: 0a={outcomes['0a']}; 0b={outcomes['0b_label']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
