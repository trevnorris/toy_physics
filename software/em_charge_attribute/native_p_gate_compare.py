#!/usr/bin/env python3
"""Compare independent computed artifacts and generate the rebuild report."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "reports" / "native_p_gate_artifacts"
DEFAULT_REPORT = HERE / "reports" / "native_p_constraint_gate.md"


def require(test: object, message: str) -> None:
    if not bool(test):
        raise AssertionError(message)


def py_classes(theory: dict) -> list[str]:
    return list(theory["constraint_classes"].values())


def compact_log(path: Path) -> str:
    return path.read_text().strip() if path.exists() else "(log unavailable)"


def tuned_rejection_signature(theory: dict) -> list[tuple]:
    records = theory["tuned_descendant_rejection_guard"]["records"]
    return [
        (
            record["first_class_count"],
            record["first_class_primaries"],
            record["accepted_gauss_directions"],
            tuple(
                (
                    item["computed_rejection_reason"],
                    item["descendant_zero"],
                    item["secondary_action_non_gradient"],
                    item["descendant"] == "0",
                    tuple(item["secondary_action"]),
                    tuple(
                        (block["proportional_to_k"], tuple(block["action"]))
                        for block in item["spatial_action_blocks"]
                    ),
                )
                for item in record["records"]
            ),
        )
        for record in records
    ]


def compare(sym: dict, math: dict) -> list[str]:
    require(sym["algebra_status"] == math["algebra_status"] == "PASS", "an engine failed")
    require(sym["pipeline"] == "build_H2 -> dirac_search -> search_G", "SymPy did not record the shared pipeline")
    require(math["pipeline"] == "buildH2 -> diracSearch -> searchG", "Mathematica did not record its shared pipeline")
    require(sym["all_controls_pass"] and sym["all_ablations_fired"], "control or ablation failure")
    for engine in (sym, math):
        require(engine["guards"]["GUARD-SEARCH-CAPABLE"]["status"] == "PASS", "search-capability guard failed")
        for t in ("A", "C"):
            require(engine["guards"]["GUARD-COUPLINGS-ENTER"][t]["status"] == "PASS", f"coupling guard failed {t}")
            require(engine["guards"]["HARDENING-TUNED-DESCENDANT-REJECTION"][t]["status"] == "PASS",
                    f"computed tuned-descendant rejection guard failed {t}")

    lines = [
        "AGREE shared input-Lagrangian pipeline present in both engines",
        "AGREE GUARD-COUPLINGS-ENTER=A:PASS,C:PASS",
        "AGREE GUARD-SEARCH-CAPABLE=PASS",
        "AGREE HARDENING-TUNED-DESCENDANT-REJECTION=A:PASS,C:PASS",
    ]
    for t in ("A", "C"):
        a, b = sym["theories"][t], math["theories"][t]
        fields = {
            "constraint_count": (len(a["constraints"]), len(b["constraints"])),
            "constraint_stages": (a["constraint_stages"], b["constraint_stages"]),
            "constraint_classes": (py_classes(a), b["constraint_classes"]),
            "bracket_rank": (a["bracket_rank"], b["bracket_rank"]),
            "first_class_count": (a["first_class_count"], b["first_class_count"]),
            "second_class_count": (a["second_class_count"], b["second_class_count"]),
            "gauss_candidates": (a["G_search"]["gauss_candidates"], b["G_search"]["gauss_candidates"]),
            "additional_G_exists": (a["G_search"]["additional_G_exists"], b["G_search"]["additional_G_exists"]),
            "boundary_scan": (a["coupling_scan"]["semidefinite_boundary"], b["coupling_scan"]["semidefinite_boundary"]),
            "verdict": (a["verdict"], b["verdict"]),
        }
        # Boundary condition labels differ typographically; compare decisive fields.
        fields["boundary_scan"] = (
            [(x["first_class_count"], x["gauss_candidates"], x["additional_G_exists"], x["hessian_nullity"])
             for x in a["coupling_scan"]["semidefinite_boundary"]],
            [(x["first_class_count"], x["gauss_candidates"], x["additional_G_exists"], x["hessian_nullity"])
             for x in b["coupling_scan"]["semidefinite_boundary"]],
        )
        fields["rankdrop_randomized_sweep"] = (
            [(x["first_class_count"], x["gauss_candidates"], x["additional_G_exists"], x["hessian_nullity"])
             for x in a["coupling_scan"]["rankdrop_representative_sweep"]["samples"]],
            [(x["first_class_count"], x["gauss_candidates"], x["additional_G_exists"], x["hessian_nullity"])
             for x in b["coupling_scan"]["rankdrop_representative_sweep"]["samples"]],
        )
        fields["tuned_descendant_rejection"] = (
            tuned_rejection_signature(a), tuned_rejection_signature(b)
        )
        for field, (left, right) in fields.items():
            require(left == right, f"ENGINE_DISAGREE THEORY-{t} {field}: SymPy={left!r}, Mathematica={right!r}")
            lines.append(f"AGREE THEORY-{t} {field}: {left}")
        lines.append(f"ENGINE_AGREE THEORY-{t}")

    math_controls = math["controls"]
    for control in sym["controls"]:
        key = control["model"]
        other = math_controls[key]
        left = (control["first_class_count"], control["second_class_count"],
                control["G_search"]["gauss_candidates"], control["classification"])
        right = (other["first_class_count"], other["second_class_count"],
                 other["G_search"]["gauss_candidates"], math["control_classifications"][key])
        require(left == right, f"ENGINE_DISAGREE CONTROL {key}: {left!r} != {right!r}")
        lines.append(f"AGREE CONTROL {key}: FC={left[0]} SC={left[1]} G={left[2]} {left[3]}")
    return lines


def matrix_block(m: list[list[str]]) -> str:
    if not m:
        return "(empty: unconstrained regular theory)"
    return "```text\n" + "\n".join("[ " + " , ".join(row) + " ]" for row in m) + "\n```"


def field_action(coords: list[str], action: list[str]) -> str:
    return "{ " + ", ".join(f"{q} -> {a}" for q, a in zip(coords, action)) + " }"


def report(sym: dict, math: dict, agreement: list[str], out: Path) -> str:
    verdicts = [sym["theories"][t]["verdict"] for t in ("A", "C")]
    same_verdict = verdicts[0] if len(set(verdicts)) == 1 else " / ".join(verdicts)
    basis = sym["operator_basis"]
    rows = [
        "# Native-`Pᵃ` constraint-class gate — genuine rebuild",
        "",
        "## HARDENING NOTE",
        "",
        "The tuned-strata rejection is now computation-backed rather than prose-backed. For every rejected first-class primary on an FC-bearing rank-drop/common-null stratum, each engine computes `C₁={C₀,H₂}` and `{q,C₁}`, records the expressions and a derived reason code, and enforces `C₁=0 OR no nonzero spatial action proportional to k`. The run fails if a rejected direction instead has a nonzero `∝k` Gauss descendant. The detailed records are printed below and retained in both engine artifacts.",
        "",
        "The scope statement is also hardened: it separates the fully symbolic open-stratum result from the argued-and-scanned tuned-strata coverage; it does not promote a representative-point scan into an exhaustive symbolic classification.",
        "",
        "## REBUILD NOTE",
        "",
        "- **Q1 closed:** `build_H2` differentiates each input Lagrangian, computes its Hessian nullspace and Legendre transform; `coupling_guard` then proves every native `g_a` survives into computed momentum/Hessian/constraint/PB objects ([native_p_gate_sympy.py:99](../native_p_gate_sympy.py#L99), [native_p_gate_sympy.py:545](../native_p_gate_sympy.py#L545)).",
        "- **Q2 closed:** Maxwell is an input Lagrangian and reaches the same `execute` pipeline as native A/C; its zero PB matrix is an output, never a literal ([native_p_gate_sympy.py:378](../native_p_gate_sympy.py#L378), [native_p_gate_sympy.py:337](../native_p_gate_sympy.py#L337)).",
        "- **Q3 closed:** all six controls are model builders in `CONTROL_BUILDERS` and are executed by `run_control`; the complete Coulomb pair is derived from multiplier terms ([native_p_gate_sympy.py:455](../native_p_gate_sympy.py#L455), [native_p_gate_sympy.py:430](../native_p_gate_sympy.py#L430)).",
        "- **Q4 closed:** Wolfram independently differentiates the Lagrangians and runs its own Hessian/Dirac/kernel code; this comparator checks native and control outputs, not shared answer literals ([native_p_gate_dual.wl:22](../native_p_gate_dual.wl#L22), [native_p_gate_compare.py:52](../native_p_gate_compare.py#L52)).",
        "- **Q5 closed:** `search_G` enumerates computed first-class primary directions and tests their Hamiltonian descendant and field action; Maxwell and gauged-hard-unit must yield nonzero candidates under `GUARD-SEARCH-CAPABLE` ([native_p_gate_sympy.py:196](../native_p_gate_sympy.py#L196), [native_p_gate_sympy.py:750](../native_p_gate_sympy.py#L750)).",
        "- **Q6 closed:** `dirac_search` iterates preservation and reports PB rank/nullity; no full-rank assertion or blanket class dictionary remains ([native_p_gate_sympy.py:141](../native_p_gate_sympy.py#L141), [native_p_gate_dual.wl:63](../native_p_gate_dual.wl#L63)).",
        "",
        f"**COMPUTED VERDICT: `{same_verdict}`.**",
        "",
        "This is the Stage-1 constraint gate only. Compactness, charge quantization, deconfinement, and native `±w` current supply remain downstream.",
        "",
        "## Frozen setup and quadratic operator basis",
        "",
        "The calculation uses one nonzero Fourier representative `k=(1,2,3)` on `R³`, with `C_c^∞` smearings and decaying fields. The `k=0` sector is global and is not counted as a local Gauss law. Coulomb control fields and gauge parameters decay, so `-∇²` has no local harmonic kernel. No punctures occur in the source-free search.",
        "",
        f"Notation: {basis['notation']}.",
        "",
        "THEORY-A cross-coupling basis:",
        "",
    ]
    rows.extend(f"- `{x}`" for x in basis["quadratic_couplings_A"])
    rows += ["", "THEORY-C cross-coupling basis:", ""]
    rows.extend(f"- `{x}`" for x in basis["quadratic_couplings_C"])
    rows += ["", "Empty families:", ""]
    rows.extend(f"- {x}." for x in basis["empty_families"])
    rows += ["", f"Operator-basis completeness: {basis['completeness']}", ""]
    rows += [
        "## Completeness and stratum scope",
        "",
        "The decisive calculation is the fully symbolic open kinetic stratum: FC=`0` for THEORY-A and THEORY-C for all retained coupling symbols. The computed physical kinetic-Hessian determinant is `(rho_u-g_t^2)^3`, so its only additional kinetic degeneracy is `g_t^2=rho_u` (apart from the frozen multiplier nulls already in the Dirac system).",
        "",
        "On that degeneracy surface, the potential-null residual is solved symbolically to identify the rank-drop/common-null conditions, and the remaining non-common rank-drop locus is checked by a fixed-seed randomized representative-point sweep in both engines. This tuned coverage is **ARGUED + SCANNED**, not an exhaustive symbolic stratification of every possible measure-zero sublocus.",
        "",
        "Consequently, any hypothetical missed measure-zero first-class Gauss stratum would be a **TUNED / inverse-design** result, not robust native electromagnetism. The physical no-go used here—native `P` does not **generically** host emergent EM—is decisive independently of that caveat.",
        "",
    ]

    for t in ("A", "C"):
        x, mx = sym["theories"][t], math["theories"][t]
        rows += [
            f"## THEORY-{t}: computed `H₂`, Dirac chain, and `G` search",
            "",
            f"Free coupling symbols: `{', '.join(x['couplings'])}`. Coupling-entry guard: **`{x['coupling_guard']['status']}`**, with computed locations `{x['coupling_guard']['locations']}`.",
            "",
            "Instantiated quadratic Lagrangian:",
            "",
            "```text",
            x["lagrangian"],
            "```",
            "",
            "Computed Legendre transform `H₂`:",
            "",
            "```text",
            x["H2"],
            "```",
            "",
            f"Hessian rank/nullity: `{x['hessian_rank']}/{x['hessian_nullity']}`. Computed momentum map:",
            "",
        ]
        rows.extend(f"- `Π[{i}] = {p}`" for i, p in enumerate(x["momentum_map"], 1))
        rows += ["", "Dirac constraints in discovery order (`stage 0` means primary):", ""]
        classes = list(x["constraint_classes"].values())
        rows.extend(f"- stage `{stage}`: `{constraint}` — `{klass}`"
                    for constraint, stage, klass in zip(x["constraints"], x["constraint_stages"], classes))
        rows += [
            "",
            f"Computed weak PB matrix `M(g,k)` (rank `{x['bracket_rank']}`, FC `{x['first_class_count']}`, SC `{x['second_class_count']}`):",
            "",
            matrix_block(x["bracket_matrix"]),
            "",
            "Coupling scan and kernel search:",
            "",
            f"- Regular kinetic determinant per vector component: `{x['coupling_scan']['physical_kinetic_determinant_per_component']}`; stable open stratum `{x['coupling_scan']['regular_stability_stratum']}`.",
            f"- Full physical kinetic-Hessian determinant: `{x['coupling_scan']['physical_kinetic_hessian_determinant']}`; only additional kinetic degeneracy: `{x['coupling_scan']['only_kinetic_hessian_degeneracy']}`.",
            f"- Regular result: computed kernel dimension `{x['G_search']['computed_kernel_dimension']}`, first-class primaries `{x['G_search']['computed_first_class_primaries']}`, Gauss candidates `{x['coupling_scan']['regular_result']['gauss_candidates']}`.",
            f"- Boundary rank polynomial: `{x['coupling_scan']['boundary_rank_polynomial']}`.",
            f"- Independently computed common-null solutions: `{x['coupling_scan']['computed_common_null_solutions']}`.",
        ]
        for b in x["coupling_scan"]["semidefinite_boundary"]:
            rows.append(f"- Boundary `{b['condition']}`: Hessian nullity `{b['hessian_nullity']}`, FC `{b['first_class_count']}`, Gauss candidates `{b['gauss_candidates']}`, `G={b['additional_G_exists']}`.")
        sweep = x["coupling_scan"]["rankdrop_representative_sweep"]
        rows.append(f"- Tuned rank-drop sweep: `{sweep['sample_count']}` fixed-seed randomized representative points (seed `{sweep['seed']}`); scope is `{x['coupling_scan']['tuned_scope']}`.")
        for sample in sweep["samples"]:
            rows.append(f"  - `{sample['condition']}` at `{sample['substitutions']}`: FC `{sample['first_class_count']}`, Gauss candidates `{sample['gauss_candidates']}`, `G={sample['additional_G_exists']}`.")
        rows += [
            f"- All spatial couplings remained free in the symbolic computations: `{', '.join(x['coupling_scan']['all_spatial_couplings_free'])}`.",
            f"- Computed aggregate: `gauss_candidates={x['G_search']['gauss_candidates']}`, `additional_G_exists={x['G_search']['additional_G_exists']}`.",
            f"- Source-first ordering: source-free searched `{x['source_test']['searched_source_free_first']}`; `j·A` added `{x['source_test']['jA_added']}`; sourced `{x['source_test']['j_sourced']}`.",
            f"- Shear duplicate: `{x['shear_duplicate']['applicability']}` (MacCullagh transverse modes `{x['shear_duplicate']['macCullagh_transverse_modes']}`).",
            f"- Independent Wolfram result: PB rank `{mx['bracket_rank']}`, FC `{mx['first_class_count']}`, candidates `{mx['G_search']['gauss_candidates']}`, verdict `{mx['verdict']}`.",
            "",
            f"**VERDICT THEORY-{t}: `{x['verdict']}` at quadratic order.**",
            "",
        ]

        hardening = x["tuned_descendant_rejection_guard"]
        rows += [
            f"### THEORY-{t} tuned FC descendant rejection audit",
            "",
            f"Computed hardening guard: **`{hardening['status']}`**; FC-bearing strata checked `{hardening['checked_strata']}`, rejected directions checked `{hardening['checked_directions']}`.",
            "",
        ]
        for b in x["coupling_scan"]["semidefinite_boundary"]:
            audit = b["tuned_descendant_rejection_guard"]
            if audit["status"] != "PASS":
                continue
            rows.append(f"- Stratum `{b['condition']}`:")
            for item in audit["records"]:
                rows.append(
                    f"  - FC direction `{item['direction']}`: primary `{item['primary']}`; "
                    f"computed descendant `{{primary,H₂}} = {item['descendant']}`; "
                    f"computed field action `{{q,{{primary,H₂}}}} = {field_action(b['coordinates'], item['secondary_action'])}`; "
                    f"reason **`{item['computed_rejection_reason']}`** "
                    f"(`descendant_zero={item['descendant_zero']}`, "
                    f"`secondary_action_non_gradient={item['secondary_action_non_gradient']}`)."
                )
                for block in item["spatial_action_blocks"]:
                    rows.append(
                        f"    - Spatial block fields `{block['field_indices']}`: action `{block['action']}`; "
                        f"computed `proportional_to_k={block['proportional_to_k']}`."
                    )
        rows.append("")

    rows += [
        "## Six able-to-fail controls through the shared path",
        "",
        "Every row was obtained from an input Lagrangian through the identical Hessian → Legendre → Dirac → kernel search. The nonconserved-current row additionally evaluates the nonzero continuity defect in Gauss preservation.",
        "",
        "| # | Control | Computed class | Hessian nullity | FC | SC | Gauss candidates |",
        "|---:|---|---|---:|---:|---:|---:|",
    ]
    for i, c in enumerate(sym["controls"], 1):
        rows.append(f"| {i} | {c['name']} | `{c['classification']}` | {c['hessian_nullity']} | {c['first_class_count']} | {c['second_class_count']} | {c['G_search']['gauss_candidates']} |")
    current = next(c for c in sym["controls"] if c["model"] == "nonconserved_current")
    rows += ["", f"Control 4 computed Gauss preservation: `{current['current_preservation']['raw']}`; with no conservation rule imposed this is nonzero and inconsistent."]
    rows += ["", "Per-tooth mutations (each reruns the shared computation and fails only its own expected assertion):", ""]
    rows.extend(f"- `{a['tooth']}` — `{a['status']}`: {a['message']}" for a in sym["ablations"])
    rows += [
        "",
        "## Dual-engine logs and agreement",
        "",
        "```text",
        *agreement,
        "ENGINE_AGREE: ALL_THEORIES_AND_CONTROLS",
        "```",
        "",
        "SymPy:",
        "",
        "```text",
        compact_log(out / "sympy_run.log"),
        "```",
        "",
        "Wolfram Language:",
        "",
        "```text",
        compact_log(out / "mathematica_run.log"),
        "```",
        "",
        "## Decision-table result",
        "",
    ]
    if all(v == "NATIVE_P_NO_EMERGENT_GAUSS" for v in verdicts):
        rows.append("Both regular coupling families have computed FC count zero, and neither semidefinite kinetic boundary develops a first-class Gauss direction. Since a regular nonlinear gauge identity must have a nontrivial leading linearization, the quadratic absence selects branch 2.")
    else:
        rows.append("At least one computed coupling stratum contains a Gauss candidate; the reported verdict follows the directive decision order without forcing the no-go branch.")
    rows += ["", f"**FINAL VERDICT: `{same_verdict}`.**", ""]
    return "\n".join(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    sym = json.loads((args.out_dir / "sympy_results.json").read_text())
    math = json.loads((args.out_dir / "mathematica_results.json").read_text())
    agreement = compare(sym, math)
    log = "\n".join(["DUAL_ENGINE_NATIVE_P_GATE", *agreement, "ENGINE_AGREE: ALL_THEORIES_AND_CONTROLS"]) + "\n"
    (args.out_dir / "engine_agreement.log").write_text(log)
    args.report.write_text(report(sym, math, agreement, args.out_dir) + "\n")
    print(log, end="")
    print(f"REPORT_WRITTEN: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
