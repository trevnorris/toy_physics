#!/usr/bin/env python3
"""Compare the independent Phase-A engines and write the report/results."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


HERE = Path(__file__).resolve().parent
ARTIFACTS = HERE / "reports/u1_body_dynamics_artifacts"
DEFAULT_REPORT = HERE / "reports/u1_body_dynamics.md"
DEFAULT_RESULTS = HERE / "reports/u1_body_dynamics_results.yaml"
DEFAULT_FIXTURES = HERE / "u1_body_dynamics_fixtures.yaml"
ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
COMPARATOR_TEETH = ("ENGINE_CANONICAL", "ENGINE_DEPENDENCIES")


class CompareFailure(AssertionError):
    def __init__(self, tooth: str, detail: str):
        super().__init__(f"ASSERT_FAIL:{tooth}:{detail}")
        self.tooth = tooth


def require(test: Any, tooth: str, detail: str) -> None:
    if not bool(test):
        raise CompareFailure(tooth, detail)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def powers_key(powers: dict[str, Any]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted((str(k), int(v)) for k, v in powers.items()))


def combine_terms(rows: list[dict[str, Any]]) -> dict[tuple[tuple[str, int], ...], float]:
    output: dict[tuple[tuple[str, int], ...], float] = {}
    for row in rows:
        key = powers_key(row["powers"])
        output[key] = output.get(key, 0.0) + float(row["coefficient"])
    return {key: value for key, value in output.items() if abs(value) > 1e-11}


def compare_terms(a: list[dict[str, Any]], b: list[dict[str, Any]], label: str) -> None:
    aa, bb = combine_terms(a), combine_terms(b)
    require(set(aa) == set(bb), "ENGINE_DEPENDENCIES", f"{label}:monomials {set(aa) ^ set(bb)}")
    for key in aa:
        require(math.isclose(aa[key], bb[key], rel_tol=3e-10, abs_tol=3e-10),
                "ENGINE_CANONICAL", f"{label}:{key}:{aa[key]}!={bb[key]}")


def numeric_string(value: str) -> float:
    text = value.replace("`", "")
    try:
        return float(sp.N(sp.sympify(text)))
    except Exception as exc:  # pragma: no cover - comparator failure path
        raise CompareFailure("ENGINE_CANONICAL", f"not numeric: {value}") from exc


def verdict_class(verdict: str) -> str:
    return verdict.split("(", 1)[0]


def expression_dependencies_from_terms(rows: list[dict[str, Any]]) -> set[str]:
    # Dependencies are extracted from the canonical monomials themselves.
    # This is deliberately independent of either engine's dependency summary.
    return {name for row in rows for name, power in row["powers"].items() if int(power) != 0}


def compare(sym: dict[str, Any], math_artifact: dict[str, Any]) -> list[str]:
    checks: list[str] = []
    require(sym["schema_version"] == math_artifact["schema_version"] == "U1_PHASE_A_FIELD_INPUTS_V2",
            "ENGINE_CANONICAL", "schema")
    require(sym["axis1"] == math_artifact["axis1"] == "COMPUTATION_VALID", "ENGINE_CANONICAL", "axis1")
    require(sym["axis2"] == math_artifact["axis2"], "ENGINE_CANONICAL", "axis2")
    checks.append("axis verdicts derived identically")

    sym_source = {(row["id"], row["source_file"], row["source_fragment"])
                  for row in sym["source_action_completeness"]["matched_assembled_terms"]}
    math_source = {(row["id"], row["source_file"], row["source_fragment"])
                   for row in math_artifact["source_action_completeness"]["matched_assembled_terms"]}
    require(sym_source == math_source, "ENGINE_DEPENDENCIES", "source-derived action surface")
    require(len(sym["source_action_completeness"]["source_derived_polar_monomials"]) == 3 ==
            len(math_artifact["source_action_completeness"]["source_derived_polar_monomials"]),
            "ENGINE_CANONICAL", "polar source manifest")
    checks.append("stage-note action manifests and assembled coverage")

    s_action, m_action = sym["assembled_action"], math_artifact["assembled_action"]
    require(set(s_action["term_expressions"]) == set(m_action["term_expressions"]),
            "ENGINE_DEPENDENCIES", "parsed action term set")
    require({key: set(value) for key, value in s_action["term_dependencies"].items()} ==
            {key: set(value) for key, value in m_action["term_dependencies"].items()},
            "ENGINE_DEPENDENCIES", "parsed action free-symbol dependencies")
    require(len(sym["source_action_completeness"]["source_derived_mandatory_g0_records"]) ==
            len(math_artifact["source_action_completeness"]["source_derived_mandatory_g0_records"]) == 3,
            "ENGINE_CANONICAL", "R1 immutable source manifest")
    for term_id in ("quantum_pressure", "Pu_coupling", "brane_shear_gradient"):
        sr = sym["action_term_removal_probes"][term_id]
        mr = math_artifact["action_term_removal_probes"][term_id]
        require(sr["operator_entry"] == mr["operator_entry"] and sr["operator_changed"] and mr["operator_changed"],
                "ENGINE_DEPENDENCIES", f"action removal reaches operator {term_id}")
    checks.append("parsed action expressions and three remove/re-derive source probes")

    sdim = {row["expression"]: tuple(row["computed_dimensions_LTM"]) for row in sym["dimensional_firewall"]}
    mdim = {row["expression"]: tuple(row["computed_dimensions_LTM"]) for row in math_artifact["dimensional_firewall"]}
    require(sdim == mdim, "ENGINE_CANONICAL", "inline dimensions")
    checks.append(f"{len(sdim)} inline dimensional constructions")

    require(sym["co_moving_laws"]["continuity_residual"] == math_artifact["co_moving_laws"]["continuity_residual"] == "0",
            "ENGINE_CANONICAL", "continuity change of variables")
    require(sym["co_moving_laws"]["momentum_residual"] == math_artifact["co_moving_laws"]["momentum_residual"] == "0",
            "ENGINE_CANONICAL", "momentum change of variables")
    require(sym["co_moving_laws"]["scalar_chain_rule_residual"] ==
            math_artifact["co_moving_laws"]["scalar_chain_rule_residual"] == "0",
            "ENGINE_CANONICAL", "composite-field chain rule")
    require(sym["co_moving_laws"]["declared_flux_residual"] ==
            math_artifact["co_moving_laws"]["declared_flux_residual"] == "0",
            "ENGINE_CANONICAL", "declared co-moving momentum flux")
    checks.append("continuity and momentum co-moving reductions")

    stails = {row["id"]: row for row in sym["tail_channels"]}
    mtails = {row["id"]: row for row in math_artifact["tail_channels"]}
    require(set(stails) == set(mtails), "ENGINE_DEPENDENCIES", "tail channel set")
    for name in stails:
        srow, mrow = stails[name], mtails[name]
        for key in ("classification", "decay_exponent", "gap_squared", "normalizable"):
            require(srow[key] == mrow[key], "ENGINE_CANONICAL", f"tail {name}:{key}")
        require(srow["solution_residual"] == mrow["solution_residual"] == "0",
                "ENGINE_CANONICAL", f"tail residual {name}")
        if srow["zero_mode_norm_value"] is not None:
            require(math.isclose(float(srow["zero_mode_norm_value"]), float(mrow["zero_mode_norm_value"]),
                                 rel_tol=3e-10, abs_tol=3e-10),
                    "ENGINE_CANONICAL", f"tail norm {name}")
        require(set(srow["dependencies"]) == set(mrow["dependencies"]),
                "ENGINE_DEPENDENCIES", f"tail deps {name}")
    checks.append(f"{len(stails)} solved radial channels, residuals, and translated norms")

    sop, mop = sym["linearized_channel_operator"], math_artifact["linearized_channel_operator"]
    require(set(sop["entries"]) == set(mop["entries"]), "ENGINE_DEPENDENCIES", "operator entry set")
    for key in ("density_gradient", "density_EOS_curvature", "phase_gradient", "wall_gradient",
                "wall_well_curvature", "polar_tangent_gradient", "polar_radial_curvature",
                "shear_curl", "Pu_cross", "uw_curvature"):
        require(math.isclose(numeric_string(sop["entries"][key]), numeric_string(mop["entries"][key]),
                             rel_tol=1e-12, abs_tol=1e-12), "ENGINE_CANONICAL", f"operator entry {key}")
    sc, mc = sym["coupled_indicial_analysis"], math_artifact["coupled_indicial_analysis"]
    require(sc["Pu"]["classification"] == mc["Pu"]["classification"] and
            sc["changes_scalar_channel_verdict"] == mc["changes_scalar_channel_verdict"],
            "ENGINE_CANONICAL", "coupled P-u class")
    require(math.isclose(numeric_string(sc["Pu"]["witness_k"]), numeric_string(mc["Pu"]["witness_k"]),
                         rel_tol=1e-12, abs_tol=1e-12) and
            math.isclose(numeric_string(sc["Pu"]["witness_determinant"]), numeric_string(mc["Pu"]["witness_determinant"]),
                         rel_tol=1e-12, abs_tol=1e-12), "ENGINE_CANONICAL", "coupled P-u witness")
    require(sc["density_phase"]["degree_difference"] == mc["density_phase"]["degree_difference"] == "-2",
            "ENGINE_CANONICAL", "drain density-phase degree analysis")
    checks.append("action-derived operator Hessians and coupled indicial/degree verdict")

    require(sym["phase_flux_normalization"]["normalization_residual"] ==
            math_artifact["phase_flux_normalization"]["normalization_residual"] == "0",
            "ENGINE_CANONICAL", "mdot flux normalization")
    require(sym["linearized_force_balance"]["substitution_residual"] ==
            math_artifact["linearized_force_balance"]["substitution_residual"] == "0",
            "ENGINE_CANONICAL", "D_E Z")
    require(sym["zero_mode_quotient"]["Q_times_Z"] == math_artifact["zero_mode_quotient"]["Q_times_Z"] == ["0", "0"],
            "ENGINE_CANONICAL", "projector")
    checks.append("force-balance operator, translated source block, Gram/projector algebra")

    for ep in ENDPOINTS:
        sr, mr = sym["endpoint_responses"][ep], math_artifact["endpoint_responses"][ep]
        require(sr["condition"] == mr["condition"] and sr["variational_class"] == mr["variational_class"],
                "ENGINE_CANONICAL", f"endpoint metadata {ep}")
        for component in ("normal", "tangent"):
            require(math.isclose(numeric_string(sr["fluid_coefficients"][component]),
                                 numeric_string(mr["fluid_coefficients"][component]), rel_tol=1e-12, abs_tol=1e-12),
                    "ENGINE_CANONICAL", f"endpoint response {ep}:{component}")
        require(math.isclose(numeric_string(sr["shear_coefficient"]), numeric_string(mr["shear_coefficient"]),
                             rel_tol=1e-12, abs_tol=1e-12), "ENGINE_CANONICAL", f"endpoint shear {ep}")
        require(sr["fluid_residual"] == mr["fluid_residual"] == ["0", "0"],
                "ENGINE_CANONICAL", f"endpoint residual {ep}")
        require(sr["shear_residual"] == mr["shear_residual"] == ["0"],
                "ENGINE_CANONICAL", f"endpoint shear residual {ep}")
        require(sr["solve_method"] == mr["solve_method"] == "exact rank-3 boundary/constraint solve",
                "ENGINE_CANONICAL", f"endpoint solve kind {ep}")
        require(sr["declared_matrix"] != [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]],
                "ENGINE_CANONICAL", f"identity endpoint solve {ep}")
        smom, mmom = sym["evaluated_moments"][ep], math_artifact["evaluated_moments"][ep]
        require(set(smom) == set(mmom), "ENGINE_DEPENDENCIES", f"moment set {ep}")
        for moment in smom:
            require(math.isclose(float(smom[moment]["production_value"]), float(mmom[moment]["production_value"]),
                                 rel_tol=3e-10, abs_tol=3e-10),
                    "ENGINE_CANONICAL", f"moment {ep}:{moment}")
            require(set(smom[moment]["dependencies"]) == set(mmom[moment]["dependencies"]),
                    "ENGINE_DEPENDENCIES", f"moment dependencies {ep}:{moment}")
        sa, ma = sym["endpoint_effective_actions"][ep], math_artifact["endpoint_effective_actions"][ep]
        compare_terms(sa["canonical_terms"], ma["canonical_terms"], f"{ep}:L_eff")
        for coefficient in sa["coefficients"]:
            sct = sa["coefficients"][coefficient]["canonical_terms"]
            mct = ma["coefficients"][coefficient]["canonical_terms"]
            compare_terms(sct, mct, f"{ep}:{coefficient}")
            extracted_s = expression_dependencies_from_terms(sct)
            extracted_m = expression_dependencies_from_terms(mct)
            require(extracted_s == set(sa["coefficients"][coefficient]["dependencies"]),
                    "ENGINE_DEPENDENCIES", f"SymPy declared deps {ep}:{coefficient}")
            require(extracted_m == set(ma["coefficients"][coefficient]["dependencies"]),
                    "ENGINE_DEPENDENCIES", f"Wolfram declared deps {ep}:{coefficient}")
    checks.append("five solved endpoint responses and all evaluated reduced moments")
    checks.append("canonical L_eff coefficient monomials and expression-derived dependency sets")

    require(sym["reconstruction_residuals"] == math_artifact["reconstruction_residuals"] == {ep: "0" for ep in ENDPOINTS},
            "ENGINE_CANONICAL", "reconstruction")
    require(sym["channel_assignments"] == math_artifact["channel_assignments"],
            "ENGINE_CANONICAL", "channel assignments")
    checks.append("substitution reconstruction and no-double-count channel partition")

    sanc = {(row["endpoint"], row["structure"], row["ancestor"], bool(row["open_remainder"]))
            for row in sym["per_structure_ancestry"]}
    manc = {(row["endpoint"], row["structure"], row["ancestor"], bool(row["open_remainder"]))
            for row in math_artifact["per_structure_ancestry"]}
    require(sanc == manc, "ENGINE_DEPENDENCIES", "per-structure ancestry")
    require(all(row["open_remainder"] or (row["classification_before"] == "NONZERO_NATIVE_STRUCTURE" and
                                          row["classification_after_removal"] == "ABSENT" and
                                          row["after_removal_monomials"] == [])
                for row in sym["per_structure_ancestry"]),
            "ENGINE_DEPENDENCIES", "actual native removal classification")
    snodes = {(row["id"], row["type"]) for row in sym["provenance_graph"]["nodes"]}
    mnodes = {(row["id"], row["type"]) for row in math_artifact["provenance_graph"]["nodes"]}
    require(snodes == mnodes and {tuple(edge) for edge in sym["provenance_graph"]["edges"]} ==
            {tuple(edge) for edge in math_artifact["provenance_graph"]["edges"]},
            "ENGINE_DEPENDENCIES", "expression-derived provenance traversal")
    checks.append(f"{len(sanc)} per-structure ancestry/native-padding ablations")

    require(sym["parity"] == math_artifact["parity"], "ENGINE_CANONICAL", "parity")
    require({name: verdict_class(value) for name, value in sym["verdict_reachability"].items()} ==
            {name: verdict_class(value) for name, value in math_artifact["verdict_reachability"].items()},
            "ENGINE_CANONICAL", "outcome reachability")
    checks.append("three-way parity tags and four data-driven outcome classes")
    return checks


def load_ablation_logs(artifacts: Path, teeth: list[str]) -> dict[str, dict[str, str]]:
    output: dict[str, dict[str, str]] = {}
    for tooth in teeth:
        output[tooth] = {}
        for engine, folder in (("SymPy", "sympy_ablation_logs"), ("Mathematica", "math_ablation_logs")):
            path = artifacts / folder / f"{tooth}.log"
            require(path.exists(), "ENGINE_CANONICAL", f"missing ablation log {path}")
            records = [line.split(":", 2) for line in path.read_text().splitlines() if line.startswith("ASSERT_FAIL:")]
            require(not any(len(record) >= 2 and record[1] == "MUTATION_NOOP" for record in records),
                    "ENGINE_CANONICAL", f"mutation noop sentinel {engine}:{tooth}")
            require(any(len(record) == 3 and record[0] == "ASSERT_FAIL" and record[1] == tooth for record in records),
                    "ENGINE_CANONICAL", f"wrong ablation {engine}:{tooth}")
            output[tooth][engine] = f"EXIT_1@ASSERT_FAIL:{tooth}"
    return output


def table(headers: list[str], rows: list[list[Any]]) -> list[str]:
    rendered = [[str(x).replace("|", "\\|").replace("\n", " ") for x in row] for row in rows]
    return ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"] + [
        "| " + " | ".join(row) + " |" for row in rendered]


def build_report(sym: dict[str, Any], agreement: list[str], ablations: dict[str, dict[str, str]], input_sha: str) -> str:
    lines = [
        "# U1 throat-body dynamics — Phase A remediation",
        "",
        f"Input SHA-256: `{input_sha}`. Phase A Axis 1 is **`{sym['axis1']}`** and Axis 2 is **`{sym['axis2']}`**.",
        "This report halts after spec 7.0/7.1. Phase B and Phase C are `NOT_RUN(upstream)` by process. The R1 coupled analysis changes the former scalar-only physics verdict, and that changed verdict is reported without repair toward the earlier target.",
        "",
        "## Frozen setup and honest scope",
        "",
        "The medium-rest lab frame, co-moving steady family, co-moving `Omega_c`, ambient-subtracted exterior-ball IR scheme, fixed `mdot`, and future `C_mdot` ownership of acceleration-like outer flux were declared before residual evaluation. The family is an action-derived exterior solution joined at `r=a` to field-level core traces through the typed throat surface functional. This is a collective-coordinate/effective-action computation, not a nonlinear throat simulation.",
        "",
        "No production input contains a tail exponent, eigenvalue sign, class boolean, or verdict. `Mh` and `cE` remain symbolic positive OPEN action coefficients; no collective-coordinate or multipole functional is an OPEN leaf.",
        "",
        "## 7.0 — declared-input ledger",
        "",
    ]
    ledger_rows = [[row["id"], row["status"], row["root_type"], row["domain"], tuple(row["dimensions"]),
                    ", ".join(row["arguments"]) or "none", row["symmetry_class"], row["source"]]
                   for row in sym["declared_inputs"]]
    lines += table(["root", "status", "type", "domain", "[L,T,M]", "arguments", "symmetry", "source"], ledger_rows)
    lines += ["", "### Action source completeness", "",
              "Both engines loaded the cited stage-note files, parsed every executable action expression, and parsed the T0/G0 action blocks. The assembled source cover contains the GNLS Berry/flow/EOS terms, Madelung quantum pressure, wall well/gradient/shear/anisotropy and typed throat/mix functionals, all three T0 polar terms, all MacCullagh/P–u/`u_w` terms, and the OPEN `h` kinetic/gradient normalization.", ""]
    source_rows = [[row["id"], row["source_file"], f"`{row['source_fragment']}`"]
                   for row in sym["source_action_completeness"]["matched_assembled_terms"]]
    lines += table(["assembled term", "loaded source", "parsed source fragment"], source_rows)
    lines += ["", "The mandatory pieces omitted by the first build are present. The immutable source manifest independently discovers quantum pressure, the MacCullagh shear gradient, and `L_Pu`; deleting any of the three changes its derived operator entry, and their joint SOURCE_COMPLETENESS ablation fails before classification.",
              "", "### Inline dimensions, IR, and G9 declaration", ""]
    lines += table(["constructed expression", "computed [L,T,M]"],
                   [[row["expression"], tuple(row["computed_dimensions_LTM"])] for row in sym["dimensional_firewall"]])
    tol = sym["ir_scheme"]["g9_tolerance"]
    lines += ["", f"Predeclared G9 force scale: `{tol['force_scale']}`; relative bound: `{tol['relative_terms']}` with `epsilon_rigid={tol['epsilon_rigid']}` and `epsilon_quad={tol['epsilon_quad']}`. It was not fitted to a measured residual.",
              "", "### Computed co-moving laws", "", "```text",
              sym["co_moving_laws"]["continuity_native"], sym["co_moving_laws"]["continuity_comoving"],
              sym["co_moving_laws"]["momentum_native"], sym["co_moving_laws"]["momentum_comoving"],
              sym["co_moving_laws"]["surface_force"], "```", "",
              f"The composite field `{sym['co_moving_laws']['composite_field']}` is differentiated before either balance law is reduced. The scalar chain-rule, continuity, momentum, and declared-flux residuals are all zero. The surface force is a native control-volume balance; no particle `F=dP/dt` form is imported.",
              "", "## 7.1 — far-field solve that decides G1", "",
              "For a bulk channel of radial dimension `d` and harmonic number `ell`, both engines construct and solve",
              "", "```text", "A_c [f''+(d-1)f'/r-ell(ell+d-2)f/r^2] - B_c f = 0.", "```", "",
              "`B_density=e''(rho_inf)` is differentiated in-engine from `K rho^5/4`, giving `5 K rho_inf^3`; `B_chi` and the radial-polar Hessian are likewise differentiated from their potentials. The translated norm uses the full measure: `S_(d-1) integral r^(d-1) |f'|^2 dr`. For the brane shear channel the four-dimensional measure factorizes as `d^3x g_ell(w)^2 dw`, so its reported radial dimension is three while the `w` integral is finite.", ""]
    tail_rows = []
    for row in sym["tail_channels"]:
        tail_rows.append([row["id"], row["radial_dimension"], row["equation"], row["classification"],
                          row["decay_exponent"] if row["decay_exponent"] is not None else row.get("gap", "algebraic"),
                          row["zero_mode_norm_integral"], f"{row['zero_mode_norm_value']:.12g}" if row["zero_mode_norm_value"] is not None else "n/a",
                          row["normalizable"]])
    lines += table(["channel", "d", "solved radial equation", "tail class", "nu or gap", "translated norm", "value at R=a=1", "normalizable"], tail_rows)
    coupled = sym["coupled_indicial_analysis"]
    lines += ["", "### R1 coupled operator result", "",
              "The displayed scalar channel rows remain the diagonal solves, but the verdict is taken from the full action-derived block. The drain background produces a density–phase mixed Hessian; its radial degree is two powers below the phase diagonal in `d=4`, so the computed indicial degree test leaves the scalar decay exponents unchanged.", "",
              "For the localized P–u term, each engine independently solves the bulk polar half-space profile and substitutes it back into the quadratic action. The resulting surface Hessian is", "", "```text",
              str(coupled["Pu"]["hessian"]),
              f"det = {coupled['Pu']['determinant']}",
              f"witness k = {coupled['Pu']['witness_k']}; det(witness) = {coupled['Pu']['witness_determinant']}", "```", "",
              f"Its computed class is **`{coupled['Pu']['classification']}`**. Thus the normalizable translational scalar tails still exist, but the declared homogeneous base has a negative long-wavelength coupled mode and Axis 2 changes honestly to **`{sym['axis2']}`**.", "",
              f"The bound-flow phase is selected from the solved continuity family by the `mdot` surface normalization; its computed normalization residual is `{sym['phase_flux_normalization']['normalization_residual']}`. No exponent is supplied in YAML.",
              "", "### Force-balance operator and quotient", "", "```text",
              sym["linearized_force_balance"]["operator_object"],
              "Z_A=-partial_A Phi_0 including the translated throat-source block",
              f"D_E Z residual = {sym['linearized_force_balance']['substitution_residual']}",
              f"Gram = {sym['zero_mode_quotient']['gram']}",
              f"Q Z = {sym['zero_mode_quotient']['Q_times_Z']}", "```", "",
              "The base source is constructed from the declared profile and the native force-balance equation; changing the field after constructing the source makes `BASE_BALANCE` fail. Translation is then checked by Cartesian substitution into `D_E`, not by a generator-scale identity.",
              "", "### Endpoint response solves", ""]
    endpoint_rows = []
    for ep in ENDPOINTS:
        row = sym["endpoint_responses"][ep]
        endpoint_rows.append([ep, row["condition"], row["variational_class"], row["fluid_coefficients"]["normal"],
                              row["fluid_coefficients"]["tangent"], row["shear_coefficient"]])
    lines += table(["endpoint", "field condition", "class", "normal response", "tangent response", "shear response"], endpoint_rows)
    lines += ["", "E1–E5 are five distinct exact rank-3 boundary/constraint solves rather than identity systems: no-slip, free-slip, permeable texture, nonholonomic shear-lock, and Robin/Rayleigh respectively. E4 additionally constructs `delta W=lambda_A(delta V_A-C_A[delta udot_T])`, whose allowed-variation residual is zero.",
              "", "### Evaluated reduced moments and L_eff", "",
              "All angular factors and radial integrals are executed. The following production values show the computed endpoint dependence; the machine result retains every symbolic field-trace expression and dependency set.", ""]
    moment_names = list(sym["evaluated_moments"]["E1"])
    lines += table(["moment"] + list(ENDPOINTS), [[name] + [f"{sym['evaluated_moments'][ep][name]['production_value']:.12g}" for ep in ENDPOINTS] for name in moment_names])
    lines += ["", "For each endpoint the directly embedded action reduces to", "", "```text",
              "L_eff = A_X V + A_p pdot + C_Xp p V + G_VV V^2/2 + G_Vp V pdot + G_pp pdot^2/2 - K_pp p^2/2.", "```", ""]
    for ep in ENDPOINTS:
        coeff = sym["endpoint_effective_actions"][ep]["coefficients"]
        lines.append(f"- `{ep}`: `G_VV={coeff['GVV']['expression']}`; `G_Vp={coeff['GVP']['expression']}`; `G_pp={coeff['GPP']['expression']}`; `K_pp={coeff['KPP']['expression']}`.")
    lines += ["", "Direct differentiation computes both canonical momenta and `Q_p`. Reconstruction independently iterates over the parsed action expressions, substitutes the rigid field embedding, performs the moment reductions, and compares that sum with the claimed coefficient form; the residual is zero for E1–E5. `K_pp` now includes the action-derived EOS Hessian, wall-well Hessian, and P–u cross term. E4 multiplier reaction and E5 Rayleigh loss remain outside `L_eff` in their typed channels.",
              "", "## Gates and able-to-fail evidence", ""]
    lines += table(["gate", "production result"], [[key, value] for key, value in sym["gates"].items()])
    lines += ["", "G1 passes because every diagonal translated tail is normalizable, while G2 passes because the displayed norms and reduced coefficients are finite in the predeclared ambient-subtracted exterior scheme. G3 is `CLASSIFIED_BY_AXIS2`: the computed P–u determinant supplies the negative internal-mode witness. G5 passes with body-only, symmetric-postulate, and one-sided-asymmetry tags kept distinct.",
              "", "### Per-tooth object ablations", ""]
    lines += table(["tooth", "SymPy", "Mathematica"], [[tooth, row["SymPy"], row["Mathematica"]] for tooth, row in ablations.items()])
    lines += ["", "Comparator mutations independently perturb one canonical coefficient and inject a symbol into one canonical monomial; they fail `ENGINE_CANONICAL` and `ENGINE_DEPENDENCIES` respectively.",
              "", "### Outcome reachability through field data", ""]
    lines += table(["fixture", "computed outcome"], [[key, value] for key, value in sym["verdict_reachability"].items()])
    lines += ["", "These controls use the same radial solver/classifier. The fat-tail fixture changes the spatial field domain, the unstable fixture changes the EOS action coefficient, and the unresolved fixture removes the positivity stratum of the OPEN `h` coefficient; none supplies an outcome token.",
              "", "## Dual-engine agreement", ""]
    lines += [f"- `{item}`" for item in agreement]
    lines += ["", "Agreement is on source-derived action coverage, solved ODE residuals, norm values, endpoint BVP coefficients, evaluated moments, canonical additive monomials, and dependency sets extracted from those monomials—not on copied verdict strings.",
              "", "## Provenance, ancestry, and parity", "", "```text",
              "ACTION + primitive field traces -> balanced exterior family -> D_E -> translated zero mode / quotient",
              "ACTION + primitive field traces + endpoint field BVP -> evaluated moments -> L_eff -> momenta / Q_var",
              "BALANCE + CONSTRAINT(E4) + RAYLEIGH(E5) + RETURN -> total generalized force (not L_eff)", "```", "",
              f"The machine graph contains `{len(sym['provenance_graph']['nodes'])}` typed nodes and `{len(sym['provenance_graph']['edges'])}` edges, all produced by traversal of the parsed quadratic and rigid-embedding expressions. Forbidden import injection is detected by node-type traversal. `{len(sym['per_structure_ancestry'])}` nonzero additive structures are actually re-derived with their ancestor removed: certified native structures change from `NONZERO_NATIVE_STRUCTURE` to `ABSENT`, while OPEN `h` remainders stay labeled and receive no native certification.",
              "", f"Parity tags: embedding `{sym['parity']['embedding_tag']}`, symmetric background `{sym['parity']['symmetric_tag']}`, one-sided control `{sym['parity']['asymmetric_control_tag']}`.",
              "", "## Verdict and Phase-A halt", "",
              f"- Axis 1: **`{sym['axis1']}`**.", f"- Axis 2: **`{sym['axis2']}`** in all E1–E5 Phase-A cells.",
              "- The result is conditional on the declared positive coefficient stratum and the field-level core/surface family; the controls show the classifier does not force OK.",
              "- Phase B: `NOT_RUN(upstream)`.", "- Phase C: `NOT_RUN(upstream)`.", "",
              "## Proposed parameter-register rows/edges (not applied)", "",
              "Proposed OPEN rows: `Mh`, `cE`, `mdot`, `gammaSigma`, `tangentDtN`, `sleeve_core_trace`, `geon_core_bundle`, `throat_surface_functional`, `outer_surface_functional`, `E4_shear_lock`, `E5_rayleigh`, and `return_closure`, with the domains/dimensions/arguments/symmetries in the ledger above.", "",
              "Proposed edges: `(hbar,m,rho_inf,K_EOS)->density D_E`; `(aB,kappaB)->chiB D_E`; `(m,rho_inf,K_EOS,a)->P D_E`; `(rhoBr,muR,ellg)->u D_E`; `(Mh,cE)->h D_E`; `core traces+surface functional->balanced exterior family`; `endpoint trace systems->N_**,U_**`; `ACTION moments->L_eff`; `E4_shear_lock->F_constraint`; `E5_rayleigh->F_Rayleigh`; `return_closure->F_flux`.", "",
              "**HALT: Phase A is complete; no Phase-B computation is included.**", ""]
    return "\n".join(lines)


def build_results(sym: dict[str, Any], agreement: list[str], ablations: dict[str, dict[str, str]], input_sha: str) -> dict[str, Any]:
    tail_classes: dict[str, int] = {}
    for row in sym["tail_channels"]:
        tail_classes[row["classification"]] = tail_classes.get(row["classification"], 0) + 1
    response_signatures = {
        (row["fluid_coefficients"]["normal"], row["fluid_coefficients"]["tangent"], row["shear_coefficient"])
        for row in sym["endpoint_responses"].values()
    }
    summary_lines = [
        f"U1 PHASE {sym['phase']} SUMMARY — HALT AFTER THIS PHASE",
        f"Axis 1: {sym['axis1']}",
        f"Axis 2: {sym['axis2']}",
        "Scalar tails: " + ", ".join(f"{key}={value}" for key, value in sorted(tail_classes.items())),
        f"R1 coupled P-u: {sym['coupled_indicial_analysis']['Pu']['classification']} at det={sym['coupled_indicial_analysis']['Pu']['witness_determinant']}",
        f"R1 drain cross: degree shift {sym['coupled_indicial_analysis']['density_phase']['degree_difference']} (subleading)",
        f"Endpoint response map: {len(response_signatures)}/{len(sym['endpoint_responses'])} distinct solved signatures",
        f"Dual engine: ENGINE_AGREE checks={len(agreement)}",
        f"Ablations: {len(ablations)} per engine + {len(COMPARATOR_TEETH)} comparator teeth",
        f"Downstream: Phase B {sym['downstream']['phase_B']}; Phase C {sym['downstream']['phase_C']}",
    ]
    return {
        "schema_version": "U1_BODY_DYNAMICS_RESULTS_V3", "phase": "A", "halt_after_phase": True,
        "input_sha256": input_sha, "axis_1": sym["axis1"], "axis_2": sym["axis2"],
        "cells": sym["cells"], "gates": sym["gates"], "base_configuration": sym["base_configuration"],
        "declared_inputs": sym["declared_inputs"], "source_action_completeness": sym["source_action_completeness"],
        "dimensional_firewall": sym["dimensional_firewall"], "ir_scheme": sym["ir_scheme"],
        "assembled_action": sym["assembled_action"], "linearized_channel_operator": sym["linearized_channel_operator"],
        "action_term_removal_probes": sym["action_term_removal_probes"],
        "coupled_indicial_analysis": sym["coupled_indicial_analysis"],
        "typed_co_moving_laws": sym["co_moving_laws"], "tail_channels": sym["tail_channels"],
        "phase_flux_normalization": sym["phase_flux_normalization"], "linearized_force_balance": sym["linearized_force_balance"],
        "zero_mode_quotient": sym["zero_mode_quotient"], "endpoint_responses": sym["endpoint_responses"],
        "evaluated_moments": sym["evaluated_moments"], "endpoint_effective_actions": sym["endpoint_effective_actions"],
        "reconstruction_residuals": sym["reconstruction_residuals"], "channel_assignments": sym["channel_assignments"],
        "E4_reduced_virtual_work": sym["E4_reduced_virtual_work"], "per_structure_ancestry": sym["per_structure_ancestry"],
        "provenance_graph": sym["provenance_graph"], "parity": sym["parity"],
        "outcome_reachability": sym["verdict_reachability"], "dual_engine": {"status": "ENGINE_AGREE", "checks": agreement},
        "ablations": ablations, "partition": sym["partition"], "downstream": sym["downstream"], "honest_scope": sym["honest_scope"],
        "summary_lines": summary_lines,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifacts", type=Path, default=ARTIFACTS)
    parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--mutation", choices=COMPARATOR_TEETH)
    parser.add_argument("--print-summary", action="store_true")
    args = parser.parse_args()
    if args.print_summary:
        result = yaml.safe_load(args.results.read_text())
        lines = result.get("summary_lines", [])
        if len(lines) != 10:
            print(f"ASSERT_FAIL:ENGINE_CANONICAL:summary line count {len(lines)}", file=sys.stderr)
            return 1
        print("\n".join(lines))
        return 0
    try:
        sym = load_json(args.artifacts / "sympy_phase_a.json")
        math_artifact = load_json(args.artifacts / "mathematica_phase_a.json")
        fixtures = yaml.safe_load(args.fixtures.read_text())
        attacks = {row["target"]: row["value"] for row in fixtures["mutations"].get(args.mutation, []) if row["op"] == "derived_attack"}
        if attacks.get("comparator_canonical_term") == "perturb_E1_GVV":
            math_artifact = copy.deepcopy(math_artifact)
            math_artifact["endpoint_effective_actions"]["E1"]["coefficients"]["GVV"]["canonical_terms"][0]["coefficient"] += 1.0
        if attacks.get("comparator_dependency_term") == "inject_forbidden_symbol":
            math_artifact = copy.deepcopy(math_artifact)
            math_artifact["endpoint_effective_actions"]["E1"]["coefficients"]["GVV"]["canonical_terms"][0]["powers"]["point_current"] = 1
        agreement = compare(sym, math_artifact)
        if args.mutation:
            print(f"ASSERT_FAIL:MUTATION_NOOP:{args.mutation}:mutation did not reach its own assert", file=sys.stderr)
            return 1
        ablations = load_ablation_logs(args.artifacts, sym["teeth"])
        input_path = HERE / "u1_body_dynamics_inputs.yaml"
        input_sha = hashlib.sha256(input_path.read_bytes()).hexdigest()
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(build_report(sym, agreement, ablations, input_sha) + "\n")
        args.results.write_text(yaml.safe_dump(build_results(sym, agreement, ablations, input_sha),
                                               sort_keys=False, allow_unicode=True, width=140))
        (args.artifacts / "engine_agreement.log").write_text("\n".join(agreement) + "\nENGINE_AGREE\n")
        print(f"COMPARATOR: ENGINE_AGREE checks={len(agreement)}")
        print(f"AXIS1: {sym['axis1']}")
        print(f"AXIS2: {sym['axis2']}")
        return 0
    except CompareFailure as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
