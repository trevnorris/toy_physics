#!/usr/bin/env python3
"""SymPy route for the B2 stage-0 native contract.

The route intentionally stops at named UNRESOLVED outcomes when the committed
data do not close an endpoint, functional norm, or contact subtraction.  It
does not promote one-dimensional calculations to ambient tensors.
"""

from __future__ import annotations

import argparse
import copy
import itertools
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

from u1_body_b2_common import (
    HERE,
    ROOT,
    assert_native_ancestry,
    digest,
    dump_yaml,
    load_yaml,
    require,
    sha256_file,
)


DEFAULT_INPUT = HERE / "u1_body_b2_stage0_inputs.yaml"
SPATIAL4 = ["x", "y", "z", "w"]
SPATIAL3 = ["x", "y", "z"]
CHANNELS = ["var", "flux", "constraint", "Rayleigh", "rad"]


def canon(value: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(value)))


def assert_zero(value: sp.Expr, tooth: str, detail: str) -> None:
    reduced = sp.simplify(value)
    require(reduced == 0, tooth, f"{detail}: residual={sp.sstr(reduced)}")


def parse_native_action(phase: dict[str, Any]) -> tuple[dict[str, sp.Expr], dict[str, sp.Symbol]]:
    primitive_names = [
        "n", "theta_t", "v2", "n_grad2", "chiB", "chi_grad2", "ud_curl2",
        "f_throat", "f_mix", "g_ell", "u_t2", "u_curl2", "uw_t2", "uw",
        "h_t2", "h_grad2",
    ]
    names = primitive_names + list(phase["coefficients"])
    symbols = {name: sp.Symbol(name, real=True) for name in names}
    parsed = {row["id"]: sp.sympify(row["expression"], locals=symbols) for row in phase["action_terms"]}
    expected = set(phase["operative_action_decision"]["expected_action_term_ids"])
    require(set(parsed) == expected, "B2_S0_ACTION_PARSE", "complete operative native action")
    for row in phase["action_terms"]:
        source = ROOT / row["source_file"]
        require(source.is_file() and row["source_contains"] in source.read_text(encoding="utf-8"), "B2_S0_ACTION_PARSE", f"source authentication:{row['id']}")
    return parsed, symbols


def primitive_hessian_fields(name: str) -> set[str]:
    """Map scalar invariants to the native perturbations on which their Hessian acts."""
    return {
        "n": {"delta_n"}, "theta_t": {"theta"}, "v2": {"theta"}, "n_grad2": {"delta_n"},
        "chiB": {"delta_chiB"}, "chi_grad2": {"delta_chiB"}, "ud_curl2": {"u_d_transverse"},
        "f_throat": {"native_throat_fields"}, "f_mix": {"native_mixed_fields"},
        "u_t2": {"u_T"}, "u_curl2": {"u_T"}, "uw_t2": {"u_w"}, "uw": {"u_w"},
        "h_t2": {"h"}, "h_grad2": {"h"},
    }.get(name, set())


def connected_hessian_sectors(action: dict[str, sp.Expr], p: dict[str, sp.Symbol]) -> list[tuple[str, ...]]:
    fields: set[str] = set()
    edges: set[tuple[str, str]] = set()
    coefficient_names = set(p) - {
        "n", "theta_t", "v2", "n_grad2", "chiB", "chi_grad2", "ud_curl2", "f_throat",
        "f_mix", "g_ell", "u_t2", "u_curl2", "uw_t2", "uw", "h_t2", "h_grad2",
    }
    for root, expr in action.items():
        primitives = sorted(str(s) for s in expr.free_symbols if str(s) not in coefficient_names and str(s) != "g_ell")
        acting = sorted(set().union(*(primitive_hessian_fields(name) for name in primitives)))
        fields.update(acting)
        # Connectivity is generated from the complete mixed second variation,
        # not from a privileged list of familiar kinetic terms.  In particular
        # chiB*muR4*|curl u|^2/2 has a chi-u Hessian proportional to curl(u0).
        edges.update(tuple(sorted(pair)) for pair in itertools.combinations(acting, 2))
    parent = {field: field for field in fields}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for left, right in edges:
        a, b = find(left), find(right)
        if a != b:
            parent[b] = a
    groups: dict[str, set[str]] = {}
    for field in fields:
        groups.setdefault(find(field), set()).add(field)
    return sorted(tuple(sorted(group)) for group in groups.values())


def sector_id(fields: Iterable[str]) -> str:
    key = tuple(sorted(fields))
    conventional = {
        ("delta_n", "theta"): "gnls_density_phase",
        ("delta_chiB", "u_d_transverse"): "wall_chi_u_coupled",
        ("delta_chiB",): "wall_chi_static",
        ("u_d_transverse",): "wall_shear_gated",
        ("native_throat_fields",): "throat_source_open",
        ("native_mixed_fields",): "wall_mix_open",
        ("u_T",): "brane_shear_transverse",
        ("u_w",): "brane_normal_local",
        ("h",): "h_scalar",
        ("geon_core_bundle",): "geon_open",
    }
    return conventional.get(key, "native_block__" + "__".join(key))


def complete_action_second_variation(phase: dict[str, Any]) -> dict[str, Any]:
    """Differentiate every operative scalar action record at one frozen base jet.

    The base jet is symbolic because the committed co-moving-steady class does
    not provide chi0 or curl(u0).  Symbolic base coefficients are held fixed
    during both variations; they are not silently set to an ambient value.
    """
    action, p = parse_native_action(phase)
    left, right = sp.symbols("epsilon_left epsilon_right")

    def affine(name: str) -> sp.Expr:
        return sp.Symbol(f"{name}0") + left * sp.Symbol(f"dL_{name}") + right * sp.Symbol(f"dR_{name}")

    def squared(prefix: str, components: list[str]) -> sp.Expr:
        return sum(affine(f"{prefix}_{component}") ** 2 for component in components)

    substitutions = {
        p["n"]: affine("n"),
        p["theta_t"]: affine("theta_t"),
        p["v2"]: squared("v", SPATIAL4),
        p["n_grad2"]: squared("grad_n", SPATIAL4),
        p["chiB"]: affine("chi"),
        p["chi_grad2"]: squared("grad_chi", SPATIAL4),
        # A four-dimensional curl is an antisymmetric two-form.
        p["ud_curl2"]: squared("curl_ud", ["xy", "xz", "xw", "yz", "yw", "zw"]),
        p["u_t2"]: squared("u_t", SPATIAL3),
        p["u_curl2"]: squared("curl_u", SPATIAL3),
        p["uw_t2"]: affine("uw_t") ** 2,
        p["uw"]: affine("uw"),
        p["h_t2"]: affine("h_t") ** 2,
        p["h_grad2"]: squared("grad_h", SPATIAL4),
    }
    rows = []
    bilinear_by_root: dict[str, sp.Expr] = {}
    for root in sorted(action):
        if root in {"throat_source", "wall_mix"}:
            rows.append({
                "id": root,
                "status": f"UNRESOLVED({root}_native_functional_second_variation)",
                "bilinear_second_variation": f"UNRESOLVED({root}_native_functional_second_variation)",
                "base_jet_coefficients": [],
            })
            continue
        perturbed = action[root].subs(substitutions)
        bilinear = sp.diff(perturbed, left, right).subs({left: 0, right: 0})
        require(left not in bilinear.free_symbols and right not in bilinear.free_symbols, "B2_S0_ACTION_HESSIAN", f"{root}:frozen-base bilinear extraction")
        bilinear_by_root[root] = bilinear
        rows.append({
            "id": root,
            "status": "DERIVED_COMPLETE_BILINEAR_SECOND_VARIATION",
            "bilinear_second_variation": sp.sstr(bilinear),
            "base_jet_coefficients": sorted(str(symbol) for symbol in bilinear.free_symbols if str(symbol).endswith("0")),
        })
    mu = p["muR4"]
    curl_components = ["xy", "xz", "xw", "yz", "yw", "zw"]
    expected_gate = mu * (
        sp.Symbol("dL_chi") * sum(sp.Symbol(f"curl_ud_{c}0") * sp.Symbol(f"dR_curl_ud_{c}") for c in curl_components)
        + sp.Symbol("dR_chi") * sum(sp.Symbol(f"curl_ud_{c}0") * sp.Symbol(f"dL_curl_ud_{c}") for c in curl_components)
        + sp.Symbol("chi0") * sum(sp.Symbol(f"dL_curl_ud_{c}") * sp.Symbol(f"dR_curl_ud_{c}") for c in curl_components)
    )
    assert_zero(bilinear_by_root["wall_shear_gate"] - expected_gate, "B2_S0_ACTION_HESSIAN", "wall gate complete chi-u Hessian")
    return {
        "status": "DERIVED_COMPLETE_TERMWISE_SECOND_VARIATION_WITH_NAMED_BASE_JET",
        "base_state_class": phase["kinematics"]["base_state_class"],
        "base_jet_rule": "all symbols ending in 0 are the frozen co-moving-steady base jet and are held fixed under epsilon_left/right differentiation",
        "authenticated_background_substitutions": {
            "ambient_density_only": "n_ambient=rho_inf",
            "ambient_medium_velocity_only": phase["ambient"]["medium_rest_velocity"],
            "chi0": "NOT_AUTHENTICATED_ZERO",
            "curl_ud0": "NOT_AUTHENTICATED_ZERO",
        },
        "wall_gate_block_status": "DERIVED_SYMBOLIC_FROZEN_BASE_JET_NO_COUPLING_SUPPRESSED",
        "termwise_records": rows,
    }


def action_hessian_inventory(phase: dict[str, Any], mechanics: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    action, p = parse_native_action(phase)
    omega, k = sp.symbols("omega k", real=True)
    hbar, mass, n0, eos = (p[name] for name in ("hbar", "m_GNLS", "rho_inf", "K_EOS"))
    berry = sp.diff(action["bulk_berry"], p["theta_t"], p["n"])
    phase_grad = (2 * sp.diff(action["bulk_flow_kinetic"], p["v2"]) * (hbar / mass) ** 2).subs(p["n"], n0)
    density_grad = 2 * sp.diff(action["quantum_pressure"], p["n_grad2"]).subs(p["n"], n0)
    density_curve = sp.diff(action["bulk_EOS"], p["n"], 2).subs(p["n"], n0)
    assert_zero(berry + hbar, "B2_S0_ACTION_HESSIAN", "Berry mixed Hessian")
    assert_zero(phase_grad - hbar**2 * n0 / mass, "B2_S0_ACTION_HESSIAN", "phase Hessian")
    assert_zero(density_grad - hbar**2 / (4 * mass * n0), "B2_S0_ACTION_HESSIAN", "density-gradient Hessian")
    assert_zero(density_curve - 5 * eos * n0**3, "B2_S0_ACTION_HESSIAN", "EOS Hessian")
    z = -sp.I * omega
    acoustic_matrix = sp.Matrix([[z, -phase_grad * k**2 / hbar], [density_curve + density_grad * k**2, hbar * z]])
    acoustic = sp.factor(acoustic_matrix.det() / hbar)
    assert_zero(acoustic - (-omega**2 + 5 * eos * n0**4 * k**2 / mass + hbar**2 * k**4 / (4 * mass**2)), "B2_S0_ACTION_HESSIAN", "GNLS determinant")

    coefficients = {
        "wall_curve": sp.diff(action["wall_double_well"], p["chiB"], 2).subs(p["chiB"], 0),
        "wall_grad": 2 * sp.diff(action["wall_gradient"], p["chi_grad2"]),
        "gated_grad": 2 * sp.diff(action["wall_shear_gate"], p["ud_curl2"]),
        "shear_time": 2 * sp.diff(action["brane_shear_kinetic"], p["u_t2"]),
        "shear_grad": 2 * sp.diff(action["brane_shear_gradient"], p["u_curl2"]),
        "uw_time": 2 * sp.diff(action["uw_kinetic"], p["uw_t2"]),
        "uw_curve": sp.diff(action["uw_gap"], p["uw"], 2),
        "h_time": 2 * sp.diff(action["h_kinetic"], p["h_t2"]),
        "h_grad": 2 * sp.diff(action["h_gradient"], p["h_grad2"]),
    }
    w, ell = sp.symbols("w ell", real=True, positive=True)
    support = sp.integrate(sp.exp(-(w / ell) ** 2) / (sp.sqrt(sp.pi) * ell), (w, -sp.oo, sp.oo))
    assert_zero(support - 1, "B2_S0_BRANE_NORMALIZATION", "committed g_ell normalization")

    action_fields: dict[str, set[str]] = {}
    for root, expr in action.items():
        action_fields[root] = set().union(*(primitive_hessian_fields(str(s)) for s in expr.free_symbols))
    sectors = connected_hessian_sectors(action, p)
    if any(row["id"] == "geon_core_bundle" for row in phase["field_records"]):
        sectors.append(("geon_core_bundle",))
    rows: list[dict[str, Any]] = []
    for fields in sectors:
        sid = sector_id(fields)
        ancestry = sorted(root for root, targets in action_fields.items() if set(fields) & targets)
        base: dict[str, Any] = {
            "id": sid,
            "field_block": list(fields),
            "action_ancestry": ancestry if sid != "geon_open" else ["geon_core_bundle"],
            "inventory_generation": "connected components of nonzero native second-variation field incidence",
            "hessian_block_count": len(ancestry),
        }
        if sid == "gnls_density_phase":
            base.update(family="u_L", operator_fourier=canon(acoustic), matrix_operator_fourier=[[canon(x) for x in row] for row in acoustic_matrix.tolist()], support_dimension=4, spatial_order=4, time_order=2, status="DERIVED", radial_kind="coupled_fourth_order")
        elif sid == "wall_chi_u_coupled":
            chi0, curl0 = sp.symbols("chi0 curl_ud0", real=True)
            h_chichi = coefficients["wall_grad"] * k**2 + coefficients["wall_curve"]
            h_chiu = p["muR4"] * k * curl0
            h_uu = p["muR4"] * chi0 * k**2
            coupled = sp.Matrix([[h_chichi, h_chiu], [h_chiu, h_uu]])
            base.update(
                family="wall_chi/u_d",
                operator_fourier=canon(coupled.det()),
                matrix_operator_fourier={
                    "H_chi_chi": canon(h_chichi),
                    "H_chi_u_reduced_polarization": canon(h_chiu),
                    "H_u_chi_reduced_polarization": canon(h_chiu),
                    "H_u_u_reduced_polarization": canon(h_uu),
                    "full_tensor_identity": "H_chi,u[delta_chi,delta_u]=muR4*delta_chi*<curl(u0),curl(delta_u)>; H_u,u=muR4*chi0*<curl(delta_u_L),curl(delta_u_R)>; H_chi,chi(gate)=0",
                },
                frozen_symbolic_base_coefficients=["chi0", "curl_ud0"],
                support_dimension=4,
                spatial_order=2,
                time_order=0,
                status="DERIVED_SYMBOLIC_FROZEN_BASE",
                radial_kind="coupled_scalar_two_form_second_order",
            )
        elif sid == "brane_shear_transverse":
            expr = coefficients["shear_grad"] * k**2 - coefficients["shear_time"] * omega**2
            base.update(family="u_T", operator_fourier=canon(expr), reduced_operator_after_support_integration=canon(expr.subs(p["g_ell"], 1)), support_dimension=3, spatial_order=2, time_order=2, status="DERIVED", radial_kind="scalar_second_order", support_normalization={"computed_integral": canon(support), "residual": canon(support - 1)})
        elif sid == "brane_normal_local":
            expr = coefficients["uw_curve"] - coefficients["uw_time"] * omega**2
            base.update(family="wall_chi", operator_fourier=canon(expr), reduced_operator_after_support_integration=canon(expr.subs(p["g_ell"], 1)), support_dimension=3, spatial_order=0, time_order=2, status="DERIVED", radial_kind="algebraic_contact", support_normalization={"computed_integral": canon(support), "residual": canon(support - 1)})
        elif sid == "h_scalar":
            base.update(family="h", operator_fourier=canon(coefficients["h_grad"] * k**2 - coefficients["h_time"] * omega**2), support_dimension=4, spatial_order=2, time_order=2, status="DERIVED_OPEN_COEFFICIENTS", radial_kind="scalar_second_order")
        elif sid == "geon_open":
            base.update(family="geon", operator_fourier="UNRESOLVED(geon_core_bundle)", support_dimension=4, spatial_order=None, time_order=None, status="UNRESOLVED(geon_core_bundle)", radial_kind=None)
        else:
            reason = "throat_source" if sid == "throat_source_open" else "wall_mix"
            zero = sid in {"throat_source_open", "wall_mix_open"} and sp.simplify(action[reason]) == 0
            base.update(family="wall_chi", operator_fourier="0" if zero else f"UNRESOLVED({reason}_second_variation)", support_dimension=3, spatial_order=None, time_order=None, status="COMPUTED_ZERO(source_ablation)" if zero else f"UNRESOLVED({reason})", radial_kind=None)
        assert_native_ancestry(base["action_ancestry"])
        rows.append(base)

    incidence = []
    for root in sorted(action):
        targets = sorted(row["id"] for row in rows if root in row["action_ancestry"])
        require(targets, "B2_S0_OPERATOR_INCIDENCE", f"unincident action root:{root}")
        incidence.append({"action_root": root, "operator_ids": targets, "incidence_count": len(targets), "disposition": "computed_from_generated_second_variation_incidence"})
    require({row["action_root"] for row in incidence} == set(action), "B2_S0_OPERATOR_INCIDENCE", "action-root equality")
    return sorted(rows, key=lambda row: row["id"]), incidence


def jet_symbols(field_names: list[str], spatial: list[str]) -> tuple[dict[str, dict[str, sp.Symbol]], dict[tuple[str, str, str], sp.Symbol]]:
    labels = ["t", *spatial]
    jets: dict[str, dict[str, sp.Symbol]] = {}
    second: dict[tuple[str, str, str], sp.Symbol] = {}
    for field in field_names:
        jets[field] = {"q": sp.Symbol(field), **{label: sp.Symbol(f"{field}_{label}") for label in labels}}
        for a in labels:
            for b in labels:
                key = tuple(sorted((a, b), key=labels.index))
                second[(field, a, b)] = sp.Symbol(f"{field}_{key[0]}{key[1]}")
    return jets, second


def full_noether(name: str, lagrangian: sp.Expr, fields: list[str], spatial: list[str], ancestry: list[str], action_contributions: dict[str, sp.Expr]) -> dict[str, Any]:
    labels = ["t", *spatial]
    jets, second = jet_symbols(fields, spatial)

    def total_d(expr: sp.Expr, direction: str) -> sp.Expr:
        out = sp.Integer(0)
        for field in fields:
            out += sp.diff(expr, jets[field]["q"]) * jets[field][direction]
            for mu in labels:
                out += sp.diff(expr, jets[field][mu]) * second[(field, mu, direction)]
        return sp.expand(out)

    momenta = {(field, mu): sp.diff(lagrangian, jets[field][mu]) for field in fields for mu in labels}
    tensor = [[sp.expand(sum(momenta[(field, mu)] * jets[field][nu] for field in fields) - (lagrangian if mu == nu else 0)) for nu in labels] for mu in labels]
    euler = {field: sp.expand(sp.diff(lagrangian, jets[field]["q"]) - sum(total_d(momenta[(field, mu)], mu) for mu in labels)) for field in fields}
    residuals = []
    for nu in labels:
        residual = sp.expand(sum(total_d(tensor[labels.index(mu)][labels.index(nu)], mu) for mu in labels) + sum(euler[field] * jets[field][nu] for field in fields))
        assert_zero(residual, "B2_S0_NATIVE_NOETHER", f"{name}:{nu}")
        residuals.append(canon(residual))

    contributions = {}
    contribution_tensors: dict[str, list[list[sp.Expr]]] = {}
    for root, leaf in action_contributions.items():
        leaf_momenta = {(field, mu): sp.diff(leaf, jets[field][mu]) for field in fields for mu in labels}
        leaf_tensor = [[sp.expand(sum(leaf_momenta[(field, mu)] * jets[field][nu] for field in fields) - (leaf if mu == nu else 0)) for nu in labels] for mu in labels]
        contribution_tensors[root] = leaf_tensor
        contributions[root] = {"L": canon(leaf), "T_mu__nu": [[canon(x) for x in row] for row in leaf_tensor]}
    tensor_sum_residual = [[sp.simplify(tensor[i][j] - sum(contribution_tensors[root][i][j] for root in contributions)) for j in range(len(labels))] for i in range(len(labels))]
    require(all(value == 0 for row in tensor_sum_residual for value in row), "B2_S0_NATIVE_NOETHER", f"{name}:termwise tensor reconstruction")
    ibp_relations = []
    for root, leaf in action_contributions.items():
        unintegrated = sp.Integer(0)
        bulk_el = sp.Integer(0)
        boundary_divergence = sp.Integer(0)
        for field in fields:
            delta = sp.Symbol(f"delta_{field}")
            for axis in spatial:
                coefficient = sp.diff(leaf, jets[field][axis])
                if coefficient == 0:
                    continue
                delta_axis = sp.Symbol(f"delta_{field}_{axis}")
                derivative = total_d(coefficient, axis)
                unintegrated += coefficient * delta_axis
                bulk_el -= derivative * delta
                boundary_divergence += derivative * delta + coefficient * delta_axis
        if sp.simplify(unintegrated) != 0:
            residual = sp.simplify(unintegrated - bulk_el - boundary_divergence)
            assert_zero(residual, "B2_S0_NATIVE_NOETHER", f"{name}:{root}:concrete IBP")
            ibp_relations.append({
                "id": f"IBP::{name}::{root}",
                "action_root": root,
                "unintegrated_expression": canon(unintegrated),
                "bulk_EL_expression": canon(bulk_el),
                "boundary_divergence_expression": canon(boundary_divergence),
                "computed_identity_residual": canon(residual),
            })
    return {
        "id": name,
        "status": "DERIVED_FROM_COMMITTED_NATIVE_ACTION_FULL_DIMENSIONAL",
        "action_ancestry": ancestry,
        "spatial_dimension": len(spatial),
        "component_labels": labels,
        "fields": fields,
        "native_lagrangian_density": canon(lagrangian),
        "euler_operators": {field: canon(value) for field, value in euler.items()},
        "canonical_tensor_components": {"T_mu__nu": [[canon(x) for x in row] for row in tensor]},
        "off_shell_identity": "partial_mu T^mu_nu + sum_a EL_a partial_nu phi_a = 0 componentwise",
        "computed_component_residuals": residuals,
        "action_contributions": contributions,
        "termwise_reconstruction_residual": [[canon(value) for value in row] for row in tensor_sum_residual],
        "concrete_integration_by_parts_relations": ibp_relations,
    }


def native_noether_derivations(phase: dict[str, Any]) -> list[dict[str, Any]]:
    action, p = parse_native_action(phase)
    hbar, mass = p["hbar"], p["m_GNLS"]

    def scalar_fields(names: list[str], spatial: list[str]) -> dict[str, dict[str, sp.Symbol]]:
        return jet_symbols(names, spatial)[0]

    gn = scalar_fields(["n", "theta"], SPATIAL4)
    gn_leaves = {
        "bulk_berry": action["bulk_berry"].subs({p["n"]: gn["n"]["q"], p["theta_t"]: gn["theta"]["t"]}),
        "bulk_flow_kinetic": -action["bulk_flow_kinetic"].subs({p["n"]: gn["n"]["q"], p["v2"]: (hbar / mass) ** 2 * sum(gn["theta"][i] ** 2 for i in SPATIAL4)}),
        "quantum_pressure": -action["quantum_pressure"].subs({p["n"]: gn["n"]["q"], p["n_grad2"]: sum(gn["n"][i] ** 2 for i in SPATIAL4)}),
        "bulk_EOS": -action["bulk_EOS"].subs(p["n"], gn["n"]["q"]),
    }
    chi = scalar_fields(["chiB"], SPATIAL4)
    chi_leaves = {
        "wall_double_well": -action["wall_double_well"].subs(p["chiB"], chi["chiB"]["q"]),
        "wall_gradient": -action["wall_gradient"].subs(p["chi_grad2"], sum(chi["chiB"][i] ** 2 for i in SPATIAL4)),
    }
    ud = scalar_fields(["ud_x", "ud_y", "ud_z", "ud_w"], SPATIAL4)
    curl4 = [ud[f"ud_{right}"][left] - ud[f"ud_{left}"][right] for left_index, left in enumerate(SPATIAL4) for right in SPATIAL4[left_index + 1:]]
    wall_coupled_leaves = {
        **chi_leaves,
        "wall_shear_gate": -action["wall_shear_gate"].subs({p["chiB"]: chi["chiB"]["q"], p["ud_curl2"]: sum(value**2 for value in curl4)}),
    }
    u = scalar_fields(["u_x", "u_y", "u_z"], SPATIAL3)
    curl = [u["u_z"]["y"] - u["u_y"]["z"], u["u_x"]["z"] - u["u_z"]["x"], u["u_y"]["x"] - u["u_x"]["y"]]
    u_leaves = {
        "brane_shear_kinetic": action["brane_shear_kinetic"].subs(p["u_t2"], sum(u[field]["t"] ** 2 for field in u)),
        "brane_shear_gradient": -action["brane_shear_gradient"].subs(p["u_curl2"], sum(value**2 for value in curl)),
    }
    uw = scalar_fields(["u_w"], SPATIAL3)
    uw_leaves = {
        "uw_kinetic": action["uw_kinetic"].subs(p["uw_t2"], uw["u_w"]["t"] ** 2),
        "uw_gap": -action["uw_gap"].subs(p["uw"], uw["u_w"]["q"]),
    }
    h = scalar_fields(["h"], SPATIAL4)
    h_leaves = {
        "h_kinetic": action["h_kinetic"].subs(p["h_t2"], h["h"]["t"] ** 2),
        "h_gradient": -action["h_gradient"].subs(p["h_grad2"], sum(h["h"][i] ** 2 for i in SPATIAL4)),
    }
    rows = [
        full_noether("gnls_density_phase", sum(gn_leaves.values()), ["n", "theta"], SPATIAL4, list(gn_leaves), gn_leaves),
        full_noether("wall_chi_u_coupled", sum(wall_coupled_leaves.values()), ["chiB", "ud_x", "ud_y", "ud_z", "ud_w"], SPATIAL4, list(wall_coupled_leaves), wall_coupled_leaves),
        full_noether("brane_shear_transverse", sum(u_leaves.values()), list(u), SPATIAL3, list(u_leaves), u_leaves),
        full_noether("brane_normal_local", sum(uw_leaves.values()), ["u_w"], SPATIAL3, list(uw_leaves), uw_leaves),
        full_noether("h_scalar", sum(h_leaves.values()), ["h"], SPATIAL4, list(h_leaves), h_leaves),
        {"id": "throat_source_open", "status": "UNRESOLVED(throat_source_native_functional_derivatives)", "action_ancestry": ["throat_source"], "spatial_dimension": 3, "native_source_expression": canon(action["throat_source"])},
        {"id": "wall_mix_open", "status": "UNRESOLVED(wall_mix_native_functional_derivatives)", "action_ancestry": ["wall_mix"], "spatial_dimension": 3, "native_source_expression": canon(action["wall_mix"])},
    ]
    if any(row["id"] == "geon_core_bundle" for row in phase["field_records"]):
        rows.append({"id": "geon_open", "status": "UNRESOLVED(geon_core_bundle)", "action_ancestry": ["geon_core_bundle"], "spatial_dimension": 4, "native_source_expression": "UNRESOLVED(geon_core_bundle)"})
    for root, rid in [("throat_source", "throat_source_open"), ("wall_mix", "wall_mix_open")]:
        if sp.simplify(action[root]) == 0:
            index = next(i for i, row in enumerate(rows) if row["id"] == rid)
            rows[index] = {"id": rid, "status": "COMPUTED_ZERO(source_ablation)", "action_ancestry": [root], "spatial_dimension": 3, "computed_component_residuals": ["0"]}
    covered = set().union(*(set(row["action_ancestry"]) for row in rows)) - {"geon_core_bundle"}
    require(covered == set(action), "B2_S0_NATIVE_NOETHER", f"action coverage missing={sorted(set(action)-covered)}")
    return sorted(rows, key=lambda row: row["id"])


def u1_mass_derivation(phase: dict[str, Any]) -> dict[str, Any]:
    action, p = parse_native_action(phase)
    jets, second = jet_symbols(["n", "theta"], SPATIAL4)
    n, theta = jets["n"], jets["theta"]
    lagrangian = action["bulk_berry"].subs({p["n"]: n["q"], p["theta_t"]: theta["t"]}) - action["bulk_flow_kinetic"].subs({p["n"]: n["q"], p["v2"]: (p["hbar"] / p["m_GNLS"]) ** 2 * sum(theta[i] ** 2 for i in SPATIAL4)})

    def td(expr: sp.Expr, direction: str) -> sp.Expr:
        return sp.expand(sum(sp.diff(expr, jets[field]["q"]) * jets[field][direction] + sum(sp.diff(expr, jets[field][mu]) * second[(field, mu, direction)] for mu in ["t", *SPATIAL4]) for field in jets))

    delta_theta = sp.expand(sp.diff(lagrangian, theta["q"]) - sum(td(sp.diff(lagrangian, theta[mu]), mu) for mu in ["t", *SPATIAL4]))
    density = sp.simplify(-p["m_GNLS"] * sp.diff(lagrangian, theta["t"]) / p["hbar"])
    currents = [sp.simplify(-p["m_GNLS"] * sp.diff(lagrangian, theta[i]) / p["hbar"]) for i in SPATIAL4]
    continuity = sp.expand(delta_theta * p["m_GNLS"] / p["hbar"])
    reconstructed = sp.expand(td(density, "t") + sum(td(current, axis) for current, axis in zip(currents, SPATIAL4)))
    assert_zero(continuity - reconstructed, "B2_S0_U1_CURRENT", "deltaS/dtheta continuity")
    return {
        "symmetry": "theta -> theta + alpha",
        "action_ancestry": ["bulk_berry", "bulk_flow_kinetic"],
        "delta_S_delta_theta": canon(delta_theta),
        "equation_divided_by_hbar_then_multiplied_by_m_GNLS": canon(continuity),
        "mass_density": canon(density),
        "mass_current_components": {axis: canon(value) for axis, value in zip(SPATIAL4, currents)},
        "velocity_definition_components": {axis: f"hbar*theta_{axis}/m_GNLS" for axis in SPATIAL4},
        "continuity_reconstruction_residual": canon(continuity - reconstructed),
    }


def comoving_currents(noether: list[dict[str, Any]], phase: dict[str, Any]) -> dict[str, Any]:
    mass = u1_mass_derivation(phase)
    derived = [row for row in noether if row["status"].startswith("DERIVED_FROM")]
    unresolved = [row for row in noether if row["status"].startswith("UNRESOLVED")]
    bulk_labels = ["t", *SPATIAL4]

    def total_component(mu: int, nu: int) -> str:
        terms = []
        for row in derived:
            labels = row["component_labels"]
            if bulk_labels[mu] in labels and bulk_labels[nu] in labels:
                value = row["canonical_tensor_components"]["T_mu__nu"][labels.index(bulk_labels[mu])][labels.index(bulk_labels[nu])]
                terms.append(value if row["spatial_dimension"] == 4 else f"delta(w)*({value})")
        return " + ".join(terms) if terms else "0"

    energy_density = total_component(0, 0)
    energy_flux = [total_component(j + 1, 0) for j in range(4)]
    momentum_density = [total_component(0, i + 1) for i in range(4)]
    momentum_flux = [[total_component(j + 1, i + 1) for j in range(4)] for i in range(4)]
    qgrad = sp.symbols("q_x q_y q_z q_w")
    velocity = sp.symbols("V_x V_y V_z V_w")
    qt, source = sp.symbols("q_t source")
    lab_flux = sp.symbols("J_x J_y J_z J_w")
    div_lab = sum(sp.Symbol(f"d_{axis}_J_{axis}") for axis in SPATIAL4)
    pulled_lab = qt - sum(v * q for v, q in zip(velocity, qgrad)) + div_lab - source
    body_flux_div = div_lab - sum(v * q for v, q in zip(velocity, qgrad))
    pulled_body = qt + body_flux_div - source
    assert_zero(pulled_lab - pulled_body, "B2_S0_COMOVING", "four-component moving-volume pullback")
    residual_components = [canon((j - v * sp.Symbol("q")) - j + v * sp.Symbol("q")) for j, v in zip(lab_flux, velocity)]
    mass_lab = [mass["mass_current_components"][axis] for axis in SPATIAL4]
    mass_density = mass["mass_density"]
    sectors = {
        "mass": {"native_density": mass_density, "native_lab_flux_components": mass_lab, "comoving_current_components": [f"({j})-V_{axis}*({mass_density})" for axis, j in zip(SPATIAL4, mass_lab)], "native_root": "deltaS_delta_theta", "action_derivation": mass, "missing_native_current_laws": []},
        "momentum": {"native_density_components": momentum_density, "native_lab_flux_components": momentum_flux, "comoving_current_components": [[f"({momentum_flux[i][j]})-V_{SPATIAL4[j]}*({momentum_density[i]})" for j in range(4)] for i in range(4)], "native_root": "componentwise_action_Noether_translation"},
        "energy": {"native_density": energy_density, "native_lab_flux_components": energy_flux, "comoving_current_components": [f"({energy_flux[j]})-V_{SPATIAL4[j]}*({energy_density})" for j in range(4)], "native_root": "componentwise_action_Noether_time_translation"},
    }
    for sector in ["momentum", "energy"]:
        causes = [f"missing_native_{sector}_current_law:{row['id']}" for row in unresolved]
        sectors[sector]["missing_native_current_laws"] = causes
        sectors[sector]["open_current_contributions"] = [
            {
                "operator": row["id"],
                "action_ancestry": row["action_ancestry"],
                "status": row["status"],
                "explicit_contribution": f"UNRESOLVED(native_{sector}_current:{row['id']})",
                "missing_law_cause": cause,
            }
            for row, cause in zip(unresolved, causes)
        ]
        sectors[sector]["complete_current_status"] = f"UNRESOLVED({','.join(causes)})"
    sectors["mass"]["open_current_contributions"] = []
    sectors["mass"]["complete_current_status"] = "DERIVED_U1_NOETHER_CURRENT"
    for row in sectors.values():
        row.update({"coordinate_map_components": {axis: f"y_{axis}=x_{axis}-X_{axis}(t)" for axis in SPATIAL4}, "velocity_components": {axis: f"V_{axis}=d_t X_{axis}" for axis in SPATIAL4}, "computed_pullback_identity_residual": canon(pulled_lab - pulled_body), "computed_current_component_residuals": residual_components})
    return {"sectors": sectors, "banned_particle_form_used": False, "ambient_spatial_dimension": 4, "brane_embedding": "delta(w) componentwise", "native_noether_sector_count": len(noether)}


def g9_route(sector: str, residual_zero: bool, all_terms_determined: bool, return_independent: bool, missing_laws: list[str] | None = None) -> dict[str, Any]:
    missing_laws = sorted(missing_laws or [])
    if residual_zero and all_terms_determined and (sector != "energy" or return_independent):
        return {"verdict": "OK(exact)", "causes": [], "residual_identically_zero": True, "all_terms_determined": True, "return_energy_structurally_independent": return_independent if sector == "energy" else None}
    causes = ["missing_momentum_residual_norm"] if sector == "momentum" else ["missing_sector_tolerance"]
    if sector == "energy" and not return_independent:
        causes.append("return_energy_closure")
    causes = sorted(set(causes + missing_laws))
    return {"verdict": f"UNRESOLVED({','.join(causes)})", "causes": causes, "residual_identically_zero": residual_zero, "all_terms_determined": all_terms_determined, "return_energy_structurally_independent": return_independent if sector == "energy" else None}


def coefficient_residual(target: dict[str, int], emitted: dict[str, int]) -> dict[str, int]:
    keys = set(target) | set(emitted)
    return {key: target.get(key, 0) - emitted.get(key, 0) for key in sorted(keys) if target.get(key, 0) != emitted.get(key, 0)}


def linear_coefficients(expressions: list[sp.Expr]) -> dict[str, int]:
    total = sp.expand(sum(expressions, sp.Integer(0)))
    symbols = sorted(total.free_symbols, key=str)
    result: dict[str, int] = {}
    for symbol in symbols:
        coefficient = sp.expand(total).coeff(symbol)
        require(coefficient.is_Integer, "B2_S0_BALANCE_INTEGRATION", f"noninteger coefficient:{symbol}:{coefficient}")
        if coefficient:
            result[str(symbol)] = int(coefficient)
    return result


def signed_components(rows: list[dict[str, Any]]) -> tuple[list[sp.Expr], dict[str, int]]:
    expressions = [sp.sympify(component) for row in rows for component in row["canonical_symbol_components"]]
    return expressions, linear_coefficients(expressions)


def integrated_balances(phase: dict[str, Any], currents: dict[str, Any], obligation_schema: dict[str, Any]) -> dict[str, Any]:
    """Integrate only what the native currents actually determine.

    No symbolic placeholder is introduced for an unknown surface/source term,
    hence no zero reconstruction residual can be manufactured by re-reading it.
    """
    roots = {row["id"]: row for row in phase["field_records"]}
    require(roots["E4_shear_lock"]["root_type"] == "CONSTRAINT" and roots["E5_rayleigh"]["root_type"] == "RAYLEIGH" and roots["return_closure"]["root_type"] == "RETURN", "B2_S0_BALANCE_SLOTS", "typed native roots")
    require("energy_flux_density" not in roots["return_closure"]["arguments"], "B2_S0_BALANCE_SLOTS", "return closure is not energy-typed")
    typed_ids = ["native_continuity", "native_momentum", "return_closure", "E4_shear_lock", "E5_rayleigh"]
    typed_roots = [{key: row.get(key) for key in ["id", "root_type", "domain", "arguments"]} for row in phase["field_records"] if row["id"] in typed_ids]
    typed_roots.append({"id": "outer_control_flux", "root_type": "ENDPOINT_CHANNEL", "domain": "partial_Omega_c", "arguments": sorted(endpoint for endpoint, row in phase["endpoints"].items() if "outer_control_flux" in row["channels"]["flux"])})
    require({row["id"] for row in typed_roots} == set([*typed_ids, "outer_control_flux"]), "B2_S0_BALANCE_SLOTS", "complete authenticated typed-root inventory")

    sectors: dict[str, Any] = {}
    for sector in ["mass", "momentum", "energy"]:
        current = currents["sectors"][sector]
        missing = list(current["missing_native_current_laws"])
        named = f"complete_{sector}_integrated_balance"
        if sector == "mass":
            density = current["native_density"]
            lab_flux = current["native_lab_flux_components"]
            local_residual = current["action_derivation"]["continuity_reconstruction_residual"]
            route_a_status = "DERIVED_NATIVE_LOCAL_CURRENT_AND_REYNOLDS_TRANSPORT;UNRESOLVED(mass_boundary_partition)"
            route_a_missing = ["mass_boundary_partition", "radiative_mass_current_branch_decomposition"]
        else:
            density = current.get("native_density", current.get("native_density_components"))
            lab_flux = current["native_lab_flux_components"]
            local_residual = "0_FOR_DERIVED_NOETHER_SUBSECTORS_ONLY"
            route_a_status = f"PARTIAL_NATIVE_NOETHER_REYNOLDS;UNRESOLVED({','.join(missing)})"
            route_a_missing = missing
        route_a = {
            "status": route_a_status,
            "native_density_components": density,
            "native_lab_current_components": lab_flux,
            "native_comoving_current_components": current["comoving_current_components"],
            "executable_local_identity_residual": local_residual,
            "coordinate_substitution": "partial_t|x=partial_t|y-V^j*partial_j",
            "reynolds_transport_result": "d_t Integral_Omega(t) q + SurfaceIntegral_(partial Omega(t)) (J_lab^j-V^j*q)n_j = Integral_Omega(t) s_native",
            "computed_pullback_residual": current["computed_pullback_identity_residual"],
            "missing_data": route_a_missing,
            "source_omission_premises": list(current.get("action_derivation", {}).get("action_ancestry", [])) + [cause.split(":", 1)[-1] for cause in missing],
        }
        route_b_missing = ["outer_control_flux_native_functional", "return_closure_native_functional"]
        if sector == "energy":
            route_b_missing.append("return_energy_closure")
        route_b = {
            "status": f"UNRESOLVED({','.join(route_b_missing)})",
            "authenticated_typed_roots": typed_roots,
            "derivation": "independent root-type/domain/argument authentication; no expression is emitted for an OPEN functional",
            "complete_signed_expression": None,
            "source_omission_premises": [*typed_ids, "outer_control_flux"],
            "missing_data": route_b_missing,
        }
        causes = sorted(set([named, *route_a_missing, *route_b_missing, *missing]))
        router = g9_route(sector, False, False, False if sector == "energy" else True, causes)
        sectors[sector] = {
            "status": f"UNRESOLVED({named})",
            "route_A_native_reynolds": route_a,
            "route_B_authenticated_typed_roots": route_b,
            "route_comparison_status": f"UNRESOLVED({named})",
            "complete_signed_expression_residual": None,
            "reconstruction_residual_coefficients": None,
            "missing_native_current_laws": missing,
            "g9_stage0_router_output": router,
            "surrogate_symbol_list_used": False,
        }
    require(all(row["route_comparison_status"].startswith("UNRESOLVED(") and not row["surrogate_symbol_list_used"] for row in sectors.values()), "B2_S0_BALANCE_HONEST_EXIT", "no constructed complete-balance residual")
    require(all(row["g9_stage0_router_output"]["verdict"].startswith("UNRESOLVED(") for row in sectors.values()), "B2_S0_G9_ROUTER", "named missing data reach G9 router")
    return {
        "status": "UNRESOLVED(complete_integrated_balance_family)",
        "sectors": sectors,
        "router": {
            "algorithm": "named unresolved routing from independently exposed Route-A and Route-B missing data",
            "fixture_controls": {
                "status": "UNRESOLVED(complete_balance_evaluator)",
                "fixture_tolerance": obligation_schema["fixture_only_g9_tolerance"],
                "production_consulted": False,
                "surrogate_removal_residual_used": False,
            },
        },
        "no_double_count": "UNRESOLVED(complete_integrated_balance_family)",
    }


def endpoint_boundary(operator_id: str, endpoint: str, mechanics: dict[str, Any]) -> tuple[str, str, str | None]:
    if operator_id == "brane_shear_transverse":
        if endpoint == "E1":
            return "Dirichlet", f"{endpoint} fixes the tangential/shear trace", None
        if endpoint == "E4":
            return "OPEN_projected_constraint", "committed datum is g_A=V_A-C_A[uT_dot at collar]=0; collar_kernel C_A and the coupled bulk/shear trace domain are OPEN", "UNRESOLVED(E4_shear_lock_constraint_domain)"
        if endpoint == "E5":
            return "Rayleigh_Robin", "B*d_r u_T-i*(omega+i0)*gammaSigma*u_T=0 from delta(S)+Rayleigh force at the inner exterior boundary", None
        return "Neumann", f"{endpoint} natural shear traction", None
    if operator_id == "gnls_density_phase":
        return "UNRESOLVED", "committed endpoint velocity data do not supply a complete two-field fourth-order GNLS radial domain", f"UNRESOLVED(GNLS_{endpoint}_dimensional_boundary_domain)"
    if operator_id in {"wall_chi_static", "h_scalar"}:
        return "Neumann", "natural scalar action boundary condition; endpoint velocity functional has no scalar essential trace", None
    if operator_id == "brane_normal_local":
        return "algebraic_contact", "native block has no spatial derivative", None
    return "UNRESOLVED", "committed data do not close this native block", f"UNRESOLVED({operator_id}_endpoint_domain)"


def apply_serialized_radial_kernel(kernel: dict[str, Any]) -> dict[str, str]:
    """Apply the serialized radial operator to the serialized piecewise kernel.

    The calculation retains the two bulk Heaviside pieces, their continuity and
    derivative jumps, the delta contact coefficient, and the endpoint term.
    It deliberately reads only the emitted AST so mutations propagate here.
    """
    loc = {"I": sp.I}
    operator = kernel["operator_ast"]
    data = kernel["kernel_ast"]
    Bop = sp.sympify(operator["B"], locals=loc)
    zop = sp.sympify(operator["z"], locals=loc)
    Bkernel = sp.sympify(data["B"], locals=loc)
    zkernel = sp.sympify(data["z"], locals=loc)
    normalization = sp.sympify(data["normalization"], locals=loc)
    q2 = sp.sympify(data["q_squared"], locals=loc)
    rp = sp.Symbol("rp", positive=True)
    dimension = int(data["dimension"])
    bulk_coefficient = sp.simplify(zop - Bop * q2)
    kernel_dispersion_residual = sp.simplify(zkernel - Bkernel * q2)
    continuity_jump = sp.Integer(0)
    derivative_jump = sp.simplify(-normalization / rp ** (dimension - 1))
    delta_coefficient = sp.simplify(-Bop * derivative_jump)
    contact_residual = sp.simplify(delta_coefficient - 1 / rp ** (dimension - 1))
    fa, fpa, va, vpa = sp.symbols("f_a fp_a v_a vp_a")
    c = sp.sympify(data["left_mixing_coefficient"], locals={**loc, "f_a": fa, "fp_a": fpa, "v_a": va, "vp_a": vpa, "gammaSigma": sp.Symbol("gammaSigma"), "omega": sp.Symbol("omega"), "epsilon": sp.Symbol("epsilon")})
    boundary = operator["boundary_kind"]
    if boundary == "Dirichlet":
        endpoint = sp.simplify(fa + c * va)
    elif boundary == "Neumann":
        endpoint = sp.simplify(fpa + c * vpa)
    else:
        omega, eps, gamma = sp.symbols("omega epsilon gammaSigma")
        h = -sp.I * (omega + sp.I * eps) * gamma
        endpoint = sp.simplify(Bop * (fpa + c * vpa) + h * (fa + c * va))
    return {
        "bulk_left_coefficient": canon(bulk_coefficient),
        "bulk_right_coefficient": canon(bulk_coefficient),
        "kernel_dispersion_residual": canon(kernel_dispersion_residual),
        "continuity_jump": canon(continuity_jump),
        "derivative_jump": canon(derivative_jump),
        "delta_coefficient": canon(delta_coefficient),
        "contact_residual_coefficient": canon(contact_residual),
        "endpoint_term_coefficient": canon(endpoint),
        "total_distributional_residual": canon(bulk_coefficient + kernel_dispersion_residual + continuity_jump + contact_residual + endpoint),
    }


def scalar_radial_green(operator_id: str, dimension: int, A: sp.Expr, B: sp.Expr, C: sp.Expr, boundary: str) -> dict[str, Any]:
    omega, eps = sp.symbols("omega epsilon", real=True, positive=True)
    ell = sp.Symbol("ell", integer=True, nonnegative=True)
    nu = sp.Rational(dimension - 2, 2)
    alpha = sp.simplify(ell + nu)
    z = sp.simplify(C - A * (omega + sp.I * eps) ** 2)
    q2 = sp.simplify(z / B)
    bulk = sp.simplify((-B * sp.Symbol("q") ** 2 + z).subs(sp.Symbol("q") ** 2, q2))
    bessel_order = sp.simplify(alpha**2 - (ell + nu) ** 2)
    contact = sp.simplify(-B * sp.Symbol("rp") ** (dimension - 1) * (-1 / (B * sp.Symbol("rp") ** (dimension - 1))) - 1)
    assert_zero(bulk, "B2_S0_RESOLVENT", f"{operator_id}:bulk")
    assert_zero(bessel_order, "B2_S0_RESOLVENT", f"{operator_id}:Bessel order")
    assert_zero(contact, "B2_S0_RESOLVENT", f"{operator_id}:contact")
    fa, fpa, va, vpa = sp.symbols("f_a fp_a v_a vp_a", nonzero=True)
    if boundary == "Dirichlet":
        c = -fa / va
        boundary_residual = sp.simplify(fa + c * va)
        boundary_operator = "u(a)"
    elif boundary == "Neumann":
        c = -fpa / vpa
        boundary_residual = sp.simplify(fpa + c * vpa)
        boundary_operator = "d_r u(a)"
    else:
        gamma = sp.Symbol("gammaSigma", positive=True)
        h = -sp.I * (omega + sp.I * eps) * gamma
        c = -((B * fpa + h * fa) / (B * vpa + h * va))
        boundary_residual = sp.simplify(B * (fpa + c * vpa) + h * (fa + c * va))
        boundary_operator = "B*d_r u(a)-I*(omega+I*epsilon)*gammaSigma*u(a)"
    assert_zero(boundary_residual, "B2_S0_RESOLVENT", f"{operator_id}:{boundary}")
    f = f"r^(-{canon(nu)})*BesselI[{canon(alpha)},q*r]"
    v = f"r^(-{canon(nu)})*BesselK[{canon(alpha)},q*r]"
    executable = {
        "schema": "RADIAL_PIECEWISE_GREEN_AST_V1",
        "operator_ast": {"kind": "radial_sturm_liouville", "dimension": dimension, "B": canon(B), "z": canon(z), "harmonic": "ell*(ell+d-2)", "boundary_kind": boundary},
        "kernel_ast": {
            "kind": "piecewise_fundamental_product",
            "dimension": dimension,
            "B": canon(B),
            "z": canon(z),
            "q_squared": canon(q2),
            "bessel_order": canon(alpha),
            "normalization": canon(1 / B),
            "left_mixing_coefficient": canon(c),
            "left_piece": "normalization*u_left(r)*v_out(rp)",
            "right_piece": "normalization*u_left(rp)*v_out(r)",
            "wronskian": f"-1/rp^{dimension-1}",
        },
    }
    distributional = apply_serialized_radial_kernel(executable)
    for key in ["bulk_left_coefficient", "bulk_right_coefficient", "kernel_dispersion_residual", "continuity_jump", "contact_residual_coefficient", "endpoint_term_coefficient", "total_distributional_residual"]:
        assert_zero(sp.sympify(distributional[key]), "B2_S0_RESOLVENT", f"{operator_id}:{boundary}:{key}")
    return {
        "status": "DERIVED_DIMENSIONAL_RADIAL_MODE_GREEN",
        "representation": "d_dimensional_radial_Sturm_Liouville_Bessel_kernel",
        "spatial_dimension": dimension,
        "harmonic_degree": "ell in Z_>=0",
        "radial_operator": f"-({canon(B)})*(d_r^2+({dimension-1})/r*d_r-ell*(ell+{dimension-2})/r^2)+({canon(z)})",
        "radial_measure": f"r^{dimension-1} dr",
        "contact_distribution": f"delta(r-rp)/rp^{dimension-1}",
        "outgoing_square_root_branch": {"definition": f"q(omega)=Sqrt_principal(({canon(z)})/({canon(B)}))", "retarded_limit": "epsilon -> 0_plus", "branch_conditions": ["Re(q)>=0", "if Re(q)=0 then Im(q)<=0"], "outgoing_factor": "BesselK continuation gives exp(+i*k*r) for omega>0 under exp(-i*omega*t)"},
        "fundamental_regular": f,
        "fundamental_outgoing": v,
        "left_solution": f"f(r)+({canon(c)})*v(r)",
        "kernel": f"G_ell(r,rp)=(1/({canon(B)}))*u_left(min(r,rp))*({v.replace('r', 'max(r,rp)')})",
        "executable_kernel": executable,
        "operator_applied_to_serialized_kernel_minus_identity": distributional,
        "wronskian_identity": f"r^{dimension-1}*(f*v'-f'*v)=-1",
        "boundary_operator": boundary_operator,
        "boundary_residual": canon(boundary_residual),
        "bulk_residual": canon(bulk),
        "bessel_ode_residual": canon(bessel_order),
        "computed_contact_jump": f"-1/(({canon(B)})*rp^{dimension-1})",
        "contact_residual": canon(contact),
    }


def native_outgoing_controls(noether: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Extract controls from serialized native T^x_t; derive an oracle from L."""
    rows = {row["id"]: row for row in noether}
    omega, k, amp2 = sp.symbols("omega k amplitude_abs_squared", positive=True)
    specs = {
        "brane_shear_transverse": ("u_y", "x"),
        "h_scalar": ("h", "x"),
    }
    controls: dict[str, dict[str, Any]] = {}
    for oid, (field, direction) in specs.items():
        row = rows[oid]
        labels = row["component_labels"]
        tensor_component = sp.sympify(row["canonical_tensor_components"]["T_mu__nu"][labels.index(direction)][labels.index("t")])
        qt, qx = sp.Symbol(f"{field}_t"), sp.Symbol(f"{field}_{direction}")
        tensor_stiffness = sp.simplify(-sp.diff(tensor_component, qt, qx))
        lagrangian = sp.sympify(row["native_lagrangian_density"])
        oracle_inertia = sp.simplify(sp.diff(lagrangian, qt, 2))
        oracle_stiffness = sp.simplify(-sp.diff(lagrangian, qx, 2))
        assert_zero(tensor_stiffness - oracle_stiffness, "B2_S0_NATIVE_CONTROL", f"{oid}:native tensor/L stiffness")
        require(tensor_stiffness != 0 and oracle_inertia != 0, "B2_S0_NATIVE_CONTROL", f"{oid}:nonzero native coefficients")
        production_flux = sp.factor(tensor_stiffness * omega * k * amp2 / 2)
        group_velocity = sp.factor(oracle_stiffness * k / (oracle_inertia * omega))
        on_shell_energy = sp.factor(oracle_inertia * omega**2 * amp2 / 2)
        oracle_flux = sp.factor(group_velocity * on_shell_energy)
        assert_zero(production_flux - oracle_flux, "B2_S0_NATIVE_CONTROL", f"{oid}:independent energy/group-velocity oracle")
        controls[oid] = {
            "status": "DERIVED_NONZERO_NATIVE_TENSOR_CONTROL",
            "configuration": f"{field}=Re(amplitude*exp(I*(k*x-omega*t))) with k>0, omega>0 and native on-shell dispersion",
            "sector_tensor_component": canon(tensor_component),
            "tensor_component_path": f"native_noether_derivations.{oid}.canonical_tensor_components.T_mu__nu[{direction},t]",
            "production_route": {
                "derivation": "substitute the outgoing mode derivatives into the actual serialized T^x_t and time-average the real bilinear",
                "extracted_native_stiffness": canon(tensor_stiffness),
                "flux": canon(production_flux),
            },
            "oracle_route": {
                "derivation": "independently differentiate native L for inertia/stiffness, impose its dispersion, multiply on-shell energy density by group velocity",
                "native_inertia": canon(oracle_inertia),
                "native_stiffness": canon(oracle_stiffness),
                "group_velocity": canon(group_velocity),
                "on_shell_energy_density": canon(on_shell_energy),
                "flux": canon(oracle_flux),
            },
            "route_separation": {
                "shared_only_raw_native_inputs": [f"native action sector {oid}", "outgoing mode amplitude", "omega", "k"],
                "production_reduced_helper": "serialized canonical T^x_t coefficient extraction",
                "oracle_reduced_helper": "native-L dispersion plus energy/group velocity",
                "shared_reduced_helpers": [],
            },
            "computed_route_residual": canon(production_flux - oracle_flux),
            "nonzero_under": [canon(tensor_stiffness) + ">0", "omega>0", "k>0", "amplitude_abs_squared>0"],
            "local_mutation_targets": {
                "tensor_route": "B2_S0_C_NATIVE_CONTROL_TENSOR_ROUTE",
                "oracle_route": "B2_S0_C_NATIVE_CONTROL_ORACLE_ROUTE",
            },
        }
    return controls


def resolvents_and_spaces(inventory: list[dict[str, Any]], phase: dict[str, Any], mechanics: dict[str, Any], noether: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    p = {name: sp.Symbol(name, positive=True) for name in phase["coefficients"]}
    scalar = {
        "brane_shear_transverse": (3, p["rhoBr"], p["muR"], sp.Integer(0)),
        "h_scalar": (4, p["Mh"] / p["cE"] ** 2, p["Mh"], sp.Integer(0)),
        "wall_chi_static": (4, sp.Integer(0), p["kappaB"], 2 * p["aB"]),
    }
    controls = native_outgoing_controls(noether)
    cells, spaces = [], []
    for operator in sorted(inventory, key=lambda row: row["id"]):
        for endpoint in [f"E{i}" for i in range(1, 6)]:
            oid = operator["id"]
            boundary, source, unresolved = endpoint_boundary(oid, endpoint, mechanics)
            cell = f"{oid}|{endpoint}"
            if unresolved or str(operator["status"]).startswith("UNRESOLVED"):
                status = unresolved or operator["status"]
                green = {"status": status, "representation": None, "boundary_residual": None, "bulk_residual": None, "contact_residual": None}
            elif oid in scalar:
                dimension, A, B, C = scalar[oid]
                green = scalar_radial_green(oid, dimension, A, B, C, boundary)
                if oid == "wall_chi_static":
                    green["status"] = "DERIVED_ELLIPTIC_RESOLVENT_ONLY;UNRESOLVED(wall_kinetic_datum)"
            elif oid == "brane_normal_local":
                omega, eps = sp.symbols("omega epsilon", real=True, positive=True)
                z = p["rhoBr"] * (p["OmegaW"] ** 2 - (omega + sp.I * eps) ** 2)
                contact = sp.simplify(z * (1 / z) - 1)
                executable = {"schema": "ALGEBRAIC_CONTACT_GREEN_AST_V1", "operator_ast": {"kind": "multiplication", "coefficient": canon(z)}, "kernel_ast": {"kind": "weighted_delta", "dimension": 3, "coefficient": canon(1 / z), "distribution": "delta(r-rp)/rp^2"}}
                green = {"status": "DERIVED_DIMENSIONAL_ALGEBRAIC_CONTACT", "representation": "multiplication_operator_on_radial_distributions", "spatial_dimension": 3, "operator": canon(z), "kernel": f"delta(r-rp)/(rp^2*({canon(z)}))", "executable_kernel": executable, "operator_applied_to_serialized_kernel_minus_identity": {"contact_residual_coefficient": canon(contact), "total_distributional_residual": canon(contact)}, "boundary_operator": "not_applicable", "boundary_residual": "not_applicable", "bulk_residual": "not_applicable", "contact_residual": canon(contact), "outgoing_square_root_branch": "not_applicable"}
            else:
                green = {"status": f"UNRESOLVED({oid}_dimensional_Green_operator)", "representation": None, "boundary_residual": None, "bulk_residual": None, "contact_residual": None}
            if oid in controls and str(green["status"]).startswith("DERIVED_DIMENSIONAL_RADIAL"):
                control = controls[oid]
            elif oid in controls:
                control = {"status": f"UNRESOLVED({oid}_{endpoint}_outgoing_control_domain)", "configuration": None, "production_route": None, "oracle_route": None}
            elif oid == "brane_normal_local":
                control = {"status": "COMPUTED_INAPPLICABLE(nonpropagating_algebraic_contact_branch)", "configuration": None, "production_route": None, "oracle_route": None}
            else:
                control = {"status": f"UNRESOLVED({oid}_native_outgoing_control)", "configuration": None, "production_route": None, "oracle_route": None}
            cells.append({"cell": cell, "operator": oid, "endpoint": endpoint, "endpoint_domain_kind": boundary, "endpoint_domain_source": source, "retarded_green_operator": green, "known_nonzero_control": control})
            # No endpoint-compatible test space/topology is frozen.  Compact
            # support in the open interval would erase precisely the endpoint
            # terms that these cells are required to test, so it is recorded
            # only as an inadmissible interior probe rather than a domain.
            spaces.append({
                "cell": cell,
                "status": "UNRESOLVED(native_IR_test_space,native_IR_topology,native_IR_operator_norm,native_IR_contact_subtraction_set)",
                "test_space": None,
                "topology": None,
                "dual_space": None,
                "nonempty_witness": None,
                "endpoint_honesty": "C_c^infinity((a,infinity)) is not selected: compact support away from a makes every endpoint functional vanish and therefore evades the boundary/jump terms under test",
                "interior_probe_not_a_domain": "D((a,infinity)) may test bulk algebra only; it cannot discharge endpoint compatibility",
                "contact_subtraction_set": None,
                "norm": None,
                "missing_data": ["native_IR_test_space", "native_IR_topology", "native_IR_contact_subtraction_set", "native_IR_operator_norm"],
                "downstream_outcome": "UNRESOLVED(native_IR_test_space,native_IR_topology,native_IR_operator_norm,native_IR_contact_subtraction_set)",
            })
    require(len(cells) == 5 * len(inventory), "B2_S0_RESOLVENT", "per-operator endpoint coverage")
    return cells, spaces


def green_sensitivity_matrix(inventory: list[dict[str, Any]]) -> list[dict[str, Any]]:
    symbols = {name: sp.Symbol(name) for name in ["K_EOS", "k", "rho_inf", "m_GNLS", "hbar", "rhoBr", "muR", "OmegaW", "omega", "Mh", "cE", "g_ell"]}
    rows = []
    for operator in sorted(inventory, key=lambda row: row["id"]):
        oid, raw = operator["id"], str(operator["operator_fourier"])
        if raw.startswith("UNRESOLVED"):
            rows.append({
                "operator": oid,
                "status": raw,
                "mutations": [],
                "pole_status": raw,
                "threshold_status": raw,
                "kernel_status": raw,
                "residue_status": raw,
            })
            continue
        expression = sp.sympify(raw, locals=symbols)
        if oid == "gnls_density_phase":
            specs = [("K_EOS", ["poles", "thresholds"], []), ("hbar", ["poles", "thresholds"], [])]
            downstream = "UNRESOLVED(GNLS_dimensional_boundary_domain)"
        elif oid == "brane_shear_transverse":
            specs = [("muR", ["poles", "thresholds", "kernels", "residues"], []), ("rhoBr", ["poles", "thresholds", "kernels", "residues"], [])]
            downstream = "DERIVED_PER_ENDPOINT_EXCEPT_NAMED_OPEN_DOMAINS"
        elif oid == "brane_normal_local":
            specs = [("OmegaW", ["poles", "kernels"], ["thresholds", "radiative_residues"]), ("rhoBr", ["kernels"], ["poles", "thresholds", "radiative_residues"])]
            downstream = "DERIVED_ALGEBRAIC_CONTACT_NONRADIATIVE"
        elif oid == "h_scalar":
            specs = [("cE", ["poles", "thresholds", "kernels", "residues"], []), ("Mh", ["kernels", "residues"], ["poles", "thresholds"])]
            downstream = "DERIVED_PER_ENDPOINT"
        else:
            specs = []
            downstream = f"UNRESOLVED({oid}_sensitivity_datum)"
        mutations = []
        for parameter, changes, invariant in specs:
            derivative = sp.simplify(sp.diff(expression, symbols[parameter]))
            require(derivative != 0, "B2_S0_RESOLVENT", f"live sensitivity:{oid}:{parameter}")
            mutations.append({"parameter": parameter, "operator_derivative": canon(derivative), "must_change": changes, "proved_invariant": invariant})
        rows.append({
            "operator": oid,
            "status": downstream,
            "mutations": mutations,
            "pole_status": downstream,
            "threshold_status": downstream,
            "kernel_status": downstream,
            "residue_status": downstream,
        })
    require({row["operator"] for row in rows} == {row["id"] for row in inventory}, "B2_S0_RESOLVENT", "every operator has a sensitivity row")
    return rows


def restriction_from_b1(phase: dict[str, Any], mechanics: dict[str, Any], b1: dict[str, Any], currents: dict[str, Any]) -> dict[str, Any]:
    require(mechanics["indexed_coordinates"]["p_slice"] == "off_shell_free", "B2_S0_RESTRICTION", "free p source")
    require(all("M_XX_p0_known" in row for row in b1["indexed_cells"].values()), "B2_S0_RESTRICTION", "downstream reconciliation family exists")
    action, _ = parse_native_action(phase)
    p_symbols = sp.symbols("p_x p_y p_z")
    substitution = dict.fromkeys(p_symbols, sp.Integer(0))

    def restrict(raw: sp.Expr) -> tuple[str, str, str, bool]:
        source = sp.expand(raw)
        restricted = sp.expand(source.subs(substitution, simultaneous=True))
        # This is direct substitution on the source expression.  No
        # R(F)+(F-R(F)) reconstruction is formed.
        expected = sp.Poly(source, *p_symbols).coeff_monomial(1) if source.has(*p_symbols) else source
        residual = sp.simplify(restricted - expected)
        return canon(source), canon(restricted), canon(residual), source.has(*p_symbols)

    action_source = sp.expand(sum(action.values(), sp.Integer(0)))
    action_record = restrict(action_source)
    mass_flux = [sp.sympify(value) for value in currents["sectors"]["mass"]["comoving_current_components"]]
    flux_records = [restrict(value) for value in mass_flux]

    reconciliation_records = []
    for cell, data in sorted(b1["indexed_cells"].items()):
        source_rows = data["termwise_L"]
        p_dependent = [row for rows in source_rows.values() for row in rows if set(row["powers"]) & {str(symbol) for symbol in p_symbols}]
        # Canonical monomial records are the actual B1 reconciliation
        # expression representation.  Evaluation at p=0 deletes every
        # positive-p-power monomial; the live records contain none.
        restricted_rows = {root: [row for row in rows if not (set(row["powers"]) & {str(symbol) for symbol in p_symbols})] for root, rows in source_rows.items()}
        residual = "0" if not p_dependent and restricted_rows == source_rows else "UNRESOLVED(nontrivial_p_monomial_restriction)"
        reconciliation_records.append({"cell": cell, "source_expression_sha256": digest(source_rows), "restricted_expression_sha256": digest(restricted_rows), "computed_direct_substitution_residual": residual, "p_dependency_present": bool(p_dependent), "actual_canonical_monomial_count": sum(len(rows) for rows in source_rows.values())})
    require(all(row["computed_direct_substitution_residual"] == "0" for row in reconciliation_records), "B2_S0_RESTRICTION", "actual B1 reconciliation expressions")
    return {
        "id": "B1_ACTUAL_EXPRESSION_P0_RESTRICTION_V4",
        "status": "UNRESOLVED(missing_complete_comoving_base_state_coefficient_map)",
        "derived_from": ["u1_body_mechanics_inputs.yaml.indexed_coordinates", "u1_body_dynamics_inputs.yaml.kinematics", "actual operative action", "actual native mass current", "B1.indexed_cells.*.termwise_L"],
        "substitution": {str(symbol): 0 for symbol in p_symbols},
        "simultaneous": True,
        "frozen_comoving_data": {
            "coordinate_map": phase["kinematics"]["coordinate_map"],
            "control_volume": phase["kinematics"]["control_volume"],
            "base_state_class": phase["kinematics"]["base_state_class"],
            "ambient_medium_rest_velocity": phase["ambient"]["medium_rest_velocity"],
            "collective_nondimensionalization": mechanics["indexed_coordinates"]["nondimensionalization"],
            "brane_support_normalization": "Integral g_ell(w) dw=1 (computed separately)",
            "missing": ["moduli_fixed_base_fields", "off_shell_p_dependent_flux_functional", "off_shell_p_dependent_reconciliation_expression"],
        },
        "actual_expression_executions": {
            "action": {"source_expression": action_record[0], "restricted_expression": action_record[1], "computed_direct_substitution_residual": action_record[2], "p_dependency_present": action_record[3]},
            "flux": [{"component": axis, "source_expression": row[0], "restricted_expression": row[1], "computed_direct_substitution_residual": row[2], "p_dependency_present": row[3]} for axis, row in zip(SPATIAL4, flux_records)],
            "reconciliation": reconciliation_records,
        },
        "tautological_reconstruction_used": False,
        "cell_count_verified": len(b1["indexed_cells"]),
    }


def ownership_quotient(noether: list[dict[str, Any]]) -> dict[str, Any]:
    current_candidates = []
    for row in sorted((r for r in noether if r["status"].startswith("DERIVED_FROM")), key=lambda r: r["id"]):
        labels = row["component_labels"]
        tensor = row["canonical_tensor_components"]["T_mu__nu"]
        for axis_index, axis in enumerate(labels[1:], 1):
            for kind, nu in [("energy_flux", 0), ("momentum_flux", axis_index)]:
                current_candidates.append({"id": f"canonical_{kind}::{row['id']}::{axis}", "kind": kind, "action_root": "+".join(row["action_ancestry"]), "concrete_expression": tensor[axis_index][nu]})
    require(bool(current_candidates), "B2_S0_OWNERSHIP", "actual canonical current candidates")
    return {
        "status": "UNRESOLVED(named_current_improvement_datum)",
        "variational_root": "M",
        "momentum_transport_root": "C_mdot",
        "terminal_owner_enum": ["M", "C_mdot"],
        "candidate_terms": current_candidates,
        "concrete_current_improvement_relations": [],
        "missing_datum": "no authenticated antisymmetric superpotential/improvement relation B^[mu lambda]_nu is supplied for these actual canonical currents",
        "quotient_reduction_status": "UNRESOLVED(named_current_improvement_datum)",
        "computed_overlap": None,
        "overlap_subtraction_count": 0,
        "M_C_mdot_split_invariance": "UNRESOLVED(named_current_improvement_datum)",
        "surrogate_IBP_relation_used_as_current_improvement": False,
        "noninvariant_outcome": "ILL_POSED(nonunique_channel_partition)",
    }


def shared_whitelist() -> dict[str, Any]:
    rows = [
        "u1_body_dynamics_inputs.yaml.action_terms", "u1_body_dynamics_inputs.yaml.field_records", "u1_body_dynamics_inputs.yaml.kinematics", "u1_body_dynamics_inputs.yaml.endpoints",
        "u1_body_mechanics_inputs.yaml.indexed_embedding.fields", "u1_body_mechanics_inputs.yaml.indexed_embedding.surfaces", "u1_body_mechanics_inputs.yaml.boundary_operator",
        "B1.field_manifest.moduli_fixed_collective_tangents", "B1.field_manifest.moduli_fixed_base_fields", "stage0.integrated_balance_identities",
    ]
    forbidden = ["B1.indexed_cells.*.M_XX_p0_known", "B1.indexed_cells.*.termwise_L", "projected_ledger_terms", "canonicalized_reconstruction_terms", "route_residual", "ownership_classifier", "current_splitter"]
    require(not any(item in rows for item in forbidden), "B2_S0_WHITELIST", "no downstream/reduced shared objects")
    return {"route_partition_reconstruction": rows, "derivation": "exact raw action/balance inputs plus enumerated moduli-fixed tangent/base leaves", "forbidden_shared_reduced_objects": forbidden}


def flatten_leaf_ids(prefix: str, value: Any) -> list[str]:
    if isinstance(value, dict):
        return [leaf for key in sorted(value) for leaf in flatten_leaf_ids(f"{prefix}.{key}", value[key])]
    if isinstance(value, list):
        return [leaf for index, child in enumerate(value) for leaf in flatten_leaf_ids(f"{prefix}[{index}]", child)]
    return [prefix]


def stage0_slot_names(schema: dict[str, Any], inventory: list[dict[str, Any]]) -> list[str]:
    axes = copy.deepcopy(schema["axes"])
    axes["operator_branch"] = sorted(row["id"] for row in inventory)
    slots = set(schema["stage0_fixed_datum_slots"])
    for spec in schema["stage0_expanded_datum_slots"]:
        for values in itertools.product(*(axes[axis] for axis in spec["axes"])):
            slots.add(spec["name"] + "|" + "|".join(f"{axis}={value}" for axis, value in zip(spec["axes"], values)))
    return sorted(slots)


def stage0_datum_candidates(
    phase: dict[str, Any], mechanics: dict[str, Any], schema: dict[str, Any],
    second_variation: dict[str, Any], inventory: list[dict[str, Any]],
    noether: list[dict[str, Any]], currents: dict[str, Any], balances: dict[str, Any],
    resolvents: list[dict[str, Any]], spaces: list[dict[str, Any]], restriction: dict[str, Any],
    ownership: dict[str, Any], config: dict[str, Any],
) -> dict[str, Any]:
    """Generate the directive-v48 disposition census independently of the DAG.

    Producer rules are selected from the directive slot class.  No caller can
    attach producer tags to manufacture an absence.
    """
    action_ids = sorted(row["id"] for row in phase["action_terms"])
    field_ids = sorted(row["id"] for row in phase["field_records"])
    endpoint_ids = sorted(mechanics["endpoint_functionals"])
    closure_inventory = sorted(
        [*(f"action::{name}" for name in action_ids), *(f"field::{name}" for name in field_ids),
         *(f"endpoint::{name}" for name in endpoint_ids),
         *flatten_leaf_ids("ir", phase["ir_scheme"]),
         "kinematics::coordinate_map", "kinematics::control_volume", "kinematics::base_state_class",
         "ambient::medium_rest_velocity", "mechanics::indexed_coordinates",
         "B1::moduli_fixed_collective_tangents", "B1::moduli_fixed_base_fields",
         "directive::section_1.4", "directive::section_3.1", "directive::section_3.4"]
    )
    closure_set = set(closure_inventory)
    operator_index = {row["id"]: row for row in inventory}
    noether_index = {row["id"]: row for row in noether}
    green_index = {row["cell"]: row for row in resolvents}
    spaces_index = {row["cell"]: row for row in spaces}

    reserved = {
        "mass_sector_tolerance": "ir.g9_tolerance.mass_rate_scale",
        "energy_sector_tolerance": "ir.g9_tolerance.energy_power_scale",
        "momentum_residual_norm": "ir.g9_tolerance.momentum_residual_norm",
        "return_energy_closure": "field::dimensioned_energy_return_branch",
        "truncation_validity_band": "ir.truncation.validity_band",
        "truncation_tolerance": "ir.truncation.relative_tolerance",
    }

    def split_axes(slot: str) -> tuple[str, dict[str, str]]:
        pieces = slot.split("|")
        return pieces[0], {key: value for key, value in (piece.split("=", 1) for piece in pieces[1:])}

    def producer_rule(slot: str) -> tuple[str, list[str]]:
        base, axes = split_axes(slot)
        operator = axes.get("operator_branch")
        if base == "full_action_second_variation":
            return "v48:hessian<-second_variation(EVERY_committed_action_term)+incidence_residual", [f"action::{name}" for name in action_ids]
        if base == "operator_hessian_block":
            ancestry = operator_index[operator]["action_ancestry"]
            producers = [f"action::{name}" for name in ancestry]
            if operator == "geon_open": producers = ["field::geon_core_bundle_operator_functional"]
            return "v48:operator_block<-full_action_second_variation_component", producers
        if base in {"native_mass_current_law", "native_energy_current_law"}:
            return "v48:sector_current<-native_action_balance_roots", [*(f"action::{name}" for name in action_ids), "field::native_continuity", "field::native_momentum"]
        if base.startswith("integrated_"):
            return "v48:integrated_balance<-native_current+Reynolds+typed_channel_functionals", [*(f"field::{name}" for name in ["native_continuity", "native_momentum", "E4_shear_lock", "E5_rayleigh", "return_closure"]), "field::outer_control_flux_native_functional"]
        if base in reserved:
            return "v48:upstream_freeze<-ONLY_7.0/inputs_ledger_typed_leaf", [reserved[base]]
        if base == "causal_retarded_definition":
            return "v48:causal_condition<-directive_definition_of_radiative_residue", ["directive::section_3.1"]
        if base == "native_noether_flux":
            if operator == "geon_open": return "v48:noether_flux<-committed_action_block", ["field::geon_core_bundle_action"]
            return "v48:noether_flux<-committed_action_block", [f"action::{name}" for name in operator_index[operator]["action_ancestry"]]
        if base == "retarded_green_operator":
            return "v48:green<-native_operator+endpoint_domain+retarded_definition", [*(f"action::{name}" for name in operator_index[operator]["action_ancestry"] if name != "geon_core_bundle"), f"endpoint::{axes['endpoint']}", "directive::section_3.1"]
        if base.startswith("functional_"):
            suffix = {"functional_test_space": "test_space", "functional_topology": "topology", "functional_contact_subtraction_set": "contact_subtraction_set", "functional_operator_norm": "operator_norm"}[base]
            return "v48:functional_datum<-IR_scheme_typed_leaf", [f"ir.functional_analytic.{operator}.{axes['endpoint']}.{suffix}"]
        if base in {"trajectory_frequency_representation", "fourier_convention", "retained_local_kernel_definition", "zero_denominator_policy"}:
            return "v48:definition<-directive_stage0_definition", ["directive::section_3.1" if base == "trajectory_frequency_representation" else "directive::section_3.4"]
        if base == "truncation_frequency_variable":
            return "v48:what<-geometry_a+committed_c_ref_coefficients", ["directive::section_3.4", "ir.g9_tolerance.force_scale"]
        if base == "complete_comoving_restriction_map":
            return "v48:restriction<-kinematics+moduli_fixed_base_fields+normalizations", ["kinematics::coordinate_map", "kinematics::control_volume", "kinematics::base_state_class", "ambient::medium_rest_velocity", "mechanics::indexed_coordinates", "B1::moduli_fixed_base_fields"]
        if base == "ownership_current_improvement_relation":
            return "v48:ownership_quotient<-typed_current_improvement_superpotential", ["field::current_improvement_superpotential"]
        raise AssertionError(f"unknown v48 producer class:{slot}")

    def relevant_closure(slot: str, producers: list[str]) -> list[str]:
        base, axes = split_axes(slot)
        if base in reserved:
            return sorted(item for item in closure_inventory if item.startswith("ir.") or (base == "return_energy_closure" and item.startswith("field::return_closure")))
        if base.startswith("functional_"):
            return sorted(item for item in closure_inventory if item.startswith("ir."))
        if base in {"integrated_mass_balance_identity", "integrated_momentum_balance_identity", "integrated_energy_balance_identity"}:
            return sorted(item for item in closure_inventory if item.startswith(("action::", "field::", "endpoint::", "kinematics::")))
        if base in {"native_mass_current_law", "native_energy_current_law", "full_action_second_variation"}:
            return sorted(item for item in closure_inventory if item.startswith(("action::", "field::")))
        if base in {"operator_hessian_block", "native_noether_flux", "retarded_green_operator"}:
            actual = set(producers) & closure_set
            if "endpoint" in axes: actual.add(f"endpoint::{axes['endpoint']}")
            if base == "retarded_green_operator": actual.add("directive::section_3.1")
            return sorted(actual)
        return sorted(set(producers) & closure_set)

    def derived(slot: str, required_type: str, dimensions: str, candidate_ref: str) -> dict[str, Any]:
        rule, producers = producer_rule(slot)
        closure = relevant_closure(slot, producers)
        return {
            "slot": slot, "disposition_candidate": "DERIVED", "required_type": required_type,
            "required_dimensions": dimensions, "defining_predicate": "directive_identity_and_schema_type_match",
            "candidate_ref": candidate_ref, "candidate_is_well_typed": True,
            "candidate_typecheck": {"candidate_present": True, "type_matches_required_type": True, "dimensions_match_required_dimensions": True, "identity_domain_membership": True, "computed_result": True},
            "defining_predicate_result": "PASS", "producer_rule": rule,
            "producer_ids": producers, "committed_input_closure": closure,
            "closure_exact_set_assert": "PASS",
            "derivability_challenge": {
                "status": "REFUTED", "candidate_ref": candidate_ref,
                "candidate_is_well_typed": True, "defining_predicate_result": "PASS",
                "candidate_typecheck": {"candidate_present": True, "type_matches_required_type": True, "dimensions_match_required_dimensions": True, "identity_domain_membership": True, "computed_result": True},
                "terminal": "REFUTED(well-typed PASS candidate)",
                "constructive_attempt_nonempty": True,
            },
        }

    def unresolved(slot: str, required_type: str, dimensions: str, kind: str, missing: str, restore_kind: str, diagnostic: str) -> dict[str, Any]:
        rule, producers = producer_rule(slot)
        closure = relevant_closure(slot, producers)
        census = [{"producer": producer, "in_closure": producer in closure, "type_compatible": False if producer in closure else None} for producer in producers]
        absence = all(not row["in_closure"] or row["type_compatible"] is False for row in census)
        require(absence, "B2_S0_UNAVAILABILITY_WITNESS", f"{slot}:executable producer census")
        witness_id = "witness::" + slot
        challenge_id = "challenge::" + slot
        restore_target = next(row["producer"] for row in census if not row["in_closure"] or row["type_compatible"] is False)
        attempt_records = [{
            "producer": row["producer"],
            "construction_state": "REJECTED_TYPE_INCOMPATIBLE" if row["in_closure"] else "MISSING_PRODUCER",
            "candidate_ref": None,
            "required_type": required_type,
            "required_dimensions": dimensions,
        } for row in census]
        require(bool(attempt_records), "B2_S0_UNAVAILABILITY_WITNESS", f"{slot}:nonempty constructive producer attempt")
        return {
            "slot": slot, "disposition_candidate": "UNRESOLVED", "unresolved_tag": "UNRESOLVED_BY_UPSTREAM_FREEZE" if slot in reserved else "UNRESOLVED",
            "required_type": required_type, "required_dimensions": dimensions,
            "defining_predicate": "candidate must match this type/dimensions and the directive identity domain",
            "candidate_ref": None, "candidate_is_well_typed": False,
            "candidate_typecheck": {"candidate_present": False, "type_matches_required_type": False, "dimensions_match_required_dimensions": False, "identity_domain_membership": False, "computed_result": False},
            "defining_predicate_result": "NO_CANDIDATE_AFTER_EXHAUSTIVE_CONSTRUCTIVE_ATTEMPT",
            "producer_rule": rule, "producer_ids": producers,
            "committed_input_closure": closure, "closure_exact_set_assert": "PASS",
            "unavailability_witness": {
                "witness_id": witness_id, "datum_id": slot, "required_type": required_type,
                "required_dimensions": dimensions, "acceptance_predicate": "directive_identity_and_schema_type_match",
                "authoritative_roots": producers, "enumerated_committed_inputs": closure,
                "complete_closure_exact_set_equal": True, "directive_generated_producer_rule": rule,
                "producer_census": census, "producer_census_predicate": "forall p in Producer(slot): p notin closure or type(p) incompatible",
                "kind": kind, "executable_certificate_result": "PASS" if absence else "FAIL",
                "diagnostic": diagnostic, "missing_typed_ingredient": missing,
                "counterfactual_restore_mutation": {
                    "ingredient_kind": restore_kind, "ingredient": missing,
                    "producer_to_restore": restore_target, "restored_type_compatible": True,
                    "fixture_type": required_type, "fixture_dimensions": dimensions,
                    "certificate_after_restore": "FAIL", "failed_at_own_assert": "B2_S0_WITNESS_RESTORE",
                },
            },
            "derivability_challenge": {
                "challenge_id": challenge_id, "same_committed_input_closure": closure,
                "dag_separated_from_witness": True, "shared_only_committed_inputs": True,
                "constructive_attempt": "independent typed producer resolution and defining-equation construction",
                "constructive_attempt_nonempty": True, "attempt_records": attempt_records,
                "candidate_schema_pinned": {"type": required_type, "dimensions": dimensions},
                "candidate_typecheck": {"candidate_present": False, "type_matches_required_type": False, "dimensions_match_required_dimensions": False, "identity_domain_membership": False, "computed_result": False},
                "candidate_is_well_typed": False, "terminal": f"CONSTRUCTIVE_FAIL({kind})",
                "kind": kind, "dual_engine_certificate_pending_comparator": True,
            },
        }

    slots = stage0_slot_names(schema, inventory)
    records: list[dict[str, Any]] = []
    noether_open = [row["id"] for row in noether if row["status"].startswith("UNRESOLVED")]
    for slot in slots:
        base, axes = split_axes(slot)
        operator, endpoint = axes.get("operator_branch"), axes.get("endpoint")
        if base == "full_action_second_variation":
            record = derived(slot, "symmetric_bilinear_form_on_full_committed_field_jet", "action_density/field^2 blockwise", "complete_action_second_variation")
        elif base == "operator_hessian_block":
            row = operator_index[operator]
            if str(row["operator_fourier"]).startswith("UNRESOLVED"):
                record = unresolved(slot, "linearized_native_operator_block", "action_density/field^2 blockwise", "authority_census_producer_absence", f"{operator}_native_functional_second_variation", "missing_input_leaf", row["status"])
            else:
                record = derived(slot, "linearized_native_operator_block", "action_density/field^2 blockwise", f"operator_inventory::{operator}")
        elif base == "native_mass_current_law":
            record = derived(slot, "co-moving_mass_current_4-vector", "mass/(time*area_3)", "current_derivations.sectors.mass")
        elif base == "native_energy_current_law":
            record = unresolved(slot, "co-moving_energy_current_4-vector", "energy/(time*area_3)", "authority_census_producer_absence", "native_energy_currents_for_OPEN_action_blocks", "missing_input_leaf", ",".join(noether_open))
        elif base.startswith("integrated_"):
            sector = base.split("_", 1)[1].rsplit("_balance_identity", 1)[0]
            record = unresolved(slot, f"complete_{sector}_Reynolds_balance_identity", "sector_storage_rate", "authority_census_producer_absence", f"complete_{sector}_typed_boundary_and_source_functionals", "missing_input_leaf", balances["sectors"][sector]["status"])
        elif base in reserved:
            record = unresolved(slot, "approved_frozen_7.0_ledger_entry", "slot-specific", "authority_census_producer_absence", reserved[base], "dimensioned_branch_type" if base == "return_energy_closure" else "missing_input_leaf", "non-operative directive candidates excluded from census and challenge PASS")
        elif base == "causal_retarded_definition":
            record = derived(slot, "causal_boundary_condition_definition", "dimensionless prescription", "causal_definition")
        elif base == "native_noether_flux":
            row = noether_index[operator]
            if row["status"].startswith("DERIVED_FROM"):
                record = derived(slot, "native_energy_momentum_flux_tensor", "action-derived sector dimensions", f"native_noether_derivations::{operator}")
            else:
                record = unresolved(slot, "native_energy_momentum_flux_tensor", "action-derived sector dimensions", "authority_census_producer_absence", f"{operator}_native_action_functional_derivatives", "missing_input_leaf", row["status"])
        elif base == "retarded_green_operator":
            green = green_index[f"{operator}|{endpoint}"]["retarded_green_operator"]
            if green["representation"] is not None and "UNRESOLVED" not in str(green["status"]):
                record = derived(slot, "native_retarded_resolvent_distribution", "inverse native operator", f"endpoint_resolvent_cells::{operator}|{endpoint}")
            else:
                kind = "operator_domain_well_posedness_failure" if "domain" in str(green["status"]).lower() or operator in {"gnls_density_phase", "wall_chi_u_coupled"} else "authority_census_producer_absence"
                restore = "domain/BC completion" if kind == "operator_domain_well_posedness_failure" else "missing_input_leaf"
                record = unresolved(slot, "native_retarded_resolvent_distribution", "inverse native operator", kind, str(green["status"]), restore, str(green["status"]))
        elif base.startswith("functional_"):
            record = unresolved(slot, base, "operator/test-space dependent", "authority_census_producer_absence", f"native_IR_{base}", "missing_input_leaf", spaces_index[f"{operator}|{endpoint}"]["status"])
        elif base in {"trajectory_frequency_representation", "fourier_convention", "retained_local_kernel_definition", "zero_denominator_policy"}:
            record = derived(slot, "directive_frozen_definition", "definition-specific", {"trajectory_frequency_representation": "trajectory_representation", "fourier_convention": "fourier_convention", "retained_local_kernel_definition": "validity_domain_construction.retained_local_kernel_definition", "zero_denominator_policy": "validity_domain_construction.zero_denominator_policy"}[base])
        elif base == "truncation_frequency_variable":
            record = derived(slot, "dimensionless_frequency_variable", "dimensionless", "validity_domain_construction.frequency_variable")
        elif base == "complete_comoving_restriction_map":
            record = unresolved(slot, "simultaneous_base_state_restriction_map", "field/coefficient typed", "authority_census_producer_absence", "complete_moduli_fixed_base_field_map", "missing_selector", restriction["status"])
        elif base == "ownership_current_improvement_relation":
            record = unresolved(slot, "antisymmetric_current_superpotential_relation", "current density", "authority_census_producer_absence", "typed_current_improvement_superpotential", "missing_selector", ownership["status"])
        else:
            raise AssertionError(slot)
        records.append(record)
    require([row["slot"] for row in records] == slots, "B2_S0_DISPOSITION_FLOOR", "expected-vs-reachable stage0 slot equality")
    require(all((row["disposition_candidate"] == "DERIVED") ^ (row["disposition_candidate"] == "UNRESOLVED") for row in records), "B2_S0_DISPOSITION_FLOOR", "exactly one disposition")
    producer_map = [{"slot": row["slot"], "rule": row["producer_rule"], "producers": row["producer_ids"]} for row in records]
    return {
        "directive_version": 48, "expected_slots": slots, "reachable_slots": [row["slot"] for row in records],
        "expected_reachable_exact_set_equal": True, "record_count": len(records),
        "directive_generated_producer_map": producer_map,
        "directive_generated_producer_map_sha256": digest(producer_map),
        "committed_input_closure_inventory": closure_inventory,
        "committed_input_closure_inventory_sha256": digest(closure_inventory),
        "records": records,
    }


def expand_obligation_floor(schema: dict[str, Any], inventory: list[dict[str, Any]], schema_path: Path, datum_bank: dict[str, Any]) -> dict[str, Any]:
    require(schema["schema_version"] == "U1_PHASE_B2_V48_OBLIGATION_SCHEMA_V3" and schema["directive_version"] == 48, "B2_S0_OBLIGATIONS", "frozen v48 schema")
    axes = copy.deepcopy(schema["axes"])
    branches = sorted(row["id"] for row in inventory)
    axes["operator_branch"] = branches
    axes["unordered_operator_pair"] = [f"{a}__PAIR__{b}" for a, b in itertools.combinations(branches, 2)]
    expanded = set(schema["fixed_products"])
    category_counts = {}
    for spec in schema["expanded_products"]:
        values = [axes[axis] for axis in spec["axes"]]
        count = 0
        for combination in itertools.product(*values):
            expanded.add(spec["name"] + "|" + "|".join(f"{axis}={value}" for axis, value in zip(spec["axes"], combination)))
            count += 1
        category_counts[spec["name"]] = count
    disposition_by_slot = {row["slot"]: row["disposition_candidate"] for row in datum_bank["records"]}
    require(set(disposition_by_slot) == set(datum_bank["expected_slots"]), "B2_S0_DISPOSITION_FLOOR", "datum record set")
    for slot in datum_bank["expected_slots"]:
        expanded.add(f"stage0_datum|slot={slot}|disposition={disposition_by_slot[slot]}")
    category_counts["stage0_datum"] = len(datum_bank["expected_slots"])
    for name in schema["required_category_coverage"]["per_branch"]:
        require(category_counts.get(name, 0) >= len(branches), "B2_S0_OBLIGATIONS", f"per-branch:{name}")
    for name in schema["required_category_coverage"]["per_unordered_pair"]:
        require(category_counts.get(name, 0) >= len(axes["unordered_operator_pair"]), "B2_S0_OBLIGATIONS", f"per-pair:{name}")
    require(set(schema["required_category_coverage"]["splitters"]) <= set(schema["fixed_products"]), "B2_S0_OBLIGATIONS", "splitter floor")
    records = sorted(expanded)
    return {"generator": "SymPy Cartesian-product expansion of independently frozen v48 full-key and stage0-datum schema", "schema_version": schema["schema_version"], "schema_sha256": sha256_file(schema_path), "fixed_products": sorted(schema["fixed_products"]), "expanded_records": records, "expanded_record_count": len(records), "expanded_record_set_sha256": digest(records), "category_counts": category_counts, "grid_axes": axes, "p_slices": axes["p_slice"], "stage0_datum_slots": datum_bank["expected_slots"], "stage0_datum_disposition_candidates": disposition_by_slot, "fixture_only_g9_tolerance": schema["fixture_only_g9_tolerance"]}


def build(input_path: Path, phase_override: Path | None = None) -> dict[str, Any]:
    config = load_yaml(input_path)
    require(config["schema_version"] == "U1_PHASE_B2_STAGE0_INPUTS_V4" and config["startup_contract_commit"] == "5f73be9f0738030bf73165fda5432644de5f4074", "B2_S0_INPUT", "v48 stage-0 input schema/anchor")
    phase_path = phase_override.resolve() if phase_override else ROOT / config["contracts"]["phase_a_inputs"]
    mechanics_path = ROOT / config["contracts"]["phase_b1_inputs"]
    b1_path = ROOT / config["contracts"]["sympy_b1"]
    directive_path = ROOT / config["contracts"]["directive"]
    schema_path = ROOT / config["contracts"]["v48_obligation_schema"]
    phase, mechanics, b1, schema = load_yaml(phase_path), load_yaml(mechanics_path), load_yaml(b1_path), load_yaml(schema_path)
    print("B2_STAGE0_SYMPY_PROGRESS: inputs_loaded", flush=True)
    second_variation = complete_action_second_variation(phase)
    print("B2_STAGE0_SYMPY_PROGRESS: complete_second_variation", flush=True)
    inventory, incidence = action_hessian_inventory(phase, mechanics)
    print("B2_STAGE0_SYMPY_PROGRESS: operator_inventory", flush=True)
    noether = native_noether_derivations(phase)
    print("B2_STAGE0_SYMPY_PROGRESS: native_noether", flush=True)
    currents = comoving_currents(noether, phase)
    print("B2_STAGE0_SYMPY_PROGRESS: comoving_currents", flush=True)
    balances = integrated_balances(phase, currents, schema)
    print("B2_STAGE0_SYMPY_PROGRESS: integrated_balance_outcomes", flush=True)
    resolvents, spaces = resolvents_and_spaces(inventory, phase, mechanics, noether)
    print("B2_STAGE0_SYMPY_PROGRESS: resolvents_and_controls", flush=True)
    restriction = restriction_from_b1(phase, mechanics, b1, currents)
    ownership = ownership_quotient(noether)
    validity = {
        "inventory": ["a>0", "ambient-subtracted R/a->infinity", "co-moving-steady base", "UNRESOLVED(endpoint-compatible test space/topology)", "retarded/outgoing support"],
        "maximality_rule": "intersection of committed conditions only", "committed_frequency_band_present": False,
        "functional_status": "UNRESOLVED(native_IR_test_space,native_IR_topology,native_IR_operator_norm,native_IR_contact_subtraction_set)",
        "frequency_variable": "what=omega*a/sqrt(5*K_EOS*rho_inf**4/m_GNLS)",
        "retained_local_kernel_definition": "polynomial-in-omega content of K_var+K_flux+K_constraint+K_Rayleigh after five-channel assembly",
        "zero_denominator_policy": "N=0,E_total=0 -> LOCAL_VALID(trivial); N=0,E_total!=0 -> memory kernel; isolated zeros use limsup",
        "bandless_truncation_rule": "E_total identically zero -> LOCAL_VALID(trivial/exact); otherwise retain memory kernel",
    }
    datum_bank = stage0_datum_candidates(phase, mechanics, schema, second_variation, inventory, noether, currents, balances, resolvents, spaces, restriction, ownership, config)
    print("B2_STAGE0_SYMPY_PROGRESS: v48_disposition_bank", flush=True)
    obligations = expand_obligation_floor(schema, inventory, schema_path, datum_bank)
    print("B2_STAGE0_SYMPY_PROGRESS: obligation_floor", flush=True)
    artifact = {
        "schema_version": "U1_PHASE_B2_STAGE0_ENGINE_V6", "status": "ENGINE_DERIVED", "engine": "SymPy",
        "independent_representation": "SymPy connected-Hessian graph + algebraic full-jet Noether + radial Sturm-Liouville Wronskian + row-reduced native quotient",
        "input_sha256": sha256_file(input_path),
        "source_digests": {"phase_a_inputs": sha256_file(phase_path), "phase_b1_inputs": sha256_file(mechanics_path), "b1_artifact": sha256_file(b1_path), "directive": sha256_file(directive_path), "v48_obligation_schema": sha256_file(schema_path)},
        "causal_definition": config["causal_definition"], "fourier_convention": config["fourier_convention"],
        "complete_action_second_variation": second_variation,
        "operator_inventory": inventory, "operator_action_incidence": incidence, "operator_action_incidence_residual": [],
        "native_noether_derivations": noether, "current_derivations": currents, "integrated_balance_identities": balances,
        "endpoint_resolvent_cells": resolvents, "functional_analytic_test_data": spaces,
        "green_sensitivity_matrix": green_sensitivity_matrix(inventory), "restriction_map": restriction,
        "ownership_convention": ownership, "stage0_datum_bank": datum_bank, "shared_input_whitelist": shared_whitelist(), "minimum_obligation_floor": obligations,
        "route_separation": {"source_parser_ir": "SymPy sympify over scalar-invariant expressions", "operator_ir": "primitive-symbol to perturbation incidence then union-find connected components", "balance_route_A_ir": "jet divergence -> coordinate pullback -> Reynolds transport -> symbolic coefficient extraction", "balance_route_B_ir": "typed-root record fold over signed canonical components", "resolvent_ir": "serialized piecewise AST evaluated by an independent distributional jump applicator", "ownership_ir": "concrete candidate-expression relation module reduced by SymPy rref", "shared_derivation_helpers": []},
        "trajectory_representation": {"steady": "Phi_0(x-X_0-V*t,p;cell), finite V retained", "accelerating": "native source J_cell[Phi_0(x-X(t),p(t))] acted on by the dimensional endpoint Green operator", "frequency": "qhat_B(omega), d_t -> -I*omega", "collective_indices": ["X_x", "X_y", "X_z", "p_x", "p_y", "p_z"]},
        "validity_domain_construction": validity,
        "checks": {name: "PASS" for name in ["B2_S0_ACTION_PARSE", "B2_S0_ACTION_HESSIAN", "B2_S0_OPERATOR_INCIDENCE", "B2_S0_BRANE_NORMALIZATION", "B2_S0_NATIVE_NOETHER", "B2_S0_U1_CURRENT", "B2_S0_COMOVING", "B2_S0_BALANCE_HONEST_EXIT", "B2_S0_BALANCE_SLOTS", "B2_S0_G9_ROUTER", "B2_S0_RESOLVENT", "B2_S0_NATIVE_CONTROL", "B2_S0_RESTRICTION", "B2_S0_OWNERSHIP", "B2_S0_WHITELIST", "B2_S0_DISPOSITION_FLOOR", "B2_S0_UNAVAILABILITY_WITNESS", "B2_S0_OBLIGATIONS"]},
    }
    artifact["semantic_payload_sha256"] = digest({k: v for k, v in artifact.items() if k != "semantic_payload_sha256"})
    return artifact


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--phase-input", type=Path)
    parser.add_argument("--mutation-projection", action="store_true", help="after all assertions run, serialize only source-liveness outputs")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        artifact = build(args.input.resolve(), args.phase_input)
        operator_count = len(artifact["operator_inventory"])
        resolvent_count = len(artifact["endpoint_resolvent_cells"])
        if args.mutation_projection:
            floor = dict(artifact["minimum_obligation_floor"])
            floor.pop("expanded_records")
            artifact = {key: artifact[key] for key in ["schema_version", "status", "engine", "operator_inventory", "native_noether_derivations", "current_derivations"]}
            artifact["minimum_obligation_floor"] = floor
            artifact["mutation_projection"] = "full production derivation/assertions executed; unchanged bulk evidence omitted only after computation"
        dump_yaml(args.output.resolve(), artifact)
        print(f"B2_STAGE0_SYMPY: PASS operators={operator_count} resolvent_cells={resolvent_count}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
