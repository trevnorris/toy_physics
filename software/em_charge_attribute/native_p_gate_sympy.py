#!/usr/bin/env python3
"""Genuine SymPy engine for the native-P constraint-class gate.

Every model, including all controls, is supplied as a quadratic Lagrangian to
the same build_H2 -> dirac_search -> search_G pipeline.  No constraint matrix,
constraint class, or G verdict is supplied by a model builder.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import sympy as sp


HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "reports" / "native_p_gate_artifacts"
SCHEMA = "native_p_constraint_gate/sympy/rebuild-v1"
K = sp.Matrix([sp.Integer(1), sp.Integer(2), sp.Integer(3)])
K2 = int(K.dot(K))
TOOTH_CODES = {
    "maxwell": 41,
    "gauged_hard_unit": 42,
    "bare_sigma": 43,
    "nonconserved_current": 44,
    "gauge_fixed_maxwell": 45,
    "global_only": 46,
}


class GateFailure(RuntimeError):
    pass


class ControlFailure(GateFailure):
    def __init__(self, tooth: str, message: str):
        super().__init__(f"CONTROL_ASSERT[{tooth}]: {message}")
        self.tooth = tooth


def require(test: object, message: str) -> None:
    if not bool(test):
        raise GateFailure(message)


def control_require(tooth: str, test: object, message: str) -> None:
    if not bool(test):
        raise ControlFailure(tooth, message)


def pb(f: sp.Expr, g: sp.Expr, pairs: list[tuple[sp.Symbol, sp.Symbol]]) -> sp.Expr:
    return sp.factor(sum(sp.diff(f, q) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, q) for q, p in pairs))


def expr_strings(items: list[sp.Expr]) -> list[str]:
    return [str(sp.factor(x)) for x in items]


def matrix_strings(m: sp.Matrix) -> list[list[str]]:
    return [[str(sp.factor(m[i, j])) for j in range(m.cols)] for i in range(m.rows)]


def independent_rows(m: sp.Matrix) -> list[int]:
    return list(m.T.rref()[1])


@dataclass
class Model:
    key: str
    label: str
    q: list[sp.Symbol]
    v: list[sp.Symbol]
    lagrangian: sp.Expr
    couplings: list[sp.Symbol]
    scalar_connection: list[int]
    spatial_connections: list[list[int]]
    multiplier_fields: list[int]
    current_data: dict[str, object] | None = None


@dataclass
class HamiltonianBuild:
    model: Model
    momenta: list[sp.Symbol]
    momentum_map: list[sp.Expr]
    hessian: sp.Matrix
    hessian_rank: int
    hessian_nullspace: list[sp.Matrix]
    primary_constraints: list[sp.Expr]
    hamiltonian: sp.Expr
    pairs: list[tuple[sp.Symbol, sp.Symbol]]


def build_H2(model: Model) -> HamiltonianBuild:
    """Legendre-transform an input Lagrangian; primaries come only from Hessian nulls."""
    n = len(model.q)
    require(len(model.v) == n, f"{model.key}: q/v size mismatch")
    momenta = list(sp.symbols(" ".join(f"Pi_{model.key}_{i}" for i in range(n)), real=True))
    momentum_map = [sp.diff(model.lagrangian, vi) for vi in model.v]
    hessian = sp.Matrix([[sp.diff(momentum_map[i], model.v[j]) for j in range(n)] for i in range(n)])
    nulls = hessian.T.nullspace()
    residual = sp.Matrix(momenta) - sp.Matrix(momentum_map).subs({vi: 0 for vi in model.v})
    primaries = [sp.factor((nu.T * residual)[0]) for nu in nulls]

    # Select an invertible row/column minor and set the arbitrary null velocities
    # to zero.  Their omitted terms are restored as total-Hamiltonian primary
    # multipliers by the Dirac algorithm.
    rank = int(hessian.rank())
    row_ids = independent_rows(hessian)
    col_ids = list(hessian.rref()[1])
    require(len(row_ids) == len(col_ids) == rank, f"{model.key}: Hessian minor selection failed")
    solution = {vi: sp.Integer(0) for vi in model.v}
    if rank:
        minor = hessian.extract(row_ids, col_ids)
        rhs = sp.Matrix([residual[i] for i in row_ids])
        solved = minor.inv() * rhs
        solution.update({model.v[j]: sp.factor(solved[a]) for a, j in enumerate(col_ids)})
    hc = sp.factor((sp.Matrix(momenta).dot(sp.Matrix(model.v)) - model.lagrangian).subs(solution))
    # A valid representative cannot retain an arbitrary velocity.
    require(not any(hc.has(vi) for vi in model.v), f"{model.key}: Legendre transform retained a velocity")
    return HamiltonianBuild(model, momenta, momentum_map, hessian, rank, nulls, primaries, hc,
                            list(zip(model.q, momenta)))


def linear_row(expr: sp.Expr, z: list[sp.Symbol]) -> sp.Matrix:
    return sp.Matrix([[sp.diff(expr, zi) for zi in z]])


def adds_independent(expr: sp.Expr, constraints: list[sp.Expr], z: list[sp.Symbol]) -> bool:
    row = linear_row(expr, z)
    if row == sp.zeros(1, len(z)):
        return False
    old = sp.Matrix.vstack(*(linear_row(c, z) for c in constraints)) if constraints else sp.zeros(0, len(z))
    return int(sp.Matrix.vstack(old, row).rank()) > int(old.rank())


def dirac_search(build: HamiltonianBuild) -> dict[str, object]:
    """Iterate all primary-preservation compatibility conditions to closure."""
    constraints = list(build.primary_constraints)
    stages = [0] * len(constraints)
    primary_count = len(constraints)
    z = build.model.q + build.momenta
    for stage in range(1, 2 * len(z) + 1):
        if not constraints:
            break
        multiplier_matrix = sp.Matrix([[pb(c, p, build.pairs) for p in build.primary_constraints]
                                       for c in constraints])
        hflow = sp.Matrix([pb(c, build.hamiltonian, build.pairs) for c in constraints])
        compat = multiplier_matrix.T.nullspace()
        new: list[sp.Expr] = []
        for left in compat:
            candidate = sp.factor((left.T * hflow)[0])
            if adds_independent(candidate, constraints + new, z):
                new.append(candidate)
        if not new:
            break
        constraints.extend(new)
        stages.extend([stage] * len(new))
    else:
        raise GateFailure(f"{build.model.key}: Dirac iteration did not close")

    m = sp.Matrix([[pb(a, b, build.pairs) for b in constraints] for a in constraints])
    rank = int(m.rank()) if constraints else 0
    kernel = m.nullspace() if constraints else []
    first = len(kernel)
    require(first == len(constraints) - rank, f"{build.model.key}: PB rank/nullity mismatch")
    classes: dict[str, str] = {}
    for i, c in enumerate(constraints):
        individually_fc = all(sp.simplify(m[i, j]) == 0 for j in range(len(constraints)))
        classes[f"C{i + 1}:{sp.factor(c)}"] = "FIRST_CLASS" if individually_fc else "SECOND_CLASS_COMPONENT"
    return {
        "constraints_expr": constraints,
        "constraint_stages": stages,
        "primary_count": primary_count,
        "constraint_count": len(constraints),
        "bracket_matrix_expr": m,
        "bracket_rank": rank,
        "first_class_count": first,
        "second_class_count": rank,
        "kernel_expr": kernel,
        "constraint_classes": classes,
    }


def proportional_to_k(action: list[sp.Expr]) -> bool:
    if len(action) != 3 or all(sp.simplify(x) == 0 for x in action):
        return False
    return all(sp.simplify(action[i] * K[j] - action[j] * K[i]) == 0 for i in range(3) for j in range(i + 1, 3))


def search_G(build: HamiltonianBuild, d: dict[str, object]) -> dict[str, object]:
    """Search the computed PB kernel for a first-class primary -> Gauss chain."""
    constraints: list[sp.Expr] = d["constraints_expr"]  # type: ignore[assignment]
    m: sp.Matrix = d["bracket_matrix_expr"]  # type: ignore[assignment]
    r = int(d["primary_count"])
    if not constraints or not r:
        primary_fc = []
    else:
        primary_fc = m[:, :r].nullspace()
    candidates: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    for coeff in primary_fc:
        primary = sp.factor(sum(coeff[i] * constraints[i] for i in range(r)))
        descendant = sp.factor(pb(primary, build.hamiltonian, build.pairs))
        paction = [sp.factor(pb(qi, primary, build.pairs)) for qi in build.model.q]
        saction = [sp.factor(pb(qi, descendant, build.pairs)) for qi in build.model.q]
        first_class_descendant = all(sp.simplify(pb(descendant, c, build.pairs)) == 0 for c in constraints)
        scalar_primary = any(sp.simplify(paction[i]) != 0 for i in build.model.scalar_connection)
        primary_only_scalar = scalar_primary and all(
            sp.simplify(a) == 0 for i, a in enumerate(paction)
            if i not in build.model.scalar_connection and i not in build.model.multiplier_fields
        )
        spatial_blocks = [
            {
                "field_indices": [i + 1 for i in block],
                "action": expr_strings([saction[i] for i in block]),
                "proportional_to_k": proportional_to_k([saction[i] for i in block]),
            }
            for block in build.model.spatial_connections
        ]
        spatial_gauss = any(bool(block["proportional_to_k"]) for block in spatial_blocks)
        nontrivial = any(sp.simplify(x) != 0 for x in saction)
        accepted = bool(first_class_descendant and primary_only_scalar and spatial_gauss and nontrivial)
        descendant_zero = sp.simplify(descendant) == 0
        secondary_action_non_gradient = not spatial_gauss
        descendant_rejection_certified = bool(descendant_zero or secondary_action_non_gradient)
        rejection_reason = (
            "DESCENDANT_ZERO" if descendant_zero else
            "SECONDARY_ACTION_NON_GRADIENT" if secondary_action_non_gradient else
            "OTHER_GAUSS_CRITERION"
        )
        item = {
            "primary": str(primary),
            "descendant": str(descendant),
            "primary_action": expr_strings(paction),
            "secondary_action": expr_strings(saction),
            "spatial_action_blocks": spatial_blocks,
            "first_class_descendant": first_class_descendant,
            "scalar_primary_end": primary_only_scalar,
            "spatial_gradient_action": spatial_gauss,
            "descendant_zero": descendant_zero,
            "secondary_action_non_gradient": secondary_action_non_gradient,
            "descendant_rejection_certified": descendant_rejection_certified,
            "computed_rejection_reason": rejection_reason,
        }
        (candidates if accepted else rejected).append(item)
    gauss_candidates = len(candidates)
    additional = gauss_candidates > 0
    return {
        "computed_kernel_dimension": len(d["kernel_expr"]),
        "computed_first_class_primaries": len(primary_fc),
        "kernel_basis": [expr_strings(list(x)) for x in d["kernel_expr"]],
        "gauss_candidates": gauss_candidates,
        "additional_G_exists": additional,
        "candidates": candidates,
        "rejected_first_class_primaries": rejected,
    }


def tuned_rejection_guard(
    theory: str,
    condition: str,
    d: dict[str, object],
    g: dict[str, object],
) -> dict[str, object]:
    """Assert that every rejected tuned FC primary fails by its computed descendant."""
    fc = int(d["first_class_count"])
    if fc == 0:
        return {"status": "NOT_APPLICABLE_FC0", "computed": True, "checked_directions": 0}

    rejected: list[dict[str, object]] = g["rejected_first_class_primaries"]  # type: ignore[assignment]
    candidates: list[dict[str, object]] = g["candidates"]  # type: ignore[assignment]
    primary_fc = int(g["computed_first_class_primaries"])
    require(primary_fc == len(rejected) + len(candidates),
            f"HARDENING-TUNED-DESCENDANT-REJECTION THEORY-{theory} {condition}: "
            "first-class primary accounting mismatch")
    records: list[dict[str, object]] = []
    for index, item in enumerate(rejected, 1):
        # This is deliberately stronger than the full G verdict: a rejected
        # tuned direction is not allowed to hide a nonzero k-gradient descendant
        # behind failure of some other (for example scalar-primary) criterion.
        certified = bool(item["descendant_rejection_certified"])
        require(certified,
                f"HARDENING-TUNED-DESCENDANT-REJECTION THEORY-{theory} {condition} "
                f"direction {index}: rejected primary has a nonzero proportional-to-k descendant")
        records.append({
            "direction": index,
            "primary": item["primary"],
            "descendant": item["descendant"],
            "secondary_action": item["secondary_action"],
            "spatial_action_blocks": item["spatial_action_blocks"],
            "descendant_zero": item["descendant_zero"],
            "secondary_action_non_gradient": item["secondary_action_non_gradient"],
            "computed_rejection_reason": item["computed_rejection_reason"],
        })
    return {
        "status": "PASS",
        "computed": True,
        "condition": condition,
        "first_class_count": fc,
        "first_class_primaries": primary_fc,
        "accepted_gauss_directions": len(candidates),
        "checked_directions": len(records),
        "records": records,
    }


def serialize_pipeline(build: HamiltonianBuild, d: dict[str, object], g: dict[str, object]) -> dict[str, object]:
    return {
        "model": build.model.key,
        "lagrangian": str(sp.factor(build.model.lagrangian)),
        "coordinates": [str(x) for x in build.model.q],
        "velocities": [str(x) for x in build.model.v],
        "momenta": [str(x) for x in build.momenta],
        "momentum_map": expr_strings(build.momentum_map),
        "hessian": matrix_strings(build.hessian),
        "hessian_rank": build.hessian_rank,
        "hessian_nullity": len(build.hessian_nullspace),
        "primary_constraints": expr_strings(build.primary_constraints),
        "H2": str(sp.factor(build.hamiltonian)),
        "constraints": expr_strings(d["constraints_expr"]),
        "constraint_stages": d["constraint_stages"],
        "constraint_classes": d["constraint_classes"],
        "bracket_matrix": matrix_strings(d["bracket_matrix_expr"]),
        "bracket_rank": d["bracket_rank"],
        "first_class_count": d["first_class_count"],
        "second_class_count": d["second_class_count"],
        "G_search": g,
    }


def execute(model: Model) -> tuple[HamiltonianBuild, dict[str, object], dict[str, object], dict[str, object]]:
    build = build_H2(model)
    d = dirac_search(build)
    g = search_G(build, d)
    return build, d, g, serialize_pipeline(build, d, g)


def native_model(theory: str, substitutions: dict[sp.Symbol, sp.Expr] | None = None) -> Model:
    require(theory in {"A", "C"}, theory)
    p = list(sp.symbols(f"p{theory}1:4", real=True))
    u = list(sp.symbols(f"u{theory}1:4", real=True))
    s, lam, sigma, b = sp.symbols(f"s{theory} lambda{theory} sigma{theory} b{theory}", real=True)
    xi = sp.symbols(f"xi{theory}", real=True)
    q = p + u + [s, lam, sigma, b] + ([xi] if theory == "C" else [])
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    vp, vu = sp.Matrix(v[:3]), sp.Matrix(v[3:6])
    vs, vb = v[6], v[9]
    gt, gs, gd = sp.symbols(f"g_t{theory} g_s{theory} g_d{theory}", real=True)
    gb = sp.symbols(f"g_b{theory}", real=True)
    ru, ap, au, ab = sp.symbols(f"rho_u{theory} a_p{theory} a_u{theory} a_b{theory}", positive=True)
    pm, um = sp.Matrix(p), sp.Matrix(u)
    transverse = K2 * sp.eye(3) - K * K.T
    kinetic = (vp.dot(vp) + vs**2 + ru * vu.dot(vu)) / 2 + gt * vp.dot(vu)
    if theory == "A":
        kinetic += vb**2 / 2
    potential = ap * pm.dot(pm) / 2 + au * (um.T * transverse * um)[0] / 2
    potential += gs * K2 * pm.dot(um) + gd * (K.dot(pm))**2 / 2 + ab * b**2 / 2
    if theory == "A":
        potential += gb * b * K.dot(pm)
    lagrangian = kinetic - potential + lam * s + sigma * K.dot(um)
    if theory == "C":
        lagrangian += xi * b
    couplings = [gt, gs, gd] + ([gb] if theory == "A" else [])
    if substitutions:
        lagrangian = sp.factor(lagrangian.subs(substitutions))
        couplings = [c for c in couplings if c not in substitutions]
    return Model(f"native_{theory}", f"THEORY-{theory}", q, v, lagrangian, couplings,
                 scalar_connection=[], spatial_connections=[list(range(3)), list(range(3, 6))],
                 multiplier_fields=[7, 8] + ([10] if theory == "C" else []))


def maxwell_model(key: str = "maxwell", *, kinetic_a0: bool = False, defect: bool = False) -> Model:
    a0, a1, a2, a3 = sp.symbols(f"{key}_A0 {key}_A1 {key}_A2 {key}_A3", real=True)
    q = [a0, a1, a2, a3]
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    avec, va = sp.Matrix(q[1:]), sp.Matrix(v[1:])
    electric = va - K * a0
    magnetic = K2 * sp.eye(3) - K * K.T
    lagrangian = electric.dot(electric) / 2 - (avec.T * magnetic * avec)[0] / 2
    if kinetic_a0:
        lagrangian += v[0] ** 2 / 2
    current_data = None
    if key == "nonconserved_current":
        rho, rho_dot = sp.symbols("nc_rho nc_rho_dot", real=True)
        j = sp.Matrix(sp.symbols("nc_j1:4", real=True))
        lagrangian += rho * a0 + j.dot(avec)
        current_data = {"rho": rho, "rho_dot": rho_dot, "j": j, "impose_conservation": not defect}
    return Model(key, key, q, v, lagrangian, [], [0], [[1, 2, 3]], [], current_data)


def gauged_hard_unit_model(*, remove_hard: bool = False) -> Model:
    base = maxwell_model("gauged_hard")
    f1, f2, lam = sp.symbols("gh_phi1 gh_phi2 gh_lambda", real=True)
    q = base.q + [f1, f2] + ([] if remove_hard else [lam])
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    # Rebuild Maxwell with the new velocity symbols.
    avec, va = sp.Matrix(q[1:4]), sp.Matrix(v[1:4])
    electric = va - K * q[0]
    magnetic = K2 * sp.eye(3) - K * K.T
    vf = sp.Matrix(v[4:6])
    phi = sp.Matrix([f1, f2])
    rotphi = sp.Matrix([-f2, f1])
    d0phi = vf - q[0] * rotphi
    lagrangian = electric.dot(electric) / 2 - (avec.T * magnetic * avec)[0] / 2 + d0phi.dot(d0phi) / 2
    lagrangian -= phi.dot(phi) / 2  # harmless invariant spatial/potential representative
    if not remove_hard:
        lagrangian += lam * (phi.dot(phi) - 1) / 2
    return Model("gauged_hard_unit", "gauged hard unit", q, v, lagrangian, [], [0], [[1, 2, 3]],
                 [] if remove_hard else [6])


def bare_sigma_model(*, remove_hard: bool = False) -> Model:
    phi = list(sp.symbols("bs_phi1:5", real=True))
    lam = sp.symbols("bs_lambda", real=True)
    q = phi + ([] if remove_hard else [lam])
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    lagrangian = sum(x * x for x in v[:4]) / 2 - sum(x * x for x in phi) / 2
    if not remove_hard:
        lagrangian += lam * (sum(x * x for x in phi) - 1) / 2
    return Model("bare_sigma", "bare O(4) hard sigma", q, v, lagrangian, [], [], [],
                 [] if remove_hard else [4])


def gauge_fixed_maxwell_model(*, remove_coulomb: bool = False) -> Model:
    a0, a1, a2, a3, eta = sp.symbols("gf_A0 gf_A1 gf_A2 gf_A3 gf_eta", real=True)
    zeta = sp.symbols("gf_zeta", real=True)
    q = [a0, a1, a2, a3, eta] + ([] if remove_coulomb else [zeta])
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    avec, va = sp.Matrix(q[1:4]), sp.Matrix(v[1:4])
    electric = va - K * a0
    magnetic = K2 * sp.eye(3) - K * K.T
    lagrangian = electric.dot(electric) / 2 - (avec.T * magnetic * avec)[0] / 2 + eta * a0
    if not remove_coulomb:
        lagrangian += zeta * K.dot(avec)
    return Model("gauge_fixed_maxwell", "Coulomb gauge Maxwell", q, v, lagrangian, [], [0], [[1, 2, 3]],
                 [4] + ([] if remove_coulomb else [5]))


def global_only_model(*, add_maxwell: bool = False) -> Model:
    if add_maxwell:
        return maxwell_model("global_only_ablation")
    x, y = sp.symbols("global_x global_y", real=True)
    q = [x, y]
    v = list(sp.symbols("d_global_x d_global_y", real=True))
    lagrangian = (v[0] ** 2 + v[1] ** 2) / 2 - (x**2 + y**2) / 2
    return Model("global_only", "global U(1) complex scalar", q, v, lagrangian, [], [], [], [])


CONTROL_BUILDERS: dict[str, Callable[[bool], Model]] = {
    "maxwell": lambda a: maxwell_model(kinetic_a0=a),
    "gauged_hard_unit": lambda a: gauged_hard_unit_model(remove_hard=a),
    "bare_sigma": lambda a: bare_sigma_model(remove_hard=a),
    "nonconserved_current": lambda a: maxwell_model("nonconserved_current", defect=not a),
    "gauge_fixed_maxwell": lambda a: maxwell_model("gauge_fixed_ablation") if a else gauge_fixed_maxwell_model(),
    "global_only": lambda a: global_only_model(add_maxwell=a),
}


def classify_control(tooth: str, data: dict[str, object], model: Model) -> str:
    fc = int(data["first_class_count"])
    sc = int(data["second_class_count"])
    gc = int(data["G_search"]["gauss_candidates"])
    if tooth == "maxwell":
        return "FIRST_CLASS_GAUSS" if gc > 0 and sc == 0 else "NO_PRIMARY_GAUSS_CHAIN"
    if tooth == "gauged_hard_unit":
        return "MIXED" if gc > 0 and fc > 0 and sc > 0 else "NOT_MIXED"
    if tooth == "bare_sigma":
        return "SECOND_CLASS_RADIAL_NO_GAUSS" if fc == 0 and sc > 0 and gc == 0 else "BAD_SIGMA_CLASS"
    if tooth == "nonconserved_current":
        return "INCONSISTENT_PRESERVATION" if gc > 0 and data.get("current_inconsistent") is True else "CURRENT_CONSISTENT"
    if tooth == "gauge_fixed_maxwell":
        return "SECOND_CLASS_NO_LOCAL_GAUGE" if fc == 0 and sc > 0 and gc == 0 else "LOCAL_GAUGE_REMAINS"
    if tooth == "global_only":
        return "GLOBAL_CHARGE_NO_LOCAL_GAUSS" if data["hessian_nullity"] == 0 and fc == gc == 0 else "LOCAL_GAUGE_PRESENT"
    raise GateFailure(tooth)


def run_control(tooth: str, ablate: bool = False) -> dict[str, object]:
    model = CONTROL_BUILDERS[tooth](ablate)
    build, d, _, data = execute(model)
    if model.current_data is not None:
        r = int(d["primary_count"])
        m: sp.Matrix = d["bracket_matrix_expr"]  # type: ignore[assignment]
        primary_directions = m[:, :r].nullspace()
        require(bool(primary_directions), "current control has no computed first-class primary")
        coeff = primary_directions[0]
        constraints: list[sp.Expr] = d["constraints_expr"]  # type: ignore[assignment]
        primary = sp.factor(sum(coeff[i] * constraints[i] for i in range(r)))
        gauss = sp.factor(pb(primary, build.hamiltonian, build.pairs))
        rho = model.current_data["rho"]
        rho_dot = model.current_data["rho_dot"]
        preservation = sp.factor(sp.diff(gauss, rho) * rho_dot + pb(gauss, build.hamiltonian, build.pairs))
        require(preservation.has(rho_dot), "current preservation lost explicit rho_dot")
        tested = preservation
        conservation_rule = None
        if model.current_data["impose_conservation"]:
            solutions = sp.solve(sp.Eq(preservation, 0), rho_dot)
            require(bool(solutions), "could not solve the computed continuity equation")
            conservation_rule = {rho_dot: solutions[0]}
            tested = sp.factor(preservation.subs(conservation_rule))
        data["current_preservation"] = {
            "gauss": str(gauss), "raw": str(preservation),
            "conservation_rule": None if conservation_rule is None else {str(k): str(v) for k, v in conservation_rule.items()},
            "tested": str(tested),
        }
        data["current_inconsistent"] = sp.simplify(tested) != 0
    classification = classify_control(tooth, data, model)
    data["name"] = model.label
    data["classification"] = classification
    if tooth == "maxwell":
        control_require(tooth, classification == "FIRST_CLASS_GAUSS", "computed Maxwell Gauss chain missing")
    elif tooth == "gauged_hard_unit":
        control_require(tooth, classification == "MIXED", "computed first/second-class mixture missing")
    elif tooth == "bare_sigma":
        control_require(tooth, classification == "SECOND_CLASS_RADIAL_NO_GAUSS", "computed radial sector is not purely second class")
    elif tooth == "nonconserved_current":
        control_require(tooth, classification == "INCONSISTENT_PRESERVATION", "continuity defect did not spoil Gauss preservation")
    elif tooth == "gauge_fixed_maxwell":
        control_require(tooth, classification == "SECOND_CLASS_NO_LOCAL_GAUGE", "complete Coulomb fixing left a local G")
    elif tooth == "global_only":
        control_require(tooth, classification == "GLOBAL_CHARGE_NO_LOCAL_GAUSS", "global charge was mistaken for a local constraint")
    return data


def run_controls() -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    controls = [run_control(t) for t in CONTROL_BUILDERS]
    ablations: list[dict[str, str]] = []
    for tooth in CONTROL_BUILDERS:
        try:
            run_control(tooth, True)
        except ControlFailure as exc:
            require(exc.tooth == tooth, f"ablation {tooth} fired at {exc.tooth}")
            ablations.append({"tooth": tooth, "status": "FIRED_AT_OWN_ASSERT", "message": str(exc)})
        else:
            raise GateFailure(f"ablation {tooth} did not fire")
    return controls, ablations


def coupling_guard(build: HamiltonianBuild, d: dict[str, object]) -> dict[str, object]:
    """GUARD-COUPLINGS-ENTER: inspect computed Legendre/constraint/PB objects."""
    inspected = list(build.momentum_map) + list(build.hessian) + list(d["constraints_expr"]) + list(d["bracket_matrix_expr"])
    locations: dict[str, list[str]] = {}
    for coupling in build.model.couplings:
        hits: list[str] = []
        if any(x.has(coupling) for x in build.momentum_map): hits.append("momentum_map")
        if any(x.has(coupling) for x in build.hessian): hits.append("hessian")
        if any(x.has(coupling) for x in d["constraints_expr"]): hits.append("constraints")
        if any(x.has(coupling) for x in d["bracket_matrix_expr"]): hits.append("bracket_matrix")
        require(bool(hits), f"GUARD-COUPLINGS-ENTER failed: {coupling} dropped out")
        locations[str(coupling)] = hits
    return {"status": "PASS", "computed": True, "locations": locations,
            "inspected_expression_count": len(inspected)}


def operator_basis() -> dict[str, object]:
    return {
        "decision_order": 2,
        "notation": "p_i=P_parallel,i; b=chi-chi_vac; k=(1,2,3) is a nonzero Fourier representative",
        "quadratic_couplings_A": ["g_t dot(p).dot(u)", "-g_s k^2 p.u", "-(g_d/2)(k.p)^2", "-g_b b(k.p)"],
        "quadratic_couplings_C": ["g_t dot(p).dot(u)", "-g_s k^2 p.u", "-(g_d/2)(k.p)^2"],
        "empty_families": [
            "one-time/one-space bilinears are T-odd for T-even configuration fields",
            "epsilon bilinears are parity-odd or integration-by-parts zero at order two",
            "undifferentiated u is excluded by u-shift",
            "in C, b=p^2 starts at amplitude two, so b-cross terms start above quadratic order",
        ],
        "completeness": "At quadratic order, SO(3), parity, T, u-shift, <=2 derivatives, IBP and the holonomic identities leave exactly the displayed four (A) or three (C) cross terms. Fixed higher orders are finite associated-graded quotients; quadratic absence is decisive by the v6 linearization lever.",
    }


def native_analysis(theory: str) -> dict[str, object]:
    model = native_model(theory)
    build, d, g, data = execute(model)
    guard = coupling_guard(build, d)
    gt = model.couplings[0]
    ru = next(x for x in model.lagrangian.free_symbols if str(x) == f"rho_u{theory}")
    determinant = sp.factor(build.hessian.extract(range(6), range(6)).det())
    determinant_per_component = sp.factor(build.hessian.extract([0, 3], [0, 3]).det())
    require(sp.simplify(determinant - (ru - gt**2) ** 3) == 0,
            f"THEORY-{theory}: unexpected physical kinetic-Hessian degeneracy")
    require(sp.simplify(determinant_per_component - (ru - gt**2)) == 0,
            f"THEORY-{theory}: kinetic component determinant mismatch")

    # The open positive kinetic stratum is the admissible regular phase-space
    # stratum.  Boundary builds are independently run to ensure a tuned Hessian
    # null does not masquerade as a scalar-primary U(1) chain.
    boundary: list[dict[str, object]] = []
    rankdrop_sweep: list[dict[str, object]] = []
    common_null_solutions: dict[str, object] = {}
    sweep_seed = 260712 + (1 if theory == "A" else 3)
    sweep_rng = random.Random(sweep_seed)
    for sign in (-1, 1):
        free = {str(x): x for x in model.lagrangian.free_symbols}
        ap, au, gs = free[f"a_p{theory}"], free[f"a_u{theory}"], model.couplings[1]
        gd = free[f"g_d{theory}"]
        gb = free.get(f"g_b{theory}")
        # Compute the exact potential action on a transverse physical Hessian
        # null.  Its simultaneous zero is the only locus on which the rank-drop
        # chain can stop and become first class; solve it before running Dirac.
        physical_q = model.q[:6]
        potential = -model.lagrangian.subs({vi: 0 for vi in model.v})
        potential_matrix = sp.hessian(potential, physical_q)
        transverse_test = sp.Matrix([2, -1, 0])  # K.dot(test)=0
        null_test = sp.Matrix(list(-sign * transverse_test) + list(transverse_test))
        null_residual = [sp.factor(x) for x in potential_matrix * null_test]
        solved_common = sp.solve(null_residual, [ap, gs], dict=True)
        require(bool(solved_common), f"THEORY-{theory} sign {sign}: common-null solve failed")
        common_null_solutions[str(sign)] = {
            "potential_null_residual": expr_strings(null_residual),
            "solutions": [{str(k): str(v) for k, v in sol.items()} for sol in solved_common],
        }
        representative = {
            gt: sp.Integer(sign), ru: sp.Integer(1), au: sp.Integer(1), ap: sp.Integer(1),
            gs: sp.Rational(15 * sign, 28), gd: sp.Integer(0),
        }
        if gb is not None:
            representative[gb] = sp.Integer(0)
        scans = [
            ("generic semidefinite boundary", {gt: sp.Integer(sign), ru: sp.Integer(1)}),
            ("first rank-drop locus, non-common-null representative", representative),
            ("common Hessian/potential null locus",
             {gt: sp.Integer(sign), ru: sp.Integer(1), ap: 14 * au, gs: sign * au}),
        ]
        for label, subs in scans:
            bm = native_model(theory, subs)
            # Boundary artifacts consume only decisive ranks/search results.
            # The Hamiltonian, Dirac, and G-search calculations are unchanged.
            bb = build_H2(bm)
            bd = dirac_search(bb)
            bg = search_G(bb, bd)
            condition = f"rho_u=1, g_t={sign}; {label}"
            rejection_guard = tuned_rejection_guard(theory, condition, bd, bg)
            boundary.append({
                "condition": condition,
                "substitutions": {str(k): str(v) for k, v in subs.items()},
                "coordinates": [str(q) for q in bm.q],
                "first_class_count": bd["first_class_count"],
                "first_class_primaries": bg["computed_first_class_primaries"],
                "gauss_candidates": bg["gauss_candidates"],
                "additional_G_exists": bg["additional_G_exists"],
                "hessian_nullity": len(bb.hessian_nullspace),
                "rejected_first_class_primaries": bg["rejected_first_class_primaries"],
                "tuned_descendant_rejection_guard": rejection_guard,
            })

        # A fixed-seed randomized sweep samples the non-common portion of the
        # first rank-drop locus.  It is evidence about tuned strata, not a claim
        # of exhaustive symbolic stratification.
        for sample in range(1, 4):
            au_value = sweep_rng.randint(1, 5)
            ap_value = sweep_rng.randint(1, 40)
            while ap_value == 14 * au_value:
                ap_value = sweep_rng.randint(1, 40)
            sweep_subs = {
                gt: sp.Integer(sign),
                ru: sp.Integer(1),
                au: sp.Integer(au_value),
                ap: sp.Integer(ap_value),
                gs: sp.Rational(sign * (ap_value + 14 * au_value), 28),
                gd: sp.Integer(0),
            }
            if gb is not None:
                sweep_subs[gb] = sp.Integer(0)
            sm = native_model(theory, sweep_subs)
            sb = build_H2(sm)
            sd = dirac_search(sb)
            sg = search_G(sb, sd)
            condition = f"rank-drop randomized sample sign={sign} index={sample}"
            rejection_guard = tuned_rejection_guard(theory, condition, sd, sg)
            rankdrop_sweep.append({
                "condition": condition,
                "substitutions": {str(k): str(v) for k, v in sweep_subs.items()},
                "first_class_count": sd["first_class_count"],
                "first_class_primaries": sg["computed_first_class_primaries"],
                "gauss_candidates": sg["gauss_candidates"],
                "additional_G_exists": sg["additional_G_exists"],
                "hessian_nullity": len(sb.hessian_nullspace),
                "tuned_descendant_rejection_guard": rejection_guard,
            })

    scanned = boundary + rankdrop_sweep
    all_scans = [g] + [{"gauss_candidates": x["gauss_candidates"], "additional_G_exists": x["additional_G_exists"]} for x in scanned]
    any_g = any(bool(x["additional_G_exists"]) for x in all_scans)
    tuned_guards = [x["tuned_descendant_rejection_guard"] for x in scanned
                    if x["tuned_descendant_rejection_guard"]["status"] == "PASS"]
    require(bool(tuned_guards), f"THEORY-{theory}: tuned FC strata disappeared from hardening audit")
    tuned_guard = {
        "status": "PASS",
        "computed": True,
        "checked_strata": len(tuned_guards),
        "checked_directions": sum(int(x["checked_directions"]) for x in tuned_guards),
        "records": tuned_guards,
    }
    data["couplings"] = [str(x) for x in model.couplings]
    data["coupling_guard"] = guard
    data["tuned_descendant_rejection_guard"] = tuned_guard
    data["coupling_scan"] = {
        "physical_kinetic_hessian_determinant": str(determinant),
        "physical_kinetic_determinant_per_component": str(determinant_per_component),
        "only_kinetic_hessian_degeneracy": f"{determinant_per_component}=0",
        "regular_stability_stratum": f"rho_u{theory}>0 and {determinant_per_component}>0",
        "regular_result": {"first_class_count": d["first_class_count"], **g},
        "semidefinite_boundary": boundary,
        "rankdrop_representative_sweep": {
            "method": "fixed-seed randomized representative points on the solved non-common rank-drop locus",
            "seed": sweep_seed,
            "sample_count": len(rankdrop_sweep),
            "samples": rankdrop_sweep,
        },
        "boundary_rank_polynomial": "a_p + 14*a_u - 28*sign(g_t)*g_s; its zero locus is rerun, followed by its common-null sublocus a_p=14*a_u, g_s=sign(g_t)*a_u",
        "computed_common_null_solutions": common_null_solutions,
        "all_spatial_couplings_free": [str(x) for x in model.couplings[1:]],
        "tuned_scope": "residual solved symbolically and non-common rank-drop locus sampled at fixed-seed randomized representative points; not an exhaustive symbolic stratification of every tuned measure-zero sublocus",
    }
    data["G_search"]["additional_G_exists"] = any_g
    data["G_search"]["gauss_candidates"] = sum(int(x["gauss_candidates"]) for x in all_scans)
    data["source_test"] = {"searched_source_free_first": True, "jA_added": any_g, "j_sourced": any_g}
    data["shear_duplicate"] = {"applicability": "REQUIRES_G" if any_g else "NOT_APPLICABLE_NO_G", "macCullagh_transverse_modes": 2}
    if any_g:
        data["verdict"] = "FIRST_CLASS_TUNED_INVERSE_DESIGN" if not g["additional_G_exists"] else "FIRST_CLASS_GENERIC_EM_CANDIDATE"
    else:
        data["verdict"] = "NATIVE_P_NO_EMERGENT_GAUSS"
    data["decision_order"] = "quadratic (field-amplitude order 2)"
    return data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--ablate-tooth", choices=sorted(CONTROL_BUILDERS))
    args = parser.parse_args()
    if args.ablate_tooth:
        try:
            run_control(args.ablate_tooth, True)
        except ControlFailure as exc:
            print(f"ABLATION_FIRED_AT_OWN_ASSERT:{exc.tooth}")
            print(exc)
            return TOOTH_CODES[exc.tooth]
        print(f"ABLATION_DID_NOT_FIRE:{args.ablate_tooth}")
        return 90

    args.out_dir.mkdir(parents=True, exist_ok=True)
    controls, ablations = run_controls()
    # GUARD-SEARCH-CAPABLE consumes search output from the same pipeline.
    by_model = {x["model"]: x for x in controls}
    search_capable = (by_model["maxwell"]["G_search"]["gauss_candidates"] > 0 and
                      by_model["gauged_hard_unit"]["G_search"]["gauss_candidates"] > 0)
    require(search_capable, "GUARD-SEARCH-CAPABLE failed on Maxwell/gauged-hard-unit")
    theories = {t: native_analysis(t) for t in ("A", "C")}
    result = {
        "schema": SCHEMA,
        "engine": "SymPy",
        "engine_version": sp.__version__,
        "pipeline": "build_H2 -> dirac_search -> search_G",
        "phase_space": {"topology": "R^3", "smearings": "C_c^infinity(R^3)", "boundary": "decay at infinity", "zero_modes": "k=0 global sector excluded from local G"},
        "operator_basis": operator_basis(),
        "theories": theories,
        "controls": controls,
        "ablations": ablations,
        "guards": {
            "GUARD-COUPLINGS-ENTER": {t: theories[t]["coupling_guard"] for t in theories},
            "GUARD-SEARCH-CAPABLE": {"status": "PASS", "computed": True,
                "maxwell_gauss_candidates": by_model["maxwell"]["G_search"]["gauss_candidates"],
                "gauged_hard_unit_gauss_candidates": by_model["gauged_hard_unit"]["G_search"]["gauss_candidates"]},
            "HARDENING-TUNED-DESCENDANT-REJECTION": {
                t: theories[t]["tuned_descendant_rejection_guard"] for t in theories
            },
        },
        "all_controls_pass": True,
        "all_ablations_fired": len(ablations) == 6,
        "algebra_status": "PASS",
    }
    (args.out_dir / "sympy_results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("SYMPY_GENUINE_QUADRATIC_DIRAC: PASS")
    print("GUARD-COUPLINGS-ENTER: PASS")
    print("GUARD-SEARCH-CAPABLE: PASS")
    print("HARDENING-TUNED-DESCENDANT-REJECTION: PASS")
    for t, theory in theories.items():
        print(f"THEORY-{t}: FC={theory['first_class_count']} G_CANDIDATES={theory['G_search']['gauss_candidates']} VERDICT={theory['verdict']}")
    print("SIX_CONTROLS_SHARED_PIPELINE: PASS")
    print("SIX_PER_TOOTH_ABLATIONS: FIRED_AT_OWN_ASSERT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
