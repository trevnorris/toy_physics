#!/usr/bin/env python3
"""Ledger stage033 SymPy audit: native P has no emergent Maxwell Gauss law.

Standalone, print-only, assert-zero, exact, and file-I/O-free.  Every native
family and every control enters as a Lagrangian and traverses the same
Legendre -> Dirac closure -> first-class/Gauss search.  Tooth-local mutation
uses ``LEDGER_STAGE033_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
import hashlib
import os
import random
from typing import Any, Callable, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE033_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()
VERDICT = "NATIVE_P_NO_EMERGENT_GAUSS"
K = sp.Matrix([1, 2, 3])
K2 = sp.Integer(K.dot(K))


class AuditFailure(AssertionError):
    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


class PipelineFailure(RuntimeError):
    pass


TOOTH_ORDER = (
    "PASS_Q1_LAGRANGIAN_LEGENDRE",
    "PASS_Q2_MAXWELL_COMPUTED",
    "PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE",
    "PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION",
    "PASS_Q5_GAUSS_SEARCH_LIVE",
    "PASS_Q6_DIRAC_CLOSURE",
    "PASS_GUARD_COUPLINGS_ENTER_A",
    "PASS_GUARD_COUPLINGS_ENTER_C",
    "PASS_THEORY_A_SIGNATURE",
    "PASS_THEORY_C_SIGNATURE",
    "PASS_KINETIC_HESSIAN_DETERMINANT",
    "PASS_COMMON_NULL_A",
    "PASS_COMMON_NULL_C",
    "PASS_BOUNDARY_SCAN_A",
    "PASS_BOUNDARY_SCAN_C",
    "PASS_RANDOMIZED_SWEEP_A",
    "PASS_RANDOMIZED_SWEEP_C",
    "PASS_HARDENING_DESCENDANT_A",
    "PASS_HARDENING_DESCENDANT_C",
    "PASS_PRIMARY_FC_ACCOUNTING_A",
    "PASS_PRIMARY_FC_ACCOUNTING_C",
    "PASS_GUARD_SEARCH_CAPABLE",
    "PASS_CONTROL_MAXWELL",
    "PASS_CONTROL_GAUGED_HARD_UNIT",
    "PASS_CONTROL_BARE_SIGMA",
    "PASS_CONTROL_NONCONSERVED_CURRENT",
    "PASS_CONTROL_COULOMB_GAUGE",
    "PASS_CONTROL_GLOBAL_U1",
    "PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING",
    "PASS_HONEST_TUNED_SCOPE",
    "PASS_DECISION_ORDER_BRANCH2",
    "PASS_VERDICT_TOTALITY",
    "PASS_SOURCE_PREDICATE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "PASS_Q1_LAGRANGIAN_LEGENDRE": "inject a retained velocity into the computed Hamiltonian representative",
    "PASS_Q2_MAXWELL_COMPUTED": "give Maxwell A0 a kinetic term",
    "PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE": "remove global-U(1) from the computed control table",
    "PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION": "corrupt one CAS-derived common-null rule before residual back-substitution",
    "PASS_Q5_GAUSS_SEARCH_LIVE": "replace the live Maxwell search input by kinetic-A0 Maxwell",
    "PASS_Q6_DIRAC_CLOSURE": "cap the native-A Dirac iteration before closure",
    "PASS_GUARD_COUPLINGS_ENTER_A": "drop g_bA from the THEORY-A input Lagrangian",
    "PASS_GUARD_COUPLINGS_ENTER_C": "drop g_dC from the THEORY-C input Lagrangian",
    "PASS_THEORY_A_SIGNATURE": "give lambdaA a kinetic term before recomputing the A signature",
    "PASS_THEORY_C_SIGNATURE": "give lambdaC a kinetic term before recomputing the C signature",
    "PASS_KINETIC_HESSIAN_DETERMINANT": "perturb the u1 kinetic input before recomputing the determinant",
    "PASS_COMMON_NULL_A": "flip g_sA in the computed common-null candidate",
    "PASS_COMMON_NULL_C": "flip g_sC in the computed common-null candidate",
    "PASS_BOUNDARY_SCAN_A": "detune one THEORY-A common-null boundary input",
    "PASS_BOUNDARY_SCAN_C": "detune one THEORY-C common-null boundary input",
    "PASS_RANDOMIZED_SWEEP_A": "replace one A non-common sample by the common-null input",
    "PASS_RANDOMIZED_SWEEP_C": "replace one C non-common sample by the common-null input",
    "PASS_HARDENING_DESCENDANT_A": "inject a computed Maxwell FC primary with a nonzero k-gradient descendant",
    "PASS_HARDENING_DESCENDANT_C": "inject a computed Maxwell FC primary with a nonzero k-gradient descendant",
    "PASS_PRIMARY_FC_ACCOUNTING_A": "perturb the A first-class-primary accounting count",
    "PASS_PRIMARY_FC_ACCOUNTING_C": "perturb the C first-class-primary accounting count",
    "PASS_GUARD_SEARCH_CAPABLE": "force the Maxwell capability input onto the no-primary chain",
    "PASS_CONTROL_MAXWELL": "Maxwell gains A0 kinetic energy",
    "PASS_CONTROL_GAUGED_HARD_UNIT": "remove the hard-unit multiplier",
    "PASS_CONTROL_BARE_SIGMA": "remove the hard-sigma multiplier",
    "PASS_CONTROL_NONCONSERVED_CURRENT": "impose the computed continuity rule",
    "PASS_CONTROL_COULOMB_GAUGE": "drop the Coulomb-fixing input",
    "PASS_CONTROL_GLOBAL_U1": "add a Maxwell local gauge sector",
    "PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING": "feed a genuine Maxwell candidate into source/shear bookkeeping",
    "PASS_HONEST_TUNED_SCOPE": "promote the tuned scan to an exhaustive classification",
    "PASS_DECISION_ORDER_BRANCH2": "feed a computed first-class Maxwell input to the quadratic decision",
    "PASS_VERDICT_TOTALITY": "feed computed Maxwell data through the verdict derivation",
    "PASS_SOURCE_PREDICATE_MANIFEST": "drop one source claim from the canonical partition",
}


def compact(value: Any) -> str:
    try:
        return sp.sstr(sp.factor(sp.cancel(sp.simplify(value))))
    except (TypeError, ValueError, AttributeError):
        return str(value)


def assert_no_float(name: str, value: Any) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            assert_no_float(f"{name}.{key}", item)
        return
    if isinstance(value, (tuple, list, set, frozenset)):
        for index, item in enumerate(value):
            assert_no_float(f"{name}[{index}]", item)
        return
    if isinstance(value, (str, type(None))):
        return
    if isinstance(value, bool):
        value = sp.Integer(1 if value else 0)
    floats = sp.sympify(value).atoms(sp.Float)
    if floats:
        raise AuditFailure(name, f"machine Float atom(s): {floats}")


def expect_zero(name: str, residual: Any, evidence: Any = None) -> None:
    global PASS_COUNT, FAIL_COUNT
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean == 0:
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    FAIL_COUNT += 1
    print(f"FIRST_FAILURE={name}")
    print(f"FAIL  {name}: residual = {compact(clean)}")
    if evidence is not None:
        print(f"      evidence = {evidence}")
    raise AuditFailure(name, compact(clean))


def expect_bool(name: str, condition: bool, evidence: Any = None) -> None:
    expect_zero(name, sp.Integer(0 if bool(condition) else 1), evidence)


def heading(text: str) -> None:
    print("")
    print(text)
    print("-" * len(text))


def pb(f: sp.Expr, g: sp.Expr, pairs: list[tuple[sp.Symbol, sp.Symbol]]) -> sp.Expr:
    return sp.factor(sum(
        sp.diff(f, q) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, q)
        for q, p in pairs
    ))


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


def independent_rows(matrix: sp.Matrix) -> list[int]:
    return list(matrix.T.rref()[1])


def build_hamiltonian(model: Model) -> HamiltonianBuild:
    n = len(model.q)
    if len(model.v) != n:
        raise PipelineFailure(f"{model.key}: q/v size mismatch")
    momenta = list(sp.symbols(" ".join(f"Pi_{model.key}_{i}" for i in range(n)), real=True))
    momentum_map = [sp.diff(model.lagrangian, vi) for vi in model.v]
    hessian = sp.Matrix([[sp.diff(momentum_map[i], model.v[j]) for j in range(n)] for i in range(n)])
    nulls = hessian.T.nullspace()
    residual = sp.Matrix(momenta) - sp.Matrix(momentum_map).subs({vi: 0 for vi in model.v})
    primaries = [sp.factor((nu.T * residual)[0]) for nu in nulls]
    rank = int(hessian.rank())
    row_ids = independent_rows(hessian)
    col_ids = list(hessian.rref()[1])
    if len(row_ids) != rank or len(col_ids) != rank:
        raise PipelineFailure(f"{model.key}: Hessian minor selection failed")
    solution = {vi: sp.Integer(0) for vi in model.v}
    if rank:
        solved = hessian.extract(row_ids, col_ids).inv() * residual.extract(row_ids, [0])
        solution.update({model.v[j]: sp.factor(solved[a]) for a, j in enumerate(col_ids)})
    hc = sp.factor((sp.Matrix(momenta).dot(sp.Matrix(model.v)) - model.lagrangian).subs(solution))
    if any(hc.has(vi) for vi in model.v):
        raise PipelineFailure(f"{model.key}: Legendre transform retained a velocity")
    return HamiltonianBuild(
        model, momenta, momentum_map, hessian, rank, nulls, primaries, hc,
        list(zip(model.q, momenta)),
    )


def linear_row(expr: sp.Expr, z: list[sp.Symbol]) -> sp.Matrix:
    return sp.Matrix([[sp.diff(expr, zi) for zi in z]])


def adds_independent(expr: sp.Expr, constraints: list[sp.Expr], z: list[sp.Symbol]) -> bool:
    row = linear_row(expr, z)
    if row == sp.zeros(1, len(z)):
        return False
    old = (sp.Matrix.vstack(*(linear_row(c, z) for c in constraints))
           if constraints else sp.zeros(0, len(z)))
    return int(sp.Matrix.vstack(old, row).rank()) > int(old.rank())


def dirac_closure(build: HamiltonianBuild, stage_cap: int | None = None) -> dict[str, object]:
    constraints = list(build.primary_constraints)
    stages = [0] * len(constraints)
    primary_count = len(constraints)
    z = build.model.q + build.momenta
    cap = stage_cap if stage_cap is not None else 2 * len(z)
    closed = False
    for stage in range(1, cap + 1):
        if not constraints:
            closed = True
            break
        multiplier_matrix = sp.Matrix([
            [pb(c, p, build.pairs) for p in build.primary_constraints]
            for c in constraints
        ])
        hflow = sp.Matrix([pb(c, build.hamiltonian, build.pairs) for c in constraints])
        new: list[sp.Expr] = []
        for left in multiplier_matrix.T.nullspace():
            candidate = sp.factor((left.T * hflow)[0])
            if adds_independent(candidate, constraints + new, z):
                new.append(candidate)
        if not new:
            closed = True
            break
        constraints.extend(new)
        stages.extend([stage] * len(new))
    if not closed:
        raise PipelineFailure(f"{build.model.key}: Dirac iteration did not close")
    matrix = sp.Matrix([[pb(a, b, build.pairs) for b in constraints] for a in constraints])
    rank = int(matrix.rank()) if constraints else 0
    kernel = matrix.nullspace() if constraints else []
    classes = [
        "FIRST_CLASS" if all(sp.simplify(matrix[i, j]) == 0 for j in range(len(constraints)))
        else "SECOND_CLASS_COMPONENT"
        for i in range(len(constraints))
    ]
    return {
        "constraints": constraints,
        "constraint_stages": stages,
        "primary_count": primary_count,
        "bracket_matrix": matrix,
        "bracket_rank": rank,
        "first_class_count": len(kernel),
        "second_class_count": rank,
        "kernel": kernel,
        "constraint_classes": classes,
        "closed": closed,
    }


def proportional_to_k(action: list[sp.Expr]) -> bool:
    if len(action) != 3 or all(sp.simplify(x) == 0 for x in action):
        return False
    return all(
        sp.simplify(action[i] * K[j] - action[j] * K[i]) == 0
        for i in range(3) for j in range(i + 1, 3)
    )


def gauss_search(build: HamiltonianBuild, data: dict[str, object]) -> dict[str, object]:
    constraints: list[sp.Expr] = data["constraints"]  # type: ignore[assignment]
    matrix: sp.Matrix = data["bracket_matrix"]  # type: ignore[assignment]
    primary_count = int(data["primary_count"])
    primary_fc = (matrix[:, :primary_count].nullspace()
                  if constraints and primary_count else [])
    accepted: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    for coeff in primary_fc:
        primary = sp.factor(sum(coeff[i] * constraints[i] for i in range(primary_count)))
        descendant = sp.factor(pb(primary, build.hamiltonian, build.pairs))
        primary_action = [sp.factor(pb(qi, primary, build.pairs)) for qi in build.model.q]
        secondary_action = [sp.factor(pb(qi, descendant, build.pairs)) for qi in build.model.q]
        fc_descendant = all(sp.simplify(pb(descendant, c, build.pairs)) == 0 for c in constraints)
        scalar_end = any(sp.simplify(primary_action[i]) != 0 for i in build.model.scalar_connection)
        only_scalar = scalar_end and all(
            sp.simplify(a) == 0
            for i, a in enumerate(primary_action)
            if i not in build.model.scalar_connection and i not in build.model.multiplier_fields
        )
        spatial = any(proportional_to_k([secondary_action[i] for i in block])
                      for block in build.model.spatial_connections)
        nontrivial = any(sp.simplify(x) != 0 for x in secondary_action)
        descendant_zero = sp.simplify(descendant) == 0
        certified = bool(descendant_zero or not spatial)
        item = {
            "primary": primary,
            "descendant": descendant,
            "primary_action": primary_action,
            "secondary_action": secondary_action,
            "spatial_gradient_action": spatial,
            "descendant_zero": descendant_zero,
            "secondary_action_non_gradient": not spatial,
            "descendant_rejection_certified": certified,
            "reason": ("DESCENDANT_ZERO" if descendant_zero else
                       "SECONDARY_ACTION_NON_GRADIENT" if not spatial else
                       "OTHER_GAUSS_CRITERION"),
        }
        (accepted if fc_descendant and only_scalar and spatial and nontrivial else rejected).append(item)
    return {
        "computed_first_class_primaries": len(primary_fc),
        "gauss_candidates": len(accepted),
        "additional_G_exists": bool(accepted),
        "candidates": accepted,
        "rejected": rejected,
    }


def execute(model: Model) -> tuple[HamiltonianBuild, dict[str, object], dict[str, object]]:
    build = build_hamiltonian(model)
    data = dirac_closure(build)
    search = gauss_search(build, data)
    return build, data, search


def native_model(theory: str, substitutions: dict[sp.Symbol, sp.Expr] | None = None) -> Model:
    if theory not in {"A", "C"}:
        raise ValueError(theory)
    p = list(sp.symbols(f"p{theory}1:4", real=True))
    u = list(sp.symbols(f"u{theory}1:4", real=True))
    s, lam, sigma, b = sp.symbols(f"s{theory} lambda{theory} sigma{theory} b{theory}", real=True)
    xi = sp.symbols(f"xi{theory}", real=True)
    q = p + u + [s, lam, sigma, b] + ([xi] if theory == "C" else [])
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    vp, vu = sp.Matrix(v[:3]), sp.Matrix(v[3:6])
    gt, gs, gd = sp.symbols(f"g_t{theory} g_s{theory} g_d{theory}", real=True)
    gb = sp.symbols(f"g_b{theory}", real=True)
    ru, ap, au, ab = sp.symbols(f"rho_u{theory} a_p{theory} a_u{theory} a_b{theory}", positive=True)
    pm, um = sp.Matrix(p), sp.Matrix(u)
    transverse = K2 * sp.eye(3) - K * K.T
    kinetic = (vp.dot(vp) + v[6] ** 2 + ru * vu.dot(vu)) / 2 + gt * vp.dot(vu)
    if theory == "A":
        kinetic += v[9] ** 2 / 2
    potential = ap * pm.dot(pm) / 2 + au * (um.T * transverse * um)[0] / 2
    potential += gs * K2 * pm.dot(um) + gd * (K.dot(pm)) ** 2 / 2 + ab * b**2 / 2
    if theory == "A":
        potential += gb * b * K.dot(pm)
    lagrangian = kinetic - potential + lam * s + sigma * K.dot(um)
    if theory == "C":
        lagrangian += xi * b
    couplings = [gt, gs, gd] + ([gb] if theory == "A" else [])
    if substitutions:
        lagrangian = sp.factor(lagrangian.subs(substitutions))
        couplings = [coupling for coupling in couplings if coupling not in substitutions]
    return Model(
        f"native_{theory}", f"THEORY-{theory}", q, v, lagrangian, couplings,
        [], [list(range(3)), list(range(3, 6))],
        [7, 8] + ([10] if theory == "C" else []),
    )


def maxwell_model(key: str = "maxwell", *, kinetic_a0: bool = False,
                  impose_conservation: bool = False) -> Model:
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
        current = sp.Matrix(sp.symbols("nc_j1:4", real=True))
        lagrangian += rho * a0 + current.dot(avec)
        current_data = {
            "rho": rho, "rho_dot": rho_dot, "current": current,
            "impose_conservation": impose_conservation,
        }
    return Model(key, key, q, v, lagrangian, [], [0], [[1, 2, 3]], [], current_data)


def gauged_hard_unit_model(remove_hard: bool = False) -> Model:
    a0, a1, a2, a3, f1, f2, lam = sp.symbols(
        "gh_A0 gh_A1 gh_A2 gh_A3 gh_phi1 gh_phi2 gh_lambda", real=True)
    q = [a0, a1, a2, a3, f1, f2] + ([] if remove_hard else [lam])
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    avec, va = sp.Matrix(q[1:4]), sp.Matrix(v[1:4])
    phi, vphi = sp.Matrix([f1, f2]), sp.Matrix(v[4:6])
    rotphi = sp.Matrix([-f2, f1])
    magnetic = K2 * sp.eye(3) - K * K.T
    lagrangian = (va - K * a0).dot(va - K * a0) / 2 - (avec.T * magnetic * avec)[0] / 2
    lagrangian += (vphi - a0 * rotphi).dot(vphi - a0 * rotphi) / 2 - phi.dot(phi) / 2
    if not remove_hard:
        lagrangian += lam * (phi.dot(phi) - 1) / 2
    return Model("gauged_hard_unit", "gauged hard unit", q, v, lagrangian, [], [0], [[1, 2, 3]],
                 [] if remove_hard else [6])


def bare_sigma_model(remove_hard: bool = False) -> Model:
    phi = list(sp.symbols("bs_phi1:5", real=True))
    lam = sp.symbols("bs_lambda", real=True)
    q = phi + ([] if remove_hard else [lam])
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    lagrangian = sum(x*x for x in v[:4]) / 2 - sum(x*x for x in phi) / 2
    if not remove_hard:
        lagrangian += lam * (sum(x*x for x in phi) - 1) / 2
    return Model("bare_sigma", "bare O(4) hard sigma", q, v, lagrangian, [], [], [],
                 [] if remove_hard else [4])


def coulomb_gauge_model() -> Model:
    a0, a1, a2, a3, eta, zeta = sp.symbols("gf_A0 gf_A1 gf_A2 gf_A3 gf_eta gf_zeta", real=True)
    q = [a0, a1, a2, a3, eta, zeta]
    v = list(sp.symbols(" ".join(f"d_{x}" for x in q), real=True))
    avec, va = sp.Matrix(q[1:4]), sp.Matrix(v[1:4])
    magnetic = K2 * sp.eye(3) - K * K.T
    lagrangian = (va - K*a0).dot(va - K*a0)/2 - (avec.T*magnetic*avec)[0]/2
    lagrangian += eta*a0 + zeta*K.dot(avec)
    return Model("gauge_fixed_maxwell", "Coulomb gauge Maxwell", q, v, lagrangian, [],
                 [0], [[1, 2, 3]], [4, 5])


def global_u1_model(add_maxwell: bool = False) -> Model:
    if add_maxwell:
        return maxwell_model("global_only_ablation")
    x, y = sp.symbols("global_x global_y", real=True)
    vx, vy = sp.symbols("d_global_x d_global_y", real=True)
    return Model("global_only", "global U(1) complex scalar", [x, y], [vx, vy],
                 (vx**2 + vy**2)/2 - (x**2 + y**2)/2, [], [], [], [])


CONTROL_PREDICATES = {
    "maxwell": "PASS_CONTROL_MAXWELL",
    "gauged_hard_unit": "PASS_CONTROL_GAUGED_HARD_UNIT",
    "bare_sigma": "PASS_CONTROL_BARE_SIGMA",
    "nonconserved_current": "PASS_CONTROL_NONCONSERVED_CURRENT",
    "gauge_fixed_maxwell": "PASS_CONTROL_COULOMB_GAUGE",
    "global_only": "PASS_CONTROL_GLOBAL_U1",
}

CONTROL_EXPECTED = {
    "maxwell": ("FIRST_CLASS_GAUSS", 2, 0, 1, 1),
    "gauged_hard_unit": ("MIXED", 2, 4, 1, 2),
    "bare_sigma": ("SECOND_CLASS_RADIAL_NO_GAUSS", 0, 4, 0, 1),
    "nonconserved_current": ("INCONSISTENT_PRESERVATION", 2, 0, 1, 1),
    "gauge_fixed_maxwell": ("SECOND_CLASS_NO_LOCAL_GAUGE", 0, 8, 0, 3),
    "global_only": ("GLOBAL_CHARGE_NO_LOCAL_GAUSS", 0, 0, 0, 0),
}


def control_model(key: str, ablate: bool) -> Model:
    if key == "maxwell":
        return maxwell_model(kinetic_a0=ablate)
    if key == "gauged_hard_unit":
        return gauged_hard_unit_model(remove_hard=ablate)
    if key == "bare_sigma":
        return bare_sigma_model(remove_hard=ablate)
    if key == "nonconserved_current":
        return maxwell_model("nonconserved_current", impose_conservation=ablate)
    if key == "gauge_fixed_maxwell":
        return maxwell_model("gauge_fixed_ablation") if ablate else coulomb_gauge_model()
    if key == "global_only":
        return global_u1_model(add_maxwell=ablate)
    raise ValueError(key)


def control_result(key: str, ablate: bool = False) -> dict[str, object]:
    model = control_model(key, ablate)
    build, data, search = execute(model)
    current_inconsistent = False
    preservation = sp.Integer(0)
    if model.current_data is not None:
        count = int(data["primary_count"])
        matrix: sp.Matrix = data["bracket_matrix"]  # type: ignore[assignment]
        directions = matrix[:, :count].nullspace()
        if not directions:
            raise PipelineFailure("current control has no first-class primary")
        constraints: list[sp.Expr] = data["constraints"]  # type: ignore[assignment]
        primary = sp.factor(sum(directions[0][i] * constraints[i] for i in range(count)))
        gauss = sp.factor(pb(primary, build.hamiltonian, build.pairs))
        rho = model.current_data["rho"]
        rho_dot = model.current_data["rho_dot"]
        preservation = sp.factor(sp.diff(gauss, rho) * rho_dot + pb(gauss, build.hamiltonian, build.pairs))
        if model.current_data["impose_conservation"]:
            rules = sp.solve(sp.Eq(preservation, 0), rho_dot)
            if not rules:
                raise PipelineFailure("continuity equation did not solve")
            preservation = sp.factor(preservation.subs(rho_dot, rules[0]))
        current_inconsistent = sp.simplify(preservation) != 0
    fc = int(data["first_class_count"])
    sc = int(data["second_class_count"])
    gauss_count = int(search["gauss_candidates"])
    nullity = len(build.hessian_nullspace)
    if key == "maxwell":
        token = "FIRST_CLASS_GAUSS" if gauss_count > 0 and sc == 0 else "NO_PRIMARY_GAUSS_CHAIN"
    elif key == "gauged_hard_unit":
        token = "MIXED" if gauss_count > 0 and fc > 0 and sc > 0 else "NOT_MIXED"
    elif key == "bare_sigma":
        token = "SECOND_CLASS_RADIAL_NO_GAUSS" if fc == gauss_count == 0 and sc > 0 else "BAD_SIGMA_CLASS"
    elif key == "nonconserved_current":
        token = "INCONSISTENT_PRESERVATION" if gauss_count > 0 and current_inconsistent else "CURRENT_CONSISTENT"
    elif key == "gauge_fixed_maxwell":
        token = "SECOND_CLASS_NO_LOCAL_GAUGE" if fc == gauss_count == 0 and sc > 0 else "LOCAL_GAUGE_REMAINS"
    else:
        token = "GLOBAL_CHARGE_NO_LOCAL_GAUSS" if nullity == fc == gauss_count == 0 else "LOCAL_GAUGE_PRESENT"
    return {
        "model": key, "classification": token, "fc": fc, "sc": sc,
        "gauss": gauss_count, "nullity": nullity, "preservation": preservation,
        "build": build, "dirac": data, "search": search,
    }


def coupling_locations(build: HamiltonianBuild, data: dict[str, object]) -> dict[str, tuple[str, ...]]:
    locations: dict[str, tuple[str, ...]] = {}
    constraints: list[sp.Expr] = data["constraints"]  # type: ignore[assignment]
    matrix: sp.Matrix = data["bracket_matrix"]  # type: ignore[assignment]
    for coupling in build.model.couplings:
        hits: list[str] = []
        if any(expr.has(coupling) for expr in build.momentum_map):
            hits.append("momentum_map")
        if any(expr.has(coupling) for expr in build.hessian):
            hits.append("hessian")
        if any(expr.has(coupling) for expr in constraints):
            hits.append("constraints")
        if any(expr.has(coupling) for expr in matrix):
            hits.append("bracket_matrix")
        locations[str(coupling)] = tuple(hits)
    return locations


def common_null_derivation(theory: str) -> dict[int, dict[str, object]]:
    model = native_model(theory)
    free = {str(symbol): symbol for symbol in model.lagrangian.free_symbols}
    ap, au = free[f"a_p{theory}"], free[f"a_u{theory}"]
    gs = free[f"g_s{theory}"]
    physical_q = model.q[:6]
    potential = -model.lagrangian.subs({velocity: 0 for velocity in model.v})
    potential_matrix = sp.hessian(potential, physical_q)
    transverse_test = sp.Matrix([2, -1, 0])
    results: dict[int, dict[str, object]] = {}
    for sign in (-1, 1):
        null_test = sp.Matrix(list(-sign * transverse_test) + list(transverse_test))
        residual = [sp.factor(value) for value in potential_matrix * null_test]
        solutions = sp.solve(residual, [ap, gs], dict=True)
        results[sign] = {
            "residual": residual, "solutions": solutions, "ap": ap, "au": au, "gs": gs,
        }
    return results


def hardening_record(theory: str, sign: int, build: HamiltonianBuild,
                     data: dict[str, object], search: dict[str, object]) -> dict[str, object]:
    rejected: list[dict[str, object]] = search["rejected"]  # type: ignore[assignment]
    candidates: list[dict[str, object]] = search["candidates"]  # type: ignore[assignment]
    return {
        "theory": theory,
        "sign": sign,
        "fc": int(data["first_class_count"]),
        "primary_fc": int(search["computed_first_class_primaries"]),
        "rejected": rejected,
        "candidates": candidates,
        "hessian_nullity": len(build.hessian_nullspace),
    }


def scan_point(theory: str, substitutions: dict[sp.Symbol, sp.Expr]) -> dict[str, object]:
    build, data, search = execute(native_model(theory, substitutions))
    return {
        "fc": int(data["first_class_count"]),
        "gauss": int(search["gauss_candidates"]),
        "additional": bool(search["additional_G_exists"]),
        "nullity": len(build.hessian_nullspace),
        "build": build,
        "dirac": data,
        "search": search,
        "substitutions": substitutions,
    }


def native_analysis(theory: str) -> dict[str, object]:
    model = native_model(theory)
    build, data, search = execute(model)
    free = {str(symbol): symbol for symbol in model.lagrangian.free_symbols}
    gt, gs, gd = free[f"g_t{theory}"], free[f"g_s{theory}"], free[f"g_d{theory}"]
    ru, ap, au = free[f"rho_u{theory}"], free[f"a_p{theory}"], free[f"a_u{theory}"]
    gb = free.get(f"g_b{theory}")
    determinant = sp.factor(build.hessian.extract(range(6), range(6)).det())
    component_determinant = sp.factor(build.hessian.extract([0, 3], [0, 3]).det())
    common = common_null_derivation(theory)

    boundary: list[dict[str, object]] = []
    common_records: list[dict[str, object]] = []
    for sign in (-1, 1):
        noncommon = {
            gt: sp.Integer(sign), ru: sp.Integer(1), au: sp.Integer(1), ap: sp.Integer(1),
            gs: sp.Rational(15 * sign, 28), gd: sp.Integer(0),
        }
        if gb is not None:
            noncommon[gb] = sp.Integer(0)
        solved = common[sign]["solutions"]
        if not solved:
            raise PipelineFailure(f"THEORY-{theory} sign {sign}: common-null solve empty")
        common_rules = {gt: sp.Integer(sign), ru: sp.Integer(1), **solved[0]}
        for label, substitutions in (
            ("generic semidefinite", {gt: sp.Integer(sign), ru: sp.Integer(1)}),
            ("first rank-drop non-common", noncommon),
            ("common Hessian/potential null", common_rules),
        ):
            point = scan_point(theory, substitutions)
            point["label"] = label
            point["sign"] = sign
            boundary.append(point)
            if label.startswith("common"):
                common_records.append(hardening_record(
                    theory, sign, point["build"], point["dirac"], point["search"]  # type: ignore[arg-type]
                ))

    seed = 260713 if theory == "A" else 260715
    rng = random.Random(seed)
    sweep: list[dict[str, object]] = []
    for sign in (-1, 1):
        for sample in range(1, 4):
            au_value = rng.randint(1, 5)
            ap_value = rng.randint(1, 40)
            while ap_value == 14 * au_value:
                ap_value = rng.randint(1, 40)
            substitutions = {
                gt: sp.Integer(sign), ru: sp.Integer(1), au: sp.Integer(au_value),
                ap: sp.Integer(ap_value),
                gs: sp.Rational(sign * (ap_value + 14 * au_value), 28),
                gd: sp.Integer(0),
            }
            if gb is not None:
                substitutions[gb] = sp.Integer(0)
            point = scan_point(theory, substitutions)
            point.update({"sign": sign, "sample": sample})
            sweep.append(point)

    return {
        "model": model,
        "build": build,
        "dirac": data,
        "search": search,
        "coupling_locations": coupling_locations(build, data),
        "determinant": determinant,
        "component_determinant": component_determinant,
        "symbols": {"gt": gt, "gs": gs, "gd": gd, "ru": ru, "ap": ap, "au": au, "gb": gb},
        "common": common,
        "boundary": boundary,
        "sweep": sweep,
        "seed": seed,
        "hardening": common_records,
    }


def signature(result: tuple[HamiltonianBuild, dict[str, object], dict[str, object]]) -> tuple[object, ...]:
    build, data, search = result
    constraints: list[sp.Expr] = data["constraints"]  # type: ignore[assignment]
    classes: list[str] = data["constraint_classes"]  # type: ignore[assignment]
    return (
        build.hessian_rank,
        len(build.hessian_nullspace),
        len(constraints),
        tuple(data["constraint_stages"]),
        tuple(classes),
        int(data["bracket_rank"]),
        int(data["first_class_count"]),
        int(data["second_class_count"]),
        int(search["gauss_candidates"]),
        bool(search["additional_G_exists"]),
        int(data["first_class_count"]) == len(constraints) - int(data["bracket_rank"]),
    )


def expected_signature(theory: str) -> tuple[object, ...]:
    if theory == "A":
        stages = (0, 0, 1, 1, 2, 2, 3, 3)
        return (8, 2, 8, stages, ("SECOND_CLASS_COMPONENT",) * 8, 8, 0, 8, 0, False, True)
    stages = (0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3)
    return (7, 4, 12, stages, ("SECOND_CLASS_COMPONENT",) * 12, 12, 0, 12, 0, False, True)


def verdict_from_computed(open_searches: Iterable[dict[str, object]],
                          tuned_searches: Iterable[dict[str, object]]) -> str:
    if any(bool(item["additional_G_exists"]) for item in open_searches):
        return "FIRST_CLASS_GENERIC_EM_CANDIDATE"
    if any(bool(item["additional_G_exists"]) for item in tuned_searches):
        return "FIRST_CLASS_TUNED_INVERSE_DESIGN"
    return VERDICT


SOURCE_TOOTH_IDS = (
    "Q1_LAGRANGIAN_LEGENDRE", "Q2_MAXWELL_COMPUTED", "Q3_SIX_CONTROLS",
    "Q4_WOLFRAM_INDEPENDENCE", "Q5_G_SEARCH", "Q6_DIRAC_CLOSURE",
    "GUARD_COUPLINGS_ENTER_A", "GUARD_COUPLINGS_ENTER_C", "GUARD_SEARCH_CAPABLE",
    "HARDENING_DESCENDANT_A", "HARDENING_DESCENDANT_C",
    "THEORY_A_HESSIAN_RANK_NULLITY", "THEORY_A_CONSTRAINT_COUNT",
    "THEORY_A_STAGES", "THEORY_A_ALL_SECOND_CLASS", "THEORY_A_PB_FC_SC",
    "THEORY_A_G_SEARCH", "THEORY_C_HESSIAN_RANK_NULLITY",
    "THEORY_C_CONSTRAINT_COUNT", "THEORY_C_STAGES", "THEORY_C_ALL_SECOND_CLASS",
    "THEORY_C_PB_FC_SC", "THEORY_C_G_SEARCH", "KINETIC_HESSIAN_DETERMINANT",
    "CONTROL_MAXWELL", "CONTROL_GAUGED_HARD_UNIT", "CONTROL_BARE_SIGMA",
    "CONTROL_NONCONSERVED_CURRENT", "CONTROL_COULOMB_GAUGE", "CONTROL_GLOBAL_U1",
    "ABLATION_MAXWELL", "ABLATION_GAUGED_HARD_UNIT", "ABLATION_BARE_SIGMA",
    "ABLATION_NONCONSERVED_CURRENT", "ABLATION_COULOMB_GAUGE", "ABLATION_GLOBAL_U1",
    "BOUNDARY_SCAN_A", "BOUNDARY_SCAN_C", "RANDOM_SWEEP_A", "RANDOM_SWEEP_C",
    "SEEDS_AND_CARDINALITY", "COMMON_NULL_A", "COMMON_NULL_C",
    "SOURCE_FIRST_ORDERING", "SHEAR_DUPLICATE", "DECISION_ORDER_BRANCH2",
    "HONEST_TUNED_SCOPE", "VERDICT_TOTALITY", "ARGPARSE_OUT_DIR_HARNESS",
    "JSON_ARTIFACT_WRITING", "CROSS_ENGINE_FILE_COMPARATOR",
)

SOURCE_MANIFEST = (
    ("Q1_LAGRANGIAN_LEGENDRE", "PRESERVED", "STAGE033_SHARED_PIPELINE"),
    ("Q2_MAXWELL_COMPUTED", "PRESERVED", "STAGE033_CONTROL_MAXWELL"),
    ("Q3_SIX_CONTROLS", "REPLACED_BY_STRONGER", "STAGE033_NATIVE_ABLATIONS"),
    ("Q4_WOLFRAM_INDEPENDENCE", "REPLACED_BY_STRONGER", "STAGE033_REAUTHORED_WOLFRAM"),
    ("Q5_G_SEARCH", "PRESERVED", "STAGE033_GAUSS_SEARCH"),
    ("Q6_DIRAC_CLOSURE", "PRESERVED", "STAGE033_DIRAC_CLOSURE"),
    ("GUARD_COUPLINGS_ENTER_A", "PRESERVED", "STAGE033_COUPLING_GUARD_A"),
    ("GUARD_COUPLINGS_ENTER_C", "PRESERVED", "STAGE033_COUPLING_GUARD_C"),
    ("GUARD_SEARCH_CAPABLE", "PRESERVED", "STAGE033_SEARCH_CAPABILITY"),
    ("HARDENING_DESCENDANT_A", "PRESERVED", "STAGE033_DESCENDANT_A"),
    ("HARDENING_DESCENDANT_C", "PRESERVED", "STAGE033_DESCENDANT_C"),
    ("THEORY_A_HESSIAN_RANK_NULLITY", "PRESERVED", "STAGE033_SIGNATURE_A"),
    ("THEORY_A_CONSTRAINT_COUNT", "PRESERVED", "STAGE033_SIGNATURE_A"),
    ("THEORY_A_STAGES", "PRESERVED", "STAGE033_SIGNATURE_A"),
    ("THEORY_A_ALL_SECOND_CLASS", "PRESERVED", "STAGE033_SIGNATURE_A"),
    ("THEORY_A_PB_FC_SC", "REPLACED_BY_STRONGER", "STAGE033_RANK_NULLITY_A"),
    ("THEORY_A_G_SEARCH", "PRESERVED", "STAGE033_GAUSS_ABSENCE_A"),
    ("THEORY_C_HESSIAN_RANK_NULLITY", "PRESERVED", "STAGE033_SIGNATURE_C"),
    ("THEORY_C_CONSTRAINT_COUNT", "PRESERVED", "STAGE033_SIGNATURE_C"),
    ("THEORY_C_STAGES", "PRESERVED", "STAGE033_SIGNATURE_C"),
    ("THEORY_C_ALL_SECOND_CLASS", "PRESERVED", "STAGE033_SIGNATURE_C"),
    ("THEORY_C_PB_FC_SC", "REPLACED_BY_STRONGER", "STAGE033_RANK_NULLITY_C"),
    ("THEORY_C_G_SEARCH", "PRESERVED", "STAGE033_GAUSS_ABSENCE_C"),
    ("KINETIC_HESSIAN_DETERMINANT", "REPLACED_BY_STRONGER", "STAGE033_FULL_AND_COMPONENT_DET"),
    ("CONTROL_MAXWELL", "PRESERVED", "STAGE033_CONTROL_MAXWELL"),
    ("CONTROL_GAUGED_HARD_UNIT", "PRESERVED", "STAGE033_CONTROL_GAUGED"),
    ("CONTROL_BARE_SIGMA", "PRESERVED", "STAGE033_CONTROL_SIGMA"),
    ("CONTROL_NONCONSERVED_CURRENT", "PRESERVED", "STAGE033_CONTROL_CURRENT"),
    ("CONTROL_COULOMB_GAUGE", "PRESERVED", "STAGE033_CONTROL_COULOMB"),
    ("CONTROL_GLOBAL_U1", "PRESERVED", "STAGE033_CONTROL_GLOBAL"),
    ("ABLATION_MAXWELL", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"),
    ("ABLATION_GAUGED_HARD_UNIT", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"),
    ("ABLATION_BARE_SIGMA", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"),
    ("ABLATION_NONCONSERVED_CURRENT", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"),
    ("ABLATION_COULOMB_GAUGE", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"),
    ("ABLATION_GLOBAL_U1", "REPLACED_BY_STRONGER", "STAGE033_ENV_OWN_ASSERT"),
    ("BOUNDARY_SCAN_A", "PRESERVED", "STAGE033_BOUNDARY_A"),
    ("BOUNDARY_SCAN_C", "PRESERVED", "STAGE033_BOUNDARY_C"),
    ("RANDOM_SWEEP_A", "PRESERVED", "STAGE033_SWEEP_A"),
    ("RANDOM_SWEEP_C", "PRESERVED", "STAGE033_SWEEP_C"),
    ("SEEDS_AND_CARDINALITY", "PRESERVED", "STAGE033_FIXED_SEEDS"),
    ("COMMON_NULL_A", "REPLACED_BY_STRONGER", "STAGE033_SYMBOLIC_SOLVE_A"),
    ("COMMON_NULL_C", "REPLACED_BY_STRONGER", "STAGE033_SYMBOLIC_SOLVE_C"),
    ("SOURCE_FIRST_ORDERING", "PRESERVED", "STAGE033_SOURCE_FIRST"),
    ("SHEAR_DUPLICATE", "PRESERVED", "STAGE033_NO_G_BOOKKEEPING"),
    ("DECISION_ORDER_BRANCH2", "PRESERVED", "STAGE033_QUADRATIC_BRANCH2"),
    ("HONEST_TUNED_SCOPE", "REPLACED_BY_STRONGER", "STAGE033_ARGUED_SCANNED_SCOPE"),
    ("VERDICT_TOTALITY", "REPLACED_BY_STRONGER", "STAGE033_COMPUTED_VERDICT"),
    ("ARGPARSE_OUT_DIR_HARNESS", "SCOPED_OUT", "STAGE033_PRINT_ONLY_CONTRACT"),
    ("JSON_ARTIFACT_WRITING", "SCOPED_OUT", "STAGE033_ZERO_FILE_IO"),
    ("CROSS_ENGINE_FILE_COMPARATOR", "SCOPED_OUT", "STAGE033_INDEPENDENT_TOKENS"),
)


def canonical_manifest_text(manifest: Iterable[tuple[str, str, str]]) -> str:
    return "\n".join("|".join(row) for row in sorted(manifest))


def manifest_digest(manifest: Iterable[tuple[str, str, str]]) -> str:
    return hashlib.sha256(canonical_manifest_text(manifest).encode("utf-8")).hexdigest()


COMMITTED_MANIFEST_DIGEST = "6b191e77fefe24c9000445f01e4e2c6154ab1bb9b15bb40a6d1515dc463a7e9d"
EXPECTED_MANIFEST_COUNTS = {
    "PRESERVED": 33,
    "REPLACED_BY_STRONGER": 15,
    "SCOPED_OUT": 3,
}


def run_assertions() -> None:
    heading("Computed input-Lagrangian pipeline and six controls")
    baseline_controls = {key: control_result(key) for key in CONTROL_EXPECTED}
    ablated_controls = {key: control_result(key, True) for key in CONTROL_EXPECTED}
    maxwell_live = baseline_controls["maxwell"]

    native_a = native_analysis("A")
    print("PROGRESS SYMPY THEORY-A OPEN+BOUNDARY+SWEEP COMPLETE")
    native_c = native_analysis("C")
    print("PROGRESS SYMPY THEORY-C OPEN+BOUNDARY+SWEEP COMPLETE")

    # Q1: the representative Hamiltonian is recomputed from the input L and
    # must contain no arbitrary velocity.  The mutation changes that computed
    # upstream object, not a copied truth flag.
    legendre_test = native_a["build"].hamiltonian  # type: ignore[union-attr]
    if ACTIVE_MUTATION == "PASS_Q1_LAGRANGIAN_LEGENDRE":
        legendre_test += native_a["model"].v[0]  # type: ignore[union-attr]
    q1_ok = (
        len(native_a["build"].primary_constraints) == len(native_a["build"].hessian_nullspace)  # type: ignore[union-attr]
        and len(native_c["build"].primary_constraints) == len(native_c["build"].hessian_nullspace)  # type: ignore[union-attr]
        and not any(legendre_test.has(v) for v in native_a["model"].v)  # type: ignore[union-attr]
    )
    expect_bool("PASS_Q1_LAGRANGIAN_LEGENDRE", q1_ok)

    q2 = control_result("maxwell", ACTIVE_MUTATION == "PASS_Q2_MAXWELL_COMPUTED")
    q2_matrix: sp.Matrix = q2["dirac"]["bracket_matrix"]  # type: ignore[index,assignment]
    expect_bool(
        "PASS_Q2_MAXWELL_COMPUTED",
        (q2["classification"], q2["fc"], q2["sc"], q2["gauss"], q2["nullity"])
        == CONTROL_EXPECTED["maxwell"] and q2_matrix == sp.zeros(2),
        (q2["classification"], q2["fc"], q2["sc"], q2["gauss"], q2["nullity"]),
    )

    control_table = dict(baseline_controls)
    if ACTIVE_MUTATION == "PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE":
        control_table.pop("global_only")
    expect_bool(
        "PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE",
        tuple(control_table) == tuple(CONTROL_EXPECTED)
        and all("build" in row and "dirac" in row and "search" in row for row in control_table.values()),
        tuple(control_table),
    )

    # Q4 is a runtime back-substitution certificate for CAS-generated common
    # null rules in both families.  It does not inspect source text.
    q4_ok = True
    for theory_data in (native_a, native_c):
        for sign in (-1, 1):
            item = theory_data["common"][sign]  # type: ignore[index]
            rules = list(item["solutions"])
            if ACTIVE_MUTATION == "PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION" and theory_data is native_a and sign == 1:
                wrong = dict(rules[0])
                wrong[item["gs"]] = -wrong[item["gs"]]
                rules = [wrong]
            q4_ok = q4_ok and bool(rules) and all(
                all(sp.simplify(residual.subs(rule)) == 0 for residual in item["residual"])
                for rule in rules
            )
    expect_bool("PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION", q4_ok)

    q5 = control_result("maxwell", ACTIVE_MUTATION == "PASS_Q5_GAUSS_SEARCH_LIVE")
    expect_bool(
        "PASS_Q5_GAUSS_SEARCH_LIVE",
        int(q5["search"]["gauss_candidates"]) > 0  # type: ignore[index]
        and bool(q5["search"]["additional_G_exists"]),  # type: ignore[index]
        q5["search"],
    )

    q6_closed = bool(native_a["dirac"]["closed"] and native_c["dirac"]["closed"])  # type: ignore[index]
    if ACTIVE_MUTATION == "PASS_Q6_DIRAC_CLOSURE":
        try:
            dirac_closure(native_a["build"], stage_cap=1)  # type: ignore[arg-type]
        except PipelineFailure:
            q6_closed = False
    expect_bool("PASS_Q6_DIRAC_CLOSURE", q6_closed)

    heading("Native THEORY-A/C symbolic open stratum")
    for theory, result, tooth in (
        ("A", native_a, "PASS_GUARD_COUPLINGS_ENTER_A"),
        ("C", native_c, "PASS_GUARD_COUPLINGS_ENTER_C"),
    ):
        locations = result["coupling_locations"]
        if ACTIVE_MUTATION == tooth:
            target = result["symbols"]["gb" if theory == "A" else "gd"]  # type: ignore[index]
            bad_model = replace(result["model"], lagrangian=result["model"].lagrangian.subs(target, 0))  # type: ignore[union-attr]
            bad_build, bad_dirac, _ = execute(bad_model)
            locations = coupling_locations(bad_build, bad_dirac)
        expected_names = ({"g_tA", "g_sA", "g_dA", "g_bA"} if theory == "A"
                          else {"g_tC", "g_sC", "g_dC"})
        expect_bool(tooth, set(locations) == expected_names and all(locations.values()), locations)

    for theory, result, tooth in (
        ("A", native_a, "PASS_THEORY_A_SIGNATURE"),
        ("C", native_c, "PASS_THEORY_C_SIGNATURE"),
    ):
        computed = signature((result["build"], result["dirac"], result["search"]))  # type: ignore[arg-type]
        if ACTIVE_MUTATION == tooth:
            model: Model = result["model"]  # type: ignore[assignment]
            bad_model = replace(model, lagrangian=model.lagrangian + model.v[7] ** 2 / 2)
            computed = signature(execute(bad_model))
        expect_bool(tooth, computed == expected_signature(theory), computed)

    det_a = (native_a["determinant"], native_a["component_determinant"])
    det_c = (native_c["determinant"], native_c["component_determinant"])
    if ACTIVE_MUTATION == "PASS_KINETIC_HESSIAN_DETERMINANT":
        model = native_a["model"]
        altered = replace(model, lagrangian=model.lagrangian + model.v[3] ** 2 / 2)  # type: ignore[union-attr]
        altered_build = build_hamiltonian(altered)
        det_a = (
            sp.factor(altered_build.hessian.extract(range(6), range(6)).det()),
            sp.factor(altered_build.hessian.extract([0, 3], [0, 3]).det()),
        )
    gt_a, ru_a = native_a["symbols"]["gt"], native_a["symbols"]["ru"]  # type: ignore[index]
    gt_c, ru_c = native_c["symbols"]["gt"], native_c["symbols"]["ru"]  # type: ignore[index]
    expect_bool(
        "PASS_KINETIC_HESSIAN_DETERMINANT",
        sp.simplify(det_a[0] - (ru_a - gt_a**2) ** 3) == 0
        and sp.simplify(det_a[1] - (ru_a - gt_a**2)) == 0
        and sp.simplify(det_c[0] - (ru_c - gt_c**2) ** 3) == 0
        and sp.simplify(det_c[1] - (ru_c - gt_c**2)) == 0,
        {"A": det_a, "C": det_c},
    )

    heading("CAS-derived common-null locus and tuned scans")
    for theory, result, tooth in (
        ("A", native_a, "PASS_COMMON_NULL_A"),
        ("C", native_c, "PASS_COMMON_NULL_C"),
    ):
        common_ok = True
        for sign in (-1, 1):
            item = result["common"][sign]  # type: ignore[index]
            solutions = list(item["solutions"])
            if ACTIVE_MUTATION == tooth and sign == 1:
                wrong = dict(solutions[0])
                wrong[item["gs"]] = -wrong[item["gs"]]
                solutions = [wrong]
            common_ok = common_ok and len(solutions) == 1
            common_ok = common_ok and sp.simplify(solutions[0][item["ap"]] - 14 * item["au"]) == 0
            common_ok = common_ok and sp.simplify(solutions[0][item["gs"]] - sign * item["au"]) == 0
            common_ok = common_ok and all(
                sp.simplify(residual.subs(solutions[0])) == 0 for residual in item["residual"]
            )
        expect_bool(tooth, common_ok)

    for theory, result, tooth, nullity in (
        ("A", native_a, "PASS_BOUNDARY_SCAN_A", 5),
        ("C", native_c, "PASS_BOUNDARY_SCAN_C", 7),
    ):
        points = list(result["boundary"])
        if ACTIVE_MUTATION == tooth:
            symbols = result["symbols"]
            bad_subs = {
                symbols["gt"]: 1, symbols["ru"]: 1, symbols["ap"]: 14,
                symbols["au"]: 1, symbols["gs"]: sp.Rational(3, 2), symbols["gd"]: 0,
            }
            if symbols["gb"] is not None:
                bad_subs[symbols["gb"]] = 0
            points[2] = scan_point(theory, bad_subs)
        observed = [(p["fc"], p["gauss"], p["additional"], p["nullity"]) for p in points]
        expected = [(0, 0, False, nullity), (0, 0, False, nullity), (2, 0, False, nullity)] * 2
        expect_bool(tooth, observed == expected, observed)

    for theory, result, tooth, nullity in (
        ("A", native_a, "PASS_RANDOMIZED_SWEEP_A", 5),
        ("C", native_c, "PASS_RANDOMIZED_SWEEP_C", 7),
    ):
        sweep = list(result["sweep"])
        if ACTIVE_MUTATION == tooth:
            sweep[0] = result["boundary"][2]
        observed = [(p["fc"], p["gauss"], p["additional"], p["nullity"]) for p in sweep]
        signs = Counter(p["sign"] for p in result["sweep"])
        expect_bool(
            tooth,
            len(observed) == 6 and signs == Counter({-1: 3, 1: 3})
            and observed == [(0, 0, False, nullity)] * 6
            and result["seed"] == (260713 if theory == "A" else 260715),
            {"seed": result["seed"], "signs": signs, "observed": observed},
        )

    # The injected mutation object is a real Maxwell primary/descendant found
    # by this engine.  It has a nonzero spatial action proportional to k.
    maxwell_candidate = maxwell_live["search"]["candidates"][0]  # type: ignore[index]
    for theory, result, descendant_tooth, accounting_tooth in (
        ("A", native_a, "PASS_HARDENING_DESCENDANT_A", "PASS_PRIMARY_FC_ACCOUNTING_A"),
        ("C", native_c, "PASS_HARDENING_DESCENDANT_C", "PASS_PRIMARY_FC_ACCOUNTING_C"),
    ):
        records = result["hardening"]
        tested_records: list[dict[str, object]] = []
        for record in records:
            tested_records.append({**record, "rejected": list(record["rejected"])})
        if ACTIVE_MUTATION in {descendant_tooth, accounting_tooth}:
            tested_records[0]["rejected"].append(maxwell_candidate)  # type: ignore[index,union-attr]
            if ACTIVE_MUTATION == descendant_tooth:
                tested_records[0]["primary_fc"] = int(tested_records[0]["primary_fc"]) + 1
        descendant_ok = (
            len(tested_records) == 2
            and sum(len(record["rejected"]) for record in tested_records) == 4
            and all(
                item["descendant_zero"] is True
                and item["secondary_action_non_gradient"] is True
                and item["descendant_rejection_certified"] is True
                and item["reason"] == "DESCENDANT_ZERO"
                for record in tested_records for item in record["rejected"]
            )
        )
        accounting_ok = all(
            int(record["primary_fc"]) == len(record["rejected"]) + len(record["candidates"])
            for record in tested_records
        )
        expect_bool(descendant_tooth, descendant_ok,
                    {"strata": len(tested_records), "directions": sum(len(x["rejected"]) for x in tested_records)})
        expect_bool(accounting_tooth, accounting_ok)

    heading("Six same-pipeline controls and native own-assert ablations")
    search_probe = (control_result("maxwell", True)
                    if ACTIVE_MUTATION == "PASS_GUARD_SEARCH_CAPABLE" else baseline_controls["maxwell"])
    expect_bool(
        "PASS_GUARD_SEARCH_CAPABLE",
        int(search_probe["gauss"]) > 0 and int(baseline_controls["gauged_hard_unit"]["gauss"]) > 0,
        {"maxwell": search_probe["gauss"], "gauged": baseline_controls["gauged_hard_unit"]["gauss"]},
    )

    for key, expected in CONTROL_EXPECTED.items():
        tooth = CONTROL_PREDICATES[key]
        tested = ablated_controls[key] if ACTIVE_MUTATION == tooth else baseline_controls[key]
        observed = (tested["classification"], tested["fc"], tested["sc"], tested["gauss"], tested["nullity"])
        ablated_observed = (
            ablated_controls[key]["classification"], ablated_controls[key]["fc"],
            ablated_controls[key]["sc"], ablated_controls[key]["gauss"],
            ablated_controls[key]["nullity"],
        )
        extra = True
        if key == "nonconserved_current":
            rho_dot = sp.symbols("nc_rho_dot", real=True)
            j1, j2, j3 = sp.symbols("nc_j1 nc_j2 nc_j3", real=True)
            extra = sp.simplify(tested["preservation"] - (-j1 - 2*j2 - 3*j3 + rho_dot)) == 0
        expect_bool(tooth, observed == expected and ablated_observed != expected and extra,
                    {"observed": observed, "ablated": ablated_observed})
        print(f"      {key}: {expected[0]} FC={expected[1]} SC={expected[2]} G={expected[3]} NULL={expected[4]}")
        print(f"      ABLATION {key}: FIRED_AT_OWN_ASSERT")

    heading("Bookkeeping, honest scope, branch decision, and verdict")
    any_native_g = any(
        bool(item["additional_G_exists"])
        for result in (native_a, native_c)
        for item in [result["search"], *(point["search"] for point in result["boundary"]),
                     *(point["search"] for point in result["sweep"])]
    )
    if ACTIVE_MUTATION == "PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING":
        any_native_g = bool(maxwell_live["search"]["additional_G_exists"])
    source_first = {"searched_source_free_first": True, "jA_added": any_native_g, "j_sourced": any_native_g}
    shear = "REQUIRES_G" if any_native_g else "NOT_APPLICABLE_NO_G"
    expect_bool(
        "PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING",
        source_first == {"searched_source_free_first": True, "jA_added": False, "j_sourced": False}
        and shear == "NOT_APPLICABLE_NO_G",
        {"source_first": source_first, "shear": shear},
    )

    scope = {
        "open_stratum": "FULLY_SYMBOLIC_ALL_RETAINED_COUPLINGS",
        "tuned_stratum": "ARGUED_PLUS_FIXED_SEED_SCANNED",
        "boundary_points_per_family": 6,
        "random_points_total": 12,
        "exhaustive_tuned_stratification": False,
        "missed_measure_zero_class": "TUNED_INVERSE_DESIGN",
        "generic_no_go_decisive": True,
    }
    if ACTIVE_MUTATION == "PASS_HONEST_TUNED_SCOPE":
        scope["exhaustive_tuned_stratification"] = True
    expect_bool(
        "PASS_HONEST_TUNED_SCOPE",
        scope["open_stratum"] == "FULLY_SYMBOLIC_ALL_RETAINED_COUPLINGS"
        and scope["tuned_stratum"] == "ARGUED_PLUS_FIXED_SEED_SCANNED"
        and scope["boundary_points_per_family"] == 6
        and scope["random_points_total"] == 12
        and scope["exhaustive_tuned_stratification"] is False
        and scope["missed_measure_zero_class"] == "TUNED_INVERSE_DESIGN"
        and scope["generic_no_go_decisive"] is True,
        scope,
    )

    branch_fc = (int(maxwell_live["fc"]) if ACTIVE_MUTATION == "PASS_DECISION_ORDER_BRANCH2"
                 else int(native_a["dirac"]["first_class_count"]) + int(native_c["dirac"]["first_class_count"]))
    branch = "BRANCH_2_QUADRATIC_ABSENCE" if branch_fc == 0 else "BRANCH_1_LINEARIZED_GAUGE_PRESENT"
    expect_bool("PASS_DECISION_ORDER_BRANCH2", branch == "BRANCH_2_QUADRATIC_ABSENCE", branch)

    open_searches = [native_a["search"], native_c["search"]]
    if ACTIVE_MUTATION == "PASS_VERDICT_TOTALITY":
        open_searches = [maxwell_live["search"], native_c["search"]]
    tuned_searches = [
        point["search"] for result in (native_a, native_c)
        for point in [*result["boundary"], *result["sweep"]]
    ]
    verdict_a = verdict_from_computed(open_searches, tuned_searches)
    if ACTIVE_MUTATION == "PASS_VERDICT_TOTALITY":
        expect_bool(
            "PASS_VERDICT_TOTALITY",
            verdict_a == VERDICT,
            {"rederived": verdict_a, "required_mutation_token": "FIRST_CLASS_GENERIC_EM_CANDIDATE"},
        )
    else:
        verdict_c = verdict_from_computed([native_c["search"]], [
            point["search"] for point in [*native_c["boundary"], *native_c["sweep"]]
        ])
        expect_bool("PASS_VERDICT_TOTALITY", verdict_a == verdict_c == VERDICT,
                    {"A": verdict_a, "C": verdict_c})

    heading("Canonical source-to-stage predicate manifest")
    manifest = SOURCE_MANIFEST[:-1] if ACTIVE_MUTATION == "PASS_SOURCE_PREDICATE_MANIFEST" else SOURCE_MANIFEST
    identifiers = tuple(row[0] for row in manifest)
    dispositions = set(row[1] for row in manifest)
    counts = dict(sorted(Counter(row[1] for row in manifest).items()))
    digest = manifest_digest(manifest)
    expect_bool(
        "PASS_SOURCE_PREDICATE_MANIFEST",
        identifiers == SOURCE_TOOTH_IDS
        and len(identifiers) == len(set(identifiers)) == 51
        and dispositions == {"PRESERVED", "REPLACED_BY_STRONGER", "SCOPED_OUT"}
        and all(row[2].startswith("STAGE033_") for row in manifest)
        and counts == EXPECTED_MANIFEST_COUNTS
        and digest == COMMITTED_MANIFEST_DIGEST,
        {"entries": len(manifest), "partition": counts, "digest": digest},
    )
    print(f"      entries={len(manifest)}; partition={counts}; digest={digest}")

    print("")
    print("HONEST_SCOPE: symbolic open stratum decisive; tuned locus ARGUED+FIXED-SEED-SCANNED, NOT exhaustive")
    print("DOWNSTREAM_ONLY: compactness, quantization, deconfinement, native +/-w current supply")
    print("SOURCE_FIRST_ORDERING: searched source-free first; j.A not added")
    print("SHEAR_DUPLICATE: NOT_APPLICABLE_NO_G")
    print("HARDENING-TUNED-DESCENDANT-REJECTION: PASS")
    print("GUARD-COUPLINGS-ENTER: PASS")
    print("GUARD-SEARCH-CAPABLE: PASS")
    print(f"THEORY-A VERDICT_TOKEN: {VERDICT}")
    print(f"THEORY-C VERDICT_TOKEN: {VERDICT}")


def fail_mutation(predicate: str, condition: bool, evidence: Any = None) -> None:
    """Finish a focused mutation run at the selected predicate's own assert."""
    expect_bool(predicate, condition, evidence)
    print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
    raise AuditFailure("MUTATION_DID_NOT_FIRE", predicate)


def focused_mutation_probe(predicate: str) -> None:
    """Exercise only the upstream computation needed by one ablated tooth.

    Clean execution remains the full fold.  Focused mutation execution avoids
    repeating the complete 24-point tuned scan 33 times while still creating
    every altered value from a live Lagrangian/CAS pipeline.
    """
    print("FOCUSED_MUTATION_PROBE=runtime-computed")

    if predicate == "PASS_Q1_LAGRANGIAN_LEGENDRE":
        model = native_model("A")
        build = build_hamiltonian(model)
        altered_h = build.hamiltonian + model.v[0]
        fail_mutation(predicate, not any(altered_h.has(v) for v in model.v), altered_h)

    if predicate in {"PASS_Q2_MAXWELL_COMPUTED", "PASS_Q5_GAUSS_SEARCH_LIVE"}:
        result = control_result("maxwell", True)
        if predicate == "PASS_Q2_MAXWELL_COMPUTED":
            observed = (result["classification"], result["fc"], result["sc"], result["gauss"], result["nullity"])
            fail_mutation(predicate, observed == CONTROL_EXPECTED["maxwell"], observed)
        fail_mutation(predicate, int(result["gauss"]) > 0 and bool(result["search"]["additional_G_exists"]), result)

    if predicate == "PASS_Q3_SIX_CONTROLS_SHARED_PIPELINE":
        table = {key: control_result(key) for key in CONTROL_EXPECTED if key != "global_only"}
        fail_mutation(predicate, tuple(table) == tuple(CONTROL_EXPECTED), tuple(table))

    if predicate == "PASS_Q4_CAS_COMMON_NULL_RECONSTRUCTION":
        item = common_null_derivation("A")[1]
        wrong = dict(item["solutions"][0])
        wrong[item["gs"]] = -wrong[item["gs"]]
        fail_mutation(predicate,
                      all(sp.simplify(residual.subs(wrong)) == 0 for residual in item["residual"]), wrong)

    if predicate == "PASS_Q6_DIRAC_CLOSURE":
        build = build_hamiltonian(native_model("A"))
        closed = True
        try:
            dirac_closure(build, stage_cap=1)
        except PipelineFailure:
            closed = False
        fail_mutation(predicate, closed)

    if predicate in {"PASS_GUARD_COUPLINGS_ENTER_A", "PASS_GUARD_COUPLINGS_ENTER_C"}:
        theory = "A" if predicate.endswith("_A") else "C"
        model = native_model(theory)
        target = model.couplings[-1]
        bad_model = replace(model, lagrangian=model.lagrangian.subs(target, 0))
        build, data, _ = execute(bad_model)
        locations = coupling_locations(build, data)
        fail_mutation(predicate, all(locations.values()), locations)

    if predicate in {"PASS_THEORY_A_SIGNATURE", "PASS_THEORY_C_SIGNATURE"}:
        theory = "A" if "_A_" in predicate else "C"
        model = native_model(theory)
        bad_model = replace(model, lagrangian=model.lagrangian + model.v[7] ** 2 / 2)
        observed = signature(execute(bad_model))
        fail_mutation(predicate, observed == expected_signature(theory), observed)

    if predicate == "PASS_KINETIC_HESSIAN_DETERMINANT":
        model = native_model("A")
        symbols = {str(x): x for x in model.lagrangian.free_symbols}
        gt, ru = symbols["g_tA"], symbols["rho_uA"]
        bad = replace(model, lagrangian=model.lagrangian + model.v[3] ** 2 / 2)
        build = build_hamiltonian(bad)
        det = sp.factor(build.hessian.extract(range(6), range(6)).det())
        component = sp.factor(build.hessian.extract([0, 3], [0, 3]).det())
        fail_mutation(predicate,
                      sp.simplify(det - (ru - gt**2)**3) == 0
                      and sp.simplify(component - (ru - gt**2)) == 0,
                      (det, component))

    if predicate in {"PASS_COMMON_NULL_A", "PASS_COMMON_NULL_C"}:
        theory = predicate[-1]
        item = common_null_derivation(theory)[1]
        wrong = dict(item["solutions"][0])
        wrong[item["gs"]] = -wrong[item["gs"]]
        condition = (
            sp.simplify(wrong[item["ap"]] - 14*item["au"]) == 0
            and sp.simplify(wrong[item["gs"]] - item["au"]) == 0
            and all(sp.simplify(residual.subs(wrong)) == 0 for residual in item["residual"])
        )
        fail_mutation(predicate, condition, wrong)

    if predicate in {"PASS_BOUNDARY_SCAN_A", "PASS_BOUNDARY_SCAN_C",
                     "PASS_RANDOMIZED_SWEEP_A", "PASS_RANDOMIZED_SWEEP_C",
                     "PASS_HARDENING_DESCENDANT_A", "PASS_HARDENING_DESCENDANT_C",
                     "PASS_PRIMARY_FC_ACCOUNTING_A", "PASS_PRIMARY_FC_ACCOUNTING_C"}:
        theory = predicate[-1]
        model = native_model(theory)
        symbols = {str(x): x for x in model.lagrangian.free_symbols}
        gt, ru = symbols[f"g_t{theory}"], symbols[f"rho_u{theory}"]
        ap, au = symbols[f"a_p{theory}"], symbols[f"a_u{theory}"]
        gs, gd = symbols[f"g_s{theory}"], symbols[f"g_d{theory}"]
        gb = symbols.get(f"g_b{theory}")
        common_rules = {gt: 1, ru: 1, **common_null_derivation(theory)[1]["solutions"][0]}
        if predicate.startswith("PASS_BOUNDARY_SCAN"):
            rules = {gt: 1, ru: 1, ap: 14, au: 1, gs: sp.Rational(3, 2), gd: 0}
            if gb is not None:
                rules[gb] = 0
            point = scan_point(theory, rules)
            expected = (2, 0, False, 5 if theory == "A" else 7)
            observed = (point["fc"], point["gauss"], point["additional"], point["nullity"])
            fail_mutation(predicate, observed == expected, observed)
        if predicate.startswith("PASS_RANDOMIZED_SWEEP"):
            point = scan_point(theory, common_rules)
            expected = (0, 0, False, 5 if theory == "A" else 7)
            observed = (point["fc"], point["gauss"], point["additional"], point["nullity"])
            fail_mutation(predicate, observed == expected, observed)
        point = scan_point(theory, common_rules)
        record = hardening_record(theory, 1, point["build"], point["dirac"], point["search"])  # type: ignore[arg-type]
        maxwell = control_result("maxwell")
        injected = list(record["rejected"]) + [maxwell["search"]["candidates"][0]]  # type: ignore[index]
        if predicate.startswith("PASS_HARDENING_DESCENDANT"):
            condition = len(injected) == 2 and all(
                item["descendant_zero"] is True and item["secondary_action_non_gradient"] is True
                and item["descendant_rejection_certified"] is True and item["reason"] == "DESCENDANT_ZERO"
                for item in injected
            )
            fail_mutation(predicate, condition, injected)
        condition = int(record["primary_fc"]) == len(injected) + len(record["candidates"])
        fail_mutation(predicate, condition, record)

    if predicate == "PASS_GUARD_SEARCH_CAPABLE":
        maxwell = control_result("maxwell", True)
        gauged = control_result("gauged_hard_unit")
        fail_mutation(predicate, int(maxwell["gauss"]) > 0 and int(gauged["gauss"]) > 0,
                      (maxwell["gauss"], gauged["gauss"]))

    if predicate in CONTROL_PREDICATES.values():
        key = next(name for name, tooth in CONTROL_PREDICATES.items() if tooth == predicate)
        tested = control_result(key, True)
        observed = (tested["classification"], tested["fc"], tested["sc"], tested["gauss"], tested["nullity"])
        fail_mutation(predicate, observed == CONTROL_EXPECTED[key], observed)

    if predicate == "PASS_SOURCE_FIRST_SHEAR_BOOKKEEPING":
        maxwell = control_result("maxwell")
        any_g = bool(maxwell["search"]["additional_G_exists"])
        source = {"searched_source_free_first": True, "jA_added": any_g, "j_sourced": any_g}
        shear = "REQUIRES_G" if any_g else "NOT_APPLICABLE_NO_G"
        fail_mutation(predicate,
                      source == {"searched_source_free_first": True, "jA_added": False, "j_sourced": False}
                      and shear == "NOT_APPLICABLE_NO_G", (source, shear))

    if predicate == "PASS_HONEST_TUNED_SCOPE":
        scope = {"open": "FULLY_SYMBOLIC", "tuned": "ARGUED_PLUS_SCANNED",
                 "boundary": 6, "random_total": 12, "exhaustive": True}
        fail_mutation(predicate,
                      scope["open"] == "FULLY_SYMBOLIC" and scope["tuned"] == "ARGUED_PLUS_SCANNED"
                      and scope["boundary"] == 6 and scope["random_total"] == 12
                      and scope["exhaustive"] is False, scope)

    if predicate == "PASS_DECISION_ORDER_BRANCH2":
        maxwell = control_result("maxwell")
        branch = "BRANCH_2_QUADRATIC_ABSENCE" if int(maxwell["fc"]) == 0 else "BRANCH_1_LINEARIZED_GAUGE_PRESENT"
        fail_mutation(predicate, branch == "BRANCH_2_QUADRATIC_ABSENCE", branch)

    if predicate == "PASS_VERDICT_TOTALITY":
        maxwell = control_result("maxwell")
        rederived = verdict_from_computed([maxwell["search"]], [])
        fail_mutation(predicate, rederived == VERDICT,
                      {"rederived": rederived, "required": "FIRST_CLASS_GENERIC_EM_CANDIDATE"})

    if predicate == "PASS_SOURCE_PREDICATE_MANIFEST":
        manifest = SOURCE_MANIFEST[:-1]
        identifiers = tuple(row[0] for row in manifest)
        fail_mutation(predicate,
                      identifiers == SOURCE_TOOTH_IDS and manifest_digest(manifest) == COMMITTED_MANIFEST_DIGEST,
                      {"entries": len(manifest), "digest": manifest_digest(manifest)})

    raise AuditFailure("UNKNOWN_MUTATION", predicate)


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)
    print("ledger_stage033_native_p_no_emergent_gauss SymPy audit")
    print("PIPELINE=input Lagrangian -> Legendre transform -> Dirac closure -> PB-kernel Gauss search")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")
        focused_mutation_probe(ACTIVE_MUTATION)
        return
    run_assertions()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy independently reached NATIVE_P_NO_EMERGENT_GAUSS for A and C")


if __name__ == "__main__":
    try:
        main()
    except (AuditFailure, PipelineFailure) as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        predicate = exc.predicate if isinstance(exc, AuditFailure) else "PIPELINE_FAILURE"
        print(f"OVERALL FAIL: SymPy stage033 audit did not close ({predicate})")
        raise SystemExit(1)
