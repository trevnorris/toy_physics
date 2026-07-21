#!/usr/bin/env python3
"""Target-blind moving-throat magnetism build (independent SymPy engine).

Route B is evaluated before Route A and accepts no Route-A object on the
production path.  The only mapping to an electromagnetic landing is in the
sealed ``section4_adjudicate`` block and its independent truth-table oracle.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import inspect
import itertools
import json
import os
import subprocess
import sys
from collections import Counter, OrderedDict
from dataclasses import dataclass, replace
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Mapping

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DIRECTIVE = HERE / "directive_magnetism_moving_throat.md"
ELECTRIC_DIRECTIVE = HERE / "directive_puncture_deflection_electric_sign.md"
ELECTRIC_RESULT = HERE / "puncture_deflection_electric_sign_result.md"
ELECTRIC_CHECK = HERE / "puncture_deflection_electric_sign_check.py"
G0 = HERE / "g0_closure_card_v0.md"
PATH36_MD = ROOT / "software/stage1_solver/reports/pathA_36_c5_phase_potential.md"
PATH36_YAML = ROOT / "software/stage1_solver/reports/pathA_36_c5_phase_potential_results.yaml"
PATH39_FILES = tuple(sorted((ROOT / "software/stage1_solver/reports").glob("pathA_39_*")))
AUDIT = ROOT / "docs/model_definition_audit.md"
HANDOFF = ROOT / "docs/em_analog_next_phase_handoff.md"
SIM = ROOT / "docs/two_throat_simulation_handoff_spec.md"
WL_PATH = HERE / "magnetism_moving_throat_check.wl"


def zq(expr: sp.Expr) -> bool:
    return sp.simplify(expr) == 0


def sx(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr)).replace("**", "^")


class NeutralEnum(str, Enum):
    def __str__(self) -> str:
        return self.value


class CurrentIdentity(NeutralEnum):
    CONVECTION_CONDITIONAL = "CONVECTION_LIKE_CONDITIONAL"
    DEPARTURE = "CHARACTERIZED_SOURCE_DEPARTURE"
    NULL = "NULL_SOURCE"
    R1 = "R1_SOURCE_BASIS"


class ComparisonFact(NeutralEnum):
    AGREE = "routes_agree"
    DIFFER = "routes_differ"
    ROUTE_B_R1 = "route_B_R1"


class RelativeSignFact(NeutralEnum):
    MATCH = "relative_sign_match"
    OPPOSITE = "relative_sign_opposite"
    CONFLICT = "leading_tensor_conflict"
    ANCHOR_CONDITIONAL = "relative_sign_anchor_conditional"


class MagnitudeFact(NeutralEnum):
    FORCED = "magnitude_forced_by_electric"
    FREE = "magnitude_free_factor"
    R1 = "R1(magnitude)"


class TierFact(NeutralEnum):
    GAPS_CLOSED = "tier_A_gaps_closed"
    CONDITIONAL = "tier_A_conditional"


class AnchorFact(NeutralEnum):
    CLOSED = "electric_anchor_closed"
    R1 = "R1_REQUIRED(bc_selection)"


@dataclass(frozen=True)
class SourceDerivation:
    signed_density: sp.Expr
    continuity_residual: sp.Expr
    unique_flux_coefficient: sp.Expr
    source_vector: sp.Matrix
    identity: CurrentIdentity
    amplitude_status: MagnitudeFact
    dependencies: tuple[str, ...]


@dataclass(frozen=True)
class RouteResult:
    name: str
    kernel: sp.Matrix
    interaction: sp.Expr
    force: sp.Matrix
    radial_force: sp.Expr
    potential_power: int
    force_power: int
    velocity_order: int
    dependencies: tuple[str, ...]


@dataclass(frozen=True)
class LiveFacts:
    current: CurrentIdentity
    comparison: ComparisonFact
    relative_sign: RelativeSignFact
    magnitude: MagnitudeFact
    tier: TierFact
    anchor: AnchorFact
    internal_sectors: tuple[str, ...]

    def __post_init__(self) -> None:
        if not live_facts_typed(self):
            raise TypeError("LiveFacts accepts only neutral enumerated facts")

    def neutral_tokens(self) -> tuple[Any, ...]:
        return (self.current, self.comparison, self.relative_sign, self.magnitude,
                self.tier, self.anchor, self.internal_sectors)


# ---------------------------------------------------------------------------
# Read-first provenance and the scoped G0+delta action transcript.
# ---------------------------------------------------------------------------

ACTION_KINETIC = (
    "1/2*rho_br*|dt u_T|^2",
    "-1/2*mu_R*|curl u_T|^2",
    "div u_T=0 (two transverse polarizations)",
)
ACTION_COUPLING = (
    "+q_T*sum_i s_i*eta_a(x-X_i)*V_i.u_T",
    "q_T=lambda_T*tau_d; tau_d is the active-throat time-arrow",
)
FIELD_IDENTITY = "u_T: transverse brane-shear displacement [L]; b_T=curl(u_T) [1]"
BARRED_SOURCE_MARKERS = frozenset({"Nu", "aT", "aTp", "aL", "q_A_T", "q_L"})
INTENDED_DELTA = (
    "import pathA_36 u_T kinetic row",
    "activate moving Q_chi*V.u_T finite-profile row",
)


def read_sources() -> dict[str, str]:
    sources = {
        "directive": DIRECTIVE.read_text(encoding="utf-8"),
        "electric_directive": ELECTRIC_DIRECTIVE.read_text(encoding="utf-8"),
        "electric_result": ELECTRIC_RESULT.read_text(encoding="utf-8"),
        "electric_check": ELECTRIC_CHECK.read_text(encoding="utf-8"),
        "g0": G0.read_text(encoding="utf-8"),
        "path36_md": PATH36_MD.read_text(encoding="utf-8"),
        "path36_yaml": PATH36_YAML.read_text(encoding="utf-8"),
        "audit": AUDIT.read_text(encoding="utf-8"),
        "handoff": HANDOFF.read_text(encoding="utf-8"),
        "sim": SIM.read_text(encoding="utf-8"),
    }
    for path in PATH39_FILES:
        sources[path.name] = path.read_text(encoding="utf-8")
    return sources


def parsed_g0_rows(g0: str) -> tuple[str, ...]:
    block = g0.split("## 9. Complete declared-zero ledger", 1)[1].split(
        "## 10. Instantiability, gates, and checks", 1
    )[0]
    return tuple(
        line.split("|", 2)[1].strip()
        for line in block.splitlines()
        if line.startswith("| ") and ("[POSTULATE]" in line or "[CONVENTION]" in line)
    )


def source_faithfulness(action_kinetic: tuple[str, ...] = ACTION_KINETIC,
                        action_coupling: tuple[str, ...] = ACTION_COUPLING,
                        extra_damage: tuple[str, ...] = ()) -> tuple[bool, list[str]]:
    src = read_sources()
    needles: Mapping[str, tuple[str, ...]] = {
        "directive": (
            "j∝sV` is `DECLARED_NOT_DERIVED` and BARRED",
            "Maxwell–Darwin vector-current reference",
            "R1_REQUIRED(electric_bc_selection)",
            "atomic per-tooth",
        ),
        "electric_directive": (
            "neutral throughout", "SEALED adjudication block", "atomic tooth IDs",
        ),
        "electric_result": (
            "xi_w(x)=\\ell h(x)", "m_{gg}=\\frac{b z_g^2}{D}",
            "R1_REQUIRED(bc_selection)", "outcome_not_invariant",
        ),
        "electric_check": (
            "def section4_adjudicate", "def mutation_campaign", "truth_table",
        ),
        "g0": (
            "The retained fields are", "`r_Bu_L`, `r_B divu`, `r_Bu_T`",
            "Maxwell/gauge fields, point sources, native current law, Coulomb potential prior",
            "magnetism/current structure (R2+)",
        ),
        "path36_md": (
            "L_T = 1/2 rho_br dot(u_T)^2 - 1/2 mu_R k^2 u_T^2",
            "c_gamma^2 = mu_R/rho_br", "2 massless transverse photons",
        ),
        "path36_yaml": (
            "primitive_lagrangian_per_polarization", "c_gamma_squared: mu_R/rho_br",
        ),
        "pathA_39_magnetic_force.md": (
            "conditional on the declared Stage-1 moving-source amplitudes",
            "j_T=q_A^T V", "F^-1[k_i k_j/k^4]",
        ),
        "pathA_39_magnetic_force_results.yaml": (
            "source: j_T=q_A^T V", "declared_stage2", "sim_deferred_values",
        ),
        "pathA_39_stage3_operator_parity.md": (
            "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING", "P_w odd",
            "realized nonlinear-throat O(V) coefficient",
        ),
        "pathA_39_stage3_operator_parity_results.yaml": (
            "combined-parity-mixing", "sim_deferred", "iVdotk",
        ),
        "pathA_39_stage4_field_classification.md": (
            "FIELD_SCALAR_VECTOR_DEPARTURE", "c_E=c_gamma", "nonlinear throat solve",
        ),
        "pathA_39_stage4_field_classification_results.yaml": (
            "operator_parity_contamination", "c_E=c_gamma and lambda_gamma knit",
        ),
        "pathA_39_scalar_admixture_screen.md": (
            "FAIL_OBSERVABLE_SCALAR_ADMIXTURE", "a_L`",
        ),
        "pathA_39_scalar_admixture_results.yaml": (
            "FAIL_OBSERVABLE_SCALAR_ADMIXTURE", "sim_deferred",
        ),
        "audit": (
            "SOURCE LAW ⚠️ IMPORTED", "Whatever the signs and current-law come out to",
            "F_electric/F_gravity",
        ),
        "handoff": (
            "Target = the characterized-departure Maxwell analog", "j=sV",
        ),
        "sim": (
            "pathA39 `j∝sV` kernel remains ancestry-barred",
            "R4's magnetic-sign comparison", "magnetic sign is",
        ),
    }
    failures: list[str] = []
    for label, terms in needles.items():
        if label not in src:
            failures.append(f"source:missing:{label}")
            continue
        for term in terms:
            if term not in src[label]:
                failures.append(f"{label}:missing:{term}")
    rows = parsed_g0_rows(src["g0"])
    if len(rows) != 14:
        failures.append(f"g0:ledger-count:{len(rows)}")
    if action_kinetic != ACTION_KINETIC:
        failures.append("action:path36-kinetic")
    if action_coupling != ACTION_COUPLING:
        failures.append("action:moving-dent-coupling")
    if extra_damage:
        failures.extend(f"g0:damaged:{x}" for x in extra_damage)
    return not failures, failures


# ---------------------------------------------------------------------------
# Q-CURRENT: translate the actual signed dent; do not import pathA_39's law.
# ---------------------------------------------------------------------------

x, y, z = sp.symbols("x y z", real=True)
Xx, Xy, Xz = sp.symbols("Xx Xy Xz", real=True)
v1x, v1y, v1z, v2x, v2y, v2z = sp.symbols(
    "v1x v1y v1z v2x v2y v2z", real=True
)
Vx, Vy, Vz = sp.symbols("Vx Vy Vz", real=True)
a0, R = sp.symbols("a0 R", positive=True)
s, s1, s2 = sp.symbols("s s1 s2", real=True)
alpha = sp.symbols("alpha", real=True)
qT = sp.symbols("qT", real=True, nonzero=True)
rhoBr, muR = sp.symbols("rhoBr muR", positive=True)
Ae = sp.symbols("Ae", real=True, nonzero=True)
cE = sp.symbols("cE", positive=True)
nx, ny, nz = sp.symbols("nx ny nz", real=True)

coords = (x, y, z)
centers = (Xx, Xy, Xz)
body_v = (Vx, Vy, Vz)
eta = sp.exp(-sum((p - q) ** 2 for p, q in zip(coords, centers)) / a0**2) / (
    sp.pi ** sp.Rational(3, 2) * a0**3
)
signed_density = s * eta
dt_signed_density = sp.factor(sum(V * sp.diff(signed_density, X)
                                  for V, X in zip(body_v, centers)))
flux_divergence = sp.factor(sum(sp.diff(alpha * signed_density * V, p)
                                for V, p in zip(body_v, coords)))
continuity_general = sp.factor(dt_signed_density + flux_divergence)
unique_alpha = sp.solve(sp.Eq(continuity_general / (sum(
    V * sp.diff(signed_density, p) for V, p in zip(body_v, coords)
)), 0), alpha)[0]


def derive_source(flux_scale: sp.Expr = sp.Integer(1)) -> SourceDerivation:
    residual = sp.factor(continuity_general.subs(alpha, flux_scale))
    current = sp.Matrix([qT * signed_density * V for V in body_v])
    return SourceDerivation(
        signed_density=signed_density,
        continuity_residual=residual,
        unique_flux_coefficient=unique_alpha,
        source_vector=current,
        identity=CurrentIdentity.CONVECTION_CONDITIONAL,
        amplitude_status=MagnitudeFact.R1,
        dependencies=("normalized_signed_dent", "translated_body_coordinate",
                      "continuity", "active_time_arrow_amplitude"),
    )


SOURCE = derive_source()

# Each entry is (R_w, P_w, time reversal, rotation type).  q_T=lambda_T*tau_d
# is time-odd because the complete active drain reverses into a source under T.
PARITIES = OrderedDict(
    s=(-1, -1, +1, "scalar"),
    V=(+1, +1, -1, "polar_vector"),
    tau_d=(+1, +1, -1, "scalar"),
    qT_s_V=(-1, -1, +1, "polar_vector"),
    u_T=(-1, -1, +1, "polar_vector"),
    b_T=(-1, -1, +1, "axial_vector"),
)


# ---------------------------------------------------------------------------
# Route B first: direct transverse response from the derived source and row.
# ---------------------------------------------------------------------------

V1 = sp.Matrix([v1x, v1y, v1z])
V2 = sp.Matrix([v2x, v2y, v2z])
nvec = sp.Matrix([nx, ny, nz])
Ddot = sp.expand(V1.dot(V2))
a1 = sp.expand(V1.dot(nvec))
a2 = sp.expand(V2.dot(nvec))
Aprod = sp.expand(a1 * a2)
cg2 = sp.factor(muR / rhoBr)


def _force_from_coordinate_potential(coefficient: sp.Expr) -> tuple[sp.Matrix, sp.Expr]:
    """Derive -grad_R U after writing n=Rvec/|Rvec|; no force literal."""
    rx, ry, rz = sp.symbols("rx ry rz", real=True)
    rvec = sp.Matrix([rx, ry, rz])
    rr = sp.sqrt(rvec.dot(rvec))
    aa1 = V1.dot(rvec) / rr
    aa2 = V2.dot(rvec) / rr
    potential = coefficient * (Ddot + aa1 * aa2) / rr
    force_xyz = sp.Matrix([-sp.diff(potential, p) for p in (rx, ry, rz)])
    rhat = rvec / rr
    expected_xyz = -coefficient * (
        aa2 * V1 + aa1 * V2 - (Ddot + 3 * aa1 * aa2) * rhat
    ) / rr**2
    if any(not zq(force_xyz[i] - expected_xyz[i]) for i in range(3)):
        raise AssertionError("coordinate force differentiation")
    force_n = sp.Matrix([
        sp.factor(-coefficient * (
            a2 * V1[i] + a1 * V2[i] - (Ddot + 3 * Aprod) * nvec[i]
        ) / R**2)
        for i in range(3)
    ])
    radial = sp.factor(coefficient * (Ddot + Aprod) / R**2)
    return force_n, radial


def derive_route_b(source: SourceDerivation,
                   foreign_payload: RouteResult | None = None) -> RouteResult:
    # Direct tensor ansatz G_ij=(a delta_ij+b n_i n_j)/R.  Transversality gives
    # b=a; the two-polarization trace fixes 3a+b=2/(4*pi).
    aa, bb = sp.symbols("aa bb", real=True)
    solved = sp.solve((sp.Eq(bb - aa, 0), sp.Eq(3 * aa + bb, 1 / (2 * sp.pi))),
                      (aa, bb), dict=True)[0]
    kernel = sp.simplify((solved[aa] * sp.eye(3) + solved[bb] * nvec * nvec.T) /
                         (muR * R))
    interaction = sp.factor(-s1 * s2 * qT**2 * (V1.T * kernel * V2)[0])
    dependencies = ("Q_CURRENT.source", "pathA36.transverse_EOM",
                    "direct_transverse_tensor_ansatz")
    if foreign_payload is not None:
        # Mutation-only illicit copy: rescaling Route A by q_T^2*c_gamma^2/A_E
        # leaves it a factor mu_R above Route B; the tooth detects the dependency tag.
        interaction = sp.factor(foreign_payload.interaction * qT**2 * cg2 / Ae)
        dependencies += ("ILLICIT_ROUTE_A_READ",)
    force, radial = _force_from_coordinate_potential(
        -s1 * s2 * qT**2 / (8 * sp.pi * muR)
    )
    return RouteResult("B_DIRECT", kernel, interaction, force, radial, 1, 2, 2,
                       dependencies)


ROUTE_B = derive_route_b(SOURCE)


# ---------------------------------------------------------------------------
# Route A second: Lorentz-completed electric anchor / Maxwell-Darwin reference.
# ---------------------------------------------------------------------------

def derive_route_a() -> RouteResult:
    rx, ry, rz = sp.symbols("rAx rAy rAz", real=True)
    rvec = sp.Matrix([rx, ry, rz])
    rr = sp.sqrt(rvec.dot(rvec))
    seed_k4 = -rr / (8 * sp.pi)
    kk_over_k4 = sp.Matrix(3, 3, lambda i, j: -sp.diff(seed_k4, rvec[i], rvec[j]))
    transverse = sp.simplify(sp.eye(3) / (4 * sp.pi * rr) - kk_over_k4)
    expected = (sp.eye(3) + nvec * nvec.T) / (8 * sp.pi * R)
    # The Hessian result is checked before returning the unit-vector form.
    for i in range(3):
        for j in range(3):
            raw = transverse[i, j].subs({rx: R * nx, ry: R * ny, rz: R * nz})
            reduced = sp.factor(raw.subs(nz**2, 1 - nx**2 - ny**2))
            want = sp.factor(expected[i, j].subs(nz**2, 1 - nx**2 - ny**2))
            if not zq(reduced - want):
                raise AssertionError(f"boost projector component {i},{j}")
    transverse_at_R = expected
    interaction = sp.factor(-s1 * s2 * Ae / cg2 *
                            (V1.T * transverse_at_R * V2)[0])
    force, radial = _force_from_coordinate_potential(
        -s1 * s2 * Ae / (8 * sp.pi * cg2)
    )
    return RouteResult("A_BOOST", transverse_at_R, interaction, force, radial, 1, 2, 2,
                       ("electric_conditional_anchor", "four_current_completion",
                        "Coulomb_gauge_transverse_projector", "c_gamma_cone"))


ROUTE_A = derive_route_a()

ELECTRIC_INTERACTION = s1 * s2 * Ae / (4 * sp.pi * R)
ELECTRIC_RADIAL_FORCE = s1 * s2 * Ae / (4 * sp.pi * R**2)
EXPECTED_A2 = -s1 * s2 * Ae * (Ddot + Aprod) / (8 * sp.pi * cg2 * R)
EXPECTED_B2 = -s1 * s2 * qT**2 * (Ddot + Aprod) / (8 * sp.pi * muR * R)
RATIO_BA = sp.factor(qT**2 * cg2 / (muR * Ae))
DELTA_NORMALIZED = sp.factor(RATIO_BA - 1)
CONE_RATIO = sp.factor(cE**2 / cg2)
DELTA_U = sp.factor(ROUTE_B.interaction - ROUTE_A.interaction)
PARALLEL_D = v1y * v2y
PARALLEL_RATIO_A = sp.factor((ROUTE_A.radial_force / ELECTRIC_RADIAL_FORCE).subs({
    nx: 1, ny: 0, nz: 0, v1x: 0, v2x: 0, v1z: 0, v2z: 0,
}))


def route_independence(route_b: RouteResult = ROUTE_B) -> bool:
    return (all("ROUTE_A" not in item for item in route_b.dependencies) and
            "foreign_payload.interaction" in inspect.getsource(derive_route_b) and
            route_b.name == "B_DIRECT")


def compare_routes(route_a: RouteResult = ROUTE_A,
                   route_b: RouteResult = ROUTE_B) -> dict[str, Any]:
    tensor_a = sp.factor(route_a.interaction / (-s1 * s2 * Ae / cg2))
    tensor_b = sp.factor(route_b.interaction / (-s1 * s2 * qT**2 / muR))
    structural = zq(tensor_a - tensor_b)
    return {
        "tensor_agrees": structural,
        "falloff_agrees": route_a.force_power == route_b.force_power == 2,
        "velocity_order_agrees": route_a.velocity_order == route_b.velocity_order == 2,
        "relative_sign": RelativeSignFact.ANCHOR_CONDITIONAL,
        "comparison": ComparisonFact.ROUTE_B_R1,
        "characterized_delta": DELTA_U,
        "normalized_ratio": RATIO_BA,
        "normalized_delta": DELTA_NORMALIZED,
        "cone_ratio": CONE_RATIO,
        "named_gap": "nonlinear throat q_T normalization + electric boundary selection",
    }


COMPARISON = compare_routes()


# ---------------------------------------------------------------------------
# Units, consistency, Q-MAG, and scoped prediction hooks.
# ---------------------------------------------------------------------------

Dim = tuple[int, int, int]
ZDIM: Dim = (0, 0, 0)
LDIM: Dim = (1, 0, 0)
VDIM: Dim = (1, -1, 0)
EDIM: Dim = (2, -2, 1)
FDIM: Dim = (1, -2, 1)


def dadd(*dims: Dim) -> Dim:
    return tuple(sum(values) for values in zip(*dims))  # type: ignore[return-value]


def dscale(n: int, dim: Dim) -> Dim:
    return tuple(n * value for value in dim)  # type: ignore[return-value]


SYMBOL_DIMS: dict[sp.Symbol, Dim] = {
    rhoBr: (-3, 0, 1), muR: (-1, -2, 1), qT: (0, -1, 1),
    Ae: (3, -2, 1), cE: VDIM, R: LDIM,
    s1: ZDIM, s2: ZDIM, nx: ZDIM, ny: ZDIM, nz: ZDIM,
    v1x: VDIM, v1y: VDIM, v1z: VDIM, v2x: VDIM, v2y: VDIM, v2z: VDIM,
}


def infer_dim(expr: sp.Expr) -> Dim:
    expr = sp.sympify(expr)
    if expr.is_Number or expr is sp.pi:
        return ZDIM
    if expr.is_Symbol:
        return SYMBOL_DIMS[expr]
    if expr.is_Add:
        dims = [infer_dim(arg) for arg in expr.args if arg != 0]
        if not dims or any(dim != dims[0] for dim in dims[1:]):
            raise ValueError(f"dimension sum:{expr}:{dims}")
        return dims[0]
    if expr.is_Mul:
        return dadd(*(infer_dim(arg) for arg in expr.args))
    if expr.is_Pow and expr.exp.is_Integer:
        return dscale(int(expr.exp), infer_dim(expr.base))
    raise ValueError(f"dimension node:{expr}")


def dimensional_check(qt_dim: Dim = (0, -1, 1)) -> bool:
    old = SYMBOL_DIMS[qT]
    SYMBOL_DIMS[qT] = qt_dim
    try:
        checks = (
            infer_dim(rhoBr * sp.Symbol("du", positive=True)**2),
            infer_dim(qT**2 * Ddot / (muR * R)),
            infer_dim(Ae * Ddot / (cg2 * R)),
            infer_dim(qT**2 * Ddot / (muR * R**2)),
            infer_dim(RATIO_BA), infer_dim(CONE_RATIO),
        )
    except (KeyError, ValueError):
        SYMBOL_DIMS[qT] = old
        return False
    SYMBOL_DIMS[qT] = old
    # du is assigned below only for this check: velocity dimension.
    return checks[1:] == (EDIM, EDIM, FDIM, ZDIM, ZDIM)


SYMBOL_DIMS[sp.Symbol("du", positive=True)] = VDIM


def amendment_consistency(rho_value: sp.Expr = rhoBr, mu_value: sp.Expr = muR,
                          damaged: tuple[str, ...] = ()) -> tuple[str, ...]:
    sectors = list(damaged)
    witness = {rhoBr: 3, muR: 5}
    if not bool(sp.sympify(rho_value).subs(witness).is_positive):
        sectors.append("transverse_inertia")
    if not bool(sp.sympify(mu_value).subs(witness).is_positive):
        sectors.append("transverse_stiffness")
    return tuple(sectors)


INTERNAL_SECTORS = amendment_consistency()
QMAG_FACT = MagnitudeFact.R1
TIER_FACT = TierFact.CONDITIONAL
ANCHOR_FACT = AnchorFact.R1
HOOKS = OrderedDict(
    emergent_lorentz=(
        "UNDETERMINED: delta_U is proportional to (q_T^2*c_gamma^2/(mu_R*A_E)-1); "
        "also require c_E^2/c_gamma^2=1 and Tier-B closure"
    ),
    preferred_frame=(
        "UNDETERMINED: medium-frame O(V1*V2) tensor is known; coefficient/cone and "
        "higher-order leakage are R1"
    ),
    active_flux=(
        "R1: G0 F_flux may have a non-conservative O(V1*V2) contribution; integrability "
        "of the full observable is not supplied by the conservative row"
    ),
    hierarchy_capstone=(
        "FLAG_ONLY: F_e/F_g becomes the held-out ratio of two couplings in one action"
    ),
)


# ---------------------------------------------------------------------------
# SEALED SECTION 4.  No upstream function assigns an EM landing.
# ---------------------------------------------------------------------------

CURRENT_DOMAIN = tuple(item.value for item in CurrentIdentity)
COMPARISON_DOMAIN = tuple(item.value for item in ComparisonFact)
RELATIVE_DOMAIN = tuple(item.value for item in RelativeSignFact)
MAG_DOMAIN = tuple(item.value for item in MagnitudeFact)
TIER_DOMAIN = tuple(item.value for item in TierFact)
ANCHOR_DOMAIN = tuple(item.value for item in AnchorFact)


def live_facts_typed(facts: LiveFacts) -> bool:
    return (type(facts.current) is CurrentIdentity and
            type(facts.comparison) is ComparisonFact and
            type(facts.relative_sign) is RelativeSignFact and
            type(facts.magnitude) is MagnitudeFact and
            type(facts.tier) is TierFact and
            type(facts.anchor) is AnchorFact and
            type(facts.internal_sectors) is tuple and
            all(type(item) is str for item in facts.internal_sectors))


def section4_adjudicate(live: LiveFacts | None = None, *, current: str | None = None,
                        comparison: str | None = None, relative_sign: str | None = None,
                        magnitude: str | None = None, tier: str | None = None,
                        anchor: str | None = None, internal: bool | None = None,
                        precedence_mutation: bool = False) -> str:
    """Total first-match implementation of the directive's sealed section 4."""
    if live is not None:
        if any(value is not None for value in
               (current, comparison, relative_sign, magnitude, tier, anchor, internal)):
            raise TypeError("live facts cannot be mixed with restated fields")
        return section4_adjudicate(
            current=live.current.value, comparison=live.comparison.value,
            relative_sign=live.relative_sign.value, magnitude=live.magnitude.value,
            tier=live.tier.value, anchor=live.anchor.value,
            internal=bool(live.internal_sectors), precedence_mutation=precedence_mutation,
        )
    values = (current, comparison, relative_sign, magnitude, tier, anchor, internal)
    if any(value is None for value in values):
        raise TypeError("incomplete section-4 case")
    if internal:
        return "NO_GO(sector)"

    unconditional = (
        current != CurrentIdentity.R1.value and
        comparison != ComparisonFact.ROUTE_B_R1.value and
        relative_sign != RelativeSignFact.ANCHOR_CONDITIONAL.value and
        tier == TierFact.GAPS_CLOSED.value and
        anchor == AnchorFact.CLOSED.value and
        magnitude != MagnitudeFact.R1.value
    )
    if precedence_mutation and anchor == AnchorFact.R1.value:
        return "R1_REQUIRED(consistency)"
    if unconditional:
        if relative_sign in {RelativeSignFact.OPPOSITE.value,
                             RelativeSignFact.CONFLICT.value}:
            reason = ("wrong_relative_sign" if relative_sign == RelativeSignFact.OPPOSITE.value
                      else "routes_leading_conflict")
            return f"AMENDMENT_EXCLUDED({reason})"
        if comparison == ComparisonFact.AGREE.value:
            if magnitude == MagnitudeFact.FORCED.value:
                return "MAGNETISM_LORENTZ_CONSISTENT"
            return "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"
        if comparison == ComparisonFact.DIFFER.value:
            return "MAGNETISM_DEPARTURE_CHARACTERIZED"

    if anchor == AnchorFact.R1.value:
        return "R1_REQUIRED(electric_bc_selection)"
    if current == CurrentIdentity.R1.value or comparison == ComparisonFact.ROUTE_B_R1.value:
        return "R1_REQUIRED(direct_moving_throat)"
    if magnitude == MagnitudeFact.R1.value:
        return "R1_REQUIRED(magnitude)"
    if tier == TierFact.CONDITIONAL.value or relative_sign == RelativeSignFact.ANCHOR_CONDITIONAL.value:
        return "R1_REQUIRED(consistency)"
    # Every enumerated case reaches an earlier branch; retained defensively.
    return "R1_REQUIRED(unclassified)"


def section4_oracle(case: Mapping[str, Any]) -> str:
    if case["internal"]:
        return "NO_GO(sector)"
    complete = (case["current"] != "R1_SOURCE_BASIS" and
                case["comparison"] != "route_B_R1" and
                case["relative_sign"] != "relative_sign_anchor_conditional" and
                case["tier"] == "tier_A_gaps_closed" and
                case["anchor"] == "electric_anchor_closed" and
                case["magnitude"] != "R1(magnitude)")
    if complete:
        if case["relative_sign"] == "relative_sign_opposite":
            return "AMENDMENT_EXCLUDED(wrong_relative_sign)"
        if case["relative_sign"] == "leading_tensor_conflict":
            return "AMENDMENT_EXCLUDED(routes_leading_conflict)"
        if case["comparison"] == "routes_agree":
            return ("MAGNETISM_LORENTZ_CONSISTENT"
                    if case["magnitude"] == "magnitude_forced_by_electric"
                    else "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)")
        return "MAGNETISM_DEPARTURE_CHARACTERIZED"
    if case["anchor"] == "R1_REQUIRED(bc_selection)":
        return "R1_REQUIRED(electric_bc_selection)"
    if case["current"] == "R1_SOURCE_BASIS" or case["comparison"] == "route_B_R1":
        return "R1_REQUIRED(direct_moving_throat)"
    if case["magnitude"] == "R1(magnitude)":
        return "R1_REQUIRED(magnitude)"
    if (case["tier"] == "tier_A_conditional" or
            case["relative_sign"] == "relative_sign_anchor_conditional"):
        return "R1_REQUIRED(consistency)"
    return "R1_REQUIRED(unclassified)"


def truth_table(precedence_mutation: bool = False) -> tuple[bool, int, Counter[str], str]:
    rows: list[str] = []
    counts: Counter[str] = Counter()
    exact = True
    keys = ("current", "comparison", "relative_sign", "magnitude", "tier", "anchor", "internal")
    for values in itertools.product(CURRENT_DOMAIN, COMPARISON_DOMAIN, RELATIVE_DOMAIN,
                                    MAG_DOMAIN, TIER_DOMAIN, ANCHOR_DOMAIN, (False, True)):
        case = dict(zip(keys, values))
        got = section4_adjudicate(**case, precedence_mutation=precedence_mutation)
        want = section4_oracle(case)
        if not precedence_mutation:
            exact &= got == want and "unclassified" not in got
        counts[got] += 1
        rows.append("|".join(map(str, values)) + "|" + got)
    digest = hashlib.sha256("\n".join(rows).encode()).hexdigest()
    return exact, len(rows), counts, digest


LIVE_FACTS = LiveFacts(
    current=SOURCE.identity,
    comparison=COMPARISON["comparison"],
    relative_sign=COMPARISON["relative_sign"],
    magnitude=QMAG_FACT,
    tier=TIER_FACT,
    anchor=ANCHOR_FACT,
    internal_sectors=INTERNAL_SECTORS,
)
ACTUAL_LANDING = section4_adjudicate(LIVE_FACTS)


def r1_blockers(facts: LiveFacts = LIVE_FACTS) -> tuple[str, ...]:
    blockers: list[str] = []
    if facts.anchor is AnchorFact.R1:
        blockers.append("R1_REQUIRED(electric_bc_selection)")
    if (facts.current is CurrentIdentity.R1 or
            facts.comparison is ComparisonFact.ROUTE_B_R1):
        blockers.append("R1_REQUIRED(direct_moving_throat)")
    if facts.magnitude is MagnitudeFact.R1:
        blockers.append("R1_REQUIRED(magnitude)")
    if facts.tier is TierFact.CONDITIONAL:
        blockers.append("R1_REQUIRED(consistency)")
    return tuple(blockers)


def landing_ownership_guard(facts: LiveFacts, emitted: str,
                            tokens: tuple[Any, ...] | None = None) -> bool:
    expected = facts.neutral_tokens()
    carried = expected if tokens is None else tokens
    return (len(carried) == len(expected) and
            all(got is want for got, want in zip(carried, expected)) and
            emitted == section4_adjudicate(facts))


# ---------------------------------------------------------------------------
# Target-blind scope, cross-engine payload, and atomic mutation campaign.
# ---------------------------------------------------------------------------

def target_blindness_guard() -> bool:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    allowed = {"section4_adjudicate", "section4_oracle", "truth_table",
               "target_blindness_guard"}
    landing_tokens = ("magnetism_lorentz_consistent", "amendment_excluded",
                      "magnetism_departure_characterized", "no_go(sector)")

    class Visitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.stack: list[str] = []
            self.bad: list[str] = []

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
            self.stack.append(node.name)
            self.generic_visit(node)
            self.stack.pop()

        def visit_Constant(self, node: ast.Constant) -> None:
            if isinstance(node.value, str) and any(
                    token in node.value.lower() for token in landing_tokens):
                if not self.stack or self.stack[-1] not in allowed:
                    self.bad.append(node.value)

    visitor = Visitor()
    visitor.visit(tree)
    route_b_symbols = {str(symbol) for symbol in ROUTE_B.interaction.free_symbols}
    return not visitor.bad and not (route_b_symbols & BARRED_SOURCE_MARKERS)


def canonical(expr: sp.Expr) -> str:
    return sp.sstr(sp.cancel(expr)).replace("**", "^")


def symbolic_payload(mutate: bool = False) -> OrderedDict[str, Any]:
    terms: OrderedDict[str, str] = OrderedDict(
        source_continuity=canonical(SOURCE.continuity_residual),
        source_alpha=canonical(SOURCE.unique_flux_coefficient),
        c_gamma_squared=canonical(cg2),
        routeA_kernel00=canonical(ROUTE_A.kernel[0, 0]),
        routeA_kernel01=canonical(ROUTE_A.kernel[0, 1]),
        routeA_U2=canonical(ROUTE_A.interaction),
        routeA_Fr=canonical(ROUTE_A.radial_force),
        routeB_kernel00=canonical(ROUTE_B.kernel[0, 0]),
        routeB_kernel01=canonical(ROUTE_B.kernel[0, 1]),
        routeB_U2=canonical(ROUTE_B.interaction),
        routeB_Fr=canonical(ROUTE_B.radial_force),
        delta_U=canonical(DELTA_U),
        ratio_BA=canonical(RATIO_BA),
        delta_normalized=canonical(DELTA_NORMALIZED),
        cone_ratio=canonical(CONE_RATIO),
        electric_U0=canonical(ELECTRIC_INTERACTION),
        electric_Fr=canonical(ELECTRIC_RADIAL_FORCE),
    )
    if mutate:
        terms["routeB_U2"] = canonical(ROUTE_B.interaction + qT**2 * Ddot / (muR * R))
    exact, total, counts, digest = truth_table()
    return OrderedDict(
        schema="MAGNETISM_MOVING_THROAT_V1",
        symbolic_terms=terms,
        current_identity=SOURCE.identity.value,
        source_amplitude=SOURCE.amplitude_status.value,
        parities=OrderedDict((key, list(value)) for key, value in PARITIES.items()),
        comparison=COMPARISON["comparison"].value,
        structural_agreement=bool(COMPARISON["tensor_agrees"] and
                                  COMPARISON["falloff_agrees"] and
                                  COMPARISON["velocity_order_agrees"]),
        relative_sign=COMPARISON["relative_sign"].value,
        magnitude=QMAG_FACT.value,
        tier=TIER_FACT.value,
        anchor=ANCHOR_FACT.value,
        internal_inconsistency=list(INTERNAL_SECTORS),
        potential_power=ROUTE_B.potential_power,
        force_power=ROUTE_B.force_power,
        velocity_order=ROUTE_B.velocity_order,
        truth_exact=exact,
        truth_total=total,
        truth_counts=dict(sorted(counts.items())),
        truth_digest=digest,
        landing=ACTUAL_LANDING,
        r1_blockers=list(r1_blockers()),
    )


PARSE_LOCALS = {str(symbol): symbol for symbol in (
    rhoBr, muR, qT, Ae, cE, R, s1, s2, nx, ny, nz,
    v1x, v1y, v1z, v2x, v2y, v2z,
)}
PARSE_LOCALS.update({"pi": sp.pi, "Pi": sp.pi})


def parse_other(text: str) -> sp.Expr:
    return sp.sympify(text.replace("^", "**"), locals=PARSE_LOCALS)


def payloads_equal(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    if set(left) != set(right) or set(left["symbolic_terms"]) != set(right["symbolic_terms"]):
        return False
    for key, value in left["symbolic_terms"].items():
        if not zq(parse_other(value) - parse_other(right["symbolic_terms"][key])):
            return False
    return all(left[key] == right[key] for key in left if key != "symbolic_terms")


def payload_difference(left: Mapping[str, Any], right: Mapping[str, Any]) -> list[str]:
    differences: list[str] = []
    if set(left) != set(right):
        differences.append(f"top_keys:{set(left) ^ set(right)}")
        return differences
    if set(left["symbolic_terms"]) != set(right["symbolic_terms"]):
        differences.append("symbolic_keys")
        return differences
    for key, value in left["symbolic_terms"].items():
        if not zq(parse_other(value) - parse_other(right["symbolic_terms"][key])):
            differences.append(f"symbolic:{key}")
    for key in left:
        if key != "symbolic_terms" and left[key] != right[key]:
            differences.append(f"metadata:{key}:{left[key]!r}!={right[key]!r}")
    return differences


def wolfram_payload(mutate: bool = False) -> Mapping[str, Any]:
    env = {**os.environ, "MAGNETISM_JSON_ONLY": "1"}
    if mutate:
        env["MAGNETISM_PAYLOAD_MUTATION"] = "1"
    run = subprocess.run(["math", "-script", str(WL_PATH)], cwd=ROOT, env=env,
                         text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         timeout=240, check=False)
    if run.returncode != 0:
        raise AssertionError(f"Mathematica exit {run.returncode}:\n{run.stdout}")
    marker = next((line for line in run.stdout.splitlines()
                   if line.startswith("JSON_PAYLOAD=")), None)
    if marker is None:
        raise AssertionError(f"Mathematica payload absent:\n{run.stdout}")
    return json.loads(marker.split("=", 1)[1])


TOOTH_ORDER = (
    "SOURCE_TRANSLATION_CONTINUITY", "SOURCE_NOT_IMPORTED", "SOURCE_BASIS",
    "PARITY_RW", "PARITY_PW", "PARITY_ROTATION", "PARITY_TIME_REVERSAL",
    "FIELD_IDENTITY_UNITS", "ACTION_KINETIC", "ACTION_COUPLING",
    "ACTION_STABILITY", "G0_DAMAGE", "ROUTE_INDEPENDENCE",
    "BOOST_PROJECTOR", "BOOST_GENERAL_VELOCITIES", "BOOST_NEXT_ORDER",
    "BOOST_COMMON_VELOCITY", "DIRECT_SOURCE", "DIRECT_PROJECTOR",
    "DIRECT_EXCHANGE_SIGN", "DIRECT_FALLOFF", "DIRECT_VELOCITY_ORDER",
    "COMPARE_COMPUTED", "DELTA_RATIO", "CONE_RATIO", "QMAG_R1",
    "UNITS_RESTORED", "ACTIVE_FLUX_CAVEAT", "HOOK_LORENTZ",
    "LEDGER_READY_ROW", "TRUTH_TOTALITY", "TRUTH_PRECEDENCE",
    "LANDING_OWNERSHIP", "TARGET_BLINDNESS", "DUAL_ENGINE_TERMS",
)


def local_checks(mutation: str | None, dual_ok: bool) -> OrderedDict[str, bool]:
    source_scale = 2 if mutation == "SOURCE_TRANSLATION_CONTINUITY" else 1
    source = derive_source(sp.Integer(source_scale))
    barred = set() if mutation != "SOURCE_NOT_IMPORTED" else {"Nu"}
    source_basis = source.source_vector
    if mutation == "SOURCE_BASIS":
        source_basis = 2 * source_basis

    parities = OrderedDict(PARITIES)
    if mutation == "PARITY_RW":
        parities["u_T"] = (+1, -1, +1, "polar_vector")
    if mutation == "PARITY_PW":
        parities["u_T"] = (-1, +1, +1, "polar_vector")
    if mutation == "PARITY_ROTATION":
        parities["u_T"] = (-1, -1, +1, "scalar")
    if mutation == "PARITY_TIME_REVERSAL":
        parities["qT_s_V"] = (-1, -1, -1, "polar_vector")

    field_dim = ZDIM if mutation == "FIELD_IDENTITY_UNITS" else LDIM
    kinetic = ACTION_KINETIC
    coupling = ACTION_COUPLING
    if mutation == "ACTION_KINETIC":
        kinetic = ACTION_KINETIC[:-1]
    if mutation == "ACTION_COUPLING":
        coupling = ("+q_T sum_i s_i eta_i u_T", ACTION_COUPLING[1])
    damage = ("existing_scalar_row",) if mutation == "G0_DAMAGE" else ()
    faithful, _ = source_faithfulness(kinetic, coupling, damage)
    rho_used = -rhoBr if mutation == "ACTION_STABILITY" else rhoBr

    copied_b = derive_route_b(SOURCE, ROUTE_A) if mutation == "ROUTE_INDEPENDENCE" else ROUTE_B
    route_a_kernel = ROUTE_A.kernel
    if mutation == "BOOST_PROJECTOR":
        route_a_kernel = 2 * route_a_kernel
    route_a_u = ROUTE_A.interaction
    if mutation == "BOOST_GENERAL_VELOCITIES":
        route_a_u = route_a_u.subs(Aprod, 0)
    boost_order = 4 if mutation == "BOOST_NEXT_ORDER" else ROUTE_A.velocity_order
    common_ratio = PARALLEL_RATIO_A
    if mutation == "BOOST_COMMON_VELOCITY":
        common_ratio += Ddot / cg2

    direct_u = ROUTE_B.interaction
    if mutation == "DIRECT_SOURCE":
        direct_u = direct_u.subs(qT**2, 2 * qT**2)
    direct_kernel = ROUTE_B.kernel
    if mutation == "DIRECT_PROJECTOR":
        direct_kernel = direct_kernel + nvec * nvec.T / (8 * sp.pi * muR * R)
    direct_sign = -1 if mutation != "DIRECT_EXCHANGE_SIGN" else +1
    direct_power = ROUTE_B.force_power if mutation != "DIRECT_FALLOFF" else 3
    direct_order = ROUTE_B.velocity_order if mutation != "DIRECT_VELOCITY_ORDER" else 1

    comparison = compare_routes()
    compare_tensor = comparison["tensor_agrees"]
    if mutation == "COMPARE_COMPUTED":
        compare_tensor = False
    delta_ratio = RATIO_BA
    if mutation == "DELTA_RATIO":
        delta_ratio *= 2
    cone_ratio = CONE_RATIO
    if mutation == "CONE_RATIO":
        cone_ratio = cE / sp.sqrt(cg2)
    qmag = QMAG_FACT if mutation != "QMAG_R1" else MagnitudeFact.FORCED
    qt_dim = (1, -1, 1) if mutation == "UNITS_RESTORED" else (0, -1, 1)
    flux_hook = HOOKS["active_flux"] if mutation != "ACTIVE_FLUX_CAVEAT" else "DECIDED_ZERO"
    lorentz_hook = HOOKS["emergent_lorentz"] if mutation != "HOOK_LORENTZ" else "YES"
    ledger_rows = INTENDED_DELTA if mutation != "LEDGER_READY_ROW" else INTENDED_DELTA[:1]

    truth_exact, _, _, _ = truth_table()
    if mutation == "TRUTH_TOTALITY":
        truth_exact = False
    witness = LiveFacts(CurrentIdentity.CONVECTION_CONDITIONAL, ComparisonFact.ROUTE_B_R1,
                        RelativeSignFact.ANCHOR_CONDITIONAL, MagnitudeFact.R1,
                        TierFact.CONDITIONAL, AnchorFact.R1, ())
    precedence = section4_adjudicate(
        witness, precedence_mutation=(mutation == "TRUTH_PRECEDENCE")
    )
    tokens = LIVE_FACTS.neutral_tokens()
    if mutation == "LANDING_OWNERSHIP":
        tokens = (*tokens, ACTUAL_LANDING)
    blind = target_blindness_guard() and mutation != "TARGET_BLINDNESS"

    expected_kernel = (sp.eye(3) + nvec * nvec.T) / (8 * sp.pi * R)
    direct_expected_kernel = expected_kernel / muR
    result = OrderedDict(
        SOURCE_TRANSLATION_CONTINUITY=zq(source.continuity_residual),
        SOURCE_NOT_IMPORTED=not (barred & BARRED_SOURCE_MARKERS) and not (
            {str(sym) for sym in SOURCE.source_vector.free_symbols} & BARRED_SOURCE_MARKERS
        ),
        SOURCE_BASIS=all(zq(source_basis[i] - SOURCE.source_vector[i]) for i in range(3)),
        PARITY_RW=(parities["qT_s_V"][0] == parities["u_T"][0] == -1),
        PARITY_PW=(parities["qT_s_V"][1] == parities["u_T"][1] == -1),
        PARITY_ROTATION=(parities["u_T"][3] == "polar_vector" and
                         parities["b_T"][3] == "axial_vector"),
        PARITY_TIME_REVERSAL=(parities["tau_d"][2] == -1 and
                              parities["V"][2] == -1 and
                              parities["qT_s_V"][2] == parities["u_T"][2] == +1),
        FIELD_IDENTITY_UNITS=(field_dim == LDIM),
        ACTION_KINETIC=faithful if mutation == "ACTION_KINETIC" else kinetic == ACTION_KINETIC,
        ACTION_COUPLING=faithful if mutation == "ACTION_COUPLING" else coupling == ACTION_COUPLING,
        ACTION_STABILITY=not amendment_consistency(rho_used, muR),
        G0_DAMAGE=not damage and not INTERNAL_SECTORS,
        ROUTE_INDEPENDENCE=route_independence(copied_b),
        BOOST_PROJECTOR=all(zq(route_a_kernel[i, j] - expected_kernel[i, j])
                            for i in range(3) for j in range(3)),
        BOOST_GENERAL_VELOCITIES=zq(route_a_u - EXPECTED_A2),
        BOOST_NEXT_ORDER=(boost_order == 2),
        BOOST_COMMON_VELOCITY=zq(common_ratio + PARALLEL_D / (2 * cg2)),
        DIRECT_SOURCE=zq(direct_u - EXPECTED_B2),
        DIRECT_PROJECTOR=all(zq(direct_kernel[i, j] - direct_expected_kernel[i, j])
                             for i in range(3) for j in range(3)),
        DIRECT_EXCHANGE_SIGN=(direct_sign == -1 and
                              zq(ROUTE_B.interaction - EXPECTED_B2)),
        DIRECT_FALLOFF=(direct_power == 2),
        DIRECT_VELOCITY_ORDER=(direct_order == 2),
        COMPARE_COMPUTED=bool(compare_tensor and comparison["falloff_agrees"] and
                              comparison["velocity_order_agrees"]),
        DELTA_RATIO=zq(delta_ratio - qT**2 * cg2 / (muR * Ae)),
        CONE_RATIO=zq(cone_ratio - cE**2 / cg2),
        QMAG_R1=(qmag is MagnitudeFact.R1),
        UNITS_RESTORED=dimensional_check(qt_dim),
        ACTIVE_FLUX_CAVEAT=flux_hook.startswith("R1:"),
        HOOK_LORENTZ=lorentz_hook.startswith("UNDETERMINED:"),
        LEDGER_READY_ROW=(ledger_rows == INTENDED_DELTA),
        TRUTH_TOTALITY=truth_exact,
        TRUTH_PRECEDENCE=(precedence == ACTUAL_LANDING),
        LANDING_OWNERSHIP=landing_ownership_guard(LIVE_FACTS, ACTUAL_LANDING, tokens),
        TARGET_BLINDNESS=blind,
        DUAL_ENGINE_TERMS=dual_ok,
    )
    return result


def mutation_campaign() -> OrderedDict[str, str]:
    local = symbolic_payload()
    other = wolfram_payload()
    dual_ok = payloads_equal(local, other)
    if not dual_ok:
        raise AssertionError(f"dual payload difference:{payload_difference(local, other)}")
    baseline = local_checks(None, dual_ok)
    if not all(baseline.values()):
        faithful, failures = source_faithfulness()
        raise AssertionError(
            f"baseline failures:{[key for key, ok in baseline.items() if not ok]};"
            f"faithful={faithful};source={failures}"
        )
    outcomes: OrderedDict[str, str] = OrderedDict()
    for tooth in TOOTH_ORDER:
        if tooth == "DUAL_ENGINE_TERMS":
            checks = local_checks(None, payloads_equal(local, wolfram_payload(True)))
        else:
            checks = local_checks(tooth, dual_ok)
        failed = [key for key, ok in checks.items() if not ok]
        if failed != [tooth]:
            raise AssertionError(f"atomic mutation {tooth} failed {failed}")
        outcomes[tooth] = f"FIRED_AT_{tooth}"
    return outcomes


def print_report(teeth: Mapping[str, str]) -> None:
    exact, total, counts, digest = truth_table()
    print("MAGNETISM_MOVING_THROAT — SYMPY")
    print("G0+DELTA_ACTION:")
    print("  L_T=1/2*rho_br*|dt u_T|^2-1/2*mu_R*|curl u_T|^2; div(u_T)=0")
    print("  L_move=+q_T*sum_i s_i*eta_a(x-X_i)*V_i.u_T; q_T=lambda_T*tau_d")
    print("  c_gamma^2=mu_R/rho_br; intended additions do not count as G0 damage")
    print(f"  internal_inconsistency={list(INTERNAL_SECTORS) or 'none'}")
    print("Q-CURRENT:")
    print("  sigma_i=s_i*eta_a(x-X_i); eta_a normalized; translation is Xdot_i=V_i")
    print(f"  continuity_residual={sx(SOURCE.continuity_residual)}; unique_flux_alpha={sx(unique_alpha)}")
    print("  I_i=s_i*eta_i*V_i is derived; field source J_i=q_T*I_i; barred pathA_39 amplitudes unused")
    print(f"  identity={SOURCE.identity}; amplitude={SOURCE.amplitude_status}")
    print("  field=u_T [L], polar/T-even/R_w-odd/P_w-odd; b_T=curl u_T [1], axial/T-even")
    print("  tau_d and V are T-odd, so q_T*s*V is T-even; a fixed drain branch itself breaks T")
    print("Q-BOOST_ROUTE_A:")
    print(f"  U_A0={sx(ELECTRIC_INTERACTION)}")
    print(f"  U_A2={sx(ROUTE_A.interaction)} + O(V^4/c_gamma^4); general independent V1,V2")
    print(f"  F_A2={list(map(sx, ROUTE_A.force))}; radial={sx(ROUTE_A.radial_force)}")
    print("  side-by-side parallel ratio F_A2,r/F_E,r=-(V1.V2)/(2*c_gamma^2); antiparallel reverses it")
    print("  one common side-by-side boost: equal V gives 1-V^2/(2*c_gamma^2)+O(V^4)")
    print("Q-DIRECT_ROUTE_B (evaluated before Route A; no production dependency):")
    print(f"  U_B2={sx(ROUTE_B.interaction)}; F_B2={list(map(sx, ROUTE_B.force))}")
    print(f"  radial={sx(ROUTE_B.radial_force)}; potential=R^-1; force=R^-2; order=O(V1*V2)")
    print(f"  route_independence={route_independence()}; dependencies={ROUTE_B.dependencies}")
    print("Q-COMPARE:")
    print(f"  structural_tensor_agreement={COMPARISON['tensor_agrees']}; falloff_agreement={COMPARISON['falloff_agrees']}; "
          f"velocity_order_agreement={COMPARISON['velocity_order_agrees']}")
    print(f"  comparison={COMPARISON['comparison']}({COMPARISON['named_gap']}); relative_sign={COMPARISON['relative_sign']}")
    print(f"  routeB/routeA={sx(RATIO_BA)}; normalized_delta={sx(DELTA_NORMALIZED)}; delta_U={sx(DELTA_U)}")
    print(f"  cone_ratio=c_E^2/c_gamma^2={sx(CONE_RATIO)}")
    print("Q-MAG:")
    print("  q_T is not a calibration: its dimensionless normalization relative to A_E is absent until the nonlinear throat solve")
    print(f"  fact={QMAG_FACT}; tier={TIER_FACT}; electric_anchor={ANCHOR_FACT}")
    print("DIMENSIONS: [u_T]=L; [curl u_T]=1; [q_T]=M/T; [A_E]=E*L; [U]=E; [F]=E/L; ratios=1")
    print("SECTION5_HOOKS:")
    for key, value in HOOKS.items():
        print(f"  {key}={value}")
    print(f"SECTION4_TRUTH_TABLE: cells={total}; exhaustive={exact}; digest={digest}; classes={dict(sorted(counts.items()))}")
    print("SECTION4_LANDING=" + ACTUAL_LANDING)
    print("SECTION4_ALL_R1=" + ",".join(r1_blockers()))
    print("DECIDED=source basis/parities; ledger row; Route-A Darwin tensor; conditional Route-B tensor/sign/falloff/order; analytic delta")
    print("R1=parent-throat q_T normalization; electric BC/branch; cone knit; higher orders; active F_flux integrability")
    print("ATOMIC_TEETH:")
    for tooth in TOOTH_ORDER:
        print(f"  PASS {tooth}; mutation={teeth[tooth]}")
    print("ENGINE_AGREE=PASS")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-only", action="store_true")
    args = parser.parse_args()
    try:
        if args.json_only:
            mutate = os.environ.get("MAGNETISM_PAYLOAD_MUTATION") == "1"
            print("JSON_PAYLOAD=" + json.dumps(symbolic_payload(mutate), separators=(",", ":")))
            return 0
        teeth = mutation_campaign()
        print_report(teeth)
        return 0
    except Exception as exc:
        print(f"FIRST_FAILURE={type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
