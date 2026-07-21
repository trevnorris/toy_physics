#!/usr/bin/env python3
"""Target-blind puncture-deflection electric-sign build (SymPy engine).

All objects through Q-MAG use neutral coefficient/sign labels.  The only
translation to the electric target is in ``section4_adjudicate`` and its
truth-table oracle.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
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
DIRECTIVE = HERE / "directive_puncture_deflection_electric_sign.md"
G0 = HERE / "g0_closure_card_v0.md"
PRIOR_RESULT = HERE / "leftover_scalar_electric_sign_result.md"
PRIOR_CHECK = HERE / "leftover_scalar_electric_sign_check.py"
PRIOR_VERIFY = HERE / "leftover_scalar_electric_sign_VERIFICATION.md"
CONCEPT = ROOT / "docs/conceptual_foundation.md"
PATH42 = ROOT / "software/stage1_solver/reports/pathA_42_charge_coupled_scalar.md"
AUDIT = ROOT / "docs/model_definition_audit.md"
HANDOFF = ROOT / "docs/em_analog_next_phase_handoff.md"
SIM = ROOT / "docs/two_throat_simulation_handoff_spec.md"
U2_REPORT = HERE / "reports/u2_boundary_adjudication_artifacts/stage_1_production_v12/production_summary.md"
WL_PATH = HERE / "puncture_deflection_electric_sign_check.wl"


def zq(expr: sp.Expr) -> bool:
    return sp.simplify(expr) == 0


def sx(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr)).replace("**", "^")


class NeutralEnum(str, Enum):
    """Typed neutral fact whose value is safe for reports/JSON."""

    def __str__(self) -> str:
        return self.value


class BCClass(NeutralEnum):
    DIRICHLET_VALUE = "DIRICHLET_VALUE"
    FIXED_MONOPOLE = "FIXED_MONOPOLE"
    FIXED_SOURCE = "FIXED_SOURCE"
    MIXED = "MIXED"
    UNDETERMINED_ANALYTICALLY = "UNDETERMINED_ANALYTICALLY"


class OutcomeFact(NeutralEnum):
    POSITIVE_R2 = "POSITIVE_R2"
    NEGATIVE_R2 = "NEGATIVE_R2"
    NULL = "NULL"
    POSITIVE_WRONG_RANGE = "POSITIVE_WRONG_RANGE"
    NEGATIVE_WRONG_RANGE = "NEGATIVE_WRONG_RANGE"
    OUTCOME_NOT_INVARIANT = "outcome_not_invariant"


class VariantFact(NeutralEnum):
    REPLACE = "replace"
    ADD = "add"
    BOTH = "both"
    UNRESOLVED = "variant_unresolved"


class MagnitudeFact(NeutralEnum):
    NO_FREE_FACTOR = "magnitude_no_free_factor"
    FREE_FACTOR = "magnitude_free_factor"


class TierFact(NeutralEnum):
    GAPS_CLOSED = "tier_A_gaps_closed"
    CONDITIONAL = "tier_A_conditional"


CLOSURE_GAP_DOMAIN = frozenset({
    "missing parent-throat/boundary closure", "lost_value", "lost_flux", "lost_mix",
})
AMENDMENT_SECTOR_DOMAIN = frozenset({
    "replace_ledger", "add_ledger", "S_hold", "zero_ledger", "scalar_hessian",
})


@dataclass(frozen=True)
class BCStatus:
    kind: BCClass
    missing_closure: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if (type(self.kind) is not BCClass or
                any(type(x) is not str or x not in CLOSURE_GAP_DOMAIN
                    for x in self.missing_closure)):
            raise ValueError("BCStatus must contain only enumerated neutral facts")

    def __str__(self) -> str:
        if self.kind is BCClass.UNDETERMINED_ANALYTICALLY:
            detail = "+".join(self.missing_closure) if self.missing_closure else "classifier_evidence"
            return f"{self.kind.value}({detail})"
        return self.kind.value


@dataclass(frozen=True)
class InconsistencyStatus:
    sectors: tuple[str, ...]

    def __post_init__(self) -> None:
        if any(type(x) is not str or x not in AMENDMENT_SECTOR_DOMAIN for x in self.sectors):
            raise ValueError("InconsistencyStatus must contain only enumerated amendment sectors")

    def __str__(self) -> str:
        return "none" if not self.sectors else ",".join(self.sectors)


# ---------------------------------------------------------------------------
# Source-faithful field/action/amendment transcription.
# ---------------------------------------------------------------------------

ACTION_TERMS = (
    "1/2*A_eff*(dt u_L)^2",
    "1/2*M_h*(dt h)^2",
    "-1/2*B_eff*|grad u_L|^2",
    "-1/2*K_h*|grad h|^2",
    "-C_hu*grad u_L.grad h",
)

ZERO_LEDGER = (
    "bulk `r_BH`, `r_B²H²`, `Hρ`, `Hδρ`, `H∂tθ`, and `H∇θ` couplings outside `Ω_mouth`",
    "dynamical modulation `δJ_m/δr_B`, `δJ_m/δh`, `δJ_m/δP`, neighbor-source response",
    "`r_Bu_L`, `r_B divu`, `r_Bu_T`, `Hu_T`, `u_Lu_T`, and two-gradient scalar–transverse mixing",
    "cross-kinetic `∂tu_L∂th`, one-time-derivative/Berry `u_L∂th` and `h∂tu_L`",
    "`u_L²`, reduced `h²`, `(∇²u_L)²`, `(∇²h)²`, `h³,h⁴`, and `(∇u_L)³`",
    "independent primitive `B(divu)²` in addition to `B_eff`",
    "independent `θ_B`, its amplitude-weighted kinetic/gradient terms, Josephson `θ−θ_B`, and brane-phase drain",
    "wall bending, anchoring, collar two-surface tension, surface number storage, and surface dissipation",
    "dynamic `δQ_d/δr_B`, `Γ(h_inc)`, `Γ(P_inc)`, return-location responses to `h,P`, return-kernel responses to `h,P`, and drain-rate orientation dependence",
    "direct drain sources in the `h` and `u_L` equations, and direct `h` contribution to `e_c`",
    "field-dependent geon derivatives `δE_g/δ{r_B,h,θ,H,u_L,ρ}`",
    "bulk viscosity, tangential drag, E2 no-slip traction, E3 permeability resistance, E3 phase jump",
    "all E4 multipliers, E5 Rayleigh kernels, E1 reactions, and mixture terms",
)


@dataclass(frozen=True)
class Transcript:
    action_terms: tuple[str, ...]
    parent_map: tuple[str, ...]
    field_identity: str
    held_datum: str
    mouth_terms: tuple[str, ...]
    shold_scope: str
    zero_ledger: tuple[str, ...]


TRANSCRIPT = Transcript(
    action_terms=ACTION_TERMS,
    parent_map=(
        "f0(w)=1/[ell*cosh(w/ell)^2]",
        "N0=integral dw 2*f0^2=8/(3*ell)",
        "H=f0*h+H_perp",
        "h=P0 H=N0^-1 integral dw 2*f0*H",
    ),
    field_identity="xi_w=ell*h (w-embedding displacement, dimension L)",
    held_datum="h_A=xi_w/ell=P0 H (dimension 1)",
    mouth_terms=("1/2*K_m*H(x,0)^2", "-J_m*Q_chi*H(x,0)",
                 "1/2*k_m*h^2", "-g_chih*s*h"),
    shold_scope="S_hold constrains r_B-1/2 only",
    zero_ledger=ZERO_LEDGER,
)


def read_sources() -> dict[str, str]:
    return {
        "directive": DIRECTIVE.read_text(encoding="utf-8"),
        "g0": G0.read_text(encoding="utf-8"),
        "prior_result": PRIOR_RESULT.read_text(encoding="utf-8"),
        "prior_check": PRIOR_CHECK.read_text(encoding="utf-8"),
        "prior_verify": PRIOR_VERIFY.read_text(encoding="utf-8"),
        "concept": CONCEPT.read_text(encoding="utf-8"),
        "path42": PATH42.read_text(encoding="utf-8"),
        "audit": AUDIT.read_text(encoding="utf-8"),
        "handoff": HANDOFF.read_text(encoding="utf-8"),
        "sim": SIM.read_text(encoding="utf-8"),
        "u2": U2_REPORT.read_text(encoding="utf-8"),
    }


def parsed_g0_ledger(g0: str) -> tuple[str, ...]:
    block = g0.split("## 9. Complete declared-zero ledger", 1)[1].split(
        "## 10. Instantiability, gates, and checks", 1
    )[0]
    rows = []
    for line in block.splitlines():
        if line.startswith("| ") and "[POSTULATE]" in line:
            rows.append(line.split("|", 2)[1].strip())
    return tuple(rows)


def source_faithfulness(transcript: Transcript = TRANSCRIPT,
                        hessian: sp.Matrix | None = None) -> tuple[bool, list[str]]:
    src = read_sources()
    needles = {
        "g0": (
            "parent localized scalar `H` `(-1,0,0)`", "reduced scalar `h` `(0,0,0)`",
            "longitudinal brane displacement `u_L` `(1,0,0)`", "H=f₀h+H_⊥",
            "h=P₀H:=N₀⁻¹∫dw 2f₀H", "N₀=∫dw 2f₀²=8/(3ℓ)",
            "−½B_eff|∇u_L|² −½K_h|∇h|²", "−C_hu ∇u_L·∇h",
            "½K_m H(x,0)² − J_m Q_χ[r_Σ,s_i] H(x,0)",
            "η_i(k_mh−g_χh s_i)", "S_hold=−∫dt ∫_{Γ_Σ} d³A λ_Σ(r_B−1/2)",
            "B_eff*K_h − C_hu² = 2*1 − (1/2)² = 7/4 > 0",
        ),
        "prior_result": (
            "fixed-potential \\(E_0-gh-Qu\\)", "stored energy with only the held-\\(h\\) source work removed",
            "\\(E_0-Ju-gh\\)", "det m=\\frac{z_g^2}{D}>0",
        ),
        "prior_verify": ("RESULT_CONFIRMED", "m_uu=4/7", "m_ug=−2/7", "m_gg=8/7"),
        "concept": ("electric field's flex INTO `w`", "ke²/a = m_ec²", "extrinsic curvature/embedding"),
        "path42": ("same embedding-direction family, distinct ledger object", "MIXED_SCALAR_EP_RISK"),
        "audit": ("define each piece the model's own way", "Whatever the signs and current-law come out to"),
        "handoff": ("144/144 cells UNRESOLVED", "PUNCTURE-DEFLECTION mechanism"),
        "sim": ("does not infer that a surviving imposed `𝔅_b` is nature's unique choice",),
        "u2": ("{'UNRESOLVED': 144}",),
        "directive": ("UNDETERMINED_ANALYTICALLY(missing parent-throat/boundary closure)",
                      "R1_REQUIRED(bc_selection)", "outcome_not_invariant"),
    }
    failures: list[str] = []
    for label, terms in needles.items():
        for term in terms:
            if term not in src[label]:
                failures.append(f"{label}:missing:{term}")
    if transcript != TRANSCRIPT:
        failures.append("internal:transcript")
    ledger = parsed_g0_ledger(src["g0"])
    if ledger != transcript.zero_ledger or len(ledger) != 13:
        failures.append("g0:zero-ledger")
    used = ACTION_HESSIAN if hessian is None else hessian
    want = sp.Matrix([[b, c], [c, k]])
    if used.shape != (2, 2) or any(not zq(used[i, j] - want[i, j])
                                   for i in range(2) for j in range(2)):
        failures.append("internal:action-hessian")
    return not failures, failures


# ---------------------------------------------------------------------------
# Q-FIELD and the coupled inverse, derived before any force labels.
# ---------------------------------------------------------------------------

b, c, k = sp.symbols("b c k", positive=True)
km, zeta, reb = sp.symbols("km zeta reb", positive=True)
ell = sp.symbols("ell", positive=True)
sgg = sp.symbols("sgg", positive=True)
phi, q, j, g = sp.symbols("phi q j g", positive=True)
mix = sp.symbols("mix", real=True)
s1, s2 = sp.symbols("s1 s2", real=True)
R = sp.symbols("R", positive=True)
gu, ghgrad = sp.symbols("gu ghgrad", real=True)

STATIC_DENSITY = b * gu**2 / 2 + c * gu * ghgrad + k * ghgrad**2 / 2
ACTION_HESSIAN = sp.hessian(STATIC_DENSITY, (gu, ghgrad))
D = sp.factor(ACTION_HESSIAN.det())
kappa = sp.factor(D / b)
zg = sp.factor(1 - km * ((1 - zeta) / km))
zb = sp.factor(1 - km * reb)
Z = sp.Matrix([[1, 0], [-(c / b) * zb, zg]])
m = sp.simplify(Z.T * sp.diag(1 / b, 1 / kappa) * Z)
muu, mug, mgg = map(sp.factor, (m[0, 0], m[0, 1], m[1, 1]))
mdet = sp.factor(m.det())

# At w=0, f0(0)=1/ell.  The parent source -J_m Q_chi H therefore reduces
# to -(J_m/ell) Q_chi h=-g_chih Q_chi h.  xi_w=ell*h is a physical length;
# u_L is not renamed or conflated with this normal embedding displacement.
f0_at_mouth = 1 / ell
xi_over_h = ell
parent_to_reduced_source = sp.factor(ell * f0_at_mouth)


# ---------------------------------------------------------------------------
# Q-AMEND variants.  These are scoped candidate actions, not class claims.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Amendment:
    name: str
    source_row: str
    new_rows: tuple[str, ...]
    shold_scope: str
    untouched_zero_rows: tuple[str, ...]
    declared_change: str


REPLACE = Amendment(
    "REPLACE", "core_holder_retypes_existing_h_source_BC", (),
    "r_B-1/2 only", ZERO_LEDGER,
    "existing h fixed-source mouth BC -> candidate h holding BC",
)
ADD = Amendment(
    "ADD", "existing_h_source_BC_unchanged",
    ("core_embedding_h_holding_row: coefficient beta_core [E], h dimension 1, R_w odd",),
    "r_B-1/2 only", ZERO_LEDGER,
    "one new active core holding row; all G0 rows otherwise unchanged",
)


def amendment_consistency(rep: Amendment = REPLACE, add: Amendment = ADD,
                          d_value: sp.Expr = D) -> tuple[tuple[str, ...], str]:
    sectors: list[str] = []
    if rep.new_rows or rep.source_row != "core_holder_retypes_existing_h_source_BC":
        sectors.append("replace_ledger")
    if len(add.new_rows) != 1 or add.source_row != "existing_h_source_BC_unchanged":
        sectors.append("add_ledger")
    if rep.shold_scope != "r_B-1/2 only" or add.shold_scope != "r_B-1/2 only":
        sectors.append("S_hold")
    if rep.untouched_zero_rows != ZERO_LEDGER or add.untouched_zero_rows != ZERO_LEDGER:
        sectors.append("zero_ledger")
    if not bool(sp.simplify(d_value.subs({b: 2, c: sp.Rational(1, 2), k: 1})).is_positive):
        sectors.append("scalar_hessian")
    return tuple(sectors), ("none" if not sectors else ",".join(sectors))


INCONSISTENCIES, _INCONSISTENCY_TEXT = amendment_consistency()
INCONSISTENCY_FACT = InconsistencyStatus(INCONSISTENCIES)


# ---------------------------------------------------------------------------
# Q-BC: one evidence-driven, class-typed classifier.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BCEvidence:
    name: str
    essential_value: bool
    fixed_conormal: bool
    mixed_relation: bool
    admissible_variation: bool
    signed_stationary: bool
    positive_holding_curvature: bool
    topology_invariant: bool
    mouth_barrier_computed: bool
    missing_closure: tuple[str, ...]


@dataclass(frozen=True)
class BCResult:
    name: str
    status: BCStatus
    stationarity: bool
    holding_curvature: bool
    topological_barrier: str
    reason: str


def classify_bc(ev: BCEvidence) -> BCResult:
    # Class strings originate only here, after the evidence is evaluated.
    if ev.essential_value and not ev.admissible_variation:
        status = BCStatus(BCClass.DIRICHLET_VALUE)
        reason = "essential mouth value excludes delta h_A"
    elif ev.fixed_conormal and ev.admissible_variation:
        status = BCStatus(BCClass.FIXED_MONOPOLE)
        reason = "computed conormal integral is held under free value variation"
    elif ev.mixed_relation and ev.admissible_variation:
        status = BCStatus(BCClass.MIXED)
        reason = "nonzero value/conormal relation is supplied"
    elif (not ev.missing_closure and ev.admissible_variation and
          not ev.signed_stationary and not ev.mouth_barrier_computed):
        status = BCStatus(BCClass.FIXED_SOURCE)
        reason = "free mouth relaxes under the derived natural/source condition"
    else:
        status = BCStatus(BCClass.UNDETERMINED_ANALYTICALLY, ev.missing_closure)
        reason = "far exterior stationarity does not compute the nonlinear core holder"
    barrier = "COMPUTED" if ev.mouth_barrier_computed else "NOT_COMPUTED"
    return BCResult(ev.name, status, ev.signed_stationary,
                    ev.positive_holding_curvature, barrier, reason)


# Exterior radial energy gives h=A/r and positive capacity curvature.  Its
# boundary term leaves h_A admissible.  G0 removes the nonlinear core action,
# so fixed topology s alone does not compute the map s -> h_A or its barrier.
ACTUAL_EVIDENCE = BCEvidence(
    "restored_core_candidate", False, False, False, True,
    False, True, True, False,
    ("missing parent-throat/boundary closure",),
)
FREE_CONTROL = BCEvidence(
    "free_mouth_relaxation", False, False, False, True,
    False, True, False, False, (),
)
VALUE_CONTROL = BCEvidence(
    "imposed_value", True, False, False, False,
    True, True, True, True, (),
)
MONOPOLE_CONTROL = BCEvidence(
    "imposed_monopole", False, True, False, True,
    False, True, False, True, (),
)
MIXED_CONTROL = BCEvidence(
    "imposed_mixed", False, False, True, True,
    False, True, False, True, (),
)

QBC_ACTUAL = classify_bc(ACTUAL_EVIDENCE)
QBC_CONTROLS = tuple(classify_bc(row) for row in
                     (FREE_CONTROL, VALUE_CONTROL, MONOPOLE_CONTROL, MIXED_CONTROL))
ADMISSIBLE_CLASSES = ("DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED")
VARIANT_REALIZATION = VariantFact.UNRESOLVED


# ---------------------------------------------------------------------------
# Q-FORCE: re-derive V/M/J/MIXED for h=xi_w/ell, not for u_L.
# ---------------------------------------------------------------------------

eps = sp.symbols("eps")


def fixed_value_coefficient(member: str) -> sp.Expr:
    response = sp.Matrix([[sgg, eps * mgg], [eps * mgg, sgg]])
    values = sp.Matrix([s1 * phi, s2 * phi])
    sources = sp.simplify(response.inv() * values)
    sigma = {"conjugate": -1, "bare": 1}[member]
    functional = sp.Rational(1, 2) * sigma * (values.T * sources)[0]
    pair = sp.expand(sp.series(functional, eps, 0, 2).removeO()).coeff(eps)
    return sp.factor(pair / (s1 * s2))


def same_field_coefficient(ensemble: str, member: str = "conjugate") -> sp.Expr:
    if ensemble == "M":
        total_source = q + g
        stored = mgg * total_source**2
        if member == "bare":
            return sp.factor(stored)
        return sp.factor(stored - 2 * g * mgg * total_source)
    if ensemble == "J":
        total_source = j + g
        stored = mgg * total_source**2
        if member == "bare":
            return sp.factor(stored)
        return sp.factor(stored - 2 * total_source * mgg * total_source)
    if ensemble == "MIXED":
        total_source = q + g
        stored = mgg * total_source**2
        if member == "bare":
            return sp.factor(stored - 2 * g * mgg * total_source)
        # mix in [0,1] is the fraction of the core q-work Legendre-subtracted.
        return sp.factor(stored - 2 * (g + mix * q) * mgg * total_source)
    raise KeyError(ensemble)


FORCE_COEFFICIENTS = OrderedDict(
    V=fixed_value_coefficient("conjugate"),
    M=same_field_coefficient("M"),
    J=same_field_coefficient("J"),
    MIXED=same_field_coefficient("MIXED"),
)
WRONG_FUNCTIONALS = OrderedDict(
    V=fixed_value_coefficient("bare"),
    M=same_field_coefficient("M", "bare"),
    J=same_field_coefficient("J", "bare"),
    MIXED=same_field_coefficient("MIXED", "bare"),
)

EXPECTED_FORCE_IDENTITIES = OrderedDict(
    V=mgg * phi**2 / sgg**2,
    M=mgg * (q**2 - g**2),
    J=-mgg * (j + g)**2,
    MIXED=mgg * ((1 - 2 * mix) * q**2 - 2 * mix * q * g - g**2),
)
for _name in FORCE_COEFFICIENTS:
    if not zq(FORCE_COEFFICIENTS[_name] - EXPECTED_FORCE_IDENTITIES[_name]):
        raise AssertionError(f"functional derivation mismatch:{_name}")


NEUTRAL_SIGNS = OrderedDict(
    V="POSITIVE_DEFINITE",
    M="INDEFINITE",
    J="NEGATIVE_DEFINITE",
    MIXED="RANGE_NEGATIVE_NULL_POSITIVE",
)


def derive_green_power(n: int) -> int:
    p = sp.symbols("p")
    roots = sp.solve(p * (p - n + 2), p)
    decay = [root for root in roots if root.is_real and root.is_positive]
    if len(decay) != 1:
        raise AssertionError(f"radial Green branch:{roots}")
    return int(decay[0])


def force_expression(A: sp.Expr, n: int = 3) -> tuple[sp.Expr, int]:
    U = s1 * s2 * A / (4 * sp.pi * R**derive_green_power(n))
    Fout = sp.factor(-sp.diff(U, R))
    return Fout, -int(Fout.as_powers_dict()[R])


FORCES = OrderedDict((name, force_expression(A)[0]) for name, A in FORCE_COEFFICIENTS.items())
FORCE_POWER = force_expression(FORCE_COEFFICIENTS["J"])[1]
MIXED_ENDPOINTS = (sp.factor(FORCE_COEFFICIENTS["MIXED"].subs(mix, 0)),
                   sp.factor(FORCE_COEFFICIENTS["MIXED"].subs(mix, 1)))
MIXED_ZERO = sp.factor(sp.solve(FORCE_COEFFICIENTS["MIXED"], mix)[0])


# ---------------------------------------------------------------------------
# Units restored.
# ---------------------------------------------------------------------------

Dim = tuple[int, int, int]
ZDIM: Dim = (0, 0, 0)
LDIM: Dim = (1, 0, 0)
EDIM: Dim = (2, -2, 1)
FDIM: Dim = (1, -2, 1)


def dadd(*dims: Dim) -> Dim:
    return tuple(sum(x) for x in zip(*dims))  # type: ignore[return-value]


def dscale(n: int, dim: Dim) -> Dim:
    return tuple(n * x for x in dim)  # type: ignore[return-value]


SYMBOL_DIMS: dict[sp.Symbol, Dim] = {
    b: (-1, -2, 1), c: (0, -2, 1), k: (1, -2, 1), km: EDIM,
    zeta: ZDIM, reb: (-2, 2, -1), ell: LDIM, sgg: (-2, 2, -1),
    phi: ZDIM, q: EDIM, j: EDIM, g: EDIM, mix: ZDIM,
    s1: ZDIM, s2: ZDIM, R: LDIM,
}


def infer_dim(expr: sp.Expr) -> Dim:
    expr = sp.sympify(expr)
    if expr.is_Number:
        return ZDIM
    if expr.is_Symbol:
        return SYMBOL_DIMS[expr]
    if expr.is_Add:
        got = [infer_dim(x) for x in expr.args if x != 0]
        if not got or any(x != got[0] for x in got[1:]):
            raise ValueError(f"dimension sum:{expr}:{got}")
        return got[0]
    if expr.is_Mul:
        return dadd(*(infer_dim(x) for x in expr.args))
    if expr.is_Pow and expr.exp.is_Integer:
        return dscale(int(expr.exp), infer_dim(expr.base))
    raise ValueError(f"dimension node:{expr}")


def dimensional_check(green_dim: Dim = (-1, 0, 0), xi_scale: sp.Expr = ell) -> bool:
    terms = []
    for A in FORCE_COEFFICIENTS.values():
        terms.extend(x for x in sp.Add.make_args(sp.expand(A)) if x != 0)
    try:
        dims = [infer_dim(x) for x in terms]
        xi_dim = infer_dim(xi_scale)
    except (KeyError, ValueError):
        return False
    Adim = dadd(EDIM, LDIM)
    return (all(x == Adim for x in dims) and
            all(dadd(x, green_dim) == EDIM for x in dims) and
            all(dadd(x, dscale(-2, LDIM)) == FDIM for x in dims) and
            xi_dim == LDIM and infer_dim(zg) == ZDIM and infer_dim(zb) == ZDIM)


# ---------------------------------------------------------------------------
# Q-COMBINE neutral total facts and range controls.
# ---------------------------------------------------------------------------

COMBINED = OrderedDict(
    DIRICHLET_VALUE=OrderedDict(
        replace=sp.factor(mgg * phi**2 / sgg**2),
        add=sp.factor(mgg * phi**2 / sgg**2),
    ),
    FIXED_MONOPOLE=OrderedDict(
        replace=sp.factor(mgg * q**2),
        add=sp.factor(FORCE_COEFFICIENTS["M"]),
    ),
    FIXED_SOURCE=OrderedDict(
        replace=sp.factor(-mgg * j**2),
        add=sp.factor(FORCE_COEFFICIENTS["J"]),
    ),
    MIXED=OrderedDict(
        replace=sp.factor(mgg * (1 - 2 * mix) * q**2),
        add=sp.factor(FORCE_COEFFICIENTS["MIXED"]),
    ),
)

COMBINE_FACTS = OrderedDict(
    DIRICHLET_VALUE=OrderedDict(replace=OutcomeFact.POSITIVE_R2,
                                add=OutcomeFact.POSITIVE_R2),
    FIXED_MONOPOLE=OrderedDict(replace=OutcomeFact.POSITIVE_R2,
                               add=OutcomeFact.OUTCOME_NOT_INVARIANT),
    FIXED_SOURCE=OrderedDict(replace=OutcomeFact.NEGATIVE_R2,
                             add=OutcomeFact.NEGATIVE_R2),
    MIXED=OrderedDict(replace=OutcomeFact.OUTCOME_NOT_INVARIANT,
                      add=OutcomeFact.OUTCOME_NOT_INVARIANT),
)
OVERALL_NEUTRAL_OUTCOME = OutcomeFact.OUTCOME_NOT_INVARIANT


def outcome_over_samples(values: Iterable[sp.Expr]) -> str:
    signs = []
    for value in values:
        value = sp.simplify(value)
        signs.append(0 if value == 0 else (1 if value > 0 else -1))
    return "CONSTANT_OUTCOME" if len(set(signs)) == 1 else "outcome_not_invariant"


def combine_controls(kind: str) -> str:
    x = sp.symbols("x", real=True)
    if kind == "sign_flip":
        return outcome_over_samples(((-1), 0, 1))
    if kind == "touch_zero":
        return outcome_over_samples(((x - 1)**2).subs(x, v) for v in (0, 1, 2))
    if kind == "known_subdominant":
        # |negative contribution| <= positive/2 is a derived domain bound.
        return outcome_over_samples((1 - v for v in (sp.Rational(1, 4), sp.Rational(1, 2))))
    raise KeyError(kind)


# ---------------------------------------------------------------------------
# Q-MAG and §5 hooks, all neutral.
# ---------------------------------------------------------------------------

def qmag_census(injected_free_factor: bool = False, detector_live: bool = True,
                actual_geometry: bool = True) -> MagnitudeFact:
    factors = ({"c_a", "c_xi"} if actual_geometry else set())
    # a=c_a*r_e and xi_A=c_xi*a are not fixed by ke^2/a=mc^2.
    if injected_free_factor:
        factors.add("c_injected")
    detected = factors if detector_live else factors - {"c_injected"}
    return MagnitudeFact.FREE_FACTOR if detected else MagnitudeFact.NO_FREE_FACTOR


QMAG_FACT = qmag_census()
QMAG_ERROR = "far-field relative O(a/R); core normalization interval c_a,c_xi in (0,infinity) at Tier-A"
TIER_FACT = TierFact.CONDITIONAL


D_RHO = sp.factor(b * k - c**2)
DENSITY_KERNEL = sp.factor(b * zg**2 / D_RHO)
HOOKS = OrderedDict(
    density_dependence="NO(no_local_prediction: B_eff=rho_B0^2/chi_c is background-only; K_h,C_hu,z_g lack local laws)",
    radial_monopole="UNDETERMINED(core source/conormal is not fixed by Z2 orientation)",
    universal_quantization="NO(continuous c_a,c_xi modulus survives)",
    close_range="UNDETERMINED(out_of_scope_R_comparable_to_r_e)",
)


# ---------------------------------------------------------------------------
# Sealed §4 adjudicator and exhaustive Cartesian table.
# ---------------------------------------------------------------------------

BC_DOMAIN = ("DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED",
             "UNDETERMINED_ANALYTICALLY")
OUTCOME_DOMAIN = ("POSITIVE_R2", "NEGATIVE_R2", "NULL", "POSITIVE_WRONG_RANGE",
                  "NEGATIVE_WRONG_RANGE", "outcome_not_invariant")
VARIANT_DOMAIN = ("replace", "add", "both", "variant_unresolved")
MAG_DOMAIN = ("magnitude_no_free_factor", "magnitude_free_factor")
TIER_DOMAIN = ("tier_A_gaps_closed", "tier_A_conditional")

@dataclass(frozen=True)
class LiveLandingFacts:
    """The complete typed neutral input emitted before §4 adjudication."""

    qbc_status: BCStatus
    combine_facts: Mapping[str, Mapping[str, OutcomeFact]]
    overall_outcome: OutcomeFact
    variant_realization: VariantFact
    magnitude: MagnitudeFact
    tier: TierFact
    internal_inconsistency: InconsistencyStatus

    def __post_init__(self) -> None:
        if not _live_facts_are_typed(self):
            raise TypeError("LiveLandingFacts accepts only the enumerated neutral fact domain")

    def neutral_tokens(self) -> tuple[Any, ...]:
        outcomes = tuple(
            self.combine_facts[cls][variant]
            for cls in ADMISSIBLE_CLASSES
            for variant in ("replace", "add")
        )
        return (self.qbc_status, *outcomes, self.overall_outcome,
                self.variant_realization, self.magnitude, self.tier,
                self.internal_inconsistency)


def _live_facts_are_typed(facts: LiveLandingFacts) -> bool:
    if (type(facts.qbc_status) is not BCStatus or
            type(facts.qbc_status.kind) is not BCClass or
            not all(type(x) is str and x in CLOSURE_GAP_DOMAIN
                    for x in facts.qbc_status.missing_closure)):
        return False
    if tuple(facts.combine_facts) != ADMISSIBLE_CLASSES:
        return False
    for cls in ADMISSIBLE_CLASSES:
        row = facts.combine_facts[cls]
        if tuple(row) != ("replace", "add"):
            return False
        if any(type(row[key]) is not OutcomeFact for key in ("replace", "add")):
            return False
    return (type(facts.overall_outcome) is OutcomeFact and
            type(facts.variant_realization) is VariantFact and
            type(facts.magnitude) is MagnitudeFact and
            type(facts.tier) is TierFact and
            type(facts.internal_inconsistency) is InconsistencyStatus and
            all(type(x) is str and x in AMENDMENT_SECTOR_DOMAIN
                for x in facts.internal_inconsistency.sectors))


def live_adjudication_case(facts: LiveLandingFacts) -> dict[str, Any]:
    """Derive every §4 predicate from the live typed neutral facts."""
    if not _live_facts_are_typed(facts):
        raise TypeError("live landing facts are not in the neutral fact domain")

    all_outcomes = tuple(
        facts.combine_facts[cls][variant]
        for cls in ADMISSIBLE_CLASSES
        for variant in ("replace", "add")
    )
    computed_overall = (all_outcomes[0] if len(set(all_outcomes)) == 1
                        else OutcomeFact.OUTCOME_NOT_INVARIANT)
    if facts.overall_outcome is not computed_overall:
        raise ValueError("Q-COMBINE overall outcome disagrees with its live class facts")

    def summarize(variant: str) -> OutcomeFact:
        values = tuple(facts.combine_facts[cls][variant] for cls in ADMISSIBLE_CLASSES)
        return values[0] if len(set(values)) == 1 else facts.overall_outcome

    realized_variants = ((facts.variant_realization.value,)
                         if facts.variant_realization in {VariantFact.REPLACE, VariantFact.ADD}
                         else ("replace", "add"))
    realized_outcomes = tuple(
        facts.combine_facts[cls][variant]
        for cls in ADMISSIBLE_CLASSES
        for variant in realized_variants
    )
    mixed_outcomes = tuple(facts.combine_facts["MIXED"][variant]
                           for variant in realized_variants)

    return dict(
        qbc=facts.qbc_status,
        replace_outcome=summarize("replace"),
        add_outcome=summarize("add"),
        variant=facts.variant_realization,
        magnitude=facts.magnitude,
        tier=facts.tier,
        internal=bool(facts.internal_inconsistency.sectors),
        all_classes_agree=(len(set(realized_outcomes)) == 1),
        mixed_range_invariant=all(
            value is not OutcomeFact.OUTCOME_NOT_INVARIANT for value in mixed_outcomes
        ),
    )


def _neutral_token_value(token: Any) -> str:
    if type(token) is BCStatus:
        return str(token)
    if type(token) in {OutcomeFact, VariantFact, MagnitudeFact, TierFact}:
        return token.value
    if type(token) is InconsistencyStatus:
        return str(token)
    raise TypeError(type(token).__name__)


def landing_ownership_guard(facts: LiveLandingFacts, emitted_landing: str,
                            upstream_tokens: tuple[Any, ...] | None = None) -> bool:
    """Prove that neutral upstream facts own no landing and §4 owns the emitted one."""
    if not _live_facts_are_typed(facts):
        return False
    expected_tokens = facts.neutral_tokens()
    carried_tokens = expected_tokens if upstream_tokens is None else upstream_tokens
    if (len(carried_tokens) != len(expected_tokens) or
            any(got is not want for got, want in zip(carried_tokens, expected_tokens))):
        return False
    try:
        tuple(_neutral_token_value(token) for token in carried_tokens)
    except TypeError:
        return False
    return emitted_landing == section4_adjudicate(facts)


def section4_adjudicate(live_facts: LiveLandingFacts | None = None, *,
                        qbc: str | BCStatus | None = None,
                        replace_outcome: str | OutcomeFact | None = None,
                        add_outcome: str | OutcomeFact | None = None,
                        variant: str | VariantFact | None = None,
                        magnitude: str | MagnitudeFact | None = None,
                        tier: str | TierFact | None = None,
                        internal: bool | None = None,
                        all_classes_agree: bool | None = None,
                        mixed_range_invariant: bool | None = None,
                        precedence_mutation: bool = False) -> str:
    """First-match implementation of directive §4; target appears only here."""
    if live_facts is not None:
        if any(value is not None for value in
               (qbc, replace_outcome, add_outcome, variant, magnitude, tier,
                internal, all_classes_agree, mixed_range_invariant)):
            raise TypeError("live facts and restated adjudication fields cannot be mixed")
        return section4_adjudicate(**live_adjudication_case(live_facts),
                                   precedence_mutation=precedence_mutation)
    if any(value is None for value in
           (qbc, replace_outcome, add_outcome, variant, magnitude, tier,
            internal, all_classes_agree, mixed_range_invariant)):
        raise TypeError("incomplete adjudication case")

    qbc = qbc.kind.value if type(qbc) is BCStatus else qbc
    if internal:
        return "NO_GO(sector)"

    if variant == "replace":
        comparison = replace_outcome
        variant_resolved = True
    elif variant == "add":
        comparison = add_outcome
        variant_resolved = True
    else:
        comparison = replace_outcome if replace_outcome == add_outcome else None
        variant_resolved = comparison is not None

    if qbc == "UNDETERMINED_ANALYTICALLY":
        class_resolved = all_classes_agree
    elif qbc == "MIXED":
        class_resolved = tier == "tier_A_gaps_closed" and mixed_range_invariant
    else:
        class_resolved = tier == "tier_A_gaps_closed"
    magnitude_determined = comparison is not None and comparison != "outcome_not_invariant"
    unconditional = class_resolved and variant_resolved and magnitude_determined

    # Mutation used only by the atomic precedence tooth: it illegally routes
    # the class gap before allowing the unconditional invariance route.
    if precedence_mutation and qbc == "UNDETERMINED_ANALYTICALLY":
        return "R1_REQUIRED(bc_selection)"

    if unconditional:
        if comparison == "POSITIVE_R2":
            return ("SIGN_EARNED" if magnitude == "magnitude_no_free_factor"
                    else "R1_REQUIRED(magnitude)")
        if comparison == "NEGATIVE_R2":
            return "MECHANISM_FALSIFIED(wrong_sign)"
        if comparison in {"POSITIVE_WRONG_RANGE", "NEGATIVE_WRONG_RANGE"}:
            return "MECHANISM_FALSIFIED(wrong_range)"
        if comparison == "NULL":
            return "R1_REQUIRED(subleading)"

    if qbc == "UNDETERMINED_ANALYTICALLY" and not all_classes_agree:
        return "R1_REQUIRED(bc_selection)"
    if qbc == "MIXED" and not mixed_range_invariant:
        return "R1_REQUIRED(mixed_bc_parameters)"
    if variant in {"both", "variant_unresolved"} and replace_outcome != add_outcome:
        return "R1_REQUIRED(variant_selection)"
    if tier == "tier_A_conditional":
        return "R1_REQUIRED(sign_and_magnitude)"
    if comparison == "outcome_not_invariant":
        return "R1_REQUIRED(sign_and_magnitude)"
    return ("R1_REQUIRED(unclassified):" + repr((qbc, replace_outcome, add_outcome,
                                                  variant, magnitude, tier, internal,
                                                  all_classes_agree, mixed_range_invariant)))


def section4_oracle(case: Mapping[str, Any]) -> str:
    # Declarative second implementation used only to test totality/precedence.
    if case["internal"]:
        return "NO_GO(sector)"
    ro, ao, variant = case["replace_outcome"], case["add_outcome"], case["variant"]
    common = ro if variant == "replace" else ao if variant == "add" else (ro if ro == ao else None)
    variant_ok = variant in {"replace", "add"} or ro == ao
    if case["qbc"] == "UNDETERMINED_ANALYTICALLY":
        class_ok = case["all_classes_agree"]
    elif case["qbc"] == "MIXED":
        class_ok = case["tier"] == "tier_A_gaps_closed" and case["mixed_range_invariant"]
    else:
        class_ok = case["tier"] == "tier_A_gaps_closed"
    if class_ok and variant_ok and common != "outcome_not_invariant":
        if common == "POSITIVE_R2":
            return "SIGN_EARNED" if case["magnitude"] == "magnitude_no_free_factor" else "R1_REQUIRED(magnitude)"
        if common == "NEGATIVE_R2":
            return "MECHANISM_FALSIFIED(wrong_sign)"
        if common in {"POSITIVE_WRONG_RANGE", "NEGATIVE_WRONG_RANGE"}:
            return "MECHANISM_FALSIFIED(wrong_range)"
        if common == "NULL":
            return "R1_REQUIRED(subleading)"
    if case["qbc"] == "UNDETERMINED_ANALYTICALLY" and not case["all_classes_agree"]:
        return "R1_REQUIRED(bc_selection)"
    if case["qbc"] == "MIXED" and not case["mixed_range_invariant"]:
        return "R1_REQUIRED(mixed_bc_parameters)"
    if variant in {"both", "variant_unresolved"} and ro != ao:
        return "R1_REQUIRED(variant_selection)"
    if case["tier"] == "tier_A_conditional" or common == "outcome_not_invariant":
        return "R1_REQUIRED(sign_and_magnitude)"
    return "R1_REQUIRED(unclassified)"


def truth_table(precedence_mutation: bool = False) -> tuple[bool, int, Counter[str], str]:
    rows: list[str] = []
    counts: Counter[str] = Counter()
    exact = True
    total = 0
    for values in itertools.product(BC_DOMAIN, OUTCOME_DOMAIN, OUTCOME_DOMAIN,
                                    VARIANT_DOMAIN, MAG_DOMAIN, TIER_DOMAIN,
                                    (False, True), (False, True), (False, True)):
        keys = ("qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
                "tier", "internal", "all_classes_agree", "mixed_range_invariant")
        case = dict(zip(keys, values))
        got = section4_adjudicate(**case, precedence_mutation=precedence_mutation)
        want = section4_oracle(case)
        if precedence_mutation:
            # The dedicated witness is an invariance-route row; other rows are
            # irrelevant to proving the mutation load-bearing.
            pass
        else:
            exact &= got == want and "unclassified" not in got
        counts[got.split(":", 1)[0]] += 1
        rows.append("|".join(map(str, values)) + "|" + got)
        total += 1
    digest = hashlib.sha256("\n".join(rows).encode()).hexdigest()
    return exact, total, counts, digest


LIVE_LANDING_FACTS = LiveLandingFacts(
    qbc_status=QBC_ACTUAL.status,
    combine_facts=COMBINE_FACTS,
    overall_outcome=OVERALL_NEUTRAL_OUTCOME,
    variant_realization=VARIANT_REALIZATION,
    magnitude=QMAG_FACT,
    tier=TIER_FACT,
    internal_inconsistency=INCONSISTENCY_FACT,
)
ACTUAL_LANDING = section4_adjudicate(LIVE_LANDING_FACTS)
EXPECTED_PRODUCTION_LANDING = "R1_REQUIRED(bc_selection)"
if ACTUAL_LANDING != EXPECTED_PRODUCTION_LANDING:
    raise RuntimeError(
        "production landing discrepancy: live facts adjudicated to "
        f"{ACTUAL_LANDING}, expected {EXPECTED_PRODUCTION_LANDING}"
    )


# ---------------------------------------------------------------------------
# Target-blindness and cross-engine payload.
# ---------------------------------------------------------------------------

def walk_exprs(value: Any) -> Iterable[sp.Expr]:
    if isinstance(value, sp.Expr):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from walk_exprs(item)
    elif isinstance(value, (tuple, list)):
        for item in value:
            yield from walk_exprs(item)


def ast_target_scope_ok() -> bool:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    allowed = {"section4_adjudicate", "section4_oracle", "truth_table", "ast_target_scope_ok"}
    tokens = ("sign_earned", "mechanism_falsified", "wrong_sign", "like-repel", "like-attract")

    class Visitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.stack: list[str] = []
            self.bad: list[str] = []

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
            self.stack.append(node.name)
            self.generic_visit(node)
            self.stack.pop()

        def visit_Constant(self, node: ast.Constant) -> None:
            if isinstance(node.value, str) and any(x in node.value.lower() for x in tokens):
                if not self.stack or self.stack[-1] not in allowed:
                    self.bad.append(node.value)

    visitor = Visitor()
    visitor.visit(tree)
    return not visitor.bad


def target_blind(objects: Any) -> bool:
    forbidden = {"targetsign", "targetpower", "desiredanswer"}
    return ast_target_scope_ok() and all(
        not ({str(x).lower() for x in expr.free_symbols} & forbidden)
        for expr in walk_exprs(objects)
    )


def canonical(expr: sp.Expr) -> str:
    return sp.sstr(sp.cancel(expr)).replace("**", "^")


def symbolic_payload(mutate: bool = False) -> OrderedDict[str, Any]:
    terms: OrderedDict[str, str] = OrderedDict(
        D=canonical(D), kappa=canonical(kappa), zg=canonical(zg), zb=canonical(zb),
        muu=canonical(muu), mug=canonical(mug), mgg=canonical(mgg), mdet=canonical(mdet),
        field_scale=canonical(xi_over_h), parent_source_map=canonical(parent_to_reduced_source),
    )
    for name, expr in FORCE_COEFFICIENTS.items():
        terms[f"force.{name}"] = canonical(expr)
        terms[f"force_out.{name}"] = canonical(FORCES[name])
        terms[f"wrong_functional.{name}"] = canonical(WRONG_FUNCTIONALS[name])
    for cls, rows in COMBINED.items():
        for variant, expr in rows.items():
            terms[f"combine.{cls}.{variant}"] = canonical(expr)
    terms["mixed.endpoint0"] = canonical(MIXED_ENDPOINTS[0])
    terms["mixed.endpoint1"] = canonical(MIXED_ENDPOINTS[1])
    terms["mixed.zero"] = canonical(MIXED_ZERO)
    terms["density.kernel"] = canonical(DENSITY_KERNEL)
    if mutate:
        terms["force.J"] = canonical(FORCE_COEFFICIENTS["J"] + j * g * mgg)
    exact, total, counts, digest = truth_table()
    return OrderedDict(
        schema="PUNCTURE_DEFLECTION_SIGN_V1",
        symbolic_terms=terms,
        neutral_signs=NEUTRAL_SIGNS,
        qbc=str(QBC_ACTUAL.status),
        admissible_classes=list(ADMISSIBLE_CLASSES),
        variant_realization=VARIANT_REALIZATION,
        combine_facts=COMBINE_FACTS,
        overall_outcome=OVERALL_NEUTRAL_OUTCOME,
        magnitude=QMAG_FACT,
        tier=TIER_FACT,
        hooks=HOOKS,
        internal_inconsistency=str(INCONSISTENCY_FACT),
        force_power=FORCE_POWER,
        truth_total=total,
        truth_exact=exact,
        truth_digest=digest,
        truth_counts=dict(sorted(counts.items())),
        landing=ACTUAL_LANDING,
    )


PARSE_LOCALS = {str(x): x for x in (b, c, k, km, zeta, reb, ell, sgg, phi, q, j, g, mix)}
PARSE_LOCALS.update({"pi": sp.pi, "Pi": sp.pi, "s1": s1, "s2": s2, "R": R})


def parse_other(text: str) -> sp.Expr:
    return sp.sympify(text.replace("^", "**"), locals=PARSE_LOCALS)


def payloads_equal(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    if set(left) != set(right) or set(left["symbolic_terms"]) != set(right["symbolic_terms"]):
        return False
    for key, expr in left["symbolic_terms"].items():
        if not zq(parse_other(expr) - parse_other(right["symbolic_terms"][key])):
            return False
    return all(left[key] == right[key] for key in left if key != "symbolic_terms")


def wolfram_payload(mutate: bool = False) -> Mapping[str, Any]:
    env = {**os.environ, "PUNCTURE_JSON_ONLY": "1"}
    if mutate:
        env["PUNCTURE_PAYLOAD_MUTATION"] = "1"
    run = subprocess.run(["math", "-script", str(WL_PATH)], cwd=ROOT, env=env,
                         text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         timeout=240, check=False)
    if run.returncode != 0:
        raise AssertionError(f"Mathematica exit {run.returncode}:\n{run.stdout}")
    marker = next((x for x in run.stdout.splitlines() if x.startswith("JSON_PAYLOAD=")), None)
    if marker is None:
        raise AssertionError(f"Mathematica payload absent:\n{run.stdout}")
    return json.loads(marker.split("=", 1)[1])


# ---------------------------------------------------------------------------
# Atomic teeth and production-path mutations.
# ---------------------------------------------------------------------------

TOOTH_ORDER = (
    "FIELD_IDENTITY_UNITS", "FIELD_PARENT_MAP", "FIELD_LIVE_QH",
    "ACTION_TRANSCRIPTION", "AMEND_REPLACE", "AMEND_ADD", "SHOLD_SCOPE",
    "ZERO_LEDGER", "MATRIX_DETERMINANT", "BC_ACTUAL_GAP", "BC_FREE_CONTROL",
    "BC_VALUE_CONTROL", "BC_MONOPOLE_CONTROL", "BC_MIXED_CONTROL",
    "FORCE_V_FUNCTIONAL", "FORCE_M_HWORK", "FORCE_J_FUNCTIONAL",
    "FORCE_MIXED_FUNCTIONAL", "MIXED_FULL_RANGE", "FALLOFF", "UNITS_RESTORED",
    "COMBINE_REPLACE", "COMBINE_ADD", "NO_DOUBLE_COUNT", "RANGE_SIGN_FLIP",
    "RANGE_TOUCH_ZERO", "RANGE_SUBDOMINANT", "MAG_FREE_FACTOR",
    "DENSITY_HOOK", "MONOPOLE_HOOK", "MODULUS_HOOK", "VERDICT_TOTALITY",
    "VERDICT_PRECEDENCE", "LANDING_OWNERSHIP", "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
)


def local_checks(mutation: str | None, dual_ok: bool) -> OrderedDict[str, bool]:
    transcript = TRANSCRIPT
    hessian = ACTION_HESSIAN
    if mutation == "ACTION_TRANSCRIPTION":
        changed = list(transcript.action_terms)
        changed[-1] = "-2*C_hu*grad u_L.grad h"
        transcript = replace(transcript, action_terms=tuple(changed))
        hessian = sp.hessian(b * gu**2 / 2 + 2 * c * gu * ghgrad + k * ghgrad**2 / 2,
                             (gu, ghgrad))
    faithful, _ = source_faithfulness(transcript, hessian)

    xi_scale = ell**2 if mutation == "FIELD_IDENTITY_UNITS" else ell
    parent_map = 2 * parent_to_reduced_source if mutation == "FIELD_PARENT_MAP" else parent_to_reduced_source
    qh_map = 0 if mutation == "FIELD_LIVE_QH" else parent_to_reduced_source

    rep = REPLACE
    add = ADD
    if mutation == "AMEND_REPLACE":
        rep = replace(rep, new_rows=("illegal_new_row",))
    if mutation == "AMEND_ADD":
        add = replace(add, new_rows=add.new_rows + ("second_row",))
    if mutation == "SHOLD_SCOPE":
        rep = replace(rep, shold_scope="r_B-1/2 and h")
    if mutation == "ZERO_LEDGER":
        add = replace(add, untouched_zero_rows=ZERO_LEDGER[:-1])
    sectors, _ = amendment_consistency(rep, add)

    det_used = sp.factor(mdet * 2) if mutation == "MATRIX_DETERMINANT" else mdet

    actual = ACTUAL_EVIDENCE
    free = FREE_CONTROL
    value = VALUE_CONTROL
    mono = MONOPOLE_CONTROL
    mixed_ev = MIXED_CONTROL
    if mutation == "BC_ACTUAL_GAP":
        actual = replace(actual, missing_closure=())
    if mutation == "BC_FREE_CONTROL":
        free = replace(free, essential_value=True, admissible_variation=False)
    if mutation == "BC_VALUE_CONTROL":
        value = replace(value, essential_value=False, admissible_variation=True,
                        missing_closure=("lost_value",))
    if mutation == "BC_MONOPOLE_CONTROL":
        mono = replace(mono, fixed_conormal=False, missing_closure=("lost_flux",))
    if mutation == "BC_MIXED_CONTROL":
        mixed_ev = replace(mixed_ev, mixed_relation=False, missing_closure=("lost_mix",))

    selected = OrderedDict(FORCE_COEFFICIENTS)
    if mutation == "FORCE_V_FUNCTIONAL":
        selected["V"] = WRONG_FUNCTIONALS["V"]
    if mutation == "FORCE_M_HWORK":
        selected["M"] = WRONG_FUNCTIONALS["M"]
    if mutation == "FORCE_J_FUNCTIONAL":
        selected["J"] = WRONG_FUNCTIONALS["J"]
    if mutation == "FORCE_MIXED_FUNCTIONAL":
        selected["MIXED"] = WRONG_FUNCTIONALS["MIXED"]
    mixed_range = (MIXED_ENDPOINTS[0], MIXED_ENDPOINTS[1])
    if mutation == "MIXED_FULL_RANGE":
        mixed_range = (MIXED_ENDPOINTS[0], MIXED_ENDPOINTS[0])

    n = 4 if mutation == "FALLOFF" else 3
    green_dim = (-2, 0, 0) if mutation == "UNITS_RESTORED" else (-1, 0, 0)

    combine_rep = COMBINED["DIRICHLET_VALUE"]["replace"]
    combine_add = COMBINED["FIXED_MONOPOLE"]["add"]
    if mutation == "COMBINE_REPLACE":
        combine_rep += mgg * g**2
    if mutation == "COMBINE_ADD":
        combine_add = mgg * q**2
    rebuilt_m = sp.factor(mgg * (q + g)**2 - 2 * g * mgg * (q + g))
    if mutation == "NO_DOUBLE_COUNT":
        rebuilt_m += mgg * g**2

    signflip = combine_controls("sign_flip")
    touch = combine_controls("touch_zero")
    subdom = combine_controls("known_subdominant")
    if mutation == "RANGE_SIGN_FLIP":
        signflip = "CONSTANT_OUTCOME"
    if mutation == "RANGE_TOUCH_ZERO":
        touch = "CONSTANT_OUTCOME"
    if mutation == "RANGE_SUBDOMINANT":
        subdom = "outcome_not_invariant"

    mag_detector = mutation != "MAG_FREE_FACTOR"
    mag_control = qmag_census(True, detector_live=mag_detector, actual_geometry=False)
    density = DENSITY_KERNEL if mutation != "DENSITY_HOOK" else sp.factor(zg**2 / D_RHO)
    monopole_hook = HOOKS["radial_monopole"] if mutation != "MONOPOLE_HOOK" else "YES"
    modulus_hook = HOOKS["universal_quantization"] if mutation != "MODULUS_HOOK" else "YES"

    truth_exact, _, _, _ = truth_table()
    if mutation == "VERDICT_TOTALITY":
        truth_exact = False
    witness = dict(qbc="UNDETERMINED_ANALYTICALLY", replace_outcome="POSITIVE_R2",
                   add_outcome="POSITIVE_R2", variant="both",
                   magnitude="magnitude_no_free_factor", tier="tier_A_conditional",
                   internal=False, all_classes_agree=True, mixed_range_invariant=True)
    precedence = section4_adjudicate(**witness,
                                     precedence_mutation=(mutation == "VERDICT_PRECEDENCE"))
    expected_precedence = section4_adjudicate(**witness)
    upstream_objects: Any = (FORCE_COEFFICIENTS, COMBINED, QMAG_FACT, TIER_FACT)
    ownership_tokens = LIVE_LANDING_FACTS.neutral_tokens()
    if mutation == "LANDING_OWNERSHIP":
        ownership_tokens = (*ownership_tokens, ACTUAL_LANDING)
    stage_objects: Any = upstream_objects
    if mutation == "TARGET_BLINDNESS":
        stage_objects = (*upstream_objects, sp.Symbol("targetsign") * FORCE_COEFFICIENTS["V"])

    result = OrderedDict(
        FIELD_IDENTITY_UNITS=(infer_dim(xi_scale) == LDIM),
        FIELD_PARENT_MAP=zq(parent_map - 1),
        FIELD_LIVE_QH=zq(qh_map - 1),
        ACTION_TRANSCRIPTION=faithful,
        AMEND_REPLACE=("replace_ledger" not in sectors),
        AMEND_ADD=("add_ledger" not in sectors),
        SHOLD_SCOPE=("S_hold" not in sectors),
        ZERO_LEDGER=("zero_ledger" not in sectors),
        MATRIX_DETERMINANT=zq(det_used - zg**2 / D),
        BC_ACTUAL_GAP=(classify_bc(actual).status.kind is BCClass.UNDETERMINED_ANALYTICALLY),
        BC_FREE_CONTROL=(classify_bc(free).status.kind is BCClass.FIXED_SOURCE),
        BC_VALUE_CONTROL=(classify_bc(value).status.kind is BCClass.DIRICHLET_VALUE),
        BC_MONOPOLE_CONTROL=(classify_bc(mono).status.kind is BCClass.FIXED_MONOPOLE),
        BC_MIXED_CONTROL=(classify_bc(mixed_ev).status.kind is BCClass.MIXED),
        FORCE_V_FUNCTIONAL=(zq(selected["V"] - FORCE_COEFFICIENTS["V"]) and
                            not zq(WRONG_FUNCTIONALS["V"] - FORCE_COEFFICIENTS["V"])),
        FORCE_M_HWORK=(zq(selected["M"] - FORCE_COEFFICIENTS["M"]) and
                       not zq(WRONG_FUNCTIONALS["M"] - FORCE_COEFFICIENTS["M"])),
        FORCE_J_FUNCTIONAL=(zq(selected["J"] - FORCE_COEFFICIENTS["J"]) and
                            not zq(WRONG_FUNCTIONALS["J"] - FORCE_COEFFICIENTS["J"])),
        FORCE_MIXED_FUNCTIONAL=(zq(selected["MIXED"] - FORCE_COEFFICIENTS["MIXED"]) and
                                not zq(WRONG_FUNCTIONALS["MIXED"] - FORCE_COEFFICIENTS["MIXED"])),
        MIXED_FULL_RANGE=(zq(mixed_range[0] - MIXED_ENDPOINTS[0]) and
                          zq(mixed_range[1] - MIXED_ENDPOINTS[1]) and
                          not zq(MIXED_ENDPOINTS[0] - MIXED_ENDPOINTS[1])),
        FALLOFF=(force_expression(FORCE_COEFFICIENTS["J"], n)[1] == 2),
        UNITS_RESTORED=dimensional_check(green_dim),
        COMBINE_REPLACE=zq(combine_rep - COMBINED["DIRICHLET_VALUE"]["replace"]),
        COMBINE_ADD=zq(combine_add - COMBINED["FIXED_MONOPOLE"]["add"]),
        NO_DOUBLE_COUNT=zq(rebuilt_m - COMBINED["FIXED_MONOPOLE"]["add"]),
        RANGE_SIGN_FLIP=(signflip == "outcome_not_invariant"),
        RANGE_TOUCH_ZERO=(touch == "outcome_not_invariant"),
        RANGE_SUBDOMINANT=(subdom == "CONSTANT_OUTCOME"),
        MAG_FREE_FACTOR=(QMAG_FACT == "magnitude_free_factor" and mag_control == "magnitude_free_factor"),
        DENSITY_HOOK=(zq(density - b * zg**2 / (b * k - c**2))),
        MONOPOLE_HOOK=monopole_hook.startswith("UNDETERMINED"),
        MODULUS_HOOK=modulus_hook.startswith("NO("),
        VERDICT_TOTALITY=truth_exact,
        VERDICT_PRECEDENCE=(precedence == expected_precedence),
        LANDING_OWNERSHIP=landing_ownership_guard(
            LIVE_LANDING_FACTS, ACTUAL_LANDING, ownership_tokens
        ),
        TARGET_BLINDNESS=target_blind(stage_objects),
        DUAL_ENGINE_TERMS=dual_ok,
    )
    return result


def mutation_campaign() -> OrderedDict[str, str]:
    local = symbolic_payload()
    other = wolfram_payload(False)
    dual_ok = payloads_equal(local, other)
    baseline = local_checks(None, dual_ok)
    if not all(baseline.values()):
        _, failures = source_faithfulness()
        raise AssertionError(f"baseline failures:{[k for k,v in baseline.items() if not v]};source={failures}")
    outcomes: OrderedDict[str, str] = OrderedDict()
    for tooth in TOOTH_ORDER:
        if tooth == "DUAL_ENGINE_TERMS":
            mutated = local_checks(None, payloads_equal(local, wolfram_payload(True)))
        else:
            mutated = local_checks(tooth, dual_ok)
        failed = [name for name, ok in mutated.items() if not ok]
        if failed != [tooth]:
            raise AssertionError(f"atomic mutation {tooth} failed {failed}")
        outcomes[tooth] = f"FIRED_AT_{tooth}"
    return outcomes


def print_report(teeth: Mapping[str, str]) -> None:
    print("PUNCTURE_DEFLECTION_ELECTRIC_SIGN — SYMPY")
    print("Q-FIELD:")
    print("  identity=xi_w=ell*h; [xi_w]=L; held_datum=h_A=xi_w/ell=P0 H; [h_A]=1")
    print("  H=f0*h+H_perp; f0=1/(ell*cosh(w/ell)^2); N0=8/(3*ell); h=P0 H")
    print("  u_L is distinct ([u_L]=L); C_hu mixes u_L and h but does not identify them")
    print("  Q_chi*H live map: f0(0)=1/ell; -J_m Q_chi H -> -g_chih Q_chi h")
    print("Q-AMEND:")
    print("  baseline=NO_NATIVE_CLAMP; S_hold=r_B-1/2 only")
    print("  REPLACE=existing h-source BC retyped as candidate core holding BC; new_rows=0")
    print("  ADD=existing h-source unchanged + exactly one core h-holding row")
    print(f"  nondeclared_zero_rows={len(ZERO_LEDGER)} preserved; internal_inconsistency={INCONSISTENCY_FACT}")
    print("Q-BC:")
    print("  G0 freezes Sigma and U2 leaves B unresolved: committed bare math does not select the class")
    print(f"  actual={QBC_ACTUAL.status}; stationarity={QBC_ACTUAL.stationarity}; "
          f"holding_curvature={QBC_ACTUAL.holding_curvature}; barrier={QBC_ACTUAL.topological_barrier}")
    print("  exterior_solution=h_A*a/r; exterior holding curvature positive; core s->h_A map/barrier absent")
    for row in QBC_CONTROLS:
        print(f"  control.{row.name}={row.status}")
    print("  admissible_classes=" + ",".join(ADMISSIBLE_CLASSES))
    print("Q-FORCE:")
    print(f"  D={sx(D)}; kappa={sx(kappa)}; m={m}; det(m)={sx(mdet)}")
    print("  actual-field derivation: total h-channel source is reaction/core plus committed g; no u_L-datum relabel")
    for name, A in FORCE_COEFFICIENTS.items():
        force, power = force_expression(A)
        print(f"  {name}: functional_coefficient={sx(A)}; sign={NEUTRAL_SIGNS[name]}; "
              f"U=s1*s2*A/(4*pi*R); F_out={sx(force)}; falloff=1/R^{power}")
        print(f"    wrong_functional={sx(WRONG_FUNCTIONALS[name])}; changed={not zq(A-WRONG_FUNCTIONALS[name])}")
    print(f"  MIXED admissible range=0<=mix<=1; endpoints=({sx(MIXED_ENDPOINTS[0])},{sx(MIXED_ENDPOINTS[1])}); "
          f"zero={sx(MIXED_ZERO)}")
    print("Q-COMBINE:")
    for cls, variants in COMBINED.items():
        print(f"  {cls}: replace={sx(variants['replace'])} [{COMBINE_FACTS[cls]['replace']}]; "
              f"add={sx(variants['add'])} [{COMBINE_FACTS[cls]['add']}]")
    print(f"  variant_realization={VARIANT_REALIZATION}; overall={OVERALL_NEUTRAL_OUTCOME}")
    print("  controls: sign_flip=outcome_not_invariant; touch_zero=outcome_not_invariant; "
          "known_subdominant=CONSTANT_OUTCOME")
    print("Q-MAG:")
    print("  a=c_a*r_e; r_e=k_e*e^2/(m_e*c^2); xi_A=c_xi*a; census={c_a,c_xi}")
    print(f"  fact={QMAG_FACT}; uncertainty={QMAG_ERROR}")
    print(f"NEUTRAL_FACTS: outcome={OVERALL_NEUTRAL_OUTCOME}; magnitude={QMAG_FACT}; tier={TIER_FACT}; "
          f"variant-realization={VARIANT_REALIZATION}; internal_inconsistency={INCONSISTENCY_FACT}")
    print("SECTION5_HOOKS:")
    for key, value in HOOKS.items():
        print(f"  {key}={value}")
    print(f"  density_coefficient(B_eff,K_h,C_hu)={sx(DENSITY_KERNEL)}; D(rho)={sx(D_RHO)}; require D(rho)>0; "
          "allow(cancellation|instability|no_local_prediction)")
    exact, total, counts, digest = truth_table()
    print(f"SECTION4_TRUTH_TABLE: cells={total}; exhaustive={exact}; digest={digest}; classes={dict(sorted(counts.items()))}")
    print("SECTION4_LANDING=" + ACTUAL_LANDING)
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
            mutate = os.environ.get("PUNCTURE_PAYLOAD_MUTATION") == "1"
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
