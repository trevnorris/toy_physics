#!/usr/bin/env python3
"""Target-blind SymPy check for the G0 leftover-longitudinal scalar.

The neutral calculation ends with mathematical quadratic-form signatures.
Only ``section5_*`` converts those signatures to the electric target language.
"""

from __future__ import annotations

import argparse
import ast
import itertools
import json
import os
import subprocess
import sys
from collections import OrderedDict
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DIRECTIVE = HERE / "directive_leftover_scalar_electric_sign.md"
G0 = HERE / "g0_closure_card_v0.md"
AUDIT = ROOT / "docs/model_definition_audit.md"
HANDOFF = ROOT / "docs/em_analog_next_phase_handoff.md"
SIM_SPEC = ROOT / "docs/two_throat_simulation_handoff_spec.md"
PATH36 = ROOT / "software/stage1_solver/reports/pathA_36_c5_phase_potential.md"
PATH38 = ROOT / "software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md"
PATH38_YAML = ROOT / "software/stage1_solver/reports/pathA_38_results.yaml"
PATH39 = ROOT / "software/stage1_solver/reports/pathA_39_stage4_field_classification.md"
PATH24 = ROOT / "software/stage1_solver/reports/pathA_24_T1_wall.md"
U2_DIRECTIVE = HERE / "directive_u2_boundary_adjudication.md"
U2_REPORT = HERE / "reports/u2_boundary_adjudication_artifacts/stage_1_production_v12/production_summary.md"
WL_PATH = HERE / "leftover_scalar_electric_sign_check.wl"


def zq(expr: sp.Expr) -> bool:
    return sp.simplify(expr) == 0


def text_expr(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr)).replace("**", "^")


# ---------------------------------------------------------------------------
# A live transcription object, checked against every authorized §0 source.
# ---------------------------------------------------------------------------

ACTION_TERMS = (
    "1/2*A_eff*(dt u_L)^2",
    "1/2*M_h*(dt h)^2",
    "-1/2*B_eff*|grad u_L|^2",
    "-1/2*K_h*|grad h|^2",
    "-C_hu*grad u_L.grad h",
)

TEST_BCS = (
    ("V", "f_V=1 on spherical A; area_average_A(f_V)=1; u_L|A=s*u0*f_V; physical Q_u reacts"),
    ("M", "f_b=1/Area(A) on A, zero off A; integral_A f_b dS=1; conormal=s*q*f_b"),
    ("J", "f_b=1/Area(A) on A, zero off A; integral_A f_b dS=1; Omega adds -s*J0*integral_A(f_b*u_L)dS"),
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
    card: tuple[tuple[str, str], ...]
    relations: tuple[str, ...]
    mouth_terms: tuple[str, ...]
    native_bcs: tuple[str, ...]
    test_bcs: tuple[tuple[str, str], ...]
    zero_ledger: tuple[str, ...]


TRANSCRIPT = Transcript(
    action_terms=ACTION_TERMS,
    card=(("A_eff", "1"), ("M_h", "1"), ("B_eff", "2"),
          ("K_h", "1"), ("C_hu", "1/2"), ("k_m", "1"), ("g_chih", "1")),
    relations=("T_L=B_eff=rho_B0^2/chi_c>0", "K_h=M_h*c_E^2"),
    mouth_terms=("1/2*k_m*h^2", "-g_chih*s_i*h"),
    native_bcs=(
        "continuity(u_L)",
        "continuity(B_eff*d_n u_L+C_hu*d_n h)",
        "continuity(h)",
        "continuity(K_h*d_n h+C_hu*d_n u_L) except mouth source",
        "u_L->0 at IR",
    ),
    test_bcs=TEST_BCS,
    zero_ledger=ZERO_LEDGER,
)


def read_sources() -> dict[str, str]:
    paths = {
        "directive": DIRECTIVE, "g0": G0, "audit": AUDIT, "handoff": HANDOFF,
        "sim": SIM_SPEC, "path36": PATH36, "path38": PATH38,
        "path38_yaml": PATH38_YAML, "path39": PATH39,
        "path24": PATH24, "u2_directive": U2_DIRECTIVE, "u2_report": U2_REPORT,
    }
    return {name: path.read_text(encoding="utf-8") for name, path in paths.items()}


def parsed_g0_ledger(g0: str) -> tuple[tuple[str, ...], bool]:
    block = g0.split("## 9. Complete declared-zero ledger", 1)[1].split(
        "## 10. Instantiability, gates, and checks", 1
    )[0]
    rows: list[str] = []
    tags_ok = True
    for line in block.splitlines():
        if line.startswith("| ") and "[POSTULATE]" in line:
            rows.append(line.split("|", 2)[1].strip())
            tags_ok &= "[POSTULATE]" in line
    return tuple(rows), tags_ok


def source_faithfulness(
    transcript: Transcript = TRANSCRIPT,
    action_hessian: sp.Matrix | None = None,
) -> tuple[bool, list[str]]:
    src = read_sources()
    g0 = src["g0"]
    failures: list[str] = []

    exact_source_needles: dict[str, tuple[str, ...]] = {
        "g0": (
            "½A_eff (∂tu_L)² + ½M_h(∂th)²",
            "−½B_eff|∇u_L|² −½K_h|∇h|²",
            "−C_hu ∇u_L·∇h",
            "| `B_eff` | `(-1,-2,1)` | `2`",
            "| `C_hu` | `(0,-2,1)` | `1/2`",
            "| `K_h` | `(1,-2,1)` | `1`",
            "| `A_eff` | `(-3,0,1)` | `1`",
            "| `M_h` | `(-1,0,1)` | `1`",
            "| reduced `k_m=K_m/ℓ²` | `E` | `1`",
            "| reduced `g_χh=j_m=J_m/ℓ` | `E` | `1`",
            "K_h=M_hc_E²",
            "½K_m H(x,0)² − J_m Q_χ[r_Σ,s_i] H(x,0)",
            "η_i(k_mh−g_χh s_i)",
            "B_eff∂_nu_L+C_hu∂_nh",
            "K_h∂_nh+C_hu∂_nu_L",
            "u_L→0",
        ),
        "directive": (
            "u_L|_mouth = s·u₀",
            "far-field monopole strength (`∮` of the conormal) `= s·q` held",
            "a `−J φ` Legendre term with `J=s·J₀` held",
        ),
        "audit": ("I1", "I2", "I3", "U1", "U2", "U3"),
        "handoff": ("144/144 cells UNRESOLVED", "R1 two-throat solve"),
        "sim": ("a supplied closure card before use", "does **not** imply that the force sign must flip"),
        "path36": ("B_eff = rho_B0^2/chi_c",),
        "path38": ("p_static=2",),
        "path38_yaml": (
            "schema: pathA_38_throat_body_electric_sympy/v2",
            "N0_norm: 8/(3*ell)",
            "q_h_plus: 2*QE*tanh(b/ell)/b",
            "q_h_minus: -2*QE*tanh(b/ell)/b",
            "p_static: 2",
        ),
        "path39": ("q_L", "mass source to u_L, recorded separately from charge residues"),
        "path24": ("DeltaE_unwind = 0",),
        "u2_directive": ("144", "UNRESOLVED", "mouth"),
        "u2_report": ("{'UNRESOLVED': 144}",),
    }
    for label, needles in exact_source_needles.items():
        for needle in needles:
            if needle not in src[label]:
                failures.append(f"{label}:missing:{needle}")

    source_ledger, tags_ok = parsed_g0_ledger(g0)
    if not tags_ok or len(source_ledger) != 13:
        failures.append("g0:zero-ledger tag/count")

    expected = TRANSCRIPT
    for field in Transcript.__dataclass_fields__:
        if getattr(transcript, field) != getattr(expected, field):
            failures.append(f"internal-transcript:{field}")
    if transcript.zero_ledger != source_ledger:
        failures.append("internal-transcript:zero-ledger-vs-source")
    # The computational operator must be the Hessian of the transcribed
    # density; this catches a mirrored cross-factor drift, not only prose.
    action_hessian = action_hessian if action_hessian is not None else globals()["ACTION_HESSIAN"]
    expected_hessian = sp.Matrix([[globals()["b"], globals()["c"]], [globals()["c"], globals()["k"]]])
    if action_hessian.shape != (2, 2) or any(
        not zq(action_hessian[i, j] - expected_hessian[i, j])
        for i in range(2) for j in range(2)
    ):
        failures.append("internal-transcript:computational-action-hessian")
    return not failures, failures


# ---------------------------------------------------------------------------
# Coupled action, exact Robin response, and thermodynamic members.
# ---------------------------------------------------------------------------

b, c, k = sp.symbols("b c k", positive=True)
km, zeta, reb = sp.symbols("km zeta reb", positive=True)
suu, sgg = sp.symbols("suu sgg", positive=True)
sug = sp.symbols("sug", real=True)
u, qf, js, gh = sp.symbols("u qf js gh", positive=True)
sone, stwo = sp.symbols("sone stwo", real=True)
rr = sp.symbols("rr", positive=True)

grad_u, grad_h = sp.symbols("grad_u grad_h", real=True)
STATIC_DENSITY = sp.Rational(1, 2) * b * grad_u**2 + c * grad_u * grad_h + sp.Rational(1, 2) * k * grad_h**2
ACTION_HESSIAN = sp.hessian(STATIC_DENSITY, (grad_u, grad_h))
b_action, c_action, k_action = ACTION_HESSIAN[0, 0], ACTION_HESSIAN[0, 1], ACTION_HESSIAN[1, 1]
d = sp.factor(ACTION_HESSIAN.det())
kappa = sp.factor(d / b_action)

# The maximum-principle/Schur escape factor zeta is strictly positive for the
# positive Robin operator.  Defining ree_full=(1-zeta)/km retains the exact
# identity z_g=1-km<eta,L_h^-1 eta>=zeta without treating a bare source as a
# dressed monopole.
ree_full = (1 - zeta) / km
zg = sp.factor(1 - km * ree_full)
zb = 1 - km * reb

# Completing the boundary shear after w=u+(c/b)h gives this exact map from
# physical (u-conormal, h-source) data to the two diagonal bulk monopoles.
vertex = sp.Matrix([[1, 0], [-(c_action / b_action) * zb, zg]])
metric = sp.diag(1 / b_action, 1 / kappa)
mfar = sp.simplify(vertex.T * metric * vertex)
muu, mug, mgg = map(sp.factor, (mfar[0, 0], mfar[0, 1], mfar[1, 1]))
mdet = sp.factor(mfar.det())
qphysical_v = sp.factor((u - sug * gh) / suu)
qdatum_v = sp.factor(u / suu)
vhh = sp.factor(mgg - 2 * (sug / suu) * mug + (sug / suu) ** 2 * muu)


def split_readouts(total: sp.Expr, datum: sp.Symbol) -> OrderedDict[str, sp.Expr]:
    total = sp.factor(total)
    h_only = sp.factor(total.subs(datum, 0))
    u_only = sp.factor(total.subs(gh, 0))
    interference = sp.factor(total - h_only - u_only)
    return OrderedDict(
        h_only=h_only, uL_only=u_only, interference=interference, total=total
    )


@lru_cache(maxsize=None)
def fixed_value_coefficient(member: str) -> sp.Expr:
    """Solve the reacting Q values and evaluate the selected functional."""
    eps = sp.symbols("eps")
    self_matrix = sp.Matrix([[suu, sug], [sug, sgg]])
    qbar = (u - sug * gh) / suu
    dq1 = -eps * stwo * (muu * qbar + mug * gh) / suu
    dq2 = -eps * sone * (muu * qbar + mug * gh) / suu
    x1 = sp.Matrix([sone * qbar + dq1, sone * gh])
    x2 = sp.Matrix([stwo * qbar + dq2, stwo * gh])
    sigma = {"conjugate": -1, "bare": 1}[member]
    functional = (
        sp.Rational(1, 2) * sigma
        * ((x1.T * self_matrix * x1)[0] + (x2.T * self_matrix * x2)[0])
        + sigma * eps * (x1.T * mfar * x2)[0]
    )
    pair = sp.expand(sp.series(functional, eps, 0, 2).removeO()).coeff(eps)
    return sp.factor(pair / (sone * stwo))


@lru_cache(maxsize=None)
def fixed_flux_or_source_coefficient(ensemble: str, member: str) -> sp.Expr:
    datum = qf if ensemble == "M" else js
    x = sp.Matrix([datum, gh])
    bare = sp.expand((x.T * mfar * x)[0])
    u_work = sp.expand(2 * datum * (mfar * x)[0])
    h_work = sp.expand(2 * gh * (mfar * x)[1])
    if member == "bare":
        return sp.factor(bare)
    if ensemble == "M":       # q held; only the h source is Legendre-dualized
        return sp.factor(bare - h_work)
    if ensemble == "J":       # both J and g are held sources
        return sp.factor(bare - u_work - h_work)
    raise KeyError(ensemble)


@lru_cache(maxsize=None)
def derive_readouts(ensemble: str, member: str = "conjugate") -> OrderedDict[str, sp.Expr]:
    if ensemble == "V":
        total = fixed_value_coefficient("conjugate" if member == "conjugate" else "bare")
        datum = u
    else:
        total = fixed_flux_or_source_coefficient(ensemble, member)
        datum = qf if ensemble == "M" else js
    return split_readouts(total, datum)


READOUTS = OrderedDict((name, derive_readouts(name)) for name in ("V", "M", "J"))

# Exact identities are assertions, not definitions of the readouts.
EXPECTED_IDENTITIES = {
    "V": muu * u**2 / suu**2 - vhh * gh**2,
    "M": muu * qf**2 - mgg * gh**2,
    "J": -(muu * js**2 + 2 * mug * js * gh + mgg * gh**2),
}
for _ens in READOUTS:
    if not zq(READOUTS[_ens]["total"] - EXPECTED_IDENTITIES[_ens]):
        raise AssertionError(f"thermodynamic derivation mismatch: {_ens}")


def positive_certificate(expr: sp.Expr) -> bool:
    expr = sp.factor(expr)
    certificates = (muu, mgg, vhh, mdet, muu / suu**2)
    stable_card = bool(d.subs({b: 2, c: sp.Rational(1, 2), k: 1}).is_positive)
    robin_nondegenerate = bool(zg.is_positive)
    return stable_card and robin_nondegenerate and any(zq(expr - candidate) for candidate in certificates)


def quadratic_signature(expr: sp.Expr, datum: sp.Symbol) -> str:
    poly = sp.Poly(sp.expand(expr), datum, gh)
    aa = sp.factor(poly.coeff_monomial(datum**2))
    bb = sp.factor(poly.coeff_monomial(datum * gh))
    cc = sp.factor(poly.coeff_monomial(gh**2))
    determinant = sp.factor(aa * cc - bb**2 / 4)
    if positive_certificate(aa) and positive_certificate(-cc) and determinant != 0:
        return "INDEFINITE"
    if positive_certificate(-aa) and positive_certificate(determinant):
        return "NEGATIVE_DEFINITE"
    if zq(expr):
        return "ZERO"
    return "UNCLASSIFIED"


DATUMS = {"V": u, "M": qf, "J": js}
NEUTRAL_SIGNATURES = OrderedDict(
    (ens, quadratic_signature(READOUTS[ens]["total"], DATUMS[ens]))
    for ens in READOUTS
)
if tuple(NEUTRAL_SIGNATURES.values()) != ("INDEFINITE", "INDEFINITE", "NEGATIVE_DEFINITE"):
    raise AssertionError(f"neutral sign classification failed: {NEUTRAL_SIGNATURES}")


# ---------------------------------------------------------------------------
# Q1 barrier evaluation: every candidate goes through the same stationary-
# point/curvature/protection evaluator.  The injection adds an action term.
# ---------------------------------------------------------------------------

xmouth, lam, mu = sp.symbols("xmouth lam mu", real=True)
sorient = sp.symbols("sorient", real=True, nonzero=True)
uheld = sp.symbols("uheld", positive=True)


@dataclass(frozen=True)
class BarrierCandidate:
    mechanism: str
    potential: sp.Expr
    constraints: tuple[sp.Expr, ...]
    constraint_is_protected: bool
    topology_components: int
    reason: str


@dataclass(frozen=True)
class BarrierResult:
    mechanism: str
    direct_datum: bool
    stationary_at_signed_datum: bool
    positive_barrier_curvature: bool
    protected: bool
    pins: bool
    reason: str


def barrier_eval(candidate: BarrierCandidate) -> BarrierResult:
    target = sorient * uheld
    derivative = sp.diff(candidate.potential, xmouth)
    curvature = sp.diff(candidate.potential, xmouth, 2)
    target_in_data = target.free_symbols <= candidate.potential.free_symbols or any(
        target.free_symbols <= constraint.free_symbols for constraint in candidate.constraints
    )
    branches = (-1, 1)
    constraint_holds_target = bool(candidate.constraints) and all(
        all(zq(constraint.subs({sorient: branch, xmouth: branch * uheld}))
            for constraint in candidate.constraints)
        for branch in branches
    )
    stationary = constraint_holds_target or all(
        zq(derivative.subs({sorient: branch, xmouth: branch * uheld}))
        for branch in branches
    )
    positive_curvature = all(
        bool(sp.simplify(curvature.subs({lam: 1, mu: sp.Rational(1, 3), sorient: branch,
                                         uheld: 1, xmouth: branch})).is_positive)
        for branch in branches
    )
    probe = {lam: 1, mu: sp.Rational(1, 3), sorient: 1, uheld: 1}
    stationary_points = sp.solve(derivative.subs(probe), xmouth)
    local_minima = [root for root in stationary_points
                    if bool(sp.simplify(curvature.subs(probe).subs(xmouth, root)).is_positive)]
    disconnected_minima = len(local_minima) > 1 and candidate.topology_components > 1
    protected = (constraint_holds_target and candidate.constraint_is_protected) or disconnected_minima
    pins = target_in_data and stationary and positive_curvature and protected
    return BarrierResult(
        candidate.mechanism, target_in_data, stationary, positive_curvature,
        protected, pins, candidate.reason,
    )


def q1_census(inject: bool = False) -> tuple[str, list[BarrierResult]]:
    src = read_sources()
    ledger, _ = parsed_g0_ledger(src["g0"])
    ql, induced = sp.symbols("qL induced", positive=True)
    b0, c0, k0 = sp.Rational(2), sp.Rational(1, 2), sp.Rational(1)
    relaxed = sp.factor(b0 - c0**2 / k0)
    self_energy = sp.Rational(1, 2) * b0 * xmouth**2
    mixed_energy = sp.Rational(1, 2) * relaxed * xmouth**2
    induced_energy = mixed_energy - sorient * induced * xmouth
    linear_source = mixed_energy - sorient * ql * xmouth
    connected_components = 1 if "DeltaE_unwind = 0" in src["path24"] else 2
    candidates = [
        BarrierCandidate("u_L self-stiffness", self_energy, (), False, 1,
                         "positive bulk Hessian and IR zero give relaxation, not a mouth barrier"),
        BarrierCandidate("C_hu mixing", mixed_energy, (), False, 1,
                         "positive Schur stiffness is convex and supplies no essential datum"),
        BarrierCandidate("h mouth source + natural u_L conormal", induced_energy, (), False, 1,
                         "mixing induces a unique relaxable response, not a held u_L datum"),
        BarrierCandidate("IR u_L->0", mixed_energy, (xmouth,), False, 1,
                         "the IR condition fixes zero, inconsistent with a signed nonzero mouth target"),
        BarrierCandidate("conditional q_L proportional to s", linear_source, (), False, 1,
                         "a linear J-type source shifts one convex minimum but supplies no protection"),
        BarrierCandidate("wall/chi_B sector", mixed_energy, (), False, 1,
                         "the parsed POSTULATE ledger contains r_B-u_L zeros"),
        BarrierCandidate("sleeve +/-w orientation", mixed_energy, (), False, connected_components,
                         "pathA_24 has connected orientation space and zero unwinding barrier"),
        BarrierCandidate("S_hold", mixed_energy, (), False, 1,
                         "S_hold fixes r_B geometry and the parsed ledger supplies no u_L coupling"),
        BarrierCandidate("geon", mixed_energy, (), False, 1,
                         "the parsed POSTULATE ledger sets the u_L geon derivative to zero"),
    ]
    if inject:
        selector = 3 * xmouth / (2 * uheld) - xmouth**3 / (2 * uheld**3)
        injected_double_well = (lam / 4) * (xmouth**2 - uheld**2) ** 2 - mu * sorient * selector
        candidates.append(BarrierCandidate(
            "injected_s_uL_hold", injected_double_well, (), False, 2,
            "0<mu<2*lam*uheld^4/3 preserves two wells while an s-odd selector picks x=s*uheld",
        ))
    results = [barrier_eval(row) for row in candidates]
    holders = [row.mechanism for row in results if row.pins]
    return (f"NATIVE_CLAMP_EXISTS({holders[0]})" if holders else "NO_NATIVE_CLAMP"), results


Q1_STATUS, Q1_ROWS = q1_census()
Q1B_STATUS = "BOLT_ON_DEFERRED_TO_R1"
Q1B_REQUIREMENT = (
    "protected/disconnected +/-w sector plus a direct s<->u_L datum coupling; "
    "supply an exact added action and coefficient domain"
)


# ---------------------------------------------------------------------------
# Q3, 3D falloff, dimensions, and the total §5 classifier.
# ---------------------------------------------------------------------------

CLOSURE_CONTENT = {
    "V": {"essential_uL_value", "charge_odd"},
    "M": {"held_uL_conormal", "charge_odd"},
    "J": {"linear_uL_source", "charge_odd"},
}
FORBIDDEN_CLOSURE_CONTENT = {"u_T_datum", "zero_ledger_coupling", "q_M"}


def tier_a(c_value: sp.Expr = sp.Rational(1, 2), additions: Iterable[str] = ()) -> OrderedDict[str, bool]:
    introduced = set().union(*CLOSURE_CONTENT.values()) | set(additions)
    margin = sp.factor(sp.Rational(2) * sp.Rational(1) - c_value**2)
    u2_text = U2_REPORT.read_text(encoding="utf-8")
    ledger, tags = parsed_g0_ledger(G0.read_text(encoding="utf-8"))
    return OrderedDict(
        scalar_hessian_positive=bool(margin.is_positive),
        transverse_decoupled=("u_T_datum" not in introduced and "u_Lu_T" in " ".join(ledger)),
        zero_ledger_preserved=(tags and len(ledger) == 13 and
                               not (introduced & {"zero_ledger_coupling"})),
        charge_perp_mass=("q_M" not in introduced),
        u2_nonselection=("{'UNRESOLVED': 144}" in u2_text),
    )


def derive_green_power(spatial_dimension: int) -> int:
    power = sp.symbols("green_power")
    radial_harmonic_coefficient = sp.factor(power * (power - spatial_dimension + 2))
    roots = sp.solve(radial_harmonic_coefficient, power)
    decaying = [root for root in roots if root.is_real and root.is_positive]
    if len(decaying) != 1:
        raise AssertionError(f"no unique decaying radial Green branch: {roots}")
    return int(decaying[0])


def interaction_energy(coefficient: sp.Expr, spatial_dimension: int = 3) -> sp.Expr:
    return sp.factor(sone * stwo * coefficient / (4 * sp.pi * rr ** derive_green_power(spatial_dimension)))


def force_and_power(coefficient: sp.Expr, spatial_dimension: int = 3) -> tuple[sp.Expr, int]:
    force = sp.factor(-sp.diff(interaction_energy(coefficient, spatial_dimension), rr))
    exponent = -int(force.as_powers_dict()[rr])
    return force, exponent


Dim = tuple[int, int, int]
ZERO_DIM: Dim = (0, 0, 0)
DIM_L: Dim = (1, 0, 0)
DIM_E: Dim = (2, -2, 1)
DIM_F: Dim = (1, -2, 1)


def dadd(*items: Dim) -> Dim:
    return tuple(sum(column) for column in zip(*items))  # type: ignore[return-value]


def dscale(n: int, item: Dim) -> Dim:
    return tuple(n * value for value in item)  # type: ignore[return-value]


SYMBOL_DIMS: dict[sp.Symbol, Dim] = {
    b: (-1, -2, 1), c: (0, -2, 1), k: (1, -2, 1),
    km: DIM_E, zeta: ZERO_DIM, reb: (-2, 2, -1),
    suu: (0, 2, -1), sug: (-1, 2, -1), sgg: (-2, 2, -1),
    u: DIM_L, qf: (1, -2, 1), js: (1, -2, 1), gh: DIM_E,
    sone: ZERO_DIM, stwo: ZERO_DIM, rr: DIM_L,
}


def infer_dim(expr: sp.Expr) -> Dim:
    expr = sp.sympify(expr)
    if expr.is_Number:
        return ZERO_DIM
    if expr.is_Symbol:
        return SYMBOL_DIMS[expr]
    if expr.is_Add:
        dims = [infer_dim(term) for term in expr.args if term != 0]
        if not dims or any(item != dims[0] for item in dims[1:]):
            raise ValueError(f"inhomogeneous sum: {expr} -> {dims}")
        return dims[0]
    if expr.is_Mul:
        return dadd(*(infer_dim(term) for term in expr.args))
    if expr.is_Pow and expr.exp.is_Integer:
        return dscale(int(expr.exp), infer_dim(expr.base))
    raise ValueError(f"unsupported dimension node: {expr}")


def dimensional_checks(green_dim: Dim = (-1, 0, 0)) -> bool:
    coefficient_terms: list[sp.Expr] = []
    for readouts in READOUTS.values():
        for name in ("h_only", "uL_only", "interference", "total"):
            expanded = sp.expand(readouts[name])
            coefficient_terms.extend(term for term in sp.Add.make_args(expanded) if term != 0)
    try:
        coefficient_dims = [infer_dim(term) for term in coefficient_terms]
    except (KeyError, ValueError):
        return False
    target_a = dadd(DIM_E, DIM_L)
    return (
        all(item == target_a for item in coefficient_dims)
        and all(dadd(item, green_dim) == DIM_E for item in coefficient_dims)
        and all(dadd(item, dscale(-2, DIM_L)) == DIM_F for item in coefficient_dims)
        and infer_dim(zg) == ZERO_DIM and infer_dim(zb) == ZERO_DIM
    )


TIER_B = (
    "assembled gravity/drain response under datum", "curved-sleeve stability",
    "momentum/return closure", "dressed non-perturbative two-body response",
)


def section5_algebraic(signatures: Mapping[str, str], force_power: int) -> dict[str, str]:
    """The first function allowed to translate neutral signs to the EM target."""
    answer: dict[str, str] = {}
    for ensemble, signature in signatures.items():
        if signature == "INDEFINITE":
            answer[ensemble] = f"range(attract_1/R^{force_power}|null_leading|repel_1/R^{force_power})"
        elif signature == "NEGATIVE_DEFINITE":
            answer[ensemble] = f"attract_1/R^{force_power}"
        elif signature == "ZERO":
            answer[ensemble] = "null_leading"
        else:
            answer[ensemble] = f"unresolved({signature})"
    return answer


def section5_verdict(*, tier_a_pass: bool, tier_a_reason: str = "tier_A",
                     native_holder: bool = False, bolt_holder: bool = False,
                     bolt_supplied: bool = False,
                     bolt_ineffective_reason: str = "does_not_hold_datum",
                     actual_force_sign: str = "POSITIVE", actual_force_power: int = 2,
                     tier_b: str = "DEFERRED", unresolved: str | None = None,
                     algebraic: Mapping[str, str] | None = None,
                     precedence_mutation: bool = False) -> str:
    """Directive §5 first-match classifier; target facts exist only here."""
    target_force_sign, target_force_power = "POSITIVE", 2

    def holder_landing() -> str | None:
        holder = "NATIVE" if native_holder else ("BOLT_ON" if bolt_holder else None)
        if holder is None:
            return None
        if actual_force_sign == "NULL":
            return "NO_GO(holder_held_variable_null)"
        if actual_force_power != target_force_power:
            return "NO_GO(holder_held_variable_wrong_range)"
        if actual_force_sign != target_force_sign:
            return "NO_GO(holder_held_variable_wrong_sign)"
        if tier_b == "PASS":
            return f"{holder}_DERIVED"
        return f"{holder}_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED"

    # The ablation swaps the first two precedence rules in the production path.
    if precedence_mutation:
        early = holder_landing()
        if early is not None:
            return early
    if not tier_a_pass:
        return f"NO_GO({tier_a_reason})"
    if not precedence_mutation:
        held = holder_landing()
        if held is not None:
            return held
    if unresolved is not None:
        return f"UNRESOLVED({unresolved})"
    alg = algebraic or {name: "UNCLASSIFIED" for name in ("V", "M", "J")}
    result = "NO_NATIVE_CLAMP; ALGEBRAIC_SIGN={" + ",".join(
        f"{name}:{alg[name]}" for name in ("V", "M", "J")
    ) + "}; BOLT_ON_DEFERRED_TO_R1"
    if bolt_supplied and not bolt_holder:
        result += f"; BOLT_ON_INEFFECTIVE({bolt_ineffective_reason})"
    return result


def verdict_oracle(case: Mapping[str, Any]) -> str:
    """Independent declarative decision-table oracle for precedence testing."""
    if not case["tier_a_pass"]:
        return f"NO_GO({case['tier_a_reason']})"
    holder = "NATIVE" if case["native_holder"] else ("BOLT_ON" if case["bolt_holder"] else None)
    if holder:
        outcome = {
            "NULL": "NO_GO(holder_held_variable_null)",
            "WRONG_RANGE": "NO_GO(holder_held_variable_wrong_range)",
            "NEGATIVE": "NO_GO(holder_held_variable_wrong_sign)",
            "PASS": f"{holder}_DERIVED" if case["tier_b"] == "PASS"
                    else f"{holder}_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED",
        }
        key = "NULL" if case["sign"] == "NULL" else (
            "WRONG_RANGE" if case["power"] != 2 else (
                "NEGATIVE" if case["sign"] != "POSITIVE" else "PASS"
            )
        )
        return outcome[key]
    if case["unresolved"]:
        return "UNRESOLVED(named_datum)"
    suffix = "; BOLT_ON_INEFFECTIVE(does_not_hold_datum)" if case["bolt_supplied"] else ""
    return "NO_NATIVE_CLAMP; ALGEBRAIC_SIGN={V:UNCLASSIFIED,M:UNCLASSIFIED,J:UNCLASSIFIED}; BOLT_ON_DEFERRED_TO_R1" + suffix


def verdict_truth_table(precedence_mutation: bool = False) -> tuple[bool, set[str], int]:
    classes: set[str] = set()
    total = 0
    exact = True
    for values in itertools.product(
        (False, True), (False, True), (False, True), (False, True),
        ("POSITIVE", "NEGATIVE", "NULL"), (2, 3), ("PASS", "DEFERRED"), (False, True),
    ):
        tier_ok, native, bolt, bolt_supplied, sign, power, tier_b, unresolved = values
        case = dict(tier_a_pass=tier_ok, tier_a_reason="tier_A", native_holder=native,
                    bolt_holder=bolt, bolt_supplied=bolt_supplied, sign=sign, power=power,
                    tier_b=tier_b, unresolved=unresolved)
        actual = section5_verdict(
            tier_a_pass=tier_ok, native_holder=native, bolt_holder=bolt,
            bolt_supplied=bolt_supplied, actual_force_sign=sign,
            actual_force_power=power, tier_b=tier_b,
            unresolved="named_datum" if unresolved else None,
            precedence_mutation=precedence_mutation,
        )
        expected = verdict_oracle(case)
        exact &= actual == expected and isinstance(actual, str) and bool(actual)
        classes.add(actual.split(";", 1)[0])
        total += 1
    required = {
        "NO_GO(tier_A)", "NO_GO(holder_held_variable_wrong_sign)",
        "NO_GO(holder_held_variable_wrong_range)", "NO_GO(holder_held_variable_null)",
        "NATIVE_DERIVED", "BOLT_ON_DERIVED",
        "NATIVE_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED",
        "BOLT_ON_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED",
        "NO_NATIVE_CLAMP", "UNRESOLVED(named_datum)",
    }
    # BOLT_ON_INEFFECTIVE is a suffix class, checked separately.
    suffix_seen = any(
        "BOLT_ON_INEFFECTIVE" in section5_verdict(tier_a_pass=True, bolt_supplied=True)
        for _ in (0,)
    )
    return exact and required <= classes and suffix_seen, classes | ({"BOLT_ON_INEFFECTIVE"} if suffix_seen else set()), total


FORCE_EXPR, FORCE_POWER = force_and_power(READOUTS["J"]["total"])
MAIN_ALGEBRAIC = section5_algebraic(NEUTRAL_SIGNATURES, FORCE_POWER)
MAIN_LANDING = section5_verdict(tier_a_pass=all(tier_a().values()), algebraic=MAIN_ALGEBRAIC)


# ---------------------------------------------------------------------------
# Semantic target-blind guard and symbolic cross-engine protocol.
# ---------------------------------------------------------------------------

def walk_expressions(value: Any) -> Iterable[sp.Expr]:
    if isinstance(value, sp.Expr):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from walk_expressions(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from walk_expressions(item)


def ast_target_scope_ok() -> bool:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    allowed = {"section5_algebraic", "section5_verdict", "verdict_oracle", "verdict_truth_table",
               "ast_target_scope_ok"}
    bad_tokens = ("repel", "attract", "coulomb", "target_force")

    class Visitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.stack: list[str] = []
            self.bad: list[str] = []

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
            self.stack.append(node.name)
            self.generic_visit(node)
            self.stack.pop()

        def visit_Constant(self, node: ast.Constant) -> None:
            if isinstance(node.value, str) and any(token in node.value.lower() for token in bad_tokens):
                if not self.stack or self.stack[-1] not in allowed:
                    self.bad.append(node.value)

    visitor = Visitor()
    visitor.visit(tree)
    return not visitor.bad


def target_blind(stage_objects: Any) -> bool:
    forbidden = {"em_target_sign", "em_target_power", "expected_formula"}
    symbolic_clean = all(
        not ({str(symbol) for symbol in expr.free_symbols} & forbidden)
        for expr in walk_expressions(stage_objects)
    )
    neutral_clean = all(value in {"INDEFINITE", "NEGATIVE_DEFINITE", "ZERO", "UNCLASSIFIED"}
                        for value in NEUTRAL_SIGNATURES.values())
    return symbolic_clean and neutral_clean and ast_target_scope_ok()


def canonical(expr: sp.Expr) -> str:
    return sp.sstr(sp.cancel(expr)).replace("**", "^")


def symbolic_payload(payload_mutation: bool = False) -> OrderedDict[str, Any]:
    terms: OrderedDict[str, str] = OrderedDict(
        d=canonical(d), kappa=canonical(kappa), zg=canonical(zg), zb=canonical(zb),
        muu=canonical(muu), mug=canonical(mug), mgg=canonical(mgg), mdet=canonical(mdet),
        vhh=canonical(vhh), qphysical_v=canonical(qphysical_v),
    )
    for ensemble, readouts in READOUTS.items():
        for name, expression in readouts.items():
            terms[f"{ensemble}.{name}"] = canonical(expression)
    if payload_mutation:
        terms["J.interference"] = canonical(READOUTS["J"]["interference"] + js * gh)
    truth, classes, total = verdict_truth_table()
    return OrderedDict(
        schema="LEFTOVER_SCALAR_SIGN_V2_SYMBOLIC",
        symbolic_terms=terms,
        neutral_signatures=dict(NEUTRAL_SIGNATURES),
        green_power=derive_green_power(3), force_power=FORCE_POWER,
        q1=Q1_STATUS, q1b=Q1B_STATUS,
        tier_a="PASS" if all(tier_a().values()) else "FAIL",
        tier_a_details=dict(tier_a()),
        tier_b="DEFERRED_R1_R3", truth_table_total=total,
        truth_table_classes=len(classes), truth_table_pass=truth,
        landing=MAIN_LANDING,
    )


PARSE_LOCALS = {str(symbol): symbol for symbol in (b, c, k, km, zeta, reb, suu, sug, sgg, u, qf, js, gh)}


def parse_other(expression: str) -> sp.Expr:
    return sp.sympify(expression.replace("^", "**"), locals=PARSE_LOCALS)


def payloads_equal(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    if set(left) != set(right) or set(left["symbolic_terms"]) != set(right["symbolic_terms"]):
        return False
    for key, expression in left["symbolic_terms"].items():
        if not zq(parse_other(expression) - parse_other(right["symbolic_terms"][key])):
            return False
    for key in left:
        if key != "symbolic_terms" and left[key] != right[key]:
            return False
    return True


def wolfram_payload(mutate: bool = False) -> Mapping[str, Any]:
    env = {**os.environ, "LEFTOVER_JSON_ONLY": "1"}
    if mutate:
        env["LEFTOVER_PAYLOAD_MUTATION"] = "1"
    run = subprocess.run(
        ["math", "-script", str(WL_PATH)], cwd=ROOT, env=env,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=180, check=False,
    )
    if run.returncode != 0:
        raise AssertionError(f"Mathematica payload exit {run.returncode}:\n{run.stdout}")
    marker = next((line for line in run.stdout.splitlines() if line.startswith("JSON_PAYLOAD=")), None)
    if marker is None:
        raise AssertionError(f"Mathematica payload marker missing:\n{run.stdout}")
    return json.loads(marker.split("=", 1)[1])


# ---------------------------------------------------------------------------
# Every §6 tooth and its own production-path mutation.
# ---------------------------------------------------------------------------

TOOTH_ORDER = (
    "FAITHFULNESS", "Q1_BARRIER", "HELD_VARIABLE_V", "HELD_VARIABLE_M",
    "HELD_VARIABLE_J", "DOUBLE_COUNT", "C_HU_STABILITY", "FALLOFF",
    "UNITS", "Q_M_GUARD", "VERDICT_CLASSIFICATION", "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
)


def local_checks(mutation: str | None, dual_ok: bool) -> OrderedDict[str, bool]:
    transcript = TRANSCRIPT
    checked_hessian = ACTION_HESSIAN
    if mutation == "FAITHFULNESS":
        changed = list(transcript.action_terms)
        changed[4] = "-2*C_hu*grad u_L.grad h"
        transcript = replace(transcript, action_terms=tuple(changed))
        mutated_density = (sp.Rational(1, 2) * b * grad_u**2 + 2 * c * grad_u * grad_h
                           + sp.Rational(1, 2) * k * grad_h**2)
        checked_hessian = sp.hessian(mutated_density, (grad_u, grad_h))
    faithful, _ = source_faithfulness(transcript, checked_hessian)
    faithful &= (
        not zq(sp.diff(mfar[0, 1], c))
        and not zq(sp.diff(mfar[0, 1], km))
        and all(zq((mfar.subs({km: 0, zeta: 1}) - ACTION_HESSIAN.inv())[i, j])
                for i in range(2) for j in range(2))
        and len(TRANSCRIPT.zero_ledger) == 13
    )

    q1_status, _ = q1_census(inject=(mutation == "Q1_BARRIER"))

    selected: dict[str, OrderedDict[str, sp.Expr]] = {}
    for ensemble in READOUTS:
        member = "bare" if mutation == f"HELD_VARIABLE_{ensemble}" else "conjugate"
        selected[ensemble] = derive_readouts(ensemble, member)

    reconstructed = {
        ensemble: sp.factor(values["h_only"] + values["uL_only"] + values["interference"])
        for ensemble, values in READOUTS.items()
    }
    if mutation == "DOUBLE_COUNT":
        reconstructed["J"] = sp.factor(READOUTS["J"]["h_only"] + READOUTS["J"]["uL_only"])

    c_value = sp.Integer(2) if mutation == "C_HU_STABILITY" else sp.Rational(1, 2)
    tier = tier_a(c_value)
    stability_landing = section5_verdict(
        tier_a_pass=all(tier.values()), tier_a_reason="scalar_unstable"
    )

    falloff_dimension = 4 if mutation == "FALLOFF" else 3
    _, inferred_power = force_and_power(READOUTS["J"]["total"], falloff_dimension)
    unit_green = (-2, 0, 0) if mutation == "UNITS" else (-1, 0, 0)
    additions = ("q_M",) if mutation == "Q_M_GUARD" else ()
    verdict_ok, _, _ = verdict_truth_table(precedence_mutation=(mutation == "VERDICT_CLASSIFICATION"))

    stage_readouts: Any = READOUTS
    if mutation == "TARGET_BLINDNESS":
        target_symbol = sp.Symbol("em_target_sign")
        stage_readouts = {key: dict(value) for key, value in READOUTS.items()}
        stage_readouts["J"]["total"] *= target_symbol
    stage_objects = (d, kappa, mfar, stage_readouts, tuple(NEUTRAL_SIGNATURES.values()), Q1_STATUS, tier_a())

    result = OrderedDict(
        FAITHFULNESS=faithful,
        Q1_BARRIER=(q1_status == "NO_NATIVE_CLAMP"),
        HELD_VARIABLE_V=(zq(selected["V"]["total"] - READOUTS["V"]["total"])
                         and not zq(derive_readouts("V", "bare")["total"] - READOUTS["V"]["total"])),
        HELD_VARIABLE_M=(zq(selected["M"]["total"] - READOUTS["M"]["total"])
                         and not zq(derive_readouts("M", "bare")["total"] - READOUTS["M"]["total"])),
        HELD_VARIABLE_J=(zq(selected["J"]["total"] - READOUTS["J"]["total"])
                         and not zq(derive_readouts("J", "bare")["total"] - READOUTS["J"]["total"])),
        DOUBLE_COUNT=all(zq(reconstructed[name] - READOUTS[name]["total"]) for name in READOUTS),
        C_HU_STABILITY=(all(tier.values()) and stability_landing != "NO_GO(scalar_unstable)"),
        FALLOFF=(inferred_power == 2),
        UNITS=dimensional_checks(unit_green),
        Q_M_GUARD=tier_a(additions=additions)["charge_perp_mass"],
        VERDICT_CLASSIFICATION=verdict_ok,
        TARGET_BLINDNESS=target_blind(stage_objects),
        DUAL_ENGINE_TERMS=dual_ok,
    )
    return result


def mutation_campaign() -> OrderedDict[str, str]:
    local_payload = symbolic_payload()
    other = wolfram_payload(False)
    dual_baseline = payloads_equal(local_payload, other)
    baseline = local_checks(None, dual_baseline)
    if not all(baseline.values()):
        _, source_failures = source_faithfulness()
        raise AssertionError(
            f"baseline failures: {[name for name, ok in baseline.items() if not ok]}; "
            f"source_details={source_failures}"
        )

    outcomes: OrderedDict[str, str] = OrderedDict()
    for tooth in TOOTH_ORDER:
        if tooth == "DUAL_ENGINE_TERMS":
            dual_mutated = payloads_equal(local_payload, wolfram_payload(True))
            mutated = local_checks(None, dual_mutated)
        else:
            mutated = local_checks(tooth, dual_baseline)
        failed = [name for name, ok in mutated.items() if not ok]
        if failed != [tooth]:
            raise AssertionError(f"ablation {tooth} failed {failed}, expected [{tooth}]")
        if tooth == "C_HU_STABILITY":
            bad = tier_a(sp.Integer(2))
            landing = section5_verdict(tier_a_pass=all(bad.values()), tier_a_reason="scalar_unstable")
            if landing != "NO_GO(scalar_unstable)":
                raise AssertionError(f"stability mutation did not land exact NO_GO: {landing}")
        outcomes[tooth] = f"FIRED_AT_{tooth}"
    return outcomes


def print_report(ablations: Mapping[str, str]) -> None:
    print("LEFTOVER_SCALAR_ELECTRIC_SIGN — SYMPY")
    print("FAITHFUL_STATIC_FUNCTIONAL=")
    print("  S_Lh=" + " + ".join(TRANSCRIPT.action_terms))
    print("  T_L=B_eff=rho_B0^2/chi_c>0; K_h=M_h*c_E^2; C_hu=1/2; k_m=1")
    print("  mouth=sum_i integral eta_i[1/2*k_m*h^2-g_chih*s_i*h]")
    print("  native BC=" + "; ".join(TRANSCRIPT.native_bcs))
    for name, bc in TRANSCRIPT.test_bcs:
        print(f"  BC_{name}={bc}")
    print(f"  ZERO_LEDGER_POSTULATE_ROWS={len(TRANSCRIPT.zero_ledger)}; exact tags/terms verified")
    print("Q1a=" + Q1_STATUS)
    for row in Q1_ROWS:
        print(f"  {row.mechanism}: direct={row.direct_datum}; stationary={row.stationary_at_signed_datum}; "
              f"curvature={row.positive_barrier_curvature}; protected={row.protected}; pins={row.pins}; {row.reason}")
    print("  injection_control=" + q1_census(True)[0])
    print("Q1b=" + Q1B_STATUS)
    print("  required_added_action_content=" + Q1B_REQUIREMENT)
    print("Q2_RESPONSE=")
    print(f"  D={text_expr(d)}; kappa={text_expr(kappa)}")
    print(f"  z_g={text_expr(zg)}; z_b={text_expr(zb)}; both k_m-live")
    print(f"  m={mfar}; det(m)={text_expr(mdet)}>0 for D>0,z_g!=0")
    print(f"  V physical reaction Q0={text_expr(qphysical_v)}; datum-only normalization qV={text_expr(qdatum_v)}")
    for ensemble, readouts in READOUTS.items():
        print(f"  {ensemble}: quadratic_signature={NEUTRAL_SIGNATURES[ensemble]}")
        for name, coefficient in readouts.items():
            print(f"    {name}: A={text_expr(coefficient)}; U=s1*s2*A/(4*pi*R); "
                  "F_out=s1*s2*A/(4*pi*R^2)")
        rebuilt = sum(readouts[key] for key in ("h_only", "uL_only", "interference"))
        print(f"    double_count={zq(rebuilt-readouts['total'])}")
    ta = tier_a()
    print("Q3_TIER_A=" + ("PASS" if all(ta.values()) else "FAIL"))
    for key, value in ta.items():
        print(f"  {key}={value}")
    print("Q3_TIER_B=DEFERRED")
    for item in TIER_B:
        print("  " + item)
    print("SECTION5_LANDING=" + MAIN_LANDING)
    truth, classes, total = verdict_truth_table()
    print(f"VERDICT_DOMAIN_TOTAL={total}; classes_exercised={len(classes)}; pass={truth}")
    print("TEETH=")
    for tooth in TOOTH_ORDER:
        print(f"  PASS {tooth}; ablation={ablations[tooth]}")
    print("ENGINE_AGREE=PASS")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-only", action="store_true")
    args = parser.parse_args()
    try:
        if args.json_only:
            mutation = os.environ.get("LEFTOVER_PAYLOAD_MUTATION") == "1"
            print("JSON_PAYLOAD=" + json.dumps(symbolic_payload(mutation), separators=(",", ":")))
            return 0
        ablations = mutation_campaign()
        print_report(ablations)
        return 0
    except Exception as exc:
        print(f"FIRST_FAILURE={type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
