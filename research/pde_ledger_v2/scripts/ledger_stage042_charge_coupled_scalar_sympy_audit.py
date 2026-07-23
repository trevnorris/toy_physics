#!/usr/bin/env python3
"""Ledger stage042 SymPy audit: charge-coupled scalar consistency map.

Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.  This
engine computes the static exchange, radiation scaling, five-channel
adjudication, Guard A, controls, deletion sensitivity, and the final
first-match verdict.  It does not resolve the SIM-gated scalar magnitude.

Tooth-local runtime ablation uses ``LEDGER_STAGE042_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction
import numbers
import os
import re
from typing import Any, Iterable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE042_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

TOOTH_ORDER = (
    "STATIC_DETERMINANT",
    "A_QQ_RESIDUAL",
    "A_QM_RESIDUAL",
    "A_MM_RESIDUAL",
    "DECOUPLED_A_QM_NO_H_MASS",
    "MIXED_QL0_EP_RISK_TERM",
    "DECOUPLED_A_MM_DENSITY",
    "QM_FLIP_A_QM_ODD",
    "A_MM_EVEN_SIGNED_FLIP",
    "H_MASS_PROJECTION",
    "RADIATION_BARE_RATIO",
    "RADIATION_BARE_EXPONENT",
    "RADIATION_PINNED_EXPONENT",
    "RADIATION_WRONG_NORM_DISCRIMINATES",
    "STABILITY_BOUND",
    "CHANNEL_H_EP",
    "CHANNEL_RADIATION",
    "CHANNEL_UNIVERSALITY",
    "CHANNEL_UL_EP",
    "CHANNEL_PREFERRED_FRAME",
    "SUBTAGS",
    "EP_ADJUDICATION",
    "PRODUCTION_VERDICT",
    "GUARD_A_PRODUCTION",
    "GUARD_A_NEGATIVE_FIXTURES",
    "GUARD_A_DIRECT_INJECTION",
    "GUARD_A_FORGED_FLAG",
    "GUARD_A_BYPASS_REGRESSIONS",
    "GUARD_A_SCOPE_DENYLIST",
    "CTRL_B",
    "CTRL_C",
    "CTRL_D",
    "CTRL_E",
    "CTRL_F",
    "CTRL_G",
    "CTRL_H",
    "CTRL_I",
    "CTRL_J",
    "CTRL_K_WITHOUT_PINNED_KH",
    "CTRL_K",
    "DELETION_H_EP_PILLAR",
    "DELETION_GATED_INDIVIDUAL_STABLE",
    "DELETION_COLLECTIVE_NATURALLY_HIDDEN",
    "REACHABLE_FALSIFIERS",
    "LORENTZ_NECESSARY_NOT_SUFFICIENT",
    "DIMENSION_HOMOGENEITY",
    "DUAL_ENGINE_TERMS",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "STATIC_DETERMINANT": "add one to the typed determinant target",
    "A_QQ_RESIDUAL": "add one to the typed A_qq closed form",
    "A_QM_RESIDUAL": "add one to the typed A_qm closed form",
    "A_MM_RESIDUAL": "add one to the typed A_mm closed form",
    "DECOUPLED_A_QM_NO_H_MASS": "add one to the decoupled A_qm target",
    "MIXED_QL0_EP_RISK_TERM": "add one to the mixed q_L=0 target",
    "DECOUPLED_A_MM_DENSITY": "add one to the decoupled A_mm target",
    "QM_FLIP_A_QM_ODD": "perturb the q_M sign-flip substitution",
    "A_MM_EVEN_SIGNED_FLIP": "perturb the q_M flip in both A_mm checks",
    "H_MASS_PROJECTION": "corrupt the fixture h-mass projection target",
    "RADIATION_BARE_RATIO": "add one to the typed bare flux-ratio target",
    "RADIATION_BARE_EXPONENT": "insert an extra c_E in the bare ratio",
    "RADIATION_PINNED_EXPONENT": "use M_h=K_h/c_E instead of K_h/c_E^2",
    "RADIATION_WRONG_NORM_DISCRIMINATES": "replace the wrong-normalization branch by the bare branch",
    "STABILITY_BOUND": "corrupt the typed strict-stability slack",
    "CHANNEL_H_EP": "break the private static-Coulomb input",
    "CHANNEL_RADIATION": "select the corrupt-speed branch",
    "CHANNEL_UNIVERSALITY": "make b/ell species indexed",
    "CHANNEL_UL_EP": "escalate the private u_L channel to NO_GO",
    "CHANNEL_PREFERRED_FRAME": "require large-c_E suppression",
    "SUBTAGS": "remove nonzero C_hu from the private subtag ledger",
    "EP_ADJUDICATION": "report unqualified decoupled-floor EP safety",
    "PRODUCTION_VERDICT": "break static-Coulomb match in the private production-verdict witness",
    "GUARD_A_PRODUCTION": "inject a numeric P_h/P_EM into the production tree",
    "GUARD_A_NEGATIVE_FIXTURES": "authorize all five fixtures through a neutralized shared guard",
    "GUARD_A_DIRECT_INJECTION": "authorize the directly injected guarded numbers",
    "GUARD_A_FORGED_FLAG": "let a forged local earned_parent_action flag authorize emission",
    "GUARD_A_BYPASS_REGRESSIONS": "disable tuple descent in the numeric-leaf walk",
    "GUARD_A_SCOPE_DENYLIST": "add diagnostic_gain to the scoped denylist",
    "CTRL_B": "neutralize the h bare-mass-residue fixture",
    "CTRL_C": "neutralize C_hu in the mixed-risk fixture",
    "CTRL_D": "neutralize the corrupt-speed selector",
    "CTRL_E": "neutralize the large-c_E preferred-frame selector",
    "CTRL_F": "neutralize the species-indexed b/ell selector",
    "CTRL_G": "neutralize the strict-stability violation",
    "CTRL_H": "restore the static-Coulomb match",
    "CTRL_I": "neutralize the q_M sign flip",
    "CTRL_J": "neutralize the wrong-normalization selector",
    "CTRL_K_WITHOUT_PINNED_KH": "treat c_E=c_gamma alone as sufficient to escalate radiation",
    "CTRL_K": "drop the pinned-K_h fact",
    "DELETION_H_EP_PILLAR": "leave h_EP unstubbed",
    "DELETION_GATED_INDIVIDUAL_STABLE": "replace one individual earned stub by NO_GO",
    "DELETION_COLLECTIVE_NATURALLY_HIDDEN": "omit universality from the collective earned stub",
    "REACHABLE_FALSIFIERS": "drop control B from the computed active falsifier map",
    "LORENTZ_NECESSARY_NOT_SUFFICIENT": "treat c_E=c_gamma alone as an earned radiation state",
    "DIMENSION_HOMOGENEITY": "change [c_E] from L T^-1 to L^2 T^-1",
    "DUAL_ENGINE_TERMS": "drop the locally computed determinant inventory term",
    "VERDICT_REDERIVATION": "neutralize the computed control-B witness while retaining its named verdict",
    "SOURCE_TO_STAGE_MANIFEST": "mis-scope one preserved source predicate as replaced-by-stronger",
}

if len(TOOTH_ORDER) != 49:
    raise RuntimeError("stage042 tooth declaration is not exactly 49")
if set(ABLATION_DESCRIPTIONS) != set(TOOTH_ORDER):
    raise RuntimeError("stage042 ablation descriptions do not cover the tooth registry")


VERDICT_MAPPED = "SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED"
VERDICT_NO_GO = "NO_GO_CONSISTENCY"
H_EP_SAFE = "EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR"
H_EP_FIFTH = "FIFTH_FORCE_TRIGGERED"
RADIATION_SIM = "SIM_GATED"
RADIATION_TENSION = "FALSIFIABLE_TENSION"
RADIATION_AUDIT = "AUDIT_ONLY_NOT_EARNED"
UNIVERSALITY_SIM = "SIM_GATED_REQUIRED_UNIVERSALITY"
UNIVERSALITY_TENSION = "FALSIFIABLE_TENSION"
UL_SIM = "SIM_GATED"
PF_SIM = "SIM_GATED"
PF_TENSION = "PREFERRED_FRAME_TENSION"
EARNED_SAFE = "EARNED_SAFE"


class AuditFailure(AssertionError):
    """A named stage predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


def expect_bool(name: str, condition: bool, evidence: Any = None) -> None:
    global PASS_COUNT, FAIL_COUNT
    if bool(condition):
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    FAIL_COUNT += 1
    print(f"FIRST_FAILURE={name}")
    if ACTIVE_MUTATION == name:
        print(f"FIRED_AT_OWN_ASSERT={name}")
    print(f"FAIL  {name}: residual = 1")
    if evidence is not None:
        print(f"      evidence = {evidence}")
    raise AuditFailure(name)


def section(text: str) -> None:
    print("")
    print(text)
    print("-" * len(text))


def zero(expr: sp.Expr) -> bool:
    return sp.factor(sp.cancel(expr)) == 0


def speed_exponent(expr: sp.Expr, c_e: sp.Symbol) -> int:
    """SymPy source-route falloff: denominator degree minus numerator degree."""
    z = sp.Symbol("z", positive=True)
    num, den = sp.fraction(sp.together(expr.subs(c_e, z)))
    return int(sp.degree(den, z) - sp.degree(num, z))


@dataclass(frozen=True)
class Core:
    B: sp.Symbol
    K: sp.Symbol
    C: sp.Symbol
    qL: sp.Symbol
    qh: sp.Symbol
    qM: sp.Symbol
    mh: sp.Symbol
    cE: sp.Symbol
    cg: sp.Symbol
    Mh: sp.Symbol
    QE: sp.Symbol
    det: sp.Expr
    Aqq: sp.Expr
    Aqm: sp.Expr
    Amm: sp.Expr
    typed_Aqq: sp.Expr
    typed_Aqm: sp.Expr
    typed_Amm: sp.Expr
    ratio_bare: sp.Expr
    ratio_pinned: sp.Expr
    ratio_corrupt: sp.Expr
    ratio_wrong: sp.Expr
    bare_exp: int
    pinned_exp: int
    corrupt_exp: int
    wrong_exp: int


def derive_core() -> Core:
    B, K, C = sp.symbols("B_eff K_h C_hu", nonzero=True)
    qL, qh, qM, mh = sp.symbols("q_L q_h q_M m_h")
    cE, cg, Mh, QE = sp.symbols("c_E c_gamma M_h Q_E", positive=True)
    omega, d, r, kappa = sp.symbols("omega d r kappa", positive=True)
    stiffness = sp.Matrix([[B, C], [C, K]])
    inv_stiffness = sp.simplify(stiffness.inv())
    det = sp.factor(stiffness.det())
    jq = sp.Matrix([qL, qh])
    jm = sp.Matrix([qM, mh])
    Aqq = sp.factor((jq.T * inv_stiffness * jq)[0])
    Aqm = sp.factor((jq.T * inv_stiffness * jm)[0])
    Amm = sp.factor((jm.T * inv_stiffness * jm)[0])
    typed_Aqq = (B * qh**2 - 2 * C * qL * qh + K * qL**2) / det
    typed_Aqm = (
        B * mh * qh - C * mh * qL - C * qM * qh + K * qL * qM
    ) / det
    typed_Amm = (B * mh**2 - 2 * C * mh * qM + K * qM**2) / det

    accel = omega**2 * d
    scalar_gradient = qh * accel / (Mh * cE**2 * r)
    scalar_power = sp.factor(Mh * cE * scalar_gradient**2 * r**2)
    em_gradient = QE * accel / (cg**2 * r)
    em_power = sp.factor(cg * em_gradient**2 * r**2)
    ratio_bare = sp.factor(scalar_power / em_power)
    ratio_pinned = sp.factor(
        ratio_bare.subs(Mh, K / cE**2) * K * kappa / qh**2
    )
    corrupt_gradient = qh * accel / (Mh * cE * r)
    corrupt_power = sp.factor(Mh * cE * corrupt_gradient**2 * r**2)
    ratio_corrupt = sp.factor(corrupt_power / em_power)
    ratio_wrong = sp.factor(
        ratio_bare.subs(Mh, K / cE**2) * K * kappa / qh**2
    )
    return Core(
        B=B,
        K=K,
        C=C,
        qL=qL,
        qh=qh,
        qM=qM,
        mh=mh,
        cE=cE,
        cg=cg,
        Mh=Mh,
        QE=QE,
        det=det,
        Aqq=Aqq,
        Aqm=Aqm,
        Amm=Amm,
        typed_Aqq=typed_Aqq,
        typed_Aqm=typed_Aqm,
        typed_Amm=typed_Amm,
        ratio_bare=ratio_bare,
        ratio_pinned=ratio_pinned,
        ratio_corrupt=ratio_corrupt,
        ratio_wrong=ratio_wrong,
        bare_exp=speed_exponent(ratio_bare, cE),
        pinned_exp=speed_exponent(ratio_pinned, cE),
        corrupt_exp=speed_exponent(ratio_corrupt, cE),
        wrong_exp=speed_exponent(ratio_wrong, cE),
    )


@dataclass(frozen=True)
class Ledger:
    static_coulomb_match: bool = True
    h_bare_mass_residue_zero: bool = True
    h_bare_mass_residue_fixture: bool = False
    pinned_kh_fact: bool = False
    commit_cE_equals_cgamma: bool = False
    radiation_derivation_ok: bool = True
    radiation_corrupt_speed: bool = False
    radiation_wrong_normalization: bool = False
    q_ratio_global_earned: bool = False
    species_indexed_b_over_ell: bool = False
    uL_no_go: bool = False
    uL_escalated: bool = False
    suppression_requires_large_cE: bool = False
    C_hu_nonzero: bool = True
    qM_nonzero: bool = True
    qM_sign: int = 1
    stability_violation: bool = False
    production_laundering_attempt: bool = False


@dataclass(frozen=True)
class Channels:
    h_EP: str
    radiation: str
    universality: str
    u_L_EP: str
    preferred_frame: str
    cherenkov_deferred: bool
    mixed_scalar_EP_risk: bool
    laundering_refused: bool
    stability_violated: bool


def compute_channels(
    ledger: Ledger,
    core: Core,
    *,
    commit_alone_sufficient: bool = False,
) -> Channels:
    if not ledger.static_coulomb_match:
        h_ep = "NO_GO"
    elif ledger.h_bare_mass_residue_fixture or not ledger.h_bare_mass_residue_zero:
        h_ep = H_EP_FIFTH
    else:
        h_ep = H_EP_SAFE

    if ledger.radiation_corrupt_speed:
        radiation = RADIATION_AUDIT
    elif ledger.pinned_kh_fact and ledger.commit_cE_equals_cgamma:
        radiation = RADIATION_TENSION
    elif commit_alone_sufficient and ledger.commit_cE_equals_cgamma:
        radiation = EARNED_SAFE
    elif ledger.radiation_derivation_ok:
        radiation = RADIATION_SIM
    else:
        radiation = RADIATION_AUDIT

    if ledger.species_indexed_b_over_ell:
        universality = UNIVERSALITY_TENSION
    elif ledger.q_ratio_global_earned:
        universality = EARNED_SAFE
    else:
        universality = UNIVERSALITY_SIM

    if ledger.uL_no_go:
        u_l_ep = "NO_GO"
    elif ledger.uL_escalated:
        u_l_ep = UNIVERSALITY_TENSION
    else:
        u_l_ep = UL_SIM

    preferred_frame = (
        PF_TENSION if ledger.suppression_requires_large_cE else PF_SIM
    )
    mixed_term = core.Aqm.subs({core.qL: 0, core.mh: 0})
    mixed_risk = bool(
        ledger.C_hu_nonzero
        and ledger.qM_nonzero
        and mixed_term != 0
    )
    return Channels(
        h_EP=h_ep,
        radiation=radiation,
        universality=universality,
        u_L_EP=u_l_ep,
        preferred_frame=preferred_frame,
        cherenkov_deferred=radiation == RADIATION_SIM,
        mixed_scalar_EP_risk=mixed_risk,
        laundering_refused=ledger.production_laundering_attempt,
        stability_violated=ledger.stability_violation,
    )


def channel_states(channels: Channels) -> dict[str, str]:
    return {
        "h_EP": channels.h_EP,
        "radiation": channels.radiation,
        "universality": channels.universality,
        "u_L_EP": channels.u_L_EP,
        "preferred_frame": channels.preferred_frame,
    }


def verdict_from_channels(channels: Channels) -> str:
    states = channel_states(channels)
    if (
        "NO_GO" in states.values()
        or channels.laundering_refused
        or channels.stability_violated
    ):
        return VERDICT_NO_GO
    tension_tests = (
        ("h_EP", channels.h_EP == H_EP_FIFTH),
        ("radiation", channels.radiation == RADIATION_TENSION),
        ("universality", channels.universality == UNIVERSALITY_TENSION),
        ("u_L_EP", channels.u_L_EP == UNIVERSALITY_TENSION),
    )
    tension_channels = [name for name, fired in tension_tests if fired]
    if tension_channels:
        return (
            "FALSIFIABLE_TENSION(channel="
            + ",".join(tension_channels)
            + ")"
        )
    gated_or_cost = (
        channels.radiation in {RADIATION_SIM, RADIATION_AUDIT}
        or channels.universality == UNIVERSALITY_SIM
        or channels.u_L_EP == UL_SIM
        or channels.preferred_frame in {PF_SIM, PF_TENSION}
    )
    if channels.h_EP == H_EP_SAFE and gated_or_cost:
        return VERDICT_MAPPED
    if (
        channels.h_EP == H_EP_SAFE
        and all(
            states[name] == EARNED_SAFE
            for name in (
                "radiation",
                "universality",
                "u_L_EP",
                "preferred_frame",
            )
        )
    ):
        return "NATURALLY_HIDDEN"
    return VERDICT_NO_GO


def ep_adjudication(channels: Channels, *, unqualified: bool = False) -> dict[str, Any]:
    full_status = (
        EARNED_SAFE
        if channels.h_EP == H_EP_SAFE
        and channels.universality == EARNED_SAFE
        else "NOT_EARNED"
    )
    return {
        "mass_channel": channels.h_EP,
        "full_decoupled_floor_EP_safety": full_status,
        "unqualified_decoupled_floor_EP_safe_reported": unqualified,
    }


NUMERIC_LITERAL_RE = re.compile(
    r"(?<![A-Za-z0-9_])[-+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][-+]?\d+)?(?![A-Za-z0-9_])"
)
GUARDED_FIELDS = {
    "m_h",
    "c_e",
    "k_h",
    "p_h_over_p_em",
    "ph_over_pem",
    "power_ratio",
    "flux_ratio",
    "ep_magnitude",
    "fifth_force_magnitude",
    "residue_floor",
}


def normalized_token(value: Any) -> str:
    return (
        str(value)
        .strip()
        .lower()
        .replace("/", "_over_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
    )


def guarded_path(path: tuple[Any, ...], extra_fields: frozenset[str]) -> bool:
    tokens = [normalized_token(part) for part in path]
    deny = GUARDED_FIELDS | set(extra_fields)
    return any(
        token in deny
        or any(
            token == prefix or token.startswith(prefix + "_")
            for prefix in ("m_h", "c_e", "k_h")
        )
        for token in tokens
    )


def numeric_leaf(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and (
            isinstance(value, numbers.Number)
            or isinstance(value, sp.Basic) and bool(value.is_number)
        )
    )


def find_guarded_numeric_paths(
    data: Any,
    path: tuple[Any, ...] = (),
    *,
    descend_sequences: bool = True,
    scan_strings: bool = True,
    extra_fields: frozenset[str] = frozenset(),
) -> tuple[str, ...]:
    findings: list[str] = []
    if isinstance(data, dict):
        for key, value in data.items():
            findings.extend(
                find_guarded_numeric_paths(
                    value,
                    path + (key,),
                    descend_sequences=descend_sequences,
                    scan_strings=scan_strings,
                    extra_fields=extra_fields,
                )
            )
    elif isinstance(data, str):
        if (
            scan_strings
            and guarded_path(path, extra_fields)
            and NUMERIC_LITERAL_RE.search(data)
        ):
            findings.append(".".join(str(part) for part in path))
    elif isinstance(data, (list, tuple)):
        if descend_sequences:
            for index, value in enumerate(data):
                findings.extend(
                    find_guarded_numeric_paths(
                        value,
                        path + (index,),
                        descend_sequences=descend_sequences,
                        scan_strings=scan_strings,
                        extra_fields=extra_fields,
                    )
                )
    elif numeric_leaf(data) and guarded_path(path, extra_fields):
        findings.append(".".join(str(part) for part in path))
    return tuple(findings)


def absent_substrate() -> dict[str, Any]:
    return {
        "facts": {"h_time_kinetic_parent_action": "ABSENT"},
        "provenance_objects": {"h_time_kinetic_parent_action": None},
    }


def earned_substrate() -> dict[str, Any]:
    return {
        "facts": {
            "h_time_kinetic_parent_action": "DERIVED_PARENT_ACTION"
        },
        "provenance_objects": {
            "h_time_kinetic_parent_action": {
                "status": "EARNED",
                "source": "substrate:h_time_kinetic_parent_action",
            }
        },
    }


def parent_action_earned(substrate: dict[str, Any]) -> bool:
    fact = substrate["facts"].get("h_time_kinetic_parent_action")
    provenance = substrate["provenance_objects"].get(
        "h_time_kinetic_parent_action"
    )
    return bool(
        fact not in {None, "ABSENT", "MISSING", "NOT_EARNED"}
        and isinstance(provenance, dict)
        and provenance.get("status") == "EARNED"
        and provenance.get("source")
    )


def guard_scan(
    payload: dict[str, Any],
    substrate: dict[str, Any],
    *,
    allow_forged_local: bool = False,
    descend_sequences: bool = True,
    scan_strings: bool = True,
    extra_fields: frozenset[str] = frozenset(),
) -> dict[str, Any]:
    forged = bool(payload.get("earned_parent_action"))
    authorized = parent_action_earned(substrate) or (
        allow_forged_local and forged
    )
    paths = (
        ()
        if authorized
        else find_guarded_numeric_paths(
            payload,
            descend_sequences=descend_sequences,
            scan_strings=scan_strings,
            extra_fields=extra_fields,
        )
    )
    if authorized:
        result = "MAGNITUDE_EMISSION_ALLOWED"
    elif paths:
        result = "LAUNDERING_REFUSED"
    else:
        result = "NO_GUARDED_MAGNITUDE_EMITTED"
    return {"result": result, "paths": paths, "authorized": authorized}


NEGATIVE_FIXTURES = (
    ("M_h_from_N0", {"assembled_results_fixture": {"M_h": 1}}),
    ("M_h_from_K_parallel", {"assembled_results_fixture": {"M_h": 2}}),
    ("c_E_from_c_gamma_Green", {"assembled_results_fixture": {"c_E": 3}}),
    ("K_h_from_N0_cgamma2", {"assembled_results_fixture": {"K_h": 4}}),
    ("residue_floor_emission", {"assembled_results_fixture": {"residue_floor": 5}}),
)


def guard_negative_vector(*, neutralized: bool = False) -> tuple[str, ...]:
    substrate = earned_substrate() if neutralized else absent_substrate()
    return tuple(
        guard_scan(payload, substrate)["result"]
        for _, payload in NEGATIVE_FIXTURES
    )


def stub_channels(channels: Channels, **updates: str) -> Channels:
    return replace(channels, **updates)


def controls(core: Core, *, mutate: str = "") -> dict[str, dict[str, Any]]:
    baseline = verdict_from_channels(compute_channels(Ledger(), core))
    result: dict[str, dict[str, Any]] = {}

    ledger_b = Ledger(
        h_bare_mass_residue_zero=mutate == "CTRL_B",
        h_bare_mass_residue_fixture=mutate != "CTRL_B",
        C_hu_nonzero=False,
    )
    ch_b = compute_channels(ledger_b, core)
    result["B"] = {"channels": ch_b, "verdict": verdict_from_channels(ch_b)}

    ch_c = compute_channels(
        Ledger(C_hu_nonzero=mutate != "CTRL_C", qM_nonzero=True), core
    )
    result["C"] = {"channels": ch_c, "verdict": verdict_from_channels(ch_c)}

    ch_d = compute_channels(
        Ledger(radiation_corrupt_speed=mutate != "CTRL_D"), core
    )
    selected_d_exp = core.corrupt_exp if mutate != "CTRL_D" else core.bare_exp
    result["D"] = {
        "channels": ch_d,
        "verdict": verdict_from_channels(ch_d),
        "baseline_exp": core.bare_exp,
        "selected_exp": selected_d_exp,
    }

    ch_e = compute_channels(
        Ledger(suppression_requires_large_cE=mutate != "CTRL_E"), core
    )
    result["E"] = {"channels": ch_e, "verdict": verdict_from_channels(ch_e)}

    ch_f = compute_channels(
        Ledger(species_indexed_b_over_ell=mutate != "CTRL_F"), core
    )
    result["F"] = {"channels": ch_f, "verdict": verdict_from_channels(ch_f)}

    ch_g = compute_channels(
        Ledger(stability_violation=mutate != "CTRL_G"), core
    )
    result["G"] = {"channels": ch_g, "verdict": verdict_from_channels(ch_g)}

    ch_h = compute_channels(
        Ledger(static_coulomb_match=mutate == "CTRL_H"), core
    )
    result["H"] = {"channels": ch_h, "verdict": verdict_from_channels(ch_h)}

    sign = 1 if mutate == "CTRL_I" else -1
    aqm_base = core.Aqm.subs(core.mh, 0)
    signed_base = core.Amm.subs(core.mh, 0) / core.qM
    aqm_case = aqm_base.subs(core.qM, sign * core.qM)
    signed_case = signed_base.subs(core.qM, sign * core.qM)
    amm_case = core.Amm.subs(
        {core.mh: 0, core.qM: sign * core.qM}
    )
    amm_base = core.Amm.subs(core.mh, 0)
    ch_i = compute_channels(Ledger(qM_sign=sign), core)
    result["I"] = {
        "channels": ch_i,
        "verdict": verdict_from_channels(ch_i),
        "aqm_flipped": zero(aqm_case + aqm_base),
        "signed_flipped": zero(signed_case + signed_base),
        "amm_even": zero(amm_case - amm_base),
    }

    selected_j_exp = core.bare_exp if mutate == "CTRL_J" else core.wrong_exp
    ch_j = compute_channels(
        Ledger(radiation_wrong_normalization=mutate != "CTRL_J"), core
    )
    result["J"] = {
        "channels": ch_j,
        "verdict": verdict_from_channels(ch_j),
        "baseline_exp": core.bare_exp,
        "selected_exp": selected_j_exp,
    }

    ledger_k0 = Ledger(commit_cE_equals_cgamma=True)
    ch_k0 = compute_channels(
        ledger_k0,
        core,
        commit_alone_sufficient=mutate == "CTRL_K_WITHOUT_PINNED_KH",
    )
    result["K_without_pinned_Kh"] = {
        "channels": ch_k0,
        "verdict": verdict_from_channels(ch_k0),
    }

    ledger_k = Ledger(
        pinned_kh_fact=mutate != "CTRL_K",
        commit_cE_equals_cgamma=True,
    )
    ch_k = compute_channels(ledger_k, core)
    result["K"] = {"channels": ch_k, "verdict": verdict_from_channels(ch_k)}
    result["_baseline"] = {"verdict": baseline}
    return result


def deletion_outcomes(
    production: Channels,
    *,
    mutate: str = "",
) -> dict[str, Any]:
    baseline = verdict_from_channels(production)
    h_stub = (
        production
        if mutate == "DELETION_H_EP_PILLAR"
        else stub_channels(production, h_EP=EARNED_SAFE)
    )
    h_verdict = verdict_from_channels(h_stub)

    individual: dict[str, tuple[str, bool]] = {}
    for index, name in enumerate(
        ("radiation", "universality", "u_L_EP", "preferred_frame")
    ):
        value = EARNED_SAFE
        if (
            mutate == "DELETION_GATED_INDIVIDUAL_STABLE"
            and index == 0
        ):
            value = "NO_GO"
        candidate = stub_channels(production, **{name: value})
        candidate_verdict = verdict_from_channels(candidate)
        individual[name] = (
            candidate_verdict,
            candidate_verdict != baseline,
        )

    collective_updates = {
        "radiation": EARNED_SAFE,
        "universality": EARNED_SAFE,
        "u_L_EP": EARNED_SAFE,
        "preferred_frame": EARNED_SAFE,
    }
    if mutate == "DELETION_COLLECTIVE_NATURALLY_HIDDEN":
        collective_updates.pop("universality")
    collective = stub_channels(production, **collective_updates)
    collective_verdict = verdict_from_channels(collective)
    return {
        "h_EP": (h_verdict, h_verdict != baseline),
        "individual": individual,
        "collective": (
            collective_verdict,
            collective_verdict != baseline,
        ),
        "collective_channels": collective,
    }


def reachable_falsifier_map(control_map: dict[str, dict[str, Any]]) -> dict[str, str]:
    baseline = control_map["_baseline"]["verdict"]
    return {
        name: control_map[name]["verdict"]
        for name in (
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K_without_pinned_Kh",
            "K",
        )
        if control_map[name]["verdict"] != baseline
    }


Dim = tuple[Fraction, Fraction, Fraction]


def dadd(*dims: Dim) -> Dim:
    return tuple(sum(items, Fraction(0)) for items in zip(*dims))  # type: ignore[return-value]


def dscale(power: Fraction, dim: Dim) -> Dim:
    return tuple(power * item for item in dim)  # type: ignore[return-value]


def dimension_guard(*, mutate: bool = False) -> dict[str, Any]:
    zero_dim: Dim = (Fraction(0), Fraction(0), Fraction(0))
    stiffness_dim: Dim = (Fraction(1), Fraction(0), Fraction(0))
    speed_dim: Dim = (Fraction(0), Fraction(1), Fraction(-1))
    length_dim: Dim = (Fraction(0), Fraction(1), Fraction(0))
    frequency_dim: Dim = (Fraction(0), Fraction(0), Fraction(-1))
    charge_dim: Dim = (
        Fraction(1, 2),
        Fraction(3, 2),
        Fraction(-1),
    )
    dims: dict[str, Dim] = {
        "B": stiffness_dim,
        "C": stiffness_dim,
        "K": stiffness_dim,
        "qL": charge_dim,
        "qh": charge_dim,
        "qM": charge_dim,
        "mh": charge_dim,
        "cE": (
            Fraction(0),
            Fraction(2 if mutate else 1),
            Fraction(-1),
        ),
        "cg": speed_dim,
        "Mh": zero_dim,
        "QE": charge_dim,
        "omega": frequency_dim,
        "d": length_dim,
        "r": length_dim,
    }
    determinant_terms = (
        dadd(dims["B"], dims["K"]),
        dscale(Fraction(2), dims["C"]),
    )
    det_dim = determinant_terms[0]
    aqq_terms = (
        dadd(dims["B"], dscale(Fraction(2), dims["qh"])),
        dadd(dims["C"], dims["qL"], dims["qh"]),
        dadd(dims["K"], dscale(Fraction(2), dims["qL"])),
    )
    aqm_terms = (
        dadd(dims["B"], dims["mh"], dims["qh"]),
        dadd(dims["C"], dims["mh"], dims["qL"]),
        dadd(dims["C"], dims["qM"], dims["qh"]),
        dadd(dims["K"], dims["qL"], dims["qM"]),
    )
    amm_terms = (
        dadd(dims["B"], dscale(Fraction(2), dims["mh"])),
        dadd(dims["C"], dims["mh"], dims["qM"]),
        dadd(dims["K"], dscale(Fraction(2), dims["qM"])),
    )
    accel_dim = dadd(dscale(Fraction(2), dims["omega"]), dims["d"])
    scalar_gradient_dim = dadd(
        dims["qh"],
        accel_dim,
        dscale(Fraction(-1), dims["Mh"]),
        dscale(Fraction(-2), dims["cE"]),
        dscale(Fraction(-1), dims["r"]),
    )
    scalar_power_dim = dadd(
        dims["Mh"],
        dims["cE"],
        dscale(Fraction(2), scalar_gradient_dim),
        dscale(Fraction(2), dims["r"]),
    )
    em_gradient_dim = dadd(
        dims["QE"],
        accel_dim,
        dscale(Fraction(-2), dims["cg"]),
        dscale(Fraction(-1), dims["r"]),
    )
    em_power_dim = dadd(
        dims["cg"],
        dscale(Fraction(2), em_gradient_dim),
        dscale(Fraction(2), dims["r"]),
    )
    ratio_dim = dadd(
        dscale(Fraction(2), dims["qh"]),
        dscale(Fraction(-1), dims["Mh"]),
        dscale(Fraction(-2), dims["QE"]),
        dscale(Fraction(3), dims["cg"]),
        dscale(Fraction(-3), dims["cE"]),
    )
    physical_power: Dim = (Fraction(1), Fraction(2), Fraction(-3))
    quadratic_target = dadd(
        dscale(Fraction(2), charge_dim),
        dscale(Fraction(-1), stiffness_dim),
    )
    homogeneous = (
        len(set(determinant_terms)) == 1
        and determinant_terms[0] == dscale(Fraction(2), stiffness_dim)
        and len(set(aqq_terms)) == 1
        and len(set(aqm_terms)) == 1
        and len(set(amm_terms)) == 1
        and dadd(aqq_terms[0], dscale(Fraction(-1), det_dim))
        == quadratic_target
        and dadd(aqm_terms[0], dscale(Fraction(-1), det_dim))
        == quadratic_target
        and dadd(amm_terms[0], dscale(Fraction(-1), det_dim))
        == quadratic_target
        and scalar_power_dim == physical_power
        and em_power_dim == physical_power
        and ratio_dim == zero_dim
    )
    return {
        "homogeneous": homogeneous,
        "determinant_terms": determinant_terms,
        "Aqq_terms": aqq_terms,
        "Aqm_terms": aqm_terms,
        "Amm_terms": amm_terms,
        "scalar_power": scalar_power_dim,
        "em_power": em_power_dim,
        "ratio": ratio_dim,
    }


SOURCE_CORE_IDS = (
    "pathA42.py::core.A_qq_residual_zero",
    "pathA42.py::core.A_qm_residual_zero",
    "pathA42.py::core.A_mm_residual_zero",
    "pathA42.py::core.A_qm_decoupled_h_mass_zero",
    "pathA42.py::core.A_qm_mixed_qL0_nonzero_term",
    "pathA42.py::core.A_mm_decoupled_no_h_mass",
    "pathA42.py::core.mass_h_projection_zero",
    "pathA42.py::core.fixture_h_projection_nonzero_symbol",
    "pathA42.py::core.qM_flip_A_qm_residual_zero",
    "pathA42.py::core.A_mm_even_in_qM",
    "pathA42.py::core.A_mm_signed_projection_flips",
    "pathA42.py::core.stability_bound_strict",
    "pathA42.py::core.stability_violation_condition",
    "pathA42.py::core.radiation_bare_flux_ratio_matches",
    "pathA42.py::core.radiation_wrong_normalization_recomputed",
    "pathA42.py::core.radiation_bare_exponent",
    "pathA42.py::core.radiation_pinned_Kh_exponent",
    "pathA42.py::core.radiation_corrupt_speed_exponent",
    "pathA42.py::core.radiation_wrong_normalization_exponent",
)
SOURCE_CHANNEL_IDS = tuple(
    "pathA42.py::channel." + name
    for name in ("h_EP", "radiation", "universality", "u_L_EP", "preferred_frame")
)
SOURCE_ADJUDICATION_IDS = (
    "pathA42.py::verdict_from_channels",
    "pathA42.py::ep_adjudication",
    "pathA42.py::guard_A_serialization_predicate",
)
SOURCE_CONTROL_IDS = tuple(
    "pathA42.py::control." + name
    for name in (
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K_without_pinned_Kh",
        "K",
    )
)
SOURCE_DELETION_IDS = tuple(
    "pathA42.py::deletion." + name
    for name in (
        "h_EP",
        "radiation",
        "universality",
        "u_L_EP",
        "preferred_frame",
        "gated_channels_collective",
    )
)
REPLACED_SOURCE_IDS = (
    "pathA42.wl::thin_mirror_coverage_gap",
)
SCOPED_SOURCE_IDS = (
    "pathA42.py::argparse_compare_harness",
    "pathA42.py::file_writing_and_report_persistence",
    "pathA42.py::main_and_run_math_script",
    "pathA42.py::comparison_payload_digest_sha_leaf_compare",
    "pathA42.py::substrate_report_reads_and_assertContains_scans",
    "pathA42.py::substrate_and_report_narration_objects",
    "pathA42.wl::Export_comparison_payload",
)
PRESERVED_SOURCE_IDS = (
    SOURCE_CORE_IDS
    + SOURCE_CHANNEL_IDS
    + SOURCE_ADJUDICATION_IDS
    + SOURCE_CONTROL_IDS
    + SOURCE_DELETION_IDS
)
SOURCE_PREDICATE_UNIVERSE = (
    PRESERVED_SOURCE_IDS + REPLACED_SOURCE_IDS + SCOPED_SOURCE_IDS
)
SOURCE_PREDICATE_TOTAL = 53
NEW_STAGE_IDS = (
    "stage042.dual::REACHABLE_FALSIFIERS",
    "stage042.dual::LORENTZ_NECESSARY_NOT_SUFFICIENT",
    "stage042.dual::GUARD_A_SCOPE_DENYLIST",
    "stage042.dual::DIMENSION_HOMOGENEITY",
    "stage042.dual::DUAL_ENGINE_TERMS",
    "stage042.dual::SOURCE_TO_STAGE_MANIFEST",
)

if len(SOURCE_PREDICATE_UNIVERSE) != SOURCE_PREDICATE_TOTAL:
    raise RuntimeError("stage042 source predicate universe count is not 53")


def ratified_category(identifier: str) -> str:
    if identifier in PRESERVED_SOURCE_IDS:
        return "preserved-folded"
    if identifier in REPLACED_SOURCE_IDS:
        return "replaced-by-stronger"
    if identifier in SCOPED_SOURCE_IDS:
        return "scoped-out"
    if identifier in NEW_STAGE_IDS:
        return "newly-added"
    raise KeyError(identifier)


RATIFIED_MANIFEST_MAP = {
    identifier: ratified_category(identifier)
    for identifier in SOURCE_PREDICATE_UNIVERSE + NEW_STAGE_IDS
}


def computed_manifest_rows(*, mutate: bool = False) -> tuple[tuple[str, str], ...]:
    rows: list[tuple[str, str]] = []
    moved = PRESERVED_SOURCE_IDS[0]
    for identifier in SOURCE_PREDICATE_UNIVERSE + NEW_STAGE_IDS:
        category = ratified_category(identifier)
        if mutate and identifier == moved:
            category = "replaced-by-stronger"
        rows.append((identifier, category))
    return tuple(rows)


def manifest_guards(*, mutate: bool = False) -> tuple[dict[str, Any], dict[str, Any]]:
    rows = computed_manifest_rows(mutate=mutate)
    identifiers = [identifier for identifier, _ in rows]
    source_identifiers = [
        identifier
        for identifier in identifiers
        if identifier in SOURCE_PREDICATE_UNIVERSE
    ]
    source_counts = Counter(source_identifiers)
    categories = {
        category: {
            identifier
            for identifier, item_category in rows
            if item_category == category
        }
        for category in (
            "preserved-folded",
            "replaced-by-stronger",
            "newly-added",
            "scoped-out",
        )
    }
    names = tuple(categories)
    disjoint = all(
        categories[left].isdisjoint(categories[right])
        for index, left in enumerate(names)
        for right in names[index + 1 :]
    )
    coverage = {
        "ok": (
            len(source_identifiers) == SOURCE_PREDICATE_TOTAL
            and set(source_identifiers) == set(SOURCE_PREDICATE_UNIVERSE)
            and all(count == 1 for count in source_counts.values())
            and set(identifiers) == set(RATIFIED_MANIFEST_MAP)
            and disjoint
        ),
        "source_count": len(source_identifiers),
        "disjoint": disjoint,
        "partition": Counter(category for _, category in rows),
    }
    category_map = dict(rows)
    category = {
        "ok": category_map == RATIFIED_MANIFEST_MAP,
        "computed": category_map,
        "ratified": RATIFIED_MANIFEST_MAP,
    }
    return coverage, category


def canonical_inventory(
    core: Core,
    production_channels: Channels,
    production_verdict: str,
    negative_guard_vector: tuple[str, ...],
    direct_guard: str,
    forged_guard: tuple[str, str],
    bypass_guard: tuple[str, str],
    falsifiers: dict[str, str],
    deletion: dict[str, Any],
    *,
    mutate: bool = False,
) -> tuple[tuple[str, Any], ...]:
    terms: list[tuple[str, Any]] = [
        ("verdict", production_verdict),
        ("channels", tuple(channel_states(production_channels).items())),
        ("det", sp.factor(core.det)),
        ("Aqq", sp.factor(core.Aqq)),
        ("Aqm", sp.factor(core.Aqm)),
        ("Amm", sp.factor(core.Amm)),
        (
            "exponents",
            (core.bare_exp, core.pinned_exp, core.corrupt_exp, core.wrong_exp),
        ),
        ("guard_negative", negative_guard_vector),
        ("guard_direct", direct_guard),
        ("guard_forged", forged_guard),
        ("guard_bypass", bypass_guard),
        ("reachable", tuple(sorted(falsifiers.items()))),
        ("deletion_h", deletion["h_EP"]),
        ("deletion_individual", tuple(deletion["individual"].items())),
        ("deletion_collective", deletion["collective"]),
    ]
    if mutate:
        terms = [term for term in terms if term[0] != "det"]
    return tuple(terms)


def expected_inventory(core: Core) -> tuple[tuple[str, Any], ...]:
    return (
        ("verdict", VERDICT_MAPPED),
        (
            "channels",
            (
                ("h_EP", H_EP_SAFE),
                ("radiation", RADIATION_SIM),
                ("universality", UNIVERSALITY_SIM),
                ("u_L_EP", UL_SIM),
                ("preferred_frame", PF_SIM),
            ),
        ),
        ("det", core.B * core.K - core.C**2),
        ("Aqq", sp.factor(core.typed_Aqq)),
        ("Aqm", sp.factor(core.typed_Aqm)),
        ("Amm", sp.factor(core.typed_Amm)),
        ("exponents", (3, 1, 1, 1)),
        ("guard_negative", ("LAUNDERING_REFUSED",) * 5),
        ("guard_direct", "LAUNDERING_REFUSED"),
        (
            "guard_forged",
            ("LAUNDERING_REFUSED", "MAGNITUDE_EMISSION_ALLOWED"),
        ),
        (
            "guard_bypass",
            ("LAUNDERING_REFUSED", "LAUNDERING_REFUSED"),
        ),
        (
            "reachable",
            tuple(
                sorted(
                    {
                        "B": "FALSIFIABLE_TENSION(channel=h_EP)",
                        "F": "FALSIFIABLE_TENSION(channel=universality)",
                        "G": VERDICT_NO_GO,
                        "H": VERDICT_NO_GO,
                        "K": "FALSIFIABLE_TENSION(channel=radiation)",
                    }.items()
                )
            ),
        ),
        ("deletion_h", (VERDICT_NO_GO, True)),
        (
            "deletion_individual",
            (
                ("radiation", (VERDICT_MAPPED, False)),
                ("universality", (VERDICT_MAPPED, False)),
                ("u_L_EP", (VERDICT_MAPPED, False)),
                ("preferred_frame", (VERDICT_MAPPED, False)),
            ),
        ),
        ("deletion_collective", ("NATURALLY_HIDDEN", True)),
    )


def inventories_equal(
    actual: tuple[tuple[str, Any], ...],
    expected: tuple[tuple[str, Any], ...],
) -> bool:
    if tuple(name for name, _ in actual) != tuple(name for name, _ in expected):
        return False
    for (_, left), (_, right) in zip(actual, expected):
        if isinstance(left, sp.Basic) or isinstance(right, sp.Basic):
            if not zero(sp.sympify(left) - sp.sympify(right)):
                return False
        elif left != right:
            return False
    return True


def verdict_witnesses(core: Core, *, mutate: bool = False) -> tuple[str, ...]:
    production = verdict_from_channels(compute_channels(Ledger(), core))
    b_ledger = Ledger(
        h_bare_mass_residue_zero=True if mutate else False,
        h_bare_mass_residue_fixture=False if mutate else True,
    )
    b = verdict_from_channels(compute_channels(b_ledger, core))
    f = verdict_from_channels(
        compute_channels(Ledger(species_indexed_b_over_ell=True), core)
    )
    k = verdict_from_channels(
        compute_channels(
            Ledger(pinned_kh_fact=True, commit_cE_equals_cgamma=True),
            core,
        )
    )
    g = verdict_from_channels(
        compute_channels(Ledger(stability_violation=True), core)
    )
    h = verdict_from_channels(
        compute_channels(Ledger(static_coulomb_match=False), core)
    )
    bg = verdict_from_channels(
        compute_channels(
            Ledger(
                h_bare_mass_residue_zero=False,
                h_bare_mass_residue_fixture=True,
                stability_violation=True,
            ),
            core,
        )
    )
    bkf = verdict_from_channels(
        compute_channels(
            Ledger(
                h_bare_mass_residue_zero=False,
                h_bare_mass_residue_fixture=True,
                pinned_kh_fact=True,
                commit_cE_equals_cgamma=True,
                species_indexed_b_over_ell=True,
            ),
            core,
        )
    )
    all_earned = stub_channels(
        compute_channels(Ledger(), core),
        radiation=EARNED_SAFE,
        universality=EARNED_SAFE,
        u_L_EP=EARNED_SAFE,
        preferred_frame=EARNED_SAFE,
    )
    natural = verdict_from_channels(all_earned)
    return (production, b, f, k, g, h, bg, bkf, natural)


def run_assertions() -> str:
    core = derive_core()
    production_channels = compute_channels(Ledger(), core)
    production_verdict = verdict_from_channels(production_channels)

    section("Static exchange and source projections")
    typed_det = core.B * core.K - core.C**2
    if ACTIVE_MUTATION == "STATIC_DETERMINANT":
        typed_det += 1
    expect_bool(
        "STATIC_DETERMINANT",
        zero(core.det - typed_det),
        sp.factor(core.det - typed_det),
    )

    typed_aqq = core.typed_Aqq + (
        1 if ACTIVE_MUTATION == "A_QQ_RESIDUAL" else 0
    )
    expect_bool(
        "A_QQ_RESIDUAL",
        zero(core.Aqq - typed_aqq),
        sp.factor(core.Aqq - typed_aqq),
    )
    typed_aqm = core.typed_Aqm + (
        1 if ACTIVE_MUTATION == "A_QM_RESIDUAL" else 0
    )
    expect_bool(
        "A_QM_RESIDUAL",
        zero(core.Aqm - typed_aqm),
        sp.factor(core.Aqm - typed_aqm),
    )
    typed_amm = core.typed_Amm + (
        1 if ACTIVE_MUTATION == "A_MM_RESIDUAL" else 0
    )
    expect_bool(
        "A_MM_RESIDUAL",
        zero(core.Amm - typed_amm),
        sp.factor(core.Amm - typed_amm),
    )

    aqm_decoupled = core.Aqm.subs({core.C: 0, core.mh: 0})
    target_decoupled = core.qL * core.qM / core.B
    if ACTIVE_MUTATION == "DECOUPLED_A_QM_NO_H_MASS":
        target_decoupled += 1
    expect_bool(
        "DECOUPLED_A_QM_NO_H_MASS",
        zero(aqm_decoupled - target_decoupled),
        sp.factor(aqm_decoupled - target_decoupled),
    )

    aqm_mixed = core.Aqm.subs({core.qL: 0, core.mh: 0})
    target_mixed = -core.C * core.qh * core.qM / core.det
    if ACTIVE_MUTATION == "MIXED_QL0_EP_RISK_TERM":
        target_mixed += 1
    expect_bool(
        "MIXED_QL0_EP_RISK_TERM",
        zero(aqm_mixed - target_mixed),
        sp.factor(aqm_mixed - target_mixed),
    )

    amm_decoupled = core.Amm.subs({core.C: 0, core.mh: 0})
    target_amm_decoupled = core.qM**2 / core.B
    if ACTIVE_MUTATION == "DECOUPLED_A_MM_DENSITY":
        target_amm_decoupled += 1
    expect_bool(
        "DECOUPLED_A_MM_DENSITY",
        zero(amm_decoupled - target_amm_decoupled),
        sp.factor(amm_decoupled - target_amm_decoupled),
    )

    aqm_mass = core.Aqm.subs(core.mh, 0)
    odd_flip = (
        -core.qM + 1
        if ACTIVE_MUTATION == "QM_FLIP_A_QM_ODD"
        else -core.qM
    )
    odd_residual = sp.factor(aqm_mass.subs(core.qM, odd_flip) + aqm_mass)
    expect_bool("QM_FLIP_A_QM_ODD", zero(odd_residual), odd_residual)

    amm_mass = core.Amm.subs(core.mh, 0)
    signed_amm = sp.factor(amm_mass / core.qM)
    even_flip = (
        -core.qM + 1
        if ACTIVE_MUTATION == "A_MM_EVEN_SIGNED_FLIP"
        else -core.qM
    )
    even_residual = sp.factor(
        amm_mass.subs(core.qM, even_flip) - amm_mass
    )
    signed_residual = sp.factor(
        signed_amm.subs(core.qM, even_flip) + signed_amm
    )
    expect_bool(
        "A_MM_EVEN_SIGNED_FLIP",
        zero(even_residual) and zero(signed_residual),
        {"even": even_residual, "signed": signed_residual},
    )

    h_unit = sp.Matrix([0, 1])
    mass_projection = (h_unit.T * sp.Matrix([core.qM, 0]))[0]
    fixture_projection = (
        h_unit.T * sp.Matrix([core.qM, core.mh])
    )[0]
    fixture_target = core.mh + (
        1 if ACTIVE_MUTATION == "H_MASS_PROJECTION" else 0
    )
    expect_bool(
        "H_MASS_PROJECTION",
        zero(mass_projection) and zero(fixture_projection - fixture_target),
        {
            "production": mass_projection,
            "fixture_residual": fixture_projection - fixture_target,
        },
    )

    section("Radiation speed scaling and strict stability")
    target_bare = (
        core.qh**2 / (core.Mh * core.QE**2)
        * (core.cg / core.cE) ** 3
    )
    if ACTIVE_MUTATION == "RADIATION_BARE_RATIO":
        target_bare += 1
    expect_bool(
        "RADIATION_BARE_RATIO",
        zero(core.ratio_bare - target_bare),
        sp.factor(core.ratio_bare - target_bare),
    )

    bare_ratio_for_tooth = core.ratio_bare * (
        core.cE if ACTIVE_MUTATION == "RADIATION_BARE_EXPONENT" else 1
    )
    bare_exp_for_tooth = speed_exponent(bare_ratio_for_tooth, core.cE)
    expect_bool(
        "RADIATION_BARE_EXPONENT",
        bare_exp_for_tooth == 3,
        bare_exp_for_tooth,
    )

    pinned_ratio_for_tooth = core.ratio_pinned
    if ACTIVE_MUTATION == "RADIATION_PINNED_EXPONENT":
        pinned_ratio_for_tooth = sp.factor(
            core.ratio_bare.subs(core.Mh, core.K / core.cE)
            * core.K
            / core.qh**2
        )
    pinned_exp_for_tooth = speed_exponent(
        pinned_ratio_for_tooth, core.cE
    )
    expect_bool(
        "RADIATION_PINNED_EXPONENT",
        pinned_exp_for_tooth == 1,
        pinned_exp_for_tooth,
    )

    wrong_ratio_for_tooth = (
        core.ratio_bare
        if ACTIVE_MUTATION == "RADIATION_WRONG_NORM_DISCRIMINATES"
        else core.ratio_wrong
    )
    wrong_exp_for_tooth = speed_exponent(
        wrong_ratio_for_tooth, core.cE
    )
    discriminator = (
        wrong_exp_for_tooth,
        core.bare_exp,
        wrong_exp_for_tooth != core.bare_exp,
    )
    expect_bool(
        "RADIATION_WRONG_NORM_DISCRIMINATES",
        discriminator == (1, 3, True),
        discriminator,
    )

    stability_target = typed_det + (
        1 if ACTIVE_MUTATION == "STABILITY_BOUND" else 0
    )
    stable_sample = core.det.subs({core.B: 2, core.K: 3, core.C: 1})
    violated_sample = core.det.subs(
        {core.B: 2, core.K: 2, core.C: 2}
    )
    stability_actual = (
        zero(core.det - stability_target),
        bool(stable_sample > 0),
        bool(violated_sample <= 0),
    )
    expect_bool(
        "STABILITY_BOUND",
        stability_actual == (True, True, True),
        stability_actual,
    )

    section("Five computed channels, subtags, and production verdict")
    h_case = compute_channels(
        Ledger(static_coulomb_match=ACTIVE_MUTATION != "CHANNEL_H_EP"),
        core,
    )
    expect_bool("CHANNEL_H_EP", h_case.h_EP == H_EP_SAFE, h_case.h_EP)

    radiation_case = compute_channels(
        Ledger(
            radiation_corrupt_speed=ACTIVE_MUTATION == "CHANNEL_RADIATION"
        ),
        core,
    )
    expect_bool(
        "CHANNEL_RADIATION",
        radiation_case.radiation == RADIATION_SIM,
        radiation_case.radiation,
    )

    universality_case = compute_channels(
        Ledger(
            species_indexed_b_over_ell=(
                ACTIVE_MUTATION == "CHANNEL_UNIVERSALITY"
            )
        ),
        core,
    )
    expect_bool(
        "CHANNEL_UNIVERSALITY",
        universality_case.universality == UNIVERSALITY_SIM,
        universality_case.universality,
    )

    ul_case = compute_channels(
        Ledger(uL_no_go=ACTIVE_MUTATION == "CHANNEL_UL_EP"), core
    )
    expect_bool(
        "CHANNEL_UL_EP",
        ul_case.u_L_EP == UL_SIM,
        ul_case.u_L_EP,
    )

    pf_case = compute_channels(
        Ledger(
            suppression_requires_large_cE=(
                ACTIVE_MUTATION == "CHANNEL_PREFERRED_FRAME"
            )
        ),
        core,
    )
    expect_bool(
        "CHANNEL_PREFERRED_FRAME",
        pf_case.preferred_frame == PF_SIM,
        pf_case.preferred_frame,
    )

    subtag_case = compute_channels(
        Ledger(C_hu_nonzero=ACTIVE_MUTATION != "SUBTAGS"), core
    )
    subtag_actual = (
        subtag_case.cherenkov_deferred,
        subtag_case.mixed_scalar_EP_risk,
    )
    expect_bool(
        "SUBTAGS",
        subtag_actual == (True, True),
        subtag_actual,
    )

    ep_case = ep_adjudication(
        production_channels,
        unqualified=ACTIVE_MUTATION == "EP_ADJUDICATION",
    )
    ep_actual = (
        ep_case["mass_channel"],
        ep_case["full_decoupled_floor_EP_safety"],
        ep_case["unqualified_decoupled_floor_EP_safe_reported"],
    )
    expect_bool(
        "EP_ADJUDICATION",
        ep_actual == (H_EP_SAFE, "NOT_EARNED", False),
        ep_actual,
    )

    verdict_ledger = (
        Ledger(static_coulomb_match=False)
        if ACTIVE_MUTATION == "PRODUCTION_VERDICT"
        else Ledger()
    )
    verdict_case = verdict_from_channels(
        compute_channels(verdict_ledger, core)
    )
    expect_bool(
        "PRODUCTION_VERDICT",
        verdict_case == VERDICT_MAPPED,
        verdict_case,
    )

    section("Guard A scoped laundering-refusal denylist")
    production_tree = {
        "radiation": {
            "magnitude": "SIM_GATED_BY_GUARD_A_NO_NUMERIC_POWER_RATIO_EMITTED",
            "attempted_guarded_magnitude_emission": False,
        }
    }
    if ACTIVE_MUTATION == "GUARD_A_PRODUCTION":
        production_tree["radiation"]["P_h/P_EM"] = 1
    production_guard = guard_scan(production_tree, absent_substrate())
    expect_bool(
        "GUARD_A_PRODUCTION",
        production_guard["result"] == "NO_GUARDED_MAGNITUDE_EMITTED",
        production_guard,
    )

    negative_vector = guard_negative_vector(
        neutralized=ACTIVE_MUTATION == "GUARD_A_NEGATIVE_FIXTURES"
    )
    expect_bool(
        "GUARD_A_NEGATIVE_FIXTURES",
        negative_vector == ("LAUNDERING_REFUSED",) * 5,
        negative_vector,
    )

    direct_substrate = (
        earned_substrate()
        if ACTIVE_MUTATION == "GUARD_A_DIRECT_INJECTION"
        else absent_substrate()
    )
    direct_guard = guard_scan(
        {"P_h/P_EM": 1, "M_h": 1}, direct_substrate
    )["result"]
    expect_bool(
        "GUARD_A_DIRECT_INJECTION",
        direct_guard == "LAUNDERING_REFUSED",
        direct_guard,
    )

    forged_payload = {
        "earned_parent_action": True,
        "P_h_over_P_EM": 1,
    }
    forged_negative = guard_scan(
        forged_payload,
        absent_substrate(),
        allow_forged_local=ACTIVE_MUTATION == "GUARD_A_FORGED_FLAG",
    )["result"]
    forged_positive = guard_scan(
        {"P_h_over_P_EM": 1}, earned_substrate()
    )["result"]
    forged_actual = (forged_negative, forged_positive)
    expect_bool(
        "GUARD_A_FORGED_FLAG",
        forged_actual
        == ("LAUNDERING_REFUSED", "MAGNITUDE_EMISSION_ALLOWED"),
        forged_actual,
    )

    tuple_guard = guard_scan(
        {"assembled": {"power_ratio": (0.137,)}},
        absent_substrate(),
        descend_sequences=(
            ACTIVE_MUTATION != "GUARD_A_BYPASS_REGRESSIONS"
        ),
    )["result"]
    string_guard = guard_scan(
        {"assembled": {"power_ratio": "P_h/P_EM = 0.137"}},
        absent_substrate(),
    )["result"]
    bypass_actual = (tuple_guard, string_guard)
    expect_bool(
        "GUARD_A_BYPASS_REGRESSIONS",
        bypass_actual
        == ("LAUNDERING_REFUSED", "LAUNDERING_REFUSED"),
        bypass_actual,
    )

    scope_extra = (
        frozenset({"diagnostic_gain"})
        if ACTIVE_MUTATION == "GUARD_A_SCOPE_DENYLIST"
        else frozenset()
    )
    scope_guard = guard_scan(
        {"diagnostic_gain": 7},
        absent_substrate(),
        extra_fields=scope_extra,
    )["result"]
    expect_bool(
        "GUARD_A_SCOPE_DENYLIST",
        scope_guard == "NO_GUARDED_MAGNITUDE_EMITTED",
        scope_guard,
    )

    section("Twelve source controls (A is Guard A; B through K recomputed)")
    normal_controls = controls(core)
    control_case = controls(core, mutate=ACTIVE_MUTATION)

    b_actual = (
        control_case["B"]["channels"].h_EP,
        control_case["B"]["verdict"],
    )
    expect_bool(
        "CTRL_B",
        b_actual
        == (H_EP_FIFTH, "FALSIFIABLE_TENSION(channel=h_EP)"),
        b_actual,
    )

    c_actual = (
        control_case["C"]["channels"].mixed_scalar_EP_risk,
        control_case["C"]["verdict"],
    )
    expect_bool(
        "CTRL_C",
        c_actual == (True, VERDICT_MAPPED),
        c_actual,
    )

    d_actual = (
        control_case["D"]["baseline_exp"],
        control_case["D"]["selected_exp"],
        control_case["D"]["channels"].radiation,
        control_case["D"]["verdict"],
    )
    expect_bool(
        "CTRL_D",
        d_actual == (3, 1, RADIATION_AUDIT, VERDICT_MAPPED),
        d_actual,
    )

    e_actual = (
        control_case["E"]["channels"].preferred_frame,
        control_case["E"]["verdict"],
    )
    expect_bool(
        "CTRL_E",
        e_actual == (PF_TENSION, VERDICT_MAPPED),
        e_actual,
    )

    f_actual = (
        control_case["F"]["channels"].universality,
        control_case["F"]["verdict"],
    )
    expect_bool(
        "CTRL_F",
        f_actual
        == (
            UNIVERSALITY_TENSION,
            "FALSIFIABLE_TENSION(channel=universality)",
        ),
        f_actual,
    )

    g_actual = (
        control_case["G"]["channels"].stability_violated,
        control_case["G"]["verdict"],
    )
    expect_bool(
        "CTRL_G",
        g_actual == (True, VERDICT_NO_GO),
        g_actual,
    )

    h_actual = (
        control_case["H"]["channels"].h_EP,
        control_case["H"]["verdict"],
    )
    expect_bool(
        "CTRL_H",
        h_actual == ("NO_GO", VERDICT_NO_GO),
        h_actual,
    )

    i_actual = (
        control_case["I"]["aqm_flipped"],
        control_case["I"]["signed_flipped"],
        control_case["I"]["amm_even"],
        control_case["I"]["verdict"],
    )
    expect_bool(
        "CTRL_I",
        i_actual == (True, True, True, VERDICT_MAPPED),
        i_actual,
    )

    j_actual = (
        control_case["J"]["baseline_exp"],
        control_case["J"]["selected_exp"],
        control_case["J"]["verdict"],
    )
    expect_bool(
        "CTRL_J",
        j_actual == (3, 1, VERDICT_MAPPED),
        j_actual,
    )

    k0_actual = (
        control_case["K_without_pinned_Kh"]["channels"].radiation,
        control_case["K_without_pinned_Kh"]["verdict"],
    )
    expect_bool(
        "CTRL_K_WITHOUT_PINNED_KH",
        k0_actual == (RADIATION_SIM, VERDICT_MAPPED),
        k0_actual,
    )

    k_actual = (
        control_case["K"]["channels"].radiation,
        control_case["K"]["verdict"],
    )
    expect_bool(
        "CTRL_K",
        k_actual
        == (
            RADIATION_TENSION,
            "FALSIFIABLE_TENSION(channel=radiation)",
        ),
        k_actual,
    )

    section("Deletion sensitivity and reachable falsifiers")
    normal_deletion = deletion_outcomes(production_channels)
    deletion_case = deletion_outcomes(
        production_channels, mutate=ACTIVE_MUTATION
    )
    expect_bool(
        "DELETION_H_EP_PILLAR",
        deletion_case["h_EP"] == (VERDICT_NO_GO, True),
        deletion_case["h_EP"],
    )
    expected_individual = {
        name: (VERDICT_MAPPED, False)
        for name in (
            "radiation",
            "universality",
            "u_L_EP",
            "preferred_frame",
        )
    }
    expect_bool(
        "DELETION_GATED_INDIVIDUAL_STABLE",
        deletion_case["individual"] == expected_individual,
        deletion_case["individual"],
    )
    expect_bool(
        "DELETION_COLLECTIVE_NATURALLY_HIDDEN",
        deletion_case["collective"] == ("NATURALLY_HIDDEN", True),
        deletion_case["collective"],
    )

    falsifiers = reachable_falsifier_map(normal_controls)
    if ACTIVE_MUTATION == "REACHABLE_FALSIFIERS":
        falsifiers = dict(falsifiers)
        falsifiers.pop("B")
    expected_falsifiers = {
        "B": "FALSIFIABLE_TENSION(channel=h_EP)",
        "F": "FALSIFIABLE_TENSION(channel=universality)",
        "G": VERDICT_NO_GO,
        "H": VERDICT_NO_GO,
        "K": "FALSIFIABLE_TENSION(channel=radiation)",
    }
    mapped_controls = {
        name: normal_controls[name]["verdict"]
        for name in (
            "C",
            "D",
            "E",
            "I",
            "J",
            "K_without_pinned_Kh",
        )
    }
    reachable_actual = (
        falsifiers,
        mapped_controls,
    )
    reachable_expected = (
        expected_falsifiers,
        {name: VERDICT_MAPPED for name in mapped_controls},
    )
    expect_bool(
        "REACHABLE_FALSIFIERS",
        reachable_actual == reachable_expected,
        reachable_actual,
    )

    commit_only = compute_channels(
        Ledger(commit_cE_equals_cgamma=True),
        core,
        commit_alone_sufficient=(
            ACTIVE_MUTATION == "LORENTZ_NECESSARY_NOT_SUFFICIENT"
        ),
    )
    lorentz_actual = (
        commit_only.radiation,
        commit_only.preferred_frame,
        verdict_from_channels(commit_only),
    )
    expect_bool(
        "LORENTZ_NECESSARY_NOT_SUFFICIENT",
        lorentz_actual == (RADIATION_SIM, PF_SIM, VERDICT_MAPPED),
        lorentz_actual,
    )

    section("Build-global guards and first-match witness table")
    dims = dimension_guard(
        mutate=ACTIVE_MUTATION == "DIMENSION_HOMOGENEITY"
    )
    expect_bool(
        "DIMENSION_HOMOGENEITY",
        dims["homogeneous"],
        dims,
    )

    forged_inventory = (
        guard_scan(forged_payload, absent_substrate())["result"],
        guard_scan(
            {"P_h_over_P_EM": 1}, earned_substrate()
        )["result"],
    )
    bypass_inventory = (
        guard_scan(
            {"assembled": {"power_ratio": (0.137,)}},
            absent_substrate(),
        )["result"],
        guard_scan(
            {"assembled": {"power_ratio": "P_h/P_EM = 0.137"}},
            absent_substrate(),
        )["result"],
    )
    inventory_actual = canonical_inventory(
        core,
        production_channels,
        production_verdict,
        guard_negative_vector(),
        guard_scan(
            {"P_h/P_EM": 1, "M_h": 1}, absent_substrate()
        )["result"],
        forged_inventory,
        bypass_inventory,
        reachable_falsifier_map(normal_controls),
        normal_deletion,
        mutate=ACTIVE_MUTATION == "DUAL_ENGINE_TERMS",
    )
    inventory_expected = expected_inventory(core)
    expect_bool(
        "DUAL_ENGINE_TERMS",
        inventories_equal(inventory_actual, inventory_expected),
        {
            "actual_names": tuple(name for name, _ in inventory_actual),
            "expected_names": tuple(name for name, _ in inventory_expected),
        },
    )

    witness_actual = verdict_witnesses(
        core, mutate=ACTIVE_MUTATION == "VERDICT_REDERIVATION"
    )
    witness_expected = (
        VERDICT_MAPPED,
        "FALSIFIABLE_TENSION(channel=h_EP)",
        "FALSIFIABLE_TENSION(channel=universality)",
        "FALSIFIABLE_TENSION(channel=radiation)",
        VERDICT_NO_GO,
        VERDICT_NO_GO,
        VERDICT_NO_GO,
        "FALSIFIABLE_TENSION(channel=h_EP,radiation,universality)",
        "NATURALLY_HIDDEN",
    )
    expect_bool(
        "VERDICT_REDERIVATION",
        witness_actual == witness_expected,
        {"actual": witness_actual, "ratified": witness_expected},
    )

    section("Bounded source-to-stage predicate manifest")
    coverage_guard, category_guard = manifest_guards(
        mutate=ACTIVE_MUTATION == "SOURCE_TO_STAGE_MANIFEST"
    )
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        coverage_guard["ok"],
        {"assertion": "coverage/count", **coverage_guard},
    )
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        category_guard["ok"],
        {"assertion": "per-predicate category map"},
    )

    print("")
    print(f"VERDICT={production_verdict}")
    for name, state in channel_states(production_channels).items():
        print(f"CHANNEL_{name}={state}")
    print(
        "EP_ADJUDICATION=mass-channel-only; "
        "full-decoupled-floor-EP=NOT_EARNED; "
        "unqualified-safe-report=False"
    )
    print(
        "RADIATION_EXPONENTS="
        f"bare:{core.bare_exp};pinned:{core.pinned_exp};"
        f"corrupt:{core.corrupt_exp};wrong_norm:{core.wrong_exp}"
    )
    print(
        "RADIATION_MAGNITUDE="
        "SIM_GATED_BY_GUARD_A_NO_NUMERIC_POWER_RATIO_EMITTED"
    )
    print(
        "REACHABLE_FALSIFIERS="
        + ",".join(
            f"{name}->{verdict}"
            for name, verdict in sorted(
                reachable_falsifier_map(normal_controls).items()
            )
        )
    )
    print(
        "NATURALLY_HIDDEN=DEFINED_BUT_UNREACHABLE_IN_PRODUCTION;"
        "collective-four-channel-earned-stub-only"
    )
    print(
        "LORENTZ=c_E=c_gamma_NECESSARY_NOT_SUFFICIENT;"
        "real-radiating-extra-scalar;magnitude-remains-SIM_GATED-without-pinned-K_h"
    )
    print(
        "GUARD_A_SCOPE=DENYLIST_ONLY:"
        "{M_h,c_E,K_h,P_h/P_EM,EP-magnitude,residue-floor};"
        "unrelated-numeric-fields-pass"
    )
    print(
        "GATE_L_CONNECTION=EARNED_CONNECTION;"
        "same-embedding-direction-family-as-light"
    )
    print(
        "HARD_WALL=deferred-throat-embedding-solve-pins-{M_h,c_E,K_h};"
        "on-pinning-check-stellar-cooling-and-BBN/CMB"
    )
    print(
        "FRAMING=FIRST-CLASS_CHARACTERIZED_DEPARTURE;"
        "MAPPED+SIM-GATED_NOT_RESOLVED;several-reachable-falsifiers"
    )
    print(f"SOURCE_PREDICATE_TOTAL={SOURCE_PREDICATE_TOTAL}")
    print("EXECUTABLE_TOOTH_TOTAL=49")
    return production_verdict


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage042_charge_coupled_scalar SymPy audit")
    print(
        "ROUTE=Matrix.inv symbolic exchange + denominator-minus-numerator "
        "speed degrees + dataclass channel cascade + recursive Guard-A walk"
    )
    print("FILE_IO=none; CROSS_ENGINE_COMPARE=none")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(
            f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}"
        )

    verdict = run_assertions()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print("")
    print(f"TOOTH_COUNT={len(TOOTH_ORDER)}")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print(f"OVERALL PASS: SymPy independently reached {verdict}")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(
            "OVERALL FAIL: SymPy stage042 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
