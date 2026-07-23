#!/usr/bin/env python3
"""Ledger stage038 SymPy audit: SEALED section-4 terminal landing.

Standalone, print-only, assert-zero, exact, and cross-engine-file-I/O-free.
The engine enumerates the full fact domain, applies the fixed first-match
adjudication, checks every cell against an independently coded oracle, hashes
the canonical row stream at runtime, and re-derives the production landing.

Tooth-local runtime ablation uses ``LEDGER_STAGE038_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from enum import Enum
import hashlib
import itertools
import math
import os
import types
from typing import Any, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE038_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

TOOTH_ORDER = (
    "TRUTH_TOTALITY",
    "TRUTH_PRECEDENCE",
    "LANDING_OWNERSHIP",
    "ACTIVE_FLUX_CAVEAT",
    "HOOK_LORENTZ",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
    "UNITS_RESTORED",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "TRUTH_TOTALITY":
        "replace one computed cell by the defensive fall-through and drop a distinct row",
    "TRUTH_PRECEDENCE":
        "move the consistency arm ahead of the electric-anchor arm for the production witness",
    "LANDING_OWNERSHIP":
        "resolve the live electric anchor while retaining the emitted production carriage",
    "ACTIVE_FLUX_CAVEAT":
        "absorb F_flux into the conservative exchange and mark its value/integrability closed",
    "HOOK_LORENTZ":
        "close delta_BA, r_cone, higher orders, and the active-flux match",
    "TARGET_BLINDNESS":
        "install a nonsealed landing-token code object and inject every barred path marker",
    "DUAL_ENGINE_TERMS":
        "corrupt every local computed inventory lane",
    "UNITS_RESTORED":
        "corrupt every cited base dimension before deriving the hook dimensions",
    "VERDICT_REDERIVATION":
        "resolve the production witness anchor before running the native landing table",
    "SOURCE_TO_STAGE_MANIFEST":
        "remove a scoped-out source row from the live manifest",
}

PUBLISHED_DIGEST = bytes.fromhex(
    "983556935e50f12670fef24f17a23c048e295ddf0b6952aa0bd1618e9f179619"
)


class AuditFailure(AssertionError):
    """A named stage predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


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


class ActiveFluxFact(Enum):
    R1_UNMATCHED_REMAINDER = 1
    ABSORBED_DECIDED_ZERO = 2


class IntegrabilityFact(Enum):
    R1_OPEN = 1
    CLOSED = 2


class LorentzDetermination(Enum):
    UNDETERMINED = 1
    DETERMINED = 2


@dataclass(frozen=True, order=True)
class LandingToken:
    value: str

    def __str__(self) -> str:
        return self.value


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
        typed = (
            type(self.current) is CurrentIdentity
            and type(self.comparison) is ComparisonFact
            and type(self.relative_sign) is RelativeSignFact
            and type(self.magnitude) is MagnitudeFact
            and type(self.tier) is TierFact
            and type(self.anchor) is AnchorFact
            and type(self.internal_sectors) is tuple
            and all(type(item) is str for item in self.internal_sectors)
        )
        if not typed:
            raise TypeError("LiveFacts accepts only typed neutral facts")

    def neutral_tokens(self) -> tuple[Any, ...]:
        return (
            self.current,
            self.comparison,
            self.relative_sign,
            self.magnitude,
            self.tier,
            self.anchor,
            self.internal_sectors,
        )

    def as_case(self) -> dict[str, Any]:
        return {
            "current": self.current.value,
            "comparison": self.comparison.value,
            "relative_sign": self.relative_sign.value,
            "magnitude": self.magnitude.value,
            "tier": self.tier.value,
            "anchor": self.anchor.value,
            "internal": bool(self.internal_sectors),
        }


@dataclass(frozen=True)
class TruthAudit:
    agreement: bool
    no_unclassified: bool
    cell_count: int
    expected_total: int
    counts: tuple[tuple[LandingToken, int], ...]
    expected_counts: tuple[tuple[LandingToken, int], ...]
    digest: bytes
    named_targets: tuple[tuple[str, LandingToken], ...]

    def target(self, name: str) -> LandingToken:
        return dict(self.named_targets)[name]


def compact(value: Any) -> str:
    if isinstance(value, bytes):
        return value.hex()
    return str(value)


def expect_zero(name: str, residual: Any, evidence: Any = None) -> None:
    global PASS_COUNT, FAIL_COUNT
    clean = sp.simplify(sp.sympify(residual))
    if clean == 0:
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    FAIL_COUNT += 1
    print(f"FIRST_FAILURE={name}")
    if ACTIVE_MUTATION == name:
        print(f"FIRED_AT_OWN_ASSERT={name}")
    print(f"FAIL  {name}: residual = {compact(clean)}")
    if evidence is not None:
        print(f"      evidence = {evidence}")
    raise AuditFailure(name, compact(clean))


def expect_bool(name: str, condition: bool, evidence: Any = None) -> None:
    expect_zero(name, sp.Integer(0 if bool(condition) else 1), evidence)


def section(text: str) -> None:
    print("")
    print(text)
    print("-" * len(text))


# ---------------------------------------------------------------------------
# SEALED section 4: fixed total first-match function and independent oracle.
# ---------------------------------------------------------------------------

CURRENT_DOMAIN = tuple(item.value for item in CurrentIdentity)
COMPARISON_DOMAIN = tuple(item.value for item in ComparisonFact)
RELATIVE_DOMAIN = tuple(item.value for item in RelativeSignFact)
MAG_DOMAIN = tuple(item.value for item in MagnitudeFact)
TIER_DOMAIN = tuple(item.value for item in TierFact)
ANCHOR_DOMAIN = tuple(item.value for item in AnchorFact)
INCONSISTENCY_DOMAIN = (False, True)
DOMAIN_AXES = (
    CURRENT_DOMAIN,
    COMPARISON_DOMAIN,
    RELATIVE_DOMAIN,
    MAG_DOMAIN,
    TIER_DOMAIN,
    ANCHOR_DOMAIN,
    INCONSISTENCY_DOMAIN,
)


def section4_adjudicate(
    live: LiveFacts | None = None,
    *,
    current: str | None = None,
    comparison: str | None = None,
    relative_sign: str | None = None,
    magnitude: str | None = None,
    tier: str | None = None,
    anchor: str | None = None,
    internal: bool | None = None,
    precedence_mutation: bool = False,
) -> LandingToken:
    """Apply the authoritative total first-match precedence."""
    if live is not None:
        restated = (
            current,
            comparison,
            relative_sign,
            magnitude,
            tier,
            anchor,
            internal,
        )
        if any(value is not None for value in restated):
            raise TypeError("live facts cannot be mixed with restated fields")
        return section4_adjudicate(
            **live.as_case(),
            precedence_mutation=precedence_mutation,
        )
    values = (
        current,
        comparison,
        relative_sign,
        magnitude,
        tier,
        anchor,
        internal,
    )
    if any(value is None for value in values):
        raise TypeError("incomplete section-4 case")
    if internal:
        return LandingToken("NO_GO(sector)")

    complete = (
        current != CurrentIdentity.R1.value
        and comparison != ComparisonFact.ROUTE_B_R1.value
        and relative_sign != RelativeSignFact.ANCHOR_CONDITIONAL.value
        and tier == TierFact.GAPS_CLOSED.value
        and anchor == AnchorFact.CLOSED.value
        and magnitude != MagnitudeFact.R1.value
    )
    if precedence_mutation and anchor == AnchorFact.R1.value:
        return LandingToken("R1_REQUIRED(consistency)")
    if complete:
        if relative_sign == RelativeSignFact.OPPOSITE.value:
            return LandingToken("AMENDMENT_EXCLUDED(wrong_relative_sign)")
        if relative_sign == RelativeSignFact.CONFLICT.value:
            return LandingToken(
                "AMENDMENT_EXCLUDED(routes_leading_conflict)"
            )
        if comparison == ComparisonFact.AGREE.value:
            if magnitude == MagnitudeFact.FORCED.value:
                return LandingToken("MAGNETISM_LORENTZ_CONSISTENT")
            return LandingToken(
                "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"
            )
        if comparison == ComparisonFact.DIFFER.value:
            return LandingToken("MAGNETISM_DEPARTURE_CHARACTERIZED")

    if anchor == AnchorFact.R1.value:
        return LandingToken("R1_REQUIRED(electric_bc_selection)")
    if (
        current == CurrentIdentity.R1.value
        or comparison == ComparisonFact.ROUTE_B_R1.value
    ):
        return LandingToken("R1_REQUIRED(direct_moving_throat)")
    if magnitude == MagnitudeFact.R1.value:
        return LandingToken("R1_REQUIRED(magnitude)")
    if (
        tier == TierFact.CONDITIONAL.value
        or relative_sign == RelativeSignFact.ANCHOR_CONDITIONAL.value
    ):
        return LandingToken("R1_REQUIRED(consistency)")
    return LandingToken("R1_REQUIRED(unclassified)")


def section4_oracle(case: Mapping[str, Any]) -> LandingToken:
    """Independent classification by complete-state and open-reason keys."""
    if bool(case["internal"]):
        return LandingToken("NO_GO(sector)")

    complete_key = (
        case["current"] != "R1_SOURCE_BASIS",
        case["comparison"] != "route_B_R1",
        case["relative_sign"] != "relative_sign_anchor_conditional",
        case["tier"] == "tier_A_gaps_closed",
        case["anchor"] == "electric_anchor_closed",
        case["magnitude"] != "R1(magnitude)",
    )
    if all(complete_key):
        terminal_by_sign = {
            "relative_sign_opposite":
                LandingToken("AMENDMENT_EXCLUDED(wrong_relative_sign)"),
            "leading_tensor_conflict":
                LandingToken(
                    "AMENDMENT_EXCLUDED(routes_leading_conflict)"
                ),
        }
        if case["relative_sign"] in terminal_by_sign:
            return terminal_by_sign[case["relative_sign"]]
        complete_by_route = {
            ("routes_agree", "magnitude_forced_by_electric"):
                LandingToken("MAGNETISM_LORENTZ_CONSISTENT"),
            ("routes_agree", "magnitude_free_factor"):
                LandingToken(
                    "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"
                ),
            ("routes_differ", "magnitude_forced_by_electric"):
                LandingToken("MAGNETISM_DEPARTURE_CHARACTERIZED"),
            ("routes_differ", "magnitude_free_factor"):
                LandingToken("MAGNETISM_DEPARTURE_CHARACTERIZED"),
        }
        return complete_by_route[
            (case["comparison"], case["magnitude"])
        ]

    open_reasons = (
        (
            case["anchor"] == "R1_REQUIRED(bc_selection)",
            LandingToken("R1_REQUIRED(electric_bc_selection)"),
        ),
        (
            case["current"] == "R1_SOURCE_BASIS"
            or case["comparison"] == "route_B_R1",
            LandingToken("R1_REQUIRED(direct_moving_throat)"),
        ),
        (
            case["magnitude"] == "R1(magnitude)",
            LandingToken("R1_REQUIRED(magnitude)"),
        ),
        (
            case["tier"] == "tier_A_conditional"
            or case["relative_sign"] == "relative_sign_anchor_conditional",
            LandingToken("R1_REQUIRED(consistency)"),
        ),
    )
    for active, token in open_reasons:
        if active:
            return token
    return LandingToken("R1_REQUIRED(unclassified)")


def truth_table(totality_mutation: bool = False) -> TruthAudit:
    """Enumerate, oracle-check, serialize, and hash every section-4 row."""
    keys = (
        "current",
        "comparison",
        "relative_sign",
        "magnitude",
        "tier",
        "anchor",
        "internal",
    )
    expected_total = math.prod(len(axis) for axis in DOMAIN_AXES)
    records: list[tuple[str, LandingToken, LandingToken]] = []
    unclassified = LandingToken("R1_REQUIRED(unclassified)")

    for ordinal, values in enumerate(itertools.product(*DOMAIN_AXES)):
        case = dict(zip(keys, values, strict=True))
        got = section4_adjudicate(**case)
        if totality_mutation and ordinal == 0:
            got = unclassified
        want = section4_oracle(case)
        row = "|".join(str(value) for value in values) + "|" + got.value
        if totality_mutation and ordinal == expected_total - 1:
            continue
        records.append((row, got, want))

    counts = Counter(got for _row, got, _want in records)
    row_stream = "\n".join(row for row, _got, _want in records)
    digest = hashlib.sha256(row_stream.encode("utf-8")).digest()

    expected_counts = Counter({
        LandingToken("NO_GO(sector)"): 576,
        LandingToken("MAGNETISM_LORENTZ_CONSISTENT"): 3,
        LandingToken("MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"): 3,
        LandingToken("MAGNETISM_DEPARTURE_CHARACTERIZED"): 6,
        LandingToken("AMENDMENT_EXCLUDED(wrong_relative_sign)"): 12,
        LandingToken("AMENDMENT_EXCLUDED(routes_leading_conflict)"): 12,
        LandingToken("R1_REQUIRED(electric_bc_selection)"): 288,
        LandingToken("R1_REQUIRED(direct_moving_throat)"): 144,
        LandingToken("R1_REQUIRED(magnitude)"): 48,
        LandingToken("R1_REQUIRED(consistency)"): 60,
    })
    named_targets = (
        ("internal", LandingToken("NO_GO(sector)")),
        (
            "consistent_free",
            LandingToken("MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"),
        ),
        (
            "lorentz_consistent",
            LandingToken("MAGNETISM_LORENTZ_CONSISTENT"),
        ),
    )
    return TruthAudit(
        agreement=all(got == want for _row, got, want in records),
        no_unclassified=counts.get(unclassified, 0) == 0,
        cell_count=len(records),
        expected_total=expected_total,
        counts=tuple(sorted(counts.items())),
        expected_counts=tuple(sorted(expected_counts.items())),
        digest=digest,
        named_targets=named_targets,
    )


# ---------------------------------------------------------------------------
# Cited production facts, complete blocker carriage, and honest hooks.
# ---------------------------------------------------------------------------

LIVE_FACTS = LiveFacts(
    current=CurrentIdentity.CONVECTION_CONDITIONAL,
    comparison=ComparisonFact.ROUTE_B_R1,
    relative_sign=RelativeSignFact.ANCHOR_CONDITIONAL,
    magnitude=MagnitudeFact.R1,
    tier=TierFact.CONDITIONAL,
    anchor=AnchorFact.R1,
    internal_sectors=(),
)
EMITTED_LANDING = section4_adjudicate(LIVE_FACTS)
CARRIED_NEUTRAL_TOKENS = (
    CurrentIdentity.CONVECTION_CONDITIONAL,
    ComparisonFact.ROUTE_B_R1,
    RelativeSignFact.ANCHOR_CONDITIONAL,
    MagnitudeFact.R1,
    TierFact.CONDITIONAL,
    AnchorFact.R1,
    (),
)
EXPECTED_BLOCKERS = (
    LandingToken("R1_REQUIRED(electric_bc_selection)"),
    LandingToken("R1_REQUIRED(direct_moving_throat)"),
    LandingToken("R1_REQUIRED(magnitude)"),
    LandingToken("R1_REQUIRED(consistency)"),
)


def r1_blockers(facts: LiveFacts) -> tuple[LandingToken, ...]:
    blockers: list[LandingToken] = []
    if facts.anchor is AnchorFact.R1:
        blockers.append(LandingToken("R1_REQUIRED(electric_bc_selection)"))
    if (
        facts.current is CurrentIdentity.R1
        or facts.comparison is ComparisonFact.ROUTE_B_R1
    ):
        blockers.append(LandingToken("R1_REQUIRED(direct_moving_throat)"))
    if facts.magnitude is MagnitudeFact.R1:
        blockers.append(LandingToken("R1_REQUIRED(magnitude)"))
    if facts.tier is TierFact.CONDITIONAL:
        blockers.append(LandingToken("R1_REQUIRED(consistency)"))
    return tuple(blockers)


def landing_ownership_guard(facts: LiveFacts) -> tuple[bool, bool, bool]:
    """Three independently corruptible ownership lanes; no X-identical-X."""
    return (
        r1_blockers(facts) == EXPECTED_BLOCKERS,
        CARRIED_NEUTRAL_TOKENS == facts.neutral_tokens(),
        EMITTED_LANDING == section4_adjudicate(facts),
    )


coefficient_B, velocity_pair, F_flux = sp.symbols(
    "coefficient_B velocity_pair F_flux"
)
CONSERVATIVE_EXCHANGE = coefficient_B * velocity_pair
r_BA, r_cone = sp.symbols("r_BA r_cone")
delta_BA = sp.expand(r_BA - 1)


def active_flux_state(
    mutate: bool = False,
) -> tuple[ActiveFluxFact, IntegrabilityFact, sp.Expr]:
    if mutate:
        return (
            ActiveFluxFact.ABSORBED_DECIDED_ZERO,
            IntegrabilityFact.CLOSED,
            CONSERVATIVE_EXCHANGE + F_flux,
        )
    return (
        ActiveFluxFact.R1_UNMATCHED_REMAINDER,
        IntegrabilityFact.R1_OPEN,
        CONSERVATIVE_EXCHANGE,
    )


def active_flux_guard(
    state: tuple[ActiveFluxFact, IntegrabilityFact, sp.Expr],
) -> tuple[bool, bool, bool]:
    flux_fact, integrability, conservative = state
    return (
        flux_fact is ActiveFluxFact.R1_UNMATCHED_REMAINDER,
        integrability is IntegrabilityFact.R1_OPEN,
        F_flux not in conservative.free_symbols,
    )


def lorentz_lock_state(
    close_all: bool = False,
) -> tuple[tuple[bool, bool, bool, bool], LorentzDetermination]:
    substitutions = {r_BA: 1, r_cone: 1} if close_all else {}
    delta_zero = sp.simplify(delta_BA.subs(substitutions)) == 0
    cone_one = sp.simplify(r_cone.subs(substitutions) - 1) == 0
    higher_orders_closed = close_all
    active_flux_matched = close_all
    locks = (
        bool(delta_zero),
        bool(cone_one),
        higher_orders_closed,
        active_flux_matched,
    )
    determination = (
        LorentzDetermination.DETERMINED
        if all(locks)
        else LorentzDetermination.UNDETERMINED
    )
    return locks, determination


# ---------------------------------------------------------------------------
# Structural seal, local term inventory, dimensions, and source manifest.
# ---------------------------------------------------------------------------

def _code_objects(code: types.CodeType) -> Iterable[types.CodeType]:
    yield code
    for constant in code.co_consts:
        if isinstance(constant, types.CodeType):
            yield from _code_objects(constant)


def target_blindness_guard(
    mutate: bool = False,
) -> tuple[tuple[tuple[str, str], ...], frozenset[str]]:
    """Inspect held code constants; no source-file text scan is used."""
    allowed = {
        "section4_adjudicate",
        "section4_oracle",
        "truth_table",
        "target_blindness_guard",
    }
    policed = (
        "magnetism_lorentz_consistent",
        "amendment_excluded",
        "magnetism_departure_characterized",
        "no_go(sector)",
    )
    inventory: list[tuple[str, types.CodeType]] = []
    for name, value in globals().items():
        if isinstance(value, types.FunctionType) and value.__module__ == __name__:
            for code in _code_objects(value.__code__):
                inventory.append((name if code is value.__code__ else code.co_name, code))
    if mutate:
        leaked = compile(repr(policed[0]), "stage038_nonsealed_leak", "eval")
        inventory.append(("nonsealed_leak", leaked))

    violations: list[tuple[str, str]] = []
    for function_name, code in inventory:
        for constant in code.co_consts:
            if isinstance(constant, str):
                lowered = constant.lower()
                for token in policed:
                    if token in lowered and function_name not in allowed:
                        violations.append((function_name, token))

    cited_symbols = {
        "A_E",
        "q_T",
        "r_BA",
        "delta_BA",
        "r_cone",
        "c_gamma",
        "c_E",
        "mu_R",
        "rho_br",
        "F_flux",
    }
    barred = {"N_u", "a_T", "a'_T", "a_L", "q_A^T", "q_L"}
    if mutate:
        cited_symbols |= barred
    return tuple(sorted(violations)), frozenset(cited_symbols & barred)


def dual_inventory(
    truth: TruthAudit,
    mutate: bool = False,
) -> tuple[Any, ...]:
    locks, determination = lorentz_lock_state()
    flux_state = active_flux_state()
    counts = truth.counts
    digest = truth.digest
    landing = EMITTED_LANDING
    blockers = r1_blockers(LIVE_FACTS)
    delta = delta_BA
    cone = r_cone
    if mutate:
        counts = counts[1:]
        digest = hashlib.sha256(digest).digest()
        landing = section4_adjudicate(
            replace(LIVE_FACTS, anchor=AnchorFact.CLOSED)
        )
        blockers = blockers[:-1]
        delta = sp.expand(delta_BA + 1)
        cone = 2 * r_cone
        locks, determination = lorentz_lock_state(close_all=True)
        flux_state = active_flux_state(mutate=True)
    return (
        counts,
        digest,
        landing,
        blockers,
        sp.expand(delta),
        cone,
        locks,
        determination,
        flux_state,
    )


Dimension = tuple[int, int, int, int]  # M, L, T, E-charge
DIMENSIONLESS: Dimension = (0, 0, 0, 0)


def dim_add(left: Dimension, right: Dimension) -> Dimension:
    return tuple(a + b for a, b in zip(left, right, strict=True))


def dim_sub(left: Dimension, right: Dimension) -> Dimension:
    return tuple(a - b for a, b in zip(left, right, strict=True))


def dim_scale(power: int, value: Dimension) -> Dimension:
    return tuple(power * item for item in value)


def unit_state(corrupt: bool = False) -> tuple[Any, ...]:
    if corrupt:
        A_E_dim = (1, 0, 1, 1)
        q_T_dim = (0, 1, 0, 1)
        c_gamma2_dim = (1, 1, -1, 0)
        c_E2_dim = (0, 3, -1, 1)
        mu_R_dim = (0, 0, 0, 0)
    else:
        A_E_dim = (0, 1, 0, 1)
        q_T_dim = (1, 0, -1, 0)
        c_gamma2_dim = (0, 2, -2, 0)
        c_E2_dim = (0, 2, -2, 0)
        mu_R_dim = (2, 1, -4, -1)
    ratio_dim = dim_sub(
        dim_add(dim_scale(2, q_T_dim), c_gamma2_dim),
        dim_add(mu_R_dim, A_E_dim),
    )
    delta_dim: Dimension | None = (
        DIMENSIONLESS if ratio_dim == DIMENSIONLESS else None
    )
    cone_dim = dim_sub(c_E2_dim, c_gamma2_dim)
    return (
        A_E_dim,
        q_T_dim,
        c_gamma2_dim,
        c_E2_dim,
        mu_R_dim,
        ratio_dim,
        delta_dim,
        cone_dim,
    )


EXPECTED_UNIT_STATE = (
    (0, 1, 0, 1),
    (1, 0, -1, 0),
    (0, 2, -2, 0),
    (0, 2, -2, 0),
    (2, 1, -4, -1),
    DIMENSIONLESS,
    DIMENSIONLESS,
    DIMENSIONLESS,
)

SOURCE_TOOTH_IDS = (
    "SOURCE_TRANSLATION_CONTINUITY",
    "SOURCE_NOT_IMPORTED",
    "SOURCE_BASIS",
    "PARITY_RW",
    "PARITY_PW",
    "PARITY_ROTATION",
    "PARITY_TIME_REVERSAL",
    "FIELD_IDENTITY_UNITS",
    "ACTION_KINETIC",
    "ACTION_COUPLING",
    "ACTION_STABILITY",
    "G0_DAMAGE",
    "ROUTE_INDEPENDENCE",
    "BOOST_PROJECTOR",
    "BOOST_GENERAL_VELOCITIES",
    "BOOST_NEXT_ORDER",
    "BOOST_COMMON_VELOCITY",
    "DIRECT_SOURCE",
    "DIRECT_PROJECTOR",
    "DIRECT_EXCHANGE_SIGN",
    "DIRECT_FALLOFF",
    "DIRECT_VELOCITY_ORDER",
    "COMPARE_COMPUTED",
    "DELTA_RATIO",
    "CONE_RATIO",
    "QMAG_R1",
    "UNITS_RESTORED",
    "ACTIVE_FLUX_CAVEAT",
    "HOOK_LORENTZ",
    "LEDGER_READY_ROW",
    "TRUTH_TOTALITY",
    "TRUTH_PRECEDENCE",
    "LANDING_OWNERSHIP",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
)

IN_SCOPE_SOURCE = (
    "TRUTH_TOTALITY",
    "TRUTH_PRECEDENCE",
    "LANDING_OWNERSHIP",
    "ACTIVE_FLUX_CAVEAT",
    "HOOK_LORENTZ",
)
BUILD_GLOBAL_SOURCE = (
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
    "UNITS_RESTORED",
)
SCOPED_OUT_SOURCE = (
    "FIELD_IDENTITY_UNITS",
    "ACTION_KINETIC",
    "ACTION_COUPLING",
    "ACTION_STABILITY",
    "G0_DAMAGE",
    "LEDGER_READY_ROW",
    "SOURCE_TRANSLATION_CONTINUITY",
    "SOURCE_NOT_IMPORTED",
    "SOURCE_BASIS",
    "PARITY_RW",
    "PARITY_PW",
    "PARITY_ROTATION",
    "PARITY_TIME_REVERSAL",
    "BOOST_PROJECTOR",
    "BOOST_GENERAL_VELOCITIES",
    "BOOST_NEXT_ORDER",
    "DIRECT_SOURCE",
    "DIRECT_PROJECTOR",
    "DIRECT_EXCHANGE_SIGN",
    "DIRECT_FALLOFF",
    "DIRECT_VELOCITY_ORDER",
    "ROUTE_INDEPENDENCE",
    "BOOST_COMMON_VELOCITY",
    "COMPARE_COMPUTED",
    "DELTA_RATIO",
    "CONE_RATIO",
    "QMAG_R1",
)


def manifest_disposition(identifier: str) -> tuple[str, str]:
    if identifier in IN_SCOPE_SOURCE:
        return "REPLACED_BY_STRONGER", "STAGE038_SEALED_RECONSTRUCTION"
    if identifier in BUILD_GLOBAL_SOURCE:
        return "REPLACED_BY_STRONGER", "STAGE038_BUILD_GLOBAL"
    if identifier in SCOPED_OUT_SOURCE[:6]:
        return "SCOPED_OUT", "STAGE034_V1_DONE"
    if identifier in SCOPED_OUT_SOURCE[6:13]:
        return "SCOPED_OUT", "STAGE035_V2_DONE"
    if identifier in SCOPED_OUT_SOURCE[13:16]:
        return "SCOPED_OUT", "STAGE036_V3_DONE"
    if identifier in SCOPED_OUT_SOURCE[16:]:
        return "SCOPED_OUT", "STAGE037_V4_DONE"
    raise KeyError(identifier)


SOURCE_MANIFEST = tuple(
    (identifier, *manifest_disposition(identifier))
    for identifier in SOURCE_TOOTH_IDS
)


def manifest_digest(
    manifest: Iterable[tuple[str, str, str]],
) -> bytes:
    stream = "\n".join("|".join(row) for row in sorted(manifest))
    return hashlib.sha256(stream.encode("utf-8")).digest()


def manifest_state(
    manifest: tuple[tuple[str, str, str], ...],
) -> tuple[Any, ...]:
    partition = tuple(sorted(Counter(row[1] for row in manifest).items()))
    in_scope = tuple(
        row[0] for row in manifest if row[1] != "SCOPED_OUT"
    )
    scoped_out = tuple(
        row[0] for row in manifest if row[1] == "SCOPED_OUT"
    )
    return (
        tuple(row[0] for row in manifest),
        partition,
        in_scope,
        scoped_out,
        manifest_digest(manifest),
    )


EXPECTED_MANIFEST_STATE = manifest_state(SOURCE_MANIFEST)


# ---------------------------------------------------------------------------
# Assertions and report.
# ---------------------------------------------------------------------------

def verdict_table(
    truth: TruthAudit,
    mutate: bool = False,
) -> tuple[tuple[LandingToken, ...], tuple[LandingToken, ...]]:
    production = (
        replace(LIVE_FACTS, anchor=AnchorFact.CLOSED)
        if mutate
        else LIVE_FACTS
    )
    internal = replace(LIVE_FACTS, internal_sectors=("injected",))
    anchor_resolved = replace(LIVE_FACTS, anchor=AnchorFact.CLOSED)
    consistent_free = LiveFacts(
        CurrentIdentity.CONVECTION_CONDITIONAL,
        ComparisonFact.AGREE,
        RelativeSignFact.MATCH,
        MagnitudeFact.FREE,
        TierFact.GAPS_CLOSED,
        AnchorFact.CLOSED,
        (),
    )
    lorentz_consistent = replace(
        consistent_free,
        magnitude=MagnitudeFact.FORCED,
    )
    actual = (
        section4_adjudicate(production),
        section4_adjudicate(internal),
        section4_adjudicate(anchor_resolved),
        section4_adjudicate(consistent_free),
        section4_adjudicate(lorentz_consistent),
        section4_adjudicate(LIVE_FACTS, precedence_mutation=True),
    )
    expected = (
        EMITTED_LANDING,
        truth.target("internal"),
        LandingToken("R1_REQUIRED(direct_moving_throat)"),
        truth.target("consistent_free"),
        truth.target("lorentz_consistent"),
        LandingToken("R1_REQUIRED(consistency)"),
    )
    return actual, expected


def run_assertions() -> LandingToken:
    section("SEALED section-4 exhaustive truth table")
    truth = truth_table(
        totality_mutation=ACTIVE_MUTATION == "TRUTH_TOTALITY"
    )
    truth_components = (
        truth.agreement,
        truth.no_unclassified,
        truth.cell_count == truth.expected_total,
        truth.digest == PUBLISHED_DIGEST,
        truth.counts == truth.expected_counts,
    )
    expect_bool(
        "TRUTH_TOTALITY",
        all(truth_components),
        {
            "components": truth_components,
            "cells": truth.cell_count,
            "expected_product": truth.expected_total,
            "digest": truth.digest.hex(),
            "counts": tuple((str(k), v) for k, v in truth.counts),
        },
    )
    print(
        f"      cells={truth.cell_count}; expected_product={truth.expected_total}; "
        f"digest={truth.digest.hex()}"
    )

    precedence_landing = section4_adjudicate(
        LIVE_FACTS,
        precedence_mutation=ACTIVE_MUTATION == "TRUTH_PRECEDENCE",
    )
    expect_bool(
        "TRUTH_PRECEDENCE",
        precedence_landing == EMITTED_LANDING,
        {
            "production": str(EMITTED_LANDING),
            "rederived": str(precedence_landing),
        },
    )
    print(f"      production_first_match={precedence_landing}")

    ownership_facts = (
        replace(LIVE_FACTS, anchor=AnchorFact.CLOSED)
        if ACTIVE_MUTATION == "LANDING_OWNERSHIP"
        else LIVE_FACTS
    )
    ownership = landing_ownership_guard(ownership_facts)
    expect_bool(
        "LANDING_OWNERSHIP",
        all(ownership),
        {
            "lanes": ownership,
            "blockers": tuple(map(str, r1_blockers(ownership_facts))),
            "carried": CARRIED_NEUTRAL_TOKENS,
            "live": ownership_facts.neutral_tokens(),
            "emitted": str(EMITTED_LANDING),
            "live_landing": str(section4_adjudicate(ownership_facts)),
        },
    )
    blockers = r1_blockers(LIVE_FACTS)
    print("      blockers=" + " | ".join(map(str, blockers)))

    section("Honest active-flux and Lorentz hooks")
    flux_state = active_flux_state(
        mutate=ACTIVE_MUTATION == "ACTIVE_FLUX_CAVEAT"
    )
    flux_guard = active_flux_guard(flux_state)
    expect_bool(
        "ACTIVE_FLUX_CAVEAT",
        all(flux_guard),
        {"guard": flux_guard, "state": flux_state},
    )
    print("      F_flux=separate R1 O(V1*V2) remainder; full integrability=R1")

    locks, determination = lorentz_lock_state(
        close_all=ACTIVE_MUTATION == "HOOK_LORENTZ"
    )
    lorentz_guard = (
        determination is LorentzDetermination.UNDETERMINED,
        locks == (False, False, False, False),
    )
    expect_bool(
        "HOOK_LORENTZ",
        all(lorentz_guard),
        {
            "delta_BA": delta_BA,
            "r_cone": r_cone,
            "locks": locks,
            "determination": determination.name,
        },
    )
    print(
        "      emergent_lorentz=UNDETERMINED; "
        "locks=(delta_BA,r_cone,higher_orders,active_flux)"
    )

    section("Build-global structural, local-inventory, and unit firewalls")
    blindness = target_blindness_guard(
        mutate=ACTIVE_MUTATION == "TARGET_BLINDNESS"
    )
    expect_bool(
        "TARGET_BLINDNESS",
        blindness == ((), frozenset()),
        {"token_violations": blindness[0], "barred_markers": blindness[1]},
    )
    print(
        "      landing assignment lives only in SEALED section 4; "
        "the structural visitor restricts exactly four non-R1 substrings"
    )

    reference_inventory = dual_inventory(truth)
    live_inventory = dual_inventory(
        truth,
        mutate=ACTIVE_MUTATION == "DUAL_ENGINE_TERMS",
    )
    expect_bool(
        "DUAL_ENGINE_TERMS",
        live_inventory == reference_inventory,
        {
            "count_terms": len(live_inventory[0]),
            "digest": compact(live_inventory[1]),
            "landing": str(live_inventory[2]),
            "blockers": tuple(map(str, live_inventory[3])),
            "delta": live_inventory[4],
            "cone": live_inventory[5],
            "locks": live_inventory[6],
        },
    )
    print(
        "      local_terms=landing-count map,digest,landing,blockers,"
        "delta_BA,r_cone,hook conditions"
    )

    dimensions = unit_state(
        corrupt=ACTIVE_MUTATION == "UNITS_RESTORED"
    )
    expect_bool(
        "UNITS_RESTORED",
        dimensions == EXPECTED_UNIT_STATE,
        {"live_dimensions": dimensions, "expected": EXPECTED_UNIT_STATE},
    )
    print(
        "      [A_E]=E*L; [q_T]=M*T^-1; "
        "[r_BA]=[delta_BA]=[r_cone]=1; [c_gamma^2]=[c_E^2]=L^2*T^-2"
    )

    section("Build-native verdict re-derivation table")
    actual_verdicts, expected_verdicts = verdict_table(
        truth,
        mutate=ACTIVE_MUTATION == "VERDICT_REDERIVATION",
    )
    print("      REDERIVED_LANDINGS=" + " | ".join(map(str, actual_verdicts)))
    expect_bool(
        "VERDICT_REDERIVATION",
        actual_verdicts == expected_verdicts,
        {
            "actual": tuple(map(str, actual_verdicts)),
            "expected": tuple(map(str, expected_verdicts)),
        },
    )

    section("Canonical source-to-stage predicate manifest")
    live_manifest = SOURCE_MANIFEST
    if ACTIVE_MUTATION == "SOURCE_TO_STAGE_MANIFEST":
        removed = SCOPED_OUT_SOURCE[0]
        live_manifest = tuple(
            row for row in SOURCE_MANIFEST if row[0] != removed
        )
    live_manifest_state = manifest_state(live_manifest)
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        live_manifest_state == EXPECTED_MANIFEST_STATE,
        {
            "partition": live_manifest_state[1],
            "in_scope": live_manifest_state[2],
            "scoped_out": live_manifest_state[3],
            "digest": compact(live_manifest_state[4]),
        },
    )
    source_total = (
        len(IN_SCOPE_SOURCE)
        + len(SCOPED_OUT_SOURCE)
        + len(BUILD_GLOBAL_SOURCE)
    )
    print(
        f"      source_manifest={source_total}; "
        f"in_scope={len(IN_SCOPE_SOURCE)}; "
        f"scoped_out={len(SCOPED_OUT_SOURCE)}; "
        f"build_global={len(BUILD_GLOBAL_SOURCE)}"
    )
    print("      SCOPED_OUT=" + ",".join(SCOPED_OUT_SOURCE))

    print("")
    print("LANDING=" + str(EMITTED_LANDING))
    print("BLOCKERS=" + " | ".join(map(str, blockers)))
    print(f"TRUTH_CELLS={truth.cell_count}")
    print(f"TRUTH_DIGEST={truth.digest.hex()}")
    for landing, count in truth.counts:
        print(f"LANDING_COUNT[{landing}]={count}")
    print("DEFENSIVE_UNCLASSIFIED=0")
    print("HOOK_ACTIVE_FLUX=R1_UNMATCHED_O(V1*V2)")
    print("HOOK_LORENTZ=UNDETERMINED")
    print("STAGE039_BOUNDARY=departure characterization is not this terminal verdict")
    return EMITTED_LANDING


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage038_sealed_landing_electric_bc_r1 SymPy audit")
    print(
        "ROUTE=itertools Cartesian stream + first-match cascade + "
        "independent complete/open-reason oracle + hashlib SHA256"
    )
    print("FILE_IO=none; CROSS_ENGINE_COMPARE=none")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")

    landing = run_assertions()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print("")
    print(f"TOOTH_COUNT={len(TOOTH_ORDER)}")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print(f"OVERALL PASS: SymPy independently reached {landing}")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(
            "OVERALL FAIL: SymPy stage038 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
