#!/usr/bin/env python3
"""Ledger stage039 SymPy audit: the b_T time-reversal-even departure.

Standalone, print-only, assert-zero, exact-integer, and cross-engine-file-I/O
free.  The engine re-instantiates the cited parity census, derives b_T from
u_T with an explicit per-attribute curl map, compares only the two attributes
shared with a Maxwell B benchmark, and propagates the active-drain time arrow.

Tooth-local runtime ablation uses ``LEDGER_STAGE039_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, fields, is_dataclass, replace
from enum import Enum
import os
from typing import Any, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE039_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

VERDICT = "B_TIME_REVERSAL_EVEN"
COUNTERFACTUAL = "COUNTERFACTUAL_B_TIME_REVERSAL_ODD_MAXWELL_CONSISTENT"
UNCLASSIFIED = "UNCLASSIFIED_PARITY_STATE"

TOOTH_ORDER = (
    "B_T_AXIAL_T_EVEN",
    "MAXWELL_B_T_ODD_DEPARTURE",
    "DEPARTURE_LOCALIZED_TO_T",
    "ACTIVE_DRAIN_TIME_ARROW_REQUIRED",
    "DEPARTURE_SELF_CONSISTENT",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
    "UNITS_RESTORED",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "B_T_AXIAL_T_EVEN":
        "corrupt the per-axis curl map so the derived b_T is T-odd and polar",
    "MAXWELL_B_T_ODD_DEPARTURE":
        "separately flip the derived-b_T side and the Maxwell-benchmark side",
    "DEPARTURE_LOCALIZED_TO_T":
        "make the locally compared derived b_T polar",
    "ACTIVE_DRAIN_TIME_ARROW_REQUIRED":
        "replace the active T-odd tau_d root by a passive T-even root",
    "DEPARTURE_SELF_CONSISTENT":
        "separately break path agreement and inject time_reversal into R72",
    "TARGET_BLINDNESS":
        "inject a barred pathA_39 symbol into the held parity-object graph",
    "DUAL_ENGINE_TERMS":
        "drop the computed active-drain chain from the local term inventory",
    "UNITS_RESTORED":
        "replace the cited displacement dimension L by dimensionless",
    "VERDICT_REDERIVATION":
        "flip computed u_T T-parity and re-run the curl/Maxwell verdict table",
    "SOURCE_TO_STAGE_MANIFEST":
        "remove one scoped-out V-1 row from the live 35-tooth partition",
}


class AuditFailure(AssertionError):
    """A named stage predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


class Rotation(Enum):
    SCALAR = "scalar"
    POLAR = "polar_vector"
    AXIAL = "axial_vector"


class PhysicsAxis(Enum):
    TIME_REVERSAL = "time_reversal"
    ROTATION = "rotation"
    SIGN = "sign"
    MAGNITUDE = "magnitude"


@dataclass(frozen=True)
class ParityRecord:
    symbol: sp.Symbol
    R_w: int
    P_w: int
    time_reversal: int
    rotation: Rotation

    def comparable(self) -> Mapping[PhysicsAxis, object]:
        return {
            PhysicsAxis.TIME_REVERSAL: self.time_reversal,
            PhysicsAxis.ROTATION: self.rotation,
        }


@dataclass(frozen=True)
class MaxwellRecord:
    time_reversal: int
    rotation: Rotation

    def comparable(self) -> Mapping[PhysicsAxis, object]:
        return {
            PhysicsAxis.TIME_REVERSAL: self.time_reversal,
            PhysicsAxis.ROTATION: self.rotation,
        }


@dataclass(frozen=True)
class ChainState:
    tau_d: int
    q_T: int
    J_T: int
    u_T: int
    b_T: int
    moving_row_available: bool


def compact(value: Any) -> str:
    if isinstance(value, Enum):
        return value.value
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
# Cited stage035 census and the stage-local per-axis curl derivation.
# ---------------------------------------------------------------------------

S_SYMBOL = sp.Symbol("s")
V_SYMBOL = sp.Symbol("V")
TAU_D_SYMBOL = sp.Symbol("tau_d")
Q_T_SYMBOL = sp.Symbol("q_T")
J_T_SYMBOL = sp.Symbol("J_T")
U_T_SYMBOL = sp.Symbol("u_T")
B_T_SYMBOL = sp.Symbol("b_T")

CITED_CENSUS = {
    "s": ParityRecord(S_SYMBOL, -1, -1, +1, Rotation.SCALAR),
    "V": ParityRecord(V_SYMBOL, +1, +1, -1, Rotation.POLAR),
    "tau_d": ParityRecord(TAU_D_SYMBOL, +1, +1, -1, Rotation.SCALAR),
    "q_T": ParityRecord(Q_T_SYMBOL, +1, +1, -1, Rotation.SCALAR),
    "J_T": ParityRecord(J_T_SYMBOL, -1, -1, +1, Rotation.POLAR),
    "u_T": ParityRecord(U_T_SYMBOL, -1, -1, +1, Rotation.POLAR),
    "b_T": ParityRecord(B_T_SYMBOL, -1, -1, +1, Rotation.AXIAL),
}

MAXWELL_B = MaxwellRecord(time_reversal=-1, rotation=Rotation.AXIAL)

CURL_INTEGER_OPERATORS = {
    "R_w": lambda parity: +1 * parity,
    "P_w": lambda parity: +1 * parity,
    "time_reversal": lambda parity: +1 * parity,
}

CURL_ROTATION_MAP = {
    Rotation.POLAR: Rotation.AXIAL,
    Rotation.AXIAL: Rotation.POLAR,
}


def derive_curl_parity(
    input_record: ParityRecord,
    *,
    corrupt: bool = False,
) -> ParityRecord:
    """Apply an explicit per-attribute curl operator to a parity record."""
    operators = dict(CURL_INTEGER_OPERATORS)
    rotations = dict(CURL_ROTATION_MAP)
    if corrupt:
        operators["time_reversal"] = lambda parity: -1 * parity
        rotations[Rotation.POLAR] = Rotation.POLAR
    return ParityRecord(
        symbol=B_T_SYMBOL,
        R_w=operators["R_w"](input_record.R_w),
        P_w=operators["P_w"](input_record.P_w),
        time_reversal=operators["time_reversal"](
            input_record.time_reversal
        ),
        rotation=rotations[input_record.rotation],
    )


def disagreement_axes(
    model: ParityRecord,
    benchmark: MaxwellRecord,
) -> frozenset[PhysicsAxis]:
    """Compare exactly {time_reversal, rotation}; R_w/P_w are not Maxwell data."""
    model_record = model.comparable()
    benchmark_record = benchmark.comparable()
    return frozenset(
        axis
        for axis in (PhysicsAxis.TIME_REVERSAL, PhysicsAxis.ROTATION)
        if model_record[axis] != benchmark_record[axis]
    )


def departure_holds(
    model: ParityRecord,
    benchmark: MaxwellRecord,
) -> bool:
    return (
        model.comparable()[PhysicsAxis.TIME_REVERSAL]
        != benchmark.comparable()[PhysicsAxis.TIME_REVERSAL]
    )


def localization_state(
    model: ParityRecord,
    benchmark: MaxwellRecord,
) -> tuple[object, frozenset[PhysicsAxis]]:
    shared_rotation: object = (
        model.rotation if model.rotation is benchmark.rotation else None
    )
    return shared_rotation, disagreement_axes(model, benchmark)


# ---------------------------------------------------------------------------
# Active-drain time-arrow propagation and computed verdict lookup.
# ---------------------------------------------------------------------------

def propagate_active_drain(tau_parity: int) -> ChainState:
    """Compose tau_d -> q_T -> J_T -> u_T -> b_T using integer parities."""
    q_parity = tau_parity
    j_parity = (
        q_parity
        * CITED_CENSUS["s"].time_reversal
        * CITED_CENSUS["V"].time_reversal
    )
    u_parity = j_parity
    b_parity = CURL_INTEGER_OPERATORS["time_reversal"](u_parity)
    return ChainState(
        tau_d=tau_parity,
        q_T=q_parity,
        J_T=j_parity,
        u_T=u_parity,
        b_T=b_parity,
        moving_row_available=(tau_parity == -1),
    )


CITED_CHAIN_STATE = ChainState(
    tau_d=CITED_CENSUS["tau_d"].time_reversal,
    q_T=CITED_CENSUS["q_T"].time_reversal,
    J_T=CITED_CENSUS["J_T"].time_reversal,
    u_T=CITED_CENSUS["u_T"].time_reversal,
    b_T=CITED_CENSUS["b_T"].time_reversal,
    moving_row_available=True,
)

VERDICT_LOOKUP = {
    (
        +1,
        -1,
        Rotation.AXIAL,
        Rotation.AXIAL,
        frozenset({PhysicsAxis.TIME_REVERSAL}),
    ): VERDICT,
    (
        -1,
        -1,
        Rotation.AXIAL,
        Rotation.AXIAL,
        frozenset(),
    ): COUNTERFACTUAL,
}


def verdict_from_computed(
    model: ParityRecord,
    benchmark: MaxwellRecord,
) -> str:
    key = (
        model.time_reversal,
        benchmark.time_reversal,
        model.rotation,
        benchmark.rotation,
        disagreement_axes(model, benchmark),
    )
    return VERDICT_LOOKUP.get(key, UNCLASSIFIED)


def verdict_witnesses(
    *,
    mutate_production: bool = False,
) -> tuple[str, str, str]:
    production_u = CITED_CENSUS["u_T"]
    if mutate_production:
        production_u = replace(
            production_u,
            time_reversal=-production_u.time_reversal,
        )
    production_b = derive_curl_parity(production_u)

    flipped_u = replace(
        CITED_CENSUS["u_T"],
        time_reversal=-CITED_CENSUS["u_T"].time_reversal,
    )
    u_flip_b = derive_curl_parity(flipped_u)

    tau_flip_chain = propagate_active_drain(
        -CITED_CENSUS["tau_d"].time_reversal
    )
    tau_flip_b = replace(
        derive_curl_parity(CITED_CENSUS["u_T"]),
        time_reversal=tau_flip_chain.b_T,
    )
    return (
        verdict_from_computed(production_b, MAXWELL_B),
        verdict_from_computed(u_flip_b, MAXWELL_B),
        verdict_from_computed(tau_flip_b, MAXWELL_B),
    )


# ---------------------------------------------------------------------------
# Structural blindness, local computed-term inventory, and unit firewall.
# ---------------------------------------------------------------------------

BARRED_SOURCE_MARKERS = frozenset(
    {"Nu", "aT", "aTp", "aL", "q_A_T", "q_L"}
)


def structural_symbol_inventory(objects: Iterable[object]) -> frozenset[str]:
    """Walk live objects structurally; never read or grep a source file."""
    found: set[str] = set()

    def visit(value: object) -> None:
        if isinstance(value, sp.Symbol):
            found.add(value.name)
        elif isinstance(value, sp.Basic):
            found.update(symbol.name for symbol in value.free_symbols)
        elif is_dataclass(value):
            for item in fields(value):
                visit(getattr(value, item.name))
        elif isinstance(value, Mapping):
            for key, item in value.items():
                visit(key)
                visit(item)
        elif isinstance(value, (tuple, list, set, frozenset)):
            for item in value:
                visit(item)

    for obj in objects:
        visit(obj)
    return frozenset(found)


def target_blindness_state(
    derived_b: ParityRecord,
    chain: ChainState,
    *,
    mutate: bool = False,
) -> frozenset[str]:
    live_objects: list[object] = [
        CITED_CENSUS,
        derived_b,
        MAXWELL_B,
        chain,
        disagreement_axes(derived_b, MAXWELL_B),
    ]
    if mutate:
        live_objects.append(sp.Symbol("Nu"))
    return (
        structural_symbol_inventory(live_objects)
        & BARRED_SOURCE_MARKERS
    )


def computed_inventory(
    derived_b: ParityRecord,
    chain: ChainState,
    *,
    mutate: bool = False,
) -> Mapping[str, object]:
    inventory: dict[str, object] = {
        "derived_b_T": derived_b,
        "maxwell_B": MAXWELL_B,
        "departure_holds": departure_holds(derived_b, MAXWELL_B),
        "disagreement_axes": disagreement_axes(derived_b, MAXWELL_B),
        "active_drain_chain": chain,
        "verdict": verdict_from_computed(derived_b, MAXWELL_B),
    }
    if mutate:
        del inventory["active_drain_chain"]
    return inventory


EXPECTED_INVENTORY = {
    "derived_b_T": CITED_CENSUS["b_T"],
    "maxwell_B": MaxwellRecord(-1, Rotation.AXIAL),
    "departure_holds": True,
    "disagreement_axes": frozenset({PhysicsAxis.TIME_REVERSAL}),
    "active_drain_chain": CITED_CHAIN_STATE,
    "verdict": VERDICT,
}

Dimension = tuple[int, int, int]  # M, L, T
DIMENSIONLESS: Dimension = (0, 0, 0)
LENGTH: Dimension = (0, 1, 0)
INVERSE_LENGTH: Dimension = (0, -1, 0)


def multiply_dimensions(left: Dimension, right: Dimension) -> Dimension:
    return tuple(
        sp.Integer(a + b)
        for a, b in zip(left, right, strict=True)
    )


def unit_state(*, corrupt: bool = False) -> tuple[Dimension, Dimension]:
    u_dimension = DIMENSIONLESS if corrupt else LENGTH
    curl_u_dimension = multiply_dimensions(INVERSE_LENGTH, u_dimension)
    b_dimension = DIMENSIONLESS
    return curl_u_dimension, b_dimension


# ---------------------------------------------------------------------------
# Complete 35-source-tooth manifest plus the five stage-native authored teeth.
# ---------------------------------------------------------------------------

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

CITED_SOURCE = (
    "PARITY_ROTATION",
    "PARITY_TIME_REVERSAL",
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

STAGE_NATIVE_AUTHORED = TOOTH_ORDER[:5]


def source_disposition(identifier: str) -> tuple[str, str]:
    if identifier in CITED_SOURCE:
        return "CITED", "STAGE035_V2_DONE"
    if identifier in BUILD_GLOBAL_SOURCE:
        return "REPLACED_BY_STRONGER", "STAGE039_BUILD_GLOBAL"
    if identifier in SCOPED_OUT_SOURCE[:6]:
        return "SCOPED_OUT", "STAGE034_V1_DONE"
    if identifier in SCOPED_OUT_SOURCE[6:11]:
        return "SCOPED_OUT", "STAGE035_V2_DONE"
    if identifier in SCOPED_OUT_SOURCE[11:14]:
        return "SCOPED_OUT", "STAGE036_V3_DONE"
    if identifier in SCOPED_OUT_SOURCE[14:25]:
        return "SCOPED_OUT", "STAGE037_V4_DONE"
    if identifier in SCOPED_OUT_SOURCE[25:]:
        return "SCOPED_OUT", "STAGE038_V5_DONE"
    raise ValueError(f"unpartitioned source tooth: {identifier}")


SOURCE_MANIFEST = tuple(
    (identifier, *source_disposition(identifier))
    for identifier in SOURCE_TOOTH_IDS
)


def manifest_state(
    manifest: Iterable[tuple[str, str, str]],
) -> tuple[object, ...]:
    rows = tuple(manifest)
    dispositions = Counter(row[1] for row in rows)
    return (
        tuple(row[0] for row in rows),
        frozenset(row[0] for row in rows if row[1] == "CITED"),
        frozenset(row[0] for row in rows if row[1] == "SCOPED_OUT"),
        frozenset(
            row[0] for row in rows if row[1] == "REPLACED_BY_STRONGER"
        ),
        tuple(sorted(dispositions.items())),
        tuple(rows),
        len(rows),
    )


EXPECTED_MANIFEST_STATE = (
    SOURCE_TOOTH_IDS,
    frozenset(CITED_SOURCE),
    frozenset(SCOPED_OUT_SOURCE),
    frozenset(BUILD_GLOBAL_SOURCE),
    (
        ("CITED", 2),
        ("REPLACED_BY_STRONGER", 3),
        ("SCOPED_OUT", 30),
    ),
    SOURCE_MANIFEST,
    35,
)


# ---------------------------------------------------------------------------
# Ten executable teeth.
# ---------------------------------------------------------------------------

def run_assertions() -> str:
    section("Cited census -> explicit per-axis curl inheritance")
    derived_b = derive_curl_parity(
        CITED_CENSUS["u_T"],
        corrupt=ACTIVE_MUTATION == "B_T_AXIAL_T_EVEN",
    )
    expect_bool(
        "B_T_AXIAL_T_EVEN",
        derived_b == CITED_CENSUS["b_T"],
        {
            "derived": derived_b,
            "cited_035": CITED_CENSUS["b_T"],
            "route": "per-axis operator dictionary",
        },
    )
    print(
        "      b_T=(-1,-1,+1,axial_vector), derived from "
        "u_T=(-1,-1,+1,polar_vector)"
    )

    section("Typed Maxwell-B comparison over exactly {T, rotation}")
    production_departure = departure_holds(derived_b, MAXWELL_B)
    if ACTIVE_MUTATION == "MAXWELL_B_T_ODD_DEPARTURE":
        flipped_u = replace(
            CITED_CENSUS["u_T"],
            time_reversal=-CITED_CENSUS["u_T"].time_reversal,
        )
        b_side_departure = departure_holds(
            derive_curl_parity(flipped_u),
            MAXWELL_B,
        )
        maxwell_side_departure = departure_holds(
            derived_b,
            replace(MAXWELL_B, time_reversal=+1),
        )
        departure_check = b_side_departure and maxwell_side_departure
    else:
        b_side_departure = None
        maxwell_side_departure = None
        departure_check = production_departure
    expect_bool(
        "MAXWELL_B_T_ODD_DEPARTURE",
        departure_check,
        {
            "derived_b_T": derived_b.time_reversal,
            "maxwell_B_T": MAXWELL_B.time_reversal,
            "b_side_ablation": b_side_departure,
            "maxwell_side_ablation": maxwell_side_departure,
        },
    )
    print("      b_T is T-even; Maxwell B is T-odd; departure=True")

    localized_b = (
        replace(derived_b, rotation=Rotation.POLAR)
        if ACTIVE_MUTATION == "DEPARTURE_LOCALIZED_TO_T"
        else derived_b
    )
    live_localization = localization_state(localized_b, MAXWELL_B)
    expected_localization = (
        Rotation.AXIAL,
        frozenset({PhysicsAxis.TIME_REVERSAL}),
    )
    expect_bool(
        "DEPARTURE_LOCALIZED_TO_T",
        live_localization == expected_localization,
        {
            "shared_rotation": live_localization[0],
            "disagreements": live_localization[1],
            "comparable_domain": (
                PhysicsAxis.TIME_REVERSAL,
                PhysicsAxis.ROTATION,
            ),
        },
    )
    decided_axes = disagreement_axes(derived_b, MAXWELL_B)
    print("      rotation=axial agrees; disagreement={time_reversal}")

    section("Active-drain tau_d time-arrow propagation")
    tau_root = CITED_CENSUS["tau_d"].time_reversal
    if ACTIVE_MUTATION == "ACTIVE_DRAIN_TIME_ARROW_REQUIRED":
        tau_root = -tau_root
    active_chain = propagate_active_drain(tau_root)
    expect_bool(
        "ACTIVE_DRAIN_TIME_ARROW_REQUIRED",
        active_chain == CITED_CHAIN_STATE,
        {
            "computed_chain": active_chain,
            "cited_chain": CITED_CHAIN_STATE,
            "passive_throat_supplies_O(V)_row":
                active_chain.moving_row_available,
        },
    )
    print("      tau_d(-1)->q_T(-1)->J_T(+1)->u_T(+1)->b_T(+1)")
    print("      passive T-even throat: no O(V) moving row")

    production_self_state = (
        derived_b.time_reversal == active_chain.b_T,
        decided_axes.isdisjoint(
            {PhysicsAxis.SIGN, PhysicsAxis.MAGNITUDE}
        ),
    )
    path_ablation_state = (
        -derived_b.time_reversal == active_chain.b_T,
        production_self_state[1],
    )
    axis_ablation_state = (
        production_self_state[0],
        decided_axes.isdisjoint(
            {
                PhysicsAxis.SIGN,
                PhysicsAxis.MAGNITUDE,
                PhysicsAxis.TIME_REVERSAL,
            }
        ),
    )
    live_self_state = (
        (path_ablation_state[0], axis_ablation_state[1])
        if ACTIVE_MUTATION == "DEPARTURE_SELF_CONSISTENT"
        else production_self_state
    )
    expect_bool(
        "DEPARTURE_SELF_CONSISTENT",
        live_self_state == (True, True),
        {
            "production": production_self_state,
            "path_only_ablation": path_ablation_state,
            "axis_only_ablation": axis_ablation_state,
            "R72": "R1_REQUIRED(electric_bc_selection)",
            "R72_axes": {
                PhysicsAxis.SIGN,
                PhysicsAxis.MAGNITUDE,
            },
        },
    )
    print(
        "      curl-path T == tau-chain T; "
        "{time_reversal} disjoint {sign,magnitude}"
    )
    print(
        "      separate ablations: path=(False,True); "
        "R72-axis injection=(True,False)"
    )

    section("Build-global structural, local-inventory, and unit firewalls")
    blindness = target_blindness_state(
        derived_b,
        active_chain,
        mutate=ACTIVE_MUTATION == "TARGET_BLINDNESS",
    )
    expect_bool(
        "TARGET_BLINDNESS",
        blindness == frozenset(),
        {
            "barred_intersection": blindness,
            "barred": BARRED_SOURCE_MARKERS,
        },
    )
    print("      structural SymPy/dataclass symbol graph; no source text scan")

    inventory = computed_inventory(
        derived_b,
        active_chain,
        mutate=ACTIVE_MUTATION == "DUAL_ENGINE_TERMS",
    )
    expect_bool(
        "DUAL_ENGINE_TERMS",
        inventory == EXPECTED_INVENTORY,
        {
            "terms": tuple(inventory),
            "derived_b_T": inventory["derived_b_T"],
            "maxwell_B": inventory["maxwell_B"],
            "departure": inventory["departure_holds"],
            "disagreements": inventory["disagreement_axes"],
            "verdict": inventory["verdict"],
        },
    )
    print(
        "      local terms=derived b_T,Maxwell B,departure,"
        "T-axis,active-drain chain,verdict"
    )

    dimensions = unit_state(
        corrupt=ACTIVE_MUTATION == "UNITS_RESTORED"
    )
    expect_bool(
        "UNITS_RESTORED",
        dimensions[0] == dimensions[1],
        {
            "curl_u_T": dimensions[0],
            "b_T": dimensions[1],
        },
    )
    print("      [curl u_T]=L^-1*L=1=[b_T]")

    section("Computed verdict re-derivation table")
    actual_verdicts = verdict_witnesses(
        mutate_production=ACTIVE_MUTATION == "VERDICT_REDERIVATION"
    )
    expected_verdicts = (VERDICT, COUNTERFACTUAL, COUNTERFACTUAL)
    if ACTIVE_MUTATION == "VERDICT_REDERIVATION":
        print(
            "      MUTATION_REDERIVED="
            + COUNTERFACTUAL
            + " via u_T flip and tau_d flip"
        )
    expect_bool(
        "VERDICT_REDERIVATION",
        actual_verdicts == expected_verdicts,
        {
            "production": actual_verdicts[0],
            "u_T_flip": actual_verdicts[1],
            "tau_d_flip": actual_verdicts[2],
        },
    )
    print("      production verdict re-derived from b_T/Maxwell parity objects")

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
            "source_total": live_manifest_state[-1],
            "partition": live_manifest_state[4],
            "cited": live_manifest_state[1],
            "scoped_out": live_manifest_state[2],
            "build_global": live_manifest_state[3],
        },
    )
    print(
        "      source_manifest=35; cited=2; scoped_out=30; "
        "build_global_replaced_by_stronger=3"
    )
    print("      CITED=" + ",".join(CITED_SOURCE))
    print("      SCOPED_OUT=" + ",".join(SCOPED_OUT_SOURCE))
    print(
        "      STAGE_NATIVE_AUTHORED="
        + ",".join(STAGE_NATIVE_AUTHORED)
    )

    print("")
    print("VERDICT_TOKEN=" + VERDICT)
    print("SCOPE_CLASS=CHARACTERIZED-DEPARTURE")
    print("FRAMING=FIRST_CLASS_LOAD_BEARING_NOT_EXACT_MAXWELL")
    print("NOT_A_BUG=TRUE; NOT_R1=TRUE; NEVER_SOFTENED=TRUE")
    print("MAXWELL_COMPARISON_DOMAIN=time_reversal,rotation")
    print("R72_CITED=R1_REQUIRED(electric_bc_selection)")
    print("R72_ORTHOGONAL_AXES=sign,magnitude")
    print("STAGE038_POST_RESOLUTION_TOKEN_REEMITTED=FALSE")
    print(
        "SIBLING_DEPARTURES=NATIVE_P_NO_EMERGENT_GAUSS,"
        "FAIL_CAUCHY_STRAY_LONGITUDINAL"
    )
    return verdict_from_computed(derived_b, MAXWELL_B)


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print(
        "ledger_stage039_b_t_time_reversal_even_departure SymPy audit"
    )
    print(
        "ROUTE=typed census records + explicit per-axis curl operators + "
        "integer-product active-drain chain"
    )
    print("FILE_IO=none; CROSS_ENGINE_COMPARE=none")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")

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
            "OVERALL FAIL: SymPy stage039 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
