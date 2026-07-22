#!/usr/bin/env python3
"""Ledger stage032 SymPy audit: BC ensembles and the sealed R1 landing.

Standalone, print-only, assert-zero, float-free, and file-I/O-free.  Stage031
is consumed as the symbolic bundle ``{m, m_gg=B_eff*z_g^2/D,
det(m)=z_g^2/D, S_gg, A=m_gg*C, NONZERO_HA_REQUIRES_CORE_HOLDER,
1/R^2}``; none of that mechanism is re-derived here.  Stage030's D and
``D*=7/4`` are consumed transitively through stage031.  Tooth-local ablation
uses ``LEDGER_STAGE032_MUTATION``.
"""

from __future__ import annotations

from collections import Counter, OrderedDict
from dataclasses import dataclass, replace
from enum import Enum
import hashlib
import itertools
import os
from typing import Any, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE032_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

VERDICT = "R1_REQUIRED(bc_selection)"
MAGNITUDE_COBLOCKER = "R1_REQUIRED(magnitude)"
COMMITTED_DIGEST = "7627417ace0f99a17187a90efa2523ca9b68df8b64f7960d38be2dc6f553ac49"
COMMITTED_MANIFEST_DIGEST = "e2cfd11b41cdd3d95111f25215b16e917b1ce0a9623929619338a1a83fa81014"
EXPECTED_TOTAL = 23040
EXPECTED_COUNTS = {
    "MECHANISM_FALSIFIED(wrong_range)": 1008,
    "MECHANISM_FALSIFIED(wrong_sign)": 504,
    "NO_GO(sector)": 11520,
    "R1_REQUIRED(bc_selection)": 1152,
    "R1_REQUIRED(magnitude)": 252,
    "R1_REQUIRED(mixed_bc_parameters)": 1152,
    "R1_REQUIRED(sign_and_magnitude)": 2856,
    "R1_REQUIRED(subleading)": 504,
    "R1_REQUIRED(variant_selection)": 3840,
    "SIGN_EARNED": 252,
}
EXPECTED_MANIFEST_COUNTS = {
    "PRESERVED": 22,
    "REPLACED_BY_STRONGER": 12,
    "SCOPED_OUT": 10,
}


class AuditFailure(AssertionError):
    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


TOOTH_ORDER = (
    "PASS_STAGE031_CONSUMED_INPUTS",
    "PASS_STAGE030_TRANSITIVE_D_WITNESS",
    "PASS_SCOPE_FIREWALL",
    "PASS_NONZERO_HA_REQUIRES_CORE_HOLDER",
    "PASS_A_V_FORMULA",
    "PASS_A_J_FORMULA",
    "PASS_A_M_FORMULA",
    "PASS_A_MIXED_FORMULA",
    "PASS_WEAK_SIGNS_GENERAL",
    "PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS",
    "PASS_DEGENERATE_Z_G_ZERO",
    "PASS_A_M_INDEFINITE",
    "PASS_MIXED_ENDPOINTS",
    "PASS_MIXED_INTERIOR_ZERO",
    "PASS_MIXED_REGIME_Q_GT_G",
    "PASS_MIXED_REGIME_Q_EQ_G",
    "PASS_MIXED_REGIME_Q_LT_G",
    "PASS_MIXED_FULL_DOMAIN_SPANS",
    "PASS_DECIDED_CONDITIONAL_TYPING",
    "PASS_AMEND_REPLACE",
    "PASS_AMEND_ADD",
    "PASS_ZERO_LEDGER_13",
    "PASS_S_HOLD_SCOPE",
    "PASS_INTERNAL_INCONSISTENCY_NONE",
    "PASS_BC_ACTUAL_CLASSIFIER",
    "PASS_BC_FREE_CONTROL",
    "PASS_BC_VALUE_CONTROL",
    "PASS_BC_MONOPOLE_CONTROL",
    "PASS_BC_MIXED_CONTROL",
    "PASS_ADMISSIBLE_CLASSES",
    "PASS_REPLACE_TOTALS",
    "PASS_ADD_TOTALS",
    "PASS_VARIANT_UNRESOLVED",
    "PASS_NO_DOUBLE_COUNT",
    "PASS_OUTCOME_NOT_INVARIANT",
    "PASS_MAGNITUDE_FREE_FACTOR",
    "PASS_QMAG_DENSITY_HOOK",
    "PASS_QMAG_MONOPOLE_HOOK",
    "PASS_QMAG_MODULUS_HOOK",
    "PASS_QMAG_CLOSE_RANGE_HOOK",
    "PASS_MAGNITUDE_COBLOCKER",
    "PASS_RANGE_SIGN_FLIP",
    "PASS_RANGE_ZERO_TOUCH",
    "PASS_RANGE_SUBDOMINANT",
    "PASS_UNITS_A",
    "PASS_UNITS_U",
    "PASS_UNITS_F",
    "PASS_TYPED_NEUTRAL_FACTS",
    "PASS_VERDICT_TOTALITY",
    "PASS_VERDICT_PRECEDENCE",
    "PASS_TRUTH_TABLE_TOTAL",
    "PASS_TRUTH_TABLE_COUNTS",
    "PASS_TRUTH_TABLE_DIGEST",
    "PASS_LANDING_OWNERSHIP",
    "PASS_PRODUCTION_LANDING",
    "PASS_TARGET_BLINDNESS",
    "PASS_SOURCE_PREDICATE_MANIFEST",
)
MUTATION_ORDER = TOOTH_ORDER + ("MANIFEST_MISDISPOSITION",)

ABLATION_DESCRIPTIONS = {
    "PASS_STAGE031_CONSUMED_INPUTS": "remove consumed S_gg from the stage031 citation bundle",
    "PASS_STAGE030_TRANSITIVE_D_WITNESS": "D*=7/4 -> 5/4",
    "PASS_SCOPE_FIREWALL": "inject response-matrix re-derivation into stage032 ownership",
    "PASS_NONZERO_HA_REQUIRES_CORE_HOLDER": "drop the named stage031 core-holder fact",
    "PASS_A_V_FORMULA": "A_V gains an extra m_gg*phi^2 term",
    "PASS_A_J_FORMULA": "A_J sign is flipped",
    "PASS_A_M_FORMULA": "omit the held-h g^2 subtraction",
    "PASS_A_MIXED_FORMULA": "omit the mixed q*g term",
    "PASS_WEAK_SIGNS_GENERAL": "assert the opposite weak sign certificate",
    "PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS": "set the nonzero phi witness to zero",
    "PASS_DEGENERATE_Z_G_ZERO": "replace the degenerate z_g=0 witness by z_g=1",
    "PASS_A_M_INDEFINITE": "replace the q<g witness by q>g",
    "PASS_MIXED_ENDPOINTS": "evaluate the second endpoint at lambda=0",
    "PASS_MIXED_INTERIOR_ZERO": "shift lambda* by 1/2",
    "PASS_MIXED_REGIME_Q_GT_G": "remove the positive q>g sample",
    "PASS_MIXED_REGIME_Q_EQ_G": "move the q=g zero from lambda=0",
    "PASS_MIXED_REGIME_Q_LT_G": "use a positive algebraic root for q<g",
    "PASS_MIXED_FULL_DOMAIN_SPANS": "claim every fixed q/g spans all signs",
    "PASS_DECIDED_CONDITIONAL_TYPING": "mis-type class selection as DECIDED",
    "PASS_AMEND_REPLACE": "give REPLACE an illegal new row",
    "PASS_AMEND_ADD": "give ADD a second new row",
    "PASS_ZERO_LEDGER_13": "drop one preserved zero row",
    "PASS_S_HOLD_SCOPE": "extend S_hold to h",
    "PASS_INTERNAL_INCONSISTENCY_NONE": "inject a Q-AMEND inconsistency sector",
    "PASS_BC_ACTUAL_CLASSIFIER": "remove the actual missing-closure premise",
    "PASS_BC_FREE_CONTROL": "turn the free control into imposed value",
    "PASS_BC_VALUE_CONTROL": "remove the imposed value premise",
    "PASS_BC_MONOPOLE_CONTROL": "remove the imposed conormal premise",
    "PASS_BC_MIXED_CONTROL": "remove the imposed mixed relation",
    "PASS_ADMISSIBLE_CLASSES": "drop MIXED from the admissible class set",
    "PASS_REPLACE_TOTALS": "add g^2 to the REPLACE monopole total",
    "PASS_ADD_TOTALS": "collapse the ADD monopole total to q^2",
    "PASS_VARIANT_UNRESOLVED": "select the replace realization",
    "PASS_NO_DOUBLE_COUNT": "double-count held-h work",
    "PASS_OUTCOME_NOT_INVARIANT": "force all class/variant outcomes positive",
    "PASS_MAGNITUDE_FREE_FACTOR": "hide c_xi from the magnitude census",
    "PASS_QMAG_DENSITY_HOOK": "claim a local density prediction",
    "PASS_QMAG_MONOPOLE_HOOK": "claim a decided radial monopole",
    "PASS_QMAG_MODULUS_HOOK": "claim universal quantization",
    "PASS_QMAG_CLOSE_RANGE_HOOK": "claim the close-range hook is decided",
    "PASS_MAGNITUDE_COBLOCKER": "replace the magnitude co-blocker",
    "PASS_RANGE_SIGN_FLIP": "classify a sign flip as constant",
    "PASS_RANGE_ZERO_TOUCH": "classify a zero touch as constant",
    "PASS_RANGE_SUBDOMINANT": "erase the derived subdominant bound",
    "PASS_UNITS_A": "[A]=E L -> E",
    "PASS_UNITS_U": "[U]=E -> E L",
    "PASS_UNITS_F": "[F]=E/L -> E",
    "PASS_TYPED_NEUTRAL_FACTS": "inject a landing token into upstream facts",
    "PASS_VERDICT_TOTALITY": "mark one exhaustive row unclassified",
    "PASS_VERDICT_PRECEDENCE": "route the class gap before unconditional invariance",
    "PASS_TRUTH_TABLE_TOTAL": "increment the exhaustive cell total",
    "PASS_TRUTH_TABLE_COUNTS": "increment the bc-selection grouped count",
    "PASS_TRUTH_TABLE_DIGEST": "corrupt the first serialized exhaustive row",
    "PASS_LANDING_OWNERSHIP": "flip upstream Q-BC to FIXED_SOURCE and freshly re-land",
    "PASS_PRODUCTION_LANDING": "flip the production Q-BC before adjudication",
    "PASS_TARGET_BLINDNESS": "inject SIGN_EARNED into an upstream token",
    "PASS_SOURCE_PREDICATE_MANIFEST": "drop one source claim from the partition",
    "MANIFEST_MISDISPOSITION": (
        "flip AMEND_REPLACE from PRESERVED to SCOPED_OUT while retaining all 44 identifiers"
    ),
}


def primitive(predicate: str, canonical: Any, mutated: Any) -> Any:
    return mutated if ACTIVE_MUTATION == predicate else canonical


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
    if isinstance(value, (str, type(None), Enum)):
        return
    if isinstance(value, bool):
        value = sp.Integer(1) if value else sp.Integer(0)
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
    expect_zero(name, sp.Integer(0) if bool(condition) else sp.Integer(1), evidence)


def heading(text: str) -> None:
    print("")
    print(text)
    print("-" * len(text))


# ---------------------------------------------------------------------------
# Consumed symbolic inputs and stage ownership firewall.
# ---------------------------------------------------------------------------

CONSUMED_STAGE031 = (
    "m", "m_gg=B_eff*z_g^2/D", "det(m)=z_g^2/D", "S_gg",
    "A=m_gg*C", "NONZERO_HA_REQUIRES_CORE_HOLDER", "s1*s2/(4*pi*R^2)",
)
TRANSITIVE_STAGE030 = ("D", "D*=7/4")
OWNED_STAGE032 = frozenset({
    "A_V", "A_J", "A_M", "A_MIXED", "Q_AMEND_consistency",
    "Q_BC_typed_facts", "Q_COMBINE", "Q_MAG", "section4_landing",
})
FORBIDDEN_REDERIVATIONS = frozenset({
    "response_matrix_derivation", "one_over_R2_derivation", "mouth_mechanism_derivation",
})


# ---------------------------------------------------------------------------
# Four DECIDED conditional ensemble coefficients, consuming m_gg and S_gg.
# ---------------------------------------------------------------------------

mgg = sp.Symbol("m_gg", nonnegative=True, real=True)
sgg = sp.Symbol("S_gg", positive=True, real=True)
phi, j = sp.symbols("phi j", real=True)
q, g = sp.symbols("q g", positive=True, real=True)
lam = sp.Symbol("lambda", real=True)
s1, s2, eps = sp.symbols("s_1 s_2 epsilon", real=True)

# Reshape the source V ensemble through the two-mouth response, without
# reconstructing stage031's response matrix or its 1/R^2 mechanism.
response_v = sp.Matrix([[sgg, eps * mgg], [eps * mgg, sgg]])
held_v = sp.Matrix([s1 * phi, s2 * phi])
reacting_sources = sp.simplify(response_v.inv() * held_v)
omega_v = -sp.Rational(1, 2) * (held_v.T * reacting_sources)[0]
derived_A_V = sp.factor(sp.expand(sp.series(omega_v, eps, 0, 2).removeO()).coeff(eps) / (s1 * s2))

# Same-field source/reservoir work is the source engine's ensemble route.
stored_M = mgg * (q + g) ** 2
derived_A_M = sp.factor(stored_M - 2 * g * mgg * (q + g))
stored_J = mgg * (j + g) ** 2
derived_A_J = sp.factor(stored_J - 2 * (j + g) * mgg * (j + g))
derived_A_MIXED = sp.factor(stored_M - 2 * (g + lam * q) * mgg * (q + g))

A_V = mgg * phi**2 / sgg**2
A_J = -mgg * (j + g) ** 2
A_M = mgg * (q**2 - g**2)
A_MIXED = mgg * ((1 - 2 * lam) * q**2 - 2 * lam * q * g - g**2)
COEFFICIENTS = OrderedDict(V=A_V, J=A_J, M=A_M, MIXED=A_MIXED)
MIXED_ROOT = sp.factor(sp.solve(sp.Eq(A_MIXED / mgg, 0), lam)[0])


def sign_label(value: sp.Expr) -> int:
    clean = sp.simplify(value)
    if clean == 0:
        return 0
    if clean.is_positive is True:
        return 1
    if clean.is_negative is True:
        return -1
    raise AuditFailure("SIGN_LABEL", compact(clean))


def strict_sign_context(zg_value: sp.Expr, mgg_value: sp.Expr, sgg_value: sp.Expr,
                        phi_value: sp.Expr, j_value: sp.Expr, g_value: sp.Expr) -> bool:
    return bool(
        zg_value != 0 and mgg_value > 0 and sgg_value > 0
        and phi_value != 0 and (j_value + g_value) != 0
    )


# ---------------------------------------------------------------------------
# Q-AMEND, Q-BC, Q-COMBINE, and Q-MAG source reshaping.
# ---------------------------------------------------------------------------

ZERO_LEDGER = (
    "bulk r_BH,r_B^2H^2,Hrho,Hdelta-rho,Hdt-theta,Hgrad-theta outside Omega_mouth",
    "dynamic J_m modulation and neighbor-source response",
    "r_B-u_L,r_B-divu,r_B-u_T,H-u_T,u_L-u_T and scalar-transverse mixing",
    "cross kinetic and one-time-derivative Berry terms",
    "u_L^2,h^2,higher gradients and cubic/quartic terms",
    "independent primitive B(divu)^2",
    "independent theta_B and phase-drain terms",
    "wall/collar/storage/dissipation terms",
    "dynamic drain and return-kernel responses",
    "direct drain sources and direct h contribution to e_c",
    "field-dependent geon derivatives",
    "viscosity/drag/E2/E3 terms",
    "E4/E5/E1 and mixture terms",
)


@dataclass(frozen=True)
class Amendment:
    source_row: str
    new_rows: tuple[str, ...]
    shold_scope: str
    zero_rows: tuple[str, ...]


REPLACE_AMENDMENT = Amendment(
    "core_holder_retypes_existing_h_source_BC", (), "r_B-1/2 only", ZERO_LEDGER,
)
ADD_AMENDMENT = Amendment(
    "existing_h_source_BC_unchanged", ("core_embedding_h_holding_row:R_w_odd",),
    "r_B-1/2 only", ZERO_LEDGER,
)


def amendment_sectors(rep: Amendment, add: Amendment) -> tuple[str, ...]:
    sectors: list[str] = []
    if rep.source_row != "core_holder_retypes_existing_h_source_BC" or rep.new_rows:
        sectors.append("replace_ledger")
    if add.source_row != "existing_h_source_BC_unchanged" or len(add.new_rows) != 1:
        sectors.append("add_ledger")
    if rep.shold_scope != "r_B-1/2 only" or add.shold_scope != "r_B-1/2 only":
        sectors.append("S_hold")
    if rep.zero_rows != ZERO_LEDGER or add.zero_rows != ZERO_LEDGER:
        sectors.append("zero_ledger")
    return tuple(sectors)


class BCClass(str, Enum):
    DIRICHLET_VALUE = "DIRICHLET_VALUE"
    FIXED_MONOPOLE = "FIXED_MONOPOLE"
    FIXED_SOURCE = "FIXED_SOURCE"
    MIXED = "MIXED"
    UNDETERMINED_ANALYTICALLY = "UNDETERMINED_ANALYTICALLY"


class Outcome(str, Enum):
    POSITIVE_R2 = "POSITIVE_R2"
    NEGATIVE_R2 = "NEGATIVE_R2"
    NULL = "NULL"
    POSITIVE_WRONG_RANGE = "POSITIVE_WRONG_RANGE"
    NEGATIVE_WRONG_RANGE = "NEGATIVE_WRONG_RANGE"
    NOT_INVARIANT = "outcome_not_invariant"


class Variant(str, Enum):
    REPLACE = "replace"
    ADD = "add"
    BOTH = "both"
    UNRESOLVED = "variant_unresolved"


class Magnitude(str, Enum):
    NO_FREE_FACTOR = "magnitude_no_free_factor"
    FREE_FACTOR = "magnitude_free_factor"


class Tier(str, Enum):
    GAPS_CLOSED = "tier_A_gaps_closed"
    CONDITIONAL = "tier_A_conditional"


@dataclass(frozen=True)
class BCEvidence:
    essential_value: bool
    fixed_conormal: bool
    mixed_relation: bool
    admissible_variation: bool
    signed_stationary: bool
    barrier_computed: bool
    missing_closure: tuple[str, ...]


@dataclass(frozen=True)
class BCStatus:
    kind: BCClass
    missing_closure: tuple[str, ...] = ()


def classify_bc(evidence: BCEvidence) -> BCStatus:
    if evidence.essential_value and not evidence.admissible_variation:
        return BCStatus(BCClass.DIRICHLET_VALUE)
    if evidence.fixed_conormal and evidence.admissible_variation:
        return BCStatus(BCClass.FIXED_MONOPOLE)
    if evidence.mixed_relation and evidence.admissible_variation:
        return BCStatus(BCClass.MIXED)
    if (not evidence.missing_closure and evidence.admissible_variation
            and not evidence.signed_stationary and not evidence.barrier_computed):
        return BCStatus(BCClass.FIXED_SOURCE)
    return BCStatus(BCClass.UNDETERMINED_ANALYTICALLY, evidence.missing_closure)


ACTUAL_EVIDENCE = BCEvidence(False, False, False, True, False, False,
                             ("missing parent-throat/boundary closure",))
FREE_CONTROL = BCEvidence(False, False, False, True, False, False, ())
VALUE_CONTROL = BCEvidence(True, False, False, False, True, True, ())
MONOPOLE_CONTROL = BCEvidence(False, True, False, True, False, True, ())
MIXED_CONTROL = BCEvidence(False, False, True, True, False, True, ())
ADMISSIBLE_CLASSES = (
    "DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED",
)

REPLACE_TOTALS = OrderedDict(
    DIRICHLET_VALUE=mgg * phi**2 / sgg**2,
    FIXED_MONOPOLE=mgg * q**2,
    FIXED_SOURCE=-mgg * j**2,
    MIXED=mgg * (1 - 2 * lam) * q**2,
)
ADD_TOTALS = OrderedDict(
    DIRICHLET_VALUE=A_V,
    FIXED_MONOPOLE=A_M,
    FIXED_SOURCE=A_J,
    MIXED=A_MIXED,
)
COMBINE_FACTS: OrderedDict[str, OrderedDict[str, Outcome]] = OrderedDict(
    DIRICHLET_VALUE=OrderedDict(replace=Outcome.POSITIVE_R2, add=Outcome.POSITIVE_R2),
    FIXED_MONOPOLE=OrderedDict(replace=Outcome.POSITIVE_R2, add=Outcome.NOT_INVARIANT),
    FIXED_SOURCE=OrderedDict(replace=Outcome.NEGATIVE_R2, add=Outcome.NEGATIVE_R2),
    MIXED=OrderedDict(replace=Outcome.NOT_INVARIANT, add=Outcome.NOT_INVARIANT),
)
VARIANT_REALIZATION = Variant.UNRESOLVED
MAGNITUDE_FACT = Magnitude.FREE_FACTOR
TIER_FACT = Tier.CONDITIONAL
HOOKS = OrderedDict(
    density="NO(no_local_prediction)",
    radial_monopole="UNDETERMINED(core source/conormal not fixed)",
    modulus="NO(continuous core modulus)",
    close_range="UNDETERMINED(out_of_scope_R_comparable_to_r_e)",
)


def classify_samples(values: Iterable[sp.Expr]) -> str:
    signs = tuple(sign_label(sp.sympify(value)) for value in values)
    return "CONSTANT_OUTCOME" if len(set(signs)) == 1 else "outcome_not_invariant"


# ---------------------------------------------------------------------------
# Typed §4 facts, fresh production tuple, and first-match adjudication.
# ---------------------------------------------------------------------------

BC_DOMAIN = (
    "DIRICHLET_VALUE", "FIXED_MONOPOLE", "FIXED_SOURCE", "MIXED",
    "UNDETERMINED_ANALYTICALLY",
)
OUTCOME_DOMAIN = (
    "POSITIVE_R2", "NEGATIVE_R2", "NULL", "POSITIVE_WRONG_RANGE",
    "NEGATIVE_WRONG_RANGE", "outcome_not_invariant",
)
VARIANT_DOMAIN = ("replace", "add", "both", "variant_unresolved")
MAG_DOMAIN = ("magnitude_no_free_factor", "magnitude_free_factor")
TIER_DOMAIN = ("tier_A_gaps_closed", "tier_A_conditional")


@dataclass(frozen=True)
class LandingFacts:
    qbc_status: BCStatus
    class_variant_outcomes: Mapping[str, Mapping[str, Outcome]]
    variant: Variant
    magnitude: Magnitude
    tier: Tier
    inconsistency_sectors: tuple[str, ...]


@dataclass(frozen=True)
class ProductionTuple:
    qbc: str
    replace_outcome: str
    add_outcome: str
    variant: str
    magnitude: str
    tier: str
    internal: bool
    all_classes_agree: bool
    mixed_range_invariant: bool
    overall_outcome: str


def typed_landing_facts(facts: LandingFacts) -> bool:
    if type(facts.qbc_status) is not BCStatus or type(facts.qbc_status.kind) is not BCClass:
        return False
    if tuple(facts.class_variant_outcomes) != ADMISSIBLE_CLASSES:
        return False
    for cls in ADMISSIBLE_CLASSES:
        row = facts.class_variant_outcomes[cls]
        if tuple(row) != ("replace", "add"):
            return False
        if any(type(row[key]) is not Outcome for key in ("replace", "add")):
            return False
    return (
        type(facts.variant) is Variant
        and type(facts.magnitude) is Magnitude
        and type(facts.tier) is Tier
        and all(type(x) is str for x in facts.inconsistency_sectors)
    )


def neutral_tokens(facts: LandingFacts) -> tuple[Any, ...]:
    return (
        facts.qbc_status,
        *(facts.class_variant_outcomes[cls][variant]
          for cls in ADMISSIBLE_CLASSES for variant in ("replace", "add")),
        facts.variant, facts.magnitude, facts.tier, *facts.inconsistency_sectors,
    )


def production_tuple(facts: LandingFacts) -> ProductionTuple:
    if not typed_landing_facts(facts):
        raise TypeError("landing accepts only typed neutral upstream facts")
    all_outcomes = tuple(
        facts.class_variant_outcomes[cls][variant]
        for cls in ADMISSIBLE_CLASSES for variant in ("replace", "add")
    )
    overall = all_outcomes[0] if len(set(all_outcomes)) == 1 else Outcome.NOT_INVARIANT

    def summarize(variant: str) -> Outcome:
        values = tuple(facts.class_variant_outcomes[cls][variant] for cls in ADMISSIBLE_CLASSES)
        return values[0] if len(set(values)) == 1 else overall

    realized = ((facts.variant.value,) if facts.variant in {Variant.REPLACE, Variant.ADD}
                else ("replace", "add"))
    realized_outcomes = tuple(
        facts.class_variant_outcomes[cls][variant]
        for cls in ADMISSIBLE_CLASSES for variant in realized
    )
    mixed = tuple(facts.class_variant_outcomes["MIXED"][variant] for variant in realized)
    return ProductionTuple(
        qbc=facts.qbc_status.kind.value,
        replace_outcome=summarize("replace").value,
        add_outcome=summarize("add").value,
        variant=facts.variant.value,
        magnitude=facts.magnitude.value,
        tier=facts.tier.value,
        internal=bool(facts.inconsistency_sectors),
        all_classes_agree=len(set(realized_outcomes)) == 1,
        mixed_range_invariant=all(value is not Outcome.NOT_INVARIANT for value in mixed),
        overall_outcome=overall.value,
    )


def adjudicate(case: ProductionTuple | Mapping[str, Any], *, illegal_precedence: bool = False) -> str:
    if isinstance(case, ProductionTuple):
        row = case.__dict__
    else:
        row = case
    qbc = row["qbc"]
    ro = row["replace_outcome"]
    ao = row["add_outcome"]
    variant = row["variant"]
    magnitude = row["magnitude"]
    tier = row["tier"]
    internal = row["internal"]
    classes_agree = row["all_classes_agree"]
    mixed_invariant = row["mixed_range_invariant"]

    if internal:
        return "NO_GO(sector)"
    if variant == "replace":
        comparison, variant_resolved = ro, True
    elif variant == "add":
        comparison, variant_resolved = ao, True
    elif ro == ao:
        comparison, variant_resolved = ro, True
    else:
        comparison, variant_resolved = None, False
    if qbc == "UNDETERMINED_ANALYTICALLY":
        class_resolved = classes_agree
    elif qbc == "MIXED":
        class_resolved = tier == "tier_A_gaps_closed" and mixed_invariant
    else:
        class_resolved = tier == "tier_A_gaps_closed"
    magnitude_determined = comparison is not None and comparison != "outcome_not_invariant"
    unconditional = class_resolved and variant_resolved and magnitude_determined

    if illegal_precedence and qbc == "UNDETERMINED_ANALYTICALLY":
        return "R1_REQUIRED(bc_selection)"
    if unconditional:
        if comparison == "POSITIVE_R2":
            return "SIGN_EARNED" if magnitude == "magnitude_no_free_factor" else "R1_REQUIRED(magnitude)"
        if comparison == "NEGATIVE_R2":
            return "MECHANISM_FALSIFIED(wrong_sign)"
        if comparison in {"POSITIVE_WRONG_RANGE", "NEGATIVE_WRONG_RANGE"}:
            return "MECHANISM_FALSIFIED(wrong_range)"
        if comparison == "NULL":
            return "R1_REQUIRED(subleading)"
    if qbc == "UNDETERMINED_ANALYTICALLY" and not classes_agree:
        return "R1_REQUIRED(bc_selection)"
    if qbc == "MIXED" and not mixed_invariant:
        return "R1_REQUIRED(mixed_bc_parameters)"
    if variant in {"both", "variant_unresolved"} and ro != ao:
        return "R1_REQUIRED(variant_selection)"
    if tier == "tier_A_conditional" or comparison == "outcome_not_invariant":
        return "R1_REQUIRED(sign_and_magnitude)"
    return "R1_REQUIRED(unclassified)"


def declarative_oracle(row: Mapping[str, Any]) -> str:
    if row["internal"]:
        return "NO_GO(sector)"
    ro, ao, variant = row["replace_outcome"], row["add_outcome"], row["variant"]
    common = ro if variant == "replace" else ao if variant == "add" else ro if ro == ao else None
    variant_ok = variant in {"replace", "add"} or ro == ao
    qbc, tier = row["qbc"], row["tier"]
    class_ok = (
        row["all_classes_agree"] if qbc == "UNDETERMINED_ANALYTICALLY"
        else tier == "tier_A_gaps_closed" and row["mixed_range_invariant"] if qbc == "MIXED"
        else tier == "tier_A_gaps_closed"
    )
    if class_ok and variant_ok and common != "outcome_not_invariant":
        return {
            "POSITIVE_R2": ("SIGN_EARNED" if row["magnitude"] == "magnitude_no_free_factor"
                            else "R1_REQUIRED(magnitude)"),
            "NEGATIVE_R2": "MECHANISM_FALSIFIED(wrong_sign)",
            "POSITIVE_WRONG_RANGE": "MECHANISM_FALSIFIED(wrong_range)",
            "NEGATIVE_WRONG_RANGE": "MECHANISM_FALSIFIED(wrong_range)",
            "NULL": "R1_REQUIRED(subleading)",
        }[common]
    predicates = (
        (qbc == "UNDETERMINED_ANALYTICALLY" and not row["all_classes_agree"],
         "R1_REQUIRED(bc_selection)"),
        (qbc == "MIXED" and not row["mixed_range_invariant"],
         "R1_REQUIRED(mixed_bc_parameters)"),
        (variant in {"both", "variant_unresolved"} and ro != ao,
         "R1_REQUIRED(variant_selection)"),
        (tier == "tier_A_conditional" or common == "outcome_not_invariant",
         "R1_REQUIRED(sign_and_magnitude)"),
    )
    return next((landing for predicate, landing in predicates if predicate), "R1_REQUIRED(unclassified)")


@dataclass(frozen=True)
class TruthTable:
    exact: bool
    total: int
    counts: Mapping[str, int]
    digest: str
    rows: tuple[str, ...]


def exhaustive_truth_table() -> TruthTable:
    keys = (
        "qbc", "replace_outcome", "add_outcome", "variant", "magnitude",
        "tier", "internal", "all_classes_agree", "mixed_range_invariant",
    )
    rows: list[str] = []
    counts: Counter[str] = Counter()
    exact = True
    for values in itertools.product(
        BC_DOMAIN, OUTCOME_DOMAIN, OUTCOME_DOMAIN, VARIANT_DOMAIN, MAG_DOMAIN,
        TIER_DOMAIN, (False, True), (False, True), (False, True),
    ):
        case = dict(zip(keys, values))
        got = adjudicate(case)
        want = declarative_oracle(case)
        exact &= got == want and "unclassified" not in got
        counts[got] += 1
        rows.append("|".join(map(str, values)) + "|" + got)
    serialized = "\n".join(rows)
    return TruthTable(exact, len(rows), dict(sorted(counts.items())),
                      hashlib.sha256(serialized.encode("utf-8")).hexdigest(), tuple(rows))


PRODUCTION_FACTS = LandingFacts(
    qbc_status=classify_bc(ACTUAL_EVIDENCE),
    class_variant_outcomes=COMBINE_FACTS,
    variant=VARIANT_REALIZATION,
    magnitude=MAGNITUDE_FACT,
    tier=TIER_FACT,
    inconsistency_sectors=amendment_sectors(REPLACE_AMENDMENT, ADD_AMENDMENT),
)
PRODUCTION_TUPLE = production_tuple(PRODUCTION_FACTS)
PRODUCTION_LANDING = adjudicate(PRODUCTION_TUPLE)
TRUTH = exhaustive_truth_table()


def landing_ownership_guard(facts: LandingFacts, emitted_landing: str,
                            carried_tokens: tuple[Any, ...] | None = None) -> bool:
    if not typed_landing_facts(facts):
        return False
    expected = neutral_tokens(facts)
    carried = expected if carried_tokens is None else carried_tokens
    if len(carried) != len(expected) or any(got is not want for got, want in zip(carried, expected)):
        return False
    if any(isinstance(token, str) and ("R1_REQUIRED" in token or "SIGN_EARNED" in token)
           for token in carried):
        return False
    fresh_tuple = production_tuple(facts)
    fresh_landing = adjudicate(fresh_tuple)
    return emitted_landing == fresh_landing


# ---------------------------------------------------------------------------
# Canonical source-to-stage manifest: all 36 source teeth plus App/Q extras.
# ---------------------------------------------------------------------------

SOURCE_TOOTH_IDS = (
    "FIELD_IDENTITY_UNITS", "FIELD_PARENT_MAP", "FIELD_LIVE_QH",
    "ACTION_TRANSCRIPTION", "AMEND_REPLACE", "AMEND_ADD", "SHOLD_SCOPE",
    "ZERO_LEDGER", "MATRIX_DETERMINANT", "BC_ACTUAL_GAP", "BC_FREE_CONTROL",
    "BC_VALUE_CONTROL", "BC_MONOPOLE_CONTROL", "BC_MIXED_CONTROL",
    "FORCE_V_FUNCTIONAL", "FORCE_M_HWORK", "FORCE_J_FUNCTIONAL",
    "FORCE_MIXED_FUNCTIONAL", "MIXED_FULL_RANGE", "FALLOFF", "UNITS_RESTORED",
    "COMBINE_REPLACE", "COMBINE_ADD", "NO_DOUBLE_COUNT", "RANGE_SIGN_FLIP",
    "RANGE_TOUCH_ZERO", "RANGE_SUBDOMINANT", "MAG_FREE_FACTOR", "DENSITY_HOOK",
    "MONOPOLE_HOOK", "MODULUS_HOOK", "VERDICT_TOTALITY", "VERDICT_PRECEDENCE",
    "LANDING_OWNERSHIP", "TARGET_BLINDNESS", "DUAL_ENGINE_TERMS",
)
EXTRA_SOURCE_CLAIMS = (
    "EXTRA_VARIANT_UNRESOLVED", "EXTRA_QMAG_CLOSE_RANGE",
    "EXTRA_SGG_SELF_RESPONSE", "EXTRA_G0_RESPONSE_WITNESS",
    "EXTRA_RESPONSE_MATRIX_ESCAPE_FACTORS", "EXTRA_EXTERIOR_CORE_GAP",
    "EXTRA_SECTION4_DIGEST_COUNTS", "EXTRA_APP_E_ACCEPTANCE",
)


def manifest_entry(identifier: str, disposition: str, owner: str) -> tuple[str, str, str]:
    return identifier, disposition, owner


SOURCE_MANIFEST = (
    manifest_entry("FIELD_IDENTITY_UNITS", "SCOPED_OUT", "STAGE031_MECHANISM"),
    manifest_entry("FIELD_PARENT_MAP", "SCOPED_OUT", "STAGE031_MECHANISM"),
    manifest_entry("FIELD_LIVE_QH", "SCOPED_OUT", "STAGE031_MECHANISM"),
    manifest_entry("ACTION_TRANSCRIPTION", "SCOPED_OUT", "STAGE030_031_CITED"),
    manifest_entry("AMEND_REPLACE", "PRESERVED", "STAGE032_Q_AMEND"),
    manifest_entry("AMEND_ADD", "PRESERVED", "STAGE032_Q_AMEND"),
    manifest_entry("SHOLD_SCOPE", "PRESERVED", "STAGE032_Q_AMEND"),
    manifest_entry("ZERO_LEDGER", "PRESERVED", "STAGE032_Q_AMEND"),
    manifest_entry("MATRIX_DETERMINANT", "SCOPED_OUT", "STAGE031_RESPONSE"),
    manifest_entry("BC_ACTUAL_GAP", "REPLACED_BY_STRONGER", "STAGE032_TYPED_Q_BC"),
    manifest_entry("BC_FREE_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"),
    manifest_entry("BC_VALUE_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"),
    manifest_entry("BC_MONOPOLE_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"),
    manifest_entry("BC_MIXED_CONTROL", "PRESERVED", "STAGE032_TYPED_Q_BC"),
    manifest_entry("FORCE_V_FUNCTIONAL", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"),
    manifest_entry("FORCE_M_HWORK", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"),
    manifest_entry("FORCE_J_FUNCTIONAL", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"),
    manifest_entry("FORCE_MIXED_FUNCTIONAL", "REPLACED_BY_STRONGER", "STAGE032_FORMULA_SIGN"),
    manifest_entry("MIXED_FULL_RANGE", "REPLACED_BY_STRONGER", "STAGE032_THREE_REGIMES"),
    manifest_entry("FALLOFF", "SCOPED_OUT", "STAGE031_NEUTRAL_SHELL"),
    manifest_entry("UNITS_RESTORED", "REPLACED_BY_STRONGER", "STAGE032_UNITS_FIREWALL"),
    manifest_entry("COMBINE_REPLACE", "PRESERVED", "STAGE032_Q_COMBINE"),
    manifest_entry("COMBINE_ADD", "PRESERVED", "STAGE032_Q_COMBINE"),
    manifest_entry("NO_DOUBLE_COUNT", "PRESERVED", "STAGE032_Q_COMBINE"),
    manifest_entry("RANGE_SIGN_FLIP", "PRESERVED", "STAGE032_RANGE"),
    manifest_entry("RANGE_TOUCH_ZERO", "PRESERVED", "STAGE032_RANGE"),
    manifest_entry("RANGE_SUBDOMINANT", "PRESERVED", "STAGE032_RANGE"),
    manifest_entry("MAG_FREE_FACTOR", "PRESERVED", "STAGE032_Q_MAG"),
    manifest_entry("DENSITY_HOOK", "PRESERVED", "STAGE032_Q_MAG"),
    manifest_entry("MONOPOLE_HOOK", "PRESERVED", "STAGE032_Q_MAG"),
    manifest_entry("MODULUS_HOOK", "PRESERVED", "STAGE032_Q_MAG"),
    manifest_entry("VERDICT_TOTALITY", "REPLACED_BY_STRONGER", "STAGE032_GRID_COUNTS_DIGEST"),
    manifest_entry("VERDICT_PRECEDENCE", "PRESERVED", "STAGE032_LADDER"),
    manifest_entry("LANDING_OWNERSHIP", "REPLACED_BY_STRONGER", "STAGE032_UPSTREAM_RELANDING"),
    manifest_entry("TARGET_BLINDNESS", "PRESERVED", "STAGE032_NEUTRAL_UPSTREAM"),
    manifest_entry("DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER", "STAGE032_INDEPENDENT_ENGINES"),
    manifest_entry("EXTRA_VARIANT_UNRESOLVED", "PRESERVED", "STAGE032_Q_COMBINE"),
    manifest_entry("EXTRA_QMAG_CLOSE_RANGE", "PRESERVED", "STAGE032_Q_MAG"),
    manifest_entry("EXTRA_SGG_SELF_RESPONSE", "SCOPED_OUT", "STAGE031_CONSUMED"),
    manifest_entry("EXTRA_G0_RESPONSE_WITNESS", "SCOPED_OUT", "STAGE031_RESPONSE"),
    manifest_entry("EXTRA_RESPONSE_MATRIX_ESCAPE_FACTORS", "SCOPED_OUT", "STAGE031_RESPONSE"),
    manifest_entry("EXTRA_EXTERIOR_CORE_GAP", "SCOPED_OUT", "STAGE031_CONSUMED_NAMED_FACT"),
    manifest_entry("EXTRA_SECTION4_DIGEST_COUNTS", "REPLACED_BY_STRONGER", "STAGE032_COMMITTED_LITERAL"),
    manifest_entry("EXTRA_APP_E_ACCEPTANCE", "REPLACED_BY_STRONGER", "STAGE032_STANDALONE_ABLATIONS"),
)


def canonical_manifest_text(manifest: Iterable[tuple[str, str, str]]) -> str:
    return "\n".join("|".join(row) for row in sorted(manifest))


def manifest_sha256(manifest: Iterable[tuple[str, str, str]]) -> str:
    return hashlib.sha256(canonical_manifest_text(manifest).encode("utf-8")).hexdigest()


def run_consumption_and_scope() -> None:
    heading("Consumed upstream bundle and ownership firewall")
    consumed = primitive(
        "PASS_STAGE031_CONSUMED_INPUTS", CONSUMED_STAGE031,
        tuple(x for x in CONSUMED_STAGE031 if x != "S_gg"),
    )
    expect_bool(
        "PASS_STAGE031_CONSUMED_INPUTS",
        consumed == CONSUMED_STAGE031 and {
            "m", "m_gg=B_eff*z_g^2/D", "det(m)=z_g^2/D", "S_gg",
            "A=m_gg*C", "NONZERO_HA_REQUIRES_CORE_HOLDER", "s1*s2/(4*pi*R^2)",
        } == set(consumed),
        consumed,
    )
    d_star = primitive("PASS_STAGE030_TRANSITIVE_D_WITNESS", sp.Rational(7, 4), sp.Rational(5, 4))
    expect_zero("PASS_STAGE030_TRANSITIVE_D_WITNESS", d_star - sp.Rational(7, 4))
    owned = primitive("PASS_SCOPE_FIREWALL", OWNED_STAGE032,
                      OWNED_STAGE032 | {"response_matrix_derivation"})
    expect_bool("PASS_SCOPE_FIREWALL", not (set(owned) & set(FORBIDDEN_REDERIVATIONS)), owned)
    core_fact = primitive("PASS_NONZERO_HA_REQUIRES_CORE_HOLDER",
                          "NONZERO_HA_REQUIRES_CORE_HOLDER", "")
    expect_bool("PASS_NONZERO_HA_REQUIRES_CORE_HOLDER",
                core_fact == "NONZERO_HA_REQUIRES_CORE_HOLDER")
    print("      stage031 mechanism cited; stage030 D*=7/4 cited transitively; no mechanism re-derived")


def run_ensembles() -> None:
    heading("Four DECIDED conditional ensembles")
    av_test = primitive("PASS_A_V_FORMULA", derived_A_V, derived_A_V + mgg * phi**2)
    expect_zero("PASS_A_V_FORMULA", av_test - A_V, {"A_V": av_test})
    aj_test = primitive("PASS_A_J_FORMULA", derived_A_J, -derived_A_J)
    expect_zero("PASS_A_J_FORMULA", aj_test - A_J, {"A_J": aj_test})
    am_test = primitive("PASS_A_M_FORMULA", derived_A_M, mgg * q**2)
    expect_zero("PASS_A_M_FORMULA", am_test - A_M, {"A_M": am_test})
    mixed_test = primitive("PASS_A_MIXED_FORMULA", derived_A_MIXED,
                           derived_A_MIXED + 2 * lam * q * g * mgg)
    expect_zero("PASS_A_MIXED_FORMULA", mixed_test - A_MIXED, {"A_MIXED": mixed_test})

    mnn = sp.Symbol("m_nonnegative", nonnegative=True, real=True)
    s_pos = sp.Symbol("S_positive", positive=True, real=True)
    p_real, j_real, g_real = sp.symbols("phi_real j_real g_real", real=True)
    weak_v = sp.ask(sp.Q.nonnegative(mnn * p_real**2 / s_pos**2)) is True
    weak_j = sp.ask(sp.Q.nonpositive(-mnn * (j_real + g_real) ** 2)) is True
    weak_ok = primitive("PASS_WEAK_SIGNS_GENERAL", weak_v and weak_j, False)
    expect_bool("PASS_WEAK_SIGNS_GENERAL", weak_ok,
                {"A_V>=0": weak_v, "A_J<=0": weak_j, "strict": False})

    strict_phi = primitive("PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS", sp.Integer(3), sp.Integer(0))
    strict_subs = {mgg: sp.Rational(8, 7), sgg: sp.Integer(2), phi: strict_phi,
                   j: sp.Integer(2), g: sp.Integer(1)}
    strict_active = strict_sign_context(1, strict_subs[mgg], strict_subs[sgg],
                                        strict_subs[phi], strict_subs[j], strict_subs[g])
    strict_v = sp.simplify(A_V.subs(strict_subs))
    strict_j = sp.simplify(A_J.subs(strict_subs))
    expect_bool("PASS_STRICT_SIGNS_NONDEGENERATE_WITNESS",
                strict_active and strict_v > 0 and strict_j < 0,
                {"m_gg": strict_subs[mgg], "S_gg": 2, "A_V": strict_v, "A_J": strict_j})

    b_eff, delta, zg = sp.symbols("B_eff D z_g", positive=True, real=True)
    degenerate_zg = primitive("PASS_DEGENERATE_Z_G_ZERO", sp.Integer(0), sp.Integer(1))
    degenerate_mgg = b_eff * degenerate_zg**2 / delta
    degenerate = tuple(sp.simplify(expr.subs(mgg, degenerate_mgg)) for expr in COEFFICIENTS.values())
    degenerate_strict = strict_sign_context(degenerate_zg, degenerate_mgg, 1, 1, 1, 1)
    expect_bool("PASS_DEGENERATE_Z_G_ZERO",
                all(value == 0 for value in degenerate) and not degenerate_strict,
                {"z_g": degenerate_zg, "A_X": degenerate, "strict_asserts_active": degenerate_strict})

    negative_q = primitive("PASS_A_M_INDEFINITE", sp.Integer(1), sp.Integer(3))
    m_pos = sp.simplify(A_M.subs({mgg: 1, q: 3, g: 1}))
    m_neg = sp.simplify(A_M.subs({mgg: 1, q: negative_q, g: 2}))
    expect_bool("PASS_A_M_INDEFINITE", m_pos > 0 and m_neg < 0,
                {"positive witness": m_pos, "negative witness": m_neg})

    endpoint_one = primitive("PASS_MIXED_ENDPOINTS", sp.Integer(1), sp.Integer(0))
    endpoints = (sp.factor(A_MIXED.subs(lam, 0)), sp.factor(A_MIXED.subs(lam, endpoint_one)))
    expect_bool("PASS_MIXED_ENDPOINTS",
                sp.simplify(endpoints[0] - A_M) == 0
                and sp.simplify(endpoints[1] + mgg * (q + g) ** 2) == 0,
                endpoints)

    root_test = primitive("PASS_MIXED_INTERIOR_ZERO", MIXED_ROOT, MIXED_ROOT + sp.Rational(1, 2))
    expect_bool("PASS_MIXED_INTERIOR_ZERO",
                sp.simplify(root_test - (q - g) / (2 * q)) == 0
                and sp.simplify(A_MIXED.subs(lam, root_test)) == 0,
                root_test)

    d = sp.Symbol("d", positive=True, real=True)
    g0 = sp.Symbol("g0", positive=True, real=True)
    m0 = sp.Symbol("m0", positive=True, real=True)
    q_gt = g0 + d
    root_gt = sp.simplify((q_gt - g0) / (2 * q_gt))
    positive_lambda = primitive("PASS_MIXED_REGIME_Q_GT_G", root_gt / 2, (1 + root_gt) / 2)
    gt_values = tuple(sp.simplify(A_MIXED.subs({mgg: m0, q: q_gt, g: g0, lam: x}))
                      for x in (positive_lambda, root_gt, (1 + root_gt) / 2))
    expect_bool("PASS_MIXED_REGIME_Q_GT_G",
                root_gt.is_positive is True and sp.factor(1 - root_gt).is_positive is True
                and sign_label(gt_values[0]) == 1 and sign_label(gt_values[1]) == 0
                and sign_label(gt_values[2]) == -1,
                {"lambda*": root_gt, "signs": tuple(map(sign_label, gt_values))})

    eq_zero_lambda = primitive("PASS_MIXED_REGIME_Q_EQ_G", sp.Integer(0), sp.Rational(1, 2))
    eq_zero = sp.simplify(A_MIXED.subs({mgg: m0, q: g0, g: g0, lam: eq_zero_lambda}))
    eq_interior = sp.simplify(A_MIXED.subs({mgg: m0, q: g0, g: g0, lam: sp.Rational(1, 2)}))
    expect_bool("PASS_MIXED_REGIME_Q_EQ_G", eq_zero == 0 and eq_interior.is_negative is True,
                {"A(0)": eq_zero, "A(1/2)": eq_interior})

    q_lt = g0
    g_lt = g0 + d
    root_lt_canonical = sp.simplify((q_lt - g_lt) / (2 * q_lt))
    root_lt = primitive("PASS_MIXED_REGIME_Q_LT_G", root_lt_canonical, -root_lt_canonical)
    lt0 = sp.factor(A_MIXED.subs({mgg: m0, q: q_lt, g: g_lt, lam: 0}))
    lt1 = sp.factor(A_MIXED.subs({mgg: m0, q: q_lt, g: g_lt, lam: 1}))
    lt_slope = sp.factor(sp.diff(A_MIXED, lam).subs({mgg: m0, q: q_lt, g: g_lt}))
    expect_bool("PASS_MIXED_REGIME_Q_LT_G",
                root_lt.is_negative is True and lt0.is_negative is True and lt1.is_negative is True
                and lt_slope.is_negative is True,
                {"lambda*": root_lt, "endpoints": (lt0, lt1)})

    full_domain_claim = primitive("PASS_MIXED_FULL_DOMAIN_SPANS",
                                  "FULL_PARAMETER_AND_MAGNITUDE_DOMAIN",
                                  "EVERY_FIXED_Q_OVER_G")
    full_samples = (
        A_MIXED.subs({mgg: 1, q: 3, g: 1, lam: 0}),
        A_MIXED.subs({mgg: 1, q: 1, g: 1, lam: 0}),
        A_MIXED.subs({mgg: 1, q: 1, g: 2, lam: 0}),
    )
    expect_bool("PASS_MIXED_FULL_DOMAIN_SPANS",
                full_domain_claim == "FULL_PARAMETER_AND_MAGNITUDE_DOMAIN"
                and tuple(map(sign_label, full_samples)) == (1, 0, -1),
                {"meaning": full_domain_claim, "cross-regime signs": tuple(map(sign_label, full_samples))})

    typing = primitive("PASS_DECIDED_CONDITIONAL_TYPING",
                       ("coefficients=DECIDED_given_class", "bc_selection=R1",
                        "magnitude=R1", "mixed_parameters=R1", "variant=R1"),
                       ("coefficients=DECIDED_given_class", "bc_selection=DECIDED",
                        "magnitude=R1", "mixed_parameters=R1", "variant=R1"))
    expect_bool("PASS_DECIDED_CONDITIONAL_TYPING",
                typing[0] == "coefficients=DECIDED_given_class"
                and typing[1:] == ("bc_selection=R1", "magnitude=R1",
                                   "mixed_parameters=R1", "variant=R1"), typing)
    print("      weak signs hold generally; strict V/J signs are witness-only; z_g=0 activates no strict assert")


def run_source_blocks() -> None:
    heading("Q-AMEND, typed Q-BC, Q-COMBINE, and Q-MAG")
    rep = primitive("PASS_AMEND_REPLACE", REPLACE_AMENDMENT,
                    replace(REPLACE_AMENDMENT, new_rows=("illegal",)))
    expect_bool("PASS_AMEND_REPLACE", "replace_ledger" not in amendment_sectors(rep, ADD_AMENDMENT), rep)
    add = primitive("PASS_AMEND_ADD", ADD_AMENDMENT,
                    replace(ADD_AMENDMENT, new_rows=ADD_AMENDMENT.new_rows + ("second",)))
    expect_bool("PASS_AMEND_ADD", "add_ledger" not in amendment_sectors(REPLACE_AMENDMENT, add), add)
    zeros = primitive("PASS_ZERO_LEDGER_13", ZERO_LEDGER, ZERO_LEDGER[:-1])
    expect_bool("PASS_ZERO_LEDGER_13", len(zeros) == 13 and zeros == ZERO_LEDGER, len(zeros))
    shold = primitive("PASS_S_HOLD_SCOPE", "r_B-1/2 only", "r_B-1/2 and h")
    expect_bool("PASS_S_HOLD_SCOPE", shold == "r_B-1/2 only", shold)
    sectors = primitive("PASS_INTERNAL_INCONSISTENCY_NONE",
                        amendment_sectors(REPLACE_AMENDMENT, ADD_AMENDMENT), ("add_ledger",))
    expect_bool("PASS_INTERNAL_INCONSISTENCY_NONE", sectors == (), {"internal_inconsistency": "none" if not sectors else sectors})

    actual = primitive("PASS_BC_ACTUAL_CLASSIFIER", ACTUAL_EVIDENCE,
                       replace(ACTUAL_EVIDENCE, missing_closure=()))
    expect_bool("PASS_BC_ACTUAL_CLASSIFIER",
                classify_bc(actual).kind is BCClass.UNDETERMINED_ANALYTICALLY,
                classify_bc(actual))
    free = primitive("PASS_BC_FREE_CONTROL", FREE_CONTROL,
                     replace(FREE_CONTROL, essential_value=True, admissible_variation=False))
    expect_bool("PASS_BC_FREE_CONTROL", classify_bc(free).kind is BCClass.FIXED_SOURCE, classify_bc(free))
    value = primitive("PASS_BC_VALUE_CONTROL", VALUE_CONTROL,
                      replace(VALUE_CONTROL, essential_value=False, admissible_variation=True,
                              missing_closure=("lost_value",)))
    expect_bool("PASS_BC_VALUE_CONTROL", classify_bc(value).kind is BCClass.DIRICHLET_VALUE, classify_bc(value))
    mono = primitive("PASS_BC_MONOPOLE_CONTROL", MONOPOLE_CONTROL,
                     replace(MONOPOLE_CONTROL, fixed_conormal=False, missing_closure=("lost_flux",)))
    expect_bool("PASS_BC_MONOPOLE_CONTROL", classify_bc(mono).kind is BCClass.FIXED_MONOPOLE, classify_bc(mono))
    mixed = primitive("PASS_BC_MIXED_CONTROL", MIXED_CONTROL,
                      replace(MIXED_CONTROL, mixed_relation=False, missing_closure=("lost_mix",)))
    expect_bool("PASS_BC_MIXED_CONTROL", classify_bc(mixed).kind is BCClass.MIXED, classify_bc(mixed))
    classes = primitive("PASS_ADMISSIBLE_CLASSES", ADMISSIBLE_CLASSES, ADMISSIBLE_CLASSES[:-1])
    expect_bool("PASS_ADMISSIBLE_CLASSES", classes == ADMISSIBLE_CLASSES and len(classes) == 4, classes)

    replace_test = OrderedDict(REPLACE_TOTALS)
    if ACTIVE_MUTATION == "PASS_REPLACE_TOTALS":
        replace_test["FIXED_MONOPOLE"] += mgg * g**2
    replace_expected = (A_V, mgg * q**2, -mgg * j**2, mgg * (1 - 2 * lam) * q**2)
    expect_bool("PASS_REPLACE_TOTALS",
                all(sp.simplify(got - want) == 0 for got, want in zip(replace_test.values(), replace_expected)),
                replace_test)
    add_test = OrderedDict(ADD_TOTALS)
    if ACTIVE_MUTATION == "PASS_ADD_TOTALS":
        add_test["FIXED_MONOPOLE"] = mgg * q**2
    expect_bool("PASS_ADD_TOTALS",
                all(sp.simplify(add_test[name] - COEFFICIENTS[key]) == 0
                    for name, key in zip(ADMISSIBLE_CLASSES, ("V", "M", "J", "MIXED"))), add_test)
    variant = primitive("PASS_VARIANT_UNRESOLVED", VARIANT_REALIZATION, Variant.REPLACE)
    expect_bool("PASS_VARIANT_UNRESOLVED", variant is Variant.UNRESOLVED, variant)
    rebuilt_m = sp.factor(stored_M - 2 * g * mgg * (q + g))
    rebuilt_m = primitive("PASS_NO_DOUBLE_COUNT", rebuilt_m, rebuilt_m + mgg * g**2)
    expect_zero("PASS_NO_DOUBLE_COUNT", rebuilt_m - A_M, rebuilt_m)
    outcome_facts = primitive("PASS_OUTCOME_NOT_INVARIANT", COMBINE_FACTS,
                              OrderedDict((cls, OrderedDict(replace=Outcome.POSITIVE_R2,
                                                           add=Outcome.POSITIVE_R2))
                                          for cls in ADMISSIBLE_CLASSES))
    flattened = tuple(outcome_facts[cls][variant] for cls in ADMISSIBLE_CLASSES
                      for variant in ("replace", "add"))
    overall = flattened[0] if len(set(flattened)) == 1 else Outcome.NOT_INVARIANT
    expect_bool("PASS_OUTCOME_NOT_INVARIANT", overall is Outcome.NOT_INVARIANT, overall)

    magnitude_factors = primitive("PASS_MAGNITUDE_FREE_FACTOR", {"c_a", "c_xi"}, {"c_a"})
    expect_bool("PASS_MAGNITUDE_FREE_FACTOR",
                magnitude_factors == {"c_a", "c_xi"} and MAGNITUDE_FACT is Magnitude.FREE_FACTOR,
                magnitude_factors)
    density = primitive("PASS_QMAG_DENSITY_HOOK", HOOKS["density"], "YES(local_prediction)")
    expect_bool("PASS_QMAG_DENSITY_HOOK", density == "NO(no_local_prediction)", density)
    monopole = primitive("PASS_QMAG_MONOPOLE_HOOK", HOOKS["radial_monopole"], "DECIDED(nonzero)")
    expect_bool("PASS_QMAG_MONOPOLE_HOOK", monopole.startswith("UNDETERMINED"), monopole)
    modulus = primitive("PASS_QMAG_MODULUS_HOOK", HOOKS["modulus"], "YES(universal)")
    expect_bool("PASS_QMAG_MODULUS_HOOK", modulus.startswith("NO("), modulus)
    close_range = primitive("PASS_QMAG_CLOSE_RANGE_HOOK", HOOKS["close_range"], "DECIDED")
    expect_bool("PASS_QMAG_CLOSE_RANGE_HOOK", close_range.startswith("UNDETERMINED"), close_range)
    coblocker = primitive("PASS_MAGNITUDE_COBLOCKER", MAGNITUDE_COBLOCKER,
                          "R1_REQUIRED(bc_selection)")
    expect_bool("PASS_MAGNITUDE_COBLOCKER", coblocker == "R1_REQUIRED(magnitude)", coblocker)


def run_ranges_and_units() -> None:
    heading("Range controls and units firewall")
    sign_flip = classify_samples((-1, 0, 1))
    sign_flip = primitive("PASS_RANGE_SIGN_FLIP", sign_flip, "CONSTANT_OUTCOME")
    expect_bool("PASS_RANGE_SIGN_FLIP", sign_flip == "outcome_not_invariant", sign_flip)
    zero_touch = classify_samples((1, 0, 1))
    zero_touch = primitive("PASS_RANGE_ZERO_TOUCH", zero_touch, "CONSTANT_OUTCOME")
    expect_bool("PASS_RANGE_ZERO_TOUCH", zero_touch == "outcome_not_invariant", zero_touch)
    subdominant = classify_samples((sp.Rational(3, 4), sp.Rational(1, 2)))
    subdominant = primitive("PASS_RANGE_SUBDOMINANT", subdominant, "outcome_not_invariant")
    expect_bool("PASS_RANGE_SUBDOMINANT", subdominant == "CONSTANT_OUTCOME", subdominant)

    dim_L = (1, 0, 0)
    dim_E = (2, -2, 1)
    dim_A = tuple(x + y for x, y in zip(dim_E, dim_L))
    a_dim = primitive("PASS_UNITS_A", dim_A, dim_E)
    expect_bool("PASS_UNITS_A", a_dim == (3, -2, 1), {"[A]": a_dim, "target": "E L=M L^3 T^-2"})
    u_dim = primitive("PASS_UNITS_U", tuple(x - y for x, y in zip(dim_A, dim_L)), dim_A)
    expect_bool("PASS_UNITS_U", u_dim == dim_E, {"[U]": u_dim, "target": "E"})
    f_dim = primitive("PASS_UNITS_F", tuple(x - 2 * y for x, y in zip(dim_A, dim_L)), dim_E)
    expect_bool("PASS_UNITS_F", f_dim == (1, -2, 1), {"[F]": f_dim, "target": "E/L"})


def run_landing() -> None:
    heading("Sealed 23040-cell section-4 ladder and landing ownership")
    typed_tokens: tuple[Any, ...] = neutral_tokens(PRODUCTION_FACTS)
    if ACTIVE_MUTATION == "PASS_TYPED_NEUTRAL_FACTS":
        typed_tokens = typed_tokens + (VERDICT,)
    expect_bool("PASS_TYPED_NEUTRAL_FACTS",
                typed_landing_facts(PRODUCTION_FACTS)
                and len(typed_tokens) == len(neutral_tokens(PRODUCTION_FACTS))
                and all(not (isinstance(x, str) and "R1_REQUIRED" in x) for x in typed_tokens),
                typed_tokens)

    truth_exact = primitive("PASS_VERDICT_TOTALITY", TRUTH.exact, False)
    expect_bool("PASS_VERDICT_TOTALITY", truth_exact and all("unclassified" not in key for key in TRUTH.counts),
                {"exact": truth_exact, "landings": tuple(TRUTH.counts)})
    precedence_witness = {
        "qbc": "UNDETERMINED_ANALYTICALLY", "replace_outcome": "POSITIVE_R2",
        "add_outcome": "POSITIVE_R2", "variant": "both",
        "magnitude": "magnitude_no_free_factor", "tier": "tier_A_conditional",
        "internal": False, "all_classes_agree": True, "mixed_range_invariant": True,
    }
    precedence = adjudicate(
        precedence_witness,
        illegal_precedence=ACTIVE_MUTATION == "PASS_VERDICT_PRECEDENCE",
    )
    expect_bool("PASS_VERDICT_PRECEDENCE",
                precedence == "SIGN_EARNED" and PRODUCTION_LANDING == VERDICT,
                {"invariance witness": precedence, "production first match": PRODUCTION_LANDING})

    total = primitive("PASS_TRUTH_TABLE_TOTAL", TRUTH.total, TRUTH.total + 1)
    expect_zero("PASS_TRUTH_TABLE_TOTAL", total - EXPECTED_TOTAL, total)
    counts = dict(TRUTH.counts)
    if ACTIVE_MUTATION == "PASS_TRUTH_TABLE_COUNTS":
        counts[VERDICT] += 1
    expect_bool("PASS_TRUTH_TABLE_COUNTS", counts == EXPECTED_COUNTS and sum(counts.values()) == EXPECTED_TOTAL,
                counts)
    digest_rows = list(TRUTH.rows)
    if ACTIVE_MUTATION == "PASS_TRUTH_TABLE_DIGEST":
        digest_rows[0] += "|MUTATED"
    digest = hashlib.sha256("\n".join(digest_rows).encode("utf-8")).hexdigest()
    expect_bool("PASS_TRUTH_TABLE_DIGEST", digest == COMMITTED_DIGEST,
                {"computed": digest, "committed": COMMITTED_DIGEST})

    alternative_facts = replace(PRODUCTION_FACTS, qbc_status=BCStatus(BCClass.FIXED_SOURCE))
    alternative_tuple = production_tuple(alternative_facts)
    alternative_landing = adjudicate(alternative_tuple)
    ownership_facts = primitive("PASS_LANDING_OWNERSHIP", PRODUCTION_FACTS, alternative_facts)
    ownership_tuple = production_tuple(ownership_facts)
    ownership_landing = adjudicate(ownership_tuple)
    injected = neutral_tokens(PRODUCTION_FACTS) + (PRODUCTION_LANDING,)
    ownership_ok = (
        landing_ownership_guard(ownership_facts, PRODUCTION_LANDING)
        and not landing_ownership_guard(PRODUCTION_FACTS, PRODUCTION_LANDING, injected)
        and ownership_tuple.overall_outcome == production_tuple(ownership_facts).overall_outcome
        and ownership_tuple.all_classes_agree == production_tuple(ownership_facts).all_classes_agree
        and alternative_landing == "R1_REQUIRED(sign_and_magnitude)"
        and alternative_landing != VERDICT
        and ownership_landing == PRODUCTION_LANDING
    )
    expect_bool("PASS_LANDING_OWNERSHIP", ownership_ok,
                {"upstream_qbc": ownership_facts.qbc_status.kind.value,
                 "recomputed_overall": ownership_tuple.overall_outcome,
                 "recomputed_all_classes_agree": ownership_tuple.all_classes_agree,
                 "fresh_landing": ownership_landing,
                 "named_different_landing": alternative_landing,
                 "injected_landing_rejected": not landing_ownership_guard(PRODUCTION_FACTS, PRODUCTION_LANDING, injected)})

    production_facts = primitive("PASS_PRODUCTION_LANDING", PRODUCTION_FACTS, alternative_facts)
    live_tuple = production_tuple(production_facts)
    live_landing = adjudicate(live_tuple)
    expect_bool("PASS_PRODUCTION_LANDING",
                live_tuple.overall_outcome == "outcome_not_invariant"
                and not live_tuple.all_classes_agree
                and live_landing == VERDICT
                and not live_tuple.internal,
                {"tuple": live_tuple, "landing": live_landing})

    upstream_vocabulary = {
        *(status.value for status in BCClass), *(outcome.value for outcome in Outcome),
        *(variant.value for variant in Variant), *(magnitude.value for magnitude in Magnitude),
        *(tier.value for tier in Tier), "missing parent-throat/boundary closure",
    }
    target_token = primitive("PASS_TARGET_BLINDNESS", "outcome_not_invariant", "SIGN_EARNED")
    expect_bool("PASS_TARGET_BLINDNESS",
                target_token in upstream_vocabulary
                and all("targetsign" not in str(expr).lower() for expr in COEFFICIENTS.values())
                and "R1_REQUIRED" not in "|".join(upstream_vocabulary),
                target_token)
    print(f"      cells={TRUTH.total}; digest={TRUTH.digest}; first production match={PRODUCTION_LANDING}")


def run_manifest() -> None:
    heading("Canonical source-to-stage predicate manifest")
    manifest = primitive("PASS_SOURCE_PREDICATE_MANIFEST", SOURCE_MANIFEST, SOURCE_MANIFEST[:-1])
    manifest = primitive(
        "MANIFEST_MISDISPOSITION",
        manifest,
        tuple(
            manifest_entry(identifier, "SCOPED_OUT", owner)
            if identifier == "AMEND_REPLACE"
            else row
            for row in SOURCE_MANIFEST
            for identifier, _disposition, owner in (row,)
        ),
    )
    identifiers = tuple(row[0] for row in manifest)
    source_part = identifiers[:len(SOURCE_TOOTH_IDS)]
    extra_part = identifiers[len(SOURCE_TOOTH_IDS):]
    dispositions = {row[1] for row in manifest}
    partition_counts = dict(sorted(Counter(row[1] for row in manifest).items()))
    manifest_digest = manifest_sha256(manifest)
    condition = (
        source_part == SOURCE_TOOTH_IDS
        and extra_part == EXTRA_SOURCE_CLAIMS
        and len(identifiers) == len(set(identifiers)) == 44
        and dispositions == {"PRESERVED", "REPLACED_BY_STRONGER", "SCOPED_OUT"}
        and all(row[2].startswith("STAGE") for row in manifest)
        and partition_counts == EXPECTED_MANIFEST_COUNTS
        and manifest_digest == COMMITTED_MANIFEST_DIGEST
    )
    expect_bool("PASS_SOURCE_PREDICATE_MANIFEST", condition,
                {"entries": len(manifest), "dispositions": sorted(dispositions),
                 "partition": partition_counts, "digest": manifest_digest})
    print(f"      entries={len(manifest)}; partition={partition_counts}; digest={manifest_digest}")


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in MUTATION_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage032_electric_sign_bc_ensembles_landing SymPy audit")
    print("CONSUMES_STAGE031={m,m_gg=B_eff*z_g^2/D,detm=z_g^2/D,S_gg,A=m_gg*C,1/R^2,NONZERO_HA_REQUIRES_CORE_HOLDER}")
    print("CONSUMES_STAGE030_TRANSITIVELY={D,D*=7/4}")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")

    run_consumption_and_scope()
    run_ensembles()
    run_source_blocks()
    run_ranges_and_units()
    run_landing()
    run_manifest()

    print("")
    print("INTERNAL_INCONSISTENCY=none")
    print("CO_BLOCKER=R1_REQUIRED(magnitude)")
    print("RESOLVER=SIM_DEFERRED(parent-throat boundary functional/barrier/map s->h_A)")
    print(f"VERDICT_TOKEN: {VERDICT}")
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified stage032 far-field ensembles and sealed R1 landing")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage032 audit did not close ({exc.predicate})")
        raise SystemExit(1)
