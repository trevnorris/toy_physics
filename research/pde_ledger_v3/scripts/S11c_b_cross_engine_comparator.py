#!/usr/bin/env python3
"""S11c-b T7 cross-engine operand comparator.

This program is a measurement instrument.  It reads the two committed S11c-b
tag streams, applies only the closed mechanical reconciliation map, and emits
both operands plus their typed recursive A-minus-B residual.  Algebraic and
structural disagreements are output data and do not change the exit status.

The SymPy engine is deliberately never imported or run.  In particular,
S11c_b_exports.py is not used as a tag stream.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# S11c-a is the verified mechanical base.  Its parsers, association/text
# containers, closed name map, capture-safe BoundIntegral canonicalizer, and
# typed recursive residual are reused; none of its family classification or
# run-status machinery is imported.
import S11c_a_cross_engine_comparator as s11ca  # noqa: E402
from S11b_cross_engine_comparator import (  # noqa: E402
    Association,
    TextAtom,
    _convert_parsed_containers,
    parse_sympy_payload,
    parse_wolfram_payload,
    residual,
)


DEFAULT_PY = SCRIPT_DIR / "out" / "S11c_b_brane_operator_sympy_audit.out"
DEFAULT_WL = (
    SCRIPT_DIR.parent
    / "mathematica"
    / "out"
    / "S11c_b_brane_operator_mathematica_audit.out"
)

PY_TAG_LINE = re.compile(r"^(PY_S11CB_[A-Z0-9_]+):\s?(.*)$")
WL_TAG_LINE = re.compile(r"^(WL_S11CB_[A-Z0-9_]+):")

AXIS_ORDER = (
    "OBJECT",
    "BRANCH",
    "DENSITY",
    "SECTOR",
    "SOURCE",
    "DIRECTION",
    "DOF",
)
Key = tuple[tuple[str, str], ...]

BRANCHES = {"LAB_HELD", "MATERIAL_ADVECTED"}
DENSITIES = {"RHO4_CONSTANT", "RHOBR_CONSTANT"}
SOURCES = {"W_BG", "MU_R_BG"}
DIRECTIONS = {"1", "2", "3"}
SECTOR_RENAME = {
    "TRANSVERSE_TO_THICKNESS": "TRANSVERSE_TO_THICKNESS",
    "TRANSVERSE_TO_THETA_EW_UL": "TRANSVERSE_TO_THICKNESS",
    "THICKNESS_TO_TRANSVERSE": "THICKNESS_TO_TRANSVERSE",
    "THETA_EW_UL_TO_TRANSVERSE": "THICKNESS_TO_TRANSVERSE",
}
MAX_RATIONAL_CANON_RAW_CHARS = 150_000
OBJECT_RENAME = {
    "ENERGY_BASIS_VARIABLE": "ENERGY_BASIS",
}

ENERGY_FAMILIES = {
    "ENERGY_BASIS_VARIABLE",
    "ENERGY_BASIS_COUNT",
    "ENERGY_BASIS_NEW_INVARIANTS",
    "ENERGY_BASIS_OMISSIONS",
}
ADMISSIBILITY_FAMILIES = {
    "ADMISSIBILITY_OPERATOR_OPERAND",
    "ADMISSIBILITY_SUPPORT_OPERAND",
    "ADMISSIBILITY_RESIDUAL",
}
NESTED_CONTROL_PREFIXES = (
    "REP_INVARIANCE_",
    "CONTROL_INDEPENDENCE_",
    "UNIFORM_LIMIT_",
)


class InputError(ValueError):
    """A committed stream does not satisfy the declared tag grammar."""


class AxisError(ValueError):
    """A case key is not a single-valued typed-axis key."""


@dataclass
class ParsedCase:
    key: Key
    value: object | None
    raw: str
    error: str | None = None
    compare_value: object | None = None
    note: str | None = None
    context: tuple[tuple[str, object], ...] = ()
    loader: Callable[[], object] | None = field(default=None, repr=False)
    reduce_divergence: bool = False


@dataclass
class FamilyCases:
    py: list[ParsedCase] = field(default_factory=list)
    wl: list[ParsedCase] = field(default_factory=list)
    extraction_notes: list[str] = field(default_factory=list)


@dataclass
class Accounting:
    join: int = 0
    py_only: int = 0
    wl_only: int = 0
    duplicate_key: int = 0
    parse_failed: int = 0
    axis_set_mismatch: int = 0
    reasons: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class RawCAS:
    """An as-emitted operand whose outer container already blocks subtraction."""

    engine: str
    raw: str


@dataclass(frozen=True)
class SymbolicDifference:
    """Exact unevaluated A-minus-B when global rational canon is oversized."""

    expression: sp.Expr
    reason: str


def load_py(path: Path) -> dict[str, str]:
    """Load the mandatory one-line PY srepr stream without running its engine."""
    if not path.is_file():
        raise InputError(f"PY input does not exist: {path}")
    output: dict[str, str] = {}
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            match = PY_TAG_LINE.fullmatch(line.rstrip("\r\n"))
            if match is None:
                raise InputError(f"{path}:{line_number}: not one PY S11c-b tagged row")
            full_tag, payload = match.groups()
            name = full_tag.removeprefix("PY_S11CB_")
            if name in output:
                raise InputError(f"{path}:{line_number}: duplicate PY tag {name}")
            output[name] = payload
    if not output:
        raise InputError(f"{path}: no PY S11c-b tags")
    return output


def load_wl(path: Path) -> dict[str, str]:
    """Load WL tags, reassembling a payload through the next tag line."""
    if not path.is_file():
        raise InputError(f"WL input does not exist: {path}")
    output: dict[str, str] = {}
    current_name: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        nonlocal current_name, chunks
        if current_name is None:
            return
        if current_name in output:
            raise InputError(f"{path}: duplicate WL tag {current_name}")
        output[current_name] = "".join(chunks).strip()

    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            clean = line.rstrip("\r\n")
            match = WL_TAG_LINE.match(clean)
            if match is not None:
                finish()
                full_tag = match.group(1)
                current_name = full_tag.removeprefix("WL_S11CB_")
                chunks = [clean.split(":", 1)[1]]
            elif current_name is None:
                if clean.strip():
                    raise InputError(f"{path}:{line_number}: content before first WL tag")
            else:
                chunks.append(clean)
    finish()
    if not output:
        raise InputError(f"{path}: no WL S11c-b tags")
    return output


# Re-express the verified delimiter-aware structural readers under this
# module's public API; tests and future family extractors use these names.
split_top = s11ca.split_top
arrow = s11ca.arrow
wl_assoc_pairs = s11ca.wl_assoc_pairs
wl_field = s11ca.wl_field
py_tuple_args = s11ca.py_tuple_args
py_pair = s11ca.py_pair
py_top_pairs = s11ca.py_top_pairs
py_text = s11ca.py_text
py_key_tokens = s11ca.py_key_tokens
py_field = s11ca.py_field


def make_key(items: Iterable[tuple[str, str]]) -> Key:
    """Build one fixed-order typed key, rejecting duplicate or unknown axes."""
    decoded: dict[str, str] = {}
    for axis, original in items:
        value = str(original)
        if axis in decoded:
            raise AxisError(f"two values for axis {axis}: {decoded[axis]} and {value}")
        if axis == "OBJECT":
            value = OBJECT_RENAME.get(value, value)
            if not value:
                raise AxisError("empty OBJECT value")
        elif axis == "BRANCH" and value not in BRANCHES:
            raise AxisError(f"invalid BRANCH value {value}")
        elif axis == "DENSITY" and value not in DENSITIES:
            raise AxisError(f"invalid DENSITY value {value}")
        elif axis == "SECTOR":
            if value not in SECTOR_RENAME:
                raise AxisError(f"invalid SECTOR value {value}")
            value = SECTOR_RENAME[value]
        elif axis == "SOURCE" and value not in SOURCES:
            raise AxisError(f"invalid SOURCE value {value}")
        elif axis == "DIRECTION" and value not in DIRECTIONS:
            raise AxisError(f"invalid DIRECTION value {value}; expected 1, 2, or 3")
        elif axis == "DOF" and not value:
            raise AxisError("empty DOF value")
        decoded[axis] = value
    unknown = set(decoded) - set(AXIS_ORDER)
    if unknown:
        raise AxisError(f"unknown axes {sorted(unknown)}")
    return tuple((axis, decoded[axis]) for axis in AXIS_ORDER if axis in decoded)


def key_dict(key: Key) -> dict[str, str]:
    return dict(key)


def key_replace(key: Key, *, remove: Sequence[str] = (), **updates: str) -> Key:
    items = {axis: value for axis, value in key if axis not in set(remove)}
    items.update(updates)
    return make_key(items.items())


def decode_py_key(raw: str, schema: Sequence[str]) -> Key:
    tokens = py_key_tokens(raw)
    if len(tokens) != len(schema):
        raise AxisError(f"PY schema {tuple(schema)} got {len(tokens)} tokens: {tokens}")
    return make_key(zip(schema, tokens))


def decode_wl_key(raw: str, schema: Sequence[str]) -> Key:
    tokens = tuple(raw.split("|"))
    if len(tokens) != len(schema):
        raise AxisError(f"WL schema {tuple(schema)} got {len(tokens)} tokens: {tokens}")
    typed: list[tuple[str, str]] = []
    for axis, token in zip(schema, tokens):
        if axis == "DIRECTION":
            if not token.startswith("DIRECTION_"):
                raise AxisError(f"untyped WL DIRECTION token {token}")
            token = token.removeprefix("DIRECTION_")
        typed.append((axis, token))
    return make_key(typed)


def py_string_map(raw: str) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    for entry in py_tuple_args(raw):
        key_raw, value_raw = py_pair(entry)
        output.append((py_text(key_raw), value_raw))
    return output


def py_residual_map(raw: str) -> list[tuple[str, str]]:
    """Decode a residual map whose emitted key is itself a key residual."""
    output: list[tuple[str, str]] = []
    for entry in py_tuple_args(raw):
        key_residual, value_raw = py_pair(entry)
        tokens = py_key_tokens(key_residual)
        if not tokens or len(set(tokens)) != 1:
            raise InputError(f"residual map key has no unique textual label: {tokens}")
        output.append((tokens[0], value_raw))
    return output


# ---------- closed name/CAS-form map ----------

# The inherited S11c-a map is applied first.  These are spelling-only S11c-b
# extensions for heads present in the measured WL payload; no energy-basis
# representative identity is encoded here.
EXTRA_HEAD = {
    "widthBackground": "W_bg",
    "modulusBackground": "mu_R_bg",
    "thetaField": "theta",
    "eWave": "e_W",
    "eWField": "e_W",
    "eWBackground": "e_W_bg",
    "virtualEw": "delta_v_e_W",
    "longitudinalTrialPotential": "phi_L",
    "longitudinalTestPotential": "psi_L_s11cb",
    "thetaTrial": "theta_probe",
    "ewTrial": "e_W_probe",
    "thetaTest": "v_theta_s11cb",
    "ewTest": "v_e_W_s11cb",
    "transverseTrialPotentialOne": "A_T_s11cb_1",
    "transverseTrialPotentialTwo": "A_T_s11cb_2",
    "transverseTrialPotentialThree": "A_T_s11cb_3",
    "transverseTestPotentialOne": "A_T_s11cb_1",
    "transverseTestPotentialTwo": "A_T_s11cb_2",
    "transverseTestPotentialThree": "A_T_s11cb_3",
    "forceHoldOne": "f_hold_u_1_0",
    "forceHoldTwo": "f_hold_u_2_0",
    "forceHoldThree": "f_hold_u_3_0",
    "forceHoldTheta": "f_hold_theta_0",
    "forceHoldEw": "f_hold_e_W_0",
    "tractionHoldUpperOne": "t_hold_plus_0_1",
    "tractionHoldUpperTwo": "t_hold_plus_0_2",
    "tractionHoldUpperThree": "t_hold_plus_0_3",
    "tractionHoldUpperW": "t_hold_plus_0_4",
    "tractionHoldLowerOne": "t_hold_minus_0_1",
    "tractionHoldLowerTwo": "t_hold_minus_0_2",
    "tractionHoldLowerThree": "t_hold_minus_0_3",
    "tractionHoldLowerW": "t_hold_minus_0_4",
}
EXTRA_SYMBOL = {
    "WZero": "W_0",
    "LWidth": "L_W",
    "frequency": "omega",
    "w1ProfileZero": "w1_profile",
    "m1ProfileZero": "m1_profile",
    "w1JetOne": "w1_profile_d1",
    "w1JetTwo": "w1_profile_d2",
    "w1JetThree": "w1_profile_d3",
    "m1JetOne": "m1_profile_d1",
    "m1JetTwo": "m1_profile_d2",
    "m1JetThree": "m1_profile_d3",
    "kappaW": "kappa_W",
    "gradientThetaEwCoefficient": "kappa_theta_W",
    "gradientThetaCoefficient": "kappa_theta",
    "modulusDivU": "B_div",
    "cCross": "C",
    "thetaDivUCoefficient": "G_theta_u",
    "ewDivUCoefficient": "G_W_u",
    "kW": "k_W",
    "muW": "mu_W",
    "modulusShear": "mu_S",
    "backgroundOrder": "__ONE__",
    "waveOrder": "__ONE__",
    "virtualOrder": "__ONE__",
    "wave": "__ONE__",
    "order": "__ONE__",
}


def _extra_basic(value: sp.Basic) -> sp.Basic:
    def mapped_head(node: AppliedUndef) -> sp.Basic:
        name = node.func.__name__
        if name in {"widthProfileJet", "modulusProfileJet"}:
            if not all(item.is_Integer for item in node.args):
                return node
            base = "w1_profile" if name == "widthProfileJet" else "m1_profile"
            suffix = "".join(
                f"d{index + 1}" * int(order)
                for index, order in enumerate(node.args)
            )
            return sp.Symbol(base + ("_" + suffix if suffix else ""))
        if "XJETX" in name:
            head, suffix = name.split("XJETX", 1)
            target = EXTRA_HEAD.get(head)
            if target is not None:
                return sp.Symbol(s11ca.canon_jet_name(target + "_" + suffix.replace("X", "_")))
        target = EXTRA_HEAD.get(name)
        if target is not None:
            return sp.Symbol(target)
        return node

    value = value.replace(
        lambda node: isinstance(node, AppliedUndef)
        and (
            node.func.__name__ in EXTRA_HEAD
            or node.func.__name__ in {"widthProfileJet", "modulusProfileJet"}
            or "XJETX" in node.func.__name__
        ),
        mapped_head,
    )
    replacements: dict[sp.Symbol, sp.Basic] = {}
    for symbol in value.atoms(sp.Symbol):
        name = symbol.name
        if name.startswith("THOLDX"):
            replacements[symbol] = sp.Symbol(name.removeprefix("THOLDX"))
            continue
        if "XJETX" in name:
            head, suffix = name.split("XJETX", 1)
            target = EXTRA_HEAD.get(head)
            if target is not None:
                replacements[symbol] = sp.Symbol(
                    s11ca.canon_jet_name(target + "_" + suffix.replace("X", "_"))
                )
                continue
        target = EXTRA_SYMBOL.get(name)
        if target == "__ONE__":
            replacements[symbol] = sp.Integer(1)
        elif target is not None:
            replacements[symbol] = sp.Symbol(target)
    return value.xreplace(replacements) if replacements else value


def _extra_value(value: object) -> object:
    if isinstance(value, Association):
        return Association(tuple((key, _extra_value(item)) for key, item in value.entries))
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(_extra_basic)
    if isinstance(value, tuple):
        return tuple(_extra_value(item) for item in value)
    if isinstance(value, Relational):
        return type(value)(
            _extra_value(value.lhs), _extra_value(value.rhs), evaluate=False
        )
    if isinstance(value, sp.Basic):
        return _extra_basic(value)
    return value


def preprocess_wl(raw: str) -> str:
    # Inactive[Div] is retained as an explicit held operator.  The adjointness
    # quotient canonicalizer below removes it only when it is the declared
    # compact-support total-divergence term.
    return s11ca.preprocess_wl(raw.replace("Inactive[Div]", "HeldDiv"))


def _parse_py_cached(raw: str, coefficient: bool, branch: str | None) -> object:
    value = parse_sympy_payload(raw)
    protected = {
        symbol: sp.Symbol("THOLDX" + symbol.name)
        for symbol in value.atoms(sp.Symbol)
        if symbol.name.startswith("t_hold_")
    }
    if protected:
        value = value.xreplace(protected)
    if coefficient:
        value = s11ca.coeff_epsilon(value)
    value = s11ca.canonical_value(value, "PY", branch)
    return _convert_parsed_containers(_extra_value(value))


def _parse_wl_cached(raw: str, branch: str | None) -> object:
    value = parse_wolfram_payload(preprocess_wl(raw))
    value = s11ca.canonical_value(value, "WL", branch)
    return _convert_parsed_containers(_extra_value(value))


def parse_py_value(raw: str, key: Key, *, coefficient: bool = True) -> object:
    return _parse_py_cached(raw, coefficient, key_dict(key).get("BRANCH"))


def parse_wl_value(raw: str, key: Key) -> object:
    return _parse_wl_cached(raw, key_dict(key).get("BRANCH"))


# Public aliases document that the S11c-a capture-safe integral canon is the
# one used here, including binder and bounds.
BOUND_INTEGRAL = s11ca.BOUND_INTEGRAL
BOUND_BINDER = s11ca.BOUND_BINDER
CAPTURED_FREE_BINDER = s11ca.CAPTURED_FREE_BINDER
canonical_integrals = s11ca.canonical_integrals


# ---------- compact-support total-divergence quotient ----------

DEPENDENT_BASES = {
    "u_T_1",
    "u_T_2",
    "u_T_3",
    "A_T_s11cb_1",
    "A_T_s11cb_2",
    "A_T_s11cb_3",
    "phi_L",
    "psi_L_s11cb",
    "theta_probe",
    "e_W_probe",
    "v_theta_s11cb",
    "v_e_W_s11cb",
}
DIFFERENTIABLE_BASES = DEPENDENT_BASES | {
    "W_bg",
    "mu_R_bg",
    "w1_profile",
    "m1_profile",
}
JET_SUFFIX = re.compile(r"^(.*)_(d[123](?:d[123])*)$")


def split_jet_name(name: str) -> tuple[str, tuple[int, ...]] | None:
    match = JET_SUFFIX.fullmatch(name)
    if match is not None:
        base, suffix = match.groups()
        if base in DIFFERENTIABLE_BASES:
            return base, tuple(int(item) for item in re.findall(r"d([123])", suffix))
    if name in DIFFERENTIABLE_BASES:
        return name, ()
    return None


def jet_symbol(base: str, directions: Sequence[int]) -> sp.Symbol:
    ordered = tuple(sorted(directions))
    suffix = "" if not ordered else "_" + "".join(f"d{item}" for item in ordered)
    return sp.Symbol(base + suffix)


def formal_total_derivative(expression: sp.Expr, direction: int) -> sp.Expr:
    terms: list[sp.Expr] = []
    for symbol in expression.free_symbols:
        decoded = split_jet_name(symbol.name)
        if decoded is None:
            continue
        base, directions = decoded
        partial = sp.diff(expression, symbol)
        if partial != 0:
            terms.append(partial * jet_symbol(base, (*directions, direction)))
    return sp.Add(*terms)


def _drop_held_divergences(expression: sp.Expr) -> sp.Expr:
    return expression.replace(
        lambda node: isinstance(node, AppliedUndef)
        and node.func.__name__ == "HeldDiv",
        lambda _node: sp.S.Zero,
    )


def modulo_total_divergence(value: object) -> object:
    """Return the spatial Euler signature, a canonical density quotient.

    The Euler operator annihilates total in-plane divergences for the declared
    compact-support trial/test fields.  Returning every dependent-field Euler
    derivative retains non-divergence density content; it does not assert
    adjointness or attach a target value.
    """
    if not isinstance(value, sp.Expr) or isinstance(value, Relational):
        return TextAtom(f"UNREDUCED_NONSCALAR_{type(value).__name__}")
    expression = _drop_held_divergences(value)
    symbols_by_base: dict[str, list[tuple[sp.Symbol, tuple[int, ...]]]] = defaultdict(list)
    for symbol in expression.free_symbols:
        decoded = split_jet_name(symbol.name)
        if decoded is not None and decoded[0] in DEPENDENT_BASES:
            symbols_by_base[decoded[0]].append((symbol, decoded[1]))
    if not symbols_by_base:
        return expression
    entries: list[tuple[str, object]] = []
    for base in sorted(symbols_by_base):
        euler = sp.S.Zero
        for symbol, directions in symbols_by_base[base]:
            term = sp.diff(expression, symbol)
            for direction in directions:
                term = -formal_total_derivative(term, direction)
            euler += term
        entries.append((base, sp.expand(euler)))
    return Association(tuple(entries))


# ---------- case construction ----------

def build_case(
    engine: str,
    key: Key,
    raw: str,
    parser: Callable[[], object],
    *,
    note: str | None = None,
    reduce_divergence: bool = False,
) -> ParsedCase:
    # Parsing is intentionally lazy.  Several measured control families hold
    # tens of megabytes of srepr in one tag; materializing all cases before the
    # first comparison needlessly multiplies peak memory.
    return ParsedCase(
        key,
        None,
        raw,
        note=note,
        loader=parser,
        reduce_divergence=reduce_divergence,
    )


def materialize(case: ParsedCase) -> bool:
    """Parse one case once; return whether this call introduced an error."""
    if case.loader is None:
        return False
    loader = case.loader
    case.loader = None
    try:
        case.value = loader()
        case.compare_value = (
            modulo_total_divergence(case.value)
            if case.reduce_divergence
            else case.value
        )
        if case.reduce_divergence:
            case.context = (("DIVERGENCE_REDUCED", case.compare_value),)
        return False
    except Exception as error:
        case.error = f"{type(error).__name__}: {error}"
        return True


def release_case(case: ParsedCase | None) -> None:
    """Drop a rendered operand before the next multi-megabyte case is parsed."""
    if case is None:
        return
    case.value = None
    case.compare_value = None
    case.context = ()


def py_case(
    key: Key,
    raw: str,
    *,
    coefficient: bool = True,
    note: str | None = None,
    reduce_divergence: bool = False,
) -> ParsedCase:
    return build_case(
        "PY",
        key,
        raw,
        lambda: parse_py_value(raw, key, coefficient=coefficient),
        note=note,
        reduce_divergence=reduce_divergence,
    )


def wl_case(
    key: Key,
    raw: str,
    *,
    note: str | None = None,
    reduce_divergence: bool = False,
) -> ParsedCase:
    return build_case(
        "WL",
        key,
        raw,
        lambda: parse_wl_value(raw, key),
        note=note,
        reduce_divergence=reduce_divergence,
    )


def py_sum_case(
    key: Key,
    raws: Sequence[str],
    *,
    note: str | None = None,
    reduce_divergence: bool = False,
) -> ParsedCase:
    rendered = " + ".join(raws)

    def parse_sum() -> object:
        values = [parse_py_value(raw, key, coefficient=True) for raw in raws]
        if not all(isinstance(value, sp.Expr) for value in values):
            raise InputError("weak-pairing sum contains a non-scalar leaf")
        return sp.Add(*values)

    return build_case(
        "PY", key, rendered, parse_sum, note=note, reduce_divergence=reduce_divergence
    )


def raw_control_case(engine: str, key: Key, raw: str, note: str | None) -> ParsedCase:
    """Surface a bespoke control container without materializing its huge leaves.

    The measured PY and WL control objects use different outer schemas (for
    example operator rows versus balance maps).  The inherited typed residual
    therefore stops at the outer structure and never reaches an algebraic
    leaf.  Recording that outer signature is mathematically equivalent to the
    recursive route while avoiding expansion of tens of megabytes that cannot
    be subtracted as scalars.
    """
    try:
        if engine == "PY":
            try:
                signature = "ASSOCIATION_KEYS:" + ",".join(
                    key_name for key_name, _value in py_string_map(raw)
                )
            except Exception:
                signature = f"TUPLE_ARITY:{len(py_tuple_args(raw))}"
        else:
            signature = "ASSOCIATION_KEYS:" + ",".join(
                key_name for key_name, _value in wl_assoc_pairs(raw)
            )
        signature = f"{engine}:{signature}"
        return ParsedCase(
            key,
            RawCAS(engine, raw),
            raw,
            compare_value=TextAtom(signature),
            note="outer_container_blocks_scalar_subtraction; raw_operand_emitted"
            + (f" | {note}" if note else ""),
        )
    except Exception as error:
        return ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}", note=note)


def extraction_error(engine: str, family: str, raw: str, error: Exception) -> ParsedCase:
    key = make_key((("OBJECT", f"{engine}_EXTRACTION_ERROR"),))
    return ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}")


# ---------- per-family extractors ----------

def extract_energy(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    quotient_note = "non_unique_energy_quotient_compared_as_emitted_no_representative_fold"
    if engine == "PY":
        for key_raw, record_raw in py_top_pairs(payload):
            key = decode_py_key(key_raw, ("BRANCH",))
            raw = py_field(record_raw)
            output.append(py_case(key, raw, coefficient=False, note=quotient_note))
        return output

    selected = {
        "ENERGY_BASIS_VARIABLE": "REPRESENTATIVE",
        "ENERGY_BASIS_COUNT": "COUNT_OPERAND",
        "ENERGY_BASIS_OMISSIONS": "OMITTED_FROM_CARRIED_INPUT",
    }
    for branch, record_raw in wl_assoc_pairs(payload):
        key = make_key((("BRANCH", branch),))
        raw = record_raw if family == "ENERGY_BASIS_NEW_INVARIANTS" else wl_field(
            record_raw, selected[family]
        )
        output.append(wl_case(key, raw, note=quotient_note))
    return output


def extract_slab(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for key_raw, record_raw in py_top_pairs(payload):
            base_key = decode_py_key(key_raw, ("BRANCH", "DENSITY"))
            values = dict(py_string_map(py_field(record_raw)))
            row_map = {
                "U_BODY_BALANCE": ("U_MOMENTUM_ROWS", "U"),
                "THETA_BALANCE": ("MU_THETA", "THETA"),
                "E_W_BALANCE": ("THICKNESS_ROW", "E_W"),
            }
            for source, (canonical, dof) in row_map.items():
                forms = dict(py_string_map(values[source]))
                for form, raw in forms.items():
                    object_name = canonical if form == "EXPANDED" else f"{source}.{form}"
                    key = key_replace(base_key, OBJECT=object_name, DOF=dof)
                    output.append(py_case(key, raw))
            key = key_replace(base_key, OBJECT="MASS_EVOLUTION_ROW", DOF="MASS")
            output.append(py_case(key, values["ADVECTIVE_MASS_OPERAND"]))
            key = key_replace(base_key, OBJECT="FACE_SHAPE_SUBSTRATE", DOF="ALL")
            output.append(py_case(key, values["FACE_FLUX_BOUNDARY_OPERANDS"]))
            key = key_replace(base_key, OBJECT="MU_THETA_FACE_BINDING", DOF="MU_THETA")
            output.append(py_case(key, values["MU_THETA_FACE_BINDING"]))
        return output

    for key_raw, record_raw in wl_assoc_pairs(payload):
        base_key = decode_wl_key(key_raw, ("BRANCH", "DENSITY"))
        for row_name, body_raw in wl_assoc_pairs(wl_field(record_raw, "ROWS")):
            dof = {
                "U_MOMENTUM_ROWS": "U",
                "THICKNESS_ROW": "E_W",
                "MASS_EVOLUTION_ROW": "MASS",
                "CENTER_FACE_GENERALIZED_ROW": "CENTER_FACE",
            }[row_name]
            key = key_replace(base_key, OBJECT=row_name, DOF=dof)
            output.append(wl_case(key, wl_field(body_raw, "EXPRESSION")))
        source = wl_field(record_raw, "DIVERGENCE_FORM_SOURCE", "EXPRESSION")
        for name, raw in wl_assoc_pairs(source):
            canonical = "MU_THETA" if name == "MU_THETA" else f"DIVERGENCE_FORM_SOURCE.{name}"
            dof = "THETA" if name == "MU_THETA" else name
            key = key_replace(base_key, OBJECT=canonical, DOF=dof)
            output.append(wl_case(key, raw))
        key = key_replace(base_key, OBJECT="VIRTUAL_CONSTRAINT", DOF="E_W")
        output.append(wl_case(key, wl_field(record_raw, "VIRTUAL_CONSTRAINT", "EXPRESSION")))
        key = key_replace(base_key, OBJECT="FACE_SHAPE_SUBSTRATE", DOF="ALL")
        output.append(wl_case(key, wl_field(record_raw, "FACE_SHAPE_SUBSTRATE")))
    return output


def extract_term_origins(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for key_raw, record_raw in py_top_pairs(payload):
            base_key = decode_py_key(key_raw, ("BRANCH", "DENSITY"))
            for origin, raw in py_string_map(py_field(record_raw)):
                key = key_replace(base_key, OBJECT=origin)
                output.append(py_case(key, raw))
        return output
    for key_raw, record_raw in wl_assoc_pairs(payload):
        base_key = decode_wl_key(key_raw, ("BRANCH", "DENSITY"))
        for origin, raw in wl_assoc_pairs(wl_field(record_raw, "EXPRESSION")):
            key = key_replace(base_key, OBJECT=origin)
            output.append(wl_case(key, raw))
    return output


def extract_mu(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for key_raw, record_raw in py_top_pairs(payload):
            key = decode_py_key(key_raw, ("BRANCH", "DENSITY"))
            items = py_tuple_args(py_field(record_raw))
            if len(items) != 2:
                raise InputError("PY MU_THETA_OPERATOR VALUE is not (reserved-name, expression)")
            output.append(py_case(key, items[1]))
        return output
    for key_raw, record_raw in wl_assoc_pairs(payload):
        key = decode_wl_key(key_raw, ("BRANCH", "DENSITY"))
        raw = wl_field(record_raw, "VARIABLE_COEFFICIENT_VARIATIONAL_OPERAND", "EXPRESSION")
        output.append(wl_case(key, raw))
    return output


def _py_pairing_raws(raw: str) -> list[str]:
    items = py_tuple_args(raw)
    if len(items) != 2:
        raise InputError("PY variational pairing entry is not (label, block)")
    return [value for _name, value in py_string_map(items[1])]


def extract_coupling(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for key_raw, record_raw in py_top_pairs(payload):
            base_key = decode_py_key(key_raw, ("BRANCH", "DENSITY"))
            values = dict(py_string_map(py_field(record_raw)))
            adjoint = py_tuple_args(values["VARIATIONAL_ADJOINTNESS"])
            if len(adjoint) != 3:
                raise InputError("PY VARIATIONAL_ADJOINTNESS is not its declared three-slot record")
            forward_raws = _py_pairing_raws(adjoint[1])
            reverse_raws = _py_pairing_raws(adjoint[2])
            for sector, raws in (
                ("TRANSVERSE_TO_THICKNESS", forward_raws),
                ("THICKNESS_TO_TRANSVERSE", reverse_raws),
            ):
                key = key_replace(base_key, SECTOR=sector)
                output.append(py_sum_case(key, raws))
            for object_name, sector, raws in (
                ("ADJOINTNESS_OPERAND_FORWARD", "TRANSVERSE_TO_THICKNESS", forward_raws),
                ("ADJOINTNESS_OPERAND_REVERSE", "THICKNESS_TO_TRANSVERSE", reverse_raws),
            ):
                key = key_replace(base_key, OBJECT=object_name, SECTOR=sector)
                output.append(py_sum_case(key, raws, reduce_divergence=True))
            residual_raws = [*forward_raws, *(f"Mul(Integer(-1), {raw})" for raw in reverse_raws)]
            key = key_replace(base_key, OBJECT="ADJOINTNESS_RESIDUAL")
            output.append(
                py_sum_case(
                    key,
                    residual_raws,
                    note="raw_density_and_spatial_Euler_signature_emitted",
                    reduce_divergence=True,
                )
            )
            key = key_replace(base_key, OBJECT="ADJOINTNESS_RELATION")
            output.append(py_case(key, adjoint[0], coefficient=False))
        return output

    for key_raw, record_raw in wl_assoc_pairs(payload):
        base_key = decode_wl_key(key_raw, ("BRANCH", "DENSITY"))
        expression = wl_field(record_raw, "EXPRESSION")
        for wl_sector in ("TRANSVERSE_TO_THETA_EW_UL", "THETA_EW_UL_TO_TRANSVERSE"):
            key = key_replace(base_key, SECTOR=wl_sector)
            raw = wl_field(
                expression,
                wl_sector,
                "WEAK_PAIRING",
                "PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP",
            )
            output.append(wl_case(key, raw))
        for object_name, wl_sector in (
            ("ADJOINTNESS_OPERAND_FORWARD", "TRANSVERSE_TO_THETA_EW_UL"),
            ("ADJOINTNESS_OPERAND_REVERSE", "THETA_EW_UL_TO_TRANSVERSE"),
        ):
            key = key_replace(base_key, OBJECT=object_name, SECTOR=wl_sector)
            raw = wl_field(
                expression,
                object_name,
                "PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP",
            )
            output.append(wl_case(key, raw, reduce_divergence=True))
        key = key_replace(base_key, OBJECT="ADJOINTNESS_RESIDUAL")
        raw = wl_field(
            expression,
            "ADJOINTNESS_RESIDUAL",
            "PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP",
        )
        output.append(
            wl_case(
                key,
                raw,
                note="raw_density_and_spatial_Euler_signature_emitted",
                reduce_divergence=True,
            )
        )
        key = key_replace(base_key, OBJECT="ADJOINTNESS_RELATION")
        output.append(wl_case(key, wl_field(expression, "ADJOINTNESS_RELATION")))
    return output


def _admissibility_py_parts(family: str, record_raw: str) -> list[tuple[str, str, str]]:
    value_raw = py_field(record_raw)
    pairs = py_residual_map(value_raw) if family == "ADMISSIBILITY_RESIDUAL" else py_string_map(value_raw)
    values = dict(pairs)
    body_pairs = (
        py_residual_map(values["BODY_FORCE"])
        if family == "ADMISSIBILITY_RESIDUAL"
        else py_string_map(values["BODY_FORCE"])
    )
    output = [("BODY_FORCE", dof, raw) for dof, raw in body_pairs]
    faces = values["PER_FACE_TRACTION"]
    if family == "ADMISSIBILITY_RESIDUAL":
        output.append(("PER_FACE_TRACTION", "ALL_FACES", faces))
    else:
        for entry in py_tuple_args(faces):
            face_raw, raw = py_pair(entry)
            face_tokens = py_key_tokens(face_raw)
            if face_tokens == ("1",):
                dof = "PLUS"
            elif face_tokens == ("-1",):
                dof = "MINUS"
            else:
                raise AxisError(f"untyped PY admissibility face key {face_tokens}")
            output.append(("PER_FACE_TRACTION", dof, raw))
    return output


def extract_admissibility(
    engine: str, family: str, payload: str
) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for key_raw, record_raw in py_top_pairs(payload):
            base_key = decode_py_key(key_raw, ("BRANCH", "DENSITY"))
            for object_name, dof, raw in _admissibility_py_parts(family, record_raw):
                key = key_replace(base_key, OBJECT=object_name, DOF=dof)
                output.append(py_case(key, raw, coefficient=False))
        return output

    root_field = {
        "ADMISSIBILITY_OPERATOR_OPERAND": "BACKGROUND_ORDER_BALANCE",
        "ADMISSIBILITY_SUPPORT_OPERAND": "DECLARED_SUPPORT_BUNDLE",
        "ADMISSIBILITY_RESIDUAL": "RESIDUAL_OPERATOR_MINUS_SUPPORT",
    }[family]
    for key_raw, record_raw in wl_assoc_pairs(payload):
        base_key = decode_wl_key(key_raw, ("BRANCH", "DENSITY"))
        bundle = wl_field(record_raw, root_field)
        body = wl_field(bundle, "BULK_DOF_BODY_FORCE")
        for dof, raw in wl_assoc_pairs(body):
            key = key_replace(base_key, OBJECT="BODY_FORCE", DOF=dof)
            output.append(wl_case(key, raw))
        faces = wl_field(bundle, "PER_FACE_TRACTION")
        if family == "ADMISSIBILITY_RESIDUAL":
            key = key_replace(base_key, OBJECT="PER_FACE_TRACTION", DOF="ALL_FACES")
            output.append(wl_case(key, faces))
        else:
            for face, raw in wl_assoc_pairs(faces):
                dof = {"UPPER": "PLUS", "LOWER": "MINUS"}.get(face)
                if dof is None:
                    raise AxisError(f"untyped WL admissibility face key {face}")
                key = key_replace(base_key, OBJECT="PER_FACE_TRACTION", DOF=dof)
                output.append(wl_case(key, raw))
    return output


def extract_control(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    is_residual = family.endswith("_RESIDUAL")
    if engine == "PY":
        schema = (
            ("OBJECT", "BRANCH", "DENSITY", "SOURCE", "DIRECTION")
            if family.startswith("CONTROL_FORM_")
            else ("OBJECT", "BRANCH", "DENSITY")
        )
        for key_raw, record_raw in py_top_pairs(payload):
            key = decode_py_key(key_raw, schema)
            raw = py_field(record_raw)
            coefficient = key_dict(key).get("OBJECT") != "ENERGY_BASIS"
            note = (
                "non_unique_energy_quotient_compared_as_emitted_no_representative_fold"
                if not coefficient
                else None
            )
            if key_dict(key).get("OBJECT") == "TRANSVERSE_DISPERSION":
                output.append(py_case(key, raw, coefficient=coefficient, note=note))
            else:
                output.append(raw_control_case("PY", key, raw, note))
        return output

    if family.startswith("CONTROL_FORM_"):
        for key_raw, record_raw in wl_assoc_pairs(payload):
            key = decode_wl_key(
                key_raw, ("OBJECT", "BRANCH", "DENSITY", "SOURCE", "DIRECTION")
            )
            note = (
                "non_unique_energy_quotient_compared_as_emitted_no_representative_fold"
                if key_dict(key).get("OBJECT") == "ENERGY_BASIS"
                else None
            )
            raw = wl_field(record_raw, "EXPRESSION")
            if key_dict(key).get("OBJECT") == "TRANSVERSE_DISPERSION":
                output.append(wl_case(key, raw, note=note))
            else:
                output.append(raw_control_case("WL", key, raw, note))
        return output

    for key_raw, record_raw in wl_assoc_pairs(payload):
        base_key = decode_wl_key(key_raw, ("BRANCH", "DENSITY"))
        object_map = wl_field(record_raw, "RESIDUAL_A_MINUS_B") if is_residual else record_raw
        for object_name, raw in wl_assoc_pairs(object_map):
            key = key_replace(base_key, OBJECT=object_name)
            if not is_residual:
                raw = wl_field(raw, "EXPRESSION")
            if key_dict(key).get("OBJECT") == "TRANSVERSE_DISPERSION":
                output.append(wl_case(key, raw))
            else:
                output.append(raw_control_case("WL", key, raw, None))
    return output


def extract_dimensions(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for object_raw, raw in py_top_pairs(payload):
            key = make_key((("OBJECT", py_text(object_raw)),))
            output.append(py_case(key, raw, coefficient=False))
        return output
    named = wl_field(payload, "NAMED_OBJECTS")
    for object_name, raw in wl_assoc_pairs(named):
        key = make_key((("OBJECT", object_name),))
        output.append(wl_case(key, raw))
    return output


def extract_homogeneity(engine: str, family: str, payload: str) -> list[ParsedCase]:
    # Both streams carry integer [L,T,M] trees, but their internal labels are
    # deliberately bespoke.  The family-level typed object keeps the complete
    # trees visible and lets the inherited recursive residual operate
    # component-wise wherever their containers align.
    key = make_key((("OBJECT", "HOMOGENEITY_PACKAGE"),))
    if engine == "PY":
        return [py_case(key, payload, coefficient=False)]
    raw = wl_field(payload, "RESIDUAL_A_MINUS_B") if family == "HOMOGENEITY_RESIDUAL" else payload
    return [wl_case(key, raw)]


def extract_generic(engine: str, payload: str) -> list[ParsedCase]:
    key = make_key((("OBJECT", "UNENUMERATED_FAMILY_PAYLOAD"),))
    return [
        py_case(key, payload, coefficient=False)
        if engine == "PY"
        else wl_case(key, payload)
    ]


def extract_engine(engine: str, family: str, payload: str) -> list[ParsedCase]:
    if family in ENERGY_FAMILIES:
        return extract_energy(engine, family, payload)
    if family == "SLAB_OPERATOR":
        return extract_slab(engine, payload)
    if family in {"SLAB_OPERATOR_TERM_ORIGINS", "COUPLING_KERNEL_TERM_ORIGINS"}:
        return extract_term_origins(engine, payload)
    if family == "MU_THETA_OPERATOR":
        return extract_mu(engine, payload)
    if family == "COUPLING_KERNEL":
        return extract_coupling(engine, payload)
    if family in ADMISSIBILITY_FAMILIES:
        return extract_admissibility(engine, family, payload)
    if family.startswith(NESTED_CONTROL_PREFIXES) or family.startswith("CONTROL_FORM_"):
        return extract_control(engine, family, payload)
    if family == "DIMENSIONS":
        return extract_dimensions(engine, payload)
    if family.startswith("HOMOGENEITY_"):
        return extract_homogeneity(engine, family, payload)
    return extract_generic(engine, payload)


def extract_family(
    family: str, py_payload: str | None, wl_payload: str | None
) -> FamilyCases:
    cases = FamilyCases()
    for engine, payload in (("PY", py_payload), ("WL", wl_payload)):
        target = cases.py if engine == "PY" else cases.wl
        if payload is None:
            cases.extraction_notes.append(f"{engine} tag missing")
            continue
        try:
            target.extend(extract_engine(engine, family, payload))
            if not target:
                error = InputError("extractor produced zero cases")
                target.append(extraction_error(engine, family, payload, error))
                cases.extraction_notes.append(f"{engine} extractor produced zero cases")
        except Exception as error:
            target.append(extraction_error(engine, family, payload, error))
            cases.extraction_notes.append(
                f"{engine} extraction failed: {type(error).__name__}: {error}"
            )
    return cases


# ---------- unconditional emission and accounting ----------

def serialise(value: object) -> str:
    class_name = type(value).__name__
    if class_name == "ResidualFailure" and hasattr(value, "reason"):
        return f"ResidualFailure(reason={getattr(value, 'reason')!r})"
    if class_name == "BooleanNotResidualable":
        return "BooleanNotResidualable()"
    if class_name == "Mismatch" and hasattr(value, "kind"):
        return (
            f"Mismatch(kind={getattr(value, 'kind')!r}, "
            f"detail={getattr(value, 'detail', None)!r})"
        )
    if class_name == "UndecidedResidual" and hasattr(value, "reason"):
        return f"UndecidedResidual(reason={getattr(value, 'reason')!r})"
    if class_name == "ResidualAssociation" and hasattr(value, "entries"):
        body = ", ".join(
            f"{key!r}: {serialise(item)}" for key, item in getattr(value, "entries")
        )
        return "ResidualAssociation({" + body + "})"
    if isinstance(value, sp.Basic):
        return sp.srepr(value)
    if isinstance(value, RawCAS):
        return value.raw
    if isinstance(value, SymbolicDifference):
        return (
            f"SymbolicDifference(reason={value.reason!r}, "
            f"expression={sp.srepr(value.expression)})"
        )
    if isinstance(value, TextAtom):
        return f"TextAtom({value.value!r})"
    if isinstance(value, Association):
        body = ", ".join(f"{key!r}: {serialise(item)}" for key, item in value.entries)
        return "Association({" + body + "})"
    if isinstance(value, tuple):
        body = ", ".join(serialise(item) for item in value)
        return "(" + body + ("," if len(value) == 1 else "") + ")"
    return repr(value)


def serialise_key(key: Key) -> str:
    return "(" + ", ".join(f"{axis}={value}" for axis, value in key) + ")"


def display_operand(case: ParsedCase | None) -> str:
    if case is None:
        return "<MISSING>"
    if case.error is not None:
        raw = case.raw
        preview = raw if len(raw) <= 320 else raw[:320] + f"...<nchar={len(raw)}>"
        return f"<PARSE_FAILED {case.error}; RAW={preview}>"
    return serialise(case.value)


def emit_case(
    family: str,
    key: Key,
    py_case_value: ParsedCase | None,
    wl_case_value: ParsedCase | None,
    difference: object,
    *,
    note: str | None = None,
) -> None:
    print(f"CASE family={family} key={serialise_key(key)}", flush=True)
    print(f"operand_A = {display_operand(py_case_value)}", flush=True)
    print(f"operand_B = {display_operand(wl_case_value)}", flush=True)
    print(f"A_minus_B = {serialise(difference)}", flush=True)
    notes = [item for item in (note, py_case_value.note if py_case_value else None, wl_case_value.note if wl_case_value else None) if item]
    if notes:
        print(f"case_note = {' | '.join(dict.fromkeys(notes))}", flush=True)
    for label, value in py_case_value.context if py_case_value is not None else ():
        print(f"context.PY_{label} = {serialise(value)}", flush=True)
    for label, value in wl_case_value.context if wl_case_value is not None else ():
        print(f"context.WL_{label} = {serialise(value)}", flush=True)


def exact_rational_residual(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.together(sp.expand(left - right)))


def replace_budget_failures(value: object, left: object, right: object) -> object:
    if (
        type(value).__name__ == "ResidualFailure"
        and getattr(value, "reason", None) == "RESIDUAL_BUDGET_EXCEEDED"
        and isinstance(left, sp.Expr)
        and isinstance(right, sp.Expr)
    ):
        return exact_rational_residual(left, right)
    if (
        type(value).__name__ == "ResidualAssociation"
        and isinstance(left, Association)
        and isinstance(right, Association)
    ):
        left_items = left.as_dict()
        right_items = right.as_dict()
        return type(value)(
            tuple(
                (key, replace_budget_failures(item, left_items[key], right_items[key]))
                for key, item in value.entries
            )
        )
    if isinstance(value, tuple):
        if isinstance(left, Relational) and isinstance(right, Relational):
            operands = ((left.lhs, right.lhs), (left.rhs, right.rhs))
        elif isinstance(left, tuple) and isinstance(right, tuple):
            operands = tuple(zip(left, right))
        else:
            return value
        if len(value) != len(operands):
            return value
        return tuple(
            replace_budget_failures(item, a, b)
            for item, (a, b) in zip(value, operands)
        )
    return value


def typed_residual(
    left: object,
    right: object,
    family: str,
    budget: float,
    raw_chars: int,
) -> object:
    if (
        isinstance(left, sp.Expr)
        and isinstance(right, sp.Expr)
        and raw_chars > MAX_RATIONAL_CANON_RAW_CHARS
    ):
        return SymbolicDifference(
            sp.Add(left, -right, evaluate=False),
            "exact_rational_canon_deferred_for_oversized_leaf",
        )
    value = residual(left, right, name=family, leaf_budget_seconds=budget)
    if raw_chars > MAX_RATIONAL_CANON_RAW_CHARS:
        # Container leaves retain the typed budget marker.  Both complete raw
        # operands and (for adjointness) both Euler signatures are emitted, so
        # this cannot turn an unresolved large leaf into apparent agreement.
        return value
    return replace_budget_failures(value, left, right)


def axis_mismatch_detail(py_key: Key, wl_key: Key) -> str | None:
    py_axes = key_dict(py_key)
    wl_axes = key_dict(wl_key)
    common = set(py_axes) & set(wl_axes)
    if not common and py_axes and wl_axes:
        return None
    if any(py_axes[axis] != wl_axes[axis] for axis in common):
        return None
    py_only_axes = sorted(set(py_axes) - set(wl_axes))
    wl_only_axes = sorted(set(wl_axes) - set(py_axes))
    if not py_only_axes and not wl_only_axes:
        return None
    details: list[str] = []
    if py_only_axes:
        details.append("WL missing " + ",".join(py_only_axes))
    if wl_only_axes:
        details.append("PY missing " + ",".join(wl_only_axes))
    return "; ".join(details)


def compare_family(family: str, cases: FamilyCases, *, budget: float) -> Accounting:
    accounting = Accounting()
    py_by_key: dict[Key, list[ParsedCase]] = defaultdict(list)
    wl_by_key: dict[Key, list[ParsedCase]] = defaultdict(list)
    for case in cases.py:
        py_by_key[case.key].append(case)
        accounting.parse_failed += int(case.error is not None)
    for case in cases.wl:
        wl_by_key[case.key].append(case)
        accounting.parse_failed += int(case.error is not None)

    duplicate_keys = {
        key for key, values in py_by_key.items() if len(values) > 1
    } | {key for key, values in wl_by_key.items() if len(values) > 1}
    accounting.duplicate_key = sum(
        len(py_by_key.get(key, ())) + len(wl_by_key.get(key, ()))
        for key in duplicate_keys
    )

    common_keys = (set(py_by_key) & set(wl_by_key)) - duplicate_keys
    for key in sorted(common_keys, key=serialise_key):
        py_item = py_by_key[key][0]
        wl_item = wl_by_key[key][0]
        accounting.parse_failed += int(materialize(py_item))
        accounting.parse_failed += int(materialize(wl_item))
        accounting.join += 1
        if py_item.error is not None or wl_item.error is not None:
            difference: object = TextAtom("UNDEFINED_PARSE_FAILED")
        else:
            left = py_item.compare_value
            right = wl_item.compare_value
            difference = typed_residual(
                left,
                right,
                family,
                budget,
                len(py_item.raw) + len(wl_item.raw),
            )
        emit_case(family, key, py_item, wl_item, difference)
        release_case(py_item)
        release_case(wl_item)

    for key in sorted(duplicate_keys, key=serialise_key):
        py_values = py_by_key.get(key, ())
        wl_values = wl_by_key.get(key, ())
        width = max(len(py_values), len(wl_values))
        for index in range(width):
            py_item = py_values[index] if index < len(py_values) else None
            wl_item = wl_values[index] if index < len(wl_values) else None
            if py_item is not None:
                accounting.parse_failed += int(materialize(py_item))
            if wl_item is not None:
                accounting.parse_failed += int(materialize(wl_item))
            emit_case(
                family,
                key,
                py_item,
                wl_item,
                TextAtom("UNDEFINED_DUPLICATE_KEY"),
                note=f"duplicate_key occurrence={index + 1}/{width}",
            )
            release_case(py_item)
            release_case(wl_item)

    py_unmatched = set(py_by_key) - set(wl_by_key) - duplicate_keys
    wl_unmatched = set(wl_by_key) - set(py_by_key) - duplicate_keys
    accounting.py_only = sum(len(py_by_key[key]) for key in py_unmatched)
    accounting.wl_only = sum(len(wl_by_key[key]) for key in wl_unmatched)

    mismatch_py: set[Key] = set()
    mismatch_wl: set[Key] = set()
    reasons_by_key: dict[tuple[str, Key], set[str]] = defaultdict(set)
    for py_key in py_unmatched:
        for wl_key in wl_unmatched:
            detail = axis_mismatch_detail(py_key, wl_key)
            if detail is None:
                continue
            mismatch_py.add(py_key)
            mismatch_wl.add(wl_key)
            reasons_by_key[("PY", py_key)].add(detail)
            reasons_by_key[("WL", wl_key)].add(detail)
            accounting.reasons.add(detail)
    accounting.axis_set_mismatch = sum(len(py_by_key[key]) for key in mismatch_py) + sum(
        len(wl_by_key[key]) for key in mismatch_wl
    )

    for engine, keys, cases_by_key in (
        ("PY", py_unmatched, py_by_key),
        ("WL", wl_unmatched, wl_by_key),
    ):
        for key in sorted(keys, key=serialise_key):
            reasons = sorted(reasons_by_key.get((engine, key), ()))
            note = "axis_set_mismatch: " + " | ".join(reasons) if reasons else f"{engine.lower()}_only"
            for case in cases_by_key[key]:
                accounting.parse_failed += int(materialize(case))
                emit_case(
                    family,
                    key,
                    case if engine == "PY" else None,
                    case if engine == "WL" else None,
                    TextAtom("UNDEFINED_UNJOINED"),
                    note=note,
                )
                release_case(case)

    accounting.reasons.update(cases.extraction_notes)
    if family in {
        "ENERGY_BASIS_VARIABLE",
        "ENERGY_BASIS_COUNT",
        "ENERGY_BASIS_NEW_INVARIANTS",
    }:
        accounting.reasons.add("energy representative quotient compared as emitted")
    reason_text = " | ".join(sorted(accounting.reasons)) if accounting.reasons else "none"
    print(
        f"ACCOUNTING {family} "
        "{"
        f"join={accounting.join}, py_only={accounting.py_only}, wl_only={accounting.wl_only}, "
        f"duplicate_key={accounting.duplicate_key}, parse_failed={accounting.parse_failed}, "
        f"axis_set_mismatch={accounting.axis_set_mismatch}"
        f"}} reasons={reason_text}",
        flush=True,
    )
    return accounting


def emit_local_inventory(py_tags: dict[str, str], wl_tags: dict[str, str]) -> None:
    py_local = sorted(name for name in py_tags if name.startswith("LOCAL_"))
    wl_local = sorted(name for name in wl_tags if name.startswith("LOCAL_"))
    print(f"LOCAL_INVENTORY engine=PY excluded={py_local}", flush=True)
    print(f"LOCAL_INVENTORY engine=WL excluded={wl_local}", flush=True)
    if "LOCAL_TAG_NAMES" in wl_tags:
        print(f"LOCAL_INVENTORY engine=WL payload={wl_tags['LOCAL_TAG_NAMES']}", flush=True)


def run(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--py", type=Path, default=DEFAULT_PY)
    parser.add_argument("--wl", type=Path, default=DEFAULT_WL)
    parser.add_argument(
        "--family",
        action="append",
        help="optional shared family name to measure (repeatable); default is every emitted family",
    )
    parser.add_argument(
        "--residual-leaf-budget",
        type=float,
        default=0.1,
        help="seconds for the optional leaf factorization before exact expand/cancel fallback",
    )
    arguments = parser.parse_args(argv)
    started = time.monotonic()
    try:
        py_tags = load_py(arguments.py)
        wl_tags = load_wl(arguments.wl)
    except Exception as error:
        print(f"OPERATIONAL_ERROR {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        return 2

    emit_local_inventory(py_tags, wl_tags)
    families = sorted(
        {
            name
            for name in (*py_tags, *wl_tags)
            if not name.startswith("LOCAL_")
        }
    )
    if arguments.family:
        families = sorted(set(arguments.family))
    results: dict[str, Accounting] = {}
    for family in families:
        cases = extract_family(family, py_tags.get(family), wl_tags.get(family))
        results[family] = compare_family(
            family, cases, budget=arguments.residual_leaf_budget
        )

    elapsed = time.monotonic() - started
    print(
        f"RUN_ACCOUNTING families={len(families)} "
        f"families_with_join={sum(item.join > 0 for item in results.values())} "
        f"families_with_unpaired={sum((item.py_only + item.wl_only) > 0 for item in results.values())} "
        f"parse_failed={sum(item.parse_failed for item in results.values())} "
        f"duplicate_key={sum(item.duplicate_key for item in results.values())} "
        f"runtime_seconds={elapsed:.3f}",
        flush=True,
    )
    print(
        "MEASUREMENT_SCOPE supplied_unfalsifiable=sections_1_to_3,"
        "admissibility_support_premise,supplied_bookkeeping residual_target=none",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
