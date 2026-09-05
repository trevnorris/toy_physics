#!/usr/bin/env python3
"""S11c-c1 T7 cross-engine operand comparator.

This is a measurement instrument.  It reads the two committed c1 tag streams,
uses axis-typed keys, applies only the closed spelling/CAS-form map, and emits
both operands and the typed recursive A-minus-B residual.  Algebraic,
structural, and coverage differences are output data and never affect the exit
code.  Neither measured engine is imported or run.
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

import S11c_a_cross_engine_comparator as s11ca  # noqa: E402
from S11b_cross_engine_comparator import (  # noqa: E402
    Association,
    BooleanNotResidualable,
    Mismatch,
    ResidualAssociation,
    ResidualFailure,
    TextAtom,
    UndecidedResidual,
    _convert_parsed_containers,
    mechanical_lower_camel,
    parse_sympy_payload,
    parse_wolfram_payload,
    residual,
)


DEFAULT_PY = SCRIPT_DIR / "out" / "S11c_c1_bulk_closure_sympy_audit.out"
DEFAULT_WL = (
    SCRIPT_DIR.parent
    / "mathematica"
    / "out"
    / "S11c_c1_bulk_closure_mathematica_audit.out"
)

PY_TAG_LINE = re.compile(r"^(PY_S11CC1_[A-Z0-9_]+):\s?(.*)$")
WL_TAG_LINE = re.compile(r"^(WL_S11CC1_[A-Z0-9_]+):")

# LEAF and CASE are structural coordinates, not physics-axis catch-alls.  Every
# physics-bearing token named by the c1 constructor has its own axis below.
AXIS_ORDER = (
    "OBJECT",
    "MUTATION",
    "BRANCH",
    "FACE",
    "DENSITY",
    "DIRECTION",
    "PARITY",
    "REGIME_OUT",
    "REGIME_IN",
    "SCENARIO",
    "LIMIT",
    "PORT",
    "CASE",
    "LEAF",
)
Key = tuple[tuple[str, str], ...]

BRANCHES = {"LAB_HELD", "MATERIAL_ADVECTED"}
FACES = {"1", "-1"}
DENSITIES = {"RHO4_CONSTANT", "RHOBR_CONSTANT"}
DIRECTIONS = {"1", "2", "3"}
PARITIES = {
    "THICKNESS",
    "CENTRE_SHIFT",
    "DELTA_W",
    "ZETA_C",
    "COMBINATION_COMPARISON",
    "PLUS_FACE_CONSTRUCTION_CONTROL",
}
PY_REGIMES = {"PROPAGATING", "EVANESCENT", "GRAZING"}
WL_REGIMES = {
    "OUTPUT_PROPAGATING",
    "OUTPUT_EVANESCENT",
    "OUTPUT_GRAZING",
    "INPUT_PROPAGATING",
    "INPUT_EVANESCENT",
    "INPUT_GRAZING",
}
REGIME_OUT_VALUES = PY_REGIMES | {value for value in WL_REGIMES if value.startswith("OUTPUT_")}
REGIME_IN_VALUES = PY_REGIMES | {value for value in WL_REGIMES if value.startswith("INPUT_")}
MUTATIONS = {
    "BASE",
    "SIGNFLIP_INPUT",
    "SIGNFLIP_OUTPUT",
    "MOMENTUMFREEZE_INPUT",
    "MOMENTUMFREEZE_OUTPUT",
}
PORTS = {"A", "V", "X"}
SCENARIOS = {
    "REAL_OMEGA_PROPAGATING_IMPERMEABLE_LAMBDA_X0_ZERO",
    "EVANESCENT_PURELY_REACTIVE",
    "PURELY_REACTIVE_BLOCK",
}
LIMITS = {"OMEGA_TAU_TO_ZERO", "OMEGA_TAU_TO_INFINITY"}
MAX_RATIONAL_CANON_RAW_CHARS = 150_000


class InputError(ValueError):
    """A committed stream or requested extraction violates its grammar."""


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
    loader: Callable[[], object] | None = field(default=None, repr=False)


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
    extracted_leaves_py: int = 0
    extracted_leaves_wl: int = 0
    reasons: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class RawCAS:
    engine: str
    raw: str


@dataclass(frozen=True)
class SymbolicDifference:
    expression: sp.Expr
    reason: str


def load_py(path: Path) -> dict[str, str]:
    """Load the mandatory single-line PY srepr stream."""
    if not path.is_file():
        raise InputError(f"PY input does not exist: {path}")
    output: dict[str, str] = {}
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            match = PY_TAG_LINE.fullmatch(line.rstrip("\r\n"))
            if match is None:
                raise InputError(f"{path}:{line_number}: not one PY S11c-c1 tagged row")
            full_tag, payload = match.groups()
            name = full_tag.removeprefix("PY_S11CC1_")
            if name in output:
                raise InputError(f"{path}:{line_number}: duplicate PY tag {name}")
            output[name] = payload
    if not output:
        raise InputError(f"{path}: no PY S11c-c1 tags")
    return output


def load_wl(path: Path) -> dict[str, str]:
    """Load WL tags, reassembling each payload through the next tag line."""
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
                current_name = match.group(1).removeprefix("WL_S11CC1_")
                chunks = [clean.split(":", 1)[1]]
            elif current_name is None:
                if clean.strip():
                    raise InputError(f"{path}:{line_number}: content before first WL tag")
            else:
                chunks.append(clean)
    finish()
    if not output:
        raise InputError(f"{path}: no WL S11c-c1 tags")
    return output


# Verified delimiter-aware structural readers.
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
    """Build a fixed-order key, rejecting duplicate and untyped axes."""
    decoded: dict[str, str] = {}
    for axis, original in items:
        value = str(original)
        if axis in decoded:
            raise AxisError(f"two values for axis {axis}: {decoded[axis]} and {value}")
        if axis not in AXIS_ORDER:
            raise AxisError(f"unknown axis {axis}")
        if axis == "OBJECT" and not value:
            raise AxisError("empty OBJECT value")
        if axis == "BRANCH" and value not in BRANCHES:
            raise AxisError(f"invalid BRANCH value {value}")
        if axis == "FACE" and value not in FACES:
            raise AxisError(f"invalid FACE value {value}; expected 1 or -1")
        if axis == "DENSITY" and value not in DENSITIES:
            raise AxisError(f"invalid DENSITY value {value}")
        if axis == "DIRECTION" and value not in DIRECTIONS:
            raise AxisError(f"invalid DIRECTION value {value}; expected 1, 2, or 3")
        if axis == "PARITY" and value not in PARITIES:
            raise AxisError(f"invalid PARITY value {value}")
        if axis == "REGIME_OUT" and value not in REGIME_OUT_VALUES:
            raise AxisError(f"invalid {axis} value {value}")
        if axis == "REGIME_IN" and value not in REGIME_IN_VALUES:
            raise AxisError(f"invalid {axis} value {value}")
        if axis == "MUTATION" and value not in MUTATIONS:
            raise AxisError(f"invalid MUTATION value {value}")
        if axis == "PORT" and value not in PORTS:
            raise AxisError(f"invalid PORT value {value}")
        if axis == "SCENARIO" and value not in SCENARIOS:
            raise AxisError(f"invalid SCENARIO value {value}")
        if axis == "LIMIT" and value not in LIMITS:
            raise AxisError(f"invalid LIMIT value {value}")
        if axis in {"CASE", "LEAF"} and not value:
            raise AxisError(f"empty {axis} value")
        decoded[axis] = value
    return tuple((axis, decoded[axis]) for axis in AXIS_ORDER if axis in decoded)


def key_dict(key: Key) -> dict[str, str]:
    return dict(key)


def key_replace(key: Key, *, remove: Sequence[str] = (), **updates: str) -> Key:
    values = {axis: value for axis, value in key if axis not in set(remove)}
    values.update(updates)
    return make_key(values.items())


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
        prefix = {
            "FACE": "FACE_",
            "DIRECTION": "DIRECTION_",
            "PARITY": "PARITY_",
        }.get(axis)
        if prefix is not None:
            if not token.startswith(prefix):
                raise AxisError(f"untyped WL {axis} token {token}")
            token = token.removeprefix(prefix)
        typed.append((axis, token))
    return make_key(typed)


# ---------- closed spelling/CAS-form map ----------

# Every active c1 addition is bare-symbol to bare-symbol.
C1_BARE_SYMBOL = {
    "cS0": "c_s0",
    "rhoM": "rho_m",
    "W0": "W_0",
    "LW": "L_W",
    "sigmaW": "sigma_W",
    "etaBg": "eta_bg",
    "epsilonShape": "epsilon_shape",
    "LambdaA0": "Lambda_A_0",
    "LambdaV0": "Lambda_V_0",
    "LambdaX0": "Lambda_X_0",
    "tauA": "tau_A",
    "tauV": "tau_V",
}
REQUIRED_C1_PY_SYMBOLS = frozenset(C1_BARE_SYMBOL.values())
ACTIVE_C1_BARE_SYMBOL = dict(C1_BARE_SYMBOL)


def checked_mechanical_symbol_map(reserved_names: Iterable[str]) -> dict[str, str]:
    """Return injective lowerCamel -> PY spellings or emit a collision finding."""
    inverse: dict[str, str] = {}
    collisions: dict[str, set[str]] = defaultdict(set)
    for name in reserved_names:
        target = mechanical_lower_camel(name)
        if target in inverse and inverse[target] != name:
            collisions[target].update((inverse[target], name))
        inverse[target] = name
    if collisions:
        detail = "; ".join(
            f"{target}<-{sorted(names)}" for target, names in sorted(collisions.items())
        )
        raise InputError(f"non-injective mechanical spelling map: {detail}")
    return inverse


def real_py_symbol_names(raw_payloads: Iterable[str]) -> set[str]:
    pattern = re.compile(r"Symbol\('([^']+)'(?:,|\))")
    return {name for raw in raw_payloads for name in pattern.findall(raw)}


def verify_active_spelling_map(py_payloads: Iterable[str]) -> tuple[dict[str, str], set[str]]:
    names = real_py_symbol_names(py_payloads)
    checked_mechanical_symbol_map(names)
    missing = REQUIRED_C1_PY_SYMBOLS - names
    active = {
        wl_name: py_name
        for wl_name, py_name in C1_BARE_SYMBOL.items()
        if py_name in names and mechanical_lower_camel(py_name) == wl_name
    }
    invalid = set(C1_BARE_SYMBOL) - set(active)
    if missing or invalid:
        raise InputError(
            f"c1 spelling fold failed real-symbol/mechanical check: "
            f"missing={sorted(missing)} invalid={sorted(invalid)}"
        )
    return active, names


def preprocess_wl(raw: str) -> str:
    """Use inherited parsing folds without geometry-axis or active-head damage."""
    protected = re.sub(
        r"\bmuThetaOperand(?=\s*\[)", "C1ProtectedmuThetaOperand", raw
    )
    # The inherited sibling map already knows some of these spellings.  Shield
    # bare occurrences so this comparator activates them only after its own
    # real-PY-Symbol and injectivity checks.  Applied heads are intentionally
    # excluded by the lookahead.
    for wl_name in sorted(C1_BARE_SYMBOL, key=len, reverse=True):
        protected = re.sub(
            rf'(?<![A-Za-z0-9_\"]){re.escape(wl_name)}(?![A-Za-z0-9_\"]|\s*\[)',
            f"C1BareProtected{wl_name}",
            protected,
        )
    # Hold unsupported multi-range integrals before the inherited one-range
    # parser sees them.  Keep every argument; nested Limit uses the held-head
    # rewrite below.  One-range Integrate and Equal retain their inherited paths.
    def hold_multi_range(match: re.Match[str]) -> str:
        body, _ = s11ca.read_bracket(protected, match.end() - 1)
        if len(split_top(body)) > 2:
            return "HeldInactiveIntegrate["
        return match.group(0)

    protected = re.sub(r"Inactive\[Integrate\]\[", hold_multi_range, protected)
    text = s11ca.preprocess_wl(protected)
    # Undo S11c-a's geometry vocabulary.  C1 numeric association rules are
    # inner FACE keys and must remain the literal typed values 1/-1.
    text = text.replace('"PLUS" ->', '"1" ->').replace('"MINUS" ->', '"-1" ->')
    # PROFILE is inherited for verified S11c-a fields, but the c1 contract
    # explicitly seals applied w1Jet/m1Jet heads.  Protect their arguments
    # from S11c-a's applied-head-to-symbol conversion.
    text = re.sub(
        r"\b((?:w1|m1)Jet[A-Za-z0-9_]*)\[",
        r"C1Protected\1[",
        text,
    )
    # Mathematica's Python parser cannot read Inactive[head][...].  Preserve
    # every non-Equal/non-Integrate head as an unevaluated applied head.
    text = re.sub(r"Inactive\[([A-Za-z][A-Za-z0-9]*)\]\[", r"HeldInactive\1[", text)
    return text


def _restore_protected_heads(value: sp.Basic) -> sp.Basic:
    value = value.replace(
        lambda node: isinstance(node, AppliedUndef)
        and node.func.__name__.startswith("C1Protected"),
        lambda node: sp.Function(node.func.__name__.removeprefix("C1Protected"))(
            *node.args
        ),
    )
    replacements: dict[sp.Symbol, sp.Symbol] = {}
    for symbol in value.atoms(sp.Symbol):
        if not symbol.name.startswith("C1BareProtected"):
            continue
        wl_name = symbol.name.removeprefix("C1BareProtected")
        final_name = ACTIVE_C1_BARE_SYMBOL.get(wl_name, wl_name)
        replacements[symbol] = sp.Symbol(
            final_name, **getattr(symbol, "_assumptions_orig", {})
        )
    return value.xreplace(replacements) if replacements else value


def _assumption_preserving_symbols(value: sp.Basic) -> sp.Basic:
    replacements: dict[sp.Symbol, sp.Basic] = {}
    for symbol in value.atoms(sp.Symbol):
        if isinstance(symbol, sp.Dummy):
            continue
        name = s11ca.canon_jet_name(s11ca.pynorm(symbol.name))
        if name == "Infinity":
            replacements[symbol] = sp.oo
        elif name == "w":
            replacements[symbol] = s11ca.W_SYMBOL
        elif name != symbol.name:
            replacements[symbol] = sp.Symbol(
                name, **getattr(symbol, "_assumptions_orig", {})
            )
    return value.xreplace(replacements) if replacements else value


def _c1_bare_symbols(value: sp.Basic) -> sp.Basic:
    replacements = {
        symbol: sp.Symbol(ACTIVE_C1_BARE_SYMBOL[symbol.name])
        for symbol in value.atoms(sp.Symbol)
        if symbol.name in ACTIVE_C1_BARE_SYMBOL
    }
    return value.xreplace(replacements) if replacements else value


def canonical_basic(value: sp.Basic, engine: str, branch: str | None) -> sp.Basic:
    if engine == "WL":
        value = s11ca.canon_wl_basic(s11ca.spint_to_integral(value), branch)
        value = _restore_protected_heads(value)
        value = _c1_bare_symbols(value)
    else:
        value = s11ca.py_to_dwin(value)
    value = _assumption_preserving_symbols(value)
    return s11ca.canonical_integrals(value)


def canonical_value(value: object, engine: str, branch: str | None) -> object:
    if isinstance(value, Association):
        return Association(
            tuple((key, canonical_value(item, engine, branch)) for key, item in value.entries)
        )
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(lambda item: canonical_basic(item, engine, branch))
    if isinstance(value, (tuple, list, sp.Tuple)):
        return tuple(canonical_value(item, engine, branch) for item in value)
    if isinstance(value, Relational):
        return type(value)(
            canonical_value(value.lhs, engine, branch),
            canonical_value(value.rhs, engine, branch),
            evaluate=False,
        )
    if isinstance(value, sp.Basic):
        return canonical_basic(value, engine, branch)
    return value


def _parse_py(raw: str, coefficient: bool, branch: str | None) -> object:
    value = parse_sympy_payload(raw)
    if coefficient:
        value = s11ca.coeff_epsilon(value)
    return _convert_parsed_containers(canonical_value(value, "PY", branch))


def _parse_wl(raw: str, branch: str | None) -> object:
    value = parse_wolfram_payload(preprocess_wl(raw))
    return _convert_parsed_containers(canonical_value(value, "WL", branch))


def parse_py_value(raw: str, key: Key, *, coefficient: bool = True) -> object:
    return _parse_py(raw, coefficient, key_dict(key).get("BRANCH"))


def parse_wl_value(raw: str, key: Key) -> object:
    return _parse_wl(raw, key_dict(key).get("BRANCH"))


BOUND_INTEGRAL = s11ca.BOUND_INTEGRAL
BOUND_BINDER = s11ca.BOUND_BINDER
CAPTURED_FREE_BINDER = s11ca.CAPTURED_FREE_BINDER
canonical_integrals = s11ca.canonical_integrals


# ---------- lazy cases ----------

def build_case(
    engine: str,
    key: Key,
    raw: str,
    parser: Callable[[], object],
    *,
    note: str | None = None,
) -> ParsedCase:
    return ParsedCase(key, None, raw, note=note, loader=parser)


def materialize(case: ParsedCase) -> bool:
    if case.loader is None:
        return False
    loader = case.loader
    case.loader = None
    try:
        case.value = loader()
        case.compare_value = case.value
        return False
    except Exception as error:
        case.error = f"{type(error).__name__}: {error}"
        return True


def release_case(case: ParsedCase | None) -> None:
    if case is None:
        return
    case.value = None
    case.compare_value = None


def py_case(key: Key, raw: str, *, coefficient: bool = False, note: str | None = None) -> ParsedCase:
    return build_case(
        "PY", key, raw, lambda: parse_py_value(raw, key, coefficient=coefficient), note=note
    )


def wl_case(key: Key, raw: str, *, note: str | None = None) -> ParsedCase:
    return build_case("WL", key, raw, lambda: parse_wl_value(raw, key), note=note)


RAW_CONTROL_WHITELIST = {
    ("DTN_OPERATOR", "WHOLE_OBJECT"),
    ("NONINVERTIBILITY_CONDITION", "OPERATOR"),
}


def raw_control_case(
    engine: str,
    key: Key,
    raw: str,
    note: str | None = None,
    *,
    family: str,
    leaf: str,
) -> ParsedCase:
    """Create a raw outer-signature case only for the sealed operator algebra."""
    if (family, leaf) not in RAW_CONTROL_WHITELIST:
        raise InputError(
            f"raw_control_case refused for {family}.{leaf}: descend to scalar leaves"
        )
    try:
        if engine == "PY":
            try:
                names = [name for name, _ in _py_named_pairs(raw)]
                signature = "ASSOCIATION_KEYS:" + ",".join(names)
            except Exception:
                try:
                    signature = f"TUPLE_ARITY:{len(py_tuple_args(raw))}"
                except Exception:
                    match = re.match(r"([A-Za-z][A-Za-z0-9_]*)\s*\(", raw.strip())
                    signature = "RAW_HEAD:" + (match.group(1) if match else "SCALAR")
        else:
            try:
                signature = "ASSOCIATION_KEYS:" + ",".join(
                    name for name, _ in wl_assoc_pairs(raw)
                )
            except Exception:
                match = re.match(r"([A-Za-z][A-Za-z0-9_]*)\s*\[", raw.strip())
                signature = "RAW_HEAD:" + (match.group(1) if match else "SCALAR")
        return ParsedCase(
            key,
            RawCAS(engine, raw),
            raw,
            compare_value=TextAtom(f"{engine}:{signature}"),
            note="outer_operator_algebra_signature; raw_operand_emitted"
            + (f" | {note}" if note else ""),
        )
    except Exception as error:
        return ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}", note=note)


def extraction_error(engine: str, family: str, raw: str, error: Exception) -> ParsedCase:
    key = make_key((("OBJECT", family), ("LEAF", f"{engine}_EXTRACTION_ERROR")))
    return ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}")


# ---------- raw structural extraction ----------

def _py_named_pairs(raw: str) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for entry in py_tuple_args(raw):
        key_raw, value_raw = py_pair(entry)
        pairs.append((py_text(key_raw), value_raw))
    return pairs


def _py_named_dict(raw: str) -> dict[str, str]:
    return dict(_py_named_pairs(raw))


def _wl_dict(raw: str) -> dict[str, str]:
    return dict(wl_assoc_pairs(raw))


def _maybe_named_pairs(engine: str, raw: str) -> list[tuple[str, str]] | None:
    try:
        pairs = _py_named_pairs(raw) if engine == "PY" else wl_assoc_pairs(raw)
    except Exception:
        return None
    return pairs


def _append_leaf(
    output: list[ParsedCase],
    engine: str,
    key: Key,
    leaf: str,
    raw: str,
    *,
    note: str | None = None,
) -> None:
    leaf_key = key_replace(key, LEAF=leaf)
    output.append(
        py_case(leaf_key, raw, coefficient=False, note=note)
        if engine == "PY"
        else wl_case(leaf_key, raw, note=note)
    )


FIELD_ALIAS = {
    "VALUE": "EXPRESSION",
    "MULTIGRADE": "MULTIGRADE_EPSILON_ETA_SIGMAW",
    "DELTA_P": "PRESSURE",
    "J": "RELATIVE_FLUX",
    "T": "TRACTION",
    "TRACTION_NORMAL_AMPLITUDE": "TRACTION",
}


def _flatten_named(
    output: list[ParsedCase],
    engine: str,
    key: Key,
    raw: str,
    prefix: str,
    *,
    aliases: dict[str, str] | None = None,
    note: str | None = None,
) -> None:
    """Descend raw named containers so mismatched shells cannot hide siblings."""
    pairs = _maybe_named_pairs(engine, raw)
    if pairs:
        rename = FIELD_ALIAS if aliases is None else {**FIELD_ALIAS, **aliases}
        for name, child in pairs:
            canonical = rename.get(name, name)
            path = f"{prefix}.{canonical}" if prefix else canonical
            _flatten_named(output, engine, key, child, path, aliases=aliases, note=note)
        return
    _append_leaf(output, engine, key, prefix or "EXPRESSION", raw, note=note)


def _py_outer(payload: str) -> list[tuple[str, str]]:
    return py_top_pairs(payload)


def _wl_outer(payload: str) -> list[tuple[str, str]]:
    return wl_assoc_pairs(payload)


def _py_face_pairs(raw: str) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    for key_raw, value_raw in py_top_pairs(raw):
        tokens = py_key_tokens(key_raw)
        if len(tokens) != 1:
            raise AxisError(f"PY inner FACE key is not one token: {key_raw[:80]}")
        output.append((make_key((("FACE", tokens[0]),))[0][1], value_raw))
    return output


def _wl_face_pairs(raw: str) -> list[tuple[str, str]]:
    return [
        (key_dict(decode_wl_key(name, ("FACE",)))["FACE"], value)
        for name, value in wl_assoc_pairs(raw)
    ]


def _extract_case_map_primary(
    engine: str,
    payload: str,
    py_schema: Sequence[str],
    wl_schema: Sequence[str],
    *,
    py_path: Sequence[str] = ("VALUE",),
    wl_path: Sequence[str] = ("EXPRESSION",),
    leaf: str = "EXPRESSION",
    extra: Iterable[tuple[str, str]] = (),
) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    for raw_key, record in outer:
        key = (
            decode_py_key(raw_key, py_schema)
            if engine == "PY"
            else decode_wl_key(raw_key, wl_schema)
        )
        if extra:
            key = make_key((*key, *extra))
        try:
            item = record
            for field_name in py_path if engine == "PY" else wl_path:
                item = (
                    _py_named_dict(item)[field_name]
                    if engine == "PY"
                    else _wl_dict(item)[field_name]
                )
            _append_leaf(output, engine, key, leaf, item)
        except Exception as error:
            output.append(extraction_error(engine, leaf, record, error))
    return output


def _mutation_for_family(family: str) -> str:
    if family == "CONTROL_BRANCH_BASE_OPERAND":
        return "BASE"
    match = re.fullmatch(r"CONTROL_BRANCH_(.+)_OPERAND", family)
    if match is None or match.group(1) not in MUTATIONS:
        raise InputError(f"cannot derive MUTATION from {family}")
    return match.group(1)


# ---------- bespoke family extractors ----------

def extract_dtn_flat(engine: str, payload: str) -> list[ParsedCase]:
    if engine == "PY":
        return _extract_case_map_primary(
            engine, payload, ("BRANCH", "FACE"), (), leaf="FLAT_SYMBOL"
        )
    output: list[ParsedCase] = []
    root = make_key(())
    fields = _wl_dict(payload)
    derivation = _wl_dict(fields["C1_INDEPENDENT_DERIVATION"])
    _append_leaf(output, engine, root, "FLAT_SYMBOL", derivation["EXPRESSION"])
    for name, raw in fields.items():
        if name == "C1_INDEPENDENT_DERIVATION":
            continue
        _flatten_named(output, engine, root, raw, name)
    return output


def extract_dtn_operator(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for raw_key, record in _py_outer(payload):
            key = decode_py_key(raw_key, ("BRANCH", "FACE"))
            raw = _py_named_dict(record)["VALUE"]
            output.append(
                raw_control_case(
                    engine, key_replace(key, LEAF="OPERATOR"), raw,
                    family="DTN_OPERATOR", leaf="WHOLE_OBJECT"
                )
            )
    else:
        for branch, record in _wl_outer(payload):
            key = make_key((("BRANCH", branch), ("LEAF", "OPERATOR")))
            raw = _wl_dict(_wl_dict(record)["COMPOSITION"])["EXPRESSION"]
            output.append(
                raw_control_case(
                    engine, key, raw, family="DTN_OPERATOR", leaf="WHOLE_OBJECT"
                )
            )
    return output


def extract_dtn_kernel(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for raw_key, record in _py_outer(payload):
            key = decode_py_key(raw_key, ("BRANCH", "FACE"))
            values = _py_named_dict(_py_named_dict(record)["VALUE"])
            _append_leaf(output, engine, key, "KERNEL_EXPRESSION", values["FIRST_SHAPE"])
            for name in ("FLAT_DIAGONAL", "ASSEMBLED"):
                _append_leaf(output, engine, key, name, values[name])
    else:
        for raw_key, record in _wl_outer(payload):
            key = decode_wl_key(raw_key, ("BRANCH", "FACE"))
            fields = _wl_dict(record)
            kernel = _wl_dict(fields["KERNEL"])
            _append_leaf(output, engine, key, "KERNEL_EXPRESSION", kernel["EXPRESSION"])
            for name in ("OUTPUT_LEG", "INPUT_LEG", "MOMENTUM_TRANSFER"):
                _append_leaf(output, engine, key, name, fields[name])
    return output


def extract_rigid(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if family.endswith("RESIDUAL"):
        if engine == "PY":
            return _extract_case_map_primary(
                engine, payload, ("BRANCH", "FACE"), (), leaf="RESIDUAL"
            )
        for branch, record in _wl_outer(payload):
            key = make_key((("BRANCH", branch),))
            _append_leaf(output, engine, key, "RESIDUAL", _wl_dict(record)["EXPRESSION"])
        return output
    if engine == "PY":
        for raw_key, record in _py_outer(payload):
            key = decode_py_key(raw_key, ("BRANCH", "FACE"))
            values = py_tuple_args(_py_named_dict(record)["VALUE"])
            names = ("CURVED_DIAGONAL_OPERAND_A", "FLAT_TRANSLATED_OPERAND_B")
            for name, raw in zip(names, values):
                _append_leaf(output, engine, key, name, raw)
    else:
        for branch, record in _wl_outer(payload):
            key = make_key((("BRANCH", branch),))
            fields = _wl_dict(record)
            for name in ("CURVED_DIAGONAL_OPERAND_A", "FLAT_TRANSLATED_OPERAND_B"):
                _append_leaf(output, engine, key, name, _wl_dict(fields[name])["EXPRESSION"])
    return output


def extract_regime_pair(engine: str, payload: str) -> list[ParsedCase]:
    if engine == "PY":
        output: list[ParsedCase] = []
        schema = ("BRANCH", "FACE", "REGIME_OUT", "REGIME_IN")
        for raw_key, record in _py_outer(payload):
            key = decode_py_key(raw_key, schema)
            fields = _py_named_dict(record)
            _append_leaf(output, engine, key, "KERNEL_EXPRESSION", fields["VALUE"])
            if "CONDITIONS" in fields:
                _append_leaf(output, engine, key, "CONDITIONS", fields["CONDITIONS"])
        return output
    output: list[ParsedCase] = []
    schema = ("BRANCH", "FACE", "REGIME_OUT", "REGIME_IN")
    for raw_key, record in _wl_outer(payload):
        key = decode_wl_key(raw_key, schema)
        fields = _wl_dict(record)
        expression = fields.get("EXPRESSION")
        for container_name in ("KERNEL_OBJECT", "KERNEL"):
            if expression is None and container_name in fields:
                expression = _wl_dict(fields[container_name])["EXPRESSION"]
        if expression is not None:
            _append_leaf(output, engine, key, "KERNEL_EXPRESSION", expression)
        for name in ("OUTPUT_CONDITION", "INPUT_CONDITION", "OUTPUT_BRANCH_VALUE", "INPUT_BRANCH_VALUE"):
            if name in fields:
                _append_leaf(output, engine, key, name, fields[name])
    return output


def extract_parity(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    for raw_key, record in outer:
        branch = py_text(raw_key) if engine == "PY" else raw_key
        key = make_key((("BRANCH", branch),))
        if engine == "PY":
            value = _py_named_dict(record)["VALUE"]
            fields: dict[str, str] = {}
            for entry in py_tuple_args(value):
                parts = py_tuple_args(entry)
                name = py_text(parts[0])
                fields[name] = parts[1] if len(parts) == 2 else "Tuple(" + ",".join(parts[1:]) + ")"
            _append_leaf(output, engine, key, "PARITY_MATRIX", fields["THICKNESS_CENTRE_BASIS"])
            _append_leaf(output, engine, key, "FACE_BASIS", fields["FACE_BASIS"])
            _append_leaf(output, engine, key, "OFF_DIAGONAL_BLOCKS", fields["OFF_DIAGONAL_BLOCKS"])
        else:
            _append_leaf(output, engine, key, "PARITY_MATRIX", _wl_dict(record)["EXPRESSION"])
    return output


def extract_hermitian_like(engine: str, family: str, payload: str) -> list[ParsedCase]:
    leaf = "HERMITIAN_PART" if family == "DTN_HERMITIAN_PART" else "REACTIVE_PART"
    if engine == "WL" or family == "DTN_REACTIVE_PART":
        return _extract_case_map_primary(
            engine, payload, ("BRANCH", "FACE"), ("BRANCH", "FACE"), leaf=leaf
        )
    output: list[ParsedCase] = []
    for raw_key, record in _py_outer(payload):
        key = decode_py_key(raw_key, ("BRANCH", "FACE"))
        values = _py_named_dict(_py_named_dict(record)["VALUE"])
        _append_leaf(output, engine, key, leaf, values["HERMITIAN_PART"])
        for name in ("FULL_BARE_DTN", "TRUE_AREA_ADJOINT", "PAIRING_MEASURE_SOURCE"):
            _append_leaf(output, engine, key, name, values[name])
    return output


def extract_inertial(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    schema = ("BRANCH", "FACE", "SCENARIO")
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    for raw_key, record in outer:
        key = decode_py_key(raw_key, schema) if engine == "PY" else decode_wl_key(raw_key, schema)
        if engine == "PY":
            fields = _py_named_dict(_py_named_dict(record)["VALUE"])
            for name, raw in fields.items():
                _append_leaf(output, engine, key, name, raw)
        else:
            fields = _wl_dict(record)
            if "PER_FACE_OBJECT" in fields:
                _append_leaf(
                    output, engine, key, "M_ADD", _wl_dict(fields["PER_FACE_OBJECT"])["EXPRESSION"]
                )
            for name in ("DEFINING_RELATION", "OUTWARD_ACCELERATION"):
                if name in fields:
                    _append_leaf(output, engine, key, name, fields[name])
    return output


def extract_grazing(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for raw_key, record in _py_outer(payload):
            key = decode_py_key(raw_key, ("BRANCH", "FACE"))
            _flatten_named(output, engine, key, _py_named_dict(record)["VALUE"], "BEHAVIOUR")
    else:
        schema = ("BRANCH", "FACE", "REGIME_OUT", "REGIME_IN")
        for raw_key, record in _wl_outer(payload):
            key = decode_wl_key(raw_key, schema)
            fields = _wl_dict(record)
            _append_leaf(output, engine, key, "BEHAVIOUR", fields["BEHAVIOUR_OBJECT"])
            for name in ("COALESCENCE_RELATIONS", "STRICT_REST_FRAME_PARAMETER"):
                if name in fields:
                    _append_leaf(output, engine, key, name, fields[name])
    return output


def extract_term_origins(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for name_raw, record in _py_outer(payload):
            name = py_text(name_raw)
            key = make_key((("OBJECT", name),))
            fields = _py_named_dict(record)
            _flatten_named(output, engine, key, fields["VALUE"], "ORIGIN")
    else:
        for branch, record in _wl_outer(payload):
            fields = _wl_dict(record)
            origins = _wl_dict(fields["ORIGINS"])
            for name, package in origins.items():
                key = make_key((("OBJECT", name), ("BRANCH", branch)))
                package_fields = _wl_dict(package)
                _append_leaf(output, engine, key, "ORIGIN", package_fields["EXPRESSION"])
            if "ORIGIN_SUM_RESIDUAL" in fields:
                key = make_key((("OBJECT", "ORIGIN_SUM"), ("BRANCH", branch)))
                _flatten_named(output, engine, key, fields["ORIGIN_SUM_RESIDUAL"], "RESIDUAL")
    return output


def extract_face_response(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        top = _py_named_dict(payload)
        opaque = top["OPAQUE_MU_THETA_OPERATOR"]
        for raw_key, record in _py_outer(opaque):
            key = decode_py_key(raw_key, ("BRANCH", "DENSITY"))
            _append_leaf(
                output, engine, key, "OPAQUE_MU_THETA_OPERATOR",
                _py_named_dict(record)["VALUE"],
                note="mu_theta_operands_preserved_without_registry_fold",
            )
        for raw_key, record in _py_outer(top["CASES"]):
            key = decode_py_key(raw_key, ("BRANCH", "FACE", "DENSITY"))
            values = _py_named_dict(_py_named_dict(record)["VALUE"])
            for name, raw in values.items():
                canonical = FIELD_ALIAS.get(name, name)
                _flatten_named(output, engine, key, raw, canonical)
    else:
        for raw_key, record in _wl_outer(payload):
            tokens = raw_key.split("|")
            if tokens[0] == "PARITY":
                if len(tokens) != 3:
                    raise AxisError(f"FACE_RESPONSE parity view key has {len(tokens)} tokens")
                key = make_key(
                    (("OBJECT", "FACE_RESPONSE"), ("BRANCH", tokens[1]),
                     ("DENSITY", tokens[2]), ("CASE", "PARITY_VIEW"))
                )
            else:
                key = decode_wl_key(raw_key, ("BRANCH", "FACE", "DENSITY"))
            for name, raw in _wl_outer(record):
                if name in {"PRESSURE", "RELATIVE_FLUX", "TRACTION"}:
                    package = _maybe_named_pairs(engine, raw)
                    package_fields = dict(package) if package else {}
                    if "EXPRESSION" in package_fields:
                        _append_leaf(output, engine, key, name, package_fields["EXPRESSION"])
                    else:
                        _flatten_named(output, engine, key, raw, name)
                else:
                    _flatten_named(output, engine, key, raw, name)
    return output


def _unwrap_wl_expression(raw: str) -> str:
    fields = _maybe_named_pairs("WL", raw)
    if fields:
        values = dict(fields)
        if "EXPRESSION" in values:
            return values["EXPRESSION"]
    return raw


def extract_face_response_coeffs(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    schema = ("BRANCH", "FACE", "DENSITY")
    for raw_key, record in outer:
        key = decode_py_key(raw_key, schema) if engine == "PY" else decode_wl_key(raw_key, schema)
        if engine == "PY":
            fields = _py_named_dict(_py_named_dict(record)["VALUE"])
            for name in ("PRESSURE", "RELATIVE_FLUX", "TRACTION_NORMAL_AMPLITUDE"):
                rows = py_tuple_args(fields[name])
                for order, row in zip(("FLAT", "FIRST_SHAPE_KERNEL"), rows):
                    values = py_tuple_args(row)
                    for coefficient, raw in zip(
                        ("FACE_VELOCITY_COEFFICIENT", "MU_THETA_COEFFICIENT"), values
                    ):
                        base = "TRACTION" if name == "TRACTION_NORMAL_AMPLITUDE" else name
                        _append_leaf(output, engine, key, f"{base}.{order}.{coefficient}", raw)
            for name in ("RESOLVENT_SYMBOL_OUTPUT", "RESOLVENT_SYMBOL_INPUT"):
                _append_leaf(output, engine, key, name, fields[name])
        else:
            fields = _wl_dict(record)
            mapping = {
                "PRESSURE_FLAT": "PRESSURE.FLAT",
                "RELATIVE_FLUX_FLAT": "RELATIVE_FLUX.FLAT",
                "TRACTION_FLAT_VECTOR": "TRACTION.FLAT",
                "PRESSURE_FIRST_SHAPE_KERNEL": "PRESSURE.FIRST_SHAPE_KERNEL",
                "RELATIVE_FLUX_FIRST_SHAPE_KERNEL": "RELATIVE_FLUX.FIRST_SHAPE_KERNEL",
                "TRACTION_FIRST_SHAPE_KERNEL_VECTOR": "TRACTION.FIRST_SHAPE_KERNEL",
            }
            for source, target in mapping.items():
                if source not in fields:
                    continue
                children = _maybe_named_pairs(engine, fields[source])
                if children:
                    for coefficient, raw in children:
                        _append_leaf(
                            output, engine, key, f"{target}.{coefficient}",
                            _unwrap_wl_expression(raw),
                        )
                else:
                    _append_leaf(output, engine, key, target, fields[source])
    return output


def _flatten_coeff_object_py(
    output: list[ParsedCase], key: Key, raw: str, prefix: str
) -> None:
    fields = _py_named_dict(raw)
    for name in ("PRESSURE", "RELATIVE_FLUX", "TRACTION_NORMAL_AMPLITUDE"):
        if name not in fields:
            continue
        rows = py_tuple_args(fields[name])
        for order, row in zip(("FLAT", "FIRST_SHAPE_KERNEL"), rows):
            # The WL control shells sometimes expose the complete response
            # expression at this address rather than its two coefficients.
            # Keep that mechanically addressable container comparison as well
            # as every scalar coefficient below it.  A tuple/scalar shape
            # difference is a measured Mismatch, not a zero-extraction path.
            base = "TRACTION" if name == "TRACTION_NORMAL_AMPLITUDE" else name
            _append_leaf(output, "PY", key, f"{prefix}.{base}.{order}", row)
            values = py_tuple_args(row)
            for coefficient, item in zip(
                ("FACE_VELOCITY_COEFFICIENT", "MU_THETA_COEFFICIENT"), values
            ):
                _append_leaf(
                    output, "PY", key,
                    f"{prefix}.{base}.{order}.{coefficient}", item,
                )
    for name in ("RESOLVENT_SYMBOL_OUTPUT", "RESOLVENT_SYMBOL_INPUT"):
        if name in fields:
            _append_leaf(output, "PY", key, f"{prefix}.{name}", fields[name])


def _flatten_coeff_object_wl(
    output: list[ParsedCase], key: Key, raw: str, prefix: str
) -> None:
    fields = _wl_dict(raw)
    mapping = {
        "PRESSURE_FLAT": "PRESSURE.FLAT",
        "RELATIVE_FLUX_FLAT": "RELATIVE_FLUX.FLAT",
        "TRACTION_FLAT_VECTOR": "TRACTION.FLAT",
        "PRESSURE_FIRST_SHAPE_KERNEL": "PRESSURE.FIRST_SHAPE_KERNEL",
        "RELATIVE_FLUX_FIRST_SHAPE_KERNEL": "RELATIVE_FLUX.FIRST_SHAPE_KERNEL",
        "TRACTION_FIRST_SHAPE_KERNEL_VECTOR": "TRACTION.FIRST_SHAPE_KERNEL",
    }
    for source, target in mapping.items():
        if source not in fields:
            continue
        children = _maybe_named_pairs("WL", fields[source])
        if children:
            for coefficient, item in children:
                _append_leaf(
                    output, "WL", key,
                    f"{prefix}.{target}.{coefficient}", _unwrap_wl_expression(item),
                )
        else:
            try:
                values = _wl_list(fields[source])
            except Exception:
                values = []
            if len(values) == 2:
                for coefficient, item in zip(
                    ("FACE_VELOCITY_COEFFICIENT", "MU_THETA_COEFFICIENT"), values
                ):
                    _append_leaf(
                        output, "WL", key, f"{prefix}.{target}.{coefficient}", item
                    )
            else:
                _append_leaf(output, "WL", key, f"{prefix}.{target}", fields[source])
    for name, item in fields.items():
        if name not in mapping:
            _flatten_named(output, "WL", key, item, f"{prefix}.{name}")


def extract_noninvertibility(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for raw_key, record in _py_outer(payload):
            key = decode_py_key(raw_key, ("BRANCH", "FACE"))
            values = _py_named_dict(_py_named_dict(record)["VALUE"])
            operator_key = key_replace(key, LEAF="OPERATOR")
            output.append(
                raw_control_case(
                    engine, operator_key, values["OPERATOR"],
                    family="NONINVERTIBILITY_CONDITION", leaf="OPERATOR",
                )
            )
            aliases = {
                "NONINVERTIBILITY": "NONINVERTIBILITY_RELATION",
                "FLAT_DIAGONAL_SYMBOL": "FLAT_DIAGONAL_SYMBOL_RELATION",
            }
            for name, raw in values.items():
                if name != "OPERATOR":
                    _append_leaf(output, engine, key, aliases.get(name, name), raw)
    else:
        for branch, record in _wl_outer(payload):
            key = make_key((("BRANCH", branch),))
            fields = _wl_dict(record)
            output.append(
                raw_control_case(
                    engine, key_replace(key, LEAF="OPERATOR"), fields["OPERATOR"],
                    family="NONINVERTIBILITY_CONDITION", leaf="OPERATOR",
                )
            )
            for name, raw in fields.items():
                if name != "OPERATOR":
                    _flatten_named(output, engine, key, raw, name)
    return output


def _wl_list(raw: str) -> list[str]:
    value = raw.strip()
    if not (value.startswith("{") and value.endswith("}")):
        raise InputError(f"expected WL list, got {value[:80]}")
    body = value[1:-1].strip()
    return [] if not body else split_top(body)


def extract_locus(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    root = make_key((("OBJECT", family),))
    suffix = family.removeprefix("DEGENERATE_LOCI_")
    if suffix == "EQUATIONS":
        if engine == "PY":
            _append_leaf(output, engine, root, "EQUATIONS", payload)
        else:
            for name, raw in _wl_dict(payload).items():
                _flatten_named(output, engine, root, raw, name)
        return output
    if suffix == "SOLUTION":
        if engine == "PY":
            variable, solution = py_tuple_args(payload)
            _append_leaf(output, engine, root, "VARIABLES", f"Tuple({variable})")
            _append_leaf(output, engine, root, "SOLUTION", solution)
        else:
            fields = _wl_dict(payload)
            for name, raw in fields.items():
                _flatten_named(output, engine, root, raw, name)
        return output
    if suffix in {"IDENTICALLY_SATISFIED", "INCONSISTENT"}:
        if engine == "PY":
            status, test, operands = py_tuple_args(payload)
            _append_leaf(output, engine, root, "STATUS_TOKEN", status)
            _append_leaf(output, engine, root, "TEST_OBJECT", test)
            _append_leaf(output, engine, root, "OPERANDS", operands)
        else:
            fields = _wl_dict(payload)
            for name, raw in fields.items():
                _flatten_named(output, engine, root, raw, name)
        return output
    if suffix == "REAL_ADMISSIBLE":
        py_records = py_tuple_args(payload) if engine == "PY" else []
        wl_records = _wl_list(payload) if engine == "WL" else []
        records = py_records if engine == "PY" else wl_records
        for index, record in enumerate(records):
            # The committed streams each carry one real-admissible case.  Do
            # not make reordered multi-case records agree by ordinal: beyond
            # the singleton shape the raw branch representation is deliberately
            # an engine-specific CASE token and therefore remains unpaired.
            if len(records) == 1:
                case_token = "SINGLE_REAL_ADMISSIBLE"
            elif engine == "PY":
                case_token = f"PY_BRANCH_{py_tuple_args(record)[0]}"
            else:
                case_token = f"WL_BRANCH_{_wl_dict(record).get('BRANCH', index)}"
            key = key_replace(root, CASE=case_token)
            if engine == "PY":
                items = py_tuple_args(record)
                if len(items) != 4:
                    raise InputError("PY real-admissible record is not four fields")
                for name, raw in zip(("BRANCH", "STATUS_TOKEN", "TEST_OBJECT", "OPERANDS"), items):
                    _append_leaf(output, engine, key, name, raw)
            else:
                fields = _wl_dict(record)
                for name, raw in fields.items():
                    _flatten_named(output, engine, key, raw, name)
        return output
    raise InputError(f"unknown locus family {family}")


def extract_port(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if family == "PERMEABLE_PORT_HERMITIAN":
        py_schema = ("BRANCH", "PARITY", "REGIME_OUT", "REGIME_IN")
        wl_schema = ("BRANCH", "DENSITY", "REGIME_OUT", "REGIME_IN", "PARITY")
    else:
        py_schema = (
            "BRANCH", "PARITY", "REGIME_OUT", "REGIME_IN", "PORT", "LIMIT"
        )
        wl_schema = ("BRANCH", "DENSITY", "REGIME_OUT", "REGIME_IN", "PARITY")
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    for raw_key, record in outer:
        key = decode_py_key(raw_key, py_schema) if engine == "PY" else decode_wl_key(raw_key, wl_schema)
        fields = _py_named_dict(record) if engine == "PY" else _wl_dict(record)
        handled: set[str] = set()
        if engine == "PY":
            value = fields.get("VALUE")
            if family == "PERMEABLE_PORT_HERMITIAN":
                fields = _py_named_dict(value) if value is not None else fields
                candidate = fields.get("BLOCK_HERMITIAN_FORM")
                if candidate is not None:
                    handled.add("BLOCK_HERMITIAN_FORM")
            else:
                candidate = value
                if value is not None:
                    handled.add("VALUE")
        else:
            candidate = fields.get("HERMITIAN_FORM", fields.get("PORT_MATRIX"))
            if "HERMITIAN_FORM" in fields:
                handled.add("HERMITIAN_FORM")
            elif "PORT_MATRIX" in fields:
                handled.add("PORT_MATRIX")
        if candidate is not None:
            _flatten_named(output, engine, key, candidate, "HERMITIAN_FORM")
        for name in ("SIGN_TEST_OBJECTS", "STATUS_TOKEN", "PORT_MATRIX"):
            if name in fields and not (name == "PORT_MATRIX" and candidate == fields[name]):
                _flatten_named(output, engine, key, fields[name], name)
                handled.add(name)
        if engine == "WL" and family == "PERMEABLE_DISSIPATION_VS_OMEGA_TAU":
            memory = fields.get("INDEPENDENT_MEMORY_LIMITS")
            if memory is not None:
                handled.add("INDEPENDENT_MEMORY_LIMITS")
                for tau_name, limits_raw in _wl_outer(memory):
                    match = re.fullmatch(r"tau([AVX])", tau_name)
                    if match is None:
                        raise AxisError(f"untyped memory port {tau_name}")
                    port_key = key_replace(key, PORT=match.group(1))
                    limits = _wl_dict(limits_raw)
                    for source, limit_name in (
                        ("ZERO_LIMIT", "OMEGA_TAU_TO_ZERO"),
                        ("INFINITE_LIMIT", "OMEGA_TAU_TO_INFINITY"),
                    ):
                        if source in limits:
                            _flatten_named(
                                output, engine, key_replace(port_key, LIMIT=limit_name),
                                limits[source], "HERMITIAN_FORM",
                            )
                    if "OMEGA_TAU_OPERAND" in limits:
                        _flatten_named(
                            output, engine, port_key, limits["OMEGA_TAU_OPERAND"],
                            "OMEGA_TAU_OPERAND",
                        )
        # Preserve every remaining construction/corruption/bookkeeping operand.
        # This is essential for the WL parity-construction controls and the PY
        # FULL_BARE_DTN/MEMORY_HERMITIAN_FORM leaves, which have no primary
        # PORT_MATRIX in several committed cases.
        for name, raw in fields.items():
            if name not in handled:
                _flatten_named(output, engine, key, raw, name)
    return output


def extract_energy(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    for raw_key, record in outer:
        key = (
            decode_py_key(raw_key, ("BRANCH", "FACE", "SCENARIO"))
            if engine == "PY"
            else decode_wl_key(raw_key, ("BRANCH", "FACE"))
        )
        if family == "ENERGY_FACE_TRACTION_OPERAND":
            if engine == "PY":
                values = _py_named_dict(_py_named_dict(record)["VALUE"])
                for name, raw in values.items():
                    _append_leaf(output, engine, key, name, raw)
            else:
                fields = _wl_dict(record)
                pairing = _wl_dict(fields["FACE_TRACTION_PAIRING"])["EXPRESSION"]
                _append_leaf(output, engine, key, "BASELINE", pairing)
                for name, raw in fields.items():
                    if name != "FACE_TRACTION_PAIRING":
                        _flatten_named(output, engine, key, raw, name)
        elif family == "ENERGY_BULK_FARFIELD_FLUX_OPERAND":
            if engine == "PY":
                values = _py_named_dict(_py_named_dict(record)["VALUE"])
                for scenario, raw in values.items():
                    parts = py_tuple_args(raw)
                    if len(parts) >= 4:
                        _append_leaf(output, engine, key, f"{scenario}.OUTGOING_FLUX", parts[1])
                        _append_leaf(
                            output, engine, key,
                            f"{scenario}.COMPARISON_ORIENTED_OPERAND", parts[3]
                        )
                    else:
                        _append_leaf(output, engine, key, scenario, raw)
            else:
                fields = _wl_dict(record)
                for name in ("OUTGOING_FLUX", "COMPARISON_ORIENTED_OPERAND"):
                    if name in fields:
                        _append_leaf(output, engine, key, f"BASELINE.{name}", _wl_dict(fields[name])["EXPRESSION"])
                for name, raw in fields.items():
                    if name not in {"OUTGOING_FLUX", "COMPARISON_ORIENTED_OPERAND"}:
                        _flatten_named(output, engine, key, raw, name)
        else:
            if engine == "PY":
                values = _py_named_dict(_py_named_dict(record)["VALUE"])
                for name, raw in values.items():
                    _append_leaf(output, engine, key, name, raw)
            else:
                fields = _wl_dict(record)
                for name in ("OPERAND_A", "OPERAND_B"):
                    if name in fields:
                        _append_leaf(output, engine, key, name, fields[name])
                residual_raw = fields["RESIDUAL_A_MINUS_B"]
                residual_fields = _maybe_named_pairs(engine, residual_raw)
                _append_leaf(
                    output, engine, key, "BASELINE",
                    dict(residual_fields)["EXPRESSION"] if residual_fields else residual_raw,
                )
                for name, raw in fields.items():
                    if name not in {"OPERAND_A", "OPERAND_B", "RESIDUAL_A_MINUS_B"}:
                        _flatten_named(output, engine, key, raw, name)

    return output


def _decode_py_object_key(raw_key: str) -> Key:
    tokens = py_key_tokens(raw_key)
    if len(tokens) == 2:
        return make_key((("OBJECT", tokens[0]), ("BRANCH", tokens[1])))
    if len(tokens) == 3:
        return make_key(
            (("OBJECT", tokens[0]), ("BRANCH", tokens[1]), ("DENSITY", tokens[2]))
        )
    raise AxisError(f"object/branch key has {len(tokens)} tokens: {tokens}")


def _flatten_py_face_object(
    output: list[ParsedCase], engine: str, key: Key, raw: str, prefix: str
) -> None:
    for face, item in _py_face_pairs(raw):
        face_key = key_replace(key, FACE=face)
        if key_dict(key).get("OBJECT") == "FACE_RESPONSE":
            _flatten_coeff_object_py(output, face_key, item, prefix)
            continue
        fields = _maybe_named_pairs("PY", item)
        values = dict(fields) if fields else {}
        if "VALUE" in values:
            _append_leaf(output, engine, face_key, prefix, values["VALUE"])
            for name, child in fields:
                if name != "VALUE":
                    _flatten_named(output, engine, face_key, child, f"{prefix}.{name}")
        else:
            _flatten_named(output, engine, face_key, item, prefix)


def _flatten_wl_face_object(
    output: list[ParsedCase], key: Key, raw: str, prefix: str
) -> bool:
    try:
        pairs = _wl_face_pairs(raw)
    except Exception:
        return False
    for face, item in pairs:
        _flatten_named(output, "WL", key_replace(key, FACE=face), item, prefix)
    return True


def extract_rep_invariance(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for raw_key, value in _py_outer(payload):
            key = _decode_py_object_key(raw_key)
            prefix = "RESIDUAL" if family.endswith("RESIDUAL") else "OPERAND"
            _flatten_py_face_object(output, engine, key, value, prefix)
        return output
    for branch, record in _wl_outer(payload):
        fields = _wl_dict(record)
        if family.endswith("RESIDUAL"):
            if "DTN_RESIDUAL_A_MINUS_B" in fields:
                key = make_key(
                    (("OBJECT", "DTN"), ("BRANCH", branch), ("FACE", "1"))
                )
                _append_leaf(output, engine, key, "RESIDUAL", fields["DTN_RESIDUAL_A_MINUS_B"])
            response = fields.get("FACE_RESPONSE_RESIDUAL_A_MINUS_B")
            if response is not None:
                for raw_key, item in _wl_outer(response):
                    key = decode_wl_key(raw_key, ("FACE", "DENSITY"))
                    key = make_key(
                        (("OBJECT", "FACE_RESPONSE"), ("BRANCH", branch), *key)
                    )
                    _flatten_coeff_object_wl(output, key, item, "RESIDUAL")
        else:
            if "DTN" in fields:
                key = make_key(
                    (("OBJECT", "DTN"), ("BRANCH", branch), ("FACE", "1"))
                )
                dtn = _wl_dict(fields["DTN"])
                _append_leaf(output, engine, key, "OPERAND", dtn.get("EXPRESSION", fields["DTN"]))
            response = fields.get("FACE_RESPONSE")
            if response is not None:
                for raw_key, item in _wl_outer(response):
                    key = decode_wl_key(raw_key, ("FACE", "DENSITY"))
                    key = make_key(
                        (("OBJECT", "FACE_RESPONSE"), ("BRANCH", branch), *key)
                    )
                    _flatten_coeff_object_wl(output, key, item, "OPERAND")
    return output


def extract_independence(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    prefix = "RESIDUAL" if family.endswith("RESIDUAL") else "OPERAND"
    if engine == "PY":
        for raw_key, value in _py_outer(payload):
            key = _decode_py_object_key(raw_key)
            _flatten_py_face_object(output, engine, key, value, prefix)
        return output
    for branch, record in _wl_outer(payload):
        fields = _wl_dict(record)
        dtn_name = "DTN_UPPER_FACE"
        if dtn_name in fields:
            key = make_key((("OBJECT", "DTN"), ("BRANCH", branch), ("FACE", "1")))
            _append_leaf(output, engine, key, prefix, fields[dtn_name])
        response = fields.get("FACE_RESPONSE")
        if response is not None:
            for density, item in _wl_outer(response):
                key = make_key(
                    (("OBJECT", "FACE_RESPONSE"), ("BRANCH", branch),
                     ("FACE", "1"), ("DENSITY", density))
                )
                _flatten_coeff_object_wl(output, key, item, prefix)
        for name, raw in fields.items():
            if name not in {dtn_name, "FACE_RESPONSE"}:
                key = make_key((("OBJECT", name), ("BRANCH", branch), ("FACE", "1")))
                _flatten_named(output, engine, key, raw, prefix)
    return output


def _control_form_wl_object(
    output: list[ParsedCase], key: Key, raw: str, prefix: str, *,
    coefficient_object: bool,
) -> None:
    fields = _maybe_named_pairs("WL", raw)
    if not fields:
        _append_leaf(output, "WL", key, prefix, raw)
        return
    names = {name for name, _ in fields}
    if names and all(name.startswith("FACE_") for name in names):
        for face, item in _wl_face_pairs(raw):
            face_key = key_replace(key, FACE=face)
            if coefficient_object:
                _flatten_coeff_object_wl(output, face_key, item, prefix)
            else:
                _flatten_named(output, "WL", face_key, item, prefix)
        return
    for name, child in fields:
        if name.startswith("FACE_"):
            face = key_dict(decode_wl_key(name, ("FACE",)))["FACE"]
            _flatten_named(output, "WL", key_replace(key, FACE=face), child, prefix)
        else:
            _flatten_named(output, "WL", key, child, f"{prefix}.{name}")


def extract_control_form(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    schema = ("OBJECT", "BRANCH", "DENSITY", "DIRECTION")
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    for raw_key, record in outer:
        key = decode_py_key(raw_key, schema) if engine == "PY" else decode_wl_key(raw_key, schema)
        if engine == "PY":
            prefix = "RESIDUAL" if family.endswith("RESIDUAL") else "OBJECT"
            if key_dict(key)["OBJECT"] == "DTN":
                # WL emits one scalar kernel here, while PY emits a FACE-keyed
                # container.  Compare those actual outer objects at the parent
                # key so residual() reports the structural disagreement; keep
                # every PY FACE scalar below as a separate unpaired finding.
                # Assigning the WL scalar to FACE=1 would pre-register the very
                # face-selection difference the frozen contract seals.
                _append_leaf(
                    output, engine, key, prefix, record,
                    note="face_keyed_container_vs_unkeyed_scalar",
                )
            for face, item in _py_face_pairs(record):
                face_key = key_replace(key, FACE=face)
                if key_dict(key)["OBJECT"] == "FACE_RESPONSE":
                    _flatten_coeff_object_py(output, face_key, item, prefix)
                else:
                    _flatten_named(output, engine, face_key, item, prefix)
            continue
        fields = _wl_dict(record)
        if family.endswith("RESIDUAL"):
            selected = (
                ("OPERAND_A", "OPERAND_A"),
                ("OPERAND_B", "OPERAND_B"),
                ("RESIDUAL_A_MINUS_B", "RESIDUAL"),
            )
            for field_name, prefix in selected:
                if field_name not in fields:
                    continue
                container = fields[field_name]
                children = _wl_dict(container)
                for name, raw in children.items():
                    target_prefix = prefix if name == "OBJECT" else f"{prefix}.{name}"
                    _control_form_wl_object(
                        output, key, raw, target_prefix,
                        coefficient_object=(
                            name == "OBJECT" and key_dict(key)["OBJECT"] == "FACE_RESPONSE"
                        ),
                    )
        else:
            for name, raw in fields.items():
                _control_form_wl_object(
                    output, key, raw, name,
                    coefficient_object=(
                        name == "OBJECT" and key_dict(key)["OBJECT"] == "FACE_RESPONSE"
                    ),
                )
    return output


def extract_uniform(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    py_schema = ("OBJECT", "REGIME_OUT", "REGIME_IN", "PARITY")
    wl_schema = ("OBJECT", "REGIME_OUT", "REGIME_IN", "PARITY", "DENSITY")
    outer = _py_outer(payload) if engine == "PY" else _wl_outer(payload)
    for raw_key, value in outer:
        key = decode_py_key(raw_key, py_schema) if engine == "PY" else decode_wl_key(raw_key, wl_schema)
        _flatten_named(output, engine, key, value, "OPERAND")
    return output


def extract_zero_jet(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    key = make_key((("OBJECT", "ZERO_JET"),))
    if engine == "PY":
        fields = _py_named_dict(payload)
        _flatten_named(output, engine, key, fields["VALUE"], "OPERAND")
        for name, raw in fields.items():
            if name != "VALUE":
                _flatten_named(output, engine, key, raw, name)
    else:
        for branch, raw in _wl_outer(payload):
            _flatten_named(output, engine, key_replace(key, BRANCH=branch), raw, "OPERAND")
    return output


def _py_branch_key(raw_key: str, mutation: str) -> Key:
    tokens = py_key_tokens(raw_key)
    if tokens[0] == "RADIATION_REAL_AXIS":
        if len(tokens) != 2:
            raise AxisError(f"radiation branch key has {len(tokens)} tokens")
        return make_key(
            (("OBJECT", tokens[0]), ("MUTATION", mutation), ("REGIME_OUT", tokens[1]))
        )
    if len(tokens) not in {2, 3}:
        raise AxisError(f"control branch key has {len(tokens)} tokens: {tokens}")
    items: list[tuple[str, str]] = [
        ("OBJECT", tokens[0]), ("MUTATION", mutation), ("BRANCH", tokens[1])
    ]
    if len(tokens) == 3:
        items.extend((("DENSITY", tokens[2]), ("FACE", "1")))
    return make_key(items)


def extract_control_branch(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        if family == "CONTROL_BRANCH_RESIDUAL":
            for raw_key, raw in _py_outer(payload):
                tokens = py_key_tokens(raw_key)
                mutation = tokens[0]
                base_raw = "Tuple(" + ", ".join(
                    f"Str('{token}')" for token in tokens[1:]
                ) + ")"
                key = _py_branch_key(base_raw, mutation)
                if key_dict(key).get("OBJECT") == "FACE_RESPONSE":
                    _flatten_coeff_object_py(output, key, raw, "RESIDUAL")
                else:
                    _flatten_named(output, engine, key, raw, "RESIDUAL")
        else:
            mutation = _mutation_for_family(family)
            for raw_key, raw in _py_outer(payload):
                key = _py_branch_key(raw_key, mutation)
                if key_dict(key).get("OBJECT") == "FACE_RESPONSE":
                    _flatten_coeff_object_py(output, key, raw, "OPERAND")
                else:
                    _flatten_named(output, engine, key, raw, "OPERAND")
        return output
    if family == "CONTROL_BRANCH_RESIDUAL":
        for mutation, branches in _wl_outer(payload):
            for branch, objects in _wl_outer(branches):
                for object_name, raw in _wl_outer(objects):
                    if object_name == "FACE_RESPONSE":
                        for raw_key, item in _wl_outer(raw):
                            face_density = decode_wl_key(raw_key, ("FACE", "DENSITY"))
                            key = make_key(
                                (("OBJECT", object_name), ("MUTATION", mutation),
                                 ("BRANCH", branch), *face_density)
                            )
                            _flatten_coeff_object_wl(output, key, item, "RESIDUAL")
                    else:
                        key = make_key(
                            (("OBJECT", object_name), ("MUTATION", mutation), ("BRANCH", branch))
                        )
                        _flatten_named(output, engine, key, raw, "RESIDUAL")
    else:
        mutation = _mutation_for_family(family)
        for branch, objects in _wl_outer(payload):
            for object_name, raw in _wl_outer(objects):
                if object_name == "FACE_RESPONSE":
                    for raw_key, item in _wl_outer(raw):
                        face_density = decode_wl_key(raw_key, ("FACE", "DENSITY"))
                        key = make_key(
                            (("OBJECT", object_name), ("MUTATION", mutation),
                             ("BRANCH", branch), *face_density)
                        )
                        _flatten_coeff_object_wl(output, key, item, "OPERAND")
                else:
                    key = make_key(
                        (("OBJECT", object_name), ("MUTATION", mutation), ("BRANCH", branch))
                    )
                    _flatten_named(output, engine, key, raw, "OPERAND")
    return output


def extract_dimensions(engine: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for name_raw, vector in _py_outer(payload):
            key = make_key((("OBJECT", py_text(name_raw)),))
            _append_leaf(output, engine, key, "DIMENSION_L_T_M", vector)
    else:
        fields = _wl_dict(payload)
        for name, vector in _wl_outer(fields["OBJECT_DIMENSIONS"]):
            key = make_key((("OBJECT", name),))
            _append_leaf(output, engine, key, "DIMENSION_L_T_M", vector)
        root = make_key((("OBJECT", "DIMENSION_DERIVATION"),))
        for name, equation in _wl_outer(fields["DERIVATION_EQUATIONS"]):
            _append_leaf(output, engine, root, name, equation)
        for name, ratio in _wl_outer(fields.get("RATIO_DEFINITIONS", "<||>")):
            key = make_key((("OBJECT", name),))
            _append_leaf(output, engine, key, "RATIO_DEFINITION", ratio)
    return output


def extract_homogeneity(engine: str, family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    if engine == "PY":
        for name_raw, raw in _py_outer(payload):
            name = py_text(name_raw)
            key = make_key((("OBJECT", name),))
            if family.endswith("RESIDUAL"):
                _append_leaf(output, engine, key, "DIMENSION_L_T_M", raw)
            else:
                expression, vector = py_tuple_args(raw)
                _append_leaf(output, engine, key, "EXPRESSION", expression)
                _append_leaf(output, engine, key, "DIMENSION_L_T_M", vector)
        return output
    root = make_key((("OBJECT", "HOMOGENEITY"),))
    fields = _wl_dict(payload)
    for name, raw in fields.items():
        _flatten_named(output, engine, root, raw, name)
    return output


ENERGY_FAMILIES = {
    "ENERGY_FACE_TRACTION_OPERAND",
    "ENERGY_BULK_FARFIELD_FLUX_OPERAND",
    "ENERGY_RESIDUAL",
}
LOCUS_FAMILIES = {
    "DEGENERATE_LOCI_EQUATIONS",
    "DEGENERATE_LOCI_SOLUTION",
    "DEGENERATE_LOCI_IDENTICALLY_SATISFIED",
    "DEGENERATE_LOCI_INCONSISTENT",
    "DEGENERATE_LOCI_REAL_ADMISSIBLE",
}


def extract_engine(engine: str, family: str, payload: str) -> list[ParsedCase]:
    if family == "DTN_FLAT_SYMBOL":
        return extract_dtn_flat(engine, payload)
    if family == "DTN_OPERATOR":
        return extract_dtn_operator(engine, payload)
    if family == "DTN_KERNEL":
        return extract_dtn_kernel(engine, payload)
    if family in {"DTN_RIGID_SHIFT_OPERAND", "DTN_RIGID_SHIFT_RESIDUAL"}:
        return extract_rigid(engine, family, payload)
    if family == "DTN_BY_REGIME_PAIR":
        return extract_regime_pair(engine, payload)
    if family == "DTN_BY_PARITY":
        return extract_parity(engine, payload)
    if family in {"DTN_HERMITIAN_PART", "DTN_REACTIVE_PART"}:
        return extract_hermitian_like(engine, family, payload)
    if family == "DTN_INERTIAL_LOADING":
        return extract_inertial(engine, payload)
    if family == "DTN_GRAZING_BEHAVIOUR":
        return extract_grazing(engine, payload)
    if family == "DTN_TERM_ORIGINS":
        return extract_term_origins(engine, payload)
    if family == "FACE_RESPONSE":
        return extract_face_response(engine, payload)
    if family == "FACE_RESPONSE_COEFFS":
        return extract_face_response_coeffs(engine, payload)
    if family == "NONINVERTIBILITY_CONDITION":
        return extract_noninvertibility(engine, payload)
    if family in LOCUS_FAMILIES:
        return extract_locus(engine, family, payload)
    if family.startswith("PERMEABLE_"):
        return extract_port(engine, family, payload)
    if family in ENERGY_FAMILIES:
        return extract_energy(engine, family, payload)
    if family.startswith("REP_INVARIANCE_"):
        return extract_rep_invariance(engine, family, payload)
    if family.startswith("CONTROL_INDEPENDENCE_"):
        return extract_independence(engine, family, payload)
    if family.startswith("CONTROL_FORM_"):
        return extract_control_form(engine, family, payload)
    if family.startswith("UNIFORM_LIMIT_"):
        return extract_uniform(engine, payload)
    if family.startswith("ZERO_JET_"):
        return extract_zero_jet(engine, payload)
    if family.startswith("CONTROL_BRANCH_"):
        return extract_control_branch(engine, family, payload)
    if family == "DIMENSIONS":
        return extract_dimensions(engine, payload)
    if family.startswith("HOMOGENEITY_"):
        return extract_homogeneity(engine, family, payload)
    # A deliberately tiny hook for synthetic mechanics tests.  Measured c1
    # families are all dispatched above; adding a measured family requires a
    # payload-specific extractor, never a generic production fallback.
    if family == "SYNTHETIC":
        key = make_key((("OBJECT", "SYNTHETIC"), ("LEAF", "EXPRESSION")))
        return [py_case(key, payload) if engine == "PY" else wl_case(key, payload)]
    raise InputError(f"no bespoke c1 extractor registered for {family}")


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
    if isinstance(value, ResidualFailure):
        return f"ResidualFailure(reason={value.reason!r})"
    if isinstance(value, BooleanNotResidualable):
        return "BooleanNotResidualable()"
    if isinstance(value, Mismatch):
        return f"Mismatch(kind={value.kind!r}, detail={value.detail!r})"
    if isinstance(value, UndecidedResidual):
        return f"UndecidedResidual(reason={value.reason!r})"
    if isinstance(value, ResidualAssociation):
        body = ", ".join(f"{key!r}: {serialise(item)}" for key, item in value.entries)
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
        preview = case.raw if len(case.raw) <= 320 else case.raw[:320] + f"...<nchar={len(case.raw)}>"
        return f"<PARSE_FAILED {case.error}; RAW={preview}>"
    return serialise(case.value)


def emit_case(
    family: str,
    key: Key,
    py_item: ParsedCase | None,
    wl_item: ParsedCase | None,
    difference: object,
    *,
    note: str | None = None,
) -> None:
    print(f"CASE family={family} key={serialise_key(key)}", flush=True)
    print(f"operand_A = {display_operand(py_item)}", flush=True)
    print(f"operand_B = {display_operand(wl_item)}", flush=True)
    print(f"A_minus_B = {serialise(difference)}", flush=True)
    notes = [
        item
        for item in (
            note,
            py_item.note if py_item else None,
            wl_item.note if wl_item else None,
        )
        if item
    ]
    if notes:
        print(f"case_note = {' | '.join(dict.fromkeys(notes))}", flush=True)


def exact_rational_residual(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.together(sp.expand(left - right)))


def _replace_budget_failures(value: object, left: object, right: object) -> object:
    if (
        isinstance(value, ResidualFailure)
        and value.reason == "RESIDUAL_BUDGET_EXCEEDED"
        and isinstance(left, sp.Expr)
        and isinstance(right, sp.Expr)
    ):
        return exact_rational_residual(left, right)
    if isinstance(value, ResidualAssociation) and isinstance(left, Association) and isinstance(right, Association):
        left_items = left.as_dict()
        right_items = right.as_dict()
        return ResidualAssociation(
            tuple(
                (key, _replace_budget_failures(item, left_items[key], right_items[key]))
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
            _replace_budget_failures(item, a, b)
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
    return value if raw_chars > MAX_RATIONAL_CANON_RAW_CHARS else _replace_budget_failures(value, left, right)


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


def axis_value_mismatch_detail(py_key: Key, wl_key: Key) -> str | None:
    py_axes = key_dict(py_key)
    wl_axes = key_dict(wl_key)
    if set(py_axes) != set(wl_axes):
        return None
    for identity_axis in ("OBJECT", "LEAF"):
        if identity_axis in py_axes and py_axes[identity_axis] != wl_axes[identity_axis]:
            return None
    differences = [
        axis
        for axis in py_axes
        if axis not in {"OBJECT", "LEAF"} and py_axes[axis] != wl_axes[axis]
    ]
    if not differences:
        return None
    return "literal axis value mismatch " + ",".join(
        f"{axis}:{py_axes[axis]}!={wl_axes[axis]}" for axis in differences
    )


def compare_family(family: str, cases: FamilyCases, *, budget: float) -> Accounting:
    accounting = Accounting(
        extracted_leaves_py=len(cases.py), extracted_leaves_wl=len(cases.wl)
    )
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
            difference = typed_residual(
                py_item.compare_value,
                wl_item.compare_value,
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
                family, key, py_item, wl_item, TextAtom("UNDEFINED_DUPLICATE_KEY"),
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
            if detail is not None:
                mismatch_py.add(py_key)
                mismatch_wl.add(wl_key)
            else:
                detail = axis_value_mismatch_detail(py_key, wl_key)
            if detail is not None:
                reasons_by_key[("PY", py_key)].add(detail)
                reasons_by_key[("WL", wl_key)].add(detail)
                accounting.reasons.add(detail)
    accounting.axis_set_mismatch = sum(len(py_by_key[key]) for key in mismatch_py) + sum(
        len(wl_by_key[key]) for key in mismatch_wl
    )

    for engine, keys, by_key in (
        ("PY", py_unmatched, py_by_key), ("WL", wl_unmatched, wl_by_key)
    ):
        for key in sorted(keys, key=serialise_key):
            reasons = sorted(reasons_by_key.get((engine, key), ()))
            note = (
                "axis_set_mismatch: " + " | ".join(reasons)
                if reasons else f"{engine.lower()}_only"
            )
            for case in by_key[key]:
                accounting.parse_failed += int(materialize(case))
                emit_case(
                    family, key,
                    case if engine == "PY" else None,
                    case if engine == "WL" else None,
                    TextAtom("UNDEFINED_UNJOINED"), note=note,
                )
                release_case(case)

    accounting.reasons.update(cases.extraction_notes)
    reason_text = " | ".join(sorted(accounting.reasons)) if accounting.reasons else "none"
    print(
        f"ACCOUNTING {family} "
        "{"
        f"join={accounting.join}, py_only={accounting.py_only}, wl_only={accounting.wl_only}, "
        f"duplicate_key={accounting.duplicate_key}, parse_failed={accounting.parse_failed}, "
        f"axis_set_mismatch={accounting.axis_set_mismatch}, "
        f"extracted_leaves_py={accounting.extracted_leaves_py}, "
        f"extracted_leaves_wl={accounting.extracted_leaves_wl}"
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
    global ACTIVE_C1_BARE_SYMBOL
    previous_active_spelling = dict(ACTIVE_C1_BARE_SYMBOL)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--py", type=Path, default=DEFAULT_PY)
    parser.add_argument("--wl", type=Path, default=DEFAULT_WL)
    parser.add_argument("--family", action="append")
    parser.add_argument("--residual-leaf-budget", type=float, default=0.1)
    arguments = parser.parse_args(argv)
    started = time.monotonic()
    try:
        py_tags = load_py(arguments.py)
        wl_tags = load_wl(arguments.wl)
        names = real_py_symbol_names(py_tags.values())
        mechanical = checked_mechanical_symbol_map(names)
        if REQUIRED_C1_PY_SYMBOLS <= names:
            active, _ = verify_active_spelling_map(py_tags.values())
        else:
            active = {}
        ACTIVE_C1_BARE_SYMBOL = dict(active)
    except Exception as error:
        print(f"OPERATIONAL_ERROR {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        return 2

    print(
        f"SPELLING_INJECTIVITY collisions=0 reserved_names={len(names)} "
        f"mechanical_spellings={len(mechanical)} active_c1_folds={sorted(active.items())}",
        flush=True,
    )
    emit_local_inventory(py_tags, wl_tags)
    families = sorted(
        name for name in set(py_tags) | set(wl_tags) if not name.startswith("LOCAL_")
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
        "MEASUREMENT_SCOPE supplied_unfalsifiable=sections_1_to_2,supplied_substrate,"
        "mu_theta_operand,supplied_bookkeeping residual_target=none",
        flush=True,
    )
    ACTIVE_C1_BARE_SYMBOL = previous_active_spelling
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
