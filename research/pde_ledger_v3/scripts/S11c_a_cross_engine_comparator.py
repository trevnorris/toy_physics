#!/usr/bin/env python3
"""S11c-a T7 cross-engine operand comparator.

This is a measurement instrument, not a verdict engine.  It reads the two
committed S11c-a tag streams, applies only the closed mechanical map declared
for S11c-a, and emits both operands and their typed recursive residual for
every joined case.  Missing, duplicate, malformed, and axis-incompatible
cases are emitted and included in per-family accounting.

The SymPy engine is deliberately not imported or run.  S11c_a_exports.py is
a ledger and is deliberately not used as a tag stream.
"""

from __future__ import annotations

import argparse
import functools
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
from sympy.core.symbol import Str

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# Only the parsers, parsed-container conversion, and typed recursive residual
# are reused.  The S11b comparison/classification/render/run path is not.
from S11b_cross_engine_comparator import (  # noqa: E402
    Association,
    TextAtom,
    _convert_parsed_containers,
    parse_sympy_payload,
    parse_wolfram_payload,
    residual,
)


DEFAULT_PY = SCRIPT_DIR / "out" / "S11c_a_interface_geometry_sympy_audit.out"
DEFAULT_WL = (
    SCRIPT_DIR.parent
    / "mathematica"
    / "out"
    / "S11c_a_interface_geometry_mathematica_audit.out"
)

PY_TAG_LINE = re.compile(r"^(PY_S11CA_[A-Z0-9_]+):\s?(.*)$")
WL_TAG_LINE = re.compile(r"^(WL_S11CA_[A-Z0-9_]+):")

FAMILIES = (
    "BACKGROUND_STATE",
    "ADMISSIBILITY_PREMISE",
    "FACE_MAP_LAB_HELD",
    "FACE_MAP_MATERIAL_ADVECTED",
    "BACKGROUND_DENSITY_MAP",
    "FACE_NORMAL",
    "CONORMAL_DERIV",
    "FACE_MEASURE_SHAPE_DERIV",
    "FACE_VELOCITY",
    "RELATIVE_FLUX",
    "KINEMATIC_BALANCE",
    "TRACTION",
    "VIRTUAL_WORK_SHAPE_DERIV",
    "FACE_SHIFT",
    "PROJECTION_SHAPE_DERIV",
    "PROJECTION_STATIC_OPERAND",
    "PROJECTION_DYNAMIC_OPERAND",
    "PROJECTION_RESIDUAL",
    "PROJECTION_TERM_ORIGINS",
    "VIRTUAL_CONSTRAINT",
    "EVOLUTION_MASS_BALANCE",
    "EVOLUTION_TERM_ORIGINS",
    "CLOSURE_SHAPE_DERIV",
    "REP_INVARIANCE_EULERIAN_OPERAND",
    "REP_INVARIANCE_MATERIAL_OPERAND",
    "REP_INVARIANCE_RESIDUAL",
    "CONTROL_INDEPENDENCE_BASE_OPERAND",
    "CONTROL_INDEPENDENCE_CORRUPTED_OPERAND",
    "CONTROL_INDEPENDENCE_RESIDUAL",
    "CONTROL_FORM_BASE_OPERAND",
    "CONTROL_FORM_ABLATED_OPERAND",
    "CONTROL_FORM_RESIDUAL",
    "UNIFORM_LIMIT_S11CA_OPERAND",
    "UNIFORM_LIMIT_S11B_OPERAND",
    "UNIFORM_LIMIT_RESIDUAL",
    "DIMENSIONS",
    "HOMOGENEITY_BASE_OPERAND",
    "HOMOGENEITY_CONTROL_OPERAND",
    "HOMOGENEITY_RESIDUAL",
)

AXIS_ORDER = (
    "OBJECT",
    "BRANCH",
    "DENSITY",
    "FACE",
    "DOF",
    "VDOF",
    "FIELD",
    "ORIGIN",
    "DIRECTION",
)
Key = tuple[tuple[str, str], ...]

BRANCHES = {"LAB_HELD", "MATERIAL_ADVECTED"}
DENSITIES = {"RHO4_CONSTANT", "RHOBR_CONSTANT"}
FACES = {"PLUS", "MINUS", "BOTH_FACES"}
DOFS = {"DELTA_W", "ZETA_C"}
DIRECTIONS = {"1", "2", "3"}

PY_SCHEMA: dict[str, tuple[str, ...]] = {
    "BACKGROUND_STATE": ("BRANCH", "DENSITY"),
    "BACKGROUND_DENSITY_MAP": ("DENSITY",),
    "FACE_NORMAL": ("BRANCH", "FACE", "DOF"),
    "CONORMAL_DERIV": ("BRANCH", "FACE", "DOF"),
    "FACE_MEASURE_SHAPE_DERIV": ("BRANCH", "FACE", "DOF"),
    "FACE_VELOCITY": ("BRANCH", "FACE", "DOF"),
    "RELATIVE_FLUX": ("BRANCH", "FACE", "DOF"),
    "KINEMATIC_BALANCE": ("BRANCH", "FACE", "DOF"),
    "TRACTION": ("BRANCH", "FACE", "DOF", "DENSITY"),
    "VIRTUAL_WORK_SHAPE_DERIV": ("BRANCH", "DOF", "VDOF", "DENSITY"),
    "FACE_SHIFT": ("BRANCH", "FACE", "DOF", "DENSITY"),
    "PROJECTION_SHAPE_DERIV": ("BRANCH", "DOF", "DENSITY"),
    "PROJECTION_STATIC_OPERAND": ("BRANCH", "DOF", "DENSITY"),
    "PROJECTION_DYNAMIC_OPERAND": ("BRANCH", "DOF", "DENSITY"),
    "PROJECTION_RESIDUAL": ("BRANCH", "DOF", "DENSITY"),
    "PROJECTION_TERM_ORIGINS": ("BRANCH", "DOF", "DENSITY"),
    "VIRTUAL_CONSTRAINT": ("BRANCH", "DOF", "DENSITY"),
    "EVOLUTION_MASS_BALANCE": ("BRANCH", "DOF", "DENSITY"),
    "EVOLUTION_TERM_ORIGINS": ("BRANCH", "DOF", "DENSITY"),
    "CLOSURE_SHAPE_DERIV": ("BRANCH", "FACE", "DOF", "DENSITY"),
}

CONTROL_FORM_SCHEMA: dict[str, tuple[str, ...]] = {
    **{
        name: ("OBJECT", "BRANCH", "FACE", "DOF", "DIRECTION")
        for name in (
            "FACE_NORMAL",
            "CONORMAL_DERIV",
            "FACE_MEASURE_SHAPE_DERIV",
            "FACE_VELOCITY",
            "RELATIVE_FLUX",
            "KINEMATIC_BALANCE",
        )
    },
    **{
        name: (
            "OBJECT",
            "BRANCH",
            "FACE",
            "DOF",
            "DENSITY",
            "DIRECTION",
        )
        for name in ("TRACTION", "FACE_SHIFT", "CLOSURE_SHAPE_DERIV")
    },
    "VIRTUAL_WORK_SHAPE_DERIV": (
        "OBJECT",
        "BRANCH",
        "FACE",
        "DOF",
        "VDOF",
        "DENSITY",
        "DIRECTION",
    ),
    **{
        name: (
            "OBJECT",
            "BRANCH",
            "FACE",
            "DOF",
            "DENSITY",
            "DIRECTION",
        )
        for name in (
            "PROJECTION_SHAPE_DERIV",
            "PROJECTION_STATIC_OPERAND",
            "PROJECTION_DYNAMIC_OPERAND",
            "PROJECTION_RESIDUAL",
            "PROJECTION_TERM_ORIGINS",
            "VIRTUAL_CONSTRAINT",
            "EVOLUTION_MASS_BALANCE",
            "EVOLUTION_TERM_ORIGINS",
        )
    },
}

OBJECT_RENAME = {
    "RELATIVE_FLUX_LAW": "RELATIVE_FLUX",
    "EVOLUTION_BALANCE": "EVOLUTION",
}

FACE_SHIFT_FIELDS = (
    "PRESSURE",
    "BULK_VELOCITY_X1",
    "BULK_VELOCITY_X2",
    "BULK_VELOCITY_X3",
    "BULK_VELOCITY_W",
    "BULK_DENSITY",
    "CURRENT_X1",
    "CURRENT_X2",
    "CURRENT_X3",
    "NORMAL_CURRENT",
)


class InputError(ValueError):
    """The committed streams do not satisfy their tag/container grammar."""


class AxisError(ValueError):
    """A case key is not a single-valued typed-axis key."""


@dataclass
class ParsedCase:
    key: Key
    value: object | None
    raw: str
    context: tuple[tuple[str, object], ...] = ()
    error: str | None = None


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


def load_py(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise InputError(f"PY input does not exist: {path}")
    output: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        match = PY_TAG_LINE.fullmatch(line)
        if match is None:
            raise InputError(f"{path}:{line_number}: not one PY S11c-a tagged row")
        full_tag, payload = match.groups()
        name = full_tag.removeprefix("PY_S11CA_")
        if name in output:
            raise InputError(f"{path}:{line_number}: duplicate PY tag {name}")
        output[name] = payload
    return output


def load_wl(path: Path) -> dict[str, str]:
    """Load and reassemble multi-line WL payloads up to the next WL tag."""
    if not path.is_file():
        raise InputError(f"WL input does not exist: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    starts = [index for index, line in enumerate(lines) if WL_TAG_LINE.match(line)]
    if not starts:
        raise InputError(f"{path}: no WL S11c-a tags")
    output: dict[str, str] = {}
    for ordinal, start in enumerate(starts):
        match = WL_TAG_LINE.match(lines[start])
        if match is None:
            raise InputError(f"{path}:{start + 1}: malformed WL tag")
        end = starts[ordinal + 1] if ordinal + 1 < len(starts) else len(lines)
        full_tag = match.group(1)
        name = full_tag.removeprefix("WL_S11CA_")
        first = lines[start].split(":", 1)[1]
        payload = (first + "".join(lines[start + 1 : end])).strip()
        if name in output:
            raise InputError(f"{path}:{start + 1}: duplicate WL tag {name}")
        output[name] = payload
    return output


def split_top(text: str, separator: str = ",") -> list[str]:
    """Split at top level across PY/WL brackets, Associations, and strings."""
    parts: list[str] = []
    depth = 0
    start = 0
    index = 0
    quote: str | None = None
    escaped = False
    while index < len(text):
        character = text[index]
        if quote is not None:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            index += 1
            continue
        if character in {"'", '"'}:
            quote = character
            index += 1
            continue
        pair = text[index : index + 2]
        if pair == "<|":
            depth += 1
            index += 2
            continue
        if pair == "|>":
            depth -= 1
            index += 2
            continue
        if character in "([{":
            depth += 1
        elif character in ")]}":
            depth -= 1
        elif depth == 0 and text.startswith(separator, index):
            parts.append(text[start:index].strip())
            start = index + len(separator)
            index += len(separator)
            continue
        index += 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def arrow(entry: str) -> tuple[str, str | None]:
    parts = split_top(entry, "->")
    if len(parts) == 2:
        return parts[0], parts[1]
    return entry.strip(), None


def py_tuple_args(raw: str) -> list[str]:
    value = raw.strip()
    if not (value.startswith("Tuple(") and value.endswith(")")):
        raise InputError(f"expected PY Tuple, got {value[:80]}")
    return split_top(value[len("Tuple(") : -1])


def py_pair(raw: str) -> tuple[str, str]:
    items = py_tuple_args(raw)
    if len(items) != 2:
        raise InputError(f"expected a two-item PY pair, got {len(items)}")
    return items[0], items[1]


def py_top_pairs(payload: str) -> list[tuple[str, str]]:
    return [py_pair(entry) for entry in py_tuple_args(payload)]


def py_text(raw: str) -> str:
    match = re.fullmatch(r"Str\('([^']*)'\)", raw.strip())
    if match is None:
        raise InputError(f"expected PY Str key, got {raw[:80]}")
    return match.group(1)


def py_key_tokens(raw: str) -> tuple[str, ...]:
    tokens = re.findall(r"Str\('([^']*)'\)|Integer\((-?\d+)\)", raw)
    return tuple(text or integer for text, integer in tokens)


def py_field(value_raw: str, name: str = "VALUE") -> str:
    for field_raw in py_tuple_args(value_raw):
        key_raw, item_raw = py_pair(field_raw)
        if py_text(key_raw) == name:
            return item_raw
    raise InputError(f"PY case has no {name} field")


def wl_assoc_pairs(raw: str) -> list[tuple[str, str]]:
    value = raw.strip()
    if not (value.startswith("<|") and value.endswith("|>")):
        raise InputError(f"expected WL Association, got {value[:80]}")
    output: list[tuple[str, str]] = []
    for entry in split_top(value[2:-2]):
        key_raw, item_raw = arrow(entry)
        if item_raw is None:
            raise InputError(f"WL Association entry has no top-level Rule: {entry[:80]}")
        match = re.fullmatch(r'"([^"]*)"', key_raw.strip())
        if match is None:
            raise InputError(f"WL Association key is not a string: {key_raw[:80]}")
        output.append((match.group(1), item_raw))
    return output


def wl_field(raw: str, *path: str) -> str:
    value = raw
    for name in path:
        fields = dict(wl_assoc_pairs(value))
        if name not in fields:
            raise InputError(f"WL Association has no {'.'.join(path)} leaf")
        value = fields[name]
    return value


def wl_rule_list_fields(raw: str) -> dict[str, str]:
    value = raw.strip()
    if value.startswith("<|") and value.endswith("|>"):
        return dict(wl_assoc_pairs(value))
    if not (value.startswith("{") and value.endswith("}")):
        raise InputError(f"expected WL Rule list, got {value[:80]}")
    entries = split_top(value[1:-1])
    if len(entries) == 2:
        labels = ("DYNAMIC_SHAPE_DERIVATIVE", "STATIC_SHAPE_DERIVATIVE")
        # CONTROL_FORM residuals are emitted as the positional difference of
        # two Rule lists.  Removing only the repeated Rule head preserves the
        # complete RHS algebra (including its sign) and restores the two typed
        # variant leaves.  The zero case is simply {0,0}.
        return {
            label: re.sub(rf'"{label}"\s*->\s*', "", entry)
            for label, entry in zip(labels, entries)
        }
    output: dict[str, str] = {}
    for entry in entries:
        key_raw, item_raw = arrow(entry)
        if item_raw is None:
            raise InputError(f"WL list entry has no Rule: {entry[:80]}")
        match = re.fullmatch(r'"([^"]*)"', key_raw.strip())
        if match is None:
            raise InputError(f"WL Rule-list key is not a string: {key_raw[:80]}")
        output[match.group(1)] = item_raw
    return output


def make_key(items: Iterable[tuple[str, str]]) -> Key:
    decoded: dict[str, str] = {}
    for axis, original in items:
        value = str(original)
        if axis in decoded:
            raise AxisError(f"two values for axis {axis}: {decoded[axis]} and {value}")
        if axis == "OBJECT":
            value = OBJECT_RENAME.get(value, value)
        elif axis == "FACE":
            value = {"1": "PLUS", "-1": "MINUS", "FACE_PLUS": "PLUS", "FACE_MINUS": "MINUS"}.get(
                value, value
            )
            if value not in FACES:
                raise AxisError(f"invalid FACE value {value}")
        elif axis == "BRANCH" and value not in BRANCHES:
            raise AxisError(f"invalid BRANCH value {value}")
        elif axis == "DENSITY" and value not in DENSITIES:
            raise AxisError(f"invalid DENSITY value {value}")
        elif axis in {"DOF", "VDOF"} and value not in DOFS:
            raise AxisError(f"invalid {axis} value {value}")
        elif axis == "DIRECTION" and value not in DIRECTIONS:
            raise AxisError(f"invalid DIRECTION value {value}; expected 1, 2, or 3")
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


def decode_py_key(raw: str, schema: Sequence[str], extra: Iterable[tuple[str, str]] = ()) -> Key:
    tokens = py_key_tokens(raw)
    if len(tokens) != len(schema):
        raise AxisError(f"PY schema {tuple(schema)} got {len(tokens)} tokens: {tokens}")
    return make_key((*extra, *zip(schema, tokens)))


def decode_wl_key(raw: str, *, object_axis: bool = False) -> Key:
    tokens = tuple(raw.split("|"))
    items: list[tuple[str, str]] = []
    for index, token in enumerate(tokens):
        if object_axis and index == 0:
            items.append(("OBJECT", token))
        elif token in BRANCHES:
            items.append(("BRANCH", token))
        elif token in DENSITIES:
            items.append(("DENSITY", token))
        elif token in {"FACE_PLUS", "FACE_MINUS", "PLUS", "MINUS", "BOTH_FACES"}:
            items.append(("FACE", token))
        elif token.startswith("VIRTUAL_DOF_"):
            items.append(("VDOF", token.removeprefix("VIRTUAL_DOF_")))
        elif token.startswith("DOF_"):
            items.append(("DOF", token.removeprefix("DOF_")))
        elif token.startswith("FIELD_"):
            items.append(("FIELD", token.removeprefix("FIELD_")))
        elif token in {"SIGMA_E_ZERO", "GRADIENT_SIGMA_E_ZERO"}:
            # Explicit object-dependent WL schema used by uniform-limit T-0.
            items.append(("FIELD", token))
        elif token.startswith("ORIGIN_"):
            items.append(("ORIGIN", token.removeprefix("ORIGIN_")))
        elif token.startswith("DIRECTION_"):
            items.append(("DIRECTION", token.removeprefix("DIRECTION_")))
        else:
            raise AxisError(f"untyped WL key token {token} in {raw}")
    return make_key(items)


# ---------- closed name/CAS-form map ----------

PARAM = {
    "W0": "W_0",
    "muR": "mu_R",
    "rhoBr": "rho_br",
    "sigmaW": "sigma_W",
    "etaBg": "eta_bg",
    "LW": "L_W",
    "rhoM": "rho_m",
    "lambdaXZero": "Lambda_X_0",
    "lambdaAZero": "Lambda_A_0",
    "lambdaVZero": "Lambda_V_0",
    "tauX": "tau_X",
    "tauA": "tau_A",
    "tauV": "tau_V",
    "omega": "omega",
}

# WL field head -> (PY symbol base, constant scale).  These are the verified
# field folds.  Perturbation-current and branch-mu heads are intentionally not
# in this stripping map.
FIELD = {
    "deltaWidth": ("e_W", "W_0"),
    "virtualDeltaWidth": ("delta_v_e_W", "W_0"),
    "zetaCenter": ("zeta_c", None),
    "virtualZetaCenter": ("delta_v_zeta_c", None),
    "uOne": ("u_1", None),
    "uTwo": ("u_2", None),
    "uThree": ("u_3", None),
    "virtualUOne": ("delta_v_u_1", None),
    "virtualUTwo": ("delta_v_u_2", None),
    "virtualUThree": ("delta_v_u_3", None),
    "thetaWave": ("theta", None),
    "virtualTheta": ("delta_v_theta", None),
    "densityPerturbation": ("delta_rho_4D_bulk", None),
}

FIELD_FACE = {
    "bulkVelocityWaveW": "delta_v_bulk_{face}_4",
    "bulkVelocityWaveOne": "delta_v_bulk_{face}_1",
    "bulkVelocityWaveTwo": "delta_v_bulk_{face}_2",
    "bulkVelocityWaveThree": "delta_v_bulk_{face}_3",
    "pressurePerturbation": "delta_p_{face}",
}

PROFILE = {"w1Jet": "w1_profile", "m1Jet": "m1_profile"}
CURRENT_HEAD = {
    "currentWPerturbation": "delta_j_bulk_4",
    "currentXPerturbation1": "delta_j_bulk_1",
    "currentXPerturbation2": "delta_j_bulk_2",
    "currentXPerturbation3": "delta_j_bulk_3",
}
MU_HEAD_BY_BRANCH = {
    "LAB_HELD": "mu_theta_L",
    "MATERIAL_ADVECTED": "mu_theta_M",
}

Dwin = sp.Function("Dwin")
O_WINDOW = sp.Function("O_window")
SPINT = sp.Function("SPINT")
BOUND_INTEGRAL = sp.Function("BoundIntegral")
BOUND_BINDER = sp.Dummy("bound_integral", dummy_index=9131101)
CAPTURED_FREE_BINDER = sp.Symbol("captured_free_bound_integral")
W_SYMBOL = sp.Symbol("w")
EPSILON = sp.Symbol("epsilon_shape")


def read_bracket(text: str, opening_index: int) -> tuple[str, int]:
    if opening_index >= len(text) or text[opening_index] != "[":
        raise InputError("expected opening [")
    depth = 0
    quote: str | None = None
    escaped = False
    index = opening_index
    while index < len(text):
        character = text[index]
        if quote is not None:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            index += 1
            continue
        if character in {"'", '"'}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[opening_index + 1 : index], index + 1
        index += 1
    raise InputError("unterminated WL bracket expression")


def jet_suffix_from(orders_raw: str, args_raw: str) -> str:
    orders = split_top(orders_raw)
    arguments = split_top(args_raw)
    codes: list[str] = []
    for index, order_raw in enumerate(orders):
        argument = arguments[index] if index < len(arguments) else ""
        if order_raw.startswith("{") and order_raw.endswith("}"):
            suborders = split_top(order_raw[1:-1])
            for subindex, suborder in enumerate(suborders):
                code = ("dw", "t")[subindex] if subindex < 2 else f"q{subindex + 1}"
                codes.extend([code] * int(suborder))
        else:
            order = int(order_raw)
            if order == 0:
                continue
            code = (
                "d1"
                if argument == "x1"
                else "d2"
                if argument == "x2"
                else "d3"
                if argument == "x3"
                else "t"
                if argument == "time"
                else "dw"
                if argument == "normalCoordinate"
                else {0: "d1", 1: "d2", 2: "d3", 3: "t", 4: "t"}.get(index, f"d{index + 1}")
            )
            codes.extend([code] * order)
    time_codes = [code for code in codes if code == "t"]
    other_codes = [code for code in codes if code != "t"]
    return "_".join((*time_codes, *other_codes))


def detect_face(args_raw: str) -> str | None:
    arguments = split_top(args_raw)
    if not arguments:
        return None
    last = arguments[-1].strip()
    if not (last.startswith("{") and last.endswith("}")):
        return None
    pair = split_top(last[1:-1])
    if not pair:
        return None
    point = pair[0].strip()
    return "minus" if point.startswith("-") else "plus"


def preprocess_wl(raw: str) -> str:
    """Apply head-only folds and parser-safe derivative/window/current tokens."""
    text = raw.replace("Inactive[Equal][", "HeldEqual[")
    # The reused Association parser requires textual keys.  These are exactly
    # the two face labels, not numeric data or a positional-axis guess.
    text = re.sub(
        r"(?<![A-Za-z0-9_])(-?1)\s*->",
        lambda match: f'"{"PLUS" if match.group(1) == "1" else "MINUS"}" ->',
        text,
    )
    output: list[str] = []
    index = 0
    while index < len(text):
        if text.startswith("Inactive[Integrate][", index):
            start = index + len("Inactive[Integrate]")
            body, end = read_bracket(text, start)
            parts = split_top(body)
            if len(parts) != 2:
                raise InputError("Inactive[Integrate] does not have integrand and range")
            range_raw = parts[1].strip()
            if not (range_raw.startswith("{") and range_raw.endswith("}")):
                raise InputError("Inactive[Integrate] range is not a list")
            limits = split_top(range_raw[1:-1])
            if len(limits) != 3:
                raise InputError("Inactive[Integrate] range is not (binder,lo,hi)")
            output.append(
                "SPINT["
                + preprocess_wl(parts[0])
                + ","
                + ",".join(preprocess_wl(item) for item in limits)
                + "]"
            )
            index = end
            continue
        if text.startswith("Derivative[", index):
            start = index + len("Derivative")
            orders_raw, after_orders = read_bracket(text, start)
            if after_orders < len(text) and text[after_orders] == "[":
                head_raw, after_head = read_bracket(text, after_orders)
                if after_head < len(text) and text[after_head] == "[":
                    args_raw, end = read_bracket(text, after_head)
                    head = head_raw.strip()
                    orders = split_top(orders_raw)
                    arguments = split_top(args_raw)
                    if head == "windowFunction" and len(orders) == 2 and len(arguments) == 2:
                        output.append(
                            f"Dwin[{int(orders[0])},{int(orders[1])},"
                            f"{preprocess_wl(arguments[0])},{preprocess_wl(arguments[1])}]"
                        )
                        index = end
                        continue
                    suffix = jet_suffix_from(orders_raw, args_raw)
                    if head in CURRENT_HEAD:
                        current_index = CURRENT_HEAD[head].rsplit("_", 1)[1]
                        token = f"CURRENT{current_index}"
                        if suffix:
                            token += "XJETX" + suffix.replace("_", "X")
                        output.append(token + "[" + preprocess_wl(args_raw) + "]")
                        index = end
                        continue
                    if head == "muThetaOperand":
                        token = "MUTHETA"
                        if suffix:
                            token += "XJETX" + suffix.replace("_", "X")
                        output.append(token + "[" + preprocess_wl(args_raw) + "]")
                        index = end
                        continue
                    if head in FIELD_FACE:
                        face = detect_face(args_raw)
                        if face is not None:
                            token = head + "XFACEX" + face
                            if suffix:
                                token += "XJETX" + suffix.replace("_", "X")
                            output.append(token)
                            index = end
                            continue
                    token = head
                    if suffix:
                        token += "XJETX" + suffix.replace("_", "X")
                    output.append(token)
                    index = end
                    continue
        matched_face = None
        for head in FIELD_FACE:
            if text.startswith(head + "[", index):
                matched_face = head
                break
        if matched_face is not None:
            start = index + len(matched_face)
            args_raw, end = read_bracket(text, start)
            face = detect_face(args_raw)
            if face is not None:
                output.append(matched_face + "XFACEX" + face)
                index = end
                continue
        matched_current = None
        for head in CURRENT_HEAD:
            if text.startswith(head + "[", index):
                matched_current = head
                break
        if matched_current is not None:
            start = index + len(matched_current)
            args_raw, end = read_bracket(text, start)
            current_index = CURRENT_HEAD[matched_current].rsplit("_", 1)[1]
            output.append(f"CURRENT{current_index}[{preprocess_wl(args_raw)}]")
            index = end
            continue
        if text.startswith("windowFunction[", index):
            start = index + len("windowFunction")
            args_raw, end = read_bracket(text, start)
            arguments = split_top(args_raw)
            if len(arguments) != 2:
                raise InputError("windowFunction is not the declared two-argument window")
            output.append(
                "Dwin[0,0,"
                + preprocess_wl(arguments[0])
                + ","
                + preprocess_wl(arguments[1])
                + "]"
            )
            index = end
            continue
        output.append(text[index])
        index += 1
    return "".join(output)


def canon_jet_name(name: str) -> str:
    parts = name.split("_")
    base_parts: list[str] = []
    derivatives: list[str] = []
    for part in parts:
        if part in {"t", "dw"} or re.fullmatch(r"(?:d\d+)+", part):
            derivatives.append(part)
        elif derivatives:
            derivatives.append(part)
        else:
            base_parts.append(part)
    if not derivatives:
        return name
    has_time = False
    directions: list[int] = []
    for token in derivatives:
        if token == "t":
            has_time = True
        elif token == "dw":
            directions.append(99)
        else:
            directions.extend(int(item) for item in re.findall(r"d(\d+)", token))
    directions.sort()
    suffix = "_t" if has_time else ""
    if directions:
        suffix += "_" + "".join("dw" if item == 99 else f"d{item}" for item in directions)
    return "_".join(base_parts) + suffix


def pynorm(name: str) -> str:
    return name[4:] + "_dw" if name.startswith("d_w_") else name


def field_symbol(base: str, suffix: str = "", scale: str | None = None) -> sp.Expr:
    name = canon_jet_name(base + ("_" + suffix if suffix else ""))
    value: sp.Expr = sp.Symbol(name)
    if scale is not None:
        value = sp.Symbol(scale) * value
    return value


def canon_wl_basic(value: sp.Basic, branch: str | None) -> sp.Basic:
    def held_equal(node: sp.Basic) -> sp.Basic:
        return sp.Eq(node.args[0], node.args[1], evaluate=False)

    value = value.replace(
        lambda node: isinstance(node, AppliedUndef)
        and node.func.__name__ == "HeldEqual"
        and len(node.args) == 2,
        held_equal,
    )

    def profile(node: AppliedUndef) -> sp.Basic:
        base = PROFILE[node.func.__name__]
        indices = [int(item) for item in node.args]
        suffix = "".join(f"d{index + 1}" * order for index, order in enumerate(indices))
        return sp.Symbol(base + ("_" + suffix if suffix else ""))

    value = value.replace(
        lambda node: isinstance(node, AppliedUndef) and node.func.__name__ in PROFILE,
        profile,
    )

    def mapped_field(node: AppliedUndef) -> sp.Basic:
        base, scale = FIELD[node.func.__name__]
        return field_symbol(base, scale=scale)

    value = value.replace(
        lambda node: isinstance(node, AppliedUndef) and node.func.__name__ in FIELD,
        mapped_field,
    )

    def mapped_current(node: AppliedUndef) -> sp.Basic:
        match = re.fullmatch(r"CURRENT([1-4])(?:XJETX(.*))?", node.func.__name__)
        if match is None:
            return node
        base = f"delta_j_bulk_{match.group(1)}"
        suffix = (match.group(2) or "").replace("X", "_")
        head = canon_jet_name(base + ("_" + suffix if suffix else ""))
        # Argument list, nesting, and arity are intentionally preserved.
        return sp.Function(head)(*node.args)

    value = value.replace(
        lambda node: isinstance(node, AppliedUndef) and node.func.__name__.startswith("CURRENT"),
        mapped_current,
    )

    def mapped_mu(node: AppliedUndef) -> sp.Basic:
        if branch not in MU_HEAD_BY_BRANCH:
            return node
        function_name = node.func.__name__
        suffix = "" if function_name == "muThetaOperand" else function_name.removeprefix("MUTHETA")
        suffix = suffix.removeprefix("XJETX").replace("X", "_")
        head = MU_HEAD_BY_BRANCH[branch]
        if suffix:
            head = canon_jet_name(head + "_" + suffix)
        return sp.Function(head)(*node.args)

    value = value.replace(
        lambda node: isinstance(node, AppliedUndef)
        and (
            node.func.__name__.startswith("MUTHETA")
            or node.func.__name__ == "muThetaOperand"
        ),
        mapped_mu,
    )

    replacements: dict[sp.Symbol, sp.Expr] = {}
    for symbol in value.free_symbols:
        name = symbol.name
        if "XFACEX" in name:
            base, rest = name.split("XFACEX", 1)
            if "XJETX" in rest:
                face, suffix = rest.split("XJETX", 1)
                suffix = suffix.replace("X", "_")
            else:
                face, suffix = rest, ""
            if base in FIELD_FACE:
                target = FIELD_FACE[base].format(face=face)
                replacements[symbol] = sp.Symbol(
                    canon_jet_name(target + ("_" + suffix if suffix else ""))
                )
                continue
        if "XJETX" in name:
            base, suffix = name.split("XJETX", 1)
            suffix = suffix.replace("X", "_")
            if base in FIELD:
                target, scale = FIELD[base]
                replacements[symbol] = field_symbol(target, suffix, scale)
                continue
        if name in PARAM:
            replacements[symbol] = sp.Symbol(PARAM[name])
        elif name in {"waveOrder", "virtualOrder"}:
            replacements[symbol] = sp.Integer(1)
        elif name == "normalCoordinate":
            replacements[symbol] = W_SYMBOL
    if replacements:
        value = value.xreplace(replacements)
    return value


def py_to_dwin(value: sp.Basic) -> sp.Basic:
    def replace_subs(node: sp.Subs) -> sp.Basic:
        inner = node.expr
        if not (isinstance(inner, sp.Derivative) and isinstance(inner.expr, AppliedUndef)):
            return node
        if inner.expr.func.__name__ != "O_window" or len(inner.expr.args) != 2:
            return node
        plus_arg, minus_arg = inner.expr.args
        counts = {symbol: order for symbol, order in inner.variable_count}
        substitutions = dict(zip(node.variables, node.point))
        return Dwin(
            counts.get(plus_arg, 0),
            counts.get(minus_arg, 0),
            substitutions.get(plus_arg, plus_arg),
            substitutions.get(minus_arg, minus_arg),
        )

    value = value.replace(lambda node: isinstance(node, sp.Subs), replace_subs)
    value = value.replace(
        lambda node: isinstance(node, AppliedUndef)
        and node.func.__name__ == "O_window"
        and len(node.args) == 2,
        lambda node: Dwin(0, 0, node.args[0], node.args[1]),
    )
    return value


def spint_to_integral(value: sp.Basic) -> sp.Basic:
    return value.replace(
        lambda node: isinstance(node, AppliedUndef)
        and node.func.__name__ == "SPINT"
        and len(node.args) == 4,
        lambda node: sp.Integral(node.args[0], (node.args[1], node.args[2], node.args[3])),
    )


def canonical_symbols(value: sp.Basic) -> sp.Basic:
    replacements: dict[sp.Basic, sp.Basic] = {}
    for symbol in value.atoms(sp.Symbol):
        if isinstance(symbol, sp.Dummy):
            continue
        name = canon_jet_name(pynorm(symbol.name))
        if name == "Infinity":
            replacements[symbol] = sp.oo
        else:
            replacements[symbol] = W_SYMBOL if name == "w" else sp.Symbol(name)
    return value.xreplace(replacements) if replacements else value


def bound_integral_node(integral: sp.Integral) -> sp.Expr:
    if len(integral.limits) != 1 or len(integral.limits[0]) != 3:
        return integral
    binder, lower, upper = integral.limits[0]
    function = integral.function
    if binder != BOUND_BINDER and BOUND_BINDER in (
        function.free_symbols | lower.free_symbols | upper.free_symbols
    ):
        function = function.xreplace({BOUND_BINDER: CAPTURED_FREE_BINDER})
        lower = lower.xreplace({BOUND_BINDER: CAPTURED_FREE_BINDER})
        upper = upper.xreplace({BOUND_BINDER: CAPTURED_FREE_BINDER})
    substitutions = {binder: BOUND_BINDER}
    function = function.xreplace(substitutions)
    lower = lower.xreplace(substitutions)
    upper = upper.xreplace(substitutions)
    independent: list[sp.Expr] = []
    dependent: list[sp.Expr] = []
    for factor in sp.Mul.make_args(function):
        (dependent if factor.has(BOUND_BINDER) else independent).append(factor)
    outside = sp.Mul(*independent)
    inside = sp.Mul(*dependent)
    return outside * BOUND_INTEGRAL(BOUND_BINDER, lower, upper, inside)


def combine_bound_integrals(value: sp.Basic) -> sp.Basic:
    """Use linearity only for identical canonical bounds."""
    expanded = sp.expand(value)
    grouped: dict[tuple[sp.Expr, sp.Expr], sp.Expr] = defaultdict(lambda: sp.S.Zero)
    remainder: list[sp.Expr] = []
    for term in sp.Add.make_args(expanded):
        integrals = [
            factor
            for factor in sp.Mul.make_args(term)
            if isinstance(factor, AppliedUndef) and factor.func == BOUND_INTEGRAL
        ]
        if len(integrals) != 1:
            remainder.append(term)
            continue
        integral = integrals[0]
        binder, lower, upper, integrand = integral.args
        coefficient = term / integral
        if coefficient.has(binder):
            remainder.append(term)
            continue
        grouped[(lower, upper)] += coefficient * integrand
    combined = list(remainder)
    for (lower, upper), integrand in grouped.items():
        factored = sp.factor_terms(sp.expand(integrand))
        independent: list[sp.Expr] = []
        dependent: list[sp.Expr] = []
        for factor in sp.Mul.make_args(factored):
            (dependent if factor.has(BOUND_BINDER) else independent).append(factor)
        combined.append(
            sp.Mul(*independent)
            * BOUND_INTEGRAL(BOUND_BINDER, lower, upper, sp.Mul(*dependent))
        )
    return sp.Add(*combined)


def canonical_integrals(value: sp.Basic) -> sp.Basic:
    if not isinstance(value, sp.Expr) or isinstance(value, Relational):
        return value
    value = value.replace(lambda node: isinstance(node, sp.Integral), bound_integral_node)
    return combine_bound_integrals(value)


def canonical_basic(value: sp.Basic, engine: str, branch: str | None) -> sp.Basic:
    if engine == "WL":
        value = canon_wl_basic(spint_to_integral(value), branch)
    else:
        value = py_to_dwin(value)
    value = canonical_symbols(value)
    return canonical_integrals(value)


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


def coeff_epsilon(value: object) -> object:
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(lambda item: coeff_epsilon(item))
    if isinstance(value, (tuple, list, sp.Tuple)):
        return tuple(coeff_epsilon(item) for item in value)
    if isinstance(value, Association):
        return Association(tuple((key, coeff_epsilon(item)) for key, item in value.entries))
    if isinstance(value, sp.Expr) and not isinstance(value, (Str, Relational)):
        expanded = sp.expand(value)
        return expanded.coeff(EPSILON, 1) if EPSILON in expanded.free_symbols else expanded
    return value


def finalize_value(value: object, engine: str, key: Key) -> object:
    branch = key_dict(key).get("BRANCH")
    return _convert_parsed_containers(canonical_value(value, engine, branch))


@functools.lru_cache(maxsize=None)
def cached_parse_sympy(raw: str) -> object:
    return parse_sympy_payload(raw)


@functools.lru_cache(maxsize=None)
def cached_parse_wolfram(raw: str) -> object:
    return parse_wolfram_payload(preprocess_wl(raw))


def parse_py_value(raw: str, key: Key, *, coefficient: bool = True) -> object:
    value = cached_parse_sympy(raw)
    if coefficient:
        value = coeff_epsilon(value)
    return finalize_value(value, "PY", key)


def parse_wl_value(raw: str, key: Key) -> object:
    value = cached_parse_wolfram(raw)
    return finalize_value(value, "WL", key)


def add_values(values: Sequence[object]) -> object:
    if not values:
        return sp.S.Zero
    first = values[0]
    if all(isinstance(value, sp.MatrixBase) for value in values):
        return sum(values[1:], first)
    if all(isinstance(value, sp.Expr) for value in values):
        return sp.Add(*values)
    raise InputError(f"partition sum requires algebraic leaves, got {type(first).__name__}")


# ---------- per-family extraction ----------

TRIPLE_OBJECTS = {"FACE_NORMAL", "CONORMAL_DERIV", "FACE_MEASURE_SHAPE_DERIV"}
PROJECTION_OBJECTS = {
    "PROJECTION_SHAPE_DERIV",
    "PROJECTION_STATIC_OPERAND",
    "PROJECTION_DYNAMIC_OPERAND",
    "PROJECTION_RESIDUAL",
}
SIMPLE_OBJECTS = {
    "FACE_VELOCITY",
    "RELATIVE_FLUX",
    "TRACTION",
    "EVOLUTION_MASS_BALANCE",
    "CLOSURE_SHAPE_DERIV",
    "VIRTUAL_CONSTRAINT",
    *PROJECTION_OBJECTS,
}


def text_label(value: object) -> str:
    if isinstance(value, TextAtom):
        return value.value
    if isinstance(value, Str):
        return str(value)
    if isinstance(value, sp.Symbol):
        return value.name
    return str(value)


def parsed_pairs(value: object) -> list[tuple[str, object]]:
    if isinstance(value, Association):
        return list(value.entries)
    if not isinstance(value, (tuple, list, sp.Tuple)):
        raise InputError(f"expected a parsed pair-map, got {type(value).__name__}")
    output: list[tuple[str, object]] = []
    for entry in value:
        if not isinstance(entry, (tuple, list, sp.Tuple)) or len(entry) != 2:
            raise InputError("parsed pair-map entry is not a pair")
        output.append((text_label(entry[0]), entry[1]))
    return output


def build_case(
    engine: str,
    key: Key,
    raw: str,
    parser: Callable[[], object],
    *,
    context: tuple[tuple[str, object], ...] = (),
) -> ParsedCase:
    try:
        return ParsedCase(key, parser(), raw, context)
    except Exception as error:
        return ParsedCase(key, None, raw, context, f"{type(error).__name__}: {error}")


def py_object_cases(object_name: str, raw: str, key: Key) -> list[ParsedCase]:
    """FAMILY_FLATTENER leaf for a direct PY object operand."""

    def parse_selected(selector: Callable[[object], object], coefficient: bool = True) -> object:
        parsed = cached_parse_sympy(raw)
        selected = selector(parsed)
        if coefficient:
            selected = coeff_epsilon(selected)
        return finalize_value(selected, "PY", key)

    if object_name in TRIPLE_OBJECTS:
        return [
            build_case(
                "PY",
                key,
                raw,
                lambda: parse_selected(
                    lambda value: value[2]
                    if isinstance(value, (tuple, list, sp.Tuple)) and len(value) == 3
                    else value
                ),
            )
        ]

    if object_name == "KINEMATIC_BALANCE":
        def kinematic() -> object:
            value = cached_parse_sympy(raw)
            try:
                fields = dict(parsed_pairs(value))
            except InputError:
                return finalize_value(coeff_epsilon(value), "PY", key)
            if not {"OPERAND_A", "OPERAND_B", "RESIDUAL"}.issubset(fields):
                return finalize_value(coeff_epsilon(value), "PY", key)
            selected = tuple(
                coeff_epsilon(fields[name]) for name in ("OPERAND_A", "OPERAND_B", "RESIDUAL")
            )
            return finalize_value(selected, "PY", key)

        return [build_case("PY", key, raw, kinematic)]

    if object_name == "VIRTUAL_WORK_SHAPE_DERIV":
        try:
            parsed = cached_parse_sympy(raw)
            if not isinstance(parsed, (tuple, list, sp.Tuple)):
                return [
                    ParsedCase(
                        key,
                        finalize_value(coeff_epsilon(parsed), "PY", key),
                        raw,
                    )
                ]
            partition_container = parsed[0] if len(parsed) == 2 else parsed
            partitions = dict(parsed_pairs(partition_container))
            if not {"UPPER", "LOWER"}.issubset(partitions):
                # Residual/control leaves can already be the face-summed scalar.
                # Only the actual UPPER/LOWER representation is summed here.
                return [
                    ParsedCase(
                        key,
                        finalize_value(coeff_epsilon(parsed), "PY", key),
                        raw,
                    )
                ]
            upper = coeff_epsilon(partitions["UPPER"])
            lower = coeff_epsilon(partitions["LOWER"])
            value = finalize_value(upper + lower, "PY", key)
            context = (
                ("PY_UPPER", finalize_value(upper, "PY", key)),
                ("PY_LOWER", finalize_value(lower, "PY", key)),
            )
            return [ParsedCase(key, value, raw, context)]
        except Exception as error:
            return [ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}")]

    if object_name == "FACE_SHIFT":
        try:
            value = cached_parse_sympy(raw)
            if not isinstance(value, (tuple, list, sp.Tuple)) or len(value) != 4:
                raise InputError("FACE_SHIFT PY operand is not the declared four-slot aggregate")
            velocity = list(value[1])
            current = list(value[3])
            if len(velocity) != 4 or len(current) != 4:
                raise InputError("FACE_SHIFT velocity/current slot does not have four components")
            leaves = [value[0], *velocity, value[2], *current]
            output: list[ParsedCase] = []
            for field_name, leaf in zip(FACE_SHIFT_FIELDS, leaves):
                leaf_key = key_replace(key, FIELD=field_name)
                selected = finalize_value(coeff_epsilon(leaf), "PY", leaf_key)
                output.append(ParsedCase(leaf_key, selected, raw, (("PY_AGGREGATE_FIELD", field_name),)))
            return output
        except Exception as error:
            return [ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}")]

    if object_name == "EVOLUTION_TERM_ORIGINS":
        try:
            partitions = parsed_pairs(cached_parse_sympy(raw))
            values = [coeff_epsilon(value) for _, value in partitions]
            selected = finalize_value(add_values(values), "PY", key)
            context = tuple((f"PY_ORIGIN_{name}", finalize_value(value, "PY", key)) for name, value in partitions)
            return [ParsedCase(key, selected, raw, context)]
        except Exception as error:
            return [ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}")]

    if object_name == "PROJECTION_TERM_ORIGINS":
        try:
            dynamic_map, static_map = cached_parse_sympy(raw)
            output: list[ParsedCase] = []
            for variant, partition_map in (("DYNAMIC", dynamic_map), ("STATIC", static_map)):
                variant_key = key_replace(key, FIELD=variant)
                partitions = parsed_pairs(partition_map)
                values = [coeff_epsilon(value) for _, value in partitions]
                selected = finalize_value(add_values(values), "PY", variant_key)
                context = tuple(
                    (f"PY_{variant}_ORIGIN_{name}", finalize_value(value, "PY", variant_key))
                    for name, value in partitions
                )
                output.append(ParsedCase(variant_key, selected, raw, context))
            return output
        except Exception as error:
            return [ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}")]

    if object_name == "BACKGROUND_DENSITY_MAP":
        try:
            value = cached_parse_sympy(raw)
            if not isinstance(value, (tuple, list, sp.Tuple)) or len(value) < 3:
                raise InputError("BACKGROUND_DENSITY_MAP PY operand has fewer than three slots")
            output = []
            for field_name, index in (("SIGMA_E_ZERO", 1), ("GRADIENT_SIGMA_E_ZERO", 2)):
                field_key = key_replace(key, FIELD=field_name)
                selected = finalize_value(value[index], "PY", field_key)
                output.append(ParsedCase(field_key, selected, raw))
            return output
        except Exception as error:
            return [ParsedCase(key, None, raw, error=f"{type(error).__name__}: {error}")]

    return [
        build_case(
            "PY",
            key,
            raw,
            lambda: parse_py_value(raw, key, coefficient=True),
        )
    ]


def extract_py_base(family: str, payload: str) -> list[ParsedCase]:
    if family == "ADMISSIBILITY_PREMISE":
        key = make_key(())
        return [
            build_case(
                "PY",
                key,
                payload,
                lambda: parse_py_value(payload, key, coefficient=False),
            )
        ]

    if family in {"FACE_MAP_LAB_HELD", "FACE_MAP_MATERIAL_ADVECTED"}:
        branch = "LAB_HELD" if family.endswith("LAB_HELD") else "MATERIAL_ADVECTED"
        output = []
        for key_raw, value_raw in py_top_pairs(payload):
            key = decode_py_key(key_raw, ("FACE",), (("BRANCH", branch),))
            output.append(
                build_case(
                    "PY",
                    key,
                    value_raw,
                    lambda raw=value_raw, case_key=key: parse_py_value(
                        raw, case_key, coefficient=False
                    ),
                )
            )
        return output

    output: list[ParsedCase] = []
    schema = PY_SCHEMA[family]
    for key_raw, record_raw in py_top_pairs(payload):
        key = decode_py_key(key_raw, schema)
        object_raw = record_raw if family == "BACKGROUND_STATE" else py_field(record_raw)
        if family == "BACKGROUND_STATE":
            output.append(
                build_case(
                    "PY",
                    key,
                    object_raw,
                    lambda raw=object_raw, case_key=key: parse_py_value(
                        raw, case_key, coefficient=False
                    ),
                )
            )
        else:
            output.extend(py_object_cases(family, object_raw, key))
    return output


def wl_plain_case(key: Key, raw: str) -> ParsedCase:
    return build_case("WL", key, raw, lambda: parse_wl_value(raw, key))


def group_wl_origin_cases(
    family: str,
    entries: Sequence[tuple[str, str]],
    *,
    object_axis: bool,
    expression_wrapped: bool,
) -> list[ParsedCase]:
    grouped: dict[tuple[Key, str], list[tuple[str, object, str]]] = defaultdict(list)
    failures: list[ParsedCase] = []
    for key_raw, body_raw in entries:
        key = decode_wl_key(key_raw, object_axis=object_axis)
        origin = key_dict(key).get("ORIGIN", "UNNAMED")
        base_key = key_replace(key, remove=("ORIGIN",))
        try:
            if expression_wrapped:
                expression_raw = wl_field(body_raw, "EXPRESSION")
                if family == "PROJECTION_TERM_ORIGINS":
                    variants = wl_rule_list_fields(expression_raw)
                    leaves = {
                        "DYNAMIC": variants["DYNAMIC_SHAPE_DERIVATIVE"],
                        "STATIC": variants["STATIC_SHAPE_DERIVATIVE"],
                    }
                else:
                    leaves = {"SUM": expression_raw}
            elif family == "PROJECTION_TERM_ORIGINS":
                leaves = {
                    "DYNAMIC": wl_field(body_raw, "DYNAMIC_SHAPE_DERIVATIVE", "EXPRESSION"),
                    "STATIC": wl_field(body_raw, "STATIC_SHAPE_DERIVATIVE", "EXPRESSION"),
                }
            else:
                leaves = {"SUM": wl_field(body_raw, "SHAPE_DERIVATIVE", "EXPRESSION")}
            for variant, leaf_raw in leaves.items():
                variant_key = key_replace(base_key, FIELD=variant) if variant != "SUM" else base_key
                parsed = parse_wl_value(leaf_raw, variant_key)
                grouped[(variant_key, variant)].append((origin, parsed, leaf_raw))
        except Exception as error:
            failures.append(
                ParsedCase(base_key, None, body_raw, error=f"{type(error).__name__}: {error}")
            )
    output = list(failures)
    for (key, variant), partitions in grouped.items():
        try:
            value = finalize_value(
                add_values([item for _, item, _ in partitions]), "WL", key
            )
            context = tuple((f"WL_{variant}_ORIGIN_{name}", item) for name, item, _ in partitions)
            output.append(ParsedCase(key, value, " + ".join(raw for _, _, raw in partitions), context))
        except Exception as error:
            output.append(
                ParsedCase(
                    key,
                    None,
                    " + ".join(raw for _, _, raw in partitions),
                    error=f"{type(error).__name__}: {error}",
                )
            )
    return output


def group_wl_virtual_work(
    entries: Sequence[tuple[str, str]], *, object_axis: bool, both_face_key: bool
) -> list[ParsedCase]:
    grouped: dict[Key, list[tuple[str, object, str]]] = defaultdict(list)
    failures: list[ParsedCase] = []
    for key_raw, body_raw in entries:
        key = decode_wl_key(key_raw, object_axis=object_axis)
        face = key_dict(key).get("FACE")
        base_key = key_replace(key, remove=("FACE",))
        if both_face_key:
            base_key = key_replace(base_key, FACE="BOTH_FACES")
        try:
            leaf_raw = wl_field(body_raw, "EXPRESSION")
            grouped[base_key].append((face or "BOTH_FACES", parse_wl_value(leaf_raw, base_key), leaf_raw))
        except Exception as error:
            failures.append(
                ParsedCase(base_key, None, body_raw, error=f"{type(error).__name__}: {error}")
            )
    output = list(failures)
    for key, faces in grouped.items():
        try:
            value = finalize_value(add_values([item for _, item, _ in faces]), "WL", key)
            context = tuple((f"WL_FACE_{face}", item) for face, item, _ in faces)
            output.append(ParsedCase(key, value, " + ".join(raw for _, _, raw in faces), context))
        except Exception as error:
            output.append(
                ParsedCase(
                    key,
                    None,
                    " + ".join(raw for _, _, raw in faces),
                    error=f"{type(error).__name__}: {error}",
                )
            )
    return output


def extract_wl_base(family: str, payload: str) -> list[ParsedCase]:
    entries = wl_assoc_pairs(payload)

    if family == "ADMISSIBILITY_PREMISE":
        output = []
        for key_raw, body_raw in entries:
            key = decode_wl_key(key_raw)
            leaf_raw = wl_field(body_raw, "EXPRESSION")
            output.append(wl_plain_case(key, leaf_raw))
        return output

    if family in {"FACE_MAP_LAB_HELD", "FACE_MAP_MATERIAL_ADVECTED"}:
        branch = "LAB_HELD" if family.endswith("LAB_HELD") else "MATERIAL_ADVECTED"
        output = []
        for key_raw, body_raw in entries:
            key = make_key((('BRANCH', branch), ('FACE', key_raw)))
            leaf_raw = wl_field(body_raw, "EXPRESSION")

            def face_map(raw: str = leaf_raw, case_key: Key = key) -> object:
                prefix = "Inactive[Equal]"
                if not raw.startswith(prefix + "["):
                    raise InputError("FACE_MAP WL EXPRESSION is not Inactive[Equal]")
                body, end = read_bracket(raw, len(prefix))
                if end != len(raw):
                    raise InputError("FACE_MAP WL Equal has trailing syntax")
                sides = split_top(body)
                if len(sides) != 2:
                    raise InputError("FACE_MAP WL Equal does not have two operands")
                return parse_wl_value(sides[1], case_key)

            try:
                prefix = "Inactive[Equal]"
                body, _ = read_bracket(leaf_raw, len(prefix))
                sides = split_top(body)
                lhs = TextAtom(sides[0])
                output.append(
                    build_case(
                        "WL",
                        key,
                        leaf_raw,
                        face_map,
                        context=(("WL_EQUAL_LHS", _convert_parsed_containers(lhs)),),
                    )
                )
            except Exception as error:
                output.append(
                    ParsedCase(key, None, leaf_raw, error=f"{type(error).__name__}: {error}")
                )
        return output

    if family == "BACKGROUND_DENSITY_MAP":
        output = []
        for key_raw, body_raw in entries:
            density_key = decode_wl_key(key_raw)
            for field_name in ("SIGMA_E_ZERO", "GRADIENT_SIGMA_E_ZERO"):
                key = key_replace(density_key, FIELD=field_name)
                output.append(wl_plain_case(key, wl_field(body_raw, field_name, "EXPRESSION")))
        return output

    if family == "KINEMATIC_BALANCE":
        output = []
        paths = (
            ("OPERAND_A_SHAPE_DERIVATIVE", "EXPRESSION"),
            ("OPERAND_B_SHAPE_DERIVATIVE", "EXPRESSION"),
            ("RESIDUAL_A_MINUS_B", "EXPRESSION"),
        )
        for key_raw, body_raw in entries:
            key = decode_wl_key(key_raw)
            raws = tuple(wl_field(body_raw, *path) for path in paths)

            def kinematic(case_key: Key = key, leaves: tuple[str, ...] = raws) -> object:
                return tuple(parse_wl_value(leaf, case_key) for leaf in leaves)

            output.append(build_case("WL", key, " | ".join(raws), kinematic))
        return output

    if family == "EVOLUTION_TERM_ORIGINS":
        return group_wl_origin_cases(family, entries, object_axis=False, expression_wrapped=False)
    if family == "PROJECTION_TERM_ORIGINS":
        return group_wl_origin_cases(family, entries, object_axis=False, expression_wrapped=False)

    leaf_paths: dict[str, tuple[str, ...]] = {
        "BACKGROUND_STATE": ("EXPRESSION",),
        "FACE_NORMAL": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "CONORMAL_DERIV": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "FACE_MEASURE_SHAPE_DERIV": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "FACE_VELOCITY": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "RELATIVE_FLUX": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "TRACTION": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "VIRTUAL_WORK_SHAPE_DERIV": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "FACE_SHIFT": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "PROJECTION_SHAPE_DERIV": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "PROJECTION_STATIC_OPERAND": ("EXPRESSION",),
        "PROJECTION_DYNAMIC_OPERAND": ("EXPRESSION",),
        "PROJECTION_RESIDUAL": ("EXPRESSION",),
        "VIRTUAL_CONSTRAINT": ("NORMALIZED_VIRTUAL_MASS_VARIATION", "EXPRESSION"),
        "EVOLUTION_MASS_BALANCE": ("SHAPE_DERIVATIVE", "EXPRESSION"),
        "CLOSURE_SHAPE_DERIV": ("SHAPE_DERIVATIVE", "EXPRESSION"),
    }
    path = leaf_paths[family]
    output = []
    for key_raw, body_raw in entries:
        key = decode_wl_key(key_raw)
        output.append(wl_plain_case(key, wl_field(body_raw, *path)))
    return output


def extract_py_nested(family: str, payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    for object_key_raw, cases_raw in py_top_pairs(payload):
        object_name = py_text(object_key_raw)
        schema = PY_SCHEMA.get(object_name)
        if schema is None:
            key = make_key((('OBJECT', object_name),))
            output.append(ParsedCase(key, None, cases_raw, error=f"no PY schema for nested object {object_name}"))
            continue
        for key_raw, value_raw in py_top_pairs(cases_raw):
            key = decode_py_key(key_raw, schema, (("OBJECT", object_name),))
            output.extend(py_object_cases(object_name, value_raw, key))
    return output


def extract_py_control_form(payload: str) -> list[ParsedCase]:
    output: list[ParsedCase] = []
    for key_raw, value_raw in py_top_pairs(payload):
        tokens = py_key_tokens(key_raw)
        object_name = tokens[0] if tokens else "UNPARSED_OBJECT"
        schema = CONTROL_FORM_SCHEMA.get(object_name)
        if schema is None:
            key = make_key((('OBJECT', object_name),))
            output.append(ParsedCase(key, None, value_raw, error=f"no CONTROL_FORM schema for {object_name}"))
            continue
        key = decode_py_key(key_raw, schema)
        output.extend(py_object_cases(object_name, value_raw, key))
    return output


def extract_wl_flat_objects(family: str, payload: str, *, mode: str) -> list[ParsedCase]:
    entries = wl_assoc_pairs(payload)
    origin_entries: list[tuple[str, str]] = []
    vw_entries: list[tuple[str, str]] = []
    output: list[ParsedCase] = []
    for key_raw, body_raw in entries:
        object_name = key_raw.split("|", 1)[0]
        if object_name in {"EVOLUTION_TERM_ORIGINS", "PROJECTION_TERM_ORIGINS"}:
            origin_entries.append((key_raw, body_raw))
            continue
        if object_name == "VIRTUAL_WORK_SHAPE_DERIV" and "|FACE_" in key_raw:
            vw_entries.append((key_raw, body_raw))
            continue
        key = decode_wl_key(key_raw, object_axis=True)
        output.append(wl_plain_case(key, wl_field(body_raw, "EXPRESSION")))
    if origin_entries:
        by_object: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for entry in origin_entries:
            by_object[entry[0].split("|", 1)[0]].append(entry)
        for object_name, object_entries in by_object.items():
            output.extend(
                group_wl_origin_cases(
                    object_name,
                    object_entries,
                    object_axis=True,
                    expression_wrapped=True,
                )
            )
    if vw_entries:
        output.extend(
            group_wl_virtual_work(
                vw_entries,
                object_axis=True,
                both_face_key=(mode == "CONTROL_FORM"),
            )
        )
    return output


def extract_py_dimensions(payload: str) -> list[ParsedCase]:
    output = []
    for object_raw, value_raw in py_top_pairs(payload):
        object_name = py_text(object_raw)
        key = make_key((('OBJECT', object_name),))
        output.append(
            build_case(
                "PY",
                key,
                value_raw,
                lambda raw=value_raw, case_key=key: parse_py_value(raw, case_key, coefficient=False),
            )
        )
    return output


def extract_wl_dimensions(payload: str) -> list[ParsedCase]:
    expression_raw = wl_field(payload, "EXPRESSION")
    output = []
    for object_name, value_raw in wl_assoc_pairs(expression_raw):
        key = make_key((('OBJECT', object_name),))
        output.append(wl_plain_case(key, value_raw))
    return output


def extract_py_homogeneity(payload: str) -> list[ParsedCase]:
    output = []
    for object_raw, value_raw in py_top_pairs(payload):
        object_name = OBJECT_RENAME.get(py_text(object_raw), py_text(object_raw))
        key = make_key((('OBJECT', object_name),))

        def dimensions(raw: str = value_raw, case_key: Key = key) -> object:
            terms = cached_parse_sympy(raw)
            selected = tuple(pair[1] for pair in terms)
            return finalize_value(selected, "PY", case_key)

        output.append(build_case("PY", key, value_raw, dimensions))
    return output


def extract_wl_homogeneity(family: str, payload: str) -> list[ParsedCase]:
    leaf = "RESIDUAL" if family == "HOMOGENEITY_RESIDUAL" else "SOURCE_TERM_DIMENSIONS"
    output = []
    for object_name, body_raw in wl_assoc_pairs(payload):
        key = make_key((('OBJECT', object_name),))
        output.append(wl_plain_case(key, wl_field(body_raw, leaf)))
    return output


def extract_family(family: str, py_payload: str | None, wl_payload: str | None) -> FamilyCases:
    cases = FamilyCases()
    try:
        if py_payload is None:
            cases.extraction_notes.append("PY tag missing")
        elif family == "DIMENSIONS":
            cases.py = extract_py_dimensions(py_payload)
        elif family.startswith("HOMOGENEITY_"):
            cases.py = extract_py_homogeneity(py_payload)
        elif family.startswith(("REP_INVARIANCE_", "CONTROL_INDEPENDENCE_", "UNIFORM_LIMIT_")):
            cases.py = extract_py_nested(family, py_payload)
        elif family.startswith("CONTROL_FORM_"):
            cases.py = extract_py_control_form(py_payload)
        else:
            cases.py = extract_py_base(family, py_payload)
    except Exception as error:
        key = make_key((('OBJECT', "PY_EXTRACTION_ERROR"),))
        cases.py = [ParsedCase(key, None, py_payload or "", error=f"{type(error).__name__}: {error}")]
        cases.extraction_notes.append(f"PY extraction failed: {type(error).__name__}: {error}")

    try:
        if wl_payload is None:
            cases.extraction_notes.append("WL tag missing")
        elif family == "DIMENSIONS":
            cases.wl = extract_wl_dimensions(wl_payload)
        elif family.startswith("HOMOGENEITY_"):
            cases.wl = extract_wl_homogeneity(family, wl_payload)
        elif family.startswith(("REP_INVARIANCE_", "CONTROL_INDEPENDENCE_")):
            cases.wl = extract_wl_flat_objects(family, wl_payload, mode="NESTED")
        elif family.startswith("UNIFORM_LIMIT_"):
            cases.wl = extract_wl_flat_objects(family, wl_payload, mode="UNIFORM")
        elif family.startswith("CONTROL_FORM_"):
            cases.wl = extract_wl_flat_objects(family, wl_payload, mode="CONTROL_FORM")
        else:
            cases.wl = extract_wl_base(family, wl_payload)
    except Exception as error:
        key = make_key((('OBJECT', "WL_EXTRACTION_ERROR"),))
        cases.wl = [ParsedCase(key, None, wl_payload or "", error=f"{type(error).__name__}: {error}")]
        cases.extraction_notes.append(f"WL extraction failed: {type(error).__name__}: {error}")
    return cases


# ---------- unconditional case emission and accounting ----------

def serialise(value: object) -> str:
    class_name = type(value).__name__
    if class_name == "ResidualFailure" and hasattr(value, "reason"):
        # The two operands were printed immediately above; do not repeat
        # multi-megabyte operands inside the typed failure wrapper.
        return f"ResidualFailure(reason={getattr(value, 'reason')!r})"
    if class_name == "BooleanNotResidualable":
        return "BooleanNotResidualable()"
    if class_name == "Mismatch" and hasattr(value, "kind"):
        detail = getattr(value, "detail", None)
        return f"Mismatch(kind={getattr(value, 'kind')!r}, detail={detail!r})"
    if class_name == "UndecidedResidual" and hasattr(value, "reason"):
        return f"UndecidedResidual(reason={getattr(value, 'reason')!r})"
    if class_name == "ResidualAssociation" and hasattr(value, "entries"):
        body = ", ".join(
            f"{key!r}: {serialise(item)}" for key, item in getattr(value, "entries")
        )
        return "ResidualAssociation({" + body + "})"
    if isinstance(value, sp.Basic):
        return sp.srepr(value)
    if isinstance(value, TextAtom):
        return f"TextAtom({value.value!r})"
    if isinstance(value, Association):
        body = ", ".join(f"{key!r}: {serialise(item)}" for key, item in value.entries)
        return "Association({" + body + "})"
    if isinstance(value, tuple):
        body = ", ".join(serialise(item) for item in value)
        if len(value) == 1:
            body += ","
        return "(" + body + ")"
    return repr(value)


def serialise_key(key: Key) -> str:
    return "(" + ", ".join(f"{axis}={value}" for axis, value in key) + ")"


def display_operand(case: ParsedCase | None) -> str:
    if case is None:
        return "<MISSING>"
    if case.error is not None:
        return f"<PARSE_FAILED {case.error}; RAW={case.raw}>"
    return serialise(case.value)


def emit_case(
    family: str,
    key: Key,
    py_case: ParsedCase | None,
    wl_case: ParsedCase | None,
    difference: object,
    *,
    note: str | None = None,
) -> None:
    print(f"CASE family={family} key={serialise_key(key)}", flush=True)
    print(f"operand_A = {display_operand(py_case)}", flush=True)
    print(f"operand_B = {display_operand(wl_case)}", flush=True)
    print(f"A_minus_B = {serialise(difference)}", flush=True)
    if note is not None:
        print(f"case_note = {note}", flush=True)
    for label, value in py_case.context if py_case is not None else ():
        print(f"context.{label} = {serialise(value)}", flush=True)
    for label, value in wl_case.context if wl_case is not None else ():
        print(f"context.{label} = {serialise(value)}", flush=True)


def exact_rational_residual(py_value: sp.Expr, wl_value: sp.Expr) -> sp.Expr:
    """The closed CAS-form route: expand first, then cancel(together(.)) exactly."""
    expanded = sp.expand(py_value - wl_value)
    return sp.cancel(sp.together(expanded))


def replace_budget_failures(
    value: object,
    py_value: object,
    wl_value: object,
) -> object:
    """Finish S11b leaves whose optional final factorization exceeds its budget."""
    if (
        type(value).__name__ == "ResidualFailure"
        and getattr(value, "reason", None) == "RESIDUAL_BUDGET_EXCEEDED"
        and isinstance(py_value, sp.Expr)
        and isinstance(wl_value, sp.Expr)
    ):
        return exact_rational_residual(py_value, wl_value)

    if type(value).__name__ == "ResidualAssociation" and isinstance(
        py_value, Association
    ) and isinstance(wl_value, Association):
        py_items = py_value.as_dict()
        wl_items = wl_value.as_dict()
        repaired = tuple(
            (
                key,
                replace_budget_failures(item, py_items[key], wl_items[key]),
            )
            for key, item in value.entries
        )
        return type(value)(repaired)

    if isinstance(value, tuple):
        if isinstance(py_value, Relational) and isinstance(wl_value, Relational):
            operands = (
                (py_value.lhs, wl_value.lhs),
                (py_value.rhs, wl_value.rhs),
            )
        elif isinstance(py_value, tuple) and isinstance(wl_value, tuple):
            operands = tuple(zip(py_value, wl_value))
        else:
            return value
        if len(value) != len(operands):
            return value
        return tuple(
            replace_budget_failures(item, left, right)
            for item, (left, right) in zip(value, operands)
        )
    return value


@functools.lru_cache(maxsize=None)
def cached_typed_residual(
    py_value: object,
    wl_value: object,
    family: str,
    leaf_budget_seconds: float,
) -> object:
    value = residual(
        py_value,
        wl_value,
        name=family,
        leaf_budget_seconds=leaf_budget_seconds,
    )
    return replace_budget_failures(value, py_value, wl_value)


def axis_mismatch_detail(py_key: Key, wl_key: Key) -> str | None:
    py_axes = key_dict(py_key)
    wl_axes = key_dict(wl_key)
    common = set(py_axes) & set(wl_axes)
    if not common and (py_axes and wl_axes):
        return None
    face_granularity = False
    for axis in common:
        if py_axes[axis] == wl_axes[axis]:
            continue
        if axis == "FACE" and "BOTH_FACES" in {py_axes[axis], wl_axes[axis]}:
            face_granularity = True
            continue
        return None
    py_only_axes = sorted(set(py_axes) - set(wl_axes))
    wl_only_axes = sorted(set(wl_axes) - set(py_axes))
    if not py_only_axes and not wl_only_axes and not face_granularity:
        return None
    details = []
    if py_only_axes:
        details.append("WL missing " + ",".join(py_only_axes))
    if wl_only_axes:
        details.append("PY missing " + ",".join(wl_only_axes))
    if face_granularity:
        details.append(
            f"FACE granularity {py_axes.get('FACE')} vs {wl_axes.get('FACE')}"
        )
    return "; ".join(details)


def compare_family(
    family: str,
    cases: FamilyCases,
    *,
    leaf_budget_seconds: float,
) -> Accounting:
    accounting = Accounting()
    py_by_key: dict[Key, list[ParsedCase]] = defaultdict(list)
    wl_by_key: dict[Key, list[ParsedCase]] = defaultdict(list)
    for case in cases.py:
        py_by_key[case.key].append(case)
        if case.error is not None:
            accounting.parse_failed += 1
    for case in cases.wl:
        wl_by_key[case.key].append(case)
        if case.error is not None:
            accounting.parse_failed += 1

    duplicate_keys = {
        key for key, values in py_by_key.items() if len(values) > 1
    } | {key for key, values in wl_by_key.items() if len(values) > 1}
    accounting.duplicate_key = sum(
        len(py_by_key.get(key, ())) + len(wl_by_key.get(key, ())) for key in duplicate_keys
    )

    common_keys = (set(py_by_key) & set(wl_by_key)) - duplicate_keys
    for key in sorted(common_keys, key=serialise_key):
        py_case = py_by_key[key][0]
        wl_case = wl_by_key[key][0]
        accounting.join += 1
        if py_case.error is not None or wl_case.error is not None:
            difference: object = TextAtom("UNDEFINED_PARSE_FAILED")
        else:
            try:
                difference = cached_typed_residual(
                    py_case.value,
                    wl_case.value,
                    family,
                    leaf_budget_seconds,
                )
            except TypeError:
                difference = residual(
                    py_case.value,
                    wl_case.value,
                    name=family,
                    leaf_budget_seconds=leaf_budget_seconds,
                )
        emit_case(family, key, py_case, wl_case, difference)

    for key in sorted(duplicate_keys, key=serialise_key):
        py_values = py_by_key.get(key, ())
        wl_values = wl_by_key.get(key, ())
        width = max(len(py_values), len(wl_values))
        for index in range(width):
            emit_case(
                family,
                key,
                py_values[index] if index < len(py_values) else None,
                wl_values[index] if index < len(wl_values) else None,
                TextAtom("UNDEFINED_DUPLICATE_KEY"),
                note=f"duplicate_key occurrence={index + 1}/{width}",
            )

    py_unmatched_keys = set(py_by_key) - set(wl_by_key) - duplicate_keys
    wl_unmatched_keys = set(wl_by_key) - set(py_by_key) - duplicate_keys
    accounting.py_only = sum(len(py_by_key[key]) for key in py_unmatched_keys)
    accounting.wl_only = sum(len(wl_by_key[key]) for key in wl_unmatched_keys)

    mismatched_py: set[Key] = set()
    mismatched_wl: set[Key] = set()
    mismatch_reasons_by_key: dict[tuple[str, Key], set[str]] = defaultdict(set)
    for py_key in py_unmatched_keys:
        for wl_key in wl_unmatched_keys:
            detail = axis_mismatch_detail(py_key, wl_key)
            if detail is None:
                continue
            mismatched_py.add(py_key)
            mismatched_wl.add(wl_key)
            mismatch_reasons_by_key[("PY", py_key)].add(detail)
            mismatch_reasons_by_key[("WL", wl_key)].add(detail)
            accounting.reasons.add(detail)
    accounting.axis_set_mismatch = sum(
        len(py_by_key[key]) for key in mismatched_py
    ) + sum(len(wl_by_key[key]) for key in mismatched_wl)

    for key in sorted(py_unmatched_keys, key=serialise_key):
        reasons = sorted(mismatch_reasons_by_key.get(("PY", key), ()))
        note = "axis_set_mismatch: " + " | ".join(reasons) if reasons else "py_only"
        for case in py_by_key[key]:
            emit_case(
                family,
                key,
                case,
                None,
                TextAtom("UNDEFINED_UNJOINED"),
                note=note,
            )
    for key in sorted(wl_unmatched_keys, key=serialise_key):
        reasons = sorted(mismatch_reasons_by_key.get(("WL", key), ()))
        note = "axis_set_mismatch: " + " | ".join(reasons) if reasons else "wl_only"
        for case in wl_by_key[key]:
            emit_case(
                family,
                key,
                None,
                case,
                TextAtom("UNDEFINED_UNJOINED"),
                note=note,
            )

    accounting.reasons.update(cases.extraction_notes)
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
    inventory = wl_tags.get("LOCAL_TAG_NAMES")
    if inventory is not None:
        print(f"LOCAL_INVENTORY engine=WL payload={inventory}", flush=True)


def run(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--py", type=Path, default=DEFAULT_PY, help="committed PY .out tag stream")
    parser.add_argument("--wl", type=Path, default=DEFAULT_WL, help="committed WL .out tag stream")
    parser.add_argument(
        "--residual-leaf-budget",
        type=float,
        default=0.1,
        help=(
            "seconds allowed to S11b's factor/cancel residual leaf before the exact "
            "expand -> cancel(together(.)) fallback"
        ),
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
    results: dict[str, Accounting] = {}
    for family in FAMILIES:
        cases = extract_family(family, py_tags.get(family), wl_tags.get(family))
        results[family] = compare_family(
            family,
            cases,
            leaf_budget_seconds=arguments.residual_leaf_budget,
        )

    elapsed = time.monotonic() - started
    print(
        f"RUN_ACCOUNTING families={len(FAMILIES)} "
        f"families_with_join={sum(item.join > 0 for item in results.values())} "
        f"families_with_unpaired={sum((item.py_only + item.wl_only) > 0 for item in results.values())} "
        f"parse_failed={sum(item.parse_failed for item in results.values())} "
        f"duplicate_key={sum(item.duplicate_key for item in results.values())} "
        f"runtime_seconds={elapsed:.3f}",
        flush=True,
    )
    print(
        "MEASUREMENT_SCOPE supplied_unfalsifiable=sections_1_to_3,admissibility_premise,"
        "supplied_bookkeeping residual_target=none",
        flush=True,
    )
    # Algebraic and structural disagreements are measured results, so the run
    # returns success after emitting them.
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
