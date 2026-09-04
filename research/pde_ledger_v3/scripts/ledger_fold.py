"""Fold F9-resolved LEDGER exports and guard consumer import manifests.

This module is physics-free infrastructure.  It treats row payloads as opaque
data, never imports an engine export on its own, and never applies F9 routing.
Callers must digest-pin this shared executable input in ``BUILD_INPUT_DIGESTS``.
"""

from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Callable, Iterable, Iterator, Mapping
import hashlib
import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import sympy as sp
from sympy.core.function import AppliedUndef, UndefinedFunction
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair


_RELATIONALS = {
    "Equality": lambda left, right: sp.Eq(left, right, evaluate=False),
    "Unequality": lambda left, right: sp.Ne(left, right, evaluate=False),
    "StrictGreaterThan": lambda left, right: sp.Gt(left, right, evaluate=False),
    "StrictLessThan": lambda left, right: sp.Lt(left, right, evaluate=False),
    "GreaterThan": lambda left, right: sp.Ge(left, right, evaluate=False),
    "LessThan": lambda left, right: sp.Le(left, right, evaluate=False),
}
_MISSING = object()


class LedgerFoldError(RuntimeError):
    """Base class for a fold or under-export contract failure."""


class ManifestError(LedgerFoldError):
    """A consumer manifest is missing, malformed, or unresolved."""


class ClosureError(LedgerFoldError):
    """A recursive symbol or dimension edge cannot be resolved."""


class AmbiguousSymbolError(ClosureError):
    """A serialized symbol edge has more than one possible write-key."""


def _restore(source: str) -> object:
    """Restore the srepr form used by generated ``*_exports.py`` modules."""

    if not isinstance(source, str):
        raise TypeError(f"srepr must be str, got {type(source).__name__}")
    return eval(  # noqa: S307 - deliberately matches the existing export format
        source,
        {
            "__builtins__": {},
            "Str": Str,
            "ExprCondPair": ExprCondPair,
            **vars(sp),
            **_RELATIONALS,
        },
    )


def _import_path(path: str | Path) -> ModuleType:
    source_path = Path(path)
    if not source_path.is_file():
        raise LedgerFoldError(f"export module does not exist: {source_path}")
    identity = hashlib.sha256(str(source_path.resolve()).encode("utf-8")).hexdigest()
    module_name = f"_ledger_fold_export_{identity}"
    spec = importlib.util.spec_from_file_location(module_name, source_path)
    if spec is None or spec.loader is None:
        raise LedgerFoldError(f"cannot import export module: {source_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _ledger_from_path(path: str | Path) -> Mapping[str, Mapping[str, object]]:
    module = _import_path(path)
    ledger = getattr(module, "LEDGER", _MISSING)
    if ledger is _MISSING:
        raise LedgerFoldError(f"export module has no LEDGER: {Path(path)}")
    if not isinstance(ledger, Mapping):
        raise LedgerFoldError(f"LEDGER is not a mapping: {Path(path)}")
    malformed = sorted(
        repr(key)
        for key, row in ledger.items()
        if not isinstance(key, str) or not isinstance(row, Mapping)
    )
    if malformed:
        raise LedgerFoldError(
            f"LEDGER has non-string keys or non-mapping rows in {Path(path)}: "
            + ", ".join(malformed)
        )
    return ledger


def load_model(
    base_path: str | Path, *delta_paths: str | Path
) -> tuple[dict[str, Mapping[str, object]], dict[str, list[object]]]:
    """Import and fold a base plus chronological deltas by exact key.

    F9 has already been applied by each writer.  This function neither compares
    row values nor changes key spellings.  An exact-key delta replacement is
    retained as an F9b overwrite audit tuple ``(key, prior_step, new_step)``.
    """

    sources = (base_path, *delta_paths)
    merged: dict[str, Mapping[str, object]] = {}
    source_row_counts: list[object] = []
    overwrites: list[object] = []

    for position, source in enumerate(sources):
        ledger = _ledger_from_path(source)
        source_row_counts.append((str(Path(source)), len(ledger)))
        for key, row in ledger.items():
            if position and key in merged:
                overwrites.append((key, merged[key].get("step"), row.get("step")))
            merged[key] = row

    audit = {
        "overwrites": overwrites,
        "source_row_counts": source_row_counts,
    }
    return merged, audit


def _manifest_keys(import_keys: Iterable[str] | None) -> frozenset[str]:
    if import_keys is None:
        raise ManifestError("IMPORT_KEYS manifest is missing")
    if isinstance(import_keys, str):
        raise ManifestError(f"IMPORT_KEYS must be an iterable of keys, got {import_keys!r}")
    try:
        keys = frozenset(import_keys)
    except TypeError as exc:
        raise ManifestError(f"IMPORT_KEYS is not a valid key iterable: {exc}") from exc
    bad = sorted(repr(key) for key in keys if not isinstance(key, str))
    if bad:
        raise ManifestError("IMPORT_KEYS contains non-string keys: " + ", ".join(bad))
    return keys


def _referenced_atoms(value: object) -> set[object]:
    """Return typed reference atoms from SymPy values and ordinary containers."""

    if isinstance(value, Mapping):
        atoms: set[object] = set()
        for key, item in value.items():
            atoms.update(_referenced_atoms(key))
            atoms.update(_referenced_atoms(item))
        return atoms
    if isinstance(value, (tuple, list, set, frozenset)):
        atoms = set()
        for item in value:
            atoms.update(_referenced_atoms(item))
        return atoms
    if isinstance(value, UndefinedFunction):
        return {value}
    if isinstance(value, sp.Basic):
        atoms = set(value.atoms(sp.Symbol))
        atoms.update(function.func for function in value.atoms(AppliedUndef))
        atoms.update(value.atoms(sp.MatrixSymbol))
        return atoms
    return set()


def _atom_identity(atom: object) -> str:
    """Serialize a reference atom's full identity, including assumptions."""

    return sp.srepr(atom)


def _producer_identities(value: object) -> set[str]:
    """Return only identities declared by a row's own value, never identities it uses."""

    if isinstance(value, (sp.Symbol, sp.MatrixSymbol, UndefinedFunction)):
        return {_atom_identity(value)}
    return set()


def _identity_rows(
    fold: Mapping[str, Mapping[str, object]],
) -> Mapping[str, frozenset[str]]:
    """Index declaration-valued rows by full atom identity."""

    index: dict[str, set[str]] = defaultdict(set)
    for key, row in fold.items():
        if not isinstance(row, Mapping):
            continue
        for identity in _producer_identities(row.get("value", _MISSING)):
            index[identity].add(key)
    return {identity: frozenset(keys) for identity, keys in index.items()}


def check_consumer(
    fold: Mapping[str, Mapping[str, object]], import_keys: Iterable[str] | None
) -> dict[str, object]:
    """Resolve an exact manifest and its identity-resolvable recursive closure.

    Each typed reference atom resolves by full ``srepr`` identity to rows whose
    own value declares that atom.  A unique declaration is followed, no
    declaration denotes a free/structural atom and is skipped, and duplicate
    declarations are a genuine ambiguity.
    """

    if not isinstance(fold, Mapping):
        raise TypeError(f"fold must be a mapping, got {type(fold).__name__}")
    manifest = _manifest_keys(import_keys)
    missing_roots = sorted(manifest.difference(fold))
    if missing_roots:
        raise ManifestError(
            "IMPORT_KEYS missing exact write-key(s): " + ", ".join(missing_roots)
        )

    identity_rows = _identity_rows(fold)
    closure = set(manifest)
    pending = deque(sorted(manifest))
    symbol_edges: list[tuple[str, str, str]] = []
    dimension_edges: list[tuple[str, str]] = []

    while pending:
        source_key = pending.popleft()
        row = fold[source_key]
        if not isinstance(row, Mapping):
            raise ClosureError(f"closure row {source_key!r} is not a mapping")
        if "value" not in row:
            raise ClosureError(f"closure row {source_key!r} has no value")

        value = row["value"]
        referenced_atoms = _referenced_atoms(value)
        if _producer_identities(value):
            referenced_atoms.discard(value)

        for atom in sorted(referenced_atoms, key=_atom_identity):
            identity = _atom_identity(atom)
            candidates = identity_rows.get(identity, ())
            if not candidates:
                continue
            if len(candidates) > 1:
                competing = ", ".join(sorted(candidates))
                raise AmbiguousSymbolError(
                    f"ambiguous atom identity {identity!r} from {source_key!r}; "
                    f"competing write-keys: {competing}"
                )
            target = next(iter(candidates))
            symbol_edges.append((source_key, identity, target))
            if target not in closure:
                closure.add(target)
                pending.append(target)

        if "dimension_key" in row:
            target = row["dimension_key"]
            if not isinstance(target, str) or target not in fold:
                raise ClosureError(
                    f"dimension edge from {source_key!r} is missing exact key {target!r}"
                )
            dimension_edges.append((source_key, target))
            if target not in closure:
                closure.add(target)
                pending.append(target)

    return {
        "import_keys": manifest,
        "resolved_imports": {key: fold[key] for key in sorted(manifest)},
        "closure": frozenset(closure),
        "symbol_edges": tuple(symbol_edges),
        "dimension_edges": tuple(dimension_edges),
    }


class _AccessRecordingLedger(Mapping[str, Mapping[str, object]]):
    """Record only mapping lookups made through this proxy.

    A consumer that captures the raw fold or reaches ``_fold`` evades the
    witness.  Such access is outside the consumer contract; callers must use
    only the proxy argument supplied to their build function.
    """

    def __init__(self, fold: Mapping[str, Mapping[str, object]]) -> None:
        self._fold = fold
        self.lookups: set[str] = set()

    def __getitem__(self, key: str) -> Mapping[str, object]:
        self.lookups.add(key)
        return self._fold[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._fold)

    def __len__(self) -> int:
        return len(self._fold)

    def __contains__(self, key: object) -> bool:
        return key in self._fold


def assert_lookups_equal_manifest(
    build_fn: Callable[[Mapping[str, Mapping[str, object]]], Any],
    fold: Mapping[str, Mapping[str, object]],
    import_keys: Iterable[str] | None,
) -> dict[str, object]:
    """Run a consumer and require its recorded ``__getitem__`` set to be exact."""

    manifest = _manifest_keys(import_keys)
    proxy = _AccessRecordingLedger(fold)
    result = build_fn(proxy)
    undeclared = sorted(proxy.lookups.difference(manifest))
    unused = sorted(manifest.difference(proxy.lookups))
    if undeclared or unused:
        parts = []
        if undeclared:
            parts.append("undeclared lookup(s): " + ", ".join(undeclared))
        if unused:
            parts.append("declared-but-unused key(s): " + ", ".join(unused))
        raise ManifestError("lookup/IMPORT_KEYS mismatch; " + "; ".join(parts))
    return {"lookups": frozenset(proxy.lookups), "result": result}


def assert_delta_is_minimal(
    delta_ledger: Mapping[str, Mapping[str, object]],
    own_bind_closure: Iterable[str],
    infra_keys: Iterable[str] = (),
) -> dict[str, frozenset[str]]:
    """Require a delta to contain exactly its own closure plus named infra."""

    actual = frozenset(delta_ledger)
    closure = _manifest_keys(own_bind_closure)
    infrastructure = _manifest_keys(infra_keys)
    allowed = closure | infrastructure
    extra = sorted(actual.difference(allowed))
    missing = sorted(closure.difference(actual))
    if extra or missing:
        parts = []
        if extra:
            parts.append("extra row(s): " + ", ".join(extra))
        if missing:
            parts.append("missing row(s): " + ", ".join(missing))
        raise LedgerFoldError("delta is not minimal; " + "; ".join(parts))
    return {
        "exported_keys": actual,
        "required_keys": closure,
        "allowed_infra_keys": infrastructure,
    }


def promotion_delta(
    row_key: str,
    srepr: str,
    cls: str,
    evidence: Mapping[str, object],
    *,
    display: str | None = None,
    value_kind: str | None = None,
    step: str | None = None,
    dimension_key: str | object = _MISSING,
    f9_operands: object = _MISSING,
    route: object = _MISSING,
    corroborated_steps: Iterable[str] | object = _MISSING,
    **metadata: object,
) -> dict[str, dict[str, object]]:
    """Construct a schema-complete one-row promotion delta.

    ``evidence`` must record ``step``, ``f9_operands``, and the writer-decided
    F9 ``route`` (explicit keyword arguments may supply or override them).  The
    returned mapping is a real delta artifact, not a manifest entry.  A Python
    consumer must never re-derive an inherited object at consume time (N1).
    """

    if not isinstance(row_key, str) or not row_key:
        raise ValueError(f"row_key must be a non-empty string, got {row_key!r}")
    if not isinstance(cls, str) or not cls:
        raise ValueError(f"cls must be a non-empty string, got {cls!r}")
    if not isinstance(evidence, Mapping):
        raise TypeError(f"evidence must be a mapping, got {type(evidence).__name__}")
    forbidden = sorted({"value", "class"}.intersection(evidence))
    if forbidden:
        raise LedgerFoldError(
            "promotion evidence cannot replace reconstructed field(s): "
            + ", ".join(forbidden)
        )

    value = _restore(srepr)
    row = dict(evidence)
    row.update(metadata)
    row["display"] = display if display is not None else row.get("display", str(value))
    row["value"] = value
    row["value_kind"] = (
        value_kind if value_kind is not None else row.get("value_kind", "COMPUTED_OBJECT")
    )
    row["class"] = cls
    if step is not None:
        row["step"] = step
    if dimension_key is not _MISSING:
        row["dimension_key"] = dimension_key
    if f9_operands is not _MISSING:
        row["f9_operands"] = f9_operands
    if route is not _MISSING:
        row["route"] = route
    if corroborated_steps is not _MISSING:
        row["corroborated_steps"] = tuple(corroborated_steps)

    missing_evidence = [
        field for field in ("step", "f9_operands", "route") if field not in row
    ]
    if missing_evidence:
        raise LedgerFoldError(
            f"promotion row {row_key!r} is missing evidence field(s): "
            + ", ".join(missing_evidence)
        )
    return {row_key: row}


__all__ = [
    "AmbiguousSymbolError",
    "ClosureError",
    "LedgerFoldError",
    "ManifestError",
    "assert_delta_is_minimal",
    "assert_lookups_equal_manifest",
    "check_consumer",
    "load_model",
    "promotion_delta",
]
