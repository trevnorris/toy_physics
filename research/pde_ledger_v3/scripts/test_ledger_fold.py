#!/usr/bin/env python3
"""Decisive self-test for ledger_fold.py; imports no real engine export."""

from __future__ import annotations

from pathlib import Path
import tempfile
import time
import traceback

import sympy as sp

from ledger_fold import (
    AmbiguousSymbolError,
    ClosureError,
    LedgerFoldError,
    ManifestError,
    assert_delta_is_minimal,
    assert_lookups_equal_manifest,
    check_consumer,
    load_model,
    promotion_delta,
)


def row(value: object, step: str, **fields: object) -> dict[str, object]:
    record = {
        "display": str(value),
        "value": value,
        "value_kind": "COMPUTED_OBJECT",
        "class": "SYNTHETIC",
        "step": step,
    }
    record.update(fields)
    return record


def emit(name: str, operands: object, observed: object) -> None:
    print(f"TEST {name} OPERANDS {operands!r}")
    print(f"TEST {name} OBSERVED {observed!r}")


def write_synthetic_exports(directory: Path) -> tuple[Path, Path]:
    base_path = directory / "base_exports.py"
    delta_path = directory / "delta_exports.py"
    base_path.write_text(
        """from types import MappingProxyType
import sympy as sp
LEDGER = MappingProxyType({
    'shared': MappingProxyType({'display': '1', 'value': sp.Integer(1), 'value_kind': 'COMPUTED_OBJECT', 'class': 'SYNTHETIC', 'step': 'BASE'}),
    'x': MappingProxyType({'display': '2', 'value': sp.Integer(2), 'value_kind': 'COMPUTED_OBJECT', 'class': 'SYNTHETIC', 'step': 'BASE'}),
    'face_response': MappingProxyType({'display': 'flat', 'value': sp.Integer(11), 'value_kind': 'COMPUTED_OBJECT', 'class': 'SYNTHETIC', 'step': 'S11b'}),
})
""",
        encoding="utf-8",
    )
    delta_path.write_text(
        """from types import MappingProxyType
import sympy as sp
LEDGER = MappingProxyType({
    'shared': MappingProxyType({'display': '3', 'value': sp.Integer(3), 'value_kind': 'COMPUTED_OBJECT', 'class': 'SYNTHETIC', 'step': 'DELTA'}),
    'added': MappingProxyType({'display': '4', 'value': sp.Integer(4), 'value_kind': 'COMPUTED_OBJECT', 'class': 'SYNTHETIC', 'step': 'DELTA'}),
    's11c_c1_x': MappingProxyType({'display': '5', 'value': sp.Integer(5), 'value_kind': 'COMPUTED_OBJECT', 'class': 'SYNTHETIC', 'step': 'S11c-c1'}),
    's11c_c1_face_response': MappingProxyType({'display': 'curved', 'value': sp.Integer(12), 'value_kind': 'COMPUTED_OBJECT', 'class': 'SYNTHETIC', 'step': 'S11c-c1'}),
})
""",
        encoding="utf-8",
    )
    return base_path, delta_path


def test_last_wins(fold: dict, audit: dict, base_path: Path, delta_path: Path) -> None:
    operands = {"base": str(base_path), "deltas": [str(delta_path)], "key": "shared"}
    observed = {
        "value": fold["shared"]["value"],
        "keys": sorted(fold),
        "audit": audit,
    }
    emit("last_wins", operands, observed)
    assert fold["shared"]["value"] == 3
    assert "added" in fold
    assert audit["overwrites"] == [("shared", "BASE", "DELTA")]
    assert audit["source_row_counts"] == [(str(base_path), 3), (str(delta_path), 4)]


def test_no_f9_reapply(fold: dict) -> None:
    operands = {"already_resolved_keys": ["x", "s11c_c1_x"]}
    observed = {key: fold[key]["value"] for key in operands["already_resolved_keys"]}
    emit("no_f9_reapply", operands, observed)
    assert "x" in fold and "s11c_c1_x" in fold
    assert fold["x"]["value"] == 2
    assert fold["s11c_c1_x"]["value"] == 5


def test_predecessor_trap(fold: dict) -> None:
    flat_report = check_consumer(fold, {"face_response"})
    curved_report = check_consumer(fold, {"s11c_c1_face_response"})
    only_flat = dict(fold)
    del only_flat["s11c_c1_face_response"]
    try:
        check_consumer(only_flat, {"s11c_c1_face_response"})
    except ManifestError as exc:
        exact_key_error = str(exc)
    else:
        exact_key_error = "NO ERROR"
    flat = flat_report["resolved_imports"]["face_response"]
    curved = curved_report["resolved_imports"]["s11c_c1_face_response"]
    operands = {"flat_manifest": ["face_response"], "curved_manifest": ["s11c_c1_face_response"]}
    observed = {
        "flat_step": flat["step"],
        "curved_step": curved["step"],
        "flat_value": flat["value"],
        "curved_value": curved["value"],
        "prefixed_missing_with_bare_present": exact_key_error,
    }
    emit("predecessor_trap", operands, observed)
    assert flat is fold["face_response"]
    assert curved is fold["s11c_c1_face_response"]
    assert flat is not curved and flat["value"] != curved["value"]
    assert flat["step"] == "S11b" and curved["step"] == "S11c-c1"
    assert "s11c_c1_face_response" in exact_key_error
    assert exact_key_error != "NO ERROR"


def test_missing_symbol_closure() -> None:
    missing_fold = {"kept": row(sp.Symbol("foo") + 1, "D")}
    try:
        check_consumer(missing_fold, {"kept"})
    except ClosureError as exc:
        missing_error = str(exc)
    else:
        missing_error = "NO ERROR"
    complete_fold = dict(missing_fold)
    complete_fold["foo"] = row(sp.Symbol("foo"), "B")
    report = check_consumer(complete_fold, {"kept"})
    operands = {"manifest": ["kept"], "missing_keys": sorted(missing_fold), "complete_keys": sorted(complete_fold)}
    observed = {"missing_error": missing_error, "closure_after_add": sorted(report["closure"])}
    emit("missing_symbol_closure", operands, observed)
    assert "foo" in missing_error and missing_error != "NO ERROR"
    assert report["closure"] == frozenset({"kept", "foo"})


def test_dimension_closure() -> None:
    missing_fold = {"kept": row(sp.Integer(1), "D", dimension_key="dim_kept")}
    try:
        check_consumer(missing_fold, {"kept"})
    except ClosureError as exc:
        missing_error = str(exc)
    else:
        missing_error = "NO ERROR"
    complete_fold = dict(missing_fold)
    complete_fold["dim_kept"] = row(sp.Integer(2), "B")
    report = check_consumer(complete_fold, {"kept"})
    operands = {"manifest": ["kept"], "dimension_key": "dim_kept"}
    observed = {"missing_error": missing_error, "closure_after_add": sorted(report["closure"])}
    emit("dimension_closure", operands, observed)
    assert "dim_kept" in missing_error and missing_error != "NO ERROR"
    assert report["closure"] == frozenset({"kept", "dim_kept"})


def test_f9c_ambiguity() -> None:
    ambiguous_fold = {
        "kept": row(sp.Symbol("a1") + 1, "D"),
        "a1": row(sp.Symbol("a1"), "BASE"),
        "s11_a1": row(sp.Symbol("a1"), "DELTA"),
    }
    try:
        check_consumer(ambiguous_fold, {"kept"})
    except AmbiguousSymbolError as exc:
        ambiguity_error = str(exc)
    else:
        ambiguity_error = "NO ERROR"
    unambiguous_fold = dict(ambiguous_fold)
    del unambiguous_fold["s11_a1"]
    report = check_consumer(unambiguous_fold, {"kept"})
    operands = {"manifest": ["kept"], "colliding_keys": ["a1", "s11_a1"]}
    observed = {"ambiguity_error": ambiguity_error, "closure_without_collision": sorted(report["closure"])}
    emit("f9c_ambiguity", operands, observed)
    assert "competing write-keys: a1, s11_a1" in ambiguity_error
    assert ambiguity_error != "NO ERROR"
    assert report["closure"] == frozenset({"kept", "a1"})


def test_bidirectional_smoke() -> None:
    fold = {"declared": row(1, "B"), "other": row(2, "B")}

    def undeclared(proxy: dict) -> None:
        proxy["declared"]
        proxy["other"]

    def omitted(proxy: dict) -> None:
        proxy["declared"]

    def exact(proxy: dict) -> tuple[object, object]:
        return proxy["declared"], proxy["other"]

    try:
        assert_lookups_equal_manifest(undeclared, fold, {"declared"})
    except ManifestError as exc:
        undeclared_error = str(exc)
    else:
        undeclared_error = "NO ERROR"
    try:
        assert_lookups_equal_manifest(omitted, fold, {"declared", "other"})
    except ManifestError as exc:
        unused_error = str(exc)
    else:
        unused_error = "NO ERROR"
    exact_report = assert_lookups_equal_manifest(exact, fold, {"declared", "other"})
    operands = {"fold_keys": sorted(fold), "manifests": [["declared"], ["declared", "other"]]}
    observed = {
        "undeclared_error": undeclared_error,
        "unused_error": unused_error,
        "exact_lookups": sorted(exact_report["lookups"]),
    }
    emit("bidirectional_smoke", operands, observed)
    assert "other" in undeclared_error and "undeclared" in undeclared_error
    assert "other" in unused_error and "declared-but-unused" in unused_error
    assert exact_report["lookups"] == frozenset({"declared", "other"})


def test_minimum_mode() -> None:
    exact_delta = {"bound": row(1, "D"), "infra": row(2, "D")}
    accumulated_delta = {**exact_delta, "unbound_history": row(3, "BASE")}
    try:
        assert_delta_is_minimal(accumulated_delta, {"bound"}, {"infra"})
    except LedgerFoldError as exc:
        extra_error = str(exc)
    else:
        extra_error = "NO ERROR"
    report = assert_delta_is_minimal(exact_delta, {"bound"}, {"infra"})
    optional_infra_report = assert_delta_is_minimal(
        {"bound": exact_delta["bound"]}, {"bound"}, {"infra"}
    )
    operands = {"own_bind_closure": ["bound"], "infra_keys": ["infra"]}
    observed = {
        "extra_error": extra_error,
        "exact_keys": sorted(report["exported_keys"]),
        "keys_with_optional_infra_absent": sorted(optional_infra_report["exported_keys"]),
    }
    emit("minimum_mode", operands, observed)
    assert "unbound_history" in extra_error and extra_error != "NO ERROR"
    assert report["exported_keys"] == frozenset({"bound", "infra"})


def test_promotion_delta() -> None:
    symbol = sp.Symbol("promoted")
    evidence = {
        "step": "PRODUCER",
        "f9_operands": sp.Tuple(symbol, sp.Tuple()),
        "route": {"decision": "F9a", "write_key": "promoted"},
    }
    delta = promotion_delta(
        "promoted",
        sp.srepr(symbol),
        "SYNTHETIC",
        evidence,
        display="promoted",
    )
    promoted = delta["promoted"]
    operands = {"row_key": "promoted", "srepr": sp.srepr(symbol), "evidence": evidence}
    observed = {"delta_keys": sorted(delta), "row_fields": sorted(promoted), "value": promoted["value"]}
    emit("promotion_delta", operands, observed)
    assert promoted["value"] == symbol
    assert promoted["route"] == evidence["route"]
    assert promoted["f9_operands"] == evidence["f9_operands"]
    assert {"display", "value", "value_kind", "class", "step"}.issubset(promoted)


def main() -> int:
    started = time.perf_counter()
    failures: list[tuple[str, BaseException]] = []
    with tempfile.TemporaryDirectory(prefix="ledger_fold_test_") as temporary:
        base_path, delta_path = write_synthetic_exports(Path(temporary))
        fold, audit = load_model(base_path, delta_path)
        tests = [
            ("last_wins", lambda: test_last_wins(fold, audit, base_path, delta_path)),
            ("no_f9_reapply", lambda: test_no_f9_reapply(fold)),
            ("predecessor_trap", lambda: test_predecessor_trap(fold)),
            ("missing_symbol_closure", test_missing_symbol_closure),
            ("dimension_closure", test_dimension_closure),
            ("f9c_ambiguity", test_f9c_ambiguity),
            ("bidirectional_smoke", test_bidirectional_smoke),
            ("minimum_mode", test_minimum_mode),
            ("promotion_delta", test_promotion_delta),
        ]
        for name, test in tests:
            try:
                test()
            except BaseException as exc:  # report every independent infrastructure test
                failures.append((name, exc))
                print(f"TEST {name} FAIL {type(exc).__name__}: {exc}")
                traceback.print_exc()
            else:
                print(f"TEST {name} PASS")

    runtime = time.perf_counter() - started
    print(f"SUMMARY tests={len(tests)} failures={len(failures)} runtime_seconds={runtime:.6f}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
