#!/usr/bin/env python3
"""Decisive synthetic tests plus a guarded frozen-ledger test for ledger_fold.py."""

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


REAL_LEDGER_ROOTS = (
    "mu_theta_operator",
    "slab_operator",
    "coupling_kernel",
    "face_normal",
    "conormal_deriv",
    "face_measure_shape_deriv",
    "face_velocity",
    "relative_flux",
    "kinematic_balance",
    "traction",
    "face_shift",
    "closure_shape_deriv",
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


def test_identity_resolution() -> None:
    omega_real = sp.Symbol("omega", real=True)
    omega_plain = sp.Symbol("omega")
    fold = {
        "kept_real": row(omega_real + 1, "ROOT"),
        "kept_plain": row(omega_plain + 1, "ROOT"),
        "omega": row(omega_real, "BARE"),
        "s11b_omega": row(omega_plain, "PREFIXED"),
    }
    reports: dict[str, object] = {}
    errors: dict[str, str] = {}
    for root in ("kept_real", "kept_plain"):
        try:
            reports[root] = check_consumer(fold, {root})
        except LedgerFoldError as exc:
            errors[root] = f"{type(exc).__name__}: {exc}"
    operands = {
        "manifest_roots": ["kept_real", "kept_plain"],
        "producer_identities": {
            "omega": sp.srepr(omega_real),
            "s11b_omega": sp.srepr(omega_plain),
        },
    }
    observed = {
        "errors": errors,
        "closures": {
            root: sorted(report["closure"])
            for root, report in reports.items()
        },
        "edges": {root: report["symbol_edges"] for root, report in reports.items()},
    }
    emit("identity_resolution", operands, observed)
    assert not errors
    assert reports["kept_real"]["closure"] == frozenset({"kept_real", "omega"})
    assert reports["kept_plain"]["closure"] == frozenset(
        {"kept_plain", "s11b_omega"}
    )


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


def test_depth_2_recursion() -> None:
    leaf = sp.Symbol("leaf", positive=True)
    fold = {
        "root": row(sp.Integer(1), "ROOT", dimension_key="mid"),
        "mid": row(leaf + 1, "MID"),
        "leaf": row(leaf, "LEAF"),
    }
    report = check_consumer(fold, {"root"})
    operands = {
        "manifest": ["root"],
        "chain": ["root", "mid", "leaf"],
        "edge_kinds": ["dimension_key", "atom_identity"],
    }
    observed = {
        "closure": sorted(report["closure"]),
        "symbol_edges": report["symbol_edges"],
        "dimension_edges": report["dimension_edges"],
    }
    emit("depth_2_recursion", operands, observed)
    assert report["closure"] == frozenset({"root", "mid", "leaf"})


def test_recursive_dimension_closure() -> None:
    missing_fold = {
        "root": row(sp.Integer(1), "ROOT", dimension_key="mid"),
        "mid": row(sp.Integer(2), "MID", dimension_key="leaf"),
    }
    try:
        check_consumer(missing_fold, {"root"})
    except ClosureError as exc:
        missing_error = str(exc)
    else:
        missing_error = "NO ERROR"
    complete_fold = dict(missing_fold)
    complete_fold["leaf"] = row(sp.Integer(3), "LEAF")
    report = check_consumer(complete_fold, {"root"})
    operands = {
        "manifest": ["root"],
        "dimension_chain": ["root", "mid", "leaf"],
        "missing_keys": sorted(missing_fold),
    }
    observed = {
        "missing_error": missing_error,
        "closure_after_add": sorted(report["closure"]),
    }
    emit("recursive_dimension_closure", operands, observed)
    assert "leaf" in missing_error and "mid" in missing_error
    assert missing_error != "NO ERROR"
    assert report["closure"] == frozenset({"root", "mid", "leaf"})


def test_structural_skip() -> None:
    x = sp.Symbol("x")
    matrix = sp.MatrixSymbol("M_structural", 2, 2)
    value = sp.Tuple(sp.Function("O_window")(x), matrix)
    fold = {"kept": row(value, "ROOT")}
    error = None
    report = None
    try:
        report = check_consumer(fold, {"kept"})
    except LedgerFoldError as exc:
        error = f"{type(exc).__name__}: {exc}"
    operands = {
        "manifest": ["kept"],
        "value": value,
        "absent_structural_atoms": ["O_window", "x", "M_structural"],
    }
    observed = {
        "error": error,
        "closure": None if report is None else sorted(report["closure"]),
        "symbol_edges": None if report is None else report["symbol_edges"],
    }
    emit("structural_skip", operands, observed)
    assert error is None
    assert report is not None and report["closure"] == frozenset({"kept"})
    assert report["symbol_edges"] == ()


def test_genuine_ambiguity() -> None:
    z = sp.Symbol("z")
    ambiguous_fold = {
        "kept": row(z + 1, "ROOT"),
        "z": row(z, "BASE"),
        "s11_z": row(z, "DELTA"),
    }
    try:
        check_consumer(ambiguous_fold, {"kept"})
    except AmbiguousSymbolError as exc:
        ambiguity_error = str(exc)
    else:
        ambiguity_error = "NO ERROR"
    operands = {
        "manifest": ["kept"],
        "identity": sp.srepr(z),
        "colliding_keys": ["z", "s11_z"],
    }
    observed = {"ambiguity_error": ambiguity_error}
    emit("genuine_ambiguity", operands, observed)
    assert "competing write-keys: s11_z, z" in ambiguity_error
    assert sp.srepr(z) in ambiguity_error
    assert ambiguity_error != "NO ERROR"


def test_user_is_not_a_producer() -> None:
    w_0 = sp.Symbol("W_0")
    fold = {
        "kept": row(2 * w_0, "ROOT"),
        "W_0": row(w_0, "DECLARATION"),
        "uses_w": row(w_0 + 1, "PREMISE"),
    }
    error = None
    report = None
    try:
        report = check_consumer(fold, {"kept"})
    except LedgerFoldError as exc:
        error = f"{type(exc).__name__}: {exc}"
    operands = {
        "manifest": ["kept"],
        "declaration_key": "W_0",
        "symbol_user_key": "uses_w",
        "identity": sp.srepr(w_0),
    }
    observed = {
        "error": error,
        "closure": None if report is None else sorted(report["closure"]),
        "symbol_edges": None if report is None else report["symbol_edges"],
    }
    emit("user_is_not_a_producer", operands, observed)
    assert error is None
    assert report is not None
    assert report["closure"] == frozenset({"kept", "W_0"})
    assert "uses_w" not in report["closure"]


def test_bidirectional_smoke() -> None:
    fold = {"declared": row(1, "B"), "other": row(2, "B")}

    def undeclared(proxy: dict) -> None:
        proxy["declared"]
        proxy["other"]

    def omitted(proxy: dict) -> None:
        proxy["declared"]

    def exact(proxy: dict) -> tuple[object, object]:
        return proxy["declared"], proxy["other"]

    def absent(proxy: dict) -> None:
        proxy["absent"]

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
    try:
        assert_lookups_equal_manifest(absent, fold, {"absent"})
    except KeyError as exc:
        absent_error = repr(exc)
    else:
        absent_error = "NO ERROR"
    exact_report = assert_lookups_equal_manifest(exact, fold, {"declared", "other"})
    operands = {
        "fold_keys": sorted(fold),
        "manifests": [["declared"], ["declared", "other"], ["absent"]],
    }
    observed = {
        "undeclared_error": undeclared_error,
        "unused_error": unused_error,
        "absent_required_error": absent_error,
        "exact_lookups": sorted(exact_report["lookups"]),
    }
    emit("bidirectional_smoke", operands, observed)
    assert "other" in undeclared_error and "undeclared" in undeclared_error
    assert "other" in unused_error and "declared-but-unused" in unused_error
    assert "absent" in absent_error and absent_error != "NO ERROR"
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


def test_minimum_mode_missing_required() -> None:
    delta = {"infra": row(1, "D")}
    try:
        assert_delta_is_minimal(delta, {"bound"}, {"infra"})
    except LedgerFoldError as exc:
        missing_error = str(exc)
    else:
        missing_error = "NO ERROR"
    operands = {
        "delta_keys": sorted(delta),
        "own_bind_closure": ["bound"],
        "infra_keys": ["infra"],
    }
    observed = {"missing_error": missing_error}
    emit("minimum_mode_missing_required", operands, observed)
    assert "missing row(s): bound" in missing_error
    assert missing_error != "NO ERROR"


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


def test_promotion_delta_guards() -> None:
    symbol = sp.Symbol("promoted")
    complete_evidence = {
        "step": "PRODUCER",
        "f9_operands": sp.Tuple(symbol, sp.Tuple()),
        "route": {"decision": "F9a", "write_key": "promoted"},
        "value": sp.Integer(0),
        "class": "FORBIDDEN",
    }
    try:
        promotion_delta("promoted", sp.srepr(symbol), "SYNTHETIC", complete_evidence)
    except LedgerFoldError as exc:
        forbidden_error = str(exc)
    else:
        forbidden_error = "NO ERROR"
    try:
        promotion_delta(
            "promoted",
            sp.srepr(symbol),
            "SYNTHETIC",
            {"step": "PRODUCER"},
        )
    except LedgerFoldError as exc:
        missing_error = str(exc)
    else:
        missing_error = "NO ERROR"
    operands = {
        "row_key": "promoted",
        "forbidden_evidence_fields": ["class", "value"],
        "incomplete_evidence_fields": ["step"],
    }
    observed = {"forbidden_error": forbidden_error, "missing_error": missing_error}
    emit("promotion_delta_guards", operands, observed)
    assert "class, value" in forbidden_error and forbidden_error != "NO ERROR"
    assert "f9_operands, route" in missing_error and missing_error != "NO ERROR"


def test_real_ledger_roots() -> None:
    ledger_path = Path(__file__).with_name("S11c_b_exports.py")
    operands = {"ledger_path": str(ledger_path), "roots": list(REAL_LEDGER_ROOTS)}
    if not ledger_path.is_file():
        emit("real_ledger_roots", operands, {"guard": "SKIPPED: ledger file absent"})
        return

    from S11c_b_exports import LEDGER as real_ledger

    closures: dict[str, list[str]] = {}
    errors: dict[str, str] = {}
    for root in REAL_LEDGER_ROOTS:
        try:
            report = check_consumer(real_ledger, {root})
        except BaseException as exc:
            errors[root] = f"{type(exc).__name__}: {exc}"
        else:
            closures[root] = sorted(report["closure"])
    observed = {
        "ledger_rows": len(real_ledger),
        "resolved_roots": sorted(closures),
        "closure_sizes": {root: len(keys) for root, keys in closures.items()},
        "errors": errors,
    }
    emit("real_ledger_roots", operands, observed)
    assert not errors
    assert set(closures) == set(REAL_LEDGER_ROOTS)


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
            ("identity_resolution", test_identity_resolution),
            ("dimension_closure", test_dimension_closure),
            ("depth_2_recursion", test_depth_2_recursion),
            ("recursive_dimension_closure", test_recursive_dimension_closure),
            ("structural_skip", test_structural_skip),
            ("genuine_ambiguity", test_genuine_ambiguity),
            ("user_is_not_a_producer", test_user_is_not_a_producer),
            ("bidirectional_smoke", test_bidirectional_smoke),
            ("minimum_mode", test_minimum_mode),
            ("minimum_mode_missing_required", test_minimum_mode_missing_required),
            ("promotion_delta", test_promotion_delta),
            ("promotion_delta_guards", test_promotion_delta_guards),
            ("real_ledger_roots", test_real_ledger_roots),
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
