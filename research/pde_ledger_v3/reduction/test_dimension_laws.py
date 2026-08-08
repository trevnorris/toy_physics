"""Regression tests for bound registry dimension laws."""

from __future__ import annotations

import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest
import sympy as sp

import dimension_law_able_to_fail as able_to_fail
import engine_dimension_pin
import w3_acceptance_ablations as acceptance_ablations
from dimensional_homogeneity_gate import DimensionAudit, HOMOGENEOUS, UNDETERMINED
from dimension_law_check import (
    inspect_declaration_forms,
    inspect_symbolic_witness_coverage,
    print_dimension_law_check,
)
from dimension_laws import BoundDimensionLaw, as_symbolic_dimension, dimension_residual
from engine_dimension_pin import (
    ENGINE_RECORDS,
    REQUIRED_QIDS,
    EngineDimensionPin,
    _normalise_engine_operand,
    inspect_engine_dimension_pin,
)
from registry_read import Registry, RegistryValidationError, load_raw_documents, load_registry


HERE = Path(__file__).resolve().parent


def _row(quantities: dict, qid: str) -> dict:
    return next(row for row in quantities["quantities"] if row["qid"] == qid)


def test_brane_constituents_are_bound_laws_and_speeds_are_constant_vectors() -> None:
    registry = load_registry()
    assert all(
        isinstance(registry.dimension_declaration(qid), BoundDimensionLaw)
        for qid in ("rho_br", "mu_R", "B_comp")
    )
    assert registry.dimension_declaration("c_gamma") == (1, -1, 0)
    assert registry.dimension_declaration("c_L") == (1, -1, 0)
    assert inspect_declaration_forms(registry).passed
    assert inspect_symbolic_witness_coverage(registry).passed


def test_legacy_quantity_dimension_is_the_checked_reference_view() -> None:
    registry = load_registry()
    for qid in ("rho_br", "mu_R", "B_comp"):
        quantity = registry.quantities[registry.resolve_qid(qid)]
        law = registry.dimension_law(qid)
        assert law is not None
        assert quantity.dimension == quantity.dimension_reference
        assert quantity.dimension == law.evaluate_reference()


@pytest.mark.parametrize("dimension", (2, 3, 4))
def test_consumer_can_request_an_evaluation_at_a_stated_dimension(
    dimension: int,
) -> None:
    registry = load_registry()
    pin = inspect_engine_dimension_pin(registry)
    assert pin.passed
    engine_operand_by_qid = {
        qid: next(
            observation.engine_operand
            for observation in pin.observations
            if observation.qid == qid
        )
        for qid in REQUIRED_QIDS
    }
    for qid, engine_operand in engine_operand_by_qid.items():
        law = registry.dimension_law(qid)
        assert law is not None
        structural_values = {binding_qid: dimension for binding_qid in law.binding_qids}
        substitutions = {
            sp.Symbol(binding_qid, integer=True): sp.Integer(value)
            for binding_qid, value in structural_values.items()
        }
        evaluated_engine_operand = tuple(
            int(component.subs(substitutions)) for component in engine_operand
        )
        assert registry.dimension_at(qid, structural_values) == evaluated_engine_operand


def test_engine_dimension_pin_rejects_empty_configured_population(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(engine_dimension_pin, "ENGINE_RECORDS", ())
    assert not EngineDimensionPin((), ()).passed


def test_engine_dimension_pin_rejects_subset_population_with_required_qids() -> None:
    registry = load_registry()
    subset = tuple(record for record in ENGINE_RECORDS if record.engine == "S11-py")
    pin = inspect_engine_dimension_pin(registry, subset)
    assert pin.covered_qids == REQUIRED_QIDS
    assert all(observation.passed for observation in pin.observations)
    assert not pin.passed


def test_engine_dimension_pin_rejects_nonempty_errors() -> None:
    pin = inspect_engine_dimension_pin(load_registry())
    assert pin.passed
    assert not replace(pin, errors=("synthetic read failure",)).passed


def test_engine_dimension_pin_rejects_unmapped_symbol() -> None:
    law = load_registry().dimension_law("Q.brane.rho_br")
    assert law is not None
    with pytest.raises(ValueError, match="unmapped engine dimension symbols"):
        _normalise_engine_operand(
            (sp.Symbol("D"), sp.Symbol("unmapped"), sp.Integer(1)),
            "D",
            law,
        )


def test_modulus_minus_density_laws_cancel_the_bound_dimension() -> None:
    registry = load_registry()
    rho = as_symbolic_dimension(registry.dimension_declaration("rho_br"))
    mu = as_symbolic_dimension(registry.dimension_declaration("mu_R"))
    bulk = as_symbolic_dimension(registry.dimension_declaration("B_comp"))
    speed_squared = tuple(
        2 * value
        for value in as_symbolic_dimension(registry.dimension_declaration("c_gamma"))
    )
    assert dimension_residual(mu, rho) == (sp.Integer(2), sp.Integer(-2), sp.Integer(0))
    assert dimension_residual(bulk, rho) == (sp.Integer(2), sp.Integer(-2), sp.Integer(0))
    assert dimension_residual(dimension_residual(mu, rho), speed_squared) == (0, 0, 0)
    assert dimension_residual(dimension_residual(bulk, rho), speed_squared) == (0, 0, 0)


def test_reference_values_must_equal_the_bound_quantities_declared_values() -> None:
    quantities, relations, schema = load_raw_documents()
    _row(quantities, "Q.brane.D_brane")["value"] = 4
    audit = DimensionAudit(quantities, relations)
    by_id = {result.relation_id: result for result in audit.classify_relations()}
    assert by_id["R4"].status == UNDETERMINED
    assert by_id["R5"].status == UNDETERMINED
    assert "reference=3 declared=4 residual=1" in by_id["R4"].detail
    with pytest.raises(
        RegistryValidationError,
        match=r"reference=3 declared=4 residual=1",
    ):
        Registry.from_documents(quantities, relations, schema)


@pytest.mark.parametrize("value", (None, ["Rat", 3, 2]))
def test_gate_is_undetermined_when_a_binding_has_no_integer_value(value: object) -> None:
    quantities, relations, schema = load_raw_documents()
    row = _row(quantities, "Q.brane.D_brane")
    if value is None:
        row.pop("value")
    else:
        row["value"] = value
    audit = DimensionAudit(quantities, relations)
    by_id = {result.relation_id: result for result in audit.classify_relations()}
    assert by_id["R4"].status == UNDETERMINED
    assert by_id["R5"].status == UNDETERMINED
    assert "cannot resolve to an integer" in by_id["R4"].detail
    with pytest.raises(RegistryValidationError, match="cannot resolve to an integer"):
        Registry.from_documents(quantities, relations, schema)


def test_retained_reference_vector_is_bound_to_the_law_evaluation() -> None:
    quantities, relations, schema = load_raw_documents()
    _row(quantities, "Q.brane.rho_br")["dimension"]["exponents"] = [-9, 0, 1]
    audit = DimensionAudit(quantities, relations)
    by_id = {result.relation_id: result for result in audit.classify_relations()}
    assert by_id["R4"].status == UNDETERMINED
    assert by_id["R5"].status == UNDETERMINED
    assert "reference evaluation differs from exponents" in by_id["R4"].detail
    with pytest.raises(
        RegistryValidationError,
        match="reference evaluation differs from exponents",
    ):
        Registry.from_documents(quantities, relations, schema)


def test_binding_target_must_be_structural() -> None:
    quantities, relations, schema = load_raw_documents()
    _row(quantities, "Q.medium.hbar")["value"] = 3
    for qid in ("Q.brane.rho_br", "Q.brane.mu_R", "Q.brane.B_comp"):
        law = _row(quantities, qid)["dimension"]["law"]
        law["bindings"]["D"] = "Q.medium.hbar"
    audit = DimensionAudit(quantities, relations)
    by_id = {result.relation_id: result for result in audit.classify_relations()}
    assert by_id["R4"].status == UNDETERMINED
    assert by_id["R5"].status == UNDETERMINED
    assert "targets non-structural quantity" in by_id["R4"].detail
    with pytest.raises(
        RegistryValidationError,
        match="targets non-structural quantity",
    ):
        Registry.from_documents(quantities, relations, schema)


@pytest.mark.parametrize("mutation", ("absent", "wrong-D-coefficient"))
def test_reduction_reports_the_honest_acceptance_surface(mutation: str) -> None:
    quantities, relations, schema = load_raw_documents()
    law_rows = [
        _row(quantities, qid)
        for qid in ("Q.brane.rho_br", "Q.brane.mu_R", "Q.brane.B_comp")
    ]
    if mutation == "absent":
        for row in law_rows:
            row["dimension"].pop("law")
    else:
        law_rows[0]["dimension"]["law"]["components"][0] = [
            "Sub",
            3,
            ["Mul", 2, ["Ref", "D"]],
        ]
        for row in law_rows[1:]:
            row["dimension"]["law"]["components"][0] = [
                "Sub",
                5,
                ["Mul", 2, ["Ref", "D"]],
            ]
    registry = Registry.from_documents(quantities, relations, schema)
    forms = inspect_declaration_forms(registry)
    coverage = inspect_symbolic_witness_coverage(registry)
    gate_results = DimensionAudit(quantities, relations).classify_relations()
    assert forms.passed is (mutation == "wrong-D-coefficient")
    assert not coverage.passed
    assert all(result.status == HOMOGENEOUS for result in gate_results)


@pytest.mark.parametrize(
    ("mutation", "detail"),
    (
        ("unbound", "unbound dimension-law symbol"),
        ("unresolvable", "is unresolvable"),
    ),
)
def test_unbound_and_unresolvable_laws_are_undetermined(
    mutation: str, detail: str
) -> None:
    quantities, relations, schema = load_raw_documents()
    law = _row(quantities, "Q.brane.rho_br")["dimension"]["law"]
    if mutation == "unbound":
        law["components"][0] = ["Neg", ["Ref", "missing_D"]]
    else:
        law["bindings"]["D"] = "Q.missing.D"
    audit = DimensionAudit(quantities, relations)
    by_id = {result.relation_id: result for result in audit.classify_relations()}
    assert by_id["R4"].status == UNDETERMINED
    assert by_id["R5"].status == UNDETERMINED
    assert detail in by_id["R4"].detail
    with pytest.raises(RegistryValidationError, match=detail):
        Registry.from_documents(quantities, relations, schema)


@pytest.mark.parametrize(
    "case_name",
    (
        "absent-laws",
        "wrong-D-coefficient",
        "missing-binding-value",
        "noninteger-binding-value",
        "reference-value-mismatch",
        "reference-vector-mismatch",
        "nonstructural-binding",
        "wrong-law",
        "wrong-binding",
        "unbound",
        "unresolvable",
    ),
)
def test_able_to_fail_dimension_law_cases(case_name: str) -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "dimension_law_able_to_fail.py"), "--case", case_name],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 1
    assert not completed.stderr
    assert f"DIMENSION_LAW_ABLE_TO_FAIL_CAUGHT: {case_name}" in completed.stdout
    assert "OPERAND" in completed.stdout
    assert "RESIDUAL" in completed.stdout


def test_dimension_law_able_to_fail_aggregate() -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "dimension_law_able_to_fail.py")],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0
    assert not completed.stderr
    assert "DIMENSION_LAW_ABLE_TO_FAIL_HARNESS: PASS" in completed.stdout


def test_duplicate_pin_ablation_reports_engine_pin_discrimination() -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "w3_duplicate_pin_ablation.py")],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0
    assert not completed.stderr
    assert "DUPLICATE_PIN_EXIT_OPERAND baseline: 0" in completed.stdout
    assert "DUPLICATE_PIN_EXIT_OPERAND changed: 1" in completed.stdout
    assert "DUPLICATE_PIN_EXIT_RESIDUAL changed-minus-baseline: 1" in completed.stdout
    assert completed.stdout.count("D_COEFFICIENT_POLICED_IN_REDUCTION: YES") == 1
    assert completed.stdout.count("D_COEFFICIENT_POLICED_IN_REDUCTION: NO") == 1
    assert "W3_DUPLICATE_PIN_ABLATION: PASS" in completed.stdout


def test_dimension_law_able_to_fail_escape_fails_the_aggregate() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(HERE / "dimension_law_able_to_fail.py"),
            "--demonstrate-escape",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 1
    assert not completed.stderr
    assert "status=ESCAPED observed_exit=0" in completed.stdout
    assert "ESCAPE_OPERAND observed-child-exit: 0" in completed.stdout
    assert "ESCAPE_RESIDUAL observed-minus-expected: -1" in completed.stdout
    assert "DIMENSION_LAW_ABLE_TO_FAIL_HARNESS: FAIL" in completed.stdout
    assert "DIMENSION_LAW_ABLE_TO_FAIL_HARNESS: PASS" not in completed.stdout


def test_demonstrate_escape_uses_the_observed_child_process(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    invoked = False

    def completed_child(*_args: object, **_kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal invoked
        invoked = True
        case_name = "wrong-law"
        return subprocess.CompletedProcess(
            args=[],
            returncode=able_to_fail.EXPECTED_FAILURE_EXIT,
            stdout=f"DIMENSION_LAW_ABLE_TO_FAIL_CAUGHT: {case_name}\n",
            stderr="",
        )

    monkeypatch.setattr(able_to_fail.subprocess, "run", completed_child)
    observed_exit = able_to_fail.demonstrate_escape()
    stdout = capsys.readouterr().out
    assert invoked
    assert observed_exit == able_to_fail.EXPECTED_FAILURE_EXIT
    assert (
        f"ESCAPE_OPERAND observed-child-exit: {able_to_fail.EXPECTED_FAILURE_EXIT}"
        in stdout
    )
    assert "ESCAPE_RESIDUAL observed-minus-expected: 0" in stdout


def test_acceptance_subverdict_is_guarded(
    capsys: pytest.CaptureFixture[str],
) -> None:
    observations = [
        acceptance_ablations._expected(case_name)
        for case_name in acceptance_ablations.CASES
    ]
    observations = [
        replace(
            observation,
            coefficient_pin_passed=not observation.coefficient_pin_passed,
        )
        if observation.case_name == "wrong-D-coefficient"
        else observation
        for observation in observations
    ]
    passed = acceptance_ablations.evaluate(observations)
    stdout = capsys.readouterr().out
    assert not passed
    assert "ACCEPTANCE_CASE_GUARD wrong-D-coefficient: FAIL" in stdout


def test_printed_dimension_law_status_is_guarded(
    capsys: pytest.CaptureFixture[str],
) -> None:
    quantities, relations, schema = load_raw_documents()
    _row(quantities, "Q.brane.rho_br")["dimension"]["law"]["components"][0] = [
        "Sub",
        3,
        ["Mul", 2, ["Ref", "D"]],
    ]
    for qid in ("Q.brane.mu_R", "Q.brane.B_comp"):
        _row(quantities, qid)["dimension"]["law"]["components"][0] = [
            "Sub",
            5,
            ["Mul", 2, ["Ref", "D"]],
        ]
    registry = Registry.from_documents(quantities, relations, schema)
    forms = inspect_declaration_forms(registry)
    coverage = inspect_symbolic_witness_coverage(registry)
    passed = print_dimension_law_check(registry, forms, coverage)
    stdout = capsys.readouterr().out
    printed = next(
        line.rsplit(": ", 1)[1]
        for line in stdout.splitlines()
        if line.startswith("DIMENSION_LAW_CHECK:")
    )
    assert printed == ("PASS" if passed else "FAIL")
    assert not passed


def test_generate_rows_is_explicitly_retired(tmp_path: Path) -> None:
    output_dir = tmp_path / "generated"
    completed = subprocess.run(
        [
            sys.executable,
            str(HERE / "generate_rows.py"),
            "--output-dir",
            str(output_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 2
    assert not completed.stdout
    assert "generate_rows.py is retired" in completed.stderr
    assert not output_dir.exists()


def test_w4_weaker_implementation_population() -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "w4_weaker_implementation_runs.py")],
        check=False,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert completed.returncode == 0
    assert not completed.stderr
    assert completed.stdout.count("FAIL_AS_REQUIRED") == 3
    assert "W4_WEAKER_IMPLEMENTATIONS: PASS" in completed.stdout


def test_shifted_registry_printer_reports_nonzero_declared_binding_value() -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "w4_shifted_registry_printer.py")],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 1
    assert not completed.stderr
    assert completed.stdout.count("S11_NUMERIC_REFERENCE_COMPARISON") == 3
    assert completed.stdout.count(
        "residual=[2 - Q.brane.D_brane,0,0]; "
        "declared-bindings=Q.brane.D_brane=3; "
        "residual-at-declared-bindings=[-1,0,0]"
    ) == 3
    assert "specialisation artefact" not in completed.stdout
    assert "physics disagreement" not in completed.stdout


def test_pin_completeness_weaker_implementation_population() -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "w4_pin_completeness_runs.py")],
        check=False,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert completed.returncode == 0
    assert not completed.stderr
    assert completed.stdout.count("PIN_WEAKER_TEST_GUARD") == 4
    assert completed.stdout.count("FAIL_AS_REQUIRED") == 4
    assert "W4_PIN_COMPLETENESS_WEAKER_IMPLEMENTATIONS: PASS" in completed.stdout
