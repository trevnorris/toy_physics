from __future__ import annotations

import tempfile
from pathlib import Path

import sympy as sp
import yaml

import registry_dimension_witness as witness


HERE = Path(__file__).resolve().parent
CONVENTION = "LTM-exponent-vector-v1"
BASES = ("L", "T", "M")


def _quantity(
    qid: str = "Q.test.x",
    exponents: list[int] | None = None,
    convention: str | None = CONVENTION,
) -> witness.RegistryQuantity:
    return witness.RegistryQuantity(
        qid,
        qid.rpartition(".")[2],
        convention,
        exponents if exponents is not None else [-3, 0, 1],
        None,
    )


def _specialisation_quantity() -> witness.RegistryQuantity:
    return witness.RegistryQuantity(
        "Q.brane.D_brane",
        "D_brane",
        CONVENTION,
        [0, 0, 0],
        3,
    )


def _measured_emission(
    raw: witness.Dimension,
    **changes: object,
) -> witness.Emission:
    values: dict[str, object] = {
        "tag": "TAG",
        "raw": raw,
        "source_class": "DERIVED",
        "axis_raw": raw,
        "axis_tag": "AXIS_TAG",
    }
    values.update(changes)
    return witness.Emission(**values)  # type: ignore[arg-type]


def test_source_bound_registry_specialisation_and_residual() -> None:
    quantity = _quantity()
    emission = _measured_emission(
        (-sp.Symbol("D"), sp.Integer(0), sp.Integer(1)),
        specialisations={"D": {"registry_value": "Q.brane.D_brane"}},
    )
    row = witness._build_row(
        "fixture",
        {"declared_emitted_convention": CONVENTION},
        quantity,
        [emission],
        {quantity.qid: quantity, "Q.brane.D_brane": _specialisation_quantity()},
        CONVENTION,
        BASES,
    )
    assert row.status == "AGREEMENT"
    assert row.emissions[0].value == (-3, 0, 1)
    assert row.emissions[0].difference == (0, 0, 0)
    assert row.emissions[0].applied == ("D=3<-registry(Q.brane.D_brane)",)


def test_branch_dimension_is_bound_to_source_and_mismatch_is_explicit() -> None:
    quantity = _quantity()
    raw = (-sp.Symbol("D"), sp.Integer(0), sp.Integer(1))
    emission = _measured_emission(
        raw,
        branch_dimension=5,
        specialisations={"D": {"branch_dimension_from": "Q.brane.D_brane"}},
    )
    row = witness._build_row(
        "fixture",
        {"declared_emitted_convention": CONVENTION},
        quantity,
        [emission],
        {quantity.qid: quantity, "Q.brane.D_brane": _specialisation_quantity()},
        CONVENTION,
        BASES,
    )
    assert row.status == "BRANCH_DIMENSION_MISMATCH"
    assert row.emissions[0].value == (-5, 0, 1)
    assert row.emissions[0].registry_branch_dimension == 3
    assert row.emissions[0].difference == (-2, 0, 0)


def test_declared_multiplier_is_applied_and_printed() -> None:
    quantity = _quantity("Q.brane.c_gamma", [1, -1, 0])
    emitted = (sp.Integer(2), sp.Integer(-2), sp.Integer(0))
    baseline = witness._build_row(
        "fixture",
        {"declared_emitted_convention": CONVENTION},
        quantity,
        [_measured_emission(emitted, multiplier=sp.Rational(2))],
        {quantity.qid: quantity},
        CONVENTION,
        BASES,
    )
    perturbed = witness._build_row(
        "fixture",
        {"declared_emitted_convention": CONVENTION},
        quantity,
        [_measured_emission(emitted, multiplier=sp.Rational(3))],
        {quantity.qid: quantity},
        CONVENTION,
        BASES,
    )
    assert baseline.status == "AGREEMENT"
    assert perturbed.status == "DISAGREEMENT"
    rendered = witness._row_line(perturbed)
    assert "|declared=[1,-1,0]|multiplier={3@TAG}" in rendered
    assert "|expected(multiplier*declared)={[3,-3,0]@TAG}" in rendered
    assert "|residual(emitted-multiplier*declared)={[-1,1,0]@TAG}" in rendered


def test_echoed_source_never_becomes_agreement() -> None:
    quantity = _quantity("Q.brane.D_brane", [0, 0, 0])
    zero = (sp.Integer(0), sp.Integer(0), sp.Integer(0))
    changed = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
    baseline = witness._build_row(
        "fixture",
        {"declared_emitted_convention": CONVENTION},
        quantity,
        [_measured_emission(zero, source_class="ECHOED")],
        {quantity.qid: quantity},
        CONVENTION,
        BASES,
    )
    perturbed = witness._build_row(
        "fixture",
        {"declared_emitted_convention": CONVENTION},
        quantity,
        [_measured_emission(changed, source_class="ECHOED")],
        {quantity.qid: quantity},
        CONVENTION,
        BASES,
    )
    assert baseline.status == "ECHOED"
    assert perturbed.status == "ECHOED_MISMATCH"


def test_convention_input_is_not_a_measurement() -> None:
    quantity = _quantity()
    emitted = (sp.Integer(-3), sp.Integer(0), sp.Integer(1))
    unlabelled = witness.Emission(
        "TAG",
        emitted,
        source_class="DERIVED",
        axis_error="NO_AXIS_LABELLED_COMPONENTS",
    )
    row = witness._build_row(
        "fixture",
        {"declared_emitted_convention": CONVENTION},
        quantity,
        [unlabelled],
        {quantity.qid: quantity},
        CONVENTION,
        BASES,
    )
    assert row.status == "UNDETERMINED"
    assert row.emissions[0].difference == (0, 0, 0)
    assert "emitted_convention_input=" in witness._row_line(row)


def test_template_expanding_to_zero_tags_produces_a_guarded_row() -> None:
    with tempfile.TemporaryDirectory(prefix="registry-witness-test-") as raw_temp:
        temporary = Path(raw_temp)
        config = {
            "dimension_sources": [
                {
                    "engine": "py",
                    "package": "MAIN",
                    "symbols": {
                        "rho_br": {
                            "tag_template": "PY_TEST_{package}_D{dimension}_DIM",
                            "dimensions": [3],
                        }
                    },
                }
            ]
        }
        output_path = temporary / "engine.out"
        output_path.write_text("PY_S9_TEST_UNUSED: 0\n", encoding="utf-8")
        config_path = temporary / "checks.yaml"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        template = "PY_TEST_{package}_D{dimension}_DIM"
        manifest = {
            "schema_version": "registry-dimension-witness-v2",
            "artifacts": [
                {
                    "id": "empty-template",
                    "engine": "py",
                    "config": str(config_path),
                    "output": str(output_path),
                    "declared_emitted_convention": CONVENTION,
                    "selection": {"packages": ["MAIN"], "dimensions": [5]},
                    "source_metadata": {template: {"source_class": "DERIVED"}},
                }
            ],
        }
        manifest_path = temporary / "manifest.yaml"
        manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
        reports = witness.build_report(manifest_path, HERE / "quantities.yaml")
    target_rows = [
        row
        for row in reports[0].rows
        if row.quantity.qid == "Q.brane.rho_br"
    ]
    assert len(target_rows) == 1
    assert target_rows[0].status == "NO_TAGS"
    assert reports[0].status == "OK"


def test_artifact_with_zero_rows_has_guarded_status() -> None:
    with tempfile.TemporaryDirectory(prefix="registry-witness-test-") as raw_temp:
        temporary = Path(raw_temp)
        output_path = temporary / "engine.out"
        output_path.write_text("PY_S9_TEST_UNUSED: 0\n", encoding="utf-8")
        config_path = temporary / "checks.yaml"
        config_path.write_text("dimension_sources: []\n", encoding="utf-8")
        manifest = {
            "schema_version": "registry-dimension-witness-v2",
            "artifacts": [
                {
                    "id": "zero-row-artifact",
                    "engine": "py",
                    "config": str(config_path),
                    "output": str(output_path),
                    "declared_emitted_convention": CONVENTION,
                }
            ],
        }
        manifest_path = temporary / "manifest.yaml"
        manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
        reports = witness.build_report(manifest_path, HERE / "quantities.yaml")
    assert reports[0].rows == ()
    assert reports[0].selected_quantities == ()
    assert reports[0].status == "NO_ROWS"
    assert "NO_ROWS=1" in witness.format_report(reports)


def test_committed_manifest_has_row_coverage_for_every_selected_pair() -> None:
    manifest = yaml.safe_load((HERE / "registry_dimension_witness.yaml").read_text(encoding="utf-8"))
    reports = witness.build_report(
        HERE / "registry_dimension_witness.yaml",
        HERE / "quantities.yaml",
    )
    manifest_artifacts = [artifact["id"] for artifact in manifest["artifacts"]]
    assert [report.artifact for report in reports] == manifest_artifacts
    for report in reports:
        selected_pairs = {(report.artifact, qid) for qid in report.selected_quantities}
        row_pairs = {(report.artifact, row.quantity.qid) for row in report.rows}
        assert selected_pairs <= row_pairs
        assert report.rows or report.status == "NO_ROWS"
    rendered = witness.format_report(reports)
    assert "scope_artifacts=[" + ",".join(manifest_artifacts) + "]" in rendered
