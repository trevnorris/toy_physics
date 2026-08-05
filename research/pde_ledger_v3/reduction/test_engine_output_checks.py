#!/usr/bin/env python3
"""Able-to-fail tests for the engine-output check harness.

No expected physics result is stored here.  The real-output mutations copy or
damage values discovered at test time, and registry targets are generated from
the loaded registry at comparison time.
"""

from __future__ import annotations

import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import pytest
import sympy as sp

import engine_output_checks as checks
from registry_read import load_registry


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
LEDGER = HERE.parent
REAL_S9_OUTPUTS = {
    "wl": LEDGER / "mathematica/out/S9_light_requires_shear_mathematica_audit.out",
    "py": LEDGER / "scripts/out/S9_light_requires_shear_sympy_audit.out",
}
CONFIG = HERE / "checks_S9.yaml"
MODULE = HERE / "engine_output_checks.py"


@pytest.fixture(scope="session", autouse=True)
def ensure_real_outputs() -> None:
    """Generate the two real engine emissions, sequentially, when absent."""
    if not REAL_S9_OUTPUTS["wl"].is_file():
        with REAL_S9_OUTPUTS["wl"].open("w", encoding="utf-8") as handle:
            subprocess.run(
                [
                    "timeout",
                    "600",
                    "math",
                    "-script",
                    "research/pde_ledger_v3/mathematica/"
                    "S9_light_requires_shear_mathematica_audit.wl",
                ],
                cwd=PROJECT,
                stdout=handle,
                check=True,
                timeout=620,
            )
    if not REAL_S9_OUTPUTS["py"].is_file():
        with REAL_S9_OUTPUTS["py"].open("w", encoding="utf-8") as handle:
            subprocess.run(
                [
                    "timeout",
                    "600",
                    sys.executable,
                    "scripts/S9_light_requires_shear_sympy_audit.py",
                ],
                cwd=LEDGER,
                stdout=handle,
                check=True,
                timeout=620,
            )


def real_text(engine: str) -> str:
    return REAL_S9_OUTPUTS[engine].read_text(encoding="utf-8")


@lru_cache(maxsize=2)
def _cached_real_raw_and_values(
    engine: str = "wl",
) -> tuple[dict[str, str], dict[str, checks.ParsedValue | checks.Unparsed]]:
    raw = checks.parse_tagged_output(real_text(engine))
    return raw, checks.normalize_tags(raw)


def real_raw_and_values(
    engine: str = "wl",
) -> tuple[dict[str, str], dict[str, checks.ParsedValue | checks.Unparsed]]:
    raw, values = _cached_real_raw_and_values(engine)
    return dict(raw), dict(values)


def render_tags(raw: dict[str, str]) -> str:
    return "\n".join(f"{tag}: {value}" for tag, value in raw.items()) + "\n"


def cli(output_paths: dict[str, Path]) -> subprocess.CompletedProcess[str]:
    arguments = [sys.executable, str(MODULE), "--config", str(CONFIG)]
    for engine, path in output_paths.items():
        arguments.extend(["--output", f"{engine}={path}"])
    return subprocess.run(
        arguments,
        cwd=LEDGER,
        check=False,
        capture_output=True,
        text=True,
        timeout=180,
    )


def test_real_s9_end_to_end_cli() -> None:
    completed = cli(REAL_S9_OUTPUTS)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "CROSS_ENGINE: agree=12" in completed.stdout
    assert "PARITY (" in completed.stdout
    assert "UNPARSED (0):" in completed.stdout
    assert "LIMITATION: triage only" in completed.stdout


def test_deleting_real_control_tag_fires_parity_missing() -> None:
    raw, values = real_raw_and_values("wl")
    control_prefix = "WL_S9_X1_"
    target = next(
        tag
        for tag in raw
        if tag.startswith(control_prefix)
        and f"WL_S9_{tag[len(control_prefix):]}" in raw
    )
    suffix = target[len(control_prefix) :]
    corrupted = dict(raw)
    del corrupted[target]

    corrupted_values = dict(values)
    del corrupted_values[target]
    parity = checks.check_tag_parity({"wl": corrupted_values})
    row = next(row for row in parity.rows if row.package == "WL_S9_X1")

    assert suffix in row.missing


def test_parity_exclusion_is_reported_and_confined_to_parity() -> None:
    outputs = {
        "wl": checks.normalize_tags(
            {
                "WL_COMMON": "1",
                "WL_LOCAL_ONLY": "local_value",
                "WL_X1_COMMON": "1",
            }
        ),
        "py": checks.normalize_tags(
            {
                "PY_COMMON": "1",
                "PY_LOCAL_ONLY": "local_value",
                "PY_X1_COMMON": "1",
            }
        ),
    }
    config = {
        "default_engine": "wl",
        "parity_exclude": ["_LOCAL_"],
        "cross_engine": [
            {
                "quantity": "engine_local_value",
                "wl": "WL_LOCAL_ONLY",
                "py": "PY_LOCAL_ONLY",
            }
        ],
        "registry_residual": [],
    }

    report = checks.run_checks(config, outputs)
    rendered = checks.format_report(report)

    assert all(
        "LOCAL_ONLY" not in row.missing and "LOCAL_ONLY" not in row.extra
        for row in report.parity.rows
    )
    assert {row.engine: len(row.tags) for row in report.parity.exclusions} == {
        "py": 1,
        "wl": 1,
    }
    assert "WL_LOCAL_ONLY" in report.controls.unpaired
    assert (
        checks.UnknownSymbol("WL_LOCAL_ONLY", "local_value")
        in report.dimensions.unknown_symbols
    )
    assert report.cross_engine.rows[0].status == "AGREE"
    assert "PARITY_EXCLUDED (2):" in rendered
    assert "py: excluded=1 by_pattern={'_LOCAL_':1}" in rendered
    assert "wl: excluded=1 by_pattern={'_LOCAL_':1}" in rendered


def test_absent_parity_exclude_preserves_report_bytes() -> None:
    outputs = {
        "wl": checks.normalize_tags(
            {
                "WL_VALUE": "1",
                "WL_LOCAL_ONLY": "2",
                "WL_X1_VALUE": "1",
            }
        )
    }
    config = {
        "default_engine": "wl",
        "cross_engine": [],
        "registry_residual": [],
    }

    rendered = checks.format_report(checks.run_checks(config, outputs))

    assert rendered == "\n".join(
        [
            "CONTROL_RESPONSE: compared=1 responsive=0 invariant=1 unparsed=0 unpaired=1",
            "DIMENSIONS: total=3 checked=3 homogeneous=3 non_homogeneous=0 unknown_symbol=0 unparsed=0",
            "NON_DIMENSIONAL: total=0 kinds=[none]",
            "CROSS_ENGINE: configured=0",
            "REGISTRY_RESIDUAL: configured=0",
            "TAG_PARITY: packages=1 gaps=1",
            "DISAGREE (0):",
            "INVARIANT (1): WL_VALUE",
            "NON_HOMOGENEOUS (0):",
            "PARITY (1):",
            "  why: present-and-identical is INVARIANT; absent is indistinguishable from never computed.",
            "  wl:packages=[WL_X1]: missing=[LOCAL_ONLY] extra=[-]",
            "UNPARSED (0):",
            "UNKNOWN_SYMBOL (0):",
            "DIMENSION_ERROR (0):",
            "LIMITATION: triage only; this run does not establish physical correctness, completeness, or derivation coverage.",
        ]
    )


def test_changing_real_mapped_operand_fires_disagree_with_both_operands() -> None:
    config = checks.load_config(CONFIG)
    row = config["cross_engine"][0]
    wl_raw, wl_values = real_raw_and_values("wl")
    py_raw, py_values = real_raw_and_values("py")
    corrupted = dict(py_raw)
    corrupted[row["py"]] = f"({py_raw[row['py']]}) + cross_engine_corruption"

    corrupted_values = dict(py_values)
    corrupted_values[row["py"]] = checks.normalize_sympy(corrupted[row["py"]])
    observed = checks.check_cross_engine(
        [row], {"wl": wl_values, "py": corrupted_values}
    ).rows[0]

    assert observed.status == "DISAGREE"
    assert {engine for engine, _operand in observed.operands} == {"wl", "py"}
    assert "cross_engine_corruption" in dict(observed.operands)["py"]
    assert wl_raw[row["wl"]]


def test_symbolic_comparison_accepts_different_algebra_and_vector_shapes() -> None:
    rows = [
        {"quantity": "algebra", "wl": "WL_EQ", "py": "PY_EQ"},
        {"quantity": "vector", "wl": "WL_VECTOR", "py": "PY_VECTOR"},
    ]
    outputs = {
        "wl": checks.normalize_tags({"WL_EQ": "a/b", "WL_VECTOR": "{1, 2}"}),
        "py": checks.normalize_tags(
            {"PY_EQ": "a*b**-1", "PY_VECTOR": "Matrix([[1], [2]])"}
        ),
    }

    report = checks.check_cross_engine(rows, outputs)

    assert [row.status for row in report.rows] == ["AGREE", "AGREE"]


def test_real_boolean_is_non_dimensional_not_unknown_or_unparsed() -> None:
    _raw, values = real_raw_and_values("wl")
    target = "WL_S9_ASSUMPTIONS"
    report = checks.check_dimensions(values, checks.load_config(CONFIG))

    assert values[target].kind == checks.ValueKind.BOOLEAN
    assert target in {row.tag for row in report.non_dimensional}
    assert target not in {row.tag for row in report.unknown_symbols}
    assert target not in report.unparsed


def test_unknown_symbol_in_real_dimensionful_expression_is_named_and_operational() -> None:
    config = checks.load_config(CONFIG)
    wl_raw, wl_values = real_raw_and_values("wl")
    _py_raw, py_values = real_raw_and_values("py")
    target = "WL_S9_LAGRANGIAN"
    corrupted = dict(wl_raw)
    corrupted[target] = f"({wl_raw[target]}) + unknown_dimension_symbol"
    corrupted_values = dict(wl_values)
    corrupted_values[target] = checks.normalize_mathematica(corrupted[target])
    report = checks.run_checks(
        config,
        {"wl": corrupted_values, "py": py_values},
    )

    assert checks.UnknownSymbol(target, "unknown_dimension_symbol") in report.dimensions.unknown_symbols
    assert report.operational_failures
    assert "unknown_dimension_symbol" in checks.format_report(report)


def test_disagree_is_a_finding_and_cli_exit_is_zero(tmp_path: Path) -> None:
    config = checks.load_config(CONFIG)
    row = config["cross_engine"][0]
    py_raw, _ = real_raw_and_values("py")
    corrupted = dict(py_raw)
    corrupted[row["py"]] = f"({py_raw[row['py']]}) + cli_disagreement"
    py_path = tmp_path / "s9_py_disagree.txt"
    py_path.write_text(render_tags(corrupted), encoding="utf-8")

    completed = cli({"wl": REAL_S9_OUTPUTS["wl"], "py": py_path})

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "DISAGREE (1):" in completed.stdout
    assert "cli_disagreement" in completed.stdout
    assert f"wl[{row['wl']}]=" in completed.stdout
    assert f"py[{row['py']}]=" in completed.stdout


def test_malformed_real_cas_is_unparsed_and_cli_exit_is_nonzero(tmp_path: Path) -> None:
    config = checks.load_config(CONFIG)
    row = config["cross_engine"][0]
    py_raw, _ = real_raw_and_values("py")
    corrupted = dict(py_raw)
    corrupted[row["py"]] = "this is not valid CAS text !"
    py_path = tmp_path / "s9_py_malformed.txt"
    py_path.write_text(render_tags(corrupted), encoding="utf-8")

    completed = cli({"wl": REAL_S9_OUTPUTS["wl"], "py": py_path})

    assert completed.returncode != 0
    assert "UNPARSED (1):" in completed.stdout
    assert "this is not valid CAS text !" in completed.stdout
    assert "OPERATIONAL_FAILURE" in completed.stderr


def test_control_response_becomes_invariant_when_sole_changed_control_is_copied() -> None:
    raw, values = real_raw_and_values()
    original = checks.check_control_response(values)
    candidate = next(
        row for row in original.responsive if len(row.differing_controls) == 1
    )

    # The compare-time main value is copied; no expected physics value is stored.
    corrupted = dict(raw)
    corrupted[candidate.differing_controls[0]] = raw[candidate.tag]
    corrupted_values = dict(values)
    corrupted_values[candidate.differing_controls[0]] = values[candidate.tag]
    observed = checks.check_control_response(corrupted_values)

    assert candidate.tag in {row.tag for row in observed.invariant}
    assert candidate.tag not in {row.tag for row in observed.responsive}


def test_mixed_dimension_sum_names_its_summands() -> None:
    raw, values = real_raw_and_values()
    config = checks.load_config(CONFIG)
    target = "WL_S9_LAGRANGIAN"
    corrupted = dict(raw)
    # omega's definitional dimension differs from the action expression's
    # engine-derived dimension; neither vector is asserted in this test.
    corrupted[target] = f"({raw[target]}) + omega"

    corrupted_values = dict(values)
    corrupted_values[target] = checks.normalize_mathematica(corrupted[target])
    report = checks.check_dimensions(corrupted_values, config)
    fired = [issue for issue in report.non_homogeneous if issue.tag == target]

    assert fired
    assert len(fired[0].summands) >= 2
    assert any(term.expression == "omega" for term in fired[0].summands)
    assert all(term.dimension is not None for term in fired[0].summands)


def test_invalid_real_value_is_unparsed_and_never_skipped() -> None:
    raw, values = real_raw_and_values()
    target = next(iter(raw))
    corrupted = dict(raw)
    corrupted[target] = "this is not valid CAS text !"

    normalized = dict(values)
    normalized[target] = checks.normalize_mathematica(corrupted[target])

    assert isinstance(normalized[target], checks.Unparsed)
    assert sum(isinstance(value, checks.Unparsed) for value in normalized.values()) == 1


def test_duplicate_real_tag_raises() -> None:
    text = real_text("wl")
    raw = checks.parse_tagged_output(text)
    tag = next(iter(raw))

    with pytest.raises(ValueError, match=f"duplicate emitted tag {tag}"):
        checks.parse_tagged_output(text + f"\n{tag}: {raw[tag]}\n")


def test_cross_engine_comparison_is_symbolic_in_both_directions() -> None:
    rows = [
        {"quantity": "algebraic_equivalence", "left": "LEFT_EQ", "right": "RIGHT_EQ"},
        {"quantity": "partial_text_overlap", "left": "LEFT_NE", "right": "RIGHT_NE"},
    ]
    outputs = {
        "left": checks.normalize_tags(
            {"LEFT_EQ": "a/b", "LEFT_NE": "a + b"}
        ),
        "right": checks.normalize_tags(
            {"RIGHT_EQ": "a*b^-1", "RIGHT_NE": "a + b + c"}
        ),
    }

    report = checks.check_cross_engine(rows, outputs)

    assert [row.status for row in report.rows] == ["AGREE", "DISAGREE"]


def test_root_multiset_comparison_is_explicitly_order_insensitive() -> None:
    left = checks.normalize_mathematica("{a/b, c}")
    right = checks.normalize_mathematica("{c, a*b^-1}")
    assert checks.symbolic_multiset_equal(left, right)


def test_registry_residual_target_is_generated_at_compare_time_and_can_fail() -> None:
    registry = load_registry()
    relation = next(
        relation
        for relation in registry.relations.values()
        if relation.rhs is not None and relation.designated_output is not None
    )
    inputs = {
        qid: sp.Symbol(f"engine_input_{index}", positive=True)
        for index, qid in enumerate(relation.input_qids)
    }
    generated_output = sp.simplify(
        relation.rhs.subs({registry.symbols[qid]: value for qid, value in inputs.items()})
    )
    all_values = {**inputs, relation.designated_output: generated_output}
    qid_tags = {qid: f"TAG_{index}" for index, qid in enumerate(all_values)}
    output = {tag: all_values[qid] for qid, tag in qid_tags.items()}
    row = {"relation_id": relation.relation_id, "engine": "engine", "qids": qid_tags}

    passing = checks.check_registry_residuals([row], {"engine": output}, registry, "engine")
    corrupted = dict(output)
    corrupted[qid_tags[relation.designated_output]] = generated_output + sp.Symbol("corruption")
    failing = checks.check_registry_residuals(
        [row], {"engine": corrupted}, registry, "engine"
    )

    assert passing.rows[0].status == "ZERO"
    assert failing.rows[0].status == "NONZERO"
    assert failing.rows[0].residual != 0
