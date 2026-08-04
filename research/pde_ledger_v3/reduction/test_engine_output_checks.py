#!/usr/bin/env python3
"""Able-to-fail tests for the engine-output check harness.

No expected physics result is stored here.  The real-output mutations copy or
damage values discovered at test time, and registry targets are generated from
the loaded registry at comparison time.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest
import sympy as sp

import engine_output_checks as checks
from registry_read import load_registry


HERE = Path(__file__).resolve().parent
REAL_S9_OUTPUT = Path("/tmp/s9.txt")
CONFIG = HERE / "checks_S9.yaml"
MODULE = HERE / "engine_output_checks.py"


def real_text() -> str:
    assert REAL_S9_OUTPUT.is_file(), (
        "real S9 output is required; create it once with the Mathematica command "
        "specified in the build request"
    )
    return REAL_S9_OUTPUT.read_text(encoding="utf-8")


def real_raw_and_values() -> tuple[dict[str, str], dict[str, object | checks.Unparsed]]:
    raw = checks.parse_tagged_output(real_text())
    return raw, checks.normalize_tags(raw)


def test_real_s9_end_to_end_cli() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(MODULE),
            "--config",
            str(CONFIG),
            "--output",
            str(REAL_S9_OUTPUT),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "UNPARSED (0):" in completed.stdout
    assert "LIMITATION: triage only" in completed.stdout


def test_control_response_becomes_invariant_when_sole_changed_control_is_copied() -> None:
    raw, values = real_raw_and_values()
    original = checks.check_control_response(values)
    candidate = next(
        row for row in original.responsive if len(row.differing_controls) == 1
    )

    # The compare-time main value is copied; no expected physics value is stored.
    corrupted = dict(raw)
    corrupted[candidate.differing_controls[0]] = raw[candidate.tag]
    observed = checks.check_control_response(checks.normalize_tags(corrupted))

    assert candidate.tag in {row.tag for row in observed.invariant}
    assert candidate.tag not in {row.tag for row in observed.responsive}


def test_mixed_dimension_sum_names_its_summands() -> None:
    raw, _ = real_raw_and_values()
    config = checks.load_config(CONFIG)
    target = "WL_S9_LAGRANGIAN"
    corrupted = dict(raw)
    # omega's definitional dimension differs from the action expression's
    # engine-derived dimension; neither vector is asserted in this test.
    corrupted[target] = f"({raw[target]}) + omega"

    report = checks.check_dimensions(checks.normalize_tags(corrupted), config)
    fired = [issue for issue in report.non_homogeneous if issue.tag == target]

    assert fired
    assert len(fired[0].summands) >= 2
    assert any(term.expression == "omega" for term in fired[0].summands)
    assert all(term.dimension is not None for term in fired[0].summands)


def test_invalid_real_value_is_unparsed_and_never_skipped() -> None:
    raw, _ = real_raw_and_values()
    target = next(iter(raw))
    corrupted = dict(raw)
    corrupted[target] = "this is not valid CAS text !"

    normalized = checks.normalize_tags(corrupted)

    assert isinstance(normalized[target], checks.Unparsed)
    assert sum(isinstance(value, checks.Unparsed) for value in normalized.values()) == 1


def test_duplicate_real_tag_raises() -> None:
    text = real_text()
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
