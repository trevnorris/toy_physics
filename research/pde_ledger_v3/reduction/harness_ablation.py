#!/usr/bin/env python3
"""Deterministic, engine-free ablation battery for engine_output_checks.py.

The instrument accepts the candidate harness, config, and committed outputs.
It mutates only in-memory mappings and prints the observed oracle for each
acceptance item.  It never invokes either symbolic engine.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import re
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence

import sympy as sp


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("candidate_harness", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import harness {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def output_arguments(specs: Sequence[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for spec in specs:
        engine, separator, raw_path = spec.partition("=")
        if not separator or not engine or not raw_path or engine in result:
            raise ValueError(f"invalid --output {spec!r}")
        result[engine] = Path(raw_path)
    return result


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def compact_coverage(coverage: object) -> str:
    return f"{coverage.numerator}/{coverage.denominator} gaps={len(coverage.gaps)}"


def parsed(harness: object, value: object):
    return harness.ParsedValue(harness.value_kind(value), value)


def action_row(config: Mapping[str, Any], package: str, dimension: int, family: str):
    return next(
        row
        for row in config["cross_engine"]
        if row.get("package") == package
        and row.get("dimension") == dimension
        and row.get("family") == family
    )


def run_battery(h, config: dict[str, Any], outputs: dict[str, Mapping[str, object]]) -> list[str]:
    registry = h.load_registry(Path(h.__file__).resolve().parent)
    lines: list[str] = []

    # 1. One complete in-memory report, used as the baseline by later items.
    baseline = h.run_checks(config, outputs)
    rendered = h.format_report(baseline)
    required_sections = [
        "CONTROL_RESPONSE[wl]", "CONTROL_RESPONSE[py]", "DIMENSIONS[wl]",
        "DIMENSIONS[py]", "CROSS_ENGINE_COVERAGE", "REGISTRY_RESIDUAL_COVERAGE",
        "LAGRANGIAN_ACTION_COVERAGE", "EULER_LAGRANGE_ACTION_COVERAGE",
        "LAGRANGIAN_ACTION_VERDICTS", "LAGRANGIAN_ACTION_GAPS",
        "EULER_LAGRANGE_ACTION_VERDICTS", "EULER_LAGRANGE_ACTION_GAPS",
        "IGNORED_LINES[wl]",
    ]
    require(all(section in rendered for section in required_sections), "complete report section absent")
    lines.append(
        "ACCEPTANCE 1 PASS complete_report_sections=" + str(len(required_sections))
        + f" operational_findings={len(baseline.operational_failures)}"
    )

    # 2. Every exact action declaration reaches a comparison verdict.
    action_statuses = [row for row in baseline.cross_engine.rows if row.family]
    require(len(action_statuses) == 26 and all(row.status in h.VERDICTS for row in action_statuses), "action row did not reach verdict")
    lines.append(
        "ACCEPTANCE 2 PASS "
        + " ".join(f"{row.quantity}={row.status}" for row in action_statuses)
    )
    for family in ("lagrangian", "euler_lagrange"):
        representative = next(row for row in action_statuses if row.family == family)
        lines.append(
            f"ACCEPTANCE 2 REPRESENTATIVE {family} {representative.quantity} "
            + " ".join(f"normalized[{engine}]={operand}" for engine, operand in representative.operands)
        )

    # 3. Positional residuals preserve a nonzero RHS, scale, and order through
    # the actual CasRelation/sp.Equality normalization branches.
    residual_row = {
        "quantity": "synthetic_residual", "family": "euler_lagrange",
        "package": "SYNTHETIC", "dimension": 2,
        "left": "L", "right": "R", "cardinality": {"kind": "sequence", "count": 1},
    }
    residual_config = {
        "cells": [{"package": "SYNTHETIC", "dimension": 2, "role": "main"}],
        "control": {
            "required_suffixes": {"left": ["L"], "right": ["R"]},
            "tag_templates": {
                "left": {"main": "{suffix}", "control": "{suffix}"},
                "right": {"main": "{suffix}", "control": "{suffix}"},
            },
        },
        "action_families": {"euler_lagrange": {"left": "L", "right": "R"}},
    }
    x, y = sp.symbols("x y")
    equal = {
        "left": {"L": parsed(h, [h.CasRelation("==", x, 1)])},
        "right": {"R": parsed(h, [sp.Eq(x, 1, evaluate=False)])},
    }
    scaled = {
        "left": equal["left"],
        "right": {"R": parsed(h, [sp.Eq(2 * x, 2, evaluate=False)])},
    }
    reordered_row = {**residual_row, "cardinality": {"kind": "sequence", "count": 2}}
    reordered = {
        "left": {"L": parsed(h, [h.CasRelation("==", x, 1), h.CasRelation("==", y, 2)])},
        "right": {"R": parsed(h, [sp.Eq(y, 2, evaluate=False), sp.Eq(x, 1, evaluate=False)])},
    }
    statuses = (
        h.check_cross_engine([residual_row], equal, config=residual_config).rows[0].status,
        h.check_cross_engine([residual_row], scaled, config=residual_config).rows[0].status,
        h.check_cross_engine([reordered_row], reordered, config=residual_config).rows[0].status,
    )
    require(statuses == ("AGREE", "DISAGREE", "DISAGREE"), "action scale/order oracle failed")
    lines.append(
        f"ACCEPTANCE 3 PASS nonzero_rhs_equations={statuses[0]} "
        f"[x-1]_vs_[2*x-2]={statuses[1]} reordered={statuses[2]} multiset=forbidden"
    )

    # 4. Delete part of one declared cell: every affected row is explicitly
    # PARTIALLY_PAIRED and only the target engine coverage moves.
    controls_by_engine = {report.engine: report for report in baseline.controls}
    damaged_wl = dict(outputs["wl"])
    target_prefix = "WL_S10_XFORM_DIVONLY_D3_"
    deleted = [tag for tag in damaged_wl if tag.startswith(target_prefix) and tag[len(target_prefix):] in config["control"]["required_suffixes"]["wl"]]
    for tag in deleted:
        del damaged_wl[tag]
    damaged_control = h.check_control_response(damaged_wl, config, engine="wl")
    unchanged_py_control = h.check_control_response(outputs["py"], config, engine="py")
    unchanged_py_dimensions = h.check_dimensions(outputs["py"], config, engine="py")
    baseline_py_dimensions = next(report for report in baseline.dimensions if report.engine == "py")
    require(damaged_control.coverage.numerator < controls_by_engine["wl"].coverage.numerator, "target control coverage did not fall")
    require(
        damaged_control.partial
        and len(damaged_control.missing_cells) > len(controls_by_engine["wl"].missing_cells)
        and unchanged_py_control == controls_by_engine["py"]
        and unchanged_py_dimensions == baseline_py_dimensions,
        "partial pairing or other engine isolation failed",
    )
    lines.append(
        f"ACCEPTANCE 4 PASS deleted={len(deleted)} wl={compact_coverage(damaged_control.coverage)} "
        f"partial_rows=[{','.join(row.cell + ':' + row.suffix for row in damaged_control.partial)}] "
        f"missing_cells=[{','.join(damaged_control.missing_cells)}] "
        f"py_control_unchanged={compact_coverage(unchanged_py_control.coverage)} "
        "py_dimension_unchanged=true"
    )

    # 5. Retain one MAIN tag while preserving the full declaration.
    sparse_wl = dict(outputs["wl"])
    main_tags = [tag for tag in sparse_wl if tag.startswith("WL_S10_MAIN_")]
    keep = "WL_S10_MAIN_D3_Q3_DETERMINANT"
    for tag in main_tags:
        if tag != keep:
            del sparse_wl[tag]
    sparse_cross = h.check_cross_engine(config["cross_engine"], {"wl": sparse_wl, "py": outputs["py"]}, config=config, registry=registry)
    require(sparse_cross.coverage.denominator == baseline.cross_engine.coverage.denominator, "cross denominator moved")
    require(sparse_cross.coverage.numerator < baseline.cross_engine.coverage.numerator, "sparse output did not lower coverage")
    require(len(sparse_cross.coverage.gaps) > len(baseline.cross_engine.coverage.gaps), "sparse output did not add named gaps")
    lines.append(
        f"ACCEPTANCE 5 PASS baseline={compact_coverage(baseline.cross_engine.coverage)} "
        f"sparse={compact_coverage(sparse_cross.coverage)} denominator_unchanged=true "
        f"formula={sparse_cross.coverage.formula!r} new_named_gaps={sparse_cross.coverage.gaps[:3]!r}"
    )

    # 6. Rename main, then control package tags.  Cell matching is declaration-driven.
    renamed_main = {tag.replace("WL_S10_MAIN_", "WL_S10_RENAMED_", 1) if tag.startswith("WL_S10_MAIN_") else tag: value for tag, value in outputs["wl"].items()}
    renamed_control = {tag.replace("WL_S10_XFORM_DIVONLY_", "WL_S10_RENAMED_CONTROL_", 1) if tag.startswith("WL_S10_XFORM_DIVONLY_") else tag: value for tag, value in outputs["wl"].items()}
    main_report = h.check_control_response(renamed_main, config, engine="wl")
    control_report = h.check_control_response(renamed_control, config, engine="wl")
    require(
        main_report.missing_cells
        and control_report.missing_cells
        and main_report.coverage.numerator < controls_by_engine["wl"].coverage.numerator
        and control_report.coverage.numerator < controls_by_engine["wl"].coverage.numerator,
        "renamed packages were inferred or retained their matched coverage",
    )
    lines.append(
        f"ACCEPTANCE 6 PASS baseline={compact_coverage(controls_by_engine['wl'].coverage)} "
        f"renamed_main={compact_coverage(main_report.coverage)} main_missing=[{','.join(main_report.missing_cells)}] "
        f"renamed_control={compact_coverage(control_report.coverage)} control_missing=[{','.join(control_report.missing_cells)}]"
    )

    # 7. Empty declarations are explicit operational declaration gaps.
    empty_cross = h.check_cross_engine([], outputs, config=config, registry=registry)
    empty_registry = h.check_registry_residuals([], outputs, registry, config["default_engine"], config=config)
    require(empty_cross.declaration_failure and empty_registry.declaration_failure, "empty declaration passed")
    lines.append(f"ACCEPTANCE 7 PASS cross={empty_cross.declaration_failure!r} registry={empty_registry.declaration_failure!r}")

    # 8. Duplicate identity and empty list/mapping/multiset fixtures.
    duplicate_rows = [
        {"quantity": "a", "wl": "W", "py": "P", "cardinality": {"kind": "scalar"}},
        {"quantity": "b", "wl": "W", "py": "P", "cardinality": {"kind": "scalar"}},
    ]
    duplicate = h.check_cross_engine(duplicate_rows, {"wl": {"W": parsed(h, 1)}, "py": {"P": parsed(h, 1)}})
    empty_statuses = []
    for mode, declaration, value in [
        ("symbolic", {"kind": "sequence", "count": 1}, []),
        ("symbolic", {"kind": "mapping", "entries": 1}, {}),
        ("multiset", {"kind": "multiset", "count": 1}, []),
    ]:
        row = {"quantity": mode + str(declaration), "wl": "W", "py": "P", "mode": mode, "cardinality": declaration}
        empty_statuses.append(h.check_cross_engine([row], {"wl": {"W": parsed(h, value)}, "py": {"P": parsed(h, value)}}).rows[0].status)
    require(duplicate.coverage.denominator == duplicate.coverage.numerator == 1 and duplicate.duplicates, "duplicate inflated coverage")
    require(empty_statuses == ["EMPTY", "EMPTY", "EMPTY"], "empty operand agreed")
    lines.append(f"ACCEPTANCE 8 PASS duplicate_denominator=1 verdicts=1 duplicate_findings={len(duplicate.duplicates)} empty_statuses={empty_statuses}")

    # 9. Named control row moves INVARIANT -> RESPONSIVE.
    small_config = {
        "cells": [{"package": "MAIN", "role": "main"}, {"package": "X", "role": "control"}],
        "control": {"required_suffixes": {"e": ["A"]}, "tag_templates": {"e": {"main": "E_MAIN_{suffix}", "control": "E_{package}_{suffix}"}}},
    }
    invariant_values = {"E_MAIN_A": parsed(h, 1), "E_X_A": parsed(h, 1)}
    responsive_values = {**invariant_values, "E_X_A": parsed(h, 2)}
    invariant = h.check_control_response(invariant_values, small_config, engine="e")
    responsive = h.check_control_response(responsive_values, small_config, engine="e")
    require(invariant.rows[0].status == "INVARIANT" and responsive.rows[0].status == "RESPONSIVE", "control response did not move")
    lines.append("ACCEPTANCE 9 PASS E_MAIN_A INVARIANT->RESPONSIVE")

    # 10. Nested-only boolean exposure.
    boolean_row = {"quantity": "nested_boolean", "wl": "W", "py": "P", "cardinality": {"kind": "mapping", "entries": 1}}
    a = sp.Symbol("a")
    boolean = h.check_cross_engine([boolean_row], {"wl": {"W": parsed(h, {a: True})}, "py": {"P": parsed(h, {a: 1})}})
    boolean_harness_report = replace(
        baseline,
        cross_engine=boolean,
        actions=h.action_family_reports(boolean),
    )
    require(
        not boolean.selected_boolean_rows
        and boolean.tree_boolean_rows == ("nested_boolean",)
        and "tree boolean exposure in nested_boolean" in boolean_harness_report.operational_failures,
        "nested boolean observer failed",
    )
    lines.append("ACCEPTANCE 10 PASS selected_boolean_rows=0 tree_boolean_rows=1 named=nested_boolean operational=true")

    # 11. Remove a source in either direction; only that engine/package table falls.
    table_baseline = {engine: h._build_dimension_tables(engine, outputs[engine], config)[1] for engine in ("wl", "py")}
    source_tags = {
        "py": "PY_S10_MAIN_D2_Q6_INERTIAL_COEFFICIENT1_Q6_DIMENSIONS",
        "wl": "WL_S10_MAIN_D2_Q6_INERTIAL_COEFFICIENT_DIMENSIONS",
    }
    dispositions = []
    for engine, source in source_tags.items():
        damaged = dict(outputs[engine]); del damaged[source]
        tables = h._build_dimension_tables(engine, damaged, config)[1]
        target = next(row for row in tables if row.package == "MAIN")
        other = "wl" if engine == "py" else "py"
        require(not target.assessable and h._build_dimension_tables(other, outputs[other], config)[1] == table_baseline[other], "per-engine source isolation failed")
        dispositions.append(f"{engine}:MAIN=UNASSESSABLE")
    lines.append("ACCEPTANCE 11 PASS " + " ".join(dispositions))

    # 12. A package's valid vector is used locally; removing it never borrows.
    local_config = {
        "default_engine": "e",
        "cells": [
            {"package": "A", "role": "main"},
            {"package": "B", "role": "main"},
            {"package": "C", "role": "main"},
        ],
        "control": {
            "required_suffixes": {"e": ["EXPR"]},
            "tag_templates": {"e": {"main": "E_S10_{package}_{suffix}", "control": "E_S10_{package}_{suffix}"}},
        },
        "primitive_dimensions": {"t": [0, 1, 0]}, "dimensionless": [],
        "dimension_sources": [
            {"engine": "e", "package": "A", "symbols": {"rho": "E_S10_A_DIM"}},
            {"engine": "e", "package": "B", "symbols": {"rho": "E_S10_B_DIM"}},
            {"engine": "e", "package": "C", "symbols": {"rho": "E_S10_C_DIM"}},
        ],
    }
    local_values = h.normalize_tags(
        {
            "E_S10_A_DIM": "{1,0,0}",
            "E_S10_B_DIM": "{2,0,0}",
            "E_S10_C_DIM": "{3,0,0}",
            "E_S10_C_EXPR": "rho + t",
        }
    )
    tables, _, disagreements = h._build_dimension_tables("e", local_values, local_config)
    original_c = tables["C"]
    baseline_c_verdict = next(
        row for row in h.check_dimensions(local_values, local_config, engine="e").statuses
        if row.tag == "E_S10_C_EXPR"
    )
    require(tables["A"]["rho"] != tables["B"]["rho"] and disagreements, "package vector disagreement hidden")
    local_values["E_S10_B_DIM"] = parsed(h, [4, 0, 0])
    changed_tables, _, _ = h._build_dimension_tables("e", local_values, local_config)
    changed_c_verdict = next(
        row for row in h.check_dimensions(local_values, local_config, engine="e").statuses
        if row.tag == "E_S10_C_EXPR"
    )
    require(changed_tables["C"] == original_c and changed_c_verdict == baseline_c_verdict, "changing B changed unrelated package C")
    del local_values["E_S10_B_DIM"]
    tables, reports, _ = h._build_dimension_tables("e", local_values, local_config)
    require(
        "B" not in tables
        and not next(row for row in reports if row.package == "B").assessable
        and tables["C"] == original_c
        and next(
            row for row in h.check_dimensions(local_values, local_config, engine="e").statuses
            if row.tag == "E_S10_C_EXPR"
        ) == baseline_c_verdict,
        "package borrowed another table or changed unrelated package C",
    )
    lines.append(
        f"ACCEPTANCE 12 PASS disagreement={disagreements[0]!r} A_used_own=true "
        f"B_after_removal=UNASSESSABLE C_verdict={baseline_c_verdict.status}_unchanged=true"
    )

    # 13. Bare integers make configured vector and family tables unassessable.
    shape_statuses = []
    for shape in ("vector", "family"):
        shape_config = {
            "primitive_dimensions": {},
            "dimensionless": [],
            "dimension_sources": [
                {
                    "engine": "e",
                    "package": shape.upper(),
                    "symbols": {"rho": {"tag": f"E_S10_{shape.upper()}_DIM", "shape": shape}},
                }
            ],
        }
        shape_values = {f"E_S10_{shape.upper()}_DIM": parsed(h, 2)}
        tables, reports, _ = h._build_dimension_tables("e", shape_values, shape_config)
        report = reports[0]
        require(not report.assessable and report.package not in tables and report.failures, f"bad {shape} source built a table")
        shape_statuses.append(f"{shape}=UNASSESSABLE({report.failures[0]})")
    require(len(shape_statuses) == 2, "bad source shape passed")
    lines.append("ACCEPTANCE 13 PASS " + " ".join(shape_statuses))

    # 14. Corrected probes name their config: S10's x lookup is unassessable;
    # literal 2 and membership are no-comparison.  The unequal add compares.
    probe_config = {
        "default_engine": "e", "cells": [{"package": "MAIN", "role": "main"}],
        "control": {"required_suffixes": {"e": ["P"]}, "tag_templates": {"e": {"main": "E_MAIN_{suffix}", "control": "E_{package}_{suffix}"}}},
        "primitive_dimensions": {"t": [0,1,0]}, "dimensionless": [],
        "dimension_sources": [{"engine": "e", "package": "MAIN", "symbols": {"rho": "E_MAIN_DIM"}}],
    }
    base_probe_values = h.normalize_tags({"E_MAIN_DIM": "{1,0,0}"})
    base_probe_report = h.check_dimensions(base_probe_values, probe_config, engine="e")
    probe_status: dict[str, str] = {}
    probe_deltas: dict[str, tuple[int, int, int, int]] = {}
    for tag, payload in [
        ("E_MAIN_TWO", "2"),
        ("E_MAIN_X", "x"),
        ("E_MAIN_ELEMENT", "Element[x, Reals]"),
    ]:
        one_probe = dict(base_probe_values)
        one_probe.update(h.normalize_tags({tag: payload}))
        observed = h.check_dimensions(one_probe, probe_config, engine="e")
        probe_status[tag] = next(row.status for row in observed.statuses if row.tag == tag)
        probe_deltas[tag] = (
            observed.no_comparison - base_probe_report.no_comparison,
            observed.compared - base_probe_report.compared,
            observed.homogeneous - base_probe_report.homogeneous,
            observed.unassessable - base_probe_report.unassessable,
        )
    add_config = copy.deepcopy(probe_config); add_config["primitive_dimensions"]["x"] = [1,0,0]
    add_values = dict(base_probe_values); add_values["E_MAIN_ADD"] = parsed(h, sp.Symbol("x") + sp.Symbol("t"))
    add_report = h.check_dimensions(add_values, add_config, engine="e")
    require(probe_status["E_MAIN_TWO"] == "no_comparison" and probe_status["E_MAIN_X"] == "unassessable" and probe_status["E_MAIN_ELEMENT"] == "no_comparison", "corrected probe partition failed")
    require(
        probe_deltas["E_MAIN_TWO"] == (1, 0, 0, 0)
        and probe_deltas["E_MAIN_ELEMENT"] == (1, 0, 0, 0)
        and probe_deltas["E_MAIN_X"] == (0, 0, 0, 1),
        "probe moved an incorrect dimension bucket",
    )
    require(
        add_report.compared == base_probe_report.compared + 1
        and any(issue.tag == "E_MAIN_ADD" for issue in add_report.non_homogeneous),
        "unequal add did not compare",
    )
    lines.append(
        f"ACCEPTANCE 14 PASS config=S10 2={probe_status['E_MAIN_TWO']} delta={probe_deltas['E_MAIN_TWO']} "
        f"x={probe_status['E_MAIN_X']} delta={probe_deltas['E_MAIN_X']} "
        f"Element={probe_status['E_MAIN_ELEMENT']} delta={probe_deltas['E_MAIN_ELEMENT']} "
        f"unequal_add=compared sites={dict(add_report.proposition_sites)}"
    )

    # 15. Mapping values are walked.
    mapping_values = dict(base_probe_values)
    mapping_values["E_MAIN_MAPPING"] = parsed(h, {"k": sp.Symbol("x") + sp.Symbol("t")})
    mapping_report = h.check_dimensions(mapping_values, add_config, engine="e")
    mapping_status = next(row.status for row in mapping_report.statuses if row.tag == "E_MAIN_MAPPING")
    arithmetic = sum([mapping_report.compared, mapping_report.no_comparison, mapping_report.not_applicable, mapping_report.unwalked, mapping_report.unassessable, mapping_report.unparsed])
    require(
        mapping_status == "compared"
        and mapping_report.compared == base_probe_report.compared + 1
        and any(issue.tag == "E_MAIN_MAPPING" for issue in mapping_report.non_homogeneous)
        and mapping_report.unwalked == 0
        and arithmetic == mapping_report.total_tags,
        "mapping walk/partition failed",
    )
    lines.append(f"ACCEPTANCE 15 PASS mapping=compared non_homogeneous={len(mapping_report.non_homogeneous)} unwalked=0 partition={arithmetic}/{mapping_report.total_tags}")

    # 16. B9 structural fixtures.
    fixtures = [
        ("derivative_arity", "Derivative[1,0][u][x,t]", "Derivative[1,0,0][u][x,t,y]"),
        ("derivative_argument_order", "Derivative[1,0][u][x,t]", "Derivative[1,0][u][t,x]"),
        ("membership", "Element[x, Integers]", "Element[x, Reals]"),
        ("quoted_member", '<|"k" -> "v"|>', '<|"k" -> "w"|>'),
        ("marker", "markerA[{1}]", "markerB[{1}]"),
        ("piecewise_branch", "Piecewise[{{1,x>0}},0]", "Piecewise[{{2,x>0}},0]"),
        ("piecewise_condition", "Piecewise[{{1,x>0}},0]", "Piecewise[{{1,x>=0}},0]"),
        ("inequality_operand", "Inequality[x,Less,y]", "Inequality[x,Less,z]"),
        ("inequality_operator", "Inequality[x,Less,y]", "Inequality[x,LessEqual,y]"),
        ("real", "60.62794", "60.62795"),
    ]
    fixture_results = []
    for name, left, right in fixtures:
        lhs, rhs = h.normalize_mathematica(left), h.normalize_mathematica(right)
        require(
            not isinstance(lhs, h.Unparsed)
            and not isinstance(rhs, h.Unparsed)
            and not h.symbolic_equal(
                h._canonicalize_derivatives(lhs), h._canonicalize_derivatives(rhs)
            ),
            f"{name} structure erased",
        )
        fixture_results.append(name + "=observable")
    lines.append("ACCEPTANCE 16 PASS " + " ".join(fixture_results))

    # 17. Diagnostics/stray lines remain outside neighboring payloads.
    diagnostic = h.parse_tagged_output("WL_S10_A: 1\nSolve::svars: diagnostic\nWL_S10_B: 2\n")
    stray = h.parse_tagged_output("WL_S10_A: 1\nstray line\nWL_S10_B: 2\n")
    expected_neighbors = {"WL_S10_A": "1", "WL_S10_B": "2"}
    require(
        diagnostic == expected_neighbors
        and stray == expected_neighbors
        and diagnostic.ignored_lines == ("Solve::svars: diagnostic",)
        and stray.ignored_lines == ("stray line",)
        and "CROSS_ENGINE_COVERAGE" in rendered,
        "ignored line absorption or report truncation",
    )
    lines.append(
        f"ACCEPTANCE 17 PASS diagnostic_ignored={diagnostic.ignored_lines!r} "
        f"stray_ignored={stray.ignored_lines!r} neighbors={dict(diagnostic)!r} rest_of_report=present"
    )

    # 18. Derivative-monomial support changes under a form change but is
    # invariant under a coefficient-only change; no config label is the oracle.
    main = action_row(config, "MAIN", 3, "lagrangian")
    naming = h._build_naming_rule(config, registry)

    def derivative_support(engine: str, row: Mapping[str, Any]) -> frozenset[tuple[int, ...]]:
        source = outputs[engine][row[engine]]
        require(isinstance(source, h.ParsedValue), f"{engine}:{row['quantity']} is unparsed")
        value, _ = h._map_tree(source.value, naming.engine_maps.get(engine, {}))
        value = h._canonicalize_derivatives(value)
        require(isinstance(value, sp.Basic), f"{engine}:{row['quantity']} is not scalar")
        derivatives = sorted(value.atoms(h.CanonicalDerivative), key=str)
        require(bool(derivatives), f"{engine}:{row['quantity']} has no derivatives")
        polynomial = sp.Poly(sp.expand(value), *derivatives)
        return frozenset(monomial for monomial, coefficient in polynomial.terms() if coefficient != 0)

    baseline_support = {
        engine: derivative_support(engine, main) for engine in ("wl", "py")
    }
    form_results = []
    for package in ("XFORM_FULLGRAD", "XFORM_DIVONLY"):
        source = action_row(config, package, 3, "lagrangian")
        changed = {
            engine: derivative_support(engine, source) != baseline_support[engine]
            for engine in ("wl", "py")
        }
        require(any(changed.values()), "form control did not change derivative support in either engine")
        form_results.append(f"{package}=FORM_CHANGE(wl={str(changed['wl']).lower()},py={str(changed['py']).lower()})")
    coefficient_results = []
    for package in ("XFORM_SIGNFLIP", "XFORM_ANISO", "XCOEF_SCALE"):
        source = action_row(config, package, 3, "lagrangian")
        changed = {
            engine: derivative_support(engine, source) != baseline_support[engine]
            for engine in ("wl", "py")
        }
        require(not any(changed.values()), "coefficient control changed derivative support")
        coefficient_results.append(f"{package}=COEFFICIENT_ONLY(wl=false,py=false)")
    lines.append(
        "ACCEPTANCE 18 PASS oracle=derivative_monomial_support selected="
        + " ".join(form_results + coefficient_results)
    )

    return lines


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harness", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, action="append", metavar="ENGINE=PATH")
    args = parser.parse_args(argv)
    try:
        harness = load_module(args.harness.resolve())
        config = harness.load_config(args.config)
        paths = output_arguments(args.output)
        outputs = {engine: harness.load_output(path, syntax=engine if engine in {"wl", "py"} else None)[1] for engine, path in paths.items()}
        for line in run_battery(harness, config, outputs):
            print(line)
        return 0
    except Exception as exc:
        print(f"ABLATION_FAILURE: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
