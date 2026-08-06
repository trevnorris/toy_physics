#!/usr/bin/env python3
"""Able-to-fail tests for the declaration-oriented engine-output harness.

All real-output fixtures are committed files.  This suite never invokes either
engine and in particular never launches Mathematica.
"""

from __future__ import annotations

import copy
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import pytest
import sympy as sp
import yaml

import engine_output_checks as checks
from registry_read import load_registry


HERE = Path(__file__).resolve().parent
LEDGER = HERE.parent
S10_CONFIG = HERE / "checks_S10.yaml"
S9_CONFIG = HERE / "checks_S9.yaml"
REAL = {
    "S10": {
        "wl": LEDGER / "mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out",
        "py": LEDGER / "scripts/out/S10_brane_mode_spectrum_sympy_audit.out",
    },
    "S9": {
        "wl": LEDGER / "mathematica/out/S9_light_requires_shear_mathematica_audit.out",
        "py": LEDGER / "scripts/out/S9_light_requires_shear_sympy_audit.out",
    },
}


@lru_cache(maxsize=4)
def real_values(stage: str, engine: str) -> checks.NormalizedRecords:
    return checks.load_output(REAL[stage][engine], syntax=engine)[1]


def parsed(value: object) -> checks.ParsedValue:
    return checks.ParsedValue(checks.value_kind(value), value)


def cross_row(**updates: object) -> dict[str, object]:
    row: dict[str, object] = {
        "quantity": "row",
        "left": "LEFT",
        "right": "RIGHT",
        "cardinality": {"kind": "scalar"},
    }
    row.update(updates)
    return row


def control_config() -> dict[str, object]:
    return {
        "cells": [
            {"package": "MAIN", "dimension": 3, "role": "main", "stiffness_form": "curl"},
            {"package": "FULL", "dimension": 3, "role": "control", "stiffness_form": "fullgrad"},
            {"package": "DIV", "dimension": 3, "role": "control", "stiffness_form": "divonly"},
        ],
        "control": {
            "required_suffixes": {"e": ["A", "B"]},
            "tag_templates": {
                "e": {
                    "main": "E_{package}_D{dimension}_{suffix}",
                    "control": "E_{package}_D{dimension}_{suffix}",
                }
            },
        },
    }


def action_config(*, dimension: int = 2) -> dict[str, object]:
    return {
        "cells": [{"package": "SYNTHETIC", "dimension": dimension, "role": "main"}],
        "control": {
            "required_suffixes": {"left": ["LEFT"], "right": ["RIGHT"]},
            "tag_templates": {
                "left": {"main": "{suffix}", "control": "{suffix}"},
                "right": {"main": "{suffix}", "control": "{suffix}"},
            },
        },
        "action_families": {
            "lagrangian": {"left": "LEFT", "right": "RIGHT"},
            "euler_lagrange": {"left": "LEFT", "right": "RIGHT"},
        },
    }


def dimension_config(*, primitive_x: bool = False) -> dict[str, object]:
    primitives = {"t": [0, 1, 0]}
    if primitive_x:
        primitives["x"] = [1, 0, 0]
    return {
        "default_engine": "e",
        "cells": [{"package": "MAIN", "role": "main"}],
        "control": {
            "required_suffixes": {"e": ["PROBE"]},
            "tag_templates": {"e": {"main": "E_MAIN_{suffix}", "control": "E_{package}_{suffix}"}},
        },
        "primitive_dimensions": primitives,
        "dimensionless": [],
        "dimension_sources": [
            {"engine": "e", "package": "MAIN", "symbols": {"rho": "E_MAIN_DIM"}}
        ],
    }


def write_cli_fixture(
    tmp_path: Path, *, malformed_py_operand: bool = False
) -> tuple[Path, dict[str, Path]]:
    """Write an engine-free CLI fixture with no operational gaps."""
    config = {
        "default_engine": "wl",
        "cells": [
            {"package": "MAIN", "role": "main"},
            {"package": "CONTROL", "role": "control"},
        ],
        "control": {
            "required_suffixes": {"wl": ["VALUE"], "py": ["VALUE"]},
            "tag_templates": {
                "wl": {"main": "WL_S10_{package}_{suffix}", "control": "WL_S10_{package}_{suffix}"},
                "py": {"main": "PY_S10_{package}_{suffix}", "control": "PY_S10_{package}_{suffix}"},
            },
        },
        "primitive_dimensions": {},
        "dimensionless": [],
        "dimension_sources": [
            {
                "engine": engine,
                "package": package,
                "symbols": {
                    "rho": f"{prefix}_S10_{package}_DIM",
                    "sigma": f"{prefix}_S10_{package}_DIM",
                },
            }
            for engine, prefix in (("wl", "WL"), ("py", "PY"))
            for package in ("MAIN", "CONTROL")
        ],
        "cross_engine": [
            {
                "quantity": "physics_disagreement",
                "wl": "WL_S10_MAIN_VALUE",
                "py": "PY_S10_MAIN_VALUE",
                "cardinality": {"kind": "scalar"},
            }
        ],
        "registry_residual": [
            {
                "relation_id": "R4",
                "engine": "wl",
                "qids": {
                    "Q.brane.rho_br": "WL_S10_REGISTRY_RHO",
                    "Q.brane.mu_R": "WL_S10_REGISTRY_MU",
                    "Q.brane.c_gamma": "WL_S10_REGISTRY_C",
                },
            }
        ],
    }
    config_path = tmp_path / "checks.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    payloads = {
        "wl": {
            "WL_S10_MAIN_VALUE": '"left"',
            "WL_S10_CONTROL_VALUE": '"left"',
            "WL_S10_MAIN_DIM": "{1,0,0}",
            "WL_S10_CONTROL_DIM": "{1,0,0}",
            "WL_S10_MAIN_DIMENSION_CHECK": "rho + sigma",
            "WL_S10_CONTROL_DIMENSION_CHECK": "rho + sigma",
            "WL_S10_REGISTRY_RHO": "1",
            "WL_S10_REGISTRY_MU": "1",
            "WL_S10_REGISTRY_C": "1",
        },
        "py": {
            "PY_S10_MAIN_VALUE": "this is not valid CAS text !" if malformed_py_operand else '"right"',
            "PY_S10_CONTROL_VALUE": '"right"',
            "PY_S10_MAIN_DIM": "(1,0,0)",
            "PY_S10_CONTROL_DIM": "(1,0,0)",
            "PY_S10_MAIN_DIMENSION_CHECK": "rho + sigma",
            "PY_S10_CONTROL_DIMENSION_CHECK": "rho + sigma",
        },
    }
    paths: dict[str, Path] = {}
    for engine, tagged in payloads.items():
        path = tmp_path / f"{engine}.out"
        path.write_text("\n".join(f"{tag}: {value}" for tag, value in tagged.items()) + "\n", encoding="utf-8")
        paths[engine] = path
    return config_path, paths


def run_harness_cli(config: Path, outputs: dict[str, Path]) -> subprocess.CompletedProcess[str]:
    arguments = [sys.executable, str(HERE / "engine_output_checks.py"), "--config", str(config)]
    for engine, path in outputs.items():
        arguments.extend(["--output", f"{engine}={path}"])
    return subprocess.run(
        arguments,
        cwd=LEDGER,
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )


def test_parser_ignores_diagnostics_and_stray_lines_without_absorbing_payload() -> None:
    raw = checks.parse_tagged_output(
        "WL_S10_A: 1\nSolve::svars: diagnostic\nstray text\nWL_S10_B: 2\n"
    )
    assert raw == {"WL_S10_A": "1", "WL_S10_B": "2"}
    assert raw.ignored_lines == ("Solve::svars: diagnostic", "stray text")


def test_duplicate_valid_tag_is_retained_as_operational_metadata() -> None:
    raw = checks.parse_tagged_output("WL_S10_A: 1\nWL_S10_A: 2\nWL_S10_B: 3\n")
    assert raw["WL_S10_A"] == "1"
    assert raw["WL_S10_B"] == "3"
    assert raw.duplicate_tags == ("WL_S10_A",)


@pytest.mark.parametrize(
    ("equal_left", "equal_right", "unequal"),
    [
        (
            "Derivative[1,0][u][x,t]",
            "Derivative[1,0][u][x,t]",
            "Derivative[0,1][u][x,t]",
        ),
        ("Element[x, Integers]", "Element[x, Integers]", "Element[x, Reals]"),
        ('<|"key" -> {"member", 2}|>', '<|"key" -> {"member", 2}|>', '<|"key" -> {"other", 2}|>'),
        ("markerHead[{1,2}]", "markerHead[{1,2}]", "otherHead[{1,2}]"),
        ("Piecewise[{{1, x > 0}}, 0]", "Piecewise[{{1, x > 0}}, 0]", "Piecewise[{{2, x > 0}}, 0]"),
        ("Inequality[x, Less, y]", "Inequality[x, Less, y]", "Inequality[x, LessEqual, y]"),
        ("60.62794", "60.62794", "60.62795"),
    ],
)
def test_live_wl_constructs_preserve_structural_equality(
    equal_left: str, equal_right: str, unequal: str
) -> None:
    left = checks.normalize_mathematica(equal_left)
    right = checks.normalize_mathematica(equal_right)
    changed = checks.normalize_mathematica(unequal)
    assert not isinstance(left, checks.Unparsed)
    assert not isinstance(right, checks.Unparsed)
    assert not isinstance(changed, checks.Unparsed)
    assert checks.symbolic_equal(left, right)
    assert not checks.symbolic_equal(left, changed)


def test_derivative_argument_order_is_identity() -> None:
    left = checks.normalize_mathematica("Derivative[1,0][u][x,t]")
    right = checks.normalize_mathematica("Derivative[1,0][u][t,x]")
    assert not checks.symbolic_equal(left, right)


@pytest.mark.parametrize(
    ("wl", "py", "expected"),
    [
        ("Derivative[1,0][u][x,t]", "Derivative(u(t,x),x)", "AGREE"),
        ("Derivative[2,0][u][x,t]", "Derivative(u(t,x),(x,2))", "AGREE"),
        ("Derivative[1,0][u][x,t]", "Derivative(v(t,x),x)", "DISAGREE"),
        ("Derivative[1,0][u][x,t]", "Derivative(u(t,x),t)", "DISAGREE"),
        ("Derivative[2,0][u][x,t]", "Derivative(u(t,x),x)", "DISAGREE"),
        ("Derivative[1,0][u][x,t]", "Derivative(u(t,x,y),x)", "DISAGREE"),
    ],
)
def test_cross_engine_derivative_identity_uses_function_arguments_variable_and_order(
    wl: str, py: str, expected: str
) -> None:
    row = {"quantity": "derivative", "wl": "W", "py": "P", "cardinality": {"kind": "scalar"}}
    outputs = {
        "wl": {"W": checks.normalize_mathematica(wl)},
        "py": {"P": checks.normalize_sympy(py)},
    }
    assert checks.check_cross_engine([row], outputs).rows[0].status == expected


def test_all_committed_wl_constructs_now_parse() -> None:
    assert not [value for value in real_values("S10", "wl").values() if isinstance(value, checks.Unparsed)]
    assert not [value for value in real_values("S9", "wl").values() if isinstance(value, checks.Unparsed)]


def test_ordered_action_residual_detects_scale_and_reordering() -> None:
    row = cross_row(
        family="euler_lagrange",
        package="SYNTHETIC",
        dimension=2,
        cardinality={"kind": "sequence", "count": 2},
    )
    x, y = sp.symbols("x y")
    equal = {
        "left": {"LEFT": parsed([checks.CasRelation("==", x, 1), checks.CasRelation("==", y, 2)])},
        "right": {"RIGHT": parsed([sp.Eq(x, 1, evaluate=False), sp.Eq(y, 2, evaluate=False)])},
    }
    scaled = copy.deepcopy(equal)
    scaled["right"]["RIGHT"] = parsed([sp.Eq(2 * x, 2, evaluate=False), sp.Eq(y, 2, evaluate=False)])
    reordered = copy.deepcopy(equal)
    reordered["right"]["RIGHT"] = parsed([sp.Eq(y, 2, evaluate=False), sp.Eq(x, 1, evaluate=False)])
    config = action_config()
    assert checks.check_cross_engine([row], equal, config=config).rows[0].status == "AGREE"
    assert checks.check_cross_engine([row], scaled, config=config).rows[0].status == "DISAGREE"
    assert checks.check_cross_engine([row], reordered, config=config).rows[0].status == "DISAGREE"


def test_action_equation_branch_retains_nonzero_right_hand_side() -> None:
    x = sp.Symbol("x")
    row = cross_row(
        family="euler_lagrange",
        package="SYNTHETIC",
        dimension=1,
        cardinality={"kind": "sequence", "count": 1},
    )
    left = {"left": {"LEFT": parsed([checks.CasRelation("==", x, 3)])}}
    equal = {**left, "right": {"RIGHT": parsed([sp.Eq(x, 3, evaluate=False)])}}
    changed_rhs = {**left, "right": {"RIGHT": parsed([sp.Eq(x, 4, evaluate=False)])}}
    config = action_config(dimension=1)
    assert checks.check_cross_engine([row], equal, config=config).rows[0].status == "AGREE"
    assert checks.check_cross_engine([row], changed_rhs, config=config).rows[0].status == "DISAGREE"


def test_action_row_forbids_multiset_mode() -> None:
    row = cross_row(
        family="euler_lagrange",
        package="SYNTHETIC",
        dimension=2,
        mode="multiset",
        cardinality={"kind": "multiset", "count": 2},
    )
    outputs = {"left": {"LEFT": parsed([1, 2])}, "right": {"RIGHT": parsed([1, 2])}}
    observed = checks.check_cross_engine([row], outputs, config=action_config()).rows[0]
    assert observed.status == "DECLARATION_ERROR"
    assert "multiset" in observed.detail


@pytest.mark.parametrize(
    ("kind", "value", "declaration", "mode"),
    [
        ("list", [], {"kind": "sequence", "count": 1}, "symbolic"),
        ("mapping", {}, {"kind": "mapping", "entries": 1}, "symbolic"),
        ("multiset", [], {"kind": "multiset", "count": 1}, "multiset"),
    ],
)
def test_empty_operands_never_agree(
    kind: str, value: object, declaration: dict[str, object], mode: str
) -> None:
    row = cross_row(cardinality=declaration, mode=mode)
    outputs = {"left": {"LEFT": parsed(value)}, "right": {"RIGHT": parsed(value)}}
    observed = checks.check_cross_engine([row], outputs).rows[0]
    assert observed.status == "EMPTY", kind


def test_duplicate_pair_does_not_increase_coverage_denominator() -> None:
    first = cross_row(quantity="first")
    second = cross_row(quantity="renamed")
    outputs = {"left": {"LEFT": parsed(1)}, "right": {"RIGHT": parsed(1)}}
    report = checks.check_cross_engine([first, second], outputs)
    assert report.coverage == checks.Coverage(1, 1, (), report.coverage.formula)
    assert report.duplicates


def test_missing_operand_reduces_fixed_cross_engine_coverage() -> None:
    rows = [cross_row(quantity="one"), cross_row(quantity="two", left="LEFT2", right="RIGHT2")]
    outputs = {"left": {"LEFT": parsed(1), "LEFT2": parsed(2)}, "right": {"RIGHT": parsed(1)}}
    report = checks.check_cross_engine(rows, outputs)
    assert (report.coverage.numerator, report.coverage.denominator) == (1, 2)
    assert report.coverage.gaps == ("two:MISSING",)


def test_nested_boolean_observer_catches_deferred_truthiness_exposure() -> None:
    row = cross_row(cardinality={"kind": "mapping", "entries": 1})
    a = sp.Symbol("a")
    outputs = {"left": {"LEFT": parsed({a: True})}, "right": {"RIGHT": parsed({a: 1})}}
    report = checks.check_cross_engine([row], outputs)
    assert checks.symbolic_equal({a: True}, {a: 1})  # deferred equality semantics
    assert report.selected_boolean_rows == ()
    assert report.tree_boolean_rows == ("row",)


def test_declared_naming_rule_maps_registry_snake_case_symmetrically() -> None:
    config = checks.load_config(S10_CONFIG)
    registry = load_registry(HERE)
    row = {
        "quantity": "row",
        "wl": "LEFT",
        "py": "RIGHT",
        "cardinality": {"kind": "scalar"},
    }
    outputs = {
        "left": {"LEFT": parsed(sp.Symbol("rhoBr") + sp.Symbol("muR"))},
        "right": {"RIGHT": parsed(sp.Symbol("rho_br") + sp.Symbol("mu_R"))},
    }
    # Reuse the configured WL/PY engine names so the declared styles apply.
    outputs = {"wl": {"LEFT": outputs["left"]["LEFT"]}, "py": {"RIGHT": outputs["right"]["RIGHT"]}}
    observed = checks.check_cross_engine([row], outputs, config=config, registry=registry).rows[0]
    assert observed.status == "AGREE"
    assert "wl:rhoBr->rho_br" in observed.naming_applied


def test_omega_algebraic_aliases_are_not_in_declared_naming_layer() -> None:
    config = checks.load_config(S9_CONFIG)
    registry = load_registry(HERE)
    row = {"quantity": "omega", "wl": "W", "py": "P", "cardinality": {"kind": "scalar"}}
    outputs = {"wl": {"W": parsed(sp.Symbol("omegaSquared"))}, "py": {"P": parsed(sp.Symbol("omega2"))}}
    report = checks.check_cross_engine([row], outputs, config=config, registry=registry)
    assert report.rows[0].status == "DISAGREE"
    assert report.rows[0].legacy_algebraic_status == "AGREE"
    assert report.naming_changed_rows == ("omega",)


def test_undeclared_mechanical_spelling_mismatch_is_operational() -> None:
    config = {
        "symbol_naming": {
            "rule": "registry_snake_case_to_lower_camel",
            "engine_styles": {"wl": "lower_camel", "py": "canonical"},
            "exceptions": [],
        }
    }
    row = {"quantity": "undeclared", "wl": "W", "py": "P", "cardinality": {"kind": "scalar"}}
    outputs = {"wl": {"W": parsed(sp.Symbol("notDeclared"))}, "py": {"P": parsed(sp.Symbol("not_declared"))}}
    observed = checks.check_cross_engine([row], outputs, config=config, registry=load_registry(HERE))
    assert observed.rows[0].status == "NAMING_MISMATCH"
    assert observed.coverage == checks.Coverage(0, 1, ("undeclared:NAMING_MISMATCH",), observed.coverage.formula)
    assert observed.rows[0].undeclared_spellings == ("not_declared<->notDeclared",)


def test_bijection_detector_does_not_hide_coefficient_changes() -> None:
    a, b = sp.symbols("a b")
    row = cross_row()
    outputs = {"left": {"LEFT": parsed(2 * a + 1)}, "right": {"RIGHT": parsed(b + 1)}}
    observed = checks.check_cross_engine([row], outputs).rows[0]
    assert observed.status == "DISAGREE"
    assert observed.undeclared_spellings == ()


def test_s10_declares_only_verified_dimension_scale_and_gradient_spellings() -> None:
    config = checks.load_config(S10_CONFIG)
    names = {
        "main_d3_q6_energy_density_dimension",
        "xcoef_scale_d3_q6_stiffness_coefficients",
        "main_d3_q2_downstream_route",
        "main_d3_q7_stiffness",
        "main_d3_root2_q5_scale_ratio",
        "xform_fullgrad_d3_q7_stiffness",
    }
    rows = [row for row in config["cross_engine"] if row["quantity"] in names]
    outputs = {"wl": real_values("S10", "wl"), "py": real_values("S10", "py")}
    observed = checks.check_cross_engine(rows, outputs, config=config, registry=load_registry(HERE))
    by_name = {row.quantity: row for row in observed.rows}
    assert by_name["main_d3_q6_energy_density_dimension"].status == "AGREE"
    assert "wl:braneDimension->D" in by_name["main_d3_q6_energy_density_dimension"].naming_applied
    assert by_name["xcoef_scale_d3_q6_stiffness_coefficients"].status == "AGREE"
    assert "py:s->coefficientScale" in by_name["xcoef_scale_d3_q6_stiffness_coefficients"].naming_applied
    assert by_name["main_d3_root2_q5_scale_ratio"].status == "AGREE"
    assert set(by_name["main_d3_root2_q5_scale_ratio"].naming_applied) == {
        "py:lambdaScale->lambda_scale",
        "wl:lambdaScale->lambda_scale",
    }
    assert by_name["main_d3_q2_downstream_route"].status == "NAMING_MISMATCH"
    assert by_name["main_d3_q2_downstream_route"].undeclared_spellings == (
        "M_B<->quadraticFormRoute",
    )
    assert by_name["main_d3_q7_stiffness"].status == "AGREE"
    assert set(by_name["main_d3_q7_stiffness"].naming_applied) == {
        "wl:g1x2->g12",
        "wl:g1x3->g13",
        "wl:g2x1->g21",
        "wl:g2x3->g23",
        "wl:g3x1->g31",
        "wl:g3x2->g32",
    }
    assert by_name["main_d3_q7_stiffness"].undeclared_spellings == ()
    assert by_name["xform_fullgrad_d3_q7_stiffness"].status == "AGREE"
    assert set(by_name["xform_fullgrad_d3_q7_stiffness"].naming_applied) == {
        "wl:g1x1->g11",
        "wl:g1x2->g12",
        "wl:g1x3->g13",
        "wl:g2x1->g21",
        "wl:g2x2->g22",
        "wl:g2x3->g23",
        "wl:g3x1->g31",
        "wl:g3x2->g32",
        "wl:g3x3->g33",
    }
    assert by_name["xform_fullgrad_d3_q7_stiffness"].undeclared_spellings == ()


def test_partial_control_pair_is_named_and_reduces_control_coverage() -> None:
    config = control_config()
    values = checks.normalize_tags(
        {
            "E_MAIN_D3_A": "1", "E_FULL_D3_A": "1", "E_DIV_D3_A": "1",
            "E_MAIN_D3_B": "2", "E_FULL_D3_B": "2", "E_DIV_D3_B": "2",
        }
    )
    baseline = checks.check_control_response(values, config, engine="e")
    damaged = dict(values)
    del damaged["E_DIV_D3_A"]
    observed = checks.check_control_response(damaged, config, engine="e")
    assert baseline.coverage.numerator == baseline.coverage.denominator == 2
    assert observed.coverage.numerator == 1
    assert observed.partial[0].missing_tags == ("E_DIV_D3_A",)


def test_control_row_moves_invariant_to_responsive() -> None:
    config = control_config()
    values = checks.normalize_tags(
        {
            "E_MAIN_D3_A": "1", "E_FULL_D3_A": "1", "E_DIV_D3_A": "1",
            "E_MAIN_D3_B": "1", "E_FULL_D3_B": "1", "E_DIV_D3_B": "1",
        }
    )
    baseline = checks.check_control_response(values, config, engine="e")
    changed = dict(values)
    changed["E_DIV_D3_A"] = parsed(2)
    observed = checks.check_control_response(changed, config, engine="e")
    assert len(baseline.invariant) == 2
    assert any(row.suffix == "A" for row in observed.responsive)


def test_dimension_shapes_include_live_py_nested_family() -> None:
    vector = checks.normalize_sympy("(-D, 0, 1)")
    family = checks.normalize_mathematica("{{-D,0,1},{-D,0,1}}")
    nested = checks.normalize_sympy("(((-D, 0, 1),),)")
    assert not isinstance(vector, checks.Unparsed)
    assert not isinstance(family, checks.Unparsed)
    assert not isinstance(nested, checks.Unparsed)
    assert checks._dimension_source_value(vector.value, "vector", "v") == (-sp.Symbol("D"), 0, 1)
    assert checks._dimension_source_value(family.value, "family", "f") == (-sp.Symbol("D"), 0, 1)
    assert checks._dimension_source_value(nested.value, "nested_family", "n") == (-sp.Symbol("D"), 0, 1)


def test_removed_py_dimension_source_invalidates_only_py_package() -> None:
    config = checks.load_config(S10_CONFIG)
    wl = dict(real_values("S10", "wl"))
    py = dict(real_values("S10", "py"))
    before_wl = checks._build_dimension_tables("wl", wl, config)[1]
    target = "PY_S10_MAIN_D2_Q6_INERTIAL_COEFFICIENT1_Q6_DIMENSIONS"
    del py[target]
    after_py = checks._build_dimension_tables("py", py, config)[1]
    after_wl = checks._build_dimension_tables("wl", wl, config)[1]
    assert not next(row for row in after_py if row.package == "MAIN").assessable
    assert after_wl == before_wl


def test_package_uses_own_vector_and_never_borrows() -> None:
    config = dimension_config()
    config["cells"].append({"package": "OTHER", "role": "main"})
    config["dimension_sources"].append(
        {"engine": "e", "package": "OTHER", "symbols": {"rho": "E_OTHER_DIM"}}
    )
    values = checks.normalize_tags({"E_MAIN_DIM": "{1,0,0}", "E_OTHER_DIM": "{2,0,0}"})
    tables, reports, disagreements = checks._build_dimension_tables("e", values, config)
    assert tables["MAIN"]["rho"] != tables["OTHER"]["rho"]
    assert disagreements
    del values["E_OTHER_DIM"]
    tables, reports, _ = checks._build_dimension_tables("e", values, config)
    assert "OTHER" not in tables
    assert not next(row for row in reports if row.package == "OTHER").assessable


@pytest.mark.parametrize("shape", ["vector", "family"])
def test_bad_dimension_source_shape_is_operational(shape: str) -> None:
    source = checks.ParsedValue(checks.ValueKind.INTEGER, sp.Integer(2))
    with pytest.raises(ValueError, match="expected|family"):
        checks._dimension_source_value(source.value, shape, "source")


def test_dimension_probe_partition_uses_corrected_config_semantics() -> None:
    config = dimension_config(primitive_x=False)
    values = checks.normalize_tags(
        {
            "E_MAIN_DIM": "{1,0,0}",
            "E_MAIN_TWO": "2",
            "E_MAIN_X": "x",
            "E_MAIN_ELEMENT": "Element[x, Reals]",
        }
    )
    report = checks.check_dimensions(values, config, engine="e")
    status = {row.tag: row.status for row in report.statuses}
    assert status["E_MAIN_TWO"] == "no_comparison"
    assert status["E_MAIN_X"] == "unassessable"  # S10-style table has no x primitive.
    assert status["E_MAIN_ELEMENT"] == "no_comparison"
    assert report.total_tags == sum(
        [report.compared, report.no_comparison, report.not_applicable, report.unwalked, report.unassessable, report.unparsed]
    )


def test_dimension_declaration_cannot_go_quiet_at_zero_comparisons() -> None:
    config = dimension_config(primitive_x=True)
    values = checks.normalize_tags({"E_MAIN_DIM": "{1,0,0}", "E_MAIN_PROBE": "x+t"})
    baseline = checks.check_dimensions(values, config, engine="e")
    assert baseline.compared == 1
    assert baseline.coverage == checks.Coverage(1, 1, (), baseline.coverage.formula)

    no_cells = copy.deepcopy(config)
    no_cells["cells"] = []
    silent = checks.check_dimensions(values, no_cells, engine="e")
    assert silent.compared == 0
    assert silent.coverage.denominator == 1
    assert silent.coverage.gaps == ("MAIN",)
    harness = checks.run_checks(
        {**no_cells, "cross_engine": [], "registry_residual": []}, {"e": values}
    )
    assert "e: zero dimension comparisons against non-empty declaration" in harness.operational_failures

    no_sources = copy.deepcopy(config)
    del no_sources["dimension_sources"]
    absent = checks.check_dimensions(values, no_sources, engine="e")
    assert absent.compared == 0
    assert absent.declaration_failure == "e: dimension_sources declaration is absent or empty"


def test_evaluable_addition_in_mapping_is_compared_not_unwalked() -> None:
    config = dimension_config(primitive_x=True)
    values = checks.normalize_tags({"E_MAIN_DIM": "{1,0,0}"})
    values["E_MAIN_MAPPING"] = parsed({"key": sp.Symbol("x") + sp.Symbol("t")})
    report = checks.check_dimensions(values, config, engine="e")
    row = next(row for row in report.statuses if row.tag == "E_MAIN_MAPPING")
    assert row.status == "compared"
    assert any(issue.tag == "E_MAIN_MAPPING" for issue in report.non_homogeneous)
    assert report.unwalked == 0
    assert dict(report.proposition_sites)["add"] >= 1


def test_registry_empty_declaration_is_not_a_pass() -> None:
    registry = load_registry(HERE)
    report = checks.check_registry_residuals([], {"e": {}}, registry, "e")
    assert report.declaration_failure
    assert report.coverage.denominator == 0


def test_registry_nonzero_is_verdict_not_operational_gap() -> None:
    registry = load_registry(HERE)
    mu, rho = sp.symbols("mu_R rho_br", positive=True)
    row = {
        "relation_id": "R4",
        "engine": "e",
        "qids": {
            "Q.brane.mu_R": "MU",
            "Q.brane.rho_br": "RHO",
            "Q.brane.c_gamma": "C",
        },
    }
    outputs = {"e": {"MU": parsed(mu), "RHO": parsed(rho), "C": parsed(mu / rho)}}
    report = checks.check_registry_residuals([row], outputs, registry, "e")
    assert report.rows[0].status == "NONZERO"
    assert report.coverage.numerator == report.coverage.denominator == 1


def test_s10_config_declares_exact_action_rows_and_cardinalities() -> None:
    config = checks.load_config(S10_CONFIG)
    lagrangians = [row for row in config["cross_engine"] if row.get("family") == "lagrangian"]
    eulers = [row for row in config["cross_engine"] if row.get("family") == "euler_lagrange"]
    assert len(lagrangians) == len(eulers) == 13
    assert all(row["cardinality"] == {"kind": "scalar"} for row in lagrangians)
    assert all(row["cardinality"] == {"kind": "sequence", "count": row["dimension"]} for row in eulers)
    assert all(row.get("mode", "symbolic") != "multiset" for row in lagrangians + eulers)


def test_action_family_is_bound_to_its_declared_tag_identity() -> None:
    config = checks.load_config(S10_CONFIG)
    row = copy.deepcopy(
        next(row for row in config["cross_engine"] if row["quantity"] == "main_d2_q1_lagrangian_action")
    )
    row["wl"] = "WL_S10_MAIN_D2_Q3_DETERMINANT"
    row["py"] = "PY_S10_MAIN_D2_Q3_DETERMINANT"
    outputs = {"wl": real_values("S10", "wl"), "py": real_values("S10", "py")}
    report = checks.check_cross_engine([row], outputs, config=config, registry=load_registry(HERE))
    assert report.rows[0].status == "DECLARATION_ERROR"
    action = checks.action_family_reports(report)[0]
    assert (action.coverage.numerator, action.coverage.denominator) == (0, 1)


def test_all_real_s10_action_rows_reach_comparison_verdicts() -> None:
    config = checks.load_config(S10_CONFIG)
    rows = [row for row in config["cross_engine"] if row.get("family")]
    outputs = {"wl": real_values("S10", "wl"), "py": real_values("S10", "py")}
    report = checks.check_cross_engine(rows, outputs, config=config, registry=load_registry(HERE))
    assert len(report.rows) == 26
    assert all(row.status in checks.VERDICTS for row in report.rows)


def test_derivative_canonicalisation_keeps_coefficient_control_nonvacuous() -> None:
    row = {
        "quantity": "coefficient",
        "left": "L",
        "right": "R",
        "cardinality": {"kind": "scalar"},
    }
    left = checks.normalize_mathematica("2*Derivative[1,0][u][x,t]^2")
    right = checks.normalize_sympy("3*Derivative(u(t,x),x)**2")
    outputs = {"left": {"L": left}, "right": {"R": right}}
    assert checks.check_cross_engine([row], outputs).rows[0].status == "DISAGREE"


def test_derivative_canonicalisation_keeps_form_control_nonvacuous() -> None:
    row = {
        "quantity": "form",
        "left": "L",
        "right": "R",
        "cardinality": {"kind": "scalar"},
    }
    left = checks.normalize_mathematica("Derivative[1,0][u][x,t]^2")
    right = checks.normalize_sympy(
        "Derivative(u(t,x),x)**2 + Derivative(u(t,x),t)**2"
    )
    outputs = {"left": {"L": left}, "right": {"R": right}}
    assert checks.check_cross_engine([row], outputs).rows[0].status == "DISAGREE"


def test_physics_disagree_does_not_itself_enter_operational_status_set() -> None:
    assert "DISAGREE" not in checks.OPERATIONAL_CROSS_STATUSES


# Updated forms of every pre-rebuild regression.  They retain their historical
# names so no legacy check silently disappears; changed assertions name behavior
# that the declaration-oriented rebuild intentionally superseded.


def test_real_s9_end_to_end_cli() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(HERE / "engine_output_checks.py"),
            "--config",
            str(S9_CONFIG),
            "--output",
            f"wl={REAL['S9']['wl']}",
            "--output",
            f"py={REAL['S9']['py']}",
        ],
        cwd=LEDGER,
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    # This is intentionally updated: real operational dimension findings now
    # make S9 nonzero, while its complete comparison report is still emitted.
    assert completed.returncode == 2
    assert "CROSS_ENGINE_COVERAGE: numerator=12 denominator=12" in completed.stdout
    assert "factored_determinant: AGREE" in completed.stdout
    assert "identities=[py:omega2->omega**2]" in completed.stdout
    assert "DIMENSIONS[wl]" in completed.stdout and "DIMENSIONS[py]" in completed.stdout
    assert "OPERATIONAL_FAILURE" in completed.stderr


def test_real_s9_naming_effect_is_measured_before_declarations() -> None:
    config = checks.load_config(S9_CONFIG)
    outputs = {"wl": real_values("S9", "wl"), "py": real_values("S9", "py")}
    report = checks.check_cross_engine(
        config["cross_engine"], outputs, config=config, registry=load_registry(HERE)
    )
    assert report.naming_before_agree == 8
    assert report.naming_after_agree == 12
    assert report.naming_changed_rows == (
        "factored_determinant",
        "full_root_set",
        "transverse_speed_squared",
        "dispersion_scaling_residual_flexural",
    )


@pytest.mark.parametrize(
    ("correct", "mistyped"),
    [("symbol_identities", "symbol_identity"), ("dimension_sources", "dimension_source")],
)
def test_mistyped_declaration_key_is_rejected(
    tmp_path: Path, correct: str, mistyped: str
) -> None:
    config = checks.load_config(S9_CONFIG)
    config[mistyped] = config.pop(correct)
    path = tmp_path / "mistyped.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with pytest.raises(checks.HarnessError, match=rf"unrecognised key.*{mistyped}"):
        checks.load_config(path)


def test_deleting_real_control_tag_fires_parity_missing() -> None:
    config = checks.load_config(S9_CONFIG)
    values = dict(real_values("S9", "wl"))
    del values["WL_S9_X1_LAGRANGIAN"]
    reports = checks.check_tag_parity({"wl": values}, config)
    assert isinstance(reports, tuple)
    row = next(row for row in reports[0].rows if row.cell == "X1")
    assert "LAGRANGIAN" in row.missing


def test_parity_exclusion_is_reported_and_confined_to_parity() -> None:
    config = control_config()
    values = checks.normalize_tags(
        {
            "E_MAIN_D3_A": "1", "E_FULL_D3_A": "1", "E_DIV_D3_A": "1",
            "E_MAIN_D3_B": "2", "E_FULL_D3_B": "2", "E_DIV_D3_B": "2",
        }
    )
    baseline = checks.check_tag_parity({"e": values}, config)
    stale = copy.deepcopy(config)
    stale["parity_exclude"] = ["_DIV_"]
    # Intentionally updated: B2 makes declared identities authoritative, so a
    # legacy regex exclusion cannot remove a declared cell from parity.
    assert checks.check_tag_parity({"e": values}, stale) == baseline


def test_absent_parity_exclude_preserves_report_bytes() -> None:
    config = dimension_config()
    config.update(
        {
            "cross_engine": [],
            "registry_residual": [],
            "cells": [
                {"package": "MAIN", "role": "main"},
                {"package": "CONTROL", "role": "control"},
            ],
            "control": {
                "required_suffixes": {"e": ["VALUE"]},
                "tag_templates": {"e": {"main": "E_{package}_{suffix}", "control": "E_{package}_{suffix}"}},
            },
            "dimension_sources": [
                {"engine": "e", "package": "MAIN", "symbols": {"rho": "E_MAIN_DIM"}},
                {"engine": "e", "package": "CONTROL", "symbols": {"rho": "E_CONTROL_DIM"}},
            ],
        }
    )
    values = checks.normalize_tags(
        {"E_MAIN_VALUE": '"same"', "E_CONTROL_VALUE": '"same"', "E_MAIN_DIM": "{1,0,0}", "E_CONTROL_DIM": "{1,0,0}"}
    )
    baseline = checks.format_report(checks.run_checks(config, {"e": values}))
    with_empty_legacy_field = copy.deepcopy(config)
    with_empty_legacy_field["parity_exclude"] = []
    assert checks.format_report(checks.run_checks(with_empty_legacy_field, {"e": values})) == baseline


def test_changing_real_mapped_operand_fires_disagree_with_both_operands() -> None:
    config = checks.load_config(S9_CONFIG)
    row = next(row for row in config["cross_engine"] if row["quantity"] == "full_root_set")
    outputs = {"wl": real_values("S9", "wl"), "py": dict(real_values("S9", "py"))}
    original = outputs["py"][row["py"]]
    assert isinstance(original, checks.ParsedValue)
    sequence = list(original.value)
    sequence[0] = sequence[0] + sp.Symbol("cross_engine_corruption")
    outputs["py"][row["py"]] = parsed(sequence)
    observed = checks.check_cross_engine(
        [row], outputs, config=config, registry=load_registry(HERE)
    ).rows[0]
    assert observed.status == "DISAGREE"
    assert {engine for engine, _ in observed.operands} == {"wl", "py"}
    assert "cross_engine_corruption" in dict(observed.operands)["py"]


def test_declared_identity_cannot_reconcile_altered_real_operand() -> None:
    config = checks.load_config(S9_CONFIG)
    row = next(row for row in config["cross_engine"] if row["quantity"] == "factored_determinant")
    outputs = {"wl": real_values("S9", "wl"), "py": dict(real_values("S9", "py"))}
    original = outputs["py"][row["py"]]
    assert isinstance(original, checks.ParsedValue)
    outputs["py"][row["py"]] = parsed(2 * original.value)
    observed = checks.check_cross_engine(
        [row], outputs, config=config, registry=load_registry(HERE)
    ).rows[0]
    assert observed.status == "DISAGREE"
    assert observed.identities_applied == ("py:omega2->omega**2",)


def test_symbolic_comparison_accepts_different_algebra_and_vector_shapes() -> None:
    rows = [
        {"quantity": "algebra", "wl": "WL", "py": "PY", "cardinality": {"kind": "scalar"}},
        {"quantity": "vector", "wl": "WLV", "py": "PYV", "cardinality": {"kind": "sequence", "count": 2}},
    ]
    outputs = {
        "wl": checks.normalize_tags({"WL": "a/b", "WLV": "{1,2}"}),
        "py": checks.normalize_tags({"PY": "a*b**-1", "PYV": "Matrix([[1],[2]])"}),
    }
    assert [row.status for row in checks.check_cross_engine(rows, outputs).rows] == ["AGREE", "AGREE"]


def test_real_boolean_is_non_dimensional_not_unknown_or_unparsed() -> None:
    config = checks.load_config(S9_CONFIG)
    values = real_values("S9", "wl")
    target = "WL_S9_ASSUMPTIONS"
    report = checks.check_dimensions(values, config, engine="wl")
    row = next(row for row in report.statuses if row.tag == target)
    assert isinstance(values[target], checks.ParsedValue)
    assert values[target].kind == checks.ValueKind.BOOLEAN
    # Intentionally updated: a failed lookup is now visibly unassessable.
    assert row.status == "unassessable"


def test_unknown_symbol_in_real_dimensionful_expression_is_named_and_operational() -> None:
    config = checks.load_config(S9_CONFIG)
    values = dict(real_values("S9", "wl"))
    target = "WL_S9_LAGRANGIAN"
    original = values[target]
    assert isinstance(original, checks.ParsedValue)
    values[target] = parsed(original.value + sp.Symbol("unknown_dimension_symbol"))
    report = checks.check_dimensions(values, config, engine="wl")
    assert (target, "unknown_dimension_symbol") in report.unknown_symbols
    row = next(row for row in report.statuses if row.tag == target)
    assert row.status == "unassessable"


def test_disagree_is_a_finding_and_cli_exit_is_zero(tmp_path: Path) -> None:
    config, outputs = write_cli_fixture(tmp_path)
    completed = run_harness_cli(config, outputs)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "CROSS_ENGINE: disagree=1" in completed.stdout
    assert "PHYSICS_DISAGREEMENTS (1):" in completed.stdout
    assert "OPERATIONAL_FAILURE" not in completed.stderr


def test_malformed_real_cas_is_unparsed_and_cli_exit_is_nonzero(tmp_path: Path) -> None:
    config, outputs = write_cli_fixture(tmp_path, malformed_py_operand=True)
    completed = run_harness_cli(config, outputs)
    assert completed.returncode == 2
    assert "CROSS_ENGINE: unparsed=1" in completed.stdout
    assert "OPERATIONAL_FAILURE" in completed.stderr


def test_control_response_becomes_invariant_when_sole_changed_control_is_copied() -> None:
    config = checks.load_config(S9_CONFIG)
    values = dict(real_values("S9", "wl"))
    original = checks.check_control_response(values, config, engine="wl")
    candidate = original.responsive[0]
    for tag in candidate.differing_controls:
        values[tag] = values[candidate.main_tag]
    observed = checks.check_control_response(values, config, engine="wl")
    row = next(row for row in observed.rows if row.cell == candidate.cell and row.suffix == candidate.suffix)
    assert row.status == "INVARIANT"


def test_mixed_dimension_sum_names_its_summands() -> None:
    config = dimension_config(primitive_x=True)
    values = checks.normalize_tags({"E_MAIN_DIM": "{1,0,0}"})
    values["E_MAIN_MIXED"] = parsed(sp.Symbol("x") + sp.Symbol("t"))
    report = checks.check_dimensions(values, config, engine="e")
    issue = next(issue for issue in report.non_homogeneous if issue.tag == "E_MAIN_MIXED")
    assert {term.expression for term in issue.summands} == {"t", "x"}
    assert all(term.dimension is not None for term in issue.summands)


def test_invalid_real_value_is_unparsed_and_never_skipped() -> None:
    values = dict(real_values("S9", "wl"))
    target = next(iter(values))
    values[target] = checks.normalize_mathematica("this is not valid CAS text !")
    assert isinstance(values[target], checks.Unparsed)
    assert sum(isinstance(value, checks.Unparsed) for value in values.values()) == 1


def test_duplicate_real_tag_raises() -> None:
    text = REAL["S9"]["wl"].read_text(encoding="utf-8")
    parsed_tags = checks.parse_tagged_output(text)
    tag = next(iter(parsed_tags))
    observed = checks.parse_tagged_output(text + f"\n{tag}: {parsed_tags[tag]}\n")
    # Intentionally updated: a duplicate is retained as operational metadata
    # and does not erase the remaining report.
    assert observed.duplicate_tags == (tag,)
    assert len(observed) == len(parsed_tags)


def test_cross_engine_comparison_is_symbolic_in_both_directions() -> None:
    rows = [
        {"quantity": "equal", "left": "LE", "right": "RE", "cardinality": {"kind": "scalar"}},
        {"quantity": "unequal", "left": "LN", "right": "RN", "cardinality": {"kind": "scalar"}},
    ]
    outputs = {
        "left": checks.normalize_tags({"LE": "a/b", "LN": "a+b"}),
        "right": checks.normalize_tags({"RE": "a*b**-1", "RN": "a+b+c"}),
    }
    assert [row.status for row in checks.check_cross_engine(rows, outputs).rows] == ["AGREE", "DISAGREE"]


def test_root_multiset_comparison_is_explicitly_order_insensitive() -> None:
    left = checks.normalize_mathematica("{a/b,c}")
    right = checks.normalize_mathematica("{c,a*b^-1}")
    assert checks.symbolic_multiset_equal(left, right)


def test_registry_residual_target_is_generated_at_compare_time_and_can_fail() -> None:
    registry = load_registry(HERE)
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
    output = {tag: parsed(all_values[qid]) for qid, tag in qid_tags.items()}
    row = {"relation_id": relation.relation_id, "engine": "engine", "qids": qid_tags}
    passing = checks.check_registry_residuals([row], {"engine": output}, registry, "engine")
    corrupted = dict(output)
    corrupted[qid_tags[relation.designated_output]] = parsed(generated_output + sp.Symbol("corruption"))
    failing = checks.check_registry_residuals([row], {"engine": corrupted}, registry, "engine")
    assert passing.rows[0].status == "ZERO"
    assert failing.rows[0].status == "NONZERO"
    assert failing.rows[0].residual != 0
