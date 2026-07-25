#!/usr/bin/env python3
"""Stage-manifest v2.1 composite integration checker.

Pipeline:
  Draft 2020-12 schema -> indexes/exact parsing -> evidence/adjudication
  integrity -> IMPORT-COMPLETENESS -> integration checks -> coverage -> C7
  mutations.

The deterministic ``--self-test`` suite builds temporary v2.1 fixtures and
proves that every review evasion is able to fail in isolation.
"""

from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import importlib.util
import itertools
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import sympy as sp
from jsonschema import Draft202012Validator
from sympy.core.function import AppliedUndef
from sympy.core.numbers import NumberSymbol
from sympy.parsing.sympy_parser import parse_expr, standard_transformations


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = Path(__file__).with_name("stage_manifest_schema_v2.json")
DEFAULT_CONFIG = Path(__file__).with_name("composite_config.json")
DEFAULT_MANIFEST_DIR = Path(__file__).with_name("stages")
DEFAULT_REPORT = Path(__file__).with_name("reports") / "composite_report.md"
DIMENSIONAL_CONSISTENCY_SCOPE = (
    "Verifies manifest-internal dimensional algebra: relation homogeneity, "
    "declared-vs-recovered agreement for symbols whose dimensions are "
    "recoverable from the stage's audit script, and cross-stage agreement of "
    "shared quantities. It does NOT independently certify that a stage's "
    "dimensions are physically correct — that is owned by the stage's "
    "dual-engine unit audit."
)
CHECK_NAMES = (
    "SCHEMA",
    "EVIDENCE",
    "ADJUDICATION",
    "IMPORT",
    "C1",
    "C2",
    "C3",
    "DIMENSIONAL_CONSISTENCY",
    "C5",
    "C6",
    "C7",
    "C8",
)
BAD_STATUSES = {"FAIL", "UNSUPPORTED"}
STATUS_PRIORITY = {"SKIPPED": -1, "PASS": 0, "PARTIAL": 1, "UNSUPPORTED": 2, "FAIL": 3}
FUNCTION_NAMES = {
    "exp": sp.exp,
    "log": sp.log,
    "sqrt": sp.sqrt,
    "sin": sp.sin,
    "cos": sp.cos,
    "tan": sp.tan,
    "sinh": sp.sinh,
    "cosh": sp.cosh,
    "tanh": sp.tanh,
    "atan": sp.atan,
    "Abs": sp.Abs,
    "Derivative": sp.Derivative,
    "Integral": sp.Integral,
    "Sum": sp.Sum,
    "Function": sp.Function,
    "Eq": sp.Eq,
    "Rational": sp.Rational,
    "oo": sp.oo,
    "pi": sp.pi,
}
LIVE_DIM_RECOVERY_CACHE: dict[
    tuple[Path, str, tuple[str, ...], tuple[str, ...]], dict[str, Any]
] = {}
BARE_TUPLE_DIM_ORDER_BY_SHA256 = {
    # stage032: dim_L=(1,0,0) pins L first; dim_A=(3,-2,1) is asserted
    # against target M L^3 T^-2, independently pinning T second and M third.
    "a6d2b2edaf497b2afa7b51ee87be57e42c7d93f21899a468dd84da77632c1804": "LTM",
    # stage038 states this order next to its four-component Dimension alias.
    "80c6d8292d8eb9b78fe7a39b4fd5ce13d7d1cc5c87de91c5b6346acc27cbfbfe": (
        "M",
        "L",
        "T",
        "E-charge",
    ),
    # stage042's independent unit guard defines the three unit vectors as
    # stiffness_dim, length_dim, and frequency_dim in this positional order.
    "201578946a2424c096c34fef5f95d963dbc1a50fa71f07dee5e81fcc54390bd1": (
        "stiffness",
        "length",
        "time",
    ),
}
DIM_ANCHOR_EXCEPTIONS: frozenset[tuple[str, str]] = frozenset()
DIM_ANCHOR_PATTERN = re.compile(
    r"^scripts/ledger_(stage[0-9]{3})_[^/]+_sympy_audit\.py$"
)
GIT_HEAD_TRACKED_CACHE: dict[Path, bool] = {}


@dataclass
class Issue:
    code: str
    detail: str

    def render(self) -> str:
        return f"{self.code}: {self.detail}"


@dataclass
class CheckResult:
    name: str
    status: str = "PASS"
    issues: list[Issue] = field(default_factory=list)

    def add(self, status: str, code: str, detail: str) -> None:
        if STATUS_PRIORITY[status] > STATUS_PRIORITY[self.status]:
            self.status = status
        self.issues.append(Issue(code, detail))


@dataclass(frozen=True)
class DimensionBasis:
    id: str
    axes: tuple[str, ...]

    def render(self) -> str:
        return f"{self.id} [{', '.join(self.axes)}]"


LEGACY_DIMENSION_BASIS = DimensionBasis("LMT", ("L", "M", "T"))


@dataclass
class Coverage:
    resolved_citations: int = 0
    total_citations: int = 0
    checked_claims: int = 0
    total_claims: int = 0
    unresolved_producers: int = 0
    causal_closure: set[str] = field(default_factory=set)
    c7_edges_run: int = 0
    c7_edges_total: int = 0
    mathematica_outputs_checked: int = 0


@dataclass
class CompositeReport:
    results: dict[str, CheckResult]
    coverage: Coverage
    dimension_bases: dict[str, DimensionBasis] = field(default_factory=dict)

    @property
    def headline(self) -> str:
        statuses = [r.status for r in self.results.values()]
        if "FAIL" in statuses:
            return "FAIL"
        if "UNSUPPORTED" in statuses:
            return "UNSUPPORTED"
        if "PARTIAL" in statuses:
            return "PARTIAL"
        if "SKIPPED" in statuses:
            return "PARTIAL"
        return "PASS"

    def matrix_line(self) -> str:
        c = self.coverage
        closure = ",".join(sorted(c.causal_closure)) or "-"
        return (
            f"citations={c.resolved_citations}/{c.total_citations} "
            f"claims={c.checked_claims}/{c.total_claims} "
            f"unresolved_producers={c.unresolved_producers} "
            f"c7_edges={c.c7_edges_run}/{c.c7_edges_total} "
            f"closure={closure} "
            f"mathematica_outputs_checked={c.mathematica_outputs_checked}"
        )


class UnsupportedExpression(Exception):
    pass


@dataclass(frozen=True)
class DimSourceSnapshot:
    data: bytes
    sha256: str
    tree: ast.Module


def read_dim_source_snapshot(path: Path) -> DimSourceSnapshot:
    """Hash and parse the same single read of a dimension source."""

    data = path.read_bytes()
    return DimSourceSnapshot(
        data=data,
        sha256=hashlib.sha256(data).hexdigest(),
        tree=ast.parse(data, filename=str(path)),
    )


def _bare_tuple_dim_recovery(
    tree: ast.Module,
    dimension_count: int,
    extra_expressions: Sequence[str],
) -> dict[str, Any]:
    """Recover exact function-local ``dim_*`` tuples without executing source.

    Positional order is deliberately absent here.  The caller supplies an
    independently digest-pinned order; this routine recovers only raw vectors
    and exact dimension algebra rooted in real AST Name assignments.
    """

    Vector = tuple[Fraction, ...]
    records: dict[str, set[Vector]] = defaultdict(set)
    binding_counts: dict[str, int] = defaultdict(int)
    binding_expressions: dict[str, list[str]] = defaultdict(list)
    unsupported: set[str] = set()
    expression_envs: list[dict[str, Vector]] = []

    def scalar(node: ast.AST, env: Mapping[str, Fraction]) -> Fraction:
        if isinstance(node, ast.Constant):
            if isinstance(node.value, int) and not isinstance(node.value, bool):
                return Fraction(node.value)
            raise UnsupportedExpression("bare tuple component is not an exact integer")
        if isinstance(node, ast.Name) and node.id in env:
            return env[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = scalar(node.operand, env)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp):
            left = scalar(node.left, env)
            right = scalar(node.right, env)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow) and right.denominator == 1:
                return left ** right.numerator
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "Fraction"
            and not node.keywords
            and len(node.args) in {1, 2}
        ):
            numerator = scalar(node.args[0], env)
            denominator = (
                scalar(node.args[1], env) if len(node.args) == 2 else Fraction(1)
            )
            return numerator / denominator
        raise UnsupportedExpression(
            f"unsupported exact bare-tuple scalar expression {ast.unparse(node)}"
        )

    def tuple_comprehension(
        node: ast.Call, env: Mapping[str, Vector]
    ) -> Vector | None:
        if (
            not isinstance(node.func, ast.Name)
            or node.func.id != "tuple"
            or len(node.args) != 1
            or node.keywords
            or not isinstance(node.args[0], ast.GeneratorExp)
        ):
            return None
        generator = node.args[0]
        if len(generator.generators) != 1:
            return None
        clause = generator.generators[0]
        if clause.ifs or clause.is_async:
            return None
        if (
            not isinstance(clause.target, (ast.Tuple, ast.List))
            or not isinstance(clause.iter, ast.Call)
            or not isinstance(clause.iter.func, ast.Name)
            or clause.iter.func.id != "zip"
            or clause.iter.keywords
        ):
            return None
        target_names = [
            item.id for item in clause.target.elts if isinstance(item, ast.Name)
        ]
        if len(target_names) != len(clause.target.elts):
            return None
        vectors: list[Vector] = []
        for item in clause.iter.args:
            if not isinstance(item, ast.Name) or item.id not in env:
                return None
            vectors.append(env[item.id])
        if len(vectors) != len(target_names):
            return None
        recovered = []
        for components in zip(*vectors):
            component_env = dict(zip(target_names, components))
            recovered.append(scalar(generator.elt, component_env))
        if len(recovered) != dimension_count:
            return None
        return tuple(recovered)

    def vector(
        node: ast.AST,
        env: Mapping[str, Vector],
        *,
        dimension_algebra: bool = False,
    ) -> Vector | None:
        if isinstance(node, ast.Name):
            return env.get(node.id)
        if (
            isinstance(node, (ast.Tuple, ast.List))
            and len(node.elts) == dimension_count
        ):
            try:
                values = tuple(scalar(item, {}) for item in node.elts)
            except UnsupportedExpression:
                return None
            return values
        if isinstance(node, ast.Call):
            try:
                return tuple_comprehension(node, env)
            except UnsupportedExpression:
                return None
        if dimension_algebra and isinstance(node, ast.BinOp):
            left = vector(node.left, env, dimension_algebra=True)
            if isinstance(node.op, ast.Pow) and left is not None:
                try:
                    power = scalar(node.right, {})
                except UnsupportedExpression:
                    return None
                return scale_dim(left, power)  # type: ignore[return-value]
            right = vector(node.right, env, dimension_algebra=True)
            if left is None or right is None:
                return None
            if isinstance(node.op, ast.Mult):
                return add_dim(left, right)  # type: ignore[return-value]
            if isinstance(node.op, ast.Div):
                return sub_dim(left, right)  # type: ignore[return-value]
        return None

    def dim_target_name(target: ast.AST) -> str | None:
        while isinstance(target, (ast.Subscript, ast.Attribute)):
            target = target.value
        if isinstance(target, ast.Name) and (
            target.id.startswith("dim_") or target.id.endswith("_dim")
        ):
            return target.id
        return None

    def assignments(
        scope_node: ast.AST,
    ) -> tuple[list[ast.Assign | ast.AnnAssign], set[str]]:
        found: list[ast.Assign | ast.AnnAssign] = []
        unsafe_names: set[str] = set()

        class UnsafeAssignmentVisitor(ast.NodeVisitor):
            def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
                return

            def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
                return

            def visit_ClassDef(self, node: ast.ClassDef) -> None:
                return

            def visit_Lambda(self, node: ast.Lambda) -> None:
                return

            def visit_Assign(self, node: ast.Assign) -> None:
                for target in node.targets:
                    name = dim_target_name(target)
                    if name is not None:
                        unsafe_names.add(name)

            def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
                name = dim_target_name(node.target)
                if name is not None:
                    unsafe_names.add(name)

            def visit_AugAssign(self, node: ast.AugAssign) -> None:
                name = dim_target_name(node.target)
                if name is not None:
                    unsafe_names.add(name)

            def visit_Delete(self, node: ast.Delete) -> None:
                for target in node.targets:
                    name = dim_target_name(target)
                    if name is not None:
                        unsafe_names.add(name)

        body = scope_node.body if isinstance(
            scope_node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef)
        ) else []
        for statement in body:
            if isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue
            if isinstance(statement, (ast.Assign, ast.AnnAssign)):
                targets = (
                    statement.targets
                    if isinstance(statement, ast.Assign)
                    else [statement.target]
                )
                if all(isinstance(target, ast.Name) for target in targets):
                    found.append(statement)
                else:
                    UnsafeAssignmentVisitor().visit(statement)
                continue
            UnsafeAssignmentVisitor().visit(statement)
        return (
            sorted(found, key=lambda item: (item.lineno, item.col_offset)),
            unsafe_names,
        )

    def target_names(node: ast.Assign | ast.AnnAssign) -> list[str]:
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        return [target.id for target in targets if isinstance(target, ast.Name)]

    def references_vector(node: ast.AST, env: Mapping[str, Vector]) -> bool:
        return any(
            isinstance(child, ast.Name) and child.id in env for child in ast.walk(node)
        )

    def evaluate_scope(
        scope_node: ast.AST, scope: str, inherited: Mapping[str, Vector]
    ) -> dict[str, Vector]:
        env = dict(inherited)
        scope_assignments, unsafe_names = assignments(scope_node)
        unsupported.update(unsafe_names)
        unsupported.update(f"{scope}.{name}" for name in unsafe_names)
        for node in scope_assignments:
            names = target_names(node)
            if not names or node.value is None:
                continue
            dim_names = [
                name
                for name in names
                if name.startswith("dim_") or name.endswith("_dim")
            ]
            if not dim_names:
                continue
            recovered = vector(node.value, env)
            if recovered is None:
                if references_vector(node.value, env):
                    unsupported.update(dim_names)
                    unsupported.update(f"{scope}.{name}" for name in dim_names)
                continue
            expression = ast.unparse(node.value)
            for name in dim_names:
                env[name] = recovered
                records[name].add(recovered)
                records[f"{scope}.{name}"].add(recovered)
                binding_counts[name] += 1
                binding_counts[f"{scope}.{name}"] += 1
                binding_expressions[name].append(expression)
                binding_expressions[f"{scope}.{name}"].append(expression)
        expression_envs.append(env)
        return env

    module_env = evaluate_scope(tree, "module", {})
    functions = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    ]
    function_name_counts: dict[str, int] = defaultdict(int)
    for node in functions:
        function_name_counts[node.name] += 1
    duplicate_function_names = {
        name for name, count in function_name_counts.items() if count > 1
    }
    for node in functions:
        if node.name in duplicate_function_names:
            duplicate_assignments, duplicate_unsafe = assignments(node)
            duplicate_targets = set(duplicate_unsafe)
            for assignment in duplicate_assignments:
                duplicate_targets.update(
                    name
                    for name in target_names(assignment)
                    if name.startswith("dim_") or name.endswith("_dim")
                )
            unsupported.update(duplicate_targets)
            unsupported.update(
                f"{node.name}.{name}" for name in duplicate_targets
            )
            for name in duplicate_targets:
                binding_counts[name] += 1
                binding_counts[f"{node.name}.{name}"] += 1
            continue
        evaluate_scope(node, node.name, module_env)

    for source in extra_expressions:
        try:
            expression = ast.parse(source, mode="eval").body
        except SyntaxError:
            unsupported.add(f"expr:{source}")
            continue
        key = f"expr:{ast.unparse(expression)}"
        if duplicate_function_names:
            unsupported.add(key)
            continue
        recovered_values = {
            recovered
            for env in expression_envs
            if references_vector(expression, env)
            and (recovered := vector(expression, env, dimension_algebra=True))
            is not None
        }
        if recovered_values:
            records[key].update(recovered_values)
        else:
            unsupported.add(key)

    resolved = {
        name: [str(value) for value in next(iter(values))]
        for name, values in records.items()
        if len(values) == 1
    }
    return {
        "tuples": resolved,
        "unsupported": sorted(unsupported),
        "ambiguous": sorted(name for name, values in records.items() if len(values) > 1),
        "binding_counts": dict(sorted(binding_counts.items())),
        "binding_expressions": {
            name: expressions
            for name, expressions in sorted(binding_expressions.items())
        },
    }


def _live_dim_recovery(
    path: Path,
    order: Sequence[str],
    extra_expressions: Sequence[str] = (),
    source_bytes: bytes | None = None,
) -> dict[str, Any]:
    """Evaluate exact source dimension algebra in a recovery subprocess.

    A live ``Dim`` source is imported after AST validation.  A registered bare
    tuple source is recovered statically and is never executed.
    """

    data = path.read_bytes() if source_bytes is None else source_bytes
    tree = ast.parse(data, filename=str(path))
    if not any(
        isinstance(node, ast.ClassDef) and node.name == "Dim"
        for node in ast.walk(tree)
    ):
        return _bare_tuple_dim_recovery(tree, len(order), extra_expressions)

    module_name = f"_composite_dim_source_{hashlib.sha256(str(path).encode()).hexdigest()[:16]}"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise UnsupportedExpression(f"cannot import dim source {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    sys.path.insert(0, str(path.parent))
    try:
        exec(compile(tree, str(path), "exec"), vars(module))
    finally:
        sys.path.pop(0)

    dim_type = getattr(module, "Dim", None)
    if not isinstance(dim_type, type):
        raise UnsupportedExpression(f"{path} does not expose a live Dim class")

    records: dict[str, set[tuple[str, ...]]] = defaultdict(set)
    binding_counts: dict[str, int] = defaultdict(int)
    binding_expressions: dict[str, list[str]] = defaultdict(list)
    unsupported: set[str] = set()
    expression_envs: list[dict[str, Any]] = []

    def vector(value: Any) -> tuple[str, ...] | None:
        if not isinstance(value, dim_type):
            return None
        try:
            raw = tuple(
                getattr(value, axis_attribute_name(axis)) for axis in order
            )
        except (AttributeError, TypeError) as exc:
            raise UnsupportedExpression(
                f"live Dim value does not expose {list(order)} components: {exc}"
            ) from exc
        if len(raw) != len(order):
            raise UnsupportedExpression(f"live Dim value has {len(raw)} components")
        rendered = tuple(str(item) for item in raw)
        try:
            tuple(Fraction(item) for item in rendered)
        except (ValueError, ZeroDivisionError) as exc:
            raise UnsupportedExpression(
                f"live Dim components are not exact rationals: {rendered}"
            ) from exc
        return rendered

    def add_tuple_record(
        name: str,
        recovered: tuple[str, ...],
        scope: str,
        expression: str,
    ) -> None:
        records[name].add(recovered)
        records[f"{scope}.{name}"].add(recovered)
        binding_counts[name] += 1
        binding_counts[f"{scope}.{name}"] += 1
        binding_expressions[name].append(expression)
        binding_expressions[f"{scope}.{name}"].append(expression)

    def add_record(name: str, value: Any, scope: str, expression: str) -> None:
        recovered = vector(value)
        if recovered is not None:
            add_tuple_record(name, recovered, scope, expression)

    def exact_integer(node: ast.AST) -> bool:
        try:
            value = ast.literal_eval(node)
        except (TypeError, ValueError):
            return False
        return isinstance(value, int) and not isinstance(value, bool)

    def admitted(node: ast.AST) -> bool:
        if isinstance(node, ast.Name):
            return True
        if isinstance(node, ast.Constant):
            return isinstance(node.value, (int, str)) and not isinstance(node.value, bool)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            return admitted(node.operand)
        if isinstance(node, ast.BinOp) and isinstance(
            node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)
        ):
            return admitted(node.left) and admitted(node.right) and (
                not isinstance(node.op, ast.Pow) or exact_integer(node.right)
            )
        if isinstance(node, ast.Call):
            return (
                isinstance(node.func, ast.Name)
                and node.func.id == "Dim"
                and not node.keywords
                and len(node.args) in {0, len(order)}
                and all(admitted(arg) for arg in node.args)
            )
        if isinstance(node, (ast.Tuple, ast.List, ast.Set)):
            return all(admitted(item) for item in node.elts)
        if isinstance(node, ast.Dict):
            return all(
                key is not None
                and admitted(key)
                and admitted(value)
                for key, value in zip(node.keys, node.values)
            )
        return False

    def references_live_dim(node: ast.AST, env: Mapping[str, Any]) -> bool:
        if any(
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Name)
            and child.func.id == "Dim"
            for child in ast.walk(node)
        ):
            return True
        return any(
            isinstance(child, ast.Name) and vector(env.get(child.id)) is not None
            for child in ast.walk(node)
        )

    def assignments(scope_node: ast.AST) -> list[ast.Assign | ast.AnnAssign]:
        found: list[ast.Assign | ast.AnnAssign] = []

        class AssignmentVisitor(ast.NodeVisitor):
            def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
                if node is scope_node:
                    self.generic_visit(node)

            def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
                if node is scope_node:
                    self.generic_visit(node)

            def visit_ClassDef(self, node: ast.ClassDef) -> None:
                return

            def visit_Lambda(self, node: ast.Lambda) -> None:
                return

            def visit_Assign(self, node: ast.Assign) -> None:
                found.append(node)

            def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
                found.append(node)

        visitor = AssignmentVisitor()
        if isinstance(scope_node, ast.Module):
            for statement in scope_node.body:
                if isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    continue
                visitor.visit(statement)
        else:
            visitor.visit(scope_node)
        return sorted(found, key=lambda item: (item.lineno, item.col_offset))

    def target_names(node: ast.Assign | ast.AnnAssign) -> list[str]:
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        return [target.id for target in targets if isinstance(target, ast.Name)]

    def evaluate_scope(
        scope_node: ast.AST,
        scope: str,
        inherited_dim_env: Mapping[str, Any],
    ) -> dict[str, Any]:
        env = dict(vars(module))
        trusted_dim_env = dict(inherited_dim_env)
        for node in assignments(scope_node):
            names = target_names(node)
            if not names or node.value is None:
                continue
            if not admitted(node.value):
                if references_live_dim(node.value, env):
                    unsupported.update(names)
                    unsupported.update(f"{scope}.{name}" for name in names)
                continue
            try:
                value = eval(
                    compile(ast.Expression(node.value), str(path), "eval"),
                    {"__builtins__": {}},
                    env,
                )
            except Exception:
                if references_live_dim(node.value, env):
                    unsupported.update(names)
                    unsupported.update(f"{scope}.{name}" for name in names)
                continue
            for child in ast.walk(node.value):
                if not admitted(child):
                    continue
                try:
                    child_value = eval(
                        compile(ast.Expression(child), str(path), "eval"),
                        {"__builtins__": {}},
                        env,
                    )
                except Exception:
                    continue
                recovered_child = vector(child_value)
                if recovered_child is not None:
                    records[f"expr:{ast.unparse(child)}"].add(recovered_child)
            expression = ast.unparse(node.value)
            for name in names:
                env[name] = value
                add_record(name, value, scope, expression)
                if vector(value) is not None:
                    trusted_dim_env[name] = value
        expression_envs.append(trusted_dim_env)
        return trusted_dim_env

    # Only AST assignments to a direct Name target create recoverable bindings.
    # Importing the module supplies the live Dim class/operator semantics, but a
    # runtime attribute or a nested container label is not dimensional algebra
    # provenance for a manifest symbol.
    module_dim_env = evaluate_scope(tree, "module", {"Dim": dim_type})
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            evaluate_scope(node, node.name, module_dim_env)
    for source in extra_expressions:
        try:
            expression = ast.parse(source, mode="eval").body
        except SyntaxError:
            unsupported.add(f"expr:{source}")
            continue
        canonical_key = f"expr:{ast.unparse(expression)}"
        if not admitted(expression) or any(
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "Dim"
            for node in ast.walk(expression)
        ):
            unsupported.add(canonical_key)
            continue
        recovered_values: set[tuple[str, ...]] = set()
        for env in expression_envs:
            if not any(
                isinstance(node, ast.Name)
                and vector(env.get(node.id)) is not None
                for node in ast.walk(expression)
            ):
                continue
            try:
                value = eval(
                    compile(ast.Expression(expression), str(path), "eval"),
                    {"__builtins__": {}},
                    env,
                )
            except Exception:
                continue
            recovered = vector(value)
            if recovered is not None:
                recovered_values.add(recovered)
        if not recovered_values:
            unsupported.add(canonical_key)
            continue
        records[canonical_key].update(recovered_values)

    resolved = {
        name: list(next(iter(values)))
        for name, values in records.items()
        if len(values) == 1
    }
    ambiguous = sorted(name for name, values in records.items() if len(values) > 1)
    return {
        "tuples": resolved,
        "unsupported": sorted(unsupported),
        "ambiguous": ambiguous,
        "binding_counts": dict(sorted(binding_counts.items())),
        "binding_expressions": {
            name: expressions
            for name, expressions in sorted(binding_expressions.items())
        },
    }


def stage_number(stage_id: str) -> int:
    return int(stage_id.removeprefix("stage"))


def digest_hex(value: str) -> tuple[str, str]:
    if ":" in value:
        algorithm, hexdigest = value.split(":", 1)
        return algorithm, hexdigest
    algorithms = {40: "sha1", 64: "sha256", 128: "sha512"}
    return algorithms[len(value)], value


def hash_file(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def basis_from_manifest(manifest: Mapping[str, Any]) -> DimensionBasis:
    raw = manifest.get("dimension_basis")
    if raw is None:
        return LEGACY_DIMENSION_BASIS
    return DimensionBasis(str(raw["id"]), tuple(str(axis) for axis in raw["axes"]))


def source_order_axes(raw: Any) -> tuple[str, ...]:
    if isinstance(raw, str):
        return tuple(raw)
    if isinstance(raw, (list, tuple)):
        return tuple(str(axis) for axis in raw)
    raise ValueError(f"invalid dim_source_order {raw!r}")


def axis_attribute_name(axis: str) -> str:
    return re.sub(r"[^a-z0-9_]", "_", axis.lower())


def normalized_dim(
    raw: Mapping[str, str] | None,
    axes: Sequence[str] = LEGACY_DIMENSION_BASIS.axes,
) -> tuple[Fraction, ...]:
    raw = raw or {}
    return tuple(Fraction(raw.get(axis, "0")) for axis in axes)


def dim_json(
    dim: Sequence[Fraction],
    axes: Sequence[str] = LEGACY_DIMENSION_BASIS.axes,
) -> dict[str, str]:
    result: dict[str, str] = {}
    for axis, value in zip(axes, dim, strict=True):
        if value:
            result[axis] = str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    return result


def add_dim(a: Sequence[Fraction], b: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(x + y for x, y in zip(a, b, strict=True))


def sub_dim(a: Sequence[Fraction], b: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(x - y for x, y in zip(a, b, strict=True))


def scale_dim(a: Sequence[Fraction], power: Fraction) -> tuple[Fraction, ...]:
    return tuple(x * power for x in a)


def safe_relative_path(value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    candidate = Path(value)
    return not candidate.is_absolute() and ".." not in candidate.parts


def normalized_relative_path(value: str) -> str:
    return Path(value).as_posix()


def git_head_tracks(path: Path) -> bool:
    """Return whether the lexical path is a regular file in its HEAD tree."""

    lexical = path.absolute()
    cached = GIT_HEAD_TRACKED_CACHE.get(lexical)
    if cached is not None:
        return cached
    try:
        worktree = subprocess.run(
            ("git", "-C", str(lexical.parent), "rev-parse", "--show-toplevel"),
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
        relative = lexical.relative_to(Path(worktree).resolve()).as_posix()
        entry = subprocess.run(
            (
                "git",
                "-C",
                worktree,
                "ls-tree",
                "--full-tree",
                "HEAD",
                "--",
                relative,
            ),
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.splitlines()
        tracked = (
            len(entry) == 1
            and entry[0].split(maxsplit=1)[0] in {"100644", "100755"}
            and entry[0].endswith("\t" + relative)
        )
    except (OSError, RuntimeError, subprocess.SubprocessError, ValueError):
        tracked = False
    GIT_HEAD_TRACKED_CACHE[lexical] = tracked
    return tracked


def admissible_dim_anchor(
    stage_id: str,
    source_value: str,
    path: Path,
    roots: Sequence[Path],
) -> bool:
    """Apply checker-owned structural identity; content integrity is separate."""

    normalized = normalized_relative_path(source_value)
    structurally_owned = (
        (match := DIM_ANCHOR_PATTERN.fullmatch(normalized)) is not None
        and match.group(1) == stage_id
    )
    excepted = (stage_id, normalized) in DIM_ANCHOR_EXCEPTIONS
    if not (structurally_owned or excepted):
        return False
    for root in roots:
        lexical = (root / normalized).absolute()
        try:
            same_file = lexical.resolve(strict=True) == path.resolve(strict=True)
        except (OSError, RuntimeError):
            continue
        if same_file and git_head_tracks(lexical):
            return True
    return False


def resolve_path(value: str, roots: Sequence[Path]) -> Path | None:
    if not safe_relative_path(value):
        return None
    candidate = Path(value)
    for root in roots:
        try:
            confined_root = root.resolve()
            resolved = (confined_root / candidate).resolve(strict=True)
            resolved.relative_to(confined_root)
        except (OSError, RuntimeError, ValueError):
            continue
        if resolved.is_file():
            return resolved
    return None


def claim_ref(stage_id: str, local_id: str) -> str:
    return f"{stage_id}/{local_id}"


def local_ref_id(value: str) -> str:
    return value.split("/", 1)[1]


def all_evidence_objects(value: Any) -> Iterator[dict[str, Any]]:
    if isinstance(value, dict):
        required = {"source_path", "source_digest"}
        if required <= value.keys():
            yield value
        for child in value.values():
            yield from all_evidence_objects(child)
    elif isinstance(value, list):
        for child in value:
            yield from all_evidence_objects(child)


def certifying_evidence_objects(value: Any) -> Iterator[dict[str, Any]]:
    if isinstance(value, dict):
        required = {
            "source_path",
            "source_digest",
            "locus",
            "engine",
            "method",
        }
        if required <= value.keys():
            yield value
        for child in value.values():
            yield from certifying_evidence_objects(child)
    elif isinstance(value, list):
        for child in value:
            yield from certifying_evidence_objects(child)


def expected_mathematica_output_path(script_path: str) -> Path:
    script = Path(script_path)
    return script.parent / "out" / f"{script.stem}.out"


def locus_dimension_expression(locus: str) -> str | None:
    match = re.search(
        r"\bdimension expression\s+(.+?)"
        r"(?=\s+in\s+[A-Za-z_][A-Za-z0-9_]*|\s*;|,\s*lines?\b|$)",
        locus,
        re.I,
    )
    if match is None:
        return None
    source = match.group(1).strip()
    try:
        return ast.unparse(ast.parse(source, mode="eval").body)
    except SyntaxError as exc:
        raise UnsupportedExpression(
            f"invalid dimension expression named by locus: {source}"
        ) from exc


def locus_assignment_anchors(locus: str) -> list[tuple[str, str | None]]:
    matches = list(
        re.finditer(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*=(?!=)", locus)
    )
    anchors: list[tuple[str, str | None]] = []
    for index, match in enumerate(matches):
        tail = locus[match.end() :]
        ends = [len(tail)]
        semicolon = tail.find(";")
        if semicolon >= 0:
            ends.append(semicolon)
        lines = re.search(r",\s*lines?\b", tail, re.I)
        if lines is not None:
            ends.append(lines.start())
        if index + 1 < len(matches):
            ends.append(matches[index + 1].start() - match.end())
        raw_expression = tail[: min(ends)].strip().rstrip(",").strip()
        try:
            expression = ast.unparse(
                ast.parse(raw_expression, mode="eval").body
            )
        except SyntaxError:
            expression = None
        anchors.append((match.group(1), expression))
    return anchors


class CompositeChecker:
    def __init__(
        self,
        manifests: Sequence[dict[str, Any]],
        config: dict[str, Any],
        *,
        schema: dict[str, Any] | None = None,
        roots: Sequence[Path] | None = None,
        closed_slice: bool | None = None,
        bare_tuple_dim_orders: Mapping[str, str | Sequence[str]] | None = None,
    ) -> None:
        self.manifests = list(manifests)
        self.config = config
        self.schema = schema or json.loads(SCHEMA_PATH.read_text())
        self.roots = list(roots or (PROJECT_ROOT, PROJECT_ROOT.parent.parent, Path.cwd()))
        config_override = True if config.get("closed_slice") is True else None
        self.closed_slice_override = config_override if closed_slice is None else closed_slice
        self.closed_slice = False
        self.bare_tuple_dim_orders = dict(BARE_TUPLE_DIM_ORDER_BY_SHA256)
        if bare_tuple_dim_orders is not None:
            self.bare_tuple_dim_orders.update(bare_tuple_dim_orders)
        self.results = {name: CheckResult(name) for name in CHECK_NAMES}
        self.coverage = Coverage()
        self.dimension_bases: dict[str, DimensionBasis] = {}
        self.stages: dict[str, dict[str, Any]] = {}
        self.claims: dict[str, dict[str, Any]] = {}
        self.exports: dict[str, dict[str, Any]] = {}
        self.symbols_by_qid: dict[str, list[tuple[str, dict[str, Any]]]] = defaultdict(list)
        self.symbols_by_stage_alias: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
        self.owners: dict[str, tuple[str, str, str]] = {}
        self.parse_cache: dict[tuple[str, str], sp.Expr] = {}
        self.dim_recovery_unsupported: dict[Path, set[str]] = {}
        self.dim_recovery_ambiguous: dict[Path, set[str]] = {}
        self.dim_recovery_binding_counts: dict[Path, dict[str, int]] = {}
        self.dim_recovery_binding_expressions: dict[
            Path, dict[str, list[str]]
        ] = {}

    def run(self) -> CompositeReport:
        self.check_schema()
        if self.results["SCHEMA"].status == "FAIL":
            # Draft validation is deliberately exhaustive. Do not interpret a
            # structurally invalid document with downstream code that assumes
            # v2.1 shapes; this is how legacy/v1 inputs fail without crashing.
            self.mark_skipped(
                CHECK_NAMES[1:],
                "SCHEMA failed; downstream checks were not run",
            )
            self.coverage.total_claims = sum(
                len(m.get("claims", [])) for m in self.manifests if isinstance(m, dict)
            )
            return CompositeReport(
                self.results, self.coverage, self.dimension_bases
            )
        self.build_indexes()
        self.closed_slice = (
            self.infer_closed_slice()
            if self.closed_slice_override is None
            else self.closed_slice_override
        )
        self.check_evidence()
        self.check_adjudications()
        self.check_import_completeness()
        self.check_c1()
        self.check_c2()
        self.check_c3()
        self.check_dimensional_consistency()
        self.check_c5()
        self.check_c6()
        # C7 results are not meaningful when the declared graph itself failed.
        if self.results["C3"].status in BAD_STATUSES:
            self.mark_skipped(("C7",), "C3 failed; C7 was not run")
        else:
            self.check_c7()
        self.check_c8()
        self.coverage.total_claims = sum(len(m.get("claims", [])) for m in self.manifests)
        self.coverage.checked_claims = self.coverage.total_claims
        return CompositeReport(self.results, self.coverage, self.dimension_bases)

    @staticmethod
    def referenced_stage(current_stage: str | None, ref: Any) -> str | None:
        if not isinstance(ref, str):
            return None
        if ref.startswith("here/"):
            return current_stage
        match = re.match(r"^(stage[0-9]{3})(?:/|$)", ref)
        return match.group(1) if match else None

    def required_producer_stages(self) -> set[str]:
        required: set[str] = set()

        def add_ref(current_stage: str | None, ref: Any) -> None:
            producer = self.referenced_stage(current_stage, ref)
            if producer is not None:
                required.add(producer)

        add_ref(None, self.config.get("range_claim_ref"))
        for stage_id, manifest in self.stages.items():
            for symbol in manifest.get("symbols", []):
                add_ref(stage_id, symbol.get("definition_ref"))
            for consume in manifest.get("consumes", []):
                add_ref(stage_id, consume.get("ref"))
                for substitution in consume.get("substitutions", []):
                    add_ref(stage_id, substitution.get("backed_by"))
            for knob in manifest.get("knobs", []):
                add_ref(stage_id, knob.get("origin"))
                add_ref(stage_id, knob.get("effective_stage"))
                add_ref(stage_id, knob.get("discharge_evidence"))
            for claim in manifest.get("claims", []):
                genesis = claim.get("genesis")
                if isinstance(genesis, dict):
                    evidence = genesis.get("evidence")
                    if isinstance(evidence, dict):
                        add_ref(stage_id, evidence.get("predates"))
                    for ref in genesis.get("coordinated_with", []):
                        add_ref(stage_id, ref)
                add_ref(stage_id, claim.get("discharged_by"))
        return required

    def infer_closed_slice(self) -> bool:
        loaded_numbers = {stage_number(stage_id) for stage_id in self.stages}
        plausible_whole = bool(loaded_numbers) and loaded_numbers == set(
            range(1, max(loaded_numbers) + 1)
        )
        return (
            plausible_whole
            and self.required_producer_stages() <= set(self.stages)
        )

    def unresolved_reference_status(self, current_stage: str | None, ref: Any) -> str:
        producer = self.referenced_stage(current_stage, ref)
        if (
            not self.closed_slice
            and producer is not None
            and producer != current_stage
            and producer not in self.stages
        ):
            return "PARTIAL"
        return "FAIL"

    def mark_skipped(self, names: Iterable[str], detail: str) -> None:
        for name in names:
            result = self.results[name]
            result.status = "SKIPPED"
            result.issues.append(Issue("PREREQUISITE", detail))

    def check_schema(self) -> None:
        validator = Draft202012Validator(self.schema)
        for index, manifest in enumerate(self.manifests):
            for error in sorted(validator.iter_errors(manifest), key=lambda e: list(e.absolute_path)):
                path = "/".join(str(part) for part in error.absolute_path) or "<root>"
                stage = manifest.get("stage_id", f"manifest[{index}]")
                self.results["SCHEMA"].add("FAIL", "SCHEMA_VALIDATION", f"{stage}:{path}: {error.message}")
        if self.results["SCHEMA"].status == "FAIL":
            return

        basis_ids: dict[str, tuple[str, ...]] = {}
        for manifest in self.manifests:
            stage_id = manifest["stage_id"]
            basis = basis_from_manifest(manifest)
            self.dimension_bases[stage_id] = basis
            previous_axes = basis_ids.setdefault(basis.id, basis.axes)
            if previous_axes != basis.axes:
                self.results["SCHEMA"].add(
                    "FAIL",
                    "DIMENSION_BASIS_ID_CONFLICT",
                    f"{stage_id}: basis {basis.id!r} declares {basis.axes}, "
                    f"previously {previous_axes}",
                )
            allowed_axes = set(basis.axes)

            def validate_dim(raw: Any, path: str) -> None:
                if not isinstance(raw, dict):
                    return
                unknown = set(raw) - allowed_axes
                if unknown:
                    self.results["SCHEMA"].add(
                        "FAIL",
                        "DIMENSION_AXIS_NOT_IN_BASIS",
                        f"{stage_id}:{path}: {sorted(unknown)} not in "
                        f"{list(basis.axes)}",
                    )

            def validate_signature(raw: Any, path: str) -> None:
                if not isinstance(raw, dict):
                    return
                for index, item in enumerate(raw.get("domain", [])):
                    validate_dim(item, f"{path}/domain/{index}")
                validate_dim(raw.get("codomain"), f"{path}/codomain")

            for index, symbol in enumerate(manifest.get("symbols", [])):
                validate_dim(symbol.get("dim"), f"symbols/{index}/dim")
                validate_signature(
                    symbol.get("signature"), f"symbols/{index}/signature"
                )
                declared_order = source_order_axes(
                    symbol.get("dim_source_order")
                )
                if (
                    len(declared_order) != len(basis.axes)
                    or set(declared_order) != allowed_axes
                ):
                    self.results["SCHEMA"].add(
                        "FAIL",
                        "DIM_SOURCE_ORDER_BASIS_MISMATCH",
                        f"{stage_id}:symbols/{index}: order "
                        f"{list(declared_order)} does not span basis "
                        f"{list(basis.axes)}",
                    )
                declared_tuple = symbol.get("dim_source_tuple", [])
                if len(declared_tuple) != len(basis.axes):
                    self.results["SCHEMA"].add(
                        "FAIL",
                        "DIM_SOURCE_TUPLE_BASIS_MISMATCH",
                        f"{stage_id}:symbols/{index}: tuple length "
                        f"{len(declared_tuple)} != basis rank "
                        f"{len(basis.axes)}",
                    )
            for index, claim in enumerate(manifest.get("claims", [])):
                validate_dim(
                    claim.get("expected_dim"), f"claims/{index}/expected_dim"
                )
                validate_signature(
                    claim.get("payload", {}).get("signature"),
                    f"claims/{index}/payload/signature",
                )
            for index, consume in enumerate(manifest.get("consumes", [])):
                validate_dim(
                    consume.get("as_consumed_dim"),
                    f"consumes/{index}/as_consumed_dim",
                )
                validate_signature(
                    consume.get("as_consumed", {})
                    .get("payload", {})
                    .get("signature"),
                    f"consumes/{index}/as_consumed/payload/signature",
                )

    def build_indexes(self) -> None:
        for manifest in self.manifests:
            stage_id = manifest.get("stage_id")
            if not isinstance(stage_id, str):
                continue
            if stage_id in self.stages:
                self.results["C1"].add("FAIL", "DUPLICATE_STAGE", stage_id)
            self.stages[stage_id] = manifest
            self.coverage.causal_closure.add(stage_id)
            for claim in manifest.get("claims", []):
                cid = claim.get("id")
                if isinstance(cid, str):
                    key = claim_ref(stage_id, cid)
                    if key in self.claims:
                        self.results["C1"].add("FAIL", "DUPLICATE_CLAIM", key)
                    self.claims[key] = claim
            for export in manifest.get("exports", []):
                cid = export.get("claim_id")
                if isinstance(cid, str):
                    self.exports[claim_ref(stage_id, cid)] = export
            for symbol in manifest.get("symbols", []):
                qid = symbol.get("quantity_id")
                alias = symbol.get("parse_alias")
                if isinstance(qid, str):
                    self.symbols_by_qid[qid].append((stage_id, symbol))
                if isinstance(alias, str):
                    if alias in self.symbols_by_stage_alias[stage_id]:
                        self.results["C1"].add("FAIL", "DUPLICATE_PARSE_ALIAS", f"{stage_id}:{alias}")
                    self.symbols_by_stage_alias[stage_id][alias] = symbol

        for qid, appearances in self.symbols_by_qid.items():
            candidates: list[tuple[int, str, str, str]] = []
            for stage_id, symbol in appearances:
                definition = symbol.get("definition_ref", "")
                if isinstance(definition, str) and definition.startswith("here/"):
                    cid = local_ref_id(definition)
                    claim = self.claims.get(claim_ref(stage_id, cid))
                    if claim and self.claim_binds_symbol(claim, symbol.get("parse_alias", ""), stage_id):
                        candidates.append((stage_number(stage_id), stage_id, cid, symbol.get("parse_alias", "")))
            if candidates:
                _, owner_stage, owner_claim, owner_alias = min(candidates)
                self.owners[qid] = (owner_stage, owner_claim, owner_alias)

    def local_dict(self, stage_id: str) -> dict[str, Any]:
        result = dict(FUNCTION_NAMES)
        for alias, symbol in self.symbols_by_stage_alias.get(stage_id, {}).items():
            if "signature" in symbol:
                result[alias] = sp.Function(alias)
            else:
                result[alias] = sp.Symbol(alias, real=True)
        return result

    def stage_basis(self, stage_id: str) -> DimensionBasis:
        return self.dimension_bases.get(stage_id, LEGACY_DIMENSION_BASIS)

    def normalized_stage_dim(
        self, stage_id: str, raw: Mapping[str, str] | None
    ) -> tuple[Fraction, ...]:
        return normalized_dim(raw, self.stage_basis(stage_id).axes)

    def zero_dim(self, stage_id: str) -> tuple[Fraction, ...]:
        return normalized_dim({}, self.stage_basis(stage_id).axes)

    def parse(self, stage_id: str, source: str) -> sp.Expr:
        key = (stage_id, source)
        if key in self.parse_cache:
            return self.parse_cache[key]
        if "__" in source or re.search(r"[\[\]{};'\"`]", source):
            raise UnsupportedExpression(f"forbidden syntax: {source}")
        identifiers = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", source))
        allowed = set(self.local_dict(stage_id))
        unknown = identifiers - allowed
        if unknown:
            raise UnsupportedExpression(f"unresolved names {sorted(unknown)} in {source}")
        try:
            expr = parse_expr(
                source,
                local_dict=self.local_dict(stage_id),
                global_dict={
                    "Integer": sp.Integer,
                    "Symbol": sp.Symbol,
                    "Rational": sp.Rational,
                    "Function": sp.Function,
                    "Add": sp.Add,
                    "Mul": sp.Mul,
                    "Pow": sp.Pow,
                },
                transformations=standard_transformations,
                evaluate=False,
            )
        except Exception as exc:
            raise UnsupportedExpression(f"parse failed for {source}: {exc}") from exc
        if expr.atoms(sp.Float):
            raise UnsupportedExpression(f"float atom in {source}")
        self.parse_cache[key] = expr
        return expr

    def claim_binds_symbol(self, claim: Mapping[str, Any], alias: str, stage_id: str) -> bool:
        kind = claim.get("payload_kind")
        payload = claim.get("payload", {})
        if kind == "relation" and payload.get("relation") in {"eq", "defines"}:
            lhs = payload.get("lhs", {}).get("sympy")
            if not isinstance(lhs, str):
                return False
            try:
                parsed = self.parse(stage_id, lhs)
            except UnsupportedExpression:
                return False
            if parsed == self.local_dict(stage_id).get(alias):
                return True
            if isinstance(parsed, AppliedUndef) and parsed.func.__name__ == alias:
                return True
            return False
        if kind == "spectrum":
            return payload.get("kernel_symbol") == alias
        if kind == "operator_identity":
            return payload.get("acts_on") == alias
        return False

    def check_evidence(self) -> None:
        seen: set[tuple[str, str]] = set()
        mathematica_output_mode = any(
            isinstance(
                manifest.get("verification", {}).get("mathematica_output"),
                dict,
            )
            for manifest in self.manifests
        )
        for manifest in self.manifests:
            stage = manifest.get("stage_id", "?")
            for evidence in all_evidence_objects(manifest):
                path_value = evidence.get("source_path")
                digest_value = evidence.get("source_digest")
                if not isinstance(path_value, str) or not isinstance(digest_value, str):
                    continue
                pair = (path_value, digest_value)
                if pair in seen:
                    continue
                seen.add(pair)
                path = resolve_path(path_value, self.roots)
                if path is None:
                    self.results["EVIDENCE"].add("FAIL", "MISSING_SOURCE", f"{stage}:{path_value}")
                    continue
                algorithm, expected = digest_hex(digest_value)
                actual = hash_file(path, algorithm)
                if actual != expected:
                    self.results["EVIDENCE"].add(
                        "FAIL", "STALE_SOURCE_DIGEST", f"{stage}:{path_value}: expected {expected}, got {actual}"
                    )
            verification = manifest.get("verification", {})
            script_value = verification.get("mathematica_script")
            output = verification.get("mathematica_output")
            if (
                mathematica_output_mode
                and isinstance(script_value, str)
                and script_value
                and not isinstance(output, dict)
            ):
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "MISSING_MATHEMATICA_OUTPUT_CITATION",
                    f"{stage}:{script_value}",
                )
                continue
            if not isinstance(output, dict):
                continue
            output_path = output.get("path")
            expected_sha256 = output.get("sha256")
            if not isinstance(script_value, str) or not script_value:
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "MATHEMATICA_OUTPUT_WITHOUT_SCRIPT",
                    f"{stage}:{output_path}",
                )
                continue
            if not safe_relative_path(script_value):
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "UNSAFE_MATHEMATICA_SCRIPT_PATH",
                    f"{stage}:{script_value}",
                )
                continue
            script_path = resolve_path(script_value, self.roots)
            if script_path is None:
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "MISSING_MATHEMATICA_SCRIPT",
                    f"{stage}:{script_value}",
                )
                continue
            if not isinstance(output_path, str) or not safe_relative_path(output_path):
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "UNSAFE_MATHEMATICA_OUTPUT_PATH",
                    f"{stage}:{output_path}",
                )
                continue
            expected_output = expected_mathematica_output_path(script_value)
            if Path(output_path) != expected_output:
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "MATHEMATICA_OUTPUT_SCRIPT_MISMATCH",
                    f"{stage}:{output_path} is not {expected_output} for {script_value}",
                )
                continue
            path = resolve_path(output_path, self.roots)
            if path is None:
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "MISSING_MATHEMATICA_OUTPUT",
                    f"{stage}:{output_path}",
                )
                continue
            data = path.read_bytes()
            if not data:
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "EMPTY_MATHEMATICA_OUTPUT",
                    f"{stage}:{output_path}",
                )
                continue
            actual_sha256 = hashlib.sha256(data).hexdigest()
            if not isinstance(expected_sha256, str) or actual_sha256 != expected_sha256:
                self.results["EVIDENCE"].add(
                    "FAIL",
                    "STALE_MATHEMATICA_OUTPUT_DIGEST",
                    f"{stage}:{output_path}: expected {expected_sha256}, got {actual_sha256}",
                )
                continue
            self.coverage.mathematica_outputs_checked += 1

    def check_adjudications(self) -> None:
        for ref, claim in self.claims.items():
            if claim.get("payload_kind") != "adjudication":
                continue
            payload = claim.get("payload", {})
            cardinality = payload.get("domain_cardinality")
            counts = payload.get("bucket_counts")
            axes = payload.get("axes")
            if isinstance(counts, dict) and sum(counts.values()) != cardinality:
                self.results["ADJUDICATION"].add(
                    "FAIL", "BUCKET_COUNT_SUM", f"{ref}: {sum(counts.values())} != {cardinality}"
                )
            if isinstance(axes, dict):
                product = 1
                for size in axes.values():
                    product *= size
                if product != cardinality:
                    self.results["ADJUDICATION"].add(
                        "FAIL", "AXIS_CARDINALITY_PRODUCT", f"{ref}: {product} != {cardinality}"
                    )

    def payload_expression_strings(self, typed: Mapping[str, Any]) -> Iterator[str]:
        kind = typed.get("payload_kind")
        payload = typed.get("payload", {})
        if kind == "relation":
            for side in ("lhs", "rhs"):
                value = payload.get(side, {}).get("sympy")
                if isinstance(value, str):
                    yield value
        elif kind == "operator_identity":
            value = payload.get("definition", {}).get("sympy")
            if isinstance(value, str):
                yield value
            for key in ("operator", "acts_on"):
                value = payload.get(key)
                if isinstance(value, str):
                    yield value
        elif kind == "spectrum":
            for key in ("operator", "kernel_symbol"):
                value = payload.get(key)
                if isinstance(value, str):
                    yield value
            for key in ("kernel", "eigenvalue", "gap"):
                value = payload.get(key, {}).get("sympy")
                if isinstance(value, str):
                    yield value

    def condition_expression_strings(self, condition: Any) -> Iterator[str]:
        if not isinstance(condition, dict):
            return
        for operand in condition.get("operands", []):
            if isinstance(operand, dict) and "variables" in operand:
                yield from (v for v in operand.get("variables", []) if isinstance(v, str))
                predicate = operand.get("predicate", {}).get("sympy")
                if isinstance(predicate, str):
                    yield predicate
            elif isinstance(operand, dict):
                yield from self.condition_expression_strings(operand)

    def expression_aliases(self, stage_id: str, source: str) -> set[str]:
        table = self.symbols_by_stage_alias.get(stage_id, {})
        if source in table:
            return {source}
        expr = self.parse(stage_id, source)
        names = {symbol.name for symbol in expr.free_symbols}
        names |= {node.func.__name__ for node in expr.atoms(AppliedUndef)}
        return names & set(table)

    def check_import_completeness(self) -> None:
        result = self.results["IMPORT"]
        # Ownership and definition binding are checked before the import walk.
        for qid, appearances in self.symbols_by_qid.items():
            owner = self.owners.get(qid)
            local_valid: list[tuple[str, dict[str, Any]]] = []
            for stage_id, symbol in appearances:
                definition = symbol.get("definition_ref", "")
                if isinstance(definition, str) and definition.startswith("here/"):
                    target = self.claims.get(claim_ref(stage_id, local_ref_id(definition)))
                    if not target or not self.claim_binds_symbol(target, symbol.get("parse_alias", ""), stage_id):
                        result.add(
                            "FAIL",
                            "DEFINITION_NOT_BINDING",
                            f"{stage_id}:{qid}:{definition} does not bind {symbol.get('parse_alias')}",
                        )
                    else:
                        local_valid.append((stage_id, symbol))
                elif isinstance(definition, str):
                    target = self.claims.get(definition)
                    if target:
                        producer_stage = definition.split("/", 1)[0]
                        producer_symbols = [
                            s for sid, s in appearances if sid == producer_stage and s.get("definition_ref", "").startswith("here/")
                        ]
                        if not producer_symbols or not any(
                            self.claim_binds_symbol(target, ps.get("parse_alias", ""), producer_stage)
                            for ps in producer_symbols
                        ):
                            result.add(
                                "FAIL",
                                "DEFINITION_NOT_BINDING",
                                f"{stage_id}:{qid}:{definition} does not bind producer quantity",
                            )
            if owner:
                owner_stage, owner_claim, _ = owner
                expected = claim_ref(owner_stage, owner_claim)
                for stage_id, symbol in appearances:
                    definition = symbol.get("definition_ref")
                    if stage_id != owner_stage and isinstance(definition, str) and definition.startswith("here/"):
                        result.add(
                            "FAIL",
                            "FALSE_LOCAL_DEFINITION",
                            f"{stage_id}:{qid}: owner is {expected}, declared {definition}",
                        )
                    elif stage_id != owner_stage and definition != expected:
                        result.add(
                            "FAIL",
                            "WRONG_OWNER_REFERENCE",
                            f"{stage_id}:{qid}: expected {expected}, got {definition}",
                        )
            elif local_valid:
                # build_indexes would have found a valid owner; defensive only.
                result.add("FAIL", "OWNER_INDEX_ERROR", qid)

        for stage_id, manifest in self.stages.items():
            used_aliases: set[str] = set()
            sources: list[str] = []
            for claim in manifest.get("claims", []):
                sources.extend(self.payload_expression_strings(claim))
                sources.extend(self.condition_expression_strings(claim.get("holds_within")))
            for consume in manifest.get("consumes", []):
                if "as_consumed" in consume:
                    sources.extend(self.payload_expression_strings(consume["as_consumed"]))
                for substitution in consume.get("substitutions", []):
                    for side in ("lhs", "rhs"):
                        expr = substitution.get(side, {}).get("sympy")
                        if isinstance(expr, str):
                            sources.append(expr)
                domain = consume.get("domain")
                if isinstance(domain, dict):
                    sources.extend(v for v in domain.get("variables", []) if isinstance(v, str))
                    predicate = domain.get("predicate", {}).get("sympy")
                    if isinstance(predicate, str):
                        sources.append(predicate)
            for source in sources:
                try:
                    used_aliases |= self.expression_aliases(stage_id, source)
                except UnsupportedExpression as exc:
                    result.add("FAIL", "UNRESOLVED_EXPRESSION", f"{stage_id}: {exc}")
            consume_refs = {consume.get("ref") for consume in manifest.get("consumes", [])}
            for alias in sorted(used_aliases):
                symbol = self.symbols_by_stage_alias[stage_id][alias]
                qid = symbol["quantity_id"]
                owner = self.owners.get(qid)
                definition = symbol.get("definition_ref")
                if owner:
                    definition = claim_ref(owner[0], owner[1])
                if (
                    isinstance(definition, str)
                    and not definition.startswith("here/")
                    and definition.split("/", 1)[0] != stage_id
                ):
                    if definition not in consume_refs:
                        result.add(
                            "FAIL",
                            "UNDECLARED_IMPORT",
                            f"{stage_id}:{qid}:{definition}",
                        )

    def check_c1(self) -> None:
        result = self.results["C1"]
        for qid, appearances in self.symbols_by_qid.items():
            baseline: tuple[Any, ...] | None = None
            baseline_stage: str | None = None
            baseline_basis: DimensionBasis | None = None
            for stage_id, symbol in appearances:
                signature = symbol.get("signature")
                current = (
                    self.normalized_stage_dim(stage_id, symbol.get("dim")),
                    tuple(sorted(symbol.get("assumptions", []))),
                    json.dumps(signature, sort_keys=True),
                )
                if baseline is None:
                    baseline = current
                    baseline_stage = stage_id
                    baseline_basis = self.stage_basis(stage_id)
                elif self.stage_basis(stage_id) != baseline_basis:
                    self.results["DIMENSIONAL_CONSISTENCY"].add(
                        "FAIL",
                        "CROSS_BASIS_DIMENSION_COMPARISON",
                        f"{qid}: {baseline_stage} uses "
                        f"{baseline_basis.render() if baseline_basis else '<unknown>'}; "
                        f"{stage_id} uses {self.stage_basis(stage_id).render()}",
                    )
                    if current[1:] != baseline[1:]:
                        result.add(
                            "FAIL",
                            "QUANTITY_ID_CONFLICT",
                            f"{qid} differs at {stage_id}",
                        )
                elif current != baseline:
                    result.add("FAIL", "QUANTITY_ID_CONFLICT", f"{qid} differs at {stage_id}")
        by_name: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for manifest in self.manifests:
            for symbol in manifest.get("symbols", []):
                by_name[symbol.get("name", "")].append(symbol)
        for name, entries in by_name.items():
            qids = {entry.get("quantity_id") for entry in entries}
            if len(qids) <= 1:
                continue
            for entry in entries:
                declaration = entry.get("alias", {})
                if set(declaration.get("quantity_ids", [])) != qids:
                    result.add("FAIL", "SAME_NAME_WITHOUT_ALIAS", f"{name}: expected {sorted(qids)}")
                    break

    def apply_substitutions(
        self, stage_id: str, expr: sp.Expr, substitutions: Sequence[Mapping[str, Any]]
    ) -> sp.Expr:
        replacements: dict[sp.Expr, sp.Expr] = {}
        for substitution in substitutions:
            lhs = self.parse(stage_id, substitution["lhs"]["sympy"])
            rhs = self.parse(stage_id, substitution["rhs"]["sympy"])
            replacements[lhs] = rhs
        return expr.xreplace(replacements)

    def relation_residual(
        self, stage_id: str, payload: Mapping[str, Any], substitutions: Sequence[Mapping[str, Any]] = ()
    ) -> sp.Expr:
        lhs = self.parse(stage_id, payload["lhs"]["sympy"])
        rhs = self.parse(stage_id, payload["rhs"]["sympy"])
        return sp.simplify(self.apply_substitutions(stage_id, lhs - rhs, substitutions))

    def expressions_equal(
        self,
        producer_stage: str,
        producer_source: str,
        consumer_stage: str,
        consumer_source: str,
        substitutions: Sequence[Mapping[str, Any]],
    ) -> bool:
        producer_expr = self.parse(producer_stage, producer_source)
        consumer_expr = self.parse(consumer_stage, consumer_source)
        # Substitution names are consumer-local; compare by stable name.
        consumer_expr = self.apply_substitutions(consumer_stage, consumer_expr, substitutions)
        producer_expr = producer_expr.xreplace(
            {s: sp.Symbol(s.name) for s in producer_expr.free_symbols}
        )
        consumer_expr = consumer_expr.xreplace(
            {s: sp.Symbol(s.name) for s in consumer_expr.free_symbols}
        )
        return sp.simplify(producer_expr - consumer_expr) == 0

    def validate_substitutions(self, stage_id: str, consume: Mapping[str, Any], result: CheckResult) -> bool:
        valid = True
        for substitution in consume.get("substitutions", []):
            backed_by = substitution.get("backed_by")
            claim = self.claims.get(backed_by)
            if (
                not claim
                or claim.get("kind") != "convention"
                or backed_by not in self.exports
                or claim.get("status") == "RETIRED"
            ):
                result.add("FAIL", "UNBACKED_SUBSTITUTION", f"{stage_id}:{backed_by}")
                valid = False
        return valid

    def check_c2(self) -> None:
        result = self.results["C2"]
        for consumer_stage, manifest in self.stages.items():
            for consume in manifest.get("consumes", []):
                ref = consume.get("ref")
                producer_claim = self.claims.get(ref)
                if not producer_claim:
                    continue
                producer_stage = ref.split("/", 1)[0]
                mode = consume.get("check")
                producer_kind = producer_claim.get("payload_kind")
                consumer_kind = consume.get("as_consumed", {}).get("payload_kind")
                if not self.validate_substitutions(consumer_stage, consume, result):
                    continue
                bound_quantities = {
                    qid
                    for qid, owner in self.owners.items()
                    if claim_ref(owner[0], owner[1]) == ref
                }
                if mode == "token_match" and (
                    consume.get("producer_quantity_id")
                    or bound_quantities
                    or any(
                        symbol.get("parse_alias") == consume.get("as_consumed", {}).get("payload", {}).get("token")
                        for symbol in manifest.get("symbols", [])
                    )
                ):
                    result.add("FAIL", "TOKEN_MATCH_QUANTITY", f"{consumer_stage}:{ref}")
                    continue
                expected_kinds = {
                    "cas_equivalence": "relation",
                    "implication": "relation",
                    "specialization": "relation",
                    "value_equal": "relation",
                    "token_match": "token",
                    "spectrum_match": "spectrum",
                    "range_match": "record_range",
                    "adjudication_match": "adjudication",
                    "set_match": "set_cardinality",
                }
                if mode in expected_kinds and (
                    producer_kind != expected_kinds[mode] or consumer_kind != expected_kinds[mode]
                ):
                    result.add(
                        "UNSUPPORTED",
                        "INELIGIBLE_PAYLOAD_PAIR",
                        f"{consumer_stage}:{ref}: {mode} on {producer_kind}->{consumer_kind}",
                    )
                    continue
                try:
                    if mode in {"cas_equivalence", "implication", "specialization", "value_equal"}:
                        producer_payload = producer_claim["payload"]
                        consumer_payload = consume["as_consumed"]["payload"]
                        p = self.relation_residual(producer_stage, producer_payload)
                        c = self.relation_residual(
                            consumer_stage, consumer_payload, consume.get("substitutions", [])
                        )
                        p = p.xreplace({s: sp.Symbol(s.name) for s in p.free_symbols})
                        c = c.xreplace({s: sp.Symbol(s.name) for s in c.free_symbols})
                        equivalent = sp.simplify(p - c) == 0 or sp.simplify(p + c) == 0
                        if mode == "value_equal" and (p.free_symbols or c.free_symbols):
                            result.add("UNSUPPORTED", "NONNUMERIC_VALUE_EQUAL", f"{consumer_stage}:{ref}")
                        elif mode == "specialization" and (
                            not consume.get("specialization") or "domain" not in consume
                        ):
                            result.add("FAIL", "SPECIALIZATION_WITHOUT_DOMAIN", f"{consumer_stage}:{ref}")
                        elif not equivalent:
                            result.add("FAIL", "CITATION_DRIFT", f"{consumer_stage}:{ref}")
                    elif mode == "dim_equal":
                        producer_basis = self.stage_basis(producer_stage)
                        consumer_basis = self.stage_basis(consumer_stage)
                        if producer_basis != consumer_basis:
                            self.results["DIMENSIONAL_CONSISTENCY"].add(
                                "FAIL",
                                "CROSS_BASIS_DIMENSION_COMPARISON",
                                f"{consumer_stage}:{ref}: producer uses "
                                f"{producer_basis.render()}; consumer uses "
                                f"{consumer_basis.render()}",
                            )
                            continue
                        producer_dims = [
                            self.normalized_stage_dim(
                                producer_stage, symbol.get("dim")
                            )
                            for symbol in self.stages[producer_stage].get("symbols", [])
                            if symbol.get("definition_ref") in {ref, f"here/{ref.split('/', 1)[1]}"}
                        ]
                        if (
                            not producer_dims
                            or self.normalized_stage_dim(
                                consumer_stage, consume["as_consumed_dim"]
                            )
                            not in producer_dims
                        ):
                            result.add("FAIL", "DIM_CITATION_DRIFT", f"{consumer_stage}:{ref}")
                    elif mode == "token_match":
                        if producer_claim["payload"] != consume["as_consumed"]["payload"]:
                            result.add("FAIL", "TOKEN_DRIFT", f"{consumer_stage}:{ref}")
                    elif mode == "spectrum_match":
                        producer = producer_claim["payload"]
                        consumer = consume["as_consumed"]["payload"]
                        scalar_fields = ("operator", "kernel_symbol", "multiplicity")
                        if any(producer.get(key) != consumer.get(key) for key in scalar_fields):
                            result.add("FAIL", "SPECTRUM_FIELD_DRIFT", f"{consumer_stage}:{ref}")
                            continue
                        for key in ("kernel", "eigenvalue", "gap"):
                            if not self.expressions_equal(
                                producer_stage,
                                producer[key]["sympy"],
                                consumer_stage,
                                consumer[key]["sympy"],
                                consume.get("substitutions", []),
                            ):
                                result.add("FAIL", "SPECTRUM_FIELD_DRIFT", f"{consumer_stage}:{ref}:{key}")
                                break
                    elif mode == "range_match":
                        if producer_claim["payload"] != consume["as_consumed"]["payload"]:
                            result.add("FAIL", "RANGE_DRIFT", f"{consumer_stage}:{ref}")
                    elif mode == "adjudication_match":
                        producer = producer_claim["payload"]
                        consumer = consume["as_consumed"]["payload"]
                        fields = ("outcome_token", "domain_cardinality", "oracle_digest")
                        if any(producer.get(key) != consumer.get(key) for key in fields):
                            result.add("FAIL", "ADJUDICATION_DRIFT", f"{consumer_stage}:{ref}")
                    elif mode == "set_match":
                        producer = producer_claim["payload"]
                        consumer = consume["as_consumed"]["payload"]
                        if (
                            producer.get("count") != consumer.get("count")
                            or set(producer.get("elements", [])) != set(consumer.get("elements", []))
                        ):
                            result.add("FAIL", "SET_DRIFT", f"{consumer_stage}:{ref}")
                    elif mode == "opaque_quantity_match":
                        export = self.exports.get(ref)
                        qid = consume.get("producer_quantity_id")
                        owner = self.owners.get(qid)
                        if (
                            not export
                            or export.get("source_digest") != consume.get("producer_source_digest")
                            or not owner
                            or claim_ref(owner[0], owner[1]) != ref
                        ):
                            result.add("FAIL", "OPAQUE_QUANTITY_PIN_DRIFT", f"{consumer_stage}:{ref}:{qid}")
                    else:
                        result.add("UNSUPPORTED", "UNKNOWN_CHECK_MODE", f"{consumer_stage}:{mode}")
                except (UnsupportedExpression, KeyError, TypeError, ValueError) as exc:
                    result.add("UNSUPPORTED", "C2_PROOF_UNSUPPORTED", f"{consumer_stage}:{ref}: {exc}")

    def check_c3(self) -> None:
        result = self.results["C3"]
        for stage_id, manifest in self.stages.items():
            local_ids = {claim.get("id") for claim in manifest.get("claims", [])}
            for export in manifest.get("exports", []):
                cid = export.get("claim_id")
                if cid not in local_ids:
                    result.add("FAIL", "UNKNOWN_EXPORT", f"{stage_id}/{cid}")
                    continue
                claim = self.claims[claim_ref(stage_id, cid)]
                if export.get("source_digest") != claim.get("evidence", {}).get("source_digest"):
                    result.add("FAIL", "EXPORT_DIGEST_MISMATCH", f"{stage_id}/{cid}")
            teeth = {tooth.get("predicate") for tooth in manifest.get("verification", {}).get("teeth", [])}
            for consume in manifest.get("consumes", []):
                ref = consume.get("ref")
                producer = self.claims.get(ref)
                self.coverage.total_citations += 1
                if not producer:
                    self.coverage.unresolved_producers += 1
                    status = self.unresolved_reference_status(stage_id, ref)
                    result.add(status, "ABSENT_PRODUCER", f"{stage_id}:{ref}")
                    continue
                self.coverage.resolved_citations += 1
                if ref not in self.exports:
                    result.add("FAIL", "NON_EXPORTED_CLAIM", f"{stage_id}:{ref}")
                if producer.get("status") == "RETIRED" or "discharged_by" in producer:
                    result.add("FAIL", "NON_OPERATIVE_CLAIM", f"{stage_id}:{ref}")
                expect = consume.get("c7_expect")
                if expect and expect.get("expected_first_failure") not in teeth:
                    result.add(
                        "FAIL",
                        "UNKNOWN_C7_TOOTH",
                        f"{stage_id}:{ref}:{expect.get('expected_first_failure')}",
                    )

    def dimension_of(self, stage_id: str, expr: sp.Expr) -> tuple[Fraction, ...]:
        zero_dim = self.zero_dim(stage_id)
        if expr.is_Number or isinstance(expr, NumberSymbol):
            return zero_dim
        if isinstance(expr, sp.Symbol):
            symbol = self.symbols_by_stage_alias.get(stage_id, {}).get(expr.name)
            if not symbol:
                raise UnsupportedExpression(f"dimension of unknown symbol {expr}")
            return self.normalized_stage_dim(stage_id, symbol.get("dim"))
        if isinstance(expr, sp.Add):
            dims = [self.dimension_of(stage_id, arg) for arg in expr.args]
            nonzero_dims = [dim for arg, dim in zip(expr.args, dims) if arg != 0]
            if nonzero_dims and any(dim != nonzero_dims[0] for dim in nonzero_dims[1:]):
                axes = self.stage_basis(stage_id).axes
                raise ValueError(
                    f"addition dimension mismatch {expr}: "
                    f"{[dim_json(dim, axes) for dim in dims]}"
                )
            return nonzero_dims[0] if nonzero_dims else zero_dim
        if isinstance(expr, sp.Mul):
            total = zero_dim
            for arg in expr.args:
                total = add_dim(total, self.dimension_of(stage_id, arg))
            return total
        if isinstance(expr, sp.Pow):
            base, power = expr.args
            if not power.is_Rational:
                raise UnsupportedExpression(f"non-rational power {expr}")
            return scale_dim(self.dimension_of(stage_id, base), Fraction(int(power.p), int(power.q)))
        if isinstance(expr, sp.Derivative):
            total = self.dimension_of(stage_id, expr.expr)
            for variable, count in expr.variable_count:
                total = sub_dim(total, scale_dim(self.dimension_of(stage_id, variable), Fraction(count)))
            return total
        if isinstance(expr, sp.Integral):
            total = self.dimension_of(stage_id, expr.function)
            for limit in expr.limits:
                total = add_dim(total, self.dimension_of(stage_id, limit[0]))
            return total
        if isinstance(expr, sp.Sum):
            return self.dimension_of(stage_id, expr.function)
        if isinstance(expr, AppliedUndef):
            alias = expr.func.__name__
            symbol = self.symbols_by_stage_alias.get(stage_id, {}).get(alias)
            signature = symbol.get("signature") if symbol else None
            if not signature:
                raise UnsupportedExpression(f"call without signature {alias}")
            actual = [self.dimension_of(stage_id, arg) for arg in expr.args]
            expected = [
                self.normalized_stage_dim(stage_id, dim)
                for dim in signature["domain"]
            ]
            if actual != expected:
                raise ValueError(f"signature domain mismatch {alias}")
            return self.normalized_stage_dim(stage_id, signature["codomain"])
        if expr.func in {sp.exp, sp.log, sp.sin, sp.cos, sp.tan, sp.sinh, sp.cosh, sp.tanh, sp.atan}:
            if any(
                self.dimension_of(stage_id, arg) != zero_dim
                for arg in expr.args
            ):
                raise ValueError(f"transcendental argument is dimensional: {expr}")
            return zero_dim
        if expr.func == sp.Abs:
            return self.dimension_of(stage_id, expr.args[0])
        raise UnsupportedExpression(f"dimension rule missing for {type(expr).__name__}: {expr}")

    def recover_dim_order_and_tuples(
        self,
        path: Path,
        basis: DimensionBasis,
        extra_expressions: Sequence[str] = (),
        snapshot: DimSourceSnapshot | None = None,
    ) -> tuple[tuple[str, ...], dict[str, tuple[Fraction, ...]]]:
        try:
            snapshot = snapshot or read_dim_source_snapshot(path)
            tree = snapshot.tree
        except Exception as exc:
            raise UnsupportedExpression(f"cannot parse dim source {path}: {exc}") from exc
        field_orders: set[tuple[str, ...]] = set()
        init_orders: set[tuple[str, ...]] = set()
        doc_orders: set[tuple[str, ...]] = set()
        axis_by_attribute = {
            axis_attribute_name(axis): axis for axis in basis.axes
        }
        axis_by_casefold = {axis.casefold(): axis for axis in basis.axes}

        def resolve_axis(token: str) -> str | None:
            stripped = token.strip()
            if stripped in basis.axes:
                return stripped
            return axis_by_casefold.get(stripped.casefold())

        for child in ast.walk(tree):
            if not isinstance(child, ast.ClassDef) or child.name != "Dim":
                continue
            fields = [
                axis_by_attribute[item.target.id]
                for item in child.body
                if isinstance(item, ast.AnnAssign)
                and isinstance(item.target, ast.Name)
                and item.target.id in axis_by_attribute
            ]
            if len(fields) >= len(basis.axes):
                field_orders.add(tuple(fields[: len(basis.axes)]))
            for item in child.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == "__init__":
                    args = [
                        axis_by_attribute[arg.arg]
                        for arg in item.args.args
                        if arg.arg != "self"
                        and arg.arg in axis_by_attribute
                    ]
                    if len(args) >= len(basis.axes):
                        init_orders.add(tuple(args[: len(basis.axes)]))
            doc = ast.get_docstring(child) or ""
            for body in re.findall(r"[\[\{\(]([^\]\}\)]+)[\]\}\)]", doc):
                resolved = tuple(
                    axis
                    for token in body.split(",")
                    if (axis := resolve_axis(token)) is not None
                )
                if len(resolved) == len(basis.axes):
                    doc_orders.add(resolved)

        structural_orders = field_orders | init_orders | doc_orders
        if len(structural_orders) > 1:
            raise UnsupportedExpression(
                "Dim field/doc/init order conflict: "
                f"fields={sorted(field_orders)}, init={sorted(init_orders)}, "
                f"doc={sorted(doc_orders)}"
            )
        order = next(iter(structural_orders), None)
        source_digest = snapshot.sha256
        if order is None:
            registered_order = self.bare_tuple_dim_orders.get(source_digest)
            if registered_order is None:
                raise UnsupportedExpression(
                    "bare-tuple Dim order is not independently registered for "
                    f"source sha256 {source_digest}"
                )
            order = source_order_axes(registered_order)
        if len(order) != len(basis.axes) or set(order) != set(basis.axes):
            raise UnsupportedExpression(
                f"recovered order {list(order)} does not span basis "
                f"{list(basis.axes)}"
            )
        canonical_expressions = tuple(sorted(set(extra_expressions)))
        cache_key = (path.resolve(), source_digest, order, canonical_expressions)
        payload = LIVE_DIM_RECOVERY_CACHE.get(cache_key)
        if payload is None:
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--recover-dims",
                str(path),
                "--dim-order",
                json.dumps(order),
                "--recover-dims-stdin",
            ]
            for expression in canonical_expressions:
                command.extend(("--dim-expression", expression))
            try:
                completed = subprocess.run(
                    command,
                    cwd=path.parent,
                    check=False,
                    input=snapshot.data,
                    capture_output=True,
                    timeout=30,
                )
            except (OSError, subprocess.TimeoutExpired) as exc:
                raise UnsupportedExpression(f"live Dim recovery failed for {path}: {exc}") from exc
            marker = "DIM_RECOVERY_JSON="
            stdout = completed.stdout.decode("utf-8", errors="replace")
            stderr = completed.stderr.decode("utf-8", errors="replace")
            payload_line = next(
                (
                    line[len(marker) :]
                    for line in reversed(stdout.splitlines())
                    if line.startswith(marker)
                ),
                None,
            )
            if completed.returncode != 0 or payload_line is None:
                detail = stderr.strip() or stdout.strip()
                raise UnsupportedExpression(
                    f"live Dim recovery failed for {path}: exit={completed.returncode}: {detail}"
                )
            try:
                payload = json.loads(payload_line)
            except json.JSONDecodeError as exc:
                raise UnsupportedExpression(
                    f"invalid live Dim recovery JSON for {path}: {exc}"
                ) from exc
            LIVE_DIM_RECOVERY_CACHE[cache_key] = payload
        try:
            tuples = {
                name: tuple(Fraction(value) for value in raw)
                for name, raw in payload["tuples"].items()
            }
            self.dim_recovery_unsupported[path] = set(payload.get("unsupported", []))
            self.dim_recovery_ambiguous[path] = set(payload.get("ambiguous", []))
            self.dim_recovery_binding_counts[path] = {
                name: int(count)
                for name, count in payload.get("binding_counts", {}).items()
            }
            self.dim_recovery_binding_expressions[path] = {
                name: [str(expression) for expression in expressions]
                for name, expressions in payload.get(
                    "binding_expressions", {}
                ).items()
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise UnsupportedExpression(f"invalid live Dim recovery payload for {path}: {exc}") from exc
        return order, tuples

    def check_dimensional_consistency(self) -> None:
        result = self.results["DIMENSIONAL_CONSISTENCY"]
        reported_anchor_failures: set[tuple[str, str]] = set()
        source_cache: dict[
            Path,
            tuple[
                DimSourceSnapshot,
                tuple[str, ...],
                dict[str, tuple[Fraction, ...]],
            ],
        ] = {}
        for stage_id, manifest in self.stages.items():
            basis = self.stage_basis(stage_id)
            for claim in manifest.get("claims", []):
                if claim.get("payload_kind") != "relation":
                    continue
                payload = claim.get("payload", {})
                try:
                    lhs = self.parse(stage_id, payload["lhs"]["sympy"])
                    rhs = self.parse(stage_id, payload["rhs"]["sympy"])
                    lhs_dim = self.dimension_of(stage_id, lhs)
                    rhs_dim = self.dimension_of(stage_id, rhs)
                    if lhs != 0 and rhs != 0 and lhs_dim != rhs_dim:
                        result.add(
                            "FAIL",
                            "DIMENSIONAL_INHOMOGENEITY",
                            f"{stage_id}/{claim.get('id')}: "
                            f"{dim_json(lhs_dim, basis.axes)} != "
                            f"{dim_json(rhs_dim, basis.axes)}",
                        )
                    if "expected_dim" in claim:
                        expected = self.normalized_stage_dim(
                            stage_id, claim["expected_dim"]
                        )
                        if (lhs != 0 and lhs_dim != expected) or (rhs != 0 and rhs_dim != expected):
                            result.add(
                                "FAIL", "EXPECTED_DIM_MISMATCH", f"{stage_id}/{claim.get('id')}"
                            )
                except ValueError as exc:
                    result.add("FAIL", "DIMENSIONAL_INHOMOGENEITY", f"{stage_id}/{claim.get('id')}: {exc}")
                except (UnsupportedExpression, KeyError, TypeError) as exc:
                    result.add("UNSUPPORTED", "DIMENSION_RULE_UNSUPPORTED", f"{stage_id}/{claim.get('id')}: {exc}")

            for symbol in manifest.get("symbols", []):
                source = symbol.get("dim_source", {})
                source_value = source.get("script_path", "")
                if not safe_relative_path(source_value):
                    result.add(
                        "FAIL",
                        "UNSAFE_DIM_SOURCE_PATH",
                        f"{stage_id}:{symbol.get('parse_alias')}:{source_value}",
                    )
                    continue
                path = resolve_path(source_value, self.roots)
                if path is None:
                    result.add(
                        "UNSUPPORTED",
                        "DIM_SOURCE_MISSING",
                        f"{stage_id}:{symbol.get('parse_alias')}:{source.get('script_path')}",
                    )
                    continue
                matching_evidence: list[dict[str, Any]] = []
                for evidence in certifying_evidence_objects(manifest):
                    evidence_value = evidence.get("source_path")
                    if not safe_relative_path(evidence_value):
                        continue
                    if resolve_path(evidence_value, self.roots) == path:
                        matching_evidence.append(evidence)
                if not matching_evidence:
                    result.add(
                        "FAIL",
                        "DIM_SOURCE_NOT_STAGE_EVIDENCE",
                        f"{stage_id}:{symbol.get('parse_alias')}:{source_value}",
                    )
                    continue
                if not admissible_dim_anchor(
                    stage_id, source_value, path, self.roots
                ):
                    anchor_key = (
                        stage_id,
                        normalized_relative_path(source_value),
                    )
                    if anchor_key not in reported_anchor_failures:
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_NOT_REGISTERED",
                            f"{stage_id}:{anchor_key[1]}",
                        )
                        reported_anchor_failures.add(anchor_key)
                    continue
                try:
                    if path not in source_cache:
                        snapshot = read_dim_source_snapshot(path)
                    else:
                        snapshot = source_cache[path][0]
                    digest_backed = False
                    for evidence in matching_evidence:
                        try:
                            algorithm, expected = digest_hex(
                                evidence.get("source_digest", "")
                            )
                            actual = hashlib.new(algorithm, snapshot.data).hexdigest()
                        except (KeyError, TypeError, ValueError):
                            continue
                        if actual == expected:
                            digest_backed = True
                            break
                    if not digest_backed:
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_EVIDENCE_DIGEST_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}:{source_value}",
                        )
                        continue
                    if path not in source_cache:
                        recovered_order, tuples = self.recover_dim_order_and_tuples(
                            path, basis, snapshot=snapshot
                        )
                        source_cache[path] = snapshot, recovered_order, tuples
                    _, recovered_order, tuples = source_cache[path]
                    declared_order = source_order_axes(
                        symbol.get("dim_source_order")
                    )
                    if recovered_order != declared_order:
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_ORDER_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}: "
                            f"recovered {list(recovered_order)}, declared "
                            f"{list(declared_order)}",
                        )
                    declared_tuple = tuple(Fraction(value) for value in symbol.get("dim_source_tuple", []))
                    locus = source.get("locus", "")
                    alias = symbol.get("parse_alias")
                    name = symbol.get("name")
                    locus_tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", locus)
                    unsupported_names = self.dim_recovery_unsupported.get(path, set())
                    ambiguous_names = self.dim_recovery_ambiguous.get(path, set())
                    binding_counts = self.dim_recovery_binding_counts.get(path, {})
                    binding_expressions = self.dim_recovery_binding_expressions.get(
                        path, {}
                    )

                    def binding_keys(target: str) -> list[str]:
                        scoped = [
                            f"{scope}.{target}"
                            for scope in locus_tokens
                            if f"{scope}.{target}" in binding_counts
                        ]
                        return list(dict.fromkeys(scoped + [target]))

                    def resolve_binding(
                        target: str, expected_expression: str | None = None
                    ) -> tuple[Fraction, ...] | None:
                        for key in binding_keys(target):
                            count = binding_counts.get(key, 0)
                            if count == 0:
                                continue
                            if count != 1:
                                raise UnsupportedExpression(
                                    f"ambiguous real Name assignments for {key}"
                                )
                            if key in unsupported_names or key in ambiguous_names:
                                raise UnsupportedExpression(
                                    f"live Dim assignment unsupported or ambiguous for {key}"
                                )
                            if expected_expression is not None and binding_expressions.get(
                                key
                            ) != [expected_expression]:
                                continue
                            if key not in tuples:
                                raise UnsupportedExpression(
                                    f"real Name assignment {key} is not Dim-valued"
                                )
                            return tuples[key]
                        return None

                    source_tuple: tuple[Fraction, ...] | None
                    expression_source = locus_dimension_expression(locus)
                    if expression_source is not None:
                        expression_tree = ast.parse(
                            expression_source, mode="eval"
                        ).body
                        contains_constructor = any(
                            isinstance(node, ast.Call)
                            and isinstance(node.func, ast.Name)
                            and node.func.id == "Dim"
                            for node in ast.walk(expression_tree)
                        )
                        rooted_names = {
                            node.id
                            for node in ast.walk(expression_tree)
                            if isinstance(node, ast.Name)
                            and any(
                                binding_counts.get(key, 0) > 0
                                for key in binding_keys(node.id)
                            )
                        }
                        if contains_constructor or not rooted_names:
                            result.add(
                                "FAIL",
                                "DIM_SOURCE_EXPRESSION_NOT_BINDING_ROOTED",
                                f"{stage_id}:{alias}:{expression_source}",
                            )
                            continue
                        expression_key = "expr:" + expression_source
                        if expression_key not in tuples:
                            _, expression_tuples = self.recover_dim_order_and_tuples(
                                path,
                                basis,
                                (expression_source,),
                                snapshot=snapshot,
                            )
                            tuples.update(expression_tuples)
                            unsupported_names = self.dim_recovery_unsupported.get(
                                path, set()
                            )
                            ambiguous_names = self.dim_recovery_ambiguous.get(
                                path, set()
                            )
                        if (
                            expression_key in unsupported_names
                            or expression_key in ambiguous_names
                            or expression_key not in tuples
                        ):
                            raise UnsupportedExpression(
                                "locus-anchored dimension expression is absent, "
                                f"unsupported, or ambiguous: {expression_source}"
                            )
                        source_tuple = tuples[expression_key]
                    else:
                        assignment_anchors = locus_assignment_anchors(locus)
                        if assignment_anchors:
                            anchored: list[tuple[Fraction, ...]] = []
                            for target_name, expected_expression in assignment_anchors:
                                if expected_expression is None:
                                    continue
                                recovered = resolve_binding(
                                    target_name, expected_expression
                                )
                                if recovered is not None:
                                    anchored.append(recovered)
                            if len(anchored) != 1:
                                raise UnsupportedExpression(
                                    "locus assignment does not resolve to exactly "
                                    f"one real Dim-valued Name binding: {locus}"
                                )
                            source_tuple = anchored[0]
                        else:
                            source_tuple = None
                            fallback_targets = [
                                target
                                for target in (alias, name)
                                if isinstance(target, str)
                            ]
                            for target_name in dict.fromkeys(fallback_targets):
                                recovered = resolve_binding(target_name)
                                if recovered is not None:
                                    source_tuple = recovered
                                    break
                            if source_tuple is None:
                                raise UnsupportedExpression(
                                    "bare alias does not resolve to a unique real "
                                    f"Dim-valued Name assignment at locus {locus}"
                                )
                    if source_tuple != declared_tuple:
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_TUPLE_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}: source {source_tuple}, declared {declared_tuple}",
                        )
                    transposed = {
                        axis: value
                        for axis, value in zip(
                            recovered_order, declared_tuple, strict=True
                        )
                    }
                    named = tuple(transposed[axis] for axis in basis.axes)
                    if named != self.normalized_stage_dim(
                        stage_id, symbol.get("dim")
                    ):
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_NAMED_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}: "
                            f"{dim_json(named, basis.axes)} != "
                            f"{symbol.get('dim')}",
                        )
                except (
                    UnsupportedExpression,
                    OSError,
                    SyntaxError,
                    UnicodeError,
                    ValueError,
                ) as exc:
                    result.add(
                        "UNSUPPORTED", "DIM_SOURCE_UNSUPPORTED", f"{stage_id}:{symbol.get('parse_alias')}: {exc}"
                    )

    @staticmethod
    def registry_owner(enters: str) -> str | None:
        explicit_owner = re.search(
            r"\bOWNED\s+by\b.*?\bstage\s*([0-9]{3})\b",
            enters,
            re.I,
        )
        if explicit_owner is not None:
            return f"stage{explicit_owner.group(1)}"
        explicit_stage = re.search(r"\bstage\s*([0-9]{3})\b", enters, re.I)
        if explicit_stage is not None:
            return f"stage{explicit_stage.group(1)}"
        parenthesized = re.search(r"\([^)]*?\b([0-9]{3})\b", enters)
        if parenthesized is not None:
            return f"stage{parenthesized.group(1)}"
        return None

    def registry_rows_and_owners(
        self, path: Path, format_name: str
    ) -> tuple[set[str], dict[str, str | None], dict[str, str | None]]:
        lines = path.read_text().splitlines()
        if format_name == "plain":
            rows = {
                line.strip()
                for line in lines
                if line.strip() and not line.lstrip().startswith("#")
            }
            return rows, {row: None for row in rows}, {row: None for row in rows}
        if format_name == "markdown_master_table":
            rows: set[str] = set()
            owners: dict[str, str | None] = {}
            row_classes: dict[str, str | None] = {}
            in_table = False
            for line in lines:
                if line.startswith("| Param |"):
                    in_table = True
                    continue
                if in_table and line.startswith("|---"):
                    continue
                if in_table and not line.startswith("|"):
                    if rows:
                        break
                    continue
                if in_table:
                    cells = line.split("|")
                    if len(cells) < 4:
                        continue
                    first = cells[1]
                    tokens = re.findall(r"`([^`]+)`", first)
                    if tokens:
                        row = tokens[0]
                        rows.add(row)
                        owners[row] = self.registry_owner(cells[3])
                        row_classes[row] = cells[4].strip() if len(cells) > 4 else None
            return rows, owners, row_classes
        raise ValueError(f"unknown parameter register format {format_name}")

    def registry_rows(self, path: Path, format_name: str) -> set[str]:
        rows, _, _ = self.registry_rows_and_owners(path, format_name)
        return rows

    @staticmethod
    def registry_row_kind(row_class: str | None) -> str:
        """Classify a register row for the knob-lifecycle census."""

        if row_class is None:
            return "knob"
        if "~~" in row_class:
            return "non_knob"
        upper = row_class.upper()
        markers: list[tuple[int, str]] = []
        for token, kind in (
            ("GAP", "pending"),
            ("DERIVED", "non_knob"),
            ("CANDIDATE", "non_knob"),
            ("CONV", "non_knob"),
            ("ACTION", "knob"),
            ("FREE-UNREDUCED", "knob"),
            ("CALIB", "knob"),
            ("IMPOSED", "knob"),
        ):
            match = re.search(rf"\b{re.escape(token)}\b", upper)
            if match is not None:
                markers.append((match.start(), kind))
        if not markers:
            return "non_knob"
        _, kind = min(markers)
        if kind == "pending":
            return "pending"
        if kind == "knob":
            return "knob"
        return "non_knob"

    def row_is_loaded_non_parameter_symbol(
        self, row: str, owner: str | None
    ) -> bool:
        if owner not in self.stages:
            return False
        non_parameter_roles = {"field", "function", "operator", "index"}
        for symbol in self.stages[owner].get("symbols", []):
            identifiers = {
                symbol.get("name"),
                symbol.get("parse_alias"),
                symbol.get("latex"),
            }
            if row in identifiers and symbol.get("role") in non_parameter_roles:
                return True
        return False

    def evaluate_record_range(self, payload: Mapping[str, Any]) -> tuple[int, int, int]:
        base_low = sum(component["low"] for component in payload.get("components", {}).values())
        base_high = sum(component["high"] for component in payload.get("components", {}).values())
        axes = payload.get("convention_axes", [])
        choices = [axis.get("choices", []) for axis in axes]
        branches = list(itertools.product(*choices)) if choices else [()]
        lows = [base_low + sum(choice["low_delta"] for choice in branch) for branch in branches]
        highs = [base_high + sum(choice["high_delta"] for choice in branch) for branch in branches]
        low = min(lows)
        high = max(highs)
        return low, high, high - low

    def check_c5(self) -> None:
        result = self.results["C5"]
        register = self.config.get("parameter_register", {})
        path = resolve_path(register.get("path", ""), self.roots)
        if path is None:
            result.add("FAIL", "REGISTER_MISSING", str(register.get("path")))
            return
        actual_digest = hash_file(path)
        if actual_digest != register.get("sha256"):
            result.add(
                "FAIL",
                "REGISTER_DIGEST_CHANGED",
                f"{path}: expected {register.get('sha256')}, got {actual_digest}",
            )
            return
        try:
            rows, row_owners, row_classes = self.registry_rows_and_owners(
                path, register.get("format", "plain")
            )
        except ValueError as exc:
            result.add("UNSUPPORTED", "REGISTER_FORMAT", str(exc))
            return
        events: list[tuple[int, str, dict[str, Any]]] = []
        categories: dict[str, set[str]] = defaultdict(set)
        for stage_id, manifest in self.stages.items():
            for knob in manifest.get("knobs", []):
                events.append((stage_number(stage_id), stage_id, knob))
                categories[knob.get("registry_row")].add(knob.get("count_category"))
        unclassified = sorted(rows - set(categories))
        if unclassified:
            pending_rows: list[str] = []
            for row in unclassified:
                owner = row_owners.get(row)
                row_kind = self.registry_row_kind(row_classes.get(row))
                if row_kind == "non_knob":
                    continue
                if row_kind == "pending":
                    pending_rows.append(row)
                elif self.row_is_loaded_non_parameter_symbol(row, owner):
                    continue
                elif row_classes.get(row) is not None and (
                    owner is None or owner not in self.stages
                ):
                    pending_rows.append(row)
                elif (
                    not self.closed_slice
                    and owner is not None
                    and owner not in self.stages
                ):
                    pending_rows.append(row)
                else:
                    detail = row if owner is None else f"{row}:owner={owner}"
                    result.add("FAIL", "UNCLASSIFIED_REGISTER_ROW", detail)
            if pending_rows:
                result.add(
                    "PARTIAL",
                    "REGISTER_COVERAGE_PARTIAL",
                    f"{len(pending_rows)} knob/GAP register rows await classification or their owning stages",
                )
        for row in sorted(set(categories) - rows):
            result.add("FAIL", "UNKNOWN_REGISTER_ROW", row)
        for row, values in categories.items():
            if len(values) != 1:
                result.add("FAIL", "MULTIPLY_CLASSIFIED_ROW", f"{row}:{sorted(values)}")

        prior_by_knob: dict[str, list[dict[str, Any]]] = defaultdict(list)
        low = 0
        high = 0
        for _, stage_id, knob in sorted(
            events,
            key=lambda item: (
                item[0],
                item[1],
                item[2].get("knob_id", ""),
                item[2].get("action", ""),
            ),
        ):
            kid = knob.get("knob_id")
            action = knob.get("action")
            if action == "inherited" and not prior_by_knob.get(kid):
                origin = knob.get("origin")
                if self.unresolved_reference_status(stage_id, origin) == "PARTIAL":
                    result.add(
                        "PARTIAL",
                        "LIFECYCLE_PRODUCER_UNEXTRACTED",
                        f"{stage_id}:{kid}:{self.referenced_stage(stage_id, origin)}",
                    )
                else:
                    result.add("FAIL", "ORPHAN_INHERIT", f"{stage_id}:{kid}")
            if action == "inherited" and knob.get("count_effect") != {"low": 0, "high": 0}:
                result.add("FAIL", "INHERIT_DOUBLE_COUNT", f"{stage_id}:{kid}")
            if knob.get("pending") and action != "inherited":
                effect = knob.get("count_effect", {})
                if effect.get("low") == 0 and effect.get("high") == 0:
                    result.add("FAIL", "PENDING_DEBT_UNDERCOUNT", f"{stage_id}:{kid}")
            if action == "discharged":
                discharge_ref = knob.get("discharge_evidence")
                discharge = self.claims.get(
                    discharge_ref.replace("here/", f"{stage_id}/") if isinstance(discharge_ref, str) else ""
                )
                evidence = discharge.get("evidence", {}) if discharge else {}
                if (
                    not discharge
                    or discharge.get("status") != "DERIVED"
                    or evidence.get("engine") not in {"sympy", "mathematica"}
                    or evidence.get("method") == "prose_only"
                ):
                    result.add("FAIL", "INVALID_DISCHARGE", f"{stage_id}:{kid}:{discharge_ref}")
            if action != "inherited":
                effect = knob.get("count_effect", {})
                low += effect.get("low", 0)
                high += effect.get("high", 0)
            prior_by_knob[kid].append(knob)

        range_ref = self.config.get("range_claim_ref")
        range_claim = self.claims.get(range_ref)
        if not range_claim:
            status = self.unresolved_reference_status(None, range_ref)
            result.add(status, "RANGE_CLAIM_UNRESOLVED", str(range_ref))
            return
        if range_claim.get("payload_kind") != "record_range":
            result.add("FAIL", "RANGE_CLAIM_WRONG_KIND", str(range_ref))
            return
        payload = range_claim["payload"]
        try:
            evaluated = self.evaluate_record_range(payload)
        except (KeyError, TypeError, ValueError) as exc:
            result.add("UNSUPPORTED", "RANGE_EVALUATION", f"{range_ref}:{exc}")
            return
        declared = (payload.get("low"), payload.get("high"), payload.get("spread"))
        if evaluated != declared:
            result.add("FAIL", "RANGE_INTERNAL_DRIFT", f"{range_ref}: evaluated {evaluated}, declared {declared}")
        if (low, high, high - low) != declared:
            result.add(
                "FAIL",
                "RANGE_ENDPOINT_DRIFT",
                f"census {(low, high, high-low)} != {range_ref} {declared}",
            )

    def check_c6(self) -> None:
        result = self.results["C6"]
        graph: dict[str, set[str]] = {stage: set() for stage in self.stages}
        for consumer, manifest in self.stages.items():
            for consume in manifest.get("consumes", []):
                producer = consume.get("ref", "").split("/", 1)[0]
                if producer in self.stages:
                    graph[producer].add(consumer)
        state: dict[str, int] = defaultdict(int)
        stack: list[str] = []

        def visit(node: str) -> bool:
            state[node] = 1
            stack.append(node)
            for child in sorted(graph[node]):
                if state[child] == 0 and visit(child):
                    return True
                if state[child] == 1:
                    start = stack.index(child)
                    cycle = stack[start:] + [child]
                    result.add("FAIL", "DEPENDENCY_CYCLE", " -> ".join(cycle))
                    return True
            stack.pop()
            state[node] = 2
            return False

        for node in sorted(graph):
            if state[node] == 0 and visit(node):
                break

    def check_c7(self) -> None:
        result = self.results["C7"]
        # TODO(C7_MUTATOR_INERT): An advisory faithfulness canary needs an
        # explicitly declared sentinel facet/baseline contract. Guessing one
        # here would falsely warn on honest mutators that ignore unknown facets.
        edges: list[tuple[str, Mapping[str, Any], Mapping[str, Any]]] = []
        for consumer_stage, manifest in self.stages.items():
            for consume in manifest.get("consumes", []):
                export = self.exports.get(consume.get("ref"))
                if export:
                    edges.append((consumer_stage, consume, export))
        self.coverage.c7_edges_total = len(edges)
        groups: dict[str, tuple[Mapping[str, Any], list[tuple[str, Mapping[str, Any]]]]] = {}
        for consumer_stage, consume, export in edges:
            binding = export.get("c7_binding")
            expect = consume.get("c7_expect")
            if not binding or not expect or expect.get("facet_used") != binding.get("exported_facet"):
                continue
            ref = consume["ref"]
            groups.setdefault(ref, (binding, []))[1].append((consumer_stage, consume))
        covered = sum(len(group_edges) for _, group_edges in groups.values())
        self.coverage.c7_edges_run = covered
        if covered != len(edges):
            result.add("PARTIAL", "C7_EDGE_UNCOVERED", f"{covered}/{len(edges)} edges have executable C7 metadata")

        for ref, (binding, declarations) in sorted(groups.items()):
            command = binding.get("mutation_command", "")
            try:
                argv = shlex.split(command)
                if not argv:
                    raise ValueError("empty mutation command")
                env = os.environ.copy()
                env[binding["mutation_env"]] = binding["exported_facet"]
                completed = subprocess.run(
                    argv,
                    cwd=self.roots[0],
                    env=env,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
                lines = [line for line in completed.stdout.splitlines() if line.strip()]
                if completed.returncode != 0 or not lines:
                    raise ValueError(
                        f"exit={completed.returncode} stdout={completed.stdout!r} stderr={completed.stderr!r}"
                    )
                observed = json.loads(lines[-1])
                if not isinstance(observed, dict):
                    raise ValueError("mutation command result is not an object")
            except (OSError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired) as exc:
                result.add("UNSUPPORTED", "C7_RUNNER_UNSUPPORTED", f"{ref}: {exc}")
                continue
            expected_stages: set[str] = set()
            for consumer_stage, consume in declarations:
                expected_stages.add(consumer_stage)
                expected = consume["c7_expect"]["expected_first_failure"]
                actual = observed.get(consumer_stage, "PASS")
                if actual == "PASS":
                    result.add(
                        "FAIL",
                        "DECORATIVE_DEPENDENCY",
                        f"{consumer_stage} stayed green under {ref}:{binding['exported_facet']}",
                    )
                elif actual != expected:
                    result.add(
                        "FAIL",
                        "WRONG_FIRST_FAILURE",
                        f"{consumer_stage}:{ref}: expected {expected}, got {actual}",
                    )
            for stage_id, actual in observed.items():
                if stage_id in self.stages and stage_id not in expected_stages and actual != "PASS":
                    result.add(
                        "FAIL",
                        "UNDECLARED_DEPENDENCY",
                        f"{stage_id} fired {actual} under undeclared {ref}:{binding['exported_facet']}",
                    )

    def check_c8(self) -> None:
        result = self.results["C8"]
        genesis_statuses = {"POSTULATED", "CONV", "CALIBRATED"}
        for stage_id, manifest in self.stages.items():
            for claim in manifest.get("claims", []):
                cid = claim.get("id")
                ref = claim_ref(stage_id, cid)
                status = claim.get("status")
                genesis = claim.get("genesis")
                if status in genesis_statuses and not isinstance(genesis, dict):
                    result.add("FAIL", "GENESIS_MISSING", ref)
                    continue
                if not isinstance(genesis, dict):
                    continue
                origin = genesis.get("origin")
                if origin == "independent":
                    evidence = genesis.get("evidence")
                    if not isinstance(evidence, dict):
                        result.add("FAIL", "INDEPENDENCE_EVIDENCE_MISSING", ref)
                        continue
                    predates = evidence.get("predates")
                    normalized_predates = (
                        predates.replace("here/", f"{stage_id}/")
                        if isinstance(predates, str)
                        else ""
                    )
                    if normalized_predates not in self.claims:
                        outcome = self.unresolved_reference_status(stage_id, predates)
                        result.add(
                            outcome,
                            "GENESIS_PREDATES_UNRESOLVED",
                            f"{ref}:{predates}",
                        )
                    record_date = evidence.get("date")
                    introduced = genesis.get("introduced")
                    if (
                        isinstance(record_date, str)
                        and isinstance(introduced, str)
                        and record_date >= introduced
                    ):
                        result.add(
                            "FAIL",
                            "GENESIS_DATE_NOT_PRIOR",
                            f"{ref}:{record_date}>={introduced}",
                        )
                if origin in {"coordinated", "target_matched"}:
                    for coordinated_ref in genesis.get("coordinated_with", []):
                        normalized = coordinated_ref.replace(
                            "here/", f"{stage_id}/"
                        )
                        if normalized not in self.claims:
                            outcome = self.unresolved_reference_status(
                                stage_id, coordinated_ref
                            )
                            result.add(
                                outcome,
                                "GENESIS_COORDINATION_UNRESOLVED",
                                f"{ref}:{coordinated_ref}",
                            )
                discharged_by = claim.get("discharged_by")
                if not isinstance(discharged_by, str):
                    continue
                discharge = self.claims.get(discharged_by)
                if discharge is None:
                    outcome = self.unresolved_reference_status(stage_id, discharged_by)
                    result.add(
                        outcome,
                        "DISCHARGE_REF_UNRESOLVED",
                        f"{ref}:{discharged_by}",
                    )
                    continue
                discharge_stage = discharged_by.split("/", 1)[0]
                if (
                    stage_number(discharge_stage) <= stage_number(stage_id)
                    or discharge.get("status") != "DERIVED"
                ):
                    result.add(
                        "FAIL",
                        "INVALID_DISCHARGE_DIRECTION",
                        f"{ref}:{discharged_by}",
                    )


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def render_report(report: CompositeReport) -> str:
    lines = [
        "# Composite manifest report",
        "",
        f"Headline: **{report.headline}**",
        "",
        f"Coverage: `{report.matrix_line()}`",
        "",
        "Dimensional consistency scope: "
        + DIMENSIONAL_CONSISTENCY_SCOPE,
        "",
        "Declared dimensional bases (legacy omission is reported as the "
        "implicit `LMT [L, M, T]` basis):",
        "",
        "| Stage | Basis |",
        "|---|---|",
    ]
    lines.extend(
        f"| {stage_id} | `{basis.render()}` |"
        for stage_id, basis in sorted(report.dimension_bases.items())
    )
    lines.extend(
        [
        "",
        "| Check | Outcome |",
        "|---|---|",
        ]
    )
    lines.extend(f"| {name} | {report.results[name].status} |" for name in CHECK_NAMES)
    for name in CHECK_NAMES:
        result = report.results[name]
        lines.extend(["", f"## {name}", ""])
        if name == "DIMENSIONAL_CONSISTENCY":
            lines.extend([DIMENSIONAL_CONSISTENCY_SCOPE, ""])
        if not result.issues:
            lines.append("No findings.")
        else:
            lines.extend(f"- `{issue.code}` — {issue.detail}" for issue in result.issues)
    return "\n".join(lines) + "\n"


def _self_test_evidence(source: Path, digest: str, locus: str) -> dict[str, Any]:
    return {
        "source_path": str(source),
        "locus": locus,
        "source_digest": digest,
        "engine": "sympy",
        "method": "script_asserted",
    }


def _self_test_symbol(
    source: Path,
    digest: str,
    *,
    name: str,
    alias: str,
    qid: str,
    definition: str,
    dim: dict[str, str],
    raw: Sequence[str],
    assumptions: Sequence[str] = ("real",),
    role: str = "derived",
) -> dict[str, Any]:
    return {
        "name": name,
        "parse_alias": alias,
        "quantity_id": qid,
        "definition_ref": definition,
        "dim": dim,
        "dim_source_order": "LTM",
        "dim_source": {"script_path": str(source), "locus": alias},
        "dim_source_tuple": list(raw),
        "assumptions": list(assumptions),
        "role": role,
        "evidence": _self_test_evidence(source, digest, alias),
    }


def _self_test_claim(
    source: Path,
    digest: str,
    claim_id: str,
    *,
    kind: str = "identity",
    status: str = "DERIVED",
    payload_kind: str = "relation",
    payload: dict[str, Any] | None = None,
    expected_dim: dict[str, str] | None = None,
) -> dict[str, Any]:
    claim: dict[str, Any] = {
        "id": claim_id,
        "kind": kind,
        "status": status,
        "payload_kind": payload_kind,
        "payload": payload
        or {
            "lhs": {"sympy": "1"},
            "rhs": {"sympy": "1"},
            "relation": "eq",
        },
        "evidence": _self_test_evidence(source, digest, claim_id),
    }
    if expected_dim is not None:
        claim["expected_dim"] = expected_dim
    return claim


def _self_test_manifest(
    source: Path, digest: str, stage_id: str, symbols: list[dict[str, Any]], claims: list[dict[str, Any]]
) -> dict[str, Any]:
    return {
        "schema_version": "2.1",
        "stage_id": stage_id,
        "title": f"Synthetic {stage_id}",
        "part": "Self-test",
        "status_token": f"{stage_id.upper()}_SELF_TEST",
        "symbols": symbols,
        "claims": claims,
        "consumes": [],
        "exports": [],
        "knobs": [],
        "departures": [],
        "verification": {
            "sympy_script": str(source),
            "mathematica_script": "",
            "teeth": [],
        },
        "extraction": {
            "report": {
                "ungrounded": [],
                "provisional_consumes": [],
                "prose_script_mismatches": [],
                "judgment_calls": [],
            },
            "source_digests": [{"source_path": str(source), "source_digest": digest}],
        },
    }


def _relation_payload(lhs: str, rhs: str, relation: str = "eq") -> dict[str, Any]:
    return {"lhs": {"sympy": lhs}, "rhs": {"sympy": rhs}, "relation": relation}


def _typed_relation(lhs: str, rhs: str, relation: str = "eq") -> dict[str, Any]:
    return {"payload_kind": "relation", "payload": _relation_payload(lhs, rhs, relation)}


def _c7_binding(command: str, facet: str) -> dict[str, Any]:
    return {
        "producing_primitive": f"primitive.{facet}",
        "mutation_env": "C7_FACET",
        "mutation_command": command,
        "exported_facet": facet,
    }


def _c7_expect(facet: str, tooth: str) -> dict[str, Any]:
    return {
        "injection_point": f"inject.{facet}",
        "facet_used": facet,
        "expected_first_failure": tooth,
    }


def _build_self_test_fixture(root: Path) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, str]]:
    scripts_dir = root / "scripts"
    scripts_dir.mkdir()
    source_file = (
        scripts_dir / "ledger_stage001_self_test_sympy_audit.py"
    )
    source_file.write_text(
        "from dataclasses import dataclass\n"
        "@dataclass(frozen=True)\n"
        "class Dim:\n"
        "    \"\"\"Exact exponent vector in [L,T,M] order.\"\"\"\n"
        "    l: int = 0\n"
        "    t: int = 0\n"
        "    m: int = 0\n"
        "    def __add__(self, other): return Dim(self.l+other.l, self.t+other.t, self.m+other.m)\n"
        "    def __sub__(self, other): return Dim(self.l-other.l, self.t-other.t, self.m-other.m)\n"
        "    def __mul__(self, other): return Dim(self.l+other.l, self.t+other.t, self.m+other.m)\n"
        "    def __truediv__(self, other): return Dim(self.l-other.l, self.t-other.t, self.m-other.m)\n"
        "    def __pow__(self, power): return Dim(power*self.l, power*self.t, power*self.m)\n"
        "ONE = Dim()\n"
        "LENGTH = Dim(1, 0, 0)\n"
        "TIME = Dim(0, 1, 0)\n"
        "MASS = Dim(0, 0, 1)\n"
        "ENERGY = MASS * LENGTH**2 / TIME**2\n"
        "x = Dim(0, 0, 0)\n"
        "x3 = Dim(1, 0, 0)\n"
        "a_B = Dim(-2, -2, 1)\n"
        "ell = Dim(1, 0, 0)\n"
        "area = Dim(2, 0, 0)\n"
        "f0 = Dim(-1, 0, 0)\n"
        "y = Dim(0, 0, 0)\n"
        "z = Dim(0, 0, 0)\n"
        "def dim_factory(value): return value\n"
        "def run_units():\n"
        "    K_m = ENERGY * LENGTH**2\n"
        "    K_decoy = Dim(4, -2, 1)\n"
        "    m_uu = LENGTH**3 / ENERGY\n"
        "    A = ENERGY * LENGTH\n"
        "    bad = dim_factory(LENGTH)\n"
        "    dims = {'PASS_UNITS_K_M': ('dim_Km', K_m, K_m)}\n"
        "def container_only():\n"
        "    dims = {'container_K_m': ('label', LENGTH)}\n"
    )
    source_digest = hash_file(source_file)
    source = source_file.relative_to(root)
    stage2_source_file = (
        scripts_dir / "ledger_stage002_self_test_sympy_audit.py"
    )
    stage2_source_file.write_bytes(source_file.read_bytes())
    stage2_source = stage2_source_file.relative_to(root)
    stage3_source_file = (
        scripts_dir / "ledger_stage003_self_test_sympy_audit.py"
    )
    stage3_source_file.write_bytes(source_file.read_bytes())
    stage3_source = stage3_source_file.relative_to(root)
    bare_source_file = (
        scripts_dir / "ledger_stage003_bare_tuple_sympy_audit.py"
    )
    bare_source_file.write_text(
        "def run_units():\n"
        "    dim_L = (1, 0, 0)\n"
        "    dim_E = (2, -2, 1)\n"
        "    dim_A = tuple(x + y for x, y in zip(dim_E, dim_L))\n"
    )
    bare_source_digest = hash_file(bare_source_file)
    if (
        bare_source_digest
        != "e313965b1fde7437467ae2060bcfd4ca8bb92245734dac860c7ac0fde3982c21"
    ):
        raise RuntimeError("deterministic bare-tuple self-test source drifted")
    bare_source = bare_source_file.relative_to(root)
    unregistered_bare_source_file = (
        scripts_dir
        / "ledger_stage003_unregistered_bare_order_sympy_audit.py"
    )
    unregistered_bare_source_file.write_text(
        "# Deliberately absent from BARE_TUPLE_DIM_ORDER_BY_SHA256.\n"
        "def run_units():\n"
        "    dim_L = (1, 0, 0)\n"
        "    dim_E = (2, -2, 1)\n"
        "    dim_A = tuple(x + y for x, y in zip(dim_E, dim_L))\n"
    )
    unregistered_bare_source_digest = hash_file(unregistered_bare_source_file)
    unregistered_bare_source = unregistered_bare_source_file.relative_to(root)
    duplicate_bare_source_file = (
        scripts_dir / "ledger_stage003_duplicate_bare_sympy_audit.py"
    )
    duplicate_bare_source_file.write_text(
        "def run_units():\n"
        "    dim_L = (1, 0, 0)\n"
        "    dim_E = (2, -2, 1)\n"
        "    dim_A = tuple(x + y for x, y in zip(dim_E, dim_L))\n"
        "def run_units():\n"
        "    pass\n"
    )
    duplicate_bare_source_digest = hash_file(duplicate_bare_source_file)
    duplicate_bare_source = duplicate_bare_source_file.relative_to(root)
    conditional_bare_source_file = (
        scripts_dir / "ledger_stage003_conditional_bare_sympy_audit.py"
    )
    conditional_bare_source_file.write_text(
        "def run_units(flag):\n"
        "    if flag:\n"
        "        dim_A = (3, -2, 1)\n"
    )
    conditional_bare_source_digest = hash_file(conditional_bare_source_file)
    conditional_bare_source = conditional_bare_source_file.relative_to(root)
    mutated_bare_source_file = (
        scripts_dir / "ledger_stage003_mutated_bare_sympy_audit.py"
    )
    mutated_bare_source_file.write_text(
        "def run_units():\n"
        "    dim_A = [3, -2, 1]\n"
        "    dim_A[0] = 4\n"
    )
    mutated_bare_source_digest = hash_file(mutated_bare_source_file)
    mutated_bare_source = mutated_bare_source_file.relative_to(root)
    augmented_bare_source_file = (
        scripts_dir / "ledger_stage003_augmented_bare_sympy_audit.py"
    )
    augmented_bare_source_file.write_text(
        "def run_units():\n"
        "    dim_A = [3, -2, 1]\n"
        "    dim_A += [1, 0, 0]\n"
    )
    augmented_bare_source_digest = hash_file(augmented_bare_source_file)
    augmented_bare_source = augmented_bare_source_file.relative_to(root)
    four_axis_source_file = (
        scripts_dir / "ledger_stage038_four_axis_self_test_sympy_audit.py"
    )
    four_axis_source_file.write_text(
        "def run_units():\n"
        "    A_E_dim = (0, 1, 0, 1)\n"
        "    q_T_dim = (1, 0, -1, 0)\n"
        "    shared_dim = (0, 1, 0, 1)\n"
    )
    four_axis_source_digest = hash_file(four_axis_source_file)
    four_axis_source = four_axis_source_file.relative_to(root)
    four_axis_peer_source_file = (
        scripts_dir / "ledger_stage039_four_axis_self_test_sympy_audit.py"
    )
    four_axis_peer_source_file.write_text(
        "def run_units():\n"
        "    shared_dim = (0, 1, 0, 1)\n"
        "    shared_drift_dim = (0, 2, 0, 1)\n"
    )
    four_axis_peer_source_digest = hash_file(four_axis_peer_source_file)
    four_axis_peer_source = four_axis_peer_source_file.relative_to(root)
    fractional_basis_source_file = (
        scripts_dir
        / "ledger_stage042_fractional_basis_self_test_sympy_audit.py"
    )
    fractional_basis_source_file.write_text(
        "from fractions import Fraction\n"
        "def dimension_guard():\n"
        "    stiffness_dim = (Fraction(1), Fraction(0), Fraction(0))\n"
        "    charge_dim = (Fraction(1, 2), Fraction(3, 2), Fraction(-1))\n"
    )
    fractional_basis_source_digest = hash_file(
        fractional_basis_source_file
    )
    fractional_basis_source = fractional_basis_source_file.relative_to(root)
    attacker_dim_source_file = (
        scripts_dir / "ledger_stage001_untracked_attacker_sympy_audit.py"
    )
    attacker_dim_source_file.write_text(
        "from dataclasses import dataclass\n"
        "@dataclass(frozen=True)\n"
        "class Dim:\n"
        "    \"\"\"Exact exponent vector in [L,M,T] order.\"\"\"\n"
        "    l: int = 0\n"
        "    m: int = 0\n"
        "    t: int = 0\n"
        "    def __mul__(self, other): return Dim(self.l+other.l, self.m+other.m, self.t+other.t)\n"
        "LENGTH = Dim(1, 0, 0)\n"
        "ENERGY = Dim(2, 1, -2)\n"
        "def run_units():\n"
        "    a_B = Dim(-2, 1, -2)\n"
    )
    attacker_dim_source_digest = hash_file(attacker_dim_source_file)
    attacker_dim_source = attacker_dim_source_file.relative_to(root)
    untracked_symlink_file = (
        scripts_dir / "ledger_stage001_untracked_symlink_sympy_audit.py"
    )
    untracked_symlink_file.symlink_to(source_file.name)
    untracked_symlink = untracked_symlink_file.relative_to(root)
    other_stage_dim_source_file = (
        scripts_dir / "ledger_stage002_other_stage_sympy_audit.py"
    )
    other_stage_dim_source_file.write_bytes(source_file.read_bytes())
    other_stage_dim_source_digest = hash_file(other_stage_dim_source_file)
    other_stage_dim_source = other_stage_dim_source_file.relative_to(root)
    non_python_dim_source_file = (
        scripts_dir / "ledger_stage001_non_python_sympy_audit.py"
    )
    non_python_dim_source_file.write_text(
        "# Dimension notes\n\n"
        "This is legitimately cited evidence — but it is not Python source.\n"
    )
    non_python_dim_source_digest = hash_file(non_python_dim_source_file)
    non_python_dim_source = non_python_dim_source_file.relative_to(root)
    midway_dim_source_file = (
        scripts_dir / "midway_knob_audit_codimension_sympy.py"
    )
    midway_dim_source_file.write_bytes(source_file.read_bytes())
    midway_dim_source = midway_dim_source_file.relative_to(root)
    mathematica_dir = root / "mathematica"
    mathematica_out_dir = mathematica_dir / "out"
    mathematica_out_dir.mkdir(parents=True)
    mathematica_script_file = mathematica_dir / "selftest_math.wl"
    mathematica_script_file.write_text("Print[\"SELF-TEST\"]\n")
    mathematica_script = Path("mathematica/selftest_math.wl")
    mathematica_output_file = mathematica_out_dir / "selftest_math.out"
    mathematica_output_file.write_text("SELF-TEST MATHEMATICA OUTPUT\nPASS\n")
    mathematica_output = Path("mathematica/out/selftest_math.out")
    mathematica_output_digest = hash_file(mathematica_output_file)
    empty_mathematica_script_file = mathematica_dir / "empty_math.wl"
    empty_mathematica_script_file.write_text("Print[\"EMPTY-TEST\"]\n")
    empty_mathematica_script = Path("mathematica/empty_math.wl")
    empty_mathematica_output_file = mathematica_out_dir / "empty_math.out"
    empty_mathematica_output_file.write_bytes(b"")
    empty_mathematica_output = Path("mathematica/out/empty_math.out")
    empty_mathematica_output_digest = hash_file(empty_mathematica_output_file)
    other_mathematica_output_file = mathematica_out_dir / "other_stage.out"
    other_mathematica_output_file.write_text("OTHER STAGE OUTPUT\n")
    other_mathematica_output = Path("mathematica/out/other_stage.out")
    other_mathematica_output_digest = hash_file(other_mathematica_output_file)
    missing_mathematica_script_file = mathematica_dir / "missing_math.wl"
    missing_mathematica_script_file.write_text("Print[\"MISSING-TEST\"]\n")
    missing_mathematica_script = Path("mathematica/missing_math.wl")
    register_file = root / "parameter_register.txt"
    register_file.write_text("row_a\nrow_b\n")
    register = Path(register_file.name)
    register_digest = hash_file(register_file)
    session = root / "c7_session.py"
    session.write_text(
        "import argparse, json, os\n"
        "p=argparse.ArgumentParser(); p.add_argument('--decorative', action='store_true'); "
        "p.add_argument('--undeclared', action='store_true'); a=p.parse_args()\n"
        "facet=os.environ.get('C7_FACET','')\n"
        "mapping={\n"
        " 'facet_x': {'stage002':'TOOTH_X'},\n"
        " 'facet_ell': {'stage002':'TOOTH_ELL'},\n"
        " 'facet_area': {'stage002':'TOOTH_AREA'},\n"
        " 'facet_spectrum': {'stage002':'TOOTH_SPECTRUM'},\n"
        " 'facet_token': {'stage002':'TOOTH_TOKEN','stage003':'TOOTH_STAGE3'},\n"
        " 'facet_y': {'stage001':'TOOTH_CYCLE'},\n"
        " 'facet_dim_citation': {'stage039':'TOOTH_DIM_CITATION'},\n"
        "}\n"
        "out=dict(mapping.get(facet,{}))\n"
        "if a.decorative and facet=='facet_spectrum': out['stage002']='PASS'\n"
        "if a.undeclared and facet=='facet_spectrum': out['stage003']='TOOTH_STAGE3'\n"
        "print(json.dumps(out, sort_keys=True))\n"
    )
    command = f"{shlex.quote(sys.executable)} {shlex.quote(str(session))}"
    ev = lambda locus: _self_test_evidence(source, source_digest, locus)

    s1_symbols = [
        _self_test_symbol(
            source,
            source_digest,
            name="x",
            alias="x",
            qid="q.core.x",
            definition="here/define_x",
            dim={},
            raw=("0", "0", "0"),
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="a_B",
            alias="a_B",
            qid="q.core.a_b",
            definition="here/define_ab",
            dim={"L": "-2", "M": "1", "T": "-2"},
            raw=("-2", "-2", "1"),
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="ell",
            alias="ell",
            qid="q.core.ell",
            definition="here/define_ell",
            dim={"L": "1"},
            raw=("1", "0", "0"),
            assumptions=("positive", "real"),
            role="primitive",
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="area",
            alias="area",
            qid="q.core.area",
            definition="here/define_area",
            dim={"L": "2"},
            raw=("2", "0", "0"),
            assumptions=("positive", "real"),
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="f0",
            alias="f0",
            qid="q.core.f0",
            definition="here/zero_mode",
            dim={"L": "-1"},
            raw=("-1", "0", "0"),
            role="field",
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="K_m",
            alias="K_m",
            qid="q.core.k_m",
            definition="here/define_km",
            dim={"L": "4", "M": "1", "T": "-2"},
            raw=("4", "-2", "1"),
        ),
    ]
    spectrum = {
        "operator": "x",
        "kernel_symbol": "f0",
        "kernel": {"sympy": "1/ell"},
        "eigenvalue": {"sympy": "0"},
        "gap": {"sympy": "4/ell**2"},
        "multiplicity": 1,
    }
    range_payload = {
        "low": 2,
        "high": 3,
        "spread": 1,
        "convention_axes": [
            {
                "axis_id": "baseline",
                "choices": [{"token": "default", "low_delta": 0, "high_delta": 0}],
            }
        ],
        "components": {
            "row_a": {"low": 1, "high": 1},
            "row_b": {"low": 1, "high": 2},
        },
    }
    s1_claims = [
        _self_test_claim(source, source_digest, "define_x", payload=_relation_payload("x", "2", "defines")),
        _self_test_claim(source, source_digest, "define_ab", payload=_relation_payload("a_B", "a_B", "defines")),
        _self_test_claim(source, source_digest, "define_ell", payload=_relation_payload("ell", "ell", "defines")),
        _self_test_claim(
            source, source_digest, "define_area", payload=_relation_payload("area", "ell**2", "defines")
        ),
        _self_test_claim(
            source, source_digest, "define_km", payload=_relation_payload("K_m", "K_m", "defines")
        ),
        _self_test_claim(
            source,
            source_digest,
            "zero_mode",
            kind="spectral",
            status="EARNED",
            payload_kind="spectrum",
            payload=spectrum,
        ),
        _self_test_claim(
            source,
            source_digest,
            "count_range",
            kind="range",
            status="EARNED",
            payload_kind="record_range",
            payload=range_payload,
        ),
        _self_test_claim(
            source,
            source_digest,
            "grid_landing",
            kind="adjudication",
            status="EARNED",
            payload_kind="adjudication",
            payload={
                "domain_cardinality": 4,
                "precedence": ["NO_GO", "PASS"],
                "outcome_token": "PASS",
                "oracle_digest": "abcd",
                "bucket_counts": {"NO_GO": 2, "PASS": 2},
                "axes": {"bc": 2, "sign": 2},
            },
        ),
        _self_test_claim(
            source,
            source_digest,
            "finite_set",
            kind="set",
            status="EARNED",
            payload_kind="set_cardinality",
            payload={"set_label": "S", "count": 2, "elements": ["a", "b"]},
        ),
        _self_test_claim(
            source,
            source_digest,
            "status_token",
            kind="convention",
            status="CONV",
            payload_kind="token",
            payload={"token": "READY", "meaning": "A true status token."},
        ),
    ]
    s1_claims[-1]["genesis"] = {
        "introduced": "2026-07-24",
        "origin": "unknown",
        "note": "Synthetic status-token fixture.",
    }
    stage1 = _self_test_manifest(source, source_digest, "stage001", s1_symbols, s1_claims)
    stage1["verification"]["teeth"] = [
        {"predicate": "TOOTH_CYCLE", "mutation": "cycle", "claim_ids": ["define_x"], "evidence": ev("TOOTH_CYCLE")}
    ]
    export_facets = {
        "define_x": "facet_x",
        "define_ell": "facet_ell",
        "define_area": "facet_area",
        "zero_mode": "facet_spectrum",
        "status_token": "facet_token",
    }
    stage1["exports"] = [
        {
            "claim_id": cid,
            "source_digest": source_digest,
            "c7_binding": _c7_binding(command, facet),
        }
        for cid, facet in export_facets.items()
    ]
    stage1["knobs"] = [
        {
            "knob_id": "k.row_a",
            "symbol": "x",
            "registry_row": "row_a",
            "action": "introduced",
            "origin": "here/define_x",
            "effective_stage": "stage001",
            "count_effect": {"low": 1, "high": 1},
            "count_category": "base",
            "pending": False,
            "evidence": ev("row_a"),
        },
        {
            "knob_id": "k.row_b",
            "symbol": "ell",
            "registry_row": "row_b",
            "action": "introduced",
            "origin": "here/define_ell",
            "effective_stage": "stage001",
            "count_effect": {"low": 1, "high": 2},
            "count_category": "route_ful_debt",
            "pending": True,
            "evidence": ev("row_b"),
        },
    ]

    s2_symbols = [
        _self_test_symbol(
            source,
            source_digest,
            name="x",
            alias="x",
            qid="q.core.x",
            definition="stage001/define_x",
            dim={},
            raw=("0", "0", "0"),
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="ell",
            alias="ell",
            qid="q.core.ell",
            definition="stage001/define_ell",
            dim={"L": "1"},
            raw=("1", "0", "0"),
            assumptions=("positive", "real"),
            role="primitive",
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="area",
            alias="area",
            qid="q.core.area",
            definition="stage001/define_area",
            dim={"L": "2"},
            raw=("2", "0", "0"),
            assumptions=("positive", "real"),
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="f0",
            alias="f0",
            qid="q.core.f0",
            definition="stage001/zero_mode",
            dim={"L": "-1"},
            raw=("-1", "0", "0"),
            role="field",
        ),
        _self_test_symbol(
            source,
            source_digest,
            name="y",
            alias="y",
            qid="q.stage2.y",
            definition="here/define_y",
            dim={},
            raw=("0", "0", "0"),
        ),
    ]
    s2_claims = [
        _self_test_claim(source, source_digest, "define_y", payload=_relation_payload("y", "x", "defines")),
        _self_test_claim(
            source,
            source_digest,
            "rational_power",
            kind="dimensional",
            status="EARNED",
            payload=_relation_payload("sqrt(area)", "ell"),
            expected_dim={"L": "1"},
        ),
        _self_test_claim(source, source_digest, "discharge_b", payload=_relation_payload("y", "x")),
    ]
    stage2 = _self_test_manifest(source, source_digest, "stage002", s2_symbols, s2_claims)
    stage2["verification"]["teeth"] = [
        {"predicate": tooth, "mutation": tooth.lower(), "claim_ids": ["define_y"], "evidence": ev(tooth)}
        for tooth in ("TOOTH_X", "TOOTH_ELL", "TOOTH_AREA", "TOOTH_SPECTRUM", "TOOTH_TOKEN")
    ]
    stage2["consumes"] = [
        {
            "ref": "stage001/define_x",
            "as_consumed": _typed_relation("x", "2", "defines"),
            "check": "cas_equivalence",
            "substitutions": [],
            "c7_expect": _c7_expect("facet_x", "TOOTH_X"),
        },
        {
            "ref": "stage001/define_ell",
            "as_consumed": _typed_relation("ell", "ell", "defines"),
            "check": "cas_equivalence",
            "substitutions": [],
            "c7_expect": _c7_expect("facet_ell", "TOOTH_ELL"),
        },
        {
            "ref": "stage001/define_area",
            "as_consumed": _typed_relation("area", "ell**2", "defines"),
            "check": "cas_equivalence",
            "substitutions": [],
            "c7_expect": _c7_expect("facet_area", "TOOTH_AREA"),
        },
        {
            "ref": "stage001/zero_mode",
            "as_consumed": {"payload_kind": "spectrum", "payload": copy.deepcopy(spectrum)},
            "check": "spectrum_match",
            "substitutions": [],
            "c7_expect": _c7_expect("facet_spectrum", "TOOTH_SPECTRUM"),
        },
        {
            "ref": "stage001/status_token",
            "as_consumed": {
                "payload_kind": "token",
                "payload": {"token": "READY", "meaning": "A true status token."},
            },
            "check": "token_match",
            "substitutions": [],
            "c7_expect": _c7_expect("facet_token", "TOOTH_TOKEN"),
        },
    ]
    stage2["exports"] = [
        {
            "claim_id": "define_y",
            "source_digest": source_digest,
            "c7_binding": _c7_binding(command, "facet_y"),
        }
    ]
    stage2["knobs"] = [
        {
            "knob_id": "k.row_a",
            "symbol": "x",
            "registry_row": "row_a",
            "action": "inherited",
            "origin": "stage001/define_x",
            "effective_stage": "stage002",
            "count_effect": {"low": 0, "high": 0},
            "count_category": "base",
            "pending": False,
            "evidence": ev("inherit_row_a"),
        },
        {
            "knob_id": "k.row_b",
            "symbol": "ell",
            "registry_row": "row_b",
            "action": "discharged",
            "origin": "stage001/define_ell",
            "effective_stage": "stage002",
            "count_effect": {"low": 0, "high": 0},
            "count_category": "route_ful_debt",
            "pending": False,
            "discharge_evidence": "here/discharge_b",
            "evidence": ev("discharge_row_b"),
        },
    ]

    stage3 = _self_test_manifest(
        source,
        source_digest,
        "stage003",
        [],
        [_self_test_claim(source, source_digest, "stage3_identity")],
    )
    stage3["verification"]["teeth"] = [
        {
            "predicate": "TOOTH_STAGE3",
            "mutation": "stage3",
            "claim_ids": ["stage3_identity"],
            "evidence": ev("TOOTH_STAGE3"),
        }
    ]
    stage3["consumes"] = [
        {
            "ref": "stage001/status_token",
            "as_consumed": {
                "payload_kind": "token",
                "payload": {"token": "READY", "meaning": "A true status token."},
            },
            "check": "token_match",
            "substitutions": [],
            "c7_expect": _c7_expect("facet_token", "TOOTH_STAGE3"),
        }
    ]

    def replace_source_path(value: Any, old: str, new: str) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                if child == old:
                    value[key] = new
                else:
                    replace_source_path(child, old, new)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                if child == old:
                    value[index] = new
                else:
                    replace_source_path(child, old, new)

    replace_source_path(stage2, str(source), str(stage2_source))
    replace_source_path(stage3, str(source), str(stage3_source))

    tracked_anchor_files = [
        source_file,
        stage2_source_file,
        stage3_source_file,
        bare_source_file,
        unregistered_bare_source_file,
        duplicate_bare_source_file,
        conditional_bare_source_file,
        mutated_bare_source_file,
        augmented_bare_source_file,
        four_axis_source_file,
        four_axis_peer_source_file,
        fractional_basis_source_file,
        other_stage_dim_source_file,
        non_python_dim_source_file,
        midway_dim_source_file,
    ]
    git_env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "Composite Self-Test",
        "GIT_AUTHOR_EMAIL": "composite-self-test@example.invalid",
        "GIT_COMMITTER_NAME": "Composite Self-Test",
        "GIT_COMMITTER_EMAIL": "composite-self-test@example.invalid",
    }
    subprocess.run(
        ("git", "init", "-q", str(root)),
        check=True,
        capture_output=True,
        env=git_env,
    )
    subprocess.run(
        (
            "git",
            "-C",
            str(root),
            "add",
            "--",
            *(str(path.relative_to(root)) for path in tracked_anchor_files),
        ),
        check=True,
        capture_output=True,
        env=git_env,
    )
    subprocess.run(
        ("git", "-C", str(root), "commit", "-q", "-m", "self-test anchors"),
        check=True,
        capture_output=True,
        env=git_env,
    )

    config = {
        "schema_version": "2.1",
        "parameter_register": {
            "path": str(register),
            "sha256": register_digest,
            "format": "plain",
        },
        "range_claim_ref": "stage001/count_range",
        "closed_slice": True,
    }
    paths = {
        "source": str(stage3_source),
        "stage1_source": str(source),
        "stage2_source": str(stage2_source),
        "stage3_source": str(stage3_source),
        "source_digest": source_digest,
        "bare_source": str(bare_source),
        "bare_source_digest": bare_source_digest,
        "unregistered_bare_source": str(unregistered_bare_source),
        "unregistered_bare_source_digest": unregistered_bare_source_digest,
        "mathematica_output": str(mathematica_output),
        "mathematica_output_absolute": str(mathematica_output_file),
        "mathematica_output_digest": mathematica_output_digest,
        "mathematica_script": str(mathematica_script),
        "empty_mathematica_script": str(empty_mathematica_script),
        "empty_mathematica_output": str(empty_mathematica_output),
        "empty_mathematica_output_digest": empty_mathematica_output_digest,
        "other_mathematica_output": str(other_mathematica_output),
        "other_mathematica_output_digest": other_mathematica_output_digest,
        "missing_mathematica_script": str(missing_mathematica_script),
        "duplicate_bare_source": str(duplicate_bare_source),
        "duplicate_bare_source_digest": duplicate_bare_source_digest,
        "conditional_bare_source": str(conditional_bare_source),
        "conditional_bare_source_digest": conditional_bare_source_digest,
        "mutated_bare_source": str(mutated_bare_source),
        "mutated_bare_source_digest": mutated_bare_source_digest,
        "augmented_bare_source": str(augmented_bare_source),
        "augmented_bare_source_digest": augmented_bare_source_digest,
        "four_axis_source": str(four_axis_source),
        "four_axis_source_digest": four_axis_source_digest,
        "four_axis_peer_source": str(four_axis_peer_source),
        "four_axis_peer_source_digest": four_axis_peer_source_digest,
        "fractional_basis_source": str(fractional_basis_source),
        "fractional_basis_source_digest": fractional_basis_source_digest,
        "attacker_dim_source": str(attacker_dim_source),
        "attacker_dim_source_absolute": str(attacker_dim_source_file),
        "attacker_dim_source_digest": attacker_dim_source_digest,
        "untracked_symlink": str(untracked_symlink),
        "untracked_symlink_digest": source_digest,
        "other_stage_dim_source": str(other_stage_dim_source),
        "other_stage_dim_source_digest": other_stage_dim_source_digest,
        "non_python_dim_source": str(non_python_dim_source),
        "non_python_dim_source_digest": non_python_dim_source_digest,
        "midway_dim_source": str(midway_dim_source),
        "midway_dim_source_digest": source_digest,
        "register": str(register),
        "session": str(session),
        "command": command,
    }
    return [stage1, stage2, stage3], config, paths


def _find_claim(manifests: list[dict[str, Any]], stage: str, cid: str) -> dict[str, Any]:
    manifest = next(item for item in manifests if item["stage_id"] == stage)
    return next(item for item in manifest["claims"] if item["id"] == cid)


def _find_consume(manifests: list[dict[str, Any]], stage: str, ref: str) -> dict[str, Any]:
    manifest = next(item for item in manifests if item["stage_id"] == stage)
    return next(item for item in manifest["consumes"] if item["ref"] == ref)


def _find_export(manifests: list[dict[str, Any]], stage: str, cid: str) -> dict[str, Any]:
    manifest = next(item for item in manifests if item["stage_id"] == stage)
    return next(item for item in manifest["exports"] if item["claim_id"] == cid)


def run_self_test() -> int:
    schema = load_json(SCHEMA_PATH)
    lines = ["SELF-TEST stage-manifest composite v2.1"]
    failures: list[str] = []
    passed = 0
    total = 0

    with tempfile.TemporaryDirectory(prefix="manifest_v21_selftest_") as tmp:
        root = Path(tmp)
        baseline_manifests, baseline_config, paths = _build_self_test_fixture(root)
        self_test_bare_orders = {
            paths["bare_source_digest"]: "LTM",
            paths["duplicate_bare_source_digest"]: "LTM",
            paths["conditional_bare_source_digest"]: "LTM",
            paths["mutated_bare_source_digest"]: "LTM",
            paths["augmented_bare_source_digest"]: "LTM",
            paths["four_axis_source_digest"]: (
                "M",
                "L",
                "T",
                "E-charge",
            ),
            paths["four_axis_peer_source_digest"]: (
                "M",
                "L",
                "T",
                "E-charge",
            ),
            paths["fractional_basis_source_digest"]: (
                "stiffness",
                "length",
                "time",
            ),
        }

        def execute(
            manifests: list[dict[str, Any]],
            config: dict[str, Any],
            *,
            closed_slice: bool | None = True,
            inject_self_test_bare_orders: bool = True,
        ) -> CompositeReport:
            return CompositeChecker(
                manifests,
                config,
                schema=schema,
                roots=(root, PROJECT_ROOT, PROJECT_ROOT.parent.parent),
                closed_slice=closed_slice,
                bare_tuple_dim_orders=(
                    self_test_bare_orders
                    if inject_self_test_bare_orders
                    else {}
                ),
            ).run()

        clean = execute(copy.deepcopy(baseline_manifests), copy.deepcopy(baseline_config))
        clean_statuses = {name: clean.results[name].status for name in CHECK_NAMES}
        total += 1
        if clean.headline == "PASS" and all(status == "PASS" for status in clean_statuses.values()):
            passed += 1
            lines.append("[PASS] clean fixture | target=ALL | observed=ALL PASS")
        else:
            failures.append(f"clean fixture: {clean_statuses}")
            lines.append(f"[FAIL] clean fixture | {clean_statuses}")

        total += 1
        production_manifests = [
            load_json(path)
            for path in sorted(DEFAULT_MANIFEST_DIR.glob("*.json"))
        ]
        production_report = CompositeChecker(
            production_manifests,
            load_json(DEFAULT_CONFIG),
            roots=(PROJECT_ROOT, PROJECT_ROOT.parent.parent, Path.cwd()),
        ).run()
        production_dimensional = production_report.results["DIMENSIONAL_CONSISTENCY"]
        production_non_dimensional = {
            name: result.status
            for name, result in production_report.results.items()
            if name != "DIMENSIONAL_CONSISTENCY"
        }
        expected_production_non_dimensional = {
            name: (
                "PARTIAL" if name in {"C3", "C5", "C7"} else "PASS"
            )
            for name in CHECK_NAMES
            if name != "DIMENSIONAL_CONSISTENCY"
        }
        if (
            production_dimensional.status == "FAIL"
            and len(production_dimensional.issues) == 1
            and production_dimensional.issues[0].code == "DIM_SOURCE_NOT_REGISTERED"
            and production_dimensional.issues[0].detail.startswith("stage032:")
            and production_non_dimensional == expected_production_non_dimensional
        ):
            passed += 1
            lines.append(
                "[PASS] DIMENSIONAL_CONSISTENCY production structural anchors isolate stage032 borrow | "
                "target=DIMENSIONAL_CONSISTENCY | clean=same-stage anchors PASS planted=FAIL | "
                "code=DIM_SOURCE_NOT_REGISTERED"
            )
        else:
            failures.append(
                "DIMENSIONAL_CONSISTENCY production structural dimension anchors: "
                f"status={production_dimensional.status}, "
                f"issues={[(issue.code, issue.detail) for issue in production_dimensional.issues]}, "
                f"extras={production_non_dimensional}"
            )
            lines.append(
                "[FAIL] DIMENSIONAL_CONSISTENCY production structural anchors isolate stage032 borrow | "
                f"observed={production_dimensional.status}"
            )

        class FlippingReadSource:
            def __init__(self) -> None:
                self.calls = 0
                self.trusted = b"dim_trusted = (1, 0, 0)\n"
                self.attacker = b"dim_attacker = (9, 9, 9)\n"

            def read_bytes(self) -> bytes:
                self.calls += 1
                return self.trusted if self.calls == 1 else self.attacker

            def __str__(self) -> str:
                return "flipping_dim_source.py"

        total += 1
        flipping_source = FlippingReadSource()
        snapshot = read_dim_source_snapshot(flipping_source)  # type: ignore[arg-type]
        snapshot_names = {
            node.id
            for node in ast.walk(snapshot.tree)
            if isinstance(node, ast.Name)
        }
        if (
            flipping_source.calls == 1
            and snapshot.sha256
            == hashlib.sha256(flipping_source.trusted).hexdigest()
            and "dim_trusted" in snapshot_names
            and "dim_attacker" not in snapshot_names
        ):
            passed += 1
            lines.append(
                "[PASS] DIMENSIONAL_CONSISTENCY dimension snapshot hashes parsed bytes | "
                "clean=PASS planted=REJECTED | code=DIM_SOURCE_SNAPSHOT_COHERENT"
            )
        else:
            failures.append(
                "DIMENSIONAL_CONSISTENCY dimension snapshot hashes parsed bytes: "
                f"calls={flipping_source.calls}, sha256={snapshot.sha256}, "
                f"names={sorted(snapshot_names)}"
            )
            lines.append(
                "[FAIL] DIMENSIONAL_CONSISTENCY dimension snapshot hashes parsed bytes | "
                "code=DIM_SOURCE_SNAPSHOT_COHERENT"
            )

        def case(
            label: str,
            target: str,
            planted_status: str,
            code: str,
            mutate: Any,
            *,
            exact_target_code: bool = False,
            inject_self_test_bare_orders: bool = True,
        ) -> None:
            nonlocal passed, total
            total += 1
            manifests = copy.deepcopy(baseline_manifests)
            config = copy.deepcopy(baseline_config)
            mutate(manifests, config)
            report = execute(
                manifests,
                config,
                inject_self_test_bare_orders=inject_self_test_bare_orders,
            )
            statuses = {name: report.results[name].status for name in CHECK_NAMES}
            issue_codes = {issue.code for issue in report.results[target].issues}
            expected_non_target = {name: "PASS" for name in CHECK_NAMES if name != target}
            if target == "SCHEMA" and planted_status == "FAIL":
                expected_non_target = {name: "SKIPPED" for name in CHECK_NAMES if name != target}
            elif target == "C3" and planted_status in BAD_STATUSES:
                expected_non_target["C7"] = "SKIPPED"
            non_target = {
                name: value
                for name, value in statuses.items()
                if name != target and value != expected_non_target[name]
            }
            target_issues = report.results[target].issues
            code_ok = code in issue_codes
            if exact_target_code:
                code_ok = (
                    len(target_issues) == 1
                    and target_issues[0].code == code
                )
            ok = statuses[target] == planted_status and code_ok and not non_target
            if ok:
                passed += 1
                lines.append(
                    f"[PASS] {label} | target={target} | clean=PASS planted={planted_status} | code={code}"
                )
            else:
                failures.append(
                    f"{label}: target={target}/{planted_status}/{code}, statuses={statuses}, "
                    f"codes={sorted(issue_codes)}, non_target={non_target}"
                )
                lines.append(
                    f"[FAIL] {label} | target={target} expected={planted_status}/{code} "
                    f"observed={statuses[target]}/{sorted(issue_codes)} extras={non_target}"
                )

        def clean_case(
            label: str,
            mutate: Any,
            verify: Any = None,
        ) -> None:
            nonlocal passed, total
            total += 1
            manifests = copy.deepcopy(baseline_manifests)
            config = copy.deepcopy(baseline_config)
            mutate(manifests, config)
            report = execute(manifests, config)
            statuses = {name: report.results[name].status for name in CHECK_NAMES}
            verified = True if verify is None else bool(verify(report))
            if report.headline == "PASS" and all(
                status == "PASS" for status in statuses.values()
            ) and verified:
                passed += 1
                lines.append(
                    f"[PASS] {label} | target=ALL | observed=ALL PASS"
                )
            else:
                failures.append(
                    f"{label}: headline={report.headline}, statuses={statuses}, "
                    f"verified={verified}, "
                    f"issues={ {name: [(issue.code, issue.detail) for issue in result.issues] for name, result in report.results.items() if result.issues} }"
                )
                lines.append(
                    f"[FAIL] {label} | headline={report.headline} statuses={statuses}"
                )

        def paired_case(
            label: str,
            target: str,
            planted_status: str,
            code: str,
            arrange_clean: Any,
            plant: Any,
            verify_clean: Any = None,
            verify_planted: Any = None,
        ) -> None:
            nonlocal passed, total
            total += 1
            clean_manifests = copy.deepcopy(baseline_manifests)
            clean_config = copy.deepcopy(baseline_config)
            arrange_clean(clean_manifests, clean_config)
            clean_report = execute(clean_manifests, clean_config)
            clean_statuses = {
                name: clean_report.results[name].status for name in CHECK_NAMES
            }
            clean_verified = (
                True
                if verify_clean is None
                else bool(verify_clean(clean_report))
            )

            planted_manifests = copy.deepcopy(clean_manifests)
            planted_config = copy.deepcopy(clean_config)
            plant(planted_manifests, planted_config)
            planted_report = execute(planted_manifests, planted_config)
            planted_statuses = {
                name: planted_report.results[name].status
                for name in CHECK_NAMES
            }
            planted_codes = {
                issue.code
                for issue in planted_report.results[target].issues
            }
            planted_verified = (
                True
                if verify_planted is None
                else bool(verify_planted(planted_report))
            )
            expected_non_target = {
                name: "PASS" for name in CHECK_NAMES if name != target
            }
            if target == "SCHEMA" and planted_status == "FAIL":
                expected_non_target = {
                    name: "SKIPPED"
                    for name in CHECK_NAMES
                    if name != target
                }
            planted_extras = {
                name: status
                for name, status in planted_statuses.items()
                if name != target
                and status != expected_non_target[name]
            }
            clean_ok = (
                clean_report.headline == "PASS"
                and all(status == "PASS" for status in clean_statuses.values())
                and clean_verified
            )
            planted_ok = (
                planted_statuses[target] == planted_status
                and code in planted_codes
                and not planted_extras
                and planted_verified
            )
            if clean_ok and planted_ok:
                passed += 1
                lines.append(
                    f"[PASS] {label} | target={target} | clean=PASS "
                    f"planted={planted_status} | code={code}"
                )
            else:
                failures.append(
                    f"{label}: clean={clean_statuses}/verified={clean_verified}, "
                    f"clean_issues="
                    f"{ {name: [(issue.code, issue.detail) for issue in result.issues] for name, result in clean_report.results.items() if result.issues} }, "
                    f"planted={planted_statuses}, codes={sorted(planted_codes)}, "
                    f"extras={planted_extras}, verified={planted_verified}"
                )
                lines.append(
                    f"[FAIL] {label} | clean={clean_statuses} "
                    f"planted={planted_statuses[target]}/"
                    f"{sorted(planted_codes)} extras={planted_extras}"
                )

        total += 1
        schema_invalid_manifests = copy.deepcopy(baseline_manifests)
        schema_invalid_manifests[0]["schema_version"] = "1.0"
        schema_invalid = execute(schema_invalid_manifests, copy.deepcopy(baseline_config))
        schema_statuses = {name: schema_invalid.results[name].status for name in CHECK_NAMES}
        downstream_statuses = {name: schema_statuses[name] for name in CHECK_NAMES if name != "SCHEMA"}
        schema_codes = {issue.code for issue in schema_invalid.results["SCHEMA"].issues}
        schema_skip_ok = (
            clean.headline == "PASS"
            and schema_invalid.headline == "FAIL"
            and schema_statuses["SCHEMA"] == "FAIL"
            and "SCHEMA_VALIDATION" in schema_codes
            and all(status == "SKIPPED" for status in downstream_statuses.values())
            and "PASS" not in downstream_statuses.values()
        )
        if schema_skip_ok:
            passed += 1
            lines.append(
                "[PASS] schema-invalid short circuit | target=SCHEMA:FAIL | "
                "downstream=ALL SKIPPED"
            )
        else:
            failures.append(
                "schema-invalid short circuit: "
                f"clean={clean.headline}, planted={schema_invalid.headline}, statuses={schema_statuses}, "
                f"codes={sorted(schema_codes)}"
            )
            lines.append(
                "[FAIL] schema-invalid short circuit | expected=SCHEMA:FAIL/downstream:ALL SKIPPED "
                f"| observed={schema_statuses}"
            )

        case(
            "UNDECLARED_IMPORT",
            "IMPORT",
            "FAIL",
            "UNDECLARED_IMPORT",
            lambda m, c: next(x for x in m if x["stage_id"] == "stage002").__setitem__(
                "consumes",
                [item for item in next(x for x in m if x["stage_id"] == "stage002")["consumes"]
                 if item["ref"] != "stage001/define_ell"],
            ),
        )

        def false_local(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            stage2 = next(x for x in m if x["stage_id"] == "stage002")
            symbol = next(x for x in stage2["symbols"] if x["quantity_id"] == "q.core.x")
            symbol["definition_ref"] = "here/fake_x"
            stage2["claims"].append(
                _self_test_claim(
                    Path(paths["source"]),
                    paths["source_digest"],
                    "fake_x",
                    payload=_relation_payload("x", "x", "defines"),
                )
            )

        case(
            "FALSE_LOCAL_DEFINITION",
            "IMPORT",
            "FAIL",
            "FALSE_LOCAL_DEFINITION",
            false_local,
        )

        case(
            "definition_ref does not bind quantity",
            "IMPORT",
            "FAIL",
            "DEFINITION_NOT_BINDING",
            lambda m, c: next(
                s
                for s in next(x for x in m if x["stage_id"] == "stage002")["symbols"]
                if s["quantity_id"] == "q.core.x"
            ).__setitem__("definition_ref", "stage001/define_ab"),
        )

        case(
            "DIMENSIONAL_CONSISTENCY order-lie LTM source declared LMT",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_ORDER_MISMATCH",
            lambda m, c: next(
                s
                for s in next(x for x in m if x["stage_id"] == "stage001")["symbols"]
                if s["parse_alias"] == "a_B"
            ).__setitem__("dim_source_order", "LMT"),
        )

        case(
            "DIMENSIONAL_CONSISTENCY function-scoped composite live Dim tuple",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_TUPLE_MISMATCH",
            lambda m, c: next(
                s
                for s in next(x for x in m if x["stage_id"] == "stage001")["symbols"]
                if s["parse_alias"] == "K_m"
            ).__setitem__("dim_source_tuple", ["4", "1", "-2"]),
        )

        def f15_unbacked_anchor(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            symbol = next(
                item
                for item in next(
                    manifest
                    for manifest in manifests
                    if manifest["stage_id"] == "stage001"
                )["symbols"]
                if item["parse_alias"] == "a_B"
            )
            symbol["dim_source"]["script_path"] = paths["attacker_dim_source"]
            symbol["dim_source_order"] = "LMT"
            symbol["dim_source_tuple"] = ["-2", "1", "-2"]

        case(
            "DIMENSIONAL_CONSISTENCY attacker-authored coherent anchor lacks stage evidence",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NOT_STAGE_EVIDENCE",
            f15_unbacked_anchor,
            exact_target_code=True,
        )

        def self_cited_unregistered_anchor(
            manifests: list[dict[str, Any]],
            _: dict[str, Any],
            source_key: str,
            digest_key: str,
            order: str,
            raw: Sequence[str],
        ) -> None:
            symbol = next(
                item
                for item in next(
                    manifest
                    for manifest in manifests
                    if manifest["stage_id"] == "stage001"
                )["symbols"]
                if item["parse_alias"] == "a_B"
            )
            source = Path(paths[source_key])
            symbol["dim_source"] = {
                "script_path": str(source),
                "locus": "a_B",
            }
            symbol["dim_source_order"] = order
            symbol["dim_source_tuple"] = list(raw)
            symbol["evidence"] = _self_test_evidence(
                source,
                paths[digest_key],
                "a_B",
            )

        case(
            "DIMENSIONAL_CONSISTENCY self-cited flipped-order anchor is not admissible",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NOT_REGISTERED",
            lambda m, c: self_cited_unregistered_anchor(
                m,
                c,
                "attacker_dim_source",
                "attacker_dim_source_digest",
                "LMT",
                ("-2", "1", "-2"),
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY self-cited correct-order copied anchor is not admissible",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NOT_REGISTERED",
            lambda m, c: self_cited_unregistered_anchor(
                m,
                c,
                "other_stage_dim_source",
                "other_stage_dim_source_digest",
                "LTM",
                ("-2", "-2", "1"),
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY untracked conventional symlink is not an admissible anchor",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NOT_REGISTERED",
            lambda m, c: self_cited_unregistered_anchor(
                m,
                c,
                "untracked_symlink",
                "untracked_symlink_digest",
                "LTM",
                ("-2", "-2", "1"),
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY cross-cutting midway audit is not a stage anchor",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NOT_REGISTERED",
            lambda m, c: self_cited_unregistered_anchor(
                m,
                c,
                "midway_dim_source",
                "midway_dim_source_digest",
                "LTM",
                ("-2", "-2", "1"),
            ),
            exact_target_code=True,
        )

        clean_case(
            "DIMENSIONAL_CONSISTENCY normalized same-stage anchor spelling certifies",
            lambda m, c: next(
                symbol
                for symbol in next(
                    manifest
                    for manifest in m
                    if manifest["stage_id"] == "stage001"
                )["symbols"]
                if symbol["parse_alias"] == "a_B"
            )["dim_source"].__setitem__(
                "script_path",
                "./scripts//ledger_stage001_self_test_sympy_audit.py",
            ),
        )

        def relabel_stage_id_only(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            stage3["consumes"] = []
            stage3["verification"]["teeth"] = []
            source_path = Path(paths["stage3_source"])
            stage3["symbols"].append(
                _self_test_symbol(
                    source_path,
                    paths["source_digest"],
                    name="z",
                    alias="z",
                    qid="q.stage3.stage_id_binding",
                    definition="here/define_stage_id_z",
                    dim={},
                    raw=("0", "0", "0"),
                )
            )
            stage3["claims"].append(
                _self_test_claim(
                    source_path,
                    paths["source_digest"],
                    "define_stage_id_z",
                    payload=_relation_payload("z", "z", "defines"),
                )
            )
            stage3["stage_id"] = "stage004"

        case(
            "DIMENSIONAL_CONSISTENCY manifest stage_id relabel cannot carry its old anchor",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NOT_REGISTERED",
            relabel_stage_id_only,
            exact_target_code=True,
        )

        total += 1
        non_python_manifests = copy.deepcopy(baseline_manifests)
        non_python_config = copy.deepcopy(baseline_config)
        self_cited_unregistered_anchor(
            non_python_manifests,
            non_python_config,
            "non_python_dim_source",
            "non_python_dim_source_digest",
            "LTM",
            ("-2", "-2", "1"),
        )
        non_python_report_path = root / "non_python_dim_source_report.md"
        try:
            non_python_report = execute(
                non_python_manifests,
                non_python_config,
            )
        except Exception as exc:
            failures.append(
                "DIMENSIONAL_CONSISTENCY cited registered non-Python anchor: "
                f"uncaught={type(exc).__name__}: {exc}"
            )
            lines.append(
                "[FAIL] DIMENSIONAL_CONSISTENCY cited registered non-Python anchor reports cleanly | "
                f"code=UNCAUGHT_{type(exc).__name__}"
            )
        else:
            non_python_report_path.write_text(render_report(non_python_report))
            non_python_issues = non_python_report.results["DIMENSIONAL_CONSISTENCY"].issues
            non_python_other_statuses = {
                name: result.status
                for name, result in non_python_report.results.items()
                if name != "DIMENSIONAL_CONSISTENCY" and result.status != "PASS"
            }
            if (
                non_python_report.results["DIMENSIONAL_CONSISTENCY"].status == "UNSUPPORTED"
                and len(non_python_issues) == 1
                and non_python_issues[0].code == "DIM_SOURCE_UNSUPPORTED"
                and non_python_report_path.is_file()
                and "DIM_SOURCE_UNSUPPORTED"
                in non_python_report_path.read_text()
                and not non_python_other_statuses
            ):
                passed += 1
                lines.append(
                    "[PASS] DIMENSIONAL_CONSISTENCY cited registered non-Python anchor reports cleanly | "
                    "target=DIMENSIONAL_CONSISTENCY | clean=PASS planted=UNSUPPORTED | "
                    "code=DIM_SOURCE_UNSUPPORTED | report=written"
                )
            else:
                failures.append(
                    "DIMENSIONAL_CONSISTENCY cited registered non-Python anchor: "
                    f"status={non_python_report.results['DIMENSIONAL_CONSISTENCY'].status}, "
                    f"issues={[(issue.code, issue.detail) for issue in non_python_issues]}, "
                    f"extras={non_python_other_statuses}, "
                    f"report_exists={non_python_report_path.is_file()}"
                )
                lines.append(
                    "[FAIL] DIMENSIONAL_CONSISTENCY cited registered non-Python anchor reports cleanly"
                )

        def cross_stage_unbacked_anchor(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage1 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )
            stage2 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage002"
            )
            stage1["symbols"][0]["evidence"] = _self_test_evidence(
                Path(paths["other_stage_dim_source"]),
                paths["other_stage_dim_source_digest"],
                "x",
            )
            stage2["symbols"][0]["dim_source"]["script_path"] = paths[
                "other_stage_dim_source"
            ]

        case(
            "DIMENSIONAL_CONSISTENCY another stage anchor without local evidence backing",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NOT_STAGE_EVIDENCE",
            cross_stage_unbacked_anchor,
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY absolute dimension anchor path rejected",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "UNSAFE_DIM_SOURCE_PATH",
            lambda m, c: next(
                item
                for item in next(
                    manifest
                    for manifest in m
                    if manifest["stage_id"] == "stage001"
                )["symbols"]
                if item["parse_alias"] == "a_B"
            )["dim_source"].__setitem__(
                "script_path", paths["attacker_dim_source_absolute"]
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY parent-escape dimension anchor path rejected",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "UNSAFE_DIM_SOURCE_PATH",
            lambda m, c: next(
                item
                for item in next(
                    manifest
                    for manifest in m
                    if manifest["stage_id"] == "stage001"
                )["symbols"]
                if item["parse_alias"] == "a_B"
            )["dim_source"].__setitem__(
                "script_path", "../attacker_dim_source.py"
            ),
            exact_target_code=True,
        )

        def pi_area(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            rhs = "4*pi*ell**2"
            _find_claim(manifests, "stage001", "define_area")["payload"]["rhs"][
                "sympy"
            ] = rhs
            _find_consume(
                manifests, "stage002", "stage001/define_area"
            )["as_consumed"]["payload"]["rhs"]["sympy"] = rhs

        clean_case(
            "DIMENSIONAL_CONSISTENCY NumberSymbol pi dimensional area",
            pi_area,
        )

        def pi_inhomogeneity(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            pi_area(manifests, config)
            rhs = "4*pi*ell + ell**2"
            _find_claim(manifests, "stage001", "define_area")["payload"]["rhs"][
                "sympy"
            ] = rhs
            _find_consume(
                manifests, "stage002", "stage001/define_area"
            )["as_consumed"]["payload"]["rhs"]["sympy"] = rhs

        case(
            "DIMENSIONAL_CONSISTENCY NumberSymbol pi does not mask dimensional break",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIMENSIONAL_INHOMOGENEITY",
            pi_inhomogeneity,
            exact_target_code=True,
        )

        def add_bare_tuple_symbol(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            bare_source = Path(paths["bare_source"])
            bare_digest = paths["bare_source_digest"]
            symbol = _self_test_symbol(
                bare_source,
                bare_digest,
                name="bare_A",
                alias="bare_A",
                qid="q.bare.a",
                definition="here/define_bare_a",
                dim={"L": "3", "M": "1", "T": "-2"},
                raw=("3", "-2", "1"),
            )
            symbol["dim_source"]["locus"] = (
                "run_units: dim_A = tuple(x + y for x, y in zip(dim_E, dim_L))"
            )
            stage3["symbols"].append(symbol)
            stage3["claims"].append(
                _self_test_claim(
                    bare_source,
                    bare_digest,
                    "define_bare_a",
                    payload=_relation_payload(
                        "bare_A", "bare_A", "defines"
                    ),
                )
            )

        clean_case(
            "DIMENSIONAL_CONSISTENCY function-local bare-tuple definition certifies",
            add_bare_tuple_symbol,
        )

        def use_bare_variant(
            manifests: list[dict[str, Any]],
            config: dict[str, Any],
            source_key: str,
            digest_key: str,
            locus: str,
        ) -> None:
            add_bare_tuple_symbol(manifests, config)
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            symbol = next(
                item
                for item in stage3["symbols"]
                if item["parse_alias"] == "bare_A"
            )
            source = Path(paths[source_key])
            digest = paths[digest_key]
            symbol["dim_source"]["script_path"] = str(source)
            symbol["dim_source"]["locus"] = locus
            symbol["evidence"] = _self_test_evidence(source, digest, "bare_A")
            claim = next(
                item
                for item in stage3["claims"]
                if item["id"] == "define_bare_a"
            )
            claim["evidence"] = _self_test_evidence(
                source, digest, "define_bare_a"
            )

        case(
            "DIMENSIONAL_CONSISTENCY shadowed duplicate function body fails closed",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            lambda m, c: use_bare_variant(
                m,
                c,
                "duplicate_bare_source",
                "duplicate_bare_source_digest",
                "run_units: dim_A = tuple(x + y for x, y in zip(dim_E, dim_L))",
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY conditional-only bare dimension fails closed",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            lambda m, c: use_bare_variant(
                m,
                c,
                "conditional_bare_source",
                "conditional_bare_source_digest",
                "run_units: dim_A = (3, -2, 1)",
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY post-assignment list dimension mutation fails closed",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            lambda m, c: use_bare_variant(
                m,
                c,
                "mutated_bare_source",
                "mutated_bare_source_digest",
                "run_units: dim_A = [3, -2, 1]",
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY augmented bare dimension mutation fails closed",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            lambda m, c: use_bare_variant(
                m,
                c,
                "augmented_bare_source",
                "augmented_bare_source_digest",
                "run_units: dim_A = [3, -2, 1]",
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY production registry rejects self-test bare digest",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            add_bare_tuple_symbol,
            exact_target_code=True,
            inject_self_test_bare_orders=False,
        )

        def bare_tuple_mismatch(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            add_bare_tuple_symbol(manifests, config)
            symbol = next(
                item
                for item in next(
                    manifest
                    for manifest in manifests
                    if manifest["stage_id"] == "stage003"
                )["symbols"]
                if item["parse_alias"] == "bare_A"
            )
            symbol["dim_source_tuple"] = ["4", "-2", "1"]
            symbol["dim"] = {"L": "4", "M": "1", "T": "-2"}

        case(
            "DIMENSIONAL_CONSISTENCY function-local bare-tuple declared tuple drift",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_TUPLE_MISMATCH",
            bare_tuple_mismatch,
            exact_target_code=True,
        )

        def bare_tuple_order_lie(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            add_bare_tuple_symbol(manifests, config)
            symbol = next(
                item
                for item in next(
                    manifest
                    for manifest in manifests
                    if manifest["stage_id"] == "stage003"
                )["symbols"]
                if item["parse_alias"] == "bare_A"
            )
            symbol["dim_source_order"] = "LMT"

        case(
            "DIMENSIONAL_CONSISTENCY function-local bare-tuple manifest order lie",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_ORDER_MISMATCH",
            bare_tuple_order_lie,
            exact_target_code=True,
        )

        def unregistered_bare_tuple_order(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            add_bare_tuple_symbol(manifests, config)
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            symbol = next(
                item
                for item in stage3["symbols"]
                if item["parse_alias"] == "bare_A"
            )
            source = Path(paths["unregistered_bare_source"])
            digest = paths["unregistered_bare_source_digest"]
            symbol["dim_source"]["script_path"] = str(source)
            symbol["evidence"] = _self_test_evidence(source, digest, "bare_A")
            claim = next(
                item
                for item in stage3["claims"]
                if item["id"] == "define_bare_a"
            )
            claim["evidence"] = _self_test_evidence(
                source, digest, "define_bare_a"
            )

        case(
            "DIMENSIONAL_CONSISTENCY unregistered bare-tuple order fails closed",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            unregistered_bare_tuple_order,
            exact_target_code=True,
        )

        def unsupported_dim_source(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage3 = next(x for x in manifests if x["stage_id"] == "stage003")
            stage3["symbols"].append(
                _self_test_symbol(
                    Path(paths["source"]),
                    paths["source_digest"],
                    name="bad",
                    alias="bad",
                    qid="q.unsupported.bad",
                    definition="here/define_bad",
                    dim={"L": "1"},
                    raw=("1", "0", "0"),
                )
            )
            stage3["claims"].append(
                _self_test_claim(
                    Path(paths["source"]),
                    paths["source_digest"],
                    "define_bad",
                    payload=_relation_payload("bad", "bad", "defines"),
                )
            )

        case(
            "DIMENSIONAL_CONSISTENCY unsupported live Dim operation raises",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            unsupported_dim_source,
        )

        def container_label_only(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            symbol = _self_test_symbol(
                Path(paths["source"]),
                paths["source_digest"],
                name="container_K_m",
                alias="container_K_m",
                qid="q.unsupported.container_k_m",
                definition="here/define_container_k_m",
                dim={"L": "1"},
                raw=("1", "0", "0"),
            )
            symbol["dim_source"]["locus"] = (
                "container_only: container label container_K_m"
            )
            stage3["symbols"].append(symbol)
            stage3["claims"].append(
                _self_test_claim(
                    Path(paths["source"]),
                    paths["source_digest"],
                    "define_container_k_m",
                    payload=_relation_payload(
                        "container_K_m", "container_K_m", "defines"
                    ),
                )
            )

        case(
            "DIMENSIONAL_CONSISTENCY container label is not a Dim binding",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            container_label_only,
        )

        def decoy_assignment(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            symbol = _self_test_symbol(
                Path(paths["source"]),
                paths["source_digest"],
                name="K_decoy",
                alias="K_decoy",
                qid="q.unsupported.k_decoy",
                definition="here/define_k_decoy",
                dim={"L": "4", "M": "1", "T": "-2"},
                raw=("4", "-2", "1"),
            )
            symbol["dim_source"]["locus"] = (
                "run_units: K_decoy = ENERGY * LENGTH**2"
            )
            stage3["symbols"].append(symbol)
            stage3["claims"].append(
                _self_test_claim(
                    Path(paths["source"]),
                    paths["source_digest"],
                    "define_k_decoy",
                    payload=_relation_payload(
                        "K_decoy", "K_decoy", "defines"
                    ),
                )
            )

        case(
            "DIMENSIONAL_CONSISTENCY decoy literal cannot impersonate locus assignment",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            decoy_assignment,
        )

        case(
            "DIMENSIONAL_CONSISTENCY failed locus expression cannot fall through to alias",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_EXPRESSION_NOT_BINDING_ROOTED",
            lambda m, c: next(
                s
                for s in next(
                    x for x in m if x["stage_id"] == "stage001"
                )["symbols"]
                if s["parse_alias"] == "K_m"
            )["dim_source"].__setitem__(
                "locus", "dimension expression MISSING"
            ),
        )

        def legitimate_composite_expression(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            symbol = _self_test_symbol(
                Path(paths["source"]),
                paths["source_digest"],
                name="inverse_area",
                alias="inverse_area",
                qid="q.legitimate.inverse_area",
                definition="here/define_inverse_area",
                dim={"L": "-2"},
                raw=("-2", "0", "0"),
            )
            symbol["dim_source"]["locus"] = (
                "dimension expression LENGTH**-2 in run_units"
            )
            stage3["symbols"].append(symbol)
            stage3["claims"].append(
                _self_test_claim(
                    Path(paths["source"]),
                    paths["source_digest"],
                    "define_inverse_area",
                    payload=_relation_payload(
                        "inverse_area", "inverse_area", "defines"
                    ),
                )
            )

        total += 1
        composite_manifests = copy.deepcopy(baseline_manifests)
        composite_config = copy.deepcopy(baseline_config)
        legitimate_composite_expression(composite_manifests, composite_config)
        composite_report = execute(composite_manifests, composite_config)
        if composite_report.headline == "PASS" and all(
            result.status == "PASS"
            for result in composite_report.results.values()
        ):
            passed += 1
            lines.append(
                "[PASS] DIMENSIONAL_CONSISTENCY legitimate composite expression without named assignment | "
                "target=DIMENSIONAL_CONSISTENCY | clean=PASS planted=PASS"
            )
        else:
            failures.append(
                "DIMENSIONAL_CONSISTENCY legitimate composite expression: "
                f"headline={composite_report.headline}, "
                f"statuses={ {name: result.status for name, result in composite_report.results.items()} }, "
                f"issues={[(issue.code, issue.detail) for issue in composite_report.results['DIMENSIONAL_CONSISTENCY'].issues]}"
            )
            lines.append(
                "[FAIL] DIMENSIONAL_CONSISTENCY legitimate composite expression without named assignment"
            )

        def manifest_dim_literal(
            manifests: list[dict[str, Any]],
            _: dict[str, Any],
            *,
            alias: str,
            qid: str,
            expression: str,
            raw: Sequence[str],
            dim: dict[str, str],
        ) -> None:
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            source_path = Path(paths["stage3_source"])
            claim_id = f"define_{alias}"
            symbol = _self_test_symbol(
                source_path,
                paths["source_digest"],
                name=alias,
                alias=alias,
                qid=qid,
                definition=f"here/{claim_id}",
                dim=dim,
                raw=raw,
            )
            symbol["dim_source"]["locus"] = (
                f"dimension expression {expression} in run_units"
            )
            stage3["symbols"].append(symbol)
            stage3["claims"].append(
                _self_test_claim(
                    source_path,
                    paths["source_digest"],
                    claim_id,
                    payload=_relation_payload(alias, alias, "defines"),
                )
            )

        case(
            "DIMENSIONAL_CONSISTENCY psi cannot be certified by manifest Dim literal",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_EXPRESSION_NOT_BINDING_ROOTED",
            lambda m, c: manifest_dim_literal(
                m,
                c,
                alias="psi",
                qid="q.attack.psi",
                expression="Dim(11, -7, 13)",
                raw=("11", "-7", "13"),
                dim={"L": "11", "M": "13", "T": "-7"},
            ),
            exact_target_code=True,
        )

        case(
            "DIMENSIONAL_CONSISTENCY z_star cannot be certified by manifest Dim literal",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_EXPRESSION_NOT_BINDING_ROOTED",
            lambda m, c: manifest_dim_literal(
                m,
                c,
                alias="z_star",
                qid="q.attack.z_star",
                expression="Dim(3, 3, 3)",
                raw=("3", "3", "3"),
                dim={"L": "3", "M": "3", "T": "3"},
            ),
            exact_target_code=True,
        )

        def backstop_free_wrong_dim(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage1 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )
            stage2 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage002"
            )
            for manifest in (stage1, stage2):
                f0 = next(
                    symbol
                    for symbol in manifest["symbols"]
                    if symbol["parse_alias"] == "f0"
                )
                f0["dim"] = {"L": "7"}

        case(
            "DIMENSIONAL_CONSISTENCY relation-free symbol dimension has a source backstop",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NAMED_MISMATCH",
            backstop_free_wrong_dim,
        )

        case(
            "DIMENSIONAL_CONSISTENCY named dimension mapping drift",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIM_SOURCE_NAMED_MISMATCH",
            lambda m, c: next(
                symbol
                for symbol in next(
                    manifest
                    for manifest in m
                    if manifest["stage_id"] == "stage001"
                )["symbols"]
                if symbol["parse_alias"] == "a_B"
            ).__setitem__("dim", {"L": "99", "M": "1", "T": "-2"}),
            exact_target_code=True,
        )

        def snapshot_exception_fixture(
            label: str, exception_type: type[Exception]
        ) -> None:
            nonlocal passed, total
            total += 1
            original_reader = globals()["read_dim_source_snapshot"]
            target_path = (root / paths["stage1_source"]).resolve()

            def planted_reader(path: Path) -> DimSourceSnapshot:
                if path.resolve() == target_path:
                    raise exception_type(f"planted {exception_type.__name__}")
                return original_reader(path)

            try:
                globals()["read_dim_source_snapshot"] = planted_reader
                report = execute(
                    copy.deepcopy(baseline_manifests),
                    copy.deepcopy(baseline_config),
                )
            finally:
                globals()["read_dim_source_snapshot"] = original_reader
            codes = {
                issue.code for issue in report.results["DIMENSIONAL_CONSISTENCY"].issues
            }
            non_target = {
                name: result.status
                for name, result in report.results.items()
                if name != "DIMENSIONAL_CONSISTENCY" and result.status != "PASS"
            }
            if (
                clean.headline == "PASS"
                and report.results["DIMENSIONAL_CONSISTENCY"].status == "UNSUPPORTED"
                and codes == {"DIM_SOURCE_UNSUPPORTED"}
                and not non_target
            ):
                passed += 1
                lines.append(
                    f"[PASS] {label} | target=DIMENSIONAL_CONSISTENCY | clean=PASS "
                    "planted=UNSUPPORTED | code=DIM_SOURCE_UNSUPPORTED"
                )
            else:
                failures.append(
                    f"{label}: DIMENSIONAL_CONSISTENCY={report.results['DIMENSIONAL_CONSISTENCY'].status}, "
                    f"codes={sorted(codes)}, extras={non_target}"
                )
                lines.append(f"[FAIL] {label}")

        snapshot_exception_fixture(
            "DIMENSIONAL_CONSISTENCY dimension source OSError is contained", OSError
        )
        snapshot_exception_fixture(
            "DIMENSIONAL_CONSISTENCY dimension source ValueError is contained", ValueError
        )
        snapshot_exception_fixture(
            "DIMENSIONAL_CONSISTENCY dimension source UnicodeError is contained", UnicodeError
        )

        def mismatch_anchor_evidence_digests(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            add_bare_tuple_symbol(manifests, config)
            stage3 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage003"
            )
            next(
                symbol
                for symbol in stage3["symbols"]
                if symbol["parse_alias"] == "bare_A"
            )["evidence"]["source_digest"] = "0" * 64
            next(
                claim
                for claim in stage3["claims"]
                if claim["id"] == "define_bare_a"
            )["evidence"]["source_digest"] = "0" * 64

        total += 1
        digest_mismatch_manifests = copy.deepcopy(baseline_manifests)
        digest_mismatch_config = copy.deepcopy(baseline_config)
        mismatch_anchor_evidence_digests(
            digest_mismatch_manifests, digest_mismatch_config
        )
        digest_mismatch_report = execute(
            digest_mismatch_manifests,
            digest_mismatch_config,
        )
        digest_dimensional_codes = {
            issue.code
            for issue in digest_mismatch_report.results["DIMENSIONAL_CONSISTENCY"].issues
        }
        digest_evidence_codes = {
            issue.code
            for issue in digest_mismatch_report.results["EVIDENCE"].issues
        }
        digest_other = {
            name: result.status
            for name, result in digest_mismatch_report.results.items()
            if name not in {"DIMENSIONAL_CONSISTENCY", "EVIDENCE"} and result.status != "PASS"
        }
        if (
            clean.headline == "PASS"
            and digest_mismatch_report.results["DIMENSIONAL_CONSISTENCY"].status == "FAIL"
            and digest_dimensional_codes
            == {"DIM_SOURCE_EVIDENCE_DIGEST_MISMATCH"}
            and digest_mismatch_report.results["EVIDENCE"].status == "FAIL"
            and "STALE_SOURCE_DIGEST" in digest_evidence_codes
            and not digest_other
        ):
            passed += 1
            lines.append(
                "[PASS] DIMENSIONAL_CONSISTENCY anchor evidence digest mismatch | target=DIMENSIONAL_CONSISTENCY+EVIDENCE | "
                "clean=PASS planted=FAIL | "
                "code=DIM_SOURCE_EVIDENCE_DIGEST_MISMATCH"
            )
        else:
            failures.append(
                "DIMENSIONAL_CONSISTENCY anchor evidence digest mismatch: "
                f"DIMENSIONAL_CONSISTENCY={digest_mismatch_report.results['DIMENSIONAL_CONSISTENCY'].status}/"
                f"{sorted(digest_dimensional_codes)}, "
                f"EVIDENCE={digest_mismatch_report.results['EVIDENCE'].status}/"
                f"{sorted(digest_evidence_codes)}, extras={digest_other}"
            )
            lines.append("[FAIL] DIMENSIONAL_CONSISTENCY anchor evidence digest mismatch")

        def duplicate_stage(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            source_path = Path(paths["stage3_source"])
            for claim_id in (
                "duplicate_stage_first",
                "duplicate_stage_second",
            ):
                manifests.append(
                    _self_test_manifest(
                        source_path,
                        paths["source_digest"],
                        "stage029",
                        [],
                        [
                            _self_test_claim(
                                source_path,
                                paths["source_digest"],
                                claim_id,
                            )
                        ],
                    )
                )

        case(
            "C1 duplicate stage id",
            "C1",
            "FAIL",
            "DUPLICATE_STAGE",
            duplicate_stage,
            exact_target_code=True,
        )

        def append_symbol_free_stage029(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            source_path = Path(paths["stage3_source"])
            manifests.append(
                _self_test_manifest(
                    source_path,
                    paths["source_digest"],
                    "stage029",
                    [],
                    [
                        _self_test_claim(
                            source_path,
                            paths["source_digest"],
                            "stage029_identity",
                        )
                    ],
                )
            )

        clean_case(
            "DIMENSIONAL_CONSISTENCY symbol-free stage029 needs no dimension anchor",
            append_symbol_free_stage029,
        )

        def stage029_missing_anchor(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            append_symbol_free_stage029(manifests, config)
            stage29 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage029"
            )
            symbol = _self_test_symbol(
                Path(paths["stage3_source"]),
                paths["source_digest"],
                name="missing_dim",
                alias="missing_dim",
                qid="q.stage029.missing_dim",
                definition="here/define_missing_dim",
                dim={},
                raw=("0", "0", "0"),
            )
            symbol["dim_source"]["script_path"] = (
                "scripts/ledger_stage029_missing_sympy_audit.py"
            )
            stage29["symbols"].append(symbol)
            stage29["claims"].append(
                _self_test_claim(
                    Path(paths["stage3_source"]),
                    paths["source_digest"],
                    "define_missing_dim",
                    payload=_relation_payload(
                        "missing_dim", "missing_dim", "defines"
                    ),
                )
            )

        case(
            "DIMENSIONAL_CONSISTENCY stage029 symbol with absent anchor is not anchor-free",
            "DIMENSIONAL_CONSISTENCY",
            "UNSUPPORTED",
            "DIM_SOURCE_MISSING",
            stage029_missing_anchor,
            exact_target_code=True,
        )

        case(
            "C2 spectrum kernel mutation",
            "C2",
            "FAIL",
            "SPECTRUM_FIELD_DRIFT",
            lambda m, c: _find_consume(m, "stage002", "stage001/zero_mode")["as_consumed"]["payload"][
                "kernel"
            ].__setitem__("sympy", "2/ell"),
        )

        case(
            "C2 cas_equivalence on spectrum",
            "C2",
            "UNSUPPORTED",
            "INELIGIBLE_PAYLOAD_PAIR",
            lambda m, c: _find_consume(m, "stage002", "stage001/zero_mode").__setitem__(
                "check", "cas_equivalence"
            ),
        )

        case(
            "token_match on quantity",
            "C2",
            "FAIL",
            "TOKEN_MATCH_QUANTITY",
            lambda m, c: _find_consume(m, "stage003", "stage001/status_token").__setitem__(
                "producer_quantity_id", "q.core.x"
            ),
        )

        def opaque_consume(m: list[dict[str, Any]], digest: str) -> None:
            consume = _find_consume(m, "stage002", "stage001/define_x")
            consume.pop("as_consumed", None)
            consume["check"] = "opaque_quantity_match"
            consume["producer_quantity_id"] = "q.core.x"
            consume["producer_source_digest"] = digest

        total += 1
        opaque_pass_manifests = copy.deepcopy(baseline_manifests)
        opaque_consume(opaque_pass_manifests, paths["source_digest"])
        opaque_pass = execute(opaque_pass_manifests, copy.deepcopy(baseline_config))
        if opaque_pass.headline == "PASS":
            passed += 1
            lines.append(
                "[PASS] opaque quantity digest pin | target=C2 | clean=PASS planted=PASS"
            )
        else:
            failures.append(
                "opaque quantity digest pin passing fixture: "
                + str({name: opaque_pass.results[name].status for name in CHECK_NAMES})
            )
            lines.append("[FAIL] opaque quantity digest pin passing fixture")

        case(
            "opaque quantity digest drift",
            "C2",
            "FAIL",
            "OPAQUE_QUANTITY_PIN_DRIFT",
            lambda m, c: opaque_consume(m, "0" * 64),
        )

        def c1_dim_conflict(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            stage3 = next(x for x in m if x["stage_id"] == "stage003")
            stage3["symbols"].append(
                _self_test_symbol(
                    Path(paths["source"]),
                    paths["source_digest"],
                    name="x",
                    alias="x3",
                    qid="q.core.x",
                    definition="stage001/define_x",
                    dim={"L": "1"},
                    raw=("1", "0", "0"),
                )
            )

        case(
            "C1 quantity dimension conflict",
            "C1",
            "FAIL",
            "QUANTITY_ID_CONFLICT",
            c1_dim_conflict,
        )

        def c1_alias_collision(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            stage3 = next(x for x in m if x["stage_id"] == "stage003")
            stage3["symbols"].append(
                _self_test_symbol(
                    Path(paths["source"]),
                    paths["source_digest"],
                    name="x",
                    alias="z",
                    qid="q.other.z",
                    definition="here/define_z",
                    dim={},
                    raw=("0", "0", "0"),
                )
            )
            stage3["claims"].append(
                _self_test_claim(
                    Path(paths["source"]),
                    paths["source_digest"],
                    "define_z",
                    payload=_relation_payload("z", "z", "defines"),
                )
            )

        case(
            "C1 same-name without alias",
            "C1",
            "FAIL",
            "SAME_NAME_WITHOUT_ALIAS",
            c1_alias_collision,
        )

        def retired_claim(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            claim = _find_claim(m, "stage001", "define_x")
            claim["status"] = "RETIRED"
            claim["retired_by"] = "synthetic retirement"

        case(
            "C3 retired claim consumed",
            "C3",
            "FAIL",
            "NON_OPERATIVE_CLAIM",
            retired_claim,
        )

        def nonexported_claim(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            stage1 = next(x for x in m if x["stage_id"] == "stage001")
            stage1["exports"] = [item for item in stage1["exports"] if item["claim_id"] != "define_x"]

        case(
            "C3 non-exported claim consumed",
            "C3",
            "FAIL",
            "NON_EXPORTED_CLAIM",
            nonexported_claim,
        )

        case(
            "C7 expectation references missing tooth",
            "C3",
            "FAIL",
            "UNKNOWN_C7_TOOTH",
            lambda m, c: _find_consume(m, "stage002", "stage001/zero_mode")["c7_expect"].__setitem__(
                "expected_first_failure", "NOT_A_TOOTH"
            ),
        )

        case(
            "DIMENSIONAL_CONSISTENCY dimensional break",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIMENSIONAL_INHOMOGENEITY",
            lambda m, c: _find_claim(m, "stage002", "define_y")["payload"]["rhs"].__setitem__(
                "sympy", "x*ell"
            ),
        )

        total += 1
        rational = execute(copy.deepcopy(baseline_manifests), copy.deepcopy(baseline_config))
        if rational.results["DIMENSIONAL_CONSISTENCY"].status == "PASS":
            passed += 1
            lines.append("[PASS] DIMENSIONAL_CONSISTENCY rational-power sqrt(area)=ell | target=DIMENSIONAL_CONSISTENCY | clean=PASS planted=PASS")
        else:
            failures.append("DIMENSIONAL_CONSISTENCY rational-power fixture did not pass")
            lines.append("[FAIL] DIMENSIONAL_CONSISTENCY rational-power sqrt(area)=ell")

        def arbitrary_basis_symbol(
            source: Path,
            digest: str,
            *,
            name: str,
            alias: str,
            qid: str,
            definition: str,
            dim: dict[str, str],
            order: Sequence[str],
            raw: Sequence[str],
            locus: str,
        ) -> dict[str, Any]:
            symbol = _self_test_symbol(
                source,
                digest,
                name=name,
                alias=alias,
                qid=qid,
                definition=definition,
                dim=dim,
                raw=raw,
            )
            symbol["dim_source_order"] = list(order)
            symbol["dim_source"]["locus"] = locus
            return symbol

        four_axes = ("M", "L", "T", "E-charge")
        four_basis = {
            "id": "electromagnetic-MLTE",
            "axes": list(four_axes),
        }

        def arrange_four_axis(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            source = Path(paths["four_axis_source"])
            digest = paths["four_axis_source_digest"]
            symbols = [
                arbitrary_basis_symbol(
                    source,
                    digest,
                    name="A_E",
                    alias="A_E",
                    qid="q.stage038.a_e",
                    definition="here/define_a_e",
                    dim={"L": "1", "E-charge": "1"},
                    order=four_axes,
                    raw=("0", "1", "0", "1"),
                    locus="A_E_dim = (0, 1, 0, 1)",
                ),
                arbitrary_basis_symbol(
                    source,
                    digest,
                    name="q_T",
                    alias="q_T",
                    qid="q.stage038.q_t",
                    definition="here/define_q_t",
                    dim={"M": "1", "T": "-1"},
                    order=four_axes,
                    raw=("1", "0", "-1", "0"),
                    locus="q_T_dim = (1, 0, -1, 0)",
                ),
            ]
            claims = [
                _self_test_claim(
                    source,
                    digest,
                    "define_a_e",
                    payload=_relation_payload("A_E", "A_E", "defines"),
                ),
                _self_test_claim(
                    source,
                    digest,
                    "define_q_t",
                    payload=_relation_payload("q_T", "q_T", "defines"),
                ),
                _self_test_claim(
                    source,
                    digest,
                    "four_axis_composition",
                    kind="dimensional",
                    status="EARNED",
                    payload=_relation_payload(
                        "A_E*q_T**2", "A_E*q_T**2"
                    ),
                    expected_dim={
                        "M": "2",
                        "L": "1",
                        "T": "-2",
                        "E-charge": "1",
                    },
                ),
            ]
            stage = _self_test_manifest(
                source, digest, "stage038", symbols, claims
            )
            stage["dimension_basis"] = copy.deepcopy(four_basis)
            manifests.append(stage)

        def plant_four_axis_algebra(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            _find_claim(
                manifests, "stage038", "four_axis_composition"
            )["expected_dim"]["E-charge"] = "2"

        paired_case(
            "arbitrary 4-axis exact composition",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "EXPECTED_DIM_MISMATCH",
            arrange_four_axis,
            plant_four_axis_algebra,
            lambda report: report.dimension_bases["stage038"]
            == DimensionBasis("electromagnetic-MLTE", four_axes),
        )

        fractional_axes = ("stiffness", "length", "time")
        fractional_basis = {
            "id": "stiffness-length-time",
            "axes": list(fractional_axes),
        }

        def arrange_fractional_basis(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            source = Path(paths["fractional_basis_source"])
            digest = paths["fractional_basis_source_digest"]
            symbols = [
                arbitrary_basis_symbol(
                    source,
                    digest,
                    name="K",
                    alias="K",
                    qid="q.stage042.stiffness",
                    definition="here/define_k",
                    dim={"stiffness": "1"},
                    order=fractional_axes,
                    raw=("1", "0", "0"),
                    locus=(
                        "stiffness_dim = "
                        "(Fraction(1), Fraction(0), Fraction(0))"
                    ),
                ),
                arbitrary_basis_symbol(
                    source,
                    digest,
                    name="q",
                    alias="q",
                    qid="q.stage042.charge",
                    definition="here/define_q",
                    dim={
                        "stiffness": "1/2",
                        "length": "3/2",
                        "time": "-1",
                    },
                    order=fractional_axes,
                    raw=("1/2", "3/2", "-1"),
                    locus=(
                        "charge_dim = (Fraction(1, 2), "
                        "Fraction(3, 2), Fraction(-1))"
                    ),
                ),
            ]
            claims = [
                _self_test_claim(
                    source,
                    digest,
                    "define_k",
                    payload=_relation_payload("K", "K", "defines"),
                ),
                _self_test_claim(
                    source,
                    digest,
                    "define_q",
                    payload=_relation_payload("q", "q", "defines"),
                ),
                _self_test_claim(
                    source,
                    digest,
                    "fractional_composition",
                    kind="dimensional",
                    status="EARNED",
                    payload=_relation_payload("q**2/K", "q**2/K"),
                    expected_dim={"length": "3", "time": "-2"},
                ),
            ]
            stage = _self_test_manifest(
                source, digest, "stage042", symbols, claims
            )
            stage["dimension_basis"] = copy.deepcopy(fractional_basis)
            manifests.append(stage)

        def plant_fractional_algebra(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            _find_claim(
                manifests, "stage042", "fractional_composition"
            )["expected_dim"]["length"] = "5/2"

        paired_case(
            "non-L/M/T fractional exact composition",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "EXPECTED_DIM_MISMATCH",
            arrange_fractional_basis,
            plant_fractional_algebra,
            lambda report: report.dimension_bases["stage042"]
            == DimensionBasis("stiffness-length-time", fractional_axes),
        )

        def plant_fractional_addition_mismatch(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            _find_claim(
                manifests, "stage042", "fractional_composition"
            )["payload"] = _relation_payload("q + K", "q + K")

        def rendered_fractional_axes(report: CompositeReport) -> bool:
            details = [
                issue.detail
                for issue in report.results[
                    "DIMENSIONAL_CONSISTENCY"
                ].issues
                if issue.code == "DIMENSIONAL_INHOMOGENEITY"
            ]
            return (
                len(details) == 1
                and "addition dimension mismatch" in details[0]
                and all(
                    f"'{axis}':" in details[0]
                    for axis in fractional_axes
                )
            )

        paired_case(
            "non-L/M/T addition mismatch renders actual axes",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "DIMENSIONAL_INHOMOGENEITY",
            arrange_fractional_basis,
            plant_fractional_addition_mismatch,
            verify_planted=rendered_fractional_axes,
        )

        def plant_float_exponent(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage = next(
                item
                for item in manifests
                if item["stage_id"] == "stage042"
            )
            next(
                symbol
                for symbol in stage["symbols"]
                if symbol["parse_alias"] == "q"
            )["dim"]["stiffness"] = 0.5

        paired_case(
            "exact rational exponent rejects float",
            "SCHEMA",
            "FAIL",
            "SCHEMA_VALIDATION",
            arrange_fractional_basis,
            plant_float_exponent,
        )

        def arrange_shared_four_axis(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            producer_source = Path(paths["four_axis_source"])
            producer_digest = paths["four_axis_source_digest"]
            producer_symbol = arbitrary_basis_symbol(
                producer_source,
                producer_digest,
                name="shared",
                alias="shared",
                qid="q.shared.four_axis",
                definition="here/define_shared",
                dim={"L": "1", "E-charge": "1"},
                order=four_axes,
                raw=("0", "1", "0", "1"),
                locus="shared_dim = (0, 1, 0, 1)",
            )
            producer = _self_test_manifest(
                producer_source,
                producer_digest,
                "stage038",
                [producer_symbol],
                [
                    _self_test_claim(
                        producer_source,
                        producer_digest,
                        "define_shared",
                        payload=_relation_payload(
                            "shared", "shared", "defines"
                        ),
                    )
                ],
            )
            producer["dimension_basis"] = copy.deepcopy(four_basis)

            consumer_source = Path(paths["four_axis_peer_source"])
            consumer_digest = paths["four_axis_peer_source_digest"]
            consumer_symbol = arbitrary_basis_symbol(
                consumer_source,
                consumer_digest,
                name="shared",
                alias="shared",
                qid="q.shared.four_axis",
                definition="stage038/define_shared",
                dim={"L": "1", "E-charge": "1"},
                order=four_axes,
                raw=("0", "1", "0", "1"),
                locus="shared_dim = (0, 1, 0, 1)",
            )
            consumer = _self_test_manifest(
                consumer_source,
                consumer_digest,
                "stage039",
                [consumer_symbol],
                [],
            )
            consumer["dimension_basis"] = copy.deepcopy(four_basis)
            manifests.extend((producer, consumer))

        def plant_cross_basis(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage = next(
                item
                for item in manifests
                if item["stage_id"] == "stage039"
            )
            stage["dimension_basis"]["id"] = (
                "distinct-electromagnetic-MLTE"
            )

        paired_case(
            "cross-stage different basis is explicit",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "CROSS_BASIS_DIMENSION_COMPARISON",
            arrange_shared_four_axis,
            plant_cross_basis,
        )

        def plant_same_basis_drift(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage = next(
                item
                for item in manifests
                if item["stage_id"] == "stage039"
            )
            symbol = stage["symbols"][0]
            symbol["dim"] = {"L": "2", "E-charge": "1"}
            symbol["dim_source_tuple"] = ["0", "2", "0", "1"]
            symbol["dim_source"]["locus"] = (
                "shared_drift_dim = (0, 2, 0, 1)"
            )

        paired_case(
            "same-basis cross-stage agreement",
            "C1",
            "FAIL",
            "QUANTITY_ID_CONFLICT",
            arrange_shared_four_axis,
            plant_same_basis_drift,
        )

        def arrange_dim_equal_citation(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            producer_source = Path(paths["four_axis_source"])
            producer_digest = paths["four_axis_source_digest"]
            producer_symbol = arbitrary_basis_symbol(
                producer_source,
                producer_digest,
                name="shared",
                alias="shared",
                qid="q.stage038.dim_citation",
                definition="here/define_dim_citation",
                dim={"L": "1", "E-charge": "1"},
                order=four_axes,
                raw=("0", "1", "0", "1"),
                locus="shared_dim = (0, 1, 0, 1)",
            )
            producer = _self_test_manifest(
                producer_source,
                producer_digest,
                "stage038",
                [producer_symbol],
                [
                    _self_test_claim(
                        producer_source,
                        producer_digest,
                        "define_dim_citation",
                        payload=_relation_payload(
                            "shared", "shared", "defines"
                        ),
                    )
                ],
            )
            producer["dimension_basis"] = copy.deepcopy(four_basis)
            producer["exports"] = [
                {
                    "claim_id": "define_dim_citation",
                    "source_digest": producer_digest,
                    "c7_binding": _c7_binding(
                        paths["command"], "facet_dim_citation"
                    ),
                }
            ]

            consumer_source = Path(paths["four_axis_peer_source"])
            consumer_digest = paths["four_axis_peer_source_digest"]
            consumer = _self_test_manifest(
                consumer_source,
                consumer_digest,
                "stage039",
                [],
                [
                    _self_test_claim(
                        consumer_source,
                        consumer_digest,
                        "dim_citation_consumer",
                    )
                ],
            )
            consumer["dimension_basis"] = copy.deepcopy(four_basis)
            consumer["verification"]["teeth"] = [
                {
                    "predicate": "TOOTH_DIM_CITATION",
                    "mutation": "dim_citation",
                    "claim_ids": ["dim_citation_consumer"],
                    "evidence": _self_test_evidence(
                        consumer_source,
                        consumer_digest,
                        "TOOTH_DIM_CITATION",
                    ),
                }
            ]
            consumer["consumes"] = [
                {
                    "ref": "stage038/define_dim_citation",
                    "as_consumed_dim": {
                        "L": "1",
                        "E-charge": "1",
                    },
                    "check": "dim_equal",
                    "substitutions": [],
                    "c7_expect": _c7_expect(
                        "facet_dim_citation", "TOOTH_DIM_CITATION"
                    ),
                }
            ]
            manifests.extend((producer, consumer))

        def plant_dim_equal_cross_basis(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage = next(
                item
                for item in manifests
                if item["stage_id"] == "stage039"
            )
            stage["dimension_basis"]["id"] = (
                "distinct-electromagnetic-MLTE"
            )

        paired_case(
            "dim_equal citation rejects cross-basis comparison",
            "DIMENSIONAL_CONSISTENCY",
            "FAIL",
            "CROSS_BASIS_DIMENSION_COMPARISON",
            arrange_dim_equal_citation,
            plant_dim_equal_cross_basis,
        )

        def plant_dim_equal_drift(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            consume = _find_consume(
                manifests,
                "stage039",
                "stage038/define_dim_citation",
            )
            consume["as_consumed_dim"]["L"] = "2"

        paired_case(
            "dim_equal citation detects dimension drift",
            "C2",
            "FAIL",
            "DIM_CITATION_DRIFT",
            arrange_dim_equal_citation,
            plant_dim_equal_drift,
        )

        def unclassified(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            for manifest in m:
                manifest["knobs"] = [k for k in manifest["knobs"] if k["registry_row"] != "row_b"]

        case(
            "C5 unclassified register row",
            "C5",
            "FAIL",
            "UNCLASSIFIED_REGISTER_ROW",
            unclassified,
        )

        def owned_register_config(
            filename: str,
            extra_row: str,
            owner: str | None,
            row_class: str = "ACTION",
        ) -> dict[str, Any]:
            register = root / filename
            enters = f"synthetic ({owner})" if owner is not None else "—"
            register.write_text(
                "| Param | `[L,T,M]` dim | Enters | Class | Depends | Route |\n"
                "|---|---|---|---|---|---|\n"
                "| `row_a` | `1` | I-1 (001) | `ACTION` | — | — |\n"
                "| `row_b` | `1` | I-2 (002) | `ACTION` | — | — |\n"
                f"| `{extra_row}` | `1` | {enters} | `{row_class}` | — | — |\n"
            )
            config = copy.deepcopy(baseline_config)
            config["parameter_register"] = {
                "path": register.name,
                "sha256": hash_file(register),
                "format": "markdown_master_table",
            }
            config["closed_slice"] = False
            return config

        total += 1
        loaded_owner_config = owned_register_config(
            "loaded_owner_register.md", "row_extra", "001"
        )
        loaded_owner_report = execute(
            copy.deepcopy(baseline_manifests),
            loaded_owner_config,
            closed_slice=False,
        )
        loaded_owner_codes = {
            issue.code for issue in loaded_owner_report.results["C5"].issues
        }
        if (
            loaded_owner_report.headline == "FAIL"
            and loaded_owner_report.results["C5"].status == "FAIL"
            and "UNCLASSIFIED_REGISTER_ROW" in loaded_owner_codes
            and "REGISTER_COVERAGE_PARTIAL" not in loaded_owner_codes
            and all(
                result.status == "PASS"
                for name, result in loaded_owner_report.results.items()
                if name != "C5"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] open internal register gap owned by loaded stage | "
                "target=C5 | clean=PASS planted=FAIL | "
                "code=UNCLASSIFIED_REGISTER_ROW"
            )
        else:
            failures.append(
                "open internal register gap: "
                f"headline={loaded_owner_report.headline}, "
                f"statuses={ {name: result.status for name, result in loaded_owner_report.results.items()} }, "
                f"codes={sorted(loaded_owner_codes)}"
            )
            lines.append("[FAIL] open internal register gap owned by loaded stage")

        total += 1
        pending_owner_config = owned_register_config(
            "pending_owner_register.md", "row_pending", "999"
        )
        pending_owner_report = execute(
            copy.deepcopy(baseline_manifests),
            pending_owner_config,
            closed_slice=False,
        )
        pending_owner_codes = {
            issue.code for issue in pending_owner_report.results["C5"].issues
        }
        if (
            pending_owner_report.headline == "PARTIAL"
            and pending_owner_report.results["C5"].status == "PARTIAL"
            and "REGISTER_COVERAGE_PARTIAL" in pending_owner_codes
            and "UNCLASSIFIED_REGISTER_ROW" not in pending_owner_codes
            and all(
                result.status == "PASS"
                for name, result in pending_owner_report.results.items()
                if name != "C5"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] open pending register row owned by unloaded stage | "
                "target=C5 | clean=PASS planted=PARTIAL | "
                "code=REGISTER_COVERAGE_PARTIAL"
            )
        else:
            failures.append(
                "open pending register row: "
                f"headline={pending_owner_report.headline}, "
                f"statuses={ {name: result.status for name, result in pending_owner_report.results.items()} }, "
                f"codes={sorted(pending_owner_codes)}"
            )
            lines.append("[FAIL] open pending register row owned by unloaded stage")

        total += 1
        derived_owner_config = owned_register_config(
            "derived_owner_register.md",
            "row_derived",
            "001",
            "DERIVED",
        )
        derived_owner_report = execute(
            copy.deepcopy(baseline_manifests),
            derived_owner_config,
            closed_slice=False,
        )
        derived_owner_codes = {
            issue.code for issue in derived_owner_report.results["C5"].issues
        }
        if (
            derived_owner_report.headline == "PASS"
            and derived_owner_report.results["C5"].status == "PASS"
            and "UNCLASSIFIED_REGISTER_ROW" not in derived_owner_codes
            and all(
                result.status == "PASS"
                for result in derived_owner_report.results.values()
            )
        ):
            passed += 1
            lines.append(
                "[PASS] derived non-knob row owned by loaded stage | "
                "target=C5 | clean=PASS planted=PASS"
            )
        else:
            failures.append(
                "derived non-knob register row: "
                f"headline={derived_owner_report.headline}, "
                f"statuses={ {name: result.status for name, result in derived_owner_report.results.items()} }, "
                f"codes={sorted(derived_owner_codes)}"
            )
            lines.append("[FAIL] derived non-knob row owned by loaded stage")

        total += 1
        gap_config = owned_register_config(
            "gap_register.md",
            "m_defect",
            None,
            "GAP",
        )
        gap_report = execute(
            copy.deepcopy(baseline_manifests),
            gap_config,
            closed_slice=False,
        )
        gap_codes = {
            issue.code for issue in gap_report.results["C5"].issues
        }
        if (
            gap_report.headline == "PARTIAL"
            and gap_report.results["C5"].status == "PARTIAL"
            and "REGISTER_COVERAGE_PARTIAL" in gap_codes
            and "UNCLASSIFIED_REGISTER_ROW" not in gap_codes
            and all(
                result.status == "PASS"
                for name, result in gap_report.results.items()
                if name != "C5"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] unattributable GAP register row stays pending | "
                "target=C5 | clean=PASS planted=PARTIAL | "
                "code=REGISTER_COVERAGE_PARTIAL"
            )
        else:
            failures.append(
                "unattributable GAP register row: "
                f"headline={gap_report.headline}, "
                f"statuses={ {name: result.status for name, result in gap_report.results.items()} }, "
                f"codes={sorted(gap_codes)}"
            )
            lines.append("[FAIL] unattributable GAP register row stays pending")

        case(
            "C5 orphan inherit",
            "C5",
            "FAIL",
            "ORPHAN_INHERIT",
            lambda m, c: next(
                k
                for k in next(x for x in m if x["stage_id"] == "stage002")["knobs"]
                if k["action"] == "inherited"
            ).__setitem__("knob_id", "k.orphan"),
        )

        case(
            "C5 invalid discharge",
            "C5",
            "FAIL",
            "INVALID_DISCHARGE",
            lambda m, c: next(
                k
                for k in next(x for x in m if x["stage_id"] == "stage002")["knobs"]
                if k["action"] == "discharged"
            ).__setitem__("discharge_evidence", "here/rational_power"),
        )

        case(
            "C5 low endpoint drift",
            "C5",
            "FAIL",
            "RANGE_INTERNAL_DRIFT",
            lambda m, c: _find_claim(m, "stage001", "count_range")["payload"].__setitem__("low", 1),
        )

        case(
            "C5 high endpoint drift",
            "C5",
            "FAIL",
            "RANGE_INTERNAL_DRIFT",
            lambda m, c: _find_claim(m, "stage001", "count_range")["payload"].__setitem__("high", 4),
        )

        case(
            "C5 pinned register digest drift",
            "C5",
            "FAIL",
            "REGISTER_DIGEST_CHANGED",
            lambda m, c: c["parameter_register"].__setitem__("sha256", "0" * 64),
        )

        def cycle(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            stage1 = next(x for x in m if x["stage_id"] == "stage001")
            stage1["symbols"].append(
                _self_test_symbol(
                    Path(paths["stage1_source"]),
                    paths["source_digest"],
                    name="y",
                    alias="y",
                    qid="q.stage2.y",
                    definition="stage002/define_y",
                    dim={},
                    raw=("0", "0", "0"),
                )
            )
            stage1["consumes"].append(
                {
                    "ref": "stage002/define_y",
                    "as_consumed": _typed_relation("y", "x", "defines"),
                    "check": "cas_equivalence",
                    "substitutions": [],
                    "c7_expect": _c7_expect("facet_y", "TOOTH_CYCLE"),
                }
            )

        case(
            "C6 two-cycle",
            "C6",
            "FAIL",
            "DEPENDENCY_CYCLE",
            cycle,
        )

        case(
            "evidence stale digest",
            "EVIDENCE",
            "FAIL",
            "STALE_SOURCE_DIGEST",
            lambda m, c: next(
                s
                for s in next(x for x in m if x["stage_id"] == "stage001")["symbols"]
                if s["parse_alias"] == "a_B"
            )["evidence"].__setitem__("source_digest", "0" * 64),
        )

        def cite_mathematica_output(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage1 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )
            stage1["verification"]["mathematica_script"] = paths[
                "mathematica_script"
            ]
            stage1["verification"]["mathematica_output"] = {
                "path": paths["mathematica_output"],
                "sha256": paths["mathematica_output_digest"],
            }

        def mathematica_output_count_is_rendered(
            report: CompositeReport,
        ) -> bool:
            return (
                report.coverage.mathematica_outputs_checked == 1
                and "mathematica_outputs_checked=1" in report.matrix_line()
                and "mathematica_outputs_checked=0" in clean.matrix_line()
                and render_report(report) != render_report(clean)
            )

        clean_case(
            "Mathematica output live digest count is rendered",
            cite_mathematica_output,
            mathematica_output_count_is_rendered,
        )

        def stale_mathematica_output(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            cite_mathematica_output(manifests, config)
            next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )["verification"]["mathematica_output"]["sha256"] = "0" * 64

        case(
            "Mathematica output stale digest",
            "EVIDENCE",
            "FAIL",
            "STALE_MATHEMATICA_OUTPUT_DIGEST",
            stale_mathematica_output,
            exact_target_code=True,
        )

        def missing_mathematica_output(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            stage1 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )
            stage1["verification"]["mathematica_script"] = paths[
                "missing_mathematica_script"
            ]
            stage1["verification"]["mathematica_output"] = {
                "path": "mathematica/out/missing_math.out",
                "sha256": "0" * 64,
            }

        case(
            "Mathematica output missing file",
            "EVIDENCE",
            "FAIL",
            "MISSING_MATHEMATICA_OUTPUT",
            missing_mathematica_output,
            exact_target_code=True,
        )

        def wrong_stage_mathematica_output(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            cite_mathematica_output(manifests, config)
            stage1 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )
            stage1["verification"]["mathematica_output"] = {
                "path": paths["other_mathematica_output"],
                "sha256": paths["other_mathematica_output_digest"],
            }

        case(
            "Mathematica output from another script is rejected",
            "EVIDENCE",
            "FAIL",
            "MATHEMATICA_OUTPUT_SCRIPT_MISMATCH",
            wrong_stage_mathematica_output,
            exact_target_code=True,
        )

        def empty_mathematica_output(
            manifests: list[dict[str, Any]], _: dict[str, Any]
        ) -> None:
            stage1 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )
            stage1["verification"]["mathematica_script"] = paths[
                "empty_mathematica_script"
            ]
            stage1["verification"]["mathematica_output"] = {
                "path": paths["empty_mathematica_output"],
                "sha256": paths["empty_mathematica_output_digest"],
            }

        case(
            "Mathematica empty output is rejected",
            "EVIDENCE",
            "FAIL",
            "EMPTY_MATHEMATICA_OUTPUT",
            empty_mathematica_output,
            exact_target_code=True,
        )

        case(
            "Mathematica absolute output path is rejected",
            "EVIDENCE",
            "FAIL",
            "UNSAFE_MATHEMATICA_OUTPUT_PATH",
            lambda m, c: (
                cite_mathematica_output(m, c),
                next(
                    manifest
                    for manifest in m
                    if manifest["stage_id"] == "stage001"
                )["verification"]["mathematica_output"].__setitem__(
                    "path", paths["mathematica_output_absolute"]
                ),
            ),
            exact_target_code=True,
        )

        case(
            "Mathematica parent-escape output path is rejected",
            "EVIDENCE",
            "FAIL",
            "UNSAFE_MATHEMATICA_OUTPUT_PATH",
            lambda m, c: (
                cite_mathematica_output(m, c),
                next(
                    manifest
                    for manifest in m
                    if manifest["stage_id"] == "stage001"
                )["verification"]["mathematica_output"].__setitem__(
                    "path", "../selftest_math.out"
                ),
            ),
            exact_target_code=True,
        )

        def missing_mathematica_output_field(
            manifests: list[dict[str, Any]], config: dict[str, Any]
        ) -> None:
            cite_mathematica_output(manifests, config)
            stage1 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage001"
            )
            stage2 = next(
                manifest
                for manifest in manifests
                if manifest["stage_id"] == "stage002"
            )
            stage2["verification"]["mathematica_script"] = paths[
                "mathematica_script"
            ]
            stage2["verification"]["mathematica_output"] = copy.deepcopy(
                stage1["verification"]["mathematica_output"]
            )
            del stage1["verification"]["mathematica_output"]

        case(
            "Mathematica output field required once citation mode is active",
            "EVIDENCE",
            "FAIL",
            "MISSING_MATHEMATICA_OUTPUT_CITATION",
            missing_mathematica_output_field,
            exact_target_code=True,
        )

        def prose_engine(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            evidence = _find_claim(m, "stage001", "grid_landing")["evidence"]
            evidence["method"] = "prose_only"
            evidence["engine"] = "sympy"

        case(
            "evidence prose_only with engine:sympy",
            "SCHEMA",
            "FAIL",
            "SCHEMA_VALIDATION",
            prose_engine,
        )

        case(
            "adjudication bucket sum",
            "ADJUDICATION",
            "FAIL",
            "BUCKET_COUNT_SUM",
            lambda m, c: _find_claim(m, "stage001", "grid_landing")["payload"]["bucket_counts"].__setitem__(
                "PASS", 1
            ),
        )

        case(
            "adjudication axes product",
            "ADJUDICATION",
            "FAIL",
            "AXIS_CARDINALITY_PRODUCT",
            lambda m, c: _find_claim(m, "stage001", "grid_landing")["payload"]["axes"].__setitem__(
                "sign", 3
            ),
        )

        def kind_payload_coupling(m: list[dict[str, Any]], _: dict[str, Any]) -> None:
            claim = _find_claim(m, "stage001", "zero_mode")
            claim["kind"] = "identity"

        case(
            "claim kind-payload coupling",
            "SCHEMA",
            "FAIL",
            "SCHEMA_VALIDATION",
            kind_payload_coupling,
        )

        case(
            "C7 absent metadata caps slice at PARTIAL",
            "C7",
            "PARTIAL",
            "C7_EDGE_UNCOVERED",
            lambda m, c: _find_consume(m, "stage002", "stage001/zero_mode").pop("c7_expect"),
        )

        case(
            "C7 decorative dependency",
            "C7",
            "FAIL",
            "DECORATIVE_DEPENDENCY",
            lambda m, c: _find_export(m, "stage001", "zero_mode")["c7_binding"].__setitem__(
                "mutation_command", paths["command"] + " --decorative"
            ),
        )

        case(
            "C7 undeclared dependency",
            "C7",
            "FAIL",
            "UNDECLARED_DEPENDENCY",
            lambda m, c: _find_export(m, "stage001", "zero_mode")["c7_binding"].__setitem__(
                "mutation_command", paths["command"] + " --undeclared"
            ),
        )

        def dangling_consumer(ref: str) -> dict[str, Any]:
            manifest = _self_test_manifest(
                Path(paths["source"]),
                paths["source_digest"],
                "stage004",
                [],
                [
                    _self_test_claim(
                        Path(paths["source"]),
                        paths["source_digest"],
                        "stage4_identity",
                    )
                ],
            )
            manifest["consumes"] = [
                {
                    "ref": ref,
                    "as_consumed": _typed_relation("1", "1"),
                    "check": "cas_equivalence",
                    "substitutions": [],
                }
            ]
            return manifest

        auto_config = copy.deepcopy(baseline_config)
        auto_config["closed_slice"] = False
        auto_clean = execute(
            copy.deepcopy(baseline_manifests),
            copy.deepcopy(auto_config),
            closed_slice=None,
        )

        total += 1
        lone_stage = copy.deepcopy(
            next(
                manifest
                for manifest in baseline_manifests
                if manifest["stage_id"] == "stage003"
            )
        )
        lone_stage["consumes"] = []
        lone_stage["claims"].append(
            _self_test_claim(
                Path(paths["source"]),
                paths["source_digest"],
                "local_count_range",
                kind="range",
                status="EARNED",
                payload_kind="record_range",
                payload={
                    "low": 0,
                    "high": 0,
                    "spread": 0,
                    "convention_axes": [
                        {
                            "axis_id": "baseline",
                            "choices": [
                                {
                                    "token": "default",
                                    "low_delta": 0,
                                    "high_delta": 0,
                                }
                            ],
                        }
                    ],
                    "components": {
                        "empty": {"low": 0, "high": 0}
                    },
                },
            )
        )
        lone_config = copy.deepcopy(pending_owner_config)
        lone_config["range_claim_ref"] = "stage003/local_count_range"
        lone_config["closed_slice"] = False
        lone_report = execute(
            [lone_stage],
            lone_config,
            closed_slice=None,
        )
        lone_codes = {
            issue.code for issue in lone_report.results["C5"].issues
        }
        if (
            lone_report.headline == "PARTIAL"
            and lone_report.results["C5"].status == "PARTIAL"
            and "REGISTER_COVERAGE_PARTIAL" in lone_codes
            and "UNCLASSIFIED_REGISTER_ROW" not in lone_codes
            and all(
                result.status == "PASS"
                for name, result in lone_report.results.items()
                if name != "C5"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] ref-empty incomplete slice stays auto-open | target=C5 | "
                "clean=PASS planted=PARTIAL | code=REGISTER_COVERAGE_PARTIAL"
            )
        else:
            failures.append(
                "ref-empty incomplete slice: "
                f"headline={lone_report.headline}, "
                f"statuses={ {name: result.status for name, result in lone_report.results.items()} }, "
                f"codes={sorted(lone_codes)}"
            )
            lines.append("[FAIL] ref-empty incomplete slice stays auto-open")

        total += 1
        explicit_open_manifests = copy.deepcopy(baseline_manifests)
        explicit_open_manifests.append(dangling_consumer("stage999/missing_claim"))
        explicit_open = execute(
            explicit_open_manifests,
            copy.deepcopy(auto_config),
            closed_slice=False,
        )
        explicit_open_codes = {
            issue.code for issue in explicit_open.results["C3"].issues
        }
        if (
            auto_clean.headline == "PASS"
            and explicit_open.headline == "PARTIAL"
            and explicit_open.results["C3"].status == "PARTIAL"
            and "ABSENT_PRODUCER" in explicit_open_codes
            and all(
                result.status == "PASS"
                for name, result in explicit_open.results.items()
                if name != "C3"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] open cross-stage dangling reference | target=C3 | "
                "clean=PASS planted=PARTIAL | code=ABSENT_PRODUCER"
            )
        else:
            failures.append(
                "open cross-stage dangling reference: "
                f"clean={auto_clean.headline}, planted={explicit_open.headline}, "
                f"statuses={{name: result.status for name, result in explicit_open.results.items()}}, "
                f"codes={sorted(explicit_open_codes)}"
            )
            lines.append("[FAIL] open cross-stage dangling reference")

        total += 1
        explicit_open_internal = copy.deepcopy(baseline_manifests)
        _find_claim(
            explicit_open_internal, "stage002", "define_y"
        )["payload"]["rhs"]["sympy"] = "x*ell"
        open_internal = execute(
            explicit_open_internal,
            copy.deepcopy(auto_config),
            closed_slice=False,
        )
        open_internal_codes = {
            issue.code for issue in open_internal.results["DIMENSIONAL_CONSISTENCY"].issues
        }
        if (
            auto_clean.headline == "PASS"
            and open_internal.headline == "FAIL"
            and open_internal.results["DIMENSIONAL_CONSISTENCY"].status == "FAIL"
            and "DIMENSIONAL_INHOMOGENEITY" in open_internal_codes
            and all(
                result.status == "PASS"
                for name, result in open_internal.results.items()
                if name != "DIMENSIONAL_CONSISTENCY"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] open internal dimensional defect | target=DIMENSIONAL_CONSISTENCY | "
                "clean=PASS planted=FAIL | code=DIMENSIONAL_INHOMOGENEITY"
            )
        else:
            failures.append(
                "open internal dimensional defect: "
                f"clean={auto_clean.headline}, planted={open_internal.headline}, "
                f"statuses={ {name: result.status for name, result in open_internal.results.items()} }, "
                f"codes={sorted(open_internal_codes)}"
            )
            lines.append("[FAIL] open internal dimensional defect")

        producer = _self_test_manifest(
            Path(paths["source"]),
            paths["source_digest"],
            "stage005",
            [],
            [
                _self_test_claim(
                    Path(paths["source"]),
                    paths["source_digest"],
                    "producer_placeholder",
                )
            ],
        )
        auto_complete_manifests = copy.deepcopy(baseline_manifests)
        auto_complete_manifests.extend(
            [dangling_consumer("stage005/missing_claim"), producer]
        )

        total += 1
        auto_complete = execute(
            copy.deepcopy(auto_complete_manifests),
            copy.deepcopy(auto_config),
            closed_slice=None,
        )
        auto_complete_codes = {
            issue.code for issue in auto_complete.results["C3"].issues
        }
        if (
            auto_clean.headline == "PASS"
            and auto_complete.headline == "FAIL"
            and auto_complete.results["C3"].status == "FAIL"
            and auto_complete.results["C7"].status == "SKIPPED"
            and "ABSENT_PRODUCER" in auto_complete_codes
            and all(
                result.status == "PASS"
                for name, result in auto_complete.results.items()
                if name not in {"C3", "C7"}
            )
        ):
            passed += 1
            lines.append(
                "[PASS] auto-closed complete producer set | target=C3 | "
                "clean=PASS planted=FAIL | code=ABSENT_PRODUCER"
            )
        else:
            failures.append(
                "auto-closed complete producer set: "
                f"clean={auto_clean.headline}, planted={auto_complete.headline}, "
                f"statuses={ {name: result.status for name, result in auto_complete.results.items()} }, "
                f"codes={sorted(auto_complete_codes)}"
            )
            lines.append("[FAIL] auto-closed complete producer set")

        total += 1
        auto_incomplete_manifests = [
            manifest
            for manifest in copy.deepcopy(auto_complete_manifests)
            if manifest["stage_id"] != "stage005"
        ]
        auto_incomplete = execute(
            auto_incomplete_manifests,
            copy.deepcopy(auto_config),
            closed_slice=None,
        )
        auto_incomplete_codes = {
            issue.code for issue in auto_incomplete.results["C3"].issues
        }
        if (
            auto_clean.headline == "PASS"
            and auto_incomplete.headline == "PARTIAL"
            and auto_incomplete.results["C3"].status == "PARTIAL"
            and "ABSENT_PRODUCER" in auto_incomplete_codes
            and all(
                result.status == "PASS"
                for name, result in auto_incomplete.results.items()
                if name != "C3"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] auto-open incomplete producer set | target=C3 | "
                "clean=PASS planted=PARTIAL | code=ABSENT_PRODUCER"
            )
        else:
            failures.append(
                "auto-open incomplete producer set: "
                f"clean={auto_clean.headline}, planted={auto_incomplete.headline}, "
                f"statuses={ {name: result.status for name, result in auto_incomplete.results.items()} }, "
                f"codes={sorted(auto_incomplete_codes)}"
            )
            lines.append("[FAIL] auto-open incomplete producer set")

        total += 1
        example = load_json(Path(__file__).with_name("examples") / "stage_manifest_v2_example.json")
        errors = list(Draft202012Validator(schema).iter_errors(example))
        if not errors:
            passed += 1
            lines.append("[PASS] v2.1 example schema validation | target=SCHEMA | observed=PASS")
        else:
            failures.append(f"example schema errors: {[error.message for error in errors]}")
            lines.append("[FAIL] v2.1 example schema validation")

    lines.append(f"SUMMARY: {passed}/{total} fixtures PASS")
    if failures:
        lines.append("FAILURES:")
        lines.extend(f"- {failure}" for failure in failures)
    print("\n".join(lines))
    return 0 if not failures else 1


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="run deterministic isolated mutation fixtures")
    parser.add_argument("--manifest-dir", type=Path, default=DEFAULT_MANIFEST_DIR)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    slice_group = parser.add_mutually_exclusive_group()
    slice_group.add_argument("--closed-slice", action="store_true")
    slice_group.add_argument(
        "--open-slice",
        action="store_true",
        help="explicitly treat unresolved external producer stages as an open slice",
    )
    parser.add_argument("--recover-dims", type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--recover-dims-stdin",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--dim-order",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--dim-expression",
        action="append",
        default=[],
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)
    if args.recover_dims is not None:
        if args.dim_order is None:
            print("--recover-dims requires --dim-order", file=sys.stderr)
            return 2
        try:
            raw_order = (
                json.loads(args.dim_order)
                if args.dim_order.lstrip().startswith("[")
                else args.dim_order
            )
            payload = _live_dim_recovery(
                args.recover_dims.resolve(),
                source_order_axes(raw_order),
                args.dim_expression,
                sys.stdin.buffer.read() if args.recover_dims_stdin else None,
            )
        except Exception as exc:
            print(f"DIM RECOVERY FAIL: {exc}", file=sys.stderr)
            return 2
        print("DIM_RECOVERY_JSON=" + json.dumps(payload, sort_keys=True))
        return 0
    if args.self_test:
        return run_self_test()

    paths = sorted(args.manifest_dir.glob("*.json"))
    manifests: list[dict[str, Any]] = []
    load_errors: list[str] = []
    for path in paths:
        try:
            manifests.append(load_json(path))
        except (OSError, json.JSONDecodeError) as exc:
            load_errors.append(f"{path}: {exc}")
    if load_errors:
        print("\n".join(f"LOAD FAIL: {error}" for error in load_errors), file=sys.stderr)
        return 1
    config = load_json(args.config)
    closed_slice_override = True if args.closed_slice else False if args.open_slice else None
    report = CompositeChecker(
        manifests,
        config,
        roots=(PROJECT_ROOT, PROJECT_ROOT.parent.parent, Path.cwd()),
        closed_slice=closed_slice_override,
    ).run()
    text = render_report(report)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(text)
    print(f"{report.headline}: {report.matrix_line()}")
    for name in CHECK_NAMES:
        print(f"{name}: {report.results[name].status}")
    print(f"REPORT: {args.report}")
    return 1 if report.headline in BAD_STATUSES else 0


if __name__ == "__main__":
    raise SystemExit(main())
