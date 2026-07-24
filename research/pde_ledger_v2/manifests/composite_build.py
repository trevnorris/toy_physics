#!/usr/bin/env python3
"""Stage-manifest v2.1 composite integration checker.

Pipeline:
  Draft 2020-12 schema -> indexes/exact parsing -> evidence/adjudication
  integrity -> IMPORT-COMPLETENESS -> C1..C6 -> coverage -> C7 mutations.

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
from sympy.parsing.sympy_parser import parse_expr, standard_transformations


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = Path(__file__).with_name("stage_manifest_schema_v2.json")
DEFAULT_CONFIG = Path(__file__).with_name("composite_config.json")
DEFAULT_MANIFEST_DIR = Path(__file__).with_name("stages")
DEFAULT_REPORT = Path(__file__).with_name("reports") / "composite_report.md"
CHECK_NAMES = (
    "SCHEMA",
    "EVIDENCE",
    "ADJUDICATION",
    "IMPORT",
    "C1",
    "C2",
    "C3",
    "C4",
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
    tuple[Path, str, str, tuple[str, ...]], dict[str, Any]
] = {}


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


@dataclass
class CompositeReport:
    results: dict[str, CheckResult]
    coverage: Coverage

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
            f"closure={closure}"
        )


class UnsupportedExpression(Exception):
    pass


def _live_dim_recovery(
    path: Path, order: str, extra_expressions: Sequence[str] = ()
) -> dict[str, Any]:
    """Import ``path`` and evaluate source dimension algebra with its live Dim.

    This function runs only in the dedicated recovery subprocess.  AST
    validation limits evaluation to names, exact integer literals, containers,
    ``Dim(...)``, and the operators whose semantics must come from the audit
    script's own objects.
    """

    tree = ast.parse(path.read_text())
    module_name = f"_composite_dim_source_{hashlib.sha256(str(path).encode()).hexdigest()[:16]}"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise UnsupportedExpression(f"cannot import dim source {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)

    dim_type = getattr(module, "Dim", None)
    if not isinstance(dim_type, type):
        raise UnsupportedExpression(f"{path} does not expose a live Dim class")

    records: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    binding_counts: dict[str, int] = defaultdict(int)
    binding_expressions: dict[str, list[str]] = defaultdict(list)
    unsupported: set[str] = set()
    expression_envs: list[dict[str, Any]] = []

    def vector(value: Any) -> tuple[str, str, str] | None:
        if not isinstance(value, dim_type):
            return None
        try:
            raw = tuple(getattr(value, axis.lower()) for axis in order)
        except (AttributeError, TypeError) as exc:
            raise UnsupportedExpression(
                f"live Dim value does not expose {order} components: {exc}"
            ) from exc
        if len(raw) != 3:
            raise UnsupportedExpression(f"live Dim value has {len(raw)} components")
        rendered = tuple(str(item) for item in raw)
        try:
            tuple(Fraction(item) for item in rendered)
        except (ValueError, ZeroDivisionError) as exc:
            raise UnsupportedExpression(
                f"live Dim components are not exact rationals: {rendered}"
            ) from exc
        return rendered  # type: ignore[return-value]

    def add_tuple_record(
        name: str,
        recovered: tuple[str, str, str],
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
                and len(node.args) in {0, 3}
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
        if not admitted(expression):
            unsupported.add(canonical_key)
            continue
        recovered_values: set[tuple[str, str, str]] = set()
        for env in expression_envs:
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


def normalized_dim(raw: Mapping[str, str] | None) -> tuple[Fraction, Fraction, Fraction]:
    raw = raw or {}
    return tuple(Fraction(raw.get(axis, "0")) for axis in "LMT")  # type: ignore[return-value]


def dim_json(dim: Sequence[Fraction]) -> dict[str, str]:
    result: dict[str, str] = {}
    for axis, value in zip("LMT", dim):
        if value:
            result[axis] = str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    return result


def add_dim(a: Sequence[Fraction], b: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(x + y for x, y in zip(a, b))


def sub_dim(a: Sequence[Fraction], b: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(x - y for x, y in zip(a, b))


def scale_dim(a: Sequence[Fraction], power: Fraction) -> tuple[Fraction, ...]:
    return tuple(x * power for x in a)


ZERO_DIM = normalized_dim({})


def resolve_path(value: str, roots: Sequence[Path]) -> Path | None:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate if candidate.exists() else None
    for root in roots:
        resolved = root / candidate
        if resolved.exists():
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
    ) -> None:
        self.manifests = list(manifests)
        self.config = config
        self.schema = schema or json.loads(SCHEMA_PATH.read_text())
        self.roots = list(roots or (PROJECT_ROOT, PROJECT_ROOT.parent.parent, Path.cwd()))
        config_override = True if config.get("closed_slice") is True else None
        self.closed_slice_override = config_override if closed_slice is None else closed_slice
        self.closed_slice = False
        self.results = {name: CheckResult(name) for name in CHECK_NAMES}
        self.coverage = Coverage()
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
            return CompositeReport(self.results, self.coverage)
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
        self.check_c4()
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
        return CompositeReport(self.results, self.coverage)

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
            for stage_id, symbol in appearances:
                signature = symbol.get("signature")
                current = (
                    normalized_dim(symbol.get("dim")),
                    tuple(sorted(symbol.get("assumptions", []))),
                    json.dumps(signature, sort_keys=True),
                )
                if baseline is None:
                    baseline = current
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
                        producer_dims = [
                            normalized_dim(symbol.get("dim"))
                            for symbol in self.stages[producer_stage].get("symbols", [])
                            if symbol.get("definition_ref") in {ref, f"here/{ref.split('/', 1)[1]}"}
                        ]
                        if not producer_dims or normalized_dim(consume["as_consumed_dim"]) not in producer_dims:
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
        if expr.is_Number:
            return ZERO_DIM
        if isinstance(expr, sp.Symbol):
            symbol = self.symbols_by_stage_alias.get(stage_id, {}).get(expr.name)
            if not symbol:
                raise UnsupportedExpression(f"dimension of unknown symbol {expr}")
            return normalized_dim(symbol.get("dim"))
        if isinstance(expr, sp.Add):
            dims = [self.dimension_of(stage_id, arg) for arg in expr.args]
            nonzero_dims = [dim for arg, dim in zip(expr.args, dims) if arg != 0]
            if nonzero_dims and any(dim != nonzero_dims[0] for dim in nonzero_dims[1:]):
                raise ValueError(f"addition dimension mismatch {expr}: {list(map(dim_json, dims))}")
            return nonzero_dims[0] if nonzero_dims else ZERO_DIM
        if isinstance(expr, sp.Mul):
            total = ZERO_DIM
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
            expected = [normalized_dim(dim) for dim in signature["domain"]]
            if actual != expected:
                raise ValueError(f"signature domain mismatch {alias}")
            return normalized_dim(signature["codomain"])
        if expr.func in {sp.exp, sp.log, sp.sin, sp.cos, sp.tan, sp.sinh, sp.cosh, sp.tanh, sp.atan}:
            if any(self.dimension_of(stage_id, arg) != ZERO_DIM for arg in expr.args):
                raise ValueError(f"transcendental argument is dimensional: {expr}")
            return ZERO_DIM
        if expr.func == sp.Abs:
            return self.dimension_of(stage_id, expr.args[0])
        raise UnsupportedExpression(f"dimension rule missing for {type(expr).__name__}: {expr}")

    def recover_dim_order_and_tuples(
        self, path: Path, extra_expressions: Sequence[str] = ()
    ) -> tuple[str, dict[str, tuple[Fraction, ...]]]:
        try:
            tree = ast.parse(path.read_text())
        except Exception as exc:
            raise UnsupportedExpression(f"cannot parse dim source {path}: {exc}") from exc
        field_orders: set[str] = set()
        init_orders: set[str] = set()
        doc_orders: set[str] = set()
        for child in ast.walk(tree):
            if not isinstance(child, ast.ClassDef) or child.name != "Dim":
                continue
            fields = [
                item.target.id.upper()
                for item in child.body
                if isinstance(item, ast.AnnAssign)
                and isinstance(item.target, ast.Name)
                and item.target.id.lower() in {"l", "m", "t"}
            ]
            if len(fields) >= 3:
                field_orders.add("".join(fields[:3]))
            for item in child.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == "__init__":
                    args = [arg.arg.upper() for arg in item.args.args if arg.arg != "self"]
                    args = [arg for arg in args if arg in {"L", "M", "T"}]
                    if len(args) >= 3:
                        init_orders.add("".join(args[:3]))
            doc = ast.get_docstring(child) or ""
            match = re.search(
                r"[\[\{\(]\s*([LMT])\s*,\s*([LMT])\s*,\s*([LMT])\s*[\]\}\)]",
                doc,
                re.I,
            )
            if match:
                doc_orders.add("".join(group.upper() for group in match.groups()))

        structural_orders = field_orders | init_orders | doc_orders
        if len(structural_orders) > 1:
            raise UnsupportedExpression(
                "Dim field/doc/init order conflict: "
                f"fields={sorted(field_orders)}, init={sorted(init_orders)}, "
                f"doc={sorted(doc_orders)}"
            )
        order = next(iter(structural_orders), None)
        if order not in {"LMT", "LTM", "MLT", "MTL", "TLM", "TML"}:
            raise UnsupportedExpression("Dim order not recoverable from source structure")
        canonical_expressions = tuple(sorted(set(extra_expressions)))
        cache_key = (path.resolve(), hash_file(path), order, canonical_expressions)
        payload = LIVE_DIM_RECOVERY_CACHE.get(cache_key)
        if payload is None:
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--recover-dims",
                str(path),
                "--dim-order",
                order,
            ]
            for expression in canonical_expressions:
                command.extend(("--dim-expression", expression))
            try:
                completed = subprocess.run(
                    command,
                    cwd=path.parent,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
            except (OSError, subprocess.TimeoutExpired) as exc:
                raise UnsupportedExpression(f"live Dim recovery failed for {path}: {exc}") from exc
            marker = "DIM_RECOVERY_JSON="
            payload_line = next(
                (
                    line[len(marker) :]
                    for line in reversed(completed.stdout.splitlines())
                    if line.startswith(marker)
                ),
                None,
            )
            if completed.returncode != 0 or payload_line is None:
                detail = completed.stderr.strip() or completed.stdout.strip()
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

    def check_c4(self) -> None:
        result = self.results["C4"]
        source_cache: dict[Path, tuple[str, dict[str, tuple[Fraction, ...]]]] = {}
        for stage_id, manifest in self.stages.items():
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
                            f"{stage_id}/{claim.get('id')}: {dim_json(lhs_dim)} != {dim_json(rhs_dim)}",
                        )
                    if "expected_dim" in claim:
                        expected = normalized_dim(claim["expected_dim"])
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
                path = resolve_path(source.get("script_path", ""), self.roots)
                if path is None:
                    result.add(
                        "UNSUPPORTED",
                        "DIM_SOURCE_MISSING",
                        f"{stage_id}:{symbol.get('parse_alias')}:{source.get('script_path')}",
                    )
                    continue
                try:
                    if path not in source_cache:
                        source_cache[path] = self.recover_dim_order_and_tuples(path)
                    recovered_order, tuples = source_cache[path]
                    declared_order = symbol.get("dim_source_order")
                    if recovered_order != declared_order:
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_ORDER_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}: recovered {recovered_order}, declared {declared_order}",
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
                        expression_key = "expr:" + expression_source
                        if expression_key not in tuples:
                            _, expression_tuples = self.recover_dim_order_and_tuples(
                                path, (expression_source,)
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
                    transposed = {axis: value for axis, value in zip(recovered_order, declared_tuple)}
                    named = tuple(transposed[axis] for axis in "LMT")
                    if named != normalized_dim(symbol.get("dim")):
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_NAMED_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}: {dim_json(named)} != {symbol.get('dim')}",
                        )
                except UnsupportedExpression as exc:
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
        "| Check | Outcome |",
        "|---|---|",
    ]
    lines.extend(f"| {name} | {report.results[name].status} |" for name in CHECK_NAMES)
    for name in CHECK_NAMES:
        result = report.results[name]
        lines.extend(["", f"## {name}", ""])
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
    source = root / "dim_source.py"
    source.write_text(
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
    source_digest = hash_file(source)
    register = root / "parameter_register.txt"
    register.write_text("row_a\nrow_b\n")
    register_digest = hash_file(register)
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
        "source": str(source),
        "source_digest": source_digest,
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

        def execute(
            manifests: list[dict[str, Any]],
            config: dict[str, Any],
            *,
            closed_slice: bool | None = True,
        ) -> CompositeReport:
            return CompositeChecker(
                manifests,
                config,
                schema=schema,
                roots=(root, PROJECT_ROOT, PROJECT_ROOT.parent.parent),
                closed_slice=closed_slice,
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

        def case(
            label: str,
            target: str,
            planted_status: str,
            code: str,
            mutate: Any,
        ) -> None:
            nonlocal passed, total
            total += 1
            manifests = copy.deepcopy(baseline_manifests)
            config = copy.deepcopy(baseline_config)
            mutate(manifests, config)
            report = execute(manifests, config)
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
            ok = statuses[target] == planted_status and code in issue_codes and not non_target
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
            "C4 order-lie LTM source declared LMT",
            "C4",
            "FAIL",
            "DIM_SOURCE_ORDER_MISMATCH",
            lambda m, c: next(
                s
                for s in next(x for x in m if x["stage_id"] == "stage001")["symbols"]
                if s["parse_alias"] == "a_B"
            ).__setitem__("dim_source_order", "LMT"),
        )

        case(
            "C4 function-scoped composite live Dim tuple",
            "C4",
            "FAIL",
            "DIM_SOURCE_TUPLE_MISMATCH",
            lambda m, c: next(
                s
                for s in next(x for x in m if x["stage_id"] == "stage001")["symbols"]
                if s["parse_alias"] == "K_m"
            ).__setitem__("dim_source_tuple", ["4", "1", "-2"]),
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
            "C4 unsupported live Dim operation raises",
            "C4",
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
            "C4 container label is not a Dim binding",
            "C4",
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
            "C4 decoy literal cannot impersonate locus assignment",
            "C4",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
            decoy_assignment,
        )

        case(
            "C4 failed locus expression cannot fall through to alias",
            "C4",
            "UNSUPPORTED",
            "DIM_SOURCE_UNSUPPORTED",
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
                "[PASS] C4 legitimate composite expression without named assignment | "
                "target=C4 | clean=PASS planted=PASS"
            )
        else:
            failures.append(
                "C4 legitimate composite expression: "
                f"headline={composite_report.headline}, "
                f"statuses={ {name: result.status for name, result in composite_report.results.items()} }, "
                f"issues={[(issue.code, issue.detail) for issue in composite_report.results['C4'].issues]}"
            )
            lines.append(
                "[FAIL] C4 legitimate composite expression without named assignment"
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
            "C4 dimensional break",
            "C4",
            "FAIL",
            "DIMENSIONAL_INHOMOGENEITY",
            lambda m, c: _find_claim(m, "stage002", "define_y")["payload"]["rhs"].__setitem__(
                "sympy", "x*ell"
            ),
        )

        total += 1
        rational = execute(copy.deepcopy(baseline_manifests), copy.deepcopy(baseline_config))
        if rational.results["C4"].status == "PASS":
            passed += 1
            lines.append("[PASS] C4 rational-power sqrt(area)=ell | target=C4 | clean=PASS planted=PASS")
        else:
            failures.append("C4 rational-power fixture did not pass")
            lines.append("[FAIL] C4 rational-power sqrt(area)=ell")

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
                "path": str(register),
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
                    Path(paths["source"]),
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
            issue.code for issue in open_internal.results["C4"].issues
        }
        if (
            auto_clean.headline == "PASS"
            and open_internal.headline == "FAIL"
            and open_internal.results["C4"].status == "FAIL"
            and "DIMENSIONAL_INHOMOGENEITY" in open_internal_codes
            and all(
                result.status == "PASS"
                for name, result in open_internal.results.items()
                if name != "C4"
            )
        ):
            passed += 1
            lines.append(
                "[PASS] open internal dimensional defect | target=C4 | "
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
        "--dim-order",
        choices=("LMT", "LTM", "MLT", "MTL", "TLM", "TML"),
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
            payload = _live_dim_recovery(
                args.recover_dims.resolve(),
                args.dim_order,
                args.dim_expression,
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
