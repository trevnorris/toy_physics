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
        self.closed_slice = config.get("closed_slice", False) if closed_slice is None else closed_slice
        self.results = {name: CheckResult(name) for name in CHECK_NAMES}
        self.coverage = Coverage()
        self.stages: dict[str, dict[str, Any]] = {}
        self.claims: dict[str, dict[str, Any]] = {}
        self.exports: dict[str, dict[str, Any]] = {}
        self.symbols_by_qid: dict[str, list[tuple[str, dict[str, Any]]]] = defaultdict(list)
        self.symbols_by_stage_alias: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
        self.owners: dict[str, tuple[str, str, str]] = {}
        self.parse_cache: dict[tuple[str, str], sp.Expr] = {}

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
        self.coverage.total_claims = sum(len(m.get("claims", [])) for m in self.manifests)
        self.coverage.checked_claims = self.coverage.total_claims
        return CompositeReport(self.results, self.coverage)

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
                    status = "FAIL" if self.closed_slice else "PARTIAL"
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

    def recover_dim_order_and_tuples(self, path: Path) -> tuple[str, dict[str, tuple[Fraction, ...]]]:
        try:
            tree = ast.parse(path.read_text())
        except Exception as exc:
            raise UnsupportedExpression(f"cannot parse dim source {path}: {exc}") from exc
        field_order: str | None = None
        doc_order: str | None = None
        tuples: dict[str, tuple[Fraction, ...]] = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == "Dim":
                fields = [
                    child.target.id.upper()
                    for child in node.body
                    if isinstance(child, ast.AnnAssign)
                    and isinstance(child.target, ast.Name)
                    and child.target.id.lower() in {"l", "m", "t"}
                ]
                if len(fields) >= 3:
                    field_order = "".join(fields[:3])
                for child in node.body:
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == "__init__":
                        args = [arg.arg.upper() for arg in child.args.args if arg.arg != "self"]
                        args = [arg for arg in args if arg in {"L", "M", "T"}]
                        if len(args) >= 3:
                            field_order = "".join(args[:3])
                doc = ast.get_docstring(node) or ""
                match = re.search(r"\[\s*([LMT])\s*,\s*([LMT])\s*,\s*([LMT])\s*\]", doc, re.I)
                if match:
                    doc_order = "".join(group.upper() for group in match.groups())
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                value = node.value
                if (
                    isinstance(value, ast.Call)
                    and isinstance(value.func, ast.Name)
                    and value.func.id == "Dim"
                    and len(value.args) == 3
                ):
                    try:
                        raw = tuple(Fraction(ast.literal_eval(arg)) for arg in value.args)
                    except Exception:
                        continue
                    for target in targets:
                        if isinstance(target, ast.Name):
                            tuples[target.id] = raw
        if field_order and doc_order and field_order != doc_order:
            raise UnsupportedExpression(f"Dim field/doc order conflict: {field_order}!={doc_order}")
        order = field_order or doc_order
        if order not in {"LMT", "LTM", "MLT", "MTL", "TLM", "TML"}:
            raise UnsupportedExpression("Dim order not recoverable from source structure")
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
                    recovered_order, tuples = source_cache.setdefault(
                        path, self.recover_dim_order_and_tuples(path)
                    )
                    declared_order = symbol.get("dim_source_order")
                    if recovered_order != declared_order:
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_ORDER_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}: recovered {recovered_order}, declared {declared_order}",
                        )
                    declared_tuple = tuple(Fraction(value) for value in symbol.get("dim_source_tuple", []))
                    locus = source.get("locus", "")
                    candidates = [
                        symbol.get("parse_alias"),
                        symbol.get("name"),
                        *re.findall(r"[A-Za-z_][A-Za-z0-9_]*", locus),
                    ]
                    target = next((name for name in candidates if name in tuples), None)
                    if not target:
                        raise UnsupportedExpression(f"raw Dim tuple not found at locus {locus}")
                    if tuples[target] != declared_tuple:
                        result.add(
                            "FAIL",
                            "DIM_SOURCE_TUPLE_MISMATCH",
                            f"{stage_id}:{symbol.get('parse_alias')}: source {tuples[target]}, declared {declared_tuple}",
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

    def registry_rows(self, path: Path, format_name: str) -> set[str]:
        lines = path.read_text().splitlines()
        if format_name == "plain":
            return {line.strip() for line in lines if line.strip() and not line.lstrip().startswith("#")}
        if format_name == "markdown_master_table":
            rows: set[str] = set()
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
                    first = line.split("|", 2)[1]
                    tokens = re.findall(r"`([^`]+)`", first)
                    if tokens:
                        rows.add(tokens[0])
            return rows
        raise ValueError(f"unknown parameter register format {format_name}")

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
            rows = self.registry_rows(path, register.get("format", "plain"))
        except ValueError as exc:
            result.add("UNSUPPORTED", "REGISTER_FORMAT", str(exc))
            return
        events: list[tuple[int, str, dict[str, Any]]] = []
        categories: dict[str, set[str]] = defaultdict(set)
        for stage_id, manifest in self.stages.items():
            for knob in manifest.get("knobs", []):
                events.append((stage_number(stage_id), stage_id, knob))
                categories[knob.get("registry_row")].add(knob.get("count_category"))
        for row in sorted(rows - set(categories)):
            result.add("FAIL", "UNCLASSIFIED_REGISTER_ROW", row)
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
            status = "FAIL" if self.closed_slice else "PARTIAL"
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
        "x = Dim(0, 0, 0)\n"
        "x3 = Dim(1, 0, 0)\n"
        "a_B = Dim(-2, -2, 1)\n"
        "ell = Dim(1, 0, 0)\n"
        "area = Dim(2, 0, 0)\n"
        "f0 = Dim(-1, 0, 0)\n"
        "y = Dim(0, 0, 0)\n"
        "z = Dim(0, 0, 0)\n"
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

        def execute(manifests: list[dict[str, Any]], config: dict[str, Any]) -> CompositeReport:
            return CompositeChecker(
                manifests,
                config,
                schema=schema,
                roots=(root, PROJECT_ROOT, PROJECT_ROOT.parent.parent),
                closed_slice=True,
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
    parser.add_argument("--closed-slice", action="store_true")
    args = parser.parse_args(argv)
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
    report = CompositeChecker(
        manifests,
        config,
        roots=(PROJECT_ROOT, PROJECT_ROOT.parent.parent, Path.cwd()),
        closed_slice=args.closed_slice or config.get("closed_slice", False),
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
