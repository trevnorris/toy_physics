#!/usr/bin/env python3
"""Validate dimension-survey reports and verifier outputs.

This validator is intentionally static.  It parses YAML and referenced report
YAML, but never imports or executes an audit script.
"""

from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import posixpath
import re
import sys
import textwrap
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import yaml
from jsonschema import Draft202012Validator


REPORT_VERSION = "dimension_survey_report_v1"
VERIFICATION_VERSION = "dimension_survey_verification_v1"
REPORT_SCHEMA_NAME = "dimension_survey_report_v1.yaml"
VERIFICATION_SCHEMA_NAME = "dimension_survey_verification_v1.yaml"
SCHEMA_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCHEMA_DIR.parents[2]

REPORT_COLLECTION_PATHS = (
    ("idiom",),
    ("static_recovery_surface",),
    ("basis", "axis_evidence"),
    ("basis", "orders"),
    ("quantities",),
    ("pass_lines", "emitter_sites"),
    ("pass_lines", "call_sites"),
    ("pass_lines", "statically_resolved_markers"),
    ("pass_lines", "unresolved_dynamic_markers"),
    ("active_mutation", "locus"),
    ("active_mutation", "variable_names"),
    ("active_mutation", "gates"),
    ("wl_evidence", "pairs"),
    ("gaps",),
    ("notes",),
)

VERIFICATION_CATEGORIES = (
    "dimension_bindings",
    "quantities",
    "contract_features_and_presentation_consumers",
    "pass_emitter_sites",
    "pass_call_sites",
)

ITEM_ID_INVENTORIES = (
    "static_recovery_surface",
    "quantities",
    "dimension_bindings",
    "dimension_use_sites",
    "contract_features",
    "presentation_consumers",
    "pass_emitter_sites",
    "pass_call_sites",
)

ALGEBRA_OPERATIONS = {
    "COMPOSE",
    "DECOMPOSE",
    "POWER",
    "EQUALITY",
    "CONSTRUCT",
    "NORMALIZE",
    "COMPONENT_ACCESS",
    "EVALUATE_EXPRESSION",
    "PROJECT",
    "RESIDUAL",
}
EXPLICIT_DIMENSION_STATUSES = {
    "EXPLICIT_PY",
    "EXPLICIT_WL",
    "EXPLICIT_PARAMETER_REGISTER",
    "EXPLICIT_STAGE_NOTE",
}
REGISTER_EVIDENCE_STATUSES = {
    "TEXT_IDENTICAL",
    "TEXT_DIFFERS",
    "ABSENT_FROM_REGISTER",
    "PENDING_IN_REGISTER",
    "NOT_CHECKED",
}
REGISTER_ENTRY_STATUSES = {
    "TEXT_IDENTICAL",
    "TEXT_DIFFERS",
    "PENDING_IN_REGISTER",
}
REGISTER_TEXT_STATUSES = {"TEXT_IDENTICAL", "TEXT_DIFFERS"}
WL_EVIDENCE_STATUSES = {
    "TEXT_IDENTICAL",
    "TEXT_DIFFERS",
    "ABSENT_FROM_WL",
    "ABSENT_FROM_PY",
    "NOT_COMPARABLE_BY_CONSTRUCTION",
    "NOT_CHECKED",
}
REGISTER_HEADER = (
    "Param",
    "[L,T,M] dim",
    "Enters",
    "Class",
    "Depends on / relation",
    "Reduction route + status",
)
SOURCE_TOKEN_PATTERN = re.compile(r"\w+|[^\w\s]", re.UNICODE)
SYMPY_NODE_KIND_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
LOCUS_PATTERN = re.compile(r"^(?P<file>[^\n]+):(?P<start>[1-9][0-9]*)(?:-(?P<end>[1-9][0-9]*))?$")
CHECKER_LOCUS_PATTERN = re.compile(r"(?:^|/)composite_build\.py:[1-9][0-9]*(?:-[1-9][0-9]*)?$")
OWNERSHIP_ANCHOR_STAGES = {"021", "038", "042"}
OWNERSHIP_ANCHOR_SPECS = {
    "021": (
        ("module", "SOURCED_N0_DIM", None, None),
        ("module", "SOURCED_D0_DIM", None, None),
        ("run_mu_free_gate", "raw_dim", None, None),
        ("run_mu_free_gate", "norm_dim", None, None),
        ("run_mu_free_gate", "p0_dim", None, None),
    ),
    "038": (
        ("unit_state", "A_E_dim", "active_mutation", "NEGATIVE_CONTROL"),
        ("unit_state", "q_T_dim", "active_mutation", "NEGATIVE_CONTROL"),
        ("unit_state", "c_gamma2_dim", "active_mutation", "NEGATIVE_CONTROL"),
        ("unit_state", "c_E2_dim", "active_mutation", "NEGATIVE_CONTROL"),
        ("unit_state", "mu_R_dim", "active_mutation", "NEGATIVE_CONTROL"),
        ("unit_state", "ratio_dim", None, "NEGATIVE_CONTROL"),
        ("unit_state", "delta_dim", None, None),
        ("unit_state", "cone_dim", None, None),
    ),
    "042": (
        ("dimension_guard", "B", None, None),
        ("dimension_guard", "C", None, None),
        ("dimension_guard", "K", None, None),
        ("dimension_guard", "cE", "active_mutation", None),
        ("dimension_guard", "scalar_power_dim", None, None),
        ("dimension_guard", "em_power_dim", None, None),
        ("dimension_guard", "ratio_dim", None, None),
    ),
}
SELECTED_REPORT_BRANCHES = {
    "representation_contract": {
        "APPLICABLE": "applicableRepresentationContract",
        "NOT_APPLICABLE": "notApplicableSection",
    },
    "dimension_bindings": {
        "APPLICABLE": "dimensionBindingSectionApplicable",
        "NOT_APPLICABLE": "notApplicableSection",
    },
    "dimension_use_sites": {
        "APPLICABLE": "dimensionUseSiteSectionApplicable",
        "NOT_APPLICABLE": "notApplicableSection",
    },
}
REJECT_FIXTURE_VERSION = "dimension_survey_reject_fixture_v1"
ACCEPT_FIXTURE_VERSION = "dimension_survey_accept_fixture_v1"
REAL_REGISTER_PATH = (SCHEMA_DIR.parent / "notes" / "parameter_register.md").resolve()


@dataclass(frozen=True)
class Problem:
    code: str
    path: str
    message: str

    def render(self) -> str:
        return f"{self.code} {self.path}: {self.message}"


@dataclass(frozen=True)
class RegisterRow:
    line: int
    param: str
    dimension_cell: str
    cells: tuple[str, ...]

    @property
    def pending_dimension(self) -> bool:
        return re.match(r"^pending(?:\s|\()", self.dimension_cell) is not None


def problem(code: str, path: str, message: str) -> Problem:
    return Problem(code=code, path=path or "$", message=message)


def yaml_load(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def fixture_path(raw: str, fixture_file: Path) -> Path:
    candidate = Path(raw)
    if candidate.is_absolute():
        return candidate
    beside_fixture = fixture_file.parent / candidate
    if beside_fixture.exists():
        return beside_fixture
    return SCHEMA_DIR / candidate


def mutation_parent(document: Any, parts: Sequence[Any]) -> tuple[Any, Any]:
    if not parts:
        raise ValueError("mutation path must not be empty")
    current = document
    for part in parts[:-1]:
        if isinstance(part, int):
            if not isinstance(current, list) or not (0 <= part < len(current)):
                raise ValueError(f"mutation index {part!r} does not resolve")
            current = current[part]
        else:
            if not isinstance(current, Mapping) or part not in current:
                raise ValueError(f"mutation key {part!r} does not resolve")
            current = current[part]
    return current, parts[-1]


def apply_fixture_mutations(document: Any, mutations: Any) -> Any:
    result = copy.deepcopy(document)
    if not isinstance(mutations, list) or not mutations:
        raise ValueError("reject fixture requires a non-empty mutations list")
    for mutation in mutations:
        if not isinstance(mutation, Mapping) or not isinstance(mutation.get("path"), list):
            raise ValueError("each mutation requires a list path")
        operation = mutation.get("op")
        parent, leaf = mutation_parent(result, mutation["path"])
        if operation in {"replace", "add"}:
            if "value" not in mutation:
                raise ValueError(f"{operation} mutation requires value")
            value = copy.deepcopy(mutation["value"])
            if isinstance(leaf, int):
                if not isinstance(parent, list):
                    raise ValueError("numeric mutation leaf requires a list parent")
                if operation == "add":
                    if not (0 <= leaf <= len(parent)):
                        raise ValueError(f"add index {leaf!r} does not resolve")
                    parent.insert(leaf, value)
                else:
                    if not (0 <= leaf < len(parent)):
                        raise ValueError(f"replace index {leaf!r} does not resolve")
                    parent[leaf] = value
            else:
                if not isinstance(parent, dict):
                    raise ValueError("string mutation leaf requires a mapping parent")
                if operation == "replace" and leaf not in parent:
                    raise ValueError(f"replace key {leaf!r} does not resolve")
                parent[leaf] = value
        elif operation == "remove":
            if isinstance(leaf, int):
                if not isinstance(parent, list) or not (0 <= leaf < len(parent)):
                    raise ValueError(f"remove index {leaf!r} does not resolve")
                parent.pop(leaf)
            else:
                if not isinstance(parent, dict) or leaf not in parent:
                    raise ValueError(f"remove key {leaf!r} does not resolve")
                del parent[leaf]
        else:
            raise ValueError(f"unsupported mutation operation {operation!r}")
    return result


def load_effective_document(path: Path) -> Any:
    document = yaml_load(path)
    if (
        not isinstance(document, Mapping)
        or document.get("schema_version") not in {REJECT_FIXTURE_VERSION, ACCEPT_FIXTURE_VERSION}
    ):
        return document
    base_raw = document.get("base_document")
    if not isinstance(base_raw, str):
        raise ValueError("reject fixture requires base_document")
    base_path = fixture_path(base_raw, path)
    base = yaml_load(base_path)
    if (
        isinstance(base, Mapping)
        and base.get("schema_version") in {REJECT_FIXTURE_VERSION, ACCEPT_FIXTURE_VERSION}
    ):
        raise ValueError("example fixtures may not extend another example fixture")
    return apply_fixture_mutations(base, document.get("mutations"))


def json_path(parts: Iterable[Any]) -> str:
    rendered = "$"
    for part in parts:
        if isinstance(part, int):
            rendered += f"[{part}]"
        else:
            rendered += f".{part}"
    return rendered


def get_path(node: Any, parts: Sequence[str]) -> Any:
    current = node
    for part in parts:
        if not isinstance(current, Mapping) or part not in current:
            return None
        current = current[part]
    return current


def fact_value(node: Any) -> Any:
    return node.get("value") if isinstance(node, Mapping) else None


def items(node: Any) -> list[Any]:
    if isinstance(node, Mapping) and isinstance(node.get("items"), list):
        return node["items"]
    return []


def item_id_set(nodes: Iterable[Mapping[str, Any]]) -> set[str]:
    return {
        node["item_id"]
        for node in nodes
        if isinstance(node, Mapping) and isinstance(node.get("item_id"), str)
    }


def classification_for_report(report: Mapping[str, Any]) -> str | None:
    content = get_path(report, ("dimensional_content", "value"))
    machinery = fact_value(report.get("machinery"))
    if content == "NONE" and machinery == "ABSENT":
        return "NONE+ABSENT"
    if content in {"REAL", "UNDETERMINED"} and machinery in {"PRESENT", "ABSENT"}:
        return f"{content}+{machinery}"
    return None


def schema_problems(schema: Mapping[str, Any], document: Any) -> list[Problem]:
    if (
        isinstance(document, Mapping)
        and document.get("schema_version") == REPORT_VERSION
        and isinstance(schema.get("$defs"), Mapping)
    ):
        selected: list[Problem] = []
        for section_name, branches in SELECTED_REPORT_BRANCHES.items():
            section = document.get(section_name)
            if not isinstance(section, Mapping):
                continue
            applicability = section.get("applicability")
            definition = branches.get(applicability)
            if definition is None:
                selected.append(
                    problem(
                        "SCHEMA_STRUCTURE",
                        f"$.{section_name}.applicability",
                        "expected APPLICABLE or NOT_APPLICABLE; "
                        f"got {applicability!r}",
                    )
                )
                continue
            branch_schema = {
                "$schema": schema.get("$schema"),
                "$defs": schema["$defs"],
                "$ref": f"#/$defs/{definition}",
            }
            branch_validator = Draft202012Validator(branch_schema)
            for error in sorted(
                branch_validator.iter_errors(section),
                key=lambda entry: tuple(entry.absolute_path),
            ):
                if (
                    section_name == "representation_contract"
                    and len(error.absolute_path) >= 2
                    and error.absolute_path[0] == "expression_evaluator"
                    and error.absolute_path[1]
                    in {
                        "numeric_zero_policy",
                        "missing_symbol_policy",
                        "non_homogeneous_add_policy",
                    }
                ):
                    continue
                suffix = json_path(error.absolute_path)
                path = f"$.{section_name}" + (suffix[1:] if suffix != "$" else "")
                selected.append(
                    problem(
                        "SCHEMA_STRUCTURE",
                        path,
                        f"selected {applicability} branch; {error.message}",
                    )
                )
        contract = document.get("representation_contract")
        evaluator = (
            contract.get("expression_evaluator")
            if isinstance(contract, Mapping) and contract.get("applicability") == "APPLICABLE"
            else None
        )
        if isinstance(evaluator, Mapping):
            for policy_name in (
                "numeric_zero_policy",
                "missing_symbol_policy",
                "non_homogeneous_add_policy",
            ):
                policy = evaluator.get(policy_name)
                if not isinstance(policy, Mapping):
                    continue
                definition = (
                    "notApplicableSection"
                    if policy.get("applicability") == "NOT_APPLICABLE"
                    else "contractFeature"
                )
                branch_schema = {
                    "$schema": schema.get("$schema"),
                    "$defs": schema["$defs"],
                    "$ref": f"#/$defs/{definition}",
                }
                branch_validator = Draft202012Validator(branch_schema)
                policy_errors = sorted(
                    branch_validator.iter_errors(policy),
                    key=lambda entry: tuple(entry.absolute_path),
                )
                if policy_errors:
                    path = "$.representation_contract.expression_evaluator." + policy_name
                    messages = "; ".join(error.message for error in policy_errors)
                    selected.append(
                        problem(
                            "SCHEMA_STRUCTURE",
                            path,
                            "selected "
                            + (
                                "NOT_APPLICABLE"
                                if definition == "notApplicableSection"
                                else "policy"
                            )
                            + f" branch; {messages}",
                        )
                    )
        if selected:
            return selected

    validator = Draft202012Validator(schema)
    result: list[Problem] = []
    for error in sorted(validator.iter_errors(document), key=lambda e: tuple(e.absolute_path)):
        result.append(
            problem(
                "SCHEMA_STRUCTURE",
                json_path(error.absolute_path),
                error.message,
            )
        )
    return result


def parsed_locus(raw: str) -> tuple[str, int, int] | None:
    match = LOCUS_PATTERN.fullmatch(raw)
    if match is None:
        return None
    file_name = posixpath.normpath(match.group("file"))
    start = int(match.group("start"))
    end = int(match.group("end") or start)
    return file_name, start, end


def overlapping_loci(left: Iterable[str], right: Iterable[str]) -> list[tuple[str, str]]:
    overlaps: list[tuple[str, str]] = []
    parsed_right = [(raw, parsed_locus(raw)) for raw in right]
    for left_raw in left:
        left_parts = parsed_locus(left_raw)
        if left_parts is None:
            continue
        left_file, left_start, left_end = left_parts
        for right_raw, right_parts in parsed_right:
            if right_parts is None:
                continue
            right_file, right_start, right_end = right_parts
            if (
                left_file == right_file
                and left_start <= right_end
                and right_start <= left_end
            ):
                overlaps.append((left_raw, right_raw))
    return sorted(set(overlaps))


def directive_anchor_refs(
    stage: str,
    bindings: Sequence[Mapping[str, Any]],
) -> tuple[set[str], list[str]]:
    expected: set[str] = set()
    errors: list[str] = []
    for enclosing_scope, name, branch_state, required_use in OWNERSHIP_ANCHOR_SPECS.get(stage, ()):
        candidates: list[str] = []
        for binding in bindings:
            if fact_value(binding.get("name")) != name:
                continue
            if binding_enclosing_scope(binding) != enclosing_scope:
                continue
            if not anchor_branch_state_matches(binding, branch_state):
                continue
            uses = {fact_value(use) for use in items(binding.get("uses"))}
            if required_use is not None and required_use not in uses:
                continue
            item_id = binding.get("item_id")
            if isinstance(item_id, str):
                candidates.append(item_id)
        if len(candidates) != 1:
            errors.append(
                "directive lexical anchor "
                f"{stage}:{enclosing_scope}:{name}:{branch_state or 'unbranched'} "
                f"resolved to {len(candidates)} report bindings"
            )
        else:
            expected.add(candidates[0])
    return expected, errors


def binding_enclosing_scope(binding: Mapping[str, Any]) -> str | None:
    raw_scope = fact_value(binding.get("qualified_scope"))
    name = fact_value(binding.get("name"))
    if not isinstance(raw_scope, str) or not isinstance(name, str):
        return None
    scope = raw_scope.strip()
    escaped_name = re.escape(name)
    for pattern in (
        rf"\.{escaped_name}(?:\[[^\]]+\])?$",
        rf"\[['\"]?{escaped_name}['\"]?\]$",
    ):
        reduced = re.sub(pattern, "", scope)
        if reduced != scope:
            scope = reduced
            break
    # r13 identifies B/C/K as bindings in dimension_guard's ``dims`` map and
    # cE as a sibling map key.  The container spelling is not a lexical scope:
    # accept both source-like ``dimension_guard.dims["B"]`` and the directive's
    # ``dimension_guard.B`` convention as the same enclosing function.
    if (
        scope == "dimension_guard.dims"
        and name in {"B", "C", "K", "cE"}
    ):
        return "dimension_guard"
    return scope


def anchor_branch_state_matches(binding: Mapping[str, Any], expected_state: str | None) -> bool:
    if expected_state is None:
        return True
    if expected_state == "active_mutation":
        return fact_value(binding.get("depends_on_active_mutation")) is True
    return False


def unknown_leaf_shape_problems(node: Any, path: str = "$") -> list[Problem]:
    result: list[Problem] = []
    if isinstance(node, Mapping):
        if node.get("value") == "UNDETERMINED" and path != "$.dimensional_content":
            reason = node.get("reason")
            searched = node.get("sources_searched")
            if (
                not isinstance(reason, str)
                or not reason.strip()
                or not isinstance(searched, list)
                or not searched
                or "loci" in node
            ):
                result.append(
                    problem(
                        "REPORT_SMALLEST_LEAF_UNKNOWN",
                        path,
                        "UNDETERMINED requires reason and non-empty sources_searched "
                        "at this leaf, and cannot carry known-fact loci",
                    )
                )
        for key, value in node.items():
            result.extend(unknown_leaf_shape_problems(value, f"{path}.{key}"))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            result.extend(unknown_leaf_shape_problems(value, f"{path}[{index}]"))
    return result


def item_id_scope_problems(report: Mapping[str, Any]) -> list[Problem]:
    allowed_nodes = {
        id(node)
        for inventory in report_item_inventories(report).values()
        for node in inventory
    }
    result: list[Problem] = []

    def walk(node: Any, path: str) -> None:
        if isinstance(node, Mapping):
            if "item_id" in node and id(node) not in allowed_nodes:
                result.append(
                    problem(
                        "REPORT_ID_SCOPE",
                        f"{path}.item_id",
                        "item_id is legal only on a §4.0 inventory item or "
                        "wl_evidence comparison record",
                    )
                )
            for key, value in node.items():
                walk(value, f"{path}.{key}")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                walk(value, f"{path}[{index}]")

    walk(report, "$")
    return result


def report_preflight(report: Any, document_path: Path | None = None) -> list[Problem]:
    if not isinstance(report, Mapping):
        return [problem("SCHEMA_STRUCTURE", "$", "document must be a mapping")]
    for parts in REPORT_COLLECTION_PATHS:
        node = get_path(report, parts)
        if node == "UNDETERMINED":
            return [
                problem(
                    "REPORT_SMALLEST_LEAF_UNKNOWN",
                    json_path(parts),
                    "a collection cannot be replaced wholesale by scalar UNDETERMINED",
                )
            ]

    result: list[Problem] = []
    result.extend(unknown_leaf_shape_problems(report))
    selected_applicabilities_valid = all(
        isinstance(report.get(section_name), Mapping)
        and report[section_name].get("applicability") in branches
        for section_name, branches in SELECTED_REPORT_BRANCHES.items()
    )
    if selected_applicabilities_valid:
        result.extend(item_id_scope_problems(report))
    idiom_nodes = items(report.get("idiom"))
    idiom_seen: dict[Any, int] = {}
    idiom_loci: dict[Any, set[str]] = {}
    for index, node in enumerate(idiom_nodes):
        value = fact_value(node)
        if value in idiom_seen:
            result.append(
                problem(
                    "REPORT_IDIOM_VALUE_UNIQUENESS",
                    f"$.idiom.items[{index}].value",
                    f"{value!r} already appears at $.idiom.items[{idiom_seen[value]}].value",
                )
            )
        elif value is not None:
            idiom_seen[value] = index
        loci = node.get("loci") if isinstance(node, Mapping) else None
        if isinstance(loci, list):
            idiom_loci.setdefault(value, set()).update(
                locus for locus in loci if isinstance(locus, str)
            )
    overlap = overlapping_loci(
        idiom_loci.get("TYPE_ALIAS_TUPLE", set()),
        idiom_loci.get("BARE_TUPLE", set()),
    )
    if overlap:
        result.append(
            problem(
                "REPORT_IDIOM_MEMBERSHIP",
                "$.idiom",
                "TYPE_ALIAS_TUPLE and BARE_TUPLE have overlapping representation "
                f"line ranges: {overlap!r}",
            )
        )

    for section_name in ("dimension_bindings", "dimension_use_sites"):
        section = report.get(section_name)
        for item_index, node in enumerate(items(section)):
            uses = items(node.get("uses")) if isinstance(node, Mapping) else []
            seen_uses: dict[Any, int] = {}
            for use_index, use in enumerate(uses):
                value = fact_value(use)
                if value in seen_uses:
                    result.append(
                        problem(
                            "REPORT_USES_VALUE_UNIQUENESS",
                            f"$.{section_name}.items[{item_index}].uses.items[{use_index}].value",
                            f"{value!r} already appears in this uses set at index "
                            f"{seen_uses[value]}",
                        )
                    )
                elif value is not None:
                    seen_uses[value] = use_index

    classification = classification_for_report(report)
    for index, quantity in enumerate(items(report.get("quantities"))):
        if not isinstance(quantity, Mapping):
            continue
        if classification == "REAL+ABSENT":
            candidate_value = fact_value(quantity.get("candidate_reason"))
            if (
                not isinstance(candidate_value, str)
                or not candidate_value.strip()
                or candidate_value in {"NOT_APPLICABLE", "UNDETERMINED"}
            ):
                result.append(
                    problem(
                        "REPORT_CANDIDATE_REASON",
                        f"$.quantities.items[{index}].candidate_reason",
                        "REAL+ABSENT requires a substantive candidate reason; "
                        f"got {candidate_value!r}",
                    )
                )
        status = fact_value(quantity.get("dimension_value_status"))
        source_text = fact_value(quantity.get("dim_source_text"))
        source_loci = [
            locus
            for locus in items(quantity.get("dimension_source_loci"))
            if isinstance(locus, str)
        ]
        if (
            status in EXPLICIT_DIMENSION_STATUSES
            and (
                not isinstance(source_text, str)
                or not source_tokens(source_text)
                or not explicit_dimension_text_matches_declared_basis(
                    quantity,
                    report,
                )
                or document_path is None
                or not any(
                    source_text_occurs_at_locus(source_text, locus, document_path)
                    for locus in source_loci
                )
            )
        ):
            result.append(
                problem(
                    "REPORT_EXPLICIT_DIM_SOURCE_TEXT",
                    f"$.quantities.items[{index}].dim_source_text",
                    f"{status} requires source-present dimension text with a valid "
                    "dim_text_form and in-report basis_locus; NAMED_AXIS uses a "
                    "standalone declared spelling (or exactly 1), while POSITIONAL "
                    "requires one component per declared axis",
                )
            )
        expected_source_field = {
            "EXPLICIT_PY": "script_path",
            "EXPLICIT_WL": "wl_path",
        }.get(status)
        if expected_source_field is not None and (
            document_path is None
            or not any(
                locus_matches_identity_source(
                    locus,
                    get_path(report.get("identity", {}), (expected_source_field, "value")),
                    document_path,
                )
                for locus in source_loci
            )
        ):
            result.append(
                problem(
                    "REPORT_EXPLICIT_DIM_SOURCE_TEXT",
                    f"$.quantities.items[{index}].dimension_source_loci",
                    f"{status} requires a locus in identity.{expected_source_field}, "
                    "not a different source kind",
                )
            )
        if status == "EXPLICIT_STAGE_NOTE" and not any(
            is_stage_note_locus(locus) for locus in source_loci
        ):
            result.append(
                problem(
                    "REPORT_EXPLICIT_DIM_SOURCE_TEXT",
                    f"$.quantities.items[{index}].dimension_source_loci",
                    "EXPLICIT_STAGE_NOTE requires a notes/stages/*.md locus",
                )
            )
        evidence = quantity.get("register_evidence")
        evidence_status = evidence.get("status") if isinstance(evidence, Mapping) else None
        if evidence_status not in REGISTER_EVIDENCE_STATUSES:
            result.append(
                problem(
                    "REPORT_EVIDENCE_RECORDS",
                    f"$.quantities.items[{index}].register_evidence",
                    "every quantity requires a register evidence status",
                )
            )
        elif evidence_status in REGISTER_ENTRY_STATUSES:
            register_locus = evidence.get("register_locus")
            row = (
                master_register_row(register_locus, document_path)
                if isinstance(register_locus, str) and document_path is not None
                else None
            )
            if row is None:
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"$.quantities.items[{index}].register_evidence.register_locus",
                        "register_locus must cite a dimension row in the master "
                        "parameter table, never an R## edge row",
                    )
                )
            elif not isinstance(evidence.get("register_cell_verbatim"), str) or (
                whitespace_normalized(evidence["register_cell_verbatim"])
                != whitespace_normalized(row.dimension_cell)
            ):
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"$.quantities.items[{index}].register_evidence."
                        "register_cell_verbatim",
                        "must preserve the whole master-table dimension cell "
                        "verbatim modulo whitespace",
                    )
                )
            elif evidence_status == "PENDING_IN_REGISTER" and not row.pending_dimension:
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"$.quantities.items[{index}].register_evidence.register_locus",
                        "PENDING_IN_REGISTER requires lowercase pending in the "
                        "master-table dimension cell; uppercase PENDING is a route flag",
                    )
                )
            elif evidence_status in REGISTER_TEXT_STATUSES:
                identical = (
                    isinstance(source_text, str)
                    and whitespace_normalized(source_text)
                    == whitespace_normalized(evidence["register_cell_verbatim"])
                )
                expected_status = "TEXT_IDENTICAL" if identical else "TEXT_DIFFERS"
                if evidence_status != expected_status:
                    result.append(
                        problem(
                            "REPORT_EVIDENCE_RECORDS",
                            f"$.quantities.items[{index}].register_evidence.status",
                            f"string comparison modulo whitespace requires {expected_status}",
                        )
                    )
        elif evidence_status == "NOT_CHECKED":
            reason = evidence.get("reason") if isinstance(evidence, Mapping) else None
            if not isinstance(reason, str) or not reason.strip():
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"$.quantities.items[{index}].register_evidence.reason",
                        "NOT_CHECKED requires a substantive reason",
                    )
                )

        if status == "EXPLICIT_PARAMETER_REGISTER":
            register_rows = [
                master_register_row(locus, document_path)
                for locus in items(quantity.get("dimension_source_loci"))
                if isinstance(locus, str) and document_path is not None
            ]
            if not any(row is not None for row in register_rows):
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"$.quantities.items[{index}].dimension_source_loci",
                        "EXPLICIT_PARAMETER_REGISTER must cite a master "
                        "parameter-table dimension row",
                    )
                )
            if evidence_status == "ABSENT_FROM_REGISTER":
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"$.quantities.items[{index}].register_evidence.status",
                        "EXPLICIT_PARAMETER_REGISTER cannot claim "
                        "ABSENT_FROM_REGISTER",
                    )
                )

        if (
            status in {"CONSTRAINED_BUT_NOT_STATED", "UNDETERMINED"}
            and evidence_status
            not in {
                "TEXT_IDENTICAL",
                "TEXT_DIFFERS",
                "PENDING_IN_REGISTER",
                "ABSENT_FROM_REGISTER",
            }
        ):
            result.append(
                problem(
                    "REPORT_EVIDENCE_RECORDS",
                    f"$.quantities.items[{index}].register_evidence",
                    f"{status} requires TEXT_IDENTICAL, TEXT_DIFFERS, "
                    "PENDING_IN_REGISTER, or ABSENT_FROM_REGISTER before it "
                    "can block readiness; NOT_CHECKED does not discharge consultation",
                )
            )

    quantity_ids = {
        quantity.get("item_id")
        for quantity in items(report.get("quantities"))
        if isinstance(quantity, Mapping) and isinstance(quantity.get("item_id"), str)
    }
    wl_evidence = report.get("wl_evidence")
    wl_pairs = items(wl_evidence.get("pairs") if isinstance(wl_evidence, Mapping) else None)
    states_dimensions = (
        fact_value(wl_evidence.get("states_dimensions"))
        if isinstance(wl_evidence, Mapping)
        else None
    )
    if states_dimensions is False and wl_pairs:
        result.append(
            problem(
                "REPORT_EVIDENCE_RECORDS",
                "$.wl_evidence.pairs",
                "states_dimensions false requires an empty reason-bearing pair list",
            )
        )
    for index, pair in enumerate(wl_pairs):
        if not isinstance(pair, Mapping):
            continue
        path = f"$.wl_evidence.pairs.items[{index}]"
        status = pair.get("status")
        quantity_ref = pair.get("quantity_ref")
        if status not in WL_EVIDENCE_STATUSES:
            result.append(
                problem(
                    "REPORT_EVIDENCE_RECORDS",
                    f"{path}.status",
                    "unknown wl_evidence status",
                )
            )
            continue
        if status == "ABSENT_FROM_PY":
            if quantity_ref != "NO_SURVEYED_QUANTITY":
                result.append(
                    problem(
                        "REPORT_BINDING_REFERENCE",
                        f"{path}.quantity_ref",
                        "ABSENT_FROM_PY requires NO_SURVEYED_QUANTITY",
                    )
                )
        elif quantity_ref not in quantity_ids:
            result.append(
                problem(
                    "REPORT_BINDING_REFERENCE",
                    f"{path}.quantity_ref",
                    "quantity_ref is not a quantities item_id in this report",
                )
            )
        for side, identity_field in (("py", "script_path"), ("wl", "wl_path")):
            text = pair.get(f"{side}_text_verbatim")
            locus = pair.get(f"{side}_locus")
            if isinstance(text, str) and isinstance(locus, str):
                if (
                    document_path is None
                    or not source_text_occurs_at_locus(text, locus, document_path)
                    or not locus_matches_identity_source(
                        locus,
                        get_path(report.get("identity", {}), (identity_field, "value")),
                        document_path,
                    )
                ):
                    result.append(
                        problem(
                            "REPORT_EVIDENCE_RECORDS",
                            f"{path}.{side}_text_verbatim",
                            f"verbatim {side} text must occur at a locus in "
                            f"identity.{identity_field}",
                        )
                    )
            elif isinstance(locus, str) and (
                document_path is None
                or not locus_matches_identity_source(
                    locus,
                    get_path(report.get("identity", {}), (identity_field, "value")),
                    document_path,
                )
            ):
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"{path}.{side}_locus",
                        f"{side}_locus must resolve within identity.{identity_field}",
                    )
                )
        if status in {"TEXT_IDENTICAL", "TEXT_DIFFERS"}:
            identical = (
                isinstance(pair.get("py_text_verbatim"), str)
                and isinstance(pair.get("wl_text_verbatim"), str)
                and whitespace_normalized(pair["py_text_verbatim"])
                == whitespace_normalized(pair["wl_text_verbatim"])
            )
            expected_status = "TEXT_IDENTICAL" if identical else "TEXT_DIFFERS"
            if status != expected_status:
                result.append(
                    problem(
                        "REPORT_EVIDENCE_RECORDS",
                        f"{path}.status",
                        f"string comparison modulo whitespace requires {expected_status}",
                    )
                )

    contract = report.get("representation_contract")
    if isinstance(contract, Mapping):
        for index, node in enumerate(items(contract.get("algebra"))):
            value = fact_value(node.get("semantic_operation")) if isinstance(node, Mapping) else None
            if value is not None and value != "UNDETERMINED" and value not in ALGEBRA_OPERATIONS:
                result.append(
                    problem(
                        "REPORT_ALGEBRA_OPERATION_VOCABULARY",
                        f"$.representation_contract.algebra.items[{index}].semantic_operation.value",
                        f"{value!r} is not an r17 algebra semantic operation",
                    )
                )
        evaluator = contract.get("expression_evaluator")
        if isinstance(evaluator, Mapping):
            for index, node in enumerate(items(evaluator.get("supported_sympy_node_kinds"))):
                value = fact_value(node.get("node_kind")) if isinstance(node, Mapping) else None
                if (
                    value is not None
                    and value != "UNDETERMINED"
                    and (
                        not isinstance(value, str)
                        or SYMPY_NODE_KIND_PATTERN.fullmatch(value) is None
                    )
                ):
                    result.append(
                        problem(
                            "REPORT_EVALUATOR_NODE_KIND_VOCABULARY",
                            "$.representation_contract.expression_evaluator."
                            f"supported_sympy_node_kinds.items[{index}].node_kind.value",
                            f"{value!r} is not an unqualified verbatim SymPy class name",
                        )
                    )
    return result


def verification_preflight(verification: Any) -> list[Problem]:
    if not isinstance(verification, Mapping):
        return [problem("SCHEMA_STRUCTURE", "$", "document must be a mapping")]

    def forbidden_name_reference(node: Any, path: str = "$") -> Problem | None:
        if isinstance(node, Mapping):
            for key, value in node.items():
                child_path = f"{path}.{key}"
                if "binding_name_ref" in str(key):
                    return problem(
                        "VER_NO_NAME_MATCHING",
                        child_path,
                        "cross-references use item_id fields, never binding-name references",
                    )
                found = forbidden_name_reference(value, child_path)
                if found is not None:
                    return found
        elif isinstance(node, list):
            for index, value in enumerate(node):
                found = forbidden_name_reference(value, f"{path}[{index}]")
                if found is not None:
                    return found
        return None

    forbidden = forbidden_name_reference(verification)
    if forbidden is not None:
        return [forbidden]
    categories = verification.get("categories")
    if not isinstance(categories, Mapping) or set(categories) != set(VERIFICATION_CATEGORIES):
        actual = sorted(categories) if isinstance(categories, Mapping) else type(categories).__name__
        return [
            problem(
                "VER_CATEGORY_SET",
                "$.categories",
                f"expected five separate categories {list(VERIFICATION_CATEGORIES)!r}; got {actual!r}",
            )
        ]
    classification = verification.get("classification")
    expected = category_applicability(classification)
    for category_name in VERIFICATION_CATEGORIES:
        category = categories.get(category_name)
        if not isinstance(category, Mapping):
            continue
        wanted = expected.get(category_name)
        applicability = category.get("applicability")
        source_count = category.get("source_count")
        report_count = category.get("report_count")
        valid = (
            wanted == "APPLICABLE"
            and applicability == "APPLICABLE"
            and type(source_count) is int
            and source_count >= 0
            and type(report_count) is int
            and report_count >= 0
        ) or (
            wanted == "NOT_APPLICABLE"
            and applicability == "NOT_APPLICABLE"
            and source_count == "N/A"
            and report_count == "N/A"
        )
        if not valid:
            return [
                problem(
                    "VER_CATEGORY_APPLICABILITY",
                    f"$.categories.{category_name}",
                    f"{classification} requires {wanted} with "
                    + ("numeric cardinalities" if wanted == "APPLICABLE" else "N/A cardinalities"),
                )
            ]
    ownership = verification.get("ownership_rederivation")
    if (
        isinstance(ownership, Mapping)
        and ownership.get("applicability") == "APPLICABLE"
        and isinstance(ownership.get("items"), list)
        and not ownership["items"]
    ):
        return [
            problem(
                "VER_OWNERSHIP_ANCHOR_COVERAGE",
                "$.ownership_rederivation.items",
                "an applicable ownership re-derivation set is never empty, "
                "including when an escape is invoked",
            )
        ]
    return []


def category_applicability(classification: Any) -> dict[str, str]:
    present = classification in {"REAL+PRESENT", "UNDETERMINED+PRESENT"}
    quantities_apply = classification != "NONE+ABSENT"
    return {
        "dimension_bindings": "APPLICABLE" if present else "NOT_APPLICABLE",
        "quantities": "APPLICABLE" if quantities_apply else "NOT_APPLICABLE",
        "contract_features_and_presentation_consumers": (
            "APPLICABLE" if present else "NOT_APPLICABLE"
        ),
        "pass_emitter_sites": "APPLICABLE",
        "pass_call_sites": "APPLICABLE",
    }


def empty_collection_problems(node: Any, path: str = "$") -> list[Problem]:
    result: list[Problem] = []
    if isinstance(node, Mapping):
        if isinstance(node.get("items"), list) and not node["items"]:
            reason = node.get("reason")
            if not isinstance(reason, str) or not reason.strip():
                result.append(
                    problem(
                        "REPORT_EMPTY_COLLECTION_REASON",
                        path,
                        "an empty items collection requires a non-empty reason",
                    )
                )
        for key, value in node.items():
            result.extend(empty_collection_problems(value, f"{path}.{key}"))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            result.extend(empty_collection_problems(value, f"{path}[{index}]"))
    return result


def report_item_inventories(report: Mapping[str, Any]) -> dict[str, list[Mapping[str, Any]]]:
    contract = report.get("representation_contract")
    contract_features: list[Mapping[str, Any]] = []
    presentations: list[Mapping[str, Any]] = []
    if isinstance(contract, Mapping) and contract.get("applicability") == "APPLICABLE":
        for key in (
            "normalization_on_construction",
            "equality_semantics",
            "hashability",
            "component_access",
        ):
            value = contract.get(key)
            if isinstance(value, Mapping):
                contract_features.append(value)
        contract_features.extend(
            node for node in items(contract.get("algebra")) if isinstance(node, Mapping)
        )
        evaluator = contract.get("expression_evaluator")
        if isinstance(evaluator, Mapping):
            contract_features.extend(
                node
                for node in items(evaluator.get("supported_sympy_node_kinds"))
                if isinstance(node, Mapping)
            )
            for key in (
                "numeric_zero_policy",
                "missing_symbol_policy",
                "non_homogeneous_add_policy",
            ):
                value = evaluator.get(key)
                if isinstance(value, Mapping) and value.get("applicability") != "NOT_APPLICABLE":
                    contract_features.append(value)
        presentations.extend(
            node for node in items(contract.get("presentation")) if isinstance(node, Mapping)
        )
    bindings = report.get("dimension_bindings")
    use_sites = report.get("dimension_use_sites")
    pass_lines = report.get("pass_lines")
    return {
        "static_recovery_surface": [
            node for node in items(report.get("static_recovery_surface")) if isinstance(node, Mapping)
        ],
        "quantities": [
            node for node in items(report.get("quantities")) if isinstance(node, Mapping)
        ],
        "dimension_bindings": (
            [
                node
                for node in bindings.get("items", [])
                if isinstance(node, Mapping)
            ]
            if isinstance(bindings, Mapping) and bindings.get("applicability") == "APPLICABLE"
            else []
        ),
        "dimension_use_sites": (
            [
                node
                for node in use_sites.get("items", [])
                if isinstance(node, Mapping)
            ]
            if isinstance(use_sites, Mapping) and use_sites.get("applicability") == "APPLICABLE"
            else []
        ),
        "contract_features": contract_features,
        "presentation_consumers": presentations,
        "pass_emitter_sites": [
            node
            for node in items(pass_lines.get("emitter_sites") if isinstance(pass_lines, Mapping) else None)
            if isinstance(node, Mapping)
        ],
        "pass_call_sites": [
            node
            for node in items(pass_lines.get("call_sites") if isinstance(pass_lines, Mapping) else None)
            if isinstance(node, Mapping)
        ],
        "wl_evidence_comparisons": [
            node
            for node in items(get_path(report, ("wl_evidence", "pairs")))
            if isinstance(node, Mapping)
        ],
    }


def report_semantic_problems(report: Mapping[str, Any]) -> list[Problem]:
    result: list[Problem] = []
    result.extend(empty_collection_problems(report))

    identity = report.get("identity", {})
    method_rules = {
        "script_path": "FILESYSTEM",
        "script_sha256": "SHA256SUM",
        "wl_path": "FILESYSTEM",
        "wl_sha256": "SHA256SUM",
        "script_line_count": "WC_L",
    }
    for field, expected in method_rules.items():
        actual = get_path(identity, (field, "method"))
        if actual != expected:
            result.append(
                problem(
                    "REPORT_WHOLE_FILE_METHODS",
                    f"$.identity.{field}.method",
                    f"expected {expected}, got {actual!r}",
                )
            )

    content = get_path(report, ("dimensional_content", "value"))
    machinery = fact_value(report.get("machinery"))
    classification = classification_for_report(report)
    if classification is None:
        result.append(
            problem(
                "REPORT_APPLICABILITY_MATRIX",
                "$.dimensional_content",
                f"invalid content/machinery combination {content!r}+{machinery!r}",
            )
        )
        return result

    idioms = [fact_value(node) for node in items(report.get("idiom"))]
    definition = report.get("definition_locus")
    recovery = report.get("static_recovery_surface")
    quantities = report.get("quantities")
    contract = report.get("representation_contract")
    bindings = report.get("dimension_bindings")
    use_sites = report.get("dimension_use_sites")
    basis = report.get("basis", {})
    wl_evidence = report.get("wl_evidence", {})

    if machinery == "ABSENT":
        if idioms != ["NONE"]:
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.idiom",
                    "ABSENT machinery requires the exclusive idiom [NONE]",
                )
            )
        if fact_value(definition) != "NOT_APPLICABLE":
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.definition_locus",
                    "ABSENT machinery requires NOT_APPLICABLE",
                )
            )
        if items(recovery):
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.static_recovery_surface",
                    "ABSENT machinery has no current recovery-surface binding",
                )
            )
        for name, section in (
            ("representation_contract", contract),
            ("dimension_bindings", bindings),
            ("dimension_use_sites", use_sites),
        ):
            if not isinstance(section, Mapping) or section.get("applicability") != "NOT_APPLICABLE":
                result.append(
                    problem(
                        "REPORT_APPLICABILITY_MATRIX",
                        f"$.{name}",
                        "ABSENT machinery requires NOT_APPLICABLE",
                    )
                )
    else:
        if not idioms or "NONE" in idioms:
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.idiom",
                    "PRESENT machinery requires at least one non-NONE idiom",
                )
            )
        if not isinstance(definition, Mapping) or fact_value(definition) in {
            "NOT_APPLICABLE",
            "UNDETERMINED",
        }:
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.definition_locus",
                    "PRESENT machinery requires a resolved definition locus",
                )
            )
        for name, section in (
            ("representation_contract", contract),
            ("dimension_bindings", bindings),
            ("dimension_use_sites", use_sites),
        ):
            if not isinstance(section, Mapping) or section.get("applicability") != "APPLICABLE":
                result.append(
                    problem(
                        "REPORT_APPLICABILITY_MATRIX",
                        f"$.{name}",
                        "PRESENT machinery requires APPLICABLE",
                    )
                )

    quantity_items = items(quantities)
    if content == "NONE":
        if quantity_items:
            result.append(
                problem(
                    "REPORT_NONE_NA_DISTINCTION",
                    "$.quantities",
                    "NONE requires a reason-bearing empty quantity inventory",
                )
            )
        if items(basis.get("axis_evidence")) or items(basis.get("orders")):
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.basis",
                    "NONE requires empty axis_evidence and orders",
                )
            )
        if fact_value(basis.get("exponent_type")) != "n/a":
            result.append(
                problem(
                    "REPORT_NONE_NA_DISTINCTION",
                    "$.basis.exponent_type",
                    "NONE requires n/a rather than a numeric/exponent assertion",
                )
            )
        if fact_value(basis.get("fractional_exponents_present")) is not False:
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.basis.fractional_exponents_present",
                    "NONE requires false",
                )
            )
        if (
            fact_value(wl_evidence.get("states_dimensions")) is not False
            or items(wl_evidence.get("pairs"))
        ):
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.wl_evidence",
                    "NONE requires states_dimensions false and an empty "
                    "reason-bearing pairs list",
                )
            )
    elif not quantity_items:
        result.append(
            problem(
                "REPORT_APPLICABILITY_MATRIX",
                "$.quantities",
                f"{classification} requires quantity leaves, with uncertainty at the smallest leaf",
            )
        )

    if classification == "REAL+PRESENT":
        if not items(basis.get("axis_evidence")) or not items(basis.get("orders")):
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.basis",
                    "REAL+PRESENT requires non-empty axis_evidence and orders",
                )
            )
        if fact_value(basis.get("exponent_type")) in {"n/a", "UNDETERMINED"}:
            result.append(
                problem(
                    "REPORT_APPLICABILITY_MATRIX",
                    "$.basis.exponent_type",
                    "REAL+PRESENT requires a resolved applicable exponent type",
                )
            )

    if classification in {"REAL+ABSENT", "UNDETERMINED+ABSENT"}:
        for index, quantity in enumerate(quantity_items):
            if items(quantity.get("binding_refs")):
                result.append(
                    problem(
                        "REPORT_BINDING_REFERENCE",
                        f"$.quantities.items[{index}].binding_refs",
                        f"{classification} future candidates require empty binding_refs",
                    )
                )

    inventories = report_item_inventories(report)
    seen: dict[str, str] = {}
    for inventory_name in ITEM_ID_INVENTORIES:
        for index, node in enumerate(inventories[inventory_name]):
            item_id = node.get("item_id")
            if not isinstance(item_id, str):
                continue
            path = f"$.{inventory_name}[{index}].item_id"
            if item_id in seen:
                result.append(
                    problem(
                        "REPORT_ID_UNIQUENESS",
                        path,
                        f"{item_id!r} already defined at {seen[item_id]}",
                    )
                )
            else:
                seen[item_id] = path

    identity_keys: list[tuple[str, tuple[Any, ...], str]] = []
    for index, node in enumerate(inventories["dimension_bindings"]):
        identity_keys.append(
            (
                "binding",
                (
                    fact_value(node.get("qualified_scope")),
                    node.get("binding_locus"),
                    fact_value(node.get("branch_or_state_predicate")),
                ),
                f"$.dimension_bindings.items[{index}]",
            )
        )
    for index, node in enumerate(inventories["dimension_use_sites"]):
        identity_keys.append(
            ("use-site", (node.get("locus"),), f"$.dimension_use_sites.items[{index}]")
        )
    for index, node in enumerate(inventories["quantities"]):
        identity_keys.append(
            (
                "quantity",
                (fact_value(node.get("name")), node.get("binding_locus")),
                f"$.quantities.items[{index}]",
            )
        )
    contract = report.get("representation_contract")
    if isinstance(contract, Mapping) and contract.get("applicability") == "APPLICABLE":
        for key in (
            "normalization_on_construction",
            "equality_semantics",
            "hashability",
            "component_access",
        ):
            node = contract.get(key, {})
            identity_keys.append(
                (
                    "contract",
                    (fact_value(node.get("qualified_name")), node.get("feature")),
                    f"$.representation_contract.{key}",
                )
            )
        for index, node in enumerate(items(contract.get("algebra"))):
            identity_keys.append(
                (
                    "contract",
                    (
                        fact_value(node.get("qualified_name")),
                        node.get("feature"),
                        fact_value(node.get("semantic_operation")),
                    ),
                    f"$.representation_contract.algebra.items[{index}]",
                )
            )
        evaluator = contract.get("expression_evaluator", {})
        for index, node in enumerate(items(evaluator.get("supported_sympy_node_kinds"))):
            identity_keys.append(
                (
                    "contract",
                    (
                        fact_value(node.get("qualified_name")),
                        node.get("feature"),
                        fact_value(node.get("node_kind")),
                    ),
                    "$.representation_contract.expression_evaluator."
                    f"supported_sympy_node_kinds.items[{index}]",
                )
            )
        for key in (
            "numeric_zero_policy",
            "missing_symbol_policy",
            "non_homogeneous_add_policy",
        ):
            node = evaluator.get(key, {})
            if isinstance(node, Mapping) and node.get("applicability") != "NOT_APPLICABLE":
                identity_keys.append(
                    (
                        "contract",
                        (fact_value(node.get("qualified_name")), node.get("feature")),
                        f"$.representation_contract.expression_evaluator.{key}",
                    )
                )
        for index, node in enumerate(items(contract.get("presentation"))):
            identity_keys.append(
                (
                    "presentation",
                    (
                        fact_value(node.get("producer_qualified_name")),
                        fact_value(node.get("form")),
                        node.get("consumer_locus"),
                    ),
                    f"$.representation_contract.presentation.items[{index}]",
                )
            )
    pass_lines = report.get("pass_lines", {})
    for section in ("emitter_sites", "call_sites"):
        for index, node in enumerate(items(pass_lines.get(section))):
            identity_keys.append(
                (
                    section,
                    (node.get("locus"),),
                    f"$.pass_lines.{section}.items[{index}]",
                )
            )
    seen_identity_keys: dict[tuple[str, tuple[Any, ...]], str] = {}
    for kind, key, path in identity_keys:
        compound = (kind, key)
        if compound in seen_identity_keys:
            result.append(
                problem(
                    "REPORT_ITEM_IDENTITY",
                    path,
                    f"duplicate §4.0 identity key {key!r}; first defined at "
                    f"{seen_identity_keys[compound]}",
                )
            )
        else:
            seen_identity_keys[compound] = path

    binding_items = inventories["dimension_bindings"]
    binding_ids = item_id_set(binding_items)
    recovery_refs: list[str] = []
    for index, recovery_item in enumerate(inventories["static_recovery_surface"]):
        ref = recovery_item.get("binding_ref")
        if ref not in binding_ids:
            result.append(
                problem(
                    "REPORT_BINDING_REFERENCE",
                    f"$.static_recovery_surface.items[{index}].binding_ref",
                    f"{ref!r} is not a dimension_bindings item_id",
                )
            )
        elif ref in recovery_refs:
            result.append(
                problem(
                    "REPORT_RECOVERY_BINDING_CORRESPONDENCE",
                    f"$.static_recovery_surface.items[{index}].binding_ref",
                    f"{ref!r} is referenced by more than one recovery item",
                )
            )
        recovery_refs.append(ref)
        visibility = recovery_item.get("current_checker_visibility")
        evidence = (
            visibility.get("loci", [])
            if isinstance(visibility, Mapping) and visibility.get("value") != "UNDETERMINED"
            else visibility.get("sources_searched", [])
            if isinstance(visibility, Mapping)
            else []
        )
        has_checker = any(CHECKER_LOCUS_PATTERN.search(str(locus)) for locus in evidence)
        has_script = any(not CHECKER_LOCUS_PATTERN.search(str(locus)) for locus in evidence)
        if not (has_checker and has_script):
            result.append(
                problem(
                    "REPORT_VISIBILITY_EVIDENCE",
                    f"$.static_recovery_surface.items[{index}].current_checker_visibility",
                    "must cite both a script binding/expression and composite_build.py",
                )
            )
    if machinery == "PRESENT" and set(recovery_refs) != binding_ids:
        result.append(
            problem(
                "REPORT_RECOVERY_BINDING_CORRESPONDENCE",
                "$.static_recovery_surface",
                "recovery binding_refs must cover every dimension-valued "
                "dimension_bindings item exactly once",
            )
        )

    for inventory_name in ("quantities", "dimension_use_sites"):
        for index, node in enumerate(inventories[inventory_name]):
            for ref_index, ref in enumerate(items(node.get("binding_refs"))):
                if ref not in binding_ids:
                    result.append(
                        problem(
                            "REPORT_BINDING_REFERENCE",
                            f"$.{inventory_name}.items[{index}].binding_refs.items[{ref_index}]",
                            f"{ref!r} is not a dimension_bindings item_id",
                        )
                    )

    ownership_map = {
        "MODULE_OWNED": "YES",
        "SCRIPT_OWNED": "NO",
        "CHECK_LOCAL": "NO",
        "UNDETERMINED": "UNDETERMINED",
    }
    negative_control_binding_ids = {
        binding.get("item_id")
        for binding in binding_items
        if "NEGATIVE_CONTROL"
        in {fact_value(use) for use in items(binding.get("uses"))}
    }
    for index, quantity in enumerate(inventories["quantities"]):
        ownership = fact_value(quantity.get("ownership"))
        declaration = quantity.get("declaration_required")
        if ownership in ownership_map and declaration != ownership_map[ownership]:
            result.append(
                problem(
                    "REPORT_DECLARATION_FROM_OWNERSHIP",
                    f"$.quantities.items[{index}].declaration_required",
                    f"{ownership} requires {ownership_map[ownership]}, got {declaration!r}",
                )
            )
        refs = items(quantity.get("binding_refs"))
        if (
            refs
            and all(ref in negative_control_binding_ids for ref in refs)
            and (ownership == "MODULE_OWNED" or declaration == "YES")
        ):
            result.append(
                problem(
                    "REPORT_QUANTITY_NEGATIVE_CONTROL_PROMOTION",
                    f"$.quantities.items[{index}]",
                    "a non-empty binding_refs set containing only NEGATIVE_CONTROL "
                    "bindings cannot become a module-owned required declaration",
                )
            )

    for index, binding in enumerate(binding_items):
        binding_id = binding.get("item_id")
        ownership = fact_value(binding.get("ownership"))
        uses = [fact_value(use) for use in items(binding.get("uses"))]
        if "NEGATIVE_CONTROL" in uses and ownership == "MODULE_OWNED":
            result.append(
                problem(
                    "REPORT_NEGATIVE_CONTROL_OWNERSHIP",
                    f"$.dimension_bindings.items[{index}].ownership",
                    "NEGATIVE_CONTROL is incompatible with MODULE_OWNED",
                )
            )
        for ref_index, ref in enumerate(items(binding.get("alias_of"))):
            ref_path = f"$.dimension_bindings.items[{index}].alias_of.items[{ref_index}]"
            if ref not in binding_ids:
                result.append(
                    problem(
                        "REPORT_ALIAS_REFERENCE",
                        ref_path,
                        f"{ref!r} is not a dimension_bindings item_id",
                    )
                )
            elif ref == binding_id:
                result.append(
                    problem(
                        "REPORT_ALIAS_REFERENCE",
                        ref_path,
                        "alias_of cannot target the binding itself",
                    )
                )

    order_pairs: set[tuple[str, tuple[str, ...]]] = set()
    interpreted_axes = {
        fact_value(axis.get("interpreted_axis_name"))
        for axis in items(basis.get("axis_evidence"))
        if isinstance(axis, Mapping)
    }
    for index, order in enumerate(items(basis.get("orders"))):
        pair = (order.get("kind"), tuple(order.get("axes", [])))
        if pair in order_pairs:
            result.append(
                problem(
                    "REPORT_ORDER_IDENTITY",
                    f"$.basis.orders.items[{index}]",
                    f"duplicate (kind, axes) order item {pair!r}",
                )
            )
        order_pairs.add(pair)
        for axis_index, axis_name in enumerate(order.get("axes", [])):
            if axis_name not in interpreted_axes:
                result.append(
                    problem(
                        "REPORT_ORDER_AXIS_REFERENCE",
                        f"$.basis.orders.items[{index}].axes[{axis_index}]",
                        f"{axis_name!r} is not an axis_evidence interpreted_axis_name",
                    )
                )
    for index, order in enumerate(items(basis.get("orders"))):
        reference = order.get("same_as_order")
        if reference is not None:
            target = (reference.get("kind"), tuple(reference.get("axes", [])))
            current = (order.get("kind"), tuple(order.get("axes", [])))
            if target not in order_pairs or target == current:
                result.append(
                    problem(
                        "REPORT_ORDER_IDENTITY",
                        f"$.basis.orders.items[{index}].same_as_order",
                        "must reference a distinct existing (kind, axes) item",
                    )
                )

    feature_slots = {
        "normalization_on_construction": "CONSTRUCTION_NORMALIZATION",
        "equality_semantics": "EQUALITY_SEMANTICS",
        "hashability": "HASHABILITY",
        "component_access": "COMPONENT_ACCESS",
    }
    if isinstance(contract, Mapping) and contract.get("applicability") == "APPLICABLE":
        for key, expected_feature in feature_slots.items():
            if get_path(contract, (key, "feature")) != expected_feature:
                result.append(
                    problem(
                        "SCHEMA_STRUCTURE",
                        f"$.representation_contract.{key}.feature",
                        f"expected {expected_feature}",
                    )
                )
        for index, algebra in enumerate(items(contract.get("algebra"))):
            if algebra.get("feature") != "ALGEBRA_OPERATION":
                result.append(
                    problem(
                        "SCHEMA_STRUCTURE",
                        f"$.representation_contract.algebra.items[{index}].feature",
                        "expected ALGEBRA_OPERATION",
                    )
                )
        evaluator = contract.get("expression_evaluator", {})
        status = fact_value(evaluator.get("status"))
        node_items = items(evaluator.get("supported_sympy_node_kinds"))
        policies = (
            ("numeric_zero_policy", "NUMERIC_ZERO_POLICY"),
            ("missing_symbol_policy", "MISSING_SYMBOL_POLICY"),
            ("non_homogeneous_add_policy", "NON_HOMOGENEOUS_ADD_POLICY"),
        )
        if status == "ABSENT":
            if node_items:
                result.append(
                    problem(
                        "SCHEMA_STRUCTURE",
                        "$.representation_contract.expression_evaluator.supported_sympy_node_kinds",
                        "ABSENT evaluator cannot have node-kind feature items",
                    )
                )
            for key, _ in policies:
                if get_path(evaluator, (key, "applicability")) != "NOT_APPLICABLE":
                    result.append(
                        problem(
                            "SCHEMA_STRUCTURE",
                            f"$.representation_contract.expression_evaluator.{key}",
                            "ABSENT evaluator requires NOT_APPLICABLE",
                        )
                    )
        elif status == "PRESENT":
            for index, node in enumerate(node_items):
                if node.get("feature") != "EVALUATOR_NODE_KIND":
                    result.append(
                        problem(
                            "SCHEMA_STRUCTURE",
                            "$.representation_contract.expression_evaluator."
                            f"supported_sympy_node_kinds.items[{index}].feature",
                            "expected EVALUATOR_NODE_KIND",
                        )
                    )
            for key, expected_feature in policies:
                if get_path(evaluator, (key, "feature")) != expected_feature:
                    result.append(
                        problem(
                            "SCHEMA_STRUCTURE",
                            f"$.representation_contract.expression_evaluator.{key}.feature",
                            f"expected {expected_feature}",
                        )
                    )

    if classification == "REAL+PRESENT":
        yes_count = sum(
            quantity.get("declaration_required") == "YES" for quantity in quantity_items
        )
        zero_argument = (
            quantities.get("zero_declaration_argument") if isinstance(quantities, Mapping) else None
        )
        has_zero_argument = (
            isinstance(zero_argument, Mapping)
            and fact_value(zero_argument) not in {None, "NOT_APPLICABLE"}
        )
        if yes_count == 0 and not has_zero_argument:
            result.append(
                problem(
                    "REPORT_ZERO_DECLARATION_ARGUMENT",
                    "$.quantities.zero_declaration_argument",
                    "REAL+PRESENT with zero declaration_required YES needs a source-cited argument",
                )
            )
    return result


def resolve_report_path(raw: str, verification_path: Path) -> Path:
    candidate = Path(raw)
    if candidate.is_absolute():
        return candidate
    repo_candidate = REPO_ROOT / candidate
    if repo_candidate.exists():
        return repo_candidate
    return verification_path.parent / candidate


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_source_identity_path(raw: str, report_path: Path) -> Path:
    candidate = Path(raw)
    if candidate.is_absolute():
        return candidate
    repo_candidate = REPO_ROOT / candidate
    if repo_candidate.exists():
        return repo_candidate
    return report_path.parent / candidate


def source_identity_failures(
    report: Mapping[str, Any],
    report_path: Path,
) -> list[tuple[str, str]]:
    identity = report.get("identity")
    if not isinstance(identity, Mapping):
        return []
    failures: list[tuple[str, str]] = []
    for path_field, hash_field in (
        ("script_path", "script_sha256"),
        ("wl_path", "wl_sha256"),
    ):
        raw_path = get_path(identity, (path_field, "value"))
        if not isinstance(raw_path, str):
            continue
        resolved = resolve_source_identity_path(raw_path, report_path)
        if not resolved.exists():
            failures.append((path_field, "FILE_NOT_FOUND"))
            continue
        if not resolved.is_file():
            failures.append((path_field, "NOT_A_FILE"))
            continue
        expected_hash = get_path(identity, (hash_field, "value"))
        if expected_hash != sha256_file(resolved):
            failures.append((hash_field, "SHA256_MISMATCH"))
        if path_field == "script_path":
            with resolved.open("r", encoding="utf-8", errors="replace") as handle:
                actual_line_count = sum(1 for _ in handle)
            expected_line_count = get_path(identity, ("script_line_count", "value"))
            if expected_line_count != actual_line_count:
                failures.append(("script_line_count", "LINE_COUNT_MISMATCH"))
    field_order = {
        "script_path": 0,
        "script_sha256": 1,
        "script_line_count": 2,
        "wl_path": 3,
        "wl_sha256": 4,
    }
    return sorted(failures, key=lambda entry: (field_order.get(entry[0], 99), entry[1]))


def local_schema_ref(root_schema: Mapping[str, Any], raw_ref: str) -> Any:
    if not raw_ref.startswith("#/"):
        return None
    current: Any = root_schema
    for raw_part in raw_ref[2:].split("/"):
        part = raw_part.replace("~1", "/").replace("~0", "~")
        if not isinstance(current, Mapping) or part not in current:
            return None
        current = current[part]
    return current


def schema_locus_strings(
    document: Any,
    schema: Mapping[str, Any],
    path: str = "$",
) -> dict[str, str]:
    """Return values reached through schema nodes explicitly typed as locus."""
    found: dict[str, str] = {}

    def walk(node: Any, node_schema: Any, node_path: str) -> None:
        if not isinstance(node_schema, Mapping):
            return
        raw_ref = node_schema.get("$ref")
        if raw_ref == "#/$defs/locus":
            if isinstance(node, str):
                found.setdefault(node, node_path)
            return
        if isinstance(raw_ref, str):
            resolved = local_schema_ref(schema, raw_ref)
            if resolved is not None:
                walk(node, resolved, node_path)
        for keyword in ("allOf", "oneOf", "anyOf"):
            branches = node_schema.get(keyword)
            if isinstance(branches, list):
                selected_branches = branches
                if keyword == "oneOf":
                    valid_branches = []
                    for branch in branches:
                        branch_schema = {
                            "$schema": schema.get("$schema"),
                            "$defs": schema.get("$defs", {}),
                            **branch,
                        }
                        if Draft202012Validator(branch_schema).is_valid(node):
                            valid_branches.append(branch)
                    if valid_branches:
                        selected_branches = valid_branches
                for branch in selected_branches:
                    walk(node, branch, node_path)
        if isinstance(node, Mapping):
            properties = node_schema.get("properties")
            if isinstance(properties, Mapping):
                for key, child_schema in properties.items():
                    if key == "locus_integrity" or key not in node:
                        continue
                    walk(node[key], child_schema, f"{node_path}.{key}")
        elif isinstance(node, list):
            item_schema = node_schema.get("items")
            if isinstance(item_schema, Mapping):
                for index, value in enumerate(node):
                    walk(value, item_schema, f"{node_path}[{index}]")

    walk(document, schema, path)
    return found


def locus_file_candidates(raw_path: str, document_path: Path) -> list[Path]:
    candidate = Path(raw_path)
    if candidate.is_absolute():
        return [candidate]
    ordered = [
        REPO_ROOT / candidate,
        SCHEMA_DIR.parent / candidate,
        SCHEMA_DIR / "examples" / candidate,
        document_path.parent / candidate,
    ]
    distinct: list[Path] = []
    for path in ordered:
        if path not in distinct:
            distinct.append(path)
    return distinct


def source_tokens(raw: str) -> tuple[str, ...]:
    return tuple(SOURCE_TOKEN_PATTERN.findall(raw))


def whitespace_normalized(raw: str) -> str:
    return re.sub(r"\s+", "", raw)


def identifier_component(character: str) -> bool:
    category = unicodedata.category(character)
    return (
        character == "_"
        or category.startswith("L")
        or category.startswith("M")
        or category in {"Nd", "Pc"}
    )


def standalone_symbol_occurs(source_text: str, spelling: str) -> bool:
    start = 0
    while True:
        index = source_text.find(spelling, start)
        if index < 0:
            return False
        before = source_text[index - 1] if index > 0 else ""
        end = index + len(spelling)
        after = source_text[end] if end < len(source_text) else ""
        if (
            (not before or not identifier_component(before))
            and (not after or not identifier_component(after))
        ):
            return True
        start = index + 1


def dimension_text_has_declared_axis(
    source_text: str,
    report: Mapping[str, Any],
) -> bool:
    if whitespace_normalized(source_text) == "1":
        return True
    for axis in items(get_path(report, ("basis", "axis_evidence"))):
        if not isinstance(axis, Mapping):
            continue
        for field in ("source_spelling", "interpreted_axis_name"):
            axis_name = fact_value(axis.get(field))
            if (
                isinstance(axis_name, str)
                and axis_name
                and axis_name != "UNDETERMINED"
                and standalone_symbol_occurs(source_text, axis_name)
            ):
                return True
    return False


def python_sequence_component_count(source_text: str) -> int | None:
    try:
        module = ast.parse(textwrap.dedent(source_text))
    except (SyntaxError, ValueError):
        return None
    if len(module.body) != 1:
        return None
    statement = module.body[0]
    value: ast.expr | None
    if isinstance(statement, (ast.Assign, ast.AnnAssign)):
        value = statement.value
    elif isinstance(statement, ast.Expr):
        value = statement.value
    else:
        return None
    if isinstance(value, (ast.Tuple, ast.List)):
        return len(value.elts)
    return None


def wl_sequence_component_count(source_text: str) -> int | None:
    candidate = source_text.strip()
    if candidate.endswith(";"):
        candidate = candidate[:-1].rstrip()
    if "=" in candidate:
        candidate = candidate.split("=", 1)[1].strip()
    elif "->" in candidate:
        candidate = candidate.split("->", 1)[1].strip()
    if len(candidate) < 2 or candidate[0] != "{" or candidate[-1] != "}":
        return None

    stack: list[str] = []
    quote: str | None = None
    escaped = False
    component_count = 0
    component_has_content = False
    pairs = {")": "(", "]": "[", "}": "{"}
    for character in candidate:
        if quote is not None:
            component_has_content = True
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            continue
        if character in {"'", '"'}:
            quote = character
            component_has_content = True
        elif character in "([{":
            stack.append(character)
            if len(stack) > 1:
                component_has_content = True
        elif character in ")]}":
            if not stack or stack[-1] != pairs[character]:
                return None
            if len(stack) == 1:
                if component_has_content:
                    component_count += 1
                elif component_count:
                    return None
            else:
                component_has_content = True
            stack.pop()
        elif character == "," and len(stack) == 1:
            if not component_has_content:
                return None
            component_count += 1
            component_has_content = False
        elif not character.isspace():
            component_has_content = True
    if stack or quote is not None:
        return None
    return component_count


def positional_component_count(source_text: str) -> int | None:
    python_count = python_sequence_component_count(source_text)
    if python_count is not None:
        return python_count
    return wl_sequence_component_count(source_text)


def declared_basis_loci(report: Mapping[str, Any]) -> set[str]:
    result: set[str] = set()
    for path in (("basis", "axis_evidence"), ("basis", "orders")):
        for entry in items(get_path(report, path)):
            locus = entry.get("locus") if isinstance(entry, Mapping) else None
            if isinstance(locus, str):
                result.add(locus)
    return result


def explicit_dimension_text_matches_declared_basis(
    quantity: Mapping[str, Any],
    report: Mapping[str, Any],
) -> bool:
    source_text = fact_value(quantity.get("dim_source_text"))
    form = quantity.get("dim_text_form")
    basis_locus = quantity.get("basis_locus")
    axes = [
        axis
        for axis in items(get_path(report, ("basis", "axis_evidence")))
        if isinstance(axis, Mapping)
    ]
    if (
        not isinstance(source_text, str)
        or not axes
        or not isinstance(basis_locus, str)
        or basis_locus not in declared_basis_loci(report)
    ):
        return False
    if form == "NAMED_AXIS":
        return dimension_text_has_declared_axis(source_text, report)
    if form == "POSITIONAL":
        return positional_component_count(source_text) == len(axes)
    return False


def locus_source_file(raw_locus: str) -> str | None:
    parsed = parsed_locus(raw_locus)
    return parsed[0] if parsed is not None else None


def locus_matches_identity_source(
    raw_locus: str,
    raw_identity_path: Any,
    document_path: Path,
) -> bool:
    raw_file = locus_source_file(raw_locus)
    if raw_file is None or not isinstance(raw_identity_path, str):
        return False
    if posixpath.normpath(raw_file) == posixpath.normpath(raw_identity_path):
        return True
    identity_path = resolve_source_identity_path(raw_identity_path, document_path)
    try:
        identity_resolved = identity_path.resolve()
    except OSError:
        identity_resolved = identity_path.absolute()
    for candidate in locus_file_candidates(raw_file, document_path):
        try:
            resolved = candidate.resolve()
        except OSError:
            resolved = candidate.absolute()
        if resolved == identity_resolved:
            return True
    return False


def is_stage_note_locus(raw_locus: str) -> bool:
    raw_file = locus_source_file(raw_locus)
    if raw_file is None:
        return False
    normalized = posixpath.normpath(raw_file)
    return (
        normalized.endswith(".md")
        and (
            "/notes/stages/" in f"/{normalized.lstrip('/')}"
            or (
                "/schemas/examples/synthetic/" in f"/{normalized.lstrip('/')}"
                and normalized.endswith("_note.md")
            )
            or (
                normalized.startswith("synthetic/")
                and normalized.endswith("_note.md")
            )
        )
        and not normalized.endswith("/notes/parameter_register.md")
    )


def source_slices_at_locus(raw_locus: str, document_path: Path) -> list[str]:
    match = LOCUS_PATTERN.fullmatch(raw_locus)
    if match is None:
        return []
    start = int(match.group("start"))
    end = int(match.group("end") or start)
    if end < start:
        return []
    result: list[str] = []
    for candidate in locus_file_candidates(match.group("file"), document_path):
        if not candidate.is_file():
            continue
        try:
            lines = candidate.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        if end <= len(lines):
            result.append("\n".join(lines[start - 1 : end]))
    return result


def token_subsequence(needle: Sequence[str], haystack: Sequence[str]) -> bool:
    if not needle or len(needle) > len(haystack):
        return False
    width = len(needle)
    return any(tuple(haystack[index : index + width]) == tuple(needle) for index in range(len(haystack) - width + 1))


def source_text_occurs_at_locus(
    source_text: str,
    raw_locus: str,
    document_path: Path,
) -> bool:
    needle = source_tokens(source_text)
    occurs = any(
        token_subsequence(needle, source_tokens(source_slice))
        for source_slice in source_slices_at_locus(raw_locus, document_path)
    )
    return bool(needle) and occurs


def markdown_table_cells(line: str) -> tuple[str, ...] | None:
    stripped = line.strip()
    if not stripped.startswith("|") or not stripped.endswith("|"):
        return None
    return tuple(cell.strip() for cell in re.split(r"(?<!\\)\|", stripped[1:-1]))


def master_register_row(
    raw_locus: str,
    document_path: Path | None,
) -> RegisterRow | None:
    if document_path is None:
        return None
    locus = parsed_locus(raw_locus)
    if locus is None:
        return None
    raw_file, start, end = locus
    if start != end:
        return None
    register_candidates = locus_file_candidates(raw_file, document_path)
    if not any(
        candidate.is_file() and candidate.resolve() == REAL_REGISTER_PATH
        for candidate in register_candidates
    ):
        return None
    try:
        lines = REAL_REGISTER_PATH.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
    except OSError:
        return None
    header_indexes = []
    for index, line in enumerate(lines):
        cells = markdown_table_cells(line)
        if cells is None:
            continue
        normalized = tuple(cell.replace("`", "") for cell in cells)
        if normalized == REGISTER_HEADER:
            header_indexes.append(index)
    for header_index in header_indexes:
        row_index = header_index + 2
        while row_index < len(lines):
            cells = markdown_table_cells(lines[row_index])
            if cells is None:
                break
            one_based = row_index + 1
            # The committed master table has six cells through row 185 and
            # five cells at rows 186--191.  Param and dimension are the stable
            # first two fields in both shapes.
            if one_based == start and len(cells) >= 2:
                return RegisterRow(
                    line=one_based,
                    param=cells[0],
                    dimension_cell=cells[1],
                    cells=cells,
                )
            row_index += 1
    return None


def locus_failure(raw_locus: str, document_path: Path) -> str | None:
    match = LOCUS_PATTERN.fullmatch(raw_locus)
    if match is None:
        return "FILE_NOT_FOUND"
    candidates = locus_file_candidates(match.group("file"), document_path)
    existing = [candidate for candidate in candidates if candidate.exists()]
    if not existing:
        return "FILE_NOT_FOUND"
    files = [candidate for candidate in existing if candidate.is_file()]
    if not files:
        return "NOT_A_FILE"
    start = int(match.group("start"))
    end = int(match.group("end") or start)
    for candidate in files:
        try:
            with candidate.open("r", encoding="utf-8", errors="replace") as handle:
                line_count = sum(1 for _ in handle)
        except OSError:
            continue
        if start <= end <= line_count:
            return None
    return "LINE_OUT_OF_RANGE"


def mechanical_locus_failures(
    report: Mapping[str, Any],
    report_path: Path,
    report_schema: Mapping[str, Any],
    verification: Mapping[str, Any],
    verification_path: Path,
    verification_schema: Mapping[str, Any],
) -> list[tuple[str, str]]:
    failures: dict[str, str] = {}
    for locus in schema_locus_strings(report, report_schema):
        failure = locus_failure(locus, report_path)
        if failure is not None:
            failures[locus] = failure
    for locus in schema_locus_strings(verification, verification_schema):
        failure = locus_failure(locus, verification_path)
        if failure is not None:
            failures[locus] = failure
    return sorted(failures.items())


def category_report_inventories(report: Mapping[str, Any]) -> dict[str, list[Mapping[str, Any]]]:
    inventories = report_item_inventories(report)
    return {
        "dimension_bindings": inventories["dimension_bindings"],
        "quantities": inventories["quantities"],
        "contract_features_and_presentation_consumers": (
            inventories["contract_features"] + inventories["presentation_consumers"]
        ),
        "pass_emitter_sites": inventories["pass_emitter_sites"],
        "pass_call_sites": inventories["pass_call_sites"],
    }


def verification_semantic_problems(
    verification: Mapping[str, Any],
    verification_path: Path,
    report_schema: Mapping[str, Any],
    verification_schema: Mapping[str, Any],
) -> list[Problem]:
    result: list[Problem] = []
    identity = verification.get("identity", {})
    raw_report_path = get_path(identity, ("report_path", "value"))
    if not isinstance(raw_report_path, str):
        return [problem("VER_REPORT_REFERENCE", "$.identity.report_path", "missing report path")]
    report_path = resolve_report_path(raw_report_path, verification_path)
    if not report_path.is_file():
        return [
            problem(
                "VER_REPORT_REFERENCE",
                "$.identity.report_path.value",
                f"referenced report does not exist: {report_path}",
            )
        ]
    expected_hash = get_path(identity, ("report_sha256", "value"))
    actual_hash = sha256_file(report_path)
    if expected_hash != actual_hash:
        result.append(
            problem(
                "VER_REPORT_REFERENCE",
                "$.identity.report_sha256.value",
                f"expected referenced report sha256 {actual_hash}, got {expected_hash!r}",
            )
        )
    try:
        report = load_effective_document(report_path)
    except (OSError, yaml.YAMLError, ValueError) as exc:
        result.append(
            problem(
                "VER_REPORT_REFERENCE",
                "$.identity.report_path.value",
                f"could not load referenced report: {exc}",
            )
        )
        return result
    report_errors = report_preflight(report, report_path)
    if not report_errors:
        report_errors = schema_problems(report_schema, report)
    if not report_errors:
        report_errors = report_semantic_problems(report)
    if report_errors:
        result.append(
            problem(
                "VER_REPORT_REFERENCE",
                "$.identity.report_path.value",
                f"referenced report is invalid ({report_errors[0].render()})",
            )
        )
        return result
    if not isinstance(report, Mapping):
        return result

    locus_check = verification.get("locus_integrity", {})
    locus_status = locus_check.get("status") if isinstance(locus_check, Mapping) else None
    if locus_status == "SYNTHETIC_EXEMPT":
        examples_root = (SCHEMA_DIR / "examples").resolve()
        try:
            synthetic_report = report_path.resolve().is_relative_to(examples_root)
        except OSError:
            synthetic_report = False
        if not synthetic_report:
            result.append(
                problem(
                    "VER_LOCUS_INTEGRITY",
                    "$.locus_integrity.status",
                    "SYNTHETIC_EXEMPT is limited to committed schema examples",
                )
            )
    elif locus_status == "CHECKED":
        actual_failures = mechanical_locus_failures(
            report,
            report_path,
            report_schema,
            verification,
            verification_path,
            verification_schema,
        )
        claimed_failures = sorted(
            (
                entry.get("locus"),
                entry.get("failure"),
            )
            for entry in items(locus_check.get("invalid_loci"))
            if isinstance(entry, Mapping)
        )
        expected_result = "FAIL" if actual_failures else "PASS"
        if locus_check.get("result") != expected_result or claimed_failures != actual_failures:
            first = actual_failures[0] if actual_failures else None
            result.append(
                problem(
                    "VER_LOCUS_INTEGRITY",
                    "$.locus_integrity",
                    f"mechanical result is {expected_result}; "
                    + (
                        f"first invalid locus is {first[0]!r} ({first[1]})"
                        if first is not None
                        else "no invalid loci were found"
                    ),
                )
            )
        if actual_failures and verification.get("overall_verdict") != "FAIL":
            result.append(
                problem(
                    "VER_LOCUS_INTEGRITY",
                    "$.overall_verdict",
                    "mechanical locus failure requires overall FAIL",
                )
            )
        actual_identity_failures = source_identity_failures(report, report_path)
        identity_check = locus_check.get("source_identity")
        claimed_identity_failures = sorted(
            (
                entry.get("field"),
                entry.get("failure"),
            )
            for entry in items(
                identity_check.get("failures")
                if isinstance(identity_check, Mapping)
                else None
            )
            if isinstance(entry, Mapping)
        )
        expected_identity_result = "FAIL" if actual_identity_failures else "PASS"
        if (
            not isinstance(identity_check, Mapping)
            or identity_check.get("result") != expected_identity_result
            or claimed_identity_failures != sorted(actual_identity_failures)
        ):
            first = actual_identity_failures[0] if actual_identity_failures else None
            result.append(
                problem(
                    "VER_IDENTITY_INTEGRITY",
                    "$.locus_integrity.source_identity",
                    f"mechanical source identity result is {expected_identity_result}; "
                    + (
                        f"first failure is {first[0]} ({first[1]})"
                        if first is not None
                        else "all source identity fields reproduce"
                    ),
                )
            )
        if actual_identity_failures and verification.get("overall_verdict") != "FAIL":
            result.append(
                problem(
                    "VER_IDENTITY_INTEGRITY",
                    "$.overall_verdict",
                    "mechanical source identity failure requires overall FAIL",
                )
            )

    report_classification = classification_for_report(report)
    classification = verification.get("classification")
    if classification != report_classification:
        result.append(
            problem(
                "VER_CLASSIFICATION_MATCH",
                "$.classification",
                f"report is {report_classification}, verification says {classification!r}",
            )
        )
    if identity.get("stage") != get_path(report, ("identity", "stage")):
        result.append(
            problem(
                "VER_REPORT_REFERENCE",
                "$.identity.stage",
                "stage does not match referenced report",
            )
        )

    report_categories = category_report_inventories(report)
    categories = verification["categories"]
    all_report_ids: set[str] = set()
    for category_name, report_nodes in report_categories.items():
        category = categories[category_name]
        if category.get("applicability") != "APPLICABLE":
            continue
        expected_count = len(report_nodes)
        if category.get("report_count") != expected_count:
            result.append(
                problem(
                    "VER_REPORT_COUNT",
                    f"$.categories.{category_name}.report_count",
                    f"referenced report contains {expected_count} items",
                )
            )
        source_unmatched = items(category.get("source_unmatched"))
        report_unmatched = items(category.get("report_unmatched"))
        if (
            category.get("source_count") - category.get("report_count")
            != len(source_unmatched) - len(report_unmatched)
        ):
            result.append(
                problem(
                    "VER_CARDINALITY_RECONCILIATION",
                    f"$.categories.{category_name}",
                    "source/report cardinality difference does not match unmatched-list difference",
                )
            )
        category_report_ids = item_id_set(report_nodes)
        all_report_ids.update(category_report_ids)
        for index, entry in enumerate(report_unmatched):
            if entry.get("item_id") not in category_report_ids:
                result.append(
                    problem(
                        "VER_REPORT_REFERENCE",
                        f"$.categories.{category_name}.report_unmatched.items[{index}].item_id",
                        "report-side unmatched item_id is not in this report category",
                    )
                )

    unexplained_paths: list[str] = []
    source_seen: dict[str, str] = {}
    report_seen: dict[str, str] = {}
    for category_name in VERIFICATION_CATEGORIES:
        category = categories[category_name]
        for side, seen in (("source_unmatched", source_seen), ("report_unmatched", report_seen)):
            for index, entry in enumerate(items(category.get(side))):
                path = f"$.categories.{category_name}.{side}.items[{index}]"
                item_id = entry.get("item_id")
                if item_id in seen:
                    result.append(
                        problem(
                            "VER_REPORT_REFERENCE",
                            f"{path}.item_id",
                            f"{item_id!r} already appears at {seen[item_id]}",
                        )
                    )
                else:
                    seen[item_id] = path
                if entry.get("explanation_status") == "UNEXPLAINED":
                    unexplained_paths.append(path)
    if unexplained_paths and verification.get("overall_verdict") != "FAIL":
        result.append(
            problem(
                "VER_UNEXPLAINED_REQUIRES_FAIL",
                "$.overall_verdict",
                f"UNEXPLAINED unmatched item at {unexplained_paths[0]} requires FAIL",
            )
        )

    none_recheck = verification.get("none_classification_recheck", {})
    if classification == "NONE+ABSENT":
        if none_recheck.get("applicability") != "APPLICABLE":
            result.append(
                problem(
                    "VER_NONE_RECHECK",
                    "$.none_classification_recheck",
                    "NONE+ABSENT requires an applicable positive no-declaration recheck",
                )
            )
        elif none_recheck.get("result") == "FAIL" and verification.get("overall_verdict") != "FAIL":
            result.append(
                problem(
                    "VER_NONE_RECHECK",
                    "$.overall_verdict",
                    "failed NONE recheck requires FAIL",
                )
            )
    elif none_recheck.get("applicability") != "NOT_APPLICABLE":
        result.append(
            problem(
                "VER_NONE_RECHECK",
                "$.none_classification_recheck",
                "only NONE+ABSENT uses the no-declaration recheck",
            )
        )

    binding_nodes = report_categories["dimension_bindings"]
    binding_by_id = {node["item_id"]: node for node in binding_nodes}
    ownership = verification.get("ownership_rederivation", {})
    present = classification in {"REAL+PRESENT", "UNDETERMINED+PRESENT"}
    if present:
        if ownership.get("applicability") != "APPLICABLE":
            result.append(
                problem(
                    "VER_OWNERSHIP_MATCH",
                    "$.ownership_rederivation",
                    "PRESENT requires applicable ownership re-derivation",
                )
            )
        else:
            seen_refs: set[str] = set()
            ownership_items = [
                entry for entry in ownership.get("items", []) if isinstance(entry, Mapping)
            ]
            stage = get_path(report, ("identity", "stage"))
            expected_case = "i" if stage in OWNERSHIP_ANCHOR_STAGES else "ii"
            coverage_case = ownership.get("coverage_case")
            anchor_record = ownership.get("applicable_ownership_anchors", {})
            sample_record = ownership.get("case_i_sample", {})
            coverage = ownership.get("coverage", {})
            escape = coverage.get("escape", {}) if isinstance(coverage, Mapping) else {}
            anchor_refs = (
                set(items(anchor_record.get("binding_refs")))
                if isinstance(anchor_record, Mapping)
                and anchor_record.get("value") == "PRESENT"
                else set()
            )
            sample_refs = (
                set(items(sample_record.get("binding_refs")))
                if isinstance(sample_record, Mapping)
                and sample_record.get("applicability") == "APPLICABLE"
                else set()
            )
            selected_by_kind = {
                "OWNERSHIP_ANCHOR": {
                    entry.get("binding_ref")
                    for entry in ownership_items
                    if entry.get("selection_kind") == "OWNERSHIP_ANCHOR"
                },
                "CASE_I_SAMPLE": {
                    entry.get("binding_ref")
                    for entry in ownership_items
                    if entry.get("selection_kind") == "CASE_I_SAMPLE"
                },
                "CASE_II_FULL_COVERAGE": {
                    entry.get("binding_ref")
                    for entry in ownership_items
                    if entry.get("selection_kind") == "CASE_II_FULL_COVERAGE"
                },
            }
            coverage_problems: list[str] = []
            if coverage_case != expected_case:
                coverage_problems.append(
                    f"stage {stage!r} requires coverage case ({expected_case})"
                )
            if coverage.get("bindings_total") != len(binding_nodes):
                coverage_problems.append(
                    f"bindings_total must equal report binding count {len(binding_nodes)}"
                )
            if coverage.get("bindings_re_derived") != len(ownership_items):
                coverage_problems.append(
                    f"bindings_re_derived must equal selected item count {len(ownership_items)}"
                )
            missing_count = len(binding_nodes) - len(ownership_items)
            if escape.get("invoked") is True:
                if escape.get("bindings_not_rederived") != missing_count:
                    coverage_problems.append(
                        f"escape count must equal unselected binding count {missing_count}"
                    )
            elif escape.get("invoked") is False and coverage_case == "ii" and missing_count:
                coverage_problems.append(
                    "case (ii) without an escape must re-derive every report binding"
                )
            if coverage_case == "i":
                if anchor_record.get("value") != "PRESENT" or not anchor_refs:
                    coverage_problems.append(
                        "case (i) requires a non-empty applicable ownership-anchor reference set"
                    )
                if locus_status == "CHECKED":
                    directive_refs, anchor_resolution_errors = directive_anchor_refs(
                        str(stage),
                        binding_nodes,
                    )
                    coverage_problems.extend(anchor_resolution_errors)
                    if anchor_refs != directive_refs:
                        coverage_problems.append(
                            "declared anchor refs do not equal the r17 "
                            "script/scope/name/branch-state property set"
                        )
                remainder = set(binding_by_id) - anchor_refs
                if remainder:
                    if (
                        sample_record.get("applicability") != "APPLICABLE"
                        or not sample_refs
                    ):
                        coverage_problems.append(
                            "case (i) requires a reproducible non-empty sample of the remainder"
                        )
                elif sample_record.get("applicability") != "NOT_APPLICABLE":
                    coverage_problems.append(
                        "case (i) with no remainder requires an inapplicable sample record"
                    )
                if anchor_refs & sample_refs:
                    coverage_problems.append("anchor and sample reference sets must be disjoint")
                if not sample_refs <= remainder:
                    coverage_problems.append(
                        "case-(i) sample references must come from the non-anchor remainder"
                    )
                if selected_by_kind["OWNERSHIP_ANCHOR"] != anchor_refs:
                    coverage_problems.append(
                        "OWNERSHIP_ANCHOR item labels must equal the declared anchor refs"
                    )
                if selected_by_kind["CASE_I_SAMPLE"] != sample_refs:
                    coverage_problems.append(
                        "CASE_I_SAMPLE item labels must equal the recorded sample refs"
                    )
                if selected_by_kind["CASE_II_FULL_COVERAGE"]:
                    coverage_problems.append(
                        "case (i) cannot use CASE_II_FULL_COVERAGE item labels"
                    )
            elif coverage_case == "ii":
                if anchor_record.get("value") != "NONE":
                    coverage_problems.append(
                        "case (ii) must positively declare applicable_ownership_anchors NONE"
                    )
                if sample_record.get("applicability") != "NOT_APPLICABLE":
                    coverage_problems.append(
                        "case (ii) cannot carry a case-(i) sample"
                    )
                if (
                    selected_by_kind["OWNERSHIP_ANCHOR"]
                    or selected_by_kind["CASE_I_SAMPLE"]
                ):
                    coverage_problems.append(
                        "case (ii) cannot relabel full-coverage bindings as anchors or samples"
                    )
                if escape.get("invoked") is False and (
                    selected_by_kind["CASE_II_FULL_COVERAGE"] != set(binding_by_id)
                ):
                    coverage_problems.append(
                        "case-(ii) full-coverage labels must cover every report binding"
                    )
            if coverage_problems:
                result.append(
                    problem(
                        "VER_OWNERSHIP_ANCHOR_COVERAGE",
                        "$.ownership_rederivation",
                        coverage_problems[0],
                    )
                )
            for index, entry in enumerate(ownership_items):
                ref = entry.get("binding_ref")
                path = f"$.ownership_rederivation.items[{index}]"
                if ref not in binding_by_id:
                    result.append(
                        problem(
                            "VER_REPORT_REFERENCE",
                            f"{path}.binding_ref",
                            "ownership reference is not a report binding item_id",
                        )
                    )
                    continue
                if ref in seen_refs:
                    result.append(
                        problem(
                            "VER_REPORT_REFERENCE",
                            f"{path}.binding_ref",
                            "binding is selected more than once",
                        )
                    )
                seen_refs.add(ref)
                actual = fact_value(binding_by_id[ref].get("ownership"))
                if entry.get("reported_ownership") != actual:
                    result.append(
                        problem(
                            "VER_OWNERSHIP_MATCH",
                            f"{path}.reported_ownership",
                            f"report records {actual!r}",
                        )
                    )
                derived_match = entry.get("reported_ownership") == entry.get("rederived_ownership")
                expected_comparison = "MATCH" if derived_match else "MISMATCH"
                if entry.get("comparison") != expected_comparison:
                    result.append(
                        problem(
                            "VER_OWNERSHIP_MATCH",
                            f"{path}.comparison",
                            f"expected {expected_comparison}",
                        )
                    )
                if not derived_match and verification.get("overall_verdict") != "FAIL":
                    result.append(
                        problem(
                            "VER_OWNERSHIP_MATCH",
                            "$.overall_verdict",
                            f"ownership mismatch at {path} requires FAIL",
                        )
                    )
    elif ownership.get("applicability") != "NOT_APPLICABLE":
        result.append(
            problem(
                "VER_OWNERSHIP_MATCH",
                "$.ownership_rederivation",
                "ABSENT machinery requires NOT_APPLICABLE",
            )
        )

    zero_check = verification.get("zero_declaration_check", {})
    if classification == "REAL+PRESENT":
        if zero_check.get("applicability") != "APPLICABLE":
            result.append(
                problem(
                    "VER_ZERO_DECLARATION",
                    "$.zero_declaration_check",
                    "REAL+PRESENT requires the zero-declaration check",
                )
            )
        else:
            quantity_nodes = report_categories["quantities"]
            actual_yes = sum(node.get("declaration_required") == "YES" for node in quantity_nodes)
            zero_argument = report.get("quantities", {}).get("zero_declaration_argument")
            argument_present = (
                isinstance(zero_argument, Mapping)
                and fact_value(zero_argument) not in {None, "NOT_APPLICABLE"}
            )
            if zero_check.get("declaration_required_yes_count") != actual_yes:
                result.append(
                    problem(
                        "VER_ZERO_DECLARATION",
                        "$.zero_declaration_check.declaration_required_yes_count",
                        f"report contains {actual_yes}",
                    )
                )
            if zero_check.get("source_cited_argument_present") is not argument_present:
                result.append(
                    problem(
                        "VER_ZERO_DECLARATION",
                        "$.zero_declaration_check.source_cited_argument_present",
                        f"report argument presence is {argument_present}",
                    )
                )
            expected_result = "PASS" if actual_yes > 0 or argument_present else "FAIL"
            if zero_check.get("result") != expected_result:
                result.append(
                    problem(
                        "VER_ZERO_DECLARATION",
                        "$.zero_declaration_check.result",
                        f"expected {expected_result}",
                    )
                )
            if expected_result == "FAIL" and verification.get("overall_verdict") != "FAIL":
                result.append(
                    problem(
                        "VER_ZERO_DECLARATION",
                        "$.overall_verdict",
                        "failed zero-declaration check requires FAIL",
                    )
                )
    elif zero_check.get("applicability") != "NOT_APPLICABLE":
        result.append(
            problem(
                "VER_ZERO_DECLARATION",
                "$.zero_declaration_check",
                "only REAL+PRESENT uses the zero-declaration check",
            )
        )

    negative_check = verification.get("negative_control_integrity", {})
    if present:
        if negative_check.get("applicability") != "APPLICABLE":
            result.append(
                problem(
                    "VER_NEGATIVE_CONTROL",
                    "$.negative_control_integrity",
                    "PRESENT requires an applicable negative-control check",
                )
            )
        else:
            for index, ref in enumerate(items(negative_check.get("binding_refs"))):
                if ref not in binding_by_id:
                    result.append(
                        problem(
                            "VER_REPORT_REFERENCE",
                            f"$.negative_control_integrity.binding_refs.items[{index}]",
                            "reference is not a report binding item_id",
                        )
                    )
            if (
                negative_check.get("result") == "FAIL"
                and verification.get("overall_verdict") != "FAIL"
            ):
                result.append(
                    problem(
                        "VER_NEGATIVE_CONTROL",
                        "$.overall_verdict",
                        "failed negative-control check requires FAIL",
                    )
                )
    elif negative_check.get("applicability") != "NOT_APPLICABLE":
        result.append(
            problem(
                "VER_NEGATIVE_CONTROL",
                "$.negative_control_integrity",
                "ABSENT machinery requires NOT_APPLICABLE",
            )
        )

    gap_check = verification.get("gap_credibility", {})
    if gap_check.get("result") == "FAIL" and verification.get("overall_verdict") != "FAIL":
        result.append(
            problem(
                "VER_GAP_CREDIBILITY",
                "$.overall_verdict",
                "failed gap-credibility check requires FAIL",
            )
        )

    return result


def validate(
    schema_path: Path,
    document_path: Path,
) -> list[Problem]:
    try:
        schema = yaml_load(schema_path)
    except (OSError, yaml.YAMLError) as exc:
        return [problem("VALIDATOR_SCHEMA_IO", "$", str(exc))]
    try:
        document = load_effective_document(document_path)
    except (OSError, yaml.YAMLError, ValueError) as exc:
        return [problem("VALIDATOR_DOCUMENT_IO", "$", str(exc))]
    if not isinstance(schema, Mapping):
        return [problem("VALIDATOR_SCHEMA_IO", "$", "schema must be a mapping")]
    schema_name = schema_path.name
    if schema_name == REPORT_SCHEMA_NAME:
        preflight = report_preflight(document, document_path)
        if preflight:
            return preflight
        structural = schema_problems(schema, document)
        if structural:
            return structural
        if not isinstance(document, Mapping):
            return [problem("SCHEMA_STRUCTURE", "$", "report must be a mapping")]
        return report_semantic_problems(document)
    if schema_name == VERIFICATION_SCHEMA_NAME:
        preflight = verification_preflight(document)
        if preflight:
            return preflight
        structural = schema_problems(schema, document)
        if structural:
            return structural
        if not isinstance(document, Mapping):
            return [problem("SCHEMA_STRUCTURE", "$", "verification must be a mapping")]
        report_schema = yaml_load(schema_path.with_name(REPORT_SCHEMA_NAME))
        return verification_semantic_problems(
            document,
            document_path,
            report_schema,
            schema,
        )
    return [
        problem(
            "VALIDATOR_UNKNOWN_SCHEMA",
            "$",
            f"schema filename must be {REPORT_SCHEMA_NAME} or {VERIFICATION_SCHEMA_NAME}",
        )
    ]


def run_examples(expectations_path: Path) -> int:
    try:
        manifest = yaml_load(expectations_path)
    except (OSError, yaml.YAMLError) as exc:
        print(f"VALIDATOR_EXAMPLE_MANIFEST $: {exc}", file=sys.stderr)
        return 2
    if not isinstance(manifest, Mapping) or not isinstance(manifest.get("cases"), list):
        print("VALIDATOR_EXAMPLE_MANIFEST $: expected cases list", file=sys.stderr)
        return 2
    failures = 0
    for case in manifest["cases"]:
        name = case.get("name", "<unnamed>")
        schema_path = SCHEMA_DIR / case["schema"]
        document_path = SCHEMA_DIR / case["document"]
        problems = validate(schema_path, document_path)
        actual_codes = {entry.code for entry in problems}
        if case.get("expect") == "ACCEPT":
            passed = not problems
            if passed:
                print(f"PASS {name}: accepted")
            else:
                failures += 1
                print(f"FAIL {name}: expected acceptance")
                for entry in problems:
                    print(f"  {entry.render()}")
        else:
            expected_codes = set(case.get("expected_error_codes", []))
            passed = bool(problems) and actual_codes == expected_codes
            if passed:
                print(
                    f"PASS {name}: rejected for asserted rule(s) "
                    + ", ".join(sorted(expected_codes))
                )
            else:
                failures += 1
                print(
                    f"FAIL {name}: expected exact error code set "
                    f"{sorted(expected_codes)}, got {sorted(actual_codes)}"
                )
                for entry in problems:
                    print(f"  {entry.render()}")
    if failures:
        print(f"EXAMPLE GATE FAIL: {failures} case(s)", file=sys.stderr)
        return 1
    print(f"EXAMPLE GATE PASS: {len(manifest['cases'])} case(s)")
    return 0


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--schema", type=Path, help="report or verification schema YAML")
    parser.add_argument("--document", type=Path, help="YAML document to validate")
    parser.add_argument(
        "--test-examples",
        type=Path,
        nargs="?",
        const=SCHEMA_DIR / "example_expectations.yaml",
        help="run the committed example gate (default manifest when no path is supplied)",
    )
    args = parser.parse_args(argv)
    if args.test_examples is None and (args.schema is None or args.document is None):
        parser.error("--schema and --document are required unless --test-examples is used")
    if args.test_examples is not None and (args.schema is not None or args.document is not None):
        parser.error("--test-examples cannot be combined with --schema/--document")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    if args.test_examples is not None:
        return run_examples(args.test_examples)
    problems = validate(args.schema.resolve(), args.document.resolve())
    if problems:
        for entry in problems:
            print(entry.render(), file=sys.stderr)
        return 1
    print(f"VALID {args.document} against {args.schema}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
