#!/usr/bin/env python3
"""Git-object-backed synchronization for the repository research memory.

The deterministic layer never summarizes research.  It resolves configured
units against one committed-tree enumeration, stages immutable inputs for an
extractor, validates generated pages, and publishes them with state written
last.  Run ``python3 memory/tools/memory.py --help`` for the CLI.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import errno
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from typing import Any, Iterable, Mapping, Sequence

import yaml

try:
    from .review_contract import (
        GROK_REVIEW_MODEL, REVIEW_CONTRACT_SHA256, REVIEW_PROMPT_SHA256, REVIEW_SCHEMA_SHA256,
    )
except ImportError:  # pragma: no cover - direct CLI execution
    from review_contract import (  # type: ignore[no-redef]
        GROK_REVIEW_MODEL, REVIEW_CONTRACT_SHA256, REVIEW_PROMPT_SHA256, REVIEW_SCHEMA_SHA256,
    )


STATE_VERSION = 1
UNIT_DIGEST_VERSION = 1
INPUT_FORMAT_VERSION = 1
DEFAULT_CONFIG = "memory/_meta/config.yaml"
DEFAULT_SCHEMA = "memory/_meta/schema.md"
DEFAULT_STATE = "memory/_meta/state.json"
RUNTIME_DIRS = ("transactions", "locks", "journals", "backups")
PAGE_ROOTS = ("memory/index.md", "memory/conflicts.md", "memory/sources", "memory/topics", "memory/scripts")
LIFECYCLES = {"current", "superseded", "deleted", "retired"}
REVIEWS = {"ai_draft", "reviewed"}
PAGE_TYPES = {"source_capsule", "topic", "script_catalog", "conflict_register", "index"}
OWNERS = {"ai_generated", "human_curated", "mixed"}
EVIDENCE = {"measured", "derived", "provisional", "open", "disputed"}
READ_MODES = {"semantic", "identity_only", "excerpt"}
HARD_EXCLUDED_SOURCE_PREFIXES = (
    "research/pde_ledger",
    "research/pde_ledger_v2",
    "research/pde_ledger_v3",
)
ID_RE = re.compile(r"^[a-z0-9][a-z0-9-]*$")
STATUS_RE = re.compile(
    r"Status:\s*`lifecycle=([^`]+)`\s*·\s*`evidence=([^`]+)`\s*·\s*`memory_review=([^`]+)`"
)
CITATION_RE = re.compile(
    r"^\s*-\s+(?:Source:\s*)?`((?:research|software)/[^`]+)`\s+—\s+(.+?)\s*$",
    re.MULTILINE,
)
LINK_RE = re.compile(r"(?<!!)\[[^\]]*\]\(([^)]+)\)")
BEGIN_RE = re.compile(r"^<!-- BEGIN GENERATED:([a-z0-9][a-z0-9-]*) -->$")
END_RE = re.compile(r"^<!-- END GENERATED:([a-z0-9][a-z0-9-]*) -->$")


class MemoryErrorBase(RuntimeError):
    """Expected user-facing error."""


class DuplicateKeyError(MemoryErrorBase):
    pass


class ConfigError(MemoryErrorBase):
    pass


class GitError(MemoryErrorBase):
    pass


class LintFailure(MemoryErrorBase):
    def __init__(self, errors: Sequence[str], warnings: Sequence[str] = ()) -> None:
        self.errors = list(errors)
        self.warnings = list(warnings)
        super().__init__("lint failed:\n" + "\n".join(f"- {x}" for x in self.errors))


class UniqueKeyLoader(yaml.SafeLoader):
    """SafeLoader that refuses YAML's otherwise silent duplicate keys."""


def _construct_unique_mapping(loader: UniqueKeyLoader, node: yaml.MappingNode, deep: bool = False) -> dict[Any, Any]:
    loader.flatten_mapping(node)
    result: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            duplicate = key in result
        except TypeError as exc:
            raise DuplicateKeyError(f"unhashable YAML mapping key at line {key_node.start_mark.line + 1}") from exc
        if duplicate:
            raise DuplicateKeyError(f"duplicate YAML key {key!r} at line {key_node.start_mark.line + 1}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def load_yaml_bytes(data: bytes, label: str) -> Any:
    try:
        return yaml.load(data.decode("utf-8"), Loader=UniqueKeyLoader)
    except (UnicodeDecodeError, yaml.YAMLError, DuplicateKeyError) as exc:
        raise ConfigError(f"invalid YAML in {label}: {exc}") from exc


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, data: bytes, mode: int | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temp_path = Path(temporary)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        if mode is not None:
            os.chmod(temp_path, mode)
        os.replace(temp_path, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def write_json(path: Path, value: Any) -> None:
    atomic_write(path, canonical_json(value), 0o644)


def run_git(repo: Path, args: Sequence[str], *, input_data: bytes | None = None, check: bool = True) -> bytes:
    proc = subprocess.run(
        ["git", "-C", str(repo), *args],
        input=input_data,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if check and proc.returncode:
        detail = proc.stderr.decode("utf-8", "replace").strip()
        raise GitError(f"git {' '.join(args)} failed: {detail}")
    return proc.stdout


def find_repo(start: Path | None = None) -> Path:
    base = (start or Path.cwd()).resolve()
    output = run_git(base, ["rev-parse", "--show-toplevel"])
    return Path(output.decode("utf-8").strip()).resolve()


def resolve_commit(repo: Path, ref: str) -> str:
    return run_git(repo, ["rev-parse", "--verify", f"{ref}^{{commit}}"]).decode("ascii").strip()


@dataclasses.dataclass(frozen=True)
class TreeEntry:
    path: str
    mode: str
    object_type: str
    oid: str
    size: int | None

    def identity(self, role: str, read_mode: str, ownership: str) -> dict[str, Any]:
        return {
            "path": self.path,
            "role": role,
            "read_mode": read_mode,
            "ownership": ownership,
            "mode": self.mode,
            "object_type": self.object_type,
            "blob_oid": self.oid,
            "blob_size": self.size,
        }


def enumerate_tree(repo: Path, commit: str, roots: Sequence[str]) -> dict[str, TreeEntry]:
    """Enumerate candidate roots exactly once, preserving NUL-safe paths."""
    raw = run_git(repo, ["ls-tree", "-rz", "-l", "--full-tree", commit, "--", *roots])
    entries: dict[str, TreeEntry] = {}
    for record in raw.split(b"\0"):
        if not record:
            continue
        try:
            metadata, raw_path = record.split(b"\t", 1)
            mode_b, type_b, oid_b, size_b = metadata.split(None, 3)
            try:
                path = raw_path.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise GitError(
                    "committed candidate path is not UTF-8 and cannot be represented safely: "
                    + raw_path.hex()
                ) from exc
            size_text = size_b.strip()
            size = None if size_text == b"-" else int(size_text)
        except (ValueError, UnicodeError) as exc:
            raise GitError(f"could not parse ls-tree record {record[:160]!r}") from exc
        if path in entries:
            raise GitError(f"duplicate path from committed tree: {path}")
        entries[path] = TreeEntry(path, mode_b.decode("ascii"), type_b.decode("ascii"), oid_b.decode("ascii"), size)
    return entries


def read_blob(repo: Path, oid: str) -> bytes:
    return run_git(repo, ["cat-file", "blob", oid])


def clean_repo_path(value: Any, label: str, *, allow_memory: bool = False) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigError(f"{label} must be a non-empty string")
    if "\\" in value or value.startswith("/") or value.endswith("/"):
        raise ConfigError(f"{label} is not a normalized repository-relative POSIX path: {value!r}")
    pure = PurePosixPath(value)
    if any(part in ("", ".", "..") for part in pure.parts) or str(pure) != value:
        raise ConfigError(f"{label} is not normalized: {value!r}")
    allowed = ("research/", "software/") + (("memory/",) if allow_memory else ())
    if not value.startswith(allowed):
        raise ConfigError(f"{label} is outside the allowed roots: {value!r}")
    return value


def segment_prefix(path: str, prefix: str) -> bool:
    return path == prefix or path.startswith(prefix + "/")


def _hard_excluded_source(path: str) -> bool:
    return any(segment_prefix(path, prefix) for prefix in HARD_EXCLUDED_SOURCE_PREFIXES)


@dataclasses.dataclass
class ResolvedUnit:
    unit: dict[str, Any]
    members: list[dict[str, Any]]
    missing_required: list[str]
    semantic_bytes: int
    digest: str

    @property
    def id(self) -> str:
        return self.unit["id"]


def _excluded_by_defaults(path: str, config: Mapping[str, Any]) -> bool:
    if _hard_excluded_source(path):
        return True
    rules = config.get("semantic_excludes") or {}
    if any(segment_prefix(path, str(p).rstrip("/")) for p in rules.get("prefixes", [])):
        return True
    parts = PurePosixPath(path).parts
    if any(segment in parts for segment in rules.get("path_segments", [])):
        return True
    name = parts[-1]
    if any(name.endswith(str(s)) for s in rules.get("basename_suffixes", [])):
        return True
    if any(path.endswith(str(s)) for s in rules.get("suffixes", [])):
        return True
    return False


def _selector_matches(path: str, selector: Mapping[str, Any]) -> bool:
    prefix = selector["prefix"]
    if not path.startswith(prefix + "/"):
        return False
    relative = path[len(prefix) + 1 :]
    if not selector["recursive"] and "/" in relative:
        return False
    if path in selector.get("exclude_paths", []):
        return False
    if any(segment_prefix(path, p) for p in selector.get("exclude_prefixes", [])):
        return False
    basename = PurePosixPath(path).name
    suffixes = selector.get("suffixes")
    if suffixes and not any(path.endswith(suffix) for suffix in suffixes):
        return False
    basenames = selector.get("basenames")
    if basenames and basename not in basenames:
        return False
    prefixes = selector.get("name_prefixes")
    if prefixes and not any(basename.startswith(item) for item in prefixes):
        return False
    return True


def validate_config(raw: Any) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise ConfigError("config must be a YAML mapping")
    config = raw
    limits = config.get("read_limits", {})
    if not isinstance(limits, dict):
        raise ConfigError("read_limits must be a mapping")
    for key in (
        "semantic_member_max_bytes", "semantic_unit_max_bytes", "derived_task_max_bytes",
        "excerpt_context_bytes", "prepare_batch_max_units", "prepare_batch_max_tasks",
        "prepare_batch_max_semantic_bytes", "review_prompt_max_bytes",
    ):
        if key in limits and (
            not isinstance(limits[key], int) or isinstance(limits[key], bool) or limits[key] <= 0
        ):
            raise ConfigError(f"read_limits.{key} must be a positive integer")
    roots = config.get("discovery", {}).get("candidate_roots")
    if not isinstance(roots, list) or not roots:
        raise ConfigError("discovery.candidate_roots must be a non-empty list")
    for root in roots:
        if root not in ("research", "software"):
            raise ConfigError(f"unsupported candidate root {root!r}")
    units = config.get("units")
    if not isinstance(units, list):
        raise ConfigError("units must be a list")
    unit_ids: set[str] = set()
    capsules: set[str] = set()
    allowed_selector = set(config.get("selector_contract", {}).get("allowed_keys", []))
    for position, unit in enumerate(units):
        label = f"units[{position}]"
        if not isinstance(unit, dict):
            raise ConfigError(f"{label} must be a mapping")
        for field in ("id", "kind", "capsule", "lifecycle"):
            if field not in unit:
                raise ConfigError(f"{label} is missing {field}")
        unit_id = unit["id"]
        if not isinstance(unit_id, str) or not ID_RE.fullmatch(unit_id):
            raise ConfigError(f"{label}.id must be stable kebab-case")
        if unit_id in unit_ids:
            raise ConfigError(f"duplicate unit id: {unit_id}")
        unit_ids.add(unit_id)
        capsule = clean_repo_path(unit["capsule"], f"{label}.capsule", allow_memory=True)
        if not capsule.endswith(".md"):
            raise ConfigError(f"{label}.capsule must be Markdown")
        if capsule in capsules:
            raise ConfigError(f"duplicate capsule target: {capsule}")
        capsules.add(capsule)
        if unit["lifecycle"] not in LIFECYCLES:
            raise ConfigError(f"{label}.lifecycle has invalid value")
        editorial_focus = unit.get("editorial_focus", [])
        if (
            not isinstance(editorial_focus, list)
            or not all(isinstance(item, str) and item.strip() for item in editorial_focus)
        ):
            raise ConfigError(f"{label}.editorial_focus must be a list of non-empty strings")
        seen_exact: set[str] = set()
        for mpos, member in enumerate(unit.get("members", [])):
            mlabel = f"{label}.members[{mpos}]"
            if not isinstance(member, dict):
                raise ConfigError(f"{mlabel} must be a mapping")
            for field in ("path", "role", "read_mode", "required"):
                if field not in member:
                    raise ConfigError(f"{mlabel} is missing {field}")
            path = clean_repo_path(member["path"], f"{mlabel}.path")
            if _hard_excluded_source(path):
                raise ConfigError(f"{mlabel}.path is in a hard-excluded source root: {path}")
            if path in seen_exact:
                raise ConfigError(f"duplicate exact member in {unit_id}: {path}")
            seen_exact.add(path)
            if member["read_mode"] not in READ_MODES:
                raise ConfigError(f"{mlabel}.read_mode is invalid")
            if not isinstance(member["required"], bool):
                raise ConfigError(f"{mlabel}.required must be boolean")
            ownership = member.get("ownership", "primary")
            if ownership not in ("primary", "shared_dependency"):
                raise ConfigError(f"{mlabel}.ownership is invalid")
            if member["read_mode"] == "excerpt" and not member.get("anchors"):
                raise ConfigError(f"{mlabel} excerpt mode requires anchors")
        for spos, selector in enumerate(unit.get("selectors", [])):
            slabel = f"{label}.selectors[{spos}]"
            if not isinstance(selector, dict):
                raise ConfigError(f"{slabel} must be a mapping")
            unknown = set(selector) - allowed_selector
            if unknown:
                raise ConfigError(f"{slabel} has unsupported keys: {sorted(unknown)}")
            for field in ("prefix", "recursive", "role", "read_mode", "ownership", "required"):
                if field not in selector:
                    raise ConfigError(f"{slabel} is missing {field}")
            clean_repo_path(selector["prefix"] + "/x", f"{slabel}.prefix")
            if _hard_excluded_source(selector["prefix"]):
                raise ConfigError(
                    f"{slabel}.prefix is in a hard-excluded source root: {selector['prefix']}"
                )
            if not isinstance(selector["recursive"], bool) or not isinstance(selector["required"], bool):
                raise ConfigError(f"{slabel} recursive/required must be boolean")
            if selector["read_mode"] not in READ_MODES:
                raise ConfigError(f"{slabel}.read_mode is invalid")
            if selector["ownership"] not in ("primary", "shared_dependency"):
                raise ConfigError(f"{slabel}.ownership is invalid")
            for key in ("suffixes", "basenames", "name_prefixes"):
                if key in selector and (not isinstance(selector[key], list) or not all(isinstance(x, str) for x in selector[key])):
                    raise ConfigError(f"{slabel}.{key} must be a string list")
            for key in ("exclude_paths",):
                for value in selector.get(key, []):
                    clean_repo_path(value, f"{slabel}.{key}")
            for key in ("exclude_prefixes",):
                for value in selector.get(key, []):
                    clean_repo_path(value + "/x", f"{slabel}.{key}")
        if not unit.get("members") and not unit.get("selectors"):
            raise ConfigError(f"{label} has no members or selectors")
    supporting_lineages = config.get("supporting_lineages", [])
    if not isinstance(supporting_lineages, list):
        raise ConfigError("supporting_lineages must be a list")
    for position, lineage in enumerate(supporting_lineages):
        label = f"supporting_lineages[{position}]"
        if not isinstance(lineage, dict) or not isinstance(lineage.get("root"), str):
            raise ConfigError(f"{label} must be a mapping with a root")
        clean_repo_path(lineage["root"] + "/x", f"{label}.root")
        if _hard_excluded_source(lineage["root"]):
            raise ConfigError(f"{label}.root is in a hard-excluded source root: {lineage['root']}")
    derived_pages = config.get("derived_pages", [])
    if not isinstance(derived_pages, list):
        raise ConfigError("derived_pages must be a list")
    derived_ids: set[str] = set()
    derived_paths: set[str] = set()
    derived_by_id: dict[str, dict[str, Any]] = {}
    allowed_derived_kinds = {"topic", "script_catalog", "conflict_register", "index"}
    budgets = config.get("output_budgets", {})
    if not isinstance(budgets, dict):
        raise ConfigError("output_budgets must be a mapping")
    budget_shapes = {
        "source_capsule": {"max_words", "max_key_statements"},
        "topic": {"max_words", "max_key_statements"},
        "script_catalog": {"max_words", "max_entries", "max_entry_words"},
        "conflict_register": {"max_words", "max_entries", "max_entry_words"},
        "index": {"max_words"},
    }
    unknown_budgets = set(budgets) - set(budget_shapes)
    if unknown_budgets:
        raise ConfigError(f"output_budgets has unsupported entries: {sorted(unknown_budgets)}")
    for budget_name, values in budgets.items():
        if not isinstance(values, dict) or not values or set(values) - budget_shapes[budget_name]:
            raise ConfigError(
                f"output_budgets.{budget_name} may contain only {sorted(budget_shapes[budget_name])}"
            )
        for key, value in values.items():
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise ConfigError(f"output_budgets.{budget_name}.{key} must be a positive integer")
    for position, item in enumerate(derived_pages):
        label = f"derived_pages[{position}]"
        if not isinstance(item, dict):
            raise ConfigError(f"{label} must be a mapping")
        expected_keys = {"id", "task_kind", "page", "region", "order", "input_units", "input_pages", "budget"}
        allowed_keys = expected_keys | {"direct_sources"}
        if not expected_keys.issubset(item) or set(item) - allowed_keys:
            raise ConfigError(
                f"{label} must contain {sorted(expected_keys)} and may additionally contain direct_sources"
            )
        derived_id = item["id"]
        if not isinstance(derived_id, str) or not ID_RE.fullmatch(derived_id):
            raise ConfigError(f"{label}.id must be stable kebab-case")
        if derived_id in derived_ids or derived_id in unit_ids:
            raise ConfigError(f"duplicate/colliding derived task id: {derived_id}")
        derived_ids.add(derived_id)
        derived_by_id[derived_id] = item
        if item["task_kind"] not in allowed_derived_kinds:
            raise ConfigError(f"{label}.task_kind is invalid")
        page = clean_repo_path(item["page"], f"{label}.page", allow_memory=True)
        if not page.endswith(".md") or page in derived_paths or page in capsules:
            raise ConfigError(f"duplicate or invalid derived page target: {page}")
        derived_paths.add(page)
        if not isinstance(item["region"], str) or not ID_RE.fullmatch(item["region"]):
            raise ConfigError(f"{label}.region must be stable kebab-case")
        if not isinstance(item["order"], int) or item["order"] < 1:
            raise ConfigError(f"{label}.order must be a positive integer")
        if not isinstance(item["input_units"], list) or not all(x in unit_ids for x in item["input_units"]):
            raise ConfigError(f"{label}.input_units contains an unknown unit")
        if len(item["input_units"]) != len(set(item["input_units"])):
            raise ConfigError(f"{label}.input_units contains duplicates")
        if not isinstance(item["input_pages"], list) or not all(isinstance(x, str) for x in item["input_pages"]):
            raise ConfigError(f"{label}.input_pages must be an ID list")
        if item["budget"] not in budgets:
            raise ConfigError(f"{label}.budget does not name output_budgets entry")
        direct_sources = item.get("direct_sources")
        if direct_sources is not None:
            if not isinstance(direct_sources, list):
                raise ConfigError(f"{label}.direct_sources must be a list")
            direct_paths: set[str] = set()
            for dpos, direct in enumerate(direct_sources):
                dlabel = f"{label}.direct_sources[{dpos}]"
                if not isinstance(direct, dict) or set(direct) - {"path", "read_mode", "anchors"}:
                    raise ConfigError(f"{dlabel} has invalid shape")
                if not {"path", "read_mode"}.issubset(direct):
                    raise ConfigError(f"{dlabel} requires path and read_mode")
                path = clean_repo_path(direct["path"], f"{dlabel}.path")
                if _hard_excluded_source(path):
                    raise ConfigError(f"{dlabel}.path is in a hard-excluded source root: {path}")
                if path in direct_paths:
                    raise ConfigError(f"{label}.direct_sources repeats {path}")
                direct_paths.add(path)
                if direct["read_mode"] not in ("semantic", "excerpt"):
                    raise ConfigError(f"{dlabel}.read_mode must be semantic or excerpt")
                if direct["read_mode"] == "excerpt" and not direct.get("anchors"):
                    raise ConfigError(f"{dlabel} excerpt mode requires anchors")
    for derived_id, item in derived_by_id.items():
        for dependency in item["input_pages"]:
            if dependency not in derived_by_id:
                raise ConfigError(f"derived task {derived_id} refers to unknown input page {dependency}")
            if derived_by_id[dependency]["order"] >= item["order"]:
                raise ConfigError(f"derived task {derived_id} dependency {dependency} must have lower order")
    return config


def load_config(repo: Path, config_path: str = DEFAULT_CONFIG) -> tuple[dict[str, Any], bytes]:
    path = repo / config_path
    try:
        data = path.read_bytes()
    except OSError as exc:
        raise ConfigError(f"cannot read config {config_path}: {exc}") from exc
    return validate_config(load_yaml_bytes(data, config_path)), data


def resolve_units(config: Mapping[str, Any], entries: Mapping[str, TreeEntry]) -> dict[str, ResolvedUnit]:
    limits = config.get("read_limits", {})
    member_limit = int(limits.get("semantic_member_max_bytes", 2_000_000))
    unit_limit = int(limits.get("semantic_unit_max_bytes", 8_000_000))
    resolved: dict[str, ResolvedUnit] = {}
    primary_owners: dict[str, str] = {}
    all_claims: dict[tuple[str, str], tuple[str, str, str]] = {}
    for unit in config["units"]:
        members: dict[str, dict[str, Any]] = {}
        missing: list[str] = []

        def add(path: str, role: str, read_mode: str, ownership: str, required: bool, origin: str, spec: Mapping[str, Any]) -> None:
            entry = entries.get(path)
            if entry is None:
                if required:
                    missing.append(path if origin == "exact" else f"selector:{origin}")
                return
            if entry.object_type != "blob" or entry.mode not in ("100644", "100755"):
                if read_mode != "identity_only":
                    raise ConfigError(
                        f"{unit['id']} member {path} has unsupported {entry.mode}/{entry.object_type}; "
                        "declare identity_only to track identity without reading it"
                    )
            identity = entry.identity(role, read_mode, ownership)
            if "anchors" in spec:
                identity["anchors"] = spec["anchors"]
            if "context_bytes" in spec:
                identity["context_bytes"] = int(spec["context_bytes"])
            if "max_bytes" in spec:
                identity["max_bytes"] = int(spec["max_bytes"])
            previous = members.get(path)
            if previous is not None and previous != identity:
                raise ConfigError(f"conflicting duplicate membership for {path} inside {unit['id']}")
            members[path] = identity

        for member in unit.get("members", []):
            add(
                member["path"], member["role"], member["read_mode"], member.get("ownership", "primary"),
                member["required"], "exact", member,
            )
        for selector in unit.get("selectors", []):
            matched = []
            for path in entries:
                if _excluded_by_defaults(path, config):
                    continue
                if _selector_matches(path, selector):
                    matched.append(path)
            if selector["required"] and not matched:
                missing.append(f"selector:{selector['prefix']}")
            for path in matched:
                add(
                    path, selector["role"], selector["read_mode"], selector["ownership"], False,
                    selector["prefix"], selector,
                )
        member_list = [members[path] for path in sorted(members)]
        semantic_bytes = 0
        for item in member_list:
            if item["read_mode"] == "semantic":
                size = item["blob_size"] or 0
                effective = int(item.get("max_bytes", member_limit))
                if effective > member_limit:
                    raise ConfigError(f"{unit['id']} member {item['path']} max_bytes exceeds global limit")
                if size > effective:
                    raise ConfigError(
                        f"semantic member {item['path']} is {size} bytes, above its {effective}-byte limit; "
                        "use identity_only or a deterministic excerpt"
                    )
                semantic_bytes += size
        if semantic_bytes > unit_limit:
            raise ConfigError(f"unit {unit['id']} has {semantic_bytes} semantic bytes, above {unit_limit}")
        digest_payload = {
            "version": UNIT_DIGEST_VERSION,
            "id": unit["id"],
            "semantic_config": unit,
            "members": member_list,
            "missing_required": sorted(set(missing)),
        }
        digest = sha256_bytes(canonical_json(digest_payload))
        resolved[unit["id"]] = ResolvedUnit(unit, member_list, sorted(set(missing)), semantic_bytes, digest)
        for item in member_list:
            key = (unit["id"], item["path"])
            all_claims[key] = (item["ownership"], item["role"], item["read_mode"])
            if item["ownership"] == "primary":
                old = primary_owners.get(item["path"])
                if old is not None and old != unit["id"]:
                    raise ConfigError(f"primary ownership overlap for {item['path']}: {old}, {unit['id']}")
                primary_owners[item["path"]] = unit["id"]
    return resolved


def digest_directory(path: Path) -> str:
    records: list[dict[str, str]] = []
    if path.exists():
        for item in sorted(p for p in path.rglob("*") if p.is_file()):
            records.append({"path": item.relative_to(path).as_posix(), "sha256": sha256_file(item)})
    return sha256_bytes(canonical_json({"files": records}))


def tooling_digest() -> str:
    tool_root = Path(__file__).resolve().parent
    records: list[dict[str, str]] = []
    for name in ("memory.py", "run_isolated.py", "review_isolated.py", "review_contract.py"):
        path = tool_root / name
        if not path.is_file():
            raise ConfigError(f"missing deterministic tool component: {path}")
        records.append({"path": name, "sha256": sha256_file(path)})
    return sha256_bytes(canonical_json({"files": records}))


def policy_digests(repo: Path, config_bytes: bytes) -> dict[str, str]:
    schema_path = repo / DEFAULT_SCHEMA
    if not schema_path.is_file():
        raise ConfigError(f"missing policy file {DEFAULT_SCHEMA}")
    pieces = {
        "config_sha256": sha256_bytes(config_bytes),
        "schema_sha256": sha256_file(schema_path),
        "prompts_sha256": digest_directory(repo / "memory/prompts"),
        "tool_sha256": tooling_digest(),
        "contract_sha256": sha256_file(repo / "memory/_meta/SCANNER_DECISIONS.md"),
        "atlas_migration_sha256": sha256_file(repo / "memory/_meta/atlas-migration.yaml"),
    }
    pieces["combined_sha256"] = sha256_bytes(canonical_json(pieces))
    return pieces


def migration_requirements(repo: Path) -> dict[str, list[dict[str, Any]]]:
    path = repo / "memory/_meta/atlas-migration.yaml"
    document = load_yaml_bytes(path.read_bytes(), path.relative_to(repo).as_posix())
    if not isinstance(document, dict):
        raise ConfigError("atlas migration inventory must be a mapping")
    by_target: dict[str, list[dict[str, Any]]] = {}

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            if value.get("disposition") == "migrate" and isinstance(value.get("id"), str):
                target = value.get("target")
                if isinstance(target, str) and target.startswith("memory/") and target.endswith(".md"):
                    originals = value.get("original_sources", [])
                    if not isinstance(originals, list) or not all(isinstance(item, str) for item in originals):
                        raise ConfigError(f"migration item {value['id']} has invalid original_sources")
                    for reference in originals:
                        source_path = reference.partition("#")[0]
                        clean_repo_path(source_path, f"migration item {value['id']} original source")
                        if _hard_excluded_source(source_path):
                            raise ConfigError(
                                f"migration item {value['id']} uses hard-excluded source {source_path}"
                            )
                    by_target.setdefault(target, []).append({
                        "legacy_id": value["id"],
                        "required_original_refs": originals,
                        "reanchor_complete": value.get("reanchor_complete"),
                    })
            for nested in value.values():
                visit(nested)
        elif isinstance(value, list):
            for nested in value:
                visit(nested)

    visit(document)
    for target in by_target:
        by_target[target].sort(key=lambda item: item["legacy_id"])
    return by_target


def empty_state() -> dict[str, Any]:
    return {
        "state_version": STATE_VERSION,
        "generation": None,
        "last_fully_processed_commit": None,
        "policy": {},
        "units": {},
        "retired_units": {},
        "pages": {},
        "derived_pages": {},
        "id_registry": {"pages": {}, "statements": {}, "conflicts": {}},
        "last_result": {"status": "uninitialized", "error": None},
    }


def load_state(repo: Path) -> tuple[dict[str, Any], bytes | None]:
    path = repo / DEFAULT_STATE
    if not path.exists():
        return empty_state(), None
    data = path.read_bytes()
    try:
        state = json.loads(data)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise MemoryErrorBase(f"invalid state JSON: {exc}") from exc
    if not isinstance(state, dict) or state.get("state_version") != STATE_VERSION:
        raise MemoryErrorBase("unsupported or invalid memory state")
    for field in ("units", "retired_units", "pages", "derived_pages"):
        if not isinstance(state.get(field), dict):
            raise MemoryErrorBase(f"state field {field} must be a mapping")
    return state, data


def state_digest(data: bytes | None) -> str:
    return sha256_bytes(data if data is not None else b"")


def dirty_tracked_paths(repo: Path, roots: Sequence[str]) -> list[str]:
    raw = run_git(repo, ["status", "--porcelain=v1", "-z", "--untracked-files=no", "--", *roots])
    records = raw.split(b"\0")
    result: list[str] = []
    index = 0
    while index < len(records):
        record = records[index]
        index += 1
        if not record:
            continue
        if len(record) < 4:
            continue
        status = record[:2]
        result.append(record[3:].decode("utf-8", "surrogateescape"))
        if b"R" in status or b"C" in status:
            if index < len(records) and records[index]:
                result.append(records[index].decode("utf-8", "surrogateescape"))
                index += 1
    return sorted(set(result))


def committed_renames(repo: Path, old_commit: str | None, target_commit: str, roots: Sequence[str]) -> list[dict[str, str]]:
    if not old_commit or old_commit == target_commit:
        return []
    raw = run_git(
        repo, ["diff", "--name-status", "-z", "--find-renames", old_commit, target_commit, "--", *roots]
    )
    fields = raw.split(b"\0")
    result: list[dict[str, str]] = []
    index = 0
    while index < len(fields):
        status = fields[index].decode("ascii", "replace")
        index += 1
        if not status:
            break
        if status.startswith("R") and index + 1 < len(fields):
            old_path = fields[index].decode("utf-8", "surrogateescape")
            new_path = fields[index + 1].decode("utf-8", "surrogateescape")
            index += 2
            result.append({"similarity": status[1:] or "unknown", "old_path": old_path, "new_path": new_path})
        else:
            index += 1
    return result


def compute_actions(
    config: Mapping[str, Any], resolved: Mapping[str, ResolvedUnit], state: Mapping[str, Any],
    target_commit: str, policy: Mapping[str, str],
) -> dict[str, dict[str, Any]]:
    actions: dict[str, dict[str, Any]] = {}
    old_units = state.get("units", {})
    policy_changed = state.get("policy", {}).get("combined_sha256") != policy["combined_sha256"]
    for unit in config["units"]:
        unit_id = unit["id"]
        item = resolved[unit_id]
        old = old_units.get(unit_id)
        if old is None:
            action = "invalid_initial" if item.missing_required else "new"
        elif not item.members:
            if (
                old.get("result") == "source_deleted"
                and item.digest == old.get("unit_digest_sha256")
                and not policy_changed
                and old.get("policy_sha256") == policy["combined_sha256"]
            ):
                action = "current" if old.get("processed_commit") == target_commit else "rebase_unchanged"
            else:
                action = "source_deleted"
        elif item.digest != old.get("unit_digest_sha256"):
            action = "changed"
        elif policy_changed or old.get("policy_sha256") != policy["combined_sha256"]:
            action = "policy_changed"
        elif old.get("processed_commit") != target_commit:
            action = "rebase_unchanged"
        elif old.get("result") != "success":
            action = "retry"
        else:
            action = "current"
        actions[unit_id] = {
            "action": action,
            "capsule": unit["capsule"],
            "unit_digest_sha256": item.digest,
            "missing_required": item.missing_required,
            "member_count": len(item.members),
            "semantic_bytes": item.semantic_bytes,
        }
    for unit_id, old in old_units.items():
        if unit_id not in resolved:
            actions[unit_id] = {
                "action": "retired_from_corpus",
                "capsule": old.get("capsule"),
                "unit_digest_sha256": old.get("unit_digest_sha256"),
                "missing_required": [],
                "member_count": len(old.get("members", [])),
                "semantic_bytes": 0,
            }
    return actions


def migration_status_report(repo: Path, config: Mapping[str, Any], state: Mapping[str, Any]) -> dict[str, Any]:
    requirements = migration_requirements(repo)
    derived_by_page = {item["page"]: item["id"] for item in config.get("derived_pages", [])}
    required_ids: list[str] = []
    completed_ids: list[str] = []
    pending: list[dict[str, str]] = []
    blocked_ids: list[str] = []
    for page, items in sorted(requirements.items()):
        derived_id = derived_by_page.get(page)
        recorded = (
            set(state.get("derived_pages", {}).get(derived_id, {}).get("migration_ids", []))
            if derived_id is not None else set()
        )
        for item in items:
            legacy_id = item["legacy_id"]
            required_ids.append(legacy_id)
            if legacy_id in recorded:
                completed_ids.append(legacy_id)
                page_path = repo / page
                if page_path.is_file():
                    try:
                        _, page_body = parse_frontmatter(page_path.read_bytes(), page)
                    except LintFailure:
                        page_body = ""
                    record_match = re.search(
                        rf"^#### Migration record — {re.escape(legacy_id)}\s*$.*?"
                        r"^- `migration_disposition`: `(migrated|blocked)`\s*$",
                        page_body, re.MULTILINE | re.DOTALL,
                    )
                    if record_match and record_match.group(1) == "blocked":
                        blocked_ids.append(legacy_id)
            else:
                pending.append({
                    "legacy_id": legacy_id,
                    "target": page,
                    "reason": "not_completed" if derived_id is not None else "target_not_managed_by_derived_task",
                })
    return {
        "required_count": len(required_ids),
        "completed_count": len(completed_ids),
        "pending": pending,
        "migration_items_accounted": not pending,
        "blocked_ids": sorted(blocked_ids),
        "all_migrated": not pending and not blocked_ids,
    }


def status_report(repo: Path, ref: str = "HEAD", config_path: str = DEFAULT_CONFIG) -> dict[str, Any]:
    config, config_bytes = load_config(repo, config_path)
    commit = resolve_commit(repo, ref)
    roots = config["discovery"]["candidate_roots"]
    entries = enumerate_tree(repo, commit, roots)
    resolved = resolve_units(config, entries)
    state, _ = load_state(repo)
    policy = policy_digests(repo, config_bytes)
    actions = compute_actions(config, resolved, state, commit, policy)
    claimed = {member["path"] for unit in resolved.values() for member in unit.members}
    old_claimed = {member["path"] for unit in state.get("units", {}).values() for member in unit.get("members", [])}
    dirty = [path for path in dirty_tracked_paths(repo, roots) if path in claimed or path in old_claimed]
    pending = sorted(unit_id for unit_id, detail in actions.items() if detail["action"] != "current")
    derived_status: dict[str, dict[str, Any]] = {}
    for item in config.get("derived_pages", []):
        record = state.get("derived_pages", {}).get(item["id"])
        page_record = state.get("pages", {}).get(item["page"])
        action = "current"
        if record is None or page_record is None:
            action = "new"
        elif record.get("policy_sha256") != policy["combined_sha256"]:
            action = "policy_changed"
        elif record.get("processed_commit") != commit:
            action = "pending_refresh"
        derived_status[item["id"]] = {
            "action": action, "page": item["page"], "order": item["order"],
            "input_units": item["input_units"], "input_pages": item["input_pages"],
        }
    pending_derived = sorted(item_id for item_id, detail in derived_status.items() if detail["action"] != "current")
    configured_derived_ids = set(derived_status)
    orphaned_derived = sorted(set(state.get("derived_pages", {})) - configured_derived_ids)
    transactions_dir = repo / "memory/transactions"
    transactions = sorted(p.name for p in transactions_dir.iterdir() if p.is_dir()) if transactions_dir.is_dir() else []
    served_pages, served_integrity_errors = served_page_bytes(repo, state)
    served_orphans = sorted(_page_candidates(repo) - set(state.get("pages", {})))
    served_integrity_errors.extend(f"unrecorded page: {path}" for path in served_orphans)
    migration_report = migration_status_report(repo, config, state)
    no_legacy_references = not any(
        re.search(r"(?:^|[\s`(])(?:atlas|graph)/", data.decode("utf-8", "replace"))
        for data in served_pages.values()
    )
    benchmark = state.get("retrieval_benchmark", {})
    benchmark_pass = benchmark.get("status") == "pass" and benchmark.get("processed_commit") == commit
    fresh_agent_pass = benchmark.get("fresh_agent_status") == "pass"
    served_lint_errors, _ = lint_served(repo)
    fully_fresh = (
        not pending and not pending_derived and not orphaned_derived
        and state.get("last_fully_processed_commit") == commit
    )
    cutover_conditions = {
        "full_freshness": fully_fresh,
        "served_lint": not served_lint_errors,
        "served_integrity": not served_integrity_errors,
        "migration_items_accounted": migration_report["migration_items_accounted"],
        "migration_all_migrated": migration_report["all_migrated"],
        "no_legacy_references": no_legacy_references,
        "retrieval_benchmark_recorded_pass": benchmark_pass,
        "fresh_agent_retrieval_pass": fresh_agent_pass,
    }
    return {
        "target_commit": commit,
        "last_fully_processed_commit": state.get("last_fully_processed_commit"),
        "fully_fresh": fully_fresh,
        "policy": policy,
        "policy_drift": state.get("policy", {}).get("combined_sha256") not in (None, policy["combined_sha256"]),
        "units": actions,
        "pending_units": pending,
        "derived_pages": derived_status,
        "pending_derived_pages": pending_derived,
        "orphaned_derived_pages": orphaned_derived,
        "dirty_tracked_members": dirty,
        "renames": committed_renames(repo, state.get("last_fully_processed_commit"), commit, roots),
        "transactions": transactions,
        "tracked_candidate_count": len(entries),
        "migration": migration_report,
        "served_integrity": {
            "ok": not served_integrity_errors,
            "recorded_pages": len(served_pages),
            "errors": served_integrity_errors,
        },
        "cutover": {"ready": all(cutover_conditions.values()), "conditions": cutover_conditions},
    }


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError as exc:
        return exc.errno == errno.EPERM
    return True


def publication_paths(repo: Path) -> tuple[Path, Path]:
    return repo / "memory/journals/publication.json", repo / "memory/locks/publication.lock"


def acquire_lock(repo: Path) -> Path:
    _, lock = publication_paths(repo)
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise MemoryErrorBase(f"memory publication lock already exists: {lock.relative_to(repo)}") from exc
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        json.dump({"pid": os.getpid(), "created_at": dt.datetime.now(dt.timezone.utc).isoformat()}, handle)
        handle.flush()
        os.fsync(handle.fileno())
    return lock


def release_lock(lock: Path) -> None:
    try:
        lock.unlink()
    except FileNotFoundError:
        pass


def _restore_publication(repo: Path, journal: Mapping[str, Any]) -> None:
    backup_root = repo / journal["backup_root"]
    for record in reversed(journal.get("targets", [])):
        target = repo / record["path"]
        if record["existed"]:
            source = backup_root / record["backup"]
            atomic_write(target, source.read_bytes())
        elif target.exists():
            target.unlink()
    state_record = journal["state"]
    state_path = repo / DEFAULT_STATE
    if state_record["existed"]:
        atomic_write(state_path, (backup_root / state_record["backup"]).read_bytes())
    elif state_path.exists():
        state_path.unlink()


def _cleanup_publication(repo: Path, journal: Mapping[str, Any]) -> None:
    backup_root = repo / journal["backup_root"]
    if backup_root.exists():
        shutil.rmtree(backup_root)
    journal_path, _ = publication_paths(repo)
    if journal_path.exists():
        journal_path.unlink()


def recover_publication(repo: Path) -> str | None:
    journal_path, lock_path = publication_paths(repo)
    if not journal_path.exists():
        if lock_path.exists():
            try:
                owner = json.loads(lock_path.read_text(encoding="utf-8"))
                pid = int(owner.get("pid", -1))
            except Exception:
                pid = -1
            if pid <= 0 or not _pid_alive(pid):
                lock_path.unlink()
        return None
    if lock_path.exists():
        try:
            owner = json.loads(lock_path.read_text(encoding="utf-8"))
            pid = int(owner.get("pid", -1))
        except Exception:
            pid = -1
        if pid > 0 and pid != os.getpid() and _pid_alive(pid):
            raise MemoryErrorBase("another process is actively publishing memory pages")
    journal = json.loads(journal_path.read_text(encoding="utf-8"))
    state_path = repo / DEFAULT_STATE
    committed = state_path.exists() and sha256_file(state_path) == journal.get("new_state_sha256")
    if committed:
        _cleanup_publication(repo, journal)
        outcome = "completed cleanup for committed publication"
    else:
        _restore_publication(repo, journal)
        _cleanup_publication(repo, journal)
        outcome = "rolled back interrupted publication"
    release_lock(lock_path)
    return outcome


def init_memory(repo: Path) -> list[str]:
    created: list[str] = []
    for relative in ("memory/sources", "memory/topics", "memory/scripts", "memory/prompts", "memory/_meta", *[f"memory/{x}" for x in RUNTIME_DIRS]):
        path = repo / relative
        if not path.exists():
            path.mkdir(parents=True)
            created.append(relative + "/")
    state_path = repo / DEFAULT_STATE
    if not state_path.exists():
        write_json(state_path, empty_state())
        created.append(DEFAULT_STATE)
    return created


def deterministic_excerpt(data: bytes, member: Mapping[str, Any], default_context: int) -> bytes:
    if b"\0" in data:
        raise MemoryErrorBase(f"cannot excerpt binary-looking blob {member['path']}")
    anchors = member.get("anchors") or []
    context = int(member.get("context_bytes", default_context))
    if context <= 0:
        raise MemoryErrorBase(f"invalid excerpt context for {member['path']}")
    text = data.decode("utf-8", "replace")
    chunks: list[str] = []
    per_anchor = max(128, context // max(1, len(anchors)))
    for anchor in anchors:
        needle = anchor.get("value") if isinstance(anchor, dict) else str(anchor)
        if not isinstance(needle, str) or not needle:
            raise MemoryErrorBase(f"invalid excerpt anchor in {member['path']}")
        position = text.find(needle)
        if position < 0:
            raise MemoryErrorBase(f"excerpt anchor {needle!r} not found in {member['path']}")
        left = max(0, position - per_anchor // 2)
        right = min(len(text), position + len(needle) + per_anchor // 2)
        chunks.append(f"--- anchor: {needle} ---\n{text[left:right]}")
    return ("\n\n".join(chunks) + "\n").encode("utf-8")


def source_task_metadata(
    config: Mapping[str, Any], unit_item: Mapping[str, Any], resolved_unit: ResolvedUnit,
    state: Mapping[str, Any], target: str,
) -> dict[str, Any]:
    kind = resolved_unit.unit["kind"]
    kind_contract = {
        "paper": ("paper", "source-capsule-paper.md"),
        "componentized_paper": ("paper", "source-capsule-paper.md"),
        "research_governance": ("governance", "source-capsule-governance.md"),
        "step_family": ("step_family", "source-capsule-step.md"),
        "in_progress_step_family": ("step_family", "source-capsule-step.md"),
        "verification_project": ("verification", "source-capsule-software.md"),
        "software_project": ("software_project", "source-capsule-software.md"),
        "software_result_family": ("result", "source-capsule-software.md"),
    }
    if kind not in kind_contract:
        raise ConfigError(f"no semantic source-task contract for unit kind {kind!r}")
    source_kind, task_prompt = kind_contract[kind]
    if unit_item.get("action") == "source_deleted":
        task_prompt = "source-capsule-lifecycle.md"
    member_paths = [member["path"] for member in unit_item.get("members", [])]
    lifecycle_members = unit_item.get("prior_members", []) if unit_item.get("action") == "source_deleted" else []
    entrypoint = resolved_unit.unit.get("entrypoint")
    if entrypoint is None:
        preferred_roles = ("primary", "result", "step_record", "governance", "entrypoint")
        for role in preferred_roles:
            candidates = [member["path"] for member in unit_item.get("members", []) if member["role"] == role]
            if candidates:
                entrypoint = sorted(candidates)[0]
                break
    if entrypoint is None and member_paths:
        entrypoint = sorted(member_paths)[0]
    old_unit = state.get("units", {}).get(unit_item["id"])
    prior_page = state.get("pages", {}).get(unit_item["capsule"])
    prior_identity: dict[str, Any] | None = None
    if old_unit is not None:
        prior_identity = {
            "processed_commit": old_unit.get("processed_commit"),
            "unit_digest_sha256": old_unit.get("unit_digest_sha256"),
            "lifecycle": old_unit.get("lifecycle"),
            "result": old_unit.get("result"),
            "members": old_unit.get("members", []),
            "page": prior_page,
        }
        if entrypoint is None and old_unit.get("members"):
            entrypoint = sorted(member["path"] for member in old_unit["members"])[0]
    allowed: list[dict[str, Any]] = []
    current_identities = {(member["path"], member.get("blob_oid")) for member in unit_item.get("members", [])}
    for member in unit_item.get("members", []):
        allowed.append({
            "path": member["path"], "availability": "current", "read_mode": member["read_mode"],
            "blob_oid": member["blob_oid"], "last_seen_commit": target,
        })
    if old_unit is not None:
        for member in old_unit.get("members", []):
            if (member["path"], member.get("blob_oid")) not in current_identities:
                allowed.append({
                    "path": member["path"], "availability": "last_seen", "read_mode": member.get("read_mode"),
                    "blob_oid": member.get("blob_oid"), "last_seen_commit": old_unit.get("processed_commit"),
                })
    related_pages: set[str] = set()
    for lineage in config.get("canonical_lineages", []):
        if unit_item["id"] in lineage.get("retrieval_order", []):
            unit_by_id = {unit["id"]: unit for unit in config["units"]}
            related_pages.update(
                unit_by_id[related]["capsule"] for related in lineage["retrieval_order"]
                if related in unit_by_id and related != unit_item["id"]
            )
    limits = config.get("read_limits", {})
    return {
        "refresh_date": dt.datetime.now().astimezone().date().isoformat(),
        "page": {
            "id": f"source-{unit_item['id']}",
            "title": str(resolved_unit.unit.get("title") or unit_item["id"].replace("-", " ").title()),
            "type": "source_capsule",
            "content_owner": "ai_generated",
            "desired_lifecycle": "deleted" if unit_item["action"] == "source_deleted" else resolved_unit.unit["lifecycle"],
            "output_repository_path": unit_item["capsule"],
        },
        "source_kind": source_kind,
        "shape": "file" if len(unit_item.get("members", []) or lifecycle_members) == 1 else "bundle",
        "entrypoint": entrypoint,
        "extractor_version": config.get("extractor_contract_version"),
        "prompt_paths": ["/packet/prompts/00-snapshot-contract.md", f"/packet/prompts/{task_prompt}"],
        "allowed_citations": sorted(allowed, key=lambda item: item["path"]),
        "related_pages": sorted(related_pages),
        "editorial_focus": list(resolved_unit.unit.get("editorial_focus", [])),
        "budgets": {
            "semantic_member_max_bytes": int(limits.get("semantic_member_max_bytes", 2_000_000)),
            "semantic_unit_max_bytes": int(limits.get("semantic_unit_max_bytes", 8_000_000)),
            "resolved_semantic_bytes": resolved_unit.semantic_bytes,
            **dict(config.get("output_budgets", {}).get("source_capsule", {})),
        },
        "prior_identity": prior_identity,
        "prior_stable_ids": {
            category: sorted(
                stable_id for stable_id, owner_path in state.get("id_registry", {}).get(category, {}).items()
                if owner_path == unit_item["capsule"]
            )
            for category in ("pages", "statements", "conflicts")
        },
        "stable_id_namespace": f"source-{unit_item['id']}--",
        "lifecycle_action": ({
            "action": "source_deleted",
            "action_commit": target,
            "last_source_commit": old_unit.get("processed_commit") if old_unit else None,
            "deleted_member_count": len(lifecycle_members),
        } if unit_item.get("action") == "source_deleted" else None),
    }


def retired_source_task_metadata(
    repo: Path, config: Mapping[str, Any], unit_item: Mapping[str, Any], state: Mapping[str, Any], target: str,
) -> dict[str, Any]:
    old = state["units"][unit_item["id"]]
    capsule = old["capsule"]
    page_record = state.get("pages", {}).get(capsule)
    page_path = repo / capsule
    if page_record is None or not page_path.is_file() or sha256_file(page_path) != page_record.get("sha256"):
        raise MemoryErrorBase(f"retired unit has no intact prior served capsule: {capsule}")
    prior_front, _ = parse_frontmatter(page_path.read_bytes(), capsule)
    kind_map = {
        "paper": "paper", "componentized_paper": "paper", "research_governance": "governance",
        "step_family": "step_family", "in_progress_step_family": "step_family",
        "verification_project": "verification", "software_project": "software_project",
        "software_result_family": "result",
    }
    members = old.get("members") or old.get("last_seen_members", [])
    entrypoint = next((m["path"] for m in members if m.get("role") == "primary"), None)
    if entrypoint is None and members:
        entrypoint = members[0]["path"]
    affected = sorted(
        derived_id for derived_id, record in state.get("derived_pages", {}).items()
        if unit_item["id"] in record.get("input_units", [])
    )
    limits = config.get("read_limits", {})
    return {
        "refresh_date": dt.datetime.now().astimezone().date().isoformat(),
        "page": {
            "id": prior_front.get("id", f"source-{unit_item['id']}"),
            "title": prior_front.get("title", unit_item["id"].replace("-", " ").title()),
            "type": "source_capsule",
            "content_owner": "ai_generated",
            "desired_lifecycle": "deleted" if old.get("result") == "source_deleted" else "retired",
            "output_repository_path": capsule,
        },
        "source_kind": kind_map.get(old.get("kind"), "result"),
        "shape": "file" if len(members) == 1 else "bundle",
        "entrypoint": entrypoint,
        "extractor_version": config.get("extractor_contract_version"),
        "prompt_paths": ["/packet/prompts/00-snapshot-contract.md", "/packet/prompts/source-capsule-lifecycle.md"],
        "allowed_citations": [
            {
                "path": member["path"], "availability": "last_seen", "read_mode": member["read_mode"],
                "blob_oid": member["blob_oid"], "last_seen_commit": old.get("processed_commit"),
            } for member in members
        ],
        "related_pages": affected,
        "affected_derived_pages": affected,
        "budgets": {
            "semantic_member_max_bytes": int(limits.get("semantic_member_max_bytes", 2_000_000)),
            "semantic_unit_max_bytes": int(limits.get("semantic_unit_max_bytes", 8_000_000)),
            **dict(config.get("output_budgets", {}).get("source_capsule", {})),
        },
        "prior_identity": {
            "processed_commit": old.get("processed_commit"), "unit_digest_sha256": old.get("unit_digest_sha256"),
            "lifecycle": old.get("lifecycle"), "result": old.get("result"), "members": members, "page": page_record,
        },
        "prior_stable_ids": {
            category: sorted(
                stable_id for stable_id, owner in state.get("id_registry", {}).get(category, {}).items()
                if owner == capsule
            ) for category in ("pages", "statements", "conflicts")
        },
        "stable_id_namespace": f"source-{unit_item['id']}--",
        "lifecycle_action": {
            "action": "retired_from_corpus", "action_commit": target,
            "last_source_commit": old.get("processed_commit"),
            "reason": "unit removed from configured corpus",
        },
    }


def derived_writer_tasks(
    repo: Path, config: Mapping[str, Any], resolved: Mapping[str, ResolvedUnit], state: Mapping[str, Any],
    selected_units: set[str], source_output_units: set[str], target: str, policy: Mapping[str, str],
    unit_actions: Mapping[str, Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    definitions = {item["id"]: item for item in config.get("derived_pages", [])}
    migrations_by_target = migration_requirements(repo)
    selected: set[str] = {
        item["id"] for item in definitions.values() if selected_units.intersection(item["input_units"])
    }
    # A unit retired from config can no longer appear in current input_units.
    # Use the last successful reverse-dependency record to force every still-
    # configured dependent through an explicit refresh/removal task.
    for prior_id, prior_record in state.get("derived_pages", {}).items():
        if prior_id in definitions and selected_units.intersection(prior_record.get("input_units", [])):
            selected.add(prior_id)
    if state.get("policy", {}).get("combined_sha256") not in (None, policy["combined_sha256"]):
        selected = set(definitions)
    changed = True
    while changed:
        changed = False
        for item in definitions.values():
            if item["id"] not in selected and selected.intersection(item["input_pages"]):
                selected.add(item["id"])
                changed = True
    tasks: list[dict[str, Any]] = []
    blocked: list[dict[str, Any]] = []
    prompt_by_kind = {
        "topic": "topic-synthesis.md",
        "script_catalog": "script-catalog.md",
        "conflict_register": "conflict-register.md",
        "index": "index-refresh.md",
    }
    state_units = state.get("units", {})
    state_pages = state.get("pages", {})
    scheduled_derived: set[str] = set()
    for derived_id in sorted(selected, key=lambda item_id: (definitions[item_id]["order"], item_id)):
        item = definitions[derived_id]
        missing_prerequisites: list[str] = []
        for unit_id in item["input_units"]:
            if unit_id in selected_units:
                continue
            record = state_units.get(unit_id)
            capsule = resolved[unit_id].unit["capsule"]
            page = state_pages.get(capsule)
            if (
                record is None or page is None or record.get("unit_digest_sha256") != resolved[unit_id].digest
                or record.get("policy_sha256") != policy["combined_sha256"]
            ):
                missing_prerequisites.append(f"unit:{unit_id}")
        for dependency_id in item["input_pages"]:
            dependency = definitions[dependency_id]
            if dependency_id in selected:
                if dependency_id not in scheduled_derived:
                    missing_prerequisites.append(f"page:{dependency_id}")
                continue
            dependency_record = state.get("derived_pages", {}).get(dependency_id)
            dependency_page = state_pages.get(dependency["page"])
            if (
                dependency_record is None or dependency_page is None
                or dependency_record.get("page_sha256") != dependency_page.get("sha256")
                or dependency_record.get("policy_sha256") != policy["combined_sha256"]
                or dependency_record.get("result") != "success"
            ):
                missing_prerequisites.append(f"page:{dependency_id}")
        if missing_prerequisites:
            blocked.append({
                "id": derived_id,
                "page": item["page"],
                "reason": "missing successful prerequisite",
                "prerequisites": sorted(missing_prerequisites),
            })
            continue
        prior_derived = state.get("derived_pages", {}).get(derived_id)
        selected_source_actions = [
            unit_actions[unit_id]["action"] for unit_id in item["input_units"] if unit_id in selected_units
        ]
        selected_page_tasks = [task for task in tasks if task["task_id"] in item["input_pages"]]
        can_acknowledge = (
            prior_derived is not None
            and prior_derived.get("result") == "success"
            and prior_derived.get("policy_sha256") == policy["combined_sha256"]
            and prior_derived.get("input_units") == sorted(item["input_units"])
            and prior_derived.get("input_pages") == item["input_pages"]
            and bool(selected_source_actions or selected_page_tasks)
            and all(action == "rebase_unchanged" for action in selected_source_actions)
            and all(task.get("action") == "derived_rebase_unchanged" for task in selected_page_tasks)
        )
        if can_acknowledge:
            page_record = state_pages.get(item["page"])
            live_page = repo / item["page"]
            if page_record is None or not live_page.is_file() or sha256_file(live_page) != page_record.get("sha256"):
                missing_prerequisites.append(f"page:{derived_id}:served-integrity")
            else:
                tasks.append({
                    "task_id": derived_id,
                    "source_unit_id": None,
                    "action": "derived_rebase_unchanged",
                    "task_kind": "deterministic_state_only",
                    "phase": item["order"],
                    "required": False,
                    "output_repository_path": item["page"],
                    "staged_output_path": None,
                    "packet_path": None,
                    "semantic_contract": {
                        "input_units": {unit_id: resolved[unit_id].members for unit_id in item["input_units"]},
                        "input_pages": item["input_pages"],
                        "page": {"desired_lifecycle": prior_derived.get("lifecycle", "current")},
                    },
                    "depends_on_tasks": [task["task_id"] for task in selected_page_tasks],
                    "generated_regions": [],
                    "delegated_frontmatter": [],
                })
                scheduled_derived.add(derived_id)
                continue
        base = repo / item["page"]
        if not base.is_file():
            blocked.append({
                "id": derived_id, "page": item["page"], "reason": "mixed base page is missing",
                "prerequisites": [item["page"]],
            })
            continue
        try:
            base_front, _ = parse_frontmatter(base.read_bytes(), item["page"])
        except LintFailure as exc:
            blocked.append({
                "id": derived_id, "page": item["page"], "reason": "mixed base page is invalid",
                "prerequisites": exc.errors,
            })
            continue
        direct_members = {
            unit_id: resolved[unit_id].members for unit_id in item["input_units"]
        }
        dependency_task_ids = [unit_id for unit_id in item["input_units"] if unit_id in source_output_units]
        dependency_task_ids.extend(dependency for dependency in item["input_pages"] if dependency in selected)
        direct_by_path: dict[str, dict[str, Any]] = {}
        direct_units_by_path: dict[str, set[str]] = {}
        read_rank = {"identity_only": 0, "excerpt": 1, "semantic": 2}
        for unit_id, members in direct_members.items():
            for member in members:
                old = direct_by_path.get(member["path"])
                if old is not None and old["blob_oid"] != member["blob_oid"]:
                    raise ConfigError(f"derived task {derived_id} sees conflicting identities for {member['path']}")
                if old is None or read_rank[member["read_mode"]] > read_rank[old["read_mode"]]:
                    direct_by_path[member["path"]] = member
                direct_units_by_path.setdefault(member["path"], set()).add(unit_id)
        all_dependency_paths = set(direct_by_path)
        selected_direct_sources = item.get("direct_sources")
        if selected_direct_sources is not None:
            selected_members: dict[str, dict[str, Any]] = {}
            missing_direct_sources: list[str] = []
            for direct in selected_direct_sources:
                path = direct["path"]
                source_member = direct_by_path.get(path)
                if source_member is None:
                    missing_direct_sources.append(path)
                    continue
                if read_rank[direct["read_mode"]] > read_rank[source_member["read_mode"]]:
                    raise ConfigError(
                        f"derived task {derived_id} direct source escalates {path} from "
                        f"{source_member['read_mode']} to {direct['read_mode']}"
                    )
                selected_member = dict(source_member)
                selected_member["read_mode"] = direct["read_mode"]
                if direct["read_mode"] == "excerpt":
                    selected_member["anchors"] = direct["anchors"]
                else:
                    selected_member.pop("anchors", None)
                    selected_member.pop("context_bytes", None)
                if selected_member["object_type"] != "blob":
                    raise ConfigError(f"derived task {derived_id} cannot read non-blob direct source {path}")
                selected_members[path] = selected_member
            if missing_direct_sources:
                blocked.append({
                    "id": derived_id,
                    "page": item["page"],
                    "reason": "selected direct source is absent from the committed input-unit membership",
                    "prerequisites": sorted(missing_direct_sources),
                })
                continue
            direct_by_path = selected_members
        semantic_inputs = [
            {
                "source_path": member["path"],
                "input_units": sorted(direct_units_by_path[member["path"]]),
                "packet_path": f"inputs/direct/{member['read_mode']}/{member['path']}",
                "sandbox_path": f"/packet/inputs/direct/{member['read_mode']}/{member['path']}",
                "read_mode": member["read_mode"],
                "blob_oid": member["blob_oid"],
                "member_identity": member,
            }
            for member in (direct_by_path[path] for path in sorted(direct_by_path))
            if member["read_mode"] in ("semantic", "excerpt")
        ]
        readable_paths = {
            member["path"] for member in direct_by_path.values()
            if member["read_mode"] in ("semantic", "excerpt")
        }
        allowed_paths = sorted(readable_paths)
        identity_only = sorted(all_dependency_paths - readable_paths)
        migration_items = migrations_by_target.get(item["page"], [])
        missing_migration_sources = sorted({
            reference.split("#", 1)[0]
            for migration in migration_items for reference in migration["required_original_refs"]
            if reference.split("#", 1)[0] not in readable_paths
        })
        if missing_migration_sources:
            blocked.append({
                "id": derived_id,
                "page": item["page"],
                "reason": "migration original refs are not selected readable direct_sources",
                "prerequisites": missing_migration_sources,
            })
            continue
        excerpt_bytes = int(config.get("read_limits", {}).get("excerpt_context_bytes", 8192))
        derived_semantic_bytes = sum(
            (member["blob_size"] or 0) if member["read_mode"] == "semantic" else excerpt_bytes
            for member in direct_by_path.values() if member["read_mode"] in ("semantic", "excerpt")
        )
        derived_limit = int(config.get("read_limits", {}).get("derived_task_max_bytes", 8_000_000))
        if derived_semantic_bytes > derived_limit:
            raise ConfigError(
                f"derived task {derived_id} has {derived_semantic_bytes} deduplicated direct semantic bytes, "
                f"above {derived_limit}"
            )
        static_memory_inputs: list[dict[str, Any]] = []
        for unit_id in item["input_units"]:
            if unit_id not in source_output_units:
                capsule = resolved[unit_id].unit["capsule"]
                if capsule in state_pages:
                    unit_record = state_units[unit_id]
                    page_record = state_pages[capsule]
                    dependency_front, _ = parse_frontmatter((repo / capsule).read_bytes(), capsule)
                    static_memory_inputs.append({
                        "id": unit_id,
                        "dependency_kind": "source_capsule",
                        "page": capsule,
                        "page_sha256": page_record["sha256"],
                        "generated_commit": dependency_front.get("generated_from_commit"),
                        "lifecycle": dependency_front.get("lifecycle"),
                        "unit_digest_sha256": unit_record.get("unit_digest_sha256"),
                        "policy_sha256": unit_record.get("policy_sha256"),
                        "policy_fresh": unit_record.get("policy_sha256") == policy["combined_sha256"],
                    })
        for page_id in item["input_pages"]:
            if page_id not in selected:
                page_path = definitions[page_id]["page"]
                derived_record = state.get("derived_pages", {})[page_id]
                page_record = state_pages[page_path]
                dependency_front, _ = parse_frontmatter((repo / page_path).read_bytes(), page_path)
                static_memory_inputs.append({
                    "id": page_id,
                    "dependency_kind": "derived_page",
                    "page": page_path,
                    "page_sha256": page_record["sha256"],
                    "generated_commit": dependency_front.get("generated_from_commit"),
                    "lifecycle": dependency_front.get("lifecycle"),
                    "unit_digest_sha256": None,
                    "policy_sha256": derived_record.get("policy_sha256"),
                    "policy_fresh": derived_record.get("policy_sha256") == policy["combined_sha256"],
                })
        dynamic_memory_inputs: list[dict[str, Any]] = []
        prior_tasks = {task["task_id"]: task for task in tasks}
        for task_id in sorted(set(dependency_task_ids)):
            if task_id in resolved:
                dependency_page = resolved[task_id].unit["capsule"]
                dependency_lifecycle = "deleted" if not resolved[task_id].members else resolved[task_id].unit["lifecycle"]
                unit_digest = resolved[task_id].digest
                dependency_kind = "source_capsule"
            else:
                dependency_task = prior_tasks[task_id]
                dependency_page = definitions[task_id]["page"]
                dependency_lifecycle = dependency_task["semantic_contract"]["page"]["desired_lifecycle"]
                unit_digest = None
                dependency_kind = "derived_page"
            dynamic_memory_inputs.append({
                "task_id": task_id,
                "dependency_kind": dependency_kind,
                "page": dependency_page,
                "sandbox_path": "/packet/memory_inputs/" + dependency_page,
                "generated_commit": target,
                "lifecycle": dependency_lifecycle,
                "unit_digest_sha256": unit_digest,
                "policy_sha256": policy["combined_sha256"],
                "policy_fresh": True,
            })
        contract = {
            "refresh_date": dt.datetime.now().astimezone().date().isoformat(),
            "page": {
                "id": base_front.get("id", derived_id),
                "title": base_front.get("title", derived_id.replace("-", " ").title()),
                "type": item["task_kind"],
                "content_owner": "mixed",
                "desired_lifecycle": base_front.get("lifecycle", "current"),
                "output_repository_path": item["page"],
            },
            "prompt_paths": ["/packet/prompts/00-snapshot-contract.md", f"/packet/prompts/{prompt_by_kind[item['task_kind']]}"],
            "allowed_citations": allowed_paths,
            "migration_requirements": migration_items,
            "input_units": direct_members,
            "input_pages": item["input_pages"],
            "static_memory_inputs": static_memory_inputs,
            "dynamic_memory_inputs": dynamic_memory_inputs,
            "budget": config["output_budgets"][item["budget"]],
            "direct_semantic_bytes": derived_semantic_bytes,
            "direct_semantic_max_bytes": derived_limit,
            "prior_stable_ids": {
                category: sorted(
                    stable_id for stable_id, page in state.get("id_registry", {}).get(category, {}).items()
                    if page == item["page"]
                ) for category in ("pages", "statements", "conflicts")
            },
            "stable_id_namespace": f"{derived_id}--",
            "generated_region": item["region"],
        }
        tasks.append({
            "task_id": derived_id,
            "source_unit_id": None,
            "action": "refresh_derived",
            "task_kind": item["task_kind"],
            "phase": item["order"],
            "required": True,
            "output_repository_path": item["page"],
            "staged_output_path": f"staged/{item['page']}",
            "packet_path": f"tasks/{derived_id}/packet",
            "sandbox_packet_path": "/packet",
            "sandbox_output_path": "/output/page.md",
            "expected_output_name": "page.md",
            "semantic_inputs": semantic_inputs,
            "identity_only_members": identity_only,
            "semantic_contract": contract,
            "depends_on_tasks": sorted(set(dependency_task_ids)),
            "generated_regions": [item["region"]],
            "delegated_frontmatter": ["sources", "last_updated", "generated_from_commit", "memory_review"],
            "base_page_sha256": sha256_file(base),
        })
        scheduled_derived.add(derived_id)
    return tasks, blocked


def _selected_unit_ids(
    requested_units: Sequence[str], requested_paths: Sequence[str], resolved: Mapping[str, ResolvedUnit],
    state: Mapping[str, Any], actions: Mapping[str, Mapping[str, Any]],
) -> list[str]:
    configured = set(resolved)
    retired = set(actions) - configured
    selected: set[str] = set()
    for unit_id in requested_units:
        if unit_id not in actions:
            raise MemoryErrorBase(f"unknown configured or retired unit: {unit_id}")
        selected.add(unit_id)
    owner: dict[str, str] = {}
    for unit_id, unit in resolved.items():
        for member in unit.members:
            if member["ownership"] == "primary":
                owner[member["path"]] = unit_id
    for unit_id, old in state.get("units", {}).items():
        for member in old.get("members", []):
            if member.get("ownership", "primary") == "primary":
                owner.setdefault(member["path"], unit_id)
    for path in requested_paths:
        normalized = clean_repo_path(path, "--paths value")
        unit_id = owner.get(normalized)
        if unit_id is None:
            raise MemoryErrorBase(f"path has no configured primary owner: {normalized}")
        selected.add(unit_id)
    if not requested_units and not requested_paths:
        selected = {unit_id for unit_id, detail in actions.items() if detail["action"] != "current"}
    return sorted(selected)


def prepare_update(
    repo: Path, ref: str, requested_units: Sequence[str], requested_paths: Sequence[str],
    config_path: str = DEFAULT_CONFIG, allow_large_batch: bool = False,
) -> Path:
    config, config_bytes = load_config(repo, config_path)
    target = resolve_commit(repo, ref)
    entries = enumerate_tree(repo, target, config["discovery"]["candidate_roots"])
    resolved = resolve_units(config, entries)
    state, state_bytes = load_state(repo)
    policy = policy_digests(repo, config_bytes)
    configured_derived_ids = {item["id"] for item in config.get("derived_pages", [])}
    retired_derived_ids = sorted(set(state.get("derived_pages", {})) - configured_derived_ids)
    if retired_derived_ids:
        raise MemoryErrorBase(
            "derived-page retirement requires an explicit lifecycle/cutover task and is not automatic: "
            + ", ".join(retired_derived_ids)
        )
    actions = compute_actions(config, resolved, state, target, policy)
    selected = _selected_unit_ids(requested_units, requested_paths, resolved, state, actions)
    if not selected:
        raise MemoryErrorBase("no pending units selected")
    invalid = [unit_id for unit_id in selected if actions[unit_id]["action"] == "invalid_initial"]
    if invalid:
        detail = "; ".join(f"{x}: {actions[x]['missing_required']}" for x in invalid)
        raise MemoryErrorBase(f"required configured inputs are missing on initial ingest: {detail}")
    limits = config.get("read_limits", {})
    batch_limits = {
        "max_units": int(limits.get("prepare_batch_max_units", 8)),
        "max_tasks": int(limits.get("prepare_batch_max_tasks", 16)),
        "max_semantic_bytes": int(limits.get("prepare_batch_max_semantic_bytes", 24_000_000)),
    }
    selected_semantic_bytes = sum(resolved[unit_id].semantic_bytes for unit_id in selected if unit_id in resolved)
    batch_excess = []
    if len(selected) > batch_limits["max_units"]:
        batch_excess.append(f"{len(selected)} units > {batch_limits['max_units']}")
    if selected_semantic_bytes > batch_limits["max_semantic_bytes"]:
        batch_excess.append(
            f"{selected_semantic_bytes} semantic bytes > {batch_limits['max_semantic_bytes']}"
        )
    if batch_excess and not allow_large_batch:
        raise MemoryErrorBase(
            "prepare batch exceeds bounded defaults (" + "; ".join(batch_excess)
            + "); select fewer --units or pass --allow-large-batch explicitly"
        )
    transaction_id = f"{target[:12]}-{dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    transaction = repo / "memory/transactions" / transaction_id
    transaction.mkdir(parents=True, exist_ok=False)
    (transaction / "snapshot/semantic").mkdir(parents=True)
    (transaction / "snapshot/excerpts").mkdir(parents=True)
    (transaction / "policy/prompts").mkdir(parents=True)
    (transaction / "staged").mkdir(parents=True)
    default_context = int(limits.get("excerpt_context_bytes", 8192))
    txn_units: list[dict[str, Any]] = []
    snapshot_hashes: dict[str, str] = {}
    try:
        # Extraction is forbidden from consulting mutable live policy just as
        # it is forbidden from consulting mutable candidate sources.  Freeze
        # every policy input named by the policy digest into the packet.
        frozen_policy_sources = {
            "policy/config.yaml": repo / config_path,
            "policy/schema.md": repo / DEFAULT_SCHEMA,
            "policy/SCANNER_DECISIONS.md": repo / "memory/_meta/SCANNER_DECISIONS.md",
            "policy/atlas-migration.yaml": repo / "memory/_meta/atlas-migration.yaml",
        }
        for relative, source in frozen_policy_sources.items():
            destination = transaction / relative
            atomic_write(destination, source.read_bytes(), 0o444)
            snapshot_hashes[relative] = sha256_file(destination)
        prompts_root = repo / "memory/prompts"
        frozen_prompts: list[str] = []
        if prompts_root.is_dir():
            for source in sorted(path for path in prompts_root.rglob("*") if path.is_file()):
                prompt_relative = source.relative_to(prompts_root).as_posix()
                relative = f"policy/prompts/{prompt_relative}"
                destination = transaction / relative
                atomic_write(destination, source.read_bytes(), 0o444)
                snapshot_hashes[relative] = sha256_file(destination)
                frozen_prompts.append(relative)
        frozen_policy_checks = {
            "config_sha256": sha256_file(transaction / "policy/config.yaml"),
            "schema_sha256": sha256_file(transaction / "policy/schema.md"),
            "contract_sha256": sha256_file(transaction / "policy/SCANNER_DECISIONS.md"),
            "atlas_migration_sha256": sha256_file(transaction / "policy/atlas-migration.yaml"),
            "prompts_sha256": digest_directory(transaction / "policy/prompts"),
        }
        for digest_name, frozen_digest in frozen_policy_checks.items():
            if frozen_digest != policy[digest_name]:
                raise MemoryErrorBase(f"live {digest_name} changed while freezing transaction policy")
        for unit_id in selected:
            action = actions[unit_id]["action"]
            if unit_id not in resolved:
                old = state["units"][unit_id]
                txn_units.append({
                    "id": unit_id, "action": action, "capsule": old.get("capsule"),
                    "unit_digest_sha256": old.get("unit_digest_sha256"),
                    "members": old.get("members") or old.get("last_seen_members", []),
                    "kind": old.get("kind"), "lifecycle": "retired",
                    "missing_required": [], "requires_staged_page": True,
                })
                continue
            unit = resolved[unit_id]
            members: list[dict[str, Any]] = []
            for member in unit.members:
                staged_member = dict(member)
                if member["read_mode"] in ("semantic", "excerpt"):
                    data = read_blob(repo, member["blob_oid"])
                    if member["blob_size"] is not None and len(data) != member["blob_size"]:
                        raise GitError(f"blob size changed while reading {member['path']}")
                    if member["read_mode"] == "semantic":
                        if b"\0" in data:
                            raise MemoryErrorBase(f"semantic blob looks binary: {member['path']}")
                        relative = PurePosixPath("snapshot/semantic") / PurePosixPath(member["path"])
                        output = data
                    else:
                        relative = PurePosixPath("snapshot/excerpts") / PurePosixPath(member["path"] + ".excerpt.txt")
                        output = deterministic_excerpt(data, member, default_context)
                    output_path = transaction / relative.as_posix()
                    atomic_write(output_path, output, 0o444)
                    staged_member["snapshot_path"] = relative.as_posix()
                    staged_member["snapshot_sha256"] = sha256_bytes(output)
                    snapshot_hashes[relative.as_posix()] = staged_member["snapshot_sha256"]
                members.append(staged_member)
            old_unit = state.get("units", {}).get(unit_id, {})
            current_identity = {member["path"]: member.get("blob_oid") for member in members}
            prior_members: list[dict[str, Any]] = []
            for prior_member in old_unit.get("last_seen_members") or old_unit.get("members", []):
                if current_identity.get(prior_member["path"]) == prior_member.get("blob_oid"):
                    continue
                staged_prior = dict(prior_member)
                if prior_member.get("read_mode") in ("semantic", "excerpt"):
                    data = read_blob(repo, prior_member["blob_oid"])
                    if prior_member["read_mode"] == "semantic":
                        if b"\0" in data:
                            raise MemoryErrorBase(f"prior semantic blob looks binary: {prior_member['path']}")
                        relative = PurePosixPath("snapshot/prior") / unit_id / PurePosixPath(prior_member["path"])
                        output = data
                    else:
                        relative = (
                            PurePosixPath("snapshot/prior") / unit_id
                            / PurePosixPath(prior_member["path"] + ".excerpt.txt")
                        )
                        output = deterministic_excerpt(data, prior_member, default_context)
                    output_path = transaction / relative.as_posix()
                    atomic_write(output_path, output, 0o444)
                    staged_prior["snapshot_path"] = relative.as_posix()
                    staged_prior["snapshot_sha256"] = sha256_bytes(output)
                    snapshot_hashes[relative.as_posix()] = staged_prior["snapshot_sha256"]
                prior_members.append(staged_prior)
            requires_page = action not in ("rebase_unchanged", "retired_from_corpus")
            txn_units.append({
                "id": unit_id,
                "action": action,
                "kind": unit.unit["kind"],
                "lifecycle": unit.unit["lifecycle"],
                "capsule": unit.unit["capsule"],
                "unit_digest_sha256": unit.digest,
                "members": members,
                "prior_members": prior_members,
                "missing_required": unit.missing_required,
                "requires_staged_page": requires_page,
            })
        writer_tasks = []
        for item in txn_units:
            required = bool(item.get("requires_staged_page"))
            semantic_contract = (
                source_task_metadata(config, item, resolved[item["id"]], state, target)
                if required and item["id"] in resolved else None
            )
            if required and item["action"] == "retired_from_corpus":
                semantic_contract = retired_source_task_metadata(repo, config, item, state, target)
            if semantic_contract is not None:
                available_prompt_paths = {"/packet/" + path.removeprefix("policy/") for path in frozen_prompts}
                absent_prompts = sorted(set(semantic_contract["prompt_paths"]) - available_prompt_paths)
                if absent_prompts:
                    raise MemoryErrorBase(
                        f"semantic prompt inputs are absent for {item['id']}: {absent_prompts}"
                    )
            writer_tasks.append({
                "task_id": item["id"],
                "source_unit_id": item["id"],
                "action": item["action"],
                "task_kind": "source_capsule" if required else "deterministic_state_only",
                "phase": 1 if required else 0,
                "required": required,
                "output_repository_path": item.get("capsule"),
                "staged_output_path": f"staged/{item['capsule']}" if required else None,
                "packet_path": f"tasks/{item['id']}/packet" if required else None,
                "sandbox_packet_path": "/packet" if required else None,
                "sandbox_output_path": "/output/page.md" if required else None,
                "expected_output_name": "page.md" if required else None,
                "semantic_inputs": [
                    {
                        "source_path": member["path"],
                        "packet_path": (
                            "inputs/" + member["snapshot_path"].removeprefix("snapshot/")
                            if member.get("snapshot_path") else f"inputs/last_seen/{member['read_mode']}/{member['path']}"
                        ),
                        "sandbox_path": (
                            "/packet/inputs/" + member["snapshot_path"].removeprefix("snapshot/")
                            if member.get("snapshot_path") else f"/packet/inputs/last_seen/{member['read_mode']}/{member['path']}"
                        ),
                    }
                    for member in item.get("members", [])
                    if member.get("read_mode") in ("semantic", "excerpt")
                ] + [
                    {
                        "source_path": member["path"],
                        "availability": "last_seen",
                        "blob_oid": member.get("blob_oid"),
                        "packet_path": "inputs/" + member["snapshot_path"].removeprefix("snapshot/"),
                        "sandbox_path": "/packet/inputs/" + member["snapshot_path"].removeprefix("snapshot/"),
                    }
                    for member in item.get("prior_members", [])
                    if member.get("read_mode") in ("semantic", "excerpt") and member.get("snapshot_path")
                ],
                "identity_only_members": [
                    member["path"] for member in item.get("members", [])
                    if member.get("read_mode") == "identity_only"
                ] + [
                    member["path"] for member in item.get("prior_members", [])
                    if member.get("read_mode") == "identity_only"
                ],
                "semantic_contract": semantic_contract,
                "depends_on_tasks": [],
                "generated_regions": [],
                "delegated_frontmatter": [],
            })
        source_output_units = {task["source_unit_id"] for task in writer_tasks if task["required"]}
        derived_tasks, blocked_derived = derived_writer_tasks(
            repo, config, resolved, state, set(selected), source_output_units, target, policy, actions
        )
        available_prompt_paths = {"/packet/" + path.removeprefix("policy/") for path in frozen_prompts}
        for task in derived_tasks:
            if not task.get("required"):
                continue
            absent_prompts = sorted(set(task["semantic_contract"]["prompt_paths"]) - available_prompt_paths)
            if absent_prompts:
                raise MemoryErrorBase(f"semantic prompt inputs are absent for {task['task_id']}: {absent_prompts}")
        writer_tasks.extend(derived_tasks)
        if len([task for task in writer_tasks if task.get("required")]) > batch_limits["max_tasks"] and not allow_large_batch:
            raise MemoryErrorBase(
                f"prepare batch has {len([task for task in writer_tasks if task.get('required')])} writer tasks, "
                f"above {batch_limits['max_tasks']}; select fewer --units or pass --allow-large-batch explicitly"
            )
        task_phases = sorted({task["phase"] for task in writer_tasks})
        manifest = {
            "input_format_version": INPUT_FORMAT_VERSION,
            "transaction_id": transaction_id,
            "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "target_commit": target,
            "prior_state_sha256": state_digest(state_bytes),
            "policy": policy,
            "candidate_roots": config["discovery"]["candidate_roots"],
            "frozen_policy": {
                "config": "policy/config.yaml",
                "schema": "policy/schema.md",
                "extraction_contract": "policy/SCANNER_DECISIONS.md",
                "atlas_migration": "policy/atlas-migration.yaml",
                "prompts": frozen_prompts,
            },
            "selection": {"units": list(requested_units), "paths": list(requested_paths)},
            "batch_bounds": {
                **batch_limits,
                "selected_units": len(selected),
                "selected_semantic_bytes": selected_semantic_bytes,
                "required_writer_tasks": len([task for task in writer_tasks if task.get("required")]),
                "explicit_large_batch_override": bool(allow_large_batch),
            },
            "units": txn_units,
            "writer_tasks": writer_tasks,
            "task_order": [
                {"phase": phase, "tasks": [task["task_id"] for task in writer_tasks if task["phase"] == phase]}
                for phase in task_phases
            ],
            "derived_task_plan": {
                "status": "complete" if not blocked_derived else "blocked_partial",
                "tasks": [task["task_id"] for task in derived_tasks],
                "blocked": blocked_derived,
            },
            "configured_derived_ids": [item["id"] for item in config.get("derived_pages", [])],
            "instructions": [
                "Read semantic/excerpt inputs only from snapshot_path values in this transaction.",
                "Do not read live research/ or software/ candidate roots while extracting this transaction.",
                "Read schema and extraction prompts only from frozen_policy paths in this transaction.",
                "Citations in pages must use the original repository-relative member path, never snapshot_path.",
                "Write each required capsule under staged/<capsule path>.",
            ],
        }
        manifest_path = transaction / "manifest.json"
        write_json(manifest_path, manifest)
        initial_task_packets: dict[str, dict[str, Any]] = {}
        initial_task_documents: dict[str, dict[str, Any]] = {}
        for writer_task in writer_tasks:
            if not writer_task["required"]:
                continue
            task_id = writer_task["task_id"]
            unit_item = next(
                (item for item in txn_units if item["id"] == writer_task["source_unit_id"]), None
            )
            packet = transaction / writer_task["packet_path"]
            output = packet.parent / "output"
            packet.mkdir(parents=True)
            output.mkdir()
            task_document = {
                "transaction_id": transaction_id,
                "target_commit": target,
                "writer_task": writer_task,
                "source_unit": unit_item,
                "semantic_contract": writer_task["semantic_contract"],
                "frozen_policy": {
                    "schema": "/packet/schema.md",
                    "config": "/packet/config.yaml",
                    "extraction_contract": "/packet/SCANNER_DECISIONS.md",
                    "atlas_migration": "/packet/atlas-migration.yaml",
                    "prompts_root": "/packet/prompts",
                },
                "base_page": None,
                "output": "/output/page.md",
                "hard_boundary": (
                    "The isolated writer can see this packet and an output directory, not the live repository. "
                    "Citations use source_path values from this document."
                ),
            }
            for source_name, packet_name in (
                ("policy/schema.md", "schema.md"),
                ("policy/config.yaml", "config.yaml"),
                ("policy/SCANNER_DECISIONS.md", "SCANNER_DECISIONS.md"),
                ("policy/atlas-migration.yaml", "atlas-migration.yaml"),
            ):
                atomic_write(packet / packet_name, (transaction / source_name).read_bytes(), 0o444)
            prompt_packet = packet / "prompts"
            prompt_packet.mkdir()
            for prompt_relative in frozen_prompts:
                destination = packet / prompt_relative.removeprefix("policy/")
                atomic_write(destination, (transaction / prompt_relative).read_bytes(), 0o444)
            if unit_item is not None:
                semantic_input_by_source = {
                    value["source_path"]: value for value in writer_task["semantic_inputs"]
                }
                for member in unit_item.get("members", []):
                    source_snapshot = member.get("snapshot_path")
                    if not source_snapshot:
                        if unit_item["action"] == "retired_from_corpus" and member["read_mode"] != "identity_only":
                            data = read_blob(repo, member["blob_oid"])
                            if member["read_mode"] == "excerpt":
                                data = deterministic_excerpt(data, member, default_context)
                            destination = packet / semantic_input_by_source[member["path"]]["packet_path"]
                            atomic_write(destination, data, 0o444)
                        continue
                    destination = packet / "inputs" / source_snapshot.removeprefix("snapshot/")
                    atomic_write(destination, (transaction / source_snapshot).read_bytes(), 0o444)
                for member in unit_item.get("prior_members", []):
                    source_snapshot = member.get("snapshot_path")
                    if source_snapshot:
                        destination = packet / "inputs" / source_snapshot.removeprefix("snapshot/")
                        atomic_write(destination, (transaction / source_snapshot).read_bytes(), 0o444)
                page_record = state.get("pages", {}).get(unit_item["capsule"])
                if page_record is not None:
                    live_page = repo / unit_item["capsule"]
                    if not live_page.is_file() or sha256_file(live_page) != page_record.get("sha256"):
                        raise MemoryErrorBase(f"successful base page is missing or modified: {unit_item['capsule']}")
                    atomic_write(packet / "base_page.md", live_page.read_bytes(), 0o444)
                    task_document["base_page"] = "/packet/base_page.md"
            else:
                contract = writer_task["semantic_contract"]
                base_page = repo / writer_task["output_repository_path"]
                if not base_page.is_file() or sha256_file(base_page) != writer_task["base_page_sha256"]:
                    raise MemoryErrorBase(f"derived mixed base page changed during prepare: {writer_task['output_repository_path']}")
                atomic_write(packet / "base_page.md", base_page.read_bytes(), 0o444)
                task_document["base_page"] = "/packet/base_page.md"
                for semantic_input in writer_task["semantic_inputs"]:
                    member = semantic_input["member_identity"]
                    data = read_blob(repo, member["blob_oid"])
                    if member["read_mode"] == "excerpt":
                        data = deterministic_excerpt(data, member, default_context)
                    destination = packet / semantic_input["packet_path"]
                    atomic_write(destination, data, 0o444)
                for memory_input in contract.get("static_memory_inputs", []):
                    page_path = memory_input["page"]
                    record = state.get("pages", {}).get(page_path)
                    live_page = repo / page_path
                    if (
                        record is None or not live_page.is_file()
                        or sha256_file(live_page) != record.get("sha256")
                        or record.get("sha256") != memory_input.get("page_sha256")
                        or memory_input.get("policy_sha256") != policy["combined_sha256"]
                        or memory_input.get("policy_fresh") is not True
                    ):
                        raise MemoryErrorBase(f"successful static memory dependency is unavailable: {page_path}")
                    dependency_front, _ = parse_frontmatter(live_page.read_bytes(), page_path)
                    if dependency_front.get("generated_from_commit") != memory_input.get("generated_commit"):
                        raise MemoryErrorBase(f"static memory dependency commit disagrees with state: {page_path}")
                    if dependency_front.get("lifecycle") != memory_input.get("lifecycle"):
                        raise MemoryErrorBase(f"static memory dependency lifecycle disagrees with state: {page_path}")
                    expected_unit_digest = memory_input.get("unit_digest_sha256")
                    if expected_unit_digest is not None:
                        source_unit = dependency_front.get("source_unit")
                        actual_digest = source_unit.get("unit_digest_sha256") if isinstance(source_unit, dict) else None
                        if actual_digest != expected_unit_digest:
                            raise MemoryErrorBase(f"static memory dependency unit digest disagrees with state: {page_path}")
                    destination = packet / "memory_inputs" / page_path
                    atomic_write(destination, live_page.read_bytes(), 0o444)
            atomic_write(packet / "task.json", canonical_json(task_document), 0o444)
            atomic_write(packet / "transaction-manifest.json", canonical_json(manifest), 0o444)
            packet_files = {
                path.relative_to(packet).as_posix(): sha256_file(path)
                for path in sorted(item for item in packet.rglob("*") if item.is_file())
            }
            packet_seal = {
                "transaction_id": transaction_id,
                "task_id": task_id,
                "source_unit_id": writer_task["source_unit_id"],
                "files": packet_files,
                "combined_sha256": sha256_bytes(canonical_json(packet_files)),
            }
            write_json(packet.parent / "packet-seal.json", packet_seal)
            initial_task_packets[task_id] = packet_seal
            initial_task_documents[task_id] = task_document
        seal = {
            "manifest_sha256": sha256_file(manifest_path),
            "snapshot_files": snapshot_hashes,
            "initial_task_packets": initial_task_packets,
            "initial_task_documents": initial_task_documents,
        }
        write_json(transaction / "seal.json", seal)
        os.chmod(manifest_path, 0o444)
        os.chmod(transaction / "seal.json", 0o444)
    except Exception:
        shutil.rmtree(transaction)
        raise
    return transaction


def packet_file_hashes(packet: Path) -> dict[str, str]:
    return {
        path.relative_to(packet).as_posix(): sha256_file(path)
        for path in sorted(item for item in packet.rglob("*") if item.is_file())
    }


def revision_task_contract(manifest: Mapping[str, Any], task: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    source_unit_id = task.get("source_unit_id")
    if source_unit_id is not None:
        unit = next(item for item in manifest["units"] if item.get("id") == source_unit_id)
        identities = {
            key: [
                {field: value for field, value in member.items() if field != "snapshot_path"}
                for member in unit.get(key, [])
            ]
            for key in ("members", "prior_members")
        }
        return ({
            "task_kind": task["task_kind"], "source_unit_id": source_unit_id,
            "source_kind": unit["kind"], "output_repository_path": task["output_repository_path"],
            "member_identities": identities,
        }, {"unit_digest_sha256": unit["unit_digest_sha256"]})
    contract = task["semantic_contract"]
    return ({
        "task_kind": task["task_kind"], "source_unit_id": None,
        "output_repository_path": task["output_repository_path"],
        "generated_regions": task.get("generated_regions", []),
        "depends_on_tasks": task.get("depends_on_tasks", []),
        "input_unit_ids": sorted(contract.get("input_units", {})),
        "input_page_ids": list(contract.get("input_pages", [])),
        "direct_sources": [
            {"source_path": item.get("source_path"), "read_mode": item.get("read_mode"),
             "member_identity": item.get("member_identity")}
            for item in task.get("semantic_inputs", [])
        ],
    }, {
        "base_page_sha256": task.get("base_page_sha256"),
        "static_memory_inputs": contract.get("static_memory_inputs", []),
        "dynamic_memory_inputs": contract.get("dynamic_memory_inputs", []),
    })


def verify_task_packet_chain(
    transaction: Path, manifest: Mapping[str, Any], transaction_seal: Mapping[str, Any], task: Mapping[str, Any],
) -> dict[str, Any]:
    task_id = task["task_id"]
    packet = transaction / task["packet_path"]
    packet_seal_path = packet.parent / "packet-seal.json"
    initial = transaction_seal.get("initial_task_packets", {}).get(task_id)
    initial_document = transaction_seal.get("initial_task_documents", {}).get(task_id)
    if not isinstance(initial, dict) or not isinstance(initial_document, dict):
        raise MemoryErrorBase(f"transaction seal does not anchor initial task packet: {task_id}")
    if not packet.is_dir() or not packet_seal_path.is_file():
        raise MemoryErrorBase(f"sealed task packet is missing: {task_id}")
    current_seal = json.loads(packet_seal_path.read_text(encoding="utf-8"))
    current_files = packet_file_hashes(packet)
    current_combined = sha256_bytes(canonical_json(current_files))
    if current_seal.get("files") != current_files or current_seal.get("combined_sha256") != current_combined:
        raise MemoryErrorBase(f"current task packet seal is invalid: {task_id}")
    if current_files == initial.get("files") and current_combined == initial.get("combined_sha256"):
        return current_seal
    dynamic = initial_document.get("semantic_contract", {}).get("dynamic_memory_inputs", [])
    revision_path = packet / "revision.json"
    has_revision = revision_path.is_file()
    revision: dict[str, Any] | None = None
    revision_basis: str | None = None
    if has_revision:
        revision = json.loads(revision_path.read_text(encoding="utf-8"))
        revision_basis = revision.get("revision_basis")
    if not dynamic and not has_revision:
        raise MemoryErrorBase(f"task packet changed without declared dynamic dependencies: {task_id}")
    allowed_additions: set[str] = set()
    hydration: dict[str, Any] | None = None
    if dynamic:
        hydration_path = packet / "hydration.json"
        if not hydration_path.is_file():
            raise MemoryErrorBase(f"dynamic task packet lacks append-only hydration record: {task_id}")
        hydration = json.loads(hydration_path.read_text(encoding="utf-8"))
        if hydration.get("initial_packet_sha256") != initial.get("combined_sha256"):
            raise MemoryErrorBase(f"dynamic task packet does not chain to its initial seal: {task_id}")
        allowed_additions |= {"hydration.json"} | {f"memory_inputs/{item['page']}" for item in dynamic}
    if has_revision:
        allowed_additions |= {
            "revision.json", "revision_candidate.md", "revision_review.md",
            "revision_review_attestation.json",
        }
        if revision_basis == "machine_lint_failure":
            allowed_additions |= {"revision_lint_record.json", "revision_lint_report.md"}
    expected_files = set(initial.get("files", {})) | allowed_additions
    if set(current_files) != expected_files:
        raise MemoryErrorBase(f"task packet has undeclared additions/removals: {task_id}")
    for relative, digest in initial.get("files", {}).items():
        if relative == "task.json":
            continue
        if current_files.get(relative) != digest:
            raise MemoryErrorBase(f"task hydration modified immutable packet input {relative}: {task_id}")
    enriched = json.loads(json.dumps(initial_document))
    if dynamic:
        assert hydration is not None
        dependency_records = hydration.get("dependencies")
        if not isinstance(dependency_records, list) or {item.get("task_id") for item in dependency_records} != {
            item["task_id"] for item in dynamic
        }:
            raise MemoryErrorBase(f"dynamic hydration dependency set mismatch: {task_id}")
        enrich_maps = [
            {item["task_id"]: item for item in enriched["semantic_contract"]["dynamic_memory_inputs"]},
            {item["task_id"]: item for item in enriched["writer_task"]["semantic_contract"]["dynamic_memory_inputs"]},
        ]
        for record in dependency_records:
            dependency_id = record["task_id"]
            attestation_path = transaction / "attestations" / f"{dependency_id}.json"
            staged_path = transaction / record["staged_output_path"]
            if (
                not attestation_path.is_file() or sha256_file(attestation_path) != record.get("attestation_sha256")
                or not staged_path.is_file() or sha256_file(staged_path) != record.get("page_sha256")
            ):
                raise MemoryErrorBase(f"dynamic hydration dependency identity changed: {dependency_id}")
            attestation = json.loads(attestation_path.read_text(encoding="utf-8"))
            if (
                attestation.get("transaction_id") != manifest["transaction_id"]
                or attestation.get("task_id") != dependency_id
                or attestation.get("output_sha256") != record.get("page_sha256")
                or attestation.get("packet_sha256") != record.get("attestation_packet_sha256")
            ):
                raise MemoryErrorBase(f"dynamic hydration attestation mismatch: {dependency_id}")
            fields = {
                "page_sha256": record["page_sha256"],
                "attestation_sha256": record["attestation_sha256"],
                "attestation_packet_sha256": record["attestation_packet_sha256"],
                "attestation_identity": {
                    "transaction_id": manifest["transaction_id"], "task_id": dependency_id,
                },
            }
            for mapping in enrich_maps:
                mapping[dependency_id].update(fields)
    current_document = json.loads((packet / "task.json").read_text(encoding="utf-8"))
    if current_document != enriched:
        raise MemoryErrorBase(f"task.json contains changes beyond verified dependency identities: {task_id}")
    if has_revision:
        assert revision is not None
        candidate_path = packet / "revision_candidate.md"
        review_path = packet / "revision_review.md"
        review_attestation_path = packet / "revision_review_attestation.json"
        lint_revision = revision_basis == "machine_lint_failure"
        lint_record_path = packet / "revision_lint_record.json"
        lint_report_path = packet / "revision_lint_report.md"
        selector, prerequisites = revision_task_contract(manifest, task)
        revision_names = {
            "revision.json", "revision_candidate.md", "revision_review.md",
            "revision_review_attestation.json",
        }
        if lint_revision:
            revision_names |= {"revision_lint_record.json", "revision_lint_report.md"}
        base_files = {
            relative: digest for relative, digest in current_files.items()
            if relative not in revision_names
        }
        expected_revision = {
            "revision_version": 1,
            "transaction_id": manifest["transaction_id"],
            "task_id": task_id,
            "source_unit_id": task.get("source_unit_id"),
            "target_commit": manifest["target_commit"],
            "output_repository_path": task["output_repository_path"],
            "base_packet_sha256": sha256_bytes(canonical_json(base_files)),
            "candidate_sha256": sha256_file(candidate_path),
            "task_kind": task["task_kind"],
            "retry_selector_sha256": sha256_bytes(canonical_json(selector)),
            "current_prerequisites_sha256": sha256_bytes(canonical_json(prerequisites)),
            "prior_review_attestation_sha256": sha256_file(review_attestation_path),
            "prior_review_report_sha256": sha256_file(review_path),
        }
        for key, value in expected_revision.items():
            if revision.get(key) != value:
                raise MemoryErrorBase(f"revision packet identity mismatch for {task_id}: {key}")
        prior_expected = {
            "prior_output_repository_path": task["output_repository_path"],
            "prior_staged_output_path": task["staged_output_path"],
        }
        for key, value in prior_expected.items():
            if revision.get(key) != value:
                raise MemoryErrorBase(f"revision packet prior identity mismatch for {task_id}: {key}")
        prior_id = revision.get("prior_transaction_id")
        if not isinstance(prior_id, str) or not re.fullmatch(r"[A-Za-z0-9._-]+", prior_id):
            raise MemoryErrorBase(f"revision packet has invalid prior transaction identity: {task_id}")
        for key in (
            "prior_attestation_sha256", "prior_packet_sha256", "prior_retry_selector_sha256",
            "prior_prerequisites_sha256", "prior_review_packet_sha256",
        ):
            if not isinstance(revision.get(key), str) or not re.fullmatch(r"[0-9a-f]{64}", revision[key]):
                raise MemoryErrorBase(f"revision packet lacks verified {key}: {task_id}")
        if revision["prior_retry_selector_sha256"] != revision["retry_selector_sha256"]:
            raise MemoryErrorBase(f"revision retry selector changed between transactions: {task_id}")
        if revision.get("prior_runtime_profile") not in ("codex", "codex-candidate-reuse"):
            raise MemoryErrorBase(f"revision packet has an ineligible prior runtime profile: {task_id}")
        prior_target = revision.get("prior_target_commit")
        if not isinstance(prior_target, str) or not re.fullmatch(r"[0-9a-f]{40}", prior_target):
            raise MemoryErrorBase(f"revision packet has invalid prior target commit: {task_id}")
        if prior_target != manifest["target_commit"]:
            if task.get("source_unit_id") is None:
                raise MemoryErrorBase(f"derived revision target commit changed: {task_id}")
            if revision["prior_prerequisites_sha256"] != revision["current_prerequisites_sha256"]:
                raise MemoryErrorBase(f"cross-commit source revision unit digest changed: {task_id}")
        prior_review = json.loads(review_attestation_path.read_text(encoding="utf-8"))
        expected_review_verdict = "PASS" if lint_revision else "FAIL"
        if (
            prior_review.get("verdict") != expected_review_verdict
            or prior_review.get("runtime_profile") != "grok-review"
        ):
            raise MemoryErrorBase(
                f"revision packet does not contain the required Grok {expected_review_verdict}: {task_id}"
            )
        allowed_revision_keys = set(expected_revision) | set(prior_expected) | {
            "prior_transaction_id", "prior_target_commit", "prior_attestation_sha256", "prior_packet_sha256",
            "prior_retry_selector_sha256", "prior_prerequisites_sha256", "prior_review_packet_sha256",
            "prior_runtime_profile",
        }
        if lint_revision:
            lint_record = json.loads(lint_record_path.read_text(encoding="utf-8"))
            lint_errors = lint_record.get("errors")
            if (
                not isinstance(lint_errors, list) or not lint_errors
                or not all(
                    isinstance(error, str) and error.startswith(f"{task['output_repository_path']}:")
                    for error in lint_errors
                )
            ):
                raise MemoryErrorBase(f"lint revision packet has ambiguous or unscoped errors: {task_id}")
            lint_expected = {
                "record_version": 1, "role": "machine_lint_failure", "verdict": "FAIL",
                "transaction_id": revision["prior_transaction_id"],
                "target_commit": revision["prior_target_commit"],
                "policy_sha256": lint_record.get("policy_sha256"),
                "policy_tool_sha256": lint_record.get("policy_tool_sha256"),
                "memory_tool_sha256": sha256_file(Path(__file__).resolve()),
                "task_id": task_id, "task_kind": task["task_kind"],
                "source_unit_id": task["source_unit_id"],
                "output_repository_path": task["output_repository_path"],
                "staged_output_path": task["staged_output_path"],
                "candidate_sha256": revision["candidate_sha256"],
                "writer_attestation_sha256": revision["prior_attestation_sha256"],
                "review_attestation_sha256": revision["prior_review_attestation_sha256"],
                "errors": lint_errors, "errors_sha256": sha256_bytes(canonical_json(lint_errors)),
                "warnings": lint_record.get("warnings"),
                "report_path": f"lint-failures/{task_id}/report.md",
                "report_sha256": sha256_file(lint_report_path),
            }
            if lint_record != lint_expected or lint_report_path.read_bytes() != _render_staged_lint_failure_report(lint_record):
                raise MemoryErrorBase(f"lint revision packet record/report identity mismatch: {task_id}")
            lint_hashes = {
                "prior_lint_record_sha256": sha256_file(lint_record_path),
                "prior_lint_report_sha256": sha256_file(lint_report_path),
                "prior_lint_errors_sha256": lint_record["errors_sha256"],
                "revision_basis": "machine_lint_failure",
            }
            for key, value in lint_hashes.items():
                if revision.get(key) != value:
                    raise MemoryErrorBase(f"lint revision packet identity mismatch for {task_id}: {key}")
            allowed_revision_keys |= set(lint_hashes)
        elif revision_basis is not None:
            raise MemoryErrorBase(f"revision packet has unknown revision basis: {task_id}")
        if set(revision) != allowed_revision_keys:
            raise MemoryErrorBase(f"revision packet has undeclared metadata: {task_id}")
    return current_seal


def load_transaction(repo: Path, transaction_arg: str) -> tuple[Path, dict[str, Any]]:
    candidate = Path(transaction_arg)
    if candidate.is_absolute():
        transaction = candidate.resolve()
    elif "/" in transaction_arg:
        transaction = (repo / candidate).resolve()
    else:
        transaction = (repo / "memory/transactions" / transaction_arg).resolve()
    root = (repo / "memory/transactions").resolve()
    if root not in transaction.parents:
        raise MemoryErrorBase("transaction must be under memory/transactions")
    manifest_path = transaction / "manifest.json"
    seal_path = transaction / "seal.json"
    if not manifest_path.is_file() or not seal_path.is_file():
        raise MemoryErrorBase(f"incomplete transaction: {transaction}")
    manifest_data = manifest_path.read_bytes()
    manifest = json.loads(manifest_data)
    seal = json.loads(seal_path.read_bytes())
    if sha256_bytes(manifest_data) != seal.get("manifest_sha256"):
        raise MemoryErrorBase("transaction manifest changed after prepare")
    expected_files = seal.get("snapshot_files", {})
    actual_files: dict[str, str] = {}
    for frozen_root_name in ("snapshot", "policy"):
        frozen_root = transaction / frozen_root_name
        for path in sorted(p for p in frozen_root.rglob("*") if p.is_file()):
            rel = path.relative_to(transaction).as_posix()
            actual_files[rel] = sha256_file(path)
    if actual_files != expected_files:
        raise MemoryErrorBase("transaction snapshot changed after prepare")
    anchored_tasks = set(seal.get("initial_task_packets", {}))
    required_tasks = {item["task_id"] for item in manifest.get("writer_tasks", []) if item.get("required")}
    if anchored_tasks != required_tasks:
        raise MemoryErrorBase("transaction task-packet anchors do not match required writer tasks")
    for task in manifest.get("writer_tasks", []):
        if task.get("required"):
            verify_task_packet_chain(transaction, manifest, seal, task)
    return transaction, manifest


def parse_frontmatter(data: bytes, label: str) -> tuple[dict[str, Any], str]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise LintFailure([f"{label}: page is not UTF-8"]) from exc
    lines = text.splitlines(keepends=True)
    if not lines or lines[0].strip() != "---":
        raise LintFailure([f"{label}: missing YAML frontmatter"])
    ending = None
    for index in range(1, len(lines)):
        if lines[index].strip() == "---":
            ending = index
            break
    if ending is None:
        raise LintFailure([f"{label}: unterminated YAML frontmatter"])
    yaml_data = "".join(lines[1:ending]).encode("utf-8")
    try:
        front = load_yaml_bytes(yaml_data, label)
    except ConfigError as exc:
        raise LintFailure([str(exc)]) from exc
    if not isinstance(front, dict):
        raise LintFailure([f"{label}: frontmatter must be a mapping"])
    return front, "".join(lines[ending + 1 :])


def _blob_at(repo: Path, commit: str, path: str) -> bytes | None:
    proc = subprocess.run(
        ["git", "-C", str(repo), "show", f"{commit}:{path}"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    return proc.stdout if proc.returncode == 0 else None


def _validate_anchor(source_data: bytes, description: str) -> bool:
    text = source_data.decode("utf-8", "replace")
    if description.strip().strip("`") == "anchor-unavailable":
        return True
    tex = re.search(r"`\\label\{([^}]+)\}`", description)
    if tex:
        return f"\\label{{{tex.group(1)}}}" in text
    heading = re.search(r"heading\s+`((?:\\`|[^`])+)`", description)
    if heading:
        sought = heading.group(1).replace("\\`", "`").strip()
        return any(line.strip() == sought for line in text.splitlines())
    explicit = re.search(r"`\{#([^}]+)\}`", description)
    if explicit:
        return f"{{#{explicit.group(1)}}}" in text
    exact_text = re.search(r"(?:exact|distinctive)\s+text\s+`((?:\\`|[^`])+)`", description)
    if exact_text:
        sought = exact_text.group(1).replace("\\`", "`")
        markdown_plain = text.replace("**", "").replace("__", "")
        return sought in text or sought in markdown_plain
    function = re.search(r"function\s+`([^`]+)`", description)
    if function:
        name = re.escape(function.group(1))
        return re.search(rf"(?:def|function)\s+{name}\b|\b{name}\s*\[", text) is not None
    quoted = re.search(r"`([^`]+)`", description)
    if quoted:
        return quoted.group(1) in text
    return False


def _migration_fragment_matches(source_data: bytes, description: str, fragment: str) -> bool:
    if fragment.casefold() in description.casefold():
        return True
    text = source_data.decode("utf-8", "replace")
    described_heading = re.search(r"heading\s+`((?:\\`|[^`])+)`", description)
    if described_heading:
        wanted = fragment.casefold()
        heading_text = described_heading.group(1).replace("\\`", "`").lstrip("#").strip()
        slug = re.sub(r"[^a-z0-9\s-]", "", heading_text.casefold())
        slug = re.sub(r"[\s-]+", "-", slug).strip("-")
        if slug == wanted:
            return True
    return False


def _validate_regions(body: str, label: str, owner: str, errors: list[str]) -> None:
    stack: str | None = None
    seen: set[str] = set()
    markers = 0
    for line_no, line in enumerate(body.splitlines(), 1):
        begin = BEGIN_RE.fullmatch(line.strip())
        end = END_RE.fullmatch(line.strip())
        if begin:
            markers += 1
            region = begin.group(1)
            if stack is not None:
                errors.append(f"{label}:{line_no}: nested generated region")
            if region in seen:
                errors.append(f"{label}:{line_no}: duplicate generated region {region}")
            seen.add(region)
            stack = region
        elif end:
            markers += 1
            if stack != end.group(1):
                errors.append(f"{label}:{line_no}: unbalanced generated region {end.group(1)}")
            stack = None
    if stack is not None:
        errors.append(f"{label}: unclosed generated region {stack}")
    if owner == "mixed" and not markers:
        errors.append(f"{label}: mixed page has no generated-region markers")
    if owner != "mixed" and markers:
        errors.append(f"{label}: generated-region markers are only valid on mixed pages")


def _mask_generated_regions(body: str, declared_regions: Sequence[str], label: str) -> tuple[str | None, list[str]]:
    errors: list[str] = []
    declared = set(declared_regions)
    seen: set[str] = set()
    active: str | None = None
    output: list[str] = []
    for line in body.splitlines(keepends=True):
        stripped = line.strip()
        begin = BEGIN_RE.fullmatch(stripped)
        end = END_RE.fullmatch(stripped)
        if begin:
            region = begin.group(1)
            if active is not None:
                errors.append(f"{label}: nested generated region {region}")
            if region not in declared:
                errors.append(f"{label}: undeclared generated region {region}")
            if region in seen:
                errors.append(f"{label}: duplicate generated region {region}")
            seen.add(region)
            active = region
            output.append(line)
            output.append(f"<GENERATED-CONTENT:{region}>\n")
        elif end:
            region = end.group(1)
            if active != region:
                errors.append(f"{label}: mismatched generated-region end {region}")
            active = None
            output.append(line)
        elif active is None:
            output.append(line)
    if active is not None:
        errors.append(f"{label}: unclosed generated region {active}")
    missing = declared - seen
    if missing:
        errors.append(f"{label}: missing declared generated regions {sorted(missing)}")
    return (None if errors else "".join(output)), errors


def _generated_region_content(body: str, region: str) -> str | None:
    """Return only a declared generated region's body, excluding marker lines."""
    active = False
    content: list[str] = []
    for line in body.splitlines(keepends=True):
        begin = BEGIN_RE.fullmatch(line.strip())
        end = END_RE.fullmatch(line.strip())
        if begin and begin.group(1) == region:
            if active:
                return None
            active = True
            continue
        if end and end.group(1) == region:
            return "".join(content) if active else None
        if active:
            content.append(line)
    return None


def _word_count(text: str) -> int:
    return len(re.findall(r"\b[\w'-]+\b", text, re.UNICODE))


def _substantive_statement_blocks(body: str) -> list[tuple[str, str]]:
    headings = list(re.finditer(r"^###\s+([a-z0-9][a-z0-9-]*)\s+—[^\n]*$", body, re.MULTILINE))
    all_boundaries = list(re.finditer(r"^#{1,3}\s+", body, re.MULTILINE))
    blocks: list[tuple[str, str]] = []
    for heading in headings:
        end = next((item.start() for item in all_boundaries if item.start() > heading.start()), len(body))
        blocks.append((heading.group(1), body[heading.start():end]))
    return blocks


def _unstructured_required_capsule_sections(body: str) -> list[str]:
    """Return capsule sections containing prose outside structured H3 blocks."""
    required = (
        "Purpose and scope",
        "Computed evidence represented by the source",
        "Assumptions, exclusions, and open questions",
        "Revision and supersession relationships",
    )
    failures: list[str] = []
    for title in required:
        heading = re.search(rf"^##\s+{re.escape(title)}\s*$", body, re.MULTILINE)
        if heading is None:
            continue
        start = heading.end()
        following = re.search(r"^##\s+", body[start:], re.MULTILINE)
        end = start + following.start() if following else len(body)
        section = body[start:end]
        structured = re.sub(
            r"^###\s+[a-z0-9][a-z0-9-]*\s+—[^\n]*\n.*?(?=^###\s+|\Z)",
            "",
            section,
            flags=re.MULTILINE | re.DOTALL,
        )
        if structured.strip():
            failures.append(title)
    return failures


def lint_effective_ids(pages: Mapping[str, bytes]) -> list[str]:
    errors: list[str] = []
    registries: dict[str, dict[str, str]] = {"page": {}, "statement": {}, "conflict": {}}
    for path, data in sorted(pages.items()):
        try:
            front, body = parse_frontmatter(data, path)
        except LintFailure:
            continue
        values: list[tuple[str, str]] = []
        if isinstance(front.get("id"), str):
            values.append(("page", front["id"]))
        values.extend(("statement", item[0]) for item in _substantive_statement_blocks(body))
        if front.get("type") == "conflict_register":
            values.extend(
                ("conflict", match.group(1))
                for match in re.finditer(r"^##\s+([a-z0-9][a-z0-9-]*)\s+—", body, re.MULTILINE)
            )
        for kind, stable_id in values:
            prior = registries[kind].get(stable_id)
            if prior is not None and prior != path:
                errors.append(f"duplicate effective {kind} id {stable_id}: {prior}, {path}")
            else:
                registries[kind][stable_id] = path
    return errors


def validate_mixed_preservation(
    base_data: bytes, output_data: bytes, *, declared_regions: Sequence[str], delegated_frontmatter: Sequence[str], label: str,
) -> list[str]:
    errors: list[str] = []
    try:
        base_front, base_body = parse_frontmatter(base_data, label + " base")
        output_front, output_body = parse_frontmatter(output_data, label + " output")
    except LintFailure as exc:
        return exc.errors
    if base_front.get("content_owner") != "mixed" or output_front.get("content_owner") != "mixed":
        errors.append(f"{label}: mixed preservation requires mixed base and output")
    delegated = set(delegated_frontmatter)
    for key in sorted(set(base_front) | set(output_front)):
        if key not in delegated and base_front.get(key) != output_front.get(key):
            errors.append(f"{label}: nondelegated frontmatter changed: {key}")
    masked_base, base_errors = _mask_generated_regions(base_body, declared_regions, label + " base")
    masked_output, output_errors = _mask_generated_regions(output_body, declared_regions, label + " output")
    errors.extend(base_errors)
    errors.extend(output_errors)
    if masked_base is not None and masked_output is not None and masked_base != masked_output:
        errors.append(f"{label}: unmarked mixed-page body bytes changed")
    return errors


def lint_pages(
    repo: Path,
    pages: Mapping[str, bytes],
    *,
    commit_for_page: Mapping[str, str],
    expected_units: Mapping[str, Mapping[str, Any]],
    known_paths: set[str],
    known_pages: Mapping[str, bytes] | None = None,
) -> tuple[list[str], list[str], dict[str, dict[str, Any]]]:
    errors: list[str] = []
    warnings: list[str] = []
    parsed: dict[str, dict[str, Any]] = {}
    page_ids: dict[str, str] = {}
    statement_ids: dict[str, str] = {}
    conflict_ids: dict[str, str] = {}
    for path, data in sorted(pages.items()):
        try:
            front, body = parse_frontmatter(data, path)
        except LintFailure as exc:
            errors.extend(exc.errors)
            continue
        parsed[path] = front
        required = ("schema_version", "id", "title", "type", "lifecycle", "memory_review", "sources", "content_owner", "last_updated")
        for field in required:
            if field not in front:
                errors.append(f"{path}: missing frontmatter field {field}")
        page_id = front.get("id")
        if not isinstance(page_id, str) or not ID_RE.fullmatch(page_id):
            errors.append(f"{path}: invalid page id")
        elif page_id in page_ids:
            errors.append(f"duplicate page id {page_id}: {page_ids[page_id]}, {path}")
        else:
            page_ids[page_id] = path
        if front.get("schema_version") != 2:
            errors.append(f"{path}: schema_version must be 2")
        if front.get("type") not in PAGE_TYPES:
            errors.append(f"{path}: invalid page type")
        if front.get("lifecycle") not in LIFECYCLES:
            errors.append(f"{path}: invalid lifecycle")
        if front.get("memory_review") not in REVIEWS:
            errors.append(f"{path}: invalid memory_review")
        if front.get("content_owner") not in OWNERS:
            errors.append(f"{path}: invalid content_owner")
        if front.get("lifecycle") == "superseded" and not front.get("superseded_by"):
            errors.append(f"{path}: superseded page lacks superseded_by")
        if front.get("memory_review") == "reviewed":
            for field in ("reviewed_by", "reviewed_at", "review_record"):
                if not front.get(field):
                    errors.append(f"{path}: reviewed page lacks {field}")
        owner = front.get("content_owner")
        if isinstance(owner, str):
            _validate_regions(body, path, owner, errors)
        if owner in ("ai_generated", "mixed") and re.search(r"`(?:atlas|graph)/[^`]+`", body):
            errors.append(f"{path}: generated content cites legacy atlas/ or graph/ material")
        sources = front.get("sources")
        if not isinstance(sources, list) or not all(isinstance(source, str) for source in sources):
            errors.append(f"{path}: sources must be a path list")
            sources = []
        elif sources != sorted(set(sources)):
            errors.append(f"{path}: sources must be sorted and duplicate-free")
        commit = commit_for_page.get(path)
        if not commit:
            errors.append(f"{path}: no applicable committed snapshot")
        for source in sources:
            try:
                clean_repo_path(source, f"{path} source")
            except ConfigError as exc:
                errors.append(str(exc))
                continue
            if commit and _blob_at(repo, commit, source) is None and front.get("lifecycle") != "deleted":
                errors.append(f"{path}: source absent at {commit[:12]}: {source}")
        if owner in ("ai_generated", "mixed"):
            if front.get("generated_from_commit") != commit:
                errors.append(f"{path}: generated_from_commit does not match transaction/state commit")
        for match in re.finditer(r"^###\s+([a-z0-9][a-z0-9-]*)\s+—", body, re.MULTILINE):
            statement = match.group(1)
            if statement in statement_ids:
                errors.append(f"duplicate statement id {statement}: {statement_ids[statement]}, {path}")
            else:
                statement_ids[statement] = path
        if front.get("type") == "conflict_register":
            for match in re.finditer(r"^##\s+([a-z0-9][a-z0-9-]*)\s+—", body, re.MULTILINE):
                conflict_id = match.group(1)
                if conflict_id in conflict_ids:
                    errors.append(f"duplicate conflict id {conflict_id}: {conflict_ids[conflict_id]}, {path}")
                else:
                    conflict_ids[conflict_id] = path
            if re.search(r"Conflict state:\s*`resolved`", body):
                bases = re.findall(r"resolution_basis:\s*`?([a-z_]+)`?", body)
                if not bases or any(basis not in ("explicit_source", "reviewed_adjudication") for basis in bases):
                    errors.append(f"{path}: resolved conflict lacks an allowed resolution_basis")
        for match in STATUS_RE.finditer(body):
            lifecycle, evidence, review = match.groups()
            if lifecycle not in LIFECYCLES:
                errors.append(f"{path}: statement has invalid lifecycle {lifecycle}")
            if evidence not in EVIDENCE:
                errors.append(f"{path}: statement has invalid evidence {evidence}")
            if review not in REVIEWS:
                errors.append(f"{path}: statement has invalid memory_review {review}")
        for source, anchor in CITATION_RE.findall(body):
            if source not in sources:
                errors.append(f"{path}: cited source missing from frontmatter: {source}")
            source_data = _blob_at(repo, commit, source) if commit else None
            if source_data is None:
                if front.get("lifecycle") != "deleted":
                    errors.append(f"{path}: citation source absent: {source}")
            elif not _validate_anchor(source_data, anchor):
                errors.append(f"{path}: citation anchor not found in {source}: {anchor}")
        for raw_target in LINK_RE.findall(body):
            target = raw_target.strip().split()[0].strip("<>")
            if target.startswith(("http://", "https://", "mailto:")):
                continue
            target, fragment_separator, fragment = target.partition("#")
            if not target and fragment_separator:
                resolved_target = path
            elif target.startswith("/"):
                resolved_target = target.lstrip("/")
            else:
                resolved_target = (PurePosixPath(path).parent / PurePosixPath(target)).as_posix()
            normalized = PurePosixPath(resolved_target)
            if ".." in normalized.parts:
                # PurePosixPath does not collapse; resolve safely through filesystem.
                candidate = (repo / path).parent.joinpath(target).resolve()
                try:
                    resolved_target = candidate.relative_to(repo).as_posix()
                except ValueError:
                    errors.append(f"{path}: internal link escapes repository: {raw_target}")
                    continue
            if resolved_target.startswith("memory/") and resolved_target not in known_paths:
                errors.append(f"{path}: unresolved internal link: {raw_target}")
            elif fragment_separator and resolved_target.startswith("memory/") and known_pages is not None:
                linked_data = known_pages.get(resolved_target)
                if linked_data is None:
                    continue
                try:
                    _, linked_body = parse_frontmatter(linked_data, resolved_target)
                except LintFailure:
                    continue
                explicit = f"{{#{fragment}}}" in linked_body
                heading_match = False
                for heading in re.finditer(r"^#{1,6}\s+(.+?)\s*$", linked_body, re.MULTILINE):
                    heading_text = re.sub(r"\s+#+\s*$", "", heading.group(1)).strip()
                    slug = re.sub(r"[^a-z0-9\s-]", "", heading_text.casefold())
                    slug = re.sub(r"[\s-]+", "-", slug).strip("-")
                    if slug == fragment.casefold():
                        heading_match = True
                        break
                if not explicit and not heading_match:
                    errors.append(f"{path}: unresolved internal link fragment: {raw_target}")
        if front.get("type") == "source_capsule":
            source_unit = front.get("source_unit")
            if not isinstance(source_unit, dict):
                errors.append(f"{path}: source capsule lacks source_unit mapping")
            else:
                unit_id = source_unit.get("id")
                expected = expected_units.get(str(unit_id))
                if expected is None:
                    errors.append(f"{path}: source unit {unit_id!r} is not expected")
                else:
                    if path != expected.get("capsule"):
                        errors.append(f"{path}: capsule target disagrees with unit {unit_id}")
                    if source_unit.get("unit_digest_sha256") != expected.get("unit_digest_sha256"):
                        errors.append(f"{path}: source unit digest disagrees with resolved membership")
                    actual_members = source_unit.get("members")
                    expected_members = []
                    for member in expected.get("members", []):
                        expected_members.append({key: member.get(key) for key in (
                            "path", "role", "read_mode", "mode", "object_type", "blob_oid", "blob_size"
                        )})
                    if actual_members != expected_members:
                        errors.append(f"{path}: source unit member identities disagree with resolved membership")
    return errors, warnings, parsed


def served_page_bytes(repo: Path, state: Mapping[str, Any]) -> tuple[dict[str, bytes], list[str]]:
    pages: dict[str, bytes] = {}
    errors: list[str] = []
    for path, record in sorted(state.get("pages", {}).items()):
        disk = repo / path
        if not disk.is_file():
            errors.append(f"served page missing: {path}")
            continue
        data = disk.read_bytes()
        digest = sha256_bytes(data)
        if digest != record.get("sha256"):
            errors.append(f"served page hash is not in successful state: {path}")
            continue
        if record.get("generation") is None:
            errors.append(f"served page has no successful generation: {path}")
            continue
        pages[path] = data
    return pages, errors


def _page_candidates(repo: Path) -> set[str]:
    result: set[str] = set()
    for root in PAGE_ROOTS:
        path = repo / root
        if path.is_file() and path.suffix == ".md":
            result.add(root)
        elif path.is_dir():
            result.update(item.relative_to(repo).as_posix() for item in path.rglob("*.md"))
    return result


def lint_served(repo: Path) -> tuple[list[str], list[str]]:
    config, _ = load_config(repo)
    state, _ = load_state(repo)
    pages, errors = served_page_bytes(repo, state)
    orphans = _page_candidates(repo) - set(state.get("pages", {}))
    errors.extend(f"page is not recorded in successful state: {path}" for path in sorted(orphans))
    configured_derived = {item["id"]: item for item in config.get("derived_pages", [])}
    removed_derived = sorted(set(state.get("derived_pages", {})) - set(configured_derived))
    for derived_id in removed_derived:
        record = state["derived_pages"][derived_id]
        errors.append(
            f"configured derived page was removed without an explicit retirement/cutover task: "
            f"{derived_id} ({record.get('page')})"
        )
    configured_paths = {item["page"] for item in configured_derived.values()}
    for path, record in state.get("pages", {}).items():
        task_id = record.get("task_id")
        if task_id is not None and task_id not in configured_derived and path not in configured_paths:
            errors.append(f"served derived page is orphaned from config: {task_id} ({path})")
    commits: dict[str, str] = {}
    expected: dict[str, Mapping[str, Any]] = {}
    for unit_id, unit in state.get("units", {}).items():
        expected[unit_id] = unit
        capsule = unit.get("capsule")
        if capsule:
            page_record = state.get("pages", {}).get(capsule, {})
            commits[capsule] = page_record.get("generated_commit") or unit.get("processed_commit")
    for unit_id, unit in state.get("retired_units", {}).items():
        expected[unit_id] = unit
        capsule = unit.get("capsule")
        if capsule:
            page_record = state.get("pages", {}).get(capsule, {})
            commits[capsule] = page_record.get("generated_commit") or unit.get("processed_commit")
    for path, record in state.get("pages", {}).items():
        commits.setdefault(
            path,
            record.get("generated_commit") or record.get("processed_commit") or state.get("last_fully_processed_commit"),
        )
    lint_errors, warnings, _ = lint_pages(
        repo, pages, commit_for_page=commits, expected_units=expected, known_paths=set(pages), known_pages=pages
    )
    errors.extend(lint_errors)
    # State coverage is exact for active/retired capsules that have succeeded.
    for unit_id, unit in {**state.get("units", {}), **state.get("retired_units", {})}.items():
        capsule = unit.get("capsule")
        if unit.get("result") == "success" and capsule and capsule not in state.get("pages", {}):
            errors.append(f"successful unit {unit_id} has no served capsule page")
    migration = migration_status_report(repo, config, state)
    for pending in migration["pending"]:
        errors.append(
            f"atlas migration incomplete: {pending['legacy_id']} has not been sealed into {pending['target']}"
        )
    return errors, warnings


def lint_staged(repo: Path, transaction_arg: str) -> tuple[list[str], list[str]]:
    transaction, manifest = load_transaction(repo, transaction_arg)
    expected = {item["id"]: item for item in manifest["units"]}
    pages: dict[str, bytes] = {}
    errors: list[str] = []
    expected_paths = {
        item["output_repository_path"] for item in manifest.get("writer_tasks", []) if item.get("required")
    }
    for path in sorted(expected_paths):
        staged = transaction / "staged" / path
        if not staged.is_file():
            errors.append(f"required staged page missing: {path}")
        else:
            pages[path] = staged.read_bytes()
    actual_paths = {
        item.relative_to(transaction / "staged").as_posix()
        for item in (transaction / "staged").rglob("*") if item.is_file()
    }
    errors.extend(f"unexpected staged file: {path}" for path in sorted(actual_paths - expected_paths))
    commits = {path: manifest["target_commit"] for path in pages}
    state, _ = load_state(repo)
    served, served_errors = served_page_bytes(repo, state)
    # Effective set supports links from staged pages to already served pages.
    effective = dict(served)
    effective.update(pages)
    errors.extend(served_errors)
    errors.extend(lint_effective_ids(effective))
    lint_errors, warnings, parsed = lint_pages(
        repo, pages, commit_for_page=commits, expected_units=expected, known_paths=set(effective), known_pages=effective
    )
    errors.extend(lint_errors)
    writer_by_path = {
        task["output_repository_path"]: task for task in manifest.get("writer_tasks", []) if task.get("required")
    }
    for path, data in pages.items():
        task = writer_by_path.get(path)
        if task is None:
            errors.append(f"{path}: no declared writer task")
            continue
        contract = task.get("semantic_contract") or {}
        page_contract = contract.get("page") or {}
        front = parsed.get(path, {})
        expected_fields = {
            "id": page_contract.get("id"),
            "title": page_contract.get("title"),
            "type": page_contract.get("type"),
            "content_owner": page_contract.get("content_owner"),
            "lifecycle": page_contract.get("desired_lifecycle"),
            "last_updated": contract.get("refresh_date"),
            "source_kind": contract.get("source_kind"),
            "extractor_version": contract.get("extractor_version"),
        }
        for field, value in expected_fields.items():
            actual = front.get(field)
            if field == "last_updated" and actual is not None:
                actual = str(actual)
            if value is not None and actual != value:
                errors.append(f"{path}: manifest-owned field disagrees with task: {field}")
        if task.get("task_kind") == "source_capsule":
            source_unit = front.get("source_unit") if isinstance(front.get("source_unit"), dict) else {}
            for field in ("shape", "entrypoint"):
                if source_unit.get(field) != contract.get(field):
                    errors.append(f"{path}: source_unit.{field} disagrees with task")
        allowed_paths = {
            item["path"] if isinstance(item, dict) else item for item in contract.get("allowed_citations", [])
        }
        if isinstance(front.get("sources"), list):
            disallowed = set(front["sources"]) - allowed_paths
            if disallowed:
                errors.append(f"{path}: sources are outside the task citation allowlist: {sorted(disallowed)}")
        citation_modes: dict[str, str] = {}
        mode_rank = {"identity_only": 0, "excerpt": 1, "semantic": 2}
        for item in contract.get("allowed_citations", []):
            if isinstance(item, dict):
                source_path = item.get("path")
                read_mode = item.get("read_mode", "identity_only")
                if isinstance(source_path, str) and mode_rank.get(read_mode, -1) > mode_rank.get(
                    citation_modes.get(source_path, "identity_only"), -1
                ):
                    citation_modes[source_path] = read_mode
        for item in task.get("semantic_inputs", []):
            source_path = item.get("source_path")
            read_mode = item.get("read_mode") or item.get("member_identity", {}).get("read_mode")
            if isinstance(source_path, str) and read_mode in ("semantic", "excerpt"):
                citation_modes[source_path] = read_mode
        _, page_body = parse_frontmatter(data, path)
        structural_body = page_body
        structural_region = None
        declared_regions = task.get("generated_regions") or []
        if page_contract.get("content_owner") == "mixed" and len(declared_regions) == 1:
            structural_region = _generated_region_content(page_body, declared_regions[0])
            if structural_region is not None:
                structural_body = structural_region
        statement_blocks = _substantive_statement_blocks(structural_body)
        statement_ids = {stable_id for stable_id, _ in statement_blocks}
        if task.get("task_kind") == "source_capsule":
            for section_title in _unstructured_required_capsule_sections(structural_body):
                errors.append(
                    f"{path}: substantive prose outside statement blocks in section: {section_title}"
                )
        for heading in re.finditer(r"^###\s+([^\n]+)$", structural_body, re.MULTILINE):
            if not re.fullmatch(r"[a-z0-9][a-z0-9-]*\s+—\s+.+", heading.group(1).strip()):
                errors.append(f"{path}: generated H3 heading is not a structured substantive statement: {heading.group(1)}")
        namespace = contract.get("stable_id_namespace")
        for stable_id, block in statement_blocks:
            if isinstance(namespace, str) and not stable_id.startswith(namespace):
                errors.append(f"{path}: statement ID is outside page namespace {namespace}: {stable_id}")
            statuses = [
                match.groups() for match in re.finditer(
                    r"^Status:\s*`lifecycle=([^`]+)`\s*·\s*`evidence=([^`]+)`\s*·\s*`memory_review=([^`]+)`\s*$",
                    block, re.MULTILINE,
                )
            ]
            if len(statuses) != 1:
                errors.append(f"{path}: substantive statement {stable_id} must have exactly one status block")
                evidence = None
            else:
                _, evidence, review = statuses[0]
                if review != "ai_draft":
                    errors.append(f"{path}: generated statement {stable_id} must remain memory_review=ai_draft")
                if evidence == "measured":
                    errors.append(
                        f"{path}: AI-generated statement {stable_id} cannot assign measured until "
                        "an explicit evidence-chain schema is implemented"
                    )
            statement_citations = CITATION_RE.findall(block)
            if not statement_citations:
                errors.append(f"{path}: substantive statement {stable_id} has no direct-source citation")
            for source_path, anchor in statement_citations:
                read_mode = citation_modes.get(source_path, "identity_only")
                identity_fact = (
                    read_mode == "identity_only"
                    and anchor.strip().strip("`") == "anchor-unavailable"
                    and evidence in ("provisional", "open")
                    and re.search(
                        r"\bidentity[- ]only\b|\bnot semantically readable\b",
                        block, re.IGNORECASE,
                    ) is not None
                )
                if read_mode not in ("semantic", "excerpt") and not identity_fact:
                    errors.append(
                        f"{path}: statement {stable_id} cites unreadable identity-only source {source_path}"
                    )
                if anchor.strip().strip("`") == "anchor-unavailable" and evidence in ("measured", "derived"):
                    errors.append(
                        f"{path}: statement {stable_id} cannot use anchor-unavailable with evidence={evidence}"
                    )
        for prior_id in contract.get("prior_stable_ids", {}).get("statements", []):
            if prior_id not in statement_ids and f"`{prior_id}`" not in structural_body:
                errors.append(f"{path}: prior statement ID is not retained or explicitly accounted for: {prior_id}")
        task_budget = contract.get("budget") or contract.get("budgets", {})
        output_limit = task_budget.get("output_max_bytes")
        if isinstance(output_limit, int) and len(data) > output_limit:
            errors.append(f"{path}: output exceeds its {output_limit}-byte task budget")
        budget_body = page_body
        regions = declared_regions
        if page_contract.get("content_owner") == "mixed" and len(regions) == 1:
            region_body = _generated_region_content(page_body, regions[0])
            if region_body is None:
                errors.append(f"{path}: cannot locate generated region for budget enforcement: {regions[0]}")
            else:
                budget_body = region_body
        max_words = task_budget.get("max_words")
        if isinstance(max_words, int):
            word_count = _word_count(budget_body)
            if word_count > max_words:
                errors.append(f"{path}: generated output has {word_count} words, above its {max_words}-word task budget")
        max_statements = task_budget.get("max_key_statements")
        if isinstance(max_statements, int):
            statement_count = len(re.findall(r"^###\s+[a-z0-9][a-z0-9-]*\s+—", budget_body, re.MULTILINE))
            if statement_count > max_statements:
                errors.append(
                    f"{path}: generated output has {statement_count} key statements, "
                    f"above its {max_statements}-statement task budget"
                )
        entry_matches = list(re.finditer(r"^##\s+.+$", budget_body, re.MULTILINE))
        max_entries = task_budget.get("max_entries")
        if isinstance(max_entries, int) and len(entry_matches) > max_entries:
            errors.append(
                f"{path}: generated output has {len(entry_matches)} entries, above its {max_entries}-entry task budget"
            )
        max_entry_words = task_budget.get("max_entry_words")
        if isinstance(max_entry_words, int):
            for position, match in enumerate(entry_matches):
                end = entry_matches[position + 1].start() if position + 1 < len(entry_matches) else len(budget_body)
                entry_words = _word_count(budget_body[match.start():end])
                if entry_words > max_entry_words:
                    heading = match.group(0).strip()
                    errors.append(
                        f"{path}: entry {heading!r} has {entry_words} words, "
                        f"above its {max_entry_words}-word entry budget"
                    )
        migration_body = structural_body
        migration_headings = list(re.finditer(r"^#### Migration record — (.+?)\s*$", migration_body, re.MULTILINE))
        migration_records: dict[str, str] = {}
        for position, heading in enumerate(migration_headings):
            end = migration_headings[position + 1].start() if position + 1 < len(migration_headings) else len(migration_body)
            migration_records[heading.group(1).strip()] = migration_body[heading.start():end]
        for migration in contract.get("migration_requirements", []):
            legacy_id = migration["legacy_id"]
            record_body = migration_records.get(legacy_id)
            if record_body is None:
                errors.append(f"{path}: migration item lacks a structured generated-region record: {legacy_id}")
                continue
            if not re.search(rf"^- `legacy_id`: `{re.escape(legacy_id)}`\s*$", record_body, re.MULTILINE):
                errors.append(f"{path}: migration record legacy_id mismatch: {legacy_id}")
            disposition_match = re.search(
                r"^- `migration_disposition`: `(migrated|blocked)`\s*$", record_body, re.MULTILINE
            )
            disposition = disposition_match.group(1) if disposition_match else None
            if disposition is None:
                errors.append(f"{path}: migration record lacks a valid disposition: {legacy_id}")
            inventory_digest = manifest["policy"].get("atlas_migration_sha256")
            if not re.search(
                rf"^- `inventory_sha256`: `{re.escape(str(inventory_digest))}`\s*$", record_body, re.MULTILINE
            ):
                errors.append(f"{path}: migration record inventory digest mismatch: {legacy_id}")
            targets_line = re.search(r"^- `target_statement_ids`: \[(.*?)\]\s*$", record_body, re.MULTILINE)
            target_ids = re.findall(r"`([a-z0-9][a-z0-9-]*)`", targets_line.group(1)) if targets_line else []
            if disposition == "migrated" and not target_ids:
                errors.append(f"{path}: migrated record has no target statement IDs: {legacy_id}")
            for target_id in target_ids:
                if target_id not in statement_ids or (isinstance(namespace, str) and not target_id.startswith(namespace)):
                    errors.append(f"{path}: migration record points to invalid target statement ID {target_id}")
            citation_pairs = CITATION_RE.findall(record_body)
            for reference in migration["required_original_refs"]:
                source_path, separator, fragment = reference.partition("#")
                matching_anchors = [anchor for source, anchor in citation_pairs if source == source_path]
                if not matching_anchors:
                    if disposition != "blocked" or reference not in record_body:
                        errors.append(
                            f"{path}: migration item {legacy_id} lacks bound original-source citation {source_path}"
                        )
                elif separator:
                    source_data = _blob_at(repo, manifest["target_commit"], source_path)
                    if source_data is not None and any(
                        _migration_fragment_matches(source_data, anchor, fragment)
                        for anchor in matching_anchors
                    ):
                        continue
                    errors.append(
                        f"{path}: migration item {legacy_id} lacks required original-source anchor {reference}"
                    )
        if page_contract.get("content_owner") == "mixed":
            base = transaction / task["packet_path"] / "base_page.md"
            if not base.is_file():
                errors.append(f"{path}: mixed task lacks sealed base_page.md")
            else:
                errors.extend(validate_mixed_preservation(
                    base.read_bytes(), data, declared_regions=regions,
                    delegated_frontmatter=task.get("delegated_frontmatter", []), label=path,
                ))
    return errors, warnings


def _render_staged_lint_failure_report(record: Mapping[str, Any]) -> bytes:
    lines = [
        "Verdict: FAIL", "", "Machine: staged memory lint", "",
        f"Transaction: `{record['transaction_id']}`", "", f"Task: `{record['task_id']}`", "",
        f"Output: `{record['output_repository_path']}`", "",
        f"Candidate SHA-256: `{record['candidate_sha256']}`", "", "## Exact errors", "",
    ]
    lines.extend(f"{index}. {error}" for index, error in enumerate(record["errors"], 1))
    return ("\n".join(lines).rstrip() + "\n").encode("utf-8")


def record_staged_lint_failure(
    repo: Path, transaction_arg: str, task_id: str | None = None,
) -> dict[str, Any]:
    """Record one task-scoped staged lint failure after an exact current Grok PASS."""
    transaction, manifest = load_transaction(repo, transaction_arg)
    _recheck_transaction(repo, manifest)
    required = [task for task in manifest.get("writer_tasks", []) if task.get("required")]
    if task_id is None:
        if len(required) != 1:
            raise MemoryErrorBase("recorded staged lint failure requires --task when multiple writers are required")
        task = required[0]
        task_id = task["task_id"]
    else:
        task = next((item for item in required if item.get("task_id") == task_id), None)
        if task is None:
            raise MemoryErrorBase(f"recorded staged lint task is not a required writer: {task_id}")
    output_path = task["output_repository_path"]
    candidate = transaction / task["staged_output_path"]
    if not candidate.is_file() or candidate.is_symlink():
        raise MemoryErrorBase(f"cannot record lint failure without a regular staged candidate: {task_id}")
    candidate_sha256 = sha256_file(candidate)
    staged_hashes: dict[str, str] = {}
    for required_task in required:
        staged = transaction / required_task["staged_output_path"]
        if not staged.is_file() or staged.is_symlink():
            raise MemoryErrorBase(
                f"cannot record lint failure while required candidate is missing: {required_task['task_id']}"
            )
        staged_hashes[required_task["output_repository_path"]] = sha256_file(staged)
    attestation_errors = verify_isolation_attestations(
        repo, transaction, manifest, staged_hashes,
    )
    if attestation_errors:
        raise MemoryErrorBase(
            "cannot record lint failure before writer provenance and current Grok PASS verify:\n"
            + "\n".join(f"- {error}" for error in attestation_errors)
        )
    errors, warnings = lint_staged(repo, str(transaction))
    if not errors:
        raise MemoryErrorBase("refusing to record staged lint failure because staged lint passed")
    prefix = f"{output_path}:"
    unscoped = [error for error in errors if not error.startswith(prefix)]
    if unscoped:
        raise MemoryErrorBase(
            "refusing to record ambiguous or unscoped staged lint errors:\n"
            + "\n".join(f"- {error}" for error in unscoped)
        )
    if sha256_file(candidate) != candidate_sha256:
        raise MemoryErrorBase("staged candidate changed while lint failure was being recorded")
    record_root = transaction / "lint-failures" / task_id
    if record_root.exists() or record_root.is_symlink():
        raise MemoryErrorBase(f"recorded staged lint failure already exists for {task_id}")
    report_relative = f"lint-failures/{task_id}/report.md"
    record = {
        "record_version": 1,
        "role": "machine_lint_failure",
        "verdict": "FAIL",
        "transaction_id": manifest["transaction_id"],
        "target_commit": manifest["target_commit"],
        "policy_sha256": manifest["policy"]["combined_sha256"],
        "policy_tool_sha256": manifest["policy"]["tool_sha256"],
        "memory_tool_sha256": sha256_file(Path(__file__).resolve()),
        "task_id": task_id,
        "task_kind": task["task_kind"],
        "source_unit_id": task["source_unit_id"],
        "output_repository_path": output_path,
        "staged_output_path": task["staged_output_path"],
        "candidate_sha256": candidate_sha256,
        "writer_attestation_sha256": sha256_file(transaction / "attestations" / f"{task_id}.json"),
        "review_attestation_sha256": sha256_file(transaction / "reviews" / task_id / "attestation.json"),
        "errors": errors,
        "errors_sha256": sha256_bytes(canonical_json(errors)),
        "warnings": warnings,
        "report_path": report_relative,
        "report_sha256": None,
    }
    report_data = _render_staged_lint_failure_report(record)
    record["report_sha256"] = sha256_bytes(report_data)
    atomic_write(transaction / report_relative, report_data, 0o444)
    write_json(record_root / "record.json", record)
    return record


def _reuse_source_identity(manifest: Mapping[str, Any], source_unit_id: str | None) -> tuple[dict[str, Any], dict[str, Any]]:
    unit = next((item for item in manifest.get("units", []) if item.get("id") == source_unit_id), None)
    if source_unit_id is None or not isinstance(unit, dict):
        raise MemoryErrorBase("reviewed reuse is only supported for source-unit tasks")

    def identities(key: str) -> list[dict[str, Any]]:
        return [
            {field: value for field, value in member.items() if field != "snapshot_path"}
            for member in unit.get(key, [])
        ]

    return unit, {"members": identities("members"), "prior_members": identities("prior_members")}


def verify_reviewed_reuse_chain(
    repo: Path,
    transaction: Path,
    manifest: Mapping[str, Any],
    task: Mapping[str, Any],
    attestation: Mapping[str, Any],
    staged_sha256: str | None,
) -> list[str]:
    """Independently verify the prior Codex -> Grok PASS chain for deterministic reuse."""
    task_id = task["task_id"]
    label = f"reviewed reuse for {task_id}"
    try:
        if task.get("source_unit_id") is None:
            raise MemoryErrorBase("reviewed reuse is only supported for source-unit tasks")
        if attestation.get("model_invoked") is not False:
            raise MemoryErrorBase("reviewed reuse must explicitly record model_invoked=false")
        if attestation.get("runtime_version") is not None or attestation.get("runtime_executable_sha256") is not None:
            raise MemoryErrorBase("reviewed reuse must not claim a runtime invocation")
        chain = attestation.get("reviewed_reuse")
        if not isinstance(chain, dict):
            raise MemoryErrorBase("reviewed reuse attestation lacks its provenance chain")
        prior_id = chain.get("prior_transaction_id")
        if not isinstance(prior_id, str) or not re.fullmatch(r"[A-Za-z0-9._-]+", prior_id):
            raise MemoryErrorBase("reviewed reuse has an invalid prior transaction identity")
        prior_transaction, prior_manifest = load_transaction(repo, prior_id)
        if prior_transaction == transaction:
            raise MemoryErrorBase("reviewed reuse must come from another transaction")
        prior_task = next(
            (item for item in prior_manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None,
        )
        if not isinstance(prior_task, dict) or not prior_task.get("required"):
            raise MemoryErrorBase("reviewed reuse prior writer task is missing")
        if (
            prior_task.get("task_kind") != task.get("task_kind")
            or prior_task.get("source_unit_id") != task.get("source_unit_id")
            or prior_task.get("output_repository_path") != task.get("output_repository_path")
            or prior_task.get("staged_output_path") != task.get("staged_output_path")
            or prior_task.get("expected_output_name") != task.get("expected_output_name")
        ):
            raise MemoryErrorBase("reviewed reuse prior task/output identity mismatch")
        current_unit, current_identity = _reuse_source_identity(manifest, task.get("source_unit_id"))
        prior_unit, prior_identity = _reuse_source_identity(prior_manifest, prior_task.get("source_unit_id"))
        if current_identity != prior_identity or current_unit.get("kind") != prior_unit.get("kind"):
            raise MemoryErrorBase("reviewed reuse source member/kind identity mismatch")
        if (
            prior_manifest.get("target_commit") != manifest.get("target_commit")
            and prior_unit.get("unit_digest_sha256") != current_unit.get("unit_digest_sha256")
        ):
            raise MemoryErrorBase("cross-commit reviewed reuse source-unit digest changed")

        prior_packet = prior_transaction / prior_task["packet_path"]
        prior_seal = json.loads((prior_packet.parent / "packet-seal.json").read_text(encoding="utf-8"))
        writer_attestation_path = prior_transaction / "attestations" / f"{task_id}.json"
        if not writer_attestation_path.is_file() or writer_attestation_path.is_symlink():
            raise MemoryErrorBase("reviewed reuse prior Codex attestation is missing")
        writer_attestation = json.loads(writer_attestation_path.read_text(encoding="utf-8"))
        expected_writer = {
            "transaction_id": prior_manifest["transaction_id"],
            "task_id": task_id,
            "source_unit_id": prior_task["source_unit_id"],
            "isolation": "bubblewrap",
            "workspace_hidden": repo.as_posix(),
            "packet_path": prior_task["packet_path"],
            "packet_sha256": prior_seal["combined_sha256"],
            "output_repository_path": prior_task["output_repository_path"],
            "staged_output_path": prior_task["staged_output_path"],
            "runtime_profile": "codex",
        }
        for key, value in expected_writer.items():
            if writer_attestation.get(key) != value:
                raise MemoryErrorBase(f"prior Codex attestation mismatch: {key}")
        for key in ("runner_sha256", "runtime_executable_sha256", "output_sha256"):
            if not isinstance(writer_attestation.get(key), str) or not re.fullmatch(
                r"[0-9a-f]{64}", writer_attestation[key]
            ):
                raise MemoryErrorBase(f"prior Codex attestation lacks {key}")
        if not isinstance(writer_attestation.get("runtime_version"), str) or not writer_attestation["runtime_version"]:
            raise MemoryErrorBase("prior Codex attestation lacks runtime version")
        candidate = prior_transaction / prior_task["staged_output_path"]
        if not candidate.is_file() or candidate.is_symlink():
            raise MemoryErrorBase("reviewed reuse prior candidate is missing")
        candidate_data = candidate.read_bytes()
        candidate_sha256 = sha256_bytes(candidate_data)
        if candidate_sha256 != writer_attestation["output_sha256"]:
            raise MemoryErrorBase("reviewed reuse prior candidate hash mismatch")

        review_root = prior_transaction / "reviews" / task_id
        review_attestation_path = review_root / "attestation.json"
        if not review_attestation_path.is_file() or review_attestation_path.is_symlink():
            raise MemoryErrorBase("reviewed reuse Grok attestation is missing")
        review = json.loads(review_attestation_path.read_text(encoding="utf-8"))
        report_relative = f"reviews/{task_id}/output/report.md"
        expected_review = {
            "role": "independent_review",
            "transaction_id": prior_manifest["transaction_id"],
            "task_id": task_id,
            "packet_sha256": prior_seal["combined_sha256"],
            "candidate_path": prior_task["staged_output_path"],
            "candidate_sha256": candidate_sha256,
            "writer_attestation_sha256": sha256_file(writer_attestation_path),
            "report_path": report_relative,
            "verdict": "PASS",
            "runtime_profile": "grok-review",
            "review_model": GROK_REVIEW_MODEL,
            "review_prompt_sha256": REVIEW_PROMPT_SHA256,
            "review_schema_sha256": REVIEW_SCHEMA_SHA256,
            "review_contract_sha256": REVIEW_CONTRACT_SHA256,
            "reviewer_sha256": sha256_file(Path(__file__).resolve().parent / "review_isolated.py"),
        }
        for key, value in expected_review.items():
            if review.get(key) != value:
                raise MemoryErrorBase(f"prior Grok PASS attestation mismatch: {key}")
        for key in ("review_packet_sha256", "report_sha256", "runtime_executable_sha256"):
            if not isinstance(review.get(key), str) or not re.fullmatch(r"[0-9a-f]{64}", review[key]):
                raise MemoryErrorBase(f"prior Grok attestation lacks {key}")
        if not isinstance(review.get("runtime_version"), str) or not review["runtime_version"]:
            raise MemoryErrorBase("prior Grok attestation lacks runtime version")
        report = prior_transaction / report_relative
        if not report.is_file() or report.is_symlink() or sha256_file(report) != review["report_sha256"]:
            raise MemoryErrorBase("reviewed reuse Grok report hash mismatch")
        review_packet = review_root / "packet"
        review_files = packet_file_hashes(review_packet)
        review_sha256 = sha256_bytes(canonical_json(review_files))
        if (
            review_files != review.get("review_packet_files")
            or review_sha256 != review["review_packet_sha256"]
            or review_files.get("candidate.md") != candidate_sha256
        ):
            raise MemoryErrorBase("reviewed reuse Grok packet hash mismatch")
        expected_review_files = dict(prior_seal["files"])
        expected_review_files.update({
            "candidate.md": candidate_sha256,
            "review-prompt.md": review_files.get("review-prompt.md"),
        })
        if review_files != expected_review_files:
            raise MemoryErrorBase("reviewed reuse Grok packet is not the prior writer packet plus candidate")

        expected_chain = {
            "prior_transaction_id": prior_manifest["transaction_id"],
            "prior_target_commit": prior_manifest["target_commit"],
            "current_target_commit": manifest["target_commit"],
            "prior_packet_sha256": prior_seal["combined_sha256"],
            "prior_candidate_sha256": candidate_sha256,
            "prior_writer_attestation_sha256": sha256_file(writer_attestation_path),
            "prior_review_attestation_sha256": sha256_file(review_attestation_path),
            "prior_review_packet_sha256": review_sha256,
            "prior_review_report_sha256": sha256_file(report),
            "prior_unit_digest_sha256": prior_unit["unit_digest_sha256"],
            "current_unit_digest_sha256": current_unit["unit_digest_sha256"],
            "member_identities_sha256": sha256_bytes(canonical_json(current_identity)),
        }
        if chain != expected_chain:
            raise MemoryErrorBase("reviewed reuse provenance record does not match the verified chain")
        current_staged = transaction / task["staged_output_path"]
        current_data = current_staged.read_bytes()
        _, candidate_body = parse_frontmatter(candidate_data, prior_task["output_repository_path"])
        _, current_body = parse_frontmatter(current_data, task["output_repository_path"])
        if current_body != candidate_body:
            raise MemoryErrorBase("reviewed reuse changed candidate body content")
        expected_raw = candidate_sha256 if staged_sha256 != candidate_sha256 else None
        if attestation.get("raw_output_sha256") != expected_raw:
            raise MemoryErrorBase("reviewed reuse normalization input hash mismatch")
        if attestation.get("normalization") != "source_capsule_frontmatter_v1":
            raise MemoryErrorBase("reviewed reuse normalization contract mismatch")
    except (MemoryErrorBase, OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return [f"{label}: {exc}"]
    return []


def verify_candidate_reuse_chain(
    repo: Path,
    transaction: Path,
    manifest: Mapping[str, Any],
    task: Mapping[str, Any],
    attestation: Mapping[str, Any],
    staged_sha256: str | None,
) -> list[str]:
    """Independently verify an original Codex artifact reused for a fresh review."""
    task_id = task["task_id"]
    label = f"candidate reuse for {task_id}"
    try:
        if task.get("source_unit_id") is None:
            raise MemoryErrorBase("candidate reuse is only supported for source-unit tasks")
        if attestation.get("model_invoked") is not False:
            raise MemoryErrorBase("candidate reuse must explicitly record model_invoked=false")
        if attestation.get("runtime_version") is not None or attestation.get("runtime_executable_sha256") is not None:
            raise MemoryErrorBase("candidate reuse must not claim a runtime invocation")
        chain = attestation.get("candidate_reuse")
        if not isinstance(chain, dict):
            raise MemoryErrorBase("candidate reuse attestation lacks its provenance chain")
        prior_id = chain.get("prior_transaction_id")
        if not isinstance(prior_id, str) or not re.fullmatch(r"[A-Za-z0-9._-]+", prior_id):
            raise MemoryErrorBase("candidate reuse has an invalid prior transaction identity")
        prior_transaction, prior_manifest = load_transaction(repo, prior_id)
        if prior_transaction == transaction:
            raise MemoryErrorBase("candidate reuse must come from another transaction")
        prior_task = next(
            (item for item in prior_manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None,
        )
        if not isinstance(prior_task, dict) or not prior_task.get("required"):
            raise MemoryErrorBase("candidate reuse prior writer task is missing")
        if (
            prior_task.get("task_kind") != task.get("task_kind")
            or prior_task.get("source_unit_id") != task.get("source_unit_id")
            or prior_task.get("output_repository_path") != task.get("output_repository_path")
            or prior_task.get("staged_output_path") != task.get("staged_output_path")
            or prior_task.get("expected_output_name") != task.get("expected_output_name")
        ):
            raise MemoryErrorBase("candidate reuse prior task/output identity mismatch")
        current_unit, current_identity = _reuse_source_identity(manifest, task.get("source_unit_id"))
        prior_unit, prior_identity = _reuse_source_identity(prior_manifest, prior_task.get("source_unit_id"))
        if current_identity != prior_identity or current_unit.get("kind") != prior_unit.get("kind"):
            raise MemoryErrorBase("candidate reuse source member/kind identity mismatch")
        if current_unit.get("unit_digest_sha256") != prior_unit.get("unit_digest_sha256"):
            raise MemoryErrorBase("candidate reuse source-unit prerequisites changed")

        prior_packet = prior_transaction / prior_task["packet_path"]
        prior_seal_path = prior_packet.parent / "packet-seal.json"
        if not prior_seal_path.is_file() or prior_seal_path.is_symlink():
            raise MemoryErrorBase("candidate reuse prior packet seal is missing")
        prior_seal = json.loads(prior_seal_path.read_text(encoding="utf-8"))
        prior_packet_files = packet_file_hashes(prior_packet)
        prior_packet_sha256 = sha256_bytes(canonical_json(prior_packet_files))
        if (
            prior_seal.get("transaction_id") != prior_manifest["transaction_id"]
            or prior_seal.get("task_id") != task_id
            or prior_seal.get("source_unit_id") != prior_task["source_unit_id"]
            or prior_seal.get("files") != prior_packet_files
            or prior_seal.get("combined_sha256") != prior_packet_sha256
        ):
            raise MemoryErrorBase("candidate reuse prior writer packet changed")
        writer_attestation_path = prior_transaction / "attestations" / f"{task_id}.json"
        if not writer_attestation_path.is_file() or writer_attestation_path.is_symlink():
            raise MemoryErrorBase("candidate reuse prior Codex attestation is missing")
        writer_attestation = json.loads(writer_attestation_path.read_text(encoding="utf-8"))
        expected_writer = {
            "transaction_id": prior_manifest["transaction_id"],
            "task_id": task_id,
            "source_unit_id": prior_task["source_unit_id"],
            "isolation": "bubblewrap",
            "workspace_hidden": repo.as_posix(),
            "packet_path": prior_task["packet_path"],
            "packet_sha256": prior_packet_sha256,
            "output_repository_path": prior_task["output_repository_path"],
            "staged_output_path": prior_task["staged_output_path"],
            "runtime_profile": "codex",
        }
        for key, value in expected_writer.items():
            if writer_attestation.get(key) != value:
                raise MemoryErrorBase(f"prior Codex attestation mismatch: {key}")
        for key in ("runner_sha256", "runtime_executable_sha256", "output_sha256"):
            if not isinstance(writer_attestation.get(key), str) or not re.fullmatch(
                r"[0-9a-f]{64}", writer_attestation[key]
            ):
                raise MemoryErrorBase(f"prior Codex attestation lacks {key}")
        if not isinstance(writer_attestation.get("runtime_version"), str) or not writer_attestation["runtime_version"]:
            raise MemoryErrorBase("prior Codex attestation lacks runtime version")
        candidate = prior_transaction / prior_task["staged_output_path"]
        if not candidate.is_file() or candidate.is_symlink():
            raise MemoryErrorBase("candidate reuse prior candidate is missing")
        candidate_data = candidate.read_bytes()
        candidate_sha256 = sha256_bytes(candidate_data)
        if candidate_sha256 != writer_attestation["output_sha256"]:
            raise MemoryErrorBase("candidate reuse prior candidate hash mismatch")

        expected_chain = {
            "prior_transaction_id": prior_manifest["transaction_id"],
            "prior_target_commit": prior_manifest["target_commit"],
            "current_target_commit": manifest["target_commit"],
            "prior_manifest_sha256": sha256_file(prior_transaction / "manifest.json"),
            "prior_packet_seal_sha256": sha256_file(prior_seal_path),
            "prior_packet_sha256": prior_packet_sha256,
            "prior_candidate_sha256": candidate_sha256,
            "prior_writer_attestation_sha256": sha256_file(writer_attestation_path),
            "prior_runner_sha256": writer_attestation["runner_sha256"],
            "prior_runtime_executable_sha256": writer_attestation["runtime_executable_sha256"],
            "prior_unit_digest_sha256": prior_unit["unit_digest_sha256"],
            "current_unit_digest_sha256": current_unit["unit_digest_sha256"],
            "member_identities_sha256": sha256_bytes(canonical_json(current_identity)),
        }
        if chain != expected_chain:
            raise MemoryErrorBase("candidate reuse provenance record does not match the verified chain")
        current_staged = transaction / task["staged_output_path"]
        current_data = current_staged.read_bytes()
        _, candidate_body = parse_frontmatter(candidate_data, prior_task["output_repository_path"])
        _, current_body = parse_frontmatter(current_data, task["output_repository_path"])
        if current_body != candidate_body:
            raise MemoryErrorBase("candidate reuse changed candidate body content")
        expected_raw = candidate_sha256 if staged_sha256 != candidate_sha256 else None
        if attestation.get("raw_output_sha256") != expected_raw:
            raise MemoryErrorBase("candidate reuse normalization input hash mismatch")
        if attestation.get("normalization") != "source_capsule_frontmatter_v1":
            raise MemoryErrorBase("candidate reuse normalization contract mismatch")
    except (MemoryErrorBase, OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return [f"{label}: {exc}"]
    return []


def verify_current_grok_review(
    repo: Path,
    transaction: Path,
    manifest: Mapping[str, Any],
    task: Mapping[str, Any],
    writer_attestation_path: Path,
    writer_attestation: Mapping[str, Any],
    staged_sha256: str | None,
) -> list[str]:
    """Require one exact current Grok PASS for a production Codex candidate."""
    task_id = task["task_id"]
    label = f"current Grok review for {task_id}"
    try:
        packet = transaction / task["packet_path"]
        seal = json.loads((packet.parent / "packet-seal.json").read_text(encoding="utf-8"))
        candidate = transaction / task["staged_output_path"]
        if not candidate.is_file() or candidate.is_symlink() or sha256_file(candidate) != staged_sha256:
            raise MemoryErrorBase("current reviewed candidate is missing or changed")
        review_root = transaction / "reviews" / task_id
        if not review_root.is_dir():
            raise MemoryErrorBase("missing required Grok PASS review")
        if {path.name for path in review_root.iterdir()} != {"packet", "output", "attestation.json"}:
            raise MemoryErrorBase("review directory contains missing or unexpected artifacts")
        output_root = review_root / "output"
        if not output_root.is_dir() or {path.name for path in output_root.iterdir()} != {"report.md"}:
            raise MemoryErrorBase("review output contains missing or unexpected artifacts")
        review_attestation_path = review_root / "attestation.json"
        if not review_attestation_path.is_file() or review_attestation_path.is_symlink():
            raise MemoryErrorBase("review attestation is missing or is not a regular file")
        review = json.loads(review_attestation_path.read_text(encoding="utf-8"))
        report_relative = f"reviews/{task_id}/output/report.md"
        expected_review = {
            "role": "independent_review",
            "transaction_id": manifest["transaction_id"],
            "task_id": task_id,
            "packet_sha256": seal["combined_sha256"],
            "candidate_path": task["staged_output_path"],
            "candidate_sha256": staged_sha256,
            "writer_attestation_sha256": sha256_file(writer_attestation_path),
            "report_path": report_relative,
            "verdict": "PASS",
            "runtime_profile": "grok-review",
            "review_model": GROK_REVIEW_MODEL,
            "review_prompt_sha256": REVIEW_PROMPT_SHA256,
            "review_schema_sha256": REVIEW_SCHEMA_SHA256,
            "review_contract_sha256": REVIEW_CONTRACT_SHA256,
            "reviewer_sha256": sha256_file(Path(__file__).resolve().parent / "review_isolated.py"),
        }
        for key, value in expected_review.items():
            if review.get(key) != value:
                raise MemoryErrorBase(f"Grok PASS attestation mismatch: {key}")
        for key in ("review_packet_sha256", "report_sha256", "runtime_executable_sha256"):
            if not isinstance(review.get(key), str) or not re.fullmatch(r"[0-9a-f]{64}", review[key]):
                raise MemoryErrorBase(f"Grok attestation lacks {key}")
        if not isinstance(review.get("runtime_version"), str) or not review["runtime_version"]:
            raise MemoryErrorBase("Grok attestation lacks runtime version")
        report = transaction / report_relative
        if not report.is_file() or report.is_symlink() or sha256_file(report) != review["report_sha256"]:
            raise MemoryErrorBase("Grok report hash mismatch")
        review_packet = review_root / "packet"
        review_files = packet_file_hashes(review_packet)
        review_sha256 = sha256_bytes(canonical_json(review_files))
        if (
            review_files != review.get("review_packet_files")
            or review_sha256 != review["review_packet_sha256"]
            or review_files.get("candidate.md") != staged_sha256
        ):
            raise MemoryErrorBase("Grok review packet or candidate hash mismatch")
        expected_review_files = dict(seal["files"])
        expected_review_files.update({
            "candidate.md": staged_sha256,
            "review-prompt.md": review_files.get("review-prompt.md"),
        })
        if review_files != expected_review_files:
            raise MemoryErrorBase("Grok review packet is not the exact writer packet plus candidate")
        if writer_attestation.get("packet_sha256") != seal["combined_sha256"]:
            raise MemoryErrorBase("writer packet identity changed after review")
    except (MemoryErrorBase, OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return [f"{label}: {exc}"]
    return []


def verify_isolation_attestations(
    repo: Path, transaction: Path, manifest: Mapping[str, Any], staged_hashes: Mapping[str, str]
) -> list[str]:
    errors: list[str] = []
    runner = Path(__file__).resolve().parent / "run_isolated.py"
    runner_hash = sha256_file(runner)
    required = [task for task in manifest.get("writer_tasks", []) if task.get("required")]
    expected_attestations = {f"{task['task_id']}.json" for task in required}
    attestations_root = transaction / "attestations"
    actual_attestations = {
        path.name for path in attestations_root.iterdir() if path.is_file()
    } if attestations_root.is_dir() else set()
    expected_current_reviews: set[str] = set()
    for extra in sorted(actual_attestations - expected_attestations):
        errors.append(f"unexpected isolation attestation: {extra}")
    for task in required:
        task_id = task["task_id"]
        attestation_path = attestations_root / f"{task_id}.json"
        if not attestation_path.is_file():
            errors.append(
                f"missing bubblewrap isolation attestation, candidate-reuse attestation, "
                f"or reviewed-reuse attestation for {task_id}"
            )
            continue
        try:
            attestation = json.loads(attestation_path.read_text(encoding="utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            errors.append(f"invalid isolation attestation for {task_id}: {exc}")
            continue
        task_root = transaction / "tasks" / task_id
        packet = task_root / "packet"
        seal_path = task_root / "packet-seal.json"
        if not seal_path.is_file():
            errors.append(f"missing sealed task packet for {task_id}")
            continue
        seal = json.loads(seal_path.read_text(encoding="utf-8"))
        packet_files = {
            path.relative_to(packet).as_posix(): sha256_file(path)
            for path in sorted(item for item in packet.rglob("*") if item.is_file())
        }
        packet_hash = sha256_bytes(canonical_json(packet_files))
        if seal.get("files") != packet_files or seal.get("combined_sha256") != packet_hash:
            errors.append(f"sealed task packet changed for {task_id}")
        runtime_profile = attestation.get("runtime_profile")
        reviewed_reuse = runtime_profile == "codex-reviewed-reuse"
        candidate_reuse = runtime_profile == "codex-candidate-reuse"
        expected = {
            "transaction_id": manifest["transaction_id"],
            "task_id": task_id,
            "source_unit_id": task.get("source_unit_id"),
            "isolation": (
                "deterministic-reviewed-reuse" if reviewed_reuse else
                "deterministic-candidate-reuse" if candidate_reuse else
                "bubblewrap"
            ),
            "workspace_hidden": repo.as_posix(),
            "packet_path": task["packet_path"],
            "packet_sha256": packet_hash,
            "output_repository_path": task["output_repository_path"],
            "staged_output_path": task["staged_output_path"],
            "output_sha256": staged_hashes.get(task["output_repository_path"]),
            "runner_sha256": runner_hash,
        }
        for key, value in expected.items():
            if attestation.get(key) != value:
                errors.append(f"isolation attestation mismatch for {task_id}: {key}")
        if reviewed_reuse:
            errors.extend(verify_reviewed_reuse_chain(
                repo, transaction, manifest, task, attestation,
                staged_hashes.get(task["output_repository_path"]),
            ))
            continue
        if candidate_reuse:
            errors.extend(verify_candidate_reuse_chain(
                repo, transaction, manifest, task, attestation,
                staged_hashes.get(task["output_repository_path"]),
            ))
            expected_current_reviews.add(task_id)
            errors.extend(verify_current_grok_review(
                repo, transaction, manifest, task, attestation_path, attestation,
                staged_hashes.get(task["output_repository_path"]),
            ))
            continue
        allowed_profiles = {"codex"}
        if os.environ.get("MEMORY_TEST_ALLOW_PROFILE") == "1":
            allowed_profiles.add("_test")
        if runtime_profile not in allowed_profiles:
            errors.append(f"isolation attestation has invalid runtime profile for {task_id}")
        if not isinstance(attestation.get("runtime_version"), str) or not attestation.get("runtime_version"):
            errors.append(f"isolation attestation lacks runtime version for {task_id}")
        if not isinstance(attestation.get("runtime_executable_sha256"), str) or not re.fullmatch(
            r"[0-9a-f]{64}", attestation.get("runtime_executable_sha256", "")
        ):
            errors.append(f"isolation attestation lacks runtime executable identity for {task_id}")
        if runtime_profile == "codex":
            expected_current_reviews.add(task_id)
            errors.extend(verify_current_grok_review(
                repo, transaction, manifest, task, attestation_path, attestation,
                staged_hashes.get(task["output_repository_path"]),
            ))
    reviews_root = transaction / "reviews"
    actual_reviews = {path.name for path in reviews_root.iterdir()} if reviews_root.is_dir() else set()
    for unexpected in sorted(actual_reviews - expected_current_reviews):
        errors.append(f"unexpected current review artifact: {unexpected}")
    return errors


def _recheck_transaction(
    repo: Path, manifest: Mapping[str, Any], config_path: str = DEFAULT_CONFIG,
) -> tuple[dict[str, ResolvedUnit], dict[str, Any], bytes | None, dict[str, str]]:
    head = resolve_commit(repo, "HEAD")
    if head != manifest["target_commit"]:
        raise MemoryErrorBase(f"HEAD moved since prepare: expected {manifest['target_commit']}, found {head}")
    config, config_bytes = load_config(repo, config_path)
    policy = policy_digests(repo, config_bytes)
    if policy != manifest["policy"]:
        raise MemoryErrorBase("policy/config/schema/prompt/tool digest drift since prepare")
    state, state_bytes = load_state(repo)
    if state_digest(state_bytes) != manifest["prior_state_sha256"]:
        raise MemoryErrorBase("successful state changed since prepare")
    entries = enumerate_tree(repo, head, config["discovery"]["candidate_roots"])
    resolved = resolve_units(config, entries)
    for item in manifest["units"]:
        if item["action"] == "retired_from_corpus":
            if item["id"] in resolved:
                raise MemoryErrorBase(f"retired unit returned to config: {item['id']}")
            continue
        current = resolved.get(item["id"])
        if current is None or current.digest != item["unit_digest_sha256"]:
            raise MemoryErrorBase(f"resolved unit drift since prepare: {item['id']}")
        stripped = [{key: member.get(key) for key in (
            "path", "role", "read_mode", "ownership", "mode", "object_type", "blob_oid", "blob_size",
            "anchors", "context_bytes", "max_bytes"
        ) if key in member} for member in item["members"]]
        if stripped != current.members:
            raise MemoryErrorBase(f"transaction member identity drift: {item['id']}")
    return resolved, state, state_bytes, policy


def build_id_registry(pages: Mapping[str, bytes]) -> dict[str, dict[str, str]]:
    registry: dict[str, dict[str, str]] = {"pages": {}, "statements": {}, "conflicts": {}}
    for path, data in sorted(pages.items()):
        front, body = parse_frontmatter(data, path)
        page_id = front.get("id")
        if isinstance(page_id, str):
            registry["pages"][page_id] = path
        for match in re.finditer(r"^###\s+([a-z0-9][a-z0-9-]*)\s+—", body, re.MULTILINE):
            registry["statements"][match.group(1)] = path
        if front.get("type") == "conflict_register":
            for match in re.finditer(r"^##\s+([a-z0-9][a-z0-9-]*)\s+—", body, re.MULTILINE):
                registry["conflicts"][match.group(1)] = path
    return registry


def _new_state_after_transaction(
    old: Mapping[str, Any], manifest: Mapping[str, Any], resolved: Mapping[str, ResolvedUnit],
    page_hashes: Mapping[str, str], id_registry: Mapping[str, Mapping[str, str]],
) -> dict[str, Any]:
    state = json.loads(json.dumps(old))
    state["state_version"] = STATE_VERSION
    state["generation"] = manifest["transaction_id"]
    state["policy"] = manifest["policy"]
    state.setdefault("units", {})
    state.setdefault("retired_units", {})
    state.setdefault("pages", {})
    state.setdefault("derived_pages", {})
    target = manifest["target_commit"]
    for item in manifest["units"]:
        unit_id = item["id"]
        action = item["action"]
        if action == "retired_from_corpus":
            old_unit = state["units"].pop(unit_id)
            prior_result = old_unit.get("result")
            old_unit["result"] = "retired_from_corpus"
            old_unit["lifecycle"] = "deleted" if prior_result == "source_deleted" else "retired"
            old_unit["prior_result"] = prior_result
            old_unit["retired_at_commit"] = target
            capsule = old_unit.get("capsule")
            if capsule in page_hashes:
                old_unit["page_sha256"] = page_hashes[capsule]
                state["pages"][capsule] = {
                    "sha256": page_hashes[capsule],
                    "generation": manifest["transaction_id"],
                    "unit_id": unit_id,
                    "processed_commit": target,
                    "generated_commit": target,
                    "lifecycle": old_unit["lifecycle"],
                    "unit_digest_sha256": old_unit.get("unit_digest_sha256"),
                    "policy_sha256": manifest["policy"]["combined_sha256"],
                    "retired_from_corpus": True,
                }
            state["retired_units"][unit_id] = old_unit
            continue
        current = resolved[unit_id]
        capsule = current.unit["capsule"]
        old_record = state["units"].get(unit_id, {})
        if action == "rebase_unchanged":
            if not old_record:
                raise MemoryErrorBase(f"cannot acknowledge unchanged unit without prior state: {unit_id}")
            old_record["processed_commit"] = target
            old_record["policy_sha256"] = manifest["policy"]["combined_sha256"]
            old_record["last_action"] = action
            old_record["last_error"] = None
            state["units"][unit_id] = old_record
            page_record = state.get("pages", {}).get(capsule)
            if page_record is not None:
                page_record["processed_commit"] = target
                # generated_commit intentionally remains the snapshot that
                # actually produced the immutable served page bytes.
            continue
        record = {
            "id": unit_id,
            "kind": current.unit["kind"],
            "lifecycle": current.unit["lifecycle"],
            "capsule": capsule,
            "processed_commit": target,
            "unit_digest_sha256": current.digest,
            "policy_sha256": manifest["policy"]["combined_sha256"],
            "members": current.members,
            "missing_required": current.missing_required,
            "result": "source_deleted" if action == "source_deleted" else "success",
            "last_action": action,
            "last_error": None,
        }
        if action == "source_deleted":
            record["last_seen_members"] = old_record.get("last_seen_members") or old_record.get("members", [])
        if item.get("requires_staged_page"):
            digest = page_hashes[capsule]
            record["page_sha256"] = digest
            state["pages"][capsule] = {
                "sha256": digest,
                "generation": manifest["transaction_id"],
                "unit_id": unit_id,
                "processed_commit": target,
                "generated_commit": target,
                "lifecycle": "deleted" if action == "source_deleted" else current.unit["lifecycle"],
                "unit_digest_sha256": current.digest,
                "policy_sha256": manifest["policy"]["combined_sha256"],
            }
        else:
            record["page_sha256"] = old_record.get("page_sha256")
            if capsule in state["pages"]:
                state["pages"][capsule]["processed_commit"] = target
        state["units"][unit_id] = record
        state["retired_units"].pop(unit_id, None)
    for task in manifest.get("writer_tasks", []):
        if task.get("action") == "derived_rebase_unchanged":
            record = state.get("derived_pages", {}).get(task["task_id"])
            if record is None:
                raise MemoryErrorBase(f"cannot acknowledge unchanged derived page without state: {task['task_id']}")
            record["processed_commit"] = target
            record["last_action"] = task["action"]
            page_record = state.get("pages", {}).get(task["output_repository_path"])
            if page_record is not None:
                page_record["processed_commit"] = target
            continue
        if not task.get("required") or task.get("source_unit_id") is not None:
            continue
        page = task["output_repository_path"]
        digest = page_hashes[page]
        state["pages"][page] = {
            "sha256": digest,
            "generation": manifest["transaction_id"],
            "task_id": task["task_id"],
            "processed_commit": target,
            "generated_commit": target,
            "lifecycle": task["semantic_contract"]["page"]["desired_lifecycle"],
            "unit_digest_sha256": None,
            "policy_sha256": manifest["policy"]["combined_sha256"],
        }
        dependency_inputs = json.loads(json.dumps(
            task["semantic_contract"].get("static_memory_inputs", [])
            + task["semantic_contract"].get("dynamic_memory_inputs", [])
        ))
        for dependency in dependency_inputs:
            dependency_page = dependency.get("page")
            if dependency_page in page_hashes:
                dependency["page_sha256"] = page_hashes[dependency_page]
        state["derived_pages"][task["task_id"]] = {
            "id": task["task_id"],
            "task_kind": task["task_kind"],
            "page": page,
            "processed_commit": target,
            "policy_sha256": manifest["policy"]["combined_sha256"],
            "page_sha256": digest,
            "input_units": sorted(task["semantic_contract"].get("input_units", {})),
            "input_pages": list(task["semantic_contract"].get("input_pages", [])),
            "migration_ids": [
                item["legacy_id"] for item in task["semantic_contract"].get("migration_requirements", [])
            ],
            "migration_requirements": task["semantic_contract"].get("migration_requirements", []),
            "dependency_inputs": dependency_inputs,
            "generated_commit": target,
            "lifecycle": task["semantic_contract"]["page"]["desired_lifecycle"],
            "result": "success",
        }
    # A full baseline is stricter than content freshness: each active unit must
    # have been acknowledged at this exact target and current policy.
    fully = True
    for unit_id, current in resolved.items():
        record = state["units"].get(unit_id)
        if not record or record.get("processed_commit") != target or record.get("unit_digest_sha256") != current.digest:
            fully = False
            break
        if record.get("policy_sha256") != manifest["policy"]["combined_sha256"]:
            fully = False
            break
    if fully:
        if manifest.get("derived_task_plan", {}).get("blocked"):
            fully = False
        else:
            for derived_id in manifest.get("configured_derived_ids", []):
                record = state.get("derived_pages", {}).get(derived_id)
                if (
                    record is None or record.get("processed_commit") != target
                    or record.get("policy_sha256") != manifest["policy"]["combined_sha256"]
                ):
                    fully = False
                    break
    state["last_fully_processed_commit"] = target if fully else old.get("last_fully_processed_commit")
    state["last_result"] = {
        "status": "success",
        "transaction_id": manifest["transaction_id"],
        "target_commit": target,
        "processed_units": [item["id"] for item in manifest["units"]],
        "error": None,
    }
    state["id_registry"] = id_registry
    return state


def finalize_update(repo: Path, transaction_arg: str, config_path: str = DEFAULT_CONFIG) -> dict[str, Any]:
    lock = acquire_lock(repo)
    journal_path, _ = publication_paths(repo)
    journal: dict[str, Any] | None = None
    transaction: Path | None = None
    manifest: dict[str, Any] | None = None
    backup_root: Path | None = None
    try:
        transaction, manifest = load_transaction(repo, transaction_arg)
        resolved, old_state, old_state_bytes, _ = _recheck_transaction(repo, manifest, config_path)
        for task in manifest.get("writer_tasks", []):
            if not task.get("required") or not task.get("generated_regions"):
                continue
            live_base = repo / task["output_repository_path"]
            if not live_base.is_file() or sha256_file(live_base) != task.get("base_page_sha256"):
                raise MemoryErrorBase(
                    f"live mixed target changed since prepare: {task['output_repository_path']}"
                )
        expected_stage_paths = {
            item["output_repository_path"] for item in manifest.get("writer_tasks", []) if item.get("required")
        }
        before_lint = {
            path: sha256_file(transaction / "staged" / path)
            for path in expected_stage_paths if (transaction / "staged" / path).is_file()
        }
        errors, warnings = lint_staged(repo, str(transaction))
        if errors:
            raise LintFailure(errors, warnings)
        staged_data = {
            path: (transaction / "staged" / path).read_bytes()
            for path in expected_stage_paths
        }
        after_lint = {path: sha256_bytes(data) for path, data in staged_data.items()}
        if before_lint != after_lint:
            raise MemoryErrorBase("staged page set changed while finalize was linting it")
        actual_stage_paths = {
            item.relative_to(transaction / "staged").as_posix()
            for item in (transaction / "staged").rglob("*") if item.is_file()
        }
        if actual_stage_paths != expected_stage_paths:
            raise MemoryErrorBase("staged page set changed after completeness validation")
        attestation_errors = verify_isolation_attestations(repo, transaction, manifest, after_lint)
        if attestation_errors:
            raise LintFailure(attestation_errors)
        page_hashes: dict[str, str] = {}
        targets: list[dict[str, Any]] = []
        backup_root = repo / "memory/backups" / manifest["transaction_id"]
        backup_root.mkdir(parents=True, exist_ok=False)
        for item in manifest.get("writer_tasks", []):
            if not item.get("required"):
                continue
            relative = item["output_repository_path"]
            page_hashes[relative] = after_lint[relative]
            target = repo / relative
            existed = target.is_file()
            record = {"path": relative, "existed": existed, "backup": relative}
            if existed:
                backup = backup_root / relative
                atomic_write(backup, target.read_bytes())
            targets.append(record)
        state_path = repo / DEFAULT_STATE
        state_existed = state_path.is_file()
        state_backup = DEFAULT_STATE
        if state_existed:
            atomic_write(backup_root / state_backup, state_path.read_bytes())
        old_served_pages, served_errors = served_page_bytes(repo, old_state)
        if served_errors:
            raise LintFailure(served_errors)
        effective_pages = dict(old_served_pages)
        effective_pages.update(staged_data)
        effective_id_errors = lint_effective_ids(effective_pages)
        if effective_id_errors:
            raise LintFailure(effective_id_errors)
        id_registry = build_id_registry(effective_pages)
        new_state = _new_state_after_transaction(old_state, manifest, resolved, page_hashes, id_registry)
        new_state_bytes = canonical_json(new_state)
        journal = {
            "journal_version": 1,
            "transaction_id": manifest["transaction_id"],
            "phase": "prepared",
            "backup_root": backup_root.relative_to(repo).as_posix(),
            "targets": targets,
            "state": {"existed": state_existed, "backup": state_backup},
            "new_state_sha256": sha256_bytes(new_state_bytes),
        }
        write_json(journal_path, journal)
        for record in targets:
            atomic_write(repo / record["path"], staged_data[record["path"]])
        journal["phase"] = "pages_published"
        write_json(journal_path, journal)
        if os.environ.get("MEMORY_TEST_FAIL_AFTER_PAGES") == "1":
            raise OSError("injected publication failure after page replacement")
        atomic_write(state_path, new_state_bytes, 0o644)
        journal["phase"] = "state_committed"
        write_json(journal_path, journal)
        _cleanup_publication(repo, journal)
        return {
            "transaction_id": manifest["transaction_id"],
            "target_commit": manifest["target_commit"],
            "processed_units": [item["id"] for item in manifest["units"]],
            "published_pages": sorted(page_hashes),
            "last_fully_processed_commit": new_state["last_fully_processed_commit"],
            "warnings": warnings,
        }
    except Exception as exc:
        if journal is not None and journal_path.exists():
            _restore_publication(repo, journal)
            _cleanup_publication(repo, journal)
        elif backup_root is not None and backup_root.exists():
            shutil.rmtree(backup_root)
        if transaction is not None:
            try:
                write_json(transaction / "failure.json", {
                    "failed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                })
            except OSError:
                pass
        raise
    finally:
        release_lock(lock)


def query_memory(repo: Path, question: str, limit: int = 8) -> dict[str, Any]:
    if not question.strip():
        raise MemoryErrorBase("query cannot be empty")
    state, _ = load_state(repo)
    pages, errors = served_page_bytes(repo, state)
    if errors:
        raise LintFailure(errors)
    report = status_report(repo)
    warnings: list[str] = []
    if report["pending_units"]:
        warnings.append(
            "Memory is stale or partial; pending units: " + ", ".join(report["pending_units"])
        )
    if report["pending_derived_pages"]:
        warnings.append(
            "Derived memory pages are stale or partial: " + ", ".join(report["pending_derived_pages"])
        )
    if report["dirty_tracked_members"]:
        warnings.append(
            "Committed memory does not include dirty tracked members: " + ", ".join(report["dirty_tracked_members"])
        )
    terms = [term for term in re.findall(r"[a-z0-9][a-z0-9_-]+", question.lower()) if len(term) > 1]
    normalized_question = question.strip().lower()
    results: list[dict[str, Any]] = []
    for path, data in pages.items():
        front, body = parse_frontmatter(data, path)
        title = str(front.get("title", ""))
        page_id = str(front.get("id", ""))
        metadata = " ".join((path, page_id, title, " ".join(map(str, front.get("sources", []))))).lower()
        full = data.decode("utf-8", "replace").lower()
        score = 0
        if normalized_question in (path.lower(), page_id.lower(), title.lower()):
            score += 2000
        elif normalized_question and normalized_question in metadata:
            score += 600
        for term in terms:
            if term == page_id.lower() or term == path.lower():
                score += 300
            score += metadata.count(term) * 50
            score += min(full.count(term), 20) * 3
        matched = score > 0
        if front.get("lifecycle") == "current":
            score += 5
        if front.get("memory_review") == "reviewed":
            score += 3
        if not matched:
            continue
        lines = body.splitlines()
        matching = [index for index, line in enumerate(lines) if any(term in line.lower() for term in terms)]
        start = max(0, (matching[0] if matching else 0) - 1)
        excerpt = "\n".join(lines[start : start + 8]).strip()
        results.append({
            "path": path, "id": page_id, "title": title, "type": front.get("type"),
            "lifecycle": front.get("lifecycle"), "memory_review": front.get("memory_review"),
            "score": score, "excerpt": excerpt,
        })
    results.sort(key=lambda item: (-item["score"], item["path"]))
    return {"question": question, "warnings": warnings, "results": results[:limit]}


def print_status_text(report: Mapping[str, Any]) -> None:
    print(f"target commit: {report['target_commit']}")
    print(f"last fully processed: {report['last_fully_processed_commit'] or 'uninitialized'}")
    print(f"fully fresh: {'yes' if report['fully_fresh'] else 'no'}")
    print(f"policy drift: {'yes' if report['policy_drift'] else 'no'}")
    counts: dict[str, int] = {}
    for detail in report["units"].values():
        counts[detail["action"]] = counts.get(detail["action"], 0) + 1
    print("actions: " + ", ".join(f"{key}={counts[key]}" for key in sorted(counts)))
    for unit_id in report["pending_units"]:
        detail = report["units"][unit_id]
        suffix = f" missing={detail['missing_required']}" if detail["missing_required"] else ""
        print(f"  {unit_id}: {detail['action']}{suffix}")
    if report["pending_derived_pages"]:
        print("pending derived pages:")
        for page_id in report["pending_derived_pages"]:
            print(f"  {page_id}: {report['derived_pages'][page_id]['action']}")
    if report["dirty_tracked_members"]:
        print("dirty tracked members:")
        for path in report["dirty_tracked_members"]:
            print(f"  {path}")
    if report["transactions"]:
        print("transactions: " + ", ".join(report["transactions"]))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, help="repository root (normally auto-detected)")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("init", help="create missing runtime directories and empty state")
    status = sub.add_parser("status", help="compare configured units with a committed snapshot")
    status.add_argument("--commit", default="HEAD")
    status.add_argument("--json", action="store_true")
    update = sub.add_parser("update", help="prepare immutable inputs or finalize staged pages")
    mode = update.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--finalize", metavar="TRANSACTION")
    update.add_argument("--commit", default="HEAD")
    update.add_argument("--units", nargs="*", default=[])
    update.add_argument("--paths", nargs="*", default=[])
    update.add_argument(
        "--allow-large-batch", action="store_true",
        help="explicitly override configured prepare batch unit/task/byte bounds",
    )
    lint = sub.add_parser("lint", help="validate served pages or a staged transaction")
    lint_mode = lint.add_mutually_exclusive_group()
    lint_mode.add_argument("--served", action="store_true")
    lint_mode.add_argument("--staged", metavar="TRANSACTION")
    lint.add_argument(
        "--record", action="store_true",
        help="record one exact task-scoped staged lint failure after a current Grok PASS",
    )
    lint.add_argument("--task", help="writer task whose exact staged lint failures should be recorded")
    lint.add_argument("--json", action="store_true")
    query = sub.add_parser("query", help="retrieve ranked local memory excerpts")
    query.add_argument("question", nargs="+")
    query.add_argument("--limit", type=int, default=8)
    query.add_argument("--json", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "lint" and args.record and not args.staged:
        parser.error("--record requires --staged TRANSACTION")
    if args.command == "lint" and args.task and not args.record:
        parser.error("--task requires --record")
    try:
        repo = args.repo.resolve() if args.repo else find_repo()
        recovered = recover_publication(repo)
        if recovered:
            print(f"recovery: {recovered}", file=sys.stderr)
        if args.command == "init":
            created = init_memory(repo)
            print("created: " + (", ".join(created) if created else "nothing (already initialized)"))
        elif args.command == "status":
            report = status_report(repo, args.commit)
            if args.json:
                print(json.dumps(report, indent=2, sort_keys=True))
            else:
                print_status_text(report)
        elif args.command == "update" and args.prepare:
            transaction = prepare_update(
                repo, args.commit, args.units, args.paths, allow_large_batch=args.allow_large_batch,
            )
            print(transaction.relative_to(repo).as_posix())
        elif args.command == "update" and args.finalize:
            result = finalize_update(repo, args.finalize)
            print(json.dumps(result, indent=2, sort_keys=True))
        elif args.command == "lint":
            if args.staged:
                if args.record:
                    record = record_staged_lint_failure(repo, args.staged, args.task)
                    errors = list(record["errors"])
                    warnings = list(record["warnings"])
                    print(f"recorded: {record['report_path']}")
                else:
                    errors, warnings = lint_staged(repo, args.staged)
                mode_name = "staged"
            else:
                errors, warnings = lint_served(repo)
                mode_name = "served"
            result = {"mode": mode_name, "errors": errors, "warnings": warnings}
            if args.json:
                print(json.dumps(result, indent=2, sort_keys=True))
            else:
                for warning in warnings:
                    print(f"warning: {warning}")
                for error in errors:
                    print(f"error: {error}")
                print(f"{mode_name} lint: {len(errors)} error(s), {len(warnings)} warning(s)")
            if errors:
                return 1
        elif args.command == "query":
            result = query_memory(repo, " ".join(args.question), args.limit)
            if args.json:
                print(json.dumps(result, indent=2, sort_keys=True))
            else:
                for warning in result["warnings"]:
                    print(f"WARNING: {warning}")
                for item in result["results"]:
                    print(f"\n{item['path']} — {item['title']} (score {item['score']})")
                    if item["excerpt"]:
                        print(item["excerpt"])
        return 0
    except LintFailure as exc:
        for warning in exc.warnings:
            print(f"warning: {warning}", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        return 1
    except (MemoryErrorBase, OSError, ValueError) as exc:
        print(f"memory: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
