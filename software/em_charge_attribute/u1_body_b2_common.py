#!/usr/bin/env python3
"""Shared strict I/O and canonicalization for the U1 Phase-B2 staged build."""

from __future__ import annotations

import hashlib
import json
import os
import re
import zlib
from pathlib import Path
from typing import Any, Iterable

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
REPORTS = HERE / "reports" / "u1_body_dynamics_artifacts"
STAGE0 = REPORTS / "stage_b2_0_intake_radiative_contract"


class UniqueKeyLoader(yaml.SafeLoader):
    pass


def _construct_mapping(loader: UniqueKeyLoader, node: yaml.MappingNode, deep: bool = False) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f"duplicate YAML key: {key!r}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


def load_yaml(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.load(stream, Loader=UniqueKeyLoader)


def load_yaml_authenticated(path: Path, expected_sha256: str, consumer: str) -> tuple[Any, dict[str, Any]]:
    """Authenticate bytes on a held descriptor, then parse those exact bytes.

    The pre/post fstat evidence makes path replacement distinguishable from an
    immutable first-use blob.  Callers retain the returned record in their own
    generated artifact so authentication is attributable per consumer.
    """
    resolved = path.resolve()
    fd = os.open(resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    try:
        before = os.fstat(fd)
        chunks: list[bytes] = []
        while True:
            chunk = os.read(fd, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    raw = b"".join(chunks)
    observed = sha256_bytes(raw)
    require(observed == expected_sha256, "B2_A1_PROTECTED_FIRST_USE", f"{consumer}:{resolved}")
    require(
        (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns),
        "B2_A1_PROTECTED_REWRITE",
        f"{consumer}:{resolved}:held descriptor changed during consumption",
    )
    value = yaml.load(raw.decode("utf-8"), Loader=UniqueKeyLoader)
    evidence = {
        "consumer": consumer,
        "path": rel_repo(resolved),
        "expected_sha256": expected_sha256,
        "consumed_sha256": observed,
        "held_descriptor": True,
        "device": before.st_dev,
        "inode": before.st_ino,
        "size": before.st_size,
        "mtime_ns_before": before.st_mtime_ns,
        "mtime_ns_after": after.st_mtime_ns,
        "descriptor_stable_during_parse": True,
    }
    return value, evidence


def read_bytes_authenticated(path: Path, expected_sha256: str, consumer: str) -> tuple[bytes, dict[str, Any]]:
    """Read a non-YAML protected input from a stable held descriptor."""
    resolved = path.resolve()
    fd = os.open(resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    try:
        before = os.fstat(fd); chunks: list[bytes] = []
        while True:
            chunk = os.read(fd, 1 << 20)
            if not chunk: break
            chunks.append(chunk)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    raw = b"".join(chunks); observed = sha256_bytes(raw)
    require(observed == expected_sha256, "B2_A1_PROTECTED_FIRST_USE", f"{consumer}:{resolved}")
    require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), "B2_A1_PROTECTED_REWRITE", f"{consumer}:{resolved}")
    return raw, {"consumer": consumer, "path": rel_repo(resolved), "expected_sha256": expected_sha256, "consumed_sha256": observed, "held_descriptor": True, "device": before.st_dev, "inode": before.st_ino, "descriptor_stable_during_parse": True}


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=220), encoding="utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str).encode()


def digest(value: Any) -> str:
    return sha256_bytes(canonical_json(value))


def repository_path(raw: str, base: Path = HERE) -> Path:
    if raw.startswith("//"):
        return (ROOT / raw[2:]).resolve()
    return (base / raw).resolve()


def certificate_path(raw: str) -> Path:
    return repository_path(raw, HERE)


def rel_repo(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def _loose_git_object(object_id: str) -> tuple[str, bytes]:
    require(re.fullmatch(r"[0-9a-f]{40}", object_id) is not None, "B2_A1_IMMUTABLE_COMMIT", f"invalid object id:{object_id}")
    raw = zlib.decompress((ROOT / ".git" / "objects" / object_id[:2] / object_id[2:]).read_bytes())
    header, payload = raw.split(b"\0", 1)
    kind, size = header.decode("ascii").split(" ", 1)
    require(int(size) == len(payload), "B2_A1_IMMUTABLE_COMMIT", f"object size:{object_id}")
    return kind, payload


def _git_tree_entry(tree_id: str, name: str) -> tuple[str, str]:
    kind, tree = _loose_git_object(tree_id)
    require(kind == "tree", "B2_A1_IMMUTABLE_COMMIT", f"expected tree:{tree_id}")
    cursor = 0
    while cursor < len(tree):
        separator = tree.index(b" ", cursor)
        terminator = tree.index(b"\0", separator)
        mode = tree[cursor:separator].decode("ascii")
        entry_name = tree[separator + 1:terminator].decode("utf-8", errors="surrogateescape")
        entry_id = tree[terminator + 1:terminator + 21].hex()
        if entry_name == name:
            return mode, entry_id
        cursor = terminator + 21
    raise FileNotFoundError(name)


def git_blob(commit: str, relative: str) -> bytes:
    kind, commit_payload = _loose_git_object(commit)
    require(kind == "commit", "B2_A1_IMMUTABLE_COMMIT", f"expected commit:{commit}")
    tree_line = next((line for line in commit_payload.splitlines() if line.startswith(b"tree ")), None)
    require(tree_line is not None, "B2_A1_IMMUTABLE_COMMIT", f"commit tree:{commit}")
    object_id = tree_line.split(b" ", 1)[1].decode("ascii")
    parts = [part for part in relative.split("/") if part]
    require(bool(parts) and all(part not in {".", ".."} for part in parts), "B2_A1_IMMUTABLE_COMMIT", f"invalid path:{relative}")
    for part in parts:
        _, object_id = _git_tree_entry(object_id, part)
    blob_kind, payload = _loose_git_object(object_id)
    require(blob_kind == "blob", "B2_A1_IMMUTABLE_COMMIT", f"expected blob:{commit}:{relative}")
    return payload


def git_object_exists(commit: str, relative: str) -> bool:
    try:
        git_blob(commit, relative)
        return True
    except (FileNotFoundError, ValueError):
        return False


def require(condition: bool, tooth: str, detail: str) -> None:
    if not condition:
        raise AssertionError(f"ASSERT_FAIL:{tooth}:{detail}")


def deep_get(value: Any, path: Iterable[Any]) -> Any:
    current = value
    for part in path:
        current = current[part]
    return current


def flatten_paths(value: Any, prefix: tuple[Any, ...] = ()) -> list[tuple[Any, ...]]:
    if isinstance(value, dict):
        return [p for key, child in value.items() for p in flatten_paths(child, prefix + (key,))]
    if isinstance(value, list):
        return [p for i, child in enumerate(value) for p in flatten_paths(child, prefix + (i,))]
    return [prefix]


FORBIDDEN_ANCESTRY = {
    "Maxwell",
    "Larmor",
    "point_current",
    "point_charge",
    "classical_particle_F_equals_dPdt",
    "j_equals_sV",
}


def assert_native_ancestry(nodes: Iterable[str], tooth: str = "B2_NATIVE_ANCESTRY") -> None:
    found = sorted(set(nodes) & FORBIDDEN_ANCESTRY)
    require(not found, tooth, f"forbidden ancestry nodes {found}")


def parse_manifest_sets(certificate: dict[str, Any]) -> dict[str, dict[str, str]]:
    return {
        "P": dict(certificate["protected_artifact_sha256"]),
        "S": dict(certificate["substrate_reference_sha256"]),
        "Dep_B1": dict(certificate["runtime_dependency_sha256"]),
    }


def manifest_absolute(raw: str) -> Path:
    return certificate_path(raw)


COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")


def validate_commit(commit: str) -> None:
    require(bool(COMMIT_RE.fullmatch(commit)), "B2_A1_STARTUP_COMMIT", "startup_contract_commit must be a full lowercase 40-hex hash")


def normalized_monomials(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    normalized = []
    for row in rows:
        normalized.append({"coefficient": str(row["coefficient"]), "powers": dict(sorted(row.get("powers", {}).items()))})
    return sorted(normalized, key=lambda row: (tuple(row["powers"].items()), row["coefficient"]))
