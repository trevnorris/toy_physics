#!/usr/bin/env python3
"""Inspect, archive, restore, and safely evict research-memory transactions.

The live transaction tree is intentionally gitignored and can grow quickly.
This utility keeps cleanup separate from publication: every mutating command is
a dry run unless ``--apply`` is supplied, and eviction is allowed only after a
lossless archive has been written and independently verified.

Archives retain the complete transaction, including writer/reviewer packets,
attestations, reports, and staged candidates required by revision, candidate
reuse, or reviewed reuse.  Transactions referenced by another transaction, matching current
served state, or named by the current state generation are never evicted.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from typing import Any, Mapping, Sequence


RECEIPT_VERSION = 1
LIFECYCLE_VERSION = 1
TRANSACTION_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
ARCHIVE_PREFIX = "_archive-"


class TransactionGCError(RuntimeError):
    """A transaction cannot be classified or changed safely."""


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, data: bytes, mode: int = 0o600) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary_path, mode)
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def write_json(path: Path, value: Any) -> None:
    atomic_write(path, canonical_json(value))


def find_repo(start: Path | None = None) -> Path:
    proc = subprocess.run(
        ["git", "-C", str((start or Path.cwd()).resolve()), "rev-parse", "--show-toplevel"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if proc.returncode:
        raise TransactionGCError("not inside a Git repository")
    return Path(proc.stdout.decode("utf-8", "strict").strip()).resolve()


def transaction_root(repo: Path) -> Path:
    return repo / "memory/transactions"


def archive_root(repo: Path) -> Path:
    # Keep archive artifacts as files directly under the already-ignored
    # transaction root. A child directory would be misreported as a live
    # transaction by older versions of ``memory.py status``.
    return transaction_root(repo)


def validate_transaction_id(transaction_id: str) -> str:
    if not TRANSACTION_ID_RE.fullmatch(transaction_id) or transaction_id.startswith("_"):
        raise TransactionGCError(f"invalid transaction ID: {transaction_id!r}")
    return transaction_id


def source_path(repo: Path, transaction_id: str) -> Path:
    return transaction_root(repo) / validate_transaction_id(transaction_id)


def archive_path(repo: Path, transaction_id: str) -> Path:
    return archive_root(repo) / f"{ARCHIVE_PREFIX}{validate_transaction_id(transaction_id)}.tar.gz"


def receipt_path(repo: Path, transaction_id: str) -> Path:
    return archive_root(repo) / f"{ARCHIVE_PREFIX}{validate_transaction_id(transaction_id)}.receipt.json"


def lifecycle_path(repo: Path, transaction_id: str) -> Path:
    return archive_root(repo) / f"{ARCHIVE_PREFIX}{validate_transaction_id(transaction_id)}.lifecycle.json"


def refuse_during_publication(repo: Path) -> None:
    lock = repo / "memory/locks/publication.lock"
    if lock.exists():
        raise TransactionGCError(
            "memory publication lock exists; recover/finish publication before changing transaction storage"
        )


def list_transaction_ids(repo: Path) -> list[str]:
    root = transaction_root(repo)
    if not root.is_dir():
        return []
    return sorted(
        path.name for path in root.iterdir()
        if path.is_dir() and not path.name.startswith("_") and TRANSACTION_ID_RE.fullmatch(path.name)
    )


def _relative_record(path: Path, root: Path) -> dict[str, Any]:
    relative = path.relative_to(root).as_posix()
    metadata = path.lstat()
    if stat.S_ISLNK(metadata.st_mode):
        raise TransactionGCError(f"transaction contains a symlink: {root.name}/{relative}")
    if stat.S_ISDIR(metadata.st_mode):
        return {"path": relative + "/", "type": "directory", "mode": stat.S_IMODE(metadata.st_mode)}
    if not stat.S_ISREG(metadata.st_mode):
        raise TransactionGCError(f"transaction contains a non-regular entry: {root.name}/{relative}")
    return {
        "path": relative,
        "type": "file",
        "mode": stat.S_IMODE(metadata.st_mode),
        "size": metadata.st_size,
        "sha256": sha256_file(path),
    }


def tree_inventory(root: Path) -> list[dict[str, Any]]:
    if not root.is_dir() or root.is_symlink():
        raise TransactionGCError(f"transaction source is not a regular directory: {root}")
    return [_relative_record(path, root) for path in sorted(root.rglob("*"))]


def inventory_summary(entries: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    files = [entry for entry in entries if entry["type"] == "file"]
    return {
        "tree_sha256": sha256_bytes(canonical_json(list(entries))),
        "entry_count": len(entries),
        "file_count": len(files),
        "byte_count": sum(int(entry["size"]) for entry in files),
    }


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise TransactionGCError(f"invalid {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise TransactionGCError(f"invalid {label} {path}: expected a JSON object")
    return value


def transaction_manifest(repo: Path, transaction_id: str) -> dict[str, Any]:
    path = source_path(repo, transaction_id) / "manifest.json"
    manifest = load_json(path, "transaction manifest")
    if manifest.get("transaction_id") != transaction_id:
        raise TransactionGCError(f"transaction manifest identity mismatch: {transaction_id}")
    return manifest


def transaction_references(path: Path) -> set[str]:
    references: set[str] = set()
    for revision_path in path.glob("tasks/*/packet/revision.json"):
        revision = load_json(revision_path, "revision provenance")
        prior = revision.get("prior_transaction_id")
        if isinstance(prior, str) and TRANSACTION_ID_RE.fullmatch(prior):
            references.add(prior)
    for attestation_path in path.glob("attestations/*.json"):
        attestation = load_json(attestation_path, "writer attestation")
        for key in ("candidate_reuse", "reviewed_reuse"):
            reuse = attestation.get(key)
            prior = reuse.get("prior_transaction_id") if isinstance(reuse, dict) else None
            if isinstance(prior, str) and TRANSACTION_ID_RE.fullmatch(prior):
                references.add(prior)
    return references


def output_hashes(path: Path, manifest: Mapping[str, Any]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for task in manifest.get("writer_tasks", []):
        if not isinstance(task, dict) or not task.get("required"):
            continue
        repository_path = task.get("output_repository_path")
        staged_path = task.get("staged_output_path")
        if not isinstance(repository_path, str) or not isinstance(staged_path, str):
            continue
        candidate = path / staged_path
        if candidate.is_file() and not candidate.is_symlink():
            hashes[repository_path] = sha256_file(candidate)
    return hashes


def read_receipt(repo: Path, transaction_id: str) -> dict[str, Any] | None:
    path = receipt_path(repo, transaction_id)
    if not path.is_file():
        return None
    receipt = load_json(path, "archive receipt")
    if receipt.get("receipt_version") != RECEIPT_VERSION or receipt.get("transaction_id") != transaction_id:
        raise TransactionGCError(f"archive receipt identity/version mismatch: {transaction_id}")
    return receipt


def read_lifecycle(repo: Path, transaction_id: str) -> dict[str, Any] | None:
    path = lifecycle_path(repo, transaction_id)
    if not path.is_file():
        return None
    receipt = load_json(path, "lifecycle receipt")
    if receipt.get("lifecycle_version") != LIFECYCLE_VERSION or receipt.get("transaction_id") != transaction_id:
        raise TransactionGCError(f"lifecycle receipt identity/version mismatch: {transaction_id}")
    return receipt


def _safe_archive_member(member: tarfile.TarInfo, transaction_id: str) -> str | None:
    pure = PurePosixPath(member.name)
    if pure.is_absolute() or ".." in pure.parts or not pure.parts or pure.parts[0] != transaction_id:
        raise TransactionGCError(f"unsafe archive member path: {member.name!r}")
    if member.issym() or member.islnk() or not (member.isdir() or member.isreg()):
        raise TransactionGCError(f"unsupported archive member type: {member.name!r}")
    if len(pure.parts) == 1:
        if not member.isdir():
            raise TransactionGCError("archive transaction root is not a directory")
        return None
    return PurePosixPath(*pure.parts[1:]).as_posix()


def archive_inventory(path: Path, transaction_id: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    try:
        archive = tarfile.open(path, "r:gz")
    except (OSError, tarfile.TarError) as exc:
        raise TransactionGCError(f"cannot read archive {path}: {exc}") from exc
    with archive:
        for member in archive:
            relative = _safe_archive_member(member, transaction_id)
            if relative is None:
                continue
            key = relative + "/" if member.isdir() else relative
            if key in seen:
                raise TransactionGCError(f"duplicate archive member: {member.name}")
            seen.add(key)
            if member.isdir():
                records.append({"path": key, "type": "directory", "mode": member.mode})
                continue
            stream = archive.extractfile(member)
            if stream is None:
                raise TransactionGCError(f"archive file has no payload: {member.name}")
            digest = hashlib.sha256()
            size = 0
            with stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
                    size += len(chunk)
            if size != member.size:
                raise TransactionGCError(f"archive member size mismatch: {member.name}")
            records.append({
                "path": relative, "type": "file", "mode": member.mode,
                "size": size, "sha256": digest.hexdigest(),
            })
    return sorted(records, key=lambda item: item["path"])


def verify_archive_file(repo: Path, transaction_id: str) -> dict[str, Any]:
    receipt = read_receipt(repo, transaction_id)
    if receipt is None:
        raise TransactionGCError(f"archive receipt is missing: {transaction_id}")
    path = archive_path(repo, transaction_id)
    if not path.is_file() or path.is_symlink():
        raise TransactionGCError(f"archive is missing or not a regular file: {transaction_id}")
    if path.stat().st_size != receipt.get("archive_bytes"):
        raise TransactionGCError(f"archive size mismatch: {transaction_id}")
    if sha256_file(path) != receipt.get("archive_sha256"):
        raise TransactionGCError(f"archive hash mismatch: {transaction_id}")
    return {"archive_sha256": receipt["archive_sha256"], "archive_bytes": path.stat().st_size}


def verify_archive(repo: Path, transaction_id: str) -> dict[str, Any]:
    archive_file = verify_archive_file(repo, transaction_id)
    receipt = read_receipt(repo, transaction_id)
    assert receipt is not None  # Established by verify_archive_file.
    path = archive_path(repo, transaction_id)
    summary = inventory_summary(archive_inventory(path, transaction_id))
    for field in ("tree_sha256", "entry_count", "file_count", "byte_count"):
        if summary[field] != receipt.get(field):
            raise TransactionGCError(f"archive tree receipt mismatch for {transaction_id}: {field}")
    return {**summary, **archive_file}


def _load_state(repo: Path) -> tuple[dict[str, Any], str]:
    path = repo / "memory/_meta/state.json"
    if not path.is_file():
        return {}, sha256_bytes(b"")
    data = path.read_bytes()
    try:
        state = json.loads(data)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise TransactionGCError(f"invalid memory state: {exc}") from exc
    if not isinstance(state, dict):
        raise TransactionGCError("invalid memory state: expected a JSON object")
    return state, sha256_bytes(data)


def _git_head(repo: Path) -> str | None:
    proc = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    return proc.stdout.decode("ascii", "replace").strip() if proc.returncode == 0 else None


def _transaction_metadata(repo: Path, transaction_id: str) -> dict[str, Any]:
    source = source_path(repo, transaction_id)
    receipt = read_receipt(repo, transaction_id)
    if source.is_dir():
        manifest = transaction_manifest(repo, transaction_id)
        references = sorted(transaction_references(source))
        outputs = output_hashes(source, manifest)
        summary = inventory_summary(tree_inventory(source))
        failure = (source / "failure.json").is_file()
    elif receipt is not None:
        manifest = {
            "transaction_id": transaction_id,
            "target_commit": receipt.get("target_commit"),
            "created_at": receipt.get("transaction_created_at"),
        }
        references = list(receipt.get("references_out", []))
        outputs = dict(receipt.get("output_hashes", {}))
        summary = {
            key: receipt.get(key) for key in ("tree_sha256", "entry_count", "file_count", "byte_count")
        }
        failure = bool(receipt.get("failure_recorded"))
    else:
        raise TransactionGCError(f"transaction and archive are both missing: {transaction_id}")
    return {
        "id": transaction_id,
        "source_exists": source.is_dir(),
        "manifest": manifest,
        "references_out": references,
        "output_hashes": outputs,
        "failure_recorded": failure,
        **summary,
    }


def all_known_transaction_ids(repo: Path) -> list[str]:
    result = set(list_transaction_ids(repo))
    root = archive_root(repo)
    if root.is_dir():
        for receipt in root.glob(f"{ARCHIVE_PREFIX}*.receipt.json"):
            transaction_id = receipt.name.removeprefix(ARCHIVE_PREFIX).removesuffix(".receipt.json")
            if TRANSACTION_ID_RE.fullmatch(transaction_id):
                result.add(transaction_id)
    return sorted(result)


def classify_transactions(repo: Path) -> dict[str, dict[str, Any]]:
    state, _ = _load_state(repo)
    head = _git_head(repo)
    metadata = {transaction_id: _transaction_metadata(repo, transaction_id) for transaction_id in all_known_transaction_ids(repo)}
    incoming: dict[str, set[str]] = {transaction_id: set() for transaction_id in metadata}
    for child_id, item in metadata.items():
        for prior_id in item["references_out"]:
            incoming.setdefault(prior_id, set()).add(child_id)
    state_pages = state.get("pages", {}) if isinstance(state.get("pages"), dict) else {}
    generation = state.get("generation")
    result: dict[str, dict[str, Any]] = {}
    for transaction_id, item in metadata.items():
        state_matches = sorted(
            path for path, digest in item["output_hashes"].items()
            if isinstance(state_pages.get(path), dict) and state_pages[path].get("sha256") == digest
        )
        protected: list[str] = []
        if transaction_id == generation:
            protected.append("current_state_generation")
        if state_matches:
            protected.append("serves_current_page_hash")
        if incoming.get(transaction_id):
            protected.append("referenced_by_transaction")
        lifecycle = read_lifecycle(repo, transaction_id)
        target = item["manifest"].get("target_commit")
        if lifecycle is not None and lifecycle.get("status") == "finalized":
            lifecycle_status = "finalized_current" if state_matches or transaction_id == generation else "finalized_superseded"
        elif item["failure_recorded"]:
            lifecycle_status = "failed"
        elif target and head and target != head:
            lifecycle_status = "inactive_target"
        else:
            lifecycle_status = "unclassified_current_target"
        archive_status = "missing"
        receipt = read_receipt(repo, transaction_id)
        if receipt is not None:
            try:
                verify_archive_file(repo, transaction_id)
                archive_status = "hash_verified"
            except TransactionGCError:
                archive_status = "invalid"
        prune_eligible = (
            item["source_exists"] and archive_status == "hash_verified" and not protected
            and lifecycle_status in {"finalized_superseded", "failed", "inactive_target"}
        )
        result[transaction_id] = {
            "transaction_id": transaction_id,
            "source_exists": item["source_exists"],
            "target_commit": target,
            "head_matches_target": bool(target and head == target),
            "tree_sha256": item["tree_sha256"],
            "file_count": item["file_count"],
            "byte_count": item["byte_count"],
            "references_out": item["references_out"],
            "references_in": sorted(incoming.get(transaction_id, set())),
            "current_state_page_matches": state_matches,
            "protected_reasons": protected,
            "lifecycle_status": lifecycle_status,
            "archive_status": archive_status,
            "prune_eligible": prune_eligible,
        }
    return result


def archive_transaction(repo: Path, transaction_id: str, *, apply: bool) -> dict[str, Any]:
    transaction_id = validate_transaction_id(transaction_id)
    source = source_path(repo, transaction_id)
    if not source.is_dir():
        raise TransactionGCError(f"transaction source is missing: {transaction_id}")
    manifest = transaction_manifest(repo, transaction_id)
    before_entries = tree_inventory(source)
    before = inventory_summary(before_entries)
    plan = {
        "action": "archive", "apply": apply, "transaction_id": transaction_id,
        "source_bytes": before["byte_count"], "source_files": before["file_count"],
        "source_tree_sha256": before["tree_sha256"],
        "archive": archive_path(repo, transaction_id).relative_to(repo).as_posix(),
    }
    if not apply:
        return plan
    refuse_during_publication(repo)
    destination = archive_path(repo, transaction_id)
    existing_receipt = read_receipt(repo, transaction_id)

    def build_receipt() -> dict[str, Any]:
        return {
            "receipt_version": RECEIPT_VERSION,
            "transaction_id": transaction_id,
            "archived_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "transaction_created_at": manifest.get("created_at"),
            "target_commit": manifest.get("target_commit"),
            "manifest_sha256": sha256_file(source / "manifest.json"),
            "seal_sha256": sha256_file(source / "seal.json") if (source / "seal.json").is_file() else None,
            **before,
            "archive_path": destination.relative_to(repo).as_posix(),
            "archive_sha256": sha256_file(destination),
            "archive_bytes": destination.stat().st_size,
            "references_out": sorted(transaction_references(source)),
            "output_hashes": output_hashes(source, manifest),
            "failure_recorded": (source / "failure.json").is_file(),
        }

    if destination.exists() or existing_receipt is not None:
        if destination.is_file() and existing_receipt is not None:
            verified = verify_archive(repo, transaction_id)
            if verified["tree_sha256"] == before["tree_sha256"]:
                return {**plan, "status": "already_verified", **verified}
        elif destination.is_file() and existing_receipt is None:
            # Recover the narrow crash window after the verified archive was
            # renamed into place but before its receipt was committed.
            archived = inventory_summary(archive_inventory(destination, transaction_id))
            if archived == before:
                write_json(receipt_path(repo, transaction_id), build_receipt())
                verified = verify_archive(repo, transaction_id)
                return {**plan, "status": "recovered_receipt_and_verified", **verified}
        raise TransactionGCError(f"archive or receipt already exists with different content: {transaction_id}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{transaction_id}.", suffix=".tar.gz", dir=destination.parent)
    os.close(descriptor)
    temporary_path = Path(temporary)
    try:
        with tarfile.open(temporary_path, "w:gz", compresslevel=6) as archive:
            archive.add(source, arcname=transaction_id, recursive=True)
        after = inventory_summary(tree_inventory(source))
        if after != before:
            raise TransactionGCError(f"transaction changed while it was being archived: {transaction_id}")
        archived = inventory_summary(archive_inventory(temporary_path, transaction_id))
        if archived != before:
            raise TransactionGCError(f"archive verification did not reproduce the source tree: {transaction_id}")
        os.replace(temporary_path, destination)
        write_json(receipt_path(repo, transaction_id), build_receipt())
        verified = verify_archive(repo, transaction_id)
        return {**plan, "status": "archived_and_verified", **verified}
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def record_finalized(repo: Path, transaction_id: str, *, apply: bool) -> dict[str, Any]:
    transaction_id = validate_transaction_id(transaction_id)
    source = source_path(repo, transaction_id)
    if not source.is_dir():
        raise TransactionGCError(f"transaction source is missing: {transaction_id}")
    state, state_sha256 = _load_state(repo)
    if state.get("generation") != transaction_id or state.get("last_result", {}).get("transaction_id") != transaction_id:
        raise TransactionGCError(
            f"cannot prove finalization: current state generation/last_result is not {transaction_id}"
        )
    manifest = transaction_manifest(repo, transaction_id)
    outputs = output_hashes(source, manifest)
    state_pages = state.get("pages", {})
    unmatched = sorted(
        path for path, digest in outputs.items()
        if not isinstance(state_pages.get(path), dict) or state_pages[path].get("sha256") != digest
    )
    if unmatched:
        raise TransactionGCError(f"cannot prove finalization; staged outputs do not match state: {unmatched}")
    tree = inventory_summary(tree_inventory(source))
    record = {
        "lifecycle_version": LIFECYCLE_VERSION,
        "transaction_id": transaction_id,
        "status": "finalized",
        "recorded_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "state_sha256": state_sha256,
        "target_commit": manifest.get("target_commit"),
        "manifest_sha256": sha256_file(source / "manifest.json"),
        "seal_sha256": sha256_file(source / "seal.json") if (source / "seal.json").is_file() else None,
        "transaction_tree_sha256": tree["tree_sha256"],
        "published_output_hashes": outputs,
    }
    if apply:
        refuse_during_publication(repo)
        path = lifecycle_path(repo, transaction_id)
        existing = read_lifecycle(repo, transaction_id)
        if existing is not None:
            stable_fields = (
                "status", "state_sha256", "target_commit", "manifest_sha256", "seal_sha256",
                "transaction_tree_sha256", "published_output_hashes",
            )
            if any(existing.get(field) != record.get(field) for field in stable_fields):
                raise TransactionGCError(f"different lifecycle receipt already exists: {transaction_id}")
            record = existing
        else:
            write_json(path, record)
    return {"action": "record-finalized", "apply": apply, **record}


def prune_transaction(repo: Path, transaction_id: str, *, apply: bool) -> dict[str, Any]:
    transaction_id = validate_transaction_id(transaction_id)
    classifications = classify_transactions(repo)
    item = classifications.get(transaction_id)
    if item is None or not item["source_exists"]:
        raise TransactionGCError(f"transaction source is missing: {transaction_id}")
    verified = verify_archive(repo, transaction_id)
    current_tree = inventory_summary(tree_inventory(source_path(repo, transaction_id)))
    if current_tree["tree_sha256"] != verified["tree_sha256"]:
        raise TransactionGCError(f"source changed after verified archive: {transaction_id}")
    if not item["prune_eligible"]:
        reasons = item["protected_reasons"] or [item["lifecycle_status"]]
        raise TransactionGCError(f"transaction is not prune-eligible: {transaction_id}: {', '.join(reasons)}")
    result = {
        "action": "prune", "apply": apply, "transaction_id": transaction_id,
        "recoverable_from": archive_path(repo, transaction_id).relative_to(repo).as_posix(),
        "reclaimed_bytes": current_tree["byte_count"],
    }
    if apply:
        refuse_during_publication(repo)
        shutil.rmtree(source_path(repo, transaction_id))
        result["status"] = "source_removed_archive_retained"
    else:
        result["status"] = "dry_run"
    return result


def restore_transaction(repo: Path, transaction_id: str, *, apply: bool) -> dict[str, Any]:
    transaction_id = validate_transaction_id(transaction_id)
    verified = verify_archive(repo, transaction_id)
    destination = source_path(repo, transaction_id)
    result = {
        "action": "restore", "apply": apply, "transaction_id": transaction_id,
        "restored_bytes": verified["byte_count"], "tree_sha256": verified["tree_sha256"],
    }
    if destination.exists():
        if destination.is_dir() and inventory_summary(tree_inventory(destination))["tree_sha256"] == verified["tree_sha256"]:
            return {**result, "status": "already_present_and_verified"}
        raise TransactionGCError(f"restore destination exists with different content: {transaction_id}")
    if not apply:
        return {**result, "status": "dry_run"}
    refuse_during_publication(repo)
    root = transaction_root(repo)
    temporary = Path(tempfile.mkdtemp(prefix=f".{transaction_id}.restore.", dir=root))
    try:
        with tarfile.open(archive_path(repo, transaction_id), "r:gz") as archive:
            members = archive.getmembers()
            for member in members:
                _safe_archive_member(member, transaction_id)
            # Every member was constrained above to the transaction prefix and
            # a regular file/directory.  Returning the validated member rather
            # than using tarfile's ``data`` filter preserves the original mode
            # exactly (the standard filter deliberately strips group-write).
            archive.extractall(temporary, members=members, filter=lambda member, _: member)
        restored = temporary / transaction_id
        if inventory_summary(tree_inventory(restored))["tree_sha256"] != verified["tree_sha256"]:
            raise TransactionGCError(f"restored transaction failed tree verification: {transaction_id}")
        os.replace(restored, destination)
        return {**result, "status": "restored_and_verified"}
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def _print_results(results: Any, as_json: bool) -> None:
    if as_json:
        print(json.dumps(results, indent=2, sort_keys=True))
        return
    rows = results.values() if isinstance(results, dict) and all(isinstance(v, dict) for v in results.values()) else results
    if isinstance(rows, dict):
        rows = [rows]
    for item in rows:
        if "lifecycle_status" in item:
            protected = ",".join(item["protected_reasons"]) or "none"
            print(
                f"{item['transaction_id']}: lifecycle={item['lifecycle_status']} "
                f"archive={item['archive_status']} protected={protected} "
                f"prune_eligible={str(item['prune_eligible']).lower()} bytes={item['byte_count']}"
            )
        else:
            print(json.dumps(item, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, help="repository root (normally auto-detected)")
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    subparsers = parser.add_subparsers(dest="command", required=True)
    status_parser = subparsers.add_parser("status", help="classify transaction lifecycle and protection")
    status_parser.add_argument("transactions", nargs="*", help="optional transaction IDs")
    for command, help_text in (
        ("archive", "create and verify a lossless archive; never removes the source"),
        ("prune", "remove an eligible source tree only after verified archival"),
        ("restore", "restore and verify a transaction from its archive"),
        ("record-finalized", "record finalization only when current state proves it"),
    ):
        child = subparsers.add_parser(command, help=help_text)
        child.add_argument("transactions", nargs="+", help="transaction IDs")
        child.add_argument("--apply", action="store_true", help="perform the operation; default is dry-run")
    verify_parser = subparsers.add_parser("verify", help="fully verify archive bytes and tree hashes")
    verify_parser.add_argument("transactions", nargs="+", help="transaction IDs")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        repo = args.repo.resolve() if args.repo else find_repo()
        if args.command == "status":
            report = classify_transactions(repo)
            if args.transactions:
                requested = {validate_transaction_id(item) for item in args.transactions}
                absent = sorted(requested - set(report))
                if absent:
                    raise TransactionGCError(f"unknown transactions: {', '.join(absent)}")
                report = {key: value for key, value in report.items() if key in requested}
            _print_results(report, args.json)
            return 0
        operation = {
            "archive": archive_transaction,
            "prune": prune_transaction,
            "restore": restore_transaction,
            "record-finalized": record_finalized,
        }.get(args.command)
        if args.command == "verify":
            results = [
                {"action": "verify", "transaction_id": item, **verify_archive(repo, item)}
                for item in args.transactions
            ]
        elif operation is not None:
            if args.command == "prune" and args.apply:
                # Validate the complete requested batch before removing the
                # first source tree. This avoids a predictable half-pruned
                # batch when a later transaction is protected or unarchived.
                for item in args.transactions:
                    prune_transaction(repo, item, apply=False)
            results = [operation(repo, item, apply=args.apply) for item in args.transactions]
        else:  # pragma: no cover - argparse owns command validation
            raise TransactionGCError(f"unknown command: {args.command}")
        _print_results(results, args.json)
        return 0
    except (TransactionGCError, OSError, tarfile.TarError, ValueError) as exc:
        print(f"transaction gc: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
