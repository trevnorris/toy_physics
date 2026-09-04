#!/usr/bin/env python3
"""Small, content-hash based freshness helper for the research memory."""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


MEMORY = Path(__file__).resolve().parent
REPO = MEMORY.parent
CATALOG_PATH = MEMORY / "catalog.json"
STATE_PATH = MEMORY / "_meta" / "state.json"


class MemoryError(Exception):
    pass


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise MemoryError(f"cannot read {path.relative_to(REPO)}: {exc}") from exc


def catalog() -> dict:
    data = read_json(CATALOG_PATH)
    if data.get("version") != 1:
        raise MemoryError("unsupported catalog version")
    return data


def state() -> dict:
    if not STATE_PATH.exists():
        return {"version": 1, "units": {}, "pages": {}}
    data = read_json(STATE_PATH)
    if data.get("version") != 1:
        raise MemoryError("unsupported state version")
    data.setdefault("units", {})
    data.setdefault("pages", {})
    return data


def write_state(data: dict) -> None:
    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = STATE_PATH.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, STATE_PATH)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def tracked_paths() -> set[str]:
    command = ["git", "ls-files", "-z", "--", "research", "software"]
    result = subprocess.run(command, cwd=REPO, check=True, stdout=subprocess.PIPE)
    return {item.decode() for item in result.stdout.split(b"\0") if item}


def is_excluded(path: str, data: dict) -> bool:
    return any(path == prefix.rstrip("/") or path.startswith(prefix) for prefix in data["excluded_prefixes"])


def resolve_patterns(patterns: list[str], tracked: set[str], data: dict) -> tuple[list[str], list[str]]:
    selected: set[str] = set()
    errors: list[str] = []
    for pattern in patterns:
        matches = sorted(path for path in tracked if fnmatch.fnmatchcase(path, pattern))
        if not matches:
            errors.append(f"no tracked source matches {pattern}")
            continue
        for path in matches:
            if is_excluded(path, data):
                errors.append(f"excluded source selected: {path}")
            else:
                selected.add(path)
    return sorted(selected), errors


def unit_snapshot(unit_id: str, spec: dict, tracked: set[str], data: dict) -> tuple[dict | None, list[str]]:
    paths, errors = resolve_patterns(spec["sources"], tracked, data)
    files: dict[str, str] = {}
    for relative in paths:
        path = REPO / relative
        if not path.is_file():
            errors.append(f"tracked source is unavailable: {relative}")
            continue
        files[relative] = sha256_file(path)
    if errors:
        return None, [f"{unit_id}: {error}" for error in errors]
    return {"content_hash": stable_hash(files), "files": files}, []


def page_hash(relative: str) -> str | None:
    path = REPO / relative
    return sha256_file(path) if path.is_file() else None


def input_hash(inputs: list[str], data: dict, snapshots: dict[str, dict | None]) -> tuple[str | None, list[str]]:
    values: dict[str, object] = {}
    errors: list[str] = []
    for item in inputs:
        if item in data["units"]:
            snapshot = snapshots.get(item)
            capsule_hash = page_hash(data["units"][item]["page"])
            if snapshot is None or capsule_hash is None:
                errors.append(item)
            else:
                values[item] = {"source": snapshot["content_hash"], "page": capsule_hash}
        elif item in data["pages"]:
            digest = page_hash(data["pages"][item]["path"])
            if digest is None:
                errors.append(item)
            else:
                values[item] = digest
        else:
            errors.append(f"unknown:{item}")
    return (stable_hash(values), []) if not errors else (None, errors)


def inspect() -> tuple[dict, list[str]]:
    data = catalog()
    saved = state()
    tracked = tracked_paths()
    errors: list[str] = []
    snapshots: dict[str, dict | None] = {}
    units: dict[str, dict] = {}

    for unit_id, spec in data["units"].items():
        snapshot, unit_errors = unit_snapshot(unit_id, spec, tracked, data)
        snapshots[unit_id] = snapshot
        errors.extend(unit_errors)
        actual_page_hash = page_hash(spec["page"])
        prior = saved["units"].get(unit_id)
        if snapshot is None or actual_page_hash is None:
            label = "missing"
        elif prior is None:
            label = "new"
        elif prior.get("content_hash") != snapshot["content_hash"]:
            label = "changed"
        elif prior.get("page_hash") != actual_page_hash:
            label = "edited"
        else:
            label = "current"
        units[unit_id] = {"status": label, "page": spec["page"], "snapshot": snapshot}

    pages: dict[str, dict] = {}
    for page_id, spec in data["pages"].items():
        digest, missing_inputs = input_hash(spec["inputs"], data, snapshots)
        actual_page_hash = page_hash(spec["path"])
        prior = saved["pages"].get(page_id)
        upstream_not_current = [
            item for item in spec["inputs"]
            if (item in units and units[item]["status"] != "current")
            or (item in pages and pages[item]["status"] != "current")
        ]
        if actual_page_hash is None:
            label = "missing"
        elif missing_inputs or upstream_not_current:
            label = "blocked"
        elif prior is None:
            label = "new"
        elif prior.get("input_hash") != digest:
            label = "stale"
        elif prior.get("page_hash") != actual_page_hash:
            label = "edited"
        else:
            label = "current"
        pages[page_id] = {
            "status": label,
            "path": spec["path"],
            "input_hash": digest,
            "blocked_by": sorted(set(missing_inputs + upstream_not_current)),
        }

    return {"units": units, "pages": pages, "state": saved, "catalog": data}, errors


def command_status(as_json: bool) -> int:
    report, errors = inspect()
    compact = {
        "units": {key: value["status"] for key, value in report["units"].items()},
        "pages": {key: value["status"] for key, value in report["pages"].items()},
        "errors": errors,
    }
    if as_json:
        print(json.dumps(compact, indent=2, sort_keys=True))
    else:
        for group in ("units", "pages"):
            counts: dict[str, int] = {}
            for label in compact[group].values():
                counts[label] = counts.get(label, 0) + 1
            print(f"{group}: " + ", ".join(f"{key}={counts[key]}" for key in sorted(counts)))
            for item, label in compact[group].items():
                if label != "current":
                    print(f"  {label:8} {item}")
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
    return 1 if errors else 0


def command_mark_units(ids: list[str]) -> int:
    data = catalog()
    saved = state()
    tracked = tracked_paths()
    for unit_id in ids:
        if unit_id not in data["units"]:
            raise MemoryError(f"unknown unit: {unit_id}")
        spec = data["units"][unit_id]
        snapshot, errors = unit_snapshot(unit_id, spec, tracked, data)
        if errors or snapshot is None:
            raise MemoryError("; ".join(errors))
        digest = page_hash(spec["page"])
        if digest is None:
            raise MemoryError(f"capsule does not exist: {spec['page']}")
        saved["units"][unit_id] = {
            "content_hash": snapshot["content_hash"],
            "page_hash": digest,
        }
    write_state(saved)
    print("marked units: " + ", ".join(ids))
    return 0


def command_mark_pages(ids: list[str]) -> int:
    for page_id in ids:
        report, errors = inspect()
        if errors:
            raise MemoryError("; ".join(errors))
        if page_id not in report["catalog"]["pages"]:
            raise MemoryError(f"unknown page: {page_id}")
        page = report["pages"][page_id]
        if page["blocked_by"]:
            raise MemoryError(f"{page_id} has stale inputs: {', '.join(page['blocked_by'])}")
        digest = page_hash(page["path"])
        if digest is None or page["input_hash"] is None:
            raise MemoryError(f"page is unavailable: {page['path']}")
        saved = report["state"]
        saved["pages"][page_id] = {
            "input_hash": page["input_hash"],
            "page_hash": digest,
        }
        write_state(saved)
    print("marked pages: " + ", ".join(ids))
    return 0


def command_lint() -> int:
    report, errors = inspect()
    data = report["catalog"]
    known = set(data["units"]) | set(data["pages"])
    for page_id, spec in data["pages"].items():
        for item in spec["inputs"]:
            if item not in known:
                errors.append(f"{page_id}: unknown input {item}")
    forbidden = ("generated_from_commit:", "blob_oid:", "processed_commit:", "target_commit:", "transaction_id:", "Git object:")
    for path in MEMORY.rglob("*.md"):
        if "transactions" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        for marker in forbidden:
            if marker in text:
                errors.append(f"{path.relative_to(REPO)}: obsolete metadata {marker}")
    if errors:
        for error in sorted(set(errors)):
            print(f"error: {error}", file=sys.stderr)
        return 1
    print(f"lint: ok ({len(data['units'])} source units, {len(data['pages'])} dependent pages)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    status_parser = subparsers.add_parser("status", help="show changed source units and affected pages")
    status_parser.add_argument("--json", action="store_true")
    mark_units = subparsers.add_parser("mark-units", help="record capsules updated from current source content")
    mark_units.add_argument("ids", nargs="+")
    mark_pages = subparsers.add_parser("mark-pages", help="record dependent pages updated from current inputs")
    mark_pages.add_argument("ids", nargs="+")
    subparsers.add_parser("lint", help="check the catalog and obsolete metadata")
    args = parser.parse_args()
    if args.command == "status":
        return command_status(args.json)
    if args.command == "mark-units":
        return command_mark_units(args.ids)
    if args.command == "mark-pages":
        return command_mark_pages(args.ids)
    return command_lint()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (MemoryError, subprocess.CalledProcessError) as exc:
        print(f"memory: {exc}", file=sys.stderr)
        raise SystemExit(1)
