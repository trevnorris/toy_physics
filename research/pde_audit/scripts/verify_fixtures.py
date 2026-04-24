#!/usr/bin/env python3
"""Verify PDE audit fixture hashes and manifest coverage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify PDE audit fixture manifest")
    parser.add_argument(
        "--fixture-dir",
        default=str(Path(__file__).resolve().parent / "fixtures"),
        help="Directory containing fixture JSON files and MANIFEST.json",
    )
    args = parser.parse_args()

    fixture_dir = Path(args.fixture_dir)
    manifest_path = fixture_dir / "MANIFEST.json"
    manifest = load_json(manifest_path)
    fixtures = manifest.get("fixtures", [])

    expected = {entry["filename"]: entry for entry in fixtures}
    actual_json = {
        path.relative_to(fixture_dir).as_posix()
        for path in fixture_dir.rglob("*.json")
        if path.name != "MANIFEST.json"
    }

    ok = True
    print("PDE audit fixture manifest verification")
    print("=" * 48)
    print(f"manifest: {manifest_path}")

    missing = sorted(set(expected) - actual_json)
    extra = sorted(actual_json - set(expected))
    if missing:
        ok = False
        print("missing fixtures:")
        for name in missing:
            print(f"  - {name}")
    if extra:
        ok = False
        print("extra unmanifested fixtures:")
        for name in extra:
            print(f"  - {name}")

    for name, entry in sorted(expected.items()):
        path = fixture_dir / name
        if not path.exists():
            continue
        actual = sha256_file(path)
        expected_hash = entry["sha256"]
        status = "PASS" if actual == expected_hash else "FAIL"
        print(f"{status}  {name}  {actual}")
        if actual != expected_hash:
            ok = False

        try:
            payload = load_json(path)
        except Exception as exc:
            print(f"FAIL  {name}  invalid JSON: {exc}")
            ok = False
            continue
        schema = payload.get("schema")
        if schema != entry.get("schema"):
            print(f"FAIL  {name}  schema {schema!r} != {entry.get('schema')!r}")
            ok = False

    print("")
    print("fixture_manifest_pass:", ok)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
