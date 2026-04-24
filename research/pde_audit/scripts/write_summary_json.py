#!/usr/bin/env python3
"""Build machine-readable summaries from captured PDE audit outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


HEADER_RE = re.compile(r"^# ([A-Za-z ]+):\s*(.*)$")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_output(path: Path) -> Dict[str, Any]:
    meta: Dict[str, Any] = {
        "name": path.stem,
        "output_file": str(path),
        "output_sha256": sha256_file(path),
    }
    exit_code_seen: Optional[int] = None
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            m = HEADER_RE.match(line)
            if m:
                key = m.group(1).strip().lower().replace(" ", "_")
                meta[key] = m.group(2).strip()
                continue
            if line.startswith("EXIT_CODE:"):
                try:
                    exit_code_seen = int(line.split(":", 1)[1].strip())
                except ValueError:
                    exit_code_seen = None
    if "exit_code" in meta:
        try:
            meta["exit_code"] = int(str(meta["exit_code"]))
        except ValueError:
            pass
    if exit_code_seen is not None:
        meta["exit_code_trailer"] = exit_code_seen
    exit_code = meta.get("exit_code", exit_code_seen)
    meta["status"] = "PASS" if exit_code == 0 else "FAIL"
    return meta


def load_manifest(path: Optional[str]) -> Optional[Dict[str, Any]]:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return None
    with p.open("r", encoding="utf-8") as f:
        manifest = json.load(f)
    manifest["manifest_path"] = str(p)
    manifest["manifest_sha256"] = sha256_file(p)
    return manifest


def hash_files(paths: Iterable[Path]) -> List[Dict[str, str]]:
    items = []
    for path in sorted(paths):
        if path.is_file():
            items.append({"path": str(path), "sha256": sha256_file(path)})
    return items


def python_package_version(name: str) -> Optional[str]:
    try:
        import importlib.metadata as metadata

        return metadata.version(name)
    except Exception:
        return None


def parse_environment_capture(path: Path) -> Dict[str, str]:
    allowed_keys = {
        "MathematicaKernelVersion",
        "MathematicaVersionNumber",
        "SystemID",
        "MachineName",
        "ProcessorType",
    }
    parsed: Dict[str, str] = {}
    if not path.exists():
        return parsed
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("EXIT_CODE:"):
                continue
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            key = key.strip()
            if key not in allowed_keys:
                continue
            parsed[key] = value.strip().strip('"')
    return parsed


def environment_summary(suite: str, output_dir: Optional[Path] = None) -> Dict[str, Any]:
    env: Dict[str, Any] = {
        "suite": suite,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    if suite == "python":
        env["packages"] = {
            "numpy": python_package_version("numpy"),
            "scipy": python_package_version("scipy"),
            "sympy": python_package_version("sympy"),
        }
    if suite == "mathematica" and output_dir is not None:
        env_path = output_dir / "_environment.txt"
        env["mathematica"] = parse_environment_capture(env_path)
        if env_path.exists():
            env["mathematica_environment_output"] = {
                "path": str(env_path),
                "sha256": sha256_file(env_path),
            }
    return env


def git_head() -> Optional[str]:
    try:
        out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return None
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Write PDE audit _summary.json")
    parser.add_argument("--suite", required=True, choices=["python", "mathematica"])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--artifact-dir", default=None)
    parser.add_argument("--filter", default="", help="Only summarize stage outputs whose stem contains this filter")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    outputs = [
        parse_output(path)
        for path in sorted(output_dir.glob("stage_v2_*.txt"))
        if path.name != "_summary.txt" and (not args.filter or args.filter in path.stem)
    ]
    total = len(outputs)
    passed = sum(1 for item in outputs if item["status"] == "PASS")
    failed = total - passed

    summary: Dict[str, Any] = {
        "schema": "pde_audit_summary/v1",
        "suite": args.suite,
        "output_dir": str(output_dir),
        "git_head": git_head(),
        "environment": environment_summary(args.suite, output_dir),
        "total": total,
        "passed": passed,
        "failed": failed,
        "outputs": outputs,
    }

    manifest = load_manifest(args.manifest)
    if manifest is not None:
        summary["fixture_manifest"] = manifest

    if args.artifact_dir:
        artifact_dir = Path(args.artifact_dir)
        summary["artifacts"] = hash_files(artifact_dir.glob("*.json"))

    out_path = output_dir / "_summary.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

    print(f"Wrote {out_path}")
    print(f"TOTAL: {total}  PASS: {passed}  FAIL: {failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
