#!/usr/bin/env python3
"""Measured, location-independent evaluated-code closure guard for U2."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


CODE_SUFFIXES = {".py", ".pyc", ".so", ".sh", ".wl", ".wls", ".m", ".mx"}


class GuardFailure(RuntimeError):
    def __init__(self, assert_id: str, detail: str):
        super().__init__(detail); self.assert_id = assert_id; self.detail = detail


def require(condition: bool, assert_id: str, detail: str) -> None:
    if not condition:
        raise GuardFailure(assert_id, detail)


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_blob(repo: Path, anchor: str, rel_path: str) -> bytes | None:
    result = subprocess.run(
        ["git", "cat-file", "blob", f"{anchor}:{rel_path}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    return result.stdout if result.returncode == 0 else None


def is_under(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root); return True
    except ValueError:
        return False


def extract_trace_paths(trace: Path, cwd: Path) -> tuple[list[Path], list[Path]]:
    paths: set[Path] = set(); direct: set[Path] = set()
    quoted = re.compile(r'"((?:[^"\\]|\\.)*)"')
    for line in trace.read_text(encoding="utf-8", errors="replace").splitlines():
        if not any(token in line for token in ("execve(", "open(", "openat(", "openat2(")):
            continue
        if "execve(" not in line and "O_WRONLY" in line and "O_RDWR" not in line:
            continue
        matches = quoted.findall(line)
        if not matches:
            continue
        raw = bytes(matches[0], "utf-8").decode("unicode_escape")
        if not raw or raw.startswith(("/proc/", "/dev/")):
            continue
        path = Path(raw)
        if not path.is_absolute(): path = cwd / path
        paths.add(path.resolve(strict=False))
        if "execve(" in line:
            for argument in matches:
                candidate = Path(bytes(argument, "utf-8").decode("unicode_escape"))
                if candidate.suffix.lower() not in CODE_SUFFIXES:
                    continue
                if not candidate.is_absolute(): candidate = cwd / candidate
                direct.add(candidate.resolve(strict=False))
    return sorted(paths), sorted(direct)


def validate_code_path(
    observed_path: Path, repo: Path, anchor: str, toolchain_roots: list[Path],
    content_path: Path | None = None, logical_repo_path: str | None = None,
) -> dict[str, Any]:
    path = observed_path.resolve(strict=False)
    source = content_path or observed_path
    content = source.read_bytes() if source.is_file() else b""
    if logical_repo_path is not None:
        blob = git_blob(repo, anchor, logical_repo_path)
        require(blob is not None, "ASSERT_EXECUTED_PATH_AT_ANCHOR", f"{logical_repo_path} absent from anchor")
        require(hashlib.sha256(content).digest() == hashlib.sha256(blob).digest(), "ASSERT_EXECUTED_BYTES_MATCH_ANCHOR", f"trace-held bytes differ from {anchor}:{logical_repo_path}")
        return {"path": str(path), "classification": "anchored_task_code", "logical_path": logical_repo_path}
    if is_under(path, repo):
        rel = path.relative_to(repo).as_posix(); blob = git_blob(repo, anchor, rel)
        require(blob is not None, "ASSERT_EVALUATED_CODE_CLOSURE", f"repository task code absent from anchor: {rel}")
        require(hashlib.sha256(content).digest() == hashlib.sha256(blob).digest(), "ASSERT_EXECUTED_BYTES_MATCH_ANCHOR", f"executed repository bytes differ from anchor: {rel}")
        return {"path": str(path), "classification": "anchored_task_code", "logical_path": rel}
    for root in toolchain_roots:
        if is_under(path, root):
            return {"path": str(path), "classification": "declared_external_toolchain", "toolchain_root": str(root), "sha256": sha256_path(path) if path.is_file() else None}
    require(False, "ASSERT_EVALUATED_CODE_CLOSURE", f"evaluated code outside anchor and toolchain roots: {path}")
    raise AssertionError("unreachable")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True); parser.add_argument("--anchor", required=True)
    parser.add_argument("--environment", required=True); parser.add_argument("--trace", action="append", default=[])
    parser.add_argument("--output"); parser.add_argument("--probe-code-path"); parser.add_argument("--probe-content-path")
    parser.add_argument("--logical-repo-path"); parser.add_argument("--stage0-precommit", action="store_true")
    args = parser.parse_args()
    repo = Path(args.repo).resolve(); scratch = repo / "software/em_charge_attribute/_scratch"
    environment = yaml.safe_load(Path(args.environment).read_text(encoding="utf-8"))
    roots = [Path(row["path"]).resolve() for row in environment["declared_read_only_toolchain_roots"]]
    try:
        require(sys.flags.isolated == 1 and sys.flags.no_user_site == 1, "ASSERT_PYTHON_STARTUP_SANITIZED", "closure guard not isolated/no-user-site")
        if args.probe_code_path:
            record = validate_code_path(
                Path(args.probe_code_path), repo, args.anchor, roots,
                Path(args.probe_content_path) if args.probe_content_path else None,
                args.logical_repo_path,
            )
            print(yaml.safe_dump({"status": "PASS", "record": record}, sort_keys=False)); return 0
        observed: set[Path] = set(); direct: set[Path] = set(); trace_records = []
        for trace_name in args.trace:
            trace = Path(trace_name).resolve(); trace_records.append({"path": str(trace), "sha256": sha256_path(trace)})
            paths, entries = extract_trace_paths(trace, repo); observed.update(paths); direct.update(entries)
        code_paths = sorted(path for path in observed if path.suffix.lower() in CODE_SUFFIXES and path.exists())
        records = []; pending = []
        for path in code_paths:
            if args.stage0_precommit and is_under(path, repo):
                rel = path.relative_to(repo).as_posix()
                if git_blob(repo, args.anchor, rel) is None and path in direct:
                    require(not is_under(path, scratch), "ASSERT_EVALUATED_CODE_CLOSURE", f"writable scratch entrypoint cannot use precommit exception: {rel}")
                    require(rel.startswith("software/em_charge_attribute/u2_") or rel == "software/em_charge_attribute/run_u2_boundary_adjudication.sh", "ASSERT_STAGE0_PRECOMMIT_SCOPE", f"non-U2 pending entrypoint {rel}")
                    pending.append(rel); continue
            records.append(validate_code_path(path, repo, args.anchor, roots))
        require(not args.stage0_precommit or bool(pending), "ASSERT_STAGE0_PRECOMMIT_SCOPE", "no pending reviewed U2 entrypoints observed")
        evidence = {
            "schema_version": "U2_EVALUATED_CODE_CLOSURE_EVIDENCE_V1",
            "status": "PASS_STAGE0_PRECOMMIT" if args.stage0_precommit else "PASS_PRODUCTION",
            "anchor": args.anchor, "trace_records": trace_records,
            "observed_path_count": len(observed), "observed_code_path_count": len(code_paths),
            "validated_records": records, "stage0_pending_orchestrator_commit_paths": sorted(pending),
            "stage0_pending_direct_entrypoint_origin_paths": sorted(str(path) for path in direct if is_under(path, repo)),
            "precommit_exception_rule": "anchor-absent direct U2 entrypoints outside writable scratch only",
            "python_isolated": sys.flags.isolated == 1, "python_no_user_site": sys.flags.no_user_site == 1,
            "writable_or_external_task_code_accepted": False,
        }
        require(bool(args.output), "ASSERT_CLOSURE_OUTPUT", "missing closure output")
        Path(args.output).write_text(yaml.safe_dump(evidence, sort_keys=False, allow_unicode=True, width=140), encoding="utf-8")
        print(f"U2_CODE_CLOSURE_PASS code_paths={len(code_paths)} precommit_pending={len(pending)}")
        return 0
    except GuardFailure as failure:
        print(f"ASSERTION_FAILED {failure.assert_id}: {failure.detail}", file=sys.stderr); return 1


if __name__ == "__main__":
    raise SystemExit(main())
