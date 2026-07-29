#!/usr/bin/env python3
"""Verify the committed external authority against the complete suite inventory."""

from __future__ import annotations

import hashlib
import pathlib
import subprocess


SUITE = pathlib.Path(__file__).resolve().parent
PROJECT_ROOT = SUITE.parents[4]
AUTHORITY = SUITE.parent / "fixtures_v4.accepted.sha256"
GOVERNING = (
    SUITE.parent / "CONTRACT.md",
    SUITE.parent / "REQUIREMENTS.md",
    PROJECT_ROOT
    / "research/pde_ledger_v2/scripts/check_ledger_dimensions_pin.py",
    PROJECT_ROOT
    / "research/pde_ledger_v2/scripts/compare_dimension_artifacts.py",
    PROJECT_ROOT
    / "research/pde_ledger_v2/scripts/ledger_dimensions.accepted.sha256",
    PROJECT_ROOT
    / "research/pde_ledger_v2/scripts/ledger_dimensions.py",
    PROJECT_ROOT
    / "research/pde_ledger_v2/scripts/"
    "ledger_stage023_nullspace_underdetermination_sympy_audit.py",
    PROJECT_ROOT
    / "research/pde_ledger_v2/scripts/"
    "ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt",
    PROJECT_ROOT
    / "research/pde_ledger_v2/mathematica/out/"
    "ledger_stage023_nullspace_underdetermination_mathematica_audit.out",
    PROJECT_ROOT / "research/pde_ledger_v2/notes/stage023_step_h_evidence/include_list.tsv",
    PROJECT_ROOT / "research/pde_ledger_v2/notes/stage023_step_h_evidence/results.tsv",
)


def relative(path: pathlib.Path) -> str:
    return path.relative_to(PROJECT_ROOT).as_posix()


def git(*arguments: str) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *arguments],
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=30,
    )


def inventory() -> list[pathlib.Path]:
    files = [
        path
        for path in SUITE.rglob("*")
        if (path.is_file() or path.is_symlink())
    ]
    files.extend(GOVERNING)
    return sorted(files, key=lambda path: relative(path).encode())


def suite_bytecode() -> list[pathlib.Path]:
    return sorted(
        (
            path
            for path in SUITE.rglob("*")
            if path.name == "__pycache__" or path.suffix == ".pyc"
        ),
        key=lambda path: relative(path).encode(),
    )


def fail(message: str) -> int:
    print(message)
    return 1


def main() -> int:
    bytecode = suite_bytecode()
    if bytecode:
        return fail(f"suite-local Python bytecode is forbidden: {relative(bytecode[0])}")
    if not AUTHORITY.is_file():
        return fail(f"PENDING-COMMIT: missing authority {relative(AUTHORITY)}")
    committed = git("show", f"HEAD:{relative(AUTHORITY)}")
    if committed.returncode != 0:
        return fail(f"PENDING-COMMIT: authority is not committed at {relative(AUTHORITY)}")
    authority_bytes = AUTHORITY.read_bytes()
    if authority_bytes != committed.stdout:
        return fail(f"freeze authority differs from HEAD: {relative(AUTHORITY)}")
    dirty = git(
        "status",
        "--porcelain",
        "--untracked-files=all",
        "--",
        relative(SUITE),
        relative(AUTHORITY),
    )
    if dirty.returncode != 0:
        return fail("cannot inspect freeze worktree state")
    if dirty.stdout:
        return fail("freeze worktree differs from HEAD")

    expected: dict[str, str] = {}
    try:
        lines = authority_bytes.decode("ascii").splitlines()
    except UnicodeError:
        return fail("freeze authority is not ASCII")
    for line in lines:
        digest, separator, name = line.partition("  ")
        if (
            not separator
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
            or not name
            or name in expected
        ):
            return fail(f"malformed freeze authority line: {line!r}")
        expected[name] = digest

    actual = inventory()
    actual_names = [relative(path) for path in actual]
    if set(actual_names) != set(expected):
        return fail("freeze inventory mismatch")
    for path in actual:
        if path.is_symlink():
            return fail(f"freeze inventory contains symlink: {relative(path)}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected[relative(path)]:
            return fail(f"freeze digest mismatch: {relative(path)}")
    print(f"freeze verified: {len(actual)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
