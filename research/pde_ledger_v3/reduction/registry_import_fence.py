#!/usr/bin/env python3
"""Discover registry-importing engines and optionally verify committed output."""

from __future__ import annotations

import argparse
import ast
import difflib
import hashlib
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


HERE = Path(__file__).resolve().parent
LEDGER = HERE.parent
WL_REGISTRY_IMPORT = re.compile(r"(?:Get|Needs)\s*\[.*registry_read|registry_read`")


def _python_imports_registry(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import) and any(
            alias.name == "registry_read" for alias in node.names
        ):
            return True
        if isinstance(node, ast.ImportFrom) and node.module == "registry_read":
            return True
    return False


def _wl_imports_registry(path: Path) -> bool:
    return WL_REGISTRY_IMPORT.search(path.read_text(encoding="utf-8")) is not None


def discover_registry_importing_engines(ledger: Path = LEDGER) -> tuple[Path, ...]:
    """Compute the consumer fence from import syntax, never from a stage list."""
    candidates = tuple(sorted((ledger / "scripts").glob("*_audit.py"))) + tuple(
        sorted((ledger / "mathematica").glob("*_audit.wl"))
    )
    selected = []
    for path in candidates:
        imports = (
            _python_imports_registry(path)
            if path.suffix == ".py"
            else _wl_imports_registry(path)
        )
        if imports:
            selected.append(path.resolve())
    return tuple(selected)


def committed_output_for(engine: Path) -> Path:
    return engine.parent / "out" / f"{engine.stem}.out"


def _normalise_runtime(engine: Path, payload: bytes) -> bytes:
    if engine.name != "S10_brane_mode_spectrum_mathematica_audit.wl":
        return payload
    return b"\n".join(
        line
        for line in payload.split(b"\n")
        if not line.startswith(b"WL_S10_RUNTIME_SECONDS:")
    )


def _digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class Verification:
    engine: Path
    returncode: int
    stdout: bytes
    stderr: bytes
    committed: bytes
    timed_out: bool = False

    @property
    def output_matches(self) -> bool:
        return _normalise_runtime(self.engine, self.stdout) == _normalise_runtime(
            self.engine, self.committed
        )

    @property
    def passed(self) -> bool:
        return not self.timed_out and self.returncode == 0 and self.output_matches


def verify(engine: Path) -> Verification:
    command = (
        [sys.executable, str(engine)]
        if engine.suffix == ".py"
        else ["timeout", "900", "wolframscript", "-file", str(engine)]
    )
    committed = committed_output_for(engine).read_bytes()
    timeout_seconds = 910 if engine.suffix == ".wl" else 3600
    try:
        completed = subprocess.run(
            command,
            cwd=engine.parent,
            check=False,
            capture_output=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        return Verification(
            engine,
            124,
            error.stdout or b"",
            error.stderr or b"",
            committed,
            True,
        )
    return Verification(
        engine,
        completed.returncode,
        completed.stdout,
        completed.stderr,
        committed,
    )


def _print_verification(result: Verification) -> None:
    print(f"REGISTRY_IMPORT_ENGINE: {result.engine}")
    print(f"ENGINE_EXIT_OPERAND observed: {result.returncode}")
    print(f"ENGINE_EXIT_RESIDUAL nonzero: {result.returncode != 0}")
    print(f"ENGINE_TIMEOUT_RESIDUAL: {result.timed_out}")
    print(f"ENGINE_OUTPUT_OPERAND regenerated-sha256: {_digest(result.stdout)}")
    print(f"ENGINE_OUTPUT_OPERAND committed-sha256: {_digest(result.committed)}")
    print(f"ENGINE_OUTPUT_RESIDUAL byte-identical: {result.output_matches}")
    stderr = result.stderr.decode("utf-8", errors="replace")
    print(f"ENGINE_STDERR_BEGIN {result.engine.name}")
    print(stderr, end="" if not stderr or stderr.endswith("\n") else "\n")
    print(f"ENGINE_STDERR_END {result.engine.name}")
    if not result.output_matches:
        generated = _normalise_runtime(result.engine, result.stdout).decode(
            "utf-8", errors="replace"
        )
        committed = _normalise_runtime(result.engine, result.committed).decode(
            "utf-8", errors="replace"
        )
        diff = difflib.unified_diff(
            committed.splitlines(),
            generated.splitlines(),
            fromfile=str(committed_output_for(result.engine)),
            tofile=f"regenerated:{result.engine}",
            lineterm="",
        )
        print("ENGINE_OUTPUT_DIFF_BEGIN")
        for line in diff:
            print(line)
        print("ENGINE_OUTPUT_DIFF_END")
    print(f"REGISTRY_IMPORT_ENGINE_GUARD: {'PASS' if result.passed else 'FAIL'}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="print only computed paths")
    parser.add_argument("--verify", action="store_true", help="run sequentially and diff stdout")
    arguments = parser.parse_args(argv)
    engines = discover_registry_importing_engines()
    if arguments.list:
        for engine in engines:
            print(engine)
        return 0
    print(f"REGISTRY_IMPORT_FENCE_OPERAND discovered-count: {len(engines)}")
    for engine in engines:
        print(f"REGISTRY_IMPORT_FENCE_MEMBER: {engine}")
    print(f"REGISTRY_IMPORT_FENCE_RESIDUAL missing-committed-output: {sum(not committed_output_for(engine).is_file() for engine in engines)}")
    if not arguments.verify:
        passed = bool(engines) and all(
            committed_output_for(engine).is_file() for engine in engines
        )
        print(f"REGISTRY_IMPORT_FENCE: {'PASS' if passed else 'FAIL'}")
        return 0 if passed else 1
    results = []
    for engine in engines:
        result = verify(engine)
        _print_verification(result)
        results.append(result)
    passed = bool(results) and all(result.passed for result in results)
    print(f"REGISTRY_IMPORT_FENCE_VERIFY: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
