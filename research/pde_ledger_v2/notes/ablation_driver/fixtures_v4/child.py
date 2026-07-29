#!/usr/bin/env python3
"""Deterministic fixture child; this is a producer/checker, never a driver."""

from __future__ import annotations

import hashlib
import os
import pathlib
import sys
import time


ARGV_PROBE = "fixture; exit 97"
EXPECTED_ENVIRONMENT = {
    "LC_ALL": "C",
    "PATH": (
        str(pathlib.Path(__file__).resolve().parent)
        + ":/usr/local/bin:/usr/bin:/bin"
    ),
    "PYTHONDONTWRITEBYTECODE": "1",
    "TZ": "UTC",
}


def read(path: str) -> str:
    return pathlib.Path(path).read_text(encoding="utf-8")


def write(path: str, value: str) -> None:
    pathlib.Path(path).write_text(value, encoding="utf-8", newline="\n")


def selected_mutant(target: str, prefix: str) -> str:
    matches = [word for word in read(target).split() if prefix in word]
    if len(matches) != 1:
        raise RuntimeError(f"one mutated token required, found {matches!r}")
    return matches[0].split("=", 1)[-1]


def producer(case: str, target: str, artifact: str, reset: str) -> int:
    observed = hashlib.sha256(pathlib.Path(target).read_bytes()).hexdigest()
    print(f"FIXTURE_TARGET_SHA256={observed}")
    if case == "a1":
        text = read(target)
        stale = pathlib.Path(reset) / "fixture.cpython.pyc"
        if "A1_MUTANT_EMIT" in text:
            write(artifact, "fresh-from-A1_MUTANT_EMIT\n")
        elif stale.exists():
            print("producer-a1-stale-bytecode")
            return 29
        stale.parent.mkdir(parents=True, exist_ok=True)
        stale.write_bytes(b"fixture stale bytecode sentinel\n")
        print("producer-a1")
        return 0
    if case == "a2":
        write(artifact, "producer-was-invoked\n")
        time.sleep(1.2)
        print("producer-a2")
        return 0
    if case == "a3":
        token = selected_mutant(target, "A3_MUTANT_")
        time.sleep(1.2)
        print(f"producer:{token}")
        return 0
    if case == "a4":
        token = selected_mutant(target, "A4_MUTANT")
        time.sleep(1.2)
        write(artifact, f"produced:{token}\n")
        print(f"producer:{token}")
        return 0
    if case == "a5":
        token = selected_mutant(target, "A5_MUTANT_")
        payload = f"artifact:{token}\n"
        write(artifact, payload)
        print(f"producer.stdout:{token}")
        print(f"producer.stderr:{token}", file=sys.stderr)
        return 0
    if case == "a8":
        token = selected_mutant(target, "A8_MUTANT")
        write(
            artifact,
            "DIM|axes=L,M,T|name=scope.record|"
            f"exponents={{{len(token)}, 1, 0}}\n",
        )
        print(f"producer:{token}")
        return 0
    if case == "a9":
        token = selected_mutant(target, "A9_MUTANT_")
        print(f"producer:{token}")
        return 17 if token.endswith("_Q") else 0
    raise RuntimeError(f"unknown producer case {case!r}")


def checker(case: str, target: str, artifact: str, _reset: str) -> int:
    if case == "a1":
        if not pathlib.Path(artifact).exists():
            print("RESULT|status=ABSENT")
            return 0
        value = read(artifact)
        if "A1_MUTANT_EMIT" in value and "A1_MUTANT_EMIT" in read(target):
            print("RESULT|status=FRESH")
            return 0
        print("RESULT|status=RESIDUAL")
        return 19
    if case == "a2":
        print("checker-a2")
        return 0
    if case in {"a3", "a4"}:
        print(f"checker:{selected_mutant(target, case.upper() + '_MUTANT')}")
        return 0
    if case == "a5":
        token = selected_mutant(target, "A5_MUTANT_")
        value = read(artifact)
        if token not in value:
            print(f"checker.stderr:artifact-mismatch:{token}", file=sys.stderr)
            return 21
        print(f"checker.stdout:{token}")
        print(f"checker.stderr:{token}", file=sys.stderr)
        return 0
    if case == "a8":
        token = selected_mutant(target, "A8_MUTANT")
        print(f"RESULT|status=OBSERVED|token={token}")
        return 0
    if case == "a9":
        token = selected_mutant(target, "A9_MUTANT_")
        print(f"checker:{token}")
        return 23 if token.endswith("_R") else 0
    raise RuntimeError(f"unknown checker case {case!r}")


def main() -> int:
    if len(sys.argv) != 8:
        print(
            "usage: child.py ROLE CASE TARGET ARTIFACT RESET POISON_NAME ARGV_PROBE",
            file=sys.stderr,
        )
        return 64
    role, case, target, artifact, reset, poison_name, argv_probe = sys.argv[1:]
    if argv_probe != ARGV_PROBE:
        print("fixture argv was reinterpreted", file=sys.stderr)
        return 64
    if poison_name in os.environ or dict(os.environ) != EXPECTED_ENVIRONMENT:
        print("fixture child environment is not the configured environment", file=sys.stderr)
        return 64
    if role == "producer":
        return producer(case, target, artifact, reset)
    if role == "checker":
        return checker(case, target, artifact, reset)
    print(f"unknown role {role!r}", file=sys.stderr)
    return 64


if __name__ == "__main__":
    raise SystemExit(main())
