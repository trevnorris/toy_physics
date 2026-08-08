#!/usr/bin/env python3
"""Prove that partial unittest discovery is reported as the wrong runner."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> int:
    completed = subprocess.run(
        [sys.executable, "-m", "unittest", "discover", "-v"],
        cwd=HERE,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    print("WRONG_RUNNER_STDOUT_BEGIN")
    print(completed.stdout, end="")
    print("WRONG_RUNNER_STDOUT_END")
    print("WRONG_RUNNER_STDERR_BEGIN")
    print(completed.stderr, end="")
    print("WRONG_RUNNER_STDERR_END")
    expected_exit_class = "nonzero"
    observed_exit_class = "nonzero" if completed.returncode != 0 else "zero"
    message_present = "WRONG_RUNNER:" in completed.stderr
    print(f"RUNNER_OPERAND expected-exit-class: {expected_exit_class}")
    print(f"RUNNER_OPERAND observed-exit-class: {observed_exit_class}")
    print(
        "RUNNER_RESIDUAL: "
        f"exit-class-equal={expected_exit_class == observed_exit_class} "
        f"message-present={message_present}"
    )
    passed = observed_exit_class == expected_exit_class and message_present
    print(f"W3_RUNNER_CHECK: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
