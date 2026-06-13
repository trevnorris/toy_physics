"""Command-line harness for Step 3 coupled branch engineering smoke."""

from __future__ import annotations

import sys

from .config import HarnessConfig
from .coupled_branch import run_step3, write_step3_report


def main() -> int:
    config = HarnessConfig(
        run_root="software/stage1_solver/runs/step3_coupled_branch_smoke",
        report_path="software/stage1_solver/reports/step3_coupled_branch_smoke.md",
    )
    results = run_step3(config)
    report_path = write_step3_report(results, config.report_path)
    print(f"wrote engineering-smoke report: {report_path}")
    if not results["passed"]:
        print("engineering-smoke gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("engineering-smoke gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

