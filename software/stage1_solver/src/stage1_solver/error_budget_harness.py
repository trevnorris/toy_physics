"""Command-line harness for Step 7 numerical error-budget composition."""

from __future__ import annotations

import sys

from .config import HarnessConfig
from .error_budget import run_step7, write_step7_report


def main() -> int:
    config = HarnessConfig(
        run_root="software/stage1_solver/runs/step7_error_budget",
        report_path="software/stage1_solver/reports/step7_error_budget.md",
    )
    results = run_step7(config)
    report_path = write_step7_report(results, config.report_path)
    print(f"wrote error-budget report: {report_path}")
    print(f"wrote machine-readable table: {results['machine_readable_table']}")
    print(f"diagnostics digest: {results['diagnostics_digest']}")
    if not results["passed"]:
        print("error-budget gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("error-budget gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
