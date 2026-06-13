"""Command-line harness for Step 4 coupled-branch grid convergence."""

from __future__ import annotations

import sys

from .config import HarnessConfig
from .convergence import run_step4, write_step4_report


def main() -> int:
    config = HarnessConfig(
        run_root="software/stage1_solver/runs/step4_convergence_study",
        report_path="software/stage1_solver/reports/step4_convergence_study.md",
    )
    results = run_step4(config)
    report_path = write_step4_report(results, config.report_path)
    print(f"wrote convergence report: {report_path}")
    print(f"wrote machine-readable table: {results['machine_readable_table']}")
    if not results["passed"]:
        print("convergence engineering gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("convergence engineering gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
