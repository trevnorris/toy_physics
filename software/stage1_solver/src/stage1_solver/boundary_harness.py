"""Command-line harness for Step 5 boundary characterization."""

from __future__ import annotations

import sys

from .boundary_characterization import run_step5, write_step5_report
from .config import HarnessConfig


def main() -> int:
    config = HarnessConfig(
        run_root="software/stage1_solver/runs/step5_boundary_characterization",
        report_path="software/stage1_solver/reports/step5_boundary_characterization.md",
    )
    results = run_step5(config)
    report_path = write_step5_report(results, config.report_path)
    print(f"wrote boundary-characterization report: {report_path}")
    print(f"wrote machine-readable table: {results['machine_readable_table']}")
    if not results["passed"]:
        print("boundary-characterization gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("boundary-characterization gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
