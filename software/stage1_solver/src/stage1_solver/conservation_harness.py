"""Command-line harness for Step 6 conservation diagnostics."""

from __future__ import annotations

import sys

from .config import HarnessConfig
from .conservation_diagnostics import run_step6, write_step6_report


def main() -> int:
    config = HarnessConfig(
        run_root="software/stage1_solver/runs/step6_conservation_diagnostics",
        report_path="software/stage1_solver/reports/step6_conservation_diagnostics.md",
    )
    results = run_step6(config)
    report_path = write_step6_report(results, config.report_path)
    print(f"wrote conservation report: {report_path}")
    print(f"wrote machine-readable table: {results['machine_readable_table']}")
    print(f"diagnostics digest: {results['diagnostics_digest']}")
    if not results["passed"]:
        print("conservation diagnostics gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("conservation diagnostics gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
