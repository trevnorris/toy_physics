"""Command-line harness for Step 8a P2 tangent operator."""

from __future__ import annotations

import sys

from .config import HarnessConfig
from .p2_tangent import run_step8a, write_step8a_report


def main() -> int:
    config = HarnessConfig(
        run_root="software/stage1_solver/runs/step8a_p2_tangent",
        report_path="software/stage1_solver/reports/step8a_p2_tangent.md",
    )
    results = run_step8a(config)
    report_path = write_step8a_report(results, config.report_path)
    print(f"wrote P2 tangent report: {report_path}")
    print(f"wrote machine-readable table: {results['machine_readable_table']}")
    print(f"diagnostics digest: {results['diagnostics_digest']}")
    if not results["passed"]:
        print("P2 tangent engineering gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("P2 tangent engineering gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
