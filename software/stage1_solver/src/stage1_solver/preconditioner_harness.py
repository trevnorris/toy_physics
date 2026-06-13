"""Command-line harness for Step 3b coupled-branch preconditioning."""

from __future__ import annotations

import sys

from .config import HarnessConfig
from .coupled_branch import run_step3b_preconditioner, write_step3b_report


def main() -> int:
    config = HarnessConfig(
        run_root="software/stage1_solver/runs/step3b_preconditioner",
        report_path="software/stage1_solver/reports/step3b_preconditioner.md",
    )
    results = run_step3b_preconditioner(config)
    report_path = write_step3b_report(results, config.report_path)
    print(f"wrote preconditioner report: {report_path}")
    if not results["passed"]:
        print("preconditioner engineering gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("preconditioner engineering gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
