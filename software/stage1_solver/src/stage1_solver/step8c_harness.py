"""Command-line harness for Step 8c conservation and surrogate response."""

from __future__ import annotations

import sys

from .p2_conservation_response import run_step8c, step8c_default_config, write_step8c_report


def main() -> int:
    config = step8c_default_config()
    results = run_step8c(config)
    report_path = write_step8c_report(results, config.report_path)
    print(f"wrote Step 8c report: {report_path}")
    print(f"wrote machine-readable table: {results['machine_readable_table']}")
    print(f"diagnostics digest: {results['diagnostics_digest']}")
    if not results["passed"]:
        print("Step 8c engineering gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("Step 8c engineering gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
