"""Command-line harness for Step 8b driven P2 tangent and absorber."""

from __future__ import annotations

import sys

from .p2_driven_absorber import run_step8b, step8b_default_config, write_step8b_report


def main() -> int:
    config = step8b_default_config()
    results = run_step8b(config)
    report_path = write_step8b_report(results, config.report_path)
    print(f"wrote Step 8b report: {report_path}")
    print(f"wrote machine-readable table: {results['machine_readable_table']}")
    print(f"diagnostics digest: {results['diagnostics_digest']}")
    if not results["passed"]:
        print("Step 8b engineering gate failed; see report diagnostics", file=sys.stderr)
        return 1
    print("Step 8b engineering gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
