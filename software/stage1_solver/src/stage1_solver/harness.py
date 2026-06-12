"""Command-line validation harness for Build Directive step 1."""

from __future__ import annotations

import sys

from .benchmarks import run_all_benchmarks
from .config import HarnessConfig
from .report import write_report


def main() -> int:
    config = HarnessConfig()
    results = run_all_benchmarks(config)
    report_path = write_report(results, config.report_path)
    print(f"wrote validation report: {report_path}")
    if not results["passed"]:
        print("validation gate failed", file=sys.stderr)
        return 1
    print("validation gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
