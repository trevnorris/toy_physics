"""Command-line validation harness for Stage-1 solver validation."""

from __future__ import annotations

import sys

from .benchmarks import run_all_benchmarks
from .config import HarnessConfig
from .mms_benchmarks import run_all_mms_benchmarks
from .report import write_report


def main() -> int:
    config = HarnessConfig()
    step1 = run_all_benchmarks(config)
    mms = run_all_mms_benchmarks(config)
    results = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "step1": step1,
        "mms": mms,
        "passed": step1["passed"] and mms["passed"],
    }
    report_path = write_report(results, config.report_path)
    print(f"wrote validation report: {report_path}")
    if not results["passed"]:
        print("validation gate failed", file=sys.stderr)
        return 1
    print("validation gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
