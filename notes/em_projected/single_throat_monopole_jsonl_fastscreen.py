#!/usr/bin/env python3
"""Partial falsification screen for single_throat_monopole.py JSONL diagnostics."""
from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any


DEFAULT_THRESHOLDS = {
    "dP_slope_target": -1.0,
    "dP_slope_tol": 0.35,
    "geff_slope_target": -2.0,
    "geff_slope_tol": 0.35,
    "mach_max": 0.6,
    "min_fit_points": 8,
}
FLOAT_EDGE_EPS = 1.0e-12


def load_events(path: pathlib.Path) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    for line_no, line in enumerate(path.read_text().splitlines(), start=1):
        line = line.strip()
        if not line:
            continue
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError as exc:
            events.append({"event": "parse_error", "line": line_no, "message": str(exc)})
    return events


def classify_events(events: list[dict[str, Any]], thresholds: dict[str, float] | None = None) -> dict[str, Any]:
    limits = dict(DEFAULT_THRESHOLDS)
    if thresholds:
        limits.update(thresholds)

    diag_events = [event for event in events if event.get("event") == "diag"]
    init_event = next((event for event in events if event.get("event") == "init"), None)
    done_event = next((event for event in reversed(events) if event.get("event") == "done"), None)
    parse_errors = [event for event in events if event.get("event") == "parse_error"]

    if not diag_events:
        return {
            "status": "INCOMPLETE",
            "failures": [],
            "warnings": ["no diagnostic events found"]
            + [f"malformed JSONL line {event['line']}: {event['message']}" for event in parse_errors],
            "latest_diag": None,
            "init": init_event,
            "done": done_event,
            "thresholds": limits,
        }

    latest = diag_events[-1]
    failures: list[str] = []
    warnings: list[str] = [
        f"malformed JSONL line {event['line']}: {event['message']}" for event in parse_errors
    ]

    fits = latest.get("fits", {})
    mach_max = float(latest.get("mach_max", float("nan")))
    dP_slope = float(fits.get("dP_slope", float("nan")))
    dP_npts = int(fits.get("dP_npts", 0))
    geff_slope = float(fits.get("geff_slope", float("nan")))
    geff_npts = int(fits.get("geff_npts", 0))

    if dP_npts < limits["min_fit_points"]:
        warnings.append(f"insufficient dP fit points: {dP_npts}")
    elif abs(dP_slope - limits["dP_slope_target"]) > limits["dP_slope_tol"] + FLOAT_EDGE_EPS:
        failures.append(f"dP slope misses -1 target: {dP_slope:.6g}")

    if geff_npts < limits["min_fit_points"]:
        warnings.append(f"insufficient geff fit points: {geff_npts}")
    elif abs(geff_slope - limits["geff_slope_target"]) > limits["geff_slope_tol"] + FLOAT_EDGE_EPS:
        failures.append(f"g_eff slope misses -2 target: {geff_slope:.6g}")

    if mach_max > limits["mach_max"] + FLOAT_EDGE_EPS:
        failures.append(f"flow is too compressible for weak-field screen: mach_max={mach_max:.6g}")

    if failures:
        status = "FAIL"
    elif warnings:
        status = "INCOMPLETE"
    else:
        status = "PASS"

    return {
        "status": status,
        "failures": failures,
        "warnings": warnings,
        "latest_diag": latest,
        "init": init_event,
        "done": done_event,
        "thresholds": limits,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("jsonl_log", help="captured stdout log from single_throat_monopole.py")
    parser.add_argument("--output-json", help="write verdict JSON to this path")
    args = parser.parse_args()

    verdict = classify_events(load_events(pathlib.Path(args.jsonl_log)))
    print("SINGLE-THROAT MONOPOLE FASTSCREEN")
    print("  status      =", verdict["status"])
    latest = verdict["latest_diag"]
    if latest is not None:
        fits = latest.get("fits", {})
        print("  mach_max    =", latest.get("mach_max"))
        print("  dP_slope    =", fits.get("dP_slope"))
        print("  geff_slope  =", fits.get("geff_slope"))
        print("  dP_npts     =", fits.get("dP_npts"))
        print("  geff_npts   =", fits.get("geff_npts"))
    if verdict["failures"]:
        print("  failures:")
        for item in verdict["failures"]:
            print("   -", item)
    if verdict["warnings"]:
        print("  warnings:")
        for item in verdict["warnings"]:
            print("   -", item)
    if args.output_json:
        pathlib.Path(args.output_json).write_text(json.dumps(verdict, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
