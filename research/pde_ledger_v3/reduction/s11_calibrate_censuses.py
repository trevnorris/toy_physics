#!/usr/bin/env python3
"""Build byte-shaped planted records and run both production census entrypoints."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys

from s11_census_common import read_record_map
from s11_census_math import split_delimited


REPO = Path("/var/projects/toy_physics")
REDUCTION = REPO / "research/pde_ledger_v3/reduction"
WL_RECORD = REPO / "research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out"
PY_RECORD = REPO / "research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out"
DEFAULT_SCRATCH = Path("/home/trevnorris/.s11_build/census_build")


def _line(data: dict[str, object], tag: str, payload: str | None = None) -> str:
    item = data[tag]
    value = item.payload if payload is None else payload  # type: ignore[attr-defined]
    return f"{tag}: {value}\n"


def _write_record(path: Path, lines: list[str]) -> None:
    path.write_text("".join(lines), encoding="utf-8")


def build_plants(scratch: Path) -> dict[str, Path]:
    wl = read_record_map(WL_RECORD)
    py = read_record_map(PY_RECORD)

    omitted_stem = "WL_S11_MAIN_D2_ROOT1_RANK_DROP_COEFF"
    omitted_payload = wl[omitted_stem + "_SOLUTION"].payload  # type: ignore[attr-defined]
    old_fragment = '"SOLUTION_SET" -> {{muR -> bComp}}'
    if omitted_payload.count(old_fragment) != 1:
        raise RuntimeError("omitted-branch source shape changed")
    omitted_payload = omitted_payload.replace(
        old_fragment, '"SOLUTION_SET" -> {}', 1
    )
    provenance_stem = "WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS"
    na_stem = "WL_S11_MAIN_D2_STRATUM1_ROOT_COUNT_ALL_CHANGE_LOCUS"
    omitted_path = scratch / "calibration_omitted_wl.out"
    _write_record(
        omitted_path,
        [
            _line(wl, omitted_stem + "_EQUATIONS"),
            _line(wl, omitted_stem + "_SOLUTION", omitted_payload),
            _line(wl, provenance_stem + "_EQUATIONS"),
            _line(wl, provenance_stem + "_SOLUTION"),
            _line(wl, na_stem + "_EQUATIONS"),
            _line(wl, na_stem + "_SOLUTION"),
        ],
    )

    spurious_stem = "PY_S11_MAIN_D2_ROOT1_RANK_DROP_COEFF"
    donor_stem = "PY_S11_MAIN_D2_ROOT1_KW_ZERO_LOCUS"
    target_payload = py[spurious_stem + "_SOLUTION"].payload  # type: ignore[attr-defined]
    donor_payload = py[donor_stem + "_SOLUTION"].payload  # type: ignore[attr-defined]
    target_branches = split_delimited(target_payload, "Tuple(", ")")
    donor_branches = split_delimited(donor_payload, "Tuple(", ")")
    if len(target_branches) != 1 or not donor_branches:
        raise RuntimeError("spurious-branch source shape changed")
    spurious_payload = "Tuple(" + ", ".join((target_branches[0], donor_branches[0])) + ")"
    spurious_path = scratch / "calibration_spurious_py.out"
    _write_record(
        spurious_path,
        [
            _line(py, spurious_stem + "_EQUATIONS"),
            _line(py, spurious_stem + "_SOLUTION", spurious_payload),
        ],
    )

    decided_tag = "PY_S11_MAIN_D2_ROOT1_RANK_DROP_K_IDENTICALLY_SATISFIED"
    decided_payload = py[decided_tag].payload  # type: ignore[attr-defined]
    old_token = "Str('PROVED_FALSE')"
    if decided_payload.count(old_token) != 1:
        raise RuntimeError("decidable-status source shape changed")
    decided_payload = decided_payload.replace(old_token, "Str('UNDECIDED')", 1)
    condition_tag = "PY_S11_MAIN_D2_ROOT1_RANK_DROP_K_INCONSISTENT"
    condition_stem = condition_tag[: -len("_INCONSISTENT")]
    probe_path = scratch / "calibration_decidable_py.out"
    _write_record(
        probe_path,
        [
            _line(py, decided_tag, decided_payload),
            _line(py, condition_stem + "_EQUATIONS"),
            _line(py, condition_tag),
        ],
    )
    return {
        "omitted": omitted_path,
        "spurious": spurious_path,
        "decidable": probe_path,
    }


def _run(command: list[str], log: Path) -> tuple[int, str]:
    completed = subprocess.run(
        command,
        cwd=str(REPO),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    log.write_text(completed.stdout, encoding="utf-8")
    return completed.returncode, completed.stdout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH)
    args = parser.parse_args()
    scratch = args.scratch_dir.resolve(strict=True)
    if not os.access(scratch, os.W_OK):
        print(
            "CALIBRATION_ENVIRONMENT "
            f"scratch={shlex.quote(str(scratch))} writable=False "
            "verdict=CALIBRATION_EXECUTION_BLOCKED",
            flush=True,
        )
        return 73
    plants = build_plants(scratch)
    cases = (
        (
            "planted_omitted_branch_wl",
            [sys.executable, str(REDUCTION / "s11_containment_census.py"), str(plants["omitted"])],
            "OMITTED_BRANCH",
        ),
        (
            "planted_spurious_branch_py",
            [sys.executable, str(REDUCTION / "s11_containment_census.py"), str(plants["spurious"])],
            "SPURIOUS_BRANCH",
        ),
        (
            "planted_decidable_undecided_py",
            [sys.executable, str(REDUCTION / "s11_undecided_probe_census.py"), str(plants["decidable"])],
            "DECIDED_UNDECIDED_RECORD",
        ),
    )
    misses = 0
    for case, command, expected in cases:
        print(
            "CALIBRATION_COMMAND "
            f"case={case} argv={shlex.join(command)} state=STARTED",
            flush=True,
        )
        returncode, stdout = _run(command, scratch / f"{case}.stdout")
        print(stdout, end="", flush=True)
        observed = stdout.count(f"verdict={expected}") + stdout.count(
            f"verdict={expected}_SAMPLED"
        )
        parse_failures = stdout.count("verdict=IN_POPULATION_PARSE_FAILURE")
        detected = returncode == 0 and observed > 0 and parse_failures == 0
        misses += int(not detected)
        verdict = "CALIBRATION_DETECTED" if detected else "CALIBRATION_MISS"
        print(
            "CALIBRATION_CASE "
            f"case={case} command={shlex.quote(shlex.join(command))} "
            f"returncode={returncode} expected_token={expected} "
            f"observed_count={observed} parse_failures={parse_failures} "
            f"verdict={verdict}",
            flush=True,
        )
        if not detected:
            break
    round_verdict = "CALIBRATION_PASS" if misses == 0 else "CALIBRATION_FAIL"
    print(
        f"CALIBRATION_SUMMARY cases={len(cases)} misses={misses} verdict={round_verdict}",
        flush=True,
    )
    return 0 if misses == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
