#!/usr/bin/env python3
"""Independent able-to-fail check for the S10 engine-output harness.

Written BEFORE the harness repair was implemented, and deliberately kept out of
the repository until afterwards, so that the thing being tested could not be
built to satisfy it.

Each case corrupts ONE payload in a COPY of an engine output (or one key in a
copy of the config) and asserts that ONE named counter in the harness report
moves in the stated direction.  A counter that does not move is a layer that is
not measuring what its name says.

The point is not that the harness reports something.  It is that the harness
reports something DIFFERENT when the input is wrong.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

LEDGER = Path("/var/projects/toy_physics/research/pde_ledger_v3")
HARNESS = LEDGER / "reduction" / "engine_output_checks.py"
CONFIG = LEDGER / "reduction" / "checks_S10.yaml"
WL_OUT = LEDGER / "mathematica" / "out" / "S10_brane_mode_spectrum_mathematica_audit.out"
PY_OUT = LEDGER / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out"

SUMMARY = re.compile(r"^(?P<layer>[A-Z_]+): (?P<body>.*)$")
FIELD = re.compile(r"(?P<key>[a-z_]+)=(?P<value>-?\d+)")


@dataclass(frozen=True)
class Report:
    exit_code: int
    counters: dict[str, int]
    stderr: str

    def get(self, name: str) -> int | None:
        return self.counters.get(name)


def run_harness(config: Path, wl: Path, py: Path) -> Report:
    proc = subprocess.run(
        [
            sys.executable,
            str(HARNESS),
            "--config",
            str(config),
            "--output",
            f"wl={wl}",
            "--output",
            f"py={py}",
        ],
        cwd=str(LEDGER),
        capture_output=True,
        text=True,
        timeout=900,
    )
    counters: dict[str, int] = {}
    for line in proc.stdout.splitlines():
        match = SUMMARY.match(line)
        if not match:
            continue
        layer = match.group("layer").lower()
        for field in FIELD.finditer(match.group("body")):
            counters[f"{layer}.{field.group('key')}"] = int(field.group("value"))
    return Report(proc.returncode, counters, proc.stderr[:400])


def find_tag(path: Path, pattern: str) -> tuple[int, str, str]:
    """First line whose tag matches `pattern`; returns index, tag, payload."""
    compiled = re.compile(pattern)
    for index, line in enumerate(path.read_text().splitlines()):
        if ":" not in line:
            continue
        tag, _, payload = line.partition(":")
        if compiled.search(tag):
            return index, tag, payload.strip()
    raise SystemExit(f"no tag matching {pattern!r} in {path}")


def rewrite_line(source: Path, target: Path, index: int, new_line: str) -> None:
    lines = source.read_text().splitlines()
    lines[index] = new_line
    target.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------
# The cases.  Each returns (label, expectation, corrupt_fn).
# corrupt_fn receives a scratch dir and returns (config, wl, py) to run.
# --------------------------------------------------------------------------


def case_cross_engine(scratch: Path):
    """A payload that is half of a configured cross-engine pair is altered."""
    index, tag, payload = find_tag(WL_OUT, r"^WL_S10_MAIN_D3_Q3_DETERMINANT$")
    wl = scratch / "wl.out"
    rewrite_line(WL_OUT, wl, index, f"{tag}: {payload} + 1")
    return CONFIG, wl, PY_OUT


def _rewrite_control_copies(scratch: Path, suffix: str, make_payload) -> Path:
    """Rewrite every control package's copy of one suffix.

    `make_payload(package, main_payload)` returns the new payload.
    """
    lines = WL_OUT.read_text().splitlines()
    main_tag = f"WL_S10_MAIN_{suffix}"
    main_payload = None
    for line in lines:
        tag, _, payload = line.partition(":")
        if tag == main_tag:
            main_payload = payload.strip()
            break
    if main_payload is None:
        raise SystemExit(f"no main tag {main_tag}")
    out = []
    touched = 0
    for line in lines:
        tag, _, _ = line.partition(":")
        if tag.startswith("WL_S10_X") and tag.endswith(f"_{suffix}"):
            package = tag[len("WL_S10_") : -len(f"_{suffix}") - 1]
            out.append(f"{tag}: {make_payload(package, main_payload)}")
            touched += 1
        else:
            out.append(line)
    if touched == 0:
        raise SystemExit(f"no control copies of {suffix}")
    path = scratch / "wl.out"
    path.write_text("\n".join(out) + "\n")
    return path


def case_control_response(_scratch: Path):
    """Differential test: force one suffix INVARIANT, then force it RESPONSIVE.

    Comparing two corrupted runs against EACH OTHER, rather than against the
    baseline, makes the result independent of what the tag happened to be.
    Under the same corruption site, `invariant` must be strictly higher in the
    all-agree run and `responsive` strictly higher in the all-differ run.
    """
    suffix = "D3_Q3_DETERMINANT"
    with tempfile.TemporaryDirectory() as raw_a, tempfile.TemporaryDirectory() as raw_b:
        agree = _rewrite_control_copies(
            Path(raw_a), suffix, lambda _package, main: main
        )
        report_agree = run_harness(CONFIG, agree, PY_OUT)
        differ = _rewrite_control_copies(
            Path(raw_b),
            suffix,
            lambda package, main: f"{main} + {abs(hash(package)) % 97 + 1}",
        )
        report_differ = run_harness(CONFIG, differ, PY_OUT)

    inv_a = report_agree.get("control_response.invariant")
    inv_b = report_differ.get("control_response.invariant")
    res_a = report_agree.get("control_response.responsive")
    res_b = report_differ.get("control_response.responsive")
    if None in (inv_a, inv_b, res_a, res_b):
        return False, f"counters absent (invariant={inv_a}/{inv_b}, responsive={res_a}/{res_b})"
    moved = inv_a > inv_b and res_b > res_a
    return moved, (
        f"all-agree invariant={inv_a} responsive={res_a}; "
        f"all-differ invariant={inv_b} responsive={res_b}"
    )


def case_tag_parity(scratch: Path):
    """A tag present in main is deleted from a control package."""
    index, _, _ = find_tag(WL_OUT, r"^WL_S10_XFORM_DIVONLY_D3_Q3_DETERMINANT$")
    lines = WL_OUT.read_text().splitlines()
    del lines[index]
    wl = scratch / "wl.out"
    wl.write_text("\n".join(lines) + "\n")
    return CONFIG, wl, PY_OUT


def case_dimensions(scratch: Path):
    """A dimensionful sum is made inhomogeneous by adding a bare number."""
    index, tag, payload = find_tag(WL_OUT, r"^WL_S10_MAIN_D3_Q3_DETERMINANT$")
    wl = scratch / "wl.out"
    rewrite_line(WL_OUT, wl, index, f"{tag}: {payload} + 1")
    return CONFIG, wl, PY_OUT


def case_declared_package_absent(scratch: Path):
    """The config declares a control package that the output does not contain.

    This is the H0b guard: a configured layer that compares nothing must fail,
    not report a quiet zero.
    """
    text = CONFIG.read_text()
    if "control_packages:" not in text:
        return None  # guard not implemented; the case reports that
    text = text.replace(
        "control_packages: [", "control_packages: [XFORM_NOT_A_REAL_PACKAGE, ", 1
    )
    config = scratch / "checks.yaml"
    config.write_text(text)
    return config, WL_OUT, PY_OUT


def case_coefficient_dimension_split(scratch: Path):
    """One entry of a per-coefficient dimension list is changed.

    The config maps the family to a single symbol, so the entries must agree.
    A disagreement must be reported, not silently resolved to the first entry.
    """
    index, tag, payload = find_tag(
        WL_OUT, r"^WL_S10_MAIN_D3_Q6_INERTIAL_COEFFICIENT_DIMENSIONS$"
    )
    # Replace the first vector's leading component with something else.
    corrupted = payload.replace("{-braneDimension, 0, 1}", "{-braneDimension, 7, 1}", 1)
    if corrupted == payload:
        corrupted = payload.replace("{", "{9 + ", 1)
    wl = scratch / "wl.out"
    rewrite_line(WL_OUT, wl, index, f"{tag}: {corrupted}")
    return CONFIG, wl, PY_OUT


CASES = [
    ("cross-engine", "cross_engine.disagree", "rises", case_cross_engine),
    ("control-response", "SELF", "differential", case_control_response),
    ("tag-parity", "tag_parity.gaps", "rises", case_tag_parity),
    ("dimensions", "dimensions.non_homogeneous", "rises", case_dimensions),
    ("declared-package-absent", "EXIT", "nonzero", case_declared_package_absent),
    ("coefficient-dimension-split", "EXIT", "nonzero", case_coefficient_dimension_split),
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", help="run one case by label")
    args = parser.parse_args()

    print("BASELINE")
    baseline = run_harness(CONFIG, WL_OUT, PY_OUT)
    for key in sorted(baseline.counters):
        print(f"  {key} = {baseline.counters[key]}")
    print(f"  EXIT = {baseline.exit_code}")
    if baseline.stderr:
        print(f"  stderr: {baseline.stderr.splitlines()[0][:160]}")
    print()

    failures = 0
    for label, counter, direction, corrupt in CASES:
        if args.only and args.only != label:
            continue
        with tempfile.TemporaryDirectory() as raw:
            scratch = Path(raw)
            prepared = corrupt(scratch)
            if prepared is None:
                print(f"SKIP  {label}: precondition absent (feature not implemented)")
                failures += 1
                continue
            if counter == "SELF":
                moved, detail = prepared
                status = "ABLE-TO-FAIL" if moved else "⛔ DEAD"
                if not moved:
                    failures += 1
                print(f"{status:14} {label:28} {detail}")
                continue
            config, wl, py = prepared
            report = run_harness(config, wl, py)

        if counter == "EXIT":
            # An exit code alone is not evidence when the baseline already
            # exits nonzero for unrelated reasons.  Require a NEW complaint:
            # some text in the corrupted run's stderr that the clean run did
            # not produce.
            baseline_words = set(re.findall(r"[A-Za-z_]{4,}", baseline.stderr))
            new_words = set(re.findall(r"[A-Za-z_]{4,}", report.stderr)) - baseline_words
            moved = report.exit_code != 0 and bool(new_words)
            detail = (
                f"exit {baseline.exit_code} -> {report.exit_code}; "
                f"new stderr terms: {sorted(new_words)[:6]}"
            )
        else:
            before, after = baseline.get(counter), report.get(counter)
            if before is None or after is None:
                moved = False
                detail = f"{counter} absent (before={before}, after={after})"
            elif direction == "rises":
                moved = after > before
                detail = f"{counter} {before} -> {after}"
            else:
                moved = after < before
                detail = f"{counter} {before} -> {after}"

        status = "ABLE-TO-FAIL" if moved else "⛔ DEAD"
        if not moved:
            failures += 1
        print(f"{status:14} {label:28} {detail}")
        if not moved and report.stderr:
            print(f"               stderr: {report.stderr.splitlines()[0][:150]}")

    print()
    print(f"{failures} case(s) did not demonstrate a live check.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
