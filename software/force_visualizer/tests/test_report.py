"""Headless numeric-report acceptance check from build spec section 6.1."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PACKAGE_ROOT.parents[1]
REPORT_PATH = PACKAGE_ROOT / "report.py"
OUTPUT_PATH = PACKAGE_ROOT / "output" / "verification_report.txt"


def test_report_runs_headlessly_and_writes_stdout_transcript() -> None:
    completed = subprocess.run(
        [sys.executable, str(REPORT_PATH)],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    assert completed.stdout == OUTPUT_PATH.read_text(encoding="utf-8")
    assert "GRAVITY" in completed.stdout
    assert "LIGHT" in completed.stdout
    assert "CHARGE" in completed.stdout
    assert "MAGNETISM" in completed.stdout
    assert "CHARACTERIZED DEPARTURES" in completed.stdout
    assert "SUMMARY: 18 checks passed of 18 total." in completed.stdout
    assert "gravity drain-field direction" in completed.stdout
    assert "electric-field directions" in completed.stdout
    assert "magnetic-field circulation" in completed.stdout
    assert "[CALIBRATED-MAGNITUDE] monopole residual amplitude=" in completed.stdout
    assert "[CALIBRATED-MAGNITUDE] dipole residual amplitude=" in completed.stdout
    assert "[CALIBRATED-MAGNITUDE] preferred-frame residual=" in completed.stdout
    assert "zero by calibrated c_E=c_gamma default" in completed.stdout
