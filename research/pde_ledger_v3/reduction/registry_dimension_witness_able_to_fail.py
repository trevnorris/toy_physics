#!/usr/bin/env python3
"""Calibrate inherited, multiplier, echoed, and branch-bound witness sources."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

import yaml


HERE = Path(__file__).resolve().parent
LEDGER = HERE.parent
CHECKER = HERE / "registry_dimension_witness.py"
MANIFEST = HERE / "registry_dimension_witness.yaml"
QUANTITIES = HERE / "quantities.yaml"
CASES = ("inherited-config", "multiplier", "echoed", "branch-dimension")
CASE_TARGETS = {
    "inherited-config": ("S9-py", "Q.brane.rho_br"),
    "multiplier": ("S9-py", "Q.brane.c_gamma"),
    "echoed": ("S10-py", "Q.brane.D_brane"),
    "branch-dimension": ("S10-py", "Q.brane.rho_br"),
}
EXPECTED_TRANSITIONS = {
    "inherited-config": ("AGREEMENT", "DISAGREEMENT"),
    "multiplier": ("AGREEMENT", "DISAGREEMENT"),
    "echoed": ("ECHOED", "ECHOED_MISMATCH"),
    "branch-dimension": ("AGREEMENT", "BRANCH_DIMENSION_MISMATCH"),
}


def _read_yaml(path: Path) -> dict[str, object]:
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise ValueError(f"{path}: top level must be a mapping")
    return document


def _artifact_manifest(case: str, branch_dimension: int = 3) -> dict[str, object]:
    target_artifact, _target_qid = CASE_TARGETS[case]
    manifest = _read_yaml(MANIFEST)
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise ValueError("manifest artifacts must be a list")
    selected = [
        dict(artifact)
        for artifact in artifacts
        if isinstance(artifact, dict) and artifact.get("id") == target_artifact
    ]
    if len(selected) != 1:
        raise ValueError(f"expected one {target_artifact} artifact")
    artifact = selected[0]
    if "config" in artifact:
        config = Path(str(artifact["config"]))
        artifact["config"] = str(config if config.is_absolute() else HERE / config)
    output = Path(str(artifact["output"]))
    artifact["output"] = str(output if output.is_absolute() else HERE / output)

    if case == "inherited-config":
        artifact["extra_dimension_sources"] = {}
    elif case == "multiplier":
        extras = artifact.get("extra_dimension_sources")
        if not isinstance(extras, dict) or "Q.brane.c_gamma" not in extras:
            raise ValueError("multiplier source is missing")
        artifact["selection"] = {"packages": []}
        artifact["extra_dimension_sources"] = {
            "Q.brane.c_gamma": extras["Q.brane.c_gamma"]
        }
    elif case == "echoed":
        extras = artifact.get("extra_dimension_sources")
        if not isinstance(extras, dict) or "Q.brane.D_brane" not in extras:
            raise ValueError("echoed source is missing")
        artifact["selection"] = {"packages": []}
        artifact["extra_dimension_sources"] = {
            "Q.brane.D_brane": extras["Q.brane.D_brane"]
        }
    elif case == "branch-dimension":
        artifact["selection"] = {
            "packages": ["MAIN"],
            "dimensions": [branch_dimension],
        }
        artifact["extra_dimension_sources"] = {}
    else:
        raise ValueError(f"unknown case {case}")
    manifest["artifacts"] = [artifact]
    return manifest


def _write_yaml(path: Path, document: object) -> None:
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")


def _mutate_declared_exponent(destination: Path) -> None:
    quantities = _read_yaml(QUANTITIES)
    rows = quantities.get("quantities")
    if not isinstance(rows, list):
        raise ValueError("quantities list is missing")
    targets = [
        row
        for row in rows
        if isinstance(row, dict) and row.get("qid") == "Q.brane.rho_br"
    ]
    if len(targets) != 1:
        raise ValueError("expected one Q.brane.rho_br declaration")
    dimension = targets[0].get("dimension")
    if not isinstance(dimension, dict) or not isinstance(dimension.get("exponents"), list):
        raise ValueError("Q.brane.rho_br has no exponent list")
    exponents = dimension["exponents"]
    if len(exponents) != 3 or not isinstance(exponents[0], int):
        raise ValueError("Q.brane.rho_br does not have an integer exponent triple")
    exponents[0] += 1
    _write_yaml(destination, quantities)


def _mutate_multiplier(manifest: dict[str, object]) -> None:
    artifact = manifest["artifacts"][0]  # type: ignore[index]
    sources = artifact["extra_dimension_sources"]["Q.brane.c_gamma"]  # type: ignore[index]
    if not isinstance(sources, list) or len(sources) != 1 or not isinstance(sources[0], dict):
        raise ValueError("expected one c_gamma multiplier source")
    if sources[0].get("multiplier") != 2:
        raise ValueError("c_gamma multiplier baseline is not 2")
    sources[0]["multiplier"] = 3


def _mutate_echoed_output(source: Path, destination: Path) -> None:
    target_tag = "PY_S10_LOCAL_REGISTRY_D_BRANE_DIMENSION"
    lines = source.read_text(encoding="utf-8").splitlines()
    matches = [index for index, line in enumerate(lines) if line.startswith(target_tag + ":")]
    if len(matches) != 1:
        raise ValueError(f"expected one {target_tag} output record")
    index = matches[0]
    mutated = lines[index].replace("(0, 0, 0)", "(1, 0, 0)", 1)
    if mutated == lines[index]:
        raise ValueError(f"{target_tag} did not contain the expected vector")
    lines[index] = mutated
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _target_line(case: str, stdout: str) -> str:
    artifact, qid = CASE_TARGETS[case]
    prefix = f"WITNESS|artifact={artifact}|quantity={qid}|"
    rows = [line for line in stdout.splitlines() if line.startswith(prefix)]
    if len(rows) != 1:
        raise ValueError(f"expected one target witness line, found {len(rows)}")
    return rows[0]


def _status(line: str) -> str:
    for field in line.split("|"):
        if field.startswith("status="):
            return field.partition("=")[2]
    raise ValueError("target witness line has no status")


def _run(manifest: Path, quantities: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(CHECKER),
            "--manifest",
            str(manifest),
            "--quantities",
            str(quantities),
        ],
        cwd=LEDGER,
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )


def run_case(case: str) -> int:
    expected_baseline, expected_perturbed = EXPECTED_TRANSITIONS[case]
    with tempfile.TemporaryDirectory(prefix="registry-dimension-witness-") as raw_temp:
        temporary = Path(raw_temp)
        baseline_manifest = _artifact_manifest(case)
        baseline_manifest_path = temporary / "baseline-manifest.yaml"
        _write_yaml(baseline_manifest_path, baseline_manifest)

        perturbed_manifest = _artifact_manifest(
            case, branch_dimension=5 if case == "branch-dimension" else 3
        )
        perturbed_manifest_path = temporary / "perturbed-manifest.yaml"
        perturbed_quantities = QUANTITIES
        if case == "inherited-config":
            perturbed_quantities = temporary / "quantities.yaml"
            _mutate_declared_exponent(perturbed_quantities)
        elif case == "multiplier":
            _mutate_multiplier(perturbed_manifest)
        elif case == "echoed":
            artifact = perturbed_manifest["artifacts"][0]  # type: ignore[index]
            source_output = Path(str(artifact["output"]))
            perturbed_output = temporary / source_output.name
            _mutate_echoed_output(source_output, perturbed_output)
            artifact["output"] = str(perturbed_output)
        _write_yaml(perturbed_manifest_path, perturbed_manifest)

        baseline = _run(baseline_manifest_path, QUANTITIES)
        perturbed = _run(perturbed_manifest_path, perturbed_quantities)
        if baseline.stderr or perturbed.stderr:
            print(
                "CALIBRATION_OPERATIONAL_STDERR|baseline="
                + (baseline.stderr.strip() or "none")
                + "|perturbed="
                + (perturbed.stderr.strip() or "none")
            )
            return 2
        baseline_line = _target_line(case, baseline.stdout)
        perturbed_line = _target_line(case, perturbed.stdout)
        baseline_status, perturbed_status = _status(baseline_line), _status(perturbed_line)
        deviations = sum(
            (
                baseline.returncode != 0,
                perturbed.returncode != 1,
                baseline_status != expected_baseline,
                perturbed_status != expected_perturbed,
            )
        )
        print(
            f"CALIBRATION_GUARD|case={case}|baseline=<{baseline_line}>"
            f"|perturbed=<{perturbed_line}>"
            f"|expected_transition={expected_baseline}->{expected_perturbed}"
            f"|observed_transition={baseline_status}->{perturbed_status}"
            f"|exit_codes={baseline.returncode}->{perturbed.returncode}"
            f"|deviation_count={deviations}"
        )
        return 1 if deviations else 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, choices=CASES)
    args = parser.parse_args(argv)
    try:
        return run_case(args.case)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"CALIBRATION_OPERATIONAL_FAILURE: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
