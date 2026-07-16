#!/usr/bin/env python3
"""Anchored, fail-closed B2 startup trust gate for traced mode."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy
import yaml

from u1_body_b2_common import (
    HERE,
    ROOT,
    UniqueKeyLoader,
    certificate_path,
    digest,
    dump_yaml,
    git_blob,
    git_object_exists,
    load_yaml,
    load_yaml_authenticated,
    manifest_absolute,
    parse_manifest_sets,
    read_bytes_authenticated,
    rel_repo,
    require,
    sha256_bytes,
    sha256_file,
    validate_commit,
)


CERT_PIN = "650656fd2ef8a87884161825d25eced2620d8602099efbb172a960073480373c"
IMMUTABLE_COMMIT = "ef934360b031bef54b37ca96d5c73cb10b0d15fd"
AMENDED_PAYLOAD = "b23993cca80dc3e6a790abcf68c1af63aa804fc47b06b153b9f224ccf27f899d"
CORE_TRACES = "7ab6a1254b1056957a6233e014d32a0d0e1ad1e2559066f7899709d07237acb0"


def manifest_git_path(raw: str) -> str:
    if raw.startswith("//"):
        return raw[2:]
    return f"software/em_charge_attribute/{raw}"


def binary_identity(command: str) -> dict[str, Any]:
    resolved = shutil.which(command)
    require(resolved is not None, "B2_A1_ENVIRONMENT", f"missing {command}")
    path = Path(resolved).resolve()
    return {"requested": command, "resolved": str(path), "sha256": sha256_file(path)}


def bootstrap_environment_context(no_network_shim: Path) -> dict[str, Any]:
    python_path = Path(sys.executable).resolve()
    wolfram = binary_identity("wolframscript")
    real_wolfram_raw = os.environ.get("WOLFRAMSCRIPT_REAL")
    require(real_wolfram_raw is not None, "B2_A1_ENVIRONMENT", "WOLFRAMSCRIPT_REAL missing")
    real_wolfram_path = Path(real_wolfram_raw).resolve()
    require(real_wolfram_path.is_file(), "B2_A1_ENVIRONMENT", "real wolframscript missing")
    version_proc = subprocess.run([wolfram["resolved"], "-version"], check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    require(version_proc.returncode == 0, "B2_A1_ENVIRONMENT", "wolframscript version probe")
    proc_version = Path("/proc/version").read_bytes()
    return {
        "record_role": "diagnostic bootstrap context only; superseded for trust by the first-use environment_closure assembled after tracing",
        "python": {
            "resolved": str(python_path),
            "sha256": sha256_file(python_path),
            "version": sys.version.replace("\n", " "),
            "implementation": platform.python_implementation(),
            "prefix": sys.prefix,
        },
        "kernel": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "proc_version_sha256": sha256_bytes(proc_version),
            "proc_version": proc_version.decode(errors="replace").strip(),
        },
        "nonstandard_package_versions_diagnostic_only": {"sympy": sympy.__version__, "pyyaml": yaml.__version__},
        "wolframscript": {**wolfram, "serialized_wrapper": True, "real_resolved": str(real_wolfram_path), "real_sha256": sha256_file(real_wolfram_path), "version_output": version_proc.stdout.strip()},
        "no_network_shim": {"path": str(no_network_shim.resolve()), "sha256": sha256_file(no_network_shim)},
    }


def read_only_root_sandbox_attestation() -> dict[str, Any]:
    require(os.environ.get("B2_READ_ONLY_ROOT_SANDBOX") == "1", "B2_A1_READ_ONLY_ROOT", "runner did not attest the v48 sandbox")
    root_rows = []
    for line in Path("/proc/self/mountinfo").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) > 6 and fields[4] == "/":
            root_rows.append({"mount_point": fields[4], "mount_options": fields[5], "read_only": "ro" in fields[5].split(",")})
    require(root_rows and all(row["read_only"] for row in root_rows), "B2_A1_READ_ONLY_ROOT", f"root mount is not read-only:{root_rows}")
    return {"status": "PASS", "mechanism": "bubblewrap --ro-bind / /", "root_mounts": root_rows, "writable_mount_policy": os.environ.get("B2_WRITABLE_MOUNT_POLICY", "missing")}


def normalized_core_traces(phase_input: dict[str, Any]) -> dict[str, Any]:
    result = {}
    for key, row in sorted(phase_input["core_traces"].items()):
        value = Fraction(str(row["value"]))
        result[key] = {"field": row["field"], "surface": row["surface"], "value": f"{float(value):.12f}"}
    return result


def verify_manifest(name: str, rows: dict[str, str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    verified, immutable = [], []
    for raw, expected in rows.items():
        path = manifest_absolute(raw)
        require(path.is_file(), "B2_A1_MANIFEST", f"{name}:missing:{raw}")
        _, authentication = read_bytes_authenticated(path, expected, f"startup:manifest:{name}")
        actual = authentication["consumed_sha256"]
        verified.append({"path": raw, "sha256": actual, "first_use_authentication": authentication})
        git_path = manifest_git_path(raw)
        if git_object_exists(IMMUTABLE_COMMIT, git_path):
            committed = sha256_bytes(git_blob(IMMUTABLE_COMMIT, git_path))
            require(committed == actual, "B2_A1_IMMUTABLE_COMMIT", git_path)
            immutable.append({"path": raw, "status": "present_and_equal", "sha256": committed})
        else:
            immutable.append({"path": raw, "status": "absent_certificate_pinned_only", "sha256": actual})
    return verified, immutable


def strict_snapshot_check(config: dict[str, Any], expected: dict[Path, str]) -> dict[str, Any]:
    results = ROOT / config["contracts"]["mutable_results"]
    result_snapshot = ROOT / config["contracts"]["b1_results_snapshot"]
    report = ROOT / config["contracts"]["mutable_report"]
    report_snapshot = ROOT / config["contracts"]["b1_report_snapshot"]
    result_live, result_live_auth = read_bytes_authenticated(results, expected[results.resolve()], "startup:mutable_results")
    result_frozen, result_frozen_auth = read_bytes_authenticated(result_snapshot, expected[result_snapshot.resolve()], "startup:results_snapshot")
    report_live, report_live_auth = read_bytes_authenticated(report, expected[report.resolve()], "startup:mutable_report")
    report_frozen, report_frozen_auth = read_bytes_authenticated(report_snapshot, expected[report_snapshot.resolve()], "startup:report_snapshot")
    results_bytes_equal = result_live == result_frozen
    report_bytes_equal = report_live == report_frozen
    results_semantic_equal = yaml.load(result_live.decode("utf-8"), Loader=UniqueKeyLoader) == yaml.load(result_frozen.decode("utf-8"), Loader=UniqueKeyLoader)
    require(results_bytes_equal or results_semantic_equal, "B2_A1_B1_CONTAINMENT", "startup flat B1 YAML snapshot")
    require(report_bytes_equal, "B2_A1_B1_CONTAINMENT", "startup B1 Markdown snapshot")
    return {
        "startup_layout": "flat_B1_pre_migration",
        "yaml_byte_equal": results_bytes_equal,
        "yaml_duplicate_key_rejected_recursive_equal": results_semantic_equal,
        "markdown_byte_equal_above_future_append_boundary": report_bytes_equal,
        "snapshot_sha256": {"yaml": result_frozen_auth["consumed_sha256"], "markdown": report_frozen_auth["consumed_sha256"]},
        "first_use_authentication": [result_live_auth, result_frozen_auth, report_live_auth, report_frozen_auth],
        "future_markdown_append_boundary": report_frozen.decode("utf-8").count("\n") + 1,
    }


def amended_digest_gate(config: dict[str, Any], certificate: dict[str, Any], expected: dict[Path, str]) -> dict[str, Any]:
    sym_path = (ROOT / config["contracts"]["sympy_b1"]).resolve(); math_path = (ROOT / config["contracts"]["mathematica_b1"]).resolve(); phase_path = (ROOT / config["contracts"]["phase_a_inputs"]).resolve()
    sym, sym_auth = load_yaml_authenticated(sym_path, expected[sym_path], "startup:amended_sympy_B1")
    math, math_auth = load_yaml_authenticated(math_path, expected[math_path], "startup:amended_wolfram_B1")
    phase_input, phase_auth = load_yaml_authenticated(phase_path, expected[phase_path], "startup:amended_phase_input")
    comparator_digest = digest(sym["amended_phase_a_payload"])
    sym_digest = sym["source_contract"]["amended_payload_sha256"]
    math_digest = math["source_contract"]["amended_payload_sha256"]
    math_payload_digest = digest(math["amended_phase_a_payload"])
    require(
        certificate["phase_a_amended_payload_sha256"] == sym_digest == math_digest == comparator_digest == math_payload_digest == AMENDED_PAYLOAD,
        "B2_A1_PHASE_A_DIGEST",
        "amended payload three-way",
    )
    core_computed = digest(normalized_core_traces(phase_input))
    require(
        certificate["core_traces_sha256"] == sym["source_contract"]["core_traces_sha256"] == math["source_contract"]["core_traces_sha256"] == core_computed == CORE_TRACES,
        "B2_A1_CORE_TRACES",
        "core traces pin",
    )
    return {
        "certificate": certificate["phase_a_amended_payload_sha256"],
        "sympy": sym_digest,
        "mathematica": math_digest,
        "comparator_recomputed_sympy_payload": comparator_digest,
        "comparator_recomputed_mathematica_payload": math_payload_digest,
        "core_traces": {"certificate": certificate["core_traces_sha256"], "comparator_recomputed": core_computed},
        "first_use_authentication": [sym_auth, math_auth, phase_auth],
    }


def build(config_path: Path, commit: str, trust_mode: str, no_network_shim: Path, pinned_environment_digest: str | None = None) -> dict[str, Any]:
    validate_commit(commit)
    require(trust_mode == "traced", "B2_A1_TRUST_MODE", "this run is required to use traced mode")
    config = load_yaml(config_path)
    require(config.get("schema_version") == "U1_PHASE_B2_STAGE0_INPUTS_V4" and config.get("directive_version") == 48 and config.get("startup_contract_commit") == commit == "5f73be9f0738030bf73165fda5432644de5f4074", "B2_A1_CONTRACT_ANCHOR", "v48 commit equation")
    sandbox = read_only_root_sandbox_attestation()
    directive = ROOT / config["contracts"]["directive"]
    certificate_path_value = ROOT / config["contracts"]["certificate"]
    certificate_commit = git_blob(commit, rel_repo(certificate_path_value))
    directive_commit = git_blob(commit, rel_repo(directive))
    certificate_working, certificate_auth = read_bytes_authenticated(certificate_path_value, CERT_PIN, "startup:approval_certificate")
    directive_expected = sha256_bytes(directive_commit)
    directive_working, directive_auth = read_bytes_authenticated(directive, directive_expected, "startup:directive")
    digests = {
        "certificate_pin": CERT_PIN,
        "certificate_commit_blob": sha256_bytes(certificate_commit),
        "certificate_working_tree": sha256_bytes(certificate_working),
        "directive_commit_blob": sha256_bytes(directive_commit),
        "directive_working_tree": sha256_bytes(directive_working),
    }
    require(digests["certificate_pin"] == digests["certificate_commit_blob"] == digests["certificate_working_tree"], "B2_A1_CONTRACT_ANCHOR", "certificate three-way")
    require(digests["directive_commit_blob"] == digests["directive_working_tree"], "B2_A1_CONTRACT_ANCHOR", "directive commit/working equality")
    certificate = yaml.load(certificate_working.decode("utf-8"), Loader=UniqueKeyLoader)
    require(certificate["status"] == "APPROVED" and certificate["git_commit"] == IMMUTABLE_COMMIT, "B2_A1_CERTIFICATE", "approval certificate status/commit")
    sets = parse_manifest_sets(certificate)
    expected_by_path = {manifest_absolute(raw).resolve(): expected for rows in sets.values() for raw, expected in rows.items()}
    results_path = (ROOT / config["contracts"]["mutable_results"]).resolve(); results_snapshot = (ROOT / config["contracts"]["b1_results_snapshot"]).resolve()
    report_path = (ROOT / config["contracts"]["mutable_report"]).resolve(); report_snapshot = (ROOT / config["contracts"]["b1_report_snapshot"]).resolve()
    expected_by_path[results_path] = expected_by_path[results_snapshot]
    expected_by_path[report_path] = expected_by_path[report_snapshot]
    manifests, immutable_rows = {}, []
    for name, rows in sets.items():
        verified, immutable = verify_manifest(name, rows)
        manifests[name] = {"count": len(verified), "entries": verified, "set_sha256": digest(verified)}
        immutable_rows.extend({"manifest": name, **row} for row in immutable)
    governing = sorted(
        raw for raw in sets["Dep_B1"]
        if raw.startswith("directive_") or raw in {"//docs/em_u1_body_definition.md", "//docs/development_pipeline.md"}
    )
    required_governing = {
        "directive_u1_body_dynamics.md", "directive_u1_body_dynamics_rebuild.md", "directive_u1_phaseA_remediation.md",
        "directive_u1_phaseB_mechanics.md", "directive_u1_phaseB_rebuild.md", "directive_u1_phaseB_remediation2.md", "directive_u1_phaseB_remediation3.md",
        "//docs/em_u1_body_definition.md", "//docs/development_pipeline.md",
    }
    require(set(governing) == required_governing, "B2_A1_GOVERNING_SET", "governing-contract enumeration")
    bootstrap_environment = bootstrap_environment_context(no_network_shim)
    if pinned_environment_digest is not None:
        require(len(pinned_environment_digest) == 64 and all(ch in "0123456789abcdef" for ch in pinned_environment_digest), "B2_A1_ENVIRONMENT_PIN", "complete closure digest syntax; full assertion is the preceding environment gate")
    return {
        "schema_version": "U1_PHASE_B2_STARTUP_GATE_V3",
        "status": "PASS",
        "startup_contract_commit": commit,
        "trust_mode": trust_mode,
        "contract_digests": digests,
        "certificate": {"status": certificate["status"], "immutable_commit": certificate["git_commit"], "schema_version": certificate["schema_version"]},
        "manifest_verification": manifests,
        "immutable_commit_comparison": {"commit": IMMUTABLE_COMMIT, "entries": immutable_rows},
        "set_definitions": {
            "P": "path set of certificate.protected_artifact_sha256",
            "S": "path set of certificate.substrate_reference_sha256",
            "Dep_B1": "path set of certificate.runtime_dependency_sha256",
            "G": governing,
        },
        "first_use_authentication": [certificate_auth, directive_auth],
        "mutable_aggregate_startup_containment": strict_snapshot_check(config, expected_by_path),
        "phase_a_amended_digest_gate": amended_digest_gate(config, certificate, expected_by_path),
        "bootstrap_environment_context": bootstrap_environment,
        "read_only_root_sandbox": sandbox,
        "environment_closure_digest_supplied": pinned_environment_digest,
        "environment_pin_assertion_owner": "u1_body_b2_environment_gate.py complete closure rehash" if pinned_environment_digest is not None else "stage0 first-use trace assembly pending",
        "completion_recheck_required_against_same_commit": True,
        "checks": {
            "B2_A1_CONTRACT_ANCHOR": "PASS", "B2_A1_CERTIFICATE": "PASS", "B2_A1_MANIFEST": "PASS",
            "B2_A1_IMMUTABLE_COMMIT": "PASS", "B2_A1_B1_CONTAINMENT": "PASS", "B2_A1_PHASE_A_DIGEST": "PASS",
            "B2_A1_CORE_TRACES": "PASS", "B2_A1_ENVIRONMENT": "PASS(first_use_audit_bootstrap)", "B2_A1_ENVIRONMENT_PIN": "PASS(complete_environment_gate_preceded_startup)" if pinned_environment_digest is not None else "PENDING_STAGE0_FIRST_USE_CLOSURE", "B2_A1_TRUST_MODE": "PASS", "B2_A1_READ_ONLY_ROOT": "PASS",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--startup-contract-commit", required=True)
    parser.add_argument("--trust-mode", required=True)
    parser.add_argument("--no-network-shim", type=Path, required=True)
    parser.add_argument("--environment-identity-digest")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        artifact = build(args.input.resolve(), args.startup_contract_commit, args.trust_mode, args.no_network_shim.resolve(), args.environment_identity_digest)
        dump_yaml(args.output.resolve(), artifact)
        print(f"B2_STARTUP_GATE: PASS commit={args.startup_contract_commit} mode={args.trust_mode}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
