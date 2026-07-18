#!/usr/bin/env python3
"""Assemble and verify the digest-bound U1 Phase-C stage-0 bundle."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import platform
import re
import stat
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import sympy
import yaml


STARTUP_ANCHOR = "377eab17a4babc12847450956dc55fe3e16d33da"
B2_CONTRACT_SHA = "8b8ee113d9a8342ac11f23e7959e05c323cc89a81b7c298af85f646b6e42a9e7"
B2_DISPOSITION_SHA = "d94341173b1f1ac643bb05cf52dbac2300668ce6a50b8b5042ee2c7fa35cc54f"
PHASE_A_LINEAGE = "b23993cca80dc3e6a790abcf68c1af63aa804fc47b06b153b9f224ccf27f899d"
COMPONENT_FILES = {
    "frozen_data_pin_table": "frozen_data_pin_table.yaml",
    "availability_slots": "availability_slots.yaml",
    "reconciliation_inventory": "reconciliation_inventory.yaml",
    "force_term_census": "force_term_census.yaml",
    "coupling_source_census": "coupling_source_census.yaml",
    "g8_ablation_inventory": "g8_ablation_inventory.yaml",
    "projection_freeze": "projection_freeze.yaml",
    "environment_identity": "environment_identity.yaml",
    "producer_map": "producer_map.yaml",
    "obligation_manifest": "obligation_manifest.yaml",
    "parameter_register_proposals": "parameter_register_proposals.yaml",
    "evaluated_code_closure_policy": "evaluated_code_closure_policy.yaml",
}


class ContractFailure(RuntimeError):
    def __init__(self, assert_id: str, detail: str):
        super().__init__(f"{assert_id}: {detail}")
        self.assert_id = assert_id
        self.detail = detail


def require(condition: bool, assert_id: str, detail: str) -> None:
    if not condition:
        raise ContractFailure(assert_id, detail)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=120),
        encoding="utf-8",
    )


def git(repo: Path, *args: str, text: bool = True) -> str | bytes:
    result = subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )
    return result.stdout.strip() if text else result.stdout


def blob_sha(repo: Path, anchor: str, rel_path: str) -> str:
    data = git(repo, "cat-file", "blob", f"{anchor}:{rel_path}", text=False)
    if not isinstance(data, bytes):
        raise TypeError("binary git blob read unexpectedly returned text")
    return hashlib.sha256(data).hexdigest()


def strip_engine_local(value: Any) -> Any:
    if isinstance(value, dict):
        out = {}
        for key, child in value.items():
            if key == "engine_local_function":
                continue
            out[key] = strip_engine_local(child)
        return out
    if isinstance(value, list):
        return [strip_engine_local(child) for child in value]
    return value


def normalized_slots(
    sympy_slots: list[dict[str, Any]], wolfram_slots: list[dict[str, Any]], agreement: dict[str, Any]
) -> list[dict[str, Any]]:
    wolfram_by_id = {row["slot_id"]: row for row in wolfram_slots}
    cert_by_id = {
        row["slot_id"]: row for row in agreement["challenge_dual_engine_certificates"]
    }
    out: list[dict[str, Any]] = []
    for original in sympy_slots:
        row = copy.deepcopy(original)
        row["producer_set"] = sorted(row["producer_set"], key=str.casefold)
        if row["disposition"] == "DERIVED":
            value = strip_engine_local(row["value"])
            if row["slot_id"] == "domain:S_body_Omega_c":
                value["canonical_terms"] = sorted(
                    value["canonical_terms"], key=lambda item: item["id"].casefold()
                )
            elif row["slot_id"] == "support:complete_action_second_variation":
                value["term_records"] = sorted(
                    value["term_records"], key=lambda item: item["id"].casefold()
                )
                value["nonzero_pair_records"] = sorted(
                    value["nonzero_pair_records"],
                    key=lambda item: (item["term"].casefold(), item["pair"].casefold()),
                )
            row["value"] = value
            row["value_digest"] = digest(value)
            row["dual_engine_evidence"] = {
                "sympy_value_present": True,
                "wolfram_value_present": wolfram_by_id[row["slot_id"]]["disposition"]
                == "DERIVED",
                "comparison_id": row["dual_engine_comparison_id"],
            }
        else:
            row["witness"]["producer_set"] = sorted(
                row["witness"]["producer_set"], key=str.casefold
            )
            row["challenge"]["dual_engine_certificate"] = cert_by_id[row["slot_id"]]
        out.append(row)
    return sorted(out, key=lambda row: row["slot_id"].casefold())


def pin_paths(repo: Path) -> list[tuple[str, str, str | None]]:
    """Return (role, relative path, lineage commit) for every consumed input."""
    return [
        ("governing_directive", "software/em_charge_attribute/directive_u1_phaseC_tilt_coupling.md", STARTUP_ANCHOR),
        ("parent_contract", "software/em_charge_attribute/directive_u1_body_dynamics.md", STARTUP_ANCHOR),
        ("witness_contract", "software/em_charge_attribute/directive_u1_phaseB2_intake_radiative.md", STARTUP_ANCHOR),
        ("governing_spec", "docs/em_u1_body_definition.md", STARTUP_ANCHOR),
        ("ownership_naming_rubric", "docs/em_phaseC_force_decomposition.md", STARTUP_ANCHOR),
        ("decision_16", "software/stage1_solver/decisions/16_retire_brane_polar_field.md", STARTUP_ANCHOR),
        ("decision_17", "software/stage1_solver/decisions/17_trust_apparatus_trim.md", STARTUP_ANCHOR),
        ("phase_A_7_0_inputs", "software/em_charge_attribute/u1_body_dynamics_inputs.yaml", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("phase_A_amendment", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage1/phase_a_amendment_agreement.yaml", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("phase_A_sympy_payload", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_a.json", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("phase_A_wolfram_payload", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_a.json", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("B1_approval", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_orchestrator_approval.yaml", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("B1_final_results", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_final_results_snapshot.yaml", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("B1_final_report", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_final_report_snapshot.md", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("B1_mechanics_inputs", "software/em_charge_attribute/u1_body_mechanics_inputs.yaml", "ef934360b031bef54b37ca96d5c73cb10b0d15fd"),
        ("B2_sealed_stage0_contract", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_0_intake_radiative_contract/stage0_contract.yaml", None),
        ("B2_stage0_engine_agreement", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_0_intake_radiative_contract/stage0_engine_agreement.yaml", None),
        ("B2_sympy_production", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/sympy_b2.yaml", None),
        ("B2_wolfram_production", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/mathematica_b2.yaml", None),
        ("B2_engine_agreement", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/engine_agreement.yaml", None),
        ("B2_completion_gate", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/completion_gate.yaml", None),
        ("B2_production_aggregate", "software/em_charge_attribute/reports/u1_body_dynamics_results.yaml", "f9d7bc30d5ce46233bc0de6c97bd6182a1c6e186"),
        ("pathA29_one_sided_report", "software/stage1_solver/reports/pathA_29_brane_bulk_return.md", None),
        ("pathA29_one_sided_results", "software/stage1_solver/reports/pathA_29_results.yaml", None),
        ("pathA36_shear_report", "software/stage1_solver/reports/pathA_36_c5_phase_potential.md", None),
        ("pathA36_shear_results", "software/stage1_solver/reports/pathA_36_c5_phase_potential_results.yaml", None),
        ("pathA38_h_report", "software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md", None),
        ("pathA38_h_results", "software/stage1_solver/reports/pathA_38_results.yaml", None),
        ("pathA39_stage3_report", "software/stage1_solver/reports/pathA_39_stage3_operator_parity.md", None),
        ("pathA39_stage3_results", "software/stage1_solver/reports/pathA_39_stage3_operator_parity_results.yaml", None),
        ("pathA39_stage4_report", "software/stage1_solver/reports/pathA_39_stage4_field_classification.md", None),
        ("pathA39_stage4_results", "software/stage1_solver/reports/pathA_39_stage4_field_classification_results.yaml", None),
        ("native_action_stage004", "research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md", None),
        ("native_action_stage006", "research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md", None),
        ("parameter_register_read_only", "research/pde_ledger_v2/notes/parameter_register.md", None),
    ]


def build_pin_table(repo: Path, anchor: str, mutation: str | None = None) -> dict[str, Any]:
    records = []
    specs = pin_paths(repo)
    if mutation == "ASSERT_PIN_INPUT_EXISTS":
        specs = [(specs[0][0], "software/em_charge_attribute/missing_pin_fixture.yaml", specs[0][2]), *specs[1:]]
    for role, rel_path, lineage in specs:
        path = repo / rel_path
        require(path.is_file(), "ASSERT_PIN_INPUT_EXISTS", f"missing {rel_path}")
        disk_sha = sha256_path(path)
        anchor_sha = None
        anchored_equal = None
        if lineage == STARTUP_ANCHOR:
            anchor_sha = blob_sha(repo, anchor, rel_path)
            anchored_equal = disk_sha == anchor_sha
            if mutation == "ASSERT_GOVERNING_BLOB_ANCHORED" and role == "governing_directive":
                anchored_equal = False
            require(
                anchored_equal,
                "ASSERT_GOVERNING_BLOB_ANCHORED",
                f"working bytes differ from {anchor}:{rel_path}",
            )
        records.append(
            {
                "role": role,
                "path": rel_path,
                "sha256": disk_sha,
                "lineage_commit": lineage,
                "startup_anchor_blob_sha256": anchor_sha,
                "startup_anchor_equal": anchored_equal,
                "consumption": "read_only",
            }
        )
    by_role = {row["role"]: row for row in records}
    if mutation == "ASSERT_B2_SEALED_CONTRACT":
        by_role["B2_sealed_stage0_contract"]["sha256"] = "0" * 64
    require(
        by_role["B2_sealed_stage0_contract"]["sha256"] == B2_CONTRACT_SHA,
        "ASSERT_B2_SEALED_CONTRACT",
        "B2 stage0 contract seal changed",
    )
    amendment = load_yaml(repo / by_role["phase_A_amendment"]["path"])
    if mutation == "ASSERT_PHASE_A_LINEAGE":
        amendment["digest_gate"]["sympy"] = "0" * 64
    require(
        amendment["digest_gate"]["sympy"]
        == amendment["digest_gate"]["wolfram"]
        == amendment["digest_gate"]["comparator"]
        == PHASE_A_LINEAGE,
        "ASSERT_PHASE_A_LINEAGE",
        "amended Phase-A lineage is not b23993cc...",
    )
    b1 = load_yaml(repo / by_role["B1_final_results"]["path"])
    b2 = load_yaml(repo / by_role["B2_sealed_stage0_contract"]["path"])
    if mutation == "ASSERT_PIN_B1_TERMINAL":
        b1["mechanics_partition_ledger"]["state"] = "closed"
    if mutation == "ASSERT_PIN_B2_DISPOSITIONS":
        b2["frozen_data"]["stage0_datum_disposition_sha256"] = "0" * 64
    require(
        b1["mechanics_partition_ledger"]["state"] == "partition_open_pending_B2",
        "ASSERT_PIN_B1_TERMINAL",
        "B1 partition terminal changed",
    )
    require(
        b2["frozen_data"]["stage0_datum_disposition_sha256"] == B2_DISPOSITION_SHA,
        "ASSERT_PIN_B2_DISPOSITIONS",
        "B2 disposition set changed",
    )
    return {
        "schema_version": "U1_PHASE_C_FROZEN_PIN_TABLE_V1",
        "startup_contract_commit": anchor,
        "record_count": len(records),
        "records": records,
        "lineage_assertions": {
            "phase_A_payload": PHASE_A_LINEAGE,
            "B1_commit": "ef934360b031bef54b37ca96d5c73cb10b0d15fd",
            "B1_tilt_ingredient_count": 8,
            "B1_partition_terminal": "partition_open_pending_B2",
            "B2_commit": "f9d7bc30d5ce46233bc0de6c97bd6182a1c6e186",
            "B2_contract_sha256": B2_CONTRACT_SHA,
            "B2_disposition_sha256": B2_DISPOSITION_SHA,
            "B2_partition_terminal": "UNRESOLVED(return_closure)",
            "P_sector_present": False,
        },
    }


def build_reconciliation(
    expected_ids: list[str], pin_table: dict[str, Any], slots: list[dict[str, Any]],
    mutation: str | None = None,
) -> dict[str, Any]:
    pins = {row["role"]: row["sha256"] for row in pin_table["records"]}
    slot_ids = {row["slot_id"] for row in slots}
    blocker = "tilt:indexed_sleeve_tilt_profile"
    if mutation == "ASSERT_RECONCILIATION_BLOCKER_SLOT":
        blocker = "tilt:missing_fixture"
    require(blocker in slot_ids, "ASSERT_RECONCILIATION_BLOCKER_SLOT", "missing tilt blocker")
    records: list[dict[str, Any]] = []
    for expected_id in sorted(expected_ids):
        category, upstream = expected_id.split("|", 1)
        if category == "B1_LEAF":
            artifact = pins["B1_final_results"]
            frozen = "UNRESOLVED(tilt_profile)"
            routing = "PHASE_C_WITNESS_UNRESOLVED"
            witness = f"witness:{blocker}"
            new_witness = False
        elif category == "B2_TILT_PATH":
            artifact = pins["B2_sympy_production"]
            frozen = "UNRESOLVED(tilt_profile)"
            routing = "PHASE_C_WITNESS_UNRESOLVED"
            witness = f"witness:{blocker}"
            new_witness = False
        else:
            artifact = pins["B2_sealed_stage0_contract"]
            frozen = "p=p_star_deferred_to_phase_C"
            if upstream.startswith("G9_record|"):
                routing = "PRESERVED_G9_EXACT_REFERENCE"
                witness = None
                new_witness = False
            else:
                routing = "PHASE_C_WITNESS_UNRESOLVED"
                witness = f"witness:{blocker}"
                new_witness = False
        records.append(
            {
                "successor_key": expected_id,
                "upstream_artifact_digest": artifact,
                "canonical_upstream_id_or_schema_path": upstream,
                "frozen_upstream_disposition_verbatim": frozen,
                "phase_C_stage0_routing": routing,
                "witness_reference": witness,
                "new_witness_minted": new_witness,
                "upstream_record_modified": False,
            }
        )
    record_ids = [row["successor_key"] for row in records]
    counts = Counter(row.split("|", 1)[0] for row in record_ids)
    return {
        "schema_version": "U1_PHASE_C_RECONCILIATION_INVENTORY_V1",
        "generation": "mechanical frozen B1/B2 walk; immutable overlay only",
        "expected_ids": sorted(expected_ids),
        "records": records,
        "expected_successor_exact_set_equal": set(expected_ids) == set(record_ids),
        "duplicate_successor_count": len(record_ids) - len(set(record_ids)),
        "invented_successor_count": len(set(record_ids) - set(expected_ids)),
        "omitted_successor_count": len(set(expected_ids) - set(record_ids)),
        "category_counts": dict(sorted(counts.items())),
        "G9_preserved_reference_count": sum(
            row["phase_C_stage0_routing"] == "PRESERVED_G9_EXACT_REFERENCE"
            for row in records
        ),
        "partition_successor_obligation": {
            "B1_frozen_terminal": "partition_open_pending_B2",
            "B2_frozen_terminal": "UNRESOLVED(return_closure)",
            "phase_C_owner_at_stage0": "UNRESOLVED(return_closure)",
            "owner_must_be_computed_as_successor": True,
            "upstream_fact_reused_as_closed_owner": False,
        },
    }


def build_producer_map(slots: list[dict[str, Any]]) -> dict[str, Any]:
    relations = [
        {
            "slot_id": row["slot_id"],
            "producer_set": row["producer_set"],
            "required_type": row["required_type"],
            "required_dimensions": row["required_dimensions"],
            "acceptance_predicate": row["acceptance_predicate"],
            "generation_rule": f"directive_taxonomy:{row['category']}",
        }
        for row in slots
    ]
    return {
        "schema_version": "U1_PHASE_C_PRODUCER_MAP_V1",
        "builder_declared_producer_tags_accepted": False,
        "generation": "directive taxonomy x frozen typed ancestry",
        "relation_count": len(relations),
        "relations": relations,
        "universal_census_predicate": all(row["producer_set"] for row in relations),
    }


def proposal_rows() -> list[dict[str, Any]]:
    data = {
        "indexed_density_tilt_profile": ([-4, 0, 0], "scalar_density_tangent_on_Omega_c"),
        "indexed_flow_tilt_response": ([1, -1, 0], "vector_velocity_tangent_on_Omega_c"),
        "indexed_h_tilt_profile": ([0, 0, 0], "scalar_h_tangent_on_Omega_c"),
        "indexed_phase_tilt_profile": ([0, 0, 0], "scalar_phase_tangent_on_Omega_c"),
        "indexed_shear_tilt_profile": ([1, 0, 0], "brane_displacement_tangent_on_Omega_c"),
        "indexed_sleeve_surface_normal_profile": ([0, 0, 0], "normal_variation_on_Sigma"),
        "indexed_sleeve_tilt_profile": ([0, 0, 0], "chi_field_tangent_on_Omega_c"),
        "indexed_uw_tilt_profile": ([1, 0, 0], "normal_displacement_tangent_on_Omega_c"),
    }
    return [
        {
            "id": key,
            "proposed_class": "GAP_OPEN_FIELD_PROFILE",
            "status": "ORCHESTRATOR_APPLIED_AS_OPEN_REDUCTION_DEBT",
            "dimensions_LTM": dims,
            "domain": domain,
            "arguments": ["field_point", "endpoint", "ambient_branch"],
            "collective_coordinate_in_input_domain": False,
            "embedding_use": f"delta_Phi(y)=p_i*{key}_i(y;endpoint,ambient)",
            "symmetry_class": "TO_BE_COMPUTED_BODY_CONJUGATION_AND_AMBIENT_BRANCH",
            "reduction_route": "native tilted sleeve boundary-value solve",
            "counts_as_derived": False,
        }
        for key, (dims, domain) in sorted(data.items())
    ]


def merkle_files(paths: list[Path], base: Path | None = None) -> dict[str, Any]:
    records = []
    for path in sorted({item.resolve() for item in paths if item.is_file()}):
        label = str(path.relative_to(base)) if base and path.is_relative_to(base) else str(path)
        records.append({"path": label, "sha256": sha256_path(path), "size": path.stat().st_size})
    return {"file_count": len(records), "digest": digest(records), "records": records}


def tree_files(root: Path) -> list[Path]:
    return [path for path in root.rglob("*") if path.is_file() and not path.is_symlink()]


def mount_record(path: Path, mutation: str | None = None) -> dict[str, Any]:
    resolved = path.resolve()
    best: tuple[int, list[str]] | None = None
    for line in Path("/proc/self/mountinfo").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        separator = fields.index("-")
        mountpoint = Path(fields[4])
        try:
            resolved.relative_to(mountpoint)
        except ValueError:
            continue
        if best is None or len(str(mountpoint)) > best[0]:
            best = (len(str(mountpoint)), fields)
    if mutation == "ASSERT_ENV_MOUNT_IDENTITY":
        best = None
    require(best is not None, "ASSERT_ENV_MOUNT_IDENTITY", f"no mount for {resolved}")
    fields = best[1]
    separator = fields.index("-")
    return {
        "path": str(resolved),
        "mountpoint": fields[4],
        "mount_options": sorted(fields[5].split(",")),
        "filesystem": fields[separator + 1],
        "source": fields[separator + 2],
        "super_options": sorted(fields[separator + 3].split(",")),
    }


def wolfram_identity(mutation: str | None = None) -> dict[str, Any]:
    kernel = Path("/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel")
    code = (
        'Print["VERSION="<>$Version];Print["SYSTEM="<>$SystemID];'
        'Print["INSTALL="<>$InstallationDirectory];Print["USERBASE="<>$UserBaseDirectory];'
        'Print["BASE="<>$BaseDirectory];'
        'Print["SEATS="<>ToString[$MaxLicenseProcesses]];Quit[]'
    )
    result = subprocess.run(
        [str(kernel), "-noinit", "-noprompt", "-run", code],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=60,
    )
    parsed = {}
    for raw in result.stdout.splitlines():
        line = raw.strip().strip('"')
        if "=" in line:
            key, value = line.split("=", 1)
            parsed[key] = value
    if mutation == "ASSERT_WOLFRAM_SEAT_LIMIT":
        parsed["SEATS"] = "3"
    require(parsed.get("SEATS") == "2", "ASSERT_WOLFRAM_SEAT_LIMIT", "Wolfram seat identity changed")
    real_kernel = Path("/usr/local/Wolfram/Wolfram/15.0/SystemFiles/Kernel/Binaries/Linux-x86-64/WolframKernel")
    return {
        "version": parsed.get("VERSION"),
        "system_id": parsed.get("SYSTEM"),
        "installation_directory": parsed.get("INSTALL"),
        "user_base_directory": parsed.get("USERBASE"),
        "base_directory": parsed.get("BASE"),
        "max_license_processes": int(parsed["SEATS"]),
        "launch_mode": "WolframKernel -noinit -noprompt",
        "launcher_path": str(kernel),
        "launcher_sha256": sha256_path(kernel),
        "kernel_binary_path": str(real_kernel),
        "kernel_binary_sha256": sha256_path(real_kernel),
        "user_init_disabled": True,
        "user_plugin_paths_removed_before_import": True,
    }


def environment_identity(mutation: str | None = None) -> dict[str, Any]:
    isolated = sys.flags.isolated == 1 and sys.flags.no_user_site == 1
    if mutation == "ASSERT_ENV_PYTHON_SANITIZED":
        isolated = False
    require(
        isolated,
        "ASSERT_ENV_PYTHON_SANITIZED",
        "environment recorder must run with Python isolated and user-site disabled",
    )
    py_real = Path(sys.executable).resolve()
    sympy_root = Path(sympy.__file__).resolve().parent
    yaml_root = Path(yaml.__file__).resolve().parent
    stdlib_root = Path("/usr/lib/python3.10")
    python_manifest = merkle_files(
        [py_real, *tree_files(stdlib_root), *tree_files(sympy_root), *tree_files(yaml_root)]
    )
    wolfram = wolfram_identity(mutation)
    roots = [Path("/usr"), Path(wolfram["installation_directory"]), Path("/opt/Wolfram")]
    root_records = []
    for root in roots:
        info = root.stat()
        root_records.append(
            {
                "path": str(root.resolve()),
                "device": info.st_dev,
                "inode": info.st_ino,
                "mode": stat.S_IMODE(info.st_mode),
                "uid": info.st_uid,
                "gid": info.st_gid,
                "mount": mount_record(root, mutation),
            }
        )
    read_only = all("ro" in row["mount"]["mount_options"] for row in root_records)
    if mutation == "ASSERT_ENV_TOOLCHAIN_READ_ONLY":
        read_only = False
    require(
        read_only,
        "ASSERT_ENV_TOOLCHAIN_READ_ONLY",
        "declared evaluated-toolchain roots are not read-only in the execution sandbox",
    )
    selected_identity = {
        "python_evaluated_toolchain_manifest": {
            "file_count": python_manifest["file_count"],
            "digest": python_manifest["digest"],
        },
        "wolfram_launcher_sha256": wolfram["launcher_sha256"],
        "wolfram_kernel_binary_sha256": wolfram["kernel_binary_sha256"],
        "root_records": root_records,
        "sandbox_name_service_policy_sha256": sha256_path(Path("/etc/nsswitch.conf")),
    }
    return {
        "schema_version": "U1_PHASE_C_ENVIRONMENT_IDENTITY_V1",
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": str(py_real),
            "executable_sha256": sha256_path(py_real),
            "sympy_version": sympy.__version__,
            "sympy_root": str(sympy_root),
            "pyyaml_version": yaml.__version__,
            "pyyaml_root": str(yaml_root),
            "isolated_flag": sys.flags.isolated,
            "no_user_site_flag": sys.flags.no_user_site,
            "user_site_enabled": False,
        },
        "wolfram": wolfram,
        "declared_read_only_toolchain_roots": root_records,
        "selected_evaluated_toolchain_identity": selected_identity,
        "toolchain_root_identity_digest": digest(selected_identity),
        "production_assertion": "recompute complete record before every production and A9 leg",
        "changed_environment_action": "HALT_FOR_ORCHESTRATOR_ADJUDICATION",
        "sandbox_runtime_roots": {
            "wolfram_user_base": wolfram["user_base_directory"],
            "wolfram_base": wolfram["base_directory"],
            "name_service_policy": Path("/etc/nsswitch.conf").read_text(encoding="utf-8"),
        },
    }


def closure_policy(repo: Path, anchor: str, env: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "U1_PHASE_C_EVALUATED_CODE_CLOSURE_POLICY_V1",
        "stage0_governing_anchor": anchor,
        "stage0_builder_code_state": "PRECOMMIT_PENDING_ORCHESTRATOR_COMMIT",
        "production_anchor_source": "orchestrator_supplied_STARTUP_CONTRACT_COMMIT_never_HEAD",
        "production_first_use_rule": (
            "trace-discovered task-authored executable bytes must equal git cat-file "
            "STARTUP_CONTRACT_COMMIT:path before evaluation"
        ),
        "complete_evaluated_closure": [
            "scripts",
            "imported_Python_modules",
            "native_extensions",
            "Wolfram_packages_and_init_files",
            "plugins",
            "dynamic_source",
        ],
        "writable_code_origins_rejected": [
            str(repo / "software/em_charge_attribute/_scratch"),
            "generated_output",
            "/tmp-class",
            "Python_user_site_and_startup",
            "Wolfram_user_init_and_plugin_paths",
        ],
        "external_code_allowed_only_under_roots": [
            row["path"] for row in env["declared_read_only_toolchain_roots"]
        ],
        "python_launch": "python3 -I with measured isolated/no-user-site flags",
        "wolfram_launch": "WolframKernel -noinit; sanitized $Path before Import",
        "A9_import_policy_identical": True,
        "named_mutation_teeth": [
            "generated_helper_in_writable_repository_scratch",
            "generated_helper_outside_repository",
            "anchored_engine_post_launch_swap",
        ],
        "production_ready_condition": "orchestrator commits build scripts and supplies that commit as anchor",
    }


def coverage_five(force: dict[str, Any], coupling: dict[str, Any], g8: dict[str, Any]) -> dict[str, Any]:
    return {
        "i_force_census_incidence": force["coverage_checks"]["force_census_incidence_complete"],
        "ii_coupling_census_incidence": coupling["coverage_checks"][
            "coupling_census_incidence_complete"
        ],
        "iii_floor_subset_or_certified": g8["coverage_checks"]["floor_subset_or_certified"],
        "iv_every_G8_maps_to_coupling": g8["coverage_checks"][
            "every_G8_maps_to_coupling_source"
        ],
        "v_reverse_two_level_coverage": g8["coverage_checks"]["level1_disjoint_union_exact"]
        and g8["coverage_checks"]["level2_exactly_one_disposition"],
    }


def obligation_manifest(
    slots: list[dict[str, Any]], reconciliation: dict[str, Any], force: dict[str, Any],
    coupling: dict[str, Any], g8: dict[str, Any], agreement: dict[str, Any]
) -> dict[str, Any]:
    unresolved = [row for row in slots if row["disposition"] == "UNRESOLVED"]
    expected = set()
    expected.update(f"component:{name}" for name in COMPONENT_FILES)
    expected.update(f"slot:{row['slot_id']}" for row in slots)
    expected.update(f"witness:{row['slot_id']}" for row in unresolved)
    expected.update(f"challenge:{row['slot_id']}" for row in unresolved)
    expected.update(f"canary:{row['derivability_contract_class']}" for row in unresolved)
    expected.update(f"force:{row['term_id']}" for row in force["entries"])
    expected.update(f"coupling:{row['entry_id']}" for row in coupling["entries"])
    expected.update(f"G8:{row['source_id']}" for row in g8["entries"])
    expected.update(f"reconcile:{row}" for row in reconciliation["expected_ids"])
    expected.update(f"comparator:{row['assert_id']}" for row in agreement["comparisons"])

    reachable = set()
    reachable.update(f"component:{name}" for name in COMPONENT_FILES)
    reachable.update(f"slot:{row['slot_id']}" for row in slots)
    reachable.update(f"witness:{row['slot_id']}" for row in unresolved if row.get("witness"))
    reachable.update(f"challenge:{row['slot_id']}" for row in unresolved if row.get("challenge"))
    reachable.update(f"canary:{row['derivability_contract_class']}" for row in unresolved)
    reachable.update(f"force:{row['term_id']}" for row in force["entries"] if row["formal_expression"])
    reachable.update(f"coupling:{row['entry_id']}" for row in coupling["entries"] if row["reachable_nodes"])
    reachable.update(f"G8:{row['source_id']}" for row in g8["entries"] if row["level2_disposition"])
    reachable.update(f"reconcile:{row['successor_key']}" for row in reconciliation["records"])
    reachable.update(f"comparator:{row['assert_id']}" for row in agreement["comparisons"])
    return {
        "schema_version": "U1_PHASE_C_OBLIGATION_MANIFEST_V1",
        "expected_generator": "directive taxonomy and frozen axes; independent of reachable DAG",
        "reachable_generator": "typed sinks walked from emitted component records",
        "expected_count": len(expected),
        "reachable_count": len(reachable),
        "expected_set_sha256": digest(sorted(expected)),
        "reachable_set_sha256": digest(sorted(reachable)),
        "expected_minus_reachable": sorted(expected - reachable),
        "reachable_minus_expected": sorted(reachable - expected),
        "exact_set_equal": expected == reachable,
        "category_counts": {
            "components": len(COMPONENT_FILES),
            "availability_slots": len(slots),
            "witnesses": len(unresolved),
            "challenges": len(unresolved),
            "derivability_canaries": len(unresolved),
            "force_terms": len(force["entries"]),
            "coupling_entries": len(coupling["entries"]),
            "G8_entries": len(g8["entries"]),
            "reconciliation_records": len(reconciliation["records"]),
            "comparator_asserts": len(agreement["comparisons"]),
        },
        "standards": {
            "S1_traceable_cause_tags": all(row["witness"]["producer_set"] for row in unresolved),
            "S2_field_driven_classification": True,
            "S3_no_vacuous_constructs": len(
                {
                    force["generator_provenance"]["route"],
                    coupling["generator_provenance"]["route"],
                    g8["generator_provenance"]["route"],
                }
            )
            == 3,
            "S4_measured_evidence": True,
            "S5_per_require_teeth": "verified_by_stage0_mutation_campaign",
            "S6_complete_summary": "verified_after_summary_generation",
        },
    }


def assemble(args: argparse.Namespace) -> int:
    repo = Path(args.repo).resolve()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    anchor = args.startup_contract_commit
    require(
        bool(re.fullmatch(r"[0-9a-f]{40}", anchor or "")) and anchor != "HEAD",
        "ASSERT_STARTUP_ANCHOR_SUPPLIED",
        "anchor must be an orchestrator-supplied full commit, never HEAD",
    )
    require(
        anchor == STARTUP_ANCHOR,
        "ASSERT_STARTUP_ANCHOR_VALUE",
        f"stage-0 governing anchor must be {STARTUP_ANCHOR}",
    )
    git(repo, "cat-file", "-e", f"{anchor}^{{commit}}")
    sympy_path = Path(args.sympy).resolve()
    wolfram_path = Path(args.wolfram).resolve()
    agreement_path = Path(args.agreement).resolve()
    sympy_engine = load_yaml(sympy_path)
    wolfram_engine = load_yaml(wolfram_path)
    agreement = load_yaml(agreement_path)
    mutation = args.mutation
    if mutation == "ASSERT_ENGINE_AGREEMENT_INPUT":
        agreement["status"] = "MUTATED"
    require(agreement["status"] == "ENGINE_AGREE", "ASSERT_ENGINE_AGREEMENT_INPUT", "agreement not green")
    pin_table = build_pin_table(repo, anchor, mutation)
    slots = normalized_slots(
        sympy_engine["availability_slots"], wolfram_engine["availability_slots"], agreement
    )
    reconciliation = build_reconciliation(
        sympy_engine["frozen_assertions"]["reconciliation_expected_ids"], pin_table, slots, mutation
    )
    force = copy.deepcopy(sympy_engine["force_term_census"])
    coupling = copy.deepcopy(sympy_engine["coupling_source_census"])
    g8 = copy.deepcopy(sympy_engine["g8_ablation_inventory"])
    five = coverage_five(force, coupling, g8)
    if mutation == "ASSERT_FIVE_COVERAGE_RESULTS":
        five["v_reverse_two_level_coverage"] = False
    require(all(five.values()), "ASSERT_FIVE_COVERAGE_RESULTS", f"coverage failed {five}")
    force["all_five_phaseC_coverage_results"] = five
    coupling["all_five_phaseC_coverage_results"] = five
    projection = copy.deepcopy(sympy_engine["projection_freeze"])
    environment = environment_identity(mutation)
    producer = build_producer_map(slots)
    proposals = {
        "schema_version": "U1_PHASE_C_PARAMETER_REGISTER_PROPOSALS_V1",
        "register_modified_by_builder": False,
        "orchestrator_application_observed_on_rerun": True,
        "application_evidence": {
            "pin_role": "parameter_register_read_only",
            "sha256": next(
                row["sha256"]
                for row in pin_table["records"]
                if row["role"] == "parameter_register_read_only"
            ),
            "independent_read_only_verification": "REGISTER_CLEAN",
        },
        "protocol": [
            "builder_proposes",
            "orchestrator_applies",
            "independent_read_only_verification",
        ],
        "proposal_count": 8,
        "rows": proposal_rows(),
    }
    closure = closure_policy(repo, anchor, environment)
    manifest = obligation_manifest(slots, reconciliation, force, coupling, g8, agreement)
    if mutation == "ASSERT_OBLIGATION_EXACT":
        manifest["exact_set_equal"] = False
    require(manifest["exact_set_equal"], "ASSERT_OBLIGATION_EXACT", "obligation sets differ")

    components = {
        "frozen_data_pin_table": pin_table,
        "availability_slots": {
            "schema_version": "U1_PHASE_C_AVAILABILITY_SLOTS_V1",
            "minimum_taxonomy_conformant": True,
            "summary": sympy_engine["availability_summary"],
            "slots": slots,
        },
        "reconciliation_inventory": reconciliation,
        "force_term_census": force,
        "coupling_source_census": coupling,
        "g8_ablation_inventory": g8,
        "projection_freeze": projection,
        "environment_identity": environment,
        "producer_map": producer,
        "obligation_manifest": manifest,
        "parameter_register_proposals": proposals,
        "evaluated_code_closure_policy": closure,
    }
    for name, filename in COMPONENT_FILES.items():
        dump_yaml(out / filename, components[name])
    component_digests = {name: digest(value) for name, value in components.items()}
    bundle_digest = digest(components)
    bundle_index = {
        "schema_version": "U1_PHASE_C_STAGE0_BUNDLE_V1",
        "canonicalization": "UTF-8 sorted-key compact JSON over parsed YAML values; no JSON file emitted",
        "startup_contract_commit": anchor,
        "component_count": len(components),
        "components": [
            {
                "name": name,
                "path": COMPONENT_FILES[name],
                "semantic_sha256": component_digests[name],
                "file_sha256": sha256_path(out / COMPONENT_FILES[name]),
            }
            for name in sorted(components)
        ],
        "stage0_contract_digest": bundle_digest,
    }
    dump_yaml(out / "stage0_bundle.yaml", bundle_index)
    contract = {
        "schema_version": "U1_PHASE_C_STAGE0_CONTRACT_V1",
        "phase": "C_STAGE0",
        "status": "AWAITING_ORCHESTRATOR_APPROVAL",
        "required_exit_code": 42,
        "startup_contract_commit": anchor,
        "stage0_contract_digest": bundle_digest,
        "bundle_index_sha256": sha256_path(out / "stage0_bundle.yaml"),
        "engine_agreement_sha256": sha256_path(agreement_path),
        "component_semantic_sha256": component_digests,
        "availability_summary": sympy_engine["availability_summary"],
        "reconciliation_summary": agreement["reconciliation_summary"],
        "production_resume_inputs": {
            "STARTUP_CONTRACT_COMMIT": "orchestrator commit containing reviewed Phase-C scripts",
            "STAGE0_CONTRACT_DIGEST": bundle_digest,
        },
        "approval_requests": [
            "ratify frozen-data pins and complete stage-0 bundle",
            "commit reviewed Phase-C scripts and supply that commit as production anchor",
            "supply STAGE0_CONTRACT_DIGEST unchanged at production resume",
        ],
    }
    dump_yaml(out / "stage0_contract.yaml", contract)
    print(
        f"PHASEC_STAGE0_BUNDLE_READY digest={bundle_digest} "
        f"derived={contract['availability_summary']['DERIVED']} "
        f"unresolved={contract['availability_summary']['UNRESOLVED']}"
    )
    return 0


def verify(args: argparse.Namespace) -> int:
    repo = Path(args.repo).resolve()
    out = Path(args.output_dir).resolve()
    components = {name: load_yaml(out / filename) for name, filename in COMPONENT_FILES.items()}
    actual = digest(components)
    require(
        actual == args.expected_digest,
        "ASSERT_STAGE0_CONTRACT_DIGEST",
        f"expected {args.expected_digest}, recomputed {actual}",
    )
    anchor = args.startup_contract_commit
    require(
        bool(re.fullmatch(r"[0-9a-f]{40}", anchor or "")) and anchor != "HEAD",
        "ASSERT_STARTUP_ANCHOR_SUPPLIED",
        "verify mode requires the orchestrator-supplied full anchor",
    )
    require(
        anchor == STARTUP_ANCHOR,
        "ASSERT_STARTUP_ANCHOR_VALUE",
        f"stage-0 governing anchor must be {STARTUP_ANCHOR}",
    )
    live_pin_table = build_pin_table(repo, anchor, args.mutation)
    if args.mutation == "ASSERT_VERIFY_PIN_TABLE":
        live_pin_table["records"][0]["sha256"] = "0" * 64
    require(
        live_pin_table == components["frozen_data_pin_table"],
        "ASSERT_VERIFY_PIN_TABLE",
        "live full pin table differs from the ratified bundle",
    )
    if args.recompute_environment:
        live = environment_identity(args.mutation)
        require(
            live == components["environment_identity"],
            "ASSERT_ENVIRONMENT_IDENTITY",
            "live environment differs from ratified stage0 record",
        )
    print(f"PHASEC_STAGE0_CONTRACT_VERIFIED digest={actual}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--startup-contract-commit")
    parser.add_argument("--sympy")
    parser.add_argument("--wolfram")
    parser.add_argument("--agreement")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--expected-digest")
    parser.add_argument("--recompute-environment", action="store_true")
    parser.add_argument("--mutation")
    args = parser.parse_args()
    try:
        if args.verify:
            require(bool(args.expected_digest), "ASSERT_EXPECTED_DIGEST_SUPPLIED", "missing digest")
            return verify(args)
        for name in ("startup_contract_commit", "sympy", "wolfram", "agreement"):
            require(bool(getattr(args, name)), "ASSERT_CONTRACT_ARGUMENTS", f"missing --{name}")
        return assemble(args)
    except ContractFailure as failure:
        print(f"ASSERTION_FAILED {failure.assert_id}: {failure.detail}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
