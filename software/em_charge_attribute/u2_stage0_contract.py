#!/usr/bin/env python3
"""Assemble and verify the canonically hashed U2 stage-0 bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any

import sympy
import yaml


DIRECTIVE_PROVENANCE_COMMIT = "aca0623f5e90e55fbdb0947f9ea722a52b0839a0"
V11_PHYSICS_RECORD_SET_DIGEST = "d70bf6bee368aa4383defa3a904af71eb083fe1eea34ab02444af73243b4de7a"
V11_PRODUCTION_RESULTS_SHA256 = "1832c490ede0c14a94963543a7db12e6e3e2a4fcc3c4f65f958436072f92f56e"
PHASEC_ANCHOR = "21ec07f7a4f2814fa7ca642ded7ebd353671578b"
PHASEC_STAGE0_DIGEST = "e632a8d6729d0a1b3a4ade883c28f6b21f7a29fea566318cdd6fefec8c15d0da"
B2_CONTRACT_DIGEST = "8b8ee113d9a8342ac11f23e7959e05c323cc89a81b7c298af85f646b6e42a9e7"
PHASE_A_LINEAGE = "b23993cca80dc3e6a790abcf68c1af63aa804fc47b06b153b9f224ccf27f899d"

COMPONENT_FILES = {
    "frozen_data_pin_table": "frozen_data_pin_table.yaml",
    "candidate_inventory": "candidate_inventory.yaml",
    "obligation_censuses": "obligation_censuses.yaml",
    "dependency_grid_inventory": "dependency_grid_inventory.yaml",
    "vocabulary_freeze": "vocabulary_freeze.yaml",
    "evidence_taxonomy": "evidence_taxonomy.yaml",
    "availability_slots": "availability_slots.yaml",
    "route_fixture_inventory": "route_fixture_inventory.yaml",
    "closure_template_contracts": "closure_template_contracts.yaml",
    "environment_identity": "environment_identity.yaml",
    "standard_bindings": "standard_bindings.yaml",
    "producer_map": "producer_map.yaml",
    "evaluated_code_closure_policy": "evaluated_code_closure_policy.yaml",
    "parameter_register_proposals": "parameter_register_proposals.yaml",
    "obligation_manifest": "obligation_manifest.yaml",
}


class ContractFailure(RuntimeError):
    def __init__(self, assert_id: str, detail: str):
        super().__init__(f"{assert_id}: {detail}")
        self.assert_id = assert_id
        self.detail = detail


def require(condition: bool, assert_id: str, detail: str) -> None:
    if not condition:
        raise ContractFailure(assert_id, detail)


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=140), encoding="utf-8")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(repo: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *args], cwd=repo, check=check, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True,
    )


def pin_specs() -> list[tuple[str, str, str | None]]:
    phasec0 = "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_0_tilt_coupling_contract"
    phasec1 = "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_1_tilt_coupling_production"
    b2_0 = "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_0_intake_radiative_contract"
    b2_1 = "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production"
    directive_lineage = DIRECTIVE_PROVENANCE_COMMIT
    return [
        ("governing_U2_directive", "software/em_charge_attribute/directive_u2_boundary_adjudication.md", directive_lineage),
        ("program_handoff_sixteen_standards", "docs/em_analog_next_phase_handoff.md", directive_lineage),
        ("governing_spec_v3_1", "docs/em_u1_body_definition.md", directive_lineage),
        ("phaseC_normative_directive_v7", "software/em_charge_attribute/directive_u1_phaseC_tilt_coupling.md", directive_lineage),
        ("B2_witness_contract_v48", "software/em_charge_attribute/directive_u1_phaseB2_intake_radiative.md", directive_lineage),
        ("parent_U1_dynamics_directive", "software/em_charge_attribute/directive_u1_body_dynamics.md", directive_lineage),
        ("decision_16", "software/stage1_solver/decisions/16_retire_brane_polar_field.md", directive_lineage),
        ("decision_17", "software/stage1_solver/decisions/17_trust_apparatus_trim.md", directive_lineage),
        ("ensemble_topology_rubric", "docs/em_phaseC_force_decomposition.md", directive_lineage),
        ("frozen_7_0_declarations", "software/em_charge_attribute/u1_body_dynamics_inputs.yaml", directive_lineage),
        ("phase_A_amendment", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage1/phase_a_amendment_agreement.yaml", directive_lineage),
        ("phase_A_sympy_payload", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_a.json", directive_lineage),
        ("phase_A_wolfram_payload", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_a.json", directive_lineage),
        ("B1_approval", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_orchestrator_approval.yaml", directive_lineage),
        ("B1_final_results", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_final_results_snapshot.yaml", directive_lineage),
        ("B1_final_report", "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_final_report_snapshot.md", directive_lineage),
        ("B1_mechanics_inputs", "software/em_charge_attribute/u1_body_mechanics_inputs.yaml", directive_lineage),
        ("B2_sealed_stage0_contract", f"{b2_0}/stage0_contract.yaml", None),
        ("B2_stage0_engine_agreement", f"{b2_0}/stage0_engine_agreement.yaml", None),
        ("B2_sympy_production", f"{b2_1}/sympy_b2.yaml", None),
        ("B2_wolfram_production", f"{b2_1}/mathematica_b2.yaml", None),
        ("B2_engine_agreement", f"{b2_1}/engine_agreement.yaml", None),
        ("B2_completion_gate", f"{b2_1}/completion_gate.yaml", None),
        ("B2_production_aggregate", "software/em_charge_attribute/reports/u1_body_dynamics_results.yaml", directive_lineage),
        ("phaseC_stage0_contract", f"{phasec0}/stage0_contract.yaml", None),
        ("phaseC_stage0_bundle", f"{phasec0}/stage0_bundle.yaml", None),
        ("phaseC_stage0_approval_halt", f"{phasec0}/orchestrator_approval_halt.yaml", None),
        ("phaseC_availability_slots", f"{phasec0}/availability_slots.yaml", None),
        ("phaseC_producer_map", f"{phasec0}/producer_map.yaml", None),
        ("phaseC_production_results", f"{phasec1}/production_results.yaml", None),
        ("phaseC_production_agreement", f"{phasec1}/production_engine_agreement.yaml", None),
        ("phaseC_production_terminal", f"{phasec1}/production_terminal.yaml", None),
        (
            "U2_v11_production_physics_record_reference",
            "software/em_charge_attribute/reports/u2_boundary_adjudication_artifacts/stage_1_production/production_results.yaml",
            None,
        ),
        ("pathA29_one_sided_report", "software/stage1_solver/reports/pathA_29_brane_bulk_return.md", directive_lineage),
        ("pathA29_one_sided_results", "software/stage1_solver/reports/pathA_29_results.yaml", directive_lineage),
        ("pathA38_h_report", "software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md", directive_lineage),
        ("pathA38_h_results", "software/stage1_solver/reports/pathA_38_results.yaml", directive_lineage),
        ("pathA39_stage3_report", "software/stage1_solver/reports/pathA_39_stage3_operator_parity.md", directive_lineage),
        ("pathA39_stage3_results", "software/stage1_solver/reports/pathA_39_stage3_operator_parity_results.yaml", directive_lineage),
        ("native_action_stage004", "research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md", directive_lineage),
        ("native_action_stage006", "research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md", directive_lineage),
        ("parameter_register_read_only", "research/pde_ledger_v2/notes/parameter_register.md", directive_lineage),
    ]


def build_pin_table(repo: Path, mutation: str | None = None) -> dict[str, Any]:
    specs = pin_specs()
    if mutation == "TOOTH_PIN_INPUT_EXISTS":
        specs[0] = (specs[0][0], "software/em_charge_attribute/missing_U2_pin.yaml", specs[0][2])
    records = []
    for role, rel_path, lineage in specs:
        path = repo / rel_path
        require(path.is_file(), "ASSERT_PIN_INPUT_EXISTS", f"missing frozen input {rel_path}")
        disk = path.read_bytes()
        anchor_blob = None
        if lineage:
            result = subprocess.run(
                ["git", "cat-file", "blob", f"{lineage}:{rel_path}"], cwd=repo,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            require(result.returncode == 0, "ASSERT_PIN_LINEAGE", f"{rel_path} absent at {lineage}")
            anchor_blob = result.stdout
            if mutation == "TOOTH_PIN_LINEAGE" and role == "governing_U2_directive":
                anchor_blob += b"\nU2_MUTATED_LINEAGE\n"
            require(hashlib.sha256(disk).digest() == hashlib.sha256(anchor_blob).digest(), "ASSERT_PIN_LINEAGE", f"{rel_path} differs from frozen lineage {lineage}")
        records.append({
            "role": role, "path": rel_path, "sha256": hashlib.sha256(disk).hexdigest(),
            "size": len(disk), "lineage_commit": lineage,
            "lineage_blob_sha256": hashlib.sha256(anchor_blob).hexdigest() if anchor_blob is not None else None,
            "read_only": True,
        })
    by_role = {row["role"]: row for row in records}
    phase_a = load_yaml(repo / by_role["phase_A_amendment"]["path"])
    b1 = load_yaml(repo / by_role["B1_final_results"]["path"])
    phasec0 = load_yaml(repo / by_role["phaseC_stage0_contract"]["path"])
    phasec1 = load_yaml(repo / by_role["phaseC_production_terminal"]["path"])
    if mutation == "TOOTH_B2_CONTRACT_PIN":
        by_role["B2_sealed_stage0_contract"]["sha256"] = "0" * 64
    require(by_role["B2_sealed_stage0_contract"]["sha256"] == B2_CONTRACT_DIGEST, "ASSERT_B2_CONTRACT_PIN", "sealed B2 contract changed")
    if mutation == "TOOTH_PHASE_A_LINEAGE":
        phase_a["digest_gate"]["sympy"] = "0" * 64
    require(phase_a["digest_gate"]["agreement"] and phase_a["digest_gate"]["sympy"] == phase_a["digest_gate"]["wolfram"] == PHASE_A_LINEAGE, "ASSERT_PHASE_A_LINEAGE", "Phase-A amendment lineage changed")
    if mutation == "TOOTH_B1_TERMINAL_PIN":
        b1["mechanics_partition_ledger"]["state"] = "MUTATED"
    require(b1["mechanics_partition_ledger"]["state"] == "partition_open_pending_B2", "ASSERT_B1_TERMINAL_PIN", "B1 terminal changed")
    if mutation == "TOOTH_PHASEC_DIGEST_PIN":
        phasec0["stage0_contract_digest"] = "0" * 64
    require(phasec0["stage0_contract_digest"] == PHASEC_STAGE0_DIGEST, "ASSERT_PHASEC_DIGEST_PIN", "Phase-C ratified digest changed")
    if mutation == "TOOTH_PHASEC_TERMINAL_PIN":
        phasec1["status"] = "MUTATED"
    require(phasec1["status"] == "SUCCESS_STOP" and phasec1["exit_code"] == 0 and phasec1["stage0_contract_digest"] == PHASEC_STAGE0_DIGEST, "ASSERT_PHASEC_TERMINAL_PIN", "Phase-C production terminal changed")
    if mutation == "TOOTH_V11_PHYSICS_RECORD_REFERENCE":
        by_role["U2_v11_production_physics_record_reference"]["sha256"] = "0" * 64
    require(
        by_role["U2_v11_production_physics_record_reference"]["sha256"] == V11_PRODUCTION_RESULTS_SHA256,
        "ASSERT_V11_PHYSICS_RECORD_REFERENCE", "v11 U2 production source changed",
    )
    return {
        "schema_version": "U2_FROZEN_DATA_PIN_TABLE_V1", "record_count": len(records),
        "records": records,
        "frozen_assertions": {
            "Phase_A_payload_lineage": PHASE_A_LINEAGE,
            "B1_commit": "ef934360b031bef54b37ca96d5c73cb10b0d15fd",
            "B1_partition_terminal": "partition_open_pending_B2",
            "B2_commit": "f9d7bc30d5ce46233bc0de6c97bd6182a1c6e186",
            "B2_contract_digest": B2_CONTRACT_DIGEST,
            "PhaseC_anchor": PHASEC_ANCHOR,
            "PhaseC_stage0_digest": PHASEC_STAGE0_DIGEST,
            "PhaseC_production_status": "SUCCESS_STOP",
            "U2_v11_physics_record_set_digest": V11_PHYSICS_RECORD_SET_DIGEST,
            "U2_v11_production_results_sha256": V11_PRODUCTION_RESULTS_SHA256,
            "U2_v11_ratified_stage0_digest": "9eff1b0c49e89007aea1008cb6712b0ea495168d101ce43ddce1cffaf68749c4",
            "U2_v11_production_anchor": "53529bf1729811f5ae9faa429cf836507469569b",
            "U2_v11_wrap_anchor": "5ceebb24",
        },
        "mutation_policy": "pinned data changes are amendments, never repins",
    }


def mount_record(path: Path, mutation: str | None = None) -> dict[str, Any]:
    resolved = path.resolve()
    best: tuple[int, list[str]] | None = None
    for line in Path("/proc/self/mountinfo").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        try:
            resolved.relative_to(Path(fields[4]))
        except ValueError:
            continue
        if best is None or len(fields[4]) > best[0]:
            best = (len(fields[4]), fields)
    if mutation == "TOOTH_ENV_MOUNT_IDENTITY":
        best = None
    require(best is not None, "ASSERT_ENV_MOUNT_IDENTITY", f"no mount record for {resolved}")
    fields = best[1]
    separator = fields.index("-")
    return {
        "path": str(resolved), "mountpoint": fields[4],
        "mount_options": sorted(fields[5].split(",")), "filesystem": fields[separator + 1],
        "source": fields[separator + 2], "super_options": sorted(fields[separator + 3].split(",")),
    }


def tree_identity(roots: list[Path]) -> dict[str, Any]:
    records = []
    for root in roots:
        for path in root.rglob("*"):
            if path.is_file() and not path.is_symlink():
                records.append({"path": str(path), "size": path.stat().st_size, "sha256": sha256_path(path)})
    records.sort(key=lambda row: row["path"])
    return {"file_count": len(records), "digest": digest(records)}


def wolfram_identity(mutation: str | None = None) -> dict[str, Any]:
    launcher = Path("/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel")
    code = (
        'Print["VERSION="<>$Version];Print["SYSTEM="<>$SystemID];'
        'Print["INSTALL="<>$InstallationDirectory];Print["USERBASE="<>$UserBaseDirectory];'
        'Print["BASE="<>$BaseDirectory];Print["SEATS="<>ToString[$MaxLicenseProcesses]];Quit[]'
    )
    result = subprocess.run(
        [str(launcher), "-noinit", "-noprompt", "-run", code], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=60,
    )
    parsed = {}
    for raw in result.stdout.splitlines():
        line = raw.strip().strip('"')
        if "=" in line:
            key, value = line.split("=", 1); parsed[key] = value
    if mutation == "TOOTH_WOLFRAM_SEAT_LIMIT":
        parsed["SEATS"] = "3"
    require(parsed.get("SEATS") == "2", "ASSERT_WOLFRAM_SEAT_LIMIT", "Wolfram seat limit changed")
    kernel = Path("/usr/local/Wolfram/Wolfram/15.0/SystemFiles/Kernel/Binaries/Linux-x86-64/WolframKernel")
    return {
        "version": parsed.get("VERSION"), "system_id": parsed.get("SYSTEM"),
        "installation_directory": parsed.get("INSTALL"), "user_base_directory": parsed.get("USERBASE"),
        "base_directory": parsed.get("BASE"), "max_license_processes": int(parsed["SEATS"]),
        "launcher_path": str(launcher), "launcher_sha256": sha256_path(launcher),
        "kernel_path": str(kernel), "kernel_sha256": sha256_path(kernel),
        "launch_mode": "WolframKernel -noinit -noprompt", "user_init_disabled": True,
        "user_plugin_paths_removed_before_import": True,
    }


def environment_identity(mutation: str | None = None) -> dict[str, Any]:
    sanitized = sys.flags.isolated == 1 and sys.flags.no_user_site == 1
    if mutation == "TOOTH_ENV_PYTHON_SANITIZED":
        sanitized = False
    require(sanitized, "ASSERT_ENV_PYTHON_SANITIZED", "Python not isolated/no-user-site")
    python = Path(sys.executable).resolve()
    sympy_root = Path(sympy.__file__).resolve().parent
    yaml_root = Path(yaml.__file__).resolve().parent
    stdlib = Path(f"/usr/lib/python{sys.version_info.major}.{sys.version_info.minor}")
    wolfram = wolfram_identity(mutation)
    roots = [Path("/usr"), Path(wolfram["installation_directory"]), Path("/opt/Wolfram")]
    root_records = []
    for root in roots:
        info = root.stat()
        root_records.append({
            "path": str(root.resolve()), "device": info.st_dev, "inode": info.st_ino,
            "mode": stat.S_IMODE(info.st_mode), "uid": info.st_uid, "gid": info.st_gid,
            "mount": mount_record(root, mutation),
        })
    read_only = all("ro" in row["mount"]["mount_options"] for row in root_records)
    if mutation == "TOOTH_ENV_TOOLCHAIN_READ_ONLY":
        read_only = False
    require(read_only, "ASSERT_ENV_TOOLCHAIN_READ_ONLY", "declared toolchain roots not read-only")
    evaluated = tree_identity([stdlib, sympy_root, yaml_root])
    selected = {
        "python_executable_sha256": sha256_path(python), "python_evaluated_tree": evaluated,
        "wolfram_launcher_sha256": wolfram["launcher_sha256"], "wolfram_kernel_sha256": wolfram["kernel_sha256"],
        "root_records": root_records, "nsswitch_sha256": sha256_path(Path("/etc/nsswitch.conf")),
    }
    return {
        "schema_version": "U2_ENVIRONMENT_IDENTITY_V1",
        "python": {
            "version": platform.python_version(), "implementation": platform.python_implementation(),
            "executable": str(python), "executable_sha256": sha256_path(python),
            "sympy_version": sympy.__version__, "sympy_root": str(sympy_root),
            "pyyaml_version": yaml.__version__, "pyyaml_root": str(yaml_root),
            "isolated_flag": sys.flags.isolated, "no_user_site_flag": sys.flags.no_user_site,
            "user_site_enabled": False,
        },
        "wolfram": wolfram, "declared_read_only_toolchain_roots": root_records,
        "selected_evaluated_toolchain_identity": selected,
        "toolchain_root_identity_digest": digest(selected),
        "changed_environment_action": "HALT_FOR_ORCHESTRATOR_ADJUDICATION",
        "production_assertion": "recompute before every evaluation and A9 leg",
    }


def closure_policy(repo: Path, environment: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "U2_EVALUATED_CODE_CLOSURE_POLICY_V1",
        "stage0_governing_anchor_source": "orchestrator_supplied_STARTUP_CONTRACT_COMMIT_never_HEAD",
        "concrete_anchor_recorded_in": [
            "stage0_bundle.yaml", "stage0_contract.yaml", "evaluated_code_closure_evidence.yaml",
            "runner_state.yaml", "stage_run_records.yaml",
        ],
        "stage0_builder_code_state": "ORCHESTRATOR_SUPPLIED_CODE_STACK_COMMIT",
        "production_anchor_source": "orchestrator_supplied_STARTUP_CONTRACT_COMMIT_never_HEAD",
        "complete_evaluated_closure": [
            "scripts", "imported_Python_modules", "native_extensions", "Wolfram_packages_and_init_files",
            "plugins", "dynamic_source",
        ],
        "first_use_rule": "trace-discovered task code bytes equal the launch-supplied anchor blob during stage0 and production",
        "stage0_precommit_exception": "inactive in sealed stage0; authoring-only anchor-absent direct-entrypoint mode retained by the guard",
        "external_code_allowed_only_under_roots": [row["path"] for row in environment["declared_read_only_toolchain_roots"]],
        "writable_code_origins_rejected": [str(repo / "software/em_charge_attribute/_scratch"), "/tmp", "user-site", "Wolfram-init"],
        "dual_placement_helper_teeth": ["generated_helper_repository_scratch", "generated_helper_external_tmp"],
        "production_ready_condition": "orchestrator supplies the commit containing the exact production code stack",
        "script_manifest_present": False, "reseal_ceremony_present": False,
    }


def producer_map(slots: list[dict[str, Any]]) -> dict[str, Any]:
    relations = [{
        "slot_id": row["slot_id"], "producer_set": row["producer_set"],
        "required_type": row["required_type"], "required_dimensions": row["required_dimensions"],
        "acceptance_predicate": row["acceptance_predicate"],
        "generation_rule": "U2_minimum_slot_taxonomy_x_external_source_census",
    } for row in slots]
    return {
        "schema_version": "U2_PRODUCER_MAP_V1", "builder_declared_producer_tags_accepted": False,
        "generation": "directive taxonomy x frozen typed ancestry", "relation_count": len(relations),
        "relations": relations,
    }


def expected_obligation_floor(
    components_without_manifest: dict[str, Any], agreement: dict[str, Any]
) -> set[str]:
    candidate = components_without_manifest["candidate_inventory"]
    obligations = components_without_manifest["obligation_censuses"]["censuses"]
    grid = components_without_manifest["dependency_grid_inventory"]
    candidate_ids = [row["candidate_id"] for row in candidate["candidate_records"]]
    mixtures = [row["family_id"] for row in candidate["mixture_generator_A"]]
    strata = grid["active_strata"]
    ambients = ("one_sided_pathA29", "two_sided_R_w_postulate")
    expected: set[str] = set()
    expected.update(f"component:{name}" for name in COMPONENT_FILES)
    expected.update(f"candidate:{value}" for value in candidate_ids)
    expected.update(f"obligation:{candidate_id}:{obligation}" for candidate_id, row in obligations.items() for obligation in row["generator_A"])
    expected_cells = {
        f"candidate={candidate_id}|ambient={ambient}|stratum={stratum}"
        for candidate_id in candidate_ids for ambient in ambients for stratum in strata
    }
    expected.update(f"cell:{cell_id}" for cell_id in expected_cells)
    expected.update(
        f"promotion:ambient={ambient}|context={stratum}"
        for ambient in ambients for stratum in strata
    )
    expected.update(f"route:route:{cell_id.removeprefix('candidate=')}" for cell_id in expected_cells)
    expected_slot_ids = {f"candidate_definition:{candidate_id}" for candidate_id in candidate_ids}
    expected_slot_ids.update(f"mixture_law:{family_id}" for family_id in mixtures)
    expected_slot_ids.add("basis_closure")
    expected_slot_ids.update(
        f"topology:{candidate_id}:{question}"
        for candidate_id in candidate_ids
        for question in ("sector_disconnection", "finite_energy_interpolation", "pair_annihilation")
    )
    expected_slot_ids.update(f"ensemble:{candidate_id}:boundary_action_variation" for candidate_id in candidate_ids)
    expected_slot_ids.update(f"host_location:{candidate_id}" for candidate_id in candidate_ids)
    expected_slot_ids.update(f"mechanical_closure:{candidate_id}" for candidate_id in candidate_ids)
    expected_slot_ids.update(f"template_free_data:{stratum}" for stratum in strata)
    expected_slot_ids.update(f"template_cell_specific:{candidate_id}" for candidate_id in candidate_ids)
    expected.update(f"slot:{slot_id}" for slot_id in expected_slot_ids)
    expected_unresolved = {
        "basis_closure", "candidate_definition:OTHER", "ensemble:E3:boundary_action_variation",
        "ensemble:OTHER:boundary_action_variation", "template_cell_specific:OTHER",
    }
    expected_unresolved.update(slot_id for slot_id in expected_slot_ids if slot_id.startswith(("topology:", "host_location:", "mechanical_closure:", "template_free_data:")))
    expected.update(f"witness:{slot_id}" for slot_id in expected_unresolved)
    expected.update(f"challenge:{slot_id}" for slot_id in expected_unresolved)
    expected.update(f"standard:S-{index}" for index in range(1, 7))
    expected.update(f"standard:P-{index}" for index in range(1, 11))
    expected.update(f"check:{row['expected_assert_id']}" for row in agreement["mutation_catalog"])
    return expected


def reachable_obligation_sinks(
    components_without_manifest: dict[str, Any], agreement: dict[str, Any]
) -> set[str]:
    """Walk the emitted typed sinks and their immediate ancestry relations."""
    reachable = {f"component:{name}" for name in components_without_manifest}
    reachable.add("component:obligation_manifest")
    candidate = components_without_manifest["candidate_inventory"]
    reachable.update(f"candidate:{candidate_id}" for candidate_id in candidate["candidate_axis"])
    obligations = components_without_manifest["obligation_censuses"]["censuses"]
    reachable.update(
        f"obligation:{candidate_id}:{obligation_id}"
        for candidate_id, row in obligations.items() for obligation_id in row["generator_B"]
    )
    grid = components_without_manifest["dependency_grid_inventory"]
    reachable.update(f"cell:{row['cell_id']}" for row in grid["grid_cells"])
    reachable.update(f"promotion:{row['promotion_key']}" for row in grid["promotion_contexts"])
    routes = components_without_manifest["route_fixture_inventory"]["route_records"]
    reachable.update(f"route:{row['route_id']}" for row in routes)
    slots = components_without_manifest["availability_slots"]["slots"]
    reachable.update(f"slot:{row['slot_id']}" for row in slots)
    reachable.update(
        f"witness:{row['witness']['datum_id']}" for row in slots if "witness" in row
    )
    reachable.update(
        f"challenge:{row['challenge_id'].removeprefix('challenge:')}"
        for row in slots if "challenge_id" in row
    )
    standards = components_without_manifest["standard_bindings"]["bindings"]
    reachable.update(f"standard:{row['standard_id']}" for row in standards)
    reachable.update(f"check:{row['assert_id']}" for row in agreement["checks"] if row.get("status") == "PASS")
    return reachable


def obligation_manifest(components_without_manifest: dict[str, Any], agreement: dict[str, Any]) -> dict[str, Any]:
    expected = expected_obligation_floor(components_without_manifest, agreement)
    reachable = reachable_obligation_sinks(components_without_manifest, agreement)
    candidate = components_without_manifest["candidate_inventory"]
    grid = components_without_manifest["dependency_grid_inventory"]
    slots = components_without_manifest["availability_slots"]["slots"]
    routes = components_without_manifest["route_fixture_inventory"]["route_records"]
    standards = components_without_manifest["standard_bindings"]["bindings"]
    return {
        "schema_version": "U2_OBLIGATION_MANIFEST_V1",
        "expected_generator": "independent directive-floor expansion over candidate records and frozen axes",
        "reachable_generator": "executed structural walk over emitted component typed sinks and ancestry IDs",
        "expected_count": len(expected), "reachable_count": len(reachable),
        "expected_set_sha256": digest(sorted(expected)), "reachable_set_sha256": digest(sorted(reachable)),
        "expected_minus_reachable": sorted(expected - reachable), "reachable_minus_expected": sorted(reachable - expected),
        "exact_set_equal": expected == reachable,
        "category_counts": {
            "components": len(COMPONENT_FILES), "candidates": len(candidate["candidate_axis"]),
            "grid_cells": len(grid["grid_cells"]), "promotion_contexts": len(grid["promotion_contexts"]),
            "routes": len(routes), "slots": len(slots),
            "witness_challenges": 2 * sum("witness" in row for row in slots),
            "standards": len(standards), "comparator_checks": len(agreement["checks"]),
        },
    }


def contract_mutation_catalog() -> list[dict[str, str]]:
    base = [
        ("TOOTH_STARTUP_ANCHOR_SUPPLIED", "ASSERT_STARTUP_ANCHOR_SUPPLIED"),
        ("TOOTH_STARTUP_ANCHOR_NEVER_HEAD", "ASSERT_STARTUP_ANCHOR_NEVER_HEAD"),
        ("TOOTH_GIT_ANCHOR_RESOLUTION", "ASSERT_GIT_ANCHOR_RESOLUTION"),
        ("TOOTH_PIN_INPUT_EXISTS", "ASSERT_PIN_INPUT_EXISTS"),
        ("TOOTH_PIN_LINEAGE", "ASSERT_PIN_LINEAGE"),
        ("TOOTH_B2_CONTRACT_PIN", "ASSERT_B2_CONTRACT_PIN"),
        ("TOOTH_PHASE_A_LINEAGE", "ASSERT_PHASE_A_LINEAGE"),
        ("TOOTH_B1_TERMINAL_PIN", "ASSERT_B1_TERMINAL_PIN"),
        ("TOOTH_PHASEC_DIGEST_PIN", "ASSERT_PHASEC_DIGEST_PIN"),
        ("TOOTH_PHASEC_TERMINAL_PIN", "ASSERT_PHASEC_TERMINAL_PIN"),
        ("TOOTH_V11_PHYSICS_RECORD_REFERENCE", "ASSERT_V11_PHYSICS_RECORD_REFERENCE"),
        ("TOOTH_ENGINE_AGREEMENT_INPUT", "ASSERT_ENGINE_AGREEMENT_INPUT"),
        ("TOOTH_PRODUCER_MAP", "ASSERT_PRODUCER_MAP"),
        ("TOOTH_OBLIGATION_MANIFEST", "ASSERT_OBLIGATION_MANIFEST"),
        ("TOOTH_ENV_PYTHON_SANITIZED", "ASSERT_ENV_PYTHON_SANITIZED"),
        ("TOOTH_WOLFRAM_SEAT_LIMIT", "ASSERT_WOLFRAM_SEAT_LIMIT"),
        ("TOOTH_ENV_TOOLCHAIN_READ_ONLY", "ASSERT_ENV_TOOLCHAIN_READ_ONLY"),
        ("TOOTH_ENV_MOUNT_IDENTITY", "ASSERT_ENV_MOUNT_IDENTITY"),
        ("TOOTH_ENVIRONMENT_IDENTITY", "ASSERT_ENVIRONMENT_IDENTITY"),
        ("TOOTH_PARAMETER_REGISTER", "ASSERT_PARAMETER_REGISTER_UNCHANGED"),
        ("TOOTH_EXPECTED_DIGEST", "ASSERT_EXPECTED_DIGEST_SUPPLIED"),
        ("TOOTH_STAGE0_CONTRACT_DIGEST", "ASSERT_STAGE0_CONTRACT_DIGEST"),
        ("TOOTH_RECORDED_STARTUP_ANCHOR", "ASSERT_RECORDED_STARTUP_ANCHOR"),
        ("TOOTH_VERIFY_PIN_TABLE", "ASSERT_VERIFY_PIN_TABLE"),
        ("TOOTH_APPROVAL_EXIT", "ASSERT_APPROVAL_EXIT_42"),
        ("TOOTH_CONTRACT_ZERO_INTEGRITY", "ASSERT_CONTRACT_ZERO_INTEGRITY"),
        ("TOOTH_CONTRACT_ARGUMENTS", "ASSERT_CONTRACT_ARGUMENTS"),
    ]
    base.extend((f"BUNDLE_COMPONENT:{name}", "ASSERT_STAGE0_CONTRACT_DIGEST") for name in COMPONENT_FILES)
    return [{"mutation_id": mutation, "expected_assert_id": assert_id} for mutation, assert_id in base]


def validate_anchor(repo: Path, anchor: str, mutation: str | None = None) -> None:
    value = "HEAD" if mutation == "TOOTH_STARTUP_ANCHOR_NEVER_HEAD" else anchor
    require(value != "HEAD", "ASSERT_STARTUP_ANCHOR_NEVER_HEAD", "HEAD is not an allowed startup anchor")
    require(bool(re.fullmatch(r"[0-9a-f]{40}", value or "")), "ASSERT_STARTUP_ANCHOR_SUPPLIED", "full orchestrator-supplied anchor required")
    resolved = git(repo, "rev-parse", f"{anchor}^{{commit}}").stdout.strip()
    if mutation == "TOOTH_GIT_ANCHOR_RESOLUTION":
        resolved = "0" * 40
    require(resolved == anchor, "ASSERT_GIT_ANCHOR_RESOLUTION", "anchor did not resolve identically")


def assemble(args: argparse.Namespace) -> int:
    repo = Path(args.repo).resolve(); out = Path(args.output_dir).resolve(); out.mkdir(parents=True, exist_ok=True)
    validate_anchor(repo, args.startup_contract_commit, args.mutation)
    sympy_engine = load_yaml(Path(args.sympy)); agreement = load_yaml(Path(args.agreement))
    if args.mutation == "TOOTH_ENGINE_AGREEMENT_INPUT": agreement["status"] = "MUTATED"
    require(agreement["status"] == "ENGINE_AGREE", "ASSERT_ENGINE_AGREEMENT_INPUT", "engine agreement not green")
    semantic = sympy_engine["semantic_view"]
    pin_table = build_pin_table(repo, args.mutation)
    environment = environment_identity(args.mutation)
    slots = semantic["availability_slots"]
    producer = producer_map(slots)
    if args.mutation == "TOOTH_PRODUCER_MAP": producer["relations"].pop()
    require(producer["relation_count"] == len(producer["relations"]) == len(slots), "ASSERT_PRODUCER_MAP", "producer map does not cover slot floor")
    parameter_record = next(row for row in pin_table["records"] if row["role"] == "parameter_register_read_only")
    proposals = {
        "schema_version": "U2_PARAMETER_REGISTER_PROPOSALS_V1", "register_modified_by_builder": False,
        "proposal_count": 0, "rows": [], "read_only_register_sha256": parameter_record["sha256"],
        "protocol": ["builder_proposes", "orchestrator_applies", "independent_read_only_verification"],
    }
    if args.mutation == "TOOTH_PARAMETER_REGISTER": proposals["register_modified_by_builder"] = True
    require(not proposals["register_modified_by_builder"], "ASSERT_PARAMETER_REGISTER_UNCHANGED", "builder modified frozen parameter register")
    components: dict[str, Any] = {
        "frozen_data_pin_table": pin_table,
        "candidate_inventory": {"schema_version": "U2_CANDIDATE_INVENTORY_V1", **semantic["candidate_inventory"]},
        "obligation_censuses": {"schema_version": "U2_OBLIGATION_CENSUSES_V1", "censuses": semantic["obligation_censuses"]},
        "dependency_grid_inventory": {
            "schema_version": "U2_DEPENDENCY_GRID_INVENTORY_V1",
            **semantic["open_dependency_relation"], **semantic["grid_inventory"],
        },
        "vocabulary_freeze": {
            "schema_version": "U2_VOCABULARY_FREEZE_V1", **semantic["vocabulary_freeze"],
            "executable_decision_and_record_guard_fixtures": {
                key: value for key, value in semantic["guard_fixtures"].items()
                if key not in {"closure", "template"}
            },
        },
        "evidence_taxonomy": {"schema_version": "U2_EVIDENCE_TAXONOMY_V1", **semantic["evidence_taxonomy"]},
        "availability_slots": {"schema_version": "U2_AVAILABILITY_SLOTS_V1", "summary": semantic["availability_summary"], "slots": slots},
        "route_fixture_inventory": {"schema_version": "U2_ROUTE_FIXTURE_INVENTORY_V1", **semantic["route_fixture_inventory"]},
        "closure_template_contracts": {
            "schema_version": "U2_CLOSURE_TEMPLATE_CONTRACTS_V1", "closure": semantic["closure_contract"],
            "template": semantic["template_contract"], "return_closure_ownership": semantic["return_closure_ownership"],
            "closure_guard_fixture": semantic["guard_fixtures"]["closure"],
            "template_guard_fixture": semantic["guard_fixtures"]["template"],
            "physics_record_invariance_contract": semantic["physics_record_invariance_contract"],
        },
        "environment_identity": environment,
        "standard_bindings": {"schema_version": "U2_STANDARD_BINDINGS_V1", "bindings": semantic["standard_bindings"]},
        "producer_map": producer,
        "evaluated_code_closure_policy": closure_policy(repo, environment),
        "parameter_register_proposals": proposals,
    }
    if args.mutation == "TOOTH_OBLIGATION_MANIFEST":
        # Mutate a real emitted typed sink; the independent floor remains intact.
        components["candidate_inventory"]["candidate_axis"].pop()
    manifest = obligation_manifest(components, agreement)
    require(manifest["exact_set_equal"], "ASSERT_OBLIGATION_MANIFEST", "expected/reachable obligation sets differ")
    components["obligation_manifest"] = manifest
    for name, filename in COMPONENT_FILES.items():
        dump_yaml(out / filename, components[name])
    component_digests = {name: digest(value) for name, value in components.items()}
    bundle_digest = digest(components)
    bundle = {
        "schema_version": "U2_STAGE0_BUNDLE_V1",
        "canonicalization": "UTF-8 sorted-key compact JSON over parsed YAML values",
        "startup_contract_commit": args.startup_contract_commit,
        "directive_provenance_commit": DIRECTIVE_PROVENANCE_COMMIT,
        "component_count": len(components),
        "components": [{
            "name": name, "path": COMPONENT_FILES[name], "semantic_sha256": component_digests[name],
            "file_sha256": sha256_path(out / COMPONENT_FILES[name]),
        } for name in sorted(components)],
        "stage0_contract_digest": bundle_digest,
    }
    dump_yaml(out / "stage0_bundle.yaml", bundle)
    integrity = semantic["availability_summary"]["integrity_failures"]
    require(integrity == 0, "ASSERT_CONTRACT_ZERO_INTEGRITY", "stage0 bundle contains integrity records")
    contract = {
        "schema_version": "U2_STAGE0_CONTRACT_V1", "phase": "U2_STAGE0",
        "status": "AWAITING_ORCHESTRATOR_APPROVAL", "required_exit_code": 42,
        "startup_contract_commit": args.startup_contract_commit,
        "directive_provenance_commit": DIRECTIVE_PROVENANCE_COMMIT,
        "stage0_contract_digest": bundle_digest, "bundle_index_sha256": sha256_path(out / "stage0_bundle.yaml"),
        "engine_agreement_sha256": sha256_path(Path(args.agreement)), "component_semantic_sha256": component_digests,
        "candidate_universe_digest": semantic["candidate_inventory"]["candidate_universe_digest"],
        "grid_summary": {
            "candidates": semantic["grid_inventory"]["candidate_count"],
            "ambient_branches": semantic["grid_inventory"]["ambient_count"],
            "active_strata": len(semantic["open_dependency_relation"]["active_strata"]),
            "raw_ragged_cardinality": semantic["grid_inventory"]["raw_ragged_cardinality"],
            "collapsed_cardinality": semantic["grid_inventory"]["collapsed_cardinality"],
            "promotion_contexts": semantic["grid_inventory"]["promotion_context_count"],
            "template_branch_equivalence_proofs": len(semantic["template_contract"]["branch_equivalence_proofs"]),
            "expected_conditional_templates_on_v11_landing": 16,
        },
        "availability_summary": semantic["availability_summary"],
        "vocabulary_summary": {
            "candidate_dispositions": 3, "promotion_physics_outcomes": 3,
            "ensemble_level_1": 4, "ensemble_level_2": 4, "topology_gate": 3,
            "integrity_states": 3, "failure_reasons": len(semantic["vocabulary_freeze"]["typed_failure_reasons"]),
        },
        "mutation_catalog": contract_mutation_catalog(),
        "production_resume_inputs": {
            "STARTUP_CONTRACT_COMMIT": "orchestrator commit containing reviewed U2 scripts",
            "STAGE0_CONTRACT_DIGEST": bundle_digest,
        },
        "approval_requests": [
            "ratify the U2 frozen-data pin table and complete stage-0 bundle",
            "ratify the v11 992-record non-template projection digest and 16 branch identical-BVP proofs",
            "adjudicate the generated three-family candidate universe and unresolved residual OTHER branch",
            "commit reviewed U2 scripts and supply the production STARTUP_CONTRACT_COMMIT",
            "supply STAGE0_CONTRACT_DIGEST unchanged at production resume",
        ],
    }
    if args.mutation == "TOOTH_APPROVAL_EXIT": contract["required_exit_code"] = 0
    require(contract["required_exit_code"] == 42, "ASSERT_APPROVAL_EXIT_42", "stage0 success-stop is not exit 42")
    dump_yaml(out / "stage0_contract.yaml", contract)
    print(f"U2_STAGE0_BUNDLE_READY digest={bundle_digest} components={len(components)} grid={contract['grid_summary']['raw_ragged_cardinality']}")
    return 0


def verify(args: argparse.Namespace) -> int:
    repo = Path(args.repo).resolve(); out = Path(args.output_dir).resolve()
    validate_anchor(repo, args.startup_contract_commit, args.mutation)
    require(bool(args.expected_digest), "ASSERT_EXPECTED_DIGEST_SUPPLIED", "expected stage0 digest missing")
    components = {name: load_yaml(out / filename) for name, filename in COMPONENT_FILES.items()}
    if args.mutation and args.mutation.startswith("BUNDLE_COMPONENT:"):
        name = args.mutation.split(":", 1)[1]
        require(name in components, "MUTATION_NOOP", f"unknown bundle component {name}")
        components[name]["_mutated_semantic_field"] = True
    actual = digest(components)
    expected = "0" * 64 if args.mutation == "TOOTH_STAGE0_CONTRACT_DIGEST" else args.expected_digest
    require(actual == expected, "ASSERT_STAGE0_CONTRACT_DIGEST", f"expected {expected}, recomputed {actual}")
    bundle_index = load_yaml(out / "stage0_bundle.yaml")
    contract = load_yaml(out / "stage0_contract.yaml")
    recorded_bundle_anchor = bundle_index.get("startup_contract_commit")
    recorded_contract_anchor = contract.get("startup_contract_commit")
    if args.mutation == "TOOTH_RECORDED_STARTUP_ANCHOR":
        recorded_contract_anchor = "0" * 40
    require(
        recorded_bundle_anchor == recorded_contract_anchor
        and bool(re.fullmatch(r"[0-9a-f]{40}", recorded_bundle_anchor or ""))
        and recorded_bundle_anchor != "HEAD",
        "ASSERT_RECORDED_STARTUP_ANCHOR", "stage0 bundle/contract do not record one full concrete startup anchor",
    )
    recorded_resolved = git(repo, "rev-parse", f"{recorded_bundle_anchor}^{{commit}}").stdout.strip()
    require(
        recorded_resolved == recorded_bundle_anchor,
        "ASSERT_RECORDED_STARTUP_ANCHOR", "recorded stage0 startup anchor did not resolve identically",
    )
    require(
        bundle_index.get("stage0_contract_digest") == contract.get("stage0_contract_digest") == actual
        and contract.get("bundle_index_sha256") == sha256_path(out / "stage0_bundle.yaml"),
        "ASSERT_STAGE0_CONTRACT_DIGEST", "stage0 index/contract binding differs from recomputed components",
    )
    live_pins = build_pin_table(repo, args.mutation)
    if args.mutation == "TOOTH_VERIFY_PIN_TABLE": live_pins["records"][0]["sha256"] = "0" * 64
    require(live_pins == components["frozen_data_pin_table"], "ASSERT_VERIFY_PIN_TABLE", "live pin table differs from bundled pins")
    if args.recompute_environment:
        live_environment = environment_identity(args.mutation)
        if args.mutation == "TOOTH_ENVIRONMENT_IDENTITY": live_environment["python"]["version"] = "mutated"
        require(live_environment == components["environment_identity"], "ASSERT_ENVIRONMENT_IDENTITY", "environment/toolchain identity changed")
    required_exit = 0 if args.mutation == "TOOTH_APPROVAL_EXIT" else contract["required_exit_code"]
    require(required_exit == 42, "ASSERT_APPROVAL_EXIT_42", "contract approval exit changed")
    integrity = 1 if args.mutation == "TOOTH_CONTRACT_ZERO_INTEGRITY" else components["availability_slots"]["summary"]["integrity_failures"]
    require(integrity == 0, "ASSERT_CONTRACT_ZERO_INTEGRITY", "bundle contains integrity failures")
    print(f"U2_STAGE0_CONTRACT_VERIFIED digest={actual}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True); parser.add_argument("--output-dir", required=True)
    parser.add_argument("--startup-contract-commit", required=True)
    parser.add_argument("--sympy"); parser.add_argument("--agreement")
    parser.add_argument("--verify", action="store_true"); parser.add_argument("--expected-digest")
    parser.add_argument("--recompute-environment", action="store_true"); parser.add_argument("--mutation")
    args = parser.parse_args()
    try:
        if args.verify:
            return verify(args)
        require(bool(args.sympy) and bool(args.agreement), "ASSERT_CONTRACT_ARGUMENTS", "assemble requires engine and agreement")
        return assemble(args)
    except ContractFailure as failure:
        print(f"ASSERTION_FAILED {failure.assert_id}: {failure.detail}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
