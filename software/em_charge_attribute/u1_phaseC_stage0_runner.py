#!/usr/bin/env python3
"""Staged, resumable U1 Phase-C stage-0 runner; exit 42 is success-stop."""

from __future__ import annotations

import argparse
import ast
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml


FULL_STAGE0_ANCHOR = "377eab17a4babc12847450956dc55fe3e16d33da"
SUCCESS_HALT = 42
ARTIFACT_REL = Path(
    "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/"
    "stage_c_0_tilt_coupling_contract"
)
SCRATCH_REL = Path("software/em_charge_attribute/_scratch/u1_phaseC_stage0")
STAGE_SEQUENCE = (
    "00_runner_shell_probe",
    "01_sympy_engine",
    "02_wolfram_engine",
    "03_engine_comparator",
    "04_bundle_contract",
    "05_preliminary_closure",
    "06_containment",
    "07_mutation_campaign",
    "07b_mutation_closure_probe",
    "08_final_closure",
    "09_summary_refresh",
    "10_contract_resume_preflight",
)


class StageFailure(RuntimeError):
    pass


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=120),
        encoding="utf-8",
    )


class Runner:
    def __init__(self, repo: Path, anchor: str):
        self.repo = repo
        self.anchor = anchor
        self.source = repo / "software/em_charge_attribute"
        self.artifact = repo / ARTIFACT_REL
        self.scratch = repo / SCRATCH_REL
        self.runtime_tmp = self.scratch / "runtime_tmp"
        self.state_path = self.scratch / "runner_state.yaml"
        self.records_path = self.scratch / "stage_run_records.yaml"
        self.artifact.mkdir(parents=True, exist_ok=True)
        self.scratch.mkdir(parents=True, exist_ok=True)
        self.runtime_tmp.mkdir(parents=True, exist_ok=True)
        (self.runtime_tmp / "wolfram_userbase/Licensing").mkdir(parents=True, exist_ok=True)
        (self.runtime_tmp / "wolfram_base").mkdir(parents=True, exist_ok=True)
        self.nsswitch = self.runtime_tmp / "phaseC_nsswitch.conf"
        self.nsswitch.write_text("hosts: files\n", encoding="utf-8")
        if self.state_path.exists():
            self.state = load_yaml(self.state_path)
            if self.state.get("startup_contract_commit") != anchor:
                raise StageFailure("resume state belongs to a different startup anchor")
            self.state["schema_version"] = "U1_PHASE_C_STAGE0_RUNNER_STATE_V2"
        else:
            self.state = {
                "schema_version": "U1_PHASE_C_STAGE0_RUNNER_STATE_V2",
                "startup_contract_commit": anchor,
                "terminal": "IN_PROGRESS",
                "completed_stages": {},
            }
            self.save_state()
        self.run_records = load_yaml(self.records_path) if self.records_path.exists() else {
            "schema_version": "U1_PHASE_C_STAGE_RUN_RECORDS_V2",
            "records": [],
        }
        self.run_records["schema_version"] = "U1_PHASE_C_STAGE_RUN_RECORDS_V2"

    def save_state(self) -> None:
        dump_yaml(self.state_path, self.state)

    def save_records(self) -> None:
        dump_yaml(self.records_path, self.run_records)

    def trace_path(self, stage: str) -> Path:
        return self.scratch / "traces" / f"{stage}.strace"

    def producer_hashes(self, producers: list[Path]) -> dict[str, str]:
        hashes: dict[str, str] = {}
        for path in producers:
            resolved = path.resolve()
            if not resolved.is_file():
                raise StageFailure(f"stage producer is missing: {resolved}")
            try:
                name = str(resolved.relative_to(self.repo))
            except ValueError:
                name = str(resolved)
            hashes[name] = sha256_path(resolved)
        return dict(sorted(hashes.items()))

    def stage_producers(self, *names: str) -> list[Path]:
        return [self.source / "u1_phaseC_stage0_runner.py", *[self.source / name for name in names]]

    def invalidate_from(self, stage: str) -> None:
        first = STAGE_SEQUENCE.index(stage)
        invalidated = [
            name
            for name in STAGE_SEQUENCE[first:]
            if name in self.state["completed_stages"]
        ]
        if not invalidated:
            return
        for name in invalidated:
            del self.state["completed_stages"][name]
        self.state["terminal"] = "IN_PROGRESS"
        self.state.pop("exit_code", None)
        self.state.pop("stage0_contract_digest", None)
        self.save_state()
        print(
            f"PHASEC_STAGE_INVALIDATE from={stage} records={','.join(invalidated)}",
            flush=True,
        )

    def completed(
        self, stage: str, required: list[Path], producer_sha256: dict[str, str]
    ) -> bool:
        record = self.state["completed_stages"].get(stage)
        if not record or record.get("status") != "PASS":
            return False
        if record.get("producer_script_sha256") != producer_sha256:
            return False
        expected = record.get("output_sha256", {})
        return all(
            path.is_file() and expected.get(str(path)) == sha256_path(path)
            for path in required
        )

    def sandbox_command(self, command: list[str]) -> list[str]:
        return [
            "/usr/bin/bwrap",
            "--ro-bind",
            "/",
            "/",
            "--dev-bind",
            "/dev",
            "/dev",
            "--proc",
            "/proc",
            "--unshare-net",
            "--bind",
            str(self.artifact),
            str(self.artifact),
            "--bind",
            str(self.scratch),
            str(self.scratch),
            "--bind",
            str(self.runtime_tmp),
            "/tmp",
            "--ro-bind",
            "/home/trevnorris/.Wolfram/Licensing",
            "/tmp/wolfram_userbase/Licensing",
            "--ro-bind",
            str(self.nsswitch),
            "/etc/nsswitch.conf",
            "--chdir",
            str(self.repo),
            "--setenv",
            "PYTHONNOUSERSITE",
            "1",
            "--setenv",
            "PYTHONDONTWRITEBYTECODE",
            "1",
            "--setenv",
            "TMPDIR",
            "/tmp",
            "--setenv",
            "WOLFRAM_USERBASE",
            "/tmp/wolfram_userbase",
            "--setenv",
            "WOLFRAM_BASE",
            "/tmp/wolfram_base",
            "--",
            *command,
        ]

    def run_stage(
        self,
        stage: str,
        command: list[str],
        required: list[Path],
        producers: list[Path],
        timeout: int = 600,
        outer_trace: bool = True,
    ) -> None:
        trace = self.trace_path(stage)
        all_required = [*required, *([trace] if outer_trace else [])]
        producer_sha256 = self.producer_hashes(producers)
        if self.completed(stage, all_required, producer_sha256):
            print(f"PHASEC_STAGE_RESUME_SKIP stage={stage}", flush=True)
            return
        self.invalidate_from(stage)
        trace.parent.mkdir(parents=True, exist_ok=True)
        launched = self.sandbox_command(command)
        if outer_trace:
            launched = [
                "/usr/bin/strace",
                "-f",
                "-qq",
                "-s",
                "4096",
                "-e",
                "trace=%file,%network",
                "-o",
                str(trace),
                *launched,
            ]
        print(f"PHASEC_STAGE_START stage={stage}", flush=True)
        started = time.monotonic()
        try:
            result = subprocess.run(
                launched,
                cwd=self.repo,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise StageFailure(f"stage {stage} exceeded {timeout}s; reformulation required") from exc
        elapsed = time.monotonic() - started
        self.run_records["records"].append(
            {
                "stage": stage,
                "return_code": result.returncode,
                "elapsed_seconds": round(elapsed, 6),
                "stdout": result.stdout,
                "stderr": result.stderr,
                "trace": str(trace) if outer_trace else "nested_campaign_process_evidence",
                "trace_sha256": sha256_path(trace) if outer_trace and trace.is_file() else None,
            }
        )
        self.save_records()
        if result.returncode != 0:
            raise StageFailure(
                f"stage {stage} failed rc={result.returncode}: {(result.stdout + result.stderr)[-1600:]}"
            )
        missing = [str(path) for path in all_required if not path.is_file()]
        if missing:
            raise StageFailure(f"stage {stage} omitted required outputs: {missing}")
        self.state["completed_stages"][stage] = {
            "status": "PASS",
            "elapsed_seconds": round(elapsed, 6),
            "producer_script_sha256": producer_sha256,
            "output_sha256": {str(path): sha256_path(path) for path in all_required},
        }
        self.save_state()
        print(f"PHASEC_STAGE_PASS stage={stage} elapsed={elapsed:.3f}s", flush=True)

    def python_command(self, script: str, *args: str) -> list[str]:
        return ["/usr/bin/python3", "-I", str(self.source / script), *args]

    def run(self) -> int:
        sympy = self.artifact / "sympy_stage0.yaml"
        wolfram = self.artifact / "wolfram_stage0.yaml"
        agreement = self.artifact / "stage0_engine_agreement.yaml"
        preliminary_closure = self.scratch / "preliminary_code_closure.yaml"
        closure = self.artifact / "evaluated_code_closure_evidence.yaml"
        containment = self.artifact / "containment_evidence.yaml"
        mutations = self.artifact / "mutation_campaign.yaml"
        summary = self.artifact / "stage0_summary.md"

        self.run_stage(
            "00_runner_shell_probe",
            ["/bin/bash", str(self.source / "run_u1_body_phaseC.sh"), "--shell-probe"],
            [self.scratch / "runner_shell_probe.yaml"],
            self.stage_producers("run_u1_body_phaseC.sh"),
        )
        self.run_stage(
            "01_sympy_engine",
            self.python_command(
                "u1_phaseC_stage0_sympy.py", "--repo", str(self.repo), "--output", str(sympy)
            ),
            [sympy],
            self.stage_producers("u1_phaseC_stage0_sympy.py"),
        )
        self.run_stage(
            "02_wolfram_engine",
            [
                "/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel",
                "-noinit",
                "-noprompt",
                "-script",
                str(self.source / "u1_phaseC_stage0_dual.wl"),
                "--repo",
                str(self.repo),
                "--output",
                str(wolfram),
            ],
            [wolfram],
            self.stage_producers("u1_phaseC_stage0_dual.wl"),
        )
        self.run_stage(
            "03_engine_comparator",
            self.python_command(
                "u1_phaseC_stage0_compare.py",
                "--repo",
                str(self.repo),
                "--sympy",
                str(sympy),
                "--wolfram",
                str(wolfram),
                "--output",
                str(agreement),
            ),
            [agreement],
            self.stage_producers("u1_phaseC_stage0_compare.py"),
        )
        self.run_stage(
            "04_bundle_contract",
            self.python_command(
                "u1_phaseC_stage0_contract.py",
                "--repo",
                str(self.repo),
                "--output-dir",
                str(self.artifact),
                "--startup-contract-commit",
                self.anchor,
                "--sympy",
                str(sympy),
                "--wolfram",
                str(wolfram),
                "--agreement",
                str(agreement),
            ),
            [
                self.artifact / "stage0_contract.yaml",
                self.artifact / "stage0_bundle.yaml",
                *[self.artifact / name for name in (
                    "frozen_data_pin_table.yaml",
                    "availability_slots.yaml",
                    "reconciliation_inventory.yaml",
                    "force_term_census.yaml",
                    "coupling_source_census.yaml",
                    "g8_ablation_inventory.yaml",
                    "projection_freeze.yaml",
                    "environment_identity.yaml",
                    "producer_map.yaml",
                    "obligation_manifest.yaml",
                    "parameter_register_proposals.yaml",
                    "evaluated_code_closure_policy.yaml",
                )],
            ],
            self.stage_producers("u1_phaseC_stage0_contract.py"),
        )
        first_traces = [
            self.trace_path(name)
            for name in (
                "00_runner_shell_probe",
                "01_sympy_engine",
                "02_wolfram_engine",
                "03_engine_comparator",
                "04_bundle_contract",
            )
        ]
        closure_args: list[str] = [
            "--repo",
            str(self.repo),
            "--anchor",
            self.anchor,
            "--environment",
            str(self.artifact / "environment_identity.yaml"),
        ]
        for trace in first_traces:
            closure_args.extend(["--trace", str(trace)])
        self.run_stage(
            "05_preliminary_closure",
            self.python_command(
                "u1_phaseC_code_closure_guard.py",
                *closure_args,
                "--stage0-precommit",
                "--output",
                str(preliminary_closure),
            ),
            [preliminary_closure],
            self.stage_producers("u1_phaseC_code_closure_guard.py"),
        )
        normal_traces = [*first_traces, self.trace_path("05_preliminary_closure")]
        containment_args: list[str] = [
            "--cwd",
            str(self.repo),
            "--allow-write-root",
            str(self.artifact),
            "--allow-write-root",
            str(self.scratch),
            "--allow-write-root",
            "/__phaseC_namespace_control__",
            "--allow-write-root",
            "/tmp/wolfram_userbase",
            "--allow-write-root",
            "/tmp/wolfram_base",
            "--mapped-write-root",
            f"/tmp={self.runtime_tmp}",
        ]
        for trace in normal_traces:
            containment_args.extend(["--trace", str(trace)])
        self.run_stage(
            "06_containment",
            self.python_command(
                "u1_phaseC_containment.py",
                *containment_args,
                "--output",
                str(containment),
            ),
            [containment],
            self.stage_producers("u1_phaseC_containment.py"),
        )
        mutation_scratch = self.scratch / "mutations"
        self.run_stage(
            "07_mutation_campaign",
            self.python_command(
                "u1_phaseC_stage0_mutations.py",
                "--repo",
                str(self.repo),
                "--scratch",
                str(mutation_scratch),
                "--sympy",
                str(sympy),
                "--wolfram",
                str(wolfram),
                "--agreement",
                str(agreement),
                "--bundle-dir",
                str(self.artifact),
                "--closure",
                str(preliminary_closure),
                "--containment",
                str(containment),
                "--output",
                str(mutations),
                "--summary-output",
                str(summary),
            ),
            [mutations],
            self.stage_producers(
                "u1_phaseC_stage0_mutations.py",
                "u1_phaseC_stage0_compare.py",
                "u1_phaseC_stage0_contract.py",
                "u1_phaseC_code_closure_guard.py",
                "u1_phaseC_containment.py",
            ),
            outer_trace=False,
        )
        self.run_stage(
            "07b_mutation_closure_probe",
            self.python_command(
                "u1_phaseC_stage0_mutations.py",
                "--summary-probe",
                str(mutation_scratch / "summary_probe.yaml"),
            ),
            [],
            self.stage_producers("u1_phaseC_stage0_mutations.py"),
        )
        final_closure_args = list(closure_args)
        for stage in ("05_preliminary_closure", "06_containment", "07b_mutation_closure_probe"):
            final_closure_args.extend(["--trace", str(self.trace_path(stage))])
        self.run_stage(
            "08_final_closure",
            self.python_command(
                "u1_phaseC_code_closure_guard.py",
                *final_closure_args,
                "--stage0-precommit",
                "--output",
                str(closure),
            ),
            [closure],
            self.stage_producers("u1_phaseC_code_closure_guard.py"),
        )
        self.run_stage(
            "09_summary_refresh",
            self.python_command(
                "u1_phaseC_stage0_mutations.py",
                "--summary-only",
                "--bundle-dir",
                str(self.artifact),
                "--closure",
                str(closure),
                "--containment",
                str(containment),
                "--output",
                str(mutations),
                "--summary-output",
                str(summary),
            ),
            [summary],
            self.stage_producers("u1_phaseC_stage0_mutations.py"),
        )
        contract = load_yaml(self.artifact / "stage0_contract.yaml")
        contract_digest = contract["stage0_contract_digest"]
        self.run_stage(
            "10_contract_resume_preflight",
            self.python_command(
                "u1_phaseC_stage0_contract.py",
                "--repo",
                str(self.repo),
                "--output-dir",
                str(self.artifact),
                "--verify",
                "--expected-digest",
                contract_digest,
                "--startup-contract-commit",
                self.anchor,
                "--recompute-environment",
            ),
            [self.artifact / "stage0_bundle.yaml"],
            self.stage_producers("u1_phaseC_stage0_contract.py"),
        )
        mutation_record = load_yaml(mutations)
        if mutation_record["status"] != "PASS" or mutation_record["vacuous_case_count"] != 0:
            raise StageFailure("mutation campaign did not close with zero vacuous cases")
        if contract["required_exit_code"] != SUCCESS_HALT:
            raise StageFailure("stage0 contract does not require success-stop exit 42")
        halt = {
            "schema_version": "U1_PHASE_C_STAGE0_HALT_V1",
            "status": "SUCCESS_STOP_AWAITING_ORCHESTRATOR_APPROVAL",
            "exit_code": SUCCESS_HALT,
            "startup_contract_commit": self.anchor,
            "stage0_contract_digest": contract_digest,
            "availability_summary": contract["availability_summary"],
            "reconciliation_summary": contract["reconciliation_summary"],
            "artifact_directory": str(self.artifact),
            "orchestrator_adjudication": contract["approval_requests"],
        }
        dump_yaml(self.artifact / "orchestrator_approval_halt.yaml", halt)
        self.state["terminal"] = "SUCCESS_STOP_AWAITING_ORCHESTRATOR_APPROVAL"
        self.state["exit_code"] = SUCCESS_HALT
        self.state["stage0_contract_digest"] = contract_digest
        self.save_state()
        print(
            f"PHASEC_STAGE0_SUCCESS_STOP exit=42 digest={contract_digest} "
            f"derived={contract['availability_summary']['DERIVED']} "
            f"unresolved={contract['availability_summary']['UNRESOLVED']}",
            flush=True,
        )
        return SUCCESS_HALT


def validate_sources(repo: Path) -> int:
    source = repo / "software/em_charge_attribute"
    python_files = sorted(source.glob("u1_phaseC_*.py"))
    for path in python_files:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    required = [
        source / "u1_phaseC_stage0_dual.wl",
        source / "run_u1_body_phaseC.sh",
        Path("/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel"),
        Path("/usr/bin/bwrap"),
        Path("/usr/bin/strace"),
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise StageFailure(f"runner prerequisites missing: {missing}")
    print(f"PHASEC_RUNNER_VALID sources={len(python_files) + 2} stages=12")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument("--startup-contract-commit")
    parser.add_argument("--stage0", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--self-probe", action="store_true")
    args = parser.parse_args()
    repo = Path(args.repo).resolve()
    if args.self_probe:
        output = repo / SCRATCH_REL / "runner_shell_probe.yaml"
        dump_yaml(
            output,
            {
                "schema_version": "U1_PHASE_C_RUNNER_SHELL_PROBE_V1",
                "python_isolated": sys.flags.isolated == 1,
                "python_no_user_site": sys.flags.no_user_site == 1,
                "runner_path": str(Path(__file__).resolve()),
            },
        )
        print("PHASEC_RUNNER_SHELL_PROBE_PASS")
        return 0
    if args.validate_only:
        return validate_sources(repo)
    if not args.stage0:
        parser.error("this build runner requires --stage0")
    if not args.startup_contract_commit:
        parser.error("--startup-contract-commit is orchestrator-supplied and required")
    resolved = subprocess.run(
        ["git", "rev-parse", f"{args.startup_contract_commit}^{{commit}}"],
        cwd=repo,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip()
    if resolved != FULL_STAGE0_ANCHOR:
        raise StageFailure(
            f"stage0 anchor mismatch: supplied {resolved}, required {FULL_STAGE0_ANCHOR}"
        )
    return Runner(repo, resolved).run()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except StageFailure as failure:
        print(f"PHASEC_STAGE0_BLOCKED {failure}", file=sys.stderr)
        raise SystemExit(1)
