#!/usr/bin/env python3
"""Staged, resumable U2 stage-0 runner; exit 42 is the success-stop."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml


SUCCESS_HALT = 42
ARTIFACT_REL = Path("software/em_charge_attribute/reports/u2_boundary_adjudication_artifacts/stage_0_contract")
SCRATCH_REL = Path("software/em_charge_attribute/_scratch/u2_stage0")
STAGE_SEQUENCE = (
    "00_runner_shell_probe", "01_sympy_engine", "02_wolfram_engine", "03_engine_comparator",
    "04_bundle_contract", "05_preliminary_closure", "06_process_tree_containment",
    "07_mutation_campaign", "08_mutation_closure_probe", "09_final_closure",
    "10_summary_refresh", "11_contract_resume_preflight",
)
COMPONENT_FILES = (
    "frozen_data_pin_table.yaml", "candidate_inventory.yaml", "obligation_censuses.yaml",
    "dependency_grid_inventory.yaml", "vocabulary_freeze.yaml", "evidence_taxonomy.yaml",
    "availability_slots.yaml", "route_fixture_inventory.yaml", "closure_template_contracts.yaml",
    "environment_identity.yaml", "standard_bindings.yaml", "producer_map.yaml",
    "evaluated_code_closure_policy.yaml", "parameter_register_proposals.yaml", "obligation_manifest.yaml",
)


class StageFailure(RuntimeError):
    pass


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=140), encoding="utf-8")


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class Runner:
    def __init__(self, repo: Path, anchor: str):
        self.repo = repo; self.anchor = anchor
        self.source = repo / "software/em_charge_attribute"
        self.artifact = repo / ARTIFACT_REL; self.scratch = repo / SCRATCH_REL
        self.runtime_tmp = self.scratch / "runtime_tmp"
        self.state_path = self.scratch / "runner_state.yaml"
        self.records_path = self.scratch / "stage_run_records.yaml"
        for path in (self.artifact, self.scratch, self.runtime_tmp, self.scratch / "traces", self.scratch / "stage_logs"):
            path.mkdir(parents=True, exist_ok=True)
        (self.runtime_tmp / "wolfram_userbase/Licensing").mkdir(parents=True, exist_ok=True)
        (self.runtime_tmp / "wolfram_base").mkdir(parents=True, exist_ok=True)
        self.nsswitch = self.runtime_tmp / "u2_nsswitch.conf"
        self.nsswitch.write_text("hosts: files\n", encoding="utf-8")
        if self.state_path.exists():
            self.state = load_yaml(self.state_path)
            if self.state.get("startup_contract_commit") != anchor:
                raise StageFailure("resume state belongs to a different startup anchor")
        else:
            self.state = {
                "schema_version": "U2_STAGE0_RUNNER_STATE_V1", "startup_contract_commit": anchor,
                "terminal": "IN_PROGRESS", "completed_stages": {},
            }
            self.save_state()
        self.run_records = load_yaml(self.records_path) if self.records_path.exists() else {
            "schema_version": "U2_STAGE_RUN_RECORDS_V1", "startup_contract_commit": anchor, "records": [],
        }

    def save_state(self) -> None:
        dump_yaml(self.state_path, self.state)

    def save_records(self) -> None:
        dump_yaml(self.records_path, self.run_records)

    def trace_path(self, stage: str) -> Path:
        return self.scratch / "traces" / f"{stage}.strace"

    def normalized_path(self, path: Path) -> str:
        resolved = path.resolve()
        try: return resolved.relative_to(self.repo).as_posix()
        except ValueError: return str(resolved)

    def file_set_digest(self, paths: list[Path], require_exists: bool = True) -> str:
        rows = []
        for path in sorted({value.resolve() for value in paths}, key=str):
            if not path.is_file():
                if require_exists: raise StageFailure(f"stage file missing: {path}")
                rows.append({"path": self.normalized_path(path), "sha256": None}); continue
            rows.append({"path": self.normalized_path(path), "sha256": sha256_path(path)})
        return canonical_digest(rows)

    def producer_record(self, producers: list[Path]) -> dict[str, Any]:
        paths = sorted(self.normalized_path(path) for path in producers)
        return {
            "producer_commit": self.anchor, "producer_paths": paths,
            "aggregate_producer_sha256": self.file_set_digest(producers),
            "per_script_content_map_present": False,
        }

    def stage_producers(self, *names: str) -> list[Path]:
        return [self.source / "u2_stage0_runner.py", *[self.source / name for name in names]]

    def invalidate_from(self, stage: str) -> None:
        first = STAGE_SEQUENCE.index(stage)
        invalidated = [value for value in STAGE_SEQUENCE[first:] if value in self.state["completed_stages"]]
        if not invalidated: return
        for value in invalidated: del self.state["completed_stages"][value]
        self.state["terminal"] = "IN_PROGRESS"
        for key in ("exit_code", "stage0_contract_digest"): self.state.pop(key, None)
        self.save_state(); print(f"U2_STAGE_INVALIDATE from={stage} stages={','.join(invalidated)}", flush=True)

    def completed(
        self, stage: str, outputs: list[Path], producer: dict[str, Any], input_digest: str,
    ) -> bool:
        row = self.state["completed_stages"].get(stage)
        if not row or row.get("status") != "PASS": return False
        if row.get("producer") != producer or row.get("input_digest") != input_digest: return False
        return row.get("output_digest") == self.file_set_digest(outputs, require_exists=False) and all(path.is_file() for path in outputs)

    def sandbox_command(self, command: list[str]) -> list[str]:
        return [
            "/usr/bin/bwrap", "--ro-bind", "/", "/", "--dev-bind", "/dev", "/dev", "--proc", "/proc",
            "--unshare-net", "--bind", str(self.artifact), str(self.artifact),
            "--bind", str(self.scratch), str(self.scratch), "--bind", str(self.runtime_tmp), "/tmp",
            "--ro-bind", "/home/trevnorris/.Wolfram/Licensing", "/tmp/wolfram_userbase/Licensing",
            "--ro-bind", str(self.nsswitch), "/etc/nsswitch.conf", "--chdir", str(self.repo),
            "--setenv", "PYTHONNOUSERSITE", "1", "--setenv", "PYTHONDONTWRITEBYTECODE", "1",
            "--setenv", "TMPDIR", "/tmp", "--setenv", "WOLFRAM_USERBASE", "/tmp/wolfram_userbase",
            "--setenv", "WOLFRAM_BASE", "/tmp/wolfram_base", "--", *command,
        ]

    def execute_with_progress(self, stage: str, command: list[str], timeout: int) -> tuple[int, str, str, float]:
        stdout_path = self.scratch / "stage_logs" / f"{stage}.stdout"
        stderr_path = self.scratch / "stage_logs" / f"{stage}.stderr"
        started = time.monotonic()
        with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
            process = subprocess.Popen(command, cwd=self.repo, stdout=stdout_handle, stderr=stderr_handle)
            next_progress = 30.0
            while process.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed >= timeout:
                    process.kill(); process.wait()
                    raise StageFailure(f"stage {stage} exceeded {timeout}s; reformulation required")
                if elapsed >= next_progress:
                    print(f"U2_STAGE_PROGRESS stage={stage} elapsed={elapsed:.1f}s", flush=True); next_progress += 30.0
                time.sleep(0.25)
        elapsed = time.monotonic() - started
        return process.returncode, stdout_path.read_text(errors="replace"), stderr_path.read_text(errors="replace"), elapsed

    def run_stage(
        self, stage: str, command: list[str], outputs: list[Path], producers: list[Path],
        inputs: list[Path] | None = None, timeout: int = 600, outer_trace: bool = True,
    ) -> None:
        trace = self.trace_path(stage); required = [*outputs, *([trace] if outer_trace else [])]
        producer = self.producer_record(producers)
        input_digest = canonical_digest({
            "startup_contract_commit": self.anchor,
            "file_set_sha256": self.file_set_digest(inputs or []),
            "producer_aggregate_sha256": producer["aggregate_producer_sha256"],
        })
        if self.completed(stage, required, producer, input_digest):
            print(f"U2_STAGE_RESUME_SKIP stage={stage}", flush=True); return
        self.invalidate_from(stage)
        launched = self.sandbox_command(command)
        if outer_trace:
            launched = [
                "/usr/bin/strace", "-f", "-qq", "-s", "4096", "-e", "trace=%file,%network",
                "-o", str(trace), *launched,
            ]
        print(f"U2_STAGE_START stage={stage}", flush=True)
        return_code, stdout, stderr, elapsed = self.execute_with_progress(stage, launched, timeout)
        record = {
            "stage": stage, "status": "PASS" if return_code == 0 else "FAIL", "return_code": return_code,
            "elapsed_seconds": round(elapsed, 6), "startup_contract_commit": self.anchor,
            "producer": producer, "input_digest": input_digest,
            "trace": str(trace) if outer_trace else "nested_campaign_evidence",
            "trace_sha256": sha256_path(trace) if outer_trace and trace.is_file() else None,
            "stdout_sha256": hashlib.sha256(stdout.encode()).hexdigest(),
            "stderr_sha256": hashlib.sha256(stderr.encode()).hexdigest(),
            "stdout_tail": stdout[-1200:], "stderr_tail": stderr[-1200:],
        }
        self.run_records["records"].append(record); self.save_records()
        if return_code != 0:
            raise StageFailure(f"stage {stage} failed rc={return_code}: {(stdout + stderr)[-2000:]}")
        missing = [str(path) for path in required if not path.is_file()]
        if missing: raise StageFailure(f"stage {stage} omitted required outputs: {missing}")
        output_digest = self.file_set_digest(required)
        self.state["completed_stages"][stage] = {
            "status": "PASS", "elapsed_seconds": round(elapsed, 6), "producer": producer,
            "input_digest": input_digest, "output_digest": output_digest,
        }
        self.save_state(); print(f"U2_STAGE_PASS stage={stage} elapsed={elapsed:.3f}s", flush=True)

    def python_command(self, script: str, *args: str) -> list[str]:
        return ["/usr/bin/python3", "-I", str(self.source / script), *args]

    def run(self) -> int:
        sympy = self.artifact / "sympy_stage0.yaml"; wolfram = self.artifact / "wolfram_stage0.yaml"
        agreement = self.artifact / "stage0_engine_agreement.yaml"
        preliminary_closure = self.scratch / "preliminary_code_closure.yaml"
        closure = self.artifact / "evaluated_code_closure_evidence.yaml"
        containment = self.artifact / "containment_evidence.yaml"
        mutations = self.artifact / "mutation_campaign.yaml"; summary = self.artifact / "stage0_summary.md"
        engine_inputs = [
            self.repo / "software/em_charge_attribute/u1_body_dynamics_inputs.yaml",
            self.repo / "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_0_tilt_coupling_contract/availability_slots.yaml",
            self.repo / "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_1_tilt_coupling_production/production_results.yaml",
        ]

        self.run_stage(
            "00_runner_shell_probe", ["/bin/bash", str(self.source / "run_u2_boundary_adjudication.sh"), "--shell-probe"],
            [self.scratch / "runner_shell_probe.yaml"], self.stage_producers("run_u2_boundary_adjudication.sh"),
        )
        self.run_stage(
            "01_sympy_engine", self.python_command("u2_stage0_sympy.py", "--repo", str(self.repo), "--output", str(sympy)),
            [sympy], self.stage_producers("u2_stage0_sympy.py"),
            inputs=engine_inputs,
        )
        self.run_stage(
            "02_wolfram_engine", [
                "/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel", "-noinit", "-noprompt", "-script",
                str(self.source / "u2_stage0_dual.wl"), "--repo", str(self.repo), "--output", str(wolfram),
            ], [wolfram], self.stage_producers("u2_stage0_dual.wl"),
            inputs=engine_inputs, timeout=1200,
        )
        self.run_stage(
            "03_engine_comparator", self.python_command(
                "u2_stage0_compare.py", "--sympy", str(sympy), "--wolfram", str(wolfram), "--output", str(agreement),
            ), [agreement], self.stage_producers("u2_stage0_compare.py"), inputs=[sympy, wolfram],
        )
        contract_outputs = [self.artifact / "stage0_contract.yaml", self.artifact / "stage0_bundle.yaml"] + [self.artifact / name for name in COMPONENT_FILES]
        self.run_stage(
            "04_bundle_contract", self.python_command(
                "u2_stage0_contract.py", "--repo", str(self.repo), "--output-dir", str(self.artifact),
                "--startup-contract-commit", self.anchor, "--sympy", str(sympy), "--agreement", str(agreement),
            ), contract_outputs, self.stage_producers("u2_stage0_contract.py"), inputs=[sympy, agreement],
        )
        first_stages = ("00_runner_shell_probe", "01_sympy_engine", "02_wolfram_engine", "03_engine_comparator", "04_bundle_contract")
        closure_base = [
            "--repo", str(self.repo), "--anchor", self.anchor,
            "--environment", str(self.artifact / "environment_identity.yaml"),
        ]
        for stage in first_stages: closure_base.extend(["--trace", str(self.trace_path(stage))])
        self.run_stage(
            "05_preliminary_closure", self.python_command(
                "u2_code_closure_guard.py", *closure_base, "--output", str(preliminary_closure),
            ), [preliminary_closure], self.stage_producers("u2_code_closure_guard.py"),
            inputs=[self.artifact / "environment_identity.yaml", *[self.trace_path(stage) for stage in first_stages]],
        )
        measured_stages = (*first_stages, "05_preliminary_closure")
        containment_args = [
            "--cwd", str(self.repo), "--allow-write-root", str(self.artifact), "--allow-write-root", str(self.scratch),
            "--allow-write-root", "/__u2_namespace_control__", "--allow-write-root", "/tmp/wolfram_userbase",
            "--allow-write-root", "/tmp/wolfram_base", "--mapped-write-root", f"/tmp={self.runtime_tmp}",
        ]
        for stage in measured_stages: containment_args.extend(["--trace", str(self.trace_path(stage))])
        mutation_scratch = self.scratch / "mutations"
        self.run_stage(
            "06_process_tree_containment", self.python_command(
                "u2_containment.py", *containment_args, "--output", str(containment),
            ), [containment], self.stage_producers("u2_containment.py"),
            inputs=[self.trace_path(stage) for stage in measured_stages],
        )
        self.run_stage(
            "07_mutation_campaign", self.python_command(
                "u2_stage0_mutations.py", "--repo", str(self.repo), "--scratch", str(mutation_scratch),
                "--startup-contract-commit", self.anchor,
                "--sympy", str(sympy), "--wolfram", str(wolfram), "--agreement", str(agreement),
                "--bundle-dir", str(self.artifact), "--closure", str(preliminary_closure),
                "--containment", str(containment),
                "--output", str(mutations), "--summary-output", str(summary),
            ), [mutations], self.stage_producers(
                "u2_stage0_mutations.py", "u2_stage0_compare.py", "u2_stage0_contract.py",
                "u2_code_closure_guard.py", "u2_containment.py",
            ), inputs=[sympy, wolfram, agreement, preliminary_closure, containment, *contract_outputs], outer_trace=False,
        )
        # Trace the campaign entrypoint itself without repeating the campaign; its summary probe is a real guarded route.
        self.run_stage(
            "08_mutation_closure_probe", self.python_command(
                "u2_stage0_mutations.py", "--summary-probe", str(mutation_scratch / "summary_probe.yaml"),
            ), [], self.stage_producers("u2_stage0_mutations.py"), inputs=[mutation_scratch / "summary_probe.yaml"],
        )
        final_closure_args = list(closure_base)
        for stage in ("05_preliminary_closure", "06_process_tree_containment", "08_mutation_closure_probe"):
            final_closure_args.extend(["--trace", str(self.trace_path(stage))])
        self.run_stage(
            "09_final_closure", self.python_command(
                "u2_code_closure_guard.py", *final_closure_args, "--output", str(closure),
            ), [closure], self.stage_producers("u2_code_closure_guard.py"),
            inputs=[self.artifact / "environment_identity.yaml", *[self.trace_path(stage) for stage in (*measured_stages, "06_process_tree_containment", "08_mutation_closure_probe")]],
        )
        self.run_stage(
            "10_summary_refresh", self.python_command(
                "u2_stage0_mutations.py", "--summary-only", "--bundle-dir", str(self.artifact),
                "--closure", str(closure), "--containment", str(containment),
                "--output", str(mutations), "--summary-output", str(summary),
            ), [summary], self.stage_producers("u2_stage0_mutations.py"),
            inputs=[mutations, closure, containment, *contract_outputs],
        )
        contract = load_yaml(self.artifact / "stage0_contract.yaml"); contract_digest = contract["stage0_contract_digest"]
        self.run_stage(
            "11_contract_resume_preflight", self.python_command(
                "u2_stage0_contract.py", "--repo", str(self.repo), "--output-dir", str(self.artifact),
                "--startup-contract-commit", self.anchor, "--verify", "--expected-digest", contract_digest,
                "--recompute-environment",
            ), [self.artifact / "stage0_bundle.yaml"], self.stage_producers("u2_stage0_contract.py"),
            inputs=contract_outputs,
        )
        campaign = load_yaml(mutations)
        if campaign["status"] != "PASS" or campaign["vacuous_case_count"] or campaign["control_failure_count"]:
            raise StageFailure("mutation campaign failed its nonvacuity/control gate")
        if contract["required_exit_code"] != SUCCESS_HALT:
            raise StageFailure("stage-0 contract does not require exit 42")
        halt = {
            "schema_version": "U2_STAGE0_APPROVAL_HALT_V1", "status": "SUCCESS_STOP_AWAITING_ORCHESTRATOR_APPROVAL",
            "exit_code": SUCCESS_HALT, "startup_contract_commit": self.anchor,
            "stage0_contract_digest": contract_digest, "grid_summary": contract["grid_summary"],
            "availability_summary": contract["availability_summary"],
            "vocabulary_summary": contract["vocabulary_summary"],
            "mutation_summary": {
                "teeth": campaign["tooth_count"], "controls": campaign["defect_absent_control_count"],
                "vacuous": campaign["vacuous_case_count"],
            },
            "artifact_directory": str(self.artifact), "orchestrator_adjudication": contract["approval_requests"],
        }
        dump_yaml(self.artifact / "orchestrator_approval_halt.yaml", halt)
        self.state.update({
            "terminal": "SUCCESS_STOP_AWAITING_ORCHESTRATOR_APPROVAL", "exit_code": SUCCESS_HALT,
            "stage0_contract_digest": contract_digest,
        }); self.save_state()
        g = contract["grid_summary"]; a = contract["availability_summary"]
        print(
            f"U2_STAGE0_SUCCESS_STOP exit=42 digest={contract_digest} candidates={g['candidates']} "
            f"ambient={g['ambient_branches']} strata={g['active_strata']} raw={g['raw_ragged_cardinality']} "
            f"collapsed={g['collapsed_cardinality']} promotion={g['promotion_contexts']} derived={a['DERIVED']} "
            f"unresolved={a['UNRESOLVED']} teeth={campaign['tooth_count']} controls={campaign['defect_absent_control_count']} "
            f"artifacts={self.artifact}", flush=True,
        )
        return SUCCESS_HALT


def validate_sources(repo: Path) -> int:
    source = repo / "software/em_charge_attribute"; python_files = sorted(source.glob("u2_*.py"))
    for path in python_files: ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    required = [
        source / "u2_stage0_sympy.py", source / "u2_stage0_dual.wl", source / "u2_stage0_compare.py",
        source / "u2_stage0_contract.py", source / "u2_stage0_mutations.py", source / "u2_code_closure_guard.py",
        source / "u2_containment.py", source / "u2_production_sympy.py", source / "u2_production_dual.wl",
        source / "u2_production_compare.py", source / "u2_production_mutations.py",
        source / "u2_production_runner.py", source / "u2_production_a9.py",
        source / "run_u2_boundary_adjudication.sh",
        Path("/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel"), Path("/usr/bin/bwrap"), Path("/usr/bin/strace"),
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing: raise StageFailure(f"runner prerequisites missing: {missing}")
    print(f"U2_RUNNER_VALID python_sources={len(python_files)} stage0_stages={len(STAGE_SEQUENCE)} production_stages=10"); return 0


def assert_production_launcher_at_anchor(repo: Path, anchor: str, path: Path) -> None:
    resolved = path.resolve(); rel = resolved.relative_to(repo).as_posix()
    result = subprocess.run(
        ["git", "cat-file", "blob", f"{anchor}:{rel}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise StageFailure(f"ASSERT_EVALUATED_CODE_CLOSURE production launcher absent from anchor: {rel}")
    if hashlib.sha256(resolved.read_bytes()).digest() != hashlib.sha256(result.stdout).digest():
        raise StageFailure(f"ASSERT_EXECUTED_BYTES_MATCH_ANCHOR production launcher differs: {rel}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument("--startup-contract-commit"); parser.add_argument("--stage0-contract-digest")
    parser.add_argument("--stage0", action="store_true"); parser.add_argument("--production", action="store_true")
    parser.add_argument(
        "--production-a9-preflight",
        choices=("arbiter", "fidelity", "adversarial_recompute", "read_only_review"),
    )
    parser.add_argument("--a9-output"); parser.add_argument("--production-self-probe", action="store_true")
    parser.add_argument("--validate-only", action="store_true"); parser.add_argument("--self-probe", action="store_true")
    args = parser.parse_args(); repo = Path(args.repo).resolve()
    source = repo / "software/em_charge_attribute"
    if args.production_self_probe:
        os.execv(
            "/usr/bin/python3",
            ["/usr/bin/python3", "-I", str(source / "u2_production_runner.py"), "--repo", str(repo), "--self-probe"],
        )
    if args.self_probe:
        output = repo / SCRATCH_REL / "runner_shell_probe.yaml"
        dump_yaml(output, {
            "schema_version": "U2_RUNNER_SHELL_PROBE_V1", "python_isolated": sys.flags.isolated == 1,
            "python_no_user_site": sys.flags.no_user_site == 1, "runner_path": str(Path(__file__).resolve()),
        }); print("U2_RUNNER_SHELL_PROBE_PASS"); return 0
    if args.validate_only: return validate_sources(repo)
    if args.production_a9_preflight:
        if not args.startup_contract_commit or not args.stage0_contract_digest:
            parser.error("A9 preflight requires both orchestrator-supplied anchor and stage0 digest")
        if not re.fullmatch(r"[0-9a-f]{40}", args.startup_contract_commit) or args.startup_contract_commit == "HEAD":
            raise StageFailure("A9 production anchor must be a full orchestrator-supplied commit, never HEAD")
        resolved = subprocess.run(
            ["git", "rev-parse", f"{args.startup_contract_commit}^{{commit}}"], cwd=repo,
            check=True, stdout=subprocess.PIPE, text=True,
        ).stdout.strip()
        if resolved != args.startup_contract_commit:
            raise StageFailure("A9 production anchor did not resolve identically")
        a9_script = source / "u2_production_a9.py"
        for path in (Path(__file__), source / "run_u2_boundary_adjudication.sh", a9_script):
            assert_production_launcher_at_anchor(repo, resolved, path)
        output = args.a9_output or str(
            source / "reports/u2_boundary_adjudication_artifacts/stage_1_production_v12"
            / f"a9_preflight_{args.production_a9_preflight}.yaml"
        )
        os.execv(
            "/usr/bin/python3",
            [
                "/usr/bin/python3", "-I", str(a9_script), "--repo", str(repo),
                "--leg", args.production_a9_preflight, "--startup-contract-commit", resolved,
                "--stage0-contract-digest", args.stage0_contract_digest, "--output", output,
            ],
        )
    if args.stage0 and args.production: parser.error("choose exactly one of --stage0 or --production")
    if args.production:
        if not args.startup_contract_commit or not args.stage0_contract_digest:
            parser.error("--production requires orchestrator-supplied anchor and stage0 digest")
        if not re.fullmatch(r"[0-9a-f]{40}", args.startup_contract_commit) or args.startup_contract_commit == "HEAD":
            raise StageFailure("production anchor must be a full orchestrator-supplied commit, never HEAD")
        resolved = subprocess.run(
            ["git", "rev-parse", f"{args.startup_contract_commit}^{{commit}}"], cwd=repo,
            check=True, stdout=subprocess.PIPE, text=True,
        ).stdout.strip()
        if resolved != args.startup_contract_commit:
            raise StageFailure("production anchor did not resolve identically")
        production_runner = source / "u2_production_runner.py"
        for path in (Path(__file__), source / "run_u2_boundary_adjudication.sh", production_runner):
            assert_production_launcher_at_anchor(repo, resolved, path)
        os.execv(
            "/usr/bin/python3",
            [
                "/usr/bin/python3", "-I", str(production_runner), "--repo", str(repo),
                "--startup-contract-commit", resolved, "--stage0-contract-digest", args.stage0_contract_digest,
            ],
        )
    if not args.stage0: parser.error("this launcher requires --stage0 or --production")
    if not args.startup_contract_commit or not re.fullmatch(r"[0-9a-f]{40}", args.startup_contract_commit) or args.startup_contract_commit == "HEAD":
        raise StageFailure("full orchestrator-supplied STARTUP_CONTRACT_COMMIT required, never HEAD")
    resolved = subprocess.run(
        ["git", "rev-parse", f"{args.startup_contract_commit}^{{commit}}"], cwd=repo,
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    ).stdout.strip()
    if resolved != args.startup_contract_commit:
        raise StageFailure(
            f"stage-0 anchor did not resolve identically supplied={args.startup_contract_commit} resolved={resolved}"
        )
    return Runner(repo, resolved).run()


if __name__ == "__main__":
    try: raise SystemExit(main())
    except StageFailure as failure:
        print(f"U2_STAGE0_BLOCKED {failure}", file=sys.stderr); raise SystemExit(1)
