#!/usr/bin/env python3
"""Staged and resumable production runner for U1 Phase C."""

from __future__ import annotations

import argparse
import ast
import hashlib
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml


RATIFIED_DIGEST = "83233baabd7f8e27c88d130b911691e76d01d5797da8eeb32c90bbae111ec95a"
STAGE0_ANCHOR = "377eab17a4babc12847450956dc55fe3e16d33da"
STAGE0_REL = Path(
    "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/"
    "stage_c_0_tilt_coupling_contract"
)
PRODUCTION_REL = Path(
    "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/"
    "stage_c_1_tilt_coupling_production"
)
SCRATCH_REL = Path("software/em_charge_attribute/_scratch/u1_phaseC_production")
STAGES = (
    "00_runner_shell_probe",
    "01_sympy_production",
    "02_wolfram_production",
    "03_dual_engine_compare",
    "04_mutation_campaign",
    "04b_mutation_closure_probe",
    "05_preliminary_code_closure",
    "06_preliminary_containment",
    "07_final_code_closure",
    "08_final_containment",
)


class RunnerFailure(RuntimeError):
    pass


class AmendmentHalt(RuntimeError):
    pass


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=140),
        encoding="utf-8",
    )


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_blob(repo: Path, anchor: str, rel: str) -> bytes:
    result = subprocess.run(
        ["git", "cat-file", "blob", f"{anchor}:{rel}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise RunnerFailure(f"ASSERT_EVALUATED_CODE_CLOSURE: {rel} absent from production anchor")
    return result.stdout


def assert_anchored(repo: Path, anchor: str, path: Path) -> None:
    resolved = path.resolve()
    try:
        rel = resolved.relative_to(repo).as_posix()
    except ValueError as failure:
        raise RunnerFailure(f"task executable is outside repository anchor: {resolved}") from failure
    blob = git_blob(repo, anchor, rel)
    if hashlib.sha256(resolved.read_bytes()).digest() != hashlib.sha256(blob).digest():
        raise RunnerFailure(f"ASSERT_EXECUTED_BYTES_MATCH_ANCHOR: {rel}")


class ProductionRunner:
    def __init__(self, repo: Path, anchor: str, stage0_digest: str):
        self.repo = repo
        self.anchor = anchor
        self.stage0_digest = stage0_digest
        self.source = repo / "software/em_charge_attribute"
        self.stage0 = repo / STAGE0_REL
        self.artifact = repo / PRODUCTION_REL
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
        assert_anchored(repo, anchor, Path(__file__))
        if self.state_path.exists():
            self.state = load_yaml(self.state_path)
            if self.state.get("startup_contract_commit") != anchor:
                raise RunnerFailure("resume state belongs to a different production anchor")
            if self.state.get("stage0_contract_digest") != stage0_digest:
                raise RunnerFailure("resume state belongs to a different stage-0 contract digest")
        else:
            self.state = {
                "schema_version": "U1_PHASE_C_PRODUCTION_RUNNER_STATE_V1",
                "startup_contract_commit": anchor,
                "stage0_contract_digest": stage0_digest,
                "terminal": "IN_PROGRESS",
                "completed_stages": {},
                "preflight_records": [],
            }
            self.save_state()
        self.records = load_yaml(self.records_path) if self.records_path.exists() else {
            "schema_version": "U1_PHASE_C_PRODUCTION_STAGE_RUN_RECORDS_V1",
            "records": [],
        }

    def save_state(self) -> None:
        dump_yaml(self.state_path, self.state)

    def save_records(self) -> None:
        dump_yaml(self.records_path, self.records)

    def trace(self, stage: str) -> Path:
        return self.scratch / "traces" / f"{stage}.strace"

    def preflight_trace(self, stage: str) -> Path:
        return self.scratch / "traces" / f"preflight_{stage}.strace"

    def stage_producers(self, *names: str) -> list[Path]:
        return [Path(__file__), *[self.source / name for name in names]]

    def producer_hashes(self, producers: list[Path]) -> dict[str, str]:
        hashes = {}
        for path in producers:
            assert_anchored(self.repo, self.anchor, path)
            rel = path.resolve().relative_to(self.repo).as_posix()
            hashes[rel] = sha256_path(path)
        return dict(sorted(hashes.items()))

    def completed(self, stage: str, required: list[Path], hashes: dict[str, str]) -> bool:
        record = self.state["completed_stages"].get(stage)
        if not record or record.get("status") != "PASS":
            return False
        if record.get("producer_script_sha256") != hashes:
            return False
        outputs = record.get("output_sha256", {})
        return all(path.is_file() and outputs.get(str(path)) == sha256_path(path) for path in required)

    def invalidate_from(self, stage: str) -> None:
        index = STAGES.index(stage)
        invalidated = [name for name in STAGES[index:] if name in self.state["completed_stages"]]
        for name in invalidated:
            del self.state["completed_stages"][name]
        if invalidated:
            self.state["terminal"] = "IN_PROGRESS"
            self.state.pop("exit_code", None)
            self.save_state()
            print(f"PHASEC_PRODUCTION_INVALIDATE from={stage} count={len(invalidated)}", flush=True)

    def sandbox(self, command: list[str]) -> list[str]:
        return [
            "/usr/bin/bwrap", "--ro-bind", "/", "/", "--dev-bind", "/dev", "/dev",
            "--proc", "/proc", "--unshare-net",
            "--bind", str(self.artifact), str(self.artifact),
            "--bind", str(self.scratch), str(self.scratch),
            "--bind", str(self.runtime_tmp), "/tmp",
            "--ro-bind", "/home/trevnorris/.Wolfram/Licensing", "/tmp/wolfram_userbase/Licensing",
            "--ro-bind", str(self.nsswitch), "/etc/nsswitch.conf",
            "--chdir", str(self.repo),
            "--setenv", "PYTHONNOUSERSITE", "1",
            "--setenv", "PYTHONDONTWRITEBYTECODE", "1",
            "--setenv", "TMPDIR", "/tmp",
            "--setenv", "WOLFRAM_USERBASE", "/tmp/wolfram_userbase",
            "--setenv", "WOLFRAM_BASE", "/tmp/wolfram_base",
            "--", *command,
        ]

    def traced(self, command: list[str], trace: Path) -> list[str]:
        trace.parent.mkdir(parents=True, exist_ok=True)
        return [
            "/usr/bin/strace", "-f", "-qq", "-s", "4096",
            "-e", "trace=%file,%network", "-o", str(trace),
            *self.sandbox(command),
        ]

    def execute(
        self, stage: str, command: list[str], trace: Path | None,
        timeout: int, live_progress: bool = False,
    ) -> tuple[int, str, float]:
        log = self.scratch / "logs" / f"{stage}.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        launched = self.traced(command, trace) if trace is not None else self.sandbox(command)
        started = time.monotonic()
        with log.open("wb") as handle:
            process = subprocess.Popen(launched, cwd=self.repo, stdout=handle, stderr=subprocess.STDOUT)
            next_marker = started + 15.0
            while process.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed > timeout:
                    process.kill()
                    process.wait()
                    raise RunnerFailure(f"stage {stage} exceeded {timeout}s; reformulation required")
                if live_progress and time.monotonic() >= next_marker:
                    print(f"PHASEC_WOLFRAM_LIVE_PROGRESS stage={stage} elapsed={elapsed:.1f}s", flush=True)
                    next_marker += 15.0
                time.sleep(0.25)
        elapsed = time.monotonic() - started
        return process.returncode, log.read_text(encoding="utf-8", errors="replace"), elapsed

    def reassert_stage0_and_environment(self, stage: str, producers: list[Path]) -> None:
        # This also asserts the working bytes of every imminent task executable
        # against the orchestrator-supplied production anchor before first use.
        self.producer_hashes([
            self.source / "u1_phaseC_stage0_contract.py", *producers,
        ])
        trace = self.preflight_trace(stage)
        command = [
            "/usr/bin/python3", "-I", str(self.source / "u1_phaseC_stage0_contract.py"),
            "--repo", str(self.repo), "--output-dir", str(self.stage0),
            "--verify", "--expected-digest", self.stage0_digest,
            "--startup-contract-commit", STAGE0_ANCHOR, "--recompute-environment",
        ]
        code, output, elapsed = self.execute(f"preflight_{stage}", command, trace, 600)
        if code != 0 or "PHASEC_STAGE0_CONTRACT_VERIFIED" not in output:
            # Environment changes are amendment-class HALTs, but the runner
            # itself does not manufacture exit 42; it records and stops for the
            # orchestrator to adjudicate.
            halt = {
                "schema_version": "U1_PHASE_C_PRODUCTION_AMENDMENT_HALT_V1",
                "status": "AMENDMENT_HALT",
                "exit_code": 42,
                "stage": stage,
                "startup_contract_commit": self.anchor,
                "stage0_contract_digest": self.stage0_digest,
                "reason": output[-2000:],
            }
            dump_yaml(self.artifact / "amendment_halt.yaml", halt)
            raise AmendmentHalt(f"preflight {stage}: {output[-800:]}")
        record = {
            "stage": stage, "status": "PASS", "elapsed_seconds": round(elapsed, 6),
            "stage0_contract_digest_recomputed": self.stage0_digest,
            "environment_identity_recomputed": True, "trace": str(trace),
            "trace_sha256": sha256_path(trace),
        }
        self.state["preflight_records"] = [
            row for row in self.state["preflight_records"] if row["stage"] != stage
        ] + [record]
        self.save_state()
        print(f"PHASEC_PRODUCTION_PREFLIGHT_PASS stage={stage}", flush=True)

    def run_stage(
        self, stage: str, command: list[str], required: list[Path], producers: list[Path],
        timeout: int = 600, evaluation: bool = False, live_progress: bool = False,
        outer_trace: bool = True,
    ) -> None:
        trace = self.trace(stage)
        all_required = [*required, *([trace] if outer_trace else [])]
        hashes = self.producer_hashes(producers)
        if self.completed(stage, all_required, hashes):
            print(f"PHASEC_PRODUCTION_RESUME_SKIP stage={stage}", flush=True)
            return
        self.invalidate_from(stage)
        if evaluation:
            self.reassert_stage0_and_environment(stage, producers)
        print(f"PHASEC_PRODUCTION_STAGE_START stage={stage}", flush=True)
        code, output, elapsed = self.execute(
            stage, command, trace if outer_trace else None, timeout, live_progress
        )
        self.records["records"].append({
            "stage": stage, "return_code": code, "elapsed_seconds": round(elapsed, 6),
            "output_tail": output[-4000:],
            "trace": str(trace) if outer_trace else None,
            "trace_sha256": sha256_path(trace) if outer_trace else None,
        })
        self.save_records()
        if code != 0:
            raise RunnerFailure(f"stage {stage} failed rc={code}: {output[-1600:]}")
        missing = [str(path) for path in all_required if not path.is_file()]
        if missing:
            raise RunnerFailure(f"stage {stage} omitted outputs: {missing}")
        self.state["completed_stages"][stage] = {
            "status": "PASS", "elapsed_seconds": round(elapsed, 6),
            "producer_script_sha256": hashes,
            "output_sha256": {str(path): sha256_path(path) for path in all_required},
        }
        self.save_state()
        print(f"PHASEC_PRODUCTION_STAGE_PASS stage={stage} elapsed={elapsed:.3f}s", flush=True)

    def python(self, script: str, *arguments: str) -> list[str]:
        return ["/usr/bin/python3", "-I", str(self.source / script), *arguments]

    def closure_arguments(self, traces: list[Path], output: Path) -> list[str]:
        arguments = [
            "--repo", str(self.repo), "--anchor", self.anchor,
            "--environment", str(self.stage0 / "environment_identity.yaml"),
        ]
        for trace in traces:
            arguments.extend(["--trace", str(trace)])
        arguments.extend(["--output", str(output)])
        return arguments

    def containment_arguments(self, traces: list[Path], output: Path) -> list[str]:
        arguments = [
            "--cwd", str(self.repo),
            "--allow-write-root", str(self.artifact),
            "--allow-write-root", str(self.scratch),
            "--allow-write-root", "/__phaseC_namespace_control__",
            "--allow-write-root", "/tmp/wolfram_userbase",
            "--allow-write-root", "/tmp/wolfram_base",
            "--mapped-write-root", f"/tmp={self.runtime_tmp}",
        ]
        for trace in traces:
            arguments.extend(["--trace", str(trace)])
        arguments.extend(["--output", str(output)])
        return arguments

    def run(self) -> int:
        sympy = self.artifact / "sympy_production.yaml"
        wolfram = self.artifact / "wolfram_production.yaml"
        agreement = self.artifact / "production_engine_agreement.yaml"
        results = self.artifact / "production_results.yaml"
        successors = self.artifact / "successor_records.yaml"
        summary = self.artifact / "production_summary.md"
        mutations = self.artifact / "mutation_campaign.yaml"
        preliminary_closure = self.scratch / "preliminary_code_closure.yaml"
        preliminary_containment = self.scratch / "preliminary_containment.yaml"
        final_closure = self.artifact / "evaluated_code_closure_evidence.yaml"
        final_containment = self.artifact / "containment_evidence.yaml"

        # §1.5 requires a fresh live assertion at every production resume,
        # independently of whether all later evaluation stages are resumable.
        self.reassert_stage0_and_environment("production_resume", [Path(__file__)])
        self.run_stage(
            "00_runner_shell_probe",
            ["/bin/bash", str(self.source / "run_u1_body_phaseC.sh"), "--production-shell-probe"],
            [self.scratch / "runner_shell_probe.yaml"],
            self.stage_producers("run_u1_body_phaseC.sh", "u1_phaseC_stage0_runner.py"),
        )
        self.run_stage(
            "01_sympy_production",
            self.python(
                "u1_phaseC_production_sympy.py", "--repo", str(self.repo),
                "--bundle-dir", str(self.stage0), "--stage0-contract-digest", self.stage0_digest,
                "--output", str(sympy),
            ),
            [sympy], self.stage_producers("u1_phaseC_production_sympy.py"), evaluation=True,
        )
        self.run_stage(
            "02_wolfram_production",
            [
                "/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel",
                "-noinit", "-noprompt", "-script", str(self.source / "u1_phaseC_production_dual.wl"),
                "--repo", str(self.repo), "--bundle-dir", str(self.stage0),
                "--stage0-contract-digest", self.stage0_digest, "--output", str(wolfram),
            ],
            [wolfram], self.stage_producers("u1_phaseC_production_dual.wl"),
            timeout=1200, evaluation=True, live_progress=True,
        )
        self.run_stage(
            "03_dual_engine_compare",
            self.python(
                "u1_phaseC_production_compare.py", "--sympy", str(sympy), "--wolfram", str(wolfram),
                "--bundle-dir", str(self.stage0), "--output", str(agreement),
                "--results-output", str(results), "--successors-output", str(successors),
                "--summary-output", str(summary),
            ),
            [agreement, results, successors, summary],
            self.stage_producers("u1_phaseC_production_compare.py"), evaluation=True,
        )
        mutation_scratch = self.scratch / "mutations"
        self.run_stage(
            "04_mutation_campaign",
            self.python(
                "u1_phaseC_production_mutations.py", "--repo", str(self.repo),
                "--sympy", str(sympy), "--wolfram", str(wolfram), "--bundle-dir", str(self.stage0),
                "--scratch", str(mutation_scratch), "--output", str(mutations),
                "--startup-contract-commit", self.anchor,
            ),
            [mutations], self.stage_producers(
                "u1_phaseC_production_mutations.py", "u1_phaseC_production_compare.py",
                "u1_phaseC_stage0_contract.py", "u1_phaseC_code_closure_guard.py",
                "u1_phaseC_containment.py",
            ), evaluation=True, outer_trace=False,
        )
        self.run_stage(
            "04b_mutation_closure_probe",
            self.python("u1_phaseC_production_mutations.py", "--closure-probe"),
            [], self.stage_producers("u1_phaseC_production_mutations.py"),
        )
        evaluation_traces = [
            self.preflight_trace("production_resume"),
            *[self.trace(name) for name in (
                "00_runner_shell_probe", "01_sympy_production", "02_wolfram_production",
                "03_dual_engine_compare", "04b_mutation_closure_probe",
            )],
            *[self.preflight_trace(name) for name in (
                "01_sympy_production", "02_wolfram_production",
                "03_dual_engine_compare", "04_mutation_campaign",
            )],
        ]
        self.run_stage(
            "05_preliminary_code_closure",
            self.python(
                "u1_phaseC_code_closure_guard.py",
                *self.closure_arguments(evaluation_traces, preliminary_closure),
            ),
            [preliminary_closure], self.stage_producers("u1_phaseC_code_closure_guard.py"),
        )
        traces_through_closure = [*evaluation_traces, self.trace("05_preliminary_code_closure")]
        self.run_stage(
            "06_preliminary_containment",
            self.python(
                "u1_phaseC_containment.py",
                *self.containment_arguments(traces_through_closure, preliminary_containment),
            ),
            [preliminary_containment], self.stage_producers("u1_phaseC_containment.py"),
        )
        traces_for_final_closure = [*traces_through_closure, self.trace("06_preliminary_containment")]
        self.run_stage(
            "07_final_code_closure",
            self.python(
                "u1_phaseC_code_closure_guard.py",
                *self.closure_arguments(traces_for_final_closure, final_closure),
            ),
            [final_closure], self.stage_producers("u1_phaseC_code_closure_guard.py"),
        )
        traces_for_final_containment = [*traces_for_final_closure, self.trace("07_final_code_closure")]
        self.run_stage(
            "08_final_containment",
            self.python(
                "u1_phaseC_containment.py",
                *self.containment_arguments(traces_for_final_containment, final_containment),
            ),
            [final_containment], self.stage_producers("u1_phaseC_containment.py"),
        )

        agreement_record = load_yaml(agreement)
        results_record = load_yaml(results)
        mutation_record = load_yaml(mutations)
        closure_record = load_yaml(final_closure)
        containment_record = load_yaml(final_containment)
        if agreement_record["status"] != "ENGINE_AGREE":
            raise RunnerFailure("engine agreement not green at terminal")
        if mutation_record["status"] != "PASS" or mutation_record["vacuous_case_count"] != 0:
            raise RunnerFailure("mutation armor not green at terminal")
        if closure_record["status"] != "PASS_PRODUCTION":
            raise RunnerFailure("evaluated-code closure not green at terminal")
        if containment_record["status"] != "PASS":
            raise RunnerFailure("process-tree containment not green at terminal")
        expected_preflights = {
            "production_resume", "01_sympy_production", "02_wolfram_production",
            "03_dual_engine_compare", "04_mutation_campaign",
        }
        observed_preflights = {row["stage"] for row in self.state["preflight_records"]}
        if observed_preflights != expected_preflights:
            raise RunnerFailure(
                f"preflight stage coverage mismatch: {sorted(observed_preflights)}"
            )
        terminal = {
            "schema_version": "U1_PHASE_C_PRODUCTION_TERMINAL_V1",
            "status": "SUCCESS_STOP",
            "exit_code": 0,
            "startup_contract_commit": self.anchor,
            "stage0_contract_digest": self.stage0_digest,
            "preflight_count": len(self.state["preflight_records"]),
            "preflights_all_green": observed_preflights == expected_preflights,
            "physics_headlines": {
                "tilt": results_record["tilt"]["headline"],
                "coupling": results_record["coupling_package"]["headline"],
            },
            "dual_engine": agreement_record["status"],
            "mutation_campaign": mutation_record["status"],
            "vacuous_mutations": mutation_record["vacuous_case_count"],
            "evaluated_code_closure": closure_record["status"],
            "containment": containment_record["status"],
            "network_attempt_count": containment_record["network_attempt_count"],
            "forbidden_write_attempt_count": containment_record["forbidden_write_attempt_count"],
            "A9_external_verification": "AWAITING_ORCHESTRATOR_EXTERNAL_LEGS_BEFORE_BANKING",
            "external_results_claimed": False,
            "artifact_sha256": {
                path.name: sha256_path(path) for path in (
                    agreement, results, successors, summary, mutations, final_closure, final_containment
                )
            },
        }
        dump_yaml(self.artifact / "production_terminal.yaml", terminal)
        self.state["terminal"] = "SUCCESS_STOP"
        self.state["exit_code"] = 0
        self.save_state()
        print(
            "PHASEC_PRODUCTION_SUCCESS_STOP exit=0 "
            f"tilt={results_record['tilt']['headline']} "
            f"coupling={results_record['coupling_package']['headline']}",
            flush=True,
        )
        return 0


def validate_sources(repo: Path) -> int:
    source = repo / "software/em_charge_attribute"
    python_names = (
        "u1_phaseC_production_sympy.py", "u1_phaseC_production_compare.py",
        "u1_phaseC_production_mutations.py", "u1_phaseC_production_runner.py",
        "u1_phaseC_production_a9.py",
        "u1_phaseC_stage0_runner.py", "u1_phaseC_stage0_contract.py",
        "u1_phaseC_code_closure_guard.py", "u1_phaseC_containment.py",
    )
    for name in python_names:
        path = source / name
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    required = [
        source / "u1_phaseC_production_dual.wl", source / "run_u1_body_phaseC.sh",
        Path("/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel"),
        Path("/usr/bin/bwrap"), Path("/usr/bin/strace"),
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RunnerFailure(f"production prerequisites missing: {missing}")
    print(f"PHASEC_PRODUCTION_RUNNER_VALID python={len(python_names)} stages={len(STAGES)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument("--startup-contract-commit")
    parser.add_argument("--stage0-contract-digest")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--self-probe", action="store_true")
    args = parser.parse_args()
    repo = Path(args.repo).resolve()
    if args.validate_only:
        return validate_sources(repo)
    if args.self_probe:
        output = repo / SCRATCH_REL / "runner_shell_probe.yaml"
        dump_yaml(output, {
            "schema_version": "U1_PHASE_C_PRODUCTION_RUNNER_SHELL_PROBE_V1",
            "python_isolated": sys.flags.isolated == 1,
            "python_no_user_site": sys.flags.no_user_site == 1,
            "runner_path": str(Path(__file__).resolve()),
        })
        print("PHASEC_PRODUCTION_RUNNER_SHELL_PROBE_PASS")
        return 0
    if not args.startup_contract_commit:
        raise RunnerFailure("STARTUP_CONTRACT_COMMIT is orchestrator-supplied and required")
    if not re.fullmatch(r"[0-9a-f]{40}", args.startup_contract_commit) or args.startup_contract_commit == "HEAD":
        raise RunnerFailure("production anchor must be a full orchestrator-supplied commit, never HEAD")
    resolved = subprocess.run(
        ["git", "rev-parse", f"{args.startup_contract_commit}^{{commit}}"], cwd=repo,
        check=True, stdout=subprocess.PIPE, text=True,
    ).stdout.strip()
    if resolved != args.startup_contract_commit:
        raise RunnerFailure("production anchor did not resolve identically")
    if args.stage0_contract_digest != RATIFIED_DIGEST:
        raise AmendmentHalt("orchestrator-supplied STAGE0_CONTRACT_DIGEST is missing or differs from ratified digest")
    return ProductionRunner(repo, resolved, args.stage0_contract_digest).run()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AmendmentHalt as failure:
        print(f"PHASEC_PRODUCTION_AMENDMENT_HALT {failure}", file=sys.stderr)
        raise SystemExit(42)
    except RunnerFailure as failure:
        print(f"PHASEC_PRODUCTION_BLOCKED {failure}", file=sys.stderr)
        raise SystemExit(1)
