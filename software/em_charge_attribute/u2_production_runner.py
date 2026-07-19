#!/usr/bin/env python3
"""Staged, resumable, anchor-authenticated U2 production runner."""

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


# The production launcher normally enters through an isolated interpreter, which
# ignores PYTHON* variables supplied by its caller.  Own the setting here too,
# and give every isolated child an explicit -B so the policy remains effective.
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True


RATIFIED_DIGEST = "b01a1821e908589c3698512bbb9aff874b721af6dcbfa1c3b8b1f8d33119b32b"
STAGE0_REL = Path("software/em_charge_attribute/reports/u2_boundary_adjudication_artifacts/stage_0_contract")
PRODUCTION_REL = Path("software/em_charge_attribute/reports/u2_boundary_adjudication_artifacts/stage_1_production_v12")
SCRATCH_REL = Path("software/em_charge_attribute/_scratch/u2_production")
STAGES = (
    "00_runner_shell_probe", "01_sympy_production", "02_wolfram_production",
    "03_dual_engine_compare", "04_mutation_campaign", "04b_mutation_closure_probe",
    "05_preliminary_code_closure", "06_preliminary_containment",
    "07_final_code_closure", "08_final_containment",
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
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=160), encoding="utf-8")


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def git_blob(repo: Path, anchor: str, rel: str) -> bytes:
    result = subprocess.run(
        ["git", "cat-file", "blob", f"{anchor}:{rel}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise RunnerFailure(f"ASSERT_EVALUATED_CODE_CLOSURE: {rel} absent from production anchor")
    return result.stdout


def assert_anchored(repo: Path, anchor: str, path: Path) -> str:
    try:
        rel = path.resolve().relative_to(repo).as_posix()
    except ValueError as failure:
        raise RunnerFailure(f"task executable outside repository anchor: {path}") from failure
    if hashlib.sha256(path.read_bytes()).digest() != hashlib.sha256(git_blob(repo, anchor, rel)).digest():
        raise RunnerFailure(f"ASSERT_EXECUTED_BYTES_MATCH_ANCHOR: {rel}")
    return rel


class ProductionRunner:
    def __init__(self, repo: Path, anchor: str, stage0_digest: str):
        self.repo = repo; self.anchor = anchor; self.stage0_digest = stage0_digest
        self.source = repo / "software/em_charge_attribute"
        self.stage0 = repo / STAGE0_REL; self.artifact = repo / PRODUCTION_REL
        self.scratch = repo / SCRATCH_REL; self.runtime_tmp = self.scratch / "runtime_tmp"
        self.state_path = self.scratch / "runner_state.yaml"
        self.records_path = self.scratch / "stage_run_records.yaml"
        for path in (
            self.artifact, self.scratch, self.runtime_tmp,
            self.scratch / "traces", self.scratch / "logs",
            self.runtime_tmp / "wolfram_userbase/Licensing", self.runtime_tmp / "wolfram_base",
        ):
            path.mkdir(parents=True, exist_ok=True)
        self.nsswitch = self.runtime_tmp / "u2_production_nsswitch.conf"
        self.nsswitch.write_text("hosts: files\n", encoding="utf-8")
        assert_anchored(repo, anchor, Path(__file__))
        if self.state_path.exists():
            self.state = load_yaml(self.state_path)
            if self.state.get("startup_contract_commit") != anchor:
                raise RunnerFailure("resume state belongs to a different production anchor")
            if self.state.get("stage0_contract_digest") != stage0_digest:
                raise AmendmentHalt("resume state belongs to a different stage-0 contract digest")
        else:
            self.state = {
                "schema_version": "U2_PRODUCTION_RUNNER_STATE_V1",
                "startup_contract_commit": anchor, "stage0_contract_digest": stage0_digest,
                "terminal": "IN_PROGRESS", "completed_stages": {}, "preflight_records": [],
            }
            self.save_state()
        self.records = load_yaml(self.records_path) if self.records_path.exists() else {
            "schema_version": "U2_PRODUCTION_STAGE_RUN_RECORDS_V1", "records": [],
        }

    def save_state(self) -> None:
        dump_yaml(self.state_path, self.state)

    def save_records(self) -> None:
        dump_yaml(self.records_path, self.records)

    def trace(self, stage: str) -> Path:
        return self.scratch / "traces" / f"{stage}.strace"

    def preflight_trace(self, stage: str) -> Path:
        return self.scratch / "traces" / f"preflight_{stage}.strace"

    def normalized_path(self, path: Path) -> str:
        resolved = path.resolve()
        try:
            return resolved.relative_to(self.repo).as_posix()
        except ValueError:
            return str(resolved)

    def file_set_digest(self, paths: list[Path], require_exists: bool = True) -> str:
        rows = []
        for path in sorted({value.resolve() for value in paths}, key=str):
            if not path.is_file():
                if require_exists:
                    raise RunnerFailure(f"stage file missing: {path}")
                rows.append({"path": self.normalized_path(path), "sha256": None})
            else:
                rows.append({"path": self.normalized_path(path), "sha256": sha256_path(path)})
        return canonical_digest(rows)

    def producer_record(self, producers: list[Path]) -> dict[str, Any]:
        paths = sorted(assert_anchored(self.repo, self.anchor, path) for path in producers)
        return {
            "producer_commit": self.anchor, "producer_paths": paths,
            "aggregate_producer_sha256": self.file_set_digest(producers),
            "per_script_content_map_present": False,
        }

    def stage_producers(self, *names: str) -> list[Path]:
        return [Path(__file__), *[self.source / name for name in names]]

    def completed(self, stage: str, required: list[Path], producer: dict[str, Any], input_digest: str) -> bool:
        row = self.state["completed_stages"].get(stage)
        if not row or row.get("status") != "PASS":
            return False
        return (
            row.get("producer") == producer and row.get("input_digest") == input_digest
            and row.get("output_digest") == self.file_set_digest(required, require_exists=False)
            and all(path.is_file() for path in required)
        )

    def invalidate_from(self, stage: str) -> None:
        index = STAGES.index(stage)
        invalidated = [name for name in STAGES[index:] if name in self.state["completed_stages"]]
        for name in invalidated:
            del self.state["completed_stages"][name]
        if invalidated:
            self.state["terminal"] = "IN_PROGRESS"
            self.state.pop("exit_code", None)
            self.save_state()
            print(f"U2_PRODUCTION_INVALIDATE from={stage} stages={','.join(invalidated)}", flush=True)

    def sandbox(self, command: list[str]) -> list[str]:
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

    def traced(self, command: list[str], trace: Path) -> list[str]:
        return [
            "/usr/bin/strace", "-f", "-qq", "-s", "4096", "-e", "trace=%file,%network",
            "-o", str(trace), *self.sandbox(command),
        ]

    def execute(self, stage: str, command: list[str], trace: Path | None, timeout: int, live: bool = False) -> tuple[int, str, float]:
        log = self.scratch / "logs" / f"{stage}.log"
        launched = self.traced(command, trace) if trace else self.sandbox(command)
        started = time.monotonic()
        with log.open("wb") as handle:
            process = subprocess.Popen(launched, cwd=self.repo, stdout=handle, stderr=subprocess.STDOUT)
            next_marker = started + 20
            while process.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed > timeout:
                    process.kill(); process.wait()
                    raise RunnerFailure(f"stage {stage} exceeded {timeout}s; reformulation required")
                if live and time.monotonic() >= next_marker:
                    print(f"U2_WOLFRAM_LIVE_PROGRESS stage={stage} elapsed={elapsed:.1f}s", flush=True)
                    next_marker += 20
                time.sleep(0.25)
        elapsed = time.monotonic() - started
        return process.returncode, log.read_text(encoding="utf-8", errors="replace"), elapsed

    def reassert_stage0_and_environment(self, stage: str, producers: list[Path]) -> None:
        preflight_producer = self.producer_record([self.source / "u2_stage0_contract.py", *producers])
        trace = self.preflight_trace(stage)
        command = [
            "/usr/bin/python3", "-I", "-B", str(self.source / "u2_stage0_contract.py"),
            "--repo", str(self.repo), "--output-dir", str(self.stage0),
            "--startup-contract-commit", self.anchor, "--verify",
            "--expected-digest", self.stage0_digest, "--recompute-environment",
        ]
        code, output, elapsed = self.execute(f"preflight_{stage}", command, trace, 600)
        if code != 0 or "U2_STAGE0_CONTRACT_VERIFIED" not in output:
            halt = {
                "schema_version": "U2_PRODUCTION_AMENDMENT_HALT_V1", "status": "AMENDMENT_HALT",
                "exit_code": 42, "stage": stage, "startup_contract_commit": self.anchor,
                "stage0_contract_digest": self.stage0_digest, "reason": output[-2400:],
            }
            dump_yaml(self.artifact / "amendment_halt.yaml", halt)
            raise AmendmentHalt(f"stage-0/environment preflight {stage} failed: {output[-1000:]}")
        row = {
            "stage": stage, "status": "PASS", "elapsed_seconds": round(elapsed, 6),
            "stage0_contract_digest_recomputed": self.stage0_digest,
            "environment_identity_recomputed": True, "trace": str(trace), "trace_sha256": sha256_path(trace),
            "producer": preflight_producer,
        }
        self.state["preflight_records"] = [value for value in self.state["preflight_records"] if value["stage"] != stage] + [row]
        self.save_state()
        print(f"U2_PRODUCTION_PREFLIGHT_PASS stage={stage}", flush=True)

    def run_stage(
        self, stage: str, command: list[str], outputs: list[Path], producers: list[Path],
        inputs: list[Path] | None = None, timeout: int = 600, evaluation: bool = True,
        live: bool = False, outer_trace: bool = True,
    ) -> None:
        trace = self.trace(stage); required = [*outputs, *([trace] if outer_trace else [])]
        producer = self.producer_record(producers)
        input_digest = canonical_digest({
            "startup_contract_commit": self.anchor, "stage0_contract_digest": self.stage0_digest,
            "file_set_sha256": self.file_set_digest(inputs or []),
            "aggregate_producer_sha256": producer["aggregate_producer_sha256"],
        })
        if self.completed(stage, required, producer, input_digest):
            print(f"U2_PRODUCTION_RESUME_SKIP stage={stage}", flush=True)
            return
        self.invalidate_from(stage)
        if evaluation:
            self.reassert_stage0_and_environment(stage, producers)
        print(f"U2_PRODUCTION_STAGE_START stage={stage}", flush=True)
        code, output, elapsed = self.execute(stage, command, trace if outer_trace else None, timeout, live)
        self.records["records"].append({
            "stage": stage, "return_code": code, "elapsed_seconds": round(elapsed, 6),
            "output_tail": output[-3000:], "trace": str(trace) if outer_trace else None,
            "trace_sha256": sha256_path(trace) if outer_trace and trace.is_file() else None,
        })
        self.save_records()
        if code != 0:
            raise RunnerFailure(f"stage {stage} failed rc={code}: {output[-1800:]}")
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise RunnerFailure(f"stage {stage} omitted outputs: {missing}")
        self.state["completed_stages"][stage] = {
            "status": "PASS", "elapsed_seconds": round(elapsed, 6), "producer": producer,
            "input_digest": input_digest, "output_digest": self.file_set_digest(required),
        }
        self.save_state()
        print(f"U2_PRODUCTION_STAGE_PASS stage={stage} elapsed={elapsed:.3f}s", flush=True)

    def python(self, script: str, *args: str) -> list[str]:
        return ["/usr/bin/python3", "-I", "-B", str(self.source / script), *args]

    def closure_args(self, traces: list[Path], output: Path) -> list[str]:
        values = [
            "--repo", str(self.repo), "--anchor", self.anchor,
            "--environment", str(self.stage0 / "environment_identity.yaml"),
        ]
        for trace in traces:
            values.extend(["--trace", str(trace)])
        return [*values, "--output", str(output)]

    def containment_args(self, traces: list[Path], output: Path) -> list[str]:
        values = [
            "--cwd", str(self.repo), "--allow-write-root", str(self.artifact),
            "--allow-write-root", str(self.scratch), "--allow-write-root", "/__u2_namespace_control__",
            "--allow-write-root", "/tmp/wolfram_userbase", "--allow-write-root", "/tmp/wolfram_base",
            "--mapped-write-root", f"/tmp={self.runtime_tmp}",
        ]
        for trace in traces:
            values.extend(["--trace", str(trace)])
        return [*values, "--output", str(output)]

    def run(self) -> int:
        sympy = self.artifact / "sympy_production.yaml"; wolfram = self.artifact / "wolfram_production.yaml"
        agreement = self.artifact / "production_engine_agreement.yaml"; results = self.artifact / "production_results.yaml"
        summary = self.artifact / "production_summary.md"; mutations = self.artifact / "mutation_campaign.yaml"
        preliminary_closure = self.scratch / "preliminary_code_closure.yaml"
        preliminary_containment = self.scratch / "preliminary_containment.yaml"
        final_closure = self.artifact / "evaluated_code_closure_evidence.yaml"
        final_containment = self.artifact / "containment_evidence.yaml"

        self.reassert_stage0_and_environment("production_resume", [Path(__file__)])
        self.run_stage(
            "00_runner_shell_probe", ["/bin/bash", str(self.source / "run_u2_boundary_adjudication.sh"), "--production-shell-probe"],
            [self.scratch / "runner_shell_probe.yaml"],
            self.stage_producers("run_u2_boundary_adjudication.sh", "u2_stage0_runner.py"), evaluation=False,
        )
        self.run_stage(
            "01_sympy_production", self.python(
                "u2_production_sympy.py", "--repo", str(self.repo), "--bundle-dir", str(self.stage0),
                "--stage0-contract-digest", self.stage0_digest, "--output", str(sympy),
            ), [sympy], self.stage_producers("u2_production_sympy.py", "u2_stage0_sympy.py"),
            inputs=[self.stage0 / "stage0_bundle.yaml"],
        )
        self.run_stage(
            "02_wolfram_production", [
                "/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel", "-noinit", "-noprompt", "-script",
                str(self.source / "u2_production_dual.wl"), "--repo", str(self.repo), "--bundle-dir", str(self.stage0),
                "--stage0-contract-digest", self.stage0_digest, "--output", str(wolfram),
            ], [wolfram], self.stage_producers("u2_production_dual.wl"),
            inputs=[self.stage0 / "stage0_bundle.yaml"], timeout=1200, live=True,
        )
        self.run_stage(
            "03_dual_engine_compare", self.python(
                "u2_production_compare.py", "--sympy", str(sympy), "--wolfram", str(wolfram),
                "--bundle-dir", str(self.stage0), "--output", str(agreement),
                "--results-output", str(results), "--summary-output", str(summary),
            ), [agreement, results, summary],
            self.stage_producers("u2_production_compare.py", "u2_production_sympy.py", "u2_stage0_sympy.py"),
            inputs=[sympy, wolfram],
        )
        mutation_scratch = self.scratch / "mutations"
        self.run_stage(
            "04_mutation_campaign", self.python(
                "u2_production_mutations.py", "--repo", str(self.repo), "--sympy", str(sympy),
                "--wolfram", str(wolfram), "--bundle-dir", str(self.stage0), "--scratch", str(mutation_scratch),
                "--output", str(mutations), "--startup-contract-commit", self.anchor,
            ), [mutations], self.stage_producers(
                "u2_production_mutations.py", "u2_production_compare.py", "u2_production_sympy.py",
                "u2_stage0_sympy.py", "u2_stage0_contract.py", "u2_code_closure_guard.py", "u2_containment.py",
            ), inputs=[sympy, wolfram, agreement], outer_trace=False,
        )
        self.run_stage(
            "04b_mutation_closure_probe", self.python("u2_production_mutations.py", "--closure-probe"),
            [], self.stage_producers("u2_production_mutations.py"), inputs=[mutations],
        )

        stage_names = (
            "00_runner_shell_probe", "01_sympy_production", "02_wolfram_production",
            "03_dual_engine_compare", "04b_mutation_closure_probe",
        )
        preflight_names = (
            "production_resume", "01_sympy_production", "02_wolfram_production",
            "03_dual_engine_compare", "04_mutation_campaign", "04b_mutation_closure_probe",
        )
        evaluation_traces = [*[self.trace(name) for name in stage_names], *[self.preflight_trace(name) for name in preflight_names]]
        self.run_stage(
            "05_preliminary_code_closure", self.python(
                "u2_code_closure_guard.py", *self.closure_args(evaluation_traces, preliminary_closure),
            ), [preliminary_closure], self.stage_producers("u2_code_closure_guard.py"), inputs=evaluation_traces,
        )
        traces_closure = [*evaluation_traces, self.preflight_trace("05_preliminary_code_closure"), self.trace("05_preliminary_code_closure")]
        self.run_stage(
            "06_preliminary_containment", self.python(
                "u2_containment.py", *self.containment_args(traces_closure, preliminary_containment),
            ), [preliminary_containment], self.stage_producers("u2_containment.py"), inputs=traces_closure,
        )
        traces_containment = [*traces_closure, self.preflight_trace("06_preliminary_containment"), self.trace("06_preliminary_containment")]
        self.run_stage(
            "07_final_code_closure", self.python(
                "u2_code_closure_guard.py", *self.closure_args(traces_containment, final_closure),
            ), [final_closure], self.stage_producers("u2_code_closure_guard.py"), inputs=traces_containment,
        )
        traces_final = [*traces_containment, self.preflight_trace("07_final_code_closure"), self.trace("07_final_code_closure")]
        self.run_stage(
            "08_final_containment", self.python(
                "u2_containment.py", *self.containment_args(traces_final, final_containment),
            ), [final_containment], self.stage_producers("u2_containment.py"), inputs=traces_final,
        )

        agreement_row = load_yaml(agreement); results_row = load_yaml(results); mutation_row = load_yaml(mutations)
        closure_row = load_yaml(final_closure); containment_row = load_yaml(final_containment)
        if agreement_row["status"] != "ENGINE_AGREE": raise RunnerFailure("dual-engine agreement is not green")
        if mutation_row["status"] != "PASS" or mutation_row["vacuous_case_count"]: raise RunnerFailure("mutation armor is not green")
        if closure_row["status"] != "PASS_PRODUCTION": raise RunnerFailure("evaluated-code closure is not green")
        if containment_row["status"] != "PASS": raise RunnerFailure("process-tree containment is not green")
        headlines = results_row["headlines"]
        if headlines["integrity_failures"]:
            raise RunnerFailure("zero-integrity banking gate failed")
        expected_preflights = {
            "production_resume", "01_sympy_production", "02_wolfram_production", "03_dual_engine_compare",
            "04_mutation_campaign", "04b_mutation_closure_probe", "05_preliminary_code_closure",
            "06_preliminary_containment", "07_final_code_closure", "08_final_containment",
        }
        observed = {row["stage"] for row in self.state["preflight_records"]}
        if observed != expected_preflights:
            raise RunnerFailure(f"preflight coverage mismatch: {sorted(observed)}")
        terminal = {
            "schema_version": "U2_PRODUCTION_TERMINAL_V1", "status": "SUCCESS_STOP", "exit_code": 0,
            "startup_contract_commit": self.anchor, "stage0_contract_digest": self.stage0_digest,
            "preflight_count": len(observed), "preflights_all_green": True,
            "headlines": headlines, "dual_engine": agreement_row["status"],
            "mutation_campaign": mutation_row["status"], "mutation_teeth": mutation_row["tooth_count"],
            "vacuous_mutations": mutation_row["vacuous_case_count"],
            "evaluated_code_closure": closure_row["status"], "containment": containment_row["status"],
            "network_attempt_count": containment_row["network_attempt_count"],
            "forbidden_write_attempt_count": containment_row["forbidden_write_attempt_count"],
            "A9_external_verification": "AWAITING_ORCHESTRATOR_FOUR_LEGS_BEFORE_BANKING",
            "external_results_claimed": False,
            "artifact_sha256": {
                path.name: sha256_path(path) for path in (
                    agreement, results, summary, mutations, final_closure, final_containment
                )
            },
        }
        dump_yaml(self.artifact / "production_terminal.yaml", terminal)
        self.state.update({"terminal": "SUCCESS_STOP", "exit_code": 0}); self.save_state()
        print(
            f"U2_PRODUCTION_SUCCESS_STOP exit=0 dispositions={headlines['dispositions']} "
            f"ensemble={headlines['ensemble_level_1_final']} topology={headlines['topology_gate']} "
            f"templates={headlines['posed_template_count']} promotions={headlines['promotion_outcomes']}",
            flush=True,
        )
        return 0


def validate_sources(repo: Path) -> int:
    source = repo / "software/em_charge_attribute"
    python_names = (
        "u2_production_sympy.py", "u2_production_compare.py", "u2_production_mutations.py",
        "u2_production_runner.py", "u2_production_a9.py", "u2_stage0_runner.py",
        "u2_stage0_contract.py", "u2_stage0_sympy.py", "u2_code_closure_guard.py", "u2_containment.py",
    )
    for name in python_names:
        path = source / name
        if not path.is_file(): raise RunnerFailure(f"production source missing: {name}")
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    required = [
        source / "u2_production_dual.wl", source / "run_u2_boundary_adjudication.sh",
        Path("/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel"),
        Path("/usr/bin/bwrap"), Path("/usr/bin/strace"),
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing: raise RunnerFailure(f"production prerequisites missing: {missing}")
    print(f"U2_PRODUCTION_RUNNER_VALID python={len(python_names)} stages={len(STAGES)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument("--startup-contract-commit"); parser.add_argument("--stage0-contract-digest")
    parser.add_argument("--validate-only", action="store_true"); parser.add_argument("--self-probe", action="store_true")
    args = parser.parse_args(); repo = Path(args.repo).resolve()
    if args.validate_only: return validate_sources(repo)
    if args.self_probe:
        output = repo / SCRATCH_REL / "runner_shell_probe.yaml"
        dump_yaml(output, {
            "schema_version": "U2_PRODUCTION_RUNNER_SHELL_PROBE_V1",
            "python_isolated": sys.flags.isolated == 1, "python_no_user_site": sys.flags.no_user_site == 1,
            "runner_path": str(Path(__file__).resolve()),
        })
        print("U2_PRODUCTION_RUNNER_SHELL_PROBE_PASS"); return 0
    anchor = args.startup_contract_commit
    if not anchor or not re.fullmatch(r"[0-9a-f]{40}", anchor) or anchor == "HEAD":
        raise RunnerFailure("full orchestrator-supplied STARTUP_CONTRACT_COMMIT required, never HEAD")
    resolved = subprocess.run(
        ["git", "rev-parse", f"{anchor}^{{commit}}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if resolved.returncode != 0 or resolved.stdout.strip() != anchor:
        raise RunnerFailure("production anchor did not resolve identically")
    if args.stage0_contract_digest != RATIFIED_DIGEST:
        raise AmendmentHalt("orchestrator-supplied STAGE0_CONTRACT_DIGEST differs from ratified digest")
    return ProductionRunner(repo, anchor, args.stage0_contract_digest).run()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AmendmentHalt as failure:
        print(f"U2_PRODUCTION_AMENDMENT_HALT {failure}", file=sys.stderr)
        raise SystemExit(42)
    except RunnerFailure as failure:
        print(f"U2_PRODUCTION_BLOCKED {failure}", file=sys.stderr)
        raise SystemExit(1)
