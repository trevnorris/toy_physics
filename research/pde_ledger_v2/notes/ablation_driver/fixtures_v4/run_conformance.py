#!/usr/bin/env python3
"""Run the frozen A1--A9 suite against the contract's public executable."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import os
import pathlib
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
from collections import Counter
from typing import Any, Iterator


sys.dont_write_bytecode = True

SUITE = pathlib.Path(__file__).resolve().parent
PROJECT_ROOT = SUITE.parents[4]
CONTRACT = SUITE.parent / "CONTRACT.md"
DRIVER_REL = pathlib.Path("research/pde_ledger_v2/scripts/ablation_driver.py")
DRIVER = PROJECT_ROOT / DRIVER_REL
CHILD_REL = pathlib.Path(
    "research/pde_ledger_v2/notes/ablation_driver/fixtures_v4/child.py"
)
LEGACY_DIR = PROJECT_ROOT / "research/pde_ledger_v2/notes/stage023_step_h_evidence"
STAGE_TARGET_REL = pathlib.Path(
    "research/pde_ledger_v2/scripts/"
    "ledger_stage023_nullspace_underdetermination_sympy_audit.py"
)
STAGE_ARTIFACT_REL = pathlib.Path(
    "research/pde_ledger_v2/scripts/"
    "ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt"
)
ALL_ITEMS = [f"A{number}" for number in range(1, 10)]
A2_COMPLETION_WAIT_SECONDS = 30.0
A2_POLL_SECONDS = 0.002
A7_PROGRESS_POLL_SECONDS = 0.2
A7_PROGRESS_GRACE_SECONDS = 60
KILL_FALLBACK_WAIT_SECONDS = 5.0
AMBIENT_POISON_NAME = "FIXTURES_V4_UNDECLARED_AMBIENT"
ARGV_PROBE = "fixture; exit 97"
INVOCATION_PROBE_COMMAND = "exec"
INVOCATION_PROBE_REL = pathlib.Path(
    "research/pde_ledger_v2/notes/ablation_driver/fixtures_v4/exec"
)
STAGE_MIRROR_INPUTS = (
    CHILD_REL,
    STAGE_TARGET_REL,
    STAGE_ARTIFACT_REL,
    pathlib.Path("research/pde_ledger_v2/scripts/check_ledger_dimensions_pin.py"),
    pathlib.Path("research/pde_ledger_v2/scripts/compare_dimension_artifacts.py"),
    pathlib.Path("research/pde_ledger_v2/scripts/ledger_dimensions.accepted.sha256"),
    pathlib.Path("research/pde_ledger_v2/scripts/ledger_dimensions.py"),
    pathlib.Path(
        "research/pde_ledger_v2/mathematica/out/"
        "ledger_stage023_nullspace_underdetermination_mathematica_audit.out"
    ),
    pathlib.Path("research/pde_ledger_v2/notes/stage023_step_h_evidence/include_list.tsv"),
)


def rel(path: pathlib.Path, project_root: pathlib.Path = PROJECT_ROOT) -> str:
    return path.resolve().relative_to(project_root.resolve()).as_posix()


def write_jcs(path: pathlib.Path, value: Any) -> None:
    path.write_text(oracle.canonical_json(value) + "\n", encoding="utf-8", newline="\n")


def driver_environment() -> dict[str, str]:
    environment = dict(os.environ)
    environment[AMBIENT_POISON_NAME] = "must-not-reach-configured-children"
    environment["HOME"] = "/fixture-ambient-home"
    environment["USER"] = "fixture-ambient-user"
    environment["LOGNAME"] = "fixture-ambient-logname"
    environment["LANG"] = "fixture-ambient-lang"
    environment["SHELL"] = "/fixture-ambient-shell"
    environment["TMPDIR"] = "/fixture-ambient-tmp"
    return environment


def committed_bytes(path: pathlib.Path) -> bytes:
    path_text = rel(path)
    process = subprocess.run(
        ["git", "show", f"HEAD:{path_text}"],
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        check=False,
    )
    oracle.require(process.returncode == 0, f"cannot read committed oracle {path_text}")
    return process.stdout


@contextlib.contextmanager
def runtime(label: str) -> Iterator[pathlib.Path]:
    scratch = PROJECT_ROOT / "research/pde_ledger_v2/_scratch"
    path = pathlib.Path(tempfile.mkdtemp(prefix=f"fixtures-v4-{label}-", dir=scratch))
    try:
        yield path
    finally:
        shutil.rmtree(path, ignore_errors=True)


def copy_fixture(item: str, destination: pathlib.Path) -> None:
    source = SUITE / "fixtures" / item
    for path in source.iterdir():
        if path.is_file():
            shutil.copyfile(path, destination / path.name)


def synthetic_config(
    run_dir: pathlib.Path,
    case: str,
    *,
    artifact: bool,
    artifact_entry: bool = False,
    artifact_format: str = "none",
    producer_format: str = "exit-only",
    checker_format: str = "exit-only",
) -> tuple[pathlib.Path, dict[str, Any]]:
    target = run_dir / "target.txt"
    include = run_dir / "include_list.tsv"
    artifact_path = run_dir / "artifact.txt"
    if artifact and not artifact_entry:
        artifact_path.unlink(missing_ok=True)
    target_rel = rel(target)
    artifact_arg = rel(artifact_path) if artifact else "-"
    reset_path = run_dir / "cache"
    child_argv_tail = [
        target_rel,
        artifact_arg,
        rel(reset_path),
        AMBIENT_POISON_NAME,
        ARGV_PROBE,
    ]
    config: dict[str, Any] = {
        "artifacts": (
            [
                {
                    "name": "probe",
                    "observation_format": artifact_format,
                    "path": rel(artifact_path),
                }
            ]
            if artifact
            else []
        ),
        "checkers": [
            {
                "argv": [
                    INVOCATION_PROBE_COMMAND,
                    "checker",
                    case,
                    *child_argv_tail,
                ],
                "name": "probe_checker",
                "report_format": checker_format,
                "timeout_seconds": 10,
            }
        ],
        "contract": "ablation-driver-v1",
        "environment": {
            "LC_ALL": "C",
            "PATH": f"{SUITE}:/usr/local/bin:/usr/bin:/bin",
            "PYTHONDONTWRITEBYTECODE": "1",
            "TZ": "UTC",
        },
        "evidence_path": rel(run_dir / "evidence"),
        "include_list": rel(include),
        "producer": {
            "argv": [
                INVOCATION_PROBE_COMMAND,
                "producer",
                case,
                *child_argv_tail,
            ],
            "report_format": producer_format,
            "timeout_seconds": 10,
        },
        "recovery_path": rel(run_dir / "recovery"),
        "reset_paths": [rel(reset_path)],
        "stable_inputs": [CHILD_REL.as_posix(), INVOCATION_PROBE_REL.as_posix()],
        "targets": [{"name": "stage", "path": target_rel}],
        "workdir": ".",
    }
    config_path = run_dir / "config.json"
    write_jcs(config_path, config)
    return config_path, config


def parse_success(process: subprocess.CompletedProcess[bytes], operation: str) -> dict[str, Any]:
    oracle.require(process.returncode == 0, f"{operation}: exit {process.returncode}: {process.stderr!r}")
    oracle.require(process.stderr == b"", f"{operation}: success wrote stderr")
    value = oracle.parse_standalone_jcs(process.stdout, f"{operation} stdout")
    oracle.require(isinstance(value, dict), f"{operation}: summary is not an object")
    return value


def invoke(
    operation: str,
    config_path: pathlib.Path,
    *,
    project_root: pathlib.Path = PROJECT_ROOT,
    timeout: float = 590,
) -> tuple[dict[str, Any], subprocess.CompletedProcess[bytes]]:
    process = subprocess.run(
        [DRIVER_REL.as_posix(), operation, "--config", rel(config_path, project_root)],
        cwd=project_root,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=driver_environment(),
        timeout=timeout,
        check=False,
    )
    return parse_success(process, operation), process


def start_run(
    config_path: pathlib.Path,
    project_root: pathlib.Path = PROJECT_ROOT,
) -> subprocess.Popen[bytes]:
    return subprocess.Popen(
        [DRIVER_REL.as_posix(), "run", "--config", rel(config_path, project_root)],
        cwd=project_root,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=driver_environment(),
        start_new_session=True,
    )


def wait_for(predicate: Any, message: str, seconds: float = 20.0) -> None:
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.01)
    raise oracle.FixtureFailure(message)


def kill_leftover_group(group: int) -> None:
    try:
        os.killpg(group, signal.SIGKILL)
    except ProcessLookupError:
        pass


def close_process_pipes(process: subprocess.Popen[bytes]) -> None:
    for stream in (process.stdout, process.stderr):
        if stream is not None:
            stream.close()


def finish_process(
    process: subprocess.Popen[bytes],
    seconds: float = 30.0,
) -> subprocess.CompletedProcess[bytes]:
    try:
        stdout, stderr = process.communicate(timeout=seconds)
    except subprocess.TimeoutExpired as wait_failure:
        kill_leftover_group(process.pid)
        try:
            stdout, stderr = process.communicate(timeout=KILL_FALLBACK_WAIT_SECONDS)
        except subprocess.TimeoutExpired as kill_failure:
            close_process_pipes(process)
            raise oracle.FixtureFailure(
                "fixture child process did not close after its bounded SIGKILL fallback"
            ) from kill_failure
        raise oracle.FixtureFailure(
            "fixture child process exceeded its bounded wait"
        ) from wait_failure
    return subprocess.CompletedProcess(process.args, process.returncode, stdout, stderr)


def capture_invocation(
    project_root: pathlib.Path,
    config_path: pathlib.Path,
    config: dict[str, Any],
) -> tuple[bytes, bytes, oracle.EntrySnapshot]:
    config_data = config_path.read_bytes()
    include_data = (project_root / config["include_list"]).read_bytes()
    entry_snapshot = oracle.capture_entry_snapshot(project_root, config)
    return config_data, include_data, entry_snapshot


def audit_child_observed_mutants(
    evidence: pathlib.Path,
    rows: list[dict[str, Any]],
) -> None:
    for row in rows:
        applied = row["applied"] if isinstance(row["applied"], bool) else row["applied"] == "true"
        if not applied:
            continue
        captures = row["captures"]
        if isinstance(captures, str):
            captures = oracle.parse_jcs_cell(captures, "fixture child captures")
        reference = next(
            (item for item in captures if item.get("role") == "producer.stdout"),
            None,
        )
        oracle.require(reference is not None, "fixture child producer stdout is absent")
        data = (evidence / reference["path"]).read_bytes()
        try:
            lines = data.decode("utf-8").splitlines()
        except UnicodeError as exc:
            raise oracle.FixtureFailure("fixture child producer stdout is not UTF-8") from exc
        observed = [
            line.removeprefix("FIXTURE_TARGET_SHA256=")
            for line in lines
            if line.startswith("FIXTURE_TARGET_SHA256=")
        ]
        oracle.require(len(observed) == 1, "fixture child target observation is not singular")
        oracle.require(
            observed[0] == row["target_mutant_sha256"],
            "fixture child observed bytes other than the exact C-3 mutant",
        )


def audit_run_summary(
    summary: dict[str, Any],
    evidence: pathlib.Path,
    rows: list[dict[str, Any]],
    project_root: pathlib.Path = PROJECT_ROOT,
) -> None:
    oracle.require(
        set(summary) == {"schema", "evidence", "rows", "outcome_counts", "restored"},
        "run summary member set",
    )
    oracle.require(summary["schema"] == "ablation-run-summary-v1", "run summary schema")
    oracle.require(summary["evidence"] == rel(evidence, project_root), "run summary evidence")
    oracle.require(summary["rows"] == len(rows), "run summary row count")
    counts = Counter(row["outcome"] for row in rows)
    oracle.require(
        summary["outcome_counts"] == {token: counts[token] for token in sorted(oracle.OUTCOMES)},
        "run summary outcome counts",
    )
    oracle.require(summary["restored"] is True, "run summary restoration")


def completed(
    project_root: pathlib.Path,
    config: dict[str, Any],
    summary: dict[str, Any],
    config_data: bytes,
    include_data: bytes,
    entry_snapshot: oracle.EntrySnapshot,
    *,
    required_restore_operation: str = "run",
) -> list[dict[str, Any]]:
    evidence = project_root / config["evidence_path"]
    rows = oracle.audit_completed_evidence(
        project_root,
        evidence,
        CONTRACT,
        entry_snapshot,
        invocation_config=config_data,
        invocation_include=include_data,
        required_restore_operation=required_restore_operation,
    )
    audit_run_summary(summary, evidence, rows, project_root)
    return rows


def test_a1() -> None:
    with runtime("a1") as run_dir:
        copy_fixture("a1", run_dir)
        config_path, config = synthetic_config(
            run_dir,
            "a1",
            artifact=True,
            checker_format="ledger-result-v1",
        )
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        summary, _ = invoke("run", config_path)
        evidence = PROJECT_ROOT / config["evidence_path"]
        rows = completed(
            PROJECT_ROOT,
            config,
            summary,
            config_data,
            include_data,
            entry_snapshot,
        )
        audit_child_observed_mutants(evidence, rows)
        oracle.grade_a1(rows)
        oracle.require(not (run_dir / "cache").exists(), "A1 reset path survived restoration")


def test_a2() -> None:
    with runtime("a2") as run_dir:
        copy_fixture("a2", run_dir)
        config_path, config = synthetic_config(run_dir, "a2", artifact=True)
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        process = start_run(config_path)
        artifact = run_dir / "artifact.txt"
        producer_seen = False
        deadline = time.monotonic() + A2_COMPLETION_WAIT_SECONDS
        while process.poll() is None:
            producer_seen = producer_seen or artifact.exists()
            if time.monotonic() >= deadline:
                try:
                    stop_process_group(process, item="A2")
                except Exception as shutdown_failure:
                    print(
                        "A2 shutdown FAIL "
                        f"{type(shutdown_failure).__name__}: {shutdown_failure}"
                    )
                raise oracle.FixtureFailure("A2 run exceeded its bounded completion wait")
            time.sleep(A2_POLL_SECONDS)
        finished = finish_process(process)
        producer_seen = producer_seen or artifact.exists()
        summary = parse_success(finished, "run")
        evidence = PROJECT_ROOT / config["evidence_path"]
        rows = completed(
            PROJECT_ROOT,
            config,
            summary,
            config_data,
            include_data,
            entry_snapshot,
        )
        oracle.grade_a2(rows)
        oracle.require(not producer_seen, "A2 producer was invoked for an unusable row")


def partial_keys(evidence: pathlib.Path) -> list[str]:
    header, body = oracle.read_tsv_bytes((evidence / "results.tsv").read_bytes(), "partial results")
    key_index = header.index("key")
    return [row[key_index] for row in body]


def banked_result_tree(evidence: pathlib.Path) -> dict[str, bytes]:
    paths = [evidence / "results.tsv"]
    captures = evidence / "captures"
    if captures.is_dir():
        paths.extend(path for path in captures.rglob("*") if path.is_file())
    return {
        path.relative_to(evidence).as_posix(): path.read_bytes()
        for path in sorted(paths, key=lambda item: item.as_posix().encode())
    }


def test_a3() -> None:
    with runtime("a3-reference") as reference_dir, runtime("a3-interrupted") as resumed_dir:
        copy_fixture("a3", reference_dir)
        reference_config_path, reference_config = synthetic_config(
            reference_dir, "a3", artifact=False
        )
        reference_config_data, reference_include_data, reference_entry = capture_invocation(
            PROJECT_ROOT,
            reference_config_path,
            reference_config,
        )
        reference_summary, _ = invoke("run", reference_config_path)
        reference_evidence = PROJECT_ROOT / reference_config["evidence_path"]
        reference_rows = completed(
            PROJECT_ROOT,
            reference_config,
            reference_summary,
            reference_config_data,
            reference_include_data,
            reference_entry,
        )
        reference_tree = banked_result_tree(reference_evidence)
        audit_child_observed_mutants(reference_evidence, reference_rows)

        copy_fixture("a3", resumed_dir)
        config_path, config = synthetic_config(resumed_dir, "a3", artifact=False)
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        target = resumed_dir / "target.txt"
        process = start_run(config_path)
        wait_for(
            lambda: b"A3_MUTANT_BOUNDARY" in target.read_bytes(),
            "A3 interruption boundary was not reached",
        )
        evidence = PROJECT_ROOT / config["evidence_path"]
        wait_for(
            lambda: (evidence / "results.tsv").is_file()
            and partial_keys(evidence) == ["first"],
            "A3 first row was not banked at the named boundary",
        )
        os.kill(process.pid, signal.SIGTERM)
        interrupted = finish_process(process)
        kill_leftover_group(process.pid)
        oracle.require(interrupted.returncode == 75, f"A3 signal exit {interrupted.returncode}")
        oracle.require(
            target.read_bytes() == entry_snapshot[("target", "stage")],
            "A3 signal did not restore target",
        )
        oracle.audit_restore(
            PROJECT_ROOT,
            evidence,
            config,
            entry_snapshot,
            required_operation="signal",
        )
        banked_prefix = partial_keys(evidence)

        observed: set[str] = set()
        stop = threading.Event()

        def watch() -> None:
            mapping = {
                b"A3_MUTANT_FIRST": "first",
                b"A3_MUTANT_BOUNDARY": "boundary",
                b"A3_MUTANT_LAST": "last",
            }
            while not stop.is_set():
                data = target.read_bytes()
                for token, key in mapping.items():
                    if token in data:
                        observed.add(key)
                time.sleep(0.005)

        watcher = threading.Thread(target=watch, daemon=True)
        watcher.start()
        resumed_process = start_run(config_path)
        resumed = finish_process(resumed_process)
        stop.set()
        watcher.join(timeout=2)
        resumed_summary = parse_success(resumed, "resumed run")
        resumed_rows = completed(
            PROJECT_ROOT,
            config,
            resumed_summary,
            config_data,
            include_data,
            entry_snapshot,
        )
        audit_child_observed_mutants(evidence, resumed_rows)
        resumed_tree = banked_result_tree(evidence)
        facts = {
            "banked_prefix": banked_prefix,
            "resume_observed": observed,
            "resumed_results": resumed_tree,
            "uninterrupted_results": reference_tree,
            "complete": len(resumed_rows) == len(reference_rows),
        }
        oracle.grade_a3(facts)


def a4_trial(mode: str) -> dict[str, Any]:
    with runtime(f"a4-{mode}") as run_dir:
        copy_fixture("a4", run_dir)
        config_path, config = synthetic_config(
            run_dir, "a4", artifact=True, artifact_entry=True
        )
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        target = run_dir / "target.txt"
        artifact = run_dir / "artifact.txt"
        evidence = PROJECT_ROOT / config["evidence_path"]
        if mode == "clean":
            process = start_run(config_path)
            wait_for(
                lambda: b"A4_MUTANT" in target.read_bytes() and not artifact.exists(),
                "A4 clean mutation boundary was not reached",
            )
            finished = finish_process(process)
            summary = parse_success(finished, "run")
            rows = completed(
                PROJECT_ROOT,
                config,
                summary,
                config_data,
                include_data,
                entry_snapshot,
            )
            audit_child_observed_mutants(evidence, rows)
        else:
            process = start_run(config_path)
            wait_for(
                lambda: b"A4_MUTANT" in target.read_bytes() and not artifact.exists(),
                f"A4 {mode} mutation boundary was not reached",
            )
            if mode == "signal":
                os.kill(process.pid, signal.SIGTERM)
                stopped = finish_process(process)
                kill_leftover_group(process.pid)
                oracle.require(stopped.returncode == 75, f"A4 signal exit {stopped.returncode}")
                oracle.audit_restore(
                    PROJECT_ROOT,
                    evidence,
                    config,
                    entry_snapshot,
                    required_operation="signal",
                )
            else:
                os.killpg(process.pid, signal.SIGKILL)
                stopped = finish_process(process)
                oracle.require(stopped.returncode < 0, "A4 uncatchable kill did not kill the run")
                repair_summary, _ = invoke("repair", config_path)
                oracle.require(
                    set(repair_summary) == {"schema", "evidence", "restore_event", "restored"},
                    "repair summary member set",
                )
                oracle.require(
                    repair_summary["schema"] == "ablation-repair-summary-v1"
                    and repair_summary["evidence"] == rel(evidence)
                    and repair_summary["restored"] is True,
                    "repair summary",
                )
                proof_rows = oracle.audit_restore(
                    PROJECT_ROOT,
                    evidence,
                    config,
                    entry_snapshot,
                    required_operation="repair",
                )
                oracle.require(
                    repair_summary["restore_event"] == int(proof_rows[-1]["event"]),
                    "repair summary event",
                )
            resumed_summary, _ = invoke("run", config_path)
            rows = completed(
                PROJECT_ROOT,
                config,
                resumed_summary,
                config_data,
                include_data,
                entry_snapshot,
            )
            audit_child_observed_mutants(evidence, rows)
        return {
            "target_restored": target.read_bytes() == entry_snapshot[("target", "stage")],
            "artifact_restored": artifact.read_bytes()
            == entry_snapshot[("artifact", "probe")],
        }


def test_a4() -> None:
    oracle.grade_a4(
        {
            "clean": a4_trial("clean"),
            "signal": a4_trial("signal"),
            "kill_repair": a4_trial("kill"),
        }
    )


def test_a5() -> None:
    with runtime("a5") as run_dir:
        copy_fixture("a5", run_dir)
        config_path, config = synthetic_config(run_dir, "a5", artifact=True)
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        summary, _ = invoke("run", config_path)
        evidence = PROJECT_ROOT / config["evidence_path"]
        rows = completed(
            PROJECT_ROOT,
            config,
            summary,
            config_data,
            include_data,
            entry_snapshot,
        )
        audit_child_observed_mutants(evidence, rows)
        oracle.grade_a5(include_data, evidence)


def a6_config(
    project_root: pathlib.Path,
    run_dir: pathlib.Path,
) -> tuple[pathlib.Path, dict[str, Any]]:
    target = project_root / STAGE_TARGET_REL
    config: dict[str, Any] = {
        "artifacts": [],
        "checkers": [
            {
                "argv": [
                    "python3",
                    CHILD_REL.as_posix(),
                    "checker",
                    "a3",
                    rel(target, project_root),
                    "-",
                    "fixture",
                ],
                "name": "parse_checker",
                "report_format": "exit-only",
                "timeout_seconds": 10,
            }
        ],
        "contract": "ablation-driver-v1",
        "environment": {
            "LC_ALL": "C",
            "PATH": "/usr/local/bin:/usr/bin:/bin",
            "PYTHONDONTWRITEBYTECODE": "1",
            "TZ": "UTC",
        },
        "evidence_path": rel(run_dir / "evidence", project_root),
        "include_list": (
            "research/pde_ledger_v2/notes/stage023_step_h_evidence/include_list.tsv"
        ),
        "producer": {
            "argv": [
                "python3",
                CHILD_REL.as_posix(),
                "producer",
                "a3",
                rel(target, project_root),
                "-",
                "fixture",
            ],
            "report_format": "exit-only",
            "timeout_seconds": 10,
        },
        "recovery_path": rel(run_dir / "recovery", project_root),
        "reset_paths": [],
        "stable_inputs": [CHILD_REL.as_posix()],
        "targets": [{"name": "stage023", "path": STAGE_TARGET_REL.as_posix()}],
        "workdir": ".",
    }
    config_path = run_dir / "config.json"
    write_jcs(config_path, config)
    return config_path, config


def test_a6() -> None:
    with stage023_mirror() as mirror:
        run_dir = mirror / "fixture-runtime/a6"
        run_dir.mkdir(parents=True)
        config_path, config = a6_config(mirror, run_dir)
        validation, _ = invoke("validate", config_path, project_root=mirror)
        oracle.require(not (run_dir / "evidence").exists(), "validate created evidence")
        oracle.require(not (run_dir / "recovery").exists(), "validate created recovery state")
        named_header, named_body = oracle.read_tsv_bytes(
            (SUITE / "fixtures/a6/named_rows.tsv").read_bytes(), "A6 named rows"
        )
        named = [
            (row[named_header.index("axis")], row[named_header.index("key")])
            for row in named_body
        ]
        include_data = committed_bytes(LEGACY_DIR / "include_list.tsv")
        mirror_include = mirror / config["include_list"]
        oracle.require(
            mirror_include.read_bytes() == include_data,
            "A6 mirror include-list differs from the committed oracle",
        )
        oracle.grade_a6(
            include_data,
            validation,
            named,
        )
        oracle.require(
            all(row["target"] == config["targets"][0]["name"] for row in validation["rows"]),
            "A6 resolved target",
        )
        metamorphic_a6_validation(mirror, run_dir, config, include_data)
        invalid_a6_validation(mirror, run_dir, config, include_data)
    audit_a6_result_contracts()


def metamorphic_a6_validation(
    project_root: pathlib.Path,
    run_dir: pathlib.Path,
    base_config: dict[str, Any],
    include_data: bytes,
) -> None:
    source_header, source_body = oracle.read_tsv_bytes(include_data, "A6 source list")
    source_rows = [dict(zip(source_header, cells, strict=True)) for cells in source_body]
    nonce = os.urandom(18).hex()
    header = ["nonce", "new_text", "target", "key", "line", "old_text", "axis", "record"]
    body: list[list[str]] = []
    for position, source in enumerate(source_rows, start=1):
        values = {
            **source,
            "line": str(position),
            "nonce": f"{nonce}.{position}",
            "target": "stage023",
        }
        body.append([values[name] for name in header])
    derivative_data = (
        "\n".join("\t".join(cells) for cells in [header, *body]) + "\n"
    ).encode("utf-8")
    derivative_path = run_dir / f"metamorphic-{nonce}.tsv"
    derivative_path.write_bytes(derivative_data)
    config = dict(base_config)
    config["include_list"] = rel(derivative_path, project_root)
    config_path = run_dir / f"metamorphic-{nonce}.json"
    write_jcs(config_path, config)
    validation, _ = invoke("validate", config_path, project_root=project_root)

    expected_rows = [
        {
            "cells": cells,
            "row_id": (
                f"stage023:{cells[header.index('axis')]}:{cells[header.index('key')]}"
            ),
            "target": "stage023",
        }
        for cells in body
    ]
    expected_counts = dict(
        Counter(cells[header.index("axis")] for cells in body)
    )
    oracle.require(
        validation
        == {
            "axis_counts": expected_counts,
            "columns": header,
            "list_sha256": hashlib.sha256(derivative_data).hexdigest(),
            "row_count": len(body),
            "rows": expected_rows,
            "schema": "ablation-list-validation-v1",
        },
        "A6 runtime-generated lossless validation",
    )


def invalid_a6_validation(
    project_root: pathlib.Path,
    run_dir: pathlib.Path,
    base_config: dict[str, Any],
    include_data: bytes,
) -> None:
    header, body = oracle.read_tsv_bytes(include_data, "A6 invalid-list source")
    invalid_body = [list(cells) for cells in body]
    invalid_body[0][header.index("new_text")] = invalid_body[0][header.index("old_text")]
    invalid_data = (
        "\n".join("\t".join(cells) for cells in [header, *invalid_body]) + "\n"
    ).encode("utf-8")
    nonce = os.urandom(18).hex()
    invalid_path = run_dir / f"invalid-{nonce}.tsv"
    invalid_path.write_bytes(invalid_data)
    config = dict(base_config)
    config["include_list"] = rel(invalid_path, project_root)
    config_path = run_dir / f"invalid-{nonce}.json"
    write_jcs(config_path, config)
    process = subprocess.run(
        [DRIVER_REL.as_posix(), "validate", "--config", rel(config_path, project_root)],
        cwd=project_root,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=driver_environment(),
        timeout=30,
        check=False,
    )
    oracle.require(process.returncode == 65, f"A6 invalid-list exit {process.returncode}")
    oracle.require(process.stderr != b"", "A6 invalid-list failure has no diagnostic")
    oracle.require(not (run_dir / "evidence").exists(), "invalid validate created evidence")
    oracle.require(not (run_dir / "recovery").exists(), "invalid validate created recovery state")


def audit_a6_result_contracts() -> None:
    with runtime("a6-result-audit") as run_dir:
        copy_fixture("a9", run_dir)
        config_path, config = synthetic_config(run_dir, "a9", artifact=False)
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        summary, _ = invoke("run", config_path)
        rows = completed(
            PROJECT_ROOT,
            config,
            summary,
            config_data,
            include_data,
            entry_snapshot,
        )
        audit_child_observed_mutants(
            PROJECT_ROOT / config["evidence_path"],
            rows,
        )


@contextlib.contextmanager
def stage023_mirror() -> Iterator[pathlib.Path]:
    mirror = pathlib.Path(tempfile.mkdtemp(prefix="fixtures-v4-a7-mirror-")).resolve()
    try:
        for relative in STAGE_MIRROR_INPUTS:
            source = PROJECT_ROOT / relative
            destination = mirror / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        driver_destination = mirror / DRIVER_REL
        driver_destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(DRIVER, driver_destination)
        yield mirror
    finally:
        shutil.rmtree(mirror, ignore_errors=True)


def a7_config(
    project_root: pathlib.Path,
    run_dir: pathlib.Path,
) -> tuple[pathlib.Path, dict[str, Any]]:
    config: dict[str, Any] = {
        "artifacts": [
            {
                "name": "stage023_sidecar",
                "observation_format": "ledger-dimension-v1",
                "path": STAGE_ARTIFACT_REL.as_posix(),
            }
        ],
        "checkers": [
            {
                "argv": ["python3", "compare_dimension_artifacts.py", "023"],
                "name": "stage023_comparator",
                "report_format": "ledger-result-v1",
                "timeout_seconds": 600,
            }
        ],
        "contract": "ablation-driver-v1",
        "environment": {
            "LC_ALL": "C",
            "PATH": "/usr/local/bin:/usr/bin:/bin",
            "PYTHONDONTWRITEBYTECODE": "1",
            "TZ": "UTC",
        },
        "evidence_path": rel(run_dir / "evidence", project_root),
        "include_list": (
            "research/pde_ledger_v2/notes/stage023_step_h_evidence/include_list.tsv"
        ),
        "producer": {
            "argv": [
                "python3",
                "ledger_stage023_nullspace_underdetermination_sympy_audit.py",
            ],
            "report_format": "ledger-tally-v1",
            "timeout_seconds": 600,
        },
        "recovery_path": rel(run_dir / "recovery", project_root),
        "reset_paths": ["research/pde_ledger_v2/scripts/__pycache__"],
        "stable_inputs": [
            "research/pde_ledger_v2/scripts/check_ledger_dimensions_pin.py",
            "research/pde_ledger_v2/scripts/compare_dimension_artifacts.py",
            "research/pde_ledger_v2/scripts/ledger_dimensions.accepted.sha256",
            "research/pde_ledger_v2/scripts/ledger_dimensions.py",
            "research/pde_ledger_v2/mathematica/out/"
            "ledger_stage023_nullspace_underdetermination_mathematica_audit.out",
        ],
        "targets": [{"name": "stage023", "path": STAGE_TARGET_REL.as_posix()}],
        "workdir": "research/pde_ledger_v2/scripts",
    }
    config_path = run_dir / "config.json"
    write_jcs(config_path, config)
    return config_path, config


def banked_row_count(
    project_root: pathlib.Path,
    evidence: pathlib.Path,
    config: dict[str, Any],
    include_data: bytes,
    entry_snapshot: oracle.EntrySnapshot,
    verified_count: int,
) -> int:
    path = evidence / "results.tsv"
    if not path.is_file():
        return 0
    complete_lines = path.read_bytes().split(b"\n")[:-1]
    if not complete_lines:
        return 0
    visible_count = len(complete_lines) - 1
    if visible_count <= verified_count:
        return visible_count
    result_data = b"\n".join(complete_lines) + b"\n"
    rows = oracle.parse_and_audit_results(
        project_root,
        evidence,
        config,
        include_data,
        entry_snapshot,
        result_data=result_data,
        require_complete=False,
    )
    return len(rows)


def stop_process_group(
    process: subprocess.Popen[bytes],
    *,
    item: str,
    grace_seconds: float = 15,
) -> subprocess.CompletedProcess[bytes]:
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        stdout, stderr = process.communicate(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        kill_leftover_group(process.pid)
        try:
            stdout, stderr = process.communicate(timeout=KILL_FALLBACK_WAIT_SECONDS)
        except subprocess.TimeoutExpired as kill_failure:
            close_process_pipes(process)
            raise oracle.FixtureFailure(
                f"{item} process group did not close after its bounded SIGKILL fallback"
            ) from kill_failure
    finally:
        kill_leftover_group(process.pid)
    return subprocess.CompletedProcess(process.args, process.returncode, stdout, stderr)


def run_a7_uninterrupted(
    project_root: pathlib.Path,
    config_path: pathlib.Path,
    config: dict[str, Any],
    include_data: bytes,
    entry_snapshot: oracle.EntrySnapshot,
    *,
    inactivity_seconds: float | None = None,
    progress: Any = None,
) -> subprocess.CompletedProcess[bytes]:
    row_allowance = config["producer"]["timeout_seconds"] + sum(
        checker["timeout_seconds"] for checker in config["checkers"]
    )
    inactivity = (
        float(row_allowance + A7_PROGRESS_GRACE_SECONDS)
        if inactivity_seconds is None
        else inactivity_seconds
    )
    evidence = project_root / config["evidence_path"]
    process = start_run(config_path, project_root)
    last_count = 0
    last_progress = time.monotonic()
    run_failure: BaseException | None = None
    shutdown_failure: BaseException | None = None
    finished: subprocess.CompletedProcess[bytes] | None = None
    try:
        while process.poll() is None:
            count = banked_row_count(
                project_root,
                evidence,
                config,
                include_data,
                entry_snapshot,
                last_count,
            )
            oracle.require(count >= last_count, "A7 banked-row count regressed")
            if count > last_count:
                for position in range(last_count + 1, count + 1):
                    if progress is None:
                        print(f"A7 progress: banked row {position}", flush=True)
                    else:
                        progress(position)
                last_count = count
                last_progress = time.monotonic()
            if time.monotonic() - last_progress > inactivity:
                raise oracle.FixtureFailure(
                    f"A7 made no banked-row progress for {int(inactivity)} seconds"
                )
            time.sleep(A7_PROGRESS_POLL_SECONDS)
    except BaseException as exc:
        run_failure = exc
    finally:
        try:
            finished = stop_process_group(process, item="A7")
        except BaseException as exc:
            shutdown_failure = exc
    if run_failure is not None:
        if shutdown_failure is not None:
            print(
                "A7 shutdown FAIL "
                f"{type(shutdown_failure).__name__}: {shutdown_failure}"
            )
        raise run_failure
    if shutdown_failure is not None:
        raise shutdown_failure
    final_count = banked_row_count(
        project_root,
        evidence,
        config,
        include_data,
        entry_snapshot,
        last_count,
    )
    if final_count > last_count:
        for position in range(last_count + 1, final_count + 1):
            if progress is None:
                print(f"A7 progress: banked row {position}", flush=True)
            else:
                progress(position)
    return finished


def path_state(path: pathlib.Path) -> tuple[bool, bytes | None]:
    return path.is_file(), path.read_bytes() if path.is_file() else None


def test_a7() -> None:
    caller_paths = [PROJECT_ROOT / STAGE_TARGET_REL, PROJECT_ROOT / STAGE_ARTIFACT_REL]
    caller_entry = {path: path_state(path) for path in caller_paths}
    try:
        with stage023_mirror() as mirror:
            run_dir = mirror / "fixture-runtime/a7"
            run_dir.mkdir(parents=True)
            config_path, config = a7_config(mirror, run_dir)
            config_data, include_data, entry_snapshot = capture_invocation(
                mirror,
                config_path,
                config,
            )
            finished = run_a7_uninterrupted(
                mirror,
                config_path,
                config,
                include_data,
                entry_snapshot,
            )
            summary = parse_success(finished, "A7 run")
            rows = completed(
                mirror,
                config,
                summary,
                config_data,
                include_data,
                entry_snapshot,
            )
            candidate = oracle.new_projection(rows)
            oracle.grade_a7(
                committed_bytes(LEGACY_DIR / "include_list.tsv"),
                committed_bytes(LEGACY_DIR / "results.tsv"),
                candidate,
            )
    finally:
        for path, state in caller_entry.items():
            oracle.require(
                path_state(path) == state,
                f"A7 changed caller checkout path {rel(path)}",
            )


@contextlib.contextmanager
def run_a8_fixture(
    label: str,
) -> Iterator[
    tuple[pathlib.Path, list[dict[str, Any]], oracle.EntrySnapshot]
]:
    with runtime(label) as run_dir:
        copy_fixture("a8", run_dir)
        config_path, config = synthetic_config(
            run_dir,
            "a8",
            artifact=True,
            artifact_entry=True,
            artifact_format="ledger-dimension-v1",
            checker_format="ledger-result-v1",
        )
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        summary, _ = invoke("run", config_path)
        rows = completed(
            PROJECT_ROOT,
            config,
            summary,
            config_data,
            include_data,
            entry_snapshot,
        )
        audit_child_observed_mutants(
            PROJECT_ROOT / config["evidence_path"],
            rows,
        )
        yield PROJECT_ROOT / config["evidence_path"], rows, entry_snapshot


def test_a8() -> None:
    with run_a8_fixture("a8") as (evidence, _initial_rows, entry_snapshot):
        config_path = evidence.parent / "config.json"
        include_path = evidence.parent / "include_list.tsv"
        config_path.unlink()
        include_path.unlink()
        rows = oracle.audit_completed_evidence(
            PROJECT_ROOT,
            evidence,
            CONTRACT,
            entry_snapshot,
            invocation_config=None,
            invocation_include=None,
            required_restore_operation="run",
        )
        config = oracle.parse_standalone_jcs((evidence / "config.json").read_bytes(), "config.json")
        restore_rows = oracle.audit_restore(
            PROJECT_ROOT,
            evidence,
            config,
            entry_snapshot,
            required_operation="run",
        )
        oracle.grade_a8_clean_restore_claims(evidence, rows, restore_rows, config)


def test_a9() -> None:
    with runtime("a9") as run_dir:
        copy_fixture("a9", run_dir)
        config_path, config = synthetic_config(run_dir, "a9", artifact=False)
        config_data, include_data, entry_snapshot = capture_invocation(
            PROJECT_ROOT,
            config_path,
            config,
        )
        summary, _ = invoke("run", config_path)
        rows = completed(
            PROJECT_ROOT,
            config,
            summary,
            config_data,
            include_data,
            entry_snapshot,
        )
        audit_child_observed_mutants(
            PROJECT_ROOT / config["evidence_path"],
            rows,
        )
        oracle.grade_a9(rows)


TESTS = {
    "A1": test_a1,
    "A2": test_a2,
    "A3": test_a3,
    "A4": test_a4,
    "A5": test_a5,
    "A6": test_a6,
    "A7": test_a7,
    "A8": test_a8,
    "A9": test_a9,
}


def verify_freeze_gate() -> bool:
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    process = subprocess.run(
        [sys.executable, str(SUITE / "verify_freeze.py")],
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=environment,
        timeout=60,
        check=False,
    )
    if process.stdout:
        print(process.stdout.decode("utf-8", errors="replace"), end="")
    if process.stderr:
        print(process.stderr.decode("utf-8", errors="replace"), end="", file=sys.stderr)
    if process.returncode != 0:
        print("conformance refused: fixture freeze verification failed")
        return False
    return True


TreeRecord = tuple[str, int, str | None, tuple[int, int, int, int, int] | None]
TreeSnapshot = dict[str, TreeRecord]


def file_sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def caller_tree_snapshot(
    root: pathlib.Path,
    previous: TreeSnapshot | None = None,
) -> TreeSnapshot:
    root = root.resolve()
    snapshot: TreeSnapshot = {}
    pending = [root]
    while pending:
        directory = pending.pop()
        entries = sorted(os.scandir(directory), key=lambda entry: os.fsencode(entry.name))
        for entry in entries:
            path = pathlib.Path(entry.path)
            relative = path.relative_to(root).as_posix()
            information = entry.stat(follow_symlinks=False)
            mode = stat.S_IMODE(information.st_mode)
            signature = (
                information.st_dev,
                information.st_ino,
                information.st_size,
                information.st_mtime_ns,
                information.st_ctime_ns,
            )
            if stat.S_ISDIR(information.st_mode):
                snapshot[relative] = ("directory", mode, None, None)
                pending.append(path)
            elif stat.S_ISREG(information.st_mode):
                prior = previous.get(relative) if previous is not None else None
                content = (
                    prior[2]
                    if prior is not None
                    and prior[0] == "file"
                    and prior[3] == signature
                    else file_sha256(path)
                )
                snapshot[relative] = ("file", mode, content, signature)
            elif stat.S_ISLNK(information.st_mode):
                snapshot[relative] = ("symlink", mode, os.readlink(path), None)
            elif stat.S_ISFIFO(information.st_mode):
                snapshot[relative] = ("fifo", mode, None, None)
            elif stat.S_ISSOCK(information.st_mode):
                snapshot[relative] = ("socket", mode, None, None)
            elif stat.S_ISBLK(information.st_mode):
                snapshot[relative] = ("block", mode, None, None)
            elif stat.S_ISCHR(information.st_mode):
                snapshot[relative] = ("character", mode, None, None)
            else:
                snapshot[relative] = ("other", mode, None, None)
    return snapshot


def require_caller_tree_unchanged(
    expected: TreeSnapshot,
    root: pathlib.Path = PROJECT_ROOT,
) -> None:
    actual = caller_tree_snapshot(root, expected)
    expected_paths = set(expected)
    actual_paths = set(actual)
    added = sorted(actual_paths - expected_paths, key=os.fsencode)
    removed = sorted(expected_paths - actual_paths, key=os.fsencode)
    if added:
        raise oracle.FixtureFailure(f"caller tree gained path {added[0]}")
    if removed:
        raise oracle.FixtureFailure(f"caller tree lost path {removed[0]}")
    for path in sorted(expected_paths, key=os.fsencode):
        if actual[path][:3] != expected[path][:3]:
            raise oracle.FixtureFailure(f"caller tree changed path {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("items", nargs="*")
    args = parser.parse_args()
    unknown = [item for item in args.items if item not in ALL_ITEMS]
    if unknown:
        parser.error(f"unknown item: {unknown[0]}")
    selected = args.items or ALL_ITEMS
    if not verify_freeze_gate():
        return 1
    global oracle
    import fixture_oracle as oracle

    if not DRIVER.is_file():
        print(f"driver missing at contract path: {DRIVER_REL}")
        return 2
    caller_entry = caller_tree_snapshot(PROJECT_ROOT)
    failures = 0
    for item in selected:
        failure: Exception | None = None
        try:
            TESTS[item]()
        except Exception as exc:
            failure = exc
        try:
            require_caller_tree_unchanged(caller_entry)
        except Exception as exc:
            if failure is None:
                failure = exc
            else:
                print(f"{item} cleanup FAIL {type(exc).__name__}: {exc}")
        if failure is not None:
            failures += 1
            print(f"{item} FAIL {type(failure).__name__}: {failure}")
        else:
            print(f"{item} PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
