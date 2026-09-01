#!/usr/bin/env python3
"""Run the Codex memory writer behind a filesystem hard boundary.

The launcher mounts one frozen task packet read-only at /packet, exposes one
persistent output directory at /output, hides the live workspace with a tmpfs,
clears the inherited environment, and imports exactly /output/page.md after a
successful run. Grok review is deliberately separate from this publication
path; Codex is the only authoring runtime.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence

import yaml

try:  # Direct script execution resolves the sibling as top-level ``memory``.
    from . import memory as mem
    from .review_contract import (
        GROK_REVIEW_MODEL, REVIEW_CONTRACT_SHA256, REVIEW_PROMPT_SHA256, REVIEW_SCHEMA_SHA256,
    )
except ImportError:  # pragma: no cover - exercised by CLI smoke usage
    import memory as mem  # type: ignore[no-redef]
    from review_contract import (  # type: ignore[no-redef]
        GROK_REVIEW_MODEL, REVIEW_CONTRACT_SHA256, REVIEW_PROMPT_SHA256, REVIEW_SCHEMA_SHA256,
    )


class IsolationError(mem.MemoryErrorBase):
    pass


# Hard wall-clock limit for one Codex authoring invocation.  A timed-out run is
# never imported into the transaction's staged memory.
CODEX_WRITER_TIMEOUT_SECONDS = 20 * 60


WRITER_PROMPT = (
    "Read /packet/task.json, every frozen prompt named there, /packet/schema.md, and only the sealed "
    "inputs made visible below /packet. Produce the complete requested Markdown page with no wrapper "
    "or commentary. Return that page as your final response; the trusted launcher writes the final "
    "response exactly to /output/page.md."
)

REVISION_PROMPT = (
    "Read /packet/revision.json, /packet/revision_candidate.md, /packet/revision_review.md, and "
    "/packet/revision_review_attestation.json, then read /packet/task.json, every "
    "frozen prompt named there, /packet/schema.md, and the current sealed source inputs below /packet. "
    "Use the prior candidate as the starting point: preserve its accurate, supported content and stable "
    "IDs, while resolving every blocking/major failed-review finding and addressing the current task's "
    "editorial_focus and frozen policy. Current sealed sources and prerequisite pages remain authoritative "
    "and may be used. Produce the complete revised Markdown page with no wrapper or "
    "commentary. Return that page as your final response; the trusted launcher writes the final response "
    "exactly to /output/page.md."
)

LINT_REVISION_PROMPT = (
    "Read /packet/revision.json, /packet/revision_candidate.md, /packet/revision_review.md, "
    "/packet/revision_review_attestation.json, /packet/revision_lint_record.json, and "
    "/packet/revision_lint_report.md, then read /packet/task.json, every frozen prompt named there, "
    "and /packet/schema.md. Use the prior candidate as the starting point and change only what is "
    "needed to resolve every exact machine-lint error. Preserve accurate supported content, stable IDs, "
    "and all content unrelated to those lint errors. Current sealed sources remain authoritative. Produce "
    "the complete corrected Markdown page with no wrapper or commentary. Return that page as your final "
    "response; the trusted launcher writes the final response exactly to /output/page.md."
)


@dataclasses.dataclass(frozen=True)
class RuntimeProfile:
    name: str
    version: str
    executable_sha256: str
    command: tuple[str, ...]
    ro_binds: tuple[tuple[str, str], ...]
    environment: tuple[tuple[str, str], ...]
    capture_stdout: bool = False
    stdin_data: bytes | None = None


def _probe_version(command: Sequence[str]) -> str:
    proc = subprocess.run(
        list(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        env={"HOME": "/nonexistent", "PATH": "/usr/bin:/bin", "LANG": "C.UTF-8"},
    )
    if proc.returncode:
        raise IsolationError(f"runtime version probe failed for {command[0]}")
    return proc.stdout.decode("utf-8", "replace").strip().splitlines()[0]


def _timeout_output(value: bytes | str | None) -> bytes:
    """Normalize partial output attached to ``subprocess.TimeoutExpired``."""
    if value is None:
        return b""
    if isinstance(value, bytes):
        return value
    return value.encode("utf-8", "replace")


def _preserve_timeout_diagnostics(output: Path, exc: subprocess.TimeoutExpired) -> None:
    """Best-effort preservation of partial writer diagnostics inside the task."""
    for name, data in (
        ("failed-stdout.bin", _timeout_output(exc.stdout)),
        ("failed-stderr.bin", _timeout_output(exc.stderr)),
    ):
        try:
            mem.atomic_write(output / name, data, 0o600)
        except OSError:
            # The timeout remains the primary failure if the diagnostic volume
            # itself became unavailable or unwritable.
            pass


def _required_file(path: Path, label: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_file():
        raise IsolationError(f"{label} is unavailable: {path}")
    return resolved


def runtime_profile(name: str) -> RuntimeProfile:
    home = Path.home()
    network_files = tuple(
        (str(path), str(path)) for path in (
            Path("/etc/resolv.conf"), Path("/etc/hosts"), Path("/etc/nsswitch.conf"),
            Path("/etc/ssl/certs/ca-certificates.crt"),
        ) if path.is_file()
    )
    if name == "codex":
        wrapper = _required_file(Path(shutil.which("codex") or ""), "codex CLI")
        openai_root = wrapper.parents[1].parent
        candidates = sorted(openai_root.glob("codex-linux-*/vendor/*/bin/codex"))
        if len(candidates) != 1:
            raise IsolationError(f"expected one installed native Codex runtime, found {len(candidates)}")
        executable = _required_file(candidates[0], "native Codex runtime")
        code_mode_host = _required_file(executable.with_name("codex-code-mode-host"), "Codex code-mode host")
        auth = _required_file(home / ".codex/auth.json", "Codex credential file")
        version = _probe_version([str(executable), "--version"])
        tool_names = ("cat", "find", "sed", "head", "tail", "wc", "grep", "sort", "cut", "tr", "dirname", "basename", "ls")
        tool_binds = tuple(
            (str(_required_file(Path(shutil.which(tool) or ""), f"Codex runtime tool {tool}")), f"/usr/bin/{tool}")
            for tool in tool_names
        )
        bash = _required_file(Path(shutil.which("bash") or ""), "Codex runtime shell")
        rg = _required_file(Path(shutil.which("rg") or ""), "Codex runtime search tool")
        runtime_libraries = tuple(
            (path, path) for path in (
                "/lib64/ld-linux-x86-64.so.2", "/lib/x86_64-linux-gnu/libc.so.6",
                "/lib/x86_64-linux-gnu/libm.so.6", "/lib/x86_64-linux-gnu/libpthread.so.0",
                "/lib/x86_64-linux-gnu/libdl.so.2", "/lib/x86_64-linux-gnu/librt.so.1",
                "/lib/x86_64-linux-gnu/libacl.so.1", "/lib/x86_64-linux-gnu/libgcc_s.so.1",
                "/lib/x86_64-linux-gnu/libpcre.so.3", "/lib/x86_64-linux-gnu/libpcre2-8.so.0",
                "/lib/x86_64-linux-gnu/libselinux.so.1", "/lib/x86_64-linux-gnu/libtinfo.so.6",
            ) if Path(path).is_file()
        )
        return RuntimeProfile(
            name="codex", version=version, executable_sha256=mem.sha256_file(executable),
            command=(
                "/runtime/codex", "exec", "--skip-git-repo-check", "--ephemeral",
                "--ignore-user-config", "--ignore-rules", "--sandbox", "danger-full-access",
                "--dangerously-bypass-approvals-and-sandbox", "-C", "/work",
                "--output-last-message", "/output/page.md", WRITER_PROMPT,
            ),
            ro_binds=(
                (str(executable), "/runtime/codex"),
                (str(code_mode_host), "/runtime/codex-code-mode-host"),
                (str(auth), "/runtime-home/.codex/auth.json"),
                (str(bash), "/bin/bash"), (str(bash), "/bin/sh"), (str(rg), "/usr/bin/rg"),
            ) + tool_binds + runtime_libraries + network_files,
            environment=(("HOME", "/runtime-home"), ("CODEX_HOME", "/runtime-home/.codex"),
                         ("PATH", "/runtime:/usr/bin:/bin"), ("LANG", "C.UTF-8"),
                         ("SSL_CERT_FILE", "/etc/ssl/certs/ca-certificates.crt")),
        )
    raise IsolationError(f"unknown writer runtime profile {name!r}; choose codex")


def packet_digest(packet: Path) -> tuple[dict[str, str], str]:
    files = {
        path.relative_to(packet).as_posix(): mem.sha256_file(path)
        for path in sorted(item for item in packet.rglob("*") if item.is_file())
    }
    return files, mem.sha256_bytes(mem.canonical_json(files))


def verify_packet(task_root: Path, transaction_id: str, task_id: str) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    packet = task_root / "packet"
    seal_path = task_root / "packet-seal.json"
    if not packet.is_dir() or not seal_path.is_file():
        raise IsolationError(f"task packet is incomplete for {task_id}")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    files, combined = packet_digest(packet)
    if seal.get("transaction_id") != transaction_id or seal.get("task_id") != task_id:
        raise IsolationError("task packet identity does not match transaction")
    if seal.get("files") != files or seal.get("combined_sha256") != combined:
        raise IsolationError("task packet changed after prepare")
    task = json.loads((packet / "task.json").read_text(encoding="utf-8"))
    if task.get("transaction_id") != transaction_id or task.get("writer_task", {}).get("task_id") != task_id:
        raise IsolationError("task document identity mismatch")
    return packet, seal, task


def bubblewrap_command(repo: Path, packet: Path, output: Path, profile: RuntimeProfile) -> list[str]:
    bwrap = shutil.which("bwrap")
    if bwrap is None:
        raise IsolationError("bubblewrap is required for the memory hard boundary but was not found")
    arguments = [
        bwrap, "--die-with-parent", "--unshare-all", "--share-net", "--new-session", "--clearenv"
    ]
    workspace = repo.as_posix()
    arguments.extend(["--tmpfs", "/tmp"])
    current = Path("/")
    for part in repo.parts[1:]:
        current = current / part
        if current.as_posix() == "/tmp":
            continue
        arguments.extend(["--dir", current.as_posix()])
    arguments.extend([
        "--tmpfs", workspace,
        "--proc", "/proc",
        "--dev", "/dev",
        "--dir", "/runtime",
        "--dir", "/runtime-home",
        "--dir", "/work",
        "--ro-bind", packet.as_posix(), "/packet",
        "--bind", output.as_posix(), "/output",
        "--chdir", "/work",
    ])
    destination_directories: set[str] = set()
    for _, destination in profile.ro_binds:
        parent = Path(destination).parent
        while parent != Path("/"):
            destination_directories.add(parent.as_posix())
            parent = parent.parent
    for directory in sorted(destination_directories, key=lambda value: (value.count("/"), value)):
        if directory not in ("/runtime", "/runtime-home", "/work", "/packet", "/output"):
            arguments.extend(["--dir", directory])
    for source, destination in profile.ro_binds:
        arguments.extend(["--ro-bind", source, destination])
    minimal_environment = dict(profile.environment)
    minimal_environment.update({
        "PWD": "/work", "GIT_CEILING_DIRECTORIES": "/", "GIT_DISCOVERY_ACROSS_FILESYSTEM": "0",
    })
    for key, value in sorted(minimal_environment.items()):
        arguments.extend(["--setenv", key, value])
    arguments.extend(["--", *profile.command])
    return arguments


def _test_profile(command: Sequence[str]) -> RuntimeProfile:
    if not command:
        raise IsolationError("test runtime command is empty")
    roots = tuple((root, root) for root in ("/usr", "/bin", "/lib", "/lib64", "/etc") if Path(root).exists())
    return RuntimeProfile(
        name="_test", version="test-runtime", executable_sha256=mem.sha256_bytes(command[0].encode()),
        command=tuple(command), ro_binds=roots,
        environment=(("HOME", "/runtime-home"), ("PATH", "/usr/bin:/bin"), ("LANG", "C.UTF-8")),
    )


def _source_unit_identity(manifest: dict[str, Any], source_unit_id: str | None) -> dict[str, Any]:
    unit = next((item for item in manifest.get("units", []) if item.get("id") == source_unit_id), None)
    if source_unit_id is None or unit is None:
        raise IsolationError("revision input is only supported for source-unit writer tasks")

    def identities(key: str) -> list[dict[str, Any]]:
        return [
            {field: value for field, value in member.items() if field != "snapshot_path"}
            for member in unit.get(key, [])
        ]

    return {
        "members": identities("members"),
        "prior_members": identities("prior_members"),
    }


def _revision_task_contract(manifest: dict[str, Any], writer: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return stable retry selection and hash-bearing prerequisite identities."""
    source_unit_id = writer.get("source_unit_id")
    if source_unit_id is not None:
        unit = next(item for item in manifest["units"] if item.get("id") == source_unit_id)
        selector = {
            "task_kind": writer["task_kind"],
            "source_unit_id": source_unit_id,
            "source_kind": unit["kind"],
            "output_repository_path": writer["output_repository_path"],
            "member_identities": _source_unit_identity(manifest, source_unit_id),
        }
        prerequisites = {"unit_digest_sha256": unit["unit_digest_sha256"]}
        return selector, prerequisites
    contract = writer["semantic_contract"]
    selector = {
        "task_kind": writer["task_kind"],
        "source_unit_id": None,
        "output_repository_path": writer["output_repository_path"],
        "generated_regions": writer.get("generated_regions", []),
        "depends_on_tasks": writer.get("depends_on_tasks", []),
        "input_unit_ids": sorted(contract.get("input_units", {})),
        "input_page_ids": list(contract.get("input_pages", [])),
        "direct_sources": [
            {
                "source_path": item.get("source_path"),
                "read_mode": item.get("read_mode"),
                "member_identity": item.get("member_identity"),
            }
            for item in writer.get("semantic_inputs", [])
        ],
    }
    prerequisites = {
        "base_page_sha256": writer.get("base_page_sha256"),
        "static_memory_inputs": contract.get("static_memory_inputs", []),
        "dynamic_memory_inputs": contract.get("dynamic_memory_inputs", []),
    }
    return selector, prerequisites


def _verified_failed_review(
    repo: Path,
    prior_transaction: Path,
    prior_manifest: dict[str, Any],
    prior_writer: dict[str, Any],
    prior_packet_seal: dict[str, Any],
    candidate_sha256: str,
    writer_attestation_path: Path,
) -> dict[str, Any]:
    task_id = prior_writer["task_id"]
    review_root = prior_transaction / "reviews" / task_id
    review_attestation_path = review_root / "attestation.json"
    if not review_attestation_path.is_file() or review_attestation_path.is_symlink():
        raise IsolationError("revision requires a regular failed Grok review attestation")
    review = json.loads(review_attestation_path.read_text(encoding="utf-8"))
    report_relative = f"reviews/{task_id}/output/report.md"
    expected = {
        "role": "independent_review",
        "transaction_id": prior_manifest["transaction_id"],
        "task_id": task_id,
        "packet_sha256": prior_packet_seal["combined_sha256"],
        "candidate_path": prior_writer["staged_output_path"],
        "candidate_sha256": candidate_sha256,
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "report_path": report_relative,
        "verdict": "FAIL",
        "runtime_profile": "grok-review",
        "review_model": GROK_REVIEW_MODEL,
        "review_prompt_sha256": REVIEW_PROMPT_SHA256,
        "review_schema_sha256": REVIEW_SCHEMA_SHA256,
        "review_contract_sha256": REVIEW_CONTRACT_SHA256,
        "reviewer_sha256": mem.sha256_file(Path(__file__).resolve().parent / "review_isolated.py"),
    }
    for key, value in expected.items():
        if review.get(key) != value:
            raise IsolationError(f"revision failed Grok attestation mismatch: {key}")
    report = prior_transaction / report_relative
    if not report.is_file() or report.is_symlink() or mem.sha256_file(report) != review.get("report_sha256"):
        raise IsolationError("revision failed Grok report is missing or changed")
    review_packet = review_root / "packet"
    review_files, review_sha256 = packet_digest(review_packet)
    expected_files = dict(prior_packet_seal["files"])
    expected_files.update({"candidate.md": candidate_sha256, "review-prompt.md": review_files.get("review-prompt.md")})
    if (
        review_files != expected_files or review_files != review.get("review_packet_files")
        or review_sha256 != review.get("review_packet_sha256")
    ):
        raise IsolationError("revision failed Grok review packet is stale or changed")
    return {
        "attestation_data": review_attestation_path.read_bytes(),
        "attestation_sha256": mem.sha256_file(review_attestation_path),
        "report_data": report.read_bytes(),
        "report_sha256": mem.sha256_file(report),
        "review_packet_sha256": review_sha256,
    }


def _verified_passed_review(
    prior_transaction: Path,
    prior_manifest: dict[str, Any],
    prior_writer: dict[str, Any],
    prior_packet_seal: dict[str, Any],
    candidate_sha256: str,
    writer_attestation_path: Path,
) -> dict[str, Any]:
    """Verify the exact current-contract Grok PASS bound to a prior candidate."""
    task_id = prior_writer["task_id"]
    review_root = prior_transaction / "reviews" / task_id
    review_attestation_path = review_root / "attestation.json"
    if not review_attestation_path.is_file() or review_attestation_path.is_symlink():
        raise IsolationError("lint revision requires a regular Grok PASS attestation")
    review = json.loads(review_attestation_path.read_text(encoding="utf-8"))
    report_relative = f"reviews/{task_id}/output/report.md"
    expected = {
        "role": "independent_review",
        "transaction_id": prior_manifest["transaction_id"],
        "task_id": task_id,
        "packet_sha256": prior_packet_seal["combined_sha256"],
        "candidate_path": prior_writer["staged_output_path"],
        "candidate_sha256": candidate_sha256,
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "report_path": report_relative,
        "verdict": "PASS",
        "runtime_profile": "grok-review",
        "review_model": GROK_REVIEW_MODEL,
        "review_prompt_sha256": REVIEW_PROMPT_SHA256,
        "review_schema_sha256": REVIEW_SCHEMA_SHA256,
        "review_contract_sha256": REVIEW_CONTRACT_SHA256,
        "reviewer_sha256": mem.sha256_file(Path(__file__).resolve().parent / "review_isolated.py"),
    }
    for key, value in expected.items():
        if review.get(key) != value:
            raise IsolationError(f"lint revision Grok PASS attestation mismatch: {key}")
    for key in ("review_packet_sha256", "report_sha256", "runtime_executable_sha256"):
        _sha256_identity(review.get(key), f"Grok {key}")
    if not isinstance(review.get("runtime_version"), str) or not review["runtime_version"]:
        raise IsolationError("lint revision Grok PASS attestation lacks a runtime version")
    report = prior_transaction / report_relative
    if not report.is_file() or report.is_symlink() or mem.sha256_file(report) != review["report_sha256"]:
        raise IsolationError("lint revision Grok PASS report is missing or changed")
    review_packet = review_root / "packet"
    review_files, review_sha256 = packet_digest(review_packet)
    expected_files = dict(prior_packet_seal["files"])
    expected_files.update({"candidate.md": candidate_sha256, "review-prompt.md": review_files.get("review-prompt.md")})
    if (
        review_files != expected_files or review_files != review.get("review_packet_files")
        or review_sha256 != review.get("review_packet_sha256")
    ):
        raise IsolationError("lint revision Grok PASS review packet is stale or changed")
    return {
        "attestation_data": review_attestation_path.read_bytes(),
        "attestation_sha256": mem.sha256_file(review_attestation_path),
        "report_data": report.read_bytes(),
        "report_sha256": mem.sha256_file(report),
        "review_packet_sha256": review_sha256,
    }


def admit_revision_candidate(
    repo: Path,
    transaction: Path,
    manifest: dict[str, Any],
    writer: dict[str, Any],
    packet: Path,
    packet_seal: dict[str, Any],
    prior_transaction_arg: str,
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    """Append one verified prior Codex candidate to the current sealed packet."""
    task_id = writer["task_id"]
    task_root = transaction / "tasks" / task_id
    if (packet / "revision.json").exists() or (packet / "revision_candidate.md").exists():
        raise IsolationError(f"revision input is already present for {task_id}")
    if list((task_root / "output").iterdir()) or (transaction / "attestations" / f"{task_id}.json").exists():
        raise IsolationError(f"current writer task already has output for {task_id}; prepare a new transaction")

    prior_transaction, prior_manifest = mem.load_transaction(repo, prior_transaction_arg)
    if prior_transaction == transaction:
        raise IsolationError("revision input must come from another transaction")
    prior_writer = next(
        (item for item in prior_manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None,
    )
    if prior_writer is None or not prior_writer.get("required"):
        raise IsolationError(f"revision transaction has no required writer task for {task_id}")
    if (
        prior_writer.get("source_unit_id") != writer.get("source_unit_id")
        or prior_writer.get("output_repository_path") != writer.get("output_repository_path")
        or prior_writer.get("staged_output_path") != writer.get("staged_output_path")
    ):
        raise IsolationError("revision writer task identity or output path does not match")
    current_selector, current_prerequisites = _revision_task_contract(manifest, writer)
    prior_selector, prior_prerequisites = _revision_task_contract(prior_manifest, prior_writer)
    if prior_selector != current_selector:
        raise IsolationError("revision retry selector/member identities do not match")
    target_changed = prior_manifest.get("target_commit") != manifest.get("target_commit")
    if target_changed:
        if writer.get("source_unit_id") is None:
            raise IsolationError("derived revision target commit does not match the current transaction")
        if prior_prerequisites != current_prerequisites:
            raise IsolationError("cross-commit source revision unit digest changed")

    prior_packet, prior_packet_seal, _ = verify_packet(
        prior_transaction / "tasks" / task_id, prior_manifest["transaction_id"], task_id,
    )
    del prior_packet  # The current packet receives only the attested candidate, never prior sources.
    prior_attestation_path = prior_transaction / "attestations" / f"{task_id}.json"
    if not prior_attestation_path.is_file() or prior_attestation_path.is_symlink():
        raise IsolationError(f"revision candidate lacks a Codex isolation attestation: {task_id}")
    prior_attestation = json.loads(prior_attestation_path.read_text(encoding="utf-8"))
    prior_runtime_profile = prior_attestation.get("runtime_profile")
    if prior_runtime_profile not in ("codex", "codex-candidate-reuse"):
        raise IsolationError("revision prior writer runtime_profile is not an eligible Codex candidate")
    required_attestation = {
        "transaction_id": prior_manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": writer.get("source_unit_id"),
        "isolation": (
            "deterministic-candidate-reuse" if prior_runtime_profile == "codex-candidate-reuse"
            else "bubblewrap"
        ),
        "workspace_hidden": repo.as_posix(),
        "packet_path": prior_writer["packet_path"],
        "packet_sha256": prior_packet_seal["combined_sha256"],
        "output_repository_path": prior_writer["output_repository_path"],
        "staged_output_path": prior_writer["staged_output_path"],
        "runtime_profile": prior_runtime_profile,
    }
    for key, value in required_attestation.items():
        if prior_attestation.get(key) != value:
            raise IsolationError(f"revision Codex attestation mismatch: {key}")
    if prior_runtime_profile == "codex":
        if not isinstance(prior_attestation.get("runtime_version"), str) or not prior_attestation["runtime_version"]:
            raise IsolationError("revision Codex attestation lacks a runtime version")
        hash_fields = ("runner_sha256", "runtime_executable_sha256", "output_sha256")
    else:
        if (
            prior_attestation.get("model_invoked") is not False
            or prior_attestation.get("runtime_version") is not None
            or prior_attestation.get("runtime_executable_sha256") is not None
        ):
            raise IsolationError("revision candidate-reuse attestation claims an invalid runtime invocation")
        hash_fields = ("runner_sha256", "output_sha256")
    for key in hash_fields:
        value = prior_attestation.get(key)
        if not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
            raise IsolationError(f"revision Codex attestation lacks a valid {key}")
    candidate = prior_transaction / prior_writer["staged_output_path"]
    if not candidate.is_file() or candidate.is_symlink():
        raise IsolationError("revision staged candidate is missing or is not a regular file")
    candidate_data = candidate.read_bytes()
    if mem.sha256_bytes(candidate_data) != prior_attestation["output_sha256"]:
        raise IsolationError("revision staged candidate changed after its Codex attestation")
    if prior_runtime_profile == "codex-candidate-reuse":
        reuse_errors = mem.verify_candidate_reuse_chain(
            repo, prior_transaction, prior_manifest, prior_writer,
            prior_attestation, prior_attestation["output_sha256"],
        )
        if reuse_errors:
            raise IsolationError(reuse_errors[0])
    failed_review = _verified_failed_review(
        repo, prior_transaction, prior_manifest, prior_writer, prior_packet_seal,
        prior_attestation["output_sha256"], prior_attestation_path,
    )

    selector_sha256 = mem.sha256_bytes(mem.canonical_json(current_selector))
    revision = {
        "revision_version": 1,
        "transaction_id": manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": writer["source_unit_id"],
        "target_commit": manifest["target_commit"],
        "output_repository_path": writer["output_repository_path"],
        "base_packet_sha256": packet_seal["combined_sha256"],
        "candidate_sha256": prior_attestation["output_sha256"],
        "task_kind": writer["task_kind"],
        "retry_selector_sha256": selector_sha256,
        "current_prerequisites_sha256": mem.sha256_bytes(mem.canonical_json(current_prerequisites)),
        "prior_transaction_id": prior_manifest["transaction_id"],
        "prior_target_commit": prior_manifest["target_commit"],
        "prior_output_repository_path": prior_writer["output_repository_path"],
        "prior_staged_output_path": prior_writer["staged_output_path"],
        "prior_runtime_profile": prior_runtime_profile,
        "prior_retry_selector_sha256": mem.sha256_bytes(mem.canonical_json(prior_selector)),
        "prior_prerequisites_sha256": mem.sha256_bytes(mem.canonical_json(prior_prerequisites)),
        "prior_attestation_sha256": mem.sha256_file(prior_attestation_path),
        "prior_packet_sha256": prior_packet_seal["combined_sha256"],
        "prior_review_attestation_sha256": failed_review["attestation_sha256"],
        "prior_review_packet_sha256": failed_review["review_packet_sha256"],
        "prior_review_report_sha256": failed_review["report_sha256"],
    }
    mem.atomic_write(packet / "revision_candidate.md", candidate_data, 0o444)
    mem.atomic_write(packet / "revision_review.md", failed_review["report_data"], 0o444)
    mem.atomic_write(packet / "revision_review_attestation.json", failed_review["attestation_data"], 0o444)
    mem.atomic_write(packet / "revision.json", mem.canonical_json(revision), 0o444)
    packet_files, combined = packet_digest(packet)
    revised_seal = {
        "transaction_id": manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": writer.get("source_unit_id"),
        "files": packet_files,
        "combined_sha256": combined,
    }
    mem.write_json(task_root / "packet-seal.json", revised_seal)
    mem.load_transaction(repo, str(transaction))
    return verify_packet(task_root, manifest["transaction_id"], task_id)


def _render_lint_failure_report(record: dict[str, Any]) -> bytes:
    lines = [
        "Verdict: FAIL", "", "Machine: staged memory lint", "",
        f"Transaction: `{record['transaction_id']}`", "", f"Task: `{record['task_id']}`", "",
        f"Output: `{record['output_repository_path']}`", "",
        f"Candidate SHA-256: `{record['candidate_sha256']}`", "", "## Exact errors", "",
    ]
    lines.extend(f"{index}. {error}" for index, error in enumerate(record["errors"], 1))
    return ("\n".join(lines).rstrip() + "\n").encode("utf-8")


def _verified_lint_revision_input(
    repo: Path, manifest: dict[str, Any], writer: dict[str, Any], prior_transaction_arg: str,
) -> dict[str, Any]:
    """Verify a source candidate, its PASS review, and its trusted machine-lint FAIL."""
    task_id = writer["task_id"]
    if writer.get("source_unit_id") is None:
        raise IsolationError("lint revision is only supported for source-unit tasks")
    prior_transaction, prior_manifest = mem.load_transaction(repo, prior_transaction_arg)
    if prior_manifest.get("transaction_id") == manifest.get("transaction_id"):
        raise IsolationError("lint revision input must come from another transaction")
    prior_writer = next(
        (item for item in prior_manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None,
    )
    if prior_writer is None or not prior_writer.get("required") or prior_writer.get("source_unit_id") is None:
        raise IsolationError(f"lint revision transaction has no source writer task for {task_id}")
    if (
        prior_writer.get("task_kind") != writer.get("task_kind")
        or prior_writer.get("source_unit_id") != writer.get("source_unit_id")
        or prior_writer.get("output_repository_path") != writer.get("output_repository_path")
        or prior_writer.get("staged_output_path") != writer.get("staged_output_path")
        or prior_writer.get("expected_output_name") != writer.get("expected_output_name")
    ):
        raise IsolationError("lint revision writer task identity or output path does not match")
    current_selector, current_prerequisites = _revision_task_contract(manifest, writer)
    prior_selector, prior_prerequisites = _revision_task_contract(prior_manifest, prior_writer)
    if prior_selector != current_selector:
        raise IsolationError("lint revision selector/member identities do not match")
    if prior_prerequisites != current_prerequisites:
        raise IsolationError("lint revision source-unit digest changed")

    _, prior_packet_seal, _ = verify_packet(
        prior_transaction / "tasks" / task_id, prior_manifest["transaction_id"], task_id,
    )
    writer_attestation_path = prior_transaction / "attestations" / f"{task_id}.json"
    if not writer_attestation_path.is_file() or writer_attestation_path.is_symlink():
        raise IsolationError("lint revision lacks a regular prior writer attestation")
    writer_attestation = json.loads(writer_attestation_path.read_text(encoding="utf-8"))
    runtime_profile_name = writer_attestation.get("runtime_profile")
    if runtime_profile_name not in ("codex", "codex-candidate-reuse"):
        raise IsolationError("lint revision prior runtime_profile is not an eligible Codex candidate")
    expected_writer = {
        "transaction_id": prior_manifest["transaction_id"], "task_id": task_id,
        "source_unit_id": prior_writer["source_unit_id"],
        "isolation": "deterministic-candidate-reuse" if runtime_profile_name == "codex-candidate-reuse" else "bubblewrap",
        "workspace_hidden": repo.as_posix(), "packet_path": prior_writer["packet_path"],
        "packet_sha256": prior_packet_seal["combined_sha256"],
        "output_repository_path": prior_writer["output_repository_path"],
        "staged_output_path": prior_writer["staged_output_path"], "runtime_profile": runtime_profile_name,
    }
    for key, value in expected_writer.items():
        if writer_attestation.get(key) != value:
            raise IsolationError(f"lint revision writer attestation mismatch: {key}")
    for key in ("runner_sha256", "output_sha256"):
        _sha256_identity(writer_attestation.get(key), f"writer {key}")
    if runtime_profile_name == "codex":
        _sha256_identity(writer_attestation.get("runtime_executable_sha256"), "writer runtime executable")
        if not isinstance(writer_attestation.get("runtime_version"), str) or not writer_attestation["runtime_version"]:
            raise IsolationError("lint revision Codex writer lacks a runtime version")
    elif (
        writer_attestation.get("model_invoked") is not False
        or writer_attestation.get("runtime_version") is not None
        or writer_attestation.get("runtime_executable_sha256") is not None
    ):
        raise IsolationError("lint revision candidate-reuse writer claims an invalid runtime invocation")
    candidate = prior_transaction / prior_writer["staged_output_path"]
    if not candidate.is_file() or candidate.is_symlink():
        raise IsolationError("lint revision candidate is missing or is not a regular file")
    candidate_data = candidate.read_bytes()
    candidate_sha256 = mem.sha256_bytes(candidate_data)
    if candidate_sha256 != writer_attestation["output_sha256"]:
        raise IsolationError("lint revision candidate changed after its writer attestation")
    if runtime_profile_name == "codex-candidate-reuse":
        reuse_errors = mem.verify_candidate_reuse_chain(
            repo, prior_transaction, prior_manifest, prior_writer, writer_attestation, candidate_sha256,
        )
        if reuse_errors:
            raise IsolationError(reuse_errors[0])

    passed_review = _verified_passed_review(
        prior_transaction, prior_manifest, prior_writer, prior_packet_seal,
        candidate_sha256, writer_attestation_path,
    )
    record_root = prior_transaction / "lint-failures" / task_id
    record_path = record_root / "record.json"
    report_path = record_root / "report.md"
    if (
        not record_path.is_file() or record_path.is_symlink()
        or not report_path.is_file() or report_path.is_symlink()
    ):
        raise IsolationError("lint revision requires a regular recorded lint failure and report")
    record = json.loads(record_path.read_text(encoding="utf-8"))
    errors = record.get("errors")
    warnings = record.get("warnings")
    if not isinstance(errors, list) or not errors or not all(isinstance(item, str) for item in errors):
        raise IsolationError("lint revision record has no exact lint errors")
    if not isinstance(warnings, list) or not all(isinstance(item, str) for item in warnings):
        raise IsolationError("lint revision record has invalid warnings")
    prefix = f"{prior_writer['output_repository_path']}:"
    if any(not error.startswith(prefix) for error in errors):
        raise IsolationError("lint revision record contains ambiguous or unscoped errors")
    report_data = _render_lint_failure_report(record)
    expected_record = {
        "record_version": 1, "role": "machine_lint_failure", "verdict": "FAIL",
        "transaction_id": prior_manifest["transaction_id"],
        "target_commit": prior_manifest["target_commit"],
        "policy_sha256": prior_manifest["policy"]["combined_sha256"],
        "policy_tool_sha256": prior_manifest["policy"]["tool_sha256"],
        "memory_tool_sha256": mem.sha256_file(Path(mem.__file__).resolve()),
        "task_id": task_id, "task_kind": prior_writer["task_kind"],
        "source_unit_id": prior_writer["source_unit_id"],
        "output_repository_path": prior_writer["output_repository_path"],
        "staged_output_path": prior_writer["staged_output_path"],
        "candidate_sha256": candidate_sha256,
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "review_attestation_sha256": passed_review["attestation_sha256"],
        "errors": errors, "errors_sha256": mem.sha256_bytes(mem.canonical_json(errors)),
        "warnings": warnings, "report_path": f"lint-failures/{task_id}/report.md",
        "report_sha256": mem.sha256_bytes(report_data),
    }
    if record != expected_record:
        raise IsolationError("lint revision record does not match trusted transaction/tool identities")
    if report_path.read_bytes() != report_data:
        raise IsolationError("lint revision report is missing, changed, or not deterministic")
    return {
        "prior_transaction": prior_transaction, "prior_manifest": prior_manifest,
        "prior_writer": prior_writer, "prior_packet_seal": prior_packet_seal,
        "writer_attestation_path": writer_attestation_path,
        "runtime_profile": runtime_profile_name, "candidate_data": candidate_data,
        "candidate_sha256": candidate_sha256, "passed_review": passed_review,
        "lint_record_data": record_path.read_bytes(), "lint_record_sha256": mem.sha256_file(record_path),
        "lint_report_data": report_data, "lint_report_sha256": mem.sha256_bytes(report_data),
        "lint_errors_sha256": record["errors_sha256"],
        "current_selector": current_selector, "prior_selector": prior_selector,
        "current_prerequisites": current_prerequisites, "prior_prerequisites": prior_prerequisites,
    }


def admit_lint_revision_candidate(
    repo: Path, transaction: Path, manifest: dict[str, Any], writer: dict[str, Any],
    packet: Path, packet_seal: dict[str, Any], prior_transaction_arg: str,
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    """Seal a reviewed candidate and trusted machine-lint failure into a source retry packet."""
    task_id = writer["task_id"]
    task_root = transaction / "tasks" / task_id
    if (packet / "revision.json").exists() or (packet / "revision_candidate.md").exists():
        raise IsolationError(f"revision input is already present for {task_id}")
    if list((task_root / "output").iterdir()) or (transaction / "attestations" / f"{task_id}.json").exists():
        raise IsolationError(f"current writer task already has output for {task_id}; prepare a new transaction")
    verified = _verified_lint_revision_input(repo, manifest, writer, prior_transaction_arg)
    selector_sha256 = mem.sha256_bytes(mem.canonical_json(verified["current_selector"]))
    revision = {
        "revision_version": 1, "revision_basis": "machine_lint_failure",
        "transaction_id": manifest["transaction_id"], "task_id": task_id,
        "source_unit_id": writer["source_unit_id"], "target_commit": manifest["target_commit"],
        "output_repository_path": writer["output_repository_path"],
        "base_packet_sha256": packet_seal["combined_sha256"],
        "candidate_sha256": verified["candidate_sha256"], "task_kind": writer["task_kind"],
        "retry_selector_sha256": selector_sha256,
        "current_prerequisites_sha256": mem.sha256_bytes(mem.canonical_json(verified["current_prerequisites"])),
        "prior_transaction_id": verified["prior_manifest"]["transaction_id"],
        "prior_target_commit": verified["prior_manifest"]["target_commit"],
        "prior_output_repository_path": verified["prior_writer"]["output_repository_path"],
        "prior_staged_output_path": verified["prior_writer"]["staged_output_path"],
        "prior_runtime_profile": verified["runtime_profile"],
        "prior_retry_selector_sha256": mem.sha256_bytes(mem.canonical_json(verified["prior_selector"])),
        "prior_prerequisites_sha256": mem.sha256_bytes(mem.canonical_json(verified["prior_prerequisites"])),
        "prior_attestation_sha256": mem.sha256_file(verified["writer_attestation_path"]),
        "prior_packet_sha256": verified["prior_packet_seal"]["combined_sha256"],
        "prior_review_attestation_sha256": verified["passed_review"]["attestation_sha256"],
        "prior_review_packet_sha256": verified["passed_review"]["review_packet_sha256"],
        "prior_review_report_sha256": verified["passed_review"]["report_sha256"],
        "prior_lint_record_sha256": verified["lint_record_sha256"],
        "prior_lint_report_sha256": verified["lint_report_sha256"],
        "prior_lint_errors_sha256": verified["lint_errors_sha256"],
    }
    files = {
        "revision_candidate.md": verified["candidate_data"],
        "revision_review.md": verified["passed_review"]["report_data"],
        "revision_review_attestation.json": verified["passed_review"]["attestation_data"],
        "revision_lint_record.json": verified["lint_record_data"],
        "revision_lint_report.md": verified["lint_report_data"],
        "revision.json": mem.canonical_json(revision),
    }
    for name, data in files.items():
        mem.atomic_write(packet / name, data, 0o444)
    packet_files, combined = packet_digest(packet)
    revised_seal = {
        "transaction_id": manifest["transaction_id"], "task_id": task_id,
        "source_unit_id": writer.get("source_unit_id"), "files": packet_files,
        "combined_sha256": combined,
    }
    mem.write_json(task_root / "packet-seal.json", revised_seal)
    mem.load_transaction(repo, str(transaction))
    return verify_packet(task_root, manifest["transaction_id"], task_id)


def _sha256_identity(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
        raise IsolationError(f"reviewed reuse chain lacks a valid {label}")
    return value


def _verified_reviewed_candidate(
    repo: Path,
    manifest: dict[str, Any],
    writer: dict[str, Any],
    prior_transaction_arg: str,
) -> dict[str, Any]:
    """Verify a prior Codex candidate and its Grok PASS review without reading served pages."""
    task_id = writer["task_id"]
    prior_transaction, prior_manifest = mem.load_transaction(repo, prior_transaction_arg)
    if prior_manifest.get("transaction_id") == manifest.get("transaction_id"):
        raise IsolationError("reviewed reuse must come from another transaction")
    prior_writer = next(
        (item for item in prior_manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None,
    )
    if prior_writer is None or not prior_writer.get("required"):
        raise IsolationError(f"reviewed transaction has no required writer task for {task_id}")
    if writer.get("source_unit_id") is None or prior_writer.get("source_unit_id") is None:
        raise IsolationError("reviewed reuse is only supported for source-unit tasks")
    if (
        prior_writer.get("task_kind") != writer.get("task_kind")
        or prior_writer.get("source_unit_id") != writer.get("source_unit_id")
        or prior_writer.get("output_repository_path") != writer.get("output_repository_path")
        or prior_writer.get("staged_output_path") != writer.get("staged_output_path")
        or prior_writer.get("expected_output_name") != writer.get("expected_output_name")
    ):
        raise IsolationError("reviewed reuse task identity or output path does not match")
    current_identity = _source_unit_identity(manifest, writer["source_unit_id"])
    prior_identity = _source_unit_identity(prior_manifest, prior_writer["source_unit_id"])
    if prior_identity != current_identity:
        raise IsolationError("reviewed reuse source-unit member identities do not match")
    current_unit = next(item for item in manifest["units"] if item.get("id") == writer["source_unit_id"])
    prior_unit = next(item for item in prior_manifest["units"] if item.get("id") == writer["source_unit_id"])
    if prior_unit.get("kind") != current_unit.get("kind"):
        raise IsolationError("reviewed reuse source-unit kind does not match")
    if (
        prior_manifest.get("target_commit") != manifest.get("target_commit")
        and prior_unit.get("unit_digest_sha256") != current_unit.get("unit_digest_sha256")
    ):
        raise IsolationError("cross-commit reviewed reuse source-unit digest changed")

    _, prior_packet_seal, _ = verify_packet(
        prior_transaction / "tasks" / task_id, prior_manifest["transaction_id"], task_id,
    )
    writer_attestation_path = prior_transaction / "attestations" / f"{task_id}.json"
    if not writer_attestation_path.is_file() or writer_attestation_path.is_symlink():
        raise IsolationError("reviewed reuse lacks a regular Codex writer attestation")
    writer_attestation = json.loads(writer_attestation_path.read_text(encoding="utf-8"))
    expected_writer = {
        "transaction_id": prior_manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": prior_writer["source_unit_id"],
        "isolation": "bubblewrap",
        "workspace_hidden": repo.as_posix(),
        "packet_path": prior_writer["packet_path"],
        "packet_sha256": prior_packet_seal["combined_sha256"],
        "output_repository_path": prior_writer["output_repository_path"],
        "staged_output_path": prior_writer["staged_output_path"],
        "runtime_profile": "codex",
    }
    for key, value in expected_writer.items():
        if writer_attestation.get(key) != value:
            raise IsolationError(f"reviewed reuse Codex attestation mismatch: {key}")
    for key in ("runner_sha256", "runtime_executable_sha256", "output_sha256"):
        _sha256_identity(writer_attestation.get(key), f"Codex {key}")
    if not isinstance(writer_attestation.get("runtime_version"), str) or not writer_attestation["runtime_version"]:
        raise IsolationError("reviewed reuse Codex attestation lacks a runtime version")
    candidate = prior_transaction / prior_writer["staged_output_path"]
    if not candidate.is_file() or candidate.is_symlink():
        raise IsolationError("reviewed reuse candidate is missing or is not a regular file")
    candidate_data = candidate.read_bytes()
    candidate_sha256 = mem.sha256_bytes(candidate_data)
    if candidate_sha256 != writer_attestation["output_sha256"]:
        raise IsolationError("reviewed reuse candidate changed after its Codex attestation")

    review_root = prior_transaction / "reviews" / task_id
    review_attestation_path = review_root / "attestation.json"
    if not review_attestation_path.is_file() or review_attestation_path.is_symlink():
        raise IsolationError("reviewed reuse lacks a regular Grok review attestation")
    review_attestation = json.loads(review_attestation_path.read_text(encoding="utf-8"))
    expected_report_path = f"reviews/{task_id}/output/report.md"
    expected_review = {
        "role": "independent_review",
        "transaction_id": prior_manifest["transaction_id"],
        "task_id": task_id,
        "packet_sha256": prior_packet_seal["combined_sha256"],
        "candidate_path": prior_writer["staged_output_path"],
        "candidate_sha256": candidate_sha256,
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "report_path": expected_report_path,
        "verdict": "PASS",
        "runtime_profile": "grok-review",
        "review_model": GROK_REVIEW_MODEL,
        "review_prompt_sha256": REVIEW_PROMPT_SHA256,
        "review_schema_sha256": REVIEW_SCHEMA_SHA256,
        "review_contract_sha256": REVIEW_CONTRACT_SHA256,
        "reviewer_sha256": mem.sha256_file(Path(__file__).resolve().parent / "review_isolated.py"),
    }
    for key, value in expected_review.items():
        if review_attestation.get(key) != value:
            raise IsolationError(f"reviewed reuse Grok attestation mismatch: {key}")
    for key in ("review_packet_sha256", "report_sha256", "runtime_executable_sha256"):
        _sha256_identity(review_attestation.get(key), f"Grok {key}")
    if not isinstance(review_attestation.get("runtime_version"), str) or not review_attestation["runtime_version"]:
        raise IsolationError("reviewed reuse Grok attestation lacks a runtime version")
    report = prior_transaction / expected_report_path
    if not report.is_file() or report.is_symlink() or mem.sha256_file(report) != review_attestation["report_sha256"]:
        raise IsolationError("reviewed reuse Grok report is missing or changed")
    review_packet = review_root / "packet"
    if not review_packet.is_dir():
        raise IsolationError("reviewed reuse Grok packet is missing")
    review_packet_files, review_packet_sha256 = packet_digest(review_packet)
    if (
        review_packet_files != review_attestation.get("review_packet_files")
        or review_packet_sha256 != review_attestation["review_packet_sha256"]
        or review_packet_files.get("candidate.md") != candidate_sha256
    ):
        raise IsolationError("reviewed reuse Grok packet or candidate identity changed")
    expected_review_files = dict(prior_packet_seal["files"])
    expected_review_files.update({
        "candidate.md": candidate_sha256,
        "review-prompt.md": review_packet_files.get("review-prompt.md"),
    })
    if review_packet_files != expected_review_files:
        raise IsolationError("reviewed reuse Grok packet is not the writer packet plus reviewed candidate")

    return {
        "prior_transaction": prior_transaction,
        "prior_manifest": prior_manifest,
        "prior_writer": prior_writer,
        "prior_unit": prior_unit,
        "current_unit": current_unit,
        "member_identity_sha256": mem.sha256_bytes(mem.canonical_json(current_identity)),
        "candidate_data": candidate_data,
        "candidate_sha256": candidate_sha256,
        "prior_packet_sha256": prior_packet_seal["combined_sha256"],
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "review_attestation_sha256": mem.sha256_file(review_attestation_path),
        "review_packet_sha256": review_packet_sha256,
        "review_report_sha256": mem.sha256_file(report),
    }


def _verified_candidate_for_recheck(
    repo: Path,
    manifest: dict[str, Any],
    writer: dict[str, Any],
    prior_transaction_arg: str,
) -> dict[str, Any]:
    """Verify an original Codex writer artifact for a fresh current-policy review."""
    task_id = writer["task_id"]
    prior_transaction, prior_manifest = mem.load_transaction(repo, prior_transaction_arg)
    if prior_manifest.get("transaction_id") == manifest.get("transaction_id"):
        raise IsolationError("candidate reuse must come from another transaction")
    prior_writer = next(
        (item for item in prior_manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None,
    )
    if prior_writer is None or not prior_writer.get("required"):
        raise IsolationError(f"candidate transaction has no required writer task for {task_id}")
    if writer.get("source_unit_id") is None or prior_writer.get("source_unit_id") is None:
        raise IsolationError("candidate reuse is only supported for source-unit tasks")
    if (
        prior_writer.get("task_kind") != writer.get("task_kind")
        or prior_writer.get("source_unit_id") != writer.get("source_unit_id")
        or prior_writer.get("output_repository_path") != writer.get("output_repository_path")
        or prior_writer.get("staged_output_path") != writer.get("staged_output_path")
        or prior_writer.get("expected_output_name") != writer.get("expected_output_name")
    ):
        raise IsolationError("candidate reuse task identity or output path does not match")
    current_identity = _source_unit_identity(manifest, writer["source_unit_id"])
    prior_identity = _source_unit_identity(prior_manifest, prior_writer["source_unit_id"])
    if prior_identity != current_identity:
        raise IsolationError("candidate reuse source-unit member identities do not match")
    current_unit = next(item for item in manifest["units"] if item.get("id") == writer["source_unit_id"])
    prior_unit = next(item for item in prior_manifest["units"] if item.get("id") == writer["source_unit_id"])
    if prior_unit.get("kind") != current_unit.get("kind"):
        raise IsolationError("candidate reuse source-unit kind does not match")
    if prior_unit.get("unit_digest_sha256") != current_unit.get("unit_digest_sha256"):
        raise IsolationError("candidate reuse source-unit prerequisites changed")

    _, prior_packet_seal, _ = verify_packet(
        prior_transaction / "tasks" / task_id, prior_manifest["transaction_id"], task_id,
    )
    prior_seal_path = prior_transaction / "tasks" / task_id / "packet-seal.json"
    writer_attestation_path = prior_transaction / "attestations" / f"{task_id}.json"
    if not writer_attestation_path.is_file() or writer_attestation_path.is_symlink():
        raise IsolationError("candidate reuse lacks a regular Codex writer attestation")
    writer_attestation = json.loads(writer_attestation_path.read_text(encoding="utf-8"))
    expected_writer = {
        "transaction_id": prior_manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": prior_writer["source_unit_id"],
        "isolation": "bubblewrap",
        "workspace_hidden": repo.as_posix(),
        "packet_path": prior_writer["packet_path"],
        "packet_sha256": prior_packet_seal["combined_sha256"],
        "output_repository_path": prior_writer["output_repository_path"],
        "staged_output_path": prior_writer["staged_output_path"],
        "runtime_profile": "codex",
    }
    for key, value in expected_writer.items():
        if writer_attestation.get(key) != value:
            raise IsolationError(f"candidate reuse Codex attestation mismatch: {key}")
    for key in ("runner_sha256", "runtime_executable_sha256", "output_sha256"):
        _sha256_identity(writer_attestation.get(key), f"Codex {key}")
    if not isinstance(writer_attestation.get("runtime_version"), str) or not writer_attestation["runtime_version"]:
        raise IsolationError("candidate reuse Codex attestation lacks a runtime version")
    candidate = prior_transaction / prior_writer["staged_output_path"]
    if not candidate.is_file() or candidate.is_symlink():
        raise IsolationError("candidate reuse candidate is missing or is not a regular file")
    candidate_data = candidate.read_bytes()
    candidate_sha256 = mem.sha256_bytes(candidate_data)
    if candidate_sha256 != writer_attestation["output_sha256"]:
        raise IsolationError("candidate reuse candidate changed after its Codex attestation")

    return {
        "prior_transaction": prior_transaction,
        "prior_manifest": prior_manifest,
        "prior_writer": prior_writer,
        "prior_unit": prior_unit,
        "current_unit": current_unit,
        "member_identity_sha256": mem.sha256_bytes(mem.canonical_json(current_identity)),
        "candidate_data": candidate_data,
        "candidate_sha256": candidate_sha256,
        "prior_manifest_sha256": mem.sha256_file(prior_transaction / "manifest.json"),
        "prior_packet_seal_sha256": mem.sha256_file(prior_seal_path),
        "prior_packet_sha256": prior_packet_seal["combined_sha256"],
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "prior_runner_sha256": writer_attestation["runner_sha256"],
        "prior_runtime_executable_sha256": writer_attestation["runtime_executable_sha256"],
    }


def reuse_candidate_for_recheck(
    repo: Path, transaction_arg: str, task_id: str, prior_transaction_arg: str,
) -> dict[str, Any]:
    """Stage an attested Codex candidate for a mandatory fresh Grok review."""
    transaction, manifest = mem.load_transaction(repo, transaction_arg)
    writer = next((item for item in manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None)
    if writer is None or not writer.get("required"):
        raise IsolationError(f"transaction has no required writer task for {task_id}")
    if writer.get("source_unit_id") is None:
        raise IsolationError("candidate reuse is only supported for source-unit tasks")
    task_root = transaction / "tasks" / task_id
    packet, packet_seal, task_document = verify_packet(task_root, manifest["transaction_id"], task_id)
    del packet
    staged = transaction / writer["staged_output_path"]
    attestation_path = transaction / "attestations" / f"{task_id}.json"
    if list((task_root / "output").iterdir()) or staged.exists() or attestation_path.exists():
        raise IsolationError(f"current task already has output for {task_id}; prepare a new transaction")
    chain = _verified_candidate_for_recheck(repo, manifest, writer, prior_transaction_arg)
    normalized, raw_output_sha256 = normalize_source_capsule_frontmatter(
        chain["candidate_data"], task_document, manifest, writer,
    )
    output_sha256 = mem.sha256_bytes(normalized)
    mem.atomic_write(staged, normalized, 0o644)
    reuse_chain = {
        "prior_transaction_id": chain["prior_manifest"]["transaction_id"],
        "prior_target_commit": chain["prior_manifest"]["target_commit"],
        "current_target_commit": manifest["target_commit"],
        "prior_manifest_sha256": chain["prior_manifest_sha256"],
        "prior_packet_seal_sha256": chain["prior_packet_seal_sha256"],
        "prior_packet_sha256": chain["prior_packet_sha256"],
        "prior_candidate_sha256": chain["candidate_sha256"],
        "prior_writer_attestation_sha256": chain["writer_attestation_sha256"],
        "prior_runner_sha256": chain["prior_runner_sha256"],
        "prior_runtime_executable_sha256": chain["prior_runtime_executable_sha256"],
        "prior_unit_digest_sha256": chain["prior_unit"]["unit_digest_sha256"],
        "current_unit_digest_sha256": chain["current_unit"]["unit_digest_sha256"],
        "member_identities_sha256": chain["member_identity_sha256"],
    }
    attestation = {
        "attestation_version": 1,
        "transaction_id": manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": writer["source_unit_id"],
        "isolation": "deterministic-candidate-reuse",
        "workspace_hidden": repo.as_posix(),
        "packet_path": writer["packet_path"],
        "packet_sha256": packet_seal["combined_sha256"],
        "output_repository_path": writer["output_repository_path"],
        "staged_output_path": writer["staged_output_path"],
        "output_sha256": output_sha256,
        "raw_output_sha256": raw_output_sha256,
        "normalization": "source_capsule_frontmatter_v1",
        "runner_sha256": mem.sha256_file(Path(__file__).resolve()),
        "runtime_profile": "codex-candidate-reuse",
        "runtime_version": None,
        "runtime_executable_sha256": None,
        "model_invoked": False,
        "candidate_reuse": reuse_chain,
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    mem.write_json(attestation_path, attestation)
    return attestation


def reuse_reviewed_candidate(
    repo: Path, transaction_arg: str, task_id: str, prior_transaction_arg: str,
) -> dict[str, Any]:
    """Deterministically normalize and stage a prior Codex+Grok-PASS candidate."""
    transaction, manifest = mem.load_transaction(repo, transaction_arg)
    writer = next((item for item in manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None)
    if writer is None or not writer.get("required"):
        raise IsolationError(f"transaction has no required writer task for {task_id}")
    if writer.get("source_unit_id") is None:
        raise IsolationError("reviewed reuse is only supported for source-unit tasks")
    task_root = transaction / "tasks" / task_id
    packet, packet_seal, task_document = verify_packet(task_root, manifest["transaction_id"], task_id)
    del packet
    staged = transaction / writer["staged_output_path"]
    attestation_path = transaction / "attestations" / f"{task_id}.json"
    if list((task_root / "output").iterdir()) or staged.exists() or attestation_path.exists():
        raise IsolationError(f"current task already has output for {task_id}; prepare a new transaction")
    chain = _verified_reviewed_candidate(repo, manifest, writer, prior_transaction_arg)
    normalized, raw_output_sha256 = normalize_source_capsule_frontmatter(
        chain["candidate_data"], task_document, manifest, writer,
    )
    output_sha256 = mem.sha256_bytes(normalized)
    mem.atomic_write(staged, normalized, 0o644)
    reuse_chain = {
        "prior_transaction_id": chain["prior_manifest"]["transaction_id"],
        "prior_target_commit": chain["prior_manifest"]["target_commit"],
        "current_target_commit": manifest["target_commit"],
        "prior_packet_sha256": chain["prior_packet_sha256"],
        "prior_candidate_sha256": chain["candidate_sha256"],
        "prior_writer_attestation_sha256": chain["writer_attestation_sha256"],
        "prior_review_attestation_sha256": chain["review_attestation_sha256"],
        "prior_review_packet_sha256": chain["review_packet_sha256"],
        "prior_review_report_sha256": chain["review_report_sha256"],
        "prior_unit_digest_sha256": chain["prior_unit"]["unit_digest_sha256"],
        "current_unit_digest_sha256": chain["current_unit"]["unit_digest_sha256"],
        "member_identities_sha256": chain["member_identity_sha256"],
    }
    attestation = {
        "attestation_version": 1,
        "transaction_id": manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": writer["source_unit_id"],
        "isolation": "deterministic-reviewed-reuse",
        "workspace_hidden": repo.as_posix(),
        "packet_path": writer["packet_path"],
        "packet_sha256": packet_seal["combined_sha256"],
        "output_repository_path": writer["output_repository_path"],
        "staged_output_path": writer["staged_output_path"],
        "output_sha256": output_sha256,
        "raw_output_sha256": raw_output_sha256,
        "normalization": "source_capsule_frontmatter_v1",
        "runner_sha256": mem.sha256_file(Path(__file__).resolve()),
        "runtime_profile": "codex-reviewed-reuse",
        "runtime_version": None,
        "runtime_executable_sha256": None,
        "model_invoked": False,
        "reviewed_reuse": reuse_chain,
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    mem.write_json(attestation_path, attestation)
    return attestation


def normalize_source_capsule_frontmatter(
    data: bytes, task_document: dict[str, Any], manifest: dict[str, Any], writer: dict[str, Any],
) -> tuple[bytes, str | None]:
    """Replace machine-owned source metadata with the frozen transaction values."""
    source_unit_id = writer.get("source_unit_id")
    if source_unit_id is None:
        return data, None
    front, body = mem.parse_frontmatter(data, writer["output_repository_path"])
    unit = task_document.get("source_unit")
    contract = task_document.get("semantic_contract")
    if not isinstance(unit, dict) or not isinstance(contract, dict):
        raise IsolationError("source capsule task lacks frozen unit or semantic contract")
    page = contract.get("page")
    if not isinstance(page, dict):
        raise IsolationError("source capsule task lacks frozen page metadata")
    members = [
        {key: member.get(key) for key in (
            "path", "role", "read_mode", "mode", "object_type", "blob_oid", "blob_size",
        )}
        for member in unit.get("members", [])
    ]
    citations = contract.get("allowed_citations", [])
    sources = sorted({item.get("path") for item in citations if isinstance(item, dict) and isinstance(item.get("path"), str)})
    normalized = dict(front)
    normalized.update({
        "schema_version": 2,
        "id": page["id"],
        "title": page["title"],
        "type": page["type"],
        "lifecycle": page["desired_lifecycle"],
        "memory_review": "ai_draft",
        "sources": sources,
        "content_owner": page["content_owner"],
        "last_updated": contract["refresh_date"],
        "generated_from_commit": manifest["target_commit"],
        "source_kind": contract["source_kind"],
        "source_unit": {
            "id": source_unit_id,
            "shape": contract["shape"],
            "entrypoint": contract["entrypoint"],
            "unit_digest_sha256": unit["unit_digest_sha256"],
            "members": members,
        },
        "extractor_version": contract["extractor_version"],
    })
    rendered = ("---\n" + yaml.safe_dump(normalized, sort_keys=False, allow_unicode=True) + "---\n" + body).encode("utf-8")
    if rendered == data:
        return data, None
    return rendered, mem.sha256_bytes(data)


def run_task(
    repo: Path, transaction_arg: str, task_id: str, profile_name: str | Sequence[str],
    revise_from: str | None = None, revise_lint_from: str | None = None,
) -> dict[str, Any]:
    if revise_from is not None and revise_lint_from is not None:
        raise IsolationError("review revision and lint revision modes are mutually exclusive")
    if (revise_from is not None or revise_lint_from is not None) and not isinstance(profile_name, str):
        raise IsolationError("revision runs require the named codex writer profile")
    profile = runtime_profile(profile_name) if isinstance(profile_name, str) else _test_profile(profile_name)
    transaction, manifest = mem.load_transaction(repo, transaction_arg)
    writer = next((item for item in manifest.get("writer_tasks", []) if item["task_id"] == task_id), None)
    if writer is None or not writer.get("required"):
        raise IsolationError(f"transaction has no required writer task for {task_id}")
    task_root = transaction / "tasks" / task_id
    packet, packet_seal, task_document = verify_packet(task_root, manifest["transaction_id"], task_id)
    dynamic_inputs = task_document.get("semantic_contract", {}).get("dynamic_memory_inputs", [])
    needs_hydration = bool(dynamic_inputs) and not (packet / "hydration.json").is_file()
    hydration_records: list[dict[str, Any]] = []
    for dependency in dynamic_inputs if needs_hydration else []:
        dependency_id = dependency["task_id"]
        dependency_writer = next(
            (item for item in manifest.get("writer_tasks", []) if item.get("task_id") == dependency_id), None
        )
        if dependency_writer is None or not dependency_writer.get("required"):
            raise IsolationError(f"prerequisite task is not declared as a required writer: {dependency_id}")
        dependency_attestation_path = transaction / "attestations" / f"{dependency_id}.json"
        if not dependency_attestation_path.is_file():
            raise IsolationError(f"prerequisite task has no isolation attestation: {dependency_id}")
        dependency_attestation = json.loads(dependency_attestation_path.read_text(encoding="utf-8"))
        if (
            dependency_attestation.get("isolation") != "bubblewrap"
            or dependency_attestation.get("transaction_id") != manifest["transaction_id"]
            or dependency_attestation.get("task_id") != dependency_id
            or dependency_attestation.get("output_repository_path") != dependency.get("page")
            or dependency_writer.get("output_repository_path") != dependency.get("page")
        ):
            raise IsolationError(f"prerequisite task was not isolated: {dependency_id}")
        dependency_task_root = transaction / "tasks" / dependency_id
        dependency_packet_seal = json.loads(
            (dependency_task_root / "packet-seal.json").read_text(encoding="utf-8")
        )
        if (
            dependency_attestation.get("packet_sha256") != dependency_packet_seal.get("combined_sha256")
            or dependency_attestation.get("runner_sha256") != mem.sha256_file(Path(__file__).resolve())
        ):
            raise IsolationError(f"prerequisite attestation identity mismatch: {dependency_id}")
        dependency_staged = transaction / dependency_attestation["staged_output_path"]
        if not dependency_staged.is_file() or mem.sha256_file(dependency_staged) != dependency_attestation.get("output_sha256"):
            raise IsolationError(f"prerequisite staged page changed: {dependency_id}")
        if dependency.get("policy_fresh") is not True or dependency.get("policy_sha256") != manifest["policy"]["combined_sha256"]:
            raise IsolationError(f"prerequisite policy metadata is stale: {dependency_id}")
        dependency_front, _ = mem.parse_frontmatter(dependency_staged.read_bytes(), dependency["page"])
        if dependency_front.get("generated_from_commit") != dependency.get("generated_commit"):
            raise IsolationError(f"prerequisite generated commit mismatch: {dependency_id}")
        if dependency_front.get("lifecycle") != dependency.get("lifecycle"):
            raise IsolationError(f"prerequisite lifecycle mismatch: {dependency_id}")
        expected_unit_digest = dependency.get("unit_digest_sha256")
        if expected_unit_digest is not None:
            source_unit = dependency_front.get("source_unit")
            actual_digest = source_unit.get("unit_digest_sha256") if isinstance(source_unit, dict) else None
            if actual_digest != expected_unit_digest:
                raise IsolationError(f"prerequisite unit digest mismatch: {dependency_id}")
        destination = packet / "memory_inputs" / dependency["page"]
        if destination.exists():
            raise IsolationError(f"dynamic memory input already exists; refusing to overwrite: {dependency['page']}")
        mem.atomic_write(destination, dependency_staged.read_bytes(), 0o444)
        attestation_sha256 = mem.sha256_file(dependency_attestation_path)
        enriched_fields = {
            "page_sha256": dependency_attestation["output_sha256"],
            "attestation_sha256": attestation_sha256,
            "attestation_packet_sha256": dependency_attestation["packet_sha256"],
            "attestation_identity": {
                "transaction_id": manifest["transaction_id"], "task_id": dependency_id,
            },
        }
        dependency.update(enriched_fields)
        duplicate_dependency = next(
            item for item in task_document["writer_task"]["semantic_contract"]["dynamic_memory_inputs"]
            if item["task_id"] == dependency_id
        )
        duplicate_dependency.update(enriched_fields)
        hydration_records.append({
            "task_id": dependency_id,
            "page": dependency["page"],
            "staged_output_path": dependency_attestation["staged_output_path"],
            "page_sha256": dependency_attestation["output_sha256"],
            "attestation_sha256": attestation_sha256,
            "attestation_packet_sha256": dependency_attestation["packet_sha256"],
        })
    if needs_hydration:
        initial_packet_sha256 = packet_seal["combined_sha256"]
        mem.atomic_write(packet / "task.json", mem.canonical_json(task_document), 0o444)
        mem.atomic_write(packet / "hydration.json", mem.canonical_json({
            "hydration_version": 1,
            "initial_packet_sha256": initial_packet_sha256,
            "dependencies": hydration_records,
        }), 0o444)
        packet_files, combined = packet_digest(packet)
        packet_seal = {
            "transaction_id": manifest["transaction_id"],
            "task_id": task_id,
            "source_unit_id": writer.get("source_unit_id"),
            "files": packet_files,
            "combined_sha256": combined,
        }
        mem.write_json(task_root / "packet-seal.json", packet_seal)
        mem.load_transaction(repo, str(transaction))
        packet, packet_seal, task_document = verify_packet(task_root, manifest["transaction_id"], task_id)
    if revise_from is not None:
        packet, packet_seal, task_document = admit_revision_candidate(
            repo, transaction, manifest, writer, packet, packet_seal, revise_from,
        )
        if not profile.command or profile.command[-1] != WRITER_PROMPT:
            raise IsolationError("codex writer command does not have the expected prompt boundary")
        profile = dataclasses.replace(profile, command=(*profile.command[:-1], REVISION_PROMPT))
    elif revise_lint_from is not None:
        packet, packet_seal, task_document = admit_lint_revision_candidate(
            repo, transaction, manifest, writer, packet, packet_seal, revise_lint_from,
        )
        if not profile.command or profile.command[-1] != WRITER_PROMPT:
            raise IsolationError("codex writer command does not have the expected prompt boundary")
        profile = dataclasses.replace(profile, command=(*profile.command[:-1], LINT_REVISION_PROMPT))
    output = task_root / "output"
    existing = list(output.iterdir())
    if existing:
        raise IsolationError(f"isolated output is not empty for {task_id}; refusing to overwrite")
    sandbox_command = bubblewrap_command(repo, packet, output, profile)
    try:
        proc = subprocess.run(
            sandbox_command, check=False, input=profile.stdin_data,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=CODEX_WRITER_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as exc:
        _preserve_timeout_diagnostics(output, exc)
        partial_stderr = _timeout_output(exc.stderr).decode("utf-8", "replace").strip()
        suffix = f": {partial_stderr[-4000:]}" if partial_stderr else ""
        raise IsolationError(
            f"isolated Codex writer timed out after {CODEX_WRITER_TIMEOUT_SECONDS} seconds{suffix}"
        ) from exc
    if proc.returncode:
        diagnostic = (proc.stderr or proc.stdout or b"").decode("utf-8", "replace").strip()
        if len(diagnostic) > 4000:
            diagnostic = diagnostic[-4000:]
        suffix = f": {diagnostic}" if diagnostic else ""
        raise IsolationError(f"isolated writer exited with status {proc.returncode}{suffix}")
    if profile.capture_stdout:
        if list(output.iterdir()):
            raise IsolationError(f"{profile.name} profile wrote undeclared output while stdout capture was active")
        mem.atomic_write(output / writer["expected_output_name"], proc.stdout or b"", 0o600)
    produced = list(output.iterdir())
    expected_name = writer["expected_output_name"]
    if len(produced) != 1 or produced[0].name != expected_name:
        names = sorted(item.name for item in produced)
        raise IsolationError(f"isolated writer must create exactly {expected_name}; found {names}")
    generated = produced[0]
    if generated.is_symlink() or not generated.is_file():
        raise IsolationError("isolated output must be one regular file, not a link or directory")
    data = generated.read_bytes()
    raw_output_sha256 = None
    if isinstance(profile_name, str):
        data, raw_output_sha256 = normalize_source_capsule_frontmatter(data, task_document, manifest, writer)
    digest = mem.sha256_bytes(data)
    staged_relative = writer["staged_output_path"]
    staged = transaction / staged_relative
    mem.atomic_write(staged, data, 0o644)
    attestation = {
        "attestation_version": 1,
        "transaction_id": manifest["transaction_id"],
        "task_id": task_id,
        "source_unit_id": writer.get("source_unit_id"),
        "isolation": "bubblewrap",
        "workspace_hidden": repo.as_posix(),
        "packet_path": writer["packet_path"],
        "packet_sha256": packet_seal["combined_sha256"],
        "output_repository_path": writer["output_repository_path"],
        "staged_output_path": staged_relative,
        "output_sha256": digest,
        "raw_output_sha256": raw_output_sha256,
        "normalization": "source_capsule_frontmatter_v1" if raw_output_sha256 else None,
        "runner_sha256": mem.sha256_file(Path(__file__).resolve()),
        "runtime_profile": profile.name,
        "runtime_version": profile.version,
        "runtime_executable_sha256": profile.executable_sha256,
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    attestation_path = transaction / "attestations" / f"{task_id}.json"
    mem.write_json(attestation_path, attestation)
    return attestation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, help="repository root (normally auto-detected)")
    parser.add_argument("transaction", help="transaction ID or path")
    parser.add_argument("--task", required=True, help="writer task ID (source capsule or derived page)")
    parser.add_argument("--profile", choices=("codex",), help="isolated Codex writer profile")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--revise-from", metavar="TRANSACTION",
        help="revise from the same task's verified Codex candidate in another transaction",
    )
    mode.add_argument(
        "--revise-lint-from", metavar="TRANSACTION",
        help="revise a Grok-PASS source candidate from its trusted recorded staged-lint failure",
    )
    mode.add_argument(
        "--reuse-reviewed-from", metavar="TRANSACTION",
        help="without a model call, normalize and reuse a prior Codex candidate with a Grok PASS review",
    )
    mode.add_argument(
        "--reuse-candidate-from", metavar="TRANSACTION",
        help="without a model call, normalize a prior Codex candidate for a fresh current Grok review",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    deterministic_mode = args.reuse_reviewed_from is not None or args.reuse_candidate_from is not None
    if deterministic_mode and args.profile is not None:
        parser.error("deterministic reuse must not be combined with --profile")
    if not deterministic_mode and args.profile is None:
        parser.error("--profile codex is required unless a deterministic reuse mode is used")
    try:
        repo = args.repo.resolve() if args.repo else mem.find_repo()
        recovered = mem.recover_publication(repo)
        if recovered:
            print(f"recovery: {recovered}", file=sys.stderr)
        if args.reuse_reviewed_from is not None:
            result = reuse_reviewed_candidate(repo, args.transaction, args.task, args.reuse_reviewed_from)
        elif args.reuse_candidate_from is not None:
            result = reuse_candidate_for_recheck(repo, args.transaction, args.task, args.reuse_candidate_from)
        else:
            result = run_task(
                repo, args.transaction, args.task, args.profile,
                args.revise_from, args.revise_lint_from,
            )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (mem.MemoryErrorBase, OSError, ValueError) as exc:
        print(f"memory isolation: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
