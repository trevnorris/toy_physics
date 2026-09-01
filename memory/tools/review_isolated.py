#!/usr/bin/env python3
"""Run Grok as a read-only checker of one Codex-generated memory page.

This command has no publication path. It mounts the original sealed writer
packet and the attested staged page read-only, hides the live repository, and
stores only a review report plus a reviewer attestation under the transaction.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Sequence

try:
    from . import memory as mem
    from . import run_isolated as writer_runner
    from .review_contract import (
        GROK_REVIEW_MODEL, REVIEW_CONTRACT_SHA256, REVIEW_PROMPT, REVIEW_PROMPT_SHA256,
        REVIEW_SCHEMA, REVIEW_SCHEMA_SHA256,
    )
except ImportError:  # pragma: no cover - direct CLI execution
    import memory as mem  # type: ignore[no-redef]
    import run_isolated as writer_runner  # type: ignore[no-redef]
    from review_contract import (  # type: ignore[no-redef]
        GROK_REVIEW_MODEL, REVIEW_CONTRACT_SHA256, REVIEW_PROMPT, REVIEW_PROMPT_SHA256,
        REVIEW_SCHEMA, REVIEW_SCHEMA_SHA256,
    )


class ReviewError(mem.MemoryErrorBase):
    pass


# Hard wall-clock limit for one Grok checking invocation.  Timeout artifacts
# remain an incomplete review and can be archived/retried with --retry-incomplete.
GROK_REVIEW_TIMEOUT_SECONDS = 15 * 60


def grok_profile() -> writer_runner.RuntimeProfile:
    home = Path.home()
    wrapper = writer_runner._required_file(Path(shutil.which("grok") or ""), "Grok CLI")
    executable = writer_runner._required_file(wrapper.resolve(), "native Grok runtime")
    auth = writer_runner._required_file(home / ".grok/auth.json", "Grok credential file")
    agent_id = writer_runner._required_file(home / ".grok/agent_id", "Grok agent identity")
    version = writer_runner._probe_version([str(executable), "--version"])
    network_files = tuple(
        (str(path), str(path)) for path in (
            Path("/etc/resolv.conf"), Path("/etc/hosts"), Path("/etc/nsswitch.conf"),
            Path("/etc/ssl/certs/ca-certificates.crt"),
        ) if path.is_file()
    )
    libraries = tuple(
        (path, path) for path in (
            "/lib64/ld-linux-x86-64.so.2", "/lib/x86_64-linux-gnu/libc.so.6",
            "/lib/x86_64-linux-gnu/libm.so.6", "/lib/x86_64-linux-gnu/libpthread.so.0",
            "/lib/x86_64-linux-gnu/libdl.so.2", "/lib/x86_64-linux-gnu/librt.so.1",
            "/lib/x86_64-linux-gnu/libgcc_s.so.1",
        ) if Path(path).is_file()
    )
    return writer_runner.RuntimeProfile(
        name="grok-review",
        version=version,
        executable_sha256=mem.sha256_file(executable),
        command=(
            "/runtime/grok", "--verbatim", "--no-subagents", "--disable-web-search",
            "--model", GROK_REVIEW_MODEL,
            "--output-format", "json", "--cwd", "/packet", "--permission-mode", "bypassPermissions",
            "--json-schema", REVIEW_SCHEMA, "--prompt-file", "/packet/review-prompt.md",
        ),
        ro_binds=(
            (str(executable), "/runtime/grok"),
            (str(auth.resolve()), "/runtime-home/.grok/auth.json"),
            (str(agent_id.resolve()), "/runtime-home/.grok/agent_id"),
        ) + libraries + network_files,
        environment=(
            ("HOME", "/runtime-home"), ("PATH", "/runtime"), ("LANG", "C.UTF-8"),
            ("SSL_CERT_FILE", "/etc/ssl/certs/ca-certificates.crt"),
        ),
        capture_stdout=True,
    )


def compose_review_prompt(review_packet: Path) -> bytes:
    """Inline all allowed evidence for Grok's single-turn, tool-free review."""
    task_path = review_packet / "task.json"
    task = json.loads(task_path.read_text(encoding="utf-8"))
    contract = task.get("semantic_contract", {})
    prompt_paths = contract.get("prompt_paths", [])
    semantic_inputs = task.get("writer_task", {}).get("semantic_inputs", [])
    sections = [REVIEW_PROMPT, "", "# Immutable task", "", task_path.read_text(encoding="utf-8")]
    sections.extend(["", "# Authoritative schema", "", (review_packet / "schema.md").read_text(encoding="utf-8")])
    for prompt_path in prompt_paths:
        if not isinstance(prompt_path, str) or not prompt_path.startswith("/packet/"):
            raise ReviewError("review task contains an invalid frozen prompt path")
        relative = Path(prompt_path.removeprefix("/packet/"))
        frozen = review_packet / relative
        if not frozen.is_file():
            raise ReviewError(f"review task frozen prompt is missing: {relative.as_posix()}")
        sections.extend(["", f"# Frozen author prompt: {relative.as_posix()}", "", frozen.read_text(encoding="utf-8")])
    for semantic in semantic_inputs:
        packet_path = semantic.get("packet_path") if isinstance(semantic, dict) else None
        source_path = semantic.get("source_path") if isinstance(semantic, dict) else None
        if not isinstance(packet_path, str) or not isinstance(source_path, str):
            raise ReviewError("review task contains an invalid semantic input record")
        frozen = review_packet / packet_path
        if not frozen.is_file():
            raise ReviewError(f"review semantic input is missing: {packet_path}")
        sections.extend(["", f"# Sealed source: {source_path}", "", frozen.read_text(encoding="utf-8")])
    sections.extend(["", "# Candidate under review", "", (review_packet / "candidate.md").read_text(encoding="utf-8")])
    rendered = ("\n".join(sections).rstrip() + "\n").encode("utf-8")
    config = mem.load_yaml_bytes((review_packet / "config.yaml").read_bytes(), "frozen review config")
    limit = int(config.get("read_limits", {}).get("review_prompt_max_bytes", 12_000_000))
    if len(rendered) > limit:
        raise ReviewError(f"review prompt is {len(rendered)} bytes, above {limit}-byte limit")
    return rendered


def _failure_detail(proc: subprocess.CompletedProcess[bytes]) -> str:
    raw = proc.stderr or proc.stdout or b""
    diagnostic = raw.decode("utf-8", "replace").strip().replace("\x00", "")
    return diagnostic[-2000:] if diagnostic else "no diagnostic output"


def _parse_structured_review(report_text: str) -> dict[str, object]:
    decoder = json.JSONDecoder()
    index = 0
    objects: list[dict[str, object]] = []
    while index < len(report_text):
        while index < len(report_text) and report_text[index].isspace():
            index += 1
        if index >= len(report_text):
            break
        try:
            value, index = decoder.raw_decode(report_text, index)
        except json.JSONDecodeError as exc:
            raise ReviewError("Grok structured review text was not JSON") from exc
        if not isinstance(value, dict):
            raise ReviewError("Grok structured review contained a non-object value")
        objects.append(value)
    if not objects:
        raise ReviewError("Grok structured review was empty")
    structured = objects[-1]
    verdict = structured.get("verdict")
    findings = structured.get("findings")
    if verdict not in ("PASS", "FAIL") or not isinstance(findings, list):
        raise ReviewError("Grok structured review did not satisfy the verdict/findings contract")
    if {item.get("verdict") for item in objects[:-1]} - {verdict}:
        raise ReviewError("Grok emitted contradictory structured verdicts")
    if verdict == "FAIL" and not any(
        finding.get("severity") in ("blocking", "major")
        for finding in findings if isinstance(finding, dict)
    ):
        raise ReviewError("Grok returned FAIL without a blocking or major actionable finding")
    if verdict == "PASS" and any(
        finding.get("severity") in ("blocking", "major")
        for finding in findings if isinstance(finding, dict)
    ):
        raise ReviewError("Grok returned PASS with a blocking or major finding")
    return structured


def _structured_from_envelope(envelope: object) -> dict[str, object]:
    if not isinstance(envelope, dict):
        raise ReviewError("Grok output envelope was not an object")
    final_structured = envelope.get("structuredOutput")
    if isinstance(final_structured, dict):
        # Grok may stream provisional schema objects through ``text`` while it
        # reasons.  ``structuredOutput`` is the CLI's final JSON-schema result
        # and is therefore the only authoritative verdict when present.
        return _parse_structured_review(json.dumps(final_structured, separators=(",", ":")))
    report_text = envelope.get("text")
    if not isinstance(report_text, str):
        raise ReviewError("Grok JSON output contained neither structuredOutput nor a text review")
    return _parse_structured_review(report_text)


def _completed_review_evidence(review_root: Path) -> str | None:
    """Return the first recognizable PASS/FAIL artifact in an existing review."""
    attestation_path = review_root / "attestation.json"
    if attestation_path.is_file():
        try:
            attestation = json.loads(attestation_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            pass
        else:
            if isinstance(attestation, dict) and attestation.get("verdict") in ("PASS", "FAIL"):
                return f"{attestation.get('verdict')} attestation"

    report_path = review_root / "output" / "report.md"
    if report_path.is_file():
        try:
            first_line = report_path.read_text(encoding="utf-8").splitlines()[0]
        except (OSError, UnicodeError, IndexError):
            pass
        else:
            if first_line in ("Verdict: PASS", "Verdict: FAIL"):
                return first_line

    rejected_path = review_root / "output" / "rejected-output.json"
    if rejected_path.is_file():
        try:
            envelope = json.loads(rejected_path.read_text(encoding="utf-8"))
            structured = _structured_from_envelope(envelope)
        except (OSError, UnicodeError, json.JSONDecodeError, ReviewError):
            pass
        else:
            return f"{structured['verdict']} structured verdict"
    return None


def _archive_incomplete_review(transaction: Path, task_id: str, review_root: Path) -> Path:
    """Atomically move an incomplete active attempt to its next stable archive name."""
    if not review_root.is_dir() or review_root.is_symlink():
        raise ReviewError("incomplete active review is not a regular directory")
    archive_root = transaction / "rejected-reviews"
    if archive_root.is_symlink() or (archive_root.exists() and not archive_root.is_dir()):
        raise ReviewError("incomplete-review archive root is not a regular directory")
    archive_root.mkdir(exist_ok=True)
    attempt = 1
    while True:
        archive = archive_root / f"{task_id}-attempt-{attempt:04d}"
        if not archive.exists() and not archive.is_symlink():
            break
        attempt += 1
    review_root.rename(archive)
    return archive


def run_review(
    repo: Path, transaction_arg: str, task_id: str, *, recover_rejected: bool = False,
    retry_incomplete: bool = False,
) -> dict[str, object]:
    if recover_rejected and retry_incomplete:
        raise ReviewError("--recover-rejected and --retry-incomplete are mutually exclusive")
    transaction, manifest = mem.load_transaction(repo, transaction_arg)
    writer = next((item for item in manifest.get("writer_tasks", []) if item.get("task_id") == task_id), None)
    if writer is None or not writer.get("required"):
        raise ReviewError(f"transaction has no required writer task for {task_id}")
    packet, packet_seal, _ = writer_runner.verify_packet(
        transaction / "tasks" / task_id, manifest["transaction_id"], task_id,
    )
    writer_attestation_path = transaction / "attestations" / f"{task_id}.json"
    if not writer_attestation_path.is_file():
        raise ReviewError(f"Codex writer attestation is missing for {task_id}")
    writer_attestation = json.loads(writer_attestation_path.read_text(encoding="utf-8"))
    writer_profile = writer_attestation.get("runtime_profile")
    if writer_profile not in ("codex", "codex-candidate-reuse") or writer_attestation.get("task_id") != task_id:
        raise ReviewError(f"candidate was not produced by the isolated Codex writer: {task_id}")
    candidate = transaction / writer["staged_output_path"]
    if not candidate.is_file() or mem.sha256_file(candidate) != writer_attestation.get("output_sha256"):
        raise ReviewError(f"attested staged candidate is missing or changed for {task_id}")
    if writer_profile == "codex-candidate-reuse":
        expected_reuse = {
            "transaction_id": manifest["transaction_id"],
            "task_id": task_id,
            "source_unit_id": writer.get("source_unit_id"),
            "isolation": "deterministic-candidate-reuse",
            "workspace_hidden": repo.as_posix(),
            "packet_path": writer["packet_path"],
            "packet_sha256": packet_seal["combined_sha256"],
            "output_repository_path": writer["output_repository_path"],
            "staged_output_path": writer["staged_output_path"],
            "model_invoked": False,
        }
        for key, value in expected_reuse.items():
            if writer_attestation.get(key) != value:
                raise ReviewError(f"candidate reuse attestation mismatch for {task_id}: {key}")
        reuse_errors = mem.verify_candidate_reuse_chain(
            repo, transaction, manifest, writer, writer_attestation, mem.sha256_file(candidate),
        )
        if reuse_errors:
            raise ReviewError(reuse_errors[0])

    review_root = transaction / "reviews" / task_id
    output = review_root / "output"
    review_packet = review_root / "packet"
    rejected_path = output / "rejected-output.json"
    retry_archive: Path | None = None
    if recover_rejected:
        if not review_packet.is_dir() or not rejected_path.is_file() or (review_root / "attestation.json").exists():
            raise ReviewError(f"no recoverable rejected Grok output exists for {task_id}")
    else:
        if retry_incomplete:
            if not review_root.exists():
                raise ReviewError(f"no incomplete Grok review exists for {task_id}")
            completed_evidence = _completed_review_evidence(review_root)
            if completed_evidence is not None:
                raise ReviewError(
                    f"review for {task_id} contains completed {completed_evidence}; refusing to overwrite"
                )
            retry_archive = _archive_incomplete_review(transaction, task_id, review_root)
        elif review_root.exists():
            raise ReviewError(f"review output already exists for {task_id}; refusing to overwrite")
        shutil.copytree(packet, review_packet, copy_function=shutil.copy2)
        shutil.copy2(candidate, review_packet / "candidate.md")
        try:
            prompt_data = compose_review_prompt(review_packet)
        except Exception:
            shutil.rmtree(review_root)
            raise
        mem.atomic_write(review_packet / "review-prompt.md", prompt_data, 0o444)
        for path in review_packet.rglob("*"):
            if path.is_file():
                path.chmod(0o444)
        output.mkdir(parents=True)
    review_packet_files, review_packet_sha256 = writer_runner.packet_digest(review_packet)
    expected_review_files, _ = writer_runner.packet_digest(packet)
    expected_review_files["candidate.md"] = mem.sha256_file(candidate)
    expected_review_files["review-prompt.md"] = mem.sha256_file(review_packet / "review-prompt.md")
    expected_review_sha256 = mem.sha256_bytes(mem.canonical_json(expected_review_files))
    if review_packet_files != expected_review_files or review_packet_sha256 != expected_review_sha256:
        raise ReviewError("derived review packet does not match the sealed writer packet plus candidate")
    profile = grok_profile()
    if recover_rejected:
        raw_stdout = rejected_path.read_bytes()
    else:
        command = writer_runner.bubblewrap_command(repo, review_packet, output, profile)
        try:
            proc = subprocess.run(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                timeout=GROK_REVIEW_TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired as exc:
            stdout = writer_runner._timeout_output(exc.stdout)
            stderr = writer_runner._timeout_output(exc.stderr)
            mem.atomic_write(output / "failed-stdout.bin", stdout, 0o600)
            mem.atomic_write(output / "failed-stderr.bin", stderr, 0o600)
            diagnostic = stderr.decode("utf-8", "replace").strip().replace("\x00", "")
            suffix = f": {diagnostic[-2000:]}" if diagnostic else ""
            raise ReviewError(
                f"isolated Grok reviewer timed out after {GROK_REVIEW_TIMEOUT_SECONDS} seconds{suffix}"
            ) from exc
        if proc.returncode:
            mem.atomic_write(output / "failed-stdout.bin", proc.stdout or b"", 0o600)
            mem.atomic_write(output / "failed-stderr.bin", proc.stderr or b"", 0o600)
            raise ReviewError(f"isolated Grok reviewer exited with status {proc.returncode}: {_failure_detail(proc)}")
        raw_stdout = proc.stdout or b""
    try:
        envelope = json.loads(raw_stdout.decode("utf-8", "strict"))
        structured = _structured_from_envelope(envelope)
    except (UnicodeError, json.JSONDecodeError, ReviewError):
        mem.atomic_write(output / "rejected-output.json", raw_stdout, 0o600)
        raise
    rendered = [f"Verdict: {structured['verdict']}", ""]
    if not structured["findings"]:
        rendered.append("No findings.")
    for index, finding in enumerate(structured["findings"], 1):
        rendered.extend([
            f"## {index}. {finding['severity'].title()}: {finding['summary']}", "",
            f"Candidate: {finding['candidate_location']}", "",
            f"Evidence: {finding['source_evidence']}", "",
        ])
    report_data = ("\n".join(rendered).rstrip() + "\n").encode("utf-8")
    report = output / "report.md"
    mem.atomic_write(report, report_data, 0o600)
    attestation: dict[str, object] = {
        "attestation_version": 1,
        "role": "independent_review",
        "transaction_id": manifest["transaction_id"],
        "task_id": task_id,
        "packet_sha256": packet_seal["combined_sha256"],
        "review_packet_sha256": review_packet_sha256,
        "review_packet_files": review_packet_files,
        "candidate_path": writer["staged_output_path"],
        "candidate_sha256": mem.sha256_file(candidate),
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "report_path": report.relative_to(transaction).as_posix(),
        "report_sha256": mem.sha256_file(report),
        "verdict": structured["verdict"],
        "runtime_profile": profile.name,
        "review_model": GROK_REVIEW_MODEL,
        "review_prompt_sha256": REVIEW_PROMPT_SHA256,
        "review_schema_sha256": REVIEW_SCHEMA_SHA256,
        "review_contract_sha256": REVIEW_CONTRACT_SHA256,
        "runtime_version": profile.version,
        "runtime_executable_sha256": profile.executable_sha256,
        "reviewer_sha256": mem.sha256_file(Path(__file__).resolve()),
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "recovered_from_rejected_output": recover_rejected,
        "retried_from_incomplete_archive": (
            retry_archive.relative_to(transaction).as_posix() if retry_archive is not None else None
        ),
    }
    mem.write_json(review_root / "attestation.json", attestation)
    return attestation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, help="repository root (normally auto-detected)")
    parser.add_argument("transaction", help="transaction ID or path")
    parser.add_argument("--task", required=True, help="attested Codex writer task ID")
    recovery = parser.add_mutually_exclusive_group()
    recovery.add_argument(
        "--recover-rejected", action="store_true",
        help="parse a preserved rejected Grok envelope without another model call",
    )
    recovery.add_argument(
        "--retry-incomplete", action="store_true",
        help="archive an incomplete infrastructure-failed attempt and start a clean review",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        repo = args.repo.resolve() if args.repo else mem.find_repo()
        result = run_review(
            repo, args.transaction, args.task,
            recover_rejected=args.recover_rejected,
            retry_incomplete=args.retry_incomplete,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (mem.MemoryErrorBase, OSError, UnicodeError, ValueError) as exc:
        print(f"memory review isolation: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
