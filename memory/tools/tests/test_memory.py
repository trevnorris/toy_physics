from __future__ import annotations

import dataclasses
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

import yaml

from memory.tools import memory as mem
from memory.tools import review_isolated
from memory.tools import run_isolated


def git(repo: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", "-C", str(repo), *args], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if proc.returncode:
        raise AssertionError(f"git {' '.join(args)} failed: {proc.stderr}")
    return proc.stdout.strip()


class MemoryRepo:
    def __init__(self, root: Path, units: list[dict] | None = None) -> None:
        self.root = root
        git(root, "init", "-q")
        git(root, "config", "user.email", "tests@example.invalid")
        git(root, "config", "user.name", "Memory Tests")
        git(root, "config", "commit.gpgsign", "false")
        self.write("memory/_meta/schema.md", "# Schema\n")
        self.write("memory/_meta/SCANNER_DECISIONS.md", "# Contract\n")
        self.write("memory/_meta/atlas-migration.yaml", "schema_version: 1\n")
        self.write("memory/prompts/00-snapshot-contract.md", "# Snapshot contract\n")
        self.write("memory/prompts/source-capsule-paper.md", "# Paper task\n")
        self.write("memory/prompts/source-capsule-lifecycle.md", "# Lifecycle task\n")
        self.write(".gitignore", "ignored/\n")
        self.write_config(units or [])

    @staticmethod
    def unit(unit_id: str, source: str, capsule: str | None = None, **member_overrides: object) -> dict:
        member = {
            "path": source,
            "role": "primary",
            "read_mode": "semantic",
            "required": True,
        }
        member.update(member_overrides)
        return {
            "id": unit_id,
            "kind": "paper",
            "capsule": capsule or f"memory/sources/{unit_id}.md",
            "lifecycle": "current",
            "members": [member],
        }

    def config(self, units: list[dict]) -> dict:
        return {
            "schema_version": 2,
            "extractor_contract_version": 1,
            "discovery": {"candidate_roots": ["research", "software"]},
            "read_limits": {
                "semantic_member_max_bytes": 100000,
                "semantic_unit_max_bytes": 200000,
                "excerpt_context_bytes": 256,
            },
            "semantic_excludes": {"prefixes": [], "path_segments": [], "basename_suffixes": [], "suffixes": []},
            "selector_contract": {
                "allowed_keys": [
                    "prefix", "recursive", "suffixes", "basenames", "name_prefixes",
                    "exclude_paths", "exclude_prefixes", "role", "read_mode", "ownership", "required",
                ]
            },
            "units": units,
        }

    def write_config(self, units: list[dict]) -> None:
        self.write("memory/_meta/config.yaml", yaml.safe_dump(self.config(units), sort_keys=False))

    def write(self, relative: str, content: str | bytes) -> Path:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(content, bytes):
            path.write_bytes(content)
        else:
            path.write_text(content, encoding="utf-8")
        return path

    def commit(self, message: str = "snapshot") -> str:
        git(self.root, "add", ".")
        git(self.root, "commit", "-qm", message)
        return git(self.root, "rev-parse", "HEAD")

    def prepare(self, units: list[str] | None = None, paths: list[str] | None = None) -> Path:
        return mem.prepare_update(self.root, "HEAD", units or [], paths or [])

    def staged_page(self, transaction: Path, unit_id: str, *, title: str = "Test", lifecycle: str = "current", bad_anchor: bool = False, bad_link: bool = False) -> Path:
        manifest = json.loads((transaction / "manifest.json").read_text())
        unit = next(item for item in manifest["units"] if item["id"] == unit_id)
        writer = next(item for item in manifest["writer_tasks"] if item["source_unit_id"] == unit_id)
        contract = writer["semantic_contract"]
        sources = sorted(item["path"] for item in unit["members"])
        members = [
            {key: item.get(key) for key in ("path", "role", "read_mode", "mode", "object_type", "blob_oid", "blob_size")}
            for item in unit["members"]
        ]
        source = sources[0] if sources else contract["allowed_citations"][0]["path"]
        anchor = "heading `## Does not exist`" if bad_anchor else "heading `## Claim`"
        body = {
            "schema_version": 2,
            "id": f"source-{unit_id}",
            "title": contract["page"]["title"],
            "type": "source_capsule",
            "lifecycle": lifecycle,
            "memory_review": "ai_draft",
            "sources": sources or [source],
            "content_owner": "ai_generated",
            "last_updated": contract["refresh_date"],
            "generated_from_commit": manifest["target_commit"],
            "source_kind": contract["source_kind"],
            "source_unit": {
                "id": unit_id,
                "shape": contract["shape"],
                "entrypoint": contract["entrypoint"],
                "unit_digest_sha256": unit["unit_digest_sha256"],
                "members": members,
            },
            "extractor_version": contract["extractor_version"],
        }
        markdown = "---\n" + yaml.safe_dump(body, sort_keys=False) + "---\n\n"
        markdown += f"### source-{unit_id}--test-statement — Claim\n\n"
        markdown += "Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`\n\n"
        if sources or contract.get("allowed_citations"):
            markdown += f"Sources:\n\n- `{source}` — {anchor}\n"
        markdown += f"\nFixture variant: {title}\n"
        if bad_link:
            markdown += "\n[missing](../topics/not-there.md)\n"
        target = transaction / "staged" / unit["capsule"]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(markdown, encoding="utf-8")
        packet_seal = json.loads((transaction / "tasks" / unit_id / "packet-seal.json").read_text())
        attestation = {
            "attestation_version": 1,
            "transaction_id": manifest["transaction_id"],
            "task_id": unit_id,
            "source_unit_id": unit_id,
            "isolation": "bubblewrap",
            "workspace_hidden": self.root.as_posix(),
            "packet_path": writer["packet_path"],
            "packet_sha256": packet_seal["combined_sha256"],
            "output_repository_path": writer["output_repository_path"],
            "staged_output_path": writer["staged_output_path"],
            "output_sha256": mem.sha256_file(target),
            "runner_sha256": mem.sha256_file(Path(run_isolated.__file__).resolve()),
            "runtime_profile": "_test",
            "runtime_version": "test-runtime",
            "runtime_executable_sha256": mem.sha256_bytes(b"test-fixture-writer"),
            "completed_at": "2026-08-24T00:00:00+00:00",
        }
        mem.write_json(transaction / "attestations" / f"{unit_id}.json", attestation)
        return target


def add_review_attestation(transaction: Path, task_id: str, verdict: str = "PASS") -> None:
    """Create a hash-complete review fixture without invoking Grok."""
    manifest = json.loads((transaction / "manifest.json").read_text())
    writer = next(item for item in manifest["writer_tasks"] if item["task_id"] == task_id)
    candidate = transaction / writer["staged_output_path"]
    writer_attestation_path = transaction / "attestations" / f"{task_id}.json"
    writer_attestation = json.loads(writer_attestation_path.read_text())
    if writer_attestation.get("runtime_profile") == "_test":
        writer_attestation["runtime_profile"] = "codex"
    mem.write_json(writer_attestation_path, writer_attestation)
    writer_packet = transaction / writer["packet_path"]
    packet_seal = json.loads((writer_packet.parent / "packet-seal.json").read_text())
    review_root = transaction / "reviews" / task_id
    review_packet = review_root / "packet"
    shutil.copytree(writer_packet, review_packet)
    shutil.copy2(candidate, review_packet / "candidate.md")
    mem.atomic_write(review_packet / "review-prompt.md", review_isolated.compose_review_prompt(review_packet), 0o444)
    review_files, review_sha256 = run_isolated.packet_digest(review_packet)
    report = review_root / "output/report.md"
    mem.atomic_write(report, f"Verdict: {verdict}\n".encode(), 0o600)
    review = {
        "attestation_version": 1,
        "role": "independent_review",
        "transaction_id": manifest["transaction_id"],
        "task_id": task_id,
        "packet_sha256": packet_seal["combined_sha256"],
        "review_packet_sha256": review_sha256,
        "review_packet_files": review_files,
        "candidate_path": writer["staged_output_path"],
        "candidate_sha256": mem.sha256_file(candidate),
        "writer_attestation_sha256": mem.sha256_file(writer_attestation_path),
        "report_path": report.relative_to(transaction).as_posix(),
        "report_sha256": mem.sha256_file(report),
        "verdict": verdict,
        "runtime_profile": "grok-review",
        "review_model": review_isolated.GROK_REVIEW_MODEL,
        "review_prompt_sha256": review_isolated.REVIEW_PROMPT_SHA256,
        "review_schema_sha256": review_isolated.REVIEW_SCHEMA_SHA256,
        "review_contract_sha256": review_isolated.REVIEW_CONTRACT_SHA256,
        "runtime_version": "grok-test-reviewer",
        "runtime_executable_sha256": mem.sha256_bytes(b"grok-test-runtime"),
        "reviewer_sha256": mem.sha256_file(Path(review_isolated.__file__).resolve()),
        "completed_at": "2026-08-25T00:00:00+00:00",
        "recovered_from_rejected_output": False,
    }
    mem.write_json(review_root / "attestation.json", review)


class MemoryToolTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.old_test_profile = os.environ.get("MEMORY_TEST_ALLOW_PROFILE")
        os.environ["MEMORY_TEST_ALLOW_PROFILE"] = "1"

    def tearDown(self) -> None:
        if self.old_test_profile is None:
            os.environ.pop("MEMORY_TEST_ALLOW_PROFILE", None)
        else:
            os.environ["MEMORY_TEST_ALLOW_PROFILE"] = self.old_test_profile
        self.temp.cleanup()

    def test_duplicate_yaml_keys_and_config_ownership_overlap_are_rejected(self) -> None:
        with self.assertRaises(mem.ConfigError):
            mem.load_yaml_bytes(b"a: 1\na: 2\n", "duplicate")
        units = [
            MemoryRepo.unit("one", "research/a.md"),
            MemoryRepo.unit("two", "research/a.md"),
        ]
        repo = MemoryRepo(self.root, units)
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        config, _ = mem.load_config(self.root)
        entries = mem.enumerate_tree(self.root, git(self.root, "rev-parse", "HEAD"), ["research", "software"])
        with self.assertRaisesRegex(mem.ConfigError, "primary ownership overlap"):
            mem.resolve_units(config, entries)

        invalid = repo.config([])
        invalid["output_budgets"] = {"topic": {"max_words": 0, "surprise": 1}}
        repo.write("memory/_meta/config.yaml", yaml.safe_dump(invalid, sort_keys=False))
        with self.assertRaisesRegex(mem.ConfigError, "may contain only"):
            mem.load_config(self.root)

    def test_selector_root_nested_filters_and_ignored_volume_independence(self) -> None:
        selector_unit = {
            "id": "selected",
            "kind": "paper",
            "capsule": "memory/sources/selected.md",
            "lifecycle": "current",
            "selectors": [{
                "prefix": "research/bundle", "recursive": False, "suffixes": [".md"],
                "name_prefixes": ["A"], "role": "paper_member", "read_mode": "semantic",
                "ownership": "primary", "required": True,
            }],
        }
        repo = MemoryRepo(self.root, [selector_unit])
        repo.write("research/bundle/A-root.md", "## Claim\nroot\n")
        repo.write("research/bundle/B-root.md", "## Claim\nwrong prefix\n")
        repo.write("research/bundle/nested/A-child.md", "## Claim\nnested\n")
        repo.commit()
        first = mem.status_report(self.root)
        for index in range(300):
            repo.write(f"ignored/noise/{index}.dat", "noise")
        second = mem.status_report(self.root)
        self.assertEqual(first["tracked_candidate_count"], second["tracked_candidate_count"])
        self.assertEqual(first["units"]["selected"]["member_count"], 1)
        config, _ = mem.load_config(self.root)
        entries = mem.enumerate_tree(self.root, first["target_commit"], ["research", "software"])
        self.assertEqual(mem.resolve_units(config, entries)["selected"].members[0]["path"], "research/bundle/A-root.md")

    def test_hard_excluded_roots_reject_unit_members_and_selectors(self) -> None:
        repo = MemoryRepo(self.root)
        invalid_member = repo.config([
            MemoryRepo.unit("ledger", "research/pde_ledger/step.md"),
        ])
        with self.assertRaisesRegex(mem.ConfigError, "hard-excluded source root"):
            mem.validate_config(invalid_member)

        invalid_selector = repo.config([{
            "id": "ledger",
            "kind": "paper",
            "capsule": "memory/sources/ledger.md",
            "lifecycle": "current",
            "selectors": [{
                "prefix": "research/pde_ledger_v3",
                "recursive": True,
                "role": "paper_member",
                "read_mode": "semantic",
                "ownership": "primary",
                "required": True,
            }],
        }])
        with self.assertRaisesRegex(mem.ConfigError, "hard-excluded source root"):
            mem.validate_config(invalid_selector)

    def test_hard_excluded_roots_reject_direct_sources_and_supporting_lineages(self) -> None:
        repo = MemoryRepo(self.root)
        invalid_direct = repo.config([
            MemoryRepo.unit("paper", "research/a.md"),
        ])
        invalid_direct["output_budgets"] = {"topic": {"max_words": 100, "max_key_statements": 2}}
        invalid_direct["derived_pages"] = [{
            "id": "topic-one",
            "task_kind": "topic",
            "page": "memory/topics/one.md",
            "region": "working-position",
            "order": 20,
            "input_units": ["paper"],
            "input_pages": [],
            "direct_sources": [{
                "path": "research/pde_ledger_v2/claim.md",
                "read_mode": "semantic",
            }],
            "budget": "topic",
        }]
        with self.assertRaisesRegex(mem.ConfigError, "hard-excluded source root"):
            mem.validate_config(invalid_direct)

        invalid_lineage = repo.config([])
        invalid_lineage["supporting_lineages"] = [{
            "id": "legacy",
            "root": "research/pde_ledger",
            "selection": "on_demand",
        }]
        with self.assertRaisesRegex(mem.ConfigError, "hard-excluded source root"):
            mem.validate_config(invalid_lineage)

    def test_hard_excluded_roots_reject_migration_original_sources(self) -> None:
        repo = MemoryRepo(self.root)
        repo.write(
            "memory/_meta/atlas-migration.yaml",
            "schema_version: 1\nitems:\n"
            "  - id: forbidden\n"
            "    disposition: migrate\n"
            "    target: memory/topics/one.md\n"
            "    original_sources:\n"
            "      - research/pde_ledger_v3/CHARTER.md#scope\n",
        )
        with self.assertRaisesRegex(mem.ConfigError, "hard-excluded source"):
            mem.migration_requirements(self.root)

    def test_hard_excluded_roots_filter_broad_selectors_and_are_segment_safe(self) -> None:
        broad = {
            "id": "selected",
            "kind": "paper",
            "capsule": "memory/sources/selected.md",
            "lifecycle": "current",
            "selectors": [{
                "prefix": "research",
                "recursive": True,
                "suffixes": [".md"],
                "role": "paper_member",
                "read_mode": "semantic",
                "ownership": "primary",
                "required": True,
            }],
        }
        repo = MemoryRepo(self.root, [broad])
        repo.write("research/ordinary.md", "## Claim\nordinary\n")
        repo.write("research/pde_ledger/hidden.md", "## Claim\nhidden\n")
        repo.write("research/pde_ledger_v2/hidden.md", "## Claim\nhidden\n")
        repo.write("research/pde_ledger_v3/hidden.md", "## Claim\nhidden\n")
        repo.write("research/pde_ledger_v30/allowed.md", "## Claim\nallowed\n")
        commit = repo.commit()
        config, _ = mem.load_config(self.root)
        entries = mem.enumerate_tree(self.root, commit, ["research", "software"])
        paths = [member["path"] for member in mem.resolve_units(config, entries)["selected"].members]
        self.assertEqual(paths, ["research/ordinary.md", "research/pde_ledger_v30/allowed.md"])

    def test_mode_identity_and_identity_only_symlink(self) -> None:
        units = [
            MemoryRepo.unit("script", "research/run.sh"),
            MemoryRepo.unit("link", "research/link", read_mode="identity_only"),
        ]
        repo = MemoryRepo(self.root, units)
        script = repo.write("research/run.sh", "## Claim\necho yes\n")
        script.chmod(0o644)
        os.symlink("run.sh", self.root / "research/link")
        repo.commit("plain")
        first = mem.status_report(self.root)
        first_digest = first["units"]["script"]["unit_digest_sha256"]
        script.chmod(0o755)
        repo.commit("executable")
        second = mem.status_report(self.root)
        self.assertNotEqual(first_digest, second["units"]["script"]["unit_digest_sha256"])
        config, _ = mem.load_config(self.root)
        entries = mem.enumerate_tree(self.root, second["target_commit"], ["research", "software"])
        resolved = mem.resolve_units(config, entries)
        self.assertEqual(resolved["script"].members[0]["mode"], "100755")
        self.assertEqual(resolved["link"].members[0]["mode"], "120000")

    def test_dirty_worktree_never_changes_prepared_blob(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\ncommitted text\n")
        repo.commit()
        repo.write("research/a.md", "## Claim\nDIRTY SECRET\n")
        report = mem.status_report(self.root)
        self.assertEqual(report["dirty_tracked_members"], ["research/a.md"])
        transaction = repo.prepare()
        snapshot = transaction / "snapshot/semantic/research/a.md"
        self.assertIn("committed text", snapshot.read_text())
        self.assertNotIn("DIRTY SECRET", snapshot.read_text())

    def test_prepare_freezes_schema_config_and_prompts(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\ncommitted\n")
        repo.write("memory/prompts/source.md", "Never read live candidates.\n")
        repo.commit()
        transaction = repo.prepare()
        manifest = json.loads((transaction / "manifest.json").read_text())
        frozen = manifest["frozen_policy"]
        self.assertEqual((transaction / frozen["schema"]).read_text(), "# Schema\n")
        self.assertIn("policy/prompts/source.md", frozen["prompts"])
        self.assertEqual(
            (transaction / "policy/prompts/source.md").read_text(),
            "Never read live candidates.\n",
        )
        task = manifest["writer_tasks"][0]
        self.assertEqual(task["source_unit_id"], "paper")
        self.assertEqual(task["staged_output_path"], "staged/memory/sources/paper.md")
        self.assertEqual(task["sandbox_output_path"], "/output/page.md")
        self.assertEqual(task["semantic_inputs"][0]["sandbox_path"], "/packet/inputs/semantic/research/a.md")
        packet = transaction / task["packet_path"]
        self.assertTrue((packet / "schema.md").is_file())
        self.assertTrue((packet / "inputs/semantic/research/a.md").is_file())
        repo.write("memory/_meta/schema.md", "# mutable live schema\n")
        self.assertEqual((transaction / frozen["schema"]).read_text(), "# Schema\n")

    def test_initial_packet_anchor_rejects_user_reseal(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        packet = transaction / "tasks/paper/packet"
        document = json.loads((packet / "task.json").read_text())
        document["untrusted_change"] = True
        mem.atomic_write(packet / "task.json", mem.canonical_json(document), 0o444)
        files, combined = run_isolated.packet_digest(packet)
        seal_path = transaction / "tasks/paper/packet-seal.json"
        seal = json.loads(seal_path.read_text())
        seal.update({"files": files, "combined_sha256": combined})
        mem.write_json(seal_path, seal)
        with self.assertRaisesRegex(mem.MemoryErrorBase, "changed without declared dynamic dependencies"):
            mem.load_transaction(self.root, str(transaction))

    def test_named_runtime_profiles_are_minimal_and_versioned(self) -> None:
        packet = self.root / "packet"
        output = self.root / "output"
        packet.mkdir()
        output.mkdir()
        for name in ("codex",):
            profile = run_isolated.runtime_profile(name)
            self.assertTrue(profile.version)
            self.assertTrue(profile.enforce_resource_limits)
            self.assertRegex(profile.executable_sha256, r"^[0-9a-f]{64}$")
            command = run_isolated.bubblewrap_command(self.root, packet, output, profile)
            self.assertIn("--clearenv", command)
            self.assertFalse(any(
                source == str(Path.home() / ".codex")
                for source, _ in profile.ro_binds
            ))
            self.assertIn(
                "/runtime/codex-code-mode-host",
                [destination for _, destination in profile.ro_binds],
            )
            version_profile = dataclasses.replace(
                profile, command=(f"/runtime/{name}", "--version"), capture_stdout=False, stdin_data=None,
            )
            probe = subprocess.run(
                run_isolated.bubblewrap_command(self.root, packet, output, version_profile),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(probe.returncode, 0, probe.stderr.decode("utf-8", "replace"))
        with self.assertRaisesRegex(run_isolated.IsolationError, "choose codex"):
            run_isolated.runtime_profile("untrusted-writer")
        with self.assertRaisesRegex(SystemExit, "2"):
            run_isolated.build_parser().parse_args([
                "transaction", "--task", "paper", "--profile", "untrusted-writer",
            ])
        self.assertIn("writer task ID", run_isolated.build_parser().format_help())

    def test_owned_process_termination_targets_only_popen_child(self) -> None:
        proc = mock.Mock()
        proc.poll.side_effect = [None, None]
        proc.wait.side_effect = [
            subprocess.TimeoutExpired(["owned-child"], run_isolated.OWNED_PROCESS_TERM_GRACE_SECONDS),
            0,
        ]

        run_isolated._terminate_owned_process(proc)

        proc.terminate.assert_called_once_with()
        proc.kill.assert_called_once_with()
        self.assertEqual(proc.wait.call_count, 2)

    def test_owned_cgroup_configuration_and_stats_are_private(self) -> None:
        parent = self.root / "delegated.slice"
        current = parent / "terminal.scope"
        current.mkdir(parents=True)
        (parent / "cgroup.subtree_control").write_text("memory pids\n", encoding="ascii")
        with (
            mock.patch.object(run_isolated, "_current_cgroup_path", return_value=current),
            mock.patch.object(run_isolated.uuid, "uuid4") as uuid_mock,
        ):
            uuid_mock.return_value.hex = "1234567890abcdef"
            owned = run_isolated._create_owned_cgroup()

        self.assertEqual(owned.parent, parent)
        self.assertIn("codex-memory-writer-", owned.name)
        self.assertEqual(
            (owned / "memory.max").read_text(encoding="ascii"),
            str(run_isolated.CODEX_WRITER_MEMORY_MAX_BYTES),
        )
        self.assertEqual(
            (owned / "memory.swap.max").read_text(encoding="ascii"),
            str(run_isolated.CODEX_WRITER_SWAP_MAX_BYTES),
        )
        self.assertEqual(
            (owned / "pids.max").read_text(encoding="ascii"),
            str(run_isolated.CODEX_WRITER_PIDS_MAX),
        )
        (owned / "memory.peak").write_text("4096\n", encoding="ascii")
        (owned / "memory.events").write_text("oom 0\noom_kill 0\n", encoding="ascii")
        stats = run_isolated._owned_cgroup_stats(owned)
        self.assertEqual(stats["memory_peak"], 4096)
        self.assertEqual(stats["memory_events"], {"oom": 0, "oom_kill": 0})

    def test_owned_process_is_stopped_before_exact_cgroup_attachment(self) -> None:
        profile = run_isolated._test_profile(("/bin/true",))
        profile = dataclasses.replace(profile, enforce_resource_limits=True)
        proc = mock.Mock(pid=424242)
        owned = self.root / "owned-cgroup"
        with (
            mock.patch.object(run_isolated, "_create_owned_cgroup", return_value=owned),
            mock.patch.object(run_isolated, "_attach_stopped_pid") as attach_mock,
            mock.patch.object(run_isolated.subprocess, "Popen", return_value=proc) as popen_mock,
        ):
            actual_proc, actual_owned = run_isolated._spawn_owned_process(
                ["/runtime/codex", "exec"], profile, subprocess.DEVNULL, subprocess.DEVNULL,
                use_cgroup=True,
            )

        self.assertIs(actual_proc, proc)
        self.assertEqual(actual_owned, owned)
        attach_mock.assert_called_once_with(owned, 424242)
        launch = popen_mock.call_args.args[0]
        self.assertEqual(
            launch[:3],
            ["/bin/sh", "-c", 'ulimit -f "$1" || exit 125; shift; kill -STOP "$$"; exec "$@"'],
        )
        self.assertEqual(launch[4], str(run_isolated.CODEX_WRITER_FILE_MAX_BYTES // 512))
        self.assertEqual(launch[-2:], ["/runtime/codex", "exec"])

    def test_unexpected_writer_interruption_terminates_exact_child(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        proc = mock.Mock()
        proc.communicate.side_effect = KeyboardInterrupt()
        with (
            mock.patch.object(run_isolated, "_spawn_owned_process", return_value=(proc, None)),
            mock.patch.object(run_isolated, "_terminate_owned_process") as terminate_mock,
        ):
            with self.assertRaises(KeyboardInterrupt):
                run_isolated.run_task(
                    self.root, str(transaction), "paper", ["/bin/true"],
                )

        terminate_mock.assert_called_once_with(proc)
        self.assertFalse((transaction / "staged/memory/sources/paper.md").exists())
        self.assertFalse((transaction / "attestations/paper.json").exists())

    def test_cgroup_attachment_interruption_terminates_exact_child(self) -> None:
        profile = dataclasses.replace(
            run_isolated._test_profile(("/bin/true",)), enforce_resource_limits=True,
        )
        proc = mock.Mock(pid=424242)
        proc.poll.return_value = None
        owned = self.root / "interrupted-cgroup"
        with (
            mock.patch.object(run_isolated, "_create_owned_cgroup", return_value=owned),
            mock.patch.object(run_isolated, "_attach_stopped_pid", side_effect=KeyboardInterrupt()),
            mock.patch.object(run_isolated.os, "kill") as kill_mock,
            mock.patch.object(run_isolated, "_terminate_owned_process") as terminate_mock,
            mock.patch.object(run_isolated.subprocess, "Popen", return_value=proc),
        ):
            with self.assertRaises(KeyboardInterrupt):
                run_isolated._spawn_owned_process(
                    ["/runtime/codex"], profile, subprocess.DEVNULL, subprocess.DEVNULL,
                    use_cgroup=True,
                )

        kill_mock.assert_called_once_with(424242, signal.SIGCONT)
        terminate_mock.assert_called_once_with(proc)

    def test_writer_output_size_limit_prevents_unbounded_parent_read(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        with mock.patch.object(run_isolated, "CODEX_WRITER_OUTPUT_MAX_BYTES", 5):
            with self.assertRaisesRegex(run_isolated.IsolationError, "exceeds 5 bytes"):
                run_isolated.run_task(
                    self.root, str(transaction), "paper",
                    ["/bin/sh", "-c", "printf 123456 > /output/page.md"],
                )

        self.assertFalse((transaction / "staged/memory/sources/paper.md").exists())
        self.assertFalse((transaction / "attestations/paper.json").exists())

    def test_failed_writer_logs_retain_only_bounded_tails(self) -> None:
        output = self.root / "output"
        output.mkdir()
        stdout_log = self.root / "stdout.bin"
        stderr_log = self.root / "stderr.bin"
        stdout_log.write_bytes(b"0123456789")
        stderr_log.write_bytes(b"abcdefghij")
        with mock.patch.object(run_isolated, "FAILED_LOG_MAX_BYTES", 4):
            run_isolated._preserve_process_logs(output, stdout_log, stderr_log)

        self.assertEqual((output / "failed-stdout.bin").read_bytes(), b"6789")
        self.assertEqual((output / "failed-stderr.bin").read_bytes(), b"ghij")
        self.assertFalse(stdout_log.exists())
        self.assertFalse(stderr_log.exists())

    def test_codex_writer_timeout_preserves_diagnostics_without_import(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        script = (
            "printf 'partial writer stdout'; "
            "printf 'writer stalled' >&2; "
            "sleep 10"
        )
        with mock.patch.object(run_isolated, "CODEX_WRITER_TIMEOUT_SECONDS", 0.2):
            with self.assertRaisesRegex(
                run_isolated.IsolationError, r"timed out after 0.2 seconds: writer stalled",
            ):
                run_isolated.run_task(
                    self.root, str(transaction), "paper", ["/bin/sh", "-c", script],
                )

        output = transaction / "tasks/paper/output"
        self.assertEqual((output / "failed-stdout.bin").read_bytes(), b"partial writer stdout")
        self.assertEqual((output / "failed-stderr.bin").read_bytes(), b"writer stalled")
        self.assertFalse((transaction / "staged/memory/sources/paper.md").exists())
        self.assertFalse((transaction / "attestations/paper.json").exists())

    def test_revision_admits_only_attested_codex_candidate_across_policy_drift(self) -> None:
        unit = MemoryRepo.unit("paper", "research/a.md")
        unit["editorial_focus"] = ["Preserve the established result."]
        repo = MemoryRepo(self.root, [unit])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()

        prior = repo.prepare()
        prior_page = repo.staged_page(prior, "paper", title="prior candidate")
        prior_attestation_path = prior / "attestations/paper.json"
        prior_attestation = json.loads(prior_attestation_path.read_text())
        prior_attestation["runtime_profile"] = "codex"
        mem.write_json(prior_attestation_path, prior_attestation)
        add_review_attestation(prior, "paper", "FAIL")

        revised_unit = MemoryRepo.unit("paper", "research/a.md")
        revised_unit["editorial_focus"] = ["Preserve the result and repair its current policy framing."]
        repo.write_config([revised_unit])  # Policy drift only; HEAD and source blobs remain identical.
        current = repo.prepare()
        prior_manifest = json.loads((prior / "manifest.json").read_text())
        current_manifest = json.loads((current / "manifest.json").read_text())
        self.assertEqual(prior_manifest["target_commit"], current_manifest["target_commit"])
        self.assertNotEqual(
            prior_manifest["units"][0]["unit_digest_sha256"],
            current_manifest["units"][0]["unit_digest_sha256"],
        )

        fake_profile = dataclasses.replace(
            run_isolated._test_profile((
                "/bin/sh", "-c", "cp /packet/revision_candidate.md /output/page.md",
                run_isolated.WRITER_PROMPT,
            )),
            name="codex",
        )
        with mock.patch.object(run_isolated, "runtime_profile", return_value=fake_profile):
            attestation = run_isolated.run_task(
                self.root, str(current), "paper", "codex", revise_from=str(prior),
            )
        packet = current / "tasks/paper/packet"
        revision = json.loads((packet / "revision.json").read_text())
        self.assertEqual((packet / "revision_candidate.md").read_bytes(), prior_page.read_bytes())
        self.assertEqual(revision["prior_runtime_profile"], "codex")
        self.assertEqual(revision["prior_transaction_id"], prior_manifest["transaction_id"])
        self.assertNotEqual(
            revision["prior_prerequisites_sha256"], revision["current_prerequisites_sha256"],
        )
        self.assertEqual(
            revision["prior_review_report_sha256"],
            mem.sha256_file(prior / "reviews/paper/output/report.md"),
        )
        self.assertIn("editorial_focus", run_isolated.REVISION_PROMPT)
        self.assertIn("every failed-review finding", run_isolated.REVISION_PROMPT)
        self.assertIn("Current sealed sources and prerequisite pages", run_isolated.REVISION_PROMPT)
        self.assertEqual(attestation["runtime_profile"], "codex")
        revised_front, _ = mem.parse_frontmatter(
            (current / "staged/memory/sources/paper.md").read_bytes(), "revised fixture",
        )
        self.assertEqual(
            revised_front["source_unit"]["unit_digest_sha256"],
            current_manifest["units"][0]["unit_digest_sha256"],
        )
        self.assertEqual(attestation["normalization"], "source_capsule_frontmatter_v1")
        self.assertEqual(attestation["packet_sha256"], json.loads(
            (current / "tasks/paper/packet-seal.json").read_text()
        )["combined_sha256"])
        self.assertFalse((current / "attestations/paper.json").is_symlink())
        mem.load_transaction(self.root, str(current))

    def test_revision_rejects_non_codex_or_unattested_imports(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        repo.staged_page(prior, "paper")
        attestation_path = prior / "attestations/paper.json"
        attestation = json.loads(attestation_path.read_text())
        attestation["runtime_profile"] = "untrusted-writer"
        mem.write_json(attestation_path, attestation)
        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        writer = next(item for item in current_manifest["writer_tasks"] if item["task_id"] == "paper")
        packet, packet_seal, _ = run_isolated.verify_packet(
            current / "tasks/paper", current_manifest["transaction_id"], "paper",
        )
        with self.assertRaisesRegex(run_isolated.IsolationError, "runtime_profile"):
            run_isolated.admit_revision_candidate(
                self.root, current, current_manifest, writer, packet, packet_seal, str(prior),
            )
        with self.assertRaisesRegex(mem.MemoryErrorBase, "under memory/transactions"):
            run_isolated.admit_revision_candidate(
                self.root, current, current_manifest, writer, packet, packet_seal, "/tmp/arbitrary.md",
            )
        parsed = run_isolated.build_parser().parse_args([
            "current", "--task", "paper", "--profile", "codex", "--revise-from", "prior",
        ])
        self.assertEqual(parsed.revise_from, "prior")

    def test_source_revision_crosses_unrelated_commit_and_finalizes(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        prior_page = repo.staged_page(prior, "paper", title="failed candidate")
        prior_attestation_path = prior / "attestations/paper.json"
        prior_attestation = json.loads(prior_attestation_path.read_text())
        prior_attestation["runtime_profile"] = "codex"
        mem.write_json(prior_attestation_path, prior_attestation)
        add_review_attestation(prior, "paper", "FAIL")

        repo.write("README.md", "unrelated committed change\n")
        git(self.root, "add", "README.md")
        git(self.root, "commit", "-qm", "unrelated change")
        current = repo.prepare()
        prior_manifest = json.loads((prior / "manifest.json").read_text())
        current_manifest = json.loads((current / "manifest.json").read_text())
        self.assertNotEqual(prior_manifest["target_commit"], current_manifest["target_commit"])

        fake_profile = dataclasses.replace(
            run_isolated._test_profile((
                "/bin/sh", "-c", "cp /packet/revision_candidate.md /output/page.md",
                run_isolated.WRITER_PROMPT,
            )), name="codex",
        )
        with mock.patch.object(run_isolated, "runtime_profile", return_value=fake_profile):
            run_isolated.run_task(
                self.root, str(current), "paper", "codex", revise_from=str(prior),
            )
        revision = json.loads((current / "tasks/paper/packet/revision.json").read_text())
        self.assertEqual(revision["prior_target_commit"], prior_manifest["target_commit"])
        self.assertEqual(revision["target_commit"], current_manifest["target_commit"])
        self.assertEqual(
            revision["prior_prerequisites_sha256"], revision["current_prerequisites_sha256"],
        )
        revised = current / "staged/memory/sources/paper.md"
        revised_front, revised_body = mem.parse_frontmatter(revised.read_bytes(), "revised")
        _, prior_body = mem.parse_frontmatter(prior_page.read_bytes(), "prior")
        self.assertEqual(revised_body, prior_body)
        self.assertEqual(revised_front["generated_from_commit"], current_manifest["target_commit"])
        add_review_attestation(current, "paper", "PASS")
        result = mem.finalize_update(self.root, str(current))
        self.assertEqual(result["transaction_id"], current_manifest["transaction_id"])

    def test_cross_commit_source_revision_rejects_changed_member(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        repo.staged_page(prior, "paper")
        prior_attestation_path = prior / "attestations/paper.json"
        prior_attestation = json.loads(prior_attestation_path.read_text())
        prior_attestation["runtime_profile"] = "codex"
        mem.write_json(prior_attestation_path, prior_attestation)
        add_review_attestation(prior, "paper", "FAIL")

        repo.write("research/a.md", "## Claim\nChanged blob\n")
        git(self.root, "add", "research/a.md")
        git(self.root, "commit", "-qm", "change source")
        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        writer = next(item for item in current_manifest["writer_tasks"] if item["task_id"] == "paper")
        packet, packet_seal, _ = run_isolated.verify_packet(
            current / "tasks/paper", current_manifest["transaction_id"], "paper",
        )
        with self.assertRaisesRegex(run_isolated.IsolationError, "selector/member identities"):
            run_isolated.admit_revision_candidate(
                self.root, current, current_manifest, writer, packet, packet_seal, str(prior),
            )

    def test_revision_from_verified_candidate_reuse_and_fail_finalizes(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        original = repo.prepare()
        repo.staged_page(original, "paper", title="original candidate")
        original_attestation_path = original / "attestations/paper.json"
        original_attestation = json.loads(original_attestation_path.read_text())
        original_attestation["runtime_profile"] = "codex"
        mem.write_json(original_attestation_path, original_attestation)

        reused = repo.prepare()
        run_isolated.reuse_candidate_for_recheck(
            self.root, str(reused), "paper", str(original),
        )
        add_review_attestation(reused, "paper", "FAIL")
        current = repo.prepare()
        fake_profile = dataclasses.replace(
            run_isolated._test_profile((
                "/bin/sh", "-c", "cp /packet/revision_candidate.md /output/page.md",
                run_isolated.WRITER_PROMPT,
            )), name="codex",
        )
        with mock.patch.object(run_isolated, "runtime_profile", return_value=fake_profile):
            run_isolated.run_task(
                self.root, str(current), "paper", "codex", revise_from=str(reused),
            )
        revision = json.loads((current / "tasks/paper/packet/revision.json").read_text())
        self.assertEqual(revision["prior_runtime_profile"], "codex-candidate-reuse")
        self.assertEqual(
            revision["prior_attestation_sha256"],
            mem.sha256_file(reused / "attestations/paper.json"),
        )
        add_review_attestation(current, "paper", "PASS")
        result = mem.finalize_update(self.root, str(current))
        current_manifest = json.loads((current / "manifest.json").read_text())
        self.assertEqual(result["transaction_id"], current_manifest["transaction_id"])

    def test_revision_from_candidate_reuse_rejects_tampered_original_chain(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        original = repo.prepare()
        original_candidate = repo.staged_page(original, "paper")
        original_attestation_path = original / "attestations/paper.json"
        original_attestation = json.loads(original_attestation_path.read_text())
        original_attestation["runtime_profile"] = "codex"
        mem.write_json(original_attestation_path, original_attestation)
        reused = repo.prepare()
        run_isolated.reuse_candidate_for_recheck(
            self.root, str(reused), "paper", str(original),
        )
        add_review_attestation(reused, "paper", "FAIL")
        original_candidate.write_text(original_candidate.read_text() + "tampered\n", encoding="utf-8")

        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        writer = next(item for item in current_manifest["writer_tasks"] if item["task_id"] == "paper")
        packet, packet_seal, _ = run_isolated.verify_packet(
            current / "tasks/paper", current_manifest["transaction_id"], "paper",
        )
        with self.assertRaisesRegex(run_isolated.IsolationError, "prior candidate hash mismatch"):
            run_isolated.admit_revision_candidate(
                self.root, current, current_manifest, writer, packet, packet_seal, str(reused),
            )

        missing_root = self.root / "missing-chain"
        missing_root.mkdir()
        missing_repo = MemoryRepo(missing_root, [MemoryRepo.unit("paper", "research/a.md")])
        missing_repo.write("research/a.md", "## Claim\nA\n")
        missing_repo.commit()
        missing_original = missing_repo.prepare()
        missing_repo.staged_page(missing_original, "paper")
        missing_original_attestation_path = missing_original / "attestations/paper.json"
        missing_original_attestation = json.loads(missing_original_attestation_path.read_text())
        missing_original_attestation["runtime_profile"] = "codex"
        mem.write_json(missing_original_attestation_path, missing_original_attestation)
        missing_reused = missing_repo.prepare()
        run_isolated.reuse_candidate_for_recheck(
            missing_root, str(missing_reused), "paper", str(missing_original),
        )
        add_review_attestation(missing_reused, "paper", "FAIL")
        missing_attestation_path = missing_reused / "attestations/paper.json"
        missing_attestation = json.loads(missing_attestation_path.read_text())
        del missing_attestation["candidate_reuse"]
        mem.write_json(missing_attestation_path, missing_attestation)
        missing_current = missing_repo.prepare()
        missing_manifest = json.loads((missing_current / "manifest.json").read_text())
        missing_writer = next(
            item for item in missing_manifest["writer_tasks"] if item["task_id"] == "paper"
        )
        missing_packet, missing_seal, _ = run_isolated.verify_packet(
            missing_current / "tasks/paper", missing_manifest["transaction_id"], "paper",
        )
        with self.assertRaisesRegex(run_isolated.IsolationError, "lacks its provenance chain"):
            run_isolated.admit_revision_candidate(
                missing_root, missing_current, missing_manifest, missing_writer,
                missing_packet, missing_seal, str(missing_reused),
            )

    def test_revision_rejects_reviewed_reuse_candidate(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        original = repo.prepare()
        repo.staged_page(original, "paper")
        original_attestation_path = original / "attestations/paper.json"
        original_attestation = json.loads(original_attestation_path.read_text())
        original_attestation["runtime_profile"] = "codex"
        mem.write_json(original_attestation_path, original_attestation)
        add_review_attestation(original, "paper", "PASS")
        reviewed = repo.prepare()
        run_isolated.reuse_reviewed_candidate(
            self.root, str(reviewed), "paper", str(original),
        )

        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        writer = next(item for item in current_manifest["writer_tasks"] if item["task_id"] == "paper")
        packet, packet_seal, _ = run_isolated.verify_packet(
            current / "tasks/paper", current_manifest["transaction_id"], "paper",
        )
        with self.assertRaisesRegex(run_isolated.IsolationError, "not an eligible Codex candidate"):
            run_isolated.admit_revision_candidate(
                self.root, current, current_manifest, writer, packet, packet_seal, str(reviewed),
            )

    def test_recorded_lint_failure_revision_repairs_and_finalizes(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        failed = repo.prepare()
        repo.staged_page(failed, "paper", bad_anchor=True)
        failed_attestation_path = failed / "attestations/paper.json"
        failed_attestation = json.loads(failed_attestation_path.read_text())
        failed_attestation["runtime_profile"] = "codex"
        mem.write_json(failed_attestation_path, failed_attestation)
        add_review_attestation(failed, "paper", "PASS")
        record = mem.record_staged_lint_failure(self.root, str(failed))
        self.assertEqual(record["role"], "machine_lint_failure")
        self.assertEqual(record["verdict"], "FAIL")
        self.assertTrue(all(
            error.startswith("memory/sources/paper.md:") for error in record["errors"]
        ))
        self.assertEqual(
            record["candidate_sha256"], mem.sha256_file(failed / "staged/memory/sources/paper.md"),
        )

        current = repo.prepare()
        fake_profile = dataclasses.replace(
            run_isolated._test_profile((
                "/bin/sh", "-c",
                "sed 's/Does not exist/Claim/' /packet/revision_candidate.md > /output/page.md",
                run_isolated.WRITER_PROMPT,
            )), name="codex",
        )
        with mock.patch.object(run_isolated, "runtime_profile", return_value=fake_profile):
            run_isolated.run_task(
                self.root, str(current), "paper", "codex", revise_lint_from=str(failed),
            )
        revision = json.loads((current / "tasks/paper/packet/revision.json").read_text())
        self.assertEqual(revision["revision_basis"], "machine_lint_failure")
        self.assertEqual(
            revision["prior_lint_record_sha256"],
            mem.sha256_file(failed / "lint-failures/paper/record.json"),
        )
        self.assertEqual(json.loads(
            (current / "tasks/paper/packet/revision_review_attestation.json").read_text()
        )["verdict"], "PASS")
        add_review_attestation(current, "paper", "PASS")
        errors, _ = mem.lint_staged(self.root, str(current))
        self.assertEqual(errors, [])
        result = mem.finalize_update(self.root, str(current))
        self.assertEqual(
            result["transaction_id"], json.loads((current / "manifest.json").read_text())["transaction_id"],
        )

    def test_recorded_lint_failure_refuses_pass_missing_review_and_unscoped_errors(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        passed = repo.prepare()
        repo.staged_page(passed, "paper")
        passed_attestation_path = passed / "attestations/paper.json"
        passed_attestation = json.loads(passed_attestation_path.read_text())
        passed_attestation["runtime_profile"] = "codex"
        mem.write_json(passed_attestation_path, passed_attestation)
        add_review_attestation(passed, "paper", "PASS")
        with self.assertRaisesRegex(mem.MemoryErrorBase, "staged lint passed"):
            mem.record_staged_lint_failure(self.root, str(passed))

        missing_review = repo.prepare()
        repo.staged_page(missing_review, "paper", bad_anchor=True)
        missing_attestation_path = missing_review / "attestations/paper.json"
        missing_attestation = json.loads(missing_attestation_path.read_text())
        missing_attestation["runtime_profile"] = "codex"
        mem.write_json(missing_attestation_path, missing_attestation)
        with self.assertRaisesRegex(mem.MemoryErrorBase, "current Grok PASS"):
            mem.record_staged_lint_failure(self.root, str(missing_review))

        # A served-state error has no unique staged-task prefix and must never be delegated.
        state_path = self.root / "memory/_meta/state.json"
        state, _ = mem.load_state(self.root)
        state["pages"]["memory/sources/missing.md"] = {
            "sha256": "0" * 64, "generation": "fixture",
        }
        mem.write_json(state_path, state)
        ambiguous = repo.prepare()
        repo.staged_page(ambiguous, "paper", bad_anchor=True, bad_link=True)
        ambiguous_attestation_path = ambiguous / "attestations/paper.json"
        ambiguous_attestation = json.loads(ambiguous_attestation_path.read_text())
        ambiguous_attestation["runtime_profile"] = "codex"
        mem.write_json(ambiguous_attestation_path, ambiguous_attestation)
        add_review_attestation(ambiguous, "paper", "PASS")
        with self.assertRaisesRegex(mem.MemoryErrorBase, "ambiguous or unscoped"):
            mem.record_staged_lint_failure(self.root, str(ambiguous))

    def test_lint_revision_rejects_tampered_or_arbitrary_record(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        failed = repo.prepare()
        repo.staged_page(failed, "paper", bad_anchor=True)
        attestation_path = failed / "attestations/paper.json"
        attestation = json.loads(attestation_path.read_text())
        attestation["runtime_profile"] = "codex"
        mem.write_json(attestation_path, attestation)
        add_review_attestation(failed, "paper", "PASS")
        mem.record_staged_lint_failure(self.root, str(failed))
        report = failed / "lint-failures/paper/report.md"
        report.chmod(0o644)
        report.write_text(report.read_text() + "tampered\n", encoding="utf-8")

        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        writer = next(item for item in current_manifest["writer_tasks"] if item["task_id"] == "paper")
        packet, packet_seal, _ = run_isolated.verify_packet(
            current / "tasks/paper", current_manifest["transaction_id"], "paper",
        )
        with self.assertRaisesRegex(run_isolated.IsolationError, "report"):
            run_isolated.admit_lint_revision_candidate(
                self.root, current, current_manifest, writer, packet, packet_seal, str(failed),
            )

        parsed = run_isolated.build_parser().parse_args([
            "current", "--task", "paper", "--profile", "codex",
            "--revise-lint-from", "failed",
        ])
        self.assertEqual(parsed.revise_lint_from, "failed")
        with self.assertRaisesRegex(SystemExit, "2"):
            run_isolated.build_parser().parse_args([
                "current", "--task", "paper", "--profile", "codex",
                "--revise-from", "failed", "--revise-lint-from", "failed",
            ])

    def test_revision_supports_exact_derived_selector_and_prerequisite_hashes(self) -> None:
        unit = MemoryRepo.unit("paper", "research/a.md")
        repo = MemoryRepo(self.root, [unit])
        repo.write("memory/prompts/topic-synthesis.md", "# Topic\n")
        config = repo.config([unit])
        config["output_budgets"] = {
            "source_capsule": {"max_words": 100}, "topic": {"max_words": 100},
        }
        config["derived_pages"] = [{
            "id": "topic-one", "task_kind": "topic", "page": "memory/topics/one.md",
            "region": "working-position", "order": 20, "input_units": ["paper"],
            "input_pages": [], "budget": "topic",
            "direct_sources": [{"path": "research/a.md", "read_mode": "semantic"}],
        }]
        repo.write("memory/_meta/config.yaml", yaml.safe_dump(config, sort_keys=False))
        repo.write("research/a.md", "## Claim\nA\n")
        repo.write(
            "memory/topics/one.md",
            "---\nschema_version: 2\nid: topic-one\ntitle: Topic\ntype: topic\n"
            "lifecycle: current\nmemory_review: ai_draft\nsources: []\ncontent_owner: mixed\n"
            "last_updated: 2026-08-25\ngenerated_from_commit: pending\n---\n\n"
            "<!-- BEGIN GENERATED:working-position -->\npending\n"
            "<!-- END GENERATED:working-position -->\n",
        )
        repo.commit()

        prior = repo.prepare()
        repo.staged_page(prior, "paper")
        prior_attestation = run_isolated.run_task(
            self.root, str(prior), "topic-one",
            ["/bin/sh", "-c", "cp /packet/base_page.md /output/page.md"],
        )
        prior_attestation["runtime_profile"] = "codex"
        mem.write_json(prior / "attestations/topic-one.json", prior_attestation)
        add_review_attestation(prior, "topic-one", "FAIL")

        current = repo.prepare()
        repo.staged_page(current, "paper")
        fake_profile = dataclasses.replace(
            run_isolated._test_profile((
                "/bin/sh", "-c", "cp /packet/revision_candidate.md /output/page.md",
                run_isolated.WRITER_PROMPT,
            )), name="codex",
        )
        with mock.patch.object(run_isolated, "runtime_profile", return_value=fake_profile):
            run_isolated.run_task(
                self.root, str(current), "topic-one", "codex", revise_from=str(prior),
            )
        packet = current / "tasks/topic-one/packet"
        revision = json.loads((packet / "revision.json").read_text())
        self.assertIsNone(revision["source_unit_id"])
        self.assertEqual(revision["task_kind"], "topic")
        self.assertRegex(revision["current_prerequisites_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(revision["retry_selector_sha256"], revision["prior_retry_selector_sha256"])
        self.assertTrue((packet / "revision_review.md").is_file())
        self.assertTrue((packet / "revision_review_attestation.json").is_file())
        mem.load_transaction(self.root, str(current))

        repo.write("README.md", "unrelated commit after derived failure\n")
        git(self.root, "add", "README.md")
        git(self.root, "commit", "-qm", "unrelated derived commit")
        cross_commit = repo.prepare()
        cross_manifest = json.loads((cross_commit / "manifest.json").read_text())
        cross_writer = next(
            item for item in cross_manifest["writer_tasks"] if item["task_id"] == "topic-one"
        )
        cross_packet, cross_seal, _ = run_isolated.verify_packet(
            cross_commit / "tasks/topic-one", cross_manifest["transaction_id"], "topic-one",
        )
        with self.assertRaisesRegex(run_isolated.IsolationError, "derived revision target commit"):
            run_isolated.admit_revision_candidate(
                self.root, cross_commit, cross_manifest, cross_writer,
                cross_packet, cross_seal, str(prior),
            )

    def test_reviewed_reuse_normalizes_and_finalizes_without_model_runtime(self) -> None:
        unit = MemoryRepo.unit("paper", "research/a.md")
        unit["editorial_focus"] = ["Original policy focus."]
        repo = MemoryRepo(self.root, [unit])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        prior_page = repo.staged_page(prior, "paper", title="reviewed body")
        add_review_attestation(prior, "paper", "PASS")

        current_unit = MemoryRepo.unit("paper", "research/a.md")
        current_unit["editorial_focus"] = ["Current tool-only policy focus."]
        repo.write_config([current_unit])
        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        attestation = run_isolated.reuse_reviewed_candidate(
            self.root, str(current), "paper", str(prior),
        )
        self.assertEqual(attestation["runtime_profile"], "codex-reviewed-reuse")
        self.assertEqual(attestation["isolation"], "deterministic-reviewed-reuse")
        self.assertIs(attestation["model_invoked"], False)
        self.assertIsNone(attestation["runtime_version"])
        self.assertIsNone(attestation["runtime_executable_sha256"])
        reused = current / "staged/memory/sources/paper.md"
        prior_front, prior_body = mem.parse_frontmatter(prior_page.read_bytes(), "prior")
        reused_front, reused_body = mem.parse_frontmatter(reused.read_bytes(), "reused")
        self.assertEqual(reused_body, prior_body)
        self.assertNotEqual(
            prior_front["source_unit"]["unit_digest_sha256"],
            reused_front["source_unit"]["unit_digest_sha256"],
        )
        self.assertEqual(
            reused_front["source_unit"]["unit_digest_sha256"],
            current_manifest["units"][0]["unit_digest_sha256"],
        )
        staged_hashes = {"memory/sources/paper.md": mem.sha256_file(reused)}
        self.assertEqual(
            mem.verify_isolation_attestations(self.root, current, current_manifest, staged_hashes), [],
        )
        result = mem.finalize_update(self.root, str(current))
        self.assertEqual(result["transaction_id"], current_manifest["transaction_id"])

    def test_reviewed_reuse_rejects_fail_and_finalizer_rechecks_tampering(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        failed = repo.prepare()
        repo.staged_page(failed, "paper")
        add_review_attestation(failed, "paper", "FAIL")
        current = repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "verdict"):
            run_isolated.reuse_reviewed_candidate(self.root, str(current), "paper", str(failed))

        passed = repo.prepare()
        repo.staged_page(passed, "paper", title="passed")
        add_review_attestation(passed, "paper", "PASS")
        reusable = repo.prepare()
        run_isolated.reuse_reviewed_candidate(self.root, str(reusable), "paper", str(passed))
        report = passed / "reviews/paper/output/report.md"
        report.write_text("tampered after reuse\n", encoding="utf-8")
        with self.assertRaisesRegex(mem.LintFailure, "Grok report hash mismatch"):
            mem.finalize_update(self.root, str(reusable))

    def test_reviewed_reuse_rejects_candidate_tampering_and_profile_argument(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        candidate = repo.staged_page(prior, "paper")
        add_review_attestation(prior, "paper", "PASS")
        candidate.write_text(candidate.read_text() + "tampered\n", encoding="utf-8")
        current = repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "candidate changed"):
            run_isolated.reuse_reviewed_candidate(self.root, str(current), "paper", str(prior))
        parsed = run_isolated.build_parser().parse_args([
            "current", "--task", "paper", "--reuse-reviewed-from", "prior",
        ])
        self.assertEqual(parsed.reuse_reviewed_from, "prior")
        with self.assertRaisesRegex(SystemExit, "2"):
            run_isolated.main([
                "current", "--task", "paper", "--profile", "codex",
                "--reuse-reviewed-from", "prior",
            ])

    def test_reviewed_reuse_crosses_unrelated_commit_and_finalizes(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        prior_page = repo.staged_page(prior, "paper", title="reviewed body")
        add_review_attestation(prior, "paper", "PASS")
        prior_manifest = json.loads((prior / "manifest.json").read_text())

        repo.write("README.md", "unrelated committed change\n")
        git(self.root, "add", "README.md")
        git(self.root, "commit", "-qm", "unrelated change")
        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        self.assertNotEqual(prior_manifest["target_commit"], current_manifest["target_commit"])
        attestation = run_isolated.reuse_reviewed_candidate(
            self.root, str(current), "paper", str(prior),
        )
        reused = current / "staged/memory/sources/paper.md"
        reused_front, reused_body = mem.parse_frontmatter(reused.read_bytes(), "reused")
        _, prior_body = mem.parse_frontmatter(prior_page.read_bytes(), "prior")
        self.assertEqual(reused_body, prior_body)
        self.assertEqual(reused_front["generated_from_commit"], current_manifest["target_commit"])
        self.assertEqual(
            attestation["reviewed_reuse"]["prior_target_commit"], prior_manifest["target_commit"],
        )
        self.assertEqual(
            attestation["reviewed_reuse"]["current_target_commit"], current_manifest["target_commit"],
        )
        self.assertFalse((current / "reviews/paper").exists())
        result = mem.finalize_update(self.root, str(current))
        self.assertEqual(result["transaction_id"], current_manifest["transaction_id"])

    def test_reviewed_reuse_cross_commit_rejects_changed_source(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        repo.staged_page(prior, "paper")
        add_review_attestation(prior, "paper", "PASS")

        repo.write("research/a.md", "## Claim\nChanged source blob\n")
        git(self.root, "add", "research/a.md")
        git(self.root, "commit", "-qm", "change source blob")
        current = repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "member identities"):
            run_isolated.reuse_reviewed_candidate(
                self.root, str(current), "paper", str(prior),
            )

    def test_candidate_reuse_requires_fresh_pass_then_finalizes(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        prior_page = repo.staged_page(prior, "paper", title="candidate to recheck")
        prior_attestation_path = prior / "attestations/paper.json"
        prior_attestation = json.loads(prior_attestation_path.read_text())
        prior_attestation["runtime_profile"] = "codex"
        mem.write_json(prior_attestation_path, prior_attestation)

        repo.write("README.md", "unrelated committed change\n")
        git(self.root, "add", "README.md")
        git(self.root, "commit", "-qm", "unrelated change")
        current = repo.prepare()
        current_manifest = json.loads((current / "manifest.json").read_text())
        prior_manifest = json.loads((prior / "manifest.json").read_text())
        self.assertNotEqual(prior_manifest["target_commit"], current_manifest["target_commit"])
        attestation = run_isolated.reuse_candidate_for_recheck(
            self.root, str(current), "paper", str(prior),
        )
        self.assertEqual(attestation["runtime_profile"], "codex-candidate-reuse")
        self.assertEqual(attestation["isolation"], "deterministic-candidate-reuse")
        self.assertIs(attestation["model_invoked"], False)
        self.assertIsNone(attestation["runtime_version"])
        self.assertIsNone(attestation["runtime_executable_sha256"])
        reused = current / "staged/memory/sources/paper.md"
        _, prior_body = mem.parse_frontmatter(prior_page.read_bytes(), "prior")
        reused_front, reused_body = mem.parse_frontmatter(reused.read_bytes(), "reused")
        self.assertEqual(reused_body, prior_body)
        self.assertEqual(reused_front["generated_from_commit"], current_manifest["target_commit"])
        self.assertEqual(
            attestation["candidate_reuse"]["prior_target_commit"], prior_manifest["target_commit"],
        )
        self.assertEqual(
            attestation["candidate_reuse"]["current_target_commit"], current_manifest["target_commit"],
        )
        self.assertRegex(attestation["candidate_reuse"]["prior_manifest_sha256"], r"^[0-9a-f]{64}$")
        with self.assertRaisesRegex(mem.LintFailure, "missing required Grok PASS review"):
            mem.finalize_update(self.root, str(current))
        review_profile = run_isolated.RuntimeProfile(
            name="grok-review", version="test-grok",
            executable_sha256=mem.sha256_bytes(b"test-grok-runtime"),
            command=("/runtime/grok",), ro_binds=(), environment=(), capture_stdout=True,
        )
        review_process = subprocess.CompletedProcess(
            ["fake-grok"], 0,
            stdout=json.dumps({
                "structuredOutput": {"verdict": "PASS", "findings": []},
            }).encode(),
            stderr=b"",
        )
        with (
            mock.patch.object(review_isolated, "grok_profile", return_value=review_profile),
            mock.patch.object(run_isolated, "bubblewrap_command", return_value=["fake-grok"]),
            mock.patch.object(review_isolated.subprocess, "run", return_value=review_process),
        ):
            review = review_isolated.run_review(self.root, str(current), "paper")
        self.assertEqual(review["verdict"], "PASS")
        result = mem.finalize_update(self.root, str(current))
        self.assertEqual(result["transaction_id"], current_manifest["transaction_id"])

    def test_candidate_reuse_rejects_tampered_prior_chain(self) -> None:
        for tamper in ("packet", "candidate", "writer"):
            with self.subTest(tamper=tamper):
                case_root = self.root / tamper
                case_root.mkdir()
                repo = MemoryRepo(case_root, [MemoryRepo.unit("paper", "research/a.md")])
                repo.write("research/a.md", "## Claim\nA\n")
                repo.commit()
                prior = repo.prepare()
                candidate = repo.staged_page(prior, "paper")
                attestation_path = prior / "attestations/paper.json"
                writer_attestation = json.loads(attestation_path.read_text())
                writer_attestation["runtime_profile"] = "codex"
                mem.write_json(attestation_path, writer_attestation)
                if tamper == "packet":
                    packet_file = prior / "tasks/paper/packet/schema.md"
                    packet_file.chmod(0o644)
                    packet_file.write_text("tampered\n", encoding="utf-8")
                elif tamper == "candidate":
                    candidate.write_text(candidate.read_text() + "tampered\n", encoding="utf-8")
                else:
                    writer_attestation["output_sha256"] = mem.sha256_bytes(b"forged")
                    mem.write_json(attestation_path, writer_attestation)
                current = repo.prepare()
                with self.assertRaises((run_isolated.IsolationError, mem.MemoryErrorBase)):
                    run_isolated.reuse_candidate_for_recheck(
                        case_root, str(current), "paper", str(prior),
                    )

    def test_candidate_reuse_rejects_source_identity_or_prerequisite_drift(self) -> None:
        unit = MemoryRepo.unit("paper", "research/a.md")
        repo = MemoryRepo(self.root, [unit])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.write("research/b.md", "## Claim\nB\n")
        repo.commit()
        prior = repo.prepare()
        repo.staged_page(prior, "paper")
        attestation_path = prior / "attestations/paper.json"
        writer_attestation = json.loads(attestation_path.read_text())
        writer_attestation["runtime_profile"] = "codex"
        mem.write_json(attestation_path, writer_attestation)

        repo.write_config([MemoryRepo.unit("paper", "research/b.md")])
        changed_members = repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "member identities"):
            run_isolated.reuse_candidate_for_recheck(
                self.root, str(changed_members), "paper", str(prior),
            )

        changed_focus = MemoryRepo.unit("paper", "research/a.md")
        changed_focus["editorial_focus"] = ["A changed source-unit prerequisite."]
        repo.write_config([changed_focus])
        changed_prerequisites = repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "prerequisites changed"):
            run_isolated.reuse_candidate_for_recheck(
                self.root, str(changed_prerequisites), "paper", str(prior),
            )

    def test_candidate_reuse_rejects_changed_blob_across_commit(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        prior = repo.prepare()
        repo.staged_page(prior, "paper")
        attestation_path = prior / "attestations/paper.json"
        writer_attestation = json.loads(attestation_path.read_text())
        writer_attestation["runtime_profile"] = "codex"
        mem.write_json(attestation_path, writer_attestation)

        repo.write("research/a.md", "## Claim\nChanged source blob\n")
        git(self.root, "add", "research/a.md")
        git(self.root, "commit", "-qm", "change source blob")
        current = repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "member identities"):
            run_isolated.reuse_candidate_for_recheck(
                self.root, str(current), "paper", str(prior),
            )

    def test_candidate_reuse_rejects_misuse_and_non_original_candidate(self) -> None:
        unit = MemoryRepo.unit("paper", "research/a.md")
        repo = MemoryRepo(self.root, [unit])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        original = repo.prepare()
        repo.staged_page(original, "paper")
        original_attestation_path = original / "attestations/paper.json"
        original_attestation = json.loads(original_attestation_path.read_text())
        original_attestation["runtime_profile"] = "codex"
        mem.write_json(original_attestation_path, original_attestation)
        add_review_attestation(original, "paper", "PASS")
        reviewed_reuse = repo.prepare()
        run_isolated.reuse_reviewed_candidate(self.root, str(reviewed_reuse), "paper", str(original))
        current = repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "attestation mismatch"):
            run_isolated.reuse_candidate_for_recheck(
                self.root, str(current), "paper", str(reviewed_reuse),
            )

        derived_root = self.root / "derived"
        derived_root.mkdir()
        derived_repo = MemoryRepo(derived_root, [unit])
        derived_repo.write("memory/prompts/topic-synthesis.md", "# Topic\n")
        derived_config = derived_repo.config([unit])
        derived_config["output_budgets"] = {
            "source_capsule": {"max_words": 100}, "topic": {"max_words": 100},
        }
        derived_config["derived_pages"] = [{
            "id": "topic-one", "task_kind": "topic", "page": "memory/topics/one.md",
            "region": "working-position", "order": 20, "input_units": ["paper"],
            "input_pages": [], "budget": "topic",
            "direct_sources": [{"path": "research/a.md", "read_mode": "semantic"}],
        }]
        derived_repo.write("memory/_meta/config.yaml", yaml.safe_dump(derived_config, sort_keys=False))
        derived_repo.write("research/a.md", "## Claim\nA\n")
        derived_repo.write(
            "memory/topics/one.md",
            "---\nschema_version: 2\nid: topic-one\ntitle: Topic\ntype: topic\n"
            "lifecycle: current\nmemory_review: ai_draft\nsources: []\ncontent_owner: mixed\n"
            "last_updated: 2026-08-25\ngenerated_from_commit: pending\n---\n\n"
            "<!-- BEGIN GENERATED:working-position -->\npending\n"
            "<!-- END GENERATED:working-position -->\n",
        )
        derived_repo.commit()
        derived_transaction = derived_repo.prepare()
        with self.assertRaisesRegex(run_isolated.IsolationError, "source-unit tasks"):
            run_isolated.reuse_candidate_for_recheck(
                derived_root, str(derived_transaction), "topic-one", str(derived_transaction),
            )

        parsed = run_isolated.build_parser().parse_args([
            "current", "--task", "paper", "--reuse-candidate-from", "prior",
        ])
        self.assertEqual(parsed.reuse_candidate_from, "prior")
        with self.assertRaisesRegex(SystemExit, "2"):
            run_isolated.main([
                "current", "--task", "paper", "--profile", "codex",
                "--reuse-candidate-from", "prior",
            ])
        with self.assertRaisesRegex(SystemExit, "2"):
            run_isolated.build_parser().parse_args([
                "current", "--task", "paper", "--revise-from", "prior",
                "--reuse-candidate-from", "other",
            ])

    def test_grok_profile_is_review_only_and_minimal(self) -> None:
        packet = self.root / "packet"
        output = self.root / "output"
        packet.mkdir()
        output.mkdir()
        profile = review_isolated.grok_profile()
        self.assertEqual(profile.name, "grok-review")
        self.assertIn("json", profile.command)
        self.assertIn("--json-schema", profile.command)
        self.assertIn("--prompt-file", profile.command)
        self.assertIn("/packet/review-prompt.md", profile.command)
        self.assertEqual(profile.command[profile.command.index("--model") + 1], review_isolated.GROK_REVIEW_MODEL)
        self.assertNotIn("--tools", profile.command)
        self.assertIn("--no-subagents", profile.command)
        self.assertIn("--disable-web-search", profile.command)
        self.assertNotIn("grok", run_isolated.build_parser().format_help())
        self.assertFalse(any(source == str(Path.home() / ".grok") for source, _ in profile.ro_binds))
        command = run_isolated.bubblewrap_command(self.root, packet, output, profile)
        self.assertIn("--clearenv", command)
        version_profile = dataclasses.replace(
            profile, command=("/runtime/grok", "--version"), capture_stdout=False,
        )
        probe = subprocess.run(
            run_isolated.bubblewrap_command(self.root, packet, output, version_profile),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        self.assertEqual(probe.returncode, 0, probe.stderr.decode("utf-8", "replace"))
        parsed = review_isolated._parse_structured_review(
            '{"findings":[],"verdict":"FAIL"}'
            '{"findings":[{"severity":"major"}],"verdict":"FAIL"}'
        )
        self.assertEqual(parsed["verdict"], "FAIL")
        self.assertEqual(len(parsed["findings"]), 1)

    def test_incomplete_grok_network_failure_retries_from_lossless_archive(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        candidate = repo.staged_page(transaction, "paper")
        writer_attestation_path = transaction / "attestations/paper.json"
        writer_attestation = json.loads(writer_attestation_path.read_text())
        writer_attestation["runtime_profile"] = "codex"
        mem.write_json(writer_attestation_path, writer_attestation)
        profile = run_isolated.RuntimeProfile(
            name="grok-review",
            version="test-grok",
            executable_sha256=mem.sha256_bytes(b"test-grok-runtime"),
            command=("/runtime/grok",),
            ro_binds=(),
            environment=(),
            capture_stdout=True,
        )
        network_failure = subprocess.CompletedProcess(
            ["fake-grok"], 75, stdout=b"partial reviewer output", stderr=b"network unavailable",
        )
        success = subprocess.CompletedProcess(
            ["fake-grok"], 0,
            stdout=json.dumps({
                "structuredOutput": {"verdict": "PASS", "findings": []},
            }).encode(),
            stderr=b"",
        )

        def inventory(root: Path) -> dict[str, tuple[str, int, bytes | None]]:
            result: dict[str, tuple[str, int, bytes | None]] = {}
            for path in sorted(root.rglob("*")):
                relative = path.relative_to(root).as_posix()
                if path.is_dir():
                    result[relative] = ("directory", path.stat().st_mode & 0o777, None)
                else:
                    result[relative] = ("file", path.stat().st_mode & 0o777, path.read_bytes())
            return result

        with (
            mock.patch.object(review_isolated, "grok_profile", return_value=profile),
            mock.patch.object(run_isolated, "bubblewrap_command", return_value=["fake-grok"]),
            mock.patch.object(review_isolated.subprocess, "run", side_effect=[network_failure, success]),
        ):
            with self.assertRaisesRegex(review_isolated.ReviewError, "network unavailable"):
                review_isolated.run_review(self.root, str(transaction), "paper")
            failed_root = transaction / "reviews/paper"
            before_archive = inventory(failed_root)
            self.assertEqual(
                (failed_root / "output/failed-stdout.bin").read_bytes(), b"partial reviewer output",
            )
            self.assertEqual(
                (failed_root / "output/failed-stderr.bin").read_bytes(), b"network unavailable",
            )
            review = review_isolated.run_review(
                self.root, str(transaction), "paper", retry_incomplete=True,
            )

        archive = transaction / "rejected-reviews/paper-attempt-0001"
        self.assertEqual(inventory(archive), before_archive)
        self.assertEqual(
            review["retried_from_incomplete_archive"],
            "rejected-reviews/paper-attempt-0001",
        )
        self.assertEqual(review["verdict"], "PASS")
        self.assertEqual(
            {path.name for path in (transaction / "reviews/paper/output").iterdir()},
            {"report.md"},
        )
        manifest = json.loads((transaction / "manifest.json").read_text())
        errors = mem.verify_isolation_attestations(
            self.root,
            transaction,
            manifest,
            {"memory/sources/paper.md": mem.sha256_file(candidate)},
        )
        self.assertEqual(errors, [], "archived incomplete reviews must be outside the active review set")
        (transaction / "reviews/unexpected").mkdir()
        errors = mem.verify_isolation_attestations(
            self.root,
            transaction,
            manifest,
            {"memory/sources/paper.md": mem.sha256_file(candidate)},
        )
        self.assertIn("unexpected current review artifact: unexpected", errors)

    def test_grok_timeout_is_lossless_and_retryable_as_incomplete(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        repo.staged_page(transaction, "paper")
        writer_attestation_path = transaction / "attestations/paper.json"
        writer_attestation = json.loads(writer_attestation_path.read_text())
        writer_attestation["runtime_profile"] = "codex"
        mem.write_json(writer_attestation_path, writer_attestation)
        profile = run_isolated.RuntimeProfile(
            name="grok-review", version="test-grok",
            executable_sha256=mem.sha256_bytes(b"test-grok-runtime"),
            command=("/runtime/grok",), ro_binds=(), environment=(), capture_stdout=True,
        )
        expired = subprocess.TimeoutExpired(
            ["fake-grok"], review_isolated.GROK_REVIEW_TIMEOUT_SECONDS,
            output="partial Grok stdout", stderr=b"review stalled",
        )
        success = subprocess.CompletedProcess(
            ["fake-grok"], 0,
            stdout=json.dumps({
                "structuredOutput": {"verdict": "PASS", "findings": []},
            }).encode(),
            stderr=b"",
        )
        with (
            mock.patch.object(review_isolated, "grok_profile", return_value=profile),
            mock.patch.object(run_isolated, "bubblewrap_command", return_value=["fake-grok"]),
            mock.patch.object(
                review_isolated.subprocess, "run", side_effect=[expired, success],
            ) as run_mock,
        ):
            with self.assertRaisesRegex(
                review_isolated.ReviewError, r"timed out after 900 seconds: review stalled",
            ):
                review_isolated.run_review(self.root, str(transaction), "paper")
            incomplete = transaction / "reviews/paper"
            self.assertEqual(
                (incomplete / "output/failed-stdout.bin").read_bytes(), b"partial Grok stdout",
            )
            self.assertEqual(
                (incomplete / "output/failed-stderr.bin").read_bytes(), b"review stalled",
            )
            self.assertIsNone(review_isolated._completed_review_evidence(incomplete))
            review = review_isolated.run_review(
                self.root, str(transaction), "paper", retry_incomplete=True,
            )

        self.assertTrue(
            (transaction / "rejected-reviews/paper-attempt-0001/output/failed-stdout.bin").is_file(),
        )
        self.assertEqual(review["verdict"], "PASS")
        self.assertEqual(
            [call.kwargs["timeout"] for call in run_mock.call_args_list],
            [review_isolated.GROK_REVIEW_TIMEOUT_SECONDS] * 2,
        )

    def test_completed_pass_and_fail_reviews_cannot_retry(self) -> None:
        for verdict in ("PASS", "FAIL"):
            with self.subTest(verdict=verdict):
                case_root = self.root / verdict.lower()
                case_root.mkdir()
                repo = MemoryRepo(case_root, [MemoryRepo.unit("paper", "research/a.md")])
                repo.write("research/a.md", "## Claim\nA\n")
                repo.commit()
                transaction = repo.prepare()
                repo.staged_page(transaction, "paper")
                add_review_attestation(transaction, "paper", verdict)
                report_before = (transaction / "reviews/paper/output/report.md").read_bytes()
                attestation_before = (transaction / "reviews/paper/attestation.json").read_bytes()
                with (
                    mock.patch.object(review_isolated, "grok_profile") as profile,
                    self.assertRaisesRegex(review_isolated.ReviewError, "completed .*refusing to overwrite"),
                ):
                    review_isolated.run_review(
                        case_root, str(transaction), "paper", retry_incomplete=True,
                    )
                profile.assert_not_called()
                self.assertEqual(
                    (transaction / "reviews/paper/output/report.md").read_bytes(), report_before,
                )
                self.assertEqual(
                    (transaction / "reviews/paper/attestation.json").read_bytes(), attestation_before,
                )
                self.assertFalse((transaction / "rejected-reviews").exists())

        parsed = review_isolated.build_parser().parse_args([
            "transaction", "--task", "paper", "--retry-incomplete",
        ])
        self.assertTrue(parsed.retry_incomplete)

    def test_prepare_batch_bounds_and_review_prompt_cap_require_explicit_control(self) -> None:
        units = [
            MemoryRepo.unit("one", "research/a.md"),
            MemoryRepo.unit("two", "research/b.md"),
        ]
        repo = MemoryRepo(self.root, units)
        config = repo.config(units)
        config["read_limits"]["prepare_batch_max_units"] = 1
        config["read_limits"]["review_prompt_max_bytes"] = 32
        repo.write("memory/_meta/config.yaml", yaml.safe_dump(config, sort_keys=False))
        repo.write("research/a.md", "## A\nA\n")
        repo.write("research/b.md", "## B\nB\n")
        repo.commit()
        with self.assertRaisesRegex(mem.MemoryErrorBase, "--allow-large-batch"):
            repo.prepare()
        transaction = mem.prepare_update(self.root, "HEAD", [], [], allow_large_batch=True)
        manifest = json.loads((transaction / "manifest.json").read_text())
        self.assertTrue(manifest["batch_bounds"]["explicit_large_batch_override"])

        review_packet = self.root / "review-packet"
        shutil.copytree(transaction / "tasks/one/packet", review_packet)
        mem.atomic_write(review_packet / "candidate.md", b"candidate\n", 0o444)
        with self.assertRaisesRegex(review_isolated.ReviewError, "review prompt is"):
            review_isolated.compose_review_prompt(review_packet)
        with self.assertRaisesRegex(review_isolated.ReviewError, "actionable finding"):
            review_isolated._parse_structured_review(
                '{"findings":[{"severity":"note"}],"verdict":"FAIL"}'
            )
        final = review_isolated._structured_from_envelope({
            "text": '{"findings":[],"verdict":"FAIL"}{"findings":[],"verdict":"PASS"}',
            "structuredOutput": {"findings": [], "verdict": "PASS"},
        })
        self.assertEqual(final["verdict"], "PASS")

    def test_capsule_required_sections_reject_unstructured_prose(self) -> None:
        body = """## Purpose and scope

Unstructured claim.

## Computed evidence represented by the source

### source-example--evidence — Evidence

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

Scoped claim.

Sources:

- `research/example.md` — heading `## Evidence`

## Assumptions, exclusions, and open questions

## Revision and supersession relationships

"""
        self.assertEqual(
            mem._unstructured_required_capsule_sections(body),
            ["Purpose and scope"],
        )

    def test_markdown_heading_anchor_allows_escaped_code_ticks(self) -> None:
        source = b"## Form control \xe2\x80\x94 stiffness `Sigma`\n"
        self.assertTrue(
            mem._validate_anchor(
                source,
                r"heading `## Form control \u2014 stiffness \`Sigma\``".replace(r"\u2014", "\u2014"),
            )
        )

    def test_semantic_limits_and_identity_only_contents(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("large", "research/large.out")])
        repo.write("research/large.out", b"x" * 100001)
        repo.commit()
        with self.assertRaisesRegex(mem.ConfigError, "above its"):
            mem.status_report(self.root)
        repo.write_config([MemoryRepo.unit("large", "research/large.out", read_mode="identity_only")])
        transaction = repo.prepare()
        manifest = json.loads((transaction / "manifest.json").read_text())
        member = manifest["units"][0]["members"][0]
        self.assertNotIn("snapshot_path", member)
        self.assertFalse(any(p.is_file() for p in (transaction / "snapshot").rglob("*")))

    def test_deterministic_excerpt_requires_and_bounds_declared_anchor(self) -> None:
        data = ("prefix " * 100 + "DISTINCT-TAG" + " suffix" * 100).encode()
        member = {"path": "research/a.md", "anchors": ["DISTINCT-TAG"], "context_bytes": 80}
        excerpt = mem.deterministic_excerpt(data, member, 8192)
        self.assertIn(b"DISTINCT-TAG", excerpt)
        self.assertLess(len(excerpt), 200)
        with self.assertRaisesRegex(mem.MemoryErrorBase, "not found"):
            mem.deterministic_excerpt(data, {**member, "anchors": ["MISSING"]}, 8192)

    def test_policy_drift_blocks_finalize(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        repo.staged_page(transaction, "paper")
        repo.write("memory/_meta/schema.md", "# Changed schema\n")
        with self.assertRaisesRegex(mem.MemoryErrorBase, "digest drift"):
            mem.finalize_update(self.root, str(transaction))

    def test_finalize_requires_isolation_attestation(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        repo.staged_page(transaction, "paper")
        (transaction / "attestations/paper.json").unlink()
        with self.assertRaisesRegex(mem.LintFailure, "missing bubblewrap isolation attestation"):
            mem.finalize_update(self.root, str(transaction))

    def test_normal_codex_finalize_requires_current_grok_pass(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        missing = repo.prepare()
        repo.staged_page(missing, "paper")
        writer_attestation_path = missing / "attestations/paper.json"
        writer_attestation = json.loads(writer_attestation_path.read_text())
        writer_attestation["runtime_profile"] = "codex"
        mem.write_json(writer_attestation_path, writer_attestation)
        with self.assertRaisesRegex(mem.LintFailure, "missing required Grok PASS review"):
            mem.finalize_update(self.root, str(missing))

        failed = repo.prepare()
        repo.staged_page(failed, "paper")
        add_review_attestation(failed, "paper", "FAIL")
        with self.assertRaisesRegex(mem.LintFailure, "verdict"):
            mem.finalize_update(self.root, str(failed))

        passed = repo.prepare()
        repo.staged_page(passed, "paper")
        add_review_attestation(passed, "paper", "PASS")
        result = mem.finalize_update(self.root, str(passed))
        self.assertEqual(result["transaction_id"], json.loads((passed / "manifest.json").read_text())["transaction_id"])

    def test_normal_review_gate_rejects_candidate_report_writer_packet_and_reviewer_tamper(self) -> None:
        for tamper in ("candidate", "report", "writer", "packet", "reviewer"):
            with self.subTest(tamper=tamper):
                case_root = self.root / tamper
                case_root.mkdir()
                repo = MemoryRepo(case_root, [MemoryRepo.unit("paper", "research/a.md")])
                repo.write("research/a.md", "## Claim\nA\n")
                repo.commit()
                transaction = repo.prepare()
                candidate = repo.staged_page(transaction, "paper")
                add_review_attestation(transaction, "paper", "PASS")
                if tamper == "candidate":
                    candidate.write_text(candidate.read_text() + "\npost-review change\n", encoding="utf-8")
                elif tamper == "report":
                    (transaction / "reviews/paper/output/report.md").write_text("changed\n", encoding="utf-8")
                elif tamper == "writer":
                    path = transaction / "attestations/paper.json"
                    value = json.loads(path.read_text())
                    value["post_review_change"] = True
                    mem.write_json(path, value)
                elif tamper == "packet":
                    packet_schema = transaction / "tasks/paper/packet/schema.md"
                    packet_schema.chmod(0o644)
                    packet_schema.write_text("changed\n", encoding="utf-8")
                else:
                    path = transaction / "reviews/paper/attestation.json"
                    value = json.loads(path.read_text())
                    value["reviewer_sha256"] = mem.sha256_bytes(b"different reviewer")
                    mem.write_json(path, value)
                with self.assertRaises(mem.MemoryErrorBase):
                    mem.finalize_update(case_root, str(transaction))

    def test_bubblewrap_runner_hides_repo_and_imports_only_declared_output(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        script = (
            f"test ! -e {self.root.as_posix()}/research && "
            "test -r /packet/schema.md && "
            "test -r /packet/inputs/semantic/research/a.md && "
            "! git rev-parse --show-toplevel >/dev/null 2>&1 && "
            "printf isolated > /output/page.md"
        )
        attestation = run_isolated.run_task(self.root, str(transaction), "paper", ["/bin/sh", "-c", script])
        self.assertEqual(attestation["isolation"], "bubblewrap")
        self.assertEqual(
            (transaction / "staged/memory/sources/paper.md").read_text(),
            "isolated",
        )

    def test_derived_tasks_are_ordered_and_hydrate_only_attested_dependencies(self) -> None:
        unit = MemoryRepo.unit("paper", "research/a.md")
        repo = MemoryRepo(self.root, [unit])
        for name in ("topic-synthesis.md", "index-refresh.md"):
            repo.write(f"memory/prompts/{name}", f"# {name}\n")
        config = repo.config([unit])
        config["output_budgets"] = {"source_capsule": {"max_words": 100}, "topic": {"max_words": 100}, "index": {"max_words": 100}}
        config["derived_pages"] = [
            {
                "id": "topic-one", "task_kind": "topic", "page": "memory/topics/one.md",
                "region": "working-position", "order": 20, "input_units": ["paper"],
                "input_pages": [], "budget": "topic",
                "direct_sources": [{"path": "research/a.md", "read_mode": "semantic"}],
            },
            {
                "id": "memory-index", "task_kind": "index", "page": "memory/index.md",
                "region": "navigation", "order": 50, "input_units": [],
                "input_pages": ["topic-one"], "budget": "index",
            },
        ]
        repo.write("memory/_meta/config.yaml", yaml.safe_dump(config, sort_keys=False))
        repo.write(
            "memory/_meta/atlas-migration.yaml",
            yaml.safe_dump({"items": [{
                "id": "legacy-topic", "disposition": "migrate",
                "target": "memory/topics/one.md", "original_sources": ["research/a.md"],
            }]}, sort_keys=False),
        )
        repo.write("research/a.md", "## Claim\nA\n")
        mixed = """---
schema_version: 2
id: {page_id}
title: {title}
type: {page_type}
lifecycle: current
memory_review: ai_draft
sources: []
content_owner: mixed
last_updated: 2026-08-24
generated_from_commit: pending
---

<!-- BEGIN GENERATED:{region} -->
pending
<!-- END GENERATED:{region} -->

## Curator notes
keep
"""
        mixed = mixed.replace("keep\n", ("keep " * 150) + "\n")
        repo.write("memory/topics/one.md", mixed.format(page_id="topic-one", title="Topic one", page_type="topic", region="working-position"))
        repo.write("memory/index.md", mixed.format(page_id="memory-index", title="Index", page_type="index", region="navigation"))
        repo.commit()
        transaction = repo.prepare()
        manifest = json.loads((transaction / "manifest.json").read_text())
        self.assertEqual([phase["phase"] for phase in manifest["task_order"]], [1, 20, 50])
        self.assertEqual(manifest["derived_task_plan"]["status"], "complete")
        topic_packet = transaction / "tasks/topic-one/packet"
        self.assertTrue((topic_packet / "base_page.md").is_file())
        self.assertTrue((topic_packet / "inputs/direct/semantic/research/a.md").is_file())
        self.assertTrue((topic_packet / "atlas-migration.yaml").is_file())
        topic_task = json.loads((topic_packet / "task.json").read_text())
        dynamic_source = topic_task["semantic_contract"]["dynamic_memory_inputs"][0]
        self.assertEqual(dynamic_source["generated_commit"], manifest["target_commit"])
        self.assertTrue(dynamic_source["policy_fresh"])
        self.assertEqual(
            topic_task["semantic_contract"]["migration_requirements"][0]["legacy_id"], "legacy-topic"
        )

        source_staged = repo.staged_page(transaction, "paper")
        source_staged.write_text(
            source_staged.read_text().replace(manifest["target_commit"], "0" * 40), encoding="utf-8"
        )
        source_attestation_path = transaction / "attestations/paper.json"
        source_attestation = json.loads(source_attestation_path.read_text())
        source_attestation["output_sha256"] = mem.sha256_file(source_staged)
        mem.write_json(source_attestation_path, source_attestation)
        with self.assertRaisesRegex(run_isolated.IsolationError, "generated commit mismatch"):
            run_isolated.run_task(
                self.root, str(transaction), "topic-one",
                ["/usr/bin/python3", "-c", "from pathlib import Path; Path('/output/page.md').write_text('x')"],
            )
        repo.staged_page(transaction, "paper")
        topic_output = (self.root / "memory/topics/one.md").read_text()
        topic_output = topic_output.replace(
            "last_updated: 2026-08-24",
            f"last_updated: {topic_task['semantic_contract']['refresh_date']}",
        )
        topic_output = topic_output.replace("sources: []", "sources:\n- research/a.md")
        topic_output = topic_output.replace("generated_from_commit: pending", f"generated_from_commit: {manifest['target_commit']}")
        topic_output = topic_output.replace(
            "pending\n<!-- END GENERATED:working-position -->",
            "### topic-one--claim — Scoped claim\n\n"
            "Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`\n\n"
            "Prepared source position.\n\n"
            "Sources:\n\n- `research/a.md` — heading `## Claim`\n\n"
            "#### Migration record — legacy-topic\n\n"
            "- `legacy_id`: `legacy-topic`\n"
            "- `migration_disposition`: `migrated`\n"
            "- `target_statement_ids`: [`topic-one--claim`]\n"
            f"- `inventory_sha256`: `{manifest['policy']['atlas_migration_sha256']}`\n\n"
            "Sources:\n\n- `research/a.md` — heading `## Claim`\n"
            "<!-- END GENERATED:working-position -->",
        )
        topic_code = f"from pathlib import Path; Path('/output/page.md').write_text({topic_output!r})"
        run_isolated.run_task(
            self.root, str(transaction), "topic-one", ["/usr/bin/python3", "-c", topic_code]
        )
        self.assertTrue((topic_packet / "memory_inputs/memory/sources/paper.md").is_file())
        hydrated_task = json.loads((topic_packet / "task.json").read_text())
        self.assertEqual(
            hydrated_task["semantic_contract"]["dynamic_memory_inputs"][0]["page_sha256"],
            mem.sha256_file(transaction / "staged/memory/sources/paper.md"),
        )
        self.assertTrue((topic_packet / "hydration.json").is_file())
        index_output = (self.root / "memory/index.md").read_text()
        index_task = json.loads((transaction / "tasks/memory-index/packet/task.json").read_text())
        index_output = index_output.replace(
            "last_updated: 2026-08-24",
            f"last_updated: {index_task['semantic_contract']['refresh_date']}",
        )
        index_output = index_output.replace("generated_from_commit: pending", f"generated_from_commit: {manifest['target_commit']}")
        index_output = index_output.replace(
            "pending\n<!-- END GENERATED:navigation -->",
            "[Topic one](topics/one.md)\n<!-- END GENERATED:navigation -->",
        )
        index_code = f"from pathlib import Path; Path('/output/page.md').write_text({index_output!r})"
        run_isolated.run_task(self.root, str(transaction), "memory-index", ["/usr/bin/python3", "-c", index_code])
        self.assertTrue((transaction / "attestations/topic-one.json").is_file())
        self.assertTrue((transaction / "attestations/memory-index.json").is_file())
        live_topic = self.root / "memory/topics/one.md"
        original_topic = live_topic.read_bytes()
        live_topic.write_bytes(original_topic + b"\nchanged after prepare\n")
        with self.assertRaisesRegex(mem.MemoryErrorBase, "live mixed target changed since prepare"):
            mem.finalize_update(self.root, str(transaction))
        live_topic.write_bytes(original_topic)
        result = mem.finalize_update(self.root, str(transaction))
        self.assertIn("memory/topics/one.md", result["published_pages"])
        state, _ = mem.load_state(self.root)
        self.assertEqual(state["derived_pages"]["topic-one"]["result"], "success")
        self.assertEqual(state["derived_pages"]["topic-one"]["migration_ids"], ["legacy-topic"])
        self.assertEqual(state["derived_pages"]["topic-one"]["dependency_inputs"][0]["page_sha256"],
                         state["pages"]["memory/sources/paper.md"]["sha256"])
        original_generated_commit = state["pages"]["memory/topics/one.md"]["generated_commit"]
        repo.write("research/unclaimed.md", "advance without semantic changes\n")
        repo.commit("rebase only")
        rebase = repo.prepare(units=["paper"])
        rebase_manifest = json.loads((rebase / "manifest.json").read_text())
        self.assertTrue(all(not task["required"] for task in rebase_manifest["writer_tasks"]))
        self.assertEqual(
            {task["action"] for task in rebase_manifest["writer_tasks"] if task["source_unit_id"] is None},
            {"derived_rebase_unchanged"},
        )
        mem.finalize_update(self.root, str(rebase))
        state, _ = mem.load_state(self.root)
        self.assertEqual(state["derived_pages"]["topic-one"]["processed_commit"], git(self.root, "rev-parse", "HEAD"))
        self.assertEqual(state["pages"]["memory/topics/one.md"]["generated_commit"], original_generated_commit)
        config["read_limits"]["derived_task_max_bytes"] = 1
        repo.write("memory/_meta/config.yaml", yaml.safe_dump(config, sort_keys=False))
        with self.assertRaisesRegex(mem.ConfigError, "deduplicated direct semantic bytes"):
            repo.prepare(paths=["research/a.md"])

    def test_direct_source_cannot_escalate_configured_excerpt(self) -> None:
        unit = MemoryRepo.unit(
            "paper", "research/a.md", read_mode="excerpt", anchors=[{"value": "## Claim"}],
        )
        repo = MemoryRepo(self.root, [unit])
        repo.write("memory/prompts/topic-synthesis.md", "# Topic\n")
        config = repo.config([unit])
        config["output_budgets"] = {"topic": {"max_words": 50}}
        config["derived_pages"] = [{
            "id": "topic-one", "task_kind": "topic", "page": "memory/topics/one.md",
            "region": "working-position", "order": 20, "input_units": ["paper"], "input_pages": [],
            "budget": "topic", "direct_sources": [{"path": "research/a.md", "read_mode": "semantic"}],
        }]
        repo.write("memory/_meta/config.yaml", yaml.safe_dump(config, sort_keys=False))
        repo.write("research/a.md", "## Claim\nA\n")
        repo.write("memory/topics/one.md", """---
schema_version: 2
id: topic-one
title: One
type: topic
lifecycle: current
memory_review: ai_draft
sources: []
content_owner: mixed
last_updated: 2026-08-24
generated_from_commit: pending
---
<!-- BEGIN GENERATED:working-position -->
pending
<!-- END GENERATED:working-position -->
""")
        repo.commit()
        with self.assertRaisesRegex(mem.ConfigError, "escalates research/a.md from excerpt to semantic"):
            repo.prepare()

    def test_mixed_generated_region_entry_budgets_are_enforced(self) -> None:
        unit = MemoryRepo.unit("paper", "research/a.md")
        repo = MemoryRepo(self.root, [unit])
        repo.write("memory/prompts/script-catalog.md", "# Script task\n")
        config = repo.config([unit])
        config["output_budgets"] = {
            "source_capsule": {"max_words": 100},
            "script_catalog": {"max_words": 100, "max_entries": 1, "max_entry_words": 3},
        }
        config["derived_pages"] = [{
            "id": "scripts-one", "task_kind": "script_catalog", "page": "memory/scripts/one.md",
            "region": "entries", "order": 30, "input_units": ["paper"], "input_pages": [],
            "budget": "script_catalog",
            "direct_sources": [{"path": "research/a.md", "read_mode": "semantic"}],
        }]
        repo.write("memory/_meta/config.yaml", yaml.safe_dump(config, sort_keys=False))
        repo.write("research/a.md", "## Claim\nA\n")
        repo.write("memory/scripts/one.md", """---
schema_version: 2
id: scripts-one
title: Scripts one
type: script_catalog
lifecycle: current
memory_review: ai_draft
sources: []
content_owner: mixed
last_updated: 2026-08-24
generated_from_commit: pending
---

<!-- BEGIN GENERATED:entries -->
pending
<!-- END GENERATED:entries -->

## Curator notes
""" + ("outside " * 200) + "\n")
        repo.commit()
        transaction = repo.prepare()
        manifest = json.loads((transaction / "manifest.json").read_text())
        repo.staged_page(transaction, "paper")
        output = (self.root / "memory/scripts/one.md").read_text()
        output = output.replace("sources: []", "sources:\n- research/a.md")
        output = output.replace("generated_from_commit: pending", f"generated_from_commit: {manifest['target_commit']}")
        output = output.replace(
            "pending\n<!-- END GENERATED:entries -->",
            "## `one`\n\nalpha beta gamma delta\n\n## `two`\n\none two three four\n"
            "<!-- END GENERATED:entries -->",
        )
        code = f"from pathlib import Path; Path('/output/page.md').write_text({output!r})"
        run_isolated.run_task(self.root, str(transaction), "scripts-one", ["/usr/bin/python3", "-c", code])
        errors, _ = mem.lint_staged(self.root, str(transaction))
        self.assertTrue(any("above its 1-entry task budget" in error for error in errors))
        self.assertTrue(any("above its 3-word entry budget" in error for error in errors))

    def test_partial_baseline_addition_deletion_and_retirement(self) -> None:
        units = [
            MemoryRepo.unit("one", "research/one.md"),
            MemoryRepo.unit("two", "research/two.md"),
        ]
        repo = MemoryRepo(self.root, units)
        repo.write("research/one.md", "## Claim\none\n")
        repo.write("research/two.md", "## Claim\ntwo\n")
        repo.commit()
        first = repo.prepare(units=["one"])
        repo.staged_page(first, "one")
        result = mem.finalize_update(self.root, str(first))
        self.assertIsNone(result["last_fully_processed_commit"])
        second = repo.prepare()
        repo.staged_page(second, "two")
        result = mem.finalize_update(self.root, str(second))
        self.assertEqual(result["last_fully_processed_commit"], git(self.root, "rev-parse", "HEAD"))

        (self.root / "research/one.md").unlink()
        repo.commit("delete one")
        deleted = repo.prepare(paths=["research/one.md"])
        self.assertTrue((deleted / "tasks/one/packet/inputs/prior/one/research/one.md").is_file())
        repo.staged_page(deleted, "one", lifecycle="deleted")
        mem.finalize_update(self.root, str(deleted))
        state, _ = mem.load_state(self.root)
        self.assertEqual(state["units"]["one"]["result"], "source_deleted")
        deletion_generated_commit = state["pages"]["memory/sources/one.md"]["generated_commit"]

        repo.write("research/unclaimed.md", "new commit without restoring one\n")
        repo.commit("advance after deletion")
        deletion_rebase = repo.prepare(units=["one"])
        rebased_manifest = json.loads((deletion_rebase / "manifest.json").read_text())
        self.assertFalse(next(task for task in rebased_manifest["writer_tasks"] if task["task_id"] == "one")["required"])
        mem.finalize_update(self.root, str(deletion_rebase))
        state, _ = mem.load_state(self.root)
        self.assertEqual(state["units"]["one"]["result"], "source_deleted")
        self.assertEqual(state["units"]["one"]["processed_commit"], git(self.root, "rev-parse", "HEAD"))
        self.assertEqual(state["pages"]["memory/sources/one.md"]["generated_commit"], deletion_generated_commit)

        repo.write_config([units[1]])
        retired = repo.prepare(units=["one"])
        retired_manifest = json.loads((retired / "manifest.json").read_text())
        desired = next(
            task["semantic_contract"]["page"]["desired_lifecycle"]
            for task in retired_manifest["writer_tasks"] if task["task_id"] == "one"
        )
        repo.staged_page(retired, "one", lifecycle=desired)
        mem.finalize_update(self.root, str(retired))
        state, _ = mem.load_state(self.root)
        self.assertNotIn("one", state["units"])
        self.assertEqual(state["retired_units"]["one"]["result"], "retired_from_corpus")

    def test_handled_failure_restores_pages_and_state(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nold\n")
        repo.commit()
        first = repo.prepare()
        repo.staged_page(first, "paper", title="Old")
        mem.finalize_update(self.root, str(first))
        old_page = (self.root / "memory/sources/paper.md").read_bytes()
        old_state = (self.root / mem.DEFAULT_STATE).read_bytes()

        repo.write("research/a.md", "## Claim\nnew\n")
        repo.commit("change")
        second = repo.prepare()
        repo.staged_page(second, "paper", title="New")
        os.environ["MEMORY_TEST_FAIL_AFTER_PAGES"] = "1"
        try:
            with self.assertRaisesRegex(OSError, "injected"):
                mem.finalize_update(self.root, str(second))
        finally:
            os.environ.pop("MEMORY_TEST_FAIL_AFTER_PAGES", None)
        self.assertEqual((self.root / "memory/sources/paper.md").read_bytes(), old_page)
        self.assertEqual((self.root / mem.DEFAULT_STATE).read_bytes(), old_state)
        self.assertFalse((self.root / "memory/journals/publication.json").exists())

    def test_interrupted_journal_recovery_rolls_back_or_cleans_up(self) -> None:
        repo = MemoryRepo(self.root, [])
        repo.commit()
        mem.init_memory(self.root)
        original_state = (self.root / mem.DEFAULT_STATE).read_bytes()
        page = repo.write("memory/sources/page.md", "new partial page\n")
        backup_root = self.root / "memory/backups/fake"
        backup_root.mkdir(parents=True)
        (backup_root / mem.DEFAULT_STATE).parent.mkdir(parents=True, exist_ok=True)
        (backup_root / mem.DEFAULT_STATE).write_bytes(original_state)
        journal = {
            "backup_root": "memory/backups/fake",
            "targets": [{"path": "memory/sources/page.md", "existed": False, "backup": "memory/sources/page.md"}],
            "state": {"existed": True, "backup": mem.DEFAULT_STATE},
            "new_state_sha256": "0" * 64,
        }
        mem.write_json(self.root / "memory/journals/publication.json", journal)
        outcome = mem.recover_publication(self.root)
        self.assertIn("rolled back", outcome or "")
        self.assertFalse(page.exists())
        self.assertEqual((self.root / mem.DEFAULT_STATE).read_bytes(), original_state)

        # If state already carries the journaled hash, recovery keeps the new
        # publication and only removes backup/journal debris.
        backup_root.mkdir(parents=True)
        current_hash = mem.sha256_file(self.root / mem.DEFAULT_STATE)
        journal["new_state_sha256"] = current_hash
        mem.write_json(self.root / "memory/journals/publication.json", journal)
        outcome = mem.recover_publication(self.root)
        self.assertIn("completed cleanup", outcome or "")
        self.assertFalse(backup_root.exists())

    def test_anchor_and_internal_link_lint(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        staged = repo.staged_page(transaction, "paper", bad_anchor=True, bad_link=True)
        staged.write_text(
            staged.read_text()
            .replace("source-paper--test-statement", "wrong-namespace")
            .replace("evidence=provisional", "evidence=measured")
            + "\n[bad local fragment](#does-not-exist)\n",
            encoding="utf-8",
        )
        errors, _ = mem.lint_staged(self.root, str(transaction))
        self.assertTrue(any("anchor not found" in error for error in errors), errors)
        self.assertTrue(any("unresolved internal link" in error for error in errors), errors)
        self.assertTrue(any("unresolved internal link fragment" in error for error in errors), errors)
        self.assertTrue(any("outside page namespace" in error for error in errors), errors)
        self.assertTrue(any("cannot assign measured" in error for error in errors), errors)

    def test_effective_ids_and_identity_only_citations_are_rejected(self) -> None:
        first_unit = MemoryRepo.unit("one", "research/one.md")
        repo = MemoryRepo(self.root, [first_unit])
        repo.write("research/one.md", "## Claim\none\n")
        repo.commit()
        first = repo.prepare()
        repo.staged_page(first, "one")
        mem.finalize_update(self.root, str(first))

        second_unit = MemoryRepo.unit("two", "research/two.md")
        repo.write_config([first_unit, second_unit])
        repo.write("research/two.md", "## Claim\ntwo\n")
        repo.commit("add second")
        second = repo.prepare(units=["two"])
        staged = repo.staged_page(second, "two")
        staged.write_text(staged.read_text().replace("id: source-two", "id: source-one"), encoding="utf-8")
        errors, _ = mem.lint_staged(self.root, str(second))
        self.assertTrue(any("duplicate effective page id source-one" in error for error in errors), errors)

        identity_repo_root = self.root / "identity"
        identity_repo_root.mkdir()
        identity_repo = MemoryRepo(
            identity_repo_root,
            [MemoryRepo.unit("identity", "research/id.md", read_mode="identity_only")],
        )
        identity_repo.write("research/id.md", "## Claim\nidentity\n")
        identity_repo.commit()
        identity_txn = identity_repo.prepare()
        identity_repo.staged_page(identity_txn, "identity")
        errors, _ = mem.lint_staged(identity_repo_root, str(identity_txn))
        self.assertTrue(any("identity-only source" in error for error in errors), errors)

        identity_page = identity_txn / "staged/memory/sources/identity.md"
        identity_page.write_text(
            identity_page.read_text(encoding="utf-8")
            .replace("Fixture variant: Test", "This artifact is identity-only; its contents were not inspected.")
            .replace("heading `## Claim`", "`anchor-unavailable`"),
            encoding="utf-8",
        )
        errors, _ = mem.lint_staged(identity_repo_root, str(identity_txn))
        self.assertFalse(any("identity-only source" in error for error in errors), errors)
        identity_page.write_text(
            identity_page.read_text(encoding="utf-8").replace(
                "This artifact is identity-only; its contents were not inspected.",
                "This artifact is not semantically readable; its contents were not inspected.",
            ),
            encoding="utf-8",
        )
        errors, _ = mem.lint_staged(identity_repo_root, str(identity_txn))
        self.assertFalse(any("identity-only source" in error for error in errors), errors)

    def test_mixed_page_preserves_unmarked_bytes_and_nondelegated_frontmatter(self) -> None:
        base = b"""---
schema_version: 2
id: topic-one
title: Curated title
type: topic
lifecycle: current
memory_review: ai_draft
sources: []
content_owner: mixed
last_updated: 2026-08-23
generated_from_commit: old
---

<!-- BEGIN GENERATED:working-position -->
old generated text
<!-- END GENERATED:working-position -->

## Curator notes

Keep this byte-for-byte.
"""
        good = base.replace(b"old generated text", b"new generated text").replace(
            b"generated_from_commit: old", b"generated_from_commit: new"
        )
        self.assertEqual(mem.validate_mixed_preservation(
            base, good, declared_regions=["working-position"],
            delegated_frontmatter=["generated_from_commit"], label="topic",
        ), [])
        bad = good.replace(b"Keep this byte-for-byte.", b"Changed curator prose.")
        errors = mem.validate_mixed_preservation(
            base, bad, declared_regions=["working-position"],
            delegated_frontmatter=["generated_from_commit"], label="topic",
        )
        self.assertTrue(any("unmarked" in error for error in errors), errors)

    def test_fallback_anchor_and_nested_conflict_source_syntax_are_recognized(self) -> None:
        self.assertTrue(mem._validate_anchor(b"anything", "`anchor-unavailable`"))
        self.assertTrue(mem._validate_anchor(
            b"**LOUD LABEL:** frozen result\n",
            "exact text `LOUD LABEL: frozen result`",
        ))
        citations = mem.CITATION_RE.findall(
            "- Position A: scoped text\n  - Source: `research/a.md` — heading `## Claim`\n"
        )
        self.assertEqual(citations, [("research/a.md", "heading `## Claim`")])

    def test_query_serves_only_recorded_pages_and_warns_on_stale_units(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nold vocabulary\n")
        repo.commit()
        transaction = repo.prepare()
        repo.staged_page(transaction, "paper", title="Old vocabulary result")
        mem.finalize_update(self.root, str(transaction))
        repo.write("research/a.md", "## Claim\nbrandnewterm\n")
        repo.commit("new source vocabulary")
        result = mem.query_memory(self.root, "brandnewterm")
        self.assertTrue(result["warnings"])
        self.assertIn("pending units: paper", result["warnings"][0])
        self.assertEqual(result["results"], [])
        (self.root / "memory/sources/paper.md").write_text("tampered", encoding="utf-8")
        status = mem.status_report(self.root)
        self.assertFalse(status["served_integrity"]["ok"])
        self.assertFalse(status["cutover"]["ready"])
        with self.assertRaises(mem.LintFailure):
            mem.query_memory(self.root, "old")

    def test_successful_state_records_stable_id_registry(self) -> None:
        repo = MemoryRepo(self.root, [MemoryRepo.unit("paper", "research/a.md")])
        repo.write("research/a.md", "## Claim\nA\n")
        repo.commit()
        transaction = repo.prepare()
        repo.staged_page(transaction, "paper")
        mem.finalize_update(self.root, str(transaction))
        state, _ = mem.load_state(self.root)
        self.assertEqual(state["id_registry"]["pages"]["source-paper"], "memory/sources/paper.md")
        self.assertEqual(
            state["id_registry"]["statements"]["source-paper--test-statement"],
            "memory/sources/paper.md",
        )


if __name__ == "__main__":
    unittest.main()
