from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from memory.tools import transaction_gc as gc


def git(repo: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", "-C", str(repo), *args], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if proc.returncode:
        raise AssertionError(f"git {' '.join(args)} failed: {proc.stderr}")
    return proc.stdout.strip()


class TransactionGCTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.repo = Path(self.temporary.name)
        git(self.repo, "init", "-q")
        git(self.repo, "config", "user.email", "tests@example.invalid")
        git(self.repo, "config", "user.name", "Transaction GC Tests")
        git(self.repo, "config", "commit.gpgsign", "false")
        (self.repo / "tracked.txt").write_text("one\n", encoding="utf-8")
        git(self.repo, "add", "tracked.txt")
        git(self.repo, "commit", "-qm", "initial")
        self.head = git(self.repo, "rev-parse", "HEAD")
        (self.repo / "memory/transactions").mkdir(parents=True)
        (self.repo / "memory/_meta").mkdir(parents=True)
        self.write_state()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_state(
        self, *, generation: str | None = None, page_path: str | None = None,
        page_sha256: str | None = None,
    ) -> None:
        state = {
            "state_version": 2,
            "generation": generation,
            "last_result": {"transaction_id": generation, "status": "success" if generation else "uninitialized"},
            "pages": {}, "units": {}, "derived_pages": {},
        }
        if page_path is not None and page_sha256 is not None:
            state["pages"][page_path] = {"sha256": page_sha256}
        (self.repo / "memory/_meta/state.json").write_text(
            json.dumps(state, sort_keys=True), encoding="utf-8",
        )

    def make_transaction(
        self, transaction_id: str, *, target: str | None = None,
        prior_transaction_id: str | None = None, failure: bool = False,
    ) -> Path:
        root = self.repo / "memory/transactions" / transaction_id
        packet = root / "tasks/paper/packet"
        packet.mkdir(parents=True)
        staged = root / "staged/memory/sources/paper.md"
        staged.parent.mkdir(parents=True)
        staged.write_text(f"candidate from {transaction_id}\n" + ("x" * 10000), encoding="utf-8")
        manifest = {
            "transaction_id": transaction_id,
            "created_at": "2026-08-25T00:00:00+00:00",
            "target_commit": target or self.head,
            "policy": {"combined_sha256": "a" * 64},
            "units": [{"id": "paper", "unit_digest_sha256": "b" * 64}],
            "writer_tasks": [{
                "task_id": "paper", "source_unit_id": "paper", "required": True,
                "output_repository_path": "memory/sources/paper.md",
                "staged_output_path": "staged/memory/sources/paper.md",
            }],
        }
        (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        (root / "seal.json").write_text(json.dumps({"transaction_id": transaction_id}), encoding="utf-8")
        (packet / "task.json").write_text(json.dumps({"transaction_id": transaction_id}), encoding="utf-8")
        if prior_transaction_id is not None:
            (packet / "revision.json").write_text(json.dumps({
                "transaction_id": transaction_id,
                "prior_transaction_id": prior_transaction_id,
            }), encoding="utf-8")
        if failure:
            (root / "failure.json").write_text(json.dumps({"error": "fixture"}), encoding="utf-8")
        return root

    def test_archive_is_dry_run_then_lossless_prune_and_restore(self) -> None:
        old_target = "0" * 40
        transaction = self.make_transaction("old", target=old_target)
        original = gc.inventory_summary(gc.tree_inventory(transaction))

        planned = gc.archive_transaction(self.repo, "old", apply=False)
        self.assertFalse((self.repo / planned["archive"]).exists())

        archived = gc.archive_transaction(self.repo, "old", apply=True)
        self.assertEqual(archived["status"], "archived_and_verified")
        self.assertEqual(archived["tree_sha256"], original["tree_sha256"])
        self.assertTrue(transaction.is_dir(), "archive must never remove its source")
        verified = gc.verify_archive(self.repo, "old")
        self.assertEqual(verified["tree_sha256"], original["tree_sha256"])
        self.assertLess(verified["archive_bytes"], original["byte_count"])

        # A crash after installing the archive but before committing its
        # receipt is recoverable without overwriting either source or archive.
        gc.receipt_path(self.repo, "old").unlink()
        recovered = gc.archive_transaction(self.repo, "old", apply=True)
        self.assertEqual(recovered["status"], "recovered_receipt_and_verified")

        plan = gc.prune_transaction(self.repo, "old", apply=False)
        self.assertEqual(plan["status"], "dry_run")
        self.assertTrue(transaction.is_dir())
        pruned = gc.prune_transaction(self.repo, "old", apply=True)
        self.assertEqual(pruned["status"], "source_removed_archive_retained")
        self.assertFalse(transaction.exists())

        restore_plan = gc.restore_transaction(self.repo, "old", apply=False)
        self.assertEqual(restore_plan["status"], "dry_run")
        restored = gc.restore_transaction(self.repo, "old", apply=True)
        self.assertEqual(restored["status"], "restored_and_verified")
        self.assertEqual(gc.inventory_summary(gc.tree_inventory(transaction)), original)

    def test_current_state_and_reference_chains_are_protected(self) -> None:
        prior = self.make_transaction("prior", target="0" * 40)
        candidate = prior / "staged/memory/sources/paper.md"
        self.write_state(generation="prior", page_path="memory/sources/paper.md", page_sha256=gc.sha256_file(candidate))
        self.make_transaction("child", prior_transaction_id="prior")
        gc.archive_transaction(self.repo, "prior", apply=True)

        status = gc.classify_transactions(self.repo)["prior"]
        self.assertIn("current_state_generation", status["protected_reasons"])
        self.assertIn("serves_current_page_hash", status["protected_reasons"])
        self.assertIn("referenced_by_transaction", status["protected_reasons"])
        self.assertEqual(status["references_in"], ["child"])
        with self.assertRaisesRegex(gc.TransactionGCError, "not prune-eligible"):
            gc.prune_transaction(self.repo, "prior", apply=True)
        self.assertTrue(prior.is_dir())

    def test_candidate_reuse_reference_chain_is_protected(self) -> None:
        prior = self.make_transaction("prior-candidate", target="0" * 40)
        child = self.make_transaction("recheck-child")
        attestation = child / "attestations/paper.json"
        attestation.parent.mkdir()
        attestation.write_text(json.dumps({
            "candidate_reuse": {"prior_transaction_id": "prior-candidate"},
        }), encoding="utf-8")
        gc.archive_transaction(self.repo, "prior-candidate", apply=True)

        status = gc.classify_transactions(self.repo)["prior-candidate"]
        self.assertIn("referenced_by_transaction", status["protected_reasons"])
        self.assertEqual(status["references_in"], ["recheck-child"])
        with self.assertRaisesRegex(gc.TransactionGCError, "not prune-eligible"):
            gc.prune_transaction(self.repo, "prior-candidate", apply=True)
        self.assertTrue(prior.is_dir())

        reviewed_prior = self.make_transaction("prior-reviewed", target="0" * 40)
        reviewed_child = self.make_transaction("reviewed-child")
        reviewed_attestation = reviewed_child / "attestations/paper.json"
        reviewed_attestation.parent.mkdir()
        reviewed_attestation.write_text(json.dumps({
            "reviewed_reuse": {"prior_transaction_id": "prior-reviewed"},
        }), encoding="utf-8")
        gc.archive_transaction(self.repo, "prior-reviewed", apply=True)
        reviewed_status = gc.classify_transactions(self.repo)["prior-reviewed"]
        self.assertIn("referenced_by_transaction", reviewed_status["protected_reasons"])
        self.assertEqual(reviewed_status["references_in"], ["reviewed-child"])
        with self.assertRaisesRegex(gc.TransactionGCError, "not prune-eligible"):
            gc.prune_transaction(self.repo, "prior-reviewed", apply=True)
        self.assertTrue(reviewed_prior.is_dir())

    def test_finalization_receipt_proves_later_supersession(self) -> None:
        transaction = self.make_transaction("published")
        candidate = transaction / "staged/memory/sources/paper.md"
        self.write_state(
            generation="published", page_path="memory/sources/paper.md",
            page_sha256=gc.sha256_file(candidate),
        )
        dry_run = gc.record_finalized(self.repo, "published", apply=False)
        self.assertEqual(dry_run["status"], "finalized")
        self.assertFalse(gc.lifecycle_path(self.repo, "published").exists())
        gc.record_finalized(self.repo, "published", apply=True)
        gc.archive_transaction(self.repo, "published", apply=True)

        self.write_state(generation="replacement", page_path="memory/sources/paper.md", page_sha256="f" * 64)
        status = gc.classify_transactions(self.repo)["published"]
        self.assertEqual(status["lifecycle_status"], "finalized_superseded")
        self.assertTrue(status["prune_eligible"])

    def test_tampered_archive_and_symlink_source_are_rejected(self) -> None:
        transaction = self.make_transaction("tamper", target="0" * 40)
        gc.archive_transaction(self.repo, "tamper", apply=True)
        with gc.archive_path(self.repo, "tamper").open("ab") as stream:
            stream.write(b"tamper")
        with self.assertRaisesRegex(gc.TransactionGCError, r"archive (?:size|hash) mismatch"):
            gc.verify_archive(self.repo, "tamper")
        with self.assertRaisesRegex(gc.TransactionGCError, r"archive (?:size|hash) mismatch"):
            gc.prune_transaction(self.repo, "tamper", apply=True)
        self.assertTrue(transaction.is_dir())

        linked = self.make_transaction("linked", target="0" * 40)
        (linked / "unsafe-link").symlink_to("manifest.json")
        with self.assertRaisesRegex(gc.TransactionGCError, "contains a symlink"):
            gc.archive_transaction(self.repo, "linked", apply=True)


if __name__ == "__main__":
    unittest.main()
