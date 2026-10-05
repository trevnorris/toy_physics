# Independent review: STATUS carry-over after the S11c cleanup fixes

## Artifact
In the repository `/var/projects/toy_physics`, Codex commit `19042324`:
- the rewritten `STATUS.md`;
- Part F of `research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md`, a carry-over table that gives every open, owed,
  deferred or debt item from the old STATUS a disposition: carried, resolved or superseded.

The old STATUS is `git show c8885c4f:STATUS.md`, 1,289 lines.

## What to check
1. **Completeness.** Is every live open, owed, deferred or debt item in the old STATUS accounted for in Part F?
   Read the old STATUS yourself and list any item Part F misses. A missed live debt is a finding.
2. **Dispositions.**
   - Each "resolved" must cite evidence that actually resolves the item. Read the cited file.
   - Each "superseded" must cite what superseded it, and the superseding text must actually replace the item.
   - Flag any debt marked resolved or superseded when the evidence only shows that work moved on, or that a
     NEXT/work-queue sentence was retired. Retiring a queue does not discharge a debt.
3. **Carried items.** Each one must appear in the new STATUS, or in a register that STATUS points to
   (`DEFECT_REGISTER.md`, `DEFERRED_HEAVY_RUNS.md`, `SUBSTRATE_REQUIREMENTS.md`, …). Check that the item is
   actually there, with an owner.
4. **New STATUS claims.** Any new sentence in STATUS stating a result or a review status must match its cited
   source.

## What you are handed
The repository and its history, including the tag `archive/pre-cleanup-2026-10-04`. The fix instructions are
`research/pde_ledger_v3/cleanup_2026_10/PHASE2_REVIEW.md`, item F6.

## Required method
This is a DOCUMENT review. Read the old STATUS before Part F. For each finding, quote both sides, with paths and
line numbers, and show the command you used (for example `grep -n` or `sed -n`) with its literal output. A finding
without these will be discarded.

## Physics filter
Report only an item that would let a live debt disappear from view, or that states a result or review status its
source does not support. No style or wording findings.

## Sandbox
Read-only. ⛔ Never modify the working tree, commit, or run CAS engines or numerical workers.

## Bounds
Write your report and exit. ⛔ Do not spawn agents. Format: a verdict line (CLEAR or NEEDS REVISION), then the
numbered findings, then what you checked and found sound.
