# Resume note 1 of the S15 Singular engine build, review pass 2: disposition

## The gate record

- **Artifact.** `research/pde_ledger_v3/directives/resume_s15_singular_engine_1.md`, revision 2, sha256 `cd284a02…`.
  It folds pass 1's four findings (`resume_s15_singular_engine_1_review_p1_disposition.md`). Revision 1 is preserved at
  `1585f67c`. The lookups are unchanged (`246b5109…`).
- **Prompt.** `/var/projects/_review_scratch/s15-singular-note1-prompt/note1_review_prompt_p2.md`, sha256 `bc82bd04…`,
  identical for both legs. It is pass 1's prompt plus a pass-2 paragraph, and it asks for a review of the whole
  revision, not only the fold.
- **Legs.** Codex (`gpt-6.1-sol`, xhigh, workspace-write) and Grok (`grok-4.7`, high, workspace sandbox), in parallel,
  2026-11-09. Reports:
  - `/var/projects/_review_scratch/s15-singular-note1-p2-codex/REVIEW.md` (sha256 `4e734151…`);
  - `/var/projects/_review_scratch/s15-singular-note1-p2-grok/REVIEW.md` (sha256 `babe3fa2…`).
- **Verdicts.** Both CLEAR, with no findings.

## What the legs established

- **Facts.** Both legs re-derived the measured account from the guard records themselves:
  - the 56 GiB cap and the start time;
  - 26,708 resource samples, from 03:42:21Z to 18:32:56Z;
  - a peak of 11.23 GiB, with every out-of-memory counter at zero;
  - the last output line and its time;
  - 144 objects in scope, with only the first reached;
  - both outcome files absent;
  - the reconciliation log, the kept stale registry, and the current empty registry.

  Neither sandbox reached the systemd user bus. The unit's state after the stop is the orchestrator's lookup L7.
- **Endings and claims.**
  - A trial stop is allowed only on what the attempt printed or computed. It is neither an ending for its object nor a
    label.
  - A stop without a full run gives an unfinished-build report, never the final report or the final transcript.
  - The full run keeps the directive's ending and final-transcript rules, including the alone, largest-cap rerun after
    an out-of-memory end at a shared cap.
  - No route closes an object other than as the directive allows.
- **Construction.** The note fixes no construction, truncation, expansion, condition or expected output.
- **Operation.**
  - Pooled admission lets overlapping jobs share the machine.
  - Stopping a unit while its launcher is alive lets the guard record the outcome and release its own reservation.
  - Waiting without repeated model calls, inside an open session, is compatible with `AGENTS.md`.

## Disposition

Nothing outstanding changes what the builder computes, what may be claimed from its output, or whether the builder
can follow the note (G4). Revision 2 is **accepted** as resume note 1 for the S15 item-4 Singular engine build.

## Next

The user resumes the builder's Codex session with this note. When the builder's report is in:
- if it is the directive's final report, the build review follows, with a fresh Claude agent and Grok;
- if it is an unfinished-build report, it goes to the user.
