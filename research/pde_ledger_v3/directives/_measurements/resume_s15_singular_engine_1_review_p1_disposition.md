# Resume note 1 of the S15 Singular engine build, review pass 1: disposition

## The gate record

- **Artifact.** `research/pde_ledger_v3/directives/resume_s15_singular_engine_1.md`, revision 1, sha256 `89f858d9…`.
  It was preserved exactly as reviewed at `1585f67c`, with its lookups (`246b5109…`).
- **Role and author.** The orchestrator wrote it. It amends how the builder of the S15 item-4 Singular engine directive
  (accepted `d4cd6d09`) runs its remaining attempts, and it records the user's 2026-11-09 decision. It is a governing
  directive, so it is reviewed until clear.
- **Prompt.** `/var/projects/_review_scratch/s15-singular-note1-prompt/note1_review_prompt_p1.md`, sha256 `97eafdf5…`,
  identical for both legs.
- **Legs.** Codex (`gpt-6.1-sol`, xhigh, workspace-write) and Grok (`grok-4.7`, high, workspace sandbox), in parallel,
  2026-11-09. Reports:
  - `/var/projects/_review_scratch/s15-singular-note1-p1-codex/REVIEW.md` (sha256 `f3bbf97c…`);
  - `/var/projects/_review_scratch/s15-singular-note1-p1-grok/REVIEW.md` (sha256 `dba6e28d…`).
- **Verdicts.** Both NEEDS REVISION.

Both legs checked the note's measured facts against the guard records, the inventory, the registry and the
reconciliation log, and found them correct. Neither found a construction, truncation, expansion or expected value
fixed by the note. Neither sandbox reached the systemd user bus, so neither re-ran lookup L7 (the unit's state).

## Findings and dispositions

Each finding was verified against the note, the directive, `AGENTS.md`, `README.md` and the guard script before it
was folded.

1. **A failed trial could close the build** (Codex 1, Grok 1). **Confirmed.**
   - The note let the builder stop after a trial and write "the directive's final report". That would leave objects in
     scope with neither their objects nor a label, against the directive's "In the final transcript, every object in
     scope ends with its item-3 object or a label".
   - The trial gate said "returns", while the directive has four ways for an attempt to end.
   - **Fold:**
     - the trial gate is "ends in one of the directive's four ways";
     - a stop without a full run gives a separate unfinished-build report, which closes no object and labels none,
       and whose trial output is not presented as the final transcript;
     - the directive's rules for the final transcript stand unchanged for the full run.
2. **The builder's stop added an ending, against "No time limits"** (Codex 2). **Confirmed.**
   - `AGENTS.md` lines 45–46: "Add no wall-clock, CPU-time or inactivity deadline unless the user asks for one for a
     specific run."
   - **Fold:**
     - the builder may stop a trial attempt only on what it printed or computed, never for elapsed time or a period
       without output;
     - a stopped trial attempt is not an ending for its object, and gives it no label;
     - the run record gives the evidence it was stopped on.
3. **The watcher ban conflicted with the standing waiting policy** (Codex 3). **Confirmed.**
   - `AGENTS.md` line 69 says "Don't poll the model while waiting", and `README.md` lines 151–162 describe a local
     completion watcher.
   - **Fold:**
     - the builder waits without repeated model calls;
     - its session stays open until its jobs have ended and its report is written;
     - whatever it uses to wait must not depend on a message reaching a session that has ended. That is the measured
       failure of its first session.
4. **No way to stop a job while its launcher waits** (Grok 2). **Confirmed.**
   - The guard launches with `systemd-run --wait` (`scripts/s11c_guarded_run.py:504`). `release_pool` refuses a release
     whose owner does not match (line 332).
   - **Fold:**
     - each job's launcher stays alive until the job ends, whatever happens to the builder's turn or session;
     - a job is stopped by stopping its unit while its launcher is alive, so the launcher records the outcome and
       releases the reservation;
     - the guard registry is never edited.

## Next

Revision 2 goes to pass 2 with the same legs, as a review of the whole revision.
