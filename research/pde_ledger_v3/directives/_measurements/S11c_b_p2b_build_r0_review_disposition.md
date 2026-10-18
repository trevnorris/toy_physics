# S11c-b P2b build, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S11c_b_p2b_build_r0_review_lookups.md` beside this file
(generator `_scratch/backlog_runs/gen_p2b_r0_lookups.sh`).

## The round
- **Artifacts.** The three S11c-b comparator scripts `S11c_b_handcoded_comparison.py`,
  `S11c_b_adjudicated_comparison.py` and `S11c_b_cross_engine_comparator.py`, and the builder's report
  `S11c_b_p2b_build_report.md`. Eight files are in the baseline `_scratch/backlog_runs/p2b_build_review_baseline_r0.sha256`,
  all OK at review.
- **Contract.** The decision list `S11c_b_p2b_gamma_bridge_directive.md` (v2, committed at `69bcc18d`), items B1–B6
  and the controls in its last paragraph.
- **Builder.** Codex `gpt-6-astra` at high effort, in a fresh session, in the build checkout
  `/var/projects/toy_physics-p2b`: an export of `da3f7947` without history, committed there as `a5da2ee`, with the
  annexed `.out` transcripts copied in as plain files. The build changed only the three scripts (435 insertions,
  52 deletions) and added the report.
- **Prompt and legs.** The prompt is `_scratch/backlog_runs/p2b_build_review_prompt_r0.md` (sha256 `ce922f61…`). The
  legs were a fresh Claude agent (opus) and Grok (`grok-4.7`, effort high), as the authorship table gives for a
  Codex-written script. Their reports are `_scratch/backlog_runs/p2b_build_review_r0_claude.md` and
  `…_r0_grok.txt`, with evidence directories beside them. Grok returned **CLEAR**; Claude returned **NEEDS REVISION**.
- **Prompt defect (mine).** The prompt named the build checkout's guard with a log directory in the main repository.
  That guard accepts only log directories under its own `_scratch/s11c/`. The Claude leg used the main repository's
  guard instead, which is the same file byte for byte. No verdict depends on it.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude F1 | **A §3a disagreement, or a count other than the typed ones, stops the whole comparator before any family is compared.** The certificate raises on a pairing degree other than one, on typed source, branch and overall counts, and on any nonzero certificate residual or false lock (`handcoded_comparison.py` 594–615, 655–656). Those raises run before the family loop (`adjudicated_comparison.py` 1427–1429, and in the handcoded `main`). The leg added one term to one Wolfram invariant: the adjudicated default run exited 2 with no line outside the P2b block. The same held for a changed Wolfram action. | **Confirmed.** Each value is printed before its raise, but the raise ends the run, so a genuine §3a disagreement hides every other family's comparison. That changes what the comparator shows when the engines disagree, which is the case it exists for (CLAUDE.md M1). The decision list's B6 says "no assert on the residual", and the repository's own test `test_disagreement_is_measurement_exit_zero` treats a printed disagreement as a measurement. The typed counts are also expected values sitting in guards. Grok's checks ran on agreeing input, where no raise fires. |
| 2 | Claude F2 | **The build's controls are a separate re-implementation, and the fold control's lock is false by its own construction.** `p2b_action_controls` (722–757) never calls the production certificate or the bridge's apply. The fold control sets the factor to one and then computes `factor_lock` as that factor's equality with the source factor (lines 740, 751). The production certificate's lock is a different predicate (line 645). | **Confirmed.** A control that sets its own outcome tests nothing (CLAUDE.md E2; corollary 3 of the build skill). Both legs ran their own mutations through the production path, and it does respond to the swap, the bare image and the fold. So the bridge is sound on this evidence, but the build's printed controls are not evidence for it. |
| 3 | Claude F3 | **The B4 fixture and the classification coverage check are built from the tables they test.** The fixture loops over `LIVE_BACKGROUND_APPLIED`, which is the comparator's own `LIVE_BACKGROUND_HEADS` (lines 198, 664). Coverage compares `APPLIED_HEAD_CLASSES` with `JET_APPLIED`, and both are unions of the same three sets (lines 197–206, 689). The leg removed `widthBackground` from the table: both `--p2b` runs exited 0, coverage printed `true`, and the nine `widthBackground` tags disappeared. | **Confirmed.** A check cannot audit its own input. The decision list asks for a fixture that shows base and shifted arguments stay distinct for `widthBackground`, `modulusBackground` and `eWBackground`. As built, dropping a head from the table drops it from the test as well. |
| 4 | Grok note, Claude note | **The builder's saved runs predate its last edits.** The handcoded and adjudicated scripts were last written at 16:05:23. The builder's guarded runs ended between 16:02 and 16:04. The Claude leg's reruns of the final code differ from the builder's output (900 added pairing lines in the adjudicated run, 6 in the comparator). | **Confirmed.** The report cites output from code that is not the final code. Both legs reran the final code themselves, so their verdicts do not depend on the stale runs. |

**Also recorded (not build defects).**
- The adjudicated default run reached `MemoryError` at 6 GiB and at 12 GiB, both at the base commit and in the build.
  That failure predates this build.
- With the bridge in place, the Claude leg sees `MU_THETA_OPERATOR` terms that remain different between the
  engines, for both sources and every key. This is a cross-engine result for the residual chain (P2a and
  adjudication), not a defect in P2b. It is recorded, not interpreted, here.

**Not clear.** Findings 1–3 change what the comparator shows or what may be claimed from the build's controls.
Finding 4 changes what the report may cite.

## What must be true after the repair
- A §3a disagreement of any kind, including an unpaired coefficient or a count that differs, is printed and stays
  visible in the comparator's output, and the comparator goes on to compare every family. Only malformed input
  stops a run. No expected count is typed into a guard.
- Each control in the decision list's last paragraph mutates the production bridge and passes through the same
  certificate computation the production run uses, so what it prints is that certificate's own residuals and lock.
  No control sets its own outcome.
- The B4 fixture names the three live-background heads independently of the table it tests, and classification
  coverage is checked against a census that is not built from the classified sets. Removing a head from a table
  changes what is printed.
- Every run the report cites is a run of the final code.

This is the builder's first repair, so the same builder session continues (G4 does not call for a new author). The
reviewed scripts and report are preserved in the main repository before the repair overwrites them.
