# Codex — answer batch I.2 paper-vs-script questions

You are being asked to act as the math authority on the moving-throat PDE project. The v2 paper-grounded red-team audit pass on batch I.2 (stages 013-023) surfaced 6 questions where the scripts and the paper disagree, and the user wants your recommendations.

## What's going on

The PDE ledger project is in a script-verification phase ("red-team audit pipeline"). Each stage NNN has:

- A paper card at `paper/stages/stage_NNN.tex`
- Optional source notes at `notes/stages/moving_throat_pde_stage{NNN}_*.md` (for the EM-projected linear stages 013-020, no per-stage notes exist — the paper card + appendix row + the original `notes/em_projected/` derivation scripts are the source)
- A SymPy audit script at `scripts/moving_throat_pde_stage{NNN}_*.py`
- A Mathematica audit script at `mathematica/moving_throat_pde_stage{NNN}_*.wl`

The v2 audit reads the paper FIRST, builds a model of what the stage is supposed to prove, then checks whether the scripts faithfully verify that. When the script's assertion and the paper's claim disagree, that's a `paper_misalignment` finding routed to the user for resolution — and you're the user's math advisor for these decisions.

Per-stage red-team audit reports for the 11 stages in batch I.2 live at `redteam/reports/stage_013.md` through `redteam/reports/stage_023.md`. The full directives (which include the questions you'll be answering, plus the auditor's reasoning) are at `redteam/directives/stage_NNN.md`. You don't have to read all of these, but they're available if you want the full context for a question.

## Batch I.2 context

Batch I.2 is the back half of the linear projected-EM core range (stages 004-021 plus 022-023 grouped-bundle bookends). Stages 013-020 are file-for-file SymPy migrations of the original `notes/em_projected/step_NN_*` derivation scripts; the paper cards under `paper/stages/` were written later as compact summaries. The audit pattern is dominated by `paper_missing_script_claim`: the scripts verify more than the paper cards declare. This is the same pattern as batch I.1.

Three exceptions worth noting:
- Stage 015 has a missing notes file (sympy docstring references `step_13_parent_throat_action_master_notes.md` which does not exist in the repository).
- Stage 018 similarly references nonexistent `step_16_parent_throat_action_bundle_master_notes.md`.
- Stage 021 is the inverse direction: the paper Output enumerates three exports, but the scripts only assert two — the third (composed wall-level odd quadrupole coefficient) is computed but not asserted against the paper's specific closed value.

## Your task

Open this questions file:

  `redteam/resolutions/batch_I2_paper_alignment.md`

There are 6 questions. For each question:

1. Read the question's file:line references and open the relevant paper card + script files.
2. Decide which option (a, b, c, or skip) is correct, based on what the math actually requires.
3. Fill in the `## Recommendation` block with:

   ```yaml
   direction: a       # or b, c, skip
   rationale: |
     1-3 sentences explaining the choice. Cite specific file:line references where helpful.
   ```

When you're done with all 6, append a `summary:` line at the end of the questions file noting `X answered, Y skipped`.

## Important rules

- **Only edit the questions file** at `redteam/resolutions/batch_I2_paper_alignment.md`. Do NOT touch any script, paper, or notes file. Recommendations are applied by a separate workflow after the user reviews them.
- **If you do not have a high-confidence answer for a question, write `direction: skip`** with a one-line note explaining why. Do NOT guess. The user will review skipped questions manually.
- **Read upstream/downstream stages if needed.** Some questions depend on conventions set elsewhere — Q1 (stage 013's `deltaP_2/deltaP_4` and sieve) depends on what stage 014 formally owns; Q3 (stage 015's wall-only/Y20/grouped blocks) needs you to check whether stages 016, 017, 020, 022, or 023 already own the same content; Q4 (stage 018's extra families) might belong in stage 021. Open the relevant paper cards as needed.
- **The earlier red-team audits already verified the algebra holds together internally** — these questions are specifically about claims that the scripts make vs. claims the paper makes, where one side may be wrong relative to the other or where the paper card may need to grow / shrink.
- **Don't edit scripts or paper.** Your recommendations will be applied by a separate Codex apply workflow after the user reviews them.
- **Iterate as needed.** If reading a stage's notes reveals you need to check a sibling stage, do it. The end goal is high-confidence recommendations or honest "skip" markers.

## Project memory worth knowing

- The Mathematica `.wl` files live in `mathematica/`, NOT in `scripts/`. SymPy `.py` files live in `scripts/`.
- The project framing is a "toy superfluid analog" aiming to unify GR-like and EM-like derivations. In private/internal context (which this is), engaging with the "unification" language directly is fine. (Public paper framing uses strict toy-analog language only.)
- The project is moving toward an end-to-end second pass, so accuracy of these recommendations matters — the first pass result is what the second pass will check against.
- The same pattern in batch I.1 was resolved as: Q1=a, Q2=b, Q3=c, Q4=a, Q5=c, Q6=a, Q7=c — a mix of expand-paper, trim-script, and add-acknowledgement directions based on what the math actually owned. Use that as a reference for the kind of judgment expected, but answer I.2 on its own merits.

## When done

Append the summary line at the bottom of the questions file (`summary: X answered, Y skipped`) and exit. Do not modify any other file.
