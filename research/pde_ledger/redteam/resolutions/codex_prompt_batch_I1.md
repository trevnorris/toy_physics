# Codex — answer batch I.1 paper-vs-script questions

You are being asked to act as the math authority on the moving-throat PDE project. The v2 paper-grounded red-team audit pass on batch I.1 (stages 001–012) surfaced 7 questions where the scripts and the paper disagree, and the user wants your recommendations.

## What's going on

The PDE ledger project is in a script-verification phase ("red-team audit pipeline"). Each stage NNN has:

- A paper card at `paper/stages/stage_NNN.tex`
- Optional source notes at `notes/stages/moving_throat_pde_stage{NNN}_*.md`
- A SymPy audit script at `scripts/moving_throat_pde_stage{NNN}_*.py`
- A Mathematica audit script at `mathematica/moving_throat_pde_stage{NNN}_*.wl`

The v2 audit reads the paper FIRST, builds a model of what the stage is supposed to prove, then checks whether the scripts faithfully verify that. When the script's assertion and the paper's claim disagree, that's a `paper_misalignment` finding routed to the user for resolution — and you're the user's math advisor for these decisions.

Per-stage red-team audit reports for the 12 stages in batch I.1 live at `redteam/reports/stage_001.md` through `redteam/reports/stage_012.md`. The full directives (which include the questions you'll be answering, plus the auditor's reasoning) are at `redteam/directives/stage_NNN.md`. You don't have to read all of these, but they're available if you want the full context for a question.

## Your task

Open this questions file:

  `redteam/resolutions/batch_I1_paper_alignment.md`

There are 7 questions, each with three options (a/b/c) and an empty `## Recommendation` block. For each question:

1. Read the question's "Files to read for context" section and open the relevant files.
2. Decide whether option (a), (b), or (c) is correct, based on what the math actually requires.
3. Fill in the `## Recommendation` block with:

   ```yaml
   direction: a       # or b, c, skip
   rationale: |
     1-3 sentences explaining the choice. Cite specific file:line references where helpful.
   ```

4. When you're done with all 7, replace the final summary comment with `summary: X answered, Y skipped`.

## Important rules

- **Only edit the questions file** at `redteam/resolutions/batch_I1_paper_alignment.md`. Do NOT touch any script, paper, or notes file. Recommendations are applied by a separate workflow after the user reviews them.
- **If you do not have a high-confidence answer for a question, write `direction: skip`** with a one-line note explaining why (e.g., "need to check upstream convention in stage X first" or "depends on a physical-interpretation question I'd want to confirm with the author"). Do NOT guess. The user will review skipped questions manually.
- **Read upstream/downstream stages if needed.** Some questions depend on conventions set elsewhere (e.g., Q1's modal-wall PDE sign affects stages 002 and 003 which reuse the PDE form; Q6 and Q7 reference downstream stages 011, 012, 013–014, 022, 023). Open the relevant paper cards / scripts as needed.
- **The earlier red-team audits already verified the algebra holds together internally** — these questions are specifically about claims that the scripts make vs. claims the paper makes, where one side may be wrong relative to the other or where the paper card may need to grow / shrink.
- **Don't edit scripts or paper.** Your recommendations will be applied by a separate workflow after the user reviews them.
- **Iterate as needed.** If reading a stage's notes reveals you need to check a sibling stage, do it. The end goal is high-confidence recommendations or honest "skip" markers.

## Project memory worth knowing

- The Mathematica `.wl` files live in `mathematica/`, NOT in `scripts/`. SymPy `.py` files live in `scripts/`.
- The project framing is a "toy superfluid analog" aiming to unify GR-like and EM-like derivations. In private/internal context (which this is), engaging with the "unification" language directly is fine. (Public paper framing uses strict toy-analog language only.)
- The project is moving toward an end-to-end second pass, so accuracy of these recommendations matters — the first pass result is what the second pass will check against.

## When done

Write the summary line at the bottom of the questions file (`summary: X answered, Y skipped`) and exit. Do not modify any other file.
