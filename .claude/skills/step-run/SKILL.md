---
name: step-run
description: Execute one PDE-ledger derivation step in the proven order from pre-registration through directive review, blind Mathematica and SymPy builds, quarantine, gates, step record, TeX card, and independent review legs. Use when running a whole ledger step rather than only building or reviewing one artifact.
allowed-tools: Bash, Read, Edit, Write, Agent
user_invocable: true
---

# Run One Ledger Step

Invoke as `/step-run <step-id>`. Follow this sequence exactly; do not substitute a check for a review.

## Runbook

1. **Walk the step with the user.** State what is being derived, what is available, what is missing, and
   where trouble is expected. Then take one move at a time with reasoning before results. Flag every
   identification before making it and state what would make it wrong; never pre-derive and present.
2. **Pre-register predictions.** Write expected script outputs to a tracked path, commit that path so
   the commit timestamp establishes priority, record the blob/commit, then move the predictions out of
   the tree. A do-not-read instruction does not survive a grep; the file must be absent.
3. **Review the build directives before any build.** Run `.claude/skills/review-legs/SKILL.md` on the
   directive packet with both directives visible as one review artifact. ⭐ The directives are
   **orchestrator-written**, so the two legs are **Codex + Grok** — ⛔ not a fresh Claude agent. Require
   the unique comparison: do the two directives specify the same physics? The directive is shared by both
   engines, so one error there can make dual-engine agreement certify wrong physics.

   ⛔⛔ **THE DO-NOT-READ LIST MEANS EVERY FILE THAT CONTAINS AN ANSWER — ⛔ NOT EVERY FILE NAMED LIKE ONE.**
   ⚠ **Measured 2026-08-03:** the list was built by *filename* (`*PREREGISTERED*`) and so omitted
   `V3_STEP_PLAN.md` and `steps/`, which the **previous step's own directive had banned**. Four separate
   loci in the plan stated expected results outright — a long-wavelength exponent, an expected
   transverse-vs-longitudinal asymmetry, a projection identity's form **plus the exact file and lines
   where it was computed**, and a validity condition verbatim. ⇒ both builders could grep, both read the
   same paragraph, and dual-engine would certify **the plan's expectation rather than a derivation**.
   ⭐ **Before writing the list, grep the live docs for the step's own tag names and expected objects, and
   quarantine what hits.** ⚠ Two of those four loci were written **that same session, by the orchestrator,
   while it was being careful about the directive** — the plan is reachable, and recording a result's
   *shape and location* leaks it just as surely as recording its value.
4. **Build the blind `.wl` first, before any `.py` exists.** Its directive carries the action and none of
   the results. Invoke `.claude/skills/build/SKILL.md`; it leak-gates before launch and starts both review
   legs itself before you inspect the artifact.
5. **Arbiter re-run.** Run the `.wl` yourself and compare its literal outputs with the pre-registration.
   Reproduction proves determinism and a match proves agreement with you; neither is review because a
   shared wrong assumption passes both.
6. **Establish a baseline before change.** Commit the generated artifact before any repair, destructive
   edit, or quarantine. A repair to an untracked file leaves no baseline and no reviewable blob.
7. **Quarantine, then build SymPy.** Move the `.wl` out of the working tree, build the SymPy audit and
   any registry insertion through `/build`, then restore the `.wl`. Verify it is byte-identical to its
   committed blob. Reviewers may read that blob with `git show <sha>:<path>` but must never restore it;
   this keeps the SymPy builder blind while `.wl` review continues.
   ⛔⛔ **QUARANTINE THE WHOLE ANSWER-BEARING SET, ⛔ not just the `.wl`** — the build directives, the
   step's `_scratch/<step>_*` files, and any raw build transcript. ⚠ `_scratch/` accumulates Codex
   transcripts holding a prior engine's **complete tag values verbatim**, and they are reachable by a
   builder that obeys every stated instruction. ⇒ `.claude/skills/review-legs/SKILL.md`
   § BLINDNESS IS ENFORCED BY ABSENCE.
8. **Run all existing gates.** Run acceptance, dim, able-to-fail, and pytest gates. For a new
   discrete row, also verify that the continuous payload is unchanged.
   ⛔⛔ **Any NEW relation also needs its own algebra checked, because no gate does it.** Add an
   assertion to the step's audit that substitutes the root the script **derived** into the registry
   residual that records it, and asserts it vanishes. Ablate it to prove it fails.
   ⚠ **Measured, not hypothetical:** `R5` was rewritten to `c_L = c_γ` — the exact claim its step
   existed to settle — and **all five gates stayed green.** Homogeneity was blind (the two moduli share
   a dimension), the acceptance fixture was blind (the designated output stays a fresh variable, so the
   constraint count is unchanged), and the Mathematica engine cannot help **by design** because it must
   not read the registry. ⇒ `research/pde_ledger_v3/DEFECT_REGISTER.md#f-r5`.
9. **Write the step record as orchestrator.** Record the walk yourself; Codex was not present for it.
10. **Have Codex write the TeX card from the record.** Invoke `/build` for the card so its own fresh-agent
    and Grok legs launch on sight, before you read it. Builder and reviewer remain separate for `.tex`.

## Finding Disposition

Before every repair or rebuild, apply the physics filter. A finding is not a mandate. Act only when it
catches a way the physics could be wrong; `Recorded, not acted on` is a complete disposition.
