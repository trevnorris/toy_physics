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
   ⛔⛔ **`_scratch/` IS GITIGNORED (`.gitignore:96`), so "commit first" silently does NOTHING there.**
   `git add` refuses the path, so the baseline step **fails open** — you believe you committed and nothing
   happened. ⚠ **Measured twice:** a SymPy audit repaired with no baseline, and a directive revision that
   overwrote its predecessor irrecoverably. ⛔ Neither was carelessness about the rule; both were the rule
   being a **no-op rather than an error** in that location.
   ⭐⭐ **THE FIX IS WHERE THINGS LIVE, ⛔ not a backup ritual** (user, 2026-08-03):
   *"`_scratch` is for things we don't care about … If it's that important it shouldn't go into `_scratch`."*

   | goes in `_scratch/` (disposable, ignored) | goes in a **TRACKED** directory |
   |---|---|
   | raw build/review transcripts, review prompts, fold notes | ⭐ **build directives** → `research/pde_ledger_v3/directives/` |

   ⇒ ⛔ **Never author a build directive under `_scratch/`.** It is the one artifact both engines share and
   the only one whose errors defeat dual-engine by construction — it belongs under version control, where
   the baseline rule actually works.
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

## ⭐⭐ STANDING CHECK before any directive ships — static or instantaneous?

Run this on every specification, and state the answer in the directive rather than leaving it implicit:

> **What timescale or rate did I just send to `0` or `∞`, and what would it have governed?**

⛔ **If the removed rate is close to what the step is trying to determine, the limit is illegitimate** —
it answers the question by assumption, and ⚠ **a closed-form result looks equally healthy either way**, so
no downstream gate will catch it.

⚠ **It never looks like removing time.** It looks like a simple constitutive law or a standard matching
condition. Observed disguises: a **memoryless** closure (removes a relaxation time); a **sharp interface**
(removes the wall's width channel); **adiabatic elimination** (removes a DOF's inertia); **linearising
about rest** when the background flows (removes convective terms); a **frozen wall**.

⭐ **Prefer generalising to taking the limit when it is cheap** — carrying one relaxation time keeps the
algebra rational and makes the instantaneous case a **reportable limit** instead of a premise.
⚠ An honestly-recorded freeze still costs: a prior step recorded *"computed with the wall width FROZEN"*
and undoing it is the entire reason its successor exists. ⇒ Recording a freeze is the minimum, ⛔ not the
fix. ⇒ `[[feedback-static-or-instantaneous-check]]`.

## ⭐⭐ SPEC AUTHORING — equations, small surface, and changing the author

⭐⭐ **Write every binding kinematic relation as an EQUATION, ⛔ never as prose.** ⚠ **Measured:** the same
coupling term fell out of a specification **four times**, every time it was **described in a sentence**.

⭐ **Keep the spec surface small.** A directive spanning many tasks and several distinct physics subsystems
will yield a finding **every round, indefinitely**. ⇒ **Split until each piece can close.**

⭐ **Prefer a re-runnable script as the derivation record** over prose an engine must interpret into algebra.

⭐⭐ **If successive revisions keep breeding defects in the material just changed, CHANGE THE AUTHOR** rather
than rewriting harder; composition still requires that **whatever writes does not review**.
⇒ `[[feedback-spec-authoring-discipline]]`.

## Finding Disposition

Before every repair or rebuild, apply the physics filter. A finding is not a mandate. Act only when it
catches a way the physics could be wrong; `Recorded, not acted on` is a complete disposition.

⛔⛔ **A FINDING ABOUT A CHECK IS NOT A FINDING ABOUT THE PHYSICS.** ⭐ **There is a dedicated red-team
phase for hardening each sector** (user, 2026-08-03) — so a weak or unfailable internal check is a
**red-team item**, ⛔ never a step blocker, and hardening it per-step is **duplicated work**.

⭐ **The question to ask when a review reports weak verification:** *do I already have better evidence for
this value than the script's own opinion of it?* If two independent derivations agree with it, **yes** —
⇒ record the limit and move on.

⚠⚠ **Measured 2026-08-03:** legs found a blind `.wl` with **correct values everywhere** but most
conclusions behind unfailable checks. The repair pass that followed added **183 lines carrying no
physics**, four banned checks-verifying-checks, and left **19 of 23 tags still unfailable** — because the
directive mandated a check reading *"the same expression that is emitted"*, which is `x === x`.

⭐⭐ **When tempted to harden one engine's internal checks, BUILD THE OTHER ENGINE INSTEAD.** The
architecture's answer to *"is this value right?"* is the blind second engine that can **disagree** — ⛔ not
a longer `allChecks` list.
