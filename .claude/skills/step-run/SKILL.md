---
name: step-run
description: Execute one PDE-ledger derivation step in the proven order from pre-registration through directive review, independently-constructed Mathematica and SymPy builds, gates, step record, TeX card, and independent review legs. Scripts must print computed objects and never state conclusions. Use when running a whole ledger step rather than only building or reviewing one artifact.
allowed-tools: Bash, Read, Edit, Write, Agent
user_invocable: true
---

# Run One Ledger Step

Invoke as `/step-run <step-id>`. Follow this sequence exactly; do not substitute a check for a review.

## ⭐⭐⭐ THE PROVEN ORDER — measured end-to-end on S9, 2026-08-04

⚠ **S9 took 4 repair rounds on the `.wl` and 3 on the `.py`.** ⛔ Almost none of that was rework: each
round surfaced a **new defect class** that the previous round's mechanisms structurally could not see.
⭐ A later step inherits all of them at once and should converge in one or two rounds. ⚠ **If it does not,
that is informative** — it means the step has a defect class S9 did not.

| # | do this | the trap it exists for |
|---|---|---|
| 1 | **Author the build directive.** Supply every equation; withhold **only** an acceptance criterion referencing an expected value. | under-specification has cost more than contamination |
| 2 | ⭐ **LEAK-GATE it.** `rg -F` every expected value, exponent, ratio, count, sign. | ⚠ **fired 4× on S9, 3× inside sentences FORBIDDING the answer**, one introduced *by the repair for another* |
| 3 | ⭐⭐ **REVIEW THE DIRECTIVE — Codex + Grok — BEFORE any build.** | the ONE artifact both engines share; an error there makes dual-engine certify wrong physics |
| 4 | **Fold, and ⛔ STOP a running build if the QUESTION is wrong.** | ⚠ S9: a build was killed mid-run over a false-negative test; the alternative was interpreting output from a bad setup |
| 5 | **Build engine 1.** ⛔ Launch in its **own tool call**, nothing chained. Prove the prompt non-empty *before*, the log >500 bytes *after*. | ⚠ 2 hangs, ~25 min; a hung Codex looks busy and **never notifies** |
| 6 | ⛔⛔ **COMMIT THE ARTIFACT** before anything repairs or ablates it. | an uncommitted baseline is what gets destroyed; a review needs a fixed target |
| 7 | ⭐ **Launch review legs ON SIGHT**, before reading results. ⚠ Serialize if both ablate Mathematica (2 seats). | a self-check discharges the felt need for an independent one |
| 8 | **Run the engine yourself and READ THE OUTPUT.** | the output is the result; ⛔ do not pre-register expected values |
| 9 | **Build engine 2** from the *same question list*, ⛔ barred from engine 1 and the registry. | a transcription agrees vacuously |
| 10 | **Run the harness over both.** ⛔ Verify it by **ablating the harness**. | a checker that mis-parses into agreement manufactures confidence |
| 11 | **Fold findings; repair; re-verify each fix against your own computation.** | ⛔ never trust a builder's report |
| 12 | **YOU write the step record.** Codex was not there for the walk. | — |
| 13 | **Codex writes the TeX card**, with its own two legs. | ⚠ a card can assert the **absence of things that now exist** |

⭐⭐ **THE RESULT THAT JUSTIFIES THE WHOLE APPARATUS, measured on S9 — ⚠ and note the two rows are
DIFFERENT KINDS of claim:**

| defect | caught by | missed by | kind |
|---|---|---|---|
| a **wrong dimension** | ⭐ the **second engine** | two review legs **and** a full ablation suite — engine 1 computed it *consistently*, ablated *correctly*, and produced a number that **looked like a signal** | ⚠ **EMPIRICAL.** A leg that happened to dimension that expression by hand *would* have caught it. ⛔ Not a structural limit — an observed miss |
| a **wrong homogeneity test** | ⭐ reviewers **deriving from scratch** | ⛔ **cross-engine comparison** | ⭐⭐ **STRUCTURAL.** The defect came from the **shared directive**, so both engines computed it the same wrong way and **agreed**. ⛔ Comparison *cannot* see this class |

⇒ ⛔ **Never treat cross-engine agreement as sufficient** — that one is provable.
⇒ ⚠ **Do not treat review as sufficient either** — that one is a measured miss, not a theorem, ⛔ but a
mechanism that missed a live defect twice is not something to lean on alone.

## Runbook

1. **Walk the step with the user.** State what is being derived, what is available, what is missing, and
   where trouble is expected. Then take one move at a time with reasoning before results. Flag every
   identification before making it and state what would make it wrong; never pre-derive and present.
2. ⛔⛔ **DO NOT PRE-REGISTER EXPECTED VALUES — RETRACTED 2026-08-04, and the reason is the whole method.**
   > ⭐⭐ **The entire point of writing the script is to GET the answer. Your job is to make sure the
   > correct QUESTION is asked.** (user, 2026-08-04)

   ⚠ **Measured:** a session wrote a hand-derived prediction table, had a review leg "independently
   derive" and agree with it, and was one tool call from writing a third check to agree again — with
   **zero scripts rebuilt** and most of a session gone. ⛔ **A committed table of expected values IS an
   acceptance criterion referencing an expected value** — reachable in the tree by a builder that
   iterates to exit 0, which converts a genuine disagreement into silent confirmation.
   ⇒ ⭐ **Spend the effort on the QUESTION instead:** is every supplied relation an equation, is every
   requested object computable from what is supplied, can each control genuinely fail, and is any
   asked-for object **ill-posed** or **tautological**?
   ⭐ **If working the setup exposes a defect, record THE DEFECT — ⛔ never the value that exposed it.**
   ⇒ `[[feedback-question-quality-not-answer-confirmation]]`.
3. **Review the build directives before any build.** Run `.claude/skills/review-legs/SKILL.md` on the
   directive packet with both directives visible as one review artifact. ⭐ The directives are
   **orchestrator-written**, so the two legs are **Codex + Grok** — ⛔ not a fresh Claude agent. Require
   the unique comparison: do the two directives specify the same physics? The directive is shared by both
   engines, so one error there can make dual-engine agreement certify wrong physics.

   ⛔⛔ **THERE IS NO DO-NOT-READ LIST — ⛔ CUT 2026-08-12.** ⚠ This step used to say to build one, and
   `CLAUDE.md` rule 12 had already cut it; the residue outlived §106's own rewrite in this same file.
   ⚠⚠ **The measurement behind it stands, and it argues the opposite of what it used to prescribe.**
   Measured 2026-08-03: the list was built by *filename* (`*PREREGISTERED*`) and so omitted
   `V3_STEP_PLAN.md` and `steps/`, which the **previous step's own directive had banned**. Four separate
   loci in the plan stated expected results outright — a long-wavelength exponent, an expected
   transverse-vs-longitudinal asymmetry, a projection identity's form **plus the exact file and lines
   where it was computed**, and a validity condition verbatim. ⇒ both builders could grep, both read the
   same paragraph, and dual-engine would certify **the plan's expectation rather than a derivation**.
   ⇒ ⭐⭐ **A list that must enumerate every file containing an answer is a denylist against a corpus that
   grows faster than the list.** ⛔ Do not grep-and-quarantine what hits.
   ⭐ **What the measurement actually obliges: STOP WRITING EXPECTED RESULTS INTO LIVE DOCS.** ⚠ Two of
   those four loci were written **that same session, by the orchestrator, while it was being careful about
   the directive** — ⚠ recording a result's *shape and location* leaks it as surely as its value.
   ⭐ The control is `CLAUDE.md` rule 5 at the moment of writing, ⛔ not a bar applied afterwards.
4. **Build the blind `.wl` first, before any `.py` exists.** Its directive carries the action and none of
   the results. Invoke `.claude/skills/build/SKILL.md`; it leak-gates before launch and starts both review
   legs itself before you inspect the artifact.
5. **Arbiter re-run.** Run the `.wl` yourself and compare its literal outputs with the pre-registration.
   Reproduction proves determinism and a match proves agreement with you; neither is review because a
   shared wrong assumption passes both.
6. **Establish a baseline before change.** Commit the generated artifact before any repair or destructive
   edit. A repair to an untracked file leaves no baseline and no reviewable blob.
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
7. **Build SymPy independently of the `.wl`.** ⚠⚠ **REWRITTEN 2026-08-04 — the quarantine ritual is CUT.**
   ⭐ **What survives, and it is the whole point: the two engines must be INDEPENDENTLY CONSTRUCTED so they
   can DISAGREE.** Write the `.wl` first, bar it from the registry, and ⛔ do not let the SymPy build be a
   transcription of it — ⭐ that is about **construction**, ⛔ not about hiding a result.
   ⛔ **CUT: moving the `.wl` out of the tree, quarantining directives, `_scratch` transcripts and the
   answer-bearing set, and the byte-identical-restore check.** ⚠ Those defended against **anchoring**;
   the measured failure was **absence of computation**, which the three clauses kill structurally and
   quarantine never touched. ⇒ `research/pde_ledger_v3/REBUILD_HANDOFF.md`.
7b. ⭐⭐ **BOTH ENGINES, EVERY STEP — and re-check a CONDITION before inheriting its conclusion.**
   ⚠ **Measured 2026-08-04:** a doc said S9 needed no SymPy engine because *"the cone is two lines of
   algebra"*, under a **conditional** rule — *"a second engine earns its place where the algebra is long
   enough that it could genuinely DISAGREE."* The condition had **expired** (26 tags had become 316) and
   the conclusion was inherited anyway. ⇒ ⭐ **A conditional rule needs its CONDITION re-checked, ⛔ not
   its conclusion quoted.**
   ⭐ **What the second engine bought, on the first try:** it caught a **wrong dimension** the `.wl` had
   been reporting — a defect that had survived **two review legs and a full ablation suite**, because the
   `.wl` computed it consistently, ablated correctly, and produced a number that *looked* like a
   meaningful signal. ⛔ Neither review nor ablation can catch that class; only a second engine computing
   the same quantity a different way.
   ⚠ **The converse, and it bounds what cross-engine buys:** a defect **shared** by both engines is
   invisible to comparison. ⇒ the reviewer deriving from scratch is still required.

7c. ⭐⭐ **RUN THE HARNESS ON BOTH OUTPUTS** — `reduction/engine_output_checks.py --config checks_<STEP>.yaml`.
   Four layers: **cross-engine agreement** (target generated at compare time from the other engine),
   **dimensional homogeneity** of every emitted expression (dimensions taken from the engine's own derived
   dimension tags), **control response** (does each tag move under some control), and **tag-set parity**
   across packages.
   ⛔⛔ **The config maps TAG NAMES, never values** — that is what keeps it leak-free, and it is not
   optional.
   ⛔ **A physics finding must EXIT 0**; only operational failure exits non-zero. Otherwise a builder
   iterating to exit 0 can make a disagreement disappear.
   ⛔ **Verify the harness by ABLATING THE HARNESS**, never by reading its self-report — a checker that
   silently mis-parses two expressions into agreement is worse than no checker, because it manufactures
   confidence.
   ⚠ **Its limits, and ⛔ do not quote it as coverage:** `INVARIANT` and `PARITY` gaps are **triage
   lists**, not failures — a quantity can be legitimately invariant or legitimately absent from a control.

7d. ⭐ **DIMENSIONS — two rules that were both learned the hard way:**
   ⭐ **Dimension the WHOLE expression by walking its tree.** ⛔ Reading the exponents of the coefficient
   symbols silently drops every other dimensionful factor — a wavevector, for instance — and reports a
   dimension short by it. ⚠ That produced a **wrong emitted value** that looked like a meaningful signal.
   ⭐ **Count FIELD FACTORS, not derivative atoms.** ⚠ A per-term analysis that sums over `Derivative`
   nodes gives **no contribution** for a bare underived field, so a gap or mass term is silently
   mis-dimensioned — with a clean exit and every gate green.
   ⚠ Homogeneity is **blind to a wrong dimensionless coefficient**. ⇒ it is a layer, ⛔ not the answer.

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

⛔⛔ **AND MAKE THE SCRIPT INCAPABLE OF STATING A CONCLUSION.** Every script directive carries the three
clauses in `.claude/skills/build/SKILL.md`: **print computed objects, never prose · print the residual,
emit BOTH operands and the residual before guarding · interpretation belongs to the step record.**
⚠ Measured 2026-08-04: named tags at named lines in **three independently-built steps** are typed prose
with no CAS object. ⛔ Do not quote a fraction. ⇒ `REBUILD_HANDOFF.md`.

⛔⛔ **DO NOT BLIND THE INPUTS — supply every equation.** ⭐ Withhold exactly one thing: an **acceptance
criterion that references an expected value**, because Codex iterates to exit 0 and will otherwise fix the
computation until it matches. ⇒ The builder's job **ends at compute-and-print**; the diff happens on the
orchestrator's side, where a mismatch is a **finding**, ⛔ not a build failure.

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
