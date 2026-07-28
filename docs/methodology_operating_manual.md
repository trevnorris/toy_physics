# AI-Assisted Research — Operating Manual (drop-in for an agent session)

**What this file is.** A self-contained methodology for doing AI-assisted derivation / modeling / analysis work **where correctness
matters and there is no external peer review**. It is written to be *followed by an AI coding agent* (e.g. a Claude Code session), not
just read by a person. It is portable: it names abstract *roles*, not specific tools — §2 tells you how to map the roles onto whatever
you have.

**How to deploy it.**
- *Human:* drop this file in your repo and **reference it from your `CLAUDE.md`** (or paste §1 + §6 at the start of a session). Tell the
  agent: *"Follow `methodology_operating_manual.md` for this work. State how you're instantiating the roles, then proceed one gate at a
  time."*
- *Agent:* treat §1 as binding for the whole session. Before starting real work, (a) state how you're instantiating the roles (§2) with
  the tools available here, and (b) confirm you will run the loop in §3 one gate at a time with a human gate between gates.

**The one idea.** An LLM will produce a plausible, internally-consistent derivation for *almost any* premise. So AI-assisted results
have no weight unless the pass/fail criteria were fixed **before** the output existed, every check **could have failed**, and the work
was verified by something **independent of what produced it**. Your job is to be that discipline, not just to produce the answer.

---

## 1. Prime directives (non-negotiable — hold these for the whole session)

1. **Separate the maker from the judge.** Whatever writes the code/derivation must **not** be the sole judge of whether it's correct.
   If you are orchestrating, do not also hand-write the computation and then "review" your own work.
2. **Every verdict-bearing check must be able to fail, and must *compute* from inputs.** No result may be a value you typed in, or a
   check that passes by construction. If a check would print the same answer whether or not the claim held, it is worthless. **A clean
   "it all works" is the *suspicious* outcome**, not the reassuring one.
3. **Truth lives in the output artifact, not in prose.** Neither you nor any model *asserting* a fact makes it true. Judge results by
   what the computation emitted (numbers, artifacts, which code path ran) — never by the confidence of the write-up.
4. **Independence is the only real verification.** Verify on a **fresh reviewer with no stake and no hint of the expected answer.** A
   reviewer that watched the work get built normalizes its errors.
5. **Pre-commit, then don't move the goalposts.** Write down pass/fail/"trustworthy-miss" criteria before running. **Never mutate a
   test after seeing its results** (don't retarget, redefine, add a "forgot-a-factor" correction, or switch to a friendlier metric). A
   miss is a *result*, not a prompt to search for a pass.
6. **A negative / "can't-be-done" verdict gets *harder* scrutiny than a positive** — it is the cheapest exit and the easiest way to
   fool yourself, especially when it matches what you already believed.
7. **The human owns scope and the final verdict. When unsure whether a deviation is "a call you can make," it isn't — halt and ask.**
   Do not silently rewrite the process to get unstuck; a later summary claiming "we now skip X" is *suspect*, not authority.

---

## 2. Roles — and how to instantiate them with your tools

The process needs several **independent viewpoints**, and — wherever you can — those viewpoints should be **different LLMs, not just
different sessions of the same one.** Heterogeneity is the point: two instances of one model share the same blind spots, so one
reviewing the other's work normalizes the *same* mistakes; different models (different training, different biases) fail *differently*,
and their disagreements are the highest-signal check you have. A "full and clear view" is precisely what you get from routing the same
artifact past several engines that go wrong in different ways. Adapt the *separation* to your tools, but preserve the
*different-engine* character wherever possible.

| Role | Responsibility |
|---|---|
| **Orchestrator** | Owns the *specification* (the test spec, the review prompts, the status notes), reviews everything, gates, records results. **Writes no computation code** — not even throwaway scratch. |
| **Executor** | Designs, writes, and **runs** the computation; iterates to a clean exit; applies all fixes. Does not pronounce the verdict. |
| **Tertiary reviewer** | *One* fresh, independent pass over each spec. Not iterated. |
| **Clean agents** | Each verification leg on its *own* fresh agent, given the artifact + the standard but no hint of the answer. |
| **Human** | Scope, pre-registration, gate approvals, framing, "too good to be true," final verdict. |

**Reference instantiation (how this project actually splits the work — three *different* LLMs + fresh subagents):**
- **Orchestrator = Claude Code.** Owns the directive/spec, the review prompts, and the status notes; reviews everything; runs the
  arbiter re-run of existing scripts unchanged; spawns the fresh clean-agent subagents for the verification legs; **writes no
  computation code**.
- **Executor = Codex, at its highest reasoning tier.** Designs, writes, and runs all the computation, on two engines (a primary CAS +
  a cross-check CAS); iterates to a clean exit; applies every fix. It never pronounces the verdict.
- **Tertiary reviewer = GLM — a third, *different* model.** One fresh independent pass over each spec, and the tie-breaker on
  load-bearing "decisive test" / reframe judgments. Deliberately a *third* model so that an orchestrator+executor agreement can't
  quietly rubber-stamp itself.
- **Clean agents = fresh Claude subagents.** Each verification leg (fidelity read, adversarial-with-ablation) on its own fresh subagent
  with no memory of how the work was built.

**The hand-off ordering matters** (this is the part that's easy to lose): iterate the **executor to green** → **one** **tertiary
(third-model)** pass → orchestrator **folds** the findings → **executor confirms to green again.** The executor *bookends* the
third-model pass (so a problem the fold introduced can't slip through), and the third model's suggestions get re-checked by the
executor because a general-purpose model tends to **over-import standard/textbook assumptions** that may not hold in a non-standard
framework. Run these **in sequence, never in parallel** — parallel reviews duplicate findings and waste context, and you lose the
fold-between-passes that makes each pass build on the last.

**Why three *different* models, concretely.** On one load-bearing "decisive test," a proposed A-vs-B discriminator was actually
*guaranteed* to pass by a mean-value-theorem identity — a tautology dressed up as a test. **Both the orchestrator (Claude) and the
executor (Codex) reviewed it SOUND.** Only the *third, independent* model (GLM) caught it. Two same-family reviews would not have.
That single episode is the whole argument for splitting the roles across heterogeneous engines rather than iterating one model against
itself.

**Instantiating with other setups:**
- **You have one agent framework that can spawn subagents (e.g. Claude Code alone):** run each **verification
  leg** as a *separate fresh subagent* — a subagent has no memory of your reasoning, which is exactly the independence you need. Keep the
  orchestrator (you) out of writing code. ⛔ **Do not delegate coding to a subagent** — an agent is a REVIEW instrument, never a coder
  (canonical rule: `docs/development_pipeline.md` Roles table). Coding stays with the **Executor** role (in this project: Codex). If the
  only tools you have are review agents, **halt and get a human-gated executor path** rather than promoting a review agent to coder.
  Spawn a *different* fresh subagent for the fidelity read and *another* for the adversarial-with-ablation.
- **You have a second model available (a different LLM):** use it for the **tertiary review** of specs and for a genuinely independent
  re-derivation. Cross-model disagreement is the highest-signal check you have.
- **You have two computation engines** (e.g. two CAS's, or a symbolic + a numeric route): run the headline both ways and require they
  **agree, each assembling the headline from its own primitives** (see §4).
- **Minimum viable version (single model, single engine):** you still get most of the value from — pre-registration, able-to-fail
  checks, truth-in-output, **fresh-subagent** verification legs, the no-mutation rule, and the human gates. What you lose is
  cross-*model* and cross-*engine* independence; compensate by making the adversarial-with-ablation leg especially aggressive and by
  routing "decisive test" judgments through a fresh subagent that never saw the build.

**Escalation rule.** A **math-correctness** question (a wrong constant, a tolerance, a suspected tautology) is resolved by
re-derivation between orchestrator and executor. A **convention-or-scope** question (which of two consistent choices; whether to drop a
claim) goes to the **human** — no engine can decide between two internally-consistent options. Escalate to the human only when a fix
changes the *conceptual nature* of the result (an effect stops holding, a sign flips, a structural claim must be abandoned).

---

## 3. The loop — run every "gate" through these phases

Decompose the work into **gates**: small, individually falsifiable questions. Do **one at a time**, with a human gate between them.

- **Phase 0 — Ground.** Read the vision/spec and prior results. Don't re-derive what's already established (re-derivation is also a
  chance to silently contradict a banked fact).
- **Phase 1 — Specify + review, before any compute.** The orchestrator writes a **directive**: the requirement, the acceptance
  criteria, and the **full verdict ladder including every way it can fail** — but **not *how* to compute it** (pre-designing the route
  biases the result and removes an independent design). Then: executor design-reviews the directive → iterate to clean → *one* tertiary
  pass → fold the findings → executor confirm-review → clean → **human gate**. *No math has run yet; most rigging is visible in the
  spec.* *(For numerical / simulation work, the directive also includes the pre-registration brief and the target-blind freeze — see
  §6.)*
- **Phase 2 — Review the run prompt too.** The per-run instructions are where subtle bias re-enters; review them before the expensive
  run.
- **Phase 3 — Execute.** The executor codes and runs it, **on two independent engines where a second is possible**, iterating to a
  clean exit. Run the **dimensional / units firewall inline** (§4). *If this gate involves a **numerical solve**, the numerical-work
  gates in §6 (pre-registration freeze, validation suite, convergence discipline) apply **before any physics claim**.*
- **Phase 4 — Tri-review (each leg independent):**
  - **(a) Re-run** the existing scripts *unchanged* to confirm reproducibility (the orchestrator may do this; it edits nothing).
  - **(b) Fidelity audit** — a *fresh* agent checks code against the intended equations **term by term**.
  - **(c) Adversarial-with-ablation** — a *fresh* agent tries to break it: mutates inputs, confirms each control actually *fires*, hunts
    for pass-by-construction.
- **Phase 5 — Claim audit (a *separate* layer from the computation checks).** A fresh reader compares every headline *claim and word*
  in the write-up ("exact," "proven," "isolated," "predicts," "derived") against what the step **actually established**. The Phase-4
  checks confirm the computation is *correct*; they are structurally blind to a correct computation being *over-stated*. Downgrade any
  claim the underlying step doesn't support. (Also confirm that what *landed* exercises the claim — a good executor sometimes corrects
  a flawed instruction; verify the destination, not the obedience.)
- **Phase 6 — Remediate (if any hole is found in Phase 4 or 5).** Fix the artifact *regardless of whether the verdict changes* (leave
  nothing a reader can call "soft"). The **executor** applies the fix; **re-verify on a NEW fresh agent** — never the one that just
  reasoned about the fix.
- **Phase 7 — Gate + record.** Update your "you are here" status note; **commit only when the human asks**; get a human gate before the
  next gate.

---

## 4. What a *real* check looks like (and the anti-patterns that fake it)

Every one of these anti-patterns is a way a check *looked* rigorous and was hollow. Self-check against them before trusting a result.

- **Able-to-fail, proven.** A control must change its verdict when you feed it a wrong input. *Anti-pattern:* a control whose FAIL
  branch is a typed-in boolean, or a "counterfactual" that was hardcoded. **Fix:** the adversarial leg must *mutate an input and re-run*
  and watch the verdict flip.
- **Derived, not emitted.** Verdicts, counts, and classifications must be *computed*. *Anti-pattern:* a "result-emitter" — a script
  that builds something adjacent to the real calculation, then prints the answer as a literal keyed on a config label. This passes
  re-runs and even cross-engine agreement (both engines transcribe the same literal). Only an independent term-by-term reader catches
  it. **Tell:** grep your own output for verdict strings / counts / matrix entries you can't trace back to the inputs.
- **Independent agreement, not an echo.** Two engines agree only if **each assembled the headline from its own primitives**.
  *Anti-pattern:* "agreement" between two byte-identical copies of the same computation (vacuously `0 == 0`).
- **A units check that can fail.** Check dimensions on the *actual* expressions with units restored, as you build them, and include an
  ablation that *should* break homogeneity and confirm it does. *Anti-pattern:* a check that forces itself to pass by *back-solving a
  free parameter* — then it's comparing a number to itself. **Fix:** put the teeth on a quantity fully fixed by sourced inputs;
  corrupt a *sourced input* and confirm the check fails.
- **Fidelity ≠ consistency.** A manufactured-solution / self-consistency test derives its own target from the same operator, so it
  proves the code is *self-consistent*, never that the operator is the *intended* one. A faithful-but-wrong operator (wrong sign,
  wrong derivative order, wrong measure factor, cubic where it should be quartic) sails through. **Fix:** a separate term-by-term
  fidelity read against the source equations.
- **Read the output on a negative verdict.** For a "can't-be-done" result, inspect the artifact: what initialized, how many
  iterations/substeps actually ran, which code path fired. *Anti-patterns to grep for:* a "prefer existing / use cached" flag that
  bypasses the computation; break-on-first-failure that skips the rest of the schedule; a structurally-hardcoded "PASS"; a cheap proxy
  quietly swapped in for the rigorous metric the spec demanded.
- **Contest a "wall" before believing it.** A numeric "it fails / folds here" is often an under-resourced solver, not physics. Re-run
  with full budget; a convincing "wall" can dissolve. (Separately, ask whether crossing it is even on the critical path.)

---

## 5. Honesty rules (the epistemics)

- **No mutation after results.** (Prime directive 5 — restated because it is the most-broken, highest-value rule.) If the branch
  misses, find the next branch on its *own* pre-registered terms; do not let this branch's miss send you hunting for a pass.
- **Calibrate-then-predict; score the held-out surplus.** Calibrating a few inputs on known data is fine — *then* the score is the
  *new* results you reach with **no further calibration**. The **form** of a law is the prediction; its overall **magnitude** may be a
  calibration input. Count **all** degrees of freedom honestly — including that *trying several candidate forms and keeping the best is
  itself a degree of freedom*. Forbidden: claiming one quantity as both calibrated-to and predicted; using as many knobs as predictions.
- **A miss is a missing-physics signal, never a rescue knob.** When you fall short, do **not** add a free parameter to close the gap
  (that's the "fits anything" fudge). Re-examine the equations and *derive* the missing term. Reaching for a parameter "to make it
  behave" is the tell — stop.
- **A no-go between requirements is a first-class success.** If the honest finding is that the requirements are mutually incompatible,
  that is a real result — report it plainly; never soften or rescue a result that breaks the idea.
- **Verify the claim, not just the computation.** All the checks above confirm the computation is right. They are blind to a correct
  computation being *over-stated* in the write-up ("exact," "proven," "isolated") beyond what it shows. Read each headline claim
  against what the step actually established.
- **Framing split.** Explore speculative/ambitious framing privately if you like; every *public* artifact states only what was shown.

---

## 6. Additional gates for numerical / simulation work

The core loop (§3) is enough for symbolic derivation. The moment a **numerical solver** enters, three more gates become mandatory —
this is where confirmation-bias tuning and noise-mistaken-for-signal live.

**6.1 Pre-registration brief + freeze-hash (no post-hoc refit).**
- Before running, write a brief: which branch / model class, which targets count as pass and which as fail, which observables and *how*
  extracted, the error budget the result must clear, and what a **trustworthy miss** looks like. Commit to it; do not edit it in
  response to data.
- **Freeze the model class target-blind:** forms, domains, tolerances, mesh/grid, the calibration objective, the optimizer, the
  tie-breaker, and *every* discrete family choice — then hash it, with no reference to any target value.
- Calibrate **only** the declared parameters to the stated anchor; freeze the calibrated values (hash again).
- Evaluate the held-out observables with the freeze locked. **Whatever comes out, stands** — no post-residual refit, no
  moving-boundary bailout added after seeing the residual, no reporting a fresh attempt as a "continuation" of a failed one.
- **Target-blind extraction:** the solver and the extraction pipeline do not know the targets; targets meet extracted values only at a
  single designated, frozen comparison step. Without this, the unconscious tuning loop is invisible from inside your own head.

**6.2 Validation-suite gate — no physics claim before it passes.** A solver earns trust only after it clears, on *independent*
benchmarks:
- a known-analytic / linear limit;
- manufactured-solution tests per operator;
- a published benchmark with a known answer;
- mesh/grid refinement at **≥3 levels** showing the expected convergence order;
- conservation diagnostics over the run;
- an explicitly stated **noise floor**.
*A result without this stack is not a result* — it is an interesting number. (Note: manufactured-solution tests prove self-consistency,
not fidelity — they do **not** substitute for the §4 term-by-term fidelity read.)

**6.3 Convergence discipline — characterize before concluding.** Do not call convergence *or* non-convergence from too few points. Run
the axis **≥7 points**, fit an extrapolation model (e.g. `value = limit + c/nᵖ`), require the extrapolated limit **jackknife-stable**
(drop the coarsest 1–2 points → the limit must not move), and confirm an **independent axis** agrees. A few rising points are not a
limit; a few wobbles are not "no limit." Invest in a heavier solver to *get* the points rather than over-reading sparse ones.

---

## 7. Pre-flight — the agent's self-check before reporting any result

Before you present a gate's result, confirm out loud:
1. ☐ The pass/fail criteria were fixed **before** this run, and I have **not** changed them in response to what I saw.
2. ☐ Every verdict-bearing number was **computed from inputs**, not typed in; I can trace each to its source.
3. ☐ Each control **fired** under a mutated input (I ran the ablation; I didn't just assert it).
4. ☐ Any cross-engine / cross-model agreement compared **independently assembled** results, not copies.
5. ☐ A **fresh** agent (not me, not the one that built it) did the fidelity read and the adversarial-with-ablation.
6. ☐ A **fresh** reader ran the **claim audit** — every headline word matches what the step actually established (no overclaim).
7. ☐ For any negative verdict, I read the **output artifact**, not just the prose, and confirmed the attack actually ran.
8. ☐ *If numerical:* the model class + calibration were **frozen target-blind** before the run (no post-residual refit); the
   **validation suite passed** before any physics claim; and convergence was **characterized** (≥7 pts, fit, jackknife-stable,
   cross-axis), not eyeballed.
9. ☐ I am reporting exactly what was shown — no rescue knob, no softened miss, no claim beyond the computation.
10. ☐ This is one gate; I am **stopping for a human gate** before rolling forward.

If any box is unchecked, the result is not ready — say so.

---

## 8. Kickoff snippet (paste to start a session)

> Follow `methodology_operating_manual.md` for this work. First, tell me how you'll instantiate the roles (orchestrator / executor /
> tertiary reviewer / clean agents) with the tools available in this repo, and how you'll get an *independent* verification leg. Then
> proceed **one gate at a time**: write the directive (requirement + acceptance + every failure mode; not the computation route), get
> it reviewed, and **stop for my approval before running anything**. Treat a clean "it all works" as suspicious, put the verdict in the
> computed output, and never change a test after seeing its result.

---

*This manual is deliberately tool-agnostic. A longer, human-facing companion with the full rationale and the case studies behind each
rule lives alongside it (`methodology_guide.md`).*
