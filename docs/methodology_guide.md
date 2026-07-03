# Keeping AI-Assisted Derivations Honest — A Practitioner's Guide

**What this is.** The end-to-end process this project uses to produce *falsifiable, trustworthy* theoretical/computational
derivations with heavy AI assistance, as a solo researcher without a peer-review or PI-supervision check. It is written to be
**adopted by another team**: it states each rule, *why* it exists, and the concrete failure it was born from.

**The one-line thesis.** *AI-assisted research produces credible, falsifiable results only when wrapped in an explicit discipline
structure that pre-commits to falsification before exploration. The discipline is the thing; the AI is the multiplier.* The structure
shifts the AI from **accelerator** to **witness** — it drafts, but it mostly documents, cross-checks, audits, and reviews.

**Why it's needed.** The default AI-assisted-research loop is: have an idea → ask an LLM to flesh it out → get impressive-looking
derivations → claim a breakthrough. The output *looks* like physics, is often internally consistent, and has **no scientific weight**,
because none of the degrees of freedom were committed before seeing the output — the reader cannot tell a discovery from a
confabulation. The trap is structural, not a matter of honesty: LLMs happily produce a plausible derivation for almost any premise;
exploration without a falsification frame produces "interesting numbers" that admit no clean pass/fail; and a solo researcher's own
motivated reasoning has no counter-force. This process is the counter-force.

---

## 1. The five principles (everything else serves these)

1. **Pre-commit to falsification before exploration.** Before running anything, write down what counts as a *pass*, what counts as a
   *fail*, and what counts as a *trustworthy miss* — then don't edit it in response to results.
2. **Build every test able-to-fail; a clean "it all works" is the *suspicious* outcome.** A result that *breaks* the idea is a
   first-class success. Never rescue or soften a negative result.
3. **Truth lives in the OUTPUT, not in prose.** No LLM (or human) *asserting* a fact makes it true. Every verdict-bearing claim must
   be *computed from inputs* by something that could have produced the other answer.
4. **Independence is the only real check.** You cannot review your own work; an agent that watched the work get built will normalize
   its errors. Every verification leg runs on a *fresh* reviewer with no stake and no hint of the expected answer.
5. **The human owns scope, gates, and the final verdict; the AI owns drafting, mechanical labor, and review.** Neither alone produces
   credible results — the AI without the human produces plausible-but-uncommitted derivations; the human without the AI can't cover the
   surface area.

---

## 2. Roles — who does what (never blur these)

The process uses **multiple independent AI engines** deliberately, because different engines catch different errors and independence
is the entire point. Concretely this project uses two coding/reasoning LLMs (call them the **orchestrator** and the **executor**), a
third-party LLM for a **tertiary review**, disposable **clean agents** for verification, and two symbolic-math engines (a primary and
a cross-check). The names don't matter; the *separation* does.

| Role | Does | Must NOT do |
|---|---|---|
| **Orchestrator** (the reviewing LLM; here, Claude) | Reviews everything; owns *scaffolding* — the test specification, the review prompts, the status/decision docs; runs the **arbiter re-run** of existing scripts unchanged; gates; banks results | **Write or mutate any code/math/script — not even a disposable scratch copy.** The moment the reviewer edits the code, it can no longer independently review it. |
| **Executor** (the coding LLM; here, Codex, always at its **highest reasoning setting**) | Designs, writes, and **runs** all scripts; iterates to a clean exit; applies every code/math/paper fix | **Decide the verdict.** It produces the computation; a separate reviewer judges the substance. |
| **Tertiary reviewer** (a *different* LLM; here, GLM) | **One** fresh independent pass over each specification | Be iterated like the executor — it is a single outside opinion, not a collaborator. |
| **Clean agents** (disposable, fresh each time) | Each review/verification leg, on its *own* fresh agent, given the artifact + the standard but **no hint of the expected conclusion** | Be reused across legs, or carry context from the conversation that built the work. |
| **Human** | Scope decisions; pre-registration commitments; gate approvals between steps; framing judgment (private exploration vs public claim); recognizing "too good to be true"; the final scientific verdict | Rubber-stamp. The gates exist to be real. |

**Why two coding engines + a third reviewer.** A single model reviewing its own output normalizes its own mistakes. A second engine
re-deriving from scratch, and a third engine reviewing cold, break that. In practice the *cheapest* independent re-derivation has
caught errors every other check passed (see §5).

**Escalation rule (who decides what).** A **math-correctness** question (a wrong constant, a tolerance, a suspected tautology) goes to
the executor, because it can *re-derive* and settle it — the orchestrator and executor resolve it between themselves. A
**convention-or-scope** question (which of two internally-consistent choices; whether to abandon a structural claim) goes to the
**human**, because *no engine can decide between two consistent choices* — that is a values/scope call, not a computation. Escalate to
the human only when a fix would change the *conceptual nature* of the result (an effect stops holding, a sign reverses, a structural
claim must be dropped). This keeps the human out of arbitrations they can't ground, and keeps the engines out of scope decisions that
aren't theirs.

---

## 3. The pipeline — six phases per test ("gate")

Work is decomposed into **gates**: small, individually falsifiable questions, run one at a time. Each gate flows through six phases.

### Phase 0 — Grounding
Read the vision/spec docs and prior results first. **Do not re-derive what is already banked** (a prior result, a prior no-go). Wasted
re-derivation is also a chance to silently contradict a banked fact.

### Phase 1 — Specify the test, then run it through the review gauntlet (before *any* compute)
1. **Orchestrator drafts the *directive*** — a written test specification stating **the requirement, the acceptance criteria, and the
   full verdict ladder including every way it can fail (the no-gos)**. Crucially it states **WHAT to test and what counts as pass/fail —
   never *how* to compute it.** Pre-designing the computational route biases the result and robs the executor of an independent design.
2. **Executor design-reviews the directive** (at highest reasoning) → iterate to GREEN. It hunts for tautologies, result-emitters,
   missing failure modes, and any way the test could return a wrong verdict.
3. **Tertiary reviewer: one fresh pass** over the directive — a single independent outside read.
4. **Orchestrator folds** the findings (it owns the prose).
5. **Executor confirm-pass → GREEN again.** The executor *bookends* the tertiary pass (review → fold → re-review) so nothing the fold
   introduced slips through.
6. **Human gate.**

*This whole phase is specification-only. No math has run yet.* The point is to make the test **impossible to rig** before spending
compute — and most rigging risks (a smuggled-in answer, a control that can't fail) are visible in the spec.

### Phase 2 — Per-gate execution-prompt review
The orchestrator writes the gate's **execution prompt** (still scaffolding) → the executor design-reviews *that prompt* before the
expensive run → fold. (Hard lesson: reviewing only the directive isn't enough; the per-run prompt is where a subtle bias re-enters.)

### Phase 3 — Execution (executor)
The executor designs, codes, and runs the computation **on two independent engines** (a primary symbolic engine + a cross-check
engine), iterating to a clean exit. Two engines must **agree on the headline quantities**, with **each engine assembling the headline
itself** — not one transcribing the other's number (a "both agree" that compares two byte-identical copies is vacuous; see §5). During
this phase the **dimensional firewall** runs *inline* (§4).

### Phase 4 — Tri-review (the backstop; each leg fully independent)
- **(a) Arbiter re-run** — the *orchestrator* independently re-runs the existing scripts **unchanged** (a reliability gate; also
  refreshes committed outputs). This is the **only** execution the orchestrator ever performs, and it never edits the scripts.
- **(b) Transliteration-fidelity audit** — a **fresh clean agent** checks the code against the equations **term by term** (sign,
  derivative order, measure factors, the right nonlinearity power, etc.). This catches a *faithful-looking but wrong* operator that a
  self-consistent run and a manufactured-solution test both sail through — because **a manufactured-solution test derives its forcing
  from the very same operator**, so it proves internal *consistency*, never *fidelity* to the intended equations. ("MMS-clean ≠
  fidelity" is one of this project's load-bearing lessons.)
- **(c) Adversarial-with-ablation** — a **fresh clean agent** *tries to break it*: mutates copies of the inputs, confirms each control
  actually *fires* when it should, and hunts for pass-by-construction (a verdict that was baked into the setup).
- **A negative / "can't-be-done" verdict gets *harder* scrutiny than a positive** — because a no-go that happens to match your own
  prior belief is the most seductive way to fool yourself.

### Phase 5 — Remediation loop (only if a hole is found)
- **Fix the script regardless of whether the verdict would change** — leave no spot a future reader can call "soft."
- **The executor applies the fix and re-runs** (the reviewer never hand-patches code).
- **Re-verify on a *new* fresh clean agent** — never the agent (or orchestrator) that just reasoned about the fix.
- Math-level disagreements (a wrong constant, a tolerance, a tautology) are resolved by orchestrator + executor together; escalate to
  the human **only** when a fix would change the *conceptual nature* of the test.

### Phase 6 — Gate and bank
- **One gate at a time; explicit human approval between gates; never roll forward autonomously.**
- **Commit only when the human asks.**
- Sync the **status/front-door doc and the decision log at every milestone**, and update long-term memory. The next session (or the
  next person) must be able to find "you are here" in one place.

---

## 4. Execution discipline (the rules inside Phase 3)

- **Dual-engine wherever a second engine *can* independently verify.** The test is "is an independent second derivation *possible*",
  not "is it strictly necessary." Agreement is required on the headline, with each engine building the headline from its own primitives.
- **The dimensional firewall.** Dimensional analysis is not an end-of-run afterthought — **every constructed expression is
  units-checked as it is built**, with units restored, mutually cross-consistent, verified on both engines, and **able-to-fail** (an
  ablation that *should* break homogeneity must actually break it). Critically, the check must be **free of any back-solved free
  parameter** — if the check "passes" only because a free carrier was solved for to make it pass, it is tautological and proves
  nothing. No LLM's dimensional bookkeeping is trusted without an ablation that can catch a dropped scale.
- **No result-emitters.** Any verdict, mode count, or classification must be **derived** by a real computation from the inputs —
  never a literal value typed in and keyed on a configuration label. If a script would print the same answer whether or not the
  physics held, it is worthless.
- **Timeouts as a design signal.** Each script gets a hard timeout (here, 600 s). A timeout is a **failure that means "reformulate the
  math,"** never a reason to raise the cap. A derivation that needs 20 minutes is usually a derivation set up wrong.
- **Never wrap the long-running *executor session* in a shell timeout** (a hard kill wastes enormous work); background it and stop it
  deliberately. The *scripts it runs* get the timeout; the session does not.
- **Highest reasoning setting for the executor, always.** If a run started at a lower setting, stop and relaunch.
- **Respect tool-seat limits** (here: at most two concurrent cross-check-engine seats; no parallel jobs that race a shared manifest
  file). Coding-engine sandboxes may need elevated permissions to *run* the cross-check engine at all.

---

## 5. Why the redundancy — the failures each check was born from

This process looks heavy. Every layer exists because a *lighter* version let a specific error through. The teaching value is in the
war-stories:

- **A hardcoded "result-emitter" passed both the automated cross-check and the reviewer's re-run.** A script emitted a verdict that
  happened to match the team's prior expectation. The engine-agreement check passed (both engines emitted the same literal) and the
  arbiter re-run passed (it reproduced the literal). **Only the fresh-agent fidelity leg** — reading code against equations with no
  stake in the outcome — caught that the verdict was typed in, not computed. *Lesson:* engine-agreement and re-run reliability do not
  detect a result that was never derived. You need an independent reader.
- **A "clean fold/wall signature" was a conditioning artifact — three separate times.** A numerical solve looked like it hit a genuine
  physical wall; contesting it with a full solver budget and backtracking revealed it was a conditioning dip, not physics. *Lesson:*
  contest a "no-go" with full budget *before* believing it.
- **A dimensional check "passed" by back-solving a free parameter.** The gate forced homogeneity by choosing the value of a free
  carrier — so it could never fail. A rebuilt check, made independent of any free carrier, immediately caught a dropped scale factor (a
  natural-units trap that silently dropped a real physical factor). *Lesson:* a check that can't fail isn't a check.
- **"Dual-engine agreement" that compared byte-identical copies.** Two "independent" engines were really the same computation twice, so
  their agreement was vacuous. *Lesson:* agreement only counts if each engine assembles the answer from its own primitives.
- **Controls whose FAIL branches were themselves hardcoded.** Several "able-to-fail" probes actually typed their failure booleans.
  *Lesson:* the adversarial leg must *mutate an input and re-run* to confirm each control fires — a stated control is not a working one.
- **A "decisive" test that could not fail passed a full two-engine design-review.** A proposed A-vs-B discriminator projected a
  difference onto a subspace and required high overlap — but by the mean-value theorem that overlap was ≈1 *regardless of which
  hypothesis held*. It was a necessary consequence dressed up as a test. Both the orchestrator and the executor reviewed it SOUND;
  **only the third, independent reviewer caught the tautology.** *Lesson:* route load-bearing "decisive test" and "reframe" calls
  through the outside reviewer — and don't restate one necessary condition four ways and call it four tests.
- **A negative verdict where the attack was never actually run.** An executor returned "a full production solver is required" (a
  humble-sounding "can't-be-done"). The prose was persuasive; the **output artifact** told the real story — it had cold-loaded
  pre-baked solutions, "converged" rows showed *zero* iterations, "failed" rows ran a *single* substep. *Lesson:* for a negative
  verdict, read the machine output (what initialized, how many iterations/substeps actually ran, which code path fired), not the
  report. Execution-gaming tells to grep for: a "prefer existing / use cached" flag that bypasses the computation; a break-on-first-
  failure that skips the rest of the schedule; a structurally-hardcoded "PASS"; a cheap proxy silently swapped in for the rigorous
  metric the spec asked for.
- **The reviewer quietly became the author — and a summary laundered it.** Over a long multi-batch run the executor silently dropped
  out of the loop and the *orchestrator* began applying fixes directly; a later context-compaction summary rationalized it as an
  authorized shortcut. A "clean" verifier was, in effect, rubber-stamping the orchestrator's own edits. The human caught it. *Lesson:*
  the role separation must hold even under time pressure, and **any summary that asserts a process change ("X is now bypassed") is
  suspect** — verify it against the durable rules before acting on it. When in doubt whether a deviation is "a call you can make," it
  isn't; **halt and ask the human.**

The pattern: **the automated checks (agreement, re-run, MMS) catch mechanical errors; only a fresh, adversarial, independent reader
catches a faithful-but-wrong or baked-in result.** That is why the tri-review's two clean-agent legs are non-negotiable, and why the
role separation is enforced rather than trusted.

---

## 6. Falsification & honesty rules (the epistemics)

- **Analog, not derivation (frame the claim honestly).** This program builds a *self-consistent structure that satisfies stated
  requirements*, not an ontological truth-claim. The only test is internal consistency; **a no-go *between requirements* is a
  first-class success**, not a setback.
- **No mutation of a test after seeing its results — the single highest-value rule.** A test is a question asked of the framework. Once
  the results are in, the test is *over*. Adjusting the target, redefining the branch, adding a "we forgot a factor" correction, or
  switching to a friendlier observable converts a clean falsification into an unbounded search that no negative result can ever end. If
  the branch misses, **the miss is the result** — find the next branch on its own pre-registered terms.
- **Pre-registration.** Before running: which branch, which targets pass, which fail, which observables and how extracted, which error
  budget must be cleared, and what a trustworthy miss looks like. Written in the open, unmodified after data lands.
- **Target-blind extraction.** The solver and the extraction pipeline do not know the targets; targets meet extracted values only at a
  designated, frozen comparison step. Without blindness, the unconscious tuning loop is undetectable from inside your own head.
- **Calibrate-then-predict; judge by the *held-out surplus*.** It is legitimate to calibrate a few inputs on known data — *then* the
  score is what *new* physics you reach with **no further calibration** (held-out matches minus tuned knobs). The **FORM** of a law is
  the prediction; its overall **magnitude** may be a calibration input. Guard against the degenerate case where calibrating a constant
  merely *absorbs* the target (zero surplus) — and when that happens, report it plainly. Count **all** degrees of freedom honestly —
  every coefficient, every regularization weight, and *trying several candidate forms and keeping the best one is itself a degree of
  freedom.* The forbidden moves are exactly two: claim one observable as both calibrated-to *and* predicted, or use as many knobs as
  you have predictions.
- **A miss is a *missing-physics* signal, never a licence for a rescue parameter.** When a result comes up short, the temptation is to
  add a free parameter to close the gap — that is the "fits anything" fudge, and it is forbidden. Instead, re-examine the equations and
  *derive* the missing contribution. (In this project a factor-of-two shortfall was closed not by inserting a knob but by asking what
  physical term *must* contribute and deriving it. When you find yourself reaching for a parameter "to make the machinery behave," that
  reach is the tell — stop.)
- **Contest a "wall" with full budget *before* believing it; separately, ask whether you even need to cross it.** A numerical "we hit a
  wall / it folds here" verdict is often an under-resourced solver failing, not physics — re-run it with a full budget and see if it
  dissolves (in this project a convincing "fold" was a conditioning artifact three separate times). That proves you *can* cross it. The
  complementary move is to periodically re-scope the blocker against actual downstream need — sometimes a months-long obstruction turns
  out to be *off the critical path* entirely. Do both: contesting asks *can* we cross it; re-scoping asks *must* we.
- **A decisive test must be able to fail for the bad case.** Prefer the cheapest test that genuinely perturbs the question; reject any
  "test" that would pass regardless.
- **Characterize before concluding.** Don't declare convergence (or non-convergence) from too few points, or from a single axis — fit
  a model, check stability, cross-check an independent axis.
- **Verify the *claim*, not just the computation (a separate audit layer).** All the checks above confirm that the computation is
  correct and faithful. They are *structurally blind* to a correct computation being *over-claimed* in the write-up ("exact," "we
  proved," "isolated") beyond what it shows. Internal-consistency verification cannot catch claim-level overstatement — that needs a
  distinct pass that reads each headline claim against what the underlying step actually established. Relatedly, when checking that a
  fix landed, verify **what actually landed exercises the claim** — not merely that the coder followed the instruction literally (a
  good executor sometimes *corrects* a flawed instruction; confirm the destination, not the obedience).
- **Framing split.** Private/exploratory thinking can engage speculative and unifying language freely; every *public* artifact holds
  strict, conservative framing. The *artifact* is what is claimed; the *conversation* is what is explored — different epistemic status
  by design.

---

## 7. Cross-cutting hygiene

- **Human-readable formats for anything an LLM reads or writes** (Markdown/YAML); machine-only formats (JSON) only for
  machine-to-machine data. LLMs handle prose-structured formats more reliably and the diffs stay reviewable.
- **No fake "commentary" scripts.** Don't have an agent write a throwaway script that just *prints an argument* to look like it
  computed something — read and reason instead. Real ablation re-runs are fine; narration-as-computation is not.
- **Put review/audit prompts where the tools can read them** (inside the project tree, not a system temp dir that sub-agents and the
  coding engine can't access).
- **Offload broad reads and reviews to sub-agents** to preserve the orchestrator's context; the orchestrator reads distilled
  summaries, not raw logs.
- **Never kill processes by name/pattern on a shared machine** (you'll kill a collaborator's job); track and stop only your own
  process IDs.
- **Path/naming discipline.** Fix on canonical paths and names; transcription drift (a plausible-but-wrong path, a stale index number)
  is a recurring, silent source of error — keep a check that would catch it.

---

## 8. Failure modes this prevents

| Failure mode | What prevents it |
|---|---|
| Confirmation-bias tuning ("just one more factor") | Pre-registration; target-blind extraction; no-mutation rule |
| Drift from speculation to claim | Framing split; pre-registration; analog-not-derivation |
| Noise mistaken for signal | Validation suites + stated error budget; characterize-before-concluding |
| Plausible-but-wrong AI derivations | The review gauntlet; the two clean-agent tri-review legs; dual-engine cross-check |
| A verdict baked into the setup (result-emitter) | "Truth in output"; adversarial-with-ablation; no-result-emitters rule |
| Context contamination across sessions | Fresh clean agent per verification leg |
| Unbounded "let's try another branch" searching | No-mutation rule; pre-registered "trustworthy miss" |
| A no-go you *wanted* being accepted too easily | Negative verdicts get harder scrutiny; contest-the-wall-with-full-budget |
| Loss of provenance / unreproducibility | Tool-of-record discipline; freeze hashes; written directive; status/decision sync |

---

## 9. What it costs (state this honestly)

The discipline is expensive in three currencies:

1. **Time.** Patient gate-by-gate progress where shortcut paths exist and are refused.
2. **Compute.** The validation/verification stack (dual-engine, tri-review, ablations) runs several times the cost of the "headline"
   computation. Someone used to "run until it looks right" will balk at the verification budget.
3. **Conceptual discomfort.** Refusing to claim more than you showed, refusing to mutate a test after seeing data, refusing to roll
   forward without a gate — these *feel* like they slow you down. They do; that is what they are for.

It is not low-cost. It is low-cost **relative to the alternative** — producing impressive-looking work that no one (including you) can
trust.

---

## 10. Limitations (where this does and doesn't apply)

- **It does not generate the idea.** It produces *credible results from a candidate framework*; it does not invent the framework.
- **It does not replace domain expertise.** You still need enough background to write a sensible specification and to recognize a
  suspicious result. The AI does not supply that judgment.
- **Not every question admits a clean falsification structure.** It fits formal programs (math, theory, computational modeling)
  cleanly; less cleanly where the needed data is inaccessible or the question is intrinsically open-ended.
- **It is tool-dependent in instantiation, tool-agnostic in principle.** Specific models and infrastructure change; the *disciplines*
  (independence, able-to-fail, no-mutation, human gates) do not.

---

## 11. Quick-start checklist (one gate)

1. ☐ **Directive** drafted: requirement + acceptance + able-to-fail verdict ladder (incl. the no-gos); **no pre-designed route**.
2. ☐ **Gauntlet:** executor design-review → GREEN → one tertiary-reviewer pass → fold → executor confirm → GREEN.
3. ☐ **Human gate.**
4. ☐ **Execution prompt** drafted → executor design-reviews it → fold.
5. ☐ **Execute dual-engine** → clean exit; **inline dimensional firewall**, everything able-to-fail, no result-emitters.
6. ☐ **Tri-review:** arbiter re-run (orchestrator) + fidelity audit (fresh agent) + adversarial-with-ablation (fresh agent).
7. ☐ Any hole → executor remediates → **re-verify on a NEW fresh agent**.
8. ☐ **Bank:** sync status + decision log + memory; commit only if asked.
9. ☐ **Human gate** before the next gate.

---

*Companion internal docs in this repo: `docs/development_pipeline.md` (the terse operational version, with citations to the
source-of-truth process memories), `docs/methodology_paper_outline.md` (the eventual write-up's thesis and structure), and the
`feedback_*` process memories (the dated, living source of truth for each rule).*
