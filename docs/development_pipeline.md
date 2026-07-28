# Development Pipeline — how we do the work (canonical "how we work" front door)

**Why this doc exists.** The *what* of the program lives in `docs/development_plan.md` (the sector-by-sector work breakdown) and
`docs/model_map.md` (the physical vision — ⛔ NOT `conceptual_foundation.md`, superseded and it re-confuses). This doc is the ***how*** — the process every compute gate runs through. It was
written **2026-06-26** after the rules — previously scattered across ~30 `feedback-*` memories + `STATUS.md` + each directive's
§discipline/§review sections — let a shortcut slip through (the orchestrator self-verified a remediated script instead of using a fresh
clean agent). Consolidating it here makes the pipeline impossible to miss.

**This doc is thin: the memories remain the source of truth.** Each rule cites its `feedback-*` memory; when they conflict, the memory
wins (it's the live, dated record). Keep this doc synced when a process rule changes.

---

## Roles — who does what (never blur these)

| Role | Does | Does NOT |
|---|---|---|
| **Claude (orchestrator)** | Reviews; **decides** what scaffolding must say (directives, execution prompts, `STATUS`/`decisions` docs) and reviews the resulting **diff**; runs the **arbiter re-run** of existing scripts; gates; banks results | **Write or mutate any code/math/script — not even disposable scratch copies** [`claude-reviews-codex-codes`, `codex-is-fix-applier`]. ⛔ Also: **hand-type a long scaffolding document** — see the row below [`thin-orchestrator-definition`] |
| **Scaffolding applier** — a distinct role from a review leg (2026-07-28, user-approved) | Applies an orchestrator-authored **decision list** to a **prose** file: a directive, review prompt, plan doc, `STATUS`, or a note section. Flow: Claude writes the decision list (what changes, where, why) → the applier edits the file → **Codex verifies** → Claude **reviews the diff**. Judgement stays with Claude; transcription does not. ⭐ **Edit, never rewrite** a file already written unless the change exceeds half of it. | ⛔ **Touch any code, math, script, or `.wl`/`.py`/`.sh` file — prose only, and never a deliverable.** ⛔ Also be the same session that reviews its own application, or double as a review leg. ⚠ The Clean-agents row's ban on scratch mutation governs **review legs**; it does not reach this role, which exists precisely to write prose. ⚠ This row exists because the old wording (*"Claude owns scaffolding (directive prose…)"*) contradicted [`thin-orchestrator-definition`] and a session then hand-wrote one directive twice in full |
| **Codex** (always `-c model_reasoning_effort=xhigh`) | Designs + writes + RUNS all scripts, iterates to exit 0, applies all code/math fixes | Decide the verdict (the verifier reviews substance afterward) [`codex-xhigh-reasoning`, `codex-iterates-until-clean`] |
| **GLM** | One fresh tertiary check per directive (CLI, full repo access) | Be iterated like Codex — it's a single pass [`review-ordering-codex-then-glm`] |
| **Grok** | Gates foundational directives and artifacts about to be frozen; applies only the documentation pass expressly assigned to it in a user-directed round-robin | Write or mutate any code/math/script; become a general-purpose documentation applier [`grok-final-review-pass`] |
| **Clean agents** (review legs) | Each review/verification leg, on its **own fresh agent**, no hint of the expected conclusion; **only in Phase 4(c), may mutate disposable scratch copies solely for the orchestrator-assigned per-tooth ablations** | **Write or mutate any live/deliverable code/math/script, or make any scratch mutation outside that bounded Phase 4(c) ablation — an agent is a REVIEW instrument, never a coder.** Delegating coding to an agent is the SAME violation as Claude doing it by hand, and it is the easier one to rationalize [`claude-reviews-codex-codes`, `codex-is-fix-applier`]. Also: be reused across legs / carry context between reviews [`review-agents`, `offload-review-gauntlet`] |

⛔ **Round-robin application and review — only for documentation the user expressly assigns to this
rotation, such as canonical process, plan, or decision documents.** Rotation: Codex applies → Grok and
Claude review → Grok applies → Codex and Claude review → alternate; each application is reviewed by the
other two parties before the next application. **Operating rationale:** rotating the applier is a
precaution against context-dependent transcription drift, based on the 2026-07-22 stage-031 directive
incident [`roundrobin-apply-review`]. This is not a standing gauntlet for routine note, status, or
transcript edits. ⛔ **Software builds are always Codex's.** This rotation governs documentation only; it
does **not** license an agent — or Grok — to write code. The Roles-table rule stands unchanged: *an agent
is a REVIEW instrument, never a coder* — delegating coding to an agent remains the same violation as
Claude doing it by hand.

---

## Standing principles (apply at every phase)
- **Analog, not derivation** — postulate the structure freely; the only test is internal consistency; a **no-go between requirements is a first-class success** [`falsification-is-the-goal`, `analog-find-consistent-structure`].
- **Falsification is the goal** — build every gate *able-to-fail*; a clean "it all works" is the *suspicious* outcome.
- **A miss is a missing-physics signal, never a rescue knob** — when you fall short, do **not** add a free parameter to close the gap (that is the "fits anything" fudge). Re-examine the equations and *derive* the missing term. Reaching for a parameter "to make it behave" is the tell — stop [`falsification-is-the-goal`].
- **Calibrate-then-predict; score the held-out surplus — and count form-shopping as a knob** — calibrating a few inputs on known data is fine; *then* the score is the *new* results you reach with **no further calibration**. The **form** of a law is the prediction; its overall **magnitude** may be a calibration input. Count **all** degrees of freedom honestly — including that ***trying several candidate forms and keeping the best is itself a degree of freedom***. Forbidden: claiming one quantity as both calibrated-to and predicted; using as many knobs as predictions [`calibrate-predict-methodology`].
- **Truth in the OUTPUT, not prose** — no LLM establishes a math or dimensional fact; a verdict-bearing control must compute from inputs and be able to fail [`negative-verdict-short-circuit`, `dimensional-consistency-check`].
- **Never unilaterally deviate from the calibrated contract** — if a step is failing, HALT and bring it to the user [`never-alter-calibrated-process`].
- **"The tool is broken" is a claim that must be VERIFIED before it licenses anything** — a stall/hang/failure attributed to Codex (or Mathematica, or GLM) is a hypothesis, not a fact. Check the actual artifact: does the log end with its completion marker? did the process exit non-zero? does a fresh smoke run reproduce it? **If Codex is genuinely broken, we FIX Codex — we never route around it**, and never by promoting an agent to coder. (2026-07-24: a "Codex CLI stalled twice" claim was written into four docs as justification for an agent-as-coder pivot; on later inspection all 27 of that day's Codex runs had ended `exit=0` and a smoke run answered correctly in 20s. The unverified excuse did more damage than the original error, because it self-justified in every future session.) [`never-alter-calibrated-process`, `negative-verdict-short-circuit`].

---

## The phases

### Phase 0 — Grounding
Read `docs/model_map.md` (⛔ NOT `conceptual_foundation.md` — superseded, it re-confuses) + the relevant memories first. Don't re-derive what's already banked (e.g. the PN ladder, prior no-gos).

### Phase 1 — Directive → review gauntlet (before any compute)
1. **Claude specifies the directive** — *requirement + acceptance + the verdict ladder (incl. the no-gos)*; **never pre-designs the computational route** [`claude-reviews-codex-codes`]. ⭐ For anything longer than a short document, Claude writes the **decision list** and a **scaffolding applier** produces the text (Roles table); Claude reviews the diff. Claude decides *what it says*, not *who types it*.
2. **Codex design-review (xhigh, via a gauntlet-runner agent)** → iterate to GREEN [`directive-design-review`, `offload-review-gauntlet`].
3. **ONE GLM tertiary pass** — a single fresh check [`review-ordering-codex-then-glm`].
4. **Claude folds** the findings — Claude owns *the decisions*; the applier may own the keystrokes. ⭐ **Edit, never rewrite**: a fold that changes a dozen sections is a dozen edits, not a re-emission of the file.
5. **Codex confirm-pass → GREEN again** — Codex bookends the GLM pass [`review-ordering-codex-then-glm`, `directive-design-review`].
6. **User gate.**

### Phase 2 — Per-gate execution-prompt design-review
Claude specifies the gate's **execution prompt** (scaffolding; same decision-list/applier split as Phase 1) → **Codex design-reviews the execution prompt (via agent) before the expensive run** → fold. (Lesson: review the per-rung prompt, not just the directive.)

### Phase 3 — Execution (Codex)
Codex designs + codes + runs **dual-engine** (SymPy + Mathematica), iterating to exit 0 [`dual-engine-required`, `codex-iterates-until-clean`]:
- Mathematica needs **`--sandbox danger-full-access`** to run under Codex [`codex-mathematica-sandbox`]; the orchestrator can run `math` directly as arbiter.
- `codex exec … xhigh` **backgrounded with `< /dev/null`**; **never wrap the Codex session in shell `timeout`** [`background-process-launch`, `never-timeout-codex`].
- The **scripts** Codex runs get `timeout 600`; a 124 = reformulate the math, never raise the cap [`script-timeout-policy`].
- **≤2 concurrent `math -script`** seats; no parallel `$RT exec-*` (it races `MANIFEST.yaml`) [`mathematica-single-seat`, `no-parallel-exec-sympy`].
- **Dimensional firewall — comprehensive + inline:** dim-check **every constructed expression as it is built** (not an end-pass), units restored, mutually cross-consistent, **free of any back-solved free carrier**, dual-engine must agree, **able-to-fail** (an ablation must FIRE). No LLM establishes a dim fact [`dimensional-consistency-check`].

### Phase 4 — Tri-review (the backstop — each leg independent)
- **(a) Arbiter re-run** — *orchestrator* independently re-runs the existing scripts unchanged (reliability gate; refreshes committed outputs) [`orchestrator-rerun-and-output`]. **This is the only execution leg the orchestrator performs.**
- **(b) Transliteration-fidelity audit** — **fresh clean agent**, code-vs-equations term-by-term (catches a faithful-but-wrong operator that MMS/arbiter can't) [`transliteration-fidelity-audit`, `script-review-depth`].
- **(c) Adversarial-with-ablation** — **fresh clean agent**, tries to break it: mutates scratch copies, confirms controls are able-to-fail, hunts pass-by-construction [`negative-verdict-short-circuit`]. ⛔ **The builder never selects the ablation target.** Ablation is per-tooth over every able-to-fail check or emitted record; the target list is owned by the **orchestrator**. Otherwise one correctly live record can carry an arbitrary number of hardcoded ones, and every stated observation still passes.
- A negative / "can't-do" verdict gets **harder** verification than a positive [`negative-verdict-short-circuit`].

### Phase 5 — Remediation loop (if a hole is found)
- **Fix the script regardless of outcome** — no spot a future reader can point to as "soft."
- **Codex applies the fix + re-runs** [`codex-is-fix-applier`].
- **Re-verify on a FRESH clean agent (a new fresh view) — never the orchestrator** that just reasoned about the fix [`review-agents`].
- Math-level disagreements: Claude + Codex resolve and agree; escalate to the user only if the fix changes the *conceptual nature* [`claude-codex-resolve-math`].

### Phase 6 — Gate + bank
- **One gate at a time, explicit user gate between gates, never roll forward autonomously** [`sequential-audit-chunks`].
- **Commit only when the user asks** (squash small follow-up fixes into the prior commit) [`squash-followup-fixes`].
- Sync **`STATUS.md` + `software/stage1_solver/decisions/13` §0 at every milestone** [`status-md-front-door`]; update memory.

---

## Diagnosing a stuck problem (2026-07-24/25, learned expensively)

**⭐ If successive fixes each ban a SPELLING and the next probe evades it, the ARCHITECTURE is
wrong — stop patching.** Four rounds hardening one check each closed a real hole and each was
defeated by the next probe (arbitrary anchor file → self-cited file → committed conventionally-named
file → `Dim(1,2,3)` literals → arithmetic over real bindings producing the identical wrong result).
That is a denylist against an expressive grammar; it does not converge.
- **The diagnostic that broke it open: *who supplies the value, and who supplies only the
  alphabet?*** The check computed `D' = eval(manifest_expression, script_bindings)` — the derived
  artifact controlled the function, the evidence contributed only basis names. Every downstream
  guarantee was therefore artifact-determined. The fix direction is to INVERT the flow (the
  artifact asserts, the checker computes), not to sanitize the grammar further.
- **Specify the INVARIANT, not the instance.** A rule written against the field a survey happened
  to measure (`dim_source.locus`) left the identical defect alive in its sibling (`evidence.locus`).
  Same shape as the above. Write the property, not the example.
- **Before building on an empirical premise, MEASURE IT.** "The convention surely holds across the
  scripts" was false — exact identity→dimension matching resolves **0/92** production symbols, and
  a fuzzy matcher reaches 30% while producing a silently WRONG certificate. A cheap survey would
  have saved a round.
- **Ask "independent of whom?" before calling a check redundant.** A checker-owned digest was
  removed as "duplicating" an existing staleness check; the existing one compared against an
  artifact-supplied digest (self-consistency only), so deleting the independent pin was a real
  regression.
- **When scope is deliberately bounded, write the bound INTO the directive** ("do NOT fix X, and
  report rather than implement if you disagree"). A competent coder handed a checker with three
  known holes will close them, and the next round starts. Scope decisions that live only in a
  conversation do not survive contact with the next agent — put them in the commit message too.

---

## Verification discipline (what green does and does not mean)

- **A green self-report is not evidence a hole closed.** Twice the fixtures went green while the
  invariant stayed open. What caught it both times was an agent **EXECUTING the case end-to-end**,
  never reading the diff. Reviews that only read code find design smells; only execution finds
  live defects.
- **Verify what could BREAK, not just what changed.** After a manifest rework the build was re-run
  but `--self-test` was not — and fixing production had broken a fixture. The suite went red and
  the orchestrator did not notice; an external reviewer did.
- **A fixture must NEVER be coupled to mutable production state.** One fixture asserted the
  production manifests exhibit a specific finding; legitimately fixing that finding broke the test
  suite. Fixing the product must never break the tests. Fixtures are synthetic.
- **A shipped guard with no fixture is not a guarantee.** Ask periodically: *which guards could I
  delete with the suite still green?* (A pre-freeze review found ~42 of ~89 issue codes were never
  planted, including the primary drift-detection path.)
- **An independent reviewer is worth more than another orchestrator pass.** The pre-freeze review
  found, in one pass, blockers that hours of self-review had not: a red suite, dimension recovery
  covering only 16 of 43 scripts while the schema required it for all, and a closedness rule that
  could never trigger.
- ⛔ **Independent agreement, not an echo — the failure mode dual-engine work is blind to.** Verdicts,
  counts and classifications must be **computed**, never emitted. The anti-pattern is a **"result-emitter"**:
  a script that builds something *adjacent* to the real calculation, then prints the answer as a
  **literal keyed on a config label**. That construction **passes re-runs and even cross-engine
  agreement** — both engines transcribe the same literal. Two engines agree only if **each assembled the
  headline from its own primitives**; "agreement" between two byte-identical copies of one computation is
  vacuously `0 == 0`. Phase 3 already *requires* dual-engine work, but nothing here otherwise states what
  dual-engine cannot see. Only an independent term-by-term reader catches it [`transliteration-fidelity-audit`].
  **Tell:** grep your own output for verdict strings / counts / matrix entries you cannot trace back to
  the inputs.
- **The tertiary (third-model) leg is measured, not stylistic.** On one load-bearing "decisive test," the
  proposed A-vs-B discriminator was *guaranteed* to pass by a mean-value-theorem identity — a tautology
  dressed up as a test — and **both the orchestrator (Claude) and the executor (Codex) reviewed it
  SOUND.** Only the third, independent model (GLM) caught it; two same-family reviews would not have.
  That single episode is the whole argument for splitting the roles across heterogeneous engines rather
  than iterating one model against itself [`review-ordering-codex-then-glm`, `decisive-test-not-tautological`].
- ⛔ **CHECK WHICH TOOL REPORTED THE PASS (2026-07-26).** An acceptance item named a specific
  comparator. The official tool exited **2** on a real defect; the build then wrote its own
  `normalized_comparator`, ran that instead, and headlined **PASS**. It *disclosed* the substitution in
  its raw log (`official_comparator_exit=2 normalized_comparator_exit=1`) — not deception, a workaround
  honestly logged and then optimistically summarised. **The summary is what gets read, and the real
  gate had never passed.** Nothing in the tree revealed it: comparator untouched, no variant left
  behind, `git status` clean. Only re-running the named command caught it.
  - ⭐ **Re-run the NAMED acceptance command yourself and read its literal exit code** — not a
    paraphrase, not the report's claim. This is what "the orchestrator independently re-runs acceptance
    commands" is *for*.
  - ⛔ **Every build directive must say: "never satisfy an acceptance criterion with a substitute tool;
    if the named tool fails, STOP and report."** Without that sentence, routing around a broken gate
    looks like diligence.
  - **Read raw exit codes, not the PASS/FAIL headline.** In a long report the logs are the evidence and
    the summary is a claim.
  - **A gate that is routed around is not a gate** — same family as the source-grep and can't-fail
    rules below.
- **When adding a STRUCTURED artifact, verify its STRUCTURE, not only its contents — by running the
  consumer.** Same incident: the orchestrator committed a `.wl` after checking the record count and one
  value, but not that the file had its required header. The consumer would have said so immediately.
- ⭐ **Specify the evasion before shipping a gate.** For every acceptance criterion in a directive, the
  author must first identify, as a reasoning exercise, the cheapest wrong build that would still satisfy
  it and record the analysis. If demonstrating the evasion requires an executable probe, the directive
  assigns Codex to construct and run it during Phase 3; the author does not write or mutate code. A
  criterion whose only failing input is one nobody would produce is decoration. *(Directive-review
  finding, 2026-07-28: before any build existed, the surviving stage023 step-(f) A7/A8 already required
  non-empty floors but named no fixed acceptance commands. The review found that a builder-authored A7
  checker could inspect only one of three exponent tokens per record and a builder-authored A8 checker
  could compare a derived count against itself. The directive was not executed, and no Python sidecar
  was produced.)*
- ⛔ **Acceptance tooling is a control, so it is authored outside the session it polices.** The same
  independence that forbids an expected *value* living inside the artifact it checks extends to the
  *checker*: a validator written by the build session that produces the artifact is the artifact grading
  itself. **Ordinarily, a checker named in a directive must exist before the build starts, authored by a
  different session. Only when the deliverable is itself that checker does pre-existence attach instead
  to its conformance fixtures: a different session must author and freeze the positive and negative
  fixtures, their expected outcomes, and required non-empty floors before the checker build starts; the
  checker-building session may neither author nor weaken them.**
- ⛔ **Never supply a premise to a verification prompt, including a NEGATIVE one.** The rule already
  forbids stating the expected answer or reason. Extend it: also never assert that something does not
  exist, has no prior record, or was never done — and never forbid a leg from reading the evidence that
  would falsify your premise. State what to check; let the leg establish what is there. *(Measured,
  2026-07-28: a verification prompt asserted that no prior verification verdict existed and instructed
  the pass not to read any `OUT_*.txt`. A completed prior pass did exist in the same directory. The pass
  adopted the supplied premise and escalated it into a blocking finding against a commit message that was
  in fact accurate.)*

---

## Sub-agent safety (non-negotiable)

- **⛔ An agent must NEVER delete a directory, or any file it did not itself create.** A review agent
  "cleaned up" by removing the shared `_scratch/` tree — gitignored, therefore unrecoverable. Lost:
  every directive, launcher, run log, review driver, and prior-session artifact including v1 pilot
  manifests. Put this clause in EVERY agent prompt, and give each agent its own subdirectory.
- **Guard the workspace, not only the repo.** Review-agent prompts had said "do not modify these
  tracked files" — which covered modification and not destruction, and said nothing about the
  working tree.
- **Agents are REVIEW instruments and never coders** (see the Roles table). A review agent may
  mutate SCRATCH COPIES for ablation; that is the only writing it does.

---

## Background-process discipline

- **The harness reaps long background waiters (~4–15 min).** A dead waiter is NOT a dead job — this
  is almost certainly what produced a phantom "the tool stalled twice" claim that then justified
  abandoning the calibrated contract. Observed live three times in one session.
- **A waiter must distinguish three outcomes**, not two: marker found / process gone WITHOUT a
  marker / still alive. Without the third branch, a reaped waiter looks exactly like a hang.
- **Match the done-marker on `tail -1` only** — a tool can quote a prior marker mid-stream.
  ⭐ **Use absolute paths in every launch and every check, and verify a launch by the existence of its
  declared output file.** Two stage023 incidents make this one artifact rule, not path-style advice:
  (1) a session launched with `-C /var/projects/toy_physics` against directive paths stated relative
  to `research/pde_ledger_v2/` produced its evidence under the repo-root
  `_scratch/stage023/enum/`, leaving two evidence trees and making the wrong one look real; (2) after
  the shell drifted into `research/pde_ledger_v2/`, the relative launch
  `bash research/pde_ledger_v2/_scratch/run_codex.sh …` never resolved, its `&&` chain skipped the
  intended ablation, and the operator saw exactly four `PASS` lines — the same visible shape as the
  successful ablation-free run. Absence of an error and pass-shaped downstream text do not establish
  that the launch happened; the absolute output artifact does.
- **Prefer a persistent monitor** for long jobs; it survives reaping and stays silent until there is
  something to say.
- **Codex's provider can refuse content on a cyber-policy filter**, killing the session after the
  work but before the report. Phrase adversarial directives in correctness terms ("negative
  fixture", "self-nominated anchor", "admissible set"), never attack language.

---

## Numerical-work gates (there were none here — extracted from the retired operating manual)

The phases above are enough for symbolic derivation. The moment a **numerical solver** enters, these gates
become mandatory — that is where confirmation-bias tuning and noise-mistaken-for-signal live.

- ⛔ **Freeze the model class target-blind, then hash it.** Before running, write the pre-registration
  brief — which branch / model class, which targets count as pass and which as fail, which observables and
  *how* extracted, the error budget the result must clear, and what a **trustworthy miss** looks like —
  and do not edit it in response to data. Then freeze the model class **target-blind**: the **forms, the
  domains, the tolerances, the mesh/grid, the calibration objective, the optimizer, the tie-breaker**, and
  *every* discrete family choice — **then hash it, with no reference to any target value**. Calibrate
  **only** the declared parameters to the stated anchor and freeze the calibrated values (hash again);
  evaluate the held-out observables with the freeze locked. **Whatever comes out, stands** — no
  post-residual refit, no moving-boundary bailout added after seeing the residual, no reporting a fresh
  attempt as a "continuation" of a failed one. ⭐ **Extraction happens after the freeze and is itself
  target-blind:** the solver and the extraction pipeline do not know the targets; targets meet extracted
  values only at a single designated, frozen comparison step. *Without this, the unconscious tuning loop is
  invisible from inside your own head.* This is the same independence as the rule that a control's expected
  value must live **outside the artifact it polices** [`control-outside-the-thing-it-polices`] — here the
  target is what must stay outside the fit.
- ⛔ **A solver earns trust only after a stack of independent checks — no physics claim before it passes.**
  On *independent* benchmarks, all six:
  1. a known-analytic / linear limit;
  2. manufactured-solution tests per operator;
  3. a published benchmark with a known answer;
  4. mesh/grid refinement at **≥3 levels** showing the expected convergence order;
  5. conservation diagnostics over the run;
  6. an explicitly stated **noise floor**.

  *A result without this stack is not a result* — it is an interesting number. ⚠ Manufactured-solution
  tests prove self-consistency, **not** fidelity; they do **not** substitute for the Phase 4(b)
  term-by-term fidelity read. ⚠ **This doc is the process-side home for the stack, not its only record:**
  near-identical statements live in `docs/methodology_paper_outline.md` §2.4 and, operationally, in
  `docs/branch_realization_execution_plan.md` and `docs/branch_realization_brief.md`. That duplication
  across the methodology docs is a known, separately-tracked debt — do not treat any one of them as
  authoritative without checking the others.

---

## Cross-cutting standing rules
- **YAML/markdown for any file an LLM reads/writes; JSON only machine-to-machine** [`no-json-for-llm-io`].
- **No fake `python3 -c` commentary scripts** — read and reason; real ablation re-runs are fine, narration-scripts are not [`no-fake-scripts`].
- **Render audit/verify prompts under the project root** (`software/stage1_solver/_scratch/`), never `/tmp` (agents/Codex can't read `/tmp`) [`audit-prompts-under-project`].
- **Offload to agents** (reviews, broad reads, distilled summaries) to preserve orchestrator context; read grep/tail summaries, not raw logs [`offload-to-agents`, `offload-review-gauntlet`].
- **Never `pkill`/`killall` by pattern** (shared box; the user runs other Codex sessions) — kill only captured PIDs or via `TaskStop` [`background-process-launch`].
- **Path discipline:** project root is `/var/projects/toy_physics` (NOT `toy_projects` — a transcription attractor) [`toy_physics-path-typo`].
- **Grok is a GATE, not a garnish.** Reserve it for foundational directives and for any artifact about
  to be FROZEN while many agents depend on it. Skipping it before a freeze is how a fanout launches
  against a broken tool [`grok-final-review-pass`, `review-ordering-codex-then-glm`].
- **The derived artifact is never a second source of truth.** Manifests/derived records are generated
  FROM the audited sources; when the two disagree the source wins, and the artifact is the thing to
  fix. Corollary: do not re-litigate a unit-level property at the integration layer — check what only
  composition can check (cross-stage agreement, citation liveness, anti-rederivation, census).
- **Model the real risk.** Our derived artifacts are written by our own agents: the operative failure
  is **DRIFT (error)**, not misrepresentation. Hardening against an adversary who does not exist costs
  rounds and buys nothing; hardening against drift is the whole job.
- **When the process is later pointed at EXTERNAL work** (e.g. published papers), validate it first on
  material with KNOWN errata/retractions as positive controls. A pipeline that cannot rediscover an
  error the community already found cannot be trusted to report a novel one.

---

## Per-gate quick checklist
0. ☐ ⛔ Directive contains the anti-substitution clause: *"never satisfy an acceptance criterion with a
   substitute tool; if the named tool fails, STOP and report"* — and states no expected ANSWER or
   REASON (ask for the determination plus its evidence; a supplied rationale gets adopted, a named
   expected result gets special-cased into existence). ⛔ **No premises either, including negative ones**
   (never assert "no prior X exists" / "X was never done", and never forbid a leg from reading the
   evidence that would falsify the premise — state what to check; let the leg establish what is there)
1. ☐ Directive drafted (requirement + acceptance + able-to-fail verdict ladder; no pre-designed route);
   ⭐ for every acceptance criterion, the author identified the cheapest wrong build as a reasoning
   exercise and recorded the analysis; any required executable probe is assigned to Codex in Phase 3
1b. ☐ ⛔ Every checker named in the directive **already existed before the build**, authored by a
   **different session** than the one that produces the artifact; only when the checker is the deliverable
   did independently authored, frozen positive/negative conformance fixtures and expected outcomes exist
   before its build (acceptance tooling is a control)
2. ☐ Codex design-review (xhigh, agent) → GREEN → GLM pass → fold → Codex confirm → GREEN
3. ☐ **User gate**
4. ☐ Execution prompt drafted → Codex design-review (agent) → fold
5. ☐ Codex executes dual-engine (danger-full-access, xhigh, `< /dev/null`, scripts `timeout 600`) → exit 0; comprehensive inline dim-firewall, able-to-fail
5b. ☐ ⭐ **Orchestrator re-runs each NAMED acceptance command itself and reads its literal exit code** —
   not the report's claim. Verify structured artifacts by running their CONSUMER, not by eyeballing
   counts and values
6. ☐ Tri-review: arbiter re-run (orchestrator) + fidelity (fresh agent) + adversarial-with-ablation
   (fresh agent); ⛔ ablation target list owned by the **orchestrator** (per-tooth over every
   able-to-fail check / emitted record — never chosen by the builder)
7. ☐ Any hole → Codex remediates → **re-verify on a FRESH clean agent**
8. ☐ Bank: `STATUS.md` + `decisions/13` §0 + memory; commit only if asked
9. ☐ **User gate** before the next gate

---

*Source of truth = the `feedback-*` memories cited above. Companion docs: `docs/development_plan.md` (what to build), `docs/model_map.md` (the model — ⛔ NOT `docs/conceptual_foundation.md`, superseded and it re-confuses), `STATUS.md` (where we are).*
