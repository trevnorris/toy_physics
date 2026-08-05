# Development Pipeline — how we do the work (canonical "how we work" front door)

**Why this doc exists.** The *what* of the program lives in `docs/development_plan.md` (the sector-by-sector work breakdown) and
`docs/model_map.md` (the physical vision — ⛔ NOT `conceptual_foundation.md`, superseded and it re-confuses). This doc is the ***how*** — the process every compute gate runs through. It was
written **2026-06-26** after the rules — previously scattered across ~30 `feedback-*` memories + `STATUS.md` + each directive's
§discipline/§review sections — let a shortcut slip through (the orchestrator self-verified a remediated script instead of using a fresh
clean agent). Consolidating it here makes the pipeline impossible to miss.

⭐⭐ **THIS DOC IS THE OPERATIVE PROCESS DOCUMENT.** What an agent or a session is to *do* is decided here. The `feedback-*` memories
hold the **why** — the incident each rule came from, the measurement behind it, the project facts — and each rule below cites its
memory so the reasoning is one hop away. ⛔ **There is no precedence question and no second process:** a conflict between this doc
and a memory is a **defect to fix**, in whichever of the two is wrong. Keep both synced when a process rule changes. *(This replaces
"the memories remain the source of truth; when they conflict, the memory wins", which handed agents two processes to choose between.)*

---

## ⭐⭐ THE POSTURE, AND THE TEST EVERY CHECK MUST PASS (user decision, 2026-07-29/30)

**This is a two-person toy-physics self-consistency project.** Checks exist to catch a **wrong
derivation**. ⛔ **Immutability is not a discipline here:** files are freely editable, and nothing is
frozen, byte-perfect, or under a custody rule. ⚠ **One carve-out, still undecided** — the numerical
solver's freeze hash, because there it *is* the mechanism of target-blindness; see the *Numerical-work
gates* section below, which is the only place a freeze survives.

⭐ **The governing test — apply it to every check, gate and process rule, including the ones below:**
- **Does it catch a way the PHYSICS could be wrong?** → **keep.**
- **Does it only catch a way the TOOLING could be wrong, or defend against a motivated adversary?** → **cut.**

This promotes *Model the real risk* (below) from a buried line to the governing principle. Our derived
artifacts are written by our own agents, so the operative failure is **drift, an honest error, or an LLM
shortcut that resembles a pass** — never misrepresentation. The apparatus this replaced had grown a
review instrument, an acceptance gate for it, a byte-authority over that gate and a documented procedure
over the authority; a full session was spent at those layers verifying no physics. Details and the
scale-back it authorises: memory [`physics-not-ceremony`].

---

## ⭐⭐⭐ THE METHOD — ⛔ THE CAS ANSWERS, ⛔ NOT THE ORCHESTRATOR (2026-08-03/05, the script rebuild)

⚠⚠ **This section is newer than everything below it and it OVERRIDES anything below that conflicts.** It is
the fold of the script rebuild, in which a session was lost to each of the failures named here.

**1. ⭐⭐ TWO ENGINES EXIST SO THEY CAN DISAGREE.** Independent **construction** — the `.wl` written first,
barred from the registry, ⛔ never a transcription of the `.py`. ⭐ **The disagreement IS the measurement.**
⛔ That is about construction, ⛔ **not** about hiding answers.

**2. ⭐⭐ A SCRIPT PRINTS COMPUTED OBJECTS; ⛔ IT NEVER STATES CONCLUSIONS.** An `emit`/`Print` payload must
be a CAS object. ⭐ Emit **both operands and the residual**, *then* guard — ⛔ a residual asserted zero
always prints `0` and carries no information. ⭐ Interpretation belongs to the **step record**.
⚠ **Measured:** named tags at named lines in **three independently-built steps** were typed prose with no
CAS object behind them, and **eight fidelity review legs missed it** — because *"does it say the right
thing"* and *"does it depend on anything"* are different questions. ⛔ **Do not quote a fraction of the
corpus**; two legs rejected that as unmeasurable.

**3. ⭐⭐⭐ NAME THE OBJECT. ⛔ DO NOT SPECIFY THE RECIPE.**
⚠ **Measured 2026-08-05, and it cost five review rounds.** A spec section specified *how* to obtain a
compared quantity — zero the velocities, divide by a weight, prove the weight nonzero, guard the division.
⇒ **Every subsequent round argued about the recipe**: is the weight unique, is the quotient well-defined,
is the residual a tautology. ⛔ **None of that was a question about the physics.** ⭐ **Every one of those
questions was manufactured by specifying a derivation path the spec never needed to specify** — and the
object was already returned by both engines' constructors.
⇒ ⭐ **Ask for the object. ⭐ Let the engine hand over what it built.**

**4. ⭐⭐⭐ IF YOU ARE DECIDING IN PROSE WHAT THE ENGINES SHOULD COMPUTE, THE METHOD IS INVERTED.**
⚠ **The tell:** many turns of reasoning toward an answer a script would settle in one.
⚠ **The cause, and it is the useful part: the instrument was broken.** The output reader could not parse a
few hundred payloads, so the cross-engine layer could not see the object under dispute; with no measurement
available the orchestrator tried to reason the answer out inside the specification.
⇒ ⭐⭐ **When the measurement is unavailable, FIX THE INSTRUMENT. ⛔ Do not reason around it** — ⭐ a broken
instrument silently promotes prose over evidence.

**5. ⭐ A DISAGREEMENT IS A FINDING, ⛔ not a defect to be prevented by better wording.** ⛔ Do not try to
make divergence impossible in prose; ⭐ name the object, emit it from both engines, **compare afterwards.**

**6. ⛔⛔ DO NOT BUILD BLINDNESS APPARATUS.** ⚠ Two rounds of it were designed and discarded in one session,
and a third was built and had to be torn out. ⭐ The measured failure is **absence of computation**, ⛔ not
anchoring on a known answer — ⇒ quarantine never touched it, and the three clauses kill it structurally.
⛔ **A do-not-read list is a denylist, and a denylist means the architecture is wrong.**
⚠ **Measured:** a denylist was written, and then a scratch-copy instruction materialised every denied path
past it. ⭐ **Withhold exactly ONE thing: an acceptance criterion that references an expected value**,
because a builder iterating to exit 0 converges on any target it can see.

**7. ⭐ PRIOR ART IS AN ORACLE, ⛔ NEVER A PREMISE.** ⭐ Check a computed result against the literature where
they should coincide — ⚠ that is a free external check, better than a review leg. ⛔ Never assume the
literature's result **for our object**: our Lagrangian is not theirs, and ⚠ **their conditions may not be
ours.** ⇒ re-check a **condition** before inheriting its **conclusion**.

### ⛔⛔ AND THE GATES THAT GO WITH IT

**8. ⛔⛔ NO COMMIT BEFORE BOTH LEGS REPORT.** ⭐ The commit is the **last** step. ⛔ Reviewing the
**directive** does not pay the tax for the **build**. ⛔ **The orchestrator's own verification is not a
leg** — ⚠ and it is **most** convincing when it finds things.
⚠ **Measured 2026-08-05:** two builds were committed with zero legs, both were reverted, and a day was lost.

**9. ⭐⭐ THE STOPPING RULE — ⛔ "both legs green" is NOT it.**
> ⭐ **Ship when no outstanding finding changes what the engines COMPUTE, or what we are ENTITLED TO CLAIM.**
> ⚠ Everything else is a **NOTE**: recorded, and shipped with.

⛔ *"A leg that finds nothing is weak evidence"* is an **orchestrator's prior** — ⭐ keep it, ⛔ never treat
a clean report as strong confirmation. ⛔ **It is NOT a leg's instruction**: in a prompt it is a quota, and
a quota is filled. ⇒ ⭐ **put the BLOCKING/NOTE buckets in the prompt, ⛔ not the pressure**, and say that
*"no blocking findings"* is a valid answer.

**10. ⭐ A SHARED SPEC IS PHYSICS-BEARING ⇒ TWO LEGS.** ⚠ The *Roles* table's "prose and process gets one
leg" ⛔ **does not** cover a build directive for an artifact both engines read: ⭐ an error there makes both
engines compute the same wrong thing and **agree**, which is precisely the class cross-checking cannot
catch. ⚠ **Measured:** on such a directive the two legs **disagreed** on the verdict four times, and the
disagreement was the most valuable output each time.

**11. ⭐ IF SUCCESSIVE REVISIONS KEEP BREEDING DEFECTS IN THE MATERIAL JUST CHANGED → CHANGE THE AUTHOR.**
⛔ Do not fold a fourth time. ⚠ **Measured:** three folds each closed the reported defects and introduced
new ones of the same severity; changing the author broke the cycle, and the **defect class** changed
immediately. ⭐ The settled **physics** is not reopened by a new author — ⛔ it is the **drafting** that
fails.

---

## Roles — who does what (never blur these)

| Role | Does | Does NOT |
|---|---|---|
| **Claude (orchestrator)** | Reviews; **decides** what scaffolding must say (directives, execution prompts, `STATUS`/`decisions` docs) and reviews the resulting **diff**; re-runs each **named acceptance command** itself; ⭐ **performs ABLATIONS**; gates; banks results | ⛔ **Review its own build.** Whatever Claude authored or repaired goes to the fresh reviewer, never back to Claude [`review-agents`]. ⛔ Also: **hand-type a long scaffolding document** — see the Scaffolding-applier row [`thin-orchestrator-definition`] |
| **Ablations** — a review instrument | A mutation whose only purpose is to make a gate fire, reverted inside the same session. Restore by `cp`. ⚠ **Revert it before the session ends** — code that survives the session, or that another party could later adopt, is authorship, and calling it "temporary" does not make it an ablation. | ⛔ **REPLACE the fresh reviewer's per-tooth ablation leg (Phase 4(b), checklist item 6).** An orchestrator ablation is **supplementary**; the fresh leg still owns the verdict. ⭐ **Measured 2026-07-28:** a fresh leg reviewing an orchestrator-run ablation found a real defect in it — the harness restored the mutated **source** between iterations but not the **emitted artifact**, so a stale sidecar decided the runs. ⚠⚠ **The "16 of 22" figure is SELF-REPORTED and repository-unverifiable:** the defective run was gitignored and is absent from the tree (`notes/stage023_step_h_evidence/` holds only `results.tsv` — the *corrected* run — plus `include_list.tsv` and the summary; the per-mutant `cmp_*.txt` are not in git). The defect class is real; the count cannot be checked. ⚠ Two further defects were claimed and are **NOT artifact-supported** (per-mutant captures were in fact retained for all 22 rows; the source restore was present) — the corrected count is **one**, and it is stated here rather than the inflated three because a rule justified by an overstated measurement is the failure this doc exists to prevent. **An ablation is not self-validating: one defect found by an independent reader is the whole point.** ⛔ Never `git checkout`/`restore`/`stash` to undo one on uncommitted work — that restores from HEAD and destroys it [`never-checkout-to-restore-uncommitted`]. ⛔ The target list is the **orchestrator's**, never the builder's [`per-tooth-ablation`] |
| **Scaffolding applier** — a distinct role from a review leg | Applies an orchestrator-authored **decision list** to a **prose** file: a directive, review prompt, plan doc, `STATUS`, or a note section. Flow: Claude writes the decision list (what changes, where, why) → the applier edits the file → Claude **reviews the diff**. Judgement stays with Claude; transcription does not. ⭐ **Edit, never rewrite** a file already written unless the change exceeds half of it. | ⛔ **Touch any code, math, script, or `.wl`/`.py`/`.sh` file — prose only, and never a deliverable.** ⛔ Also be the same session that reviews its own application, or double as a review leg. ⚠ The Clean-agents row's ban on scratch mutation governs **review legs**; it does not reach this role, which exists precisely to write prose. ⚠ This row exists because the old wording (*"Claude owns scaffolding (directive prose…)"*) contradicted [`thin-orchestrator-definition`] and a session then hand-wrote one directive twice in full |
| **Codex** (always `-c model_reasoning_effort=xhigh`) | The default builder: designs + writes + RUNS scripts, iterates to exit 0, applies code/math fixes | Decide the verdict (the fresh reviewer reviews substance afterward) [`codex-xhigh-reasoning`, `codex-iterates-until-clean`] |
| **Clean agents** (the review leg) | The review/verification leg, on its **own fresh agent**, no hint of the expected conclusion; may mutate disposable scratch copies for the orchestrator-assigned per-tooth ablations | **Write or mutate any live/deliverable code/math/script, or make any scratch mutation outside that bounded ablation — an agent is a REVIEW instrument, never a coder.** Also: be reused across legs / carry context between reviews [`review-agents`] |

⭐ **ONE BUILDER; ONE FRESH REVIEWER FOR PROSE AND PROCESS, TWO INDEPENDENT REVIEW LEGS FOR PHYSICS**
(user decision, 2026-07-29/30, amended 2026-07-30). Codex is the default builder, but any party may build;
whoever builds does not review, and every review leg runs on a **fresh** agent with no context from the build.
⭐⭐ **How many legs is decided by the ARTIFACT, not by the phase:**
- **Physics-bearing** — a stage's math, a derivation, a dimensional determination, a verdict a physics claim
  rests on → **TWO independent review legs.** Independent of the builder *and of each other*: neither leg sees
  the other's findings, and a clean first leg does not retire the second. ⭐ Independence is what was measured to
  matter (see *Verification discipline*: a load-bearing "decisive test" was a mean-value-theorem tautology that
  both its author and its executor reviewed SOUND, and only a third independent reader caught).
- **Prose and process** — a directive, a plan doc, `STATUS`, a note section, tooling → **ONE fresh leg.**

⛔ **What this does NOT restore, and what does not come back:** the four-leg Codex→GLM→fold→Codex→Grok bookend,
the round-robin applier rotation, per-chunk user gates, three-session build/fixture/implement sequences,
"acceptance tooling is authored outside the session it polices", and "Claude may never author code". A second
leg is a **second independent reader of the same artifact**, not a second *model in a fixed sequence* and not a
new phase. ⚠ **Two things are kept and are not negotiable:** every review leg is **fresh** (cheap, and it caught
real defects), and the **physics leg stays blocking**.

---

## Standing principles (apply at every phase)
- **Analog, not derivation** — postulate the structure freely; the only test is internal consistency; a **no-go between requirements is a first-class success** [`falsification-is-the-goal`, `analog-find-consistent-structure`].
- **Falsification is the goal** — build every gate *able-to-fail*; a clean "it all works" is the *suspicious* outcome.
- **A miss is a missing-physics signal, never a rescue knob** — when you fall short, do **not** add a free parameter to close the gap (that is the "fits anything" fudge). Re-examine the equations and *derive* the missing term. Reaching for a parameter "to make it behave" is the tell — stop [`falsification-is-the-goal`].
- **Calibrate-then-predict; score the held-out surplus — and count form-shopping as a knob** — calibrating a few inputs on known data is fine; *then* the score is the *new* results you reach with **no further calibration**. The **form** of a law is the prediction; its overall **magnitude** may be a calibration input. Count **all** degrees of freedom honestly — including that ***trying several candidate forms and keeping the best is itself a degree of freedom***. Forbidden: claiming one quantity as both calibrated-to and predicted; using as many knobs as predictions [`calibrate-predict-methodology`].
- **Truth in the OUTPUT, not prose** — no LLM establishes a math or dimensional fact; a verdict-bearing control must compute from inputs and be able to fail [`negative-verdict-short-circuit`, `dimensional-consistency-check`].
- **Never unilaterally deviate from the agreed contract** — the contract is now *one builder, one fresh reviewer, physics leg blocking* (Roles table); if a step is failing, HALT and bring it to the user rather than improvising a substitute. ⚠ The rule polices **silent drift**, not scope: a user-decided scale-back is a change to the contract, not a violation of it [`never-alter-calibrated-process`].
- **"The tool is broken" is a claim that must be VERIFIED before it licenses anything** — a stall/hang/failure attributed to Codex (or Mathematica, or GLM) is a hypothesis, not a fact. Check the actual artifact: does the log end with its completion marker? did the process exit non-zero? does a fresh smoke run reproduce it? **If a tool is genuinely broken, we FIX it — we never route around it.** ⚠ *(Retired clause, 2026-07-29/30: "and never by promoting an agent to coder" — any party may now build; the rule that survives is that whoever builds does not review, and that a broken tool is diagnosed, not silently replaced.)* (2026-07-24: a "Codex CLI stalled twice" claim was written into four docs as justification for an agent-as-coder pivot; on later inspection all 27 of that day's Codex runs had ended `exit=0` and a smoke run answered correctly in 20s. The unverified excuse did more damage than the original error, because it self-justified in every future session.) [`never-alter-calibrated-process`, `negative-verdict-short-circuit`].

---

## The phases

### Phase 0 — Grounding
Read `docs/model_map.md` (⛔ NOT `conceptual_foundation.md` — superseded, it re-confuses) + the relevant memories first. Don't re-derive what's already banked (e.g. the PN ladder, prior no-gos).

### Phase 1 — Directive → ONE review pass (before any compute)
1. **Claude specifies the directive** — *requirement + acceptance + the verdict ladder (incl. the no-gos)*; **never pre-designs the computational route** [`claude-reviews-codex-codes`]. ⭐ For anything longer than a short document, Claude writes the **decision list** and a **scaffolding applier** produces the text (Roles table); Claude reviews the diff. Claude decides *what it says*, not *who types it*.
2. **ONE design-review pass on a fresh reviewer** → fold the findings and go. ⭐ **Edit, never rewrite**: a fold that changes a dozen sections is a dozen edits, not a re-emission of the file. ⛔ **No bookend, no rotation, no second-model confirm pass** — a directive that needs three review rounds to be safe is the wrong shape; decompose it instead [`decompose-before-building-gates`].

### Phase 2 — Per-gate execution prompt
Claude specifies the gate's **execution prompt** (scaffolding; same decision-list/applier split as Phase 1). Review it in the same single pass as the directive when they land together; a separate round is optional, not required.

### Phase 3 — Execution (Codex)
Codex designs + codes + runs **dual-engine** (SymPy + Mathematica), iterating to exit 0 [`dual-engine-required`, `codex-iterates-until-clean`]:
- Mathematica needs **`--sandbox danger-full-access`** to run under Codex [`codex-mathematica-sandbox`]; the orchestrator can run `math` directly as arbiter.
- `codex exec … xhigh` **backgrounded with `< /dev/null`**; **never wrap the Codex session in shell `timeout`** [`background-process-launch`, `never-timeout-codex`].
- The **scripts** Codex runs get `timeout 600`; a 124 = reformulate the math, never raise the cap [`script-timeout-policy`].
- **≤2 concurrent `math -script`** seats; no parallel `$RT exec-*` (it races `MANIFEST.yaml`) [`mathematica-single-seat`, `no-parallel-exec-sympy`].
- ⭐⭐ **CO-AUTHORSHIP GUARD ON THE TWO ENGINES — the party that wrote the `.py` must not be the party that
  adjusts the `.wl` until the comparator agrees.** Chasing green by editing whichever side disagrees collapses
  two engines into one thought, and the result is a comparator that passes because both halves were tuned to
  each other — **exactly the "LLM shortcut that resembles a pass"** *THE POSTURE* names as the operative failure
  mode. A disagreement is a **finding to report**, never a difference to adjust either side into. ⛔ This is
  **not** custody: no commit ordering, no "reference custody", no freezing one engine before the other, and
  nothing here says who commits what when. It is one rule about **who may edit which side while a comparison is
  open** [`dual-engine-required`].
- **Dimensional firewall — comprehensive + inline:** dim-check **every constructed expression as it is built** (not an end-pass), units restored, mutually cross-consistent, **free of any back-solved free carrier**, dual-engine must agree, **able-to-fail** (an ablation must FIRE). No LLM establishes a dim fact [`dimensional-consistency-check`].

### Phase 4 — The backstop: an arbiter re-run plus the review legs
- **(a) Arbiter re-run** — *orchestrator* independently re-runs the existing scripts and every **named acceptance command** unchanged, and reads the literal exit codes (reliability gate; refreshes committed outputs) [`orchestrator-rerun-and-output`, `execute-your-acceptance-commands`].
- **(b) The review leg, doing both halves** — a **fresh clean agent** that (i) reads **code-vs-equations term-by-term** (catches a faithful-but-wrong operator, and the result-emitter that dual-engine agreement is blind to) [`transliteration-fidelity-audit`, `script-review-depth`], and (ii) **tries to break it**: mutates scratch copies, confirms controls are able-to-fail, hunts pass-by-construction [`negative-verdict-short-circuit`]. ⭐ **On a physics-bearing artifact this runs TWICE, on two mutually independent fresh agents** (Roles table): the second leg is launched without the first's findings, and a clean first leg does not retire it. Prose and tooling get one. ⛔ **The builder never selects the ablation target.** Ablation is per-tooth over every able-to-fail check or emitted record; the target list is owned by the **orchestrator**. Otherwise one correctly live record can carry an arbitrary number of hardcoded ones, and every stated observation still passes [`per-tooth-ablation`].
- A negative / "can't-do" verdict gets **harder** verification than a positive [`negative-verdict-short-circuit`].

### Phase 5 — Remediation loop (if a hole is found)
- **Fix the script regardless of outcome** — no spot a future reader can point to as "soft."
- **Codex applies the fix + re-runs** [`codex-is-fix-applier`].
- **Re-verify on a FRESH clean agent (a new fresh view) — never the orchestrator** that just reasoned about the fix [`review-agents`].
- Math-level disagreements: Claude + Codex resolve and agree; escalate to the user only if the fix changes the *conceptual nature* [`claude-codex-resolve-math`].

### Phase 6 — Gate + bank
- **One gate at a time** — finish and bank a gate before starting the next, because a math error cascades and finishing gate 1 can invalidate gate 2's premise. ⚠ **Not a per-chunk user gate** (removed 2026-07-29/30): stop for the user at a *decision*, a blocking finding, or a no-go — not at every chunk boundary [`sequential-audit-chunks`].
- ⭐ **Commit whenever it makes sense — commits are cheap** (user decision, 2026-07-30). ⭐⭐ **ALWAYS commit
  BEFORE a destructive or hard-to-reverse change**: a trim, a deletion, a restructure, a mass rewrite. The
  orchestrator does not wait to be asked. ⚠ The failure this prevents is real and was measured: a ~130-line
  trim was applied to an **untracked** file under the premise "git preserves what you remove" — which was
  **false**, because nothing had ever been committed. Squash small follow-up fixes into the prior commit
  [`squash-followup-fixes`].
  ⛔ Unchanged: **a builder/sub-agent still does not commit** — the orchestrator commits explicit paths after
  reviewing the diff. That rule keeps the diff reviewable and is a different rule from this one.
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
- ⭐ **INDEPENDENCE of the reviewer is what was measured — not the number of models.** On one
  load-bearing "decisive test," the proposed A-vs-B discriminator was *guaranteed* to pass by a
  mean-value-theorem identity — a tautology dressed up as a test — and **both the author (Claude) and the
  executor (Codex) reviewed it SOUND.** A third, independent reader caught it; a review by either party
  that had already reasoned about it would not have. ⇒ The lesson is that a review leg must be
  **independent of whoever built the thing** — and, on a **physics-bearing** artifact, that **one such
  reader is not enough**: this is the measurement the second independent leg rests on (Roles table).
  ⛔ It is **not** an argument for a standing multi-model bookend, a fixed model order, or a confirm pass;
  those were cut 2026-07-29/30 and stay cut. What is restored is a **second independent reader of the same
  physics**, launched without the first's findings [`decisive-test-not-tautological`, `review-agents`].
- ⛔⛔ **CHECK WHICH TOOL REPORTED THE PASS (2026-07-26) — ⭐ A KEEPER, and it passes the governing test
  squarely.** This is not adversary-defence: it is **exactly the "LLM shortcut that resembles a pass"**
  named in *THE POSTURE* as the operative failure mode. Nobody deceived anyone — a build hit a failing
  gate, worked around it, logged the workaround honestly, and then summarised optimistically. That is an
  honest error whose signature is a green headline over a gate that never passed, and it can hide a wrong
  derivation indefinitely. An acceptance item named a specific
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
    rules below. ⚠ **This survives the 2026-07-29/30 scale-back deliberately:** it costs one command and
    catches a class of green that means nothing.
- **When adding a STRUCTURED artifact, verify its STRUCTURE, not only its contents — by running the
  consumer.** Same incident: the orchestrator committed a `.wl` after checking the record count and one
  value, but not that the file had its required header. The consumer would have said so immediately.
- ⭐ **State what would make each acceptance criterion FIRE.** For every criterion, say — in one line, as
  a reasoning exercise — what input it rejects. A criterion with no such input is decoration, and that is
  the defect worth catching: **a check that cannot fail.** *(Directive-review finding, 2026-07-28: an A7
  checker was specified that could inspect only one of three exponent tokens per record, and an A8 checker
  that could compare a derived count against itself.)* ⛔ Do not turn this into an adversary exercise —
  enumerating the cheapest wrong build a motivated builder could ship is the layer this process no longer
  buys.
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

- ⛔ **Decide the model class TARGET-BLIND, and write it down before running.** The pre-registration brief
  states which branch / model class, which targets count as pass and which as fail, which observables and
  *how* extracted, the error budget the result must clear, and what a **trustworthy miss** looks like — and
  it is not edited in response to data. Settle the **forms, domains, tolerances, mesh/grid, calibration
  objective, optimizer, tie-breaker** and *every* discrete family choice **with no reference to any target
  value**. Calibrate **only** the declared parameters to the stated anchor, record the calibrated values,
  then evaluate the held-out observables without touching any of it. **Whatever comes out, stands** — no
  post-residual refit, no moving-boundary bailout added after seeing the residual, no reporting a fresh
  attempt as a "continuation" of a failed one. ⭐ **Extraction is itself target-blind:** the solver and the
  extraction pipeline do not know the targets; targets meet extracted values at a single designated
  comparison step. *Without this, the unconscious tuning loop is invisible from inside your own head.*
  ⭐ **Target-blindness is the load-bearing property** — what must not happen is a choice made *because of*
  a target [`never-supply-the-expected-reason`, `calibrate-predict-methodology`].
  ⚠⚠ **UNRESOLVED, and deliberately not resolved here — the solver's FREEZE HASH.** The 2026-07-29/30
  decision removes immutability mandates from the *process* layer, but on the numerical branch the freeze
  hash is the stated *mechanism* of target-blindness and no-post-residual-refit, not decoration:
  `docs/methodology_paper_outline.md` §2.2/§2.6 (*"the freeze-hash machinery … enforces it mechanically"*),
  `docs/branch_realization_execution_plan.md` §7 and its §"prerequisites", and
  `research/pde_audit/simulation/NONLINEAR_PROTOCOL_V2.md`'s `freeze.target_blind` /
  `freeze.no_post_residual_refit` flags. It is also **already exercised**: the Stage-1 pre-registration
  record was frozen with user sign-off (commit `3865ec9`). ⇒ **The scale-back is NOT applied to it. A user
  decision is needed** on whether the numerical freeze hash stays as the enforcement mechanism (it plausibly
  passes the governing test, since post-hoc refitting is a way the *physics* claim goes wrong) or is
  replaced by discipline alone. Until then, follow the numerical docs as written.
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
   ⭐ for every acceptance criterion, one line saying what input it REJECTS — a criterion with no such
   input is decoration
2. ☐ ONE design-review pass on a fresh reviewer → fold → go (⛔ no bookend, no rotation, no confirm pass)
3. ☐ Builder executes dual-engine (danger-full-access, xhigh, `< /dev/null`, scripts `timeout 600`) → exit 0; comprehensive inline dim-firewall, able-to-fail
4. ☐ ⭐ **Orchestrator re-runs each NAMED acceptance command itself and reads its literal exit code** —
   not the report's claim. Verify structured artifacts by running their CONSUMER, not by eyeballing
   counts and values
5. ☐ ⭐ **The PHYSICS leg, on a fresh agent — BLOCKING, and it does not get cut.** Cross-engine agreement
   is necessary and not sufficient; deriving the result from the model is the leg that establishes
   correctness. Where a workstream defines its own form of this leg, that doc owns it — for the dimension
   rewrite see `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md` §4-(c)/(c1)
6. ☐ Fresh review leg: term-by-term fidelity **and** per-tooth ablation; ⭐ **TWO mutually independent fresh
   legs if the artifact is physics-bearing** (a stage's math, a derivation, a dimensional determination) —
   one for prose/process/tooling; ⛔ ablation target list owned by the **orchestrator** (per-tooth over every
   able-to-fail check / emitted record — never chosen by the builder)
7. ☐ Any hole → the builder remediates → **re-verify on a FRESH clean agent**
8. ☐ Bank: `STATUS.md` + `decisions/13` §0 + memory; commit only if asked

---

*This doc is the operative process; the `feedback-*` memories cited above hold the why. Companion docs: `docs/development_plan.md` (what to build), `docs/model_map.md` (the model — ⛔ NOT `docs/conceptual_foundation.md`, superseded and it re-confuses), `STATUS.md` (where we are).*
