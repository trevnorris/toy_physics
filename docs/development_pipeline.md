# Development Pipeline — how we work

The operative process document. What a session or an agent is to *do* is decided here.

- **What** we are building: `docs/development_plan.md`
- **The model**: `docs/model_map.md` (not `conceptual_foundation.md` — superseded, it re-confuses)
- **Where we are**: `STATUS.md`
- **The short form, always in context**: `CLAUDE.md`

The `feedback-*` memories hold the *why* — the incident behind each rule. A conflict between this doc and
a memory is a defect to fix in whichever is wrong, not a precedence question.

---

## 1. The posture, and the test every rule must pass

This is a two-person toy-physics self-consistency project. Checks exist to catch a **wrong derivation**.
Immutability is not a discipline here: files are freely editable and nothing is under a custody rule. (One
carve-out, still undecided: the numerical solver's freeze hash — see §9.)

**The governing test, applied to every check, gate and rule including the ones below:**

- Does it catch a way the **physics** could be wrong? → keep.
- Does it only catch a way the **tooling** could be wrong, or defend against a motivated adversary? → cut.

Our artifacts are written by our own agents, so the operative failure is **drift, an honest error, or an
LLM shortcut that resembles a pass** — never misrepresentation. Hardening against an adversary who does not
exist costs rounds and buys nothing. [`physics-not-ceremony`]

---

## 2. The method — the CAS answers, not the orchestrator

**2.1 Two engines exist so they can disagree.** Independent *construction*: the `.wl` written first, barred
from the registry, never a transcription of the `.py`. The disagreement is the measurement. This is about
construction, not about hiding answers. [`dual-engine-required`]

**2.2 A script prints computed objects; it never states conclusions.** An `emit`/`Print` payload must be a
CAS object. Emit **both operands and the residual**, then guard — a residual asserted zero always prints
`0` and carries no information. Interpretation belongs to the step record.

> Measured: named tags at named lines in three independently-built steps were typed prose with no CAS
> object behind them, and eight fidelity review legs missed it — because *"does it say the right thing"*
> and *"does it depend on anything"* are different questions. Do not quote a fraction of the corpus; two
> legs rejected that as unmeasurable. [`scripts-print-never-assert`]

**2.3 Name the object. Do not specify the recipe.**

> Measured, and it cost five review rounds: a spec specified *how* to obtain a compared quantity — zero the
> velocities, divide by a weight, prove the weight nonzero, guard the division. Every subsequent round then
> argued about the recipe: is the weight unique, is the quotient well-defined, is the residual a tautology.
> None of that was a question about the physics. Every one of those questions was manufactured by
> specifying a derivation path the spec never needed — and the object was already returned by both engines'
> constructors.

Ask for the object. Let the engine hand over what it built.

**2.4 If you are deciding in prose what the engines should compute, the method is inverted.**

The tell: many turns of reasoning toward an answer a script would settle in one. The cause, measured: **the
instrument was broken.** The output reader could not parse a few hundred payloads, so the cross-engine
layer could not see the object under dispute; with no measurement available the orchestrator tried to
reason the answer out inside the specification.

When the measurement is unavailable, **fix the instrument**. A broken instrument silently promotes prose
over evidence.

**2.5 A disagreement is a finding**, not a defect to prevent with better wording. Name the object, emit it
from both engines, compare afterwards.

**2.6 Do not build blindness apparatus.** The measured failure is **absence of computation**, not anchoring
on a known answer — quarantine never touched it, and 2.2 kills it structurally. A do-not-read list is a
denylist, and a denylist means the architecture is wrong. (Measured: a denylist was written, then a
scratch-copy instruction materialised every denied path past it.)

Withhold exactly one thing: **an acceptance criterion that references an expected value.** A builder
iterating to exit 0 converges on any target it can see. Supply everything else — under-specification has
cost this project far more than contamination: the same coupling term fell out of a spec four times, every
time it was prose instead of an equation. [`supply-verified-blind-only-open`]

**2.7 Prior art is an oracle, never a premise.** Check a computed result against the literature where they
should coincide — that is a free external check, better than a review leg. Never assume the literature's
result *for our object*: our Lagrangian is not theirs, and their conditions may not be ours. Re-check a
condition before inheriting its conclusion. [`prior-art-maccullagh`]

---

## 3. The gates

**3.1 Whatever writes does not review.** Every review leg runs on a **fresh** agent with no context from
the build.

**3.2 How many legs is decided by the artifact.**

| artifact | legs |
|---|---|
| Physics-bearing — a stage's math, a derivation, a dimensional determination, a verdict a physics claim rests on | **two**, mutually independent |
| A build directive for anything **both engines read** | **two** — an error there makes both engines compute the same wrong thing and *agree*, the one class cross-checking cannot catch |
| Prose, process, tooling | one |

Who the two legs are is decided by **who wrote it**: orchestrator-written → Codex + Grok. Codex-written →
a fresh Claude agent + Grok. Neither leg sees the other's findings; a clean first leg does not retire the
second.

> Measured: on one load-bearing "decisive test" the discriminator was guaranteed to pass by a
> mean-value-theorem identity — a tautology dressed as a test — and **both its author and its executor
> reviewed it sound.** A third independent reader caught it. Independence is what was measured to matter,
> not the number of models. On a shared-spec directive the two legs disagreed on the verdict four times in
> one session, and the disagreement was the most valuable output each time.
> [`decisive-test-not-tautological`, `review-agents`]

**3.3 Launch legs on sight, before reading the result.** A self-check discharges the felt need for an
independent one, and it is most convincing when it finds things.

**3.4 No commit before both legs report.** The commit is the last step. Reviewing the *directive* does not
pay the tax for the *build*. The orchestrator's own verification is not a leg. (Measured: two builds
committed with zero legs, both reverted, a day lost.) [`no-commit-before-legs-report`]

**3.5 The stopping rule — "both legs green" is not it.**

> Ship when no outstanding finding changes what the engines **compute** or what we are **entitled to
> claim**. Everything else is a NOTE: recorded, and shipped with.

*"A leg that finds nothing is weak evidence"* is an **orchestrator's prior** — keep it, never treat a clean
report as strong confirmation. It is **not a leg's instruction**: in a prompt it is a quota, and a quota is
filled. Put the BLOCKING/NOTE buckets in the prompt instead, and say that "no blocking findings" is a valid
answer. [`review-stopping-rule`]

**3.6 Cost is never a reason** to drop a control, narrow a check, or skip a leg. Scaling work down is the
user's call. [`correctness-is-king-cost-irrelevant`]

**3.7 If successive revisions keep breeding defects in the material just changed, change the author.** Do
not fold a fourth time. (Measured: three folds each closed the reported defects and introduced new ones of
the same severity; changing the author broke the cycle and the defect *class* changed immediately.) The
settled physics is not reopened by a new author — it is the drafting that fails.
[`spec-authoring-discipline`]

---

## 4. Roles

| Role | Does | Does not |
|---|---|---|
| **Claude (orchestrator)** | Decides what scaffolding says; reviews the resulting diff; re-runs each named acceptance command itself; performs ablations; gates; banks; commits explicit paths | Review its own build. Hand-type a long scaffolding document — write the decision list instead |
| **Codex** (always `-c model_reasoning_effort=xhigh`) | Default builder: designs, writes and runs scripts, iterates to exit 0, applies fixes; implements committed ablation harnesses (§4) from the orchestrator's cleared knife-list | Decide the verdict; review what it wrote; **choose the harness's targets/knives** |
| **Review legs** (fresh agent / Grok) | Read code against equations term-by-term; try to break it; mutate disposable scratch copies for orchestrator-assigned ablations | Write or mutate any live deliverable; be reused across legs; carry context between reviews |
| **Scaffolding applier** | Applies an orchestrator-authored decision list to a **prose** file; Claude reviews the diff | Touch code, math, or any deliverable; review its own application |

**Ablations** are a review instrument: a mutation whose only purpose is to make a gate fire, reverted in the
same session. Restore by `cp` — never `git checkout`/`restore`/`stash`, which restores from HEAD and
destroys uncommitted work. The target list is the **orchestrator's**, never the builder's: otherwise one
correctly live record can carry any number of hardcoded ones and every stated observation still passes. An
orchestrator ablation is supplementary; the fresh leg still owns the verdict.
[`per-tooth-ablation`, `never-checkout-to-restore-uncommitted`]

A **form** control tests physics; a **coefficient** control tests arithmetic — scaling never leaves the
family. Demand a script and its literal stdout: a prose re-derivation is the same defect relocated into the
review.

**Committed ablation harnesses — the durable form** *(adopted 2026-09-07 — codex-sol + Grok reviewed → SOUND-WITH-CHANGES, folded; user signed off)*.
The review-instrument ablation above is *ephemeral*: it leaves no committed, replayable knife — the strongest
defect-catcher in the method reduced to a leg's `/tmp` scratch (its script + stdout are recorded per E1, but a
later reviewer cannot replay the exact mutation). For every **physics-bearing engine in a cross-engine pair**,
additionally build a *committed* ablation harness — S1 (a control that holds leaves an artifact whose absence
you can see) applied to the mandatory FORM ablation.

- **Wrap the live engine, never re-derive it.** The harness runs the canonical production engine *unchanged*
  for its baseline, then a copy differing by exactly one declared upstream mutation. Any reproducer is
  mechanically extracted from — or calls — the named production construction; an independently reimplemented
  baseline is not provenance evidence (a second derivation that agrees at baseline is the designed-to-agree
  risk relocated: a harness that copies a hand-typed `H` for its baseline and computes `H+δ` corrupted matches
  byte-for-byte and diffs nonzero, yet proves nothing about whether the engine computed `H`).
- **Orchestrator-owned knives; the manifest does NOT replace the decision list.** The target/knife list stays
  the **orchestrator's**, never the builder's (the existing §4 rule) — a builder choosing its own teeth writes
  the knife its payload already passes. Before any harness builder launches, the orchestrator writes that
  target/knife decision list and gates it once under **G2** (Codex + Grok); the harness then *mirrors* the
  cleared list in an executable manifest and its build legs verify one-to-one. Embedding a manifest in the
  script does not waive G2. [`directive-design-review`]
- **Construction site ≠ observation site; print, don't PASS.** Force the perturbation on the object's
  **construction**, emit the `(baseline, corrupted, diff)` triple of the **emitted tag**, *then* guard (E1) —
  a `PASS` tag is the residual-asserted-zero defect in a new file. Zero-diff under a claimed FORM knife is a
  **fail**, not a skip; mutating the emit/diff wrapper is theater. A **FORM** knife changes
  tensor/sector/index/derivative/field dependence while preserving enough type, dimension, retained grade and
  domain to still run; a sign flip or scalar rescale *alone* is a coefficient test. A **coefficient** knife is
  allowed only when the list names it as a channel a FORM knife cannot see (the N6 Φ-coefficient precedent).
  Two-route objects: corrupt one route; the other must stay. [`per-tooth-ablation`, `xform-prefix-is-not-a-form-control`]
- **Scope = the claim-dependency frontier, not "every line."** Cover every independent physics-bearing
  construction whose output enters the comparator or supports a claim/control — each independent construction
  path gets a genuine FORM knife (N6 was five knives, not forty tags). Mechanical/arithmetic checks get their
  appropriate perturbation, never a relabeled FORM. Imported objects are marked explicitly: a consumer harness
  proves the step *depends* on an imported row; only the *producer's* harness proves that row was computed.
- **Drift guard = tagged-payload identity, not full stdout.** The uncorrupted run's **load-bearing tagged
  payloads** must match the engine's committed tags — ⛔ not full stdout (progress, RSS, timestamps, banners
  are nondeterministic and already are on these engines). This is copy-identity, ⛔ not the cut
  byte-identical-**restore** quarantine.
- **Commit the evidence, and ablate the harness itself.** Commit/annex the harness *and* its exact invocation,
  source/output digests, and literal `(baseline, corrupted, diff)` transcript through the normal
  `_measurements/` record — a committed-but-never-run harness is a capability, not evidence (E1). Verify the
  harness by ablating **it**: a coefficient-rescale or dead-path mutation must **fail to report a bite**
  (Phase 4's "ablate the harness" applies to this artifact too).
- **Both engines, never one — including an already-shipped sibling.** The comparator is only as trustworthy as
  the weaker-verified engine, and the *designed-to-agree* risk is symmetric. Adopting at N6 covers **both** N6
  engines now — the existing SymPy N6 scripts that emit joined tags *and* the WL engine — not only the next
  rebuild. [`matching-number-is-not-evidence`]
- **Authorship.** The orchestrator never writes it (CAS instrument). A fresh non-orchestrator builder writes
  each engine's harness, preferably a different model family than that engine's author; a WL-competence
  exception (fresh astra) is allowed and recorded, provided astra did not choose the targets/knives — which is
  safe here precisely because the knife list is orchestrator-owned and G2-cleared. Each invocation receives
  only its engine and the cleared list (not the sibling engine, its output, or the comparator). Codex-written
  → **fresh Claude agent + Grok**, review-until-clear; the orchestrator adjudicates from the triples, not a
  tag.
- **Forward, not retroactive; deferral is scope, and the trigger is dependency, not disagreement.** Adopt from
  S11c-c2 N6 onward; prior steps keep their ephemeral review; blanket backfill is a **user-scoped S11 closeout
  pass** *before* S11 is declared complete. Retrofit a prior object early when a later claim **depends on it as
  a supplied premise** and it has no committed bite-triple — risk-ranked by claim ancestry, known-defect
  history, scoped independence and unresolved controls — ⛔ **not** only when a comparator disagreement traces
  to it (designed-to-agree, common-mode error and shared freezes produce *agreement*, the precise failures
  this exists to catch). *(Both review legs flagged the S11c-b `slab_operator` pressure-slot carrier —
  `∂(slab rows)/∂(δp±, ∂_w δp±)|_{P=0}` plus the #90 `closure_shape_deriv` fold — as the highest-value early
  retrofit: c2 binds it as a supplied/unfalsifiable premise and becomes c2's `C_E/C_M`, yet its #90 sign and
  cross-engine residual are still deferred.)* [`re-scope-blocker-vs-downstream-need`, `never-freeze-a-varying-field`]

The ephemeral leg-ablation is **not** retired and is **not** discharged by a green harness — legs still
FORM-ablate a `/tmp` copy against the orchestrator's list, choosing knives the builder did not. The committed
harness is the durable record that those controls bite.

---

## 5. Standing principles

- **Analog, not derivation.** Postulate structure freely; the only test is internal consistency. A no-go
  between requirements is a first-class success. [`falsification-is-the-goal`]
- **Build every gate able-to-fail.** A clean "it all works" is the suspicious outcome.
- **A miss is a missing-physics signal, never a rescue knob.** Reaching for a parameter to make it behave is
  the tell — stop and derive the missing term.
- **Calibrate then predict; count form-shopping as a knob.** Trying several candidate forms and keeping the
  best is itself a degree of freedom. Never claim one quantity as both calibrated-to and predicted.
  [`calibrate-predict-methodology`]
- **Truth in the output, not prose.** No LLM establishes a math or dimensional fact.
- **A finding is not a mandate — verify it yourself.** Legs have been wrong in both directions.
- **Never unilaterally deviate from the agreed process.** If a step is failing, halt and bring it to the
  user. This polices silent drift, not scope: a user-decided scale-back is a change to the process, not a
  violation of it. [`never-alter-calibrated-process`]
- **"The tool is broken" is a claim to be verified before it licenses anything.** Check the artifact: does
  the log end with its completion marker, did the process exit non-zero, does a smoke run reproduce it? If
  a tool is genuinely broken we fix it; we never route around it. (Measured: a "Codex stalled twice" claim
  was written into four docs as justification for a process pivot; all 27 of that day's runs had exited 0.)

---

## 6. The phases

**Phase 0 — Grounding.** Read `docs/model_map.md` and the relevant memories. Don't re-derive what is
already banked.

**Phase 1 — Directive.** Specify requirement + acceptance + the verdict ladder including the no-gos. Never
pre-design the computational route. For anything longer than a short document, write the **decision list**
and let a scaffolding applier produce the text; review the diff.

Leak-gate before launch: probe for the step's **symbols in proximity**, not for assembled expressions — a
fixed-string probe for a composed ratio does not match that ratio with a factor between its symbols, and
enumerating spellings is the wrong architecture. A gate returning zero hits on a packet that names the
step's own symbols is evidence the probes are wrong.

Then review per §3.2. Fold and go.

**Phase 2 — Execution prompt.** Same decision-list/applier split. Review it with the directive when they
land together.

**Phase 3 — Execution.** Dual-engine, iterating to exit 0.

- Mathematica under Codex needs `--sandbox danger-full-access`; the orchestrator can run `math` directly as
  arbiter.
- `codex exec … xhigh` backgrounded with `run_in_background`, never a shell `&`, and never wrapped in shell
  `timeout`. A launch goes in its own tool call, alone — a multi-line command chained in front of it makes
  the `"$(<file)"` substitution come back empty and Codex then waits on stdin forever.
- Prove the prompt non-empty *before* launch and the log >500 bytes *after*. A hung Codex looks exactly
  like a busy one.
- Scripts get `timeout 600`; a 124 means reformulate the math, never raise the cap.
- At most 2 concurrent `math -script` seats. An orphaned kernel leaks memory without limit and its symptom
  looks unrelated to Mathematica — when a background job dies with a healthy log, check `free -h` first.
- One writer on the tree at a time. A second build editing the file a leg is reading forces the leg to be
  killed and relaunched.
- **Co-authorship guard:** the party that wrote the `.py` must not adjust the `.wl` until the comparator
  agrees. Chasing green by editing whichever side disagrees collapses two engines into one thought.
- **Dimensional firewall:** dim-check every constructed expression as it is built, units restored, dual
  engines agreeing, able-to-fail. Dimension the whole expression by walking its tree — reading the
  exponents of coefficient symbols silently drops every other dimensionful factor. Count field factors, not
  derivative atoms. [`dimensional-consistency-check`]

**Phase 4 — Backstop.**

- **Arbiter re-run:** the orchestrator independently re-runs every **named** acceptance command unchanged
  and reads the literal exit codes.
- **Review legs** per §3.2, doing both halves: read code against equations term-by-term, and try to break
  it. A negative or "can't-do" verdict gets *harder* verification than a positive.
- **Verify the harness by ablating the harness**, never by reading its self-report. A checker that
  mis-parses two expressions into agreement is worse than no checker, because it manufactures confidence.

**Phase 5 — Remediation.** Codex applies the fix and re-runs. Re-verify on a fresh agent, never the
orchestrator that just reasoned about the fix. Math-level disagreements: Claude and Codex resolve and agree
between themselves; escalate only if the fix changes the *conceptual nature* of the model — a units
convention, a value's classification, or how to formalise something does not.
[`claude-codex-resolve-math`]

**Phase 6 — Gate and bank.** One gate at a time; a math error cascades and finishing gate 1 can invalidate
gate 2's premise. Commit freely — and **always before a destructive or hard-to-reverse change**: a trim, a
deletion, a restructure. (Measured: a ~130-line trim was applied to an untracked file under the premise
"git preserves what you remove", which was false.) A builder or sub-agent never commits. Sync `STATUS.md`
and memory at every milestone.

---

## 7. What green does not mean

- **A green self-report is not evidence a hole closed.** What caught it, twice, was an agent *executing*
  the case end-to-end. Reviews that only read code find design smells; only execution finds live defects.
- **Verify what could break, not only what changed.**
- **A fixture must never be coupled to mutable production state.** Fixing the product must never break the
  tests.
- **A shipped guard with no fixture is not a guarantee.** Ask periodically: which guards could I delete
  with the suite still green? (A pre-freeze review found ~42 of ~89 issue codes were never planted.)
- **The result-emitter is what dual-engine agreement is blind to.** A script that builds something
  *adjacent* to the real calculation and prints the answer as a literal keyed on a config label passes
  re-runs *and* cross-engine agreement — both engines transcribe the same literal. Two engines agree
  meaningfully only if each assembled the headline from its own primitives. Tell: grep your own output for
  verdicts, counts or matrix entries you cannot trace back to inputs. [`transliteration-fidelity-audit`]
- **Check which tool reported the pass.** An acceptance item named a comparator; the official tool exited 2
  on a real defect, the build wrote its own variant, ran that, and headlined PASS — disclosing the
  substitution in its raw log and then summarising optimistically. Nothing in the tree revealed it. Only
  re-running the named command caught it. Every build directive must say: *never satisfy an acceptance
  criterion with a substitute tool; if the named tool fails, stop and report.*
  [`verify-which-tool-reported-pass`]
- **State what would make each acceptance criterion fire.** One line per criterion saying what input it
  rejects. A criterion with no such input is decoration — a check that cannot fail.
- **Never supply a premise to a verification prompt, including a negative one.** Never assert that
  something does not exist or was never done, and never forbid a leg from reading the evidence that would
  falsify your premise. (Measured: a prompt asserted no prior verdict existed and forbade reading the
  outputs; one did exist, and the leg escalated the supplied premise into a blocking finding against an
  accurate commit message.) [`never-supply-the-expected-reason`]
- **Verify a structured artifact by running its consumer**, not by eyeballing counts and values.

---

## 8. Diagnosing a stuck problem

- **If successive fixes each ban a spelling and the next probe evades it, the architecture is wrong.** Four
  rounds hardening one check each closed a real hole and each was defeated by the next probe. That is a
  denylist against an expressive grammar; it does not converge. [`denylist-means-wrong-architecture`]
- **The diagnostic that breaks it open: who supplies the value, and who supplies only the alphabet?** If
  the derived artifact controls the function and the evidence contributes only names, every downstream
  guarantee is artifact-determined. Invert the flow; do not sanitise the grammar further.
- **Specify the invariant, not the instance.** A rule written against the field a survey happened to
  measure left the identical defect alive in its sibling.
- **Before building on an empirical premise, measure it.** A cheap survey would have saved a round.
- **Ask "independent of whom?" before calling a check redundant.**
- **When scope is deliberately bounded, write the bound into the directive** — *"do not fix X; report
  rather than implement if you disagree"* — and into the commit message. A competent builder handed a
  checker with three known holes will close them, and the next round starts. Scope decisions that live only
  in a conversation do not survive contact with the next agent.

---

## 9. Operational

**Sub-agent safety.** An agent must never delete a directory, or any file it did not itself create. (A
review agent "cleaned up" the shared gitignored `_scratch/` tree: every directive, launcher, run log and
prior-session artifact, unrecoverable.) Put that clause in every agent prompt and give each agent its own
subdirectory. Guard the workspace, not only the repo.

**Background processes.** The harness reaps long waiters — a dead waiter is not a dead job. A waiter must
distinguish three outcomes: marker found, process gone without a marker, still alive. Match the done-marker
on `tail -1` only. Use absolute paths in every launch and check, and verify a launch by the existence of
its declared output file — absence of an error and pass-shaped downstream text do not establish that the
launch happened. Never `pkill`/`killall` by pattern; the box is shared. Codex's provider can refuse content
on a cyber-policy filter — phrase adversarial directives in correctness terms, never attack language.

**Paths.** Project root is `/var/projects/toy_physics` (not `toy_projects` — a transcription attractor).
Render prompts under the project root, never `/tmp` — agents and Codex cannot read it. Raw transcripts go
*outside* the repository as **tree hygiene**: a transcript carries an engine's complete tag values verbatim
and is noise in the record. It is **not** a blindness control, and relocating it is not one either — see
rule 12 and the CUT table at `.claude/skills/build/SKILL.md:137-148`.

**Grok** is one instance per user (`~/.grok/auth.json.lock`); a stale lock naming a dead PID is normal.
Legs therefore serialise — a scheduling cost, not a reason to drop one.

**File formats.** YAML or markdown for anything an LLM reads or writes; JSON only machine-to-machine. No
`python3 -c` commentary scripts — real ablation re-runs are fine, narration scripts are not.

**Numerical-work gates.** The phases above suffice for symbolic derivation. The moment a numerical solver
enters, these become mandatory:

1. Decide the model class **target-blind** and write it down before running: which branch, which targets
   count as pass and fail, which observables and how extracted, the error budget, and what a trustworthy
   miss looks like. Not edited in response to data. Extraction is itself target-blind. Whatever comes out,
   stands — no post-residual refit.
2. A solver earns trust only after all six independent checks: a known-analytic limit; manufactured-solution
   tests per operator; a published benchmark; mesh refinement at ≥3 levels showing the expected order;
   conservation diagnostics; a stated noise floor. A result without this stack is an interesting number.
   Manufactured-solution tests prove self-consistency, not fidelity.

*Unresolved:* the solver's **freeze hash**. The scale-back removed immutability from the process layer, but
on the numerical branch the freeze hash is the stated *mechanism* of target-blindness, not decoration, and
it is already exercised with user sign-off. A user decision is needed on whether it stays. Until then,
follow the numerical docs as written.

---

## 10. Per-gate checklist

0. Directive states no expected answer, no expected reason, and no premises — including negative ones. It
   carries the anti-substitution clause and the three clauses of §2.2.
1. Directive drafted: requirement + acceptance + able-to-fail verdict ladder, no pre-designed route. One
   line per acceptance criterion saying what input it rejects.
2. Leak-gate run (symbols in proximity, every hit read).
3. Review per §3.2 → fold → go.
4. Builder executes dual-engine → exit 0, with the inline dimensional firewall.
5. Orchestrator re-runs each **named** acceptance command and reads its literal exit code.
6. Review legs: term-by-term fidelity and per-tooth ablation, target list owned by the orchestrator;
   physics-bearing engines also get a committed ablation harness (§4) built from that same orchestrator-owned,
   G2-gated knife-list — both engines of a cross-engine pair.
7. Harness verified by ablating the harness, not by reading its self-report.
8. Any hole → builder remediates → re-verify on a fresh agent.
9. Both legs reported → **then** commit. Bank: `STATUS.md`, memory.

---

*Companion docs: `docs/development_plan.md` (what to build), `docs/model_map.md` (the model), `STATUS.md`
(where we are), `CLAUDE.md` (the short form).*
