# ▶▶ SCRIPT REBUILD — read this before anything else in v3

---

# ⭐⭐⭐ STATE, 2026-08-05 — **S10's ENGINES ARE DONE. THE HARNESS IS BLOCKED ON PLUMBING.**

## ⛔⛔ THREE BLOCKERS, ALL PLUMBING, ⛔ NONE OF THEM PHYSICS — fix these first

The harness cannot yet consume S10's output. ⭐ **Every blocker is an interface mismatch between the
engines and a checker built for S9's simpler shape.** ⛔ No physics result is in question.

1. ⛔ **Mathematica diagnostics reach stdout.** `Solve::svars` prints **10 raw lines** into engine 1's
   output; the harness parses them as tags and dies with *"duplicate emitted tag Solve"*.
   ⚠ **The engine is doing the right thing otherwise** — it already emits
   `..._SOLVE_SVARS_MESSAGE` tags, which is what the repair asked for. ⭐ **Fix: `Quiet` the message at the
   solve while keeping the captured flag and its tag.** ⛔ Do not re-suppress the *information*.
2. ⛔ **`derived_dimensions` cannot take a symbolic-in-`D` payload.** ⭐ Already worked around in
   `checks_S10.yaml` by pointing at the `_SPECIALIZED` tags.
3. ⛔ **`derived_dimensions` cannot take a LIST of coefficient dimensions.** The engine emits one entry
   per coefficient — ⭐ correct, because Q6 was amended so that *which* coefficients exist is
   package-dependent — and the layer wants a single triple. ⇒ ⭐ **Either add a selector to that config
   key, or have the engines also emit one tag per named coefficient.**

⚠ **The last harness run used a FILTERED copy** of engine 1's output (blocker 1 removed by `grep`). ⛔ That
is a diagnostic convenience, ⛔ **not** a result — the committed `out/` file still contains the 10 lines.

## ✅ What IS done

| artifact | state |
|---|---|
| `directives/S10_SHARED_PHYSICS.md` | reviewed (Codex+Grok) + **6 amendments**, all leak-gated |
| `mathematica/S10_brane_mode_spectrum_mathematica_audit.wl` | **2983 tags**, 64 s, exit 0 — built blind, reviewed ×3, repaired ×3 |
| `scripts/S10_brane_mode_spectrum_sympy_audit.py` | **4227 tags**, 150 s, exit 0, **byte-identical reruns** — reviewed ×2, repaired ×2 |
| `mathematica/out/` + `scripts/out/` | ⭐ **engine output is now COMMITTED EVIDENCE** (v2's convention, adopted in v3 at last) |
| `reduction/checks_S10.yaml` | 677 name pairs; `parity_exclude: _LOCAL_` |
| `reduction/engine_output_checks.py` | `parity_exclude` added, parity layer only |

⭐ **Both engines' computed chains are verified LIVE by form ablation on independent legs** — `N3` (stacked
rank), `N7` (two different algorithms), both matrix routes under one-sided corruption **in both
directions**, per-package re-entry at the action, the dimension tree walk, and Q7's two independent sides.

## ⛔ REMAINING, in order

1. **The three blockers above**, then run the harness for real.
2. ⛔ **Verify the harness by ABLATING THE HARNESS** — ⛔ never by reading its self-report.
3. **Write the step record** (orchestrator's job) — ⛔ it should cite harness results, which is why it is
   not written yet.
4. **TeX card** — Codex, with its own two legs.
5. ⚠ **S9 `out/` retro-fit** — S9's engines have no committed output either.

## ⚠ KNOWN LIMITS — ⛔ do not let these pass silently into the record

- ⛔ **Cross-engine coverage is a SUBSET: 677 pairs out of ~2900 candidates.** ⚠ **This is a SPEC gap of
  mine, ⛔ not an engine defect** — §8 pinned the tag **prefix** grammar and never pinned the **quantity
  names**, so both engines obeyed it and still diverged (`Q5_ORIGINAL_ROOT` vs `Q5_ROOT_ORIGINAL`; engine 2
  suffixes dimension info onto object names where engine 1 emits separate families).
  ⇒ ⭐ **Pin quantity names in the spec before S11.**
- ⚠ **Several rounds rested on ONE review leg, not two** — grok was unavailable (per-user lock, and later a
  machine-wide memory exhaustion). ⇒ [[feedback-grok-single-instance-per-user]].
- ⚠ `parity_exclude` matches by **substring**; the builder flagged that as blunt and an anchored
  engine-prefix convention would be safer. ⭐ Agreed, not yet done.
- ⚠ **Runtime margin is thin.** Modest form ablations pushed both engines past 600 s, always in
  `XFORM_ANISO`.

## ⭐⭐ SIX SPEC AMENDMENTS — the durable output of this step, ⭐ and FOUR WERE MY OWN DEFECTS

⚠ Three of the four were introduced **by repairs fixing something else**, and **every one was caught by a
review leg reading the spec end-to-end**, ⛔ never by a build failing.

1. **Phase average** — the period was written with `ω`-dependent `t` limits, which are a real period only
   if `ω` is real and nonzero, **precisely the case Q3 exists to detect**.
2. **`[u] = length`** — never stated; both engines invented it and happened to agree.
3. **Q5 contradicted corollary 4** — it told the builder to emit some tags only when a value existed.
4. **Q6d** — the homogeneity check is **structurally vacuous** and cannot be repaired by a cleverer check
   (`[u]` cancels in every coefficient ratio); ⇒ emit the equation/unknown counts and **label** it.
5. ⭐⭐ **Corollary 5 — A DECLARATION MUST BE WIRED TO WHAT IT DECLARES.** Both engines typed their premise
   tags as second literals; one kept printing its old value while **281 dependent payloads moved**.
   ⭐ **The test: perturb what the tag declares; it must move.**
6. **The run record must be OBSERVED**, not declared before the sweep.

⇒ ⛔⛔ **THE SPEC IS NOW ~490 LINES AND HAS TWICE CONTRADICTED ITSELF FASTER THAN BUILDS SURFACED IT.**
⭐ **Do a consistency pass over it as its own step before S11 inherits it.**

---

## ⭐⭐⭐ WHAT WAS FOUND

An engine script can write `emit("TAG", "prose stating a conclusion")` where **no computation in the script
produced that conclusion**. The tag looks like output. It is a typed sentence. ⚠ Cross-engine agreement on
such a tag is **vacuous** — both engines carry the same author's sentence.

Measured with `reduction/derived_or_declared.py` (built 2026-08-04, ⚠ uncommitted at time of writing):

| engine | tags | gate says DERIVED | notes |
|---|---|---|---|
| `S10_two_transverse_photons_sympy_audit.py` | — | ⛔ **unmeasured** | gate cannot parse it: **duplicate emitted tag** `S10_D_TABLE_ROW` — itself a finding |
| `S11_stray_longitudinal_sympy_audit.py` | 79 | 16 | ⚠ only **1 of 6** perturbations informative |
| `S11bA_interface_response_sympy_audit.py` | 44 | 5 | 6/6 perturbations ran |
| `S11bB_interface_assembly_sympy_audit.py` | 133 | 13 | ⛔ 3 are symbol-name echo, ⛔ and the CONSTANT side over-counts |

⛔⛔ **DO NOT QUOTE A FRACTION. An earlier version of this file said "only ~10–20% of tags depend on any
computation"; BOTH review legs rejected it and they were right.** ⚠ The `DERIVED` column is a **lower
bound** and the `CONSTANT` column is **not a count of typed prose** — the gate cannot see dimension tables,
cannot see through an `assert`, and uses only six collapse pairs.
⚠ **Verified false negative:** `Zperm_difference` at `S11bB:136-139` is **genuinely computed** and emitted
at `:342`, yet the gate calls it CONSTANT — because it is asserted zero, so it always prints `0`.

⭐ **What IS established, by SOURCE READING and confirmed independently by two legs:** specific named tags
at specific lines are typed prose with no CAS object — `S11BB_PORT_DISSIPATIVITY` through
`S11BB_COEFFICIENT_ADMISSIBILITY` (`:348-360`), and the whole transverse block (`:467-471`, where the symbol
`mu_R` appears in **no expression in the file**). ⇒ ⭐ **the defect is located, not quantified**, and it
appears in **three independently-built steps** with different review histories. ⛔ Not one bad script.

⛔⛔ **THE MATHEMATICA ENGINES ARE COMPLETELY UNMEASURED.** The gate is Python-only. S11b-B's `.wl` is
**known** to carry the identical defect (`hInt` at line 86, assigned and never referenced). ⇒ treat the
`.wl` corpus as **unmeasured, not clean.**

### ⚠ How it evaded eight review legs

Every review was a **fidelity** review — *does the artifact match its source?* Under that question a prose
tag is **perfectly faithful**: the card says what the step record says, which says what the engine printed.
⭐ **Every hop is locally correct and the chain certifies nothing.**

⇒ The only leg that found it was one instructed to **ablate** — break something and check the output moved.
⭐ *"Does it say the right thing"* and *"does it depend on anything"* are different questions, and only the
second catches this.

---

## ⭐⭐⭐ THE STANDING RULE — the whole fix in three clauses

> **1. A script may PRINT computed objects. It may NOT state conclusions.**
> `emit()`'s payload must be a CAS object — an expression, a solved root, a boolean from a symbolic test.
> ⛔ Never author prose describing a result.
>
> **2. EMIT BOTH OPERANDS AND THE RESIDUAL, then guard.**
> ⛔ **A residual that is asserted zero always prints `0` and carries no information.** ⚠ **Measured:**
> `Zperm_difference` (`S11bB:136-139`) is genuinely computed, asserted zero, and emitted at `:342` — its
> printed `0` says nothing about whether the hand-typed `Zperm_given` it was compared against is right.
> ⇒ ⭐ **Emit the INDEPENDENTLY DERIVED object itself**, alongside the candidate and their difference. The
> derived operand carries the evidence; the residual only carries the agreement.
> ⚠ A guard is still **correct** — a failed invariant must not let downstream emit garbage. ⭐ Emit first,
> **then** hard-stop. ⛔ An earlier version of this file said *"do NOT assert"*; that was wrong, both legs
> caught it, and it is **retracted**.
>
> **3. Interpretation belongs to the STEP RECORD, not the script.**
> The orchestrator writes what a printed value *means*, citing it. ⛔ The script does not editorialise.

⭐⭐ **This is a STRUCTURAL control where blindness was a behavioural one.** Remove the slot a typed answer
goes in and you no longer depend on hiding the answer.
⛔⛔ **THIS DOES NOT ABOLISH THE BLIND `.wl`, AND AN EARLIER VERSION OF THIS FILE IMPLIED IT DID.** ⭐ Two
engines exist so they can **DISAGREE**; that is about **independent construction**, ⛔ not about hiding a
result. ⇒ **Keep: `.wl` written first, barred from the registry and from the `.py`.**
⇒ **Drop: quarantining answers from the builder.** Those are different mechanisms and the old docs conflate
them.

⭐ **The good pattern already exists in our own code** — `S11bB:421-443`: type a candidate `K0_slice`,
compute the determinant **independently**, `emit` **the derived polynomial itself** *and* the difference,
hard-stop if nonzero. ⇒ ⭐ **That emitted derived object is exactly what `Zperm` lacks**, and it is why `K₀`
is the file's one solid result. **Copy this shape, including the emitted operand.**

---

## ⭐⭐ WHAT TO SUPPLY, AND THE ONE THING TO WITHHOLD

⛔⛔ **DO NOT BLIND THE INPUTS. THIS WAS MY ERROR AND THE USER CORRECTED IT** (2026-08-04):
> *"I feel like it's unreasonable to say, here is only plain english do some complex calculations from it.
> … AI sucks at doing math, even in CAS. It needs some proper hand holding."*

⭐ **Under-specification has cost this ledger far more than contamination ever has.** Measured: `∇·u` fell
out of a spec **four times**, every time it was prose instead of an equation; and S11b-A's first
whole-interface spec described a system with **no linear coupling at all**, which both engines would have
faithfully turned into *"the longitudinal does not couple."*

| ⭐ **SUPPLY, in full, as EQUATIONS** | ⛔ **WITHHOLD** |
|---|---|
| the physical setup, field content, governing equations | ⛔⛔ **any acceptance criterion that references an expected value** |
| supplied premises, explicitly flagged as unfalsifiable-by-this-build | |
| the list of quantities to compute (⭐ a **question list** — leak-free by construction) | |
| the standing rule above | |

⭐⭐ **The risk is NOT a faked computation. It is the "fix until it matches" loop.** Codex iterates to
exit 0; if *"matches the recorded value"* is the exit condition, it will get there, and a genuine
disagreement — **the most valuable output available** — becomes silent confirmation.

⇒ ⭐ **The builder's job ENDS at compute-and-print. The diff happens on our side, afterwards**, where a
mismatch is a **finding**, ⛔ not a build failure.

⛔⛔ **BE HONEST ABOUT WHAT THIS BOUNDARY IS: it is ANTI-TUNING, ⛔ NOT a seal.** Both legs made this point.
A complete equation set **largely determines the answer**, and a builder can hard-code the resulting CAS
expression with no prose and no computation. ⇒ ⭐ **withholding the acceptance criterion stops the
fix-until-it-matches loop and nothing else.** ⛔ Do not treat it as blindness.

⚠ Consequence: this runs **in-repo like any normal build** — ⛔ no quarantine of answers from the builder,
no sandbox, no outside-the-repo directory. ⭐ **The blind `.wl` is retained** (see the standing rule):
independent construction, ⛔ not answer-hiding.

---

## ▶ THE PLAN

⭐ **Budget: under two days** (user, 2026-08-04). ⚠ A review leg estimated "multi-week"; ⛔ **that is wrong
and the reason matters** — the original build's cost was **discovery**: figuring out the process, walking
the physics side by side, eleven directive revisions on S11b-B, two rejected whole-interface designs.
⭐ **All of that is already paid.** The rebuild turns existing equations into computations. ⚠ The only thing
that can blow the estimate is a genuine **disagreement**, and each one is a **finding** — ⇒ the schedule
slips only in proportion to how much is actually wrong, which is the outcome worth paying for.

### Order — ⛔ dependency order, not convenience

```
S9 (BOTH engines)  →  S10  →  S11  →  S11b-A  →  S11b-B
```

⚠⚠ **CORRECTED 2026-08-04 — this line read `S9 (.wl only)` and that was wrong.** It inherited
`LAUNCH_PROMPT.md` OWED 2b, which stamped an orchestrator's decision with the user's name. ⛔ **The user
never made it**, and has never written a line in this repo. ⇒ **S9 gets a SymPy engine like every other
step.** The old exemption's rule was *conditional* — *"a second engine earns its place where the algebra
is long enough that it could genuinely DISAGREE"* — and its condition **expired** when the 26-tag script
became a 316-tag engine. ⇒ [[feedback-attributions-are-my-paraphrase]].

⚠ S11b builds on S11, which builds on **S10's mode count**. ⛔ If S10 is wrong the rest inherits it.
⭐ Each step: **both engines**, blind `.wl` first.

⛔⛔ **THE CHAIN IS NOT CLOSED — these enter from outside it and a rebuild CANNOT falsify them.** List them
as premises in each step record: **S9** consumes the charter and Maxwell (`S9:110`); **S11** adds the
postulated `B_comp` (`S11:118`); **S11b-B** receives supplied affinity, branch prescription, mass balance,
virtual-work rule and row structure (`S11bB record:136`).

### Per step

1. **Leak-gate the existing spec** — `directives/*_SHARED_PHYSICS.md`. ⭐ Strip stated **conclusions**;
   ⛔ keep **every equation**. ⚠ S11b-B's spec went through **eleven revisions**, several written *after*
   engine results existed, so it does contain conclusions. This is a strip, ⛔ not a rewrite.
2. **Build the question list** — what to compute, ⛔ never what it equals.
3. **Build under the standing rule.** Objects printed, residuals emitted before assertion.
4. **Run `reduction/derived_or_declared.py` on the result.** ⛔ It must not be the *only* check — see its
   limits below — but a script that still shows mostly CONSTANT has not been fixed.
5. **Reviewer runs the script and reads the printed outputs.** ⭐ This is now possible because there is
   real output to read. ⛔ Mandate a **FORM ablation** in every review prompt — it is the only thing that
   caught this.
6. **Diff the printed values against the recorded step record, ourselves.** Three buckets:
   - **agrees** ⇒ genuinely verified, for the first time
   - **disagrees** ⇒ ⭐ a real physics finding — ⛔ do NOT reconcile it, record it
   - **not computable** ⇒ the claim was always a supplied premise ⇒ **declare it as one**

### ⚠ Expect the third bucket to be large, and that is an honest outcome

⛔ **Do not try to make all 256 tags derived.** Many are definitions, conventions, or supplied premises.
⭐ The deliverable is a corpus where **what is derived is derived, and what is supplied says so.**

---

## ⚠ WHAT ALREADY SURVIVES — ⛔ do not redo this

- ✅ **`K₀ = B_ρ⁽³⁾ − 2CW₀ + k_W W₀² > 0`, and the slice roots.** Genuinely derived: candidate typed,
  determinant computed independently, symbolic difference printed, hard stop if nonzero. ⭐ Verified by
  direct source reading 2026-08-04.
- ✅ **Onsager route agreement in the `.wl`** (`onsagerRoutesAgree = zeroQ[...]`, line 94) — computed and
  folded into the verdict. ⚠ The SymPy side only describes it in comments.
- ⚠ **The passive/admissibility region is CORRECT but NOT machine-derived.** ⭐ Three independent **hand**
  derivations agree (orchestrator + both review legs, 2026-08-04). For a Hermitian `[[a,b],[b*,0]]`,
  `det = −|b|² ≥ 0` forces `b = 0`, giving the region in three lines. ⇒ ⭐ **a rebuild should be able to
  reproduce this; if it cannot, that is a finding.**
- ⚠ **The transverse null's ARGUMENT is correct** — every scalar invariant contains `∇·u`, and transverse
  projection annihilates `∇·u`. ⛔ But **no engine ran it**, and `mu_R` appears in **no expression** in the
  SymPy engine at all.

---

## ⚠ THE GATE'S OWN LIMITS — established by review, ⛔ do not over-trust it

`reduction/derived_or_declared.py` classifies tags DERIVED/CONSTANT by symbol-collapse perturbation.

- **L1** ⛔ It **cannot** test dimension tags (`S11BB_DIM_*`) — their text comes from a Python int-tuple
  table containing no symbols. Nor `PASS`/`FAIL` homogeneity tags. ⇒ ~61 of S11b-B's 120 are in this class:
  **not evidence either way.**
- **L2** ⛔⛔ **The deep one.** Every check in these engines is `assert`ed *before* the first `emit`, so any
  perturbation strong enough to flip a check **kills the process** ⇒ `SKIPPED` ⇒ no information.
  ⭐ **Clause 2 of the standing rule fixes this at the source.**
- **L3** ⛔⛔ **`CONSTANT` conflates *typed prose* with *computed and provably zero* — a VERIFIED FALSE
  NEGATIVE.** `Zperm_difference` (`S11bB:136-139`) is genuinely computed, then asserted zero, then emitted
  at `:342`; it always prints `0`, so the gate calls it CONSTANT. ⇒ ⭐ **the CONSTANT column is NOT a count
  of typed prose, and must never be quoted as one.**
- **L4** `MAX_PERTURBATIONS=6`; its selection never reaches `Lambda_*`, `tau_*`, `K_L`, `D_*`, `A_*`.
  ⚠ One slot was spent on `mu_R`, a symbol used nowhere.
- ⛔ **False DERIVED is real** — a tag that merely echoes an input symbol's *name* moves under collapse.
  3 of S11b-B's 13 are this.

⇒ ⭐ **The gate is a triage tool, ⛔ not a verdict.** Its CONSTANT list is **for adjudication.**

---

## ⛔⛔ PROCESS TRAPS — all measured 2026-08-04, ⛔ do not repeat

- ⛔ **Do not build blindness apparatus.** Two rounds of it were designed and discarded in one session.
  The measured failure is **absence of computation**, ⛔ not anchoring on a known answer.
- ⚠⚠ **A DELIBERATE CHOICE, ⛔ not an oversight — and an earlier version of this file contradicted itself
  here.** Rebuilding every engine and diffing every printed value **IS** a full-corpus check. ⭐ Both review
  legs argued for narrowing to load-bearing claims only; ⛔ **the user overrode that, with data** — the
  original build took days *including* discovery, so the fuller thing is affordable and buys certainty the
  narrow version cannot. ⇒ **Do the full rebuild. Do NOT re-litigate the scope**, and ⛔ do not quietly
  scale it down after reading the review transcripts.
- ⛔ **What that does NOT license: auditing tags for their own sake.** ⭐ When triaging *existing* output,
  go at the **claims** the ledger makes — the cards, the registry,
  `model_map.md`. Most tags are internal bookkeeping nobody cites. ⚠ v2 died from checks-on-checks.
- ⛔ **A reviewer's finding is not a mandate — verify it yourself.** ⚠ Measured: a leg reported the `K₀`
  boundary as "typed and never computed." ⭐ It is verified; reading lines 421–443 settled it in one pass.
- ⛔ **Launch background jobs with `run_in_background: true`, never a shell `&`** — a `&` inside a
  foreground call leaves the job untracked and no completion notification fires.
- ⛔ **Pick acceptance tags that are genuinely derived.** ⚠ My own gate directive named
  `S11BB_STABILITY_CONDITION` as a must-be-DERIVED acceptance tag; it passes only by echoing a typed matrix.

---

## ⛔ ALSO OPEN, ⛔ not forgotten

- **The S11b retraction is committed** (`3e41b463`) and is **correct on its own terms**: the passive region
  is a classification, ⛔ not a verdict, and A's velocity leak costs a **named reservoir** rather than being
  forbidden. ⚠ Two review findings on it are **folded but unapplied**: (a) `Λ_A⁰ ≥ 0` has a
  **thermodynamics-free** grounding (a pole crosses into the upper half-plane for `Λ_A⁰ < 0`) and was
  over-demoted; (b) the `v₀` caveat is applied **selectively** — it must apply to everything in the
  evanescent regime, including the Debye peak the same card sells.
- **Per-step observational field** (user decision, 2026-08-04): every step states **what measurement bounds
  this, and does it pass.** ⛔ Not yet implemented anywhere.
- **S11b-C** — the non-uniform transverse coupling. ⚠ Its requirements are in `steps/S11b_HANDOFF.md`.
  ⭐ Its pre-registered observational falsifier: *if a slit edge converted an O(1) fraction of a photon into
  the thickness channel, diffraction gratings would not work.*
- **The charge experiments** (user's, 2026-08-04) belong to the **charge sector**, ⛔ not C.
  ⭐ Step zero is pure symmetry: under `w → −w`, `ζ` is odd, `e_W` is **even**, and charge is **odd**
  ⇒ the thickness response to charge is **even in charge sign**. ⚠ `v₀` breaks that reflection, so the
  charge-odd part is suppressed by the drain.
- ⚠ `relations.yaml` is still **empty**.
