# ▶▶ SCRIPT REBUILD — read this before anything else in v3

---

# ⭐⭐⭐ STATE, 2026-08-04 — **S9 IS DONE. S10 IS NEXT.**

⭐ **The proven per-step order is `.claude/skills/step-run/SKILL.md` § THE PROVEN ORDER.** ⛔ Do not
reconstruct it from this file; that table is the operative one and it is 13 numbered steps with the trap
each exists for. ⭐ `build` and `review-legs` carry the rest.

## What S9 produced

| artifact | state |
|---|---|
| `mathematica/S9_light_requires_shear_mathematica_audit.wl` | **1559 tags, 0 typed conclusions**, 9 actions |
| `scripts/S9_light_requires_shear_sympy_audit.py` (+ `.premises`) | **635 tags, 0 typed conclusions**, 9 actions |
| `reduction/engine_output_checks.py` + `checks_S9.yaml` | 4-layer harness; **cross-engine `agree=12, disagree=0`**, dimensions **1219/1219 homogeneous** |
| `steps/S9_light_requires_shear.md` | rewritten with a full verification section |
| `paper/steps/S9_light_requires_shear.tex` | rebuilt from the record |

⭐ **The physics did not change.** Every value the pre-rebuild engine computed is identical. ⛔ No result
was revised, no sign flipped. **What changed is that each claim now survives something able to refute it.**

## ⭐⭐ THE ONE PHYSICS RESULT THAT CHANGED WHAT S9 CLAIMS

⛔⛔ **Transverse propagation does NOT require the curl-only form.** The gradient-elastic control — an
ordinary elastic solid — carries **two transverse modes at the same `c² = μ_R/ρ_br`**; it merely *also*
propagates the longitudinal. ⇒ ⭐ **curl-only buys the ABSENCE of a propagating longitudinal mode**,
⛔ not the presence of the transverse ones. ⚠ **Say "buys", ⛔ not "uniquely buys"** — five stiffness
forms were tested, which ⛔ does not establish uniqueness across all operators. ⚠ Prose reading *"light requires shear"* as established
by this computation **overclaims** — the statement is conditional on the stiffness form.

## ⛔⛔ CARRY THESE INTO S10 — each cost a repair round on S9

1. ⛔⛔ **CHECK S10's MODE-COUNT METHOD FIRST, BEFORE ANYTHING ELSE.** `M·T = 0` is **not** a transverse
   existence test — it demands `M` annihilate the *whole* transverse subspace and returns a **false
   negative** when a mode's partner sits at a different frequency. ⭐ Use **corank of `M` stacked on `kᵀ`**.
   ⚠ **S9's original engine AND its first rebuild both used the broken test**, and ⭐⭐ **S10's entire
   result is a mode count.**
2. ⛔ **Never test homogeneity by substituting a symbol** (`k² → q`). It silently no-ops when the root is
   not a bare multiple of that sum, giving a false *positive* on a gapped root and a false **negative** on
   a genuinely dispersive one. ⭐ Test by **scaling**: `ω²(λk) − λ²ω²(k)`.
3. ⭐ **Dimension the WHOLE expression tree**, ⛔ never by reading coefficient exponents — that drops
   explicit wavevector factors. And **count FIELD FACTORS**, ⛔ not `Derivative` atoms, or a bare
   undifferentiated field loses its dimension silently.
4. ⭐ **Emit the SIGN of `ω²` per root.** Without it, exponentially growing modes carry the same
   polarisation signature as propagating waves.
5. ⛔ **Emission must never be conditional on a payload's VALUE.** ⚠ A dedup applied across control
   packages deleted correctly-invariant tags; *absent* is indistinguishable from *never computed*.
6. ⭐ **Give the dynamical matrix two independent routes**, and prove independence by **one-sided
   corruption** — break one, the other must not move.
7. ⚠ **`MatrixRank` returns the GENERIC rank.** Emit the exceptional loci explicitly and check the
   enumeration is complete.

## ▶ S10 setup notes

- ⭐ **S10 already has both engines** — ⛔ no exemption question arises. Rebuild both.
- **Author `reduction/checks_S10.yaml`** — ⛔ tag **names** only, never values.
- ⚠ **Move `_SYMBOL_ALIASES` out of `engine_output_checks.py` into the per-step config.** Symbol
  spellings are step-specific; a hardcoded module dict needs a code edit per step. It maps names to
  names, so moving it is leak-free. ⚠ It already caused S9's only **false-positive `DISAGREE`**
  (`lambdaScale` vs `lambda_scale`) ⇒ ⭐ **a `DISAGREE` needs adjudication too**, exactly like `INVARIANT`.
- ⚠ Known-open on S9, ⛔ not blockers: the assumption set is **computationally inert** (physics bounded
  instead by generic rank + a verified-complete locus enumeration); one tag is `UNPARSED` (a `Piecewise`
  sign); the `.py` was written **after** the `.wl`, so its agreement is weaker than blind-first ordering.
  ⭐ At S10 both engines predate any comparison, so that caveat disappears.

---

**Opened 2026-08-04 by user decision.** ⛔ **All new ledger physics is STOPPED**, including S11b-C, until
this closes. ⭐ This file is the plan; ⛔ it supersedes `steps/S11b_HANDOFF.md` as the read-first document.

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
