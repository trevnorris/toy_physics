# ▶▶ SCRIPT REBUILD — read this before anything else in v3

**Opened 2026-08-04 by user decision.** ⛔ **All new ledger physics is STOPPED**, including S11b-C, until
this closes. ⭐ This file is the plan; ⛔ it supersedes `steps/S11b_HANDOFF.md` as the read-first document.

---

## ⭐⭐⭐ WHAT WAS FOUND

An engine script can write `emit("TAG", "prose stating a conclusion")` where **no computation in the script
produced that conclusion**. The tag looks like output. It is a typed sentence. ⚠ Cross-engine agreement on
such a tag is **vacuous** — both engines carry the same author's sentence.

Measured with `reduction/derived_or_declared.py` (built 2026-08-04, ⚠ uncommitted at time of writing):

| engine | tags | **genuinely derived** | notes |
|---|---|---|---|
| `S10_two_transverse_photons_sympy_audit.py` | — | ⛔ **unmeasured** | gate cannot parse it: **duplicate emitted tag** `S10_D_TABLE_ROW` — itself a finding |
| `S11_stray_longitudinal_sympy_audit.py` | 79 | 16 | ⚠ only **1 of 6** perturbations informative ⇒ number unreliable |
| `S11bA_interface_response_sympy_audit.py` | 44 | 5 | 6/6 perturbations ran |
| `S11bB_interface_assembly_sympy_audit.py` | 133 | **10** | 13 flagged, ⛔ **3 are symbol-name echo, not derivation** |

⇒ ⭐ **Roughly 10–20% of emitted tags depend on any computation**, across three independently-built steps
with different review histories. ⛔ **This is not one bad script.**

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
> **2. PRINT the residual. Do NOT assert it.**
> `assert residual == 0` **is the builder writing down the expected output**, and it converts an
> informative value into a binary crash. ⭐ Compute → `emit` → *then* assert if you like. The emit must
> come first so the number is on the record instead of in a traceback.
>
> **3. Interpretation belongs to the STEP RECORD, not the script.**
> The orchestrator writes what a printed value *means*, citing it. ⛔ The script does not editorialise.

⭐⭐ **This is a STRUCTURAL control, not a behavioural one.** Blindness is behavioural — quarantine the
file, sandbox the builder, audit the log — and it has failed here repeatedly, including through `git show`.
⇒ **Remove the slot where a typed answer can go and no blindness is required.**

⭐ **The good pattern already exists in our own code, used once in 133 tags** — `S11bB` lines 421–443: type
a candidate `K0_slice`, compute the determinant **independently**, `emit` the symbolic difference, hard-stop
if nonzero. ⇒ That is why `K₀` is the one solid result in the file. **Copy this shape.**

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

⚠ Consequence: this can run **in-repo like any normal build.** ⛔ No quarantine, no sandbox, no
outside-the-repo directory, no memory games. Those defended against the wrong threat.

---

## ▶ THE PLAN

### Order — ⛔ dependency order, not convenience

```
S9 (.wl only)  →  S10  →  S11  →  S11b-A  →  S11b-B
```

⚠ S11b builds on S11, which builds on **S10's mode count**. ⛔ If S10 is wrong the rest inherits it.
⭐ Each step: **both engines**, blind `.wl` first (that rule still stands and is cheap).

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
- **L3** `CONSTANT` conflates *typed prose* with *computed and provably zero*.
- **L4** `MAX_PERTURBATIONS=6`; its selection never reaches `Lambda_*`, `tau_*`, `K_L`, `D_*`, `A_*`.
  ⚠ One slot was spent on `mu_R`, a symbol used nowhere.
- ⛔ **False DERIVED is real** — a tag that merely echoes an input symbol's *name* moves under collapse.
  3 of S11b-B's 13 are this.

⇒ ⭐ **The gate is a triage tool, ⛔ not a verdict.** Its CONSTANT list is **for adjudication.**

---

## ⛔⛔ PROCESS TRAPS — all measured 2026-08-04, ⛔ do not repeat

- ⛔ **Do not build blindness apparatus.** Two rounds of it were designed and discarded in one session.
  The measured failure is **absence of computation**, ⛔ not anchoring on a known answer.
- ⛔ **Do not audit 256 tags.** ⭐ Audit the **claims** the ledger actually makes — the cards, the registry,
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
