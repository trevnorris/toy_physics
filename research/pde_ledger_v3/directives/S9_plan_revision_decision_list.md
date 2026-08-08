# Decision list — revise `HARNESS_S9_PILOT_PLAN.md`

**You are applying decisions, ⛔ not authoring a plan.** Every item below is already settled by two
independent review legs that produced executed evidence. ⛔ Do not re-litigate them, ⛔ do not improve the
plan beyond this list, and ⛔ do not soften any retraction.

**Target file:** `/var/projects/toy_physics/research/pde_ledger_v3/HARNESS_S9_PILOT_PLAN.md`
**Touch nothing else.** ⛔ No other file may be created or modified.

**Report when done:** ≤40 lines — one line per decision (`Dn: applied at §X` / `Dn: already present`) plus
the file path. ⛔ No summary of the plan, ⛔ no commentary on its merits.

---

## ⛔⛔ Standing prohibitions — these govern every decision below

1. ⛔ **Never introduce a physics value, formula, expected verdict, cardinality, or measured outcome.** The
   revised document is read by whoever builds the comparator; rule 5 — *a builder iterating to exit 0
   converges on any target it can see.*
2. ⛔ **Where a decision says "remove a value," remove it — ⛔ do not restate it, and do not paraphrase it
   into a nearby sentence.** A prohibition leaks as surely as an assertion.
3. ⛔ **Do not add new checks, layers or criteria of your own.** If a decision seems to need one, apply the
   decision as written and say so in your report.
4. ⭐ Preserve the document's existing marker conventions (⭐ ⛔ ⚠ ⇒) and its retraction blocks verbatim
   where they already exist.

---

## Part A — verify work already applied (⭐ confirm complete, ⛔ do not redo)

Several decisions were partly applied already. For each: **check it is complete and internally consistent,
and finish it if not.** Report `already present` or `completed at §X`.

**D1 · L3 is demoted to a unit-covariance check.** Confirm every fabrication/provenance claim for L3 is
gone — including in **R4 of the risk table**, which previously named L3 a fabrication defence. Add two
counterexample classes still missing: (a) **`0` satisfies the scaling law for every declared exponent**;
(b) **a matrix or coefficient list may legitimately carry different dimensions per entry, so one scalar
exponent is the wrong covariance law**; (c) a `Piecewise`/conditional object needs branchwise dimensions
plus dimensionally valid conditions, not one global test. ⭐ Rename the layer heading to
**"unit-covariance check."**

**D2 · The L4 leak is repaired.** Confirm no physics value survives anywhere in the layer.

**D3 · Witness checks live in the harness, not the engine.** Confirm the engine emits operands and
residuals and the harness renders every verdict. ⛔ Confirm the expected multiplicity value is gone from
the witness table.

**D4 · L6 is inverted.** Confirm: symbolic residual equality is **primary** for algebraic scalars; exact
evaluation is **secondary and differential**; probe points are drawn **after** implementation and re-drawn
per adjudication run; degree/pole/domain bounds are derived and enforced **independently of the party they
gate**; and nothing weaker than the current comparator may be adopted on any row.

**D5 · A1 is replaced by per-row equivalence.** Confirm it is stated as a **runnable** criterion: per-row
verdict equivalence against the current symbolic comparator, over baseline **and** an adversarial mutant
set, with the row inventory pinned by an artifact independent of the candidate comparator.

**D6 · L8 cross-step consistency exists.** **D7 · §3c builder/oracle separation exists.**

---

## Part B — new decisions to apply

### D8 · Add a FOURTH cost layer to §2: semantic binding / domain

§2 currently names three negotiation layers (tag join · symbol naming · shape selection). Both legs found a
fourth that the plan's own L6 table already depends on. Add it as a row, comprising: **coordinate/ambient
order · coefficient field · parameter premises · exceptional rank strata · branch conditions and sheet
metadata · pole avoidance and degree bounds · root multiplicity · defining polynomial.**

⭐ State plainly: these are **not shape properties**, so L0's object `kind` does not remove them — **L0 must
also carry semantic symbol roles, domains, coordinate bases, branch data and multiplicity semantics.**

⚠ Add: the S9 `omega2 = omega**2` case is an **algebraic identity, not a spelling exception**. Deleting the
symbol-identity machinery does not remove it; it **relocates** into the symbol→value map or the object
schema. ⇒ L0 may trade selector negotiation for **kind-taxonomy** negotiation; say so, and require the
fourth layer's cost to be **measured separately** rather than assumed to vanish.

### D9 · Shared-spec blindness becomes a DECLARED STANDING HOLE, not a study outcome

Both legs walked L0–L7 against a defect introduced **identically into both engines from the shared
specification** and found **no layer fires**. Codex exhibited a passing instance: a **dimensionless
coefficient error in the shared construction** — it keeps correct kinds, shapes, cardinality and
multiplicity; is dimensionless so L2/L3 are blind; preserves the structural scaling behaviour L4 tests;
makes every emitted root a genuine root of the wrong object so L5 is self-consistent; and both engines are
identical so L6 cannot fire.

⛔ **Describe the class structurally. ⛔ Do NOT write the mutant expression or any determinant into the
document.**

Apply:
1. ⭐ State the hole **in §3, before the layer list** — ⛔ not as a possible A10 finding.
2. ⭐ Require **at least one physics-bearing oracle whose construction does not come from the shared
   spec**: independent re-derivation from the action, a separately commissioned metamorphic invariant, or
   an external published result consulted **only after** computation.
3. ⛔ **An uncaught priority shared-spec mutant blocks any coverage reduction** — it is not merely
   published.
4. ⛔ Remove any claim that the pilot "proves the harness" or closes S9's measured structural miss.

### D10 · Provenance mutation must be at SOURCE level

§3b's study mutates committed **outputs** in memory. ⛔ That representation **cannot distinguish** a value
computed from the model, the identical payload typed literally, a premise used in construction, and the
same premise merely printed afterward — all four produce identical output.

Apply: the **fabrication** and **inert-premise** classes must be tested by mutating **engine source** (or
an instrumented derivation IR), ⛔ never output mappings. Cite the existing honest instrument —
`reduction/derived_or_declared.py:318` labels a tag derived if a selected perturbation changes its text —
**and its own recorded caveat at `:332-337` that constant status is not proof of literal source.** ⭐ Note
the converse failure too: a typed **symbolic** literal can change under perturbation and be misclassified
as derived.

### D11 · Harden the §3b reduction gate

As written it is satisfiable as a formality: one mutant per layer, all mutants visible, no kill threshold,
classes permitted to map to "no contract fires," and a minimal hitting set that can overfit the seeded
catalogue.

Apply, as explicit requirements: **multiple mutants per class across different object kinds** · **coupled /
common-mode mutants** · **held-out mutants not visible to the builder** · **deliberately weaker comparator
implementations as negative controls** · a stated **kill threshold** · and ⛔ **any uncaught priority mutant
blocks reduction.** ⭐ Also require the catalogue to cover the shared-spec class (D9) and S11b-like object
kinds, or else carry an explicit **"not yet tested ⇒ no reduction for those steps."**

### D12 · Extend the A6 seed catalogue, and resolve its conflict with R3

Add to **must-fire**: dropped dimensioned factor · wrong root multiplicity · a typed bare literal carrying
a non-zero claimed dimension · one-sided fabrication against a computed counterpart · **wrong operator that
shares a kernel** (must fire the operator comparison, which L6 already requires but A6 does not test) ·
pole/branch error · compensated error · missing duplicate.

⚠ Resolve the conflict: A6 says a **relabelled symbol must not fire**, while R3 says both engines must
accept the same symbol→value map **or fail loudly**. ⭐ State that pure alpha-renaming may be transparent
**only after semantic symbol roles exist** (D8); until then a relabel legitimately fails loudly.

### D13 · Harden A2, A3, A4, A7 — and reclassify A8 and A10

Each currently has a wrong implementation that passes:
- **A2** is trivially satisfied by **deleting the naming section** — the ablator iterates only declarations
  still present (`reduction/measurements/declaration_load_ablation.py:91`), so zero declarations gives zero
  load-bearing ones. ⇒ require the measurement against the **declaration set as it stands today**.
- **A3** — nothing independently verifies that the layer a row *claims* decided it actually did. ⇒ require
  independent confirmation of the layer label.
- **A4** — "is compared" is satisfiable by an always-true or projection-only predicate. ⇒ require the
  comparison to be by the kind's declared mode.
- **A7** — "measured against 3,121" has **no semantic-preservation criterion**; a minified or row-discarding
  config passes. ⇒ require preservation of the compared row set.
- ⭐ **A8 and A10 are REPORTING requirements, ⛔ not falsification criteria.** Relabel them so they cannot be
  counted as gates. Disclosure alone must not be able to pass an adoption decision.

### D14 · Restate the pilot as PHASE 1 of two — S9 alone is a smoke test

S9 is the right cheap falsifier for the naming question and ⛔ nothing more: 12 pairs, one unparsed
conditional sign, few object kinds. The classes it cannot exercise are **live in the tree today**
(`scripts/S11bA_interface_response_sympy_audit.py:37` unevaluated integrals; `:336` transcendentals;
`scripts/S11bB_interface_assembly_sympy_audit.py:29` sheet-sensitive algebra;
`mathematica/S11bB_interface_assembly_mathematica_audit.wl:166` monodromy and non-fixed multiplicity).

Apply: **Phase 1 = S9.** **Phase 2, required before any method adoption or coverage reduction**, adds
representative rows for: an S10 **matrix**, an **action/derivative** object, a **relation**, a **mapping**,
a **multiset with multiplicity**, a **parameter-dependent subspace**; and S11b's **integral**,
**transcendental**, **pole**, **branch-condition**, **sheet** and **textual** cases. ⭐ Add the synthetic
adversaries by name: **probe-set vanisher · pole hit · branch collision · wrong operator with the same
kernel · missing duplicate · symbol relabel · typed literal · inert premise.**
⛔ **State explicitly: no method adoption and no comparison reduction follows from S9 green alone.**

### D15 · Correct a count

S9 declares **six** naming exceptions and **one** symbol identity (`reduction/checks_S9.yaml:27` and `:38`).
Any text saying "7 naming exceptions" is wrong — fix every occurrence.

### D16 · Record the authorship change in §5

Add to the sequencing table and its notes:
- ⭐ This revision was **applied by Codex from an orchestrator decision list** ⇒ per rule 7, the two review
  legs for **this document** are now a **fresh Claude agent + Grok**, ⛔ not Codex + Grok.
- ⛔ **Step 1's `HARNESS_STANDARD.md` is not authored by the orchestrator** — three consecutive
  orchestrator-written artifacts each carried a central defect (rule 15).
- ⛔ **The Codex session that edits this plan must not be the Codex session that later builds the
  comparator.**

### D17 · Elevate textual emissions from a risk row to a named finding

R8 currently notes that natural-language string emissions exist. ⭐ Promote it: **an emission that is not a
parseable CAS object is an ENGINE DEFECT**, comparable by no mode, and it **lands on the S11b rebuilds** as
work — ⛔ not on the comparator as a gap to accommodate.

---

## Acceptance — how the diff will be checked

⭐ I will read the diff, not the document. It must show: **D1–D7 confirmed or completed · D8–D17 applied ·
⛔ zero new physics values · ⛔ no other file touched · ⛔ no decision re-argued.**
