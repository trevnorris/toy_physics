# Independent review — the `C17`/`C18` spec-repair decisions, and the work list they sit in

## Artifacts under review

1. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions.md`
2. ⭐ **The work list itself** — the `WHAT S10 STILL OWES` table in
   `/var/projects/toy_physics/research/pde_ledger_v3/REBUILD_HANDOFF.md`, and the ordering that puts this
   spec repair **before** the S11 engines are built.

Both are **orchestrator-written**. Nothing is applied. You are one of two independent legs; the other is
not visible to you.

⚠ `S11_SHARED_PHYSICS.md` is a **shared spec**: an error in it makes **both** engines agree on the same
wrong thing, or **disagree for a reason that is not physics**. That is why it gets this gate.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the decision list first

1. `research/pde_ledger_v3/DEFECT_REGISTER.md`, entries `C17` and `C18` **in full**.
2. `directives/S11_SHARED_PHYSICS.md` — `§5`'s locus protocol (`:230-251`), `§Q8a`/`§Q8b` (`:560-641`),
   and `§Q3`. ⭐ Enough of `§3` to know what the allowed region is.
3. The committed outputs the register cites: `mathematica/out/S11_stray_longitudinal_mathematica_audit.out`
   and `scripts/out/S11_stray_longitudinal_sympy_audit.out`, at the `XFORM_EXTRA, D = 2` stratum ordering.
   ⭐ **Open them; ⛔ do not take the register's quotation on trust.**
4. **Write down, before step 5:** how *you* would make two independently built engines emit **comparable**
   stratum information, given that one CAS can decide real-domain questions the other cannot. ⭐ Keep it.
5. **Only now** the decision list, then the `WHAT S10 STILL OWES` table.

## ⭐⭐⭐ THE GOVERNING QUESTION THIS ROUND — ⛔ answer it for EVERY item

> ⭐ *"make sure all changes we pursue are physics related and not process related"* — **user, 2026-08-10**

⇒ For **each** decision `S1`–`S5`, and **each row** of the `S10 STILL OWES` table, say which:

- ⭐ **PHYSICS** — it catches a way the physics could be **wrong**, or a way a claim would be
  **unsupported**; or
- ⛔ **PROCESS** — it is tooling, bookkeeping, or tidiness.

⛔⛔ **Attack the classification in BOTH directions, and the second is the dangerous one:**
- ⭐ Is anything kept that is really ceremony? ⇒ it should be cut.
- ⛔⛔ **Is anything CUT that is really physics?** ⚠ The table cuts the export regeneration, `C19`'s rename
  and disclosure, `F6`, and the registers. ⭐ **For each cut, try to construct a case where that omission
  lets a WRONG NUMBER or an UNSUPPORTED CLAIM reach a step record, a paper card, or a later step.** ⛔ If
  you can, the cut is wrong — name the case.

## ⭐⭐ THE ORDERING IS ALSO IN SCOPE — ⛔ not settled background

The plan repairs `C17`/`C18` **before** either S11 engine is built. ⚠ The stated basis is `C18`'s own
register text: *"one engine may omit an allowed, physics-changing component the other retains, and both
are faithful to the words."*

⭐ **Is that right?** The alternative is to build both engines against the spec as it stands and treat any
divergence as a **finding to adjudicate** — which is what *"a disagreement is a finding"* would ordinarily
demand.
⛔ Settle it on this question: **is the measured `XFORM_EXTRA, D = 2` divergence a finding about the
PHYSICS, or about the PROSE?** ⭐ Open both outputs and say which, with the evidence.

## What to check in the decisions

### 1. ⭐⭐ Does `S1` actually name an OBJECT, or has it smuggled in a recipe?

`S1` says the stratum's `Q3`/`Q4` are *"the objects obtained under that stratum's defining equations."*
⭐ Is that an object two independently built engines will compute **the same**, or does it hide a branch
choice, a solve order, or a simplification that differs between CASes?
⛔ **Test it, do not reason about it.** Take one stratum from the committed outputs, compute its `Q3`
objects under the defining equations in **your own script**, and report what you get. ⚠ If the component is
positive-dimensional and imposing its equations still requires solving for some variables in terms of
others, ⭐ **say what that does to comparability.**

### 2. ⭐⭐ Does `S1` actually fix `C17`, or relocate it?

⭐ The register measured that two allowed points on **one** component give different root counts.
⛔ Under `S1`, does that difference **disappear**, or does it reappear as a **sub-locus** inside the
component that `S1` now fails to report at all? ⚠ The second would be a **new** blindness — ⭐ a place where
the mode count changes and no tag says so. **Construct the case.**

### 3. ⭐ Is `S2`'s fallback a loophole?

`S2` lets an engine evaluate at its point provided the tag says so. ⛔ Does that make `S1` optional in
practice — both engines take the fallback, and `C17` survives with a label on it? ⭐ If so, say what would
force the stronger object.

### 4. ⭐⭐ Is `S3`'s canonical object real?

`S3` requires one engine-independent characterisation of a locus, with variable and monomial order pinned
in the spec.
⭐ **Name a concrete candidate and test it in BOTH CASes** on the `XFORM_EXTRA, D = 2` equations. ⛔ Does it
exist for the systems this step actually produces — including the empty and identically-satisfied cases the
protocol already separates? ⚠ What happens when the equations are **not polynomial**?
⛔ If no such object survives contact with the real systems, `S3` is unbuildable and must be replaced.

### 5. ⭐ Do `S4`/`S5` close the measured divergence, or only label it?

⭐ Apply them to the measured `XFORM_EXTRA, D = 2` case. Under `S4`+`S5`, what does each engine emit, and
what does the comparator then report? ⛔ Is the outcome a **finding**, or is it *"not compared"* — and is
*"not compared"* acceptable for a region where the **mode count** may move?

### 6. ⛔ Leaks, and one open question about the register itself

⛔ Does the decision list state or imply what any computation will produce? ⚠ **The test is derivability,
not literal presence.**
⭐ **Then a question the orchestrator has NOT decided:** `DEFECT_REGISTER.md#C17` contains measured root
counts for this step's own system, and a builder can open it. ⛔ Is that a leak that matters, given that
rule 12 forbids blindness apparatus and a do-not-read list is a denylist? ⭐ Give a recommendation and the
reasoning, ⛔ not a rule.

### 7. ⭐ What is missing?

⭐ Compare your step-4 list. ⚠ Does the *"does not decide"* section exclude anything **load-bearing** for
`S1`–`S5`?

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, a way a **claim would be
unsupported**, or a way these decisions would cause **wrong or incomparable artifacts** to be built.
⛔ Not style, ⛔ not naming taste.

## Method

- ⭐⭐ **Quote both sides of every finding**: the list's text, and the source it fails against.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
- ⭐⭐ Checks 1, 2, 4 and 5 are settled by **running code**. Write the script, run it, paste the **literal
  stdout**. ⛔ A prose derivation is discarded — it is the same defect class this rebuild exists to remove,
  relocated into the review.
- ⛔ Wrap **every** Wolfram kernel run in `timeout 600`; a 600 s hit is a **failed** check — report it and
  move on. ⛔ Never raise the timeout, ⛔ never more than one kernel at a time (the licence has two seats).
  ⛔ Copy anything you execute to `/tmp` and run the copy.
- ⛔ Do **not** edit anything in the working tree. Read-only.
- ⭐ End with: your `PHYSICS`/`PROCESS` verdict per item; **whether the spec repair should precede the
  engine builds, and the evidence**; which of your step-4 items the list handles and which it misses; and
  what you checked that could have failed and did not.
