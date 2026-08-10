# Independent review — the S11 spec-repair decision list, before any repair is applied

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_spec_repair_decisions_v2.md`

It is **orchestrator-written**. It is the list a builder will apply to a **914-line specification that both
S11 engines read**. Nothing has been repaired yet. You are one of two independent legs; the other is not
visible to you.

⚠ **A spec both engines read is physics-bearing: an error in it makes both engines agree on the same wrong
thing.** A defect in this list therefore propagates into two engines, their comparator, and the step
record — which is why it is reviewed before it is applied.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the decision list first

For a document, blindness comes from **reading order**. Reading the list first anchors you to its framing,
which is the thing under test.

1. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **914 lines**, the spec to be repaired.
   Read it properly and in full; it is the primary source for everything below.
2. `research/pde_ledger_v3/DEFECT_REGISTER.md`, the `C16` entry at `:517`.
3. `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the closed record for the previous step.
4. **Write down, before step 5:** what *you* would repair in this spec before two engines are built from
   it, and how you would specify each repair so that two independently-built engines cannot diverge.
   ⭐ Keep that list — you will be asked which items the decision list handles and which it misses.
5. **Only now** read the decision list.

## What to check

### 1. ⭐⭐ IS EACH ITEM BUILDABLE AS WRITTEN? — the most important question

A repair item that reads well and leaves a choice open produces **two engines that differ while both
builders followed the words.** That has already happened in this spec: `Q3:284` asks for a solution set
*"retaining multiplicity"*, and the two engines invented different constructions.

⭐ For **each** of the five items, ask: if two builders who never speak implement this, can they produce
different objects? Name the specific freedom if so.

⚠ Item 1 is the largest. It supplies a kinetic density, a distinguished axis, a stiffness functional, a
`D` sweep, two premises, and a dimensional treatment. ⭐ **Is anything still unspecified that a builder
must choose?** Check it against how `§7` specifies the existing seven packages and against `§3:106-117`.

### 2. ⭐⭐ IS ITEM 1d's REASONING CORRECT?

Item 1d declines to declare `s_ρ` dimensionless, arguing that `§7:788-791`'s justification for
`XCOEF_BSCALE`'s `s` does not transfer because `ρ_br` also appears **unscaled** in the same package's
kinetic density, whereas `B_comp` appears **only** in the scaled term of `W_XCOEF_BSCALE`.

⭐ **Check that distinction against `§7:740-746` and `Q6:406-440`.** Is it real? Does the choice change
what `Q6` computes, and does either choice make a `Q6` output unable to fail? ⛔ This is the one place the
list makes a judgement call between two defensible options — say plainly if it chose wrong, and why.

### 3. ⭐⭐ IS ITEM 1e's PROPAGATION LIST COMPLETE?

The list names four places that assume a single isotropic kinetic term with one coefficient: `Q6:408`,
`Q6:416`, `Q1:241-242`, `§7:801-806`, plus `§9:878-885`.

⭐ **Search the whole 914 lines yourself for anything else that assumes the kinetic form is fixed,
isotropic, or package-independent** — `Q2`'s two routes, `Q4`'s matrix construction, `Q5`'s scaling, `Q9`'s
census, `Q10`, `Q11`, `§8`'s tag grammar, `§5`'s corollaries. ⛔ A missed site means the new package is
built but some question computes the wrong thing for it, or emits nothing.

### 4. ⭐ IS ITEM 2 THE RIGHT FIX, AND IS 2c SAFE?

`DEFECT_REGISTER#C16` allows two fixes; the list takes "add the stacked source" and rejects "declare `Q8`
silent on transversality", giving a reason.
⭐ Is the reason sound? ⭐ Does `2a` specify the new minor family tightly enough that both engines build the
same one — in particular the treatment of the returned rank `σ`, and the interaction with `Q8a:544-547`'s
rank-`0` convention?
⚠ **`2c` forbids in-engine deduplication.** ⭐ Does that leave `STRATUM_ORDERING` comparable across
engines, given `Q8b:576`? Or does it multiply strata in a way that makes the `Q3`/`Q4` recomputation at
`Q8b:563` blow up or become ambiguous? ⛔ Name the concrete failure if so.

### 5. ⭐ DOES ITEM 3's REPLACEMENT ACTUALLY WORK?

`Q6r:476-477` currently requires a deleted directory. The list repoints it at `scripts/S10_exports.py`.
⭐ **Open that file and check**: does every coefficient that will plausibly be in `COEFFICIENT_ORDERING`
have a resolvable `dimension_key`, and does the replacement state enough for a builder to implement the
resolution? ⭐ What happens for a coefficient the `LEDGER` does not carry — is that specified?
⚠ Only one engine imports the `LEDGER`. ⭐ Is keeping `Q6r` engine-local still coherent, and does anything
else in `§8` need to change with it?

### 6. ⛔ DOES ANY ITEM LEAK AN ANSWER A BUILDER COULD CONVERGE ON?

The builder applying this list reads it in full. ⭐ Does any item state, imply, or hint at what a package,
root, rank, locus or dimension will come out to be? ⚠ **A prohibition leaks as surely as an assertion** —
check the phrasing of the "what this list does not do" section too.
⚠ `DEFECT_REGISTER#C16` contains **measured values** from a witness. ⭐ Check that the decision list does
not carry them, and say so if the repair it commissions would put them into the spec.

### 7. ⛔ DOES ANY ITEM COMMISSION DAMAGE?

⚠ A repair that touches a **correct** part of a spec breeds new defects in the material it changes.
⭐ Is any item broader than the defect it names? ⭐ Is any item's target **not actually defective** — i.e.
would the spec be correct if that item were simply dropped? ⛔ Name what would be damaged.

### 8. ⭐ WHAT IS MISSING?

⭐ Compare your step-4 list against the five items. What would you repair that this list does not?
⚠ **The spec being closed after four rounds and eight legs is not evidence it is correct.**

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, or a way this list would cause a
**wrong spec, and therefore two wrong engines**, to be built. ⛔ Not style, not formatting, not "a builder
might misread this" absent a concrete reading that produces a wrong build.

## Method

- ⭐⭐ **Quote both sides for every finding**: the decision list's text, and the spec/register/source text it
  fails against. A finding without both quotations is not usable.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
  ⛔ A path that resolves is not a source — open it.
- ⭐ Where a question is settled by computation — item 2's minors, item 4's multiplicity constructions,
  item 5's key resolution — **write a script, run it, and paste its literal stdout.** ⛔ A prose derivation
  is not accepted for a contested claim.
- ⛔ Do **not** edit the decision list, the spec, or anything else. Read-only.
- ⛔ Do **not** apply the repair or write any engine code.
- ⭐ End with: **which items on your step-4 list this list handles, and which it misses.**
- ⭐ State explicitly **what you checked that could have failed and did not** — that is part of the result,
  not filler.
