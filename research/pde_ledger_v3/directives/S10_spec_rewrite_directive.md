# Rewrite `S10_SHARED_PHYSICS.md` — ⛔ a REWRITE, ⛔ not a seventh amendment

**File to rewrite:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
**Your input, read it in full first:** `directives/S10_SPEC_CONSISTENCY_FINDINGS.md`

⛔ **Do not commit.** ⛔ **Do not edit any other file.** ⛔ Do not touch either engine, either `.out`, or
anything under `steps/`, `paper/`, or `reduction/`.

---

## ⭐ WHY YOU AND NOT ME

⚠ **I wrote this document, and four of its six amendments were my own defects — three of them introduced
by revisions that were fixing something else.** ⇒ ⭐ **The author has to change.** Read the findings file
as evidence, ⛔ not as a patch list to apply mechanically, and ⭐ tell me where you think it is wrong.

## ⭐ WHAT THIS DOCUMENT IS FOR

It is the **single written specification** from which two engines are built **blind to each other** — one
Mathematica, one SymPy. It is therefore the one artifact both share, so ⛔ **a defect in it defeats the
two-engine design by construction: both engines implement the same mistake and then agree.**

⛔⛔ **AND IT IS NOT A HISTORICAL RECORD.** The engines that exist today were built from the current text
and their output is committed; git holds the old version. ⇒ ⭐ **You are writing the document the NEXT
build reads.** ⛔ Do not preserve a sentence because an engine once implemented it.

---

## ⛔⛔ THE ONE RULE THAT MATTERS: DELETE, DON'T APPEND

⚠ **Every amendment so far was appended BELOW the instruction list it corrected — and the instruction
list is what a builder implements from.** ⭐ Measured consequence, in the findings file: §Q7 told a reader
to emit one object at line 373 and a different object at line 379, and **the two engines each followed a
different line.**

⛔⛔ **When you fix something, the sentence that caused it must be GONE.** ⭐ A reader must not be able to
find the old instruction anywhere in the file. ⚠ If a correction is worth recording, it goes in **one
clearly-marked section at the end** — ⛔ never inline next to the instruction it contradicts.

---

## ⭐ THE WORK

### ⛔⛔ 1. §Q7 — rewrite the section outright

The findings file documents the divergence. ⭐ **Write the three compared objects as EQUATIONS over the
package's own action**, so that "which stiffness" cannot be read two ways.
⛔ **Delete the object list that names `S_curl` directly.**
⛔⛔ **And delete the sentence stating what the residual is expected to be for non-curl packages.** ⚠ That
sentence tells a builder the answer, in a document whose own opening line says it supplies none.

### ⛔⛔ 2. §8 — pin the whole tag name

⚠ **Measured: both engines matched all 13 `(package, D)` prefixes exactly, and then diverged on
everything to the right of `D<n>_`.** §8 pinned about 4% of the name.

⭐ **The findings file carries a proposed `<QUANTITY>` registry, a dimension-suffix rule, a payload symbol
lexicon, and payload shapes for structured tags. ⭐ Use it as a starting point — ⛔ not as gospel.**
⭐ **Check it against both `.out` files yourself** and fix what it gets wrong or omits.

⛔ **Requirements the registry must meet:**
- ⭐ Every quantity either has **one** registry name or is **engine-local**. ⛔ No third option.
- ⛔⛔ **Retire `Q3_DETERMINANT`.** One engine used it for the factored form and the other for the
  expanded polynomial. ⚠ That is **worse than a missing pair** — a name-matching checker reports a
  **false physics disagreement.**
- ⭐ The grammar must be able to express a root **pair**, a **stratum**, and the engine-local infix.
  ⚠ Today one real tag violates all three at once.
- ⚠ §8's words say the local infix goes *"immediately after the engine prefix"* and its example shows it
  **after the step token**. ⭐ Pick one, ⛔ and delete the other.
- ⭐ **Pin the payload serialiser.** Two incompatible choices are currently both permitted.

### ⛔ 3. §Q6 — write the unwritten equations

⚠ The coefficient dimensions are said to follow from *"requiring `L` to be an energy density on a
`D`-dimensional brane"*, and **nothing else is written**. ⇒ both engines had to assume `[energy]`, the
volume element, `[∂_i]` and `[∂_t]`. ⭐ **Write all of them as equations.**

⭐ **And define the Q6d unknown-coefficient count as an equation over the action's symbol set.** ⚠ One
engine counted per coefficient *expression* including a declared-dimensionless factor; the other per
dimensionful *symbol*. ⛔ Both are correct readings, which is why no reviewer of either engine could have
caught it. ⭐ Whichever you choose, §7's "dimensionless by declaration" must become a **computational
input**, not a remark.

### ⛔ 4. §7 — write all six actions

⚠ §7 insists **every control re-enters at the ACTION** and then gives replacement *fragments*, with one
package's row being prose about a sign. ⭐ **Write each package's Lagrangian in full.** It is six lines,
and it is the object everything else is computed from.

### ⛔ 5. §Q2 — give Route A an equation

⚠ *"Strip the common trigonometric factor and extract the coefficient matrix"* leaves the overall scale
free, while Route B is written as an equation. ⭐ **Fix Route A's normalisation in an equation.**
⛔⛔ **And delete the claim that `M_A − M_B` is structurally zero for every stiffness** — ⚠ both engines
refute it in **13 of 13** packages.

### ⛔ 6. Sweep the whole file for stated results

⚠ The findings file lists ten. ⭐ **Find them yourself as well** — read every sentence against the
opening line's promise that the document supplies no result — and ⛔ **delete every one.**
⭐ A statement of **what to compute** stays. ⭐ A statement of **what it comes out to**, of what is
"expected", of what "the control working" looks like, or of what was "measured", goes.

### ⛔ 7. Resolve the live contradictions

⚠ The findings file lists them. ⭐ Two that must not survive:
- a conditional emission that corollary 4 forbids, still present six lines above the amendment that
  forbids it;
- a section requiring **every** package to emit the **full** tag set, against another requiring one
  question at a single `D` only. ⚠ The engines resolved this differently and the file gave neither a way
  to know which was right.

---

## ⛔ WHAT NOT TO CHANGE

⛔ The physics. ⛔ The action, the ansatz, the question list's *content*, the three clauses, the five
corollaries, the structural rule that every control re-enters at the action.
⭐ **You are fixing whether the document SAYS what it means, ⛔ not what it means.**
⚠ If you believe a physics statement is wrong, ⭐ **report it, ⛔ do not act on it.**

## ⭐ ACCEPTANCE — ⛔ run these, do not assert them

1. ⭐ **Grep your registry against both `.out` files** and report, per engine: how many emitted quantity
   names are in the registry, how many are not, and the full list of those that are not.
   ⚠ **Expect a large mismatch — the engines predate the registry.** ⛔ That is data, ⛔ not a failure.
2. ⭐ **A contradiction sweep**: for every instruction that says *emit X*, confirm no other sentence in
   the file says *emit something else* under the same name. ⭐ Report the method you used.
3. ⭐ **A stated-result sweep**: list every sentence you deleted under item 6, with its old line number.
4. ⛔ **Confirm the file contains no sentence beginning "an earlier version"** anywhere except the single
   end section. ⚠ Those are the append-residue markers.

## Report back — ⛔ under 40 lines, plus the three sweeps

1. One line per item 1–7: done / partially / not.
2. The three sweep outputs (these may exceed the line budget).
3. Old and new line counts.
4. ⭐ Where you think the findings file is **wrong**, and anything you changed that I did not ask for.
5. ⛔ Do not summarise the physics.
