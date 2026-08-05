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
⛔⛔ **And delete the sentence in §Q7 that states what the residual is EXPECTED to be, and for which
packages.** ⚠ It tells a builder the answer, in a document whose own opening line says it supplies none.
⛔ **Delete the neighbouring sentence describing the MEASURED failure shape under a corrupted action** —
⚠ it names the corruption that was tried and what it did, which steers an implementer just as effectively.

#### ⭐⭐⭐ THE SETTLED §Q7 EQUATIONS — ⛔ DECIDED, ⛔ NOT YOURS TO CHOOSE

⚠⚠ **These were settled between the orchestrator and a peer consult on 2026-08-05, because two review legs
correctly said a BUILDER must not be handed this choice.** ⭐ **Write them; ⛔ do not re-open them.**
⚠ If you believe one is **wrong**, ⭐ **report it and stop** — ⛔ do not substitute your own.

**⭐ What §Q7 is FOR — state this in the section, because leaving it implicit is how it drifted:**
> The model's stiffness must be written in general `D`, since the 3-vector cross product exists only at
> `D = 3`. §Q7 is the **bridge**: it checks that the general-`D` antisymmetric-derivative form **becomes**
> the ordinary curl-squared at `D = 3`. ⭐ It is a **form-and-normalisation** check on the coefficient that
> carries the wave speed. ⛔ It is **not** evidence for the mode count, which is computed in every `D` and
> does not depend on it.

**⭐ The objects, as equations:**
```
L₀     := L |_{∂ₜu → 0,  ∂ᵢuⱼ → g_ij}      the action's static part, g_ij INDEPENDENT symbols
w      := the coefficient the ACTION CONSTRUCTOR actually placed in front of the density
S_ext  := L₀ / w                            valid only once w ≠ 0 is PROVED
c_i    := Σ_{j,k} ε_ijk g_jk                COMPUTED from Levi-Civita, ⛔ never typed out
```

⛔⛔ **`w` IS READ FROM THE ACTION. ⛔ It is NEVER a retyped `sign · coefficient / 2`.** ⚠ **This is the
subtle one and it was caught by computation, not by reading:** a package may **scale** the coefficient, and
a retyped weight leaves that factor undivided ⇒ ⛔ **a COEFFICIENT control then emits a spurious FORM
residual**, which looks exactly like a finding. ⇒ ⭐ [[feedback-per-tooth-ablation]].

⛔⛔ **PROVE `w ≠ 0` BEFORE ANY SIMPLIFICATION, and emit the status.** ⚠ **Measured:** a CAS cancels
`w·S / w → S` symbolically, so a degenerate weight yields a **finite, wrong, silent** answer ⛔ — it does
**not** raise. ⇒ If nonzeroness cannot be proved, ⭐ **emit an explicit undefined status; ⛔ do not form a
residual.**

**⭐⭐ EMIT AS TWO SEPARATE CHECKS — ⛔ do not compress them into one residual:**

| check | emit | what it catches |
|---|---|---|
| **1 · extraction** | `L₀`, `w`, the `w ≠ 0` status, `S_ext`, and the reconstruction certificate `L₀ − w·S_ext` | coefficient/sign degeneracy and provenance failure |
| **2 · reference identity** | `S_ext`, `c·c`, and `S_ext − c·c` | the form-and-normalisation bridge itself |

⭐ **Why the split is required:** a sign or coefficient control moves **`w`**, ⛔ not the form residual.
⇒ Compressed into one number, those controls either vanish or masquerade as form failures. ⚠ ⛔ **An
overall sign therefore CANNOT be policed by the density residual** — if it is to be policed, it needs its
own weight diagnostic, and the section must say so.

⭐ **`c·c` is a SUPPLIED DIAGNOSTIC REFERENCE — say so in the text.** ⚠ Otherwise §4's "only the action and
the ansatz may be combined by hand" reads as forbidding the very object §Q7 compares against. ⭐ The
Levi-Civita **definition** is supplied; ⛔ the expanded polynomial is **computed**.

⭐ **Outside `D = 3`:** emit `Q7_APPLICABLE` **and the same quantity tags at every `D`**, with an explicit
not-applicable payload. ⛔ Do **not** swap in a different conditional tag name — a tag's presence may depend
only on package and quantity, ⛔ never on a payload's value.

⭐ **Residual orientation is `package_object − reference`**, matching clause 2's `A − B`.

⛔ **Delete "curl vector of the amplitude."** ⚠ It contradicts the derivative formula on the next line, and
the independent-symbol requirement is what stops the check being one expression substituted into itself.

#### ⛔⛔ STILL OPEN, and it is YOURS — the tag-name contract contradicts itself ACROSS FILES

§Q7 **orders** the Q7 tags renamed "for what they now are", while the live repair directive **forbids**
changing tag names because they are wired into a comparison config — ⚠ and that config currently maps
**both** spellings. ⭐ **Resolve it here, in the same breath as item 2's registry**, and say which name is
canonical.
⛔⛔ **A rename the harness config does not follow SILENTLY DROPS the comparison** — ⚠ and a dropped
comparison is indistinguishable from agreement. ⇒ ⭐ **Either choose a name the config already maps, or
state the config change as a REQUIRED companion edit** and list the exact rows. ⛔ You may not edit
`reduction/` yourself; ⭐ naming what must change there is part of this deliverable.

#### ⭐⭐ THE SIGN DIAGNOSTIC — settled, ⛔ write it as a FOUR-WAY classification

⚠ One marker currently covers two different situations, so a reader cannot tell which occurred.
⭐ Under the **effective premises** `P_eff` (joint premises + solver conditions + any active stratum
conditions), classify each sign call by solving `P_eff ∧ sign_condition_σ(root)` over the reals for
`σ ∈ {positive, zero, negative}`:

| verdict | condition | emit |
|---|---|---|
| `DETERMINED(σ)` | exactly one region satisfiable, the other two **proved empty** | the region |
| `PREMISES_INSUFFICIENT` | **two or more** regions satisfiable | ⭐ the **witness models**, explicitly |
| `CAS_UNDECIDED` | the solver establishes neither | the unsettled sub-expression |
| `INCONSISTENT_PREMISES` | all three empty | the premise set |

⚠ **Two opposite-sign witnesses are sufficient but ⛔ NOT necessary** — a zero/positive or zero/negative
pair equally proves the three-way sign is not determined.
⛔⛔ **And note the consequence:** an **unevaluated** sign of a quantity the premises already declare
positive is ⛔ **a CAS/proof-route failure, ⛔ NOT premise insufficiency.** ⚠ Today those are reported the
same way. ⭐ Say which tags carry the diagnostic and whether it covers **every** sign call or only some.

### ⛔⛔ 2. §8 — pin the whole tag name

⚠ **Measured: both engines matched all 13 `(package, D)` prefixes exactly, and then diverged on
everything to the right of `D<n>_`.** §8 pinned about 4% of the name.

⛔⛔ **THERE IS NO EXISTING REGISTRY. ⛔ DO NOT GO LOOKING FOR ONE — ⚠ an earlier version of this directive
said the findings file carried one, and the findings file says it was reproduced here. ⭐ Both claims are
FALSE; neither file contains it, and the pass that produced it left no artifact.** ⚠ Corrected 2026-08-05
after two independent legs hit the circular reference.

⇒ ⭐ **You are BUILDING the registry, from the two committed `.out` files.** ⛔ Not adapting one, ⛔ not
starting from a draft. ⭐ Derive it from the tag names the engines actually emit, and say what you had to
invent.

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
⛔⛔ **And delete §Q2's stated claim about the difference between the two matrix routes** — ⚠ it is a
stated result, and it contradicts the two lines below it in its own section. ⛔ Do **not** replace it with
any other claim about that difference; ⭐ the section says what to COMPUTE and stops.

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

## ⛔⛔ WHAT I HAVE TOLD YOU MUST NOT ENTER THE REWRITTEN FILE

⭐⭐ **You are a RELAY, and that is the hazard.** This directive hands you measured facts about what the two
engines currently produce — counts, which engine did what, where they diverged — ⭐ **so that you can decide
what to write.** ⛔ **None of it may appear in `S10_SHARED_PHYSICS.md`.**

⛔⛔ **AND THE SAME APPLIES, HARDER, TO `S10_SPEC_CONSISTENCY_FINDINGS.md`, WHICH THIS DIRECTIVE REQUIRES
YOU TO READ.** ⚠ **Both review legs flagged it independently: it is UNSANITISED and it is copy-paste bait.**
It carries a matrix entry ratio, a "`N` of `N` packages" count, a homogeneity-boolean count, an
unknown-coefficient count for each engine, tag totals and intersection percentages, and the sentence
stating what a residual is **expected** to be.
⇒ ⭐ **Read it for WHERE THE DEFECTS ARE. ⛔ Copy NONE of its numbers.**

⭐⭐ **AND THE OLD SPEC IS ITSELF A LEAK SOURCE — ⛔ you are deleting those sentences, ⛔ not relocating
them.** ⚠ Legs found result-bearing text throughout: a structurally-zero matrix claim contradicted by both
engines, an expected residual, an answer-shaped nullity hint, measured payload-movement counts, a
"returns a definite answer" promise, and a prior step's tag-parity gap count. ⛔ **None of it survives the
rewrite in any form.**

⚠ **Why this is the whole risk of the task:** that file is the ONE artifact both engines read. ⛔ A measured
result copied into it is not a leak into *your* work — it is a leak into **every future build**, and
cross-engine agreement on it becomes **vacuous**, because both engines are then repeating a sentence rather
than computing.

⭐ **The test to apply to every sentence you write:** *could a builder read this and learn what an answer
comes out to, or which way a comparison goes, without computing it?* ⛔ If yes, it does not go in — **even
inside a sentence forbidding it.** ⚠ **Measured three times in one session: a prohibition leaks exactly as
well as an assertion**, and one leak was introduced **by the repair for another.**

⇒ ⭐ **`S10_SHARED_PHYSICS.md` states what to COMPUTE. It never states what anything EQUALS, what is
EXPECTED, what was MEASURED, or what "the control working" looks like.**

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
5. ⭐⭐ **THE DISAMBIGUATION TEST — ⛔ this is the one that decides whether §Q7 is actually fixed.**
   ⚠ **Strengthened 2026-08-05: both legs found the earlier version SELF-GRADED and passable by writing
   text narrow enough to describe one engine.** ⇒ it is now a **per-decision** table, ⛔ not a verdict.

   ⭐ Build a row for **each** settled decision above — the static-part extraction, the weight read from
   the constructor, the `w ≠ 0` proof, the reconstruction certificate, the reference identity, the
   computed Levi-Civita reference, residual orientation, and the non-`D=3` payload. For **each row × each
   engine**, record `satisfies / violates / cannot tell`, ⭐ **quoting the sentence of YOURS that decides
   it and the engine line that settles it.**

   ⛔⛔ **A single "cannot tell" anywhere in the table means the rewrite has FAILED — keep going.**
   ⚠ ⛔ **Do not change either engine.** ⭐ You are testing your text against them, ⛔ not them against it.
   ⭐ **Your text ruling BOTH engines out is a legitimate and expected outcome** — the decisions above were
   fixed **independently of** what either engine does. ⛔ A row that "violates" is a finding to report,
   ⛔ never a reason to soften the sentence that produced it.
6. ⭐ **Tag-name coherence:** grep whatever Q7 name you make canonical against
   `reduction/checks_S10.yaml`, and report every spelling that file currently maps. ⚠ A rename the config
   does not follow **silently drops the comparison**, which is indistinguishable from agreement.

## Report back — ⛔ under 40 lines, plus the three sweeps

1. One line per item 1–7: done / partially / not.
2. The three sweep outputs (these may exceed the line budget).
3. Old and new line counts.
4. ⭐ Where you think the findings file is **wrong**, and anything you changed that I did not ask for.
5. ⛔ Do not summarise the physics.
