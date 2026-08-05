# Repair `S10_SHARED_PHYSICS.md` — ⭐ SIX targeted edits, ⛔ NOT a rewrite

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
⛔ **Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file.**
⚠ The as-built text is pinned at tag `s10-as-built`; ⭐ the step record cites that, ⛔ not your output.

⭐⭐ **SCOPE, and it is deliberately narrow (user, 2026-08-05): keep S10 to S10.**
⛔ **A full rewrite into a schema-checked contract with a tag registry is DEFERRED** —
⇒ `S10_spec_rewrite_directive_v2.md` + `S10_contract_DEFERRED_findings.md`. ⭐ That work buys **future
cross-engine comparability**, ⛔ not S10's correctness. ⛔ **Do not do any of it here.**

⛔⛔ **EXPLICITLY OUT OF SCOPE — ⛔ do not touch, even if you see it is wrong:**
the `§8` tag-name grammar and any quantity registry · the payload serialiser · the sign-diagnostic
classifier · the `Q6d` unknown-coefficient count · `report_contract`-style machinery.
⚠ **All are known-defective and all are recorded elsewhere.** ⭐ If you spot another, ⭐ **report it, ⛔ do
not fix it.**

---

## ⛔⛔ E1 — `§Q7`: rewrite the section. ⭐ It is the reason this directive exists.

⚠ **The section currently instructs TWO DIFFERENT OBJECTS.** One sentence names a fixed curl density; a
later amendment says to use the package's **own** stiffness density. ⇒ ⛔ **the two engines each followed a
different sentence and both exited clean.** ⭐ **Delete, ⛔ do not annotate.**

**⭐ What Q7 IS — state this in the section:**
> The stiffness must be written in general `D`, because the 3-vector cross product exists only at `D = 3`.
> `§Q7` checks that the general-`D` antisymmetric-derivative form **becomes** the ordinary curl-squared at
> `D = 3`. ⭐ It is a **form-and-normalisation** check on the coefficient that carries the wave speed.
> ⛔ It is **not** evidence for the mode count, which is computed in every `D` and does not depend on it.

**⭐ The comparison, as equations:**
```
g_ij   := INDEPENDENT symbols standing for ∂_i u_j        ⛔ never a k×a amplitude curl
S_pkg  := THIS package's stiffness density, obtained from THIS package's action
          with the g_ij substituted in                    ⛔ never a re-typed curl density
c_i    := Σ_{j,k} ε_ijk g_jk                              COMPUTED by the CAS
emit:    S_pkg  ·  c·c  ·  (S_pkg − c·c)                  all three
```

⛔⛔ **`c·c` MUST BE REACHED BY COMPUTATION from the Levi-Civita definition. ⛔ Do NOT write out the
expanded polynomial.** ⭐ That computed side is the **only** reason this check can fail on physics — the
density side is typed, so a hand-expanded reference would compare a typed object against a typed object.

⭐ **`S_pkg` must come from the ACTION, so that mutating the action alone — with the package selector held
FIXED — moves it.** ⚠ **Measured:** an implementation keyed on the **selector**, taking no action object at
all, passes a *"change the form in one package and watch it move"* test. ⇒ ⛔ **re-deriving from the
selector is not re-entering at the action**, and the words must make the two distinguishable.

⛔ **DELETE from `§Q7`:** the sentence naming a fixed curl density as the compared object · the sentence
stating what the residual is **expected** to be and for which packages · the sentence describing the
**measured** failure shape under a corrupted action · the "curl vector of the **amplitude**" wording, which
contradicts the derivative formula on the next line.

⚠ ⭐ **Say plainly what this check does NOT do:** it compares a **typed** density against a **computed**
reference. ⇒ it catches a wrong **form** and a wrong **normalisation**; ⛔ it does **not** establish that
the density was assembled from the action rather than re-typed — ⭐ that is what the action-mutation
requirement above is for.

## ⛔⛔ E2 — add the DISTINCTNESS premises. ⭐ Two controls are dead without them.

⚠ **Measured:** the premise set omits the distinctness conditions for the **anisotropy scale** and the
**coefficient scale**. ⛔ With the degenerate value admitted, the anisotropic package's action is
**identical to the main package's**, and the coefficient-scale package collapses the same way.
⇒ ⛔ **two of the six controls are dead.**
⭐ **Add both distinctness premises to the supplied set.**
⚠ ⭐ **And say why they are legitimate**, because the existing instruction *"do not add a premise to force a
solver to decide"* is being read as licence to omit them: ⭐ **a premise that keeps a CONTROL DISTINCT is
not a premise that forces a decision.** ⛔ Do not conflate them.

## ⛔ E3 — `§Q2`: delete the stated claim about the two matrix routes

It asserts a result about the difference between the two routes; ⚠ **both engines refute it**, and it
contradicts the two lines below it in its own section.
⛔ **Delete it. ⛔ Do not replace it with any other claim about that difference.** ⭐ The section says what
to COMPUTE and stops.

## ⛔ E4 — `§Q6`: write the FOUR unwritten equations

The coefficient dimensions are said to follow from *"requiring `L` to be an energy density on a
`D`-dimensional brane"* — ⛔ **and nothing else is written.** ⇒ both engines had to **assume** the dimension
of energy, the volume element, and the two derivative operators. ⚠ They happened to assume the same four.
⭐ **Write all four as equations.** ⛔ Nothing else in `§Q6` changes.

## ⛔ E5 — `§7`: write all six actions IN FULL

⚠ `§7` insists **every control re-enters at the ACTION** and then gives replacement **fragments**, with one
package's row being **prose about a sign**. ⇒ a reader must reconstruct the action mentally.
⭐ **Write each package's Lagrangian in full.** ⚠ It is six lines, and it is the object everything else is
computed from.

## ⛔ E6 — sweep the file for STATED RESULTS and delete them

⭐ Read every sentence against the file's own opening promise that it supplies no results.
⭐ **A statement of what to COMPUTE stays.** ⛔ A statement of what something **equals**, what is
**expected**, what "the control working" looks like, or what was **measured**, **goes**.
⛔ **Delete them; ⛔ do NOT relocate them.** ⚠ Several sit in sections that are otherwise out of scope —
⭐ **deleting a stated result is always in scope.**

---

## ⭐ ACCEPTANCE — ⛔ run these and paste literal output

1. ⭐ **The two-reader test on `§Q7`:** read your rewritten section as an implementer and list **every**
   object it could tell you to compare. ⛔ If the list has more than one entry, ⛔ **the rewrite has
   failed.**
2. ⭐ **Grep the whole file** for sentences beginning "an earlier version", and for the words `expected`,
   `measured`, `it turns out`. ⭐ Paste every hit with its line number and say why each survivor is not a
   stated result.
3. ⭐ **Confirm the six package actions** are each written as a complete Lagrangian, and that the
   distinctness premises appear in the supplied premise set.
4. ⭐ **Old and new line counts**, and a list of every section you touched. ⛔ Any section outside `E1`–`E6`
   appearing in that list is a scope violation — ⭐ report it.

## Report back — ⛔ under 25 lines

1. One line per `E1`–`E6`: done / partially / not.
2. The acceptance output.
3. ⭐ Anything you saw that is wrong and that you did **not** fix because it was out of scope. ⭐ **This is
   wanted** — it is how the deferred list stays honest.
4. ⛔ Do not report what any engine's values came out to be.
