# Four corrections to `S10_SHARED_PHYSICS.md` — ⛔ NOT a fifth rewrite of `§Q7`

**File:** `research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
**Edit in place. Do not commit. Do not modify any other file.**

A review leg computed four defects in the repaired spec. **Three of them were transcribed from the repair
directive's own suggested prose**, so the directive's acceptance gates could not catch them —
`S10_spec_Q7_repair_directive.md:27-31` and `:80`. ⛔ **The wording below is the leg's, not the previous
author's.** Apply it; do not improve on it, and do not re-open `§Q7` beyond these four edits.

None of this changes a computed value: both engines were built from the text pinned at `s10-as-built`.
What it changes is what a rebuild would be told to compute.

---

## C1 — `:378-379`. The action-mutation requirement is unsatisfiable as written

The sentence reads *"Mutating **the action alone**, with the package selector held **fixed**, must move
the emitted `S_pkg`."* But `§2` defines `S[∂u]` as the stiffness **density**, and the stiffness sign, the
stiffness coefficient and the inertial coefficients all live **outside** it. Measured through the SymPy
engine's own `build_action`, selector held at `MAIN`: of five action mutations, three move the action and
leave `S_pkg` fixed. The same is visible in the committed outputs — `MAIN`, `XFORM_SIGNFLIP` and
`XCOEF_SCALE` emit **byte-identical** `§Q7` densities on both engines.

⇒ an implementer applying the sentence either declares a correct build non-conforming, or widens the
emitted object until it does move — which changes both `§Q7` operands and the residual.

**Scope the requirement to the object it names:**

> Mutating **the stiffness density in the action**, with the package selector held **fixed**, must move
> the emitted `S_pkg`.

⭐ Keep everything around it: the *"must be the stiffness-density object … not an equal-valued rebuild"*
sentence, and the concession that provenance rests on cross-engine comparison.

## C2 — `:361-362`. `§Q7` states the dispersion relation

The clause *"a form-and-normalisation comparison bearing on the coefficient that carries the wave speed"*
asserts, in the one file both engines read, that `ω² = c²k²` with `c² ∝ μ_R/ρ_br`. That is the answer to
`Q3`'s sign test, to `Q5`'s scaling exponent, and to `Q6`'s implied speed dimension. `§Q5` at `:273-275`
forbids exactly this shape — *"naming a power tells you the answer before you compute it."*

**Delete from `"This is a form-and-normalisation comparison…"` through `"…wave speed"`, keeping
`"This is a form-and-normalisation comparison."`** Leave the following sentence (that `§Q7` is not
evidence for the mode count) untouched.

## C3 — `:89` and `:442-443`. One of the two dimensionless declarations is derivable, and `§7`'s stated reason is false for it

`§7` justifies both declarations with *"the energy-density requirement alone fixes only the sum of a
scale factor's dimension and its coefficient's."* Built from the four supplied dimension equations, the
supplied `[u]` and the `§7` actions, with **no** dimensionless declaration made, the coefficient-dimension
system has **nullity 0** for `XFORM_ANISO` and **nullity 3** for `XCOEF_SCALE`.

⇒ the justification is true for `s` and **false** for `s_ρ`, and the asymmetry runs the **opposite way**
to the distinctness asymmetry `E2` established.

**Take the second option the leg offered, because it keeps both symbols visible and adds a check rather
than removing one:** state the dimensionless premise for `s` only, and require `[s_ρ]` to be emitted as a
**solved** value with the declaration **differenced against it** — both operands and the residual, per
clause 2. Correct `§7`'s stated reason so it describes only the case it is true of.

⚠ **This is the shared-spec cause of a live cross-engine disagreement.** Measured on the pinned build:
`Q6_UNKNOWN_COEFFICIENT_DIMENSION_COUNT` reads `6` on SymPy for **all thirteen** package-and-dimension
pairs, and `6` or `9` on Mathematica according to whether the package carries a scale factor. **Both
conform to the current text**, because it never says whether a symbol declared dimensionless is still an
unknown. ⭐ **Say which it is.**

## C4 — `:477-478, 484, 488-489`. The new stratum token describes one engine and not the other, and `E7`'s mandatory disclosure was omitted

Transcribed to a regex and matched against every emitted tag containing `STRATUM`: Mathematica conforms
on all of its stratum tags; SymPy conforms on **none** of its, every one carrying `Q8_` between `D<n>` and
`STRATUM<s>`, and a further group using a skipped-stratum scope token the grammar has no slot for. The
cross-engine intersection of stratum-scoped name suffixes is **empty**.

`E7` said plainly that what is **not** acceptable is leaving this unstated. The token was added and the
disclosure was not, so the file now presents a grammar under which the committed pair looks alignable when
nothing aligns.

**Add one sentence:** the as-built engines' stratum tag names are not aligned, and their comparison is
manual. ⛔ Do not touch anything else in `§8`.

---

## Out of scope — report, do not fix

Everything else. In particular the leg also recorded, and these are **deliberately not** being fixed here:
`§9`'s supplied-vs-tested table was not updated for the new premises · two instructions were deleted along
with the results attached to them (`§Q7`'s tag-naming instruction; `§8`'s *"follow it exactly"*) ·
`§3` and `§7` disagree about what the package-domain premise set contains · `:76`'s claim about `§6` no
longer holds now that `§6` contains `§Q7` · `:458`'s *"scaling never leaves the family"* states what a
control does · `:451-452`'s replacement discloses less than the clause it replaced · the requirement that
`c_i` be computed from the Levi-Civita definition is met by **neither** as-built engine.

⭐ If you judge any of those to belong with these four, **say so and leave it alone.**

## Acceptance — run these and paste literal stdout

1. For each of `C1`–`C4`: the line numbers touched and the literal before/after text.
2. **The two-reader test on `§Q7`**, again: read the section cold and list every object it could tell you
   to compare. More than one pair is a failure.
3. **The satisfiability test `C1` exists for:** state in writing which mutations your corrected sentence
   requires to move `S_pkg`, and confirm a conforming engine can satisfy it. ⛔ If a correct build would
   still be declared non-conforming, the edit has failed.
4. Grep the file for `wave speed`, `carries the`, and `Levi` and paste every hit with its line number.
5. Old and new line counts, and every section you touched. Anything outside `C1`–`C4` is a scope
   violation — report it.

## Constraints

- No script over 600 seconds. **Do not launch Mathematica or `wolframscript`.**
- No `git add`, no `git commit`, no other git write.

## Report back — under 20 lines

1. `C1`–`C4`: done / partially / not.
2. The acceptance output.
3. Anything you measured that contradicts this directive.
