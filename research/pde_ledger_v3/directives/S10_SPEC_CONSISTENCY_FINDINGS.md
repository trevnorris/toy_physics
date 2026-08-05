# S10_SHARED_PHYSICS.md — consistency-pass findings, 2026-08-05

⭐ Produced by an independent pass whose question was **"can two readers who never speak implement this
the same way?"** — ⛔ not "does the document match the engines". Every item below is grounded in the text
or in the two engines' committed output.

⚠ **This file is the input to the spec REWRITE.** ⛔ It is not the rewrite, and ⛔ the rewrite must not be
another appended amendment — see the mechanism below.

---

## ⭐⭐ THE MECHANISM, and it is the one this ledger keeps repeating

⛔⛔ **Every amendment was APPENDED BELOW the instruction list it was correcting, and the instruction list
is what a builder implements from.**

**Measured, §Q7:** line 373 says emit `S_curl[∂u]` at `D=3`. Lines 379–380, six lines later, say compare
**the package's OWN stiffness density**, re-entering at that package's action. ⇒ **one engine followed
line 373 and the other followed line 379.**

⛔⛔ **So the dead Q7 control was a SPEC defect, ⛔ not an engine defect.** The engine that produced a
`Q7_DIFFERENCE` of `0` in all six packages was implementing the sentence it was given.
⚠ **And the appended amendment had itself identified this exact failure** — it says in as many words that
an engine re-typing `S_curl` per package produces a residual that cannot move — and still failed to
prevent it, because it did not delete the sentence that causes it.

⇒ ⭐ **REWRITE, ⛔ never annotate.** ⇒ [[feedback-quarantine-gap-governing-prose]],
[[feedback-fold-dont-write-a-new-doc]], [[feedback-spec-authoring-discipline]].

⭐ **The same residue, not yet bitten:** §Q5 line 285 still carries *"**If** the ratio is a pure power of
`lambdaScale`, also emit that exponent"* — a conditional emission — while lines 288–296 say emit the
whole set **unconditionally** and explicitly record that an earlier version contradicted corollary 4 and
an engine followed it. ⚠ Both engines happened to obey the amendment this time.

---

## ⛔⛔ STATED RESULTS — the spec's own line 8 says it supplies none

⚠ Ten were found. ⭐ Three matter:

1. ⛔⛔ **Line 385: *"For non-curl packages the residual is expected to be nonzero; that is the control
   working."*** ⇒ **§Q7's expected outcome is stated in the specification both engines read.**
   ⭐ **This must be disclosed in the step record: Q7's verdict was supplied, not derived.** The residual
   is still a computed polynomial, ⛔ but the judgement of what it should be was given.
2. ⛔ **Line 208: *"`M_A − M_B` is structurally zero for every stiffness, including a wrong one."***
   ⚠ **Both engines refute it**: `Q2_MATRIX_ENTRY_RATIO = -2` and a nonzero `Q2_MATRIX_RESIDUAL` in
   **13 of 13** packages. ⇒ a stated result that is **wrong**, and it contradicts lines 211–212 of the
   same section, which tell the builder what to do when the two matrices differ by a scalar.
3. ⛔ **Line 404–406** asserts that some §7 action has an allowed real wavevector where the null space
   dimension differs from the generic answer — ⭐ the Q8b answer, stated before Q8b runs.

⚠ Also: line 231 (a locus's real solution set is empty), line 334 (`0 of 552` homogeneity booleans move),
line 341 (`[u]` cancels identically), line 354 (the sign of the Q6d difference), line 440 (root form),
line 447 (`XCOEF_SCALE` "tests arithmetic only", retracted two lines later), line 78 (an overall positive
constant "changes nothing requested in §6" — ⚠ §6 requests `Q2_MATRIX_ENTRY_RATIO`, which it moves).

---

## ⛔ THE UNWRITTEN EQUATIONS — the failure mode that has now recurred a fifth time

⛔⛔ **The dimension of the target is never written.** §Q6 says the coefficient dimensions are *"obtained
by requiring `L` to be an energy density on a `D`-dimensional brane"* — and **nothing else**.
⚠ `[energy]`, the volume element, `[∂_i]` and `[∂_t]` are never stated. Only `[u]` is supplied.
⇒ ⭐ **Both engines had to assume four equations, and every Q6 number rests on them.** They happened to
assume the same four.

⭐ Write them:
```
[energy] = (2, −2, 1)          [u]     = (1, 0, 0)
[L]      = [energy] · length^(−D) = (2 − D, −2, 1)
[∂_i]    = (−1, 0, 0)          [∂_t]   = (0, −1, 0)
```

⛔ **No package's `L` is ever written in full.** §7 gives replacement *fragments*, and `XFORM_SIGNFLIP`'s
row is prose about a sign, requiring the reader to reconstruct the action mentally. ⚠ §7 line 419 insists
**"EVERY control re-enters at the ACTION"** and then never writes the six actions.
⇒ ⭐ **Write all six. It is six lines, and it is the object everything else is computed from.**

⛔ **Route A has no equation while Route B does.** §Q2 line 200: *"substitute the ansatz, strip the common
trigonometric factor, and extract the coefficient matrix"* — "strip" and "extract" are unspecified and
leave the overall scale free. ⇒ ⭐ **that free scale is where the `−2` comes from.**

⛔ **The Q6d counts are nouns, not definitions.** One engine counted per coefficient *expression*
including a declared-dimensionless factor (**9**); the other counted per dimensionful *symbol* (**6**).
⚠ Both are correct readings, so ⛔ **no reviewer of either engine could have caught it.** It hides today
only because the difference is `0` on both sides.

---

## ⛔ §8 PINNED 4% OF THE NAME

⭐ **Measured: the 13 `(package, D)` prefixes match exactly — both engines obeyed §8 perfectly.**
Everything to the right of `D<n>_` was left free, and **562 of 2983 / 4227 tag names exist in both**
(18.8% / 13.3%). At the `<QUANTITY>` token level: 276 distinct in one engine, 376 in the other, **33
shared**.

⭐ **The mechanism is visible in the contrast:** `N1`–`N4` agreed because the spec gave them literal short
labels; `N5`–`N7`, described in longer prose, did not.

⛔⛔ **The single largest cause is unpinned placement of dimension data** — one engine suffixes it onto
the object's own tag (`<OBJ>_Q6_DIMENSIONS`, ~2260 tags), the other builds parallel families
(`Q6_<OBJ>_DIMENSIONS`, ~1025 tags). ⇒ **~3300 tags, none comparable.**

⛔⛔ **And one collision is worse than any gap: `Q3_DETERMINANT` holds the FACTORED form in one engine and
the EXPANDED polynomial in the other.** ⇒ a name-matching checker reports a **false physics
disagreement**. ⭐ Retire the name; use `Q3_DETERMINANT_RAW` and `Q3_DETERMINANT_FACTORED`.

⚠ **Not one `Q5` tag is comparable. Not one `Q8b` stratum tag is comparable.** §8's grammar also cannot
express three things the document requires — a root **pair**, a **stratum**, and the `LOCAL` infix — and
a single real tag violates all three at once.

⚠ **§8 line 489's words and its example disagree**: *"the infix `_LOCAL_` immediately after the engine
prefix"* would be `WL_LOCAL_S10_…`; the example is `WL_S10_LOCAL_…`. ⭐ Both engines followed the example.
⛔ The harness's package matcher was written against the words.

⭐ **The pass produced a pasteable `<QUANTITY>` registry, a dimension-suffix rule, a payload symbol
lexicon, and payload shapes for structured tags.** ⚠ It is long; it is reproduced in the rewrite
directive rather than here.

---

## ⭐ THE THREE TO FIX FIRST — the pass's own ranking, and I agree with it

1. ⭐⭐ **Rewrite §Q7 — ⛔ do not append a third amendment.** It is the only place where the two engines
   computed **different objects and both exited clean**.
2. ⭐⭐ **Add §8's `<QUANTITY>` registry and the dimension-suffix rule.** Mechanical, verifiable by
   grepping the registry against emitted names **before** the next build, and nothing else in the file
   makes the engines comparable.
3. ⭐ **Define the Q6d unknown-coefficient count as an equation over the action's symbol set** — the only
   *numerical* disagreement, and it lands on the denominator of the determinacy verdict that Q6d exists
   to stop resting on a declared value.
