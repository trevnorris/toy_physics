# DIRECTIVE — repair two unearned claims in the S11 Mathematica audit

Edit **only** `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`.
⛔ Do not commit. ⛔ Do not create any other file.

---

## ⛔⛔ THE ONE RULE THAT MATTERS MOST HERE

An independent SymPy script covering this same physics **exists in this repository**. ⛔⛔ **DO NOT READ
IT.** ⛔ Do not read `research/pde_ledger_v3/scripts/`, `research/pde_ledger_v3/reduction/`,
`research/pde_ledger_v3/V3_STEP_PLAN.md`, or anything under `research/pde_ledger_v3/steps/`.

⭐ **Why this matters more than usual:** the two scripts were written independently, and their value is
that they **could disagree**. If you repair this one by copying what the other concluded, they agree
because they were *made* to — and the check is destroyed retroactively. **You are fixing a script that
overclaims. You are NOT making it match anything.**

⇒ If a repair leaves you unable to conclude something, ⭐ **the correct output is to report what you
computed and state plainly what it does and does not establish.** ⛔ Do not reach for a conclusion to
fill the gap.

---

## THE TWO DEFECTS

Both are the same class: **the script certifies more than it computes.**

### D1 · `WL_S11_TRANSVERSE_MATCHING` emits a fixed conclusion token unconditionally

The tag prints a result token that is a **literal**, chosen when the code was written and emitted
regardless of what the surrounding computation returns. Nothing computed can change it.

Additionally, the residual it does compute is **circular**: the transverse kernel is *selected* by
testing a condition, and the "check" then re-asserts that same condition on the vectors that passed it.
It cannot fail by construction.

**Repair — pick whichever is honest, and ⛔ do not manufacture certainty:**
- **either** derive the conclusion from a quantity the script actually computes, by a test that could
  return the opposite answer;
- **or** stop emitting a conclusion, print the computed quantities, and add an explicit statement of the
  **scope** of the calculation — what class of question this algebra can and cannot settle, given what
  it was and was not given as input.

⚠ Be precise about what the calculation was handed. State any physical ingredient that would be required
to reach a stronger conclusion and that is **absent** from the setup.

### D2 · `dimensionalClosure` cannot fail

The four dimension checks are algebraic rearrangements of the definitions they check: a quantity defined
as `a - b` is "verified" by adding `b` back. They test that subtraction is invertible; they say nothing
about whether the dimensional formulas are right.

**Demonstrated:** setting `displacementDimension` to a physically absurd value makes the script print
nonsense dimension vectors and **still emit `WL_S11_VERDICT: PASS`, exit 0**. The whole block is
unguarded.

**Repair:** make the dimensional results **able to fail** — a wrong premise or a wrong derived exponent
must produce a `FAIL` and a non-zero exit.
⛔ **Do not add a table of expected dimension values.** That would be a hardcoded target and is worse
than the current state. ⭐ Guard it against something the script can establish **internally** — e.g.
consistency requirements the physics imposes among the terms of the Lagrangian, which a wrong premise
would violate.

⭐ **Prove your repair works: run an ablation yourself**, perturbing a premise, and show the script now
FAILs. Report the ablated value and the resulting output. ⛔ Restore the file afterwards.

---

## ⛔ ACCEPTANCE — this is how I will check you

**Every other `WL_S11_*` line must be UNCHANGED.** This repair touches **claims**, not computations. The
invariant counts, the dynamical matrix, all four spectra with their nullities and orientations, the
cross-sector results, the degeneracy condition, the four dimension tags, `KW_SQUARED`,
`BOUND_CONDITION`, and both controls must print exactly what they print now.

⛔⛔ **If you find yourself wanting to change a computed physics value, STOP and report it instead of
changing it.** That would mean the review missed something, and it must be surfaced, not silently
absorbed.

Script conventions are unchanged: standalone, **no arguments**, imports **no files**, print-only.
Mathematica has a **2-seat** licence — never more than 2 concurrent `math -script`. ⚠ A licence-contention
run exits **40** with no output; that is not a script failure.

## Report

Under 25 lines: what you changed for D1 and for D2, the ablation you ran for D2 and its output, and
confirmation that every other tag is unchanged (show a diff of the tag lines).
⛔ Do not summarise this directive back to me.
