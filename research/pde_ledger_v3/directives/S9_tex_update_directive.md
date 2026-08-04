# S9 TeX card — update it to match the rebuilt step

⛔ **Edit in place. Do not restructure the card and do not rewrite sections that are still correct.**

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S9_light_requires_shear.tex`

**Source of truth — read it first, in full:**
`/var/projects/toy_physics/research/pde_ledger_v3/steps/S9_light_requires_shear.md`

⛔ Do not commit. ⛔ Do not modify any other file.

⚠ **Check `paper/macros.tex` before you move anything into a field.** Some fields are **suppressed in the
default build**, so reader-critical content placed in one becomes invisible in the PDF. ⚠ The existing
card carries a comment at its Verification paragraph explaining that it deliberately uses
`\paragraph{Verification.}` rather than `\stagefield{Verification}` for exactly this reason — ⭐ **respect
that and keep the content visible.**

---

## Why this update exists

The step's two engines were rebuilt. ⭐ **The physics did not change** — every value the previous engine
computed is identical today. ⛔ **But several statements in this card are now FALSE**, and two of them
assert the absence of things that now exist.

## S1 — ⛔⛔ the closing paragraph is now FALSE and must be replaced

The card ends by stating that **no separate SymPy audit was written, deliberately**, and that a second
engine *"would be a third execution of the same two lines."*

⛔ **A SymPy engine now exists**, and the reason matters: **it found a defect the Mathematica engine had**
— a wrongly computed dimension that had survived two review legs and a full ablation suite.

⭐ Replace that paragraph with the honest account from the step record: the exemption rested on a
**conditional** rule (*a second engine earns its place where the algebra could genuinely disagree*), the
condition **expired** as the engine grew, and the conclusion was inherited without re-checking it.
⭐ State what the second engine bought.

## S2 — ⛔ the "scope of that audit" paragraph is now FALSE

It says the audit is **insensitive to a wrong overall prefactor, a flipped potential sign, and the assumed
brane dimension**, each leaving its verdict green.

⛔ **All three are now addressed** — there is a sign control, there are form controls, the dimensions are
solved **symbolically in `D`**, and there is no verdict tag at all any more (a script may not state a
conclusion). ⭐ Replace this with the limits that are **actually current**, listed in the step record's
*"what this step still does not establish"* section. ⛔ Do not soften them and ⛔ do not drop any.

## S3 — Verification: replace with what the rebuild established

⭐ **The single most important thing to convey, and the card should make it plain:** two verification
mechanisms were used, and **neither could have caught what the other did.**

- a **wrong dimension** was caught by the **second engine**, and was invisible to review and ablation
  because the first engine computed it *consistently* and produced a plausible-looking number;
- a **wrong homogeneity test** was caught by **reviewers deriving from scratch**, and was invisible to
  cross-engine comparison because it came from the **shared directive**, so both engines agreed on it.

Also record: the mode count is a **rank** computation, ⛔ not the `M·T = 0` test (which returns a false
negative under anisotropic inertia); non-dispersiveness is tested by **scaling** `k → λk`; the dynamical
matrix has **two independent routes** whose independence was proven by one-sided corruption; nine actions
are exercised; and the automated consumer reports **12 cross-engine quantities agreeing, 0 disagreeing**,
with all dimensionful expressions homogeneous.

## S4 — Dimensions: they are now solved symbolically in `D`

The card derives them at `D = 3`. ⭐ Give the closed form in `D` and the result that
`[μ_R] − [ρ_br]` is **independent of `D`**, then specialise to `D = 3` for the registry comparison.

⭐ **Say what that buys**, because it retires a recorded weakness by explaining it: the old audit's
insensitivity to the assumed brane dimension was filed as a blind spot; it is an **identity** — the
speed's dimension *cannot* see `D`.

## S5 — ⭐⭐ a physics claim needs sharpening, and this is the substantive change

The card's framing implies the curl-only form is what gives light its transverse waves.

⛔ **A control shows otherwise.** Ordinary gradient elasticity — `−½ μ_R Σ(∂_i u_j)²`, an ordinary elastic
solid — carries **two transverse modes at the same `c² = μ_R/ρ_br`**. It simply *also* propagates the
longitudinal.

⇒ ⭐ **What the curl-only form uniquely buys is the ABSENCE of a propagating longitudinal mode** —
Maxwell's third demand — ⛔ **not** the presence of the transverse ones. The `what is new` table already
says *"forced by Maxwell's no-longitudinal demand"*, so this **confirms** that row rather than
contradicting it — ⭐ but it is now **computed**, and the divergence-only control shows the converse: the
roles **swap**, the transverse pair falls to `ω² = 0` and the longitudinal propagates.

⚠ ⇒ **Adjust any prose that reads as though this computation establishes "light requires shear."** It does
not: it establishes a **conditional** statement, given the stiffness form. ⛔ The section title may stay;
the body must not overclaim.

## S6 — keep, unchanged

The status (`\StatusOpen`), the falsifier and its live status, the two-sided requirement table, the
registry additions, the provenance paragraph about moves 2 and 3 resting on **one** external source with
**no executing script**, and the regime and departure paragraphs. ⭐ **P2 is still cited-and-never-computed
and the card already says so at the right strength — ⛔ do not weaken that.**

---

## Report back — under 20 lines

1. The file you changed and whether it still compiles if a build is available (⛔ do not install anything).
2. One line per section S1–S5: what you replaced, and the paragraph it now reads as.
3. Anything in the step record you could not represent in the card, and why.
4. ⭐ Anything in the card you believe is still wrong or overclaimed that this directive did not name.
   ⭐ This is wanted — the card is what a reader actually trusts.
