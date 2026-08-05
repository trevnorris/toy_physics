# S10 SymPy engine — repair round 1

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py`
**And:** its sibling `…_sympy_audit.premises`
**Its specification, ⚠ AMENDED since you last saw it:**
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md` — ⭐ **re-read §1, §3,
Q5 and Q6 in full**; four things changed.

⛔ **Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file** except the two named above.
⛔⛔ **Do not read the sibling Mathematica engine** (`mathematica/S10_brane_mode_spectrum_*.wl`), and ⛔ do
not read `research/pde_ledger_v3/steps/` or `research/pde_ledger_v3/paper/`.

⚠ After every fix, re-run and confirm it completes under **10 minutes** and exits `0`.

---

## ⭐ SPEC AMENDMENTS — apply first; ⚠ **two of these were MY errors, ⛔ not yours**

**S1 · §1 — `[u] = (1,0,0)` is now a STATED premise.** ⭐ **Emit it as a premise tag** and add it to
`.premises`. ⚠ It is currently hardcoded in the walker and **appears nowhere in the output**.

**S2 · §3 — the period average is over the PHASE**, `⟨F⟩ ≡ (1/2π)∫₀^{2π} F dφ`. ⭐ **Your code already does
this correctly.** ⛔ **But `PREMISE_PERIOD_AVERAGE` emits the interval `(0, 2π/√ω²)`** — the ω-dependent
`t`-limits the spec now forbids — and `.premises` line 2 repeats it as a "one-period **time** average".
⇒ ⭐ **Fix the premise tag and the `.premises` line to state what the code actually does.** ⚠ Those two are
the **only** places `sqrt(omegaSquared)` appears in the whole tag stream.

**S3 · Q5 — emit the ENTIRE Q5 tag set for every root, unconditionally.** ⚠ The old text contradicted
corollary 4 and you followed the old text. ⭐ Where the ratio is undefined the payload is an **explicit
marker object**; ⛔ the tag never disappears.

**S4 · Q6d — new.** Emit the number of independent dimension equations, the number of unknown coefficient
dimensions, their difference, and whether the system was over-/exactly-/under-determined. ⭐ Emit a tag
recording that the action's homogeneity booleans are **vacuous** when that difference is `≤ 0`, and ⭐ give
homogeneity tags for **solved action terms** versus **everything else** distinguishable names.
⚠ **Measured on your engine:** perturbing `[u]` moved **528** payloads and **0 of 552** homogeneity
booleans. ⛔ Do not try to make those booleans able to fail — ⭐ **label them.**

---

## ⛔⛔ R1 — BLOCKING: the "real locus solve" is a coordinate-subspace enumeration, not a solve

`RealLocus.branches` (≈439–453) is built by trying **every subset of `k` set to zero**. The payloads of
`Q3_ROOT_COINCIDENCE_LOCI` and `Q8_RANK_DROP_LOCUS` (≈700, 875–877, 889, 899) come from that alone.

⚠ **Measured:** on a locus of the form `k_i² − k_j²` it returns **only the origin branch**, while `solve`
returns the true sheets — ⛔ and those true sheets are emitted **only under a `_LOCAL_` tag**, which §8
defines as **excluded from cross-engine parity.**

⇒ ⛔⛔ **Any stratum requiring a RELATION AMONG COMPONENTS is absent** from the locus, from the
allowed-region test, and from Q8b's re-run list. ⚠ **That is exactly the stratum class Q8b exists for.**

⭐ **Also:** `sp.solve(..., domain=sp.S.Reals)` (≈474) — ⛔ **`solve` has no `domain` argument and silently
ignores it.** The real restriction is coming only from `real=True` symbols.

⭐ **Fix:** solve the locus **properly over the reals** — `solveset` / `nonlinsolve` with an explicit real
domain, or `solve` followed by an explicit reality filter you emit. ⭐ The resulting locus must be a
**parity-visible** tag, ⛔ not `_LOCAL_`. ⭐ Emit the enumeration too if you like, but ⛔ it is not the
answer.

## ⛔⛔ R2 — BLOCKING: two guards exit non-zero on a physics outcome

≈670 (`spectrum solve produced no roots`) and ≈843 (`dimension solve returned no solution`).
⚠ Q6 **explicitly anticipates** an under-determined or inconsistent dimension system **as a finding**, and
§5 requires a physics finding to **exit 0**. ⛔ Both currently abort every remaining package.

⭐ **Fix:** emit the condition as its own tag with its operands and **continue the sweep.** ⛔ Exit non-zero
only on an exception or a genuine operational failure.

## ⛔⛔ R3 — BLOCKING: an unhandled `ComplexInfinity` kills the run

`Q2_MATRIX_ENTRY_RATIO` is `M_A[0,0]/M_B[0,0]` (≈805). When the denominator vanishes the payload is `zoo`,
and `DimensionWalker.dimension` (≈167) has no rule for it ⇒ `RuntimeError`, exit 1, **every later package
lost**. ⚠ `nan` survives because `is_Number` is true; `zoo` is not.

⭐ A vanishing matrix entry is a **physics configuration**, ⛔ not an operational failure. **Fix:** handle
`zoo`/`ComplexInfinity`/`nan` explicitly in the walker, emit an explicit marker payload, and continue.

## ⛔ R4 — six stratum tags carry a hand-typed empty payload

`emit_stratum_q3_q4` (≈732–738) emits the Q3 coincidence family as `sp.Tuple()` instead of recomputing it.
⇒ ⭐ They **cannot move under any ablation.** ⭐ **Fix:** recompute them at the stratum, exactly as the
generic path does.

## ⛔ R5 — matrix payloads span multiple lines, breaking the tag grammar

≈376 continuation lines. §8 requires **one line per tag**; a line-oriented consumer truncates them.
⭐ **Fix:** render every payload on a single line. ⚠ This is the same class as the previous version's
duplicate tag name — a defect that breaks automated consumption entirely.

## ⛔ R6 — the three-way sign test degrades to an unevaluated `sign(...)`

4 of 30 `Q3_SIGN` tags (all `XFORM_ANISO`). §3's **joint** assumption set exists to prevent exactly this.
⭐ **Fix:** apply the joint predicate via `refine`/`ask` at the sign test. ⚠ If it still cannot decide,
emit an explicit **undecided** marker — ⛔ not a bare unevaluated expression, which the sibling engine
answering definitely would turn into a false physics disagreement.

## ⛔ R7 — route A silently discards terms

`.coeff(sp.cos(theta))` (≈789) drops any Euler–Lagrange term **not** proportional to `cos θ`, with no
residual emitted. ⭐ **Fix:** emit the discarded remainder as its own tag. ⛔ A silent drop is
indistinguishable from there being nothing to drop.

## ⛔ R8 — `dimension()` on an `Add` returns the FIRST term's dimension

≈156–160. ⇒ inhomogeneity inside a denominator or a `Pow` base is **invisible**; the boolean only inspects
top-level expanded terms. ⭐ **Fix:** walk into denominators and `Pow` bases, or ⭐ emit explicitly that
those positions were **not** inspected. ⛔ Do not leave it looking like full coverage.

## ⛔ R9 — a stratum whose branch fixes every component is silently skipped

≈924, with no tag; and `allowed` is a per-locus aggregate while **all** its branches are enrolled
(≈899–900). ⭐ **Fix:** emit a tag for the skip with its reason, and make the allowed-test per-branch.

## ⛔ R10 — `N7`'s two routines are applied to the same pre-processed matrix

Both run on `rank_input` (≈534, 542) rather than on `M_r`. ⚠ That **narrows the independence** the check is
built on — a defect introduced by the shared pre-processing would move both operands together.
⭐ **Fix:** compute at least one of the two from `M_r` directly.

## ⚠ R11 — runtime headroom is thin

Baseline ≈110 s, but a single stiffness-form change pushed `XFORM_ANISO` past memory at ≈542 s.
⭐ Reduce peak cost where you can **without dropping any `(package, D)` pair**; ⛔ if a pair genuinely
cannot run, emit it in `SKIPPED_PAIRS` with a reason — ⛔ never drop it silently.

---

## ⭐ What NOT to change — confirmed live by ablation

⛔ `N3` (rank of `M_r.col_join(kᵀ)`), the separate construction of `M_A`/`M_B` (one-sided corruption moved
each alone), per-package re-entry at the action, Q5's `k→λk` scaling, the dimension **tree walk**, Q7's
independent gradient symbols, and ⭐ **the registry allowlist handling** — patching locus validation off
before `load_registry` and restoring it after is **correct**, because `registry_read` otherwise opens every
`source_locus` path including ones into `steps/`. ⛔ Do not "simplify" that away.

## Report back — ⛔ under 25 lines

1. One line per `S1`–`S4` and `R1`–`R11`: fixed / partially / not, with line numbers.
2. New tag count, wall-clock, exit code.
3. ⛔ **Do not report what any value came out to be**, and ⛔ do not say whether anything "worked".
4. ⭐ Anything in this list you believe is **wrong**, or any fix that would break something the
   "do not change" list protects. ⭐ This is wanted.
