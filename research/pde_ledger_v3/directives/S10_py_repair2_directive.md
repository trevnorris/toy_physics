# S10 SymPy engine — repair round 2 (final)

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py`
**Spec, ⚠ AMENDED:** `directives/S10_SHARED_PHYSICS.md` — ⭐ **re-read §5 corollary 5, Q6d, Q7, §7**.

⛔ Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file except the `.premises` sibling.
⛔ Do not read `steps/`, `paper/`, or the sibling Mathematica engine.
⚠ Re-run after every fix: under **10 minutes**, exit **0**.

⭐ **A review leg ablated this engine and confirmed the computed chain is live** — N3, N7, both matrix
routes, per-package action re-entry, Q7's two independent sides, and the registry allowlist (measured: zero
opens under `steps/`, `paper/`, `mathematica/`). ⛔ **Do not change any of those.**

## ⛔⛔ P1 — BLOCKING: stratum objects are dimensioned with the WAVEVECTOR ERASED

`emit_stratum_q3_q4` substitutes integer literals for `k`, then reports dimensions with the **outer**
walker, where a number returns `ZERO_DIM`. ⇒ every `Q8_STRATUM#_…_Q6_DIMENSIONS`/`_HOMOGENEITY` is **short
by the wavevector**. ⚠ **Measured: 2 of 4 stratum-root dimension tags differ from their generic
counterparts.**

⛔⛔ **This is the exact defect Q6 names, and it is the one that survived two review legs and a full
ablation suite at the previous step.** ⭐ **Fix:** dimension the stratum objects with a walker that still
carries the wavevector's dimension, or ⭐ dimension the **pre-substitution** expressions. ⛔ Do not report a
dimension computed on a numerically substituted object as though it were the object's dimension.

## ⛔⛔ P2 — BLOCKING: unrecognised function nodes return dimensionless, silently

`if expr.is_Function: return ZERO_DIM` catches everything the explicit list misses — `Piecewise`, `RootOf`,
`Min`/`Max`, `re`/`im` — all of which `sp.solve` can return, and it makes the `raise` below it
**unreachable for function nodes**. ⇒ a dimensionful node reads as dimensionless with a clean exit.
⭐ **Fix:** unhandled heads must produce an **explicit indeterminate marker** and emit the head they failed
on, exactly as the walker already does elsewhere. ⛔ Never fall through to zero.

## ⛔ P3 — homogeneity does not inspect every `Add`

`dimension` returns the **first** term's dimension for an `Add`; the compensating helper only visits `Add`s
inside a `Pow` **base**. ⇒ an inhomogeneous `Add` appearing as a `Mul` factor or a function argument is
neither dimensioned correctly nor flagged, and the boolean can read **true** for an inhomogeneous
expression. ⭐ **Fix:** walk every `Add` wherever it occurs, or ⭐ emit explicitly which positions were not
inspected. ⛔ Do not leave it looking like full coverage.

## ⛔ P4 — Q6d's unknown count is a constant (corollary 5)

`unknowns` is a fixed 6-tuple regardless of the package's action. ⚠ **Measured:** an ablation removing a
coefficient from every action moved the equation count and the difference **13/13** and moved this
**0/13** ⇒ the determinacy verdict rests on a **declared denominator**. ⭐ **Fix per the amended Q6d:**
count the unknowns **from the package's own action**.

## ⛔ P5 — the SOLVED / UNSOLVED homogeneity split is drawn at one tag

Only the expanded Lagrangian is labelled `SOLVED_ACTION`; the EL system, both matrices, the averaged
Lagrangian and the determinant are **linear consequences of the same solved coefficients** and carry the
"unsolved" name — and the vacuity tag by name qualifies only the action class. ⭐ **Fix:** label every
object whose dimension is determined by the solved coefficients as `SOLVED`, and make the vacuity tag
qualify that whole class.

## ⛔ P6 — an undecidable branch is reported as a definite `false`

In `locus_allowed_data`, a branch that fixes every component but whose point fails to materialise emits
`false` rather than the `undecided` marker used in the adjacent arm. ⇒ a materialisation failure reads as
"this stratum is disallowed" and **Q8b never re-runs there**. ⭐ **Fix:** emit the undecided marker.

## ⛔ P7 — a non-reproducible payload ordering

`nested_add_dimensions` iterates `expr.atoms(sp.Pow)`, a **set**, so ordering follows per-process hash
randomisation. ⚠ **Measured: one payload differed between two runs of identical code.** ⭐ **Fix:** sort
deterministically. ⚠ A byte-level cross-engine comparison would otherwise report a disagreement that is
not one.

## ⛔ P8 — the run record is declared, not observed

`RUN_PAIRS` is emitted from the config **before** any package runs and `SKIPPED_PAIRS` is a hardcoded empty
tuple. ⭐ **Fix per the amended §7:** accumulate after each package finishes emitting; emit after the sweep.

## ⛔ P9 — one guard precedes its emit

`deduplicate_roots` raises before the raw solution tag is printed. ⭐ **Fix:** emit the raw object first,
then guard.

## ⛔ P10 — determinant method deviates from the directive

`det(method="domain-ge")` where the directive mandates `berkowitz`. ⭐ Either conform, or ⭐ emit a tag
recording the method actually used and why. ⛔ Do not leave it silently different from what was specified.

## ⚠ P11 — two 5-second `SIGALRM` budgets make payloads machine-dependent

Both emit a route tag, so it is visible — ⚠ but the same file on a slower host emits different payloads for
a non-physics reason, and the sibling engine runs on this same host. ⭐ Raise them to a value that will not
trip under normal load, or ⭐ make the fallback deterministic rather than time-dependent.

## Report — ⛔ under 18 lines
1. One line per `P1`–`P11`: fixed / partially / not, with line numbers.
2. For **P1 and P4**, state the perturbation you ran and that the tag moved.
3. New tag count, wall-clock, exit code; and confirm two consecutive runs produce **byte-identical** output.
4. ⛔ Do not report what any value came out to be. ⭐ Anything here you believe is wrong is wanted.
