# S11c-b #89 (PY) — make the retained-grade PROJECTION fast (physics byte-identical); nothing else

## 0 · Role and single deliverable

The working-tree SymPy engine
`research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` is the #89-repaired engine (jet-tower
/ variable-coefficient repair already applied and correct). Its primaries complete and emit, but the full
run does not finish because the retained-grade projection is ~5 minutes per call and the controls invoke it
on dozens of objects. **Scope of this task: make that projection fast, with byte-identical output. Change
nothing else** — in particular, ⛔ do **not** touch operator/kernel construction (`build_operator`,
`build_kernel`, `restrict_expression`, the material pullback); their cost is a separate, deferred concern.
The goal is that the **full run (all primaries + all control tasks + the `S11c_b_exports.py` write)
completes in bounded, tractable wall time and no longer hangs** — not a specific minute target. Re-run to
completion and regenerate `.out` and `S11c_b_exports.py`. Modify the engine in place.

⛔ This is a **performance** repair of one function family. The physics is done and must not move. If any
emitted object changes, the fix is wrong.

## 1 · The measured bottleneck (SUPPLIED — cite by NAME; line numbers drift with the repair)

- The projector is `first_shape_series` (locate it by `def first_shape_series`; it is **not** at the old
  L713 — that span is now the background-jet `derivative_map` loop, i.e. **physics**, do not edit it) and
  its object-level wrapper `retained_grade` (`def retained_grade` → `map_object(value, first_shape_series)`).
- `first_shape_series` substitutes `PROFILE_GRADE_SUBS`, then runs a full `sp.series(expr, eta_bg, 0, 2)`,
  then `expand`, then a per-term grade filter (keep `eta_bg≤1 ∧ sigma_W≤1`). It also **swallows** an
  `sp.series` failure (try/except) and keeps the un-expanded expression before the filter — this fallback
  semantics is part of the output and must be preserved.
- Measured cost (this session, `~/.s11_build/S11c_b_89_profile.out`): `retained_grade(operator)` **278–382 s**
  per call; `retained_grade(kernel)` **154–186 s**. A recompute is *not* faster (382 s vs 278 s), so caching
  helps only exact repeats. The controls call the projection on many **distinct** objects, so per-call speed
  is what matters.

## 2 · Hard invariant — the projection's output is byte-identical

`first_shape_series(x)` must return the **same** expression for every input `x` the run ever projects.
Consequently every emitted object is unchanged: the corrected live `S11CB_ENERGY_BASIS_COUNT` stays whatever
the repaired construction computes (⛔ do not change/tune/special-case it — its value is withheld from you and
is not a target), and the frozen-limit (Hessian-freeze) regression still yields the public **26**. The only
target value in this build is the public frozen **26**; a faster projection that returns a different
expression is a failure.

## 3 · The levers, in priority order (each must preserve output exactly)

The two **safe** levers (A, B) keep `first_shape_series`'s computation exactly as-is and are expected to
suffice; the **risky** lever (C) is optional and only if A+B do not make the run complete.

**A · Scalar memoization (safe, do this).** `first_shape_series` takes a single SymPy expression, which is
hashable; it is a pure function. Memoize it at the **scalar** level (keyed on the expression itself). This is
output-identical by construction and deduplicates the sub-expressions that recur within and across objects.
⛔ Do **not** memoize `retained_grade` directly — it receives mutable `dict` operators (unhashable). ⛔ Do
**not** add an object-level projection cache keyed on a partial builder identity such as
`(branch, representative, route)`: the controls deliberately vary `ablation_source`, `ablation_direction`,
`corrupt_material_constraint`, `background_depth`, `full_zero_source`, `include_term_origins`, and aliasing
those would make a corrupted/ablated control compare a baseline **with itself** (residual 0, a silent false
pass). If any object-level reuse is introduced at all, it must key on the **complete** builder tuple and
reuse only exact full-key matches.

**B · Parallelize the projection across cores (safe, primary speedup — the machine has 16 logical / 8
physical cores).** `retained_grade` maps `first_shape_series` (a pure function) over many **independent**
scalars, so distribute those scalar projections over a `multiprocessing` pool. This keeps `sp.series`
**unchanged** — it is the *same* computation, only scheduled across cores — so the result is byte-identical
**provided** you collect the leaf scalars from the object, project them in the pool, and **reassemble them
into the original structure in the original order**. Requirements:
  - Determinism/ordering: reassemble strictly by position; ⛔ never let out-of-order completion change the
    emitted structure. Two runs must produce byte-identical `.out`.
  - Prefer the **fork** start method (Linux default) so workers inherit the already-imported module, the
    built operators, and the caches with no re-`import`/re-init; only the scalar and its result cross the
    process boundary. Combine with lever A (memoize inside each worker and/or the parent) so tiny/recurring
    scalars are not shipped needlessly.
  - Choose the pool size from the core count; guard against nested pools (the projection must not spawn a
    pool while already inside a worker). Keep a serial fallback path (pool size 1) that is exercised by §4.
  - ⛔ Do not parallelize anything else (operator/kernel construction stays untouched, single-process).

**C · Faster per-call projection — OPTIONAL, only if A+B do not make the run complete, and only with proven
equivalence.** Replace the heavyweight `sp.series(·, eta_bg, 0, 2)` with a computation **provably equal** to
it on the retained grade. The `eta_bg`-dependence after `PROFILE_GRADE_SUBS` comes from
`W_bg → W_0(1+eta_bg·w1_profile)` (and `mu_R_bg` likewise) at `/W_bg` and, after differentiation, at **higher
inverse powers** `/W_bg^n`, plus the reciprocal density representatives (`rho4_rhobr` etc.). ⚠ `e_W_bg` is
**not** in `PROFILE_GRADE_SUBS` — do not assume it is. A first-order expansion of each `(1+eta_bg·profile)^(−n)`
factor equals `sp.series(...,eta_bg,0,2).removeO()` on the retained grade, so it is valid **only if** you
handle **every** denominator form the projection actually sees (enumerate all `/W_bg^n` across all branches,
representatives, routes, and depths), **fall back to the original `sp.series`** (and its failure-swallow
semantics) for any form not positively recognized, and keep the original `first_shape_series` reachable as
the §4 reference. If you cannot make C provably equivalent, ship A+B alone — a bounded, completing run beats
a fast-but-wrong one.

## 4 · Required self-checks (emit computed objects; print residuals, do not assert before emitting)

- **Projection equivalence — load-bearing.** Keep the current serial `first_shape_series` as an untouched
  reference and emit `reference(x) − new(x)` (new = the parallel/memoized path, and lever C if used) for a
  set of inputs that **covers every distinct denominator/scalar form**, not one sample: at minimum, projected
  scalars drawn from the SLAB operator, the coupling kernel, and the admissibility operand, across **both
  branches, both density representatives, and both routes**, plus at least one input that drives the
  `sp.series` failure path (so the fallback semantics are exercised). Confirm every residual is `0`.
- **Determinism / ordering — load-bearing for the parallel path.** Emit a whole-object residual
  `retained_grade_serial(obj) − retained_grade_parallel(obj)` for at least one operator and one kernel
  (structural difference, must be `0`), proving the parallel reassembly preserves order and structure. The
  regenerated `.out` for a fixed input must be reproducible byte-for-byte across two runs; if you cannot
  assert that cheaply, at least confirm the parallel/serial object residual is `0`.
- **Frozen-limit regression.** Emit the Hessian-freeze-switch `S11CB_ENERGY_BASIS_COUNT` and its residual
  against `26`.
- **Timing summary.** Emit wall-seconds per primary task and per control task and the total, so the remaining
  cost is visible (kernel-build cost may still be large — that is expected and out of scope; the run must
  simply complete).

## 5 · Method and script obligations

1. Print computed objects, not conclusions; print residuals, do not `assert` before emitting.
2. No hard-coded results; no tautological residuals; controls remain able to fail; ⛔ do not reduce any
   control's cases/anchorings/representatives.
3. Run the FULL engine to completion; keep both anchorings and both density representatives; regenerate
   `.out` and `S11c_b_exports.py`. Verify both are regenerated and non-empty.
4. Run `reduction/derived_or_declared.py` and `reduction/engine_output_checks.py`; report their literal
   verdicts.

## 6 · Builder report

- The unified diff of every function changed, with a one-line rationale per change that it is
  output-preserving; confirm operator/kernel construction was **not** touched.
- The literal projection-equivalence residuals (must all be `0`, including the failure-path input) and the
  frozen-limit residual against `26`.
- The per-task timing summary and total wall time (before/after if available).
- Confirmation the full run completed and regenerated `.out` + `S11c_b_exports.py`, with byte sizes.

## 7 · Supplied vs computed

SUPPLIED: the measured bottleneck and its per-call cost (§1), the hard invariant (§2), the public frozen
`26`. COMPUTED and subject to external review: the optimized projection and every residual. A passing
frozen-limit regression confirms only the frozen limit; the §4 projection-equivalence residuals are what
certify the physics is unchanged.
