**CLEAR FOR THIS PACKET-ACTION FOURIER BUILD**

This covers source correctness and faithfulness to method.md §4 for one guarded execution of the 190 fixed requests. It does not clear the full action evaluator, and it establishes no result. I found no blockers.

## What I inspected
- **Read in full:** `build.md`, `numerical-library.py`, `worker.py`, `evidence-guide.md`, and method.md §4 (lines 207–377).
- **Read in part:**
  - `manifest.json` (about 511 KB). I read `runtimeLibrary`, the source pins, the resource block and the route keys.
  - `raw-helper.py` (Journal and decode, lines 218–292).
  - method.md line 64 (`D_j`).
- **Operands:** `saved/preflight/field-79ea…-adapter-input.json` (constant field) and `field-abd7…-adapter-input.json` (degree-4 field). I also read `saved/selected/pressure-addresses.json` (wave-multiplier entries only).
- **Pinned mpmath 1.3.0:** `ctx_mp.py` `clone`, `eigen_symmetric.py` `gauss_quadrature` and the eigenvalue sort, and `matrices.py` `__iter__`.

## Coverage limits
- I did not open most of the 826 saved files, so I did not individually check all 34 vectors, all 544 address joins, or the zero-proof files. I checked the join logic in `worker.py:121-145` and two representative operands.
- I did not examine `launcher.py`, `tooling-tests.py`, `execution-authority.json`, the shared guard, the supervisor or the hook. Their identity checks (`worker.py:50-73`) are source-read only.
- No code was run, so every numerical outcome is unknown. That includes whether the 1e-25, 1e-42 and 1e-12 checks pass.

## Checks that held
**Coefficient and zero-proof joins** (`worker.py:121-145`)
- The saved `coefficientVector` is descending in tanh. The degree-4 field reads 24/11009+…, …, 8/11009+…, and `list(reversed(pairs))` at line 133 gives the ascending order that `poly_eval` expects.
- Each coefficient joins its input `left` and the `cancelled==ZERO` return.
- Rational strings such as `-27/11009` parse through `Fraction` in `scalar` (`numerical-library.py:158`).
- The wave-jet input `left` joins `a['waveMultiplier']` (for example `I*p`, consistent with `(ik)^n`).

**Families and requests**
- `family_plan` gives 13 X families plus 3 Y plus 3 Y′, which is 19.
- 19 × 2 carriers × 5 arguments is 190 (`worker.py:161-168`).
- X carrier is p0 and Y carrier is −p0 (`worker.py:166`).
- For nu=0 the argument and carrier are computed bit-identically, so `sign(0)=0` applies exactly.

**Signs, normalization and ordering** (`numerical-library.py:160-194`)
- `nu = arg − carrier` matches method §4 lines 241-243.
- X has a 1/(2π) normalization, Y has none, and the native `D_j` matches method line 64.
- The Q recurrence `Q' + (i·carrier − w/64)Q` is applied before the coefficient, and the carrier phase stays in `exp(−i·nu·z)`.
- The argument-derivative factor is `−sign·i·x`, so X gives −ix and Y gives +ix. That matches d/dl of ∫e^{ilx}.
- The shift `z = x0 + y − i·c·sign(nu)` matches method line 245.
- The contour moves the integrand, not the coefficient field.

**Constant analytic reference** (`:241-253`)
- It is separately coded and correct. It uses the Gaussian transform `s√(2π)e^{−s²nu²/2−i·nu·x0}`, the factor `(i·arg)^n`, the argument-derivative recurrence `P'+(s²·carrier−i·x0−s²·arg)P`, and the Y sign `(−1)^m`.

**Contour tail**
- |tanh| ≤ 1 holds because 5/10 < π/4. The tail bound is the sum of absolute coefficients with the stated M0/M1/M_n recurrence and the factor `exp(c²/128−c|nu|)` (`:196-206`).
- Route A's finite-contour truncation needs no vertical-segment term. The full shifted line equals the real line by Cauchy, because the strip is analytic and the integrand decays, and the tail bound covers the rest.

**Panels**
- Cuts are at 0 and −x0, widths follow `h_A` and `h_B`, and the endpoint equality checks are exact. The radii are integers and the centers dyadic, so the checks cannot fail spuriously.

**Gauss and Kronrod rules**
- mpmath's `gauss_quadrature` returns ascending nodes, and the matrix iteration works as the code uses it. The cloned contexts are independent.
- The E8 system is correct. P7 is odd and E8 is even, so only j=1,3,5,7 are needed, with four unknowns. The moment and right-hand-side signs and the exact residual check are right.
- Interlacing brackets, bisection to 1e-48 and the 15×15 weight solve are sound. K15 is exact through degree 23, so the residual check through n=23 is the right range. The embedded G7 weights come from the same G7 nodes.

**Adaptive Route B**
- Per-panel acceptance halves the allocated target, so the sum of accepted estimates cannot exceed 1e-13. The final `require` is consistent with that.
- Stagnation raises.

**Comparisons and durable evidence**
- Absolute comparisons are strict `<1e-12`, including the three analytic differences for constants.
- Every input, panel, point, value, weight, adaptive parent and child, and comparison is persisted before the final `require`.
- Failure stops with no fallback. The SQLite file closes and hashes in `finally` (`worker.py:175-177`).
- The mpmath pins sit inside `sourcePins`, which `verify_gate` checks.

## Optional comments (not prerequisites)
1. **Disk and time.** My rough estimate is 5–6k panels per request at |nu|≈27–32, with about 30–40 KB of JSON per panel and one FULL-sync commit per panel. That gives tens of GB of SQLite and hours of fsyncs over the roughly 76 large-|nu| requests. Before launch, confirm free disk. The 4 GiB limit covers memory only, and an out-of-disk error would fail closed but waste the run.
2. **Shared code between routes.** A and B share `product()`, so a shared convention error is invisible to A-versus-B. Only the constant families have an independent reference. I inspected the signs and ordering above and found none, but passing nonconstant requests does not independently test them.
3. **Constant detection by length.** `constant_reference` keys off `len(coefficients)==1` (`:242`). Saved certificates mark constants as degree 0, so this is consistent here, but it silently skips the reference if a field ever carried trailing zeros.
4. **Extra resolution cuts.** The method's ±1/s and ±1/L cuts (lines 328-337) belong to the outer k/l and t/Q arrangements, not the y-integral. I did not treat their absence in the Fourier bank as a deviation.

## Preserved unknowns
- The numerical outcomes (none computed).
- Whether pass or fail at the 190 requests transfers to any other argument. It does not.
- The embedded estimate is empirical and not a rigorous bound. The 1e-14 tail bounds are numerical evaluations of analytic bounds, not interval certificates.
- Full action units, profile and PV, the height route, tails and controls remain open.