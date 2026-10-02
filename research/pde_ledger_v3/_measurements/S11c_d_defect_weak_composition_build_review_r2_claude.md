**CLEAR FOR THIS GLOBAL WEAK-COMPOSITION BUILD**

I found no blockers. This is a source-only review: I read the files and ran nothing, so I could not check any hash and did not re-run `tests.py`. Hash and pin identity is left to the fresh gate.

## Blockers

None.

## Non-blocking tooling notes

These are optional hardening, not scientific defects.

- **Fragile identity step.** `worker.py:398-399` checks Re β and Im β with `sp.re(...).expand(complex=True)`.
  - If sympy leaves `re()` or `im()` unevaluated, `J.zero` fails closed with preserved evidence, which would burn the single authorized run.
  - Building the identity from the explicit conjugate product would avoid that.
  - I verified the identities by hand: Re β = 30/((10+δ)²+9) and Im β = (9+δ(10+δ))/((10+δ)²+9).
- **Jet multipliers not checked.** The worker never compares `waveMultiplier` to the stored jet spec, i.e. (−3i)^{n_t}(i/5)^{n_2}(i/10)^{n_3}(ip)^{n_1}. A cheap exact check over all 13,260 addresses would close this. The Leibniz control covers the ordering semantics separately.
- **Constants not derived in code.**
  - `worker.py:406-424` emits 11, 121, 18, 36/√a_*, KJ, KD, 135/8 and 665/8 as typed constants. Only the polynomial-gap inequalities and algebraic joins are machine-checked.
  - I re-derived all of them and they hold. For example, KJ = 2·(|a||μ|²·WL/4 ≤ 2/5)·121·Cq/b³ and KD uses (WL/4)|μ| ≤ 1 with b².
  - The weight facts also hold: sup w = 4/e ≤ 2, ∫w = 10, and per-root bound 4·2+10 = 18.
  - These are labelled as assessed arguments, not machine proofs.
- **Reflected-root control point.** The point k=3/2, l=2, t=1/2 has t = l−k, so qh equals qo and qs equals qi (qh=3/2, qs=2). The control is still applicable and responsive, since qh≠qs and the movement is nonzero. A point with t≠l−k would exercise the route identity more sharply.
- **Launcher timeouts.** The 5 s arming poll and 30 s handshake are startup handshakes, not computation deadlines.

## Review coverage limit

I did not read the five 13,260-address arrays. I read these:
- the 608 representatives;
- the complete 320-cell grade census;
- all 34 field objects;
- the whole-tag and typed-direct files;
- the factor operands and returns (the factor-3 operands, and the `-return.json` files, which are all literal zero);
- the controls;
- the emitting worker (`source/composition-worker.py:779-815`).

The emitter builds every address, including the zero-status ones, with the same full key set. Per its design, `tests.py` (`test_complete_grade_routes_and_corruption`, `test_actual_metadata_contracts_no_restoration`) already walks all 13,260. The future run must still check every actual address.

## What I verified against the data

- **Fields.** The constructors are in the 34 saved field objects (8 constants including 0, 26 tanh polynomials of degree ≤5).
  - Every field reduces to a polynomial in T=tanh(x/10), the recurrence (1−T²)P′/10 is correct, and the L=10 join is present.
  - The worker's no-leftover-x, no-function and finite-constant-denominator guards are strict.
- **Grade rectangle.** 16 triples × 5 rows × 2 faces × 2 slots = 320 cells. With 17 triple-component pairs × 39 jets, that gives 13,260 addresses. Status counts are 11,232 + 1,420 + 608, and zeros are preserved.
- **Joins to the saved densities.**
  - The worker's Dformula, Jformula and Bformula all match the saved constants. Prefactor ratio (WL/4i)=−5i/2 and J constant (5/32)aω² reproduce, and β=ω/(10−iω).
  - H contact is j/4 and the subtraction is A(t)(j(Q−t)−j(Q))/(2it).
  - The direct map keeps qh=q(k+t) and qs=q(l−t) as separate symbols.
- **Global bounds.**
  - b=3000/11101 and a_*²=879/400.
  - |q|≤|p|+4, via (K+4)² − K² − 453/50 = 8K + 347/50.
  - Inequality (2) holds, since |rad| ≥ |κ²−p²| for δ>0.
  - The B and Bn variation constants hold, because B depends on qo only through 1/|qo+β|, so the Q-tail needs no growth bound.
  - The numerator-gap polynomials all have nonnegative coefficients.
  - The compact q≤5 and |Q|≤6 bounds are not reused.
- **Whole-tag counts.**
  - Direct is one tag with the exact `Dwhole_<face>` signature.
  - Mixed iteration is exactly one H plus one Jwhole (10 unique face × component keys).
  - Bare, factored and whole stay distinct, and nothing is integrated twice.
- **Normal controls.** There are exactly four names. Address 9126 is a normal-slot NATIVE_SLOPE with a nonconstant consumer, and the native normals are ±i·qo.
- **New controls.**
  - Predicates select constant source and consumer fields before evaluation. The Leibniz source is deliberately nonconstant.
  - The flat point maps both qi and qo at k=l=2.
  - Evidence is emitted before the guards.
  - Selected addresses exist for both faces, and all three controls are labelled formal.
- **Gate and launcher.**
  - The gate fixes the authority path and pin, exact scope, one authorized run and no deadline, plus the argv tail and the review fields.
  - The launcher is hook-first, the guard command is unchanged, and the worker imports only stdlib before `containment()`.
  - Strict errors, copied inputs, posthashes and 'x'-mode failure saving are present.

## Runtime obligations

- **Gate.** Build a fresh pinned gate after both literal build verdicts and write the review record. The record's worker, manifest, guard, supervisor and launcher hashes must equal the gate's. The method SHA must be 98c295d4…, the authority pin 39cc989f…, and the worker 77ec3fa0….
- **Execution.** Run the hook-first launcher once, under the 4 GiB, 1-CPU, 32-task, no-deadline guard.
- **Runtime checks.**
  - Check all 13,260 addresses and every factor join at runtime.
  - Confirm that all 34 fields reconstruct to residual 0, including the first derivative.
  - Confirm that all four reused controls validate and the six new controls give exact nonzero movements.
  - Confirm that the root-route and β/κ identities are zero.
- **Afterwards.**
  - Inspect every copied operand and posthash.
  - Read stderr.
  - Treat any failure as immutable.
  - Do no automatic retry.
  - Treat a clean exit as not scientific acceptance.

## What success would and would not establish

**Would establish:**
- every saved coefficient is a finite-constant polynomial in tanh with a bounded-derivative class;
- the 13,260-address and 320-grade census matches the saved definitions;
- the global constants and numerator envelopes are algebraically joined to the saved J, direct and H objects;
- the formal controls respond.

Together with the assessed analytic arguments, that supports a continuous bilinear weak action on S×S with polynomial degree ≤3, a height action bounded by Pk^{3/2}, and the whole direct and mixed iteration entering once.

**Would not establish:**
- a machine proof of the measure theory;
- any evaluated integral, loss or finite solve;
- a bounded inverse, scattering solution or plane-wave extension;
- differentiability at matching;
- the local slab part;
- calibration;
- scientific acceptance of results.