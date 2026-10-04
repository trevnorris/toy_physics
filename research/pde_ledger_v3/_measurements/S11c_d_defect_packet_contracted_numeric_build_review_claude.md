**CLEAR FOR THIS BOUNDED CONTRACTED NUMERICAL BUILD**

This is a source-level verdict under the coverage limits below. I executed nothing and found no blocker. The mandatory runtime checks and scope limits below stay open.

## What I read, and what I did not

**Read in full:**
- `build.md`, `method.md`, `evidence-guide.md`
- `worker.py`, `prepare.py`, `numeric.py`, `geometry.py`
- `request-index.py`, `evidence-store.py`, `launcher.py`, `restore-library.py`
- `tooling-tests.py` and both logs, `execution-authority.json`
- `original-source/S11c_d_defect_packet_geometry_lib.py`
- `numeric-source-contracts.json`, `tail-bound-derivation.json`
- `new-contraction-definitions.json`, `saved/rules/B-G7-K15.json`
- `address-8347-input.json`, `absolute-bounds/D-reflected-numerator-envelope.json`

**Read in part:**
- `input-manifest.json`: header, pins, resources, library paths, and byte-size greps.
- `tail-plan.json`: first entries.
- `address-8346-return.json`: lines 1–130.
- The 8350 address input and the field, polynomial and reconstruction files for 79ea…: targeted greps and reads only.
- The label-transport arrangement file: an excerpt.
- The original preflight and inner-lib sources: targeted greps for `exponential_moment`, `weighted_tail`, `E30`, the budget loop and the profile.

**Spot-checked operands:**
- J entries 8346, 8350, 8354, 8358, 8360, 9672, 9676, 9680, 9684, 9686: source and consumer constants, wave multiplier, status.
- D entry 8347 in full.
- 12 sampled zero entries: all have a zero source field and status `EXACT_ZERO_SOURCE_JET`.

**Not read:**
- The other nine D entries (8351, 8355, 8359, 8361, 9673, 9677, 9681, 9685, 9687) were not individually inspected.
- About 3,000 of the 3,142 runtime files were not opened. That includes the 17 factor-proof sets, the 20 template contexts and injection proofs, and most field files.
- `runtime-source/*` (guard, supervisor, containment helper, hook), `packet-action-method.md` and the source-contract JSONs were not read, so the unchanged-guard claim rests on the manifest and gate pins only.
- I could not verify SHA pins, mpmath or sympy internals, or any runtime behaviour.

## What I verified

**Constant quotient and eligibility** (`prepare.py:109-118`, `:101-107`, `:131-133`)
- Numerator and denominator are rebuilt and joined to the address field, `constantValue` and the old zero reconstruction proofs. All 34 fields are checked at `:72-76`.
- The census must be 64 applicable, 20 live, 10 J, 12 `EXACT_ZERO_SOURCE_JET` and 32 `EXACT_ZERO_CONSUMER`.
- Normal = 1, grades 00→11 and ε = 1 are required.
- `fields.json` carries `constant: true`. Wave multipliers match `P_j(ik)^n` in the data: `-composition_p²` for n=2, `-1/25`, `-1/100`, `-3i`, and `1` for e_W.

**Fourier selectors and signs** (`numeric.py:87-107`)
- X uses ν = k−p0 with center −5/2. Y uses ν = p0−l with center +5/2 and no conjugation.
- Both match the original `fourier-product` and `constant_reference` source.
- A carries k^n and B carries (ik)^n/i^n, so `alpha = b·c·P·i^n` carries i^n exactly once.
- The B recurrence (`moment_polynomial`, M1, M(r+1)) matches the identities in `method.md`. I confirmed n=2 by hand against `prepare.py:137-143`.
- Native coefficients are used for the checks only. Alpha multiplies once, after the contraction (`numeric.py:335-342`).

**Contraction assembly**
- `combine` (`numeric.py:110-116`) matches all five saved definitions. It is also joined symbolically at `prepare.py:147-151`.
- `inner_node` (`numeric.py:217-235`) has input clipping on J, Dh, Dq and the wrong-root mutant, and output clipping on Dr only.
- The wrong-root formula includes the unclipped `Y0C`/`Y1C`.
- The profile is textually identical to the original, including the sinhc series.
- The root's radicand and branch are AST-joined (`prepare.py:60-66`).

**Geometry** (`geometry.py`, with Quad/Line semantics from the geometry lib)
- Six affinity slabs, with the clip lines per slab as specified.
- The `d_g` branch table matches `method.md` slab by slab.
- Carrier and profile offsets are exact.
- All pairwise graph intersections inside each slab become cuts; nothing is pruned.
- The audit rebuilds the true max/min window at both ends and the midpoint of every slab.
- Each of the eight wing mutants sets the central window on an outer slab and must fail with exactly `'actual max/min clipping window'`.
- The contracted variables' singular loci are fully covered, including the label-transport targets (m = ±κ, ±(T−K), ±M, k/l = ±κ).

**Routes** (`numeric.py:237-401`)
- Route A jacobians: GL node to z = (n+1)/2, weight·2·half·z/2, positive.
- Route B: weight·half, with the embedded Gauss weight times jac.
- Inner A discrepancies are summed panel by panel before the products.
- B pointwise target is ε(1+1/|q|)/(16(3M+4)), with a separate weighted guard of ε/4 on K and G weights. The allocation lemma checks (π + 2acosh(M/κ) < 4 + M).
- The B outer loop is global: largest leaf contribution to the worst component, with parent and child records kept.
- No retry, escalation or timer appears anywhere.
- The B derivative mutant replaces Q by the constant (i p0)^n at `numeric.py:102`, so the moment list collapses to `[1]`. That is a real change to the recurrence input.

**Comparisons and controls** (`worker.py:132-183`)
- Every address/primitive, face/primitive, grade, J, D and J+D comparison is run per plan, for A24/A48, A48/B50 and the enlarged-window pair.
- Controls use 8347/Dr and 8350/J through all routes, windows and carriers.
- Movement must exceed 10× (errors + cross-route + enlargement).

**Tails** (`prepare.py:164-185`)
- `wt` equals the original `exponential_moment`/`weighted_tail` formula as I read them.
- The outer and middle expressions match the original `contributions`. The `F3` and `E30` values agree with the saved moments.

**Requests, journal and containment**
- Full request equality, route/purpose/precision namespace joins, PENDING refusal and the 8 MiB / 20 GiB gates all match the code.
- Failure prefixes, no resume, posthashes and the guard parameters match the code and `launcher.py`.
- The gate requires the exact literal verdict string. Mine is above.

**Inherited versus new**
- *Inherited, not re-derived:* the four contraction proofs, factor and normal joins, native injection proofs, units (including inferred gamma), the strip envelopes and the absolute numerator lemma.
- *New joins:* constant quotients, wave-on-delta, family summand, centered moments, Y phase, `combine` against the saved formulas, implemented parameters, root AST, enlarged tails, geometry and requests.

## Non-blocking deviations
These do not change what the pinned operands compute. I am not asking for a revision loop.

1. **Selector gate (`prepare.py:119`).** `argumentDerivative: 0` is written into the specs rather than read from `matchingX/YFamilyInterfaces`, so nothing refuses a non-zero value. I grepped all `accepted-units/summand-*-return.json`: no non-zero value exists, and the inputs are hash-pinned.
2. **Source joins (`prepare.py:60-66`).** Only `q` is AST-joined. The profile and `wt` are not, and `wt(·,27,·)` is not cross-checked against the saved K27 values. I verified both by reading.
3. **Label provenance (`prepare.py:189`).** The matching-K27 label file is attached to all four plans, and no per-primitive transport record is emitted. This is provenance only; the geometry is complete.
4. **Early-failure receipt (`worker.py:188`).** The `finally` selects from `contracted_leaves`, which is created only in `Evaluator.__init__`. A failure before that point loses the final receipt. The original error stays in the chained traceback and the `failure` record.
5. **Wrong-root labels (`numeric.py:339-340`).** `addressed` labels the wrong-root value as `8347/Dr`, `/Dh` and `/Dq`, so mutant D subtotals triple-count. Only `8347/Dr` is meaningful.
6. **Codec check dropped (`numeric.py:31-41`).** `decode` no longer re-checks the restored MP tuple as the original `inner_lib.py:92` did.

## Mandatory runtime checks

- **Preamble refusals.** Any `exact()` join, tail-domain check or geometry audit refusal ends the run before numerics. The K27 aggregate must pass `<1e-11` with D bounds triple-counted and the J bound carrying the H overcount.
- **First-use paths.** The tests are stdlib-only and never run `Evaluator` with mpmath. The first use of `encode`/`decode`, `Estimate` operator dispatch with `mpc`, the `mpmath.mp.clone()` contexts, the embedded-weight match (7) and the adaptive loop is in the native run. I desk-checked them.
- **Budgets.** A nested ≤ ε/4, B pointwise and weighted guards, K–G, comparisons, and all 12 control rows (`worker.py:132-183`) must pass.
- **Storage and memory.**
  - Every request, including each node's transform checks, stores the full `exactFamilyContexts` (`numeric.py:153`).
  - Each A inner node commits about 10 records and costs 10 FULL-sync transactions.
  - Source constants alone put A48 at at least thousands of outer nodes times at least about 2,000 inner nodes each, on the 20 GiB reserve.
  - A capacity refusal at an early A run is plausible. Treat that as a terminal safeguard, not a retune.
  - All parsed saved inputs (`raw`, about 100+ MB of JSON) stay resident for the whole run, so watch RSS against the 4 GiB limit.
- **After the run.** Audit the complete chain with `EvidenceReader.audit(expected_records, expected_head)`. Verify all posthashes and copies.

## Scope limits

- Route B is not an independent formula for any Y transform or any n=0 X transform, which covers 16 of the 20 entries. There it equals Route A at a different precision. Only n=2 X (4 entries) exercises the moment recurrence.
- The Gaussian identity, root, profile and geometry conventions are shared between routes.
- Indicators are empirical, not rigorous enclosures. The tails are conditional on the assessed bounds.
- The finite-window controls do not cover all-real mutated responses.
- Nothing here establishes H, flat, height/PV or slope, or any packet action, current or leakage factor. J is not full native mixed.