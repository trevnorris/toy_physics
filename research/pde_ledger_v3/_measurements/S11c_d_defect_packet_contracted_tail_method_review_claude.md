**NEEDS REVISION**

I find the proposed separation mathematically appropriate in principle. The blockers below are gaps in what the proposal specifies, not objections to the idea. They are small, and I would clear the method once they are stated.

## What I checked

**The failure was an interpretation, not a numerical or physical event.**
- `prepare.py:212-216` applies `epsilon=1/(80e11)` to each analytic row. The aggregate guard at `:217-219` was never reached.
- `new-tail-address-budget-8360-J-27.json`: the total is about 9.6e-13, against a threshold of 1.25e-13. About 9.2e-13 of it is the middle (|t|>T) term.
- The 32 earlier rows (8346 through 8359, both windows) passed individually.
- The enlarged return for 8360 totals roughly 2.4e-13 (about 5e-17 outer, 2.4e-13 middle). The K29/T124 row would therefore also have breached the per-address threshold. The refusal comes from the carrier-uniform constant `E30=3^30` and the middle-tail term, not from the K27 choice.

**The original source supports the aggregate reading.**
- `preflight.py:339-348` stops when each of the outer, heightQ and middle totals is below 1e-11/3. It never imposes a 1/80 share per address.
- `baseline-method.md` §5 says "per-address/grade budget checks" without defining them, and the epsilon/80 is introduced as a numerical budget (§5, lines 277-283). So the proposal's reading is consistent with the source and was not tuned to the failure. The 1e-11 ceiling is the original number, not a new one.

**The enlarged formulas match the original `contributions` source.**
- `prepare.py:203-208` has the same right-hand sides as `preflight.py:331-336`, with the same `E30`, `F[2]`, `F[3]`, `b*`, `Cordinary` and `weighted_tail`. The failed run's AST join (`validation.py:109-119`) enforced this.
- I hand-evaluated the 8360 and 8361 middle terms from the formulas. They reproduce the saved values: about 9.2e-13 for J, 9.0e-14 for D, and a J/D ratio of about 10.

**The triangle majorant is applied per primitive.** `new-per-primitive-absolute-tail-transport.json` shows that each of Dr, Dh and Dq uses the full positive D envelope, so three copies per D address is legitimate. The Cordinary outer bound covers every component. The J bound keeps the H overcount.

**The 20 base rows and 20 enlarged returns exist.**
- Addresses 8346-8361 and 9672-9687, for J at 10 addresses and D at 10 addresses, give 10 + 30 = 40 primitives per window.
- I confirmed 20 `new-eligible-*` files. All 20 carry the unit vector [-2,-1,1].
- Bounds depend on the carrier only through the uniform |p0|<3 constants, so the two carrier sums will be numerically identical.

**The Dr window is the original window.** In `new-per-primitive-absolute-tail-transport.json`, Dr maps t = l−m. The |t|≤T clip on I(m) therefore matches the original middle region.

## Blockers (method text must state these)

1. **Scope of the 1e-11 ceiling.** The original 1e-11 covered all live addresses, including flat, height, slope and H. The proposal now spends a full 1e-11 on J/D alone, with D counted three times. It must say the ceiling is a J/D-subset ceiling, not a share of a global 1e-11. Pending components then need their own analytic accounting, with no donation either way. Their ceiling must also be compared with the 1e-9 floor of the final tolerance. This is harmless in size but must be explicit, or the later "all-real" claims will cite a budget that was spent twice.

2. **The monotonicity check is a regression check, not a certificate.** The formulas are decreasing in K and T for fixed CX and CY, so a pass is guaranteed by construction. A failure signals an argument or provenance error, not a physical problem. The method should say so and should not present it as an independent analytic result. It must compare the saved K29 returns componentwise (outer, middle, and the J H-term) against the base values.

3. **Enlarged-window values are unverified.** The proposal says to preserve the 20 saved enlarged returns and calls them observations, not an independent proof. It also forbids replay. It must decide explicitly whether a cheap exact-rational re-evaluation of the closed-form K29 expressions from the saved constants is allowed as a separate cross-check. I recommend allowing it, labelled as an independent exact arithmetic check and not a replay of the completed operation. If that is refused, the K29 aggregate rests entirely on the failed run's receipts.

4. **Units.** The sum is only meaningful if all 20 addresses share one unit and the 1e-11 and 1e-9 floors are in that unit. I confirmed [-2,-1,1] for all 20 in the saved eligibility records. The continuation must assert this before summing, and must state that the tail bounds are in original reference-unit numeric coordinates.

5. **Guard attestation for the prefix.** Only 263 operations are listed as complete. Many `require()` guards and `J.emit` records in `prepare()` are not operations (for example the accepted-contraction checks and the outgoing-source AST join). The continuation must say which of these are re-run, since they are cheap and not scientific replays, and which are attested from the evidence chain. It must also say how a mismatch is detected.

## Obligations if the method is accepted

**Coverage**
- Exact rational arithmetic with strict `<` against 1/10^11.
- Persist the 40 per-address/primitive bounds, the per-face, per-grade and per-component sums, and the full sum for each of 2 windows × 2 carriers. All decisions come after the saved operands.
- Count D three times and J once with its H overcount. No cancellation, no deduplication of shared families, no cross-carrier donation.
- Join `T≥K+4` and the `b*` and `CX`/`CY` provenance for all 40 rows. The `baselineOnly:True` label on the K29 rows in `prepare.py:212` is wrong and must not be carried forward as provenance.
- Coverage controls must refuse a dropped or duplicated primitive and a wrong window or face assignment, using the same census with no second baseline computation.
- The wrong-root and derivative numerical controls and all A24/A48/B50 and window comparisons are unchanged.

**Runtime**
- A new run directory that copies and verifies the failed bytes.
- The aggregate runs first. A failure there is preserved and stops everything, with no enlargement or retuning.
- Only after it passes may geometry and numerical quadrature run, under the unchanged guard and supervisor limits (4 GiB, 16 GiB pool, one CPU, 32 tasks, no deadline).

## Not computed, not claimed

I did not total the 40 bounds, and I did not read the minus-face K27 rows beyond 9686, which equals 8360. Reading a few saved rationals, I see 8360 J ≈ 9.6e-13, 8361 D ≈ 1.3e-13 per copy and 8350 J ≈ 7.6e-14. That is low-1e-12 order for the largest rows, but my digit counting is error-prone and this is not a pass. Neither the aggregate nor the enlarged monotonicity has been established. The constants (`E30`, Cordinary, the 36·121/b*² and 4/5·121/b*³ factors) and the strip and envelope inequalities are inherited. I have not independently re-derived them. The verdict clears none of the implementation, runtime, or computed science.

## Files read

- Fully read: `method.md`, `evidence-guide.md`, `baseline-method.md`, `prepare.py`, `original-source/preflight.py`, `failure/failure.json`, `failure/new-tail-address-budget-8360-J-27.json`, `new-tail-address-budget-8359-Dq-29.json`, `new-enlarged-tail-8360-input/decision/return.json`, `new-enlarged-tail-8361-decision/return.json`, `restored-base-tail-8361.json`, `restored-base-tail-9686.json`, `failure/new-per-primitive-absolute-tail-transport.json`.
- Partly read: `saved/preflight/tail-plan.json` (lines 1-1518 and 2380-2470 of 2470), `saved/preflight/tail-selection-prefix.json` (lines 1-1207 of 1958), `worker.py` (lines 100-205), `validation.py` (lines 95-135), `numeric.py` (lines 40-90), `saved/preflight/Fourier-envelope-constants.json` (grep only).
- Directory listings and grep counts only: the failure directory (100 of 1571 paths), `failure/new-eligible-*` (unit-vector count only), `saved/preflight/`.
- Not read: the other 18 enlarged returns, the other 17 base-tail records, `new-eligible-*` contents, `checks.json`, `evidence-chain.jsonl`, the exact-joins directory, `saved/absolute-bounds/`, `original-source/{inner,fourier,contraction,contraction_lib}.py`, `restore.py`, `geometry.py`, `evidence-store.py`, `request-index.py`, and `review-prompt.md`.

No essential operand is missing for the method-level judgment. The 40-bound aggregate itself needs a computation that this review did not do.