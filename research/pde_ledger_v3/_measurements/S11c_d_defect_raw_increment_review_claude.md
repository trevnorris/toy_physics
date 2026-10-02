**NEEDS REVISION**

This is a source-only review. Nothing was run, so I can't say whether the worker restores or executes cleanly, or what it costs. Several of the joins the packet calls "actual" are one of three things: a hash match, a literal in the worker, or a boolean copied from the author's extraction.

What I checked and found sound:
- The upper factorization is exact. I re-derived `N = −H[qi² + k·qi(H+2k)/(qh+qi) + k·qh(2k+2Q−H)/(qo+qs)]` from `qh²−qi²=−H(H+2k)` and `qo²−qs²=−H(2k+2Q−H)`. That reproduces `Br` at `worker.py:271-273`.
- The contact term is zero at H=0, and the lower boundary equations are triangular.
- The four depths are bound to one outgoing rule.
- The `Rprod` direct-once identity is correct: `∂closed[0,2]/∂D = qi·qo/((qi+aρω)(qo+aρω))`.
- Pressure and normal-jet jets are dropped at (1,1) for the right reason. The row jet slots carry `eta_bg·w1_profile`, so the direct increment there is η²σ.

## Math / method / claim findings

1. **The lower-face geometry is a literal, not a native join.**
   - `worker.py:229-233` hard-codes `face=-1` and the lab slope `face*s`. It runs only the generic graph-normal assignment, so every `face` enters as `face²`.
   - The lower derivation reproduces the upper one by construction, and `lower-upper-outward-mirror` (`:255`) is the only external check.
   - `native['geometry']['face_normal']` and the minus `face_shift` case (index 4) are never read. `face_shift` carries `−½·W_0·d_w_delta_p_minus·eta·w1_profile`, which fixes the lab-height sign. The worker discards it with `ht0 = height.subs({eta:0,sigma:0})` at `:362`.
   - The sign of lab height is therefore only exercised by a control, not joined to the native source.
   - Smallest fix: add a zero check that the minus `height` coefficient equals `face·W_0·eta·w1_profile/2`. Add another that the native minus `face_normal[0]`, differentiated in the slope, equals `normal[0]'` at `:239`.

2. **The rows actually used are never joined to the native rows.**
   - `:386-387` joins only `old['raw']` to the sum of `selectedPressureJetChildren`.
   - `:388` then builds the increments from `old['bound']`. Nothing joins `bound` to `bind(raw)`.
   - `bound` is only the pressure and jet slot terms: 4 terms for THETA, 4 for E_W, and `0` for U0–U2. The emitted key `fullPublishedBoundRow` (`:391`) is mislabeled.
   - `native-rows/*.txt` are pinned but never read. "All pressure occurrences covered" is a boolean copied from the author's extraction.
   - Smallest fix: add `J.zero(row+'-bound-join', bind(old['raw']), base)`. Parse the pinned full-row text. Check that every `delta_p_*` and `d_w_delta_p_*` occurrence lies in the selected children and that the row is affine in those slots.
   - Smaller point: the `packet-index.json` row addresses for THETA and E_W end in `EXPANDED, 0`. The `native-selected.json` addresses end in `EXPANDED`. The worker's lookup at `:385` works only with the latter.

3. **The two zero-grade sources are identical by hash, not by a join.**
   - `plus-source-restrictions.json` and `minus-source-restrictions.json` have the same sha256 and contain no face symbols. The worker uses `S0` for the plus face and `SM` for the minus face (`:183`, `:378`).
   - `:348` joins only the native source to `input_<label>['raw']`. It never checks that `S0`/`SM` is the grade-(0,0) restriction of that face's own input. The two inputs differ in size by 4 bytes.
   - Smallest fix: add a per-face zero check that the bound `input_<label>['raw']` at eta=0, sigma=0 equals the face's source00, with `bind` and the jets zeroed. Otherwise relabel the minus source as inherited and unverified.

4. **The measure and edge-delta join is typed in, not executed.**
   - `:319-322` writes `factoredEdges: 2`, `sourceForwardNormalized: False` and `profileForwardNormalized: True` as literals next to the `profile_bindings`, `kernel_apply` and `fourier_profiles` source text.
   - The dimension sum at `:317-323` uses hard-coded hat and measure dimensions.
   - The method asked to "join the executable definitions, not an ambiguous comment", and to compare normalization to the native first-shape second slot. Neither is done.
   - Smallest fix: parse the `(2π)^-3` factors and delta counts from `profile_bindings` and `kernel_apply` and assert them. Or state in the output that this join is asserted and not verified.

5. **The emitted row density carries free depths.**
   - `rawRowDensity` (`:395`) substitutes `raw`, which still has `qh`, `qs`, `qi` and `qo` as independent symbols (and numeric ω, ρ). The sheet binding is emitted separately in `physicalRaw` at `:292`.
   - Nothing ties the row density to the shared sheet, so independent depth assignments could still be inserted.
   - Smallest fix: also emit `rawRowDensity.xreplace(physical_depth_bindings)`. Zero-check that the row density with the sheet substituted equals the row increment built with `physicalRaw`.

6. **Several controls are vacuous or sit where the k-terms vanish.**
   - `double-insertion-response` (`:406`) compares `baseline` with `2*baseline`, which is nonzero by construction.
   - `lower-reference-jet-sign` (`:416`) is a pure sign flip, which the plan admits. It tests the trace, not a consumed jet.
   - The wrong-height, slope-omission and wrong-sheet controls all use `sample` at `k=0` (`:401`), where every k-proportional term vanishes. They never touch the k-dependent contact numerator or the `qs` dependence.
   - The row controls (`:407-415`) are real, but they respond only through the pressure slot. The jet slot drops at (1,1).
   - Smallest fix: use `sample` with `k≠0` and consistent `qh`, `qs`, `qo`. Replace the hand-built double with a real doubled substitution into the row.

## Pure tooling and risk (not result-changing)

- **`:277`** reduces the numerator with `subs({qh**2:…, qs**2:…})`. This is correct only while the degrees stay ≤2. SymPy turns `x**3` into `y**(3/2)` here, which would not cancel and would fail the run. My count says the degrees stay ≤2, so this should pass. Use `sp.rem` or `sp.reduced` against both dispersion relations instead.
- **`:217`** emits `sp.Ne(qi,0)` and `sp.Ne(qo,0)` as `externalDomain`. Because the symbols carry `nonzero=True`, these evaluate to `True`. They are not provenance. The `qh+qi≠0` and `qs+qo≠0` domain exists only as a prose string at `:276`.
- **`:450`** (`posthashes` in `finally`): a missing pin raises inside `finally`, so the remaining files, the integrity flag and `checks.json` are skipped.
- `raw_lower = raw` (`:290`) rests on the coefficient mirror (`:255`). It does not rest on the native lower-face hats, which are unjoined (finding 1).

## Optional

- State explicitly in `checks.json` that the U0–U2 increments are the constant 0.
- Add the "addressed times-two (symmetrization)" control the method asks for. The current "times-two" control is only the insertion check.

Source-only method agreement here is not worker or result acceptance. Findings 1–4 are the ones that change what the output can claim.