## Verdict

**CLEAR FOR THIS BOUNDED SOURCE-COMPOSITION METHOD**

I found no blocker. This is an assessment of the method text against the operands I opened. The "Not verified" section below lists what I did not check. I ran no code and left no files in the packet.

## What I checked against the raw operands

- **Direct correction versus native c2.**
  - `source/native-c2.py:518-525` puts a literal 0 at `trace_three[0,2]`, and `:533` returns `reference_three[0,2]`.
  - `raw/plus-closure-operands.json` carries the independent `plus_raw_direct` symbol at `[0,2]`. Its linear coefficient in the closed matrix is 1/((1+a_o)(1+a_i)), with a = 3c/(10q).
  - By hand, 3c/10 = (30+9i)/109 = β = 3/(10−3i). So the coefficient is R(qo)R(qi) with R = q/(q+β), and it matches the saved `factor`.
  - The method's labelling of the direct term as augmentation, not native c2, is correct.
- **External resolvents are already in the closed density.**
  - `direct/closed-density.json` has the denominators (q_i+β)(q_o+β) in `Bc`.
  - Its prefactor is sinh(5π·td)·sinh(5π(l−td−k)), which matches the method's td-route map: height momentum k+td, slope momentum l−td.
  - So the "no second resolvent pair on Dwhole" rule is the right one.
- **Reference and jet factors.**
  - `raw/*-closed-before-cancel.json` has reference = q_i q_o(60+91i)/(…).
  - I checked that (60+91i)·β = 9+30i and (60+91i)·β² = 9i. So reference equals Rprod exactly.
  - The jet is +i·q_o·reference on the plus face and −i·q_o·reference on the minus face. That matches the method's i·f·q(l) and "put it before F".
  - The raw kernels `rawKernelPlus` and `rawKernelMinus` are textually identical in `U0-retained-increment.json`, so the face dependence enters only through f and the face-specific source and consumer. The method requires separate per-face joins and does not infer them, and that requirement is necessary.
- **Native trace and slot routing.**
  - `plus-new-native-trace.json` shows T00 = 1, height_constant = 0, and T[0,2] = 0. The method's claim that a `[0,2]`-only increment leaves T01·reference12 unchanged is therefore consistent.
  - In `build_face` the retained pressure is pressure − H·normal_jet. With normal_jet = i·f·qo·ΔP this gives retained ΔP plus an excluded (2,1) piece, as the method says.
  - `plus-final-native-slot-routing.json` has the direct placeholder at zero and `sourceCompositionPerformed: false`. It correctly does not discharge the new routing obligation, and the method does not claim it does.
- **Fourier convention.** `raw/native-fourier-contract.json` and `raw/edge-delta-reduction.json` show the profile forward transform at (2π)^-3 and the source inverse at (2π)^-3. Two edge deltas leave a normalized 1/(2π) one-dimensional factor. That is consistent with the method's b̂ = (1/2π)∫e^{-itx}b, plain dk dl, and constants giving b·δ.
- **Counts.** The triple counts are right: 9 + 3 + 3 + 1 = 16.
- **Saved retained-increment records.** They show `wholeConvolutionInserted: false` and `middleIntegralEvaluated: false`, matching the "no middle integral, no evaluated convolution" scope.
- **Compact-domain handling.** The method says explicitly that [-3,3] response certificates do not give global composition, and it marks the global-momentum obligations UNRESOLVED. It makes no grazing-limit claim and inserts no cutoff. That scoping is honest.

## Blockers

None.

## Non-blocking obligations the instrument must satisfy

These are tooling and validation items, not scientific objections.

1. **Join the minus face separately.** Do not inherit it from the plus face. The packet shows plus-only raw kernels, and the sign f enters only through the jet. The minus trace must show height_constant = 0 and T00 = 1.
2. **Map the raw symbols to the density symbols.** Join `increment_difference` to l−k and `increment_transfer` to td, and join the raw q_h and q_s to the density's grazing_qh and grazing_qs. The sinh arguments are consistent, but the map is not stated in the method.
3. **Apply the adapter to H as well as ΔP.** The (2,1) normal-slot remainder is the old unprojected record. Keep it as excluded data and do not use it as a match target.
4. **Do not treat an all-zero control as responsive.** The U0, U1 and U2 rows have no pressure slots. The addressed nonzero rows must come from THETA_BALANCE and E_W_BALANCE, and the instrument must show they are nonzero per face and slot.
5. **Record the reference-level ΔP injection.** Injecting ΔP as the physical `pressure` rather than as reference-level differs only at grade (2,1). The adapter record should say which one it uses and that the difference is excluded data.

Optional wording: none needed.

## What a successful instrument would establish

- A both-face, source-addressed inventory of the ordered noncommuting composition at all 16 grade triples.
- The direct whole-convolution tag entering once, with its augmentation label.
- Formal routing identities, and tag-sensitivity controls that are not field values.

## What stays unresolved before a finite near-unity defect pilot

- Global-momentum well-definedness, test spaces and quadrature for the profile-multiplied composition.
- Conditional convergence and PV interchange.
- Any grazing or delta→0 limit.
- The old cutoffs (4/6) and the numerical matrix-route and momentum-zero debts.
- The physical value of the direct addend, which the sensitivity controls do not establish.

## Runtime evidence needed

- The new formal `build_face` adapter join, with the actual assignment text for both faces.
- The missing source and consumer cross-grade coefficients, derived from the saved rational quotients.
- Per-face nonzero control rows.
- Exact residuals persisted before any zero is required.

## Does saved evidence support no-replay implementation?

Yes for the operands I opened, but no new scientific result exists yet.

## Not verified

I read `method.md`, the evidence guide, and only the first 758 of 1,196 lines of `packet-index.json`. I did not open the following:

- The consumer, source and census files.
- `function-routes.json`.
- The minus-face trace and slot-routing files.
- The THETA and E_W increment values.
- The density body.
- The continuation controls.

These conclusions rest on the files named above.