CLEAR FOR THIS SELECTED CENTRE-DRIVE METHOD

I found no blocker. The hand candidate checks out against the saved operands, and the saved face maps are consistent with the plan's tests. This is a source-compatibility verdict only. It is not worker, runtime or science clearance, and it settles nothing about the holding law, centre uniqueness or later work balance.

**Candidate algebra**

I expanded the saved a operand for `(LAB_HELD, DELTA_W, ZETA_C, RHO4_CONSTANT)` by hand. Here `i·Λ_X_0/(ωτ_X+i) = Λ_X` and `Λ_X = Λ_X_0/(1−iωτ_X)`.

- **UPPER face:** `−χ[P₊ + ½ηw₁W₀·DwP₊] − Λ_X·μ/ρ_br·(1−ηw₁)`
- **LOWER face:** `+χ[P₋ − ½ηw₁W₀·DwP₋] + Λ_X·μ/ρ_br·(1−ηw₁)`
- **Sum:** `χ[P₋−P₊ − ηW₀w₁(DwP₋+DwP₊)/2]`, with `χ = 1−Λ_X/ρ_m`

The μ terms cancel in the pressure term and in the normal-jet term. This follows because the single symbol `mu_theta_L` appears with opposite sign on the two faces. The stored b row (`views/selected-b-center_face_generalized_row.txt`) has no μ at all and matches term by term. So the candidate is confirmed, including both pressure differences and the laboratory-derivative sum.

The sign convention is also consistent with the saved operand:

- **Pressure:** the plus face enters as `−P₊` and the minus face as `+P₋`.
- **Units:** the row has pressure dimension `[−1,−2,1]`.
- **Source extraction:** `face_generalized_force_rows` takes `diff(center_density, delta_v_zeta_c)`, so the row has no support measure.
- **c2 restriction:** the c2 audit never references `CENTER_FACE_GENERALIZED_ROW` or `ZETA_C`, so the omitted row is preserved.

**Saved affine maps**

- **Pressure:** the plus and minus column-0 and column-1 pressure strings are textually identical.
- **Normal jet:** the outward derivatives are identical (`+iqP`) on both faces. The lab derivatives are `+iqP` on plus and `−iqP` on minus, set by `sign·iqP` at `S11c_d_first_order_receiving_blocks.py:272`. So `DwP₊+DwP₋` is zero by the saved sign map.
- **Tag:** the shared tag `receiving_mu1_column_c_hat_l_minus_p` has real shared ancestry. The plus and minus `inheritedSourceNormalization` operands are textually identical, and both are joined to one per-column `chemical-c-return`. Column 1 has envelope 0, hence no tag.
- **Baseline source zeros:** `S00 = S10 = 0` is enforced by `require` at line 254. The remaining direct source is `S01`.
- **Other zeros:** the baseline chemical zero, the published zero direct velocity, and the `face-P/mu/V-transverse` zeros are all checked in the producer.

If the new cross-face identities pass, they give a compatible zero only for this selected response class. The zero comes from the symmetric DELTA_W face dictionaries and the outgoing lab-sign map. It does not say anything about free-centre uniqueness.

**Non-blocking cautions for the worker spec**

1. **Path weighting in the tag.** `directChemicalEnvelope` is `independentSigmaCoefficient/10` (producer line 269), so the saved tag already carries the `σ=λ/10` path weight. The plan's "independent η/σ until the optional path" can only be honoured for this tag by also recording `chemical-c-return.independentSigmaCoefficient` and stating the factor 10. The row itself has no σ in its pressure or normal-jet coefficients, because σ enters only through the μ binding that cancels.
2. **Two tag families.** The affine maps use `receiving_mu1_column_c_hat_l_minus_p`. The pressure assembly (line 255) uses `source01_{face}_column_c_hat_l_minus_p`. The plan should name the join between them and the single-count flags `directPressureCount` and `sourcePressureNotAddedAgain`, so the direct chemical enters once.
3. **C10/C01 on the baseline.** In the saved data this is `S00 = S10 = 0`, plus zero `(1,0)` and `(0,1)` consumers in the five-row splits I sampled (THETA delta_p). The worker should restore those zero operands rather than assume them.
4. **Normal-jet control.** On the selected first-order maps the normal-jet sum is identically zero, so the plan is right to call that control inapplicable. The formal row-level `η·W₀·w₁·χ·DwP₋` sensitivity is live, and the one-face pressure sign mutation is live generically.

**Coverage**

- **Inspected:**
  - `plan.txt`, `guide.txt` and `review-prompt.md`
  - the selected a work, b centre row, μ binding and case-dimensions views
  - b `face_generalized_force_rows` and `bind_mu_theta_operand`
  - the c2 centre-row grep
  - affine plus/minus column-0 views (pressure, both normal derivatives, direct chemical tag, normalization), column-1 pressure and envelope, and minus-0 totalMu and plus-0 totalV
  - the plus-0 `nativeFace` key list and the short-field text lines (omitted-line fields such as `PRESSURE` were not read)
  - producer source for the receiving blocks (lines 225–283) and the first-order source (lines 140–264)
  - the THETA delta_p split files
- **Not inspected:**
  - the full `virtual_work-full-display` and the other physical ZETA_C cases
  - the raw/input/return JSON records and file hashes
  - the remaining affine views (mu00, V00, minus-1 totals)
  - the `*face-P/mu/V-transverse` and `original-lab-normal-join` saved records themselves (I read only their producer checks)
  - the saved background-force-pairing, measure/support registries and native directives
  - the dimensions of the full b operand tree beyond the stored row
  - any baseline-incident face record beyond the zeros listed above

The support measure, total `R0_c`, the THETA work pairing and the holding law remain open and are not assumed here.