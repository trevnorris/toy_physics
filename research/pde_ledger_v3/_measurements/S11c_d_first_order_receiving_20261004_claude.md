**NEEDS REVISION.** The revision is narrow. The 3×3 receiving block and the A(x) step argument are correct and well posed. The plan does not yet include the transverse sector, which is exactly resonant at the grazing points and decides whether the deficit starts at O(λ²).

## What I checked by hand against the saved operands

- **Block separation, both directions.** In `LEFT-invariant-P.json`, with t and b(l) written in uniform order:
  - The elastic block acts as `1.5K2−ω² = (3/2)(l²−p²)` on both t and b(l), and as `(17/30)K2−ω²` on k(l). So D_T is correct on both transverse columns.
  - Columns 4 and 5 are proportional to k, with coefficients `i(4n²+47)/12` and `i(20n²+31)/200`.
  - Rows 4 and 5 are `(1/5,1/10,n)·X`, also proportional to k.
  - So t and b decouple exactly in both directions for all l, ω and depth. The 3×3 block on (longitudinal, θ, e_W) is even in l.
- **The q→0 limit.** `q²=p²−l²` matches `receiving-sheet.json`. The 1/d closure terms have a removable cancellation, so the determinant at d=0 does not depend on the sheet. The sheet matters only for the √Q expansion around it.
- **A(x).** The saved `etaForce` equals `A(x)·U` componentwise for both columns. Here A = 2T²+4.5T+2.5 = 9w−6m, with w=(1+T)/2, m=(1−T²)/3 and T=tanh(x/10), as in `restored-profile-rules.json`. A(−∞)=0 and A(+∞)=9.
- **The identity.** With `F̂=∫e^{−iQx}F`, `QÂ=−iÂ'` holds, and Q·δ(Q)=0.
- **Source projections.** `k·U_A=iQ/5`, and `k·f_B=0` for all l (U_B's force stays ∝ t at every x). For U_A the longitudinal source at Q=0 is (−9i)(i/5)=9/5, which is nonzero. The pointwise zero of `k(p)·f_R` hides this. Rows 4 and 5 of `etaForce` are zero, so Â never reaches the scalar rows unmultiplied. The chemical-0 and `source01Envelope` polynomials vanish at T=±1. I did not verify the sigma-force rows or the other pressure-assembly rows directly.
- **Face algebra.** The P formula follows from the native identities at `audit.py:3547-3548`: `affinity=μ−P/ρ_m`, `j=ρ_m(v_bulk−V)=A·affinity+V_mem·V`. The flat factor at q=0 is 3/(10β)=1−3i/10, matching the saved value.

## Blocker

**The transverse sector is resonant exactly where the problem sits, and the plan excludes it.**
- D_T(±p)=0, and Â carries a 1/Q pole. So `Â/D_T` has a double pole at l=p. That is the secular wavenumber shift, and it appears for both polarizations.
- U_A has both t and b components at l=p: `t·U_A=−ip/10` and `b(p)·U_A=1.2i`. Because b(l) varies with l, the b-sector pole also gives a polarization correction. That correction is the "matched-end polarization correction" the plan asks to check, and it cannot be checked without this sector.
- Inverting the scalar D_T is not a pencil inversion, because the blocks are exactly decoupled. The revised deliverable should include:
  1. a stated radiation or limiting-absorption prescription for the real poles at l=±p;
  2. the O(λ) forward 2×2 solvability integral on the doublet, including the c_b′(p) term;
  3. the reflection weight at l=−p;
  4. a Hermiticity check in the physical flux metric.
- Step 4 decides whether the deficit is O(λ²). For U_B, the imaginary part of the relative multiplier is odd and a total derivative, `Im a ∝ −(p/15)T(1−T²)`, so its integral vanishes. I have not shown this for U_A.

## Physical input needed, separate from derivation

- **Incident state or polarization.** The saved current Gram is non-diagonal (`−1.339`), and only U_A drives escape. The leakage is a rank-1 form in the (U_A, U_B) amplitudes. A coefficient needs the incident state or a generalized eigenproblem in the flux metric. U_A is not flux-normalized.
- **Power balance.** The O(λ²) deficit equals escape plus memory dissipation only if a conservative power identity holds. The held profile may also do external work.
- **Real-axis census.** The elastic longitudinal diagonal vanishes at K2=270/17 (|l|≈3.98), outside the radiating window. Whether the couplings damp it is open. Â there is exponentially small but nonzero. Threshold regularity of the 3×3 block does not cover this.

## Cautions on the face maps (not blockers)

- ρ_m is the right density in P, j and affinity. ρ_br appears only in the chemical-driver denominator. Keep μ=driver/ρ_br separate from ρ_m.
- The finite flat factor at q=0 comes from Re β=30/109>0 (the memory kernel), not from a cancellation. Without memory it would go as 1/q.
- The pressure addition `R00·S01` is already in the forcing. Reconstructing P_total=R00[V+(A/ρ)μ_total] contains it once. Use one or the other in any power sum.
- The saved orientation control (`orientation.disagreement=true`) must stay in the controls.

## Coverage

I read these in full or in the relevant part:
- `plan.txt` and `guide.txt`
- `incident-columns.json`, both force files, both chemical returns, `receiving-sheet.json`
- `LEFT-invariant-P.json` and `actual-source-chart.json`
- the w and m rules in `restored-profile-rules.json`
- the first 945 of 4398 lines of the pressure assembly (U0 only)
- the native audit around lines 3300–3390 and 3535–3565
- the chemical-driver view, part of `native-dimensions.json`, and the `rho` hits in the governing document

I did not read:
- `LEFT-native-1-original.json`, `acoustic-face-records.json` and the memory kernels
- the local and consumer source files and `S11c_d_first_order_source.py`
- the other pressure-assembly rows and the field-argument and control files, plus the rest of the source-stage artifacts

I computed no determinant. The 3×3 block's threshold determinant is still unevaluated, as the plan intends.