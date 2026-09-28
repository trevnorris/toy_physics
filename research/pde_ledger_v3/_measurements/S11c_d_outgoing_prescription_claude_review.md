Finished reading the packet; here is the review.

# Verdict: **NEEDS REVISION**

The method is sound for this fixed slice. The pole sign, the multiplicity-two residue formula, the branch and coordinate joins, and the inherited current-orientation rule are correct. It needs no global no-pole theorem, frequency campaign, A12 witness or new census. The real-axis integral already carries the complex poles and the imaginary-axis branch cut, so justification depends only on real-axis data.

Two bounded gaps remain in that real-axis data. Both can be closed inside the worker, using only saved operands and no new roots.

## Substantive blockers

**B1. The convergence/distribution domain is not established, and the saved candidate carries no domain.**
- `S11c_d_outgoing_prescription.py:290-298` builds `pv` as a `Limit` of improper `Integral`s over `(-oo, k1-ε)…(k2+ε, oo)` of `phase*fixed_inverse/fourier`.
- Nothing computes how `fixed_inverse` behaves as |k|→∞. The pointwise improper integral exists only if every entry tends to 0.
  - An entry that tends to a constant or grows gives a δ/δ′ contact term at z=zp, and its Riemann integral does not converge even when x≠0.
  - An entry that decays like 1/|k| (plausible because q_out ~ i√b|k|) gives a logarithmic singularity at x=0.
- The x=0 exclusion exists only as a comment (`:304-305`). The saved `candidateKernel` (`:322`) is stored as if it held for all `z, zp`.

*Repair:* for all 25 bound inverse entries, compute and save the exact leading large-|k| order separately for k→+∞ and k→−∞. Substitute k=±1/t and √(a+bk²)=√(at²+b)/t, then expand in t.
- Require each entry to tend to 0, and save `Ne(z, zp)` as the candidate's pointwise domain.
- If any entry is O(1) or growing, split off its polynomial part as explicit contact distributions at x=0, using the saved Fourier mass, and apply the PV only to the remainder.
- If the leading order is 1/|k|, record the diagonal as excluded.

**B2. Transfer of real-pole coverage to the actual inverse is incomplete in three explicit links.**
The underlying argument does work. On the strictly evanescent slice (checked at `:155-156`; here q² = 0.01 − 0.05 − k² < 0), a real k forces an imaginary q. The saved axis records show only factor 1 has imaginary-axis roots, which gives disks 7 and 8, lifts 14–17, and selection of 16/17. What the worker does not check:
- **(a) Off-sheet exclusion is numerical only.** Records 14/15 are real-normal but off-sheet. They are dropped only through the tolerance-based `SHEET_MEMBERSHIP` flag (`:201-204`); the exact native join at `:240` runs only for selected records. *Repair:* for every `PROVED_REAL` record not selected, require `simplify(q_out(K) + scale·Q) == 0` exactly, using the saved exact lifts.
- **(b) Singularities from denominators are unchecked.** The coverage certifies zeros of the modal determinant's *numerator* along the curve. `DENOMINATOR_EXCLUDED` (`source-definitions.md:407-409, 496`) only rules out a shared gcd; it does not rule out real-k poles in the entries of `fixed_det` or `fixed_adj`, and those would put further singularities in the PV integrand. *Repair:* require that the denominator of `together()` for `fixed_det` and each adjugate entry contains no `kn`, apart from powers of `sqrt(-square)`. Those are nonzero on this slice.
- **(c) The two singular objects are never joined.** Residues come from `det`/`adj` (`:234`), but the PV density comes from `fixed_inverse` (`:292`). *Repair:* require `fixed_inverse*fixed_det − fixed_adj == 0`. This also removes any doubt about the orientation of the "transposed cofactors".

## Verified as correct
- **Pole sign.** Fields are built as `exp(+i kn z)` (`source-definitions.md:108`), and `fourier_mass` = 2π (`:42-46`). So G(x) = (1/2π)∫e^{ikx}M⁻¹ dk, and 1/(k−k_a−i0·s) = PV + iπsδ.
  - J is the +z flux, being the coefficient of `phase_rates[1]` (`:209, 347`).
  - `:273` and the identity at `:306-315` are right.
  - The rule matches the inherited one in `S11c_d_continuum_boundary.py:198-203`.
  - In the saved data, k=+0.785 has J>0 and k=−0.785 has J<0 (modal checkpoint `:330-349`).
- **Multiplicity-two residue.** m·adj^(m−1)/det^(m) follows from matching Taylor coefficients. adj(k_a)=0 is equivalent to rank ≤ 3, and with c invertible (`:251`) the block is semisimple. The Keldysh comparison R c⁻¹ L* at `:253-254` is correct and independent. `NORMAL_PENCIL_PLUS` is the total derivative along the branch (`:555-560`).
- **Branch and coordinate joins.** Im q_out > 0 matches decay under `exp(+i q depth)` (`:280`) and `BULK_DECAY_DISK_CERTIFIED` (`:747`). The full-symbol join is exact, and the native↔physical q join at the selected poles is exact.

## Optional observations (non-blocking)
1. `implementation.md:90-91` says the integrals use "the actual saved spectral density". The code instead rebuilds the density with a re-typed phase (`:291`) and never reads the accepted `regular-spectral-density` artifact. The phase is correct, but either join the artifact or fix the wording.
2. Three checks are wiring tests, not evidence: `outgoingCurrentSign` (`:269`), the direction mutation (`:272-275`) and the scalar identity (`:306-315`). The same applies to `projectedDerivativeJoin`, which compares two identical saved formulas.
   - A non-trivial local check from saved data: the eigenvalues of −(FLUX_LEFT* ∂ₖP FLUX_RIGHT)⁻¹ (dk/dω, because the flux-frequency normalization is I) should be real with sign equal to `direction`.
   - That would confirm agreement with local ω+i0 displacement at each real block, without claiming any global equivalence.
3. `fixedPhysicalInput` is saved with `eta_bg=1/100`, but `binding()` sets `eta_bg=sigma_W=0` (`:135-136`). Record the reference binding explicitly.
4. `candidate-census.json` is written only after the guards inside the loop (`:197-213`). Save it per record, or before the guards, to meet the durability rule.
5. The exact-zero tests use `sp.simplify` at nested-radical points (`:230-240`). Failure is safe, since it stops or is caught by the numerical residue comparison, but it may be slow within 900 s. Zero tests based on minimal polynomials would be more robust.

## What turns the candidate into a result
If B1 and B2 pass, the object is a justified spatial outgoing kernel for z≠zp at this fixed input, with any contact terms explicit. It would still not show ω+i0 equivalence, radiating-domain (A12) coverage, other inputs, or the profile response. A clean process exit does not accept the numbers; the saved residuals and domain records do.
