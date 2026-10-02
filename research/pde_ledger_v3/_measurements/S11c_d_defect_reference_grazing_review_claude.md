**CLEAR FOR THIS BOUNDED REFERENCE-GRAZING METHOD**

I found no substantive blocker. I did this by hand from source and JSON, with no code run. I re-derived each formula and re-checked each constant against the packet sources.

## What I checked

**Source identities**
- **First-shape kernel.** `dtn_first_kernel` (native-c1.py:606–618) gives iρω(h·qo·qi + h(k²−κ²) − i tilt·k)/(qo·qi). After the on-shell reduction (k²−κ² → −qi²) this is exactly (1): Zh = iμ(qo−qi)/qo·h, and Zj = μk/(qo·qi)·j with k the input momentum. The middle-leg maps at native-c2.py:404–405 give the same assignment for (l,m) and (m,k).
- **Closure.** The operator-level second-order term of (I+aZ)⁻¹Z is −a·G·Z1·G·Z1·G with G = q/(q+β). That matches (2), including both orderings and no extra symmetrization. The saved `closed[0,2]` entry (minus-closure-operands.json) has the same z01·z12 structure.
- **Trace.** `reference_pressure_kernels` (native-c2.py:509–525) builds T01 = h(l−m)·i f q_m and T12 = h(m−k)·i f q_k, so T1(l,k) = i f q(k)·h(Q), as in (3).
  - Back-substitution gives X02 = R02 − T01·R12 + T01·T12·R22. The last term is η², and `retained_shape` (native-c2.py:665–668) drops it.
  - The retained (1,1) trace term is therefore −T1·Fj.
- **Both faces.** The saved lower trace gives lab height −Wηw/2 and normal jet −i·qo, so their product is +i·q·h on both faces. `outgoing_spectral` uses the same root for both faces, and the lower slope join has no sign flip (lower-native-slope-input.json).
- **Final slot.** `build_face` (native-c2.py:589–591) sets reference pressure = physical pressure − (lab height)·(reference jet), with the jet taken from the reference kernels. That is T⁻¹ truncated at the retained grades. There is no second inversion.

**Algebra**
- (4)–(6) follow from R·Z. The reference Fh is −iμqi/(qo+β)·h, which is (5).
- Term 1 plus term 3 of (7) combine to C = −iμk·qo/[(qo+β)(qi+β)], as in (8). I verified this using β = aμ.
- The contact of term 1 vanishes on its own. The surviving contact C·(W/4)·j(Q) comes entirely from the trace term.
- (10): with qm−qi = −t(m+k)/(qm+qi) and t·h(t) = (W/2)A/i, the sign and prefactor aμ²WL/4 are correct. The height contact of term 2 is exactly zero at finite δ.
- The jet normalization j = (L/2)A is consistent. The physical profile is (1+tanh ξ)/2 with ξ = x/L, so the jet transform is L·A and σ_W = ηW/L (raw-worker.py:447).

**Constants and inequalities**
- Prefactor: |aμ²WL/4| ≤ 0.16·2.5 = 2/5.
- κ_δ: the minimum is at cs = 2, δ = 0.1, giving 879/400, so κ_δ ≥ √879/20.
- Lower bound on |q+β|: first-quadrant q and β give |q+β| ≥ |β| ≥ 0.2847 > b0 = 0.2703, valid at δ > 0 as well. Also |q| ≤ 4.25 ≤ 5 and |β| ≤ 0.3.
- 1/|q| bound: |q|² = |q²| ≥ |κ²−p²| = ab and a+b ≥ 2κ. This gives 1/|q| ≤ κ^(−1/2)(a^(−1/2) + b^(−1/2)), a sum rather than a product. The integral bound 4√(2u)/√κ is correct.
- Tail: |q| ≥ |t|/√2 is tight at |t| = 12 and holds beyond. (|t|+3)(|t|+6) ≤ 15t²/8 is also tight at 12. The A·A tail is actually about 37.5·e^{30π}, so 150 is loose but valid. The constant in (12) is about 159/b0³ ≤ 200/b0³.
- (9) pieces: the bounds |A'| ≤ L/4, L²/(16π) and L·e^(−πL|s|/2) all check. The inequality 0 ≤ x·cosh x − sinh x ≤ sinh²x is proved correctly. |C| ≤ 6/(5b0) is correct.
- (14): √7 from |2k+Q| ≤ 7, 2/b0 and 2√7/b0² are correct. The normal-jet variation constant 4√7/(5b0²) is correct.
- Normal-jet sup is ≤ 2, using |μ·qi| ≤ 0.4·4.25 < 2. The sup bound relies on |qo/(qo+β)| ≤ 1 and the jet variation on the β-difference identity; both hold.

**Limit logic**
- J: pointwise convergence off the endpoint set, plus uniform local integrability (the √u bound), plus the exponential tail, gives an L1 limit by Vitali. The qi/(qm+qi) factor is correctly kept whole.
- qi → 0 sends J to 0 in L1.
- C·H is continuous and has no endpoint singularity.
- The first height term is correctly treated as a distribution. The |Q|^(1/2) Hölder bound is integrable against A/Q even when grazing and a transfer zero coincide.

## Non-blocking points (scope and wording, not scientific blockers)

1. **The direct (1,1) kernel is not in native `second`.** `kernel_bridge` sets z_three[0,2] = 0 (native-c2.py:409), so native `second` is iteration only. The direct block enters through a tagged saved slot (`direct_slot_tag` in trace-difference-identity.json, and `D` in raw-worker.py:534). The final-slot join should therefore cover iteration plus trace only, with D added once as the tagged slot. The method's "added once, outside (2)" wording is consistent with this, but the grade census should say it explicitly.
2. **Some joins rest on source text, not saved operands.** The saved trace-domain records give the value coefficient 1, the height and the output-leg jet i f qo. They do not save the three-leg trace matrix, the middle/input normal multipliers, `reference_three`, or a generic final-slot expression. Those joins must come from reading native-c2.py:509–525 and 589–591 symbolically. The saved selected-point trace files do not cover arbitrary momenta. This is feasible without replay, but it is a derivation, not a restoration.
3. **The complex-frequency prescription is new.** `outgoing_spectral` is a real-frequency Piecewise using `roots[-1]` and sign(ω). The first-quadrant definition for δ > 0 is not in the native source. It must be verified through the unrestricted defining equations and then specialized to the saved ω = 3 coefficients. The method already says this.
4. **(13) is a fixed-k statement.** Uniformity in k, parameter derivatives and the untruncated solver are correctly disclaimed. Keep that in the result wording.
5. **Retention exclusions.** η², σ² and T1² are dropped by `retained_shape`. The instrument should record this as the retention rule, not as a regularity result.

## What this would establish if executed

- Exact source-joined formulas for F0, Fh, Fj and the reference mixed iteration, C·H plus ∫J, on both faces.
- Uniform L1 and tail certificates for J, and bounded continuity of C·H, on the compact domain at δ → 0⁺.
- A unique limiting distribution for the first-height term, and bounded normal-jet multipliers.

## What stays outside it

- The full slab operator, source-field and consumer composition, and arbitrary or lower-face Fourier source maps.
- The finite inverse, the 4/6 cutoffs, a defect sweep and physical loss.
- Differentiability in the speed, and the κ = 0 and β = 0 limits.

## Essential runtime evidence

- Symbol-level source joins for Zh, Zj, T1 and the final-slot routing on both faces, with each uncombined term of (7) kept.
- Multiplication identities Z0 = μ/q, (I+aZ)·R and T·Fref = Fphys on the unrestricted laws, followed by coefficient-wise specialization to the saved ω = 3 matrices.
- Exact certificates for the constants above, including qi = 0, qo = 0, both external matches and contact-endpoint coincidence.
- The four responsive controls the method lists, with the failed identities saved.

This is source review only. It is not worker or result clearance.