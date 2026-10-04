# Verdict: CLEAR FOR THIS FIRST-ORDER MATCHED RECEIVING AND FACE METHOD, with five binding stop-gates and partial coverage

I found no blocker in the method as written. The five gates are places where the plan's own "candidate" language has to be resolved from the saved operands, never fitted.

## What I checked and found sound

- **Transverse chart.** The columns I checked are consistent with the saved data.
  - `k(l)·t = 0` and `k(l)·b(l) = 0` for every `l`, so `E(l)` spans the transverse subspace exactly.
  - `t·b = 0`, `b·b = H2·K2` and `t·t = H2`, so the stated dual rows are exact.
  - From `incident-columns.json`, `U_B = i·t` and `U_A = i(4·b(p) − 2p·t)`. So `C = i[[−2p, 1], [4, 0]]` with `det C = 4`, and `B(p) = U` holds.
  - `k(l)·U = (l−p)·U_x` is exactly right, and `U_B` has `U_x = 0`.
- **Step decomposition.** In `full-local-force-column-1.json` the force is `U_B·(2T²+4.5T+2.5)`. That is 0 at `T=−1` and 9 at `T=+1`, which matches `A(−∞)=0, A(+∞)=9` and the saved profile limits (W: 0→1, M: 0→0).
- **Distributions.** `Q·Â = −i·(A′)^` is correct for the forward convention. The delta term is annihilated by `Q`. The half-line Abel Green integrals avoid the `δ(l−p) × singular denominator` product.
- **Radiation kernel.** `i·e^{ip|x−y|}/(3p)` is the right outgoing kernel for `−(3/2)(∂²+p²)`, and `l²−p²−i0` matches `p→p+i0`. The finite transmission correction is retained through the `B′`, `D_B′` terms and the independent scalar step solution.
- **K1.** `K1 = G1 + T1†G0 + G0T1` follows from `a†(I+λT1)†G_R(I+λT1)a`.
  - It is Hermitian by construction. The plan is right that this property cannot be used to test `T1` alone.
  - Reflection enters survival at `λ²` only.
  - The saved current Gram (`[[32.88, −1.34], [−1.34, 0.274]]`) is far from diagonal, so keeping arbitrary `a` is necessary.
- **Flux definition.** `SLAB_CURRENT_MATRIX` depends on `mu_R`, `mu_S`, `G_theta_u`, the kappas, `W_0` and both momentum legs. So `G1` is a real material-plus-momentum derivative and cannot be a coordinate Gram. The mass-rate matrix couples only through θ/e_W, which vanish on the incident doublet. So it matters only if the end eigenvector acquires O(λ) θ/e_W parts, as the plan says.
- **Face maps.** The finite flat factor is correct: `3/(10·(30+9i)/109) = 1 − 3i/10`. `beta ≠ 0`, `V_mem = X_mem = 0`, and pressure stays in the load. The plan also keeps the no-double-count rule for `R00·S01`.
- **Scope discipline.** The plan treats the RIGHT native current as opaque evidence and does not claim any result for it.

## Binding stop-gates

These are not reasons to withhold the verdict.

1. **Sign of `δp`.**
   - `δp = +9/D_T′(p) = +3/p` holds only if the saved force sits on the right-hand side.
   - If `f10 = P1·ψ0` is an operator perturbation, then `D_T(p_R) + 9λ = 0` and `δp = −3/p`.
   - Fix this from the actual end-pencil sign convention. The scalar step control must reproduce it. Never fit it to a finite1% root.
   - Check that the end force is exactly `9·I` on the doublet in the end pencil. The `U_B` data show `9U_B`. The `U_A` end value should be checked the same way.
2. **Transverse–face coupling at off-wave `l`.**
   - `b(l)` has normal component `H2`, and the saved pencil has `1/depth` entries that cancel only on the selected doublet.
   - Coupling that is O(q) in each direction gives a Schur-complement term ∝ `q² = p²−l²`. Its slope at `l=p` is `−2p`, which is comparable to `D_T′ = 3p`.
   - That would change `δp`, the pole residue, and the flux normalization. The plan says "test both directions" but gives no branch for nonzero coupling.
   - Add a rule: nonzero coupling means switching to `D_eff = D_T − C·M⁻¹·C′` and re-deriving `δp`, the residues and `G`. Do not continue with scalar `D_T`.
3. **Threshold coincidence.**
   - `q = 0` at `l = ±p`, so the exterior branch point sits on the transverse pole.
   - "Other channels" then carry `√Q` terms and algebraic tails, not exponential decay.
   - State the asymptotic projection that defines `T1`, `R1` and the flux. State which sheet the 3×3 determinant uses at `q = 0`. State that the O(λ²) pole shift is sheet-dependent.
4. **Gauge and origin invariance of K1.**
   - `K1` must be unchanged under `B′ → B′ + B·M` (with `T1`, `G1` shifting consistently) and under a phase-origin shift `x0`, which changes `T1` by `iδp·x0·I`.
   - These two controls should run through the same equations as the real calculation.
5. **Orientation and sign of `G_ref`.**
   - The saved orientation control in `incident-columns.json` records `disagreement: true` (assigned −1, reversed +1, independent group sign −1).
   - Resolve the sign of `−B(−p)†J_L(−p,−p)B(−p)` from the actual saved convention before using it. Do not rely on a label.

## Physical inputs, as distinct from derivations

- **Needed from saved operands:**
  - the force-versus-RHS sign;
  - the end-material parameter shifts entering `G1` (`mu_R`, `mu_S`, `W_0`, …);
  - the RIGHT current, restored and joined to source and grade.
- **Still open and correctly flagged by the plan:**
  - the O(λ) θ/e_W content of the end eigenvector;
  - full receiving regularity;
  - the held-profile and external-work balance. A nonzero `K1` could be dissipation or work, and that cannot be told apart without that balance.

## Coverage

- **Read in full:** `plan.txt`, `guide.txt`, `incident-columns.json`, `receiving-sheet.json`, the RIGHT opaque index and receipts, `slab-slab_current_matrix.json`, `slab-mass_rate_correction_matrix.json`.
- **Read in part:**
  - `first-order-pressure-assembly.json`: only lines 1–945 of 4398, and the `U0` rows I saw are all zero.
  - `full-local-force-column-1.json`: partial.
  - `RIGHT-raw-source-binding.json`: partial.
  - `native-views/`: grepped; the `acoustic-chemical-driver` and face-record bodies were not read.
- **Not read:**
  - `full-local-force-column-0.json`;
  - the bodies of `LEFT-invariant-P.json` and `acoustic-face-records.json`;
  - the memory-kernel details and closure joins;
  - the 781 source-stage artifacts;
  - the pickles, which were not decoded and cannot be checked from hashes.
- **Not verified by me:** that `D_T` is scalar off-wave, any determinant, and any face-map value. All of these are derivations the plan assigns to a later step.

This verdict does not authorize an unspecified worker. It covers only the physics method above.