# Verdict: CLEAR (no blockers)

I found no fidelity blocker. The Lean statements match the object you specified, the proofs follow the mathematics I recomputed by hand, and the native link is consistent with the helper source I read. Everything below comes from reading files; I compiled nothing and ran nothing (see the limits section at the end).

**Reviewed packet:** `S11 D3 bulk variation K1–K4 fidelity contract v1`, aggregate SHA-256 `bb5fe9bbfbc09f3a184bb6f256daba74f9ee7d6baa85cbb67b6e4eb2ddb83b54`. This is the value as supplied; I did not recompute it.

## Definitions, assumptions and quantifiers checked

**Density** (`S11D3Bulk/Action.lean`)
- `spatialGradient J j i = J j.succ i`, i.e. `G_ij = ∂_i u_j`.
- The trace, `tr(G²)` and `|G|²` terms are written in index form, and `lagrangian v J = -(1/2)(a·div² + b·tr G² + c·|G|²)`.
- `density_identity` *proves* this equals `-(1/2)·S11D3Invariants.invariantForm v (spatialGradient J)`; it is not a definitional equality. I checked `invariantForm_apply` (`Rotation.lean:103`) and the three trace-form lemmas.
- `all_invariant_densities` reuses `SO_unique`.
- `v : Fin 3 → ℝ` is unrestricted: no sign, nonzero or symmetry assumption, and no inertia term (the time-component momentum is 0).

**Momenta and variation**
- `momentum` is defined as the actual derivative `deriv (s ↦ L(J + s·basisJet j i)) 0`.
- By hand, `momentum_eq` gives `P_ri = -(a δ_ri div + b ∂_i u_r + c ∂_r u_i)`. I also rederived `variationDensity` and the exact quadratic expansion `lagrangian_increment`.
- `relativeAction` integrates only the density *change* under `u + s·h`, so the background's total action need not be finite.
- `relativeAction_hasDerivAt` differentiates that actual integral, via `relativeAction_expansion`.
- `integrated_variation_by_parts` uses S10's `coord_integration_by_parts` for compactly supported `g`. The sign agrees with `eulerLagrange = -Σ_j ∂_j P_ji`.
- `ActionStationary` is phrased with `deriv`, but it is not vacuous: `relativeAction_hasDerivAt` proves the derivative exists for every test field.

**Stationarity ⇔ pointwise PDE** (`actionStationary_iff_eulerLagrange`)
- Single-component test fields (`single_testField`) plus `ae_eq_zero_of_integral_contDiff_smul_eq_zero` give vanishing almost everywhere.
- Continuity and `Measure.eq_of_ae_eq` then give vanishing at every point.
- The quantifiers are: for all smooth `u`, for all smooth compactly supported `h`, if and only if for all `x`.

**Name resolution**
- S11D3Bulk reuses S10Pilot names while `open S10Pilot` is in effect. Lean 4 resolves names in the current namespace first, and the argument types differ anyway (`Vec 3` vs `ℝ`), so nothing silently binds to the S10 version.
- The control files use fully qualified names where it matters.

## Findings on the questions you asked

1. **EL = (a+b)∇(div u) + cΔu on every smooth field.** Yes: `eulerLagrange_eq` (`Bulk.lean:26`), for all smooth `u` and all `x`. By hand, `-Σ_r ∂_r P_ri = a∂_i div + b Σ_r ∂_r∂_i u_r + cΔu_i`, and the middle term needs `∂_r∂_i = ∂_i∂_r`.

2. **Mixed derivatives are commuted from smoothness, not assumed.** `partial_commute` (`Calculus.lean:64`) follows from `ContDiff ℝ ∞` through `partial_second` (which identifies the iterated partial with `D²f x (e_i)(e_j)`) and mathlib's `isSymmSndFDerivAt`. There is no extra hypothesis.

3. **The current's divergence.** `boundaryCurrent` is exactly `F_i = Σ_j (u_i ∂_j u_j − u_j ∂_j u_i)`, and `boundary_identity` proves `div F = (div u)² − tr(G²)`.
   - By hand: the terms `Σ u_i ∂_i∂_j u_j` and `Σ u_j ∂_i∂_j u_i` cancel after relabelling and commuting.
   - `null_density_is_divergence` gives `L(a,−a,0) = −(a/2)·div F`. This is a divergence, not a pointwise-zero density, and no boundary physics is discarded.

4. **Bulk equivalence holds exactly when (c, a+b) agree.** `bulkEquivalent_iff` quantifies over all smooth `u` and all `x`.
   - Necessity uses static plane waves with `k = e₁`: a transverse amplitude `e₂` separates `c`, and a longitudinal amplitude `e₁` separates `a+b+c`. Both are smooth by `smooth_planeWave` and are evaluated at `x = 0`, where cos = 1.
   - Sufficiency uses the full formula `eulerLagrange_eq`, not sampling.

5. **Variational nullness is exactly c = 0 and a+b = 0.** `variationallyNull_iff` goes through `BulkEquivalent v 0`, and `null_parameterization` gives `v = (t, −t, 0)`.
   - `v` is arbitrary, so zero and negative coefficients are included. The controls exercise `0`, `(1,−1,0)` and `(−2,3,−5)`.

6. **The kernel/range dimensions belong to the actual operator.** `nullSpace_identification` proves `ker responseMap ↔ VariationallyNull`, and `bulkEquivalent_response` proves `BulkEquivalent ↔` equal `responseMap` values. So the dimensions 1 and 2 (`null_dimension`, `bulk_response_dimension`, `response_surjective`) are tied to the variational operator, not to an unrelated coefficient map.

7. **Nonzero first variation outside the null locus.** `exists_nonzero_firstVariation` proves it by contradiction (non-constructively). This is sound: the contrapositive of the stationarity equivalence combined with the plane-wave separation. `null_nonzero_example` and the positive controls show the hypotheses can actually be met.

8. **Modal operator (K3).**
   - `momentum_contraction` evaluates the actual momenta on `modeJet`.
   - `eulerLagrange_planeWave` proves `EL(planeWave) = cos(phase)·M·a` for all `ω` and all `k`, including 0.
   - `M a = −c|k|²a − (a+b)(k·a)k`, which agrees with `∇div → −kkᵀ` and `Δ → −|k|²`.
   - `homogeneous_operator` identifies `M` with `S11Homogeneous.modalOperator 0 c (a+b+c) ω k a` for all `k`, since `μ − B = −(a+b)`.
   - `modal_longitudinal` gives stiffness `a+b+c`, `modal_transverse` gives `c`, and `modal_zero_wavevector` covers `k = 0`.
   - No spectrum claim is made.

9. **Native link** (`S11_lean_d3_bulk_source_check.py` and the helper source, lines 316–378 and 478–600).
   - **Placeholder convention:** `G_i_j ↦ Derivative(u_j, x_i)`, the same as Lean. The three forms are transpose-invariant, so ordering could not change the result anyway.
   - **Sign:** the native EL computes `+Σ_i ∂_i(∂L/∂G_ij)`, so native EL(L) = −Lean EL(L), and V5 applied to the unnormalised Q gives 2·Lean EL(L). This is confirmed by `actual_L_native_EL_equals_negative_Lean_EL` and `V5_all_basis_combination`.
   - **Basis change:** I checked `[[0,1,0],[½,−½,0],[0,−1,1]]` by hand; its transpose applied to `(a+b+c, 2a, c)` returns `(a, b, c)`. `basis_map·target = native_basis` and a nonzero determinant are asserted.
   - **Normalisation** is checked against the native `V1_POLYS`, independently of the `V1_BASIS` vectors.
   - **Every V5 basis vector:** I compared the recorded outputs to `2·Lean EL(row)`. The null basis element `½(tr G)² − ½tr G²` gives terms `∂₁∂₂u₂ − ∂₂∂₁u₂`, which vanish only after SymPy's `.doit()` canonicalises the derivative order. This is exactly the documented symbolic reduction, and the Lean commutation proof justifies it for smooth fields.
   - **Controls:** the native sign and factor mutations change line 377 inside `euler_lagrange_from_placeholders` and break V5, the sign identity and the modal identity. The wrong-B map and the null/non-null target controls also pass. Wolfram material was only inspected for source anchors, and no kernel certification is claimed.

10. **Mutations** (the diagnostics shown in `S11_lean_d3_bulk_contract_checks.json`).
   - **`density_sign` / `density_half`:** unsolved polynomial goals in `density_identity`. `lagrangian_increment` also fails, which is expected. The recorded counterexample values (+½ and −1 instead of −½) are correct.
   - **`bulk_coefficient_map`:** the goal `v1 − v0 = −v1 + −v0 ∨ …` is the false coefficient identity itself, and the `+e₁` vs `−e₁` example is correct.
   - **`boundary_current_sign`:** the log has the unused-`partial_sub` linter error (made fatal by `-DwarningAsError`) *and* an unsolved-goals error in `boundary_identity` showing a residual derivative identity. I confirmed independently that the mutated identity is false: for `u = (x₁,0,0)` the mutated current has divergence 2, while `(div u)² − tr G² = 1 − 1 = 0`. Rejecting it on mathematical grounds is justified.
   - **Eight statement mutants:** each reduces to `⊢ False` after rewriting with the proved theorem, and each true counterpart passes.
   - **Instrument failures:** no import, syntax, timeout or resource failure was counted as a rejection; the classifier's `invalid` pattern excludes these.

11. **Accounting.**
   - The checks add up: 44 logs (15 dependencies + 5 builds + 1 audit + 4 source mutations + 16 paired examples + 3 positives), 21 compiled objects, and 49 `#print axioms` roots (`S11D3Bulk.lean:4–52`).
   - The source audit found no `sorry`, `admit` or custom axioms.
   - The report, instrument and source hashes in `S11_lean_d3_bulk_validation.json` match `MANIFEST.json`.

## Optional improvements (not blockers)

- **O1.** K3 is identified with the homogeneous operator only at the plane-wave level. A full-field identity with `S11Homogeneous.eulerLagrange 0 c (a+b+c)` would also hold, because the two densities differ by a multiple of the null combination `(b+c)`, but it is not stated.
- **O2.** The mutation classifier accepts any "unsolved goals" diagnostic in the named declaration, so it depends on a human reading the goal (I did). Matching the expected residual more tightly would make it self-checking.
- **O3.** The `zero_wavevector` wrapper control only checks one component. The nonzero-first-variation positive is non-constructive; an explicit witness (a plane wave times a bump function) would be stronger evidence of non-vacuity.
- **O4.** The dependency builds were reused from a run with a different instrument hash (`6ce57a…`), with guards on sources, command and objects. This is acceptable as documented, but worth noting.
- **O5.** Units and dimensional consistency are not modelled. This is already documented.

## Limits of my independent verification

- **Not run:** I did not compile any Lean, run the Python or SymPy checks, execute Wolfram, or recompute any hash (the aggregate or per-file); I had no shell tool in this session.
- **Build and test results are the packet's own records:** compilation, the axiom audit and the mutation outcomes all come from the recorded JSON and logs.
- **What I did check:**
  - I read all five S11D3Bulk modules, the audit root and the relevant imported definitions (S10Pilot Action/PlaneWave/Analytic/Variation, S11D3Invariants Quadratic/Rotation plus the `SO_unique` signature, and S11Homogeneous Action).
  - I recomputed the key identities and counterexamples by hand.
  - I read the native helper source and both instruments.
  - I checked that the recorded hashes are consistent with each other.
- **Not reviewed:** S10Controls internals and the D3 classification proof beyond the lemmas used here. As you asked, I did not re-review the earlier H/I/E/J scope.
- **Unchanged:** I modified no files.
