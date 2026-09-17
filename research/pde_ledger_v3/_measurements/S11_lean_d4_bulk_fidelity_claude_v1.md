# Fidelity review of D4C.1–D4C.4 (full D4 bulk contract): **CLEAR**

The Lean statements prove what D4C.1–D4C.4 claim, with the stated hypotheses, quantifiers and conventions. None of my findings requires a fidelity fix; the non-blocking ones are listed under Findings. This is one of the two independent reviews the contract requires; the other is still pending.

**Packet:** "S11 full D4 bulk D4C.1–D4C.4 fidelity contract v1", author Codex, aggregate SHA256 `061d0d96f0c3b62439875877aa4a107b53b6870a206b308ee8a04056bd8e8a11`. This is the value in `MANIFEST.json`; I did not recompute it.

## What I could not verify
I could only read files; I had no shell or Lean access. So:
- **Not run by me:** I did not run Lean, Lake, SymPy, Wolfram, the two instruments or any hash check.
- **Author evidence only:**
  - the builds, the 67 axiom lists and the 36 control outcomes;
  - the source, pin and object hashes, the 327 historical files and the 18 NP/T1/VC objects;
  - the native identity results and the D4 generator `--check`.
- **Proofs:** I read the tactic proofs but checked only their statements. That they compile rests on the recorded run2 logs, which I read and which are consistent (diagnostics, axiom output, check count).
- **Reviewed dependencies:** I read the imported definitions the new statements rely on, but did not re-review previously reviewed proofs (the D4 classification, D4B odd variation, S10 analytic lemmas).
- **What I did check by hand:** the basis matrix and its determinant, the native weights, the momenta, the EL operator, the modal and homogeneous maps, the currents and every witness value.
- **Other:** both Wolfram anchor strings exist (lines 863 and 910). There are no `sorry`, `admit`, `axiom` or `native_decide` in the new modules.

## Conventions checked
- **Coordinates and gradient:** `Point 4 = Fin 5 → ℝ` with time at index 0. `Even.spatialGradient J j i = J j.succ i` gives `G_ij = ∂_i u_j`, and the odd module's `spatialGradient` is definitionally the same.
  - The native helper uses the same convention (`coordinate_substitution`: `g[i,j] ↦ ∂_{x_{i+1}} u_{j+1}`, arguments `(x1..x4, t)`).
  - All four forms are unchanged under `G ↦ Gᵀ`, so a row/column slip could not change the density family in any case.
- **Invariants:** `divergence`, `transposePair` and `gradientSquare` are tr G, tr(G²) and tr(GGᵀ). `orientation` is `F01F23 − F02F13 + F03F12` with `F = G − Gᵀ`.
- **Rotations:** `SOInvariant` is invariance under `G ↦ RGRᵀ` for det R = 1, which is the correct transformation for `∂_i u_j` under a rotated field.
- **Coefficients:** they are fixed arguments `v : Vec 4`, unrestricted, and never varied. `Even.lagrangian` uses only `v 0, v 1, v 2`; the full `lagrangian` adds `S11D4Odd.lagrangian (v 3) = −(β/2)P`.

## Findings by item

### D4C.1 — density, momenta, variation
- **Density family** (`Action.lean`)
  - `density_identity` proves `L = −½·invariantForm v (G)`.
  - `all_invariant_densities` (via `SO_unique`) shows every SO(4)-invariant quadratic form is `invariantForm v` for exactly one `v`. Together these identify the entire classified family.
- **Momenta**
  - `momentum` is defined as `deriv` in the `basisJet j i` direction and is backed by a proved `HasDerivAt`.
  - `momentum_eq` splits it into even plus odd parts.
  - `Even.momentum_eq` gives `−(a δ_ri tr G + b G_ir + c G_ri)`, which matches ∂L/∂G_ri by hand.
  - `momentum_time` proves the time row is zero.
- **Relative action** (`Variation.lean`)
  - `relativeAction` is the Bochner integral of `L(jet(u+sh)) − L(jet u)`. It is not a formal stand-in for its own derivative, and no absolute action is ever defined, so none is assumed integrable.
  - Integrability is proved for the linear term and for the compactly supported quadratic term.
  - `relativeAction_expansion` gives the exact expansion, and `relativeAction_hasDerivAt` gives a genuine derivative.
  - Integration by parts uses `coord_integration_by_parts`. The fundamental lemma (`ae_eq_zero_of_integral_contDiff_smul_eq_zero`) plus continuity upgrades vanishing almost everywhere to vanishing at every point.
  - Result: `actionStationary_iff_eulerLagrange` (stationarity iff the local EL vanishes) holds for every smooth background.
- **Sign convention:** EL is `−Σ_j ∂_j p_ji`, the first-variation convention.

### D4C.2 — bulk operator
- **Operator identity** (`Bulk.lean`)
  - `Even.eulerLagrange_eq` gives `(a+b)·gradDiv + c·laplacian`, using `coordDeriv` and `partial_commute` (derivatives commute because the field is smooth).
  - `eulerLagrange_even` removes the odd part through the reviewed `S11D4Odd.eulerLagrange_zero`.
  - `gradDiv` and `laplacian` are spatial-only double derivatives, and `gradDiv_eq` ties `gradDiv` to ∂_i(div u).
- **Equivalence** (`Census.lean`)
  - `BulkEquivalent` quantifies over all smooth fields and all points.
  - Sufficiency in `bulkEquivalent_iff` comes from the full identity.
  - Necessity uses two waves, `planeWave 0 e₀ e₁` (transverse) and `planeWave 0 e₀ e₀` (longitudinal), evaluated at x = 0. I checked they give `−c` and `−(a+b+c)`.
  - These waves are only witnesses; they do not restrict the theorem's domain.
  - No sign or genericity hypotheses appear.

### D4C.3 — null family and currents
- **Null family**
  - `VariationallyNull` requires zero first variation for every smooth background and every test field.
  - `variationallyNull_iff` gives `c = 0 ∧ a+b = 0`, and `null_parameterization` gives `∃ t β, v = (t, −t, 0, β)`, with both parameters free.
- **Response map**
  - `responseMap` is a genuine `ℝ`-linear map `v ↦ (v 2, v 0 + v 1)`.
  - It is proved surjective, with image dimension 2 and kernel dimension 2 by rank–nullity, and `nullSpace_identification` shows its kernel is exactly the null family.
  - The theorem covers only the classified quadratic family, as documented.
- **Currents**
  - Even current: `Even.boundaryCurrent` is exactly `J_i = Σ_j(u_i ∂_j u_j − u_j ∂_j u_i)`. `boundary_identity` gives `(div u)² − tr(G²)`, which I confirmed by hand (the third-derivative terms cancel by commutation).
  - Odd current: `K_i = ½ Σ_j u_j M_ij`, with divergence `P`.
  - `null_density_is_divergence` writes the null density as `−t/2 · div J − β/2 · div K`. The weights are correct.
- **Nonzero witnesses** (all values checked by hand)
  - null density at `evenJet` = −1;
  - even momentum = −2;
  - odd density = −1;
  - odd momentum = −1;
  - `div J` = 2 on the smooth affine field `affineWitness = (x₁, x₂, 0, 0)`.
- **Homogeneous map:** `S11Homogeneous.modalOperator 0 c (a+b+c)` expands to `−c|k|² a − (a+b)(k·a) k`, which matches. `modal_zero_wavevector` covers k = 0. No spectrum claim is made.

### D4C.4 — controls and native check
- **Formal controls**
  - There are 16 pairs, each with a plausible wrong value, sign, factor or proposition.
  - Every recorded mutant has exit status 1, exactly one `contract_control` diagnostic reading `unsolved goals ⊢ False`, and no warnings or instrument errors.
  - Counts reconcile: 16 + 4 = 20 positive runs, of which 18 are distinct (the even and odd momentum positives each appear twice). That gives 36 executions; adding 26 dependency builds, 6 module builds and the audit root gives 69 check records and 33 objects.
- **Native instrument**
  - `gauss_jordan_solve` expresses the native basis rows in the trace/orientation basis. The rows are `(0,1,0,0)`, `(½,−½,0,0)`, `(0,−1,1,0)`, `(0,0,0,1)`; the determinant is −½; and the solved weights are `(a+b+c, 2a, c, β)`. I confirmed that the transpose of this matrix maps those weights back to `(a,b,c,β)`.
  - The helper EL is `+Σ ∂_i(∂L/∂g_ij)`, so native EL = −(Lean EL), and V5 of Q = 2·(Lean EL of L). Both relations are consistent.
  - The two source mutations replace the actual EL-assembly line of `euler_lagrange_from_placeholders`; the recorded report shows they break the EL and V5 identities.
  - The rref-normalized `PD_POLY` has leading monomial `g01·g23` with coefficient +1, which is consistent with `P_D = P`.
  - The seven target checks are real `assert`s (instrument lines 84–91). The JSON flags recording them are hard-coded but can only be written after those asserts pass.

## Non-blocking notes (optional, not required fixes)
1. **Odd current witness is jet-level only.** `odd_current_normalization` evaluates `currentJet` at a chosen value and jet, not on a smooth field. Nonzero boundary content still follows from `odd_density_normalization` together with the odd `boundary_identity`. An affine-field witness would match the even case if the "actual smooth affine-current witnesses" wording is kept.
2. **No definition-level mutations.** `MUTATIONS = []`, so no canonical definition was mutated. The canonical normalization theorems would stop compiling under such a mutation, but that was not demonstrated. A mutation like `−(1/2) → −1` in `Even.lagrangian` would be stronger evidence.
3. **Bulk-equivalence control covers one direction.** The pair only rejects wrongly declaring an equivalent pair inequivalent (a condition that is too strong). A "same c, different a+b" non-equivalence pair would test a condition that is too weak. The longitudinal and null-sign controls cover this indirectly.
4. **Homogeneous map is modal-level only.** A field-level identity with `S11Homogeneous.eulerLagrange 0 c (a+b+c)` would also hold (the two densities differ by an even null density). This is optional.
5. **Mixed-derivative canonicalization.** The native per-basis matches for the two null basis vectors depend on SymPy canonicalizing mixed derivatives. The documentation acknowledges this, and the momentum and odd-density checks cover the normalization that those zero EL expressions cannot reveal.
