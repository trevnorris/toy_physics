**Verdict: CLEAR**

I found nothing that needs fixing for D5B.1–D5B.4. The statements match the stated contract, and the documents don't claim more than the evidence supports. I have three optional notes, listed at the end.

`_measurements/S11_lean_d5_bulk_review_packet.json` isn't in the packet, so I treated `MANIFEST.json` as the authoritative manifest.

## What this review rests on
- **Source reading only.** I read the policy, coverage, fidelity and verification documents; the five `S11D5Bulk` modules and the root; the S10Pilot Action, PlaneWave, Analytic and Variation modules; `S11D4Odd/Calculus.lean`; `S11Homogeneous/Action.lean`; the matching `S10Controls` definitions; and the Classification, Forms, Rotation and CoordinateAlgebra declarations it imports.
- **Author records I checked but did not reproduce.** I compared the native instrument and its source against its report, and read the contract instrument, the recorded run4 results and the author validation.
- **What I did not do.** I executed nothing and computed no hashes. The only hash check was that the manifest, run4 source-dependency and validation hashes agree with each other.
- **Build status.** Build, axiom and control outcomes come from the author's records; I did not reproduce them. The 42 `seed_S11D5Invariants.*` objects are validated reuse of the portable D5 replay, not fresh builds. The 18 other objects were built in bulk runs 1–3 and reused under guards in run4 (the records carry `reused_from_instrument_sha256 f1927b…`). The 54 axiom lists come from the recorded run3 root log. Only the 32 controls ran fresh in run4.

## D5B.1 — Action, momentum, first variation
- **Three coefficients vs five components.** `Coeff = Vec 3` holds (a,b,c). Fields are `Point 5 → Vec 5`, with `Point 5 = Fin 6 → ℝ`, where row 0 is time. The code never swaps coefficient-space dimension for spatial dimension.
- **Gradient convention.** `spatialGradient J j i = J j.succ i` (`Action.lean:14`), so G_ij = ∂_i u_j with the derivative index as the row.
- **Lagrangian.** `lagrangian` at `Action.lean:19` is exactly −½[a(div)² + b ΣJ_{i,j}J_{j,i} + c ΣJ_{i,j}²]. `density_identity` ties it to `invariantForm_apply` (a·tr²G + b·tr(G²) + c·tr(GGᵀ)). Uniqueness comes from `SO_unique`, and `O_classification` exists in the reused module.
- **Momentum.** It is a genuine `deriv` along `basisJet` (`Action.lean:57`). By hand I get p_{r,i} = −(a·δ_{ri}·div + b·∂_i u_r + c·∂_r u_i), which matches `momentum_eq`. The time row is zero (`momentum_time`).
- **Integral.** `relativeAction` integrates the pointwise density difference (`Variation.lean:53`). Integrability is proved from compact support of the test field (`integrable_linearDensity`, `integrable_test_lagrangian`), so the background action is never assumed integrable.
- **Differentiation and the fundamental lemma.** The derivative comes from an exact quadratic expansion (`relativeAction_expansion`, then `HasDerivAt`). Integration by parts is Mathlib's `integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable`. The fundamental lemma is `ae_eq_zero_of_integral_contDiff_smul_eq_zero`, upgraded to pointwise equality by continuity.
- **Domain.** `SmoothField` is componentwise `ContDiff ℝ ∞` on ℝ⁶; `TestField` adds compact support.

## D5B.2 — Local equation and census
- **The local equation.** `eulerLagrange = −Σ_j ∂_j p_{j,i}`. `eulerLagrange_eq` gives (a+b)·grad div + c·Δu for every smooth u. Mixed partials commute via `partial_commute`, which comes from Mathlib's `isSymmSndFDerivAt`. By hand, −∂_j p_{ji} = (a+b)∂_i div u + cΔu_i, which is the correct sign for L = −Q/2.
- **Bulk equivalence.** `BulkEquivalent` quantifies over all smooth fields and all points. `bulkEquivalent_iff` states: equal iff c agrees and a+b agrees. Sufficiency uses the full formula; necessity uses a static transverse wave and a static longitudinal wave (k = e₁, a = e₂ and a = e₁). The waves are used only for necessity, as required.
- **Null family and dimensions.** `null_parameterization` gives exactly (t,−t,0). `nullSpace_identification` ties `ker responseMap` to `VariationallyNull`. Image and kernel dimensions are 2 and 1 via `range_eq_top` and rank–nullity on ℝ³.
- **Coefficients.** Every theorem is universal in v, so zero and negative coefficients are included; the `negative_coefficients_positive` control exercises this.

## D5B.3 — Current, normalization, homogeneous map
- **Current.** `boundaryCurrent` is exactly J_i = Σ_j(u_i ∂_j u_j − u_j ∂_j u_i). `boundary_identity` gives div J = (div u)² − tr(G²), and `null_density_is_divergence` gives the factor −t/2.
- **Nonzero witnesses.**
  - `evenJet` (∂₁u₁ = ∂₂u₂ = 1) gives density −1 for (1,−1,0) (`even_null_density_nonzero`).
  - It gives momentum −2 for (1,0,0) (`even_momentum_normalization`).
  - `affineWitness = (x₁, x₂, 0, 0, 0)` is proved smooth, and its jet is proved equal to `evenJet`. It gives div J = 2 (`even_current_normalization`), which is consistent with −(1/2)·2 = −1.
- **Homogeneous map.** `homogeneous_operator` holds for every ω. It gives modal operator −c|k|²a − (a+b)(k·a)k = `S11Homogeneous.modalOperator 0 c (a+b+c)`, since (μ − B) = −(a+b). `modal_zero_wavevector` covers k = 0. `eulerLagrange_planeWave` connects the modal operator to the actual local equation.
- **Boundary effects.** No declaration claims pointwise-zero density or absent boundary effects, and the documents say so explicitly.

## D5B.4 — Native identification
- **Index convention.** Native `coordinate_substitution` maps `g[i,j]` to ∂_{x_i} u_j, the same row-derivative convention as Lean. The recorded native momenta (for example, entry (1,2) = −c·G₁₂ − b·G₂₁) match `momentum_eq`.
- **Full span.** `basis_map·target = V1_BASIS` is checked with a nonzero determinant. The rows are (0,1,0), (½,−½,0), (0,−1,1). I computed by hand: det = −½ and weights (Bᵀ)⁻¹(a,b,c) = (a+b+c, 2a, c). Both match the report.
- **Sign and factor.** `euler_lagrange_from_placeholders` computes +Σ_i ∂_i(∂L/∂g_ij), which is opposite in sign to Lean's convention. The report's `actual_L_native_EL_equals_negative_Lean_EL` and `V5 = 2·Lean EL` checks pass on each basis element.
- **Native mutations.** The sign and factor mutations of the original helper each flip the V5, EL and modal identities to false. Unrelated identities stay true, which shows each mutation hits its target.
- **Scope of the native check.** It is a tested translation outside Lean's kernel. The Wolfram file was only searched for two text anchors, not run.

## Controls and mutation coverage
- **Genuine paired controls.** Four pairs discriminate on real mathematical content:
  - `null_sum_sign`: (1,1,0) is not null.
  - `omit_laplacian`: (0,0,1) is not null.
  - `bulk_equivalence`: (2,3,7) and (8,−3,7) are equivalent, since a+b agrees while a and b differ.
  - `time_momentum`.
- **Value-pinning controls.** The other pairs fix a proved constant against a wrong one: dimensions 2/1, momentum sign and factor, current factor, and the transverse/longitudinal responses in both the e₁ and fifth (e₅) directions. They show the proved values are specific; they are not source-replacement mutants, and the documents say so.
- **Normalization despite null variation.** The density −1, momentum −2 and current 2 witnesses catch normalization errors that the zero bulk response of the null family cannot.
- **Rejection reasons.** The recorded mutant outputs show a single `unsolved goals ⊢ False` in `contract_control`, and failed runs are never counted as mathematical rejections.

## Optional notes (not required fixes)
1. **Hollow native control.** `current_factor_rejected` (`S11_lean_d5_bulk_source_check.py:84`) compares 2X against X for the same expression. It never touches a computed object, so it adds nothing beyond `current_divergence` and the Lean `even_current_factor` pair. Consider relabelling it, since it is currently counted among the "three deliberately wrong formulas".
2. **Operator-level map.** The ρ=0, μ=c, B=a+b+c map identifies the operator, not the density. The D5 density minus the S11Homogeneous density equals the null term with t = −(b+c). The documents already call it an operator/modal map; saying so explicitly would prevent it being misread as an action identity.
3. **Non-constructive witness.** `exists_nonzero_firstVariation` is proved by contradiction and does not construct an explicit field. That is acceptable for the existence claim as worded.

Scope exclusions were respected: no variable coefficients, interfaces, general-dimension theory, spectrum work, S11c or production execution. I did not contact the other reviewer and did not commit anything.
