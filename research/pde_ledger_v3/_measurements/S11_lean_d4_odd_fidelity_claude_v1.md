**Verdict: CLEAR** for D4B.1–D4B.4, within the stated scope. I found no required corrections. Four optional improvements are listed under "Optional improvements".

This review is based on reading the files. I did not build Lean, run either instrument or any mutation, or recompute any hash, and I modified nothing.

## Reviewed packet
- **Packet:** `MANIFEST.json`, revision "S11 D4 odd divergence/variation D4B.1–D4B.4 fidelity contract v1", `aggregate_sha256 = 1c2fb3840e56f82e506adbddb13463a6fbc645193294a3173531026f0ae838d5`. I used this identifier as supplied and did not recompute it.
- **Hash consistency (text comparison only):** the hashes written in the manifest agree with those recorded in `S11_lean_d4_odd_validation.json` and the two check reports. This covers the instruments, both reports, the native sources, the generator, and the proof and dependency sources.
- **Recorded counts:** 19 compiled objects (14 + 4 + root) and 43 check logs (14 + 4 + 1 + 4 + 16 + 4) are internally consistent.

## Conventions checked

- **Gradient, row = derivative index.** `Jet 4 = Fin 5 → Vec 4`, and `fieldJet u x j i = coordDeriv j (u · i) x`. `spatialGradient J i j = J i.succ j`, so G_ij = ∂_i u_j with the row as the derivative index.
  - `antisym J i j = J i.succ j − J j.succ i`, so antisym is F_ij.
  - `orientationJet_eq` holds by `rfl` and gives P = F01F23 − F02F13 + F03F12.
  - I checked the monomial indices in `orientationForm` by hand; `orientationForm_apply` also proves them.
- **Action and specification.** Spec §2 gives L = T − W, and §7 gives W_XFORM_EXTRA ⊃ (β/2)P_D. The isolated odd term is therefore −(β/2)P_D, which matches `lagrangian beta J := -(beta/2) * orientationJet J`. Beta is an unrestricted `ℝ`, with no positivity or nonzero premise anywhere.
- **Native side.**
  - `derivative_placeholders` and `coordinate_substitution` map `G_{i+1}_{j+1}` to `Derivative(u_{j+1}, x_{i+1})` and `V_j` to `∂_t u_j`, so the row is the derivative index there too.
  - Derivatives are taken with respect to named symbols, so the argument order `(x1..x4, t)` does not matter. Lean coordinate a+1 corresponds to native `x_{a+1}`. The mutation "∂₁" (`coordDeriv (1 : Fin 5)`) and the native `wrong_index` (`xs[0]`) both mean x1, so they agree.
  - The instrument rebuilds `QG_ALL`, `X_ALL` and `t` itself. These are identical to the native globals: `declared_symbol` returns `sp.Symbol(name, real=True)`.
- **P_D.** `compute_q9` sets `PD_POLY` to the sum of the actual rows of the rref'd `V6_BASIS` (`V6_DIM = 1`), as §7 requires, with no rescaling. The instrument asserts `actualP − P = 0` exactly, which fixes both sign and scale. The recorded P string matches the expansion.
- **EL sign.**
  - Native `euler_lagrange_from_placeholders` computes Σ_i ∂_i(∂L/∂G_ij) + ∂_t(∂L/∂V_j) (quadratic L, no u-dependence), which is §Q1's sign.
  - Lean's `eulerLagrange` is −Σ_j ∂_j(∂L/∂J_ji), the standard ∂L/∂u − ∂·∂L/∂(∂u).
  - The two therefore have opposite signs, as documented.
- **Excluded premises.** No symmetry of G, transverse ansatz or time independence enters: the only hypotheses are `SmoothField` and `TestField`.

## Findings by item

### D4B.1: density, classification, momenta, first variation
- **`all_odd_densities`** (`Action.lean:33`) is exactly `S11D4Invariants.odd_classification`.
  - That theorem is exhaustive over all quadratic forms on real 4×4 matrices. It goes through `SO_classification`, whose proof starts from `quadratic_representation` (every quadratic form's coefficient expansion), and then `invariantForm_odd`.
  - The group action is G ↦ RGRᵀ, the correct transformation for G_ij, and `reflection = diag(−1,1,1,1)` has determinant −1.
- **`density_identity`** equates `lagrangian beta` with −½·`invariantForm ![0,0,0,beta]` on `spatialGradient J`, which is the same family.
  - *Scope note:* the classification covers quadratic forms in G only. Velocity terms (V·G, V·V) are outside it. This is consistent with the contract's wording.
- **Momentum** is `momentum beta J j i := deriv (s ↦ L(J + s·basisJet j i)) 0`, an actual derivative. `basisJet j i` has entry 1 exactly at (j, i).
  - `lagrangian_variation` establishes the `HasDerivAt`, so the `deriv` value is not a junk value.
  - `momentum_eq` gives 0 for j = 0 (the time row) and −(β/2)·`dualCurl J r i` for j = r+1.
- **I recomputed ∂P/∂G_ri for all 16 entries** and they match M:
  - ∂P/∂F01 = F23, ∂P/∂F02 = −F13, ∂P/∂F03 = F12, ∂P/∂F12 = F03, ∂P/∂F13 = −F02, ∂P/∂F23 = F01.
  - Each is antisymmetrized through F_ab = G_ab − G_ba; the diagonal is zero.
- **`linearDensity`** is Σ_{i,j} momentum_ji · H_ji, built from the derivative-defined momenta. It is not a separately posited operator.
  - `linearDensity_eq` proves it equal to `variationDensity`.
  - `variationDensity` is the exact linear coefficient in `lagrangian_increment`: L(J+sH) = L(J) + s·B(J,H) + s²·L(H). I checked the bilinear expansion by hand.
- **Witnesses, checked by hand:**
  - `nonzero_density`: G01 = G23 = 1 gives P = 1.
  - `nonzero_momentum`: at β = 2, j = 1, i = 1, the value is −(2/2)·M01 = −F23 = −1.
  - `witness_lagrangian`: −β/2.
  - `current_normalization`: ½·2·M01 = 1.

### D4B.2: divergence
- **`dualCurl_contraction`.** Σ_{i<j} F_ij M_ij = 2(F01F23 − F02F13 + F03F12), which is 2P. This matches Euler's degree-2 identity.
- **`dualCurl_divergence_zero`** (Σ_i ∂_{i+1} M_ij = 0) is the Bianchi-type identity. I checked the j = 0 case by hand: −∂1F23 + ∂2F13 − ∂3F12 = 0 once second partials commute.
- **Mixed-partial regularity.** Commutation uses only `partial_commute`, which requires `ContDiff ℝ ∞ f` and relies on Mathlib's `isSymmSndFDerivAt` (C² suffices). It is applied to `hu r`, so the regularity is appropriate. `∞` here means C^∞, not analytic.
- **`boundary_identity`.** Using `partial_mul`, `partial_sum` and `partial_const_mul` (each with smoothness side conditions), the divergence is ∂_i(u_j M_ij) = G_ij M_ij + u_j ∂_i M_ij. So div K = ½(2P + 0) = P, and the factor ½ is forced.
- **`lagrangian_is_divergence`** gives L = −(β/2) div K pointwise. No claim is made that boundary effects are absent.

### D4B.3: bulk conclusion
- **`eulerLagrange_zero`.**
  - The j = 0 term is ∂₀ of the constant 0.
  - For j = r+1 the derivative index matches the row index of M, the same pattern as `dualCurl_divergence_zero`.
  - It holds for every β and every smooth u, with arbitrary time dependence.
- **`relativeAction`** is ∫[L(jet(u+sh)) − L(jet u)]; only the difference is integrated.
  - Integrability comes from `lagrangian_change`: the integrand equals s·linearDensity(u,h) + s²·L(h).
  - The first part is integrable as continuous × compactly supported.
  - L(h) = ½·linearDensity(h,h) (`linearDensity_self`), so it is integrable as well.
  - No integrability of the background action is required, which avoids the Bochner-integral junk-value trap.
- **`relativeAction_hasDerivAt`** differentiates the actual integral through the proved quadratic expansion.
- **`integrated_variation_by_parts`** applies S10's `coord_integration_by_parts` (Mathlib `integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable`, compact support of h) term by term.
- **`firstVariation_zero` and `every_background_stationary`** quantify over all smooth u, all smooth compactly supported h and all real β. Positive controls instantiate β = −7 and β = 0.
- **Name resolution.** S10Pilot also defines `lagrangian`, `momentum`, `eulerLagrange` and `smooth_momentum`. Lean 4 resolves names declared in the current namespace (`S11D4Odd`) before opened ones, and the argument types differ in any case, so the S11 definitions are the ones used.

### D4B.4: compact native link
- **Helpers.** Exactly the ten named helpers are extracted with `ast`; the production driver is not executed.
- **Density and momenta.** The instrument checks the actual `PD_POLY`, derivative momenta equal to −βM/2 (non-trivially nonzero), and a zero velocity momentum.
- **Divergence and EL.** It checks contraction = 2P, coordinate div M = 0, div K = actual coordinate P, and native EL = 0.
- **V5 comparison.** `gauss_jordan_solve` expresses the actual V6 row in the emitted V1 basis. The recorded weights are `[0,0,0,1]`, and the weighted native V5 combination is zero. By linearity, multiplying by −β/2 keeps it zero.
- **Controls:**
  - Wrong sign (momenta − βM/2) and wrong factor (momenta + βM) are both nonzero.
  - A doubled current fails, and P is not identically zero.
  - **Added G00² term:** with badL = −(P + G00²)/2, the native EL is nonzero and equals `[−∂₁²u₁, 0, 0, 0]`. I confirmed this by hand: ∂²(−G00²/2)/∂G00² = −1.
  - **Always-zero helper:** replacing `system.append(...)` with zero yields a zero operator, which would fail the non-null assertion. This is established by inference from the asserted zero output, not by re-running the canonical assertion under the mutation. That is logically adequate.
- **Implicit commutation in SymPy.** Cancellation relies on `doit()` canonically ordering mixed `Derivative` variables. That is an implicit smoothness assumption on the native side, which Lean handles explicitly. The recorded V5 string shows the uncanonicalized pairs that then cancel.
- **Status of the link.** It is exact, tested symbolic translation, not a kernel-certified bridge. The Wolfram "link" is only a substring check for `Transpose[actionRows]`, as disclosed.

## Mutation controls

I checked each named false identity independently and read the recorded diagnostics:

1. **Density normalization.** The error is at `density_identity` (l.28): the unsolved goal is −βP = −βP/2, false at β = 2, P = 1. The secondary failure in `lagrangian_increment` is also a genuine false identity.
2. **Orientation sign.** The error is at `momentum_eq` (l.70): `J 4 2 − J 3 3 = J 3 3 − J 4 2 ∨ beta = 0`, i.e. −F23 = F23 ∨ β = 0, which is false at the witness. The secondary `dualCurl_contraction` failure is also genuine.
3. **Mixed-derivative index.** The error is at `dualCurl_divergence_zero` (l.13): an unsolved combination of second derivatives.
   - Independently, with u = (0,0,0,x1·x2) the mutated j = 0 component is ∂_{x1}(−F23 + F13 − F12) = ∂_{x1}(x1) = 1 ≠ 0, so the statement is false.
   - The residual goal also contains unnormalized atoms (`Fin.succ 3` alongside `3`), so the Lean diagnostic alone does not prove falsity. The acceptance rests on the counterexample, which is what the record claims.
   - I disregarded the later `rewrite` and `simp` errors (l.39, l.70) as evidence.
4. **Current factor.** The error is inside `boundary_identity`, in its internal `hd` step (l.63): `1·S = ½·S` with S generic.
   - With u = (0, x1, 0, x3) and i = 0, S = 1, and div K becomes 2 while P = 1.
   - The later `change` failure in `current_normalization` (l.112) is secondary, and I did not count it.
5. **Eight paired statements.** Each is a concrete wrong value (P = 0, momentum +1, current 2, contraction 1, first variation 1 at β = 3, EL component 1, density −1 at β = −2, time momentum 1).
   - The instrument accepts a rejection only if `⊢ False` appears in `contract_control`. There are 16 occurrences, one message plus one output per mutant.
   - All eight correct partners pass.
6. **Admissible examples.** The zero test field, stationarity at β = −7 and β = 0, and smoothness of the affine background all pass.
   - The hypotheses `SmoothField` and `TestField` are satisfiable, so the universal statements are not vacuous.
7. **Acceptance logic.** The instrument rejects results containing import, syntax, heartbeat or memory errors, or timeouts; positives must also be free of warnings. The recorded outcomes are the expected ones.

## Optional improvements (non-blocking)
1. **`firstVariation_zero`** states `deriv … = 0`, which on its own would also hold for a non-differentiable function. Differentiability is proved separately in `relativeAction_hasDerivAt`, so the content is present, but a corollary `HasDerivAt (relativeAction beta u h) 0 0` would make the headline self-contained.
2. **Field-level witness.** The nonzero-density witness exists only at the jet level. `smooth_affine_background_positive` proves smoothness only. A Lean lemma `orientationJet (fieldJet (![0,x 1,0,x 3]) x) = 1` would show the same thing for an actual smooth field, which the mutation-4 note currently assumes.
3. **Mixed-index diagnostic.** The mutation-3 diagnostic could be made self-evidently false, for example by normalizing `Fin.succ k` to numerals or by adding a concrete-field false statement that reduces to `False`.
4. **Wolfram evidence** is only a substring check and should stay described that way; no change is needed.

## Limits of this review
- I did not compile anything, run the SymPy instrument, run Lean mutations, or check package pins.
- The build status, axiom-audit lines, diagnostics and the native `PASS` are the author's recorded evidence. I inspected them for consistency with the sources and checked their mathematical content by hand.
- Proof-script validity (as opposed to statement fidelity) rests on that recorded build. I found no `sorry`, `admit` or `axiom` in the four modules or the root.
- The mathematical calculations above (momentum entries, contraction, Bianchi cancellation, current factor, witnesses, counterexamples, native G00² target) were done by hand.
- This review does not cover the full D4 even-family bulk census, variable β, general null Lagrangians, D5, the `XFORM_EXTRA` spectrum, stability or interfaces, S11c, production or export work, or any systematic CAS bridge. It also does not assess boundary effects, which the packet correctly does not claim.
- This is one of the two required independent reviews; the other is still outstanding under the policy.
