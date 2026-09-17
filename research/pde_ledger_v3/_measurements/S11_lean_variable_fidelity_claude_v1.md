# Independent fidelity review: S11 VC1–VC4 (variable-coefficient bulk identities and flat-interface normal-slice/trace contract)

## Verdict: **CLEAR**

The Lean statements say what the contract says, with the stated conventions, hypotheses and limits. I found no required fidelity corrections. A few optional, non-blocking observations are at the end.

**Packet identifier (as supplied, not recomputed):** `MANIFEST.json`, revision "S11 variable-coefficient and flat-interface VC1–VC4 fidelity contract v1", aggregate_sha256 `8152df86a31ecaabe89713a3fb178c7c469885f905b9e486207bbc695c02cf45`.

## What I checked by reading, and how

I read the manifest, `FORMALIZATION_POLICY.md`, the three `VARIABLE_COEFFICIENT_*` documents, all five `S11VariableCoefficients` modules and the audit root. I also read the reused definitions and theorems:
- **S10:** `Vec`, `Point`, `Jet`, `axis`, `basisJet`, `coordDeriv`, `fieldJet`, `antisym`, `SmoothField`
- **D3:** `S11D3Bulk` Action/Bulk/Calculus/Variation
- **D4:** `S11D4Odd` Action/Calculus/Boundary/Variation

I also read the native helper subset in `scripts/S11_stray_longitudinal_sympy_audit.py`, both new instruments, and the three reports. I redid the key calculations by hand.

### Conventions
- **Coordinates:** `Point D = Fin (D+1) → ℝ`, with index 0 as time. `fieldJet u x j i = coordDeriv j (u · i) x`, so `G_ij = ∂_i u_j` sits at `J i.succ j`. The native `coordinate_substitution` uses the same rule: `G[i,j] ↦ ∂_{x_i} u_j`.
- **Time momentum:** both `momentum_eq` theorems give `Fin.cases 0 …`, so the time momentum is identically zero for any jet. Profiles are `Point → …` and may depend on time. No stationarity, sign, nonzero, transverse or background-equation premise appears anywhere.
- **Fixed profiles:** in `D3.momentum` and `D4.momentum` the profile is only evaluated at `x`, never varied. Both are genuine `deriv` definitions, and `momentum_identity` (proved by `rfl`) ties them to the reviewed constant-coefficient derivatives. `pointwise_variation` reuses the reviewed `lagrangian_variation`.

### VC1 (D3)
- **Momentum:** `S11D3Bulk.momentum_eq` gives `p_{r i} = −(a δ_{ri} div + b ∂_i u_r + c ∂_r u_i)`. This matches `∂L/∂G_{ri}` for `L = −½[a(tr G)² + b Σ G_ij G_ji + c Σ G_ij²]`.
- **EL:** `EL_j = −Σ_r ∂_r p_{rj}` works out by hand to the claimed formula. The b-gradient term is `Σ_r (∂_r b) ∂_j u_r`.
- **Transposed index:** `D3.correction` (`D3.lean:26-29`) uses `fieldJet u x j.succ i` (= `∂_j u_i`) for the b term and `fieldJet u x i.succ j` (= `∂_i u_j`) for the c term. Both agree with the actual momentum, so the index is not transposed.
- **General theorem:** `eulerLagrange_eq` is stated for the derivative-defined `eulerLagrange` and every smooth profile and field. A wrong correction definition could not have been proved.
- **Constant coefficients:** `constant_profile` (by `rfl`) plus `S11D3Bulk.eulerLagrange_eq` recover `(a+b)∇div + cΔ`.
- **Null residual:** `null_profile_residual` (a = −b, c = 0) gives `(∂_j a) div u − Σ_i (∂_i a) ∂_j u_i`, as claimed.
- **Nonzero example:** `d3_nonzero_response` uses `u = (0, x2, 0)`, `a = x1`. By hand, div u = 1 and the residual is (1, 0, 0). Lean proves component 0 = 1 for all x, and the native check asserts the full vector.

### VC2 (weighted divergence and D4)
- **Product rule:** `weighted_divergence` is dimension-generic and correct: `div(aJ) = ∇a·J + a div J`.
- **D3 weighted density:** `D3.weighted_null_density` states `L = −½ div(aJ) + ½ ∇a·J`. The current `J_i = Σ_j (u_i ∂_j u_j − u_j ∂_j u_i)` matches `boundaryCurrent`.
- **D4 weighted density:** `D4.weighted_odd_density` states `L = −½ div(βK) + ½ ∇β·K`, with `K_i = ½ Σ_j u_j M_ij` (`currentJet`). Both keep the +½ gradient/current contraction.
- **D4 momentum and EL:**
  - `momentum_eq` gives `p_{ri} = −(β/2) M_{ri}` with a zero time row. Since it is proved from the actual derivative, `dualCurl` is exactly `∂P/∂G`; I spot-checked `M_01 = F23` and `M_02 = −F13`.
  - `D4.eulerLagrange_eq` gives `EL_j = ½ Σ_i (∂_i β) M_ij`, using the reviewed `dualCurl_divergence_zero`. The index order is correct (derivative row i, component j).
  - `constant_profile` reuses `eulerLagrange_zero`.
- **D4 witness:** `u = (0,0,0,x3)`, `β = x1` gives `F23 = 1`, row `M_0· = (0,1,0,0)`, so EL = (0, ½, 0, 0). Lean proves component 1 = ½, and the native check asserts the full vector.
- **Normalization:** by hand, `Σ ε_ijkl G_ij G_kl = 2P` holds, and so does `M_ij = Σ_kl ε_ijkl G_kl`. The earlier instrument asserts `P_D = P` and `Σ G_ij M_ij = 2P`. The epsilon form itself is prose only (see observation 5).

### VC3 (interface)
- **Hypotheses of `split_integration_by_parts`** (`Interface.lean:12-25`):
  - separate `pm` and `pp` defined on all of ℝ, each with `HasDerivAt` at every point of the closed `uIcc`, so both have two-sided derivatives at c;
  - one test function `h` with the same `dh` on both intervals, which gives a common trace (and differentiability at c);
  - interval integrability of `dm`, `dp` and `dh` on the appropriate sides.
- **Conclusion:** it keeps `pp b·h b − pm a·h a + (pm c − pp c)·h c` and both interior derivative integrals. I re-derived it from two applications of the standard integration-by-parts identity; the signs are right.
- **Orientation:** oriented integrals, so a ≤ c ≤ b is not required. The documents correctly reserve the adjacent-interval reading for a ≤ c ≤ b.
- **Endpoint corollary:** `compact_endpoint_split` adds only `h a = 0` and `h b = 0`.
- **Integral witness:** `interfaceWitness_eq` evaluates an actual integral: 2·1 + 5·(−1) = −3 = (2−5)·h(0). The representatives take unequal values at c.
- **Pairing criterion:** `jumpPair_zero_iff` quantifies over every `h : Vec D` and proves "pairing zero for all h iff minus = plus".
- **Tractions:** `normalFlux n p j = Σ_i n_i p_ij`. `D3.traction_eq` and `D4.traction_eq` match the momenta. The minus-to-plus orientation and the statement that stiffness traction has the opposite sign (L = −W) are both correct.
- **Limits:** the stated exclusions (no surface action, weak-trace existence, trace surjectivity, field continuity, Fubini lifting or multidimensional transmission) are accurate. The documents do not add spurious completion requirements. Nothing formally links the scalar integral theorem to `traction`, and the documents say so.

### VC4 (controls)
- **Mutant rejections:** all 12 paired mutants in `S11_lean_variable_contract_checks.json` exit 1 with exactly one diagnostic, `unsolved goals ⊢ False` in `contract_control`. None is an instrument failure. All 16 positives exit 0 with warnings treated as errors.
- **Mutant values are meaningful:**
  - 0 is the constant-coefficient prediction;
  - a negative value is the reversed correction;
  - 1 is the dropped ½;
  - −2 is the `−βM` momentum;
  - 0 and ±3 are the omitted or reversed jump.
- **What the witnesses rest on:** they evaluate the actual derivative-defined EL, tractions and integral through the general theorems.
- **Axiom audit:** 31 entries, all `[propext, Classical.choice, Quot.sound]`. The instrument also scans for `sorry`/`admit`/`axiom`.
- **Native D3 identification:** the V1 basis is solved to the three forms, and the result is asserted equal to `a tr² + b tr(G²) + c tr(GGᵀ)`.
- **Native D4 identification:** `P_D = P` is asserted.
- **Native derivative:** the prescribed-profile divergence is computed by the new instrument itself, not by `euler_lagrange_from_placeholders`. Native profiles are spatial only; native fields depend on (x, t).
- **Index control:** with `b = x1`, `a = c = 0`, `u = (x2,0,0)`, I get EL = (0,1,0) by hand, and the transposed formula gives 0, matching the report.
- **Adequacy:** the suite is adequate for its claims. The general identities are kernel proofs; the controls show the witnesses are nonvacuous and sensitive to the signs, factors and jumps.

## Optional, non-blocking observations
1. **"Exactly one" is not enforced.** The instrument accepts a rejection when *at least one* `⊢ False` diagnostic is in `contract_control` (`bool(intended)`). "Exactly one" is true of the recorded data but is not checked by the instrument.
2. **The Lean weighted-correction control is generic, not D3/D4.** `nonzero_weighted_correction` works in D = 1 with a constant J. That the D3/D4 gradient/current contraction is nonzero is shown only natively. This is disclosed (verification item 6). A Lean D3/D4 instance would be a stronger optional addition.
3. **No factor/sign mutant for the weighted density's ±½ coefficients.** The exact theorems plus the omission controls already imply that any other coefficient on the gradient term would be detected. An explicit mutant is optional.
4. **The Lean D3 witness does not detect the index transposition.** Sensitivity comes from the proved `eulerLagrange_eq` plus the native control, as disclosed.
5. **The epsilon statement is prose only.** "Epsilon summed over all indices equals 2P" is correct but is not itself machine-checked here. The checked equivalent is `Σ G_ij M_ij = 2P`.
6. **Scratch path mismatch.** `VERIFICATION.txt` names `_scratch/S11_lean_variable/verification_run2`, while the instrument writes to `lean/s11/_scratch/variable_verification`. Probably an outer wrapper directory; this does not affect the evidence.

## Limits of this review
- **Nothing was run.** I did not run Lean, Lake, Python/SymPy or any hash computation, and I did not recompute the manifest aggregate or any file hash.
- **Author evidence only:** builds, the axiom audit, the fresh run2 controls, the reuse guards, the package pins, the Mathlib hashes and the native results are the author's reported evidence. I checked them for internal consistency: 40 reuse markers, 68 records, and diagnostics matching the stated sources.
- **Mathlib and proofs:** Mathlib isn't in the packet, so I did not inspect the exact pinned form of `integral_mul_deriv_eq_deriv_mul`. Correctness of `Interface.lean` therefore rests on the reported build plus my own derivation of the same identity. Tactic proofs were reviewed for plausibility only; kernel acceptance is the author's reported result.
- **Scope:** my hand calculations cover the statements and witnesses above, not every intermediate tactic step. This clearance is limited to the bounded D3-family and D4-odd VC1–VC4 contract. It says nothing about the full S11c operator, multidimensional transmission, or CAS-wide fidelity.
