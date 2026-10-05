# NP1–NP4 fidelity review (nonlinear-pencil finite core): **CLEAR**

No fidelity corrections are required. I found no load-bearing gap: the formal statements match the bounded claims, and the stronger analytic reading is not slipped into any formal conclusion. The optional improvements listed below are not conditions of clearance.

**Packet identifier (not recomputed):** "S11 nonlinear-pencil finite-core NP1–NP4 fidelity contract v1". `MANIFEST.json` gives `aggregate_sha256` = `9c7e8e41dd1820647c3dc787ad171043944038c70d37b873699d414f49396e57` (created 2026-09-17T03:44:41Z, author Codex). Your prompt supplied no separate identifier.

## What I checked

### NP1 — `Modal.lean`
- **Types and data:** `PairingData` has separate module types X, Y, K over a field, with maps A: X→Y, V: K→X, W: Y→K. D is an actual linear equivalence `K ≃ₗ K`, and `pairing_eq` says W(A(V k)) = D k for every k.
- **Nothing extra is assumed:**
  - The name `derivative` carries no calculus content.
  - W is not assumed to annihilate L0.
  - K may be zero-dimensional, and nothing asserts that a pole exists.
- **Definitions:** `residueCandidate` is V∘D⁻¹∘W, `fieldProjection` is R∘A and `sourceProjection` is A∘R.
- **Algebraic results:**
  - `residue_sandwich` proves RAR = R, and the two idempotency theorems follow from it.
  - The range theorems give range(V) and range(A∘V), with both inclusions proved.
  - The `finrank_*` theorems give dim K when K is finite-dimensional, using the injectivity of V and of A∘V, which follows from D being invertible.
- **`residue_unique`:** its hypotheses are exactly range R ≤ range V and W A R = W, and the proof is correct.
- **`residue_of_inverse_coefficients`:**
  - Hypotheses: range V = ker L0, W L0 = 0, L0 R = 0 and L0 H + A R = I.
  - These are exactly the z⁻¹ and z⁰ coefficients of L·L⁻¹ = I for a simple pole (the right-inverse form, matching the kernel convention).
  - The proof only needs ker L0 ≤ range V.
  - The docstring states that the expansion is supplied, not constructed. The conclusion is only R = V D⁻¹ W, with no analytic claim.
- **Singular pairings:** these are excluded structurally, because the `K ≃ₗ K` field plus `pairing_eq` forces WAV to be bijective.

### NP2 — `Moments.lean`, `Laurent.lean`
- **Definition of `moment`:** (2πi)⁻¹ • ∮ over C(0,R). Mathlib's circle map runs θ from 0 to 2π, so the orientation is positive; `integral_sub_center_inv` = 2πi confirms this.
  - Positive R serves only to keep the path off zero (`nonzero_on_circle`). Orientation does not depend on the sign of R: C(0,−1) is still counterclockwise.
- **`moment_zpow`:** covers every integer exponent. The n ≠ −1 case uses `integral_sub_zpow_of_ne`, which holds for all centres and radii.
- **`moment_finite_laurent`:** covers an arbitrary `Finset` in a complete normed ℂ-space, using `integral_fun_sum` with per-term integrability. `moment_polynomial` follows by actual integration.
- **`double_log_expansion`:** I expanded (z⁻¹C₁ + z⁻²C₂)(A + zB) by hand and got exponents [−1, 0, −2, −1] on [C₁A, C₁B, C₂A, C₂B]. The moment is C₁A + C₂B, with the order kept.
- **`double_response_expansion`:** all eight terms match my hand expansion, with exponents [−1, 0, −2, −1, 0, 1, −1, 0]. The moment O₀C₁B₀ + O₀C₂B₁ + O₁C₂B₀ is complete.
  - The maps are rectangular: O is `o×n` and B is `n×u`.
  - No remainder is dropped: equality holds on the whole circle via `moment_congr`.
- **Matrix norm:** matrices use the scoped elementwise norm. In finite dimensions this does not affect the integral.

### NP3 — `Scalar.lean`, `Jordan.lean`
- **Real calculus throughout:** derivatives are computed with `deriv`, and every identity that inverts z assumes z ≠ 0 (or z ≠ ±1).
- **Scalar z²:**
  - `square_residue` = 0.
  - `square_log_moment` = 2 and `square_log_not_idempotent`.
  - `square_higher_coefficient` = 1: this is the moment of z·z⁻², i.e. the coefficient C₋₂.
- **Scalar z² − 1:**
  - `twoRoot_zeros` gives exactly ±1.
  - The radius-2 moment uses genuine `integral_sub_inv_of_mem_ball` integrals, with explicit proofs that ±1 lie inside the circle and off it. The result is 2, and it is not idempotent.
- **Jordan pencil:**
  - The pencil is z•1 − N, and `jordan_derivative` = 1.
  - `jordanResolvent` is defined as `doublePrincipal 1 N`. It is only *proved* to be a two-sided inverse for z ≠ 0 (`jordan_inverse_left/right`, using N² = 0), which is the right way round.
  - The determinant is z².
  - `jordan_kernel`: −N v = 0 ⇔ v₁ = 0, i.e. the kernel is the first axis.
  - `jordan_chain`: N e₂ = e₁ and N e₁ = 0. This matches the chain equation L(0)v₁ + L′(0)v₀ = 0, with v₀ = N v₁.
  - The resolvent moment is I and the logarithmic moment is I, which is idempotent with trace 2. The higher coefficient is N, and N ≠ 0.
- **Transfer:** entry (0,1) is z⁻², with residue 0.
- **Affine response:** (2+5z)(1+3z)/z² = 2/z² + 11/z + 15, giving residue 11. This agrees with the generic formula: 2·0·1 + 2·1·3 + 5·1·1 = 11. The frozen value is 2·0·1 = 0.

### NP4 — controls and the native link
- **Controls:**
  - Each mutant rewrites with a proved canonical fact, so an unsolved `⊢ False` means the mutated statement is equivalent to False. That is a genuine mathematical rejection.
  - I checked the recorded diagnostics for the mutants I opened (pairing normalization, projection normalization, full rank, singular pairing): each has one diagnostic in `contract_control` and one `⊢ False`. The instrument enforces this condition in `S11_lean_pole_contract_check.py:236`.
  - The mutant values are the right ones: 1 vs ½, 2 vs 1, a reversed-order entry of 0, −1, 11 vs 0/6/5. The 6 and 5 correspond to dropping O₁C₂B₀ and O₀C₂B₁ respectively.
  - The counts add up: 17 pairs plus 4 extra positives gives 21 positives. Six builds plus the audit plus 38 controls gives 45 records. The audit root has 55 `#print axioms` lines, and the recorded output shows only propext, Classical.choice and Quot.sound, including the wrapped list from run1.
- **Native link (`S11_lean_pole_source_check.py`):**
  - It has 19 checks and 4 translation controls.
  - The probe's `affineJordan` is `[[z,−1],[0,z]]`, i.e. zI − N.
  - The native `realization()` uses a = N (companion form of z²), b = e₂, c = e₁ᵀ and (zI − a)⁻¹. This matches Lean's sign convention and the (0,1) entry, and the native forcing and observation are 1+3z and 2+5z.
  - The native "both roots" count is a sum of local residues; Lean's is an actual circle integral. They are identified only at the pencil level, which is all that is claimed.
  - Historical evidence is read, not rerun, and is labelled as such.
- **Scope wording:** the assessment and fidelity documents keep the Keldysh/Fredholm existence results, J = P_X for an isolated semisimple zero, the argument principle, Riesz realizations and the homotopy count as literature or application hypotheses. None is described as Lean-certified.

## Findings (all optional, non-blocking)
1. **Singular-pairing control is disconnected from `PairingData`.** `Controls.lean:zero_pairing_not_invertible` and its paired control only show that x ↦ 0·x is not injective; they never mention `PairingData`. The exclusion itself holds by construction (the `K ≃ₗ K` field and `derivative_right_injective`). The "singular-pairing exclusion" row in `POLE_FIDELITY.md` therefore overstates what the control tests. A stronger version would prove `¬∃ d : PairingData ℂ ℂ ℂ ℂ, d.derivative = 0`.
2. **No example shows the hypotheses of `residue_of_inverse_coefficients` / `residue_unique` can hold together.** By hand they can: take `scaledPair` with L0 = 0, R = ½, H = 0. A second example with L0 ≠ 0 is diag(z, 1) with L0 = diag(0, 1). No control drops one of these hypotheses to show it is needed. Basis covariance of R follows from `residue_unique` but is not stated.
3. **`affine_response_residue` is not derived from the generic theorem.** It is proved from its own scalar expansion rather than as a 1×1 case of `double_response_moment`. I checked by hand that the two agree.
4. **The "nonzero response" facts are point values.** `square_response_nonzero` and `jordan_transfer_nonzero` evaluate at z = 1, and the `zero_residue_response` control tests only that value. The nonzero *singular* coefficient is carried by `square_higher_coefficient` together with `jordan_transfer_exact`. A control on the moment of z·`jordanTransfer` would test it directly.
5. **Docstring wording.** `Laurent.lean:17` says "B is the second derivative of the pencil", which reads as a claim. It should be worded as the application premise it is.

## Verification limits
- **Nothing was executed.** I had file-read and search tools only. I did not run Lean, Lake, Python or SymPy.
- **No hashes were computed.** Hash agreement is textual only: the hashes in `S11_lean_pole_contract_checks.json` and `S11_lean_pole_source_checks.json` match the `MANIFEST.json` entries for all seven Lean sources and both instruments/reports.
- **Mathlib was not inspected.** It isn't in the packet, so I checked lemma names and meanings from memory of the library. In particular, I could not confirm that the `convert!` tactic used in `jordan_derivative` exists at the pinned commit. Its success rests on the recorded build log (exit 0, empty output, warnings treated as errors). The statement it proves is mathematically correct.
- **Author evidence I rely on:** the build, reuse, axiom and preservation results (run1–run3), the 40 VC + 5 T1 object hashes, the package pins and `S11_lean_pole_validation.json`.
- **Not reviewed:** the historical 89-check native repair, beyond the selected functions quoted above.
