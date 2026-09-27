# Independent fidelity review: S1–S4 conditional finite-solve sensitivity packet

## Verdict: **CLEAR** for bounded S1–S4

I found no blocking fidelity or mathematical defect. The Lean statements match the stated contract. Their hypotheses are explicit and not stronger than needed. The documents do not claim a physical inverse bound, a physical residual-to-flux certificate, conservation or convergence. The optional observations below would tighten wording or controls; none of them is needed for this verdict.

## Obligation-by-obligation assessment

**1. Residual, inverse bound, scaling and perturbation (`Residual.lean`)**
- `residual A b x := A x - b`. `residual_identity` gives `x - A.symm b = A.symm (residual A b x)`. `residual_bound` then gives `‖x - A.symm b‖ ≤ κ‖A x - b‖`, using `‖A.symm‖ ≤ κ` as a supplied premise. Invertibility is a supplied `X ≃L[ℂ] Y`; nothing infers it from floating-point rank, singular values, `lstsq` output or a condition number.
- `balanced A R D = D.symm.trans (A.trans R)`, so `z ↦ R(A(D⁻¹z))`. That is K = R M D⁻¹. `balanced_residual` shows `K(Dx) − Rb = R(Ax − b)`. `unscale_error` supplies the ‖D⁻¹‖ factor for converting an error in z into an error in c.
  - This matches the native code at `S11c_d_finite_scattering.py:371–377`: R = diag(1/row_scale), D = diag(column_scale), `lstsq` on the balanced matrix with right-hand side `rhs/row_scale`, then division by column_scale. The code's `scaled_residual = residual/row_scale` equals R(Mc − b) = Kz − g.
- `perturbed_residual_bound` correctly adds two reviewed T2 results: `residual_bound` for B (with `inverse_norm_bound` giving ‖B⁻¹‖ ≤ κ/(1−κε)) and `solution_error_bound`. The result is ‖x − A⁻¹b‖ ≤ κ/(1−κε)·(‖B x − (b+db)‖ + ‖db‖ + ε‖A⁻¹b‖).
  - B is supplied, together with the premise B = A + E.
  - κ ≥ 0 and ε ≥ 0 follow from the norm premises, and the theorem requires κε < 1.
  - No `CompleteSpace` premise is needed or silently discharged, because only T2's existence theorem uses completeness.

**2. Affine observation (`Current.lean`, `Pipeline.lean`)**
- `amplitude C d x = C_i x + d_i`, and `amplitude_difference` cancels d.
- `amplitude_error_bound` uses the explicit ℓ¹ mass and the bound Σᵢ‖Cᵢ‖.
- `flux_residual_bound` expands the flux about the *full* reference amplitude `amplitude C d (A.symm b)`, which includes d, as required.
- The native modal extraction (`trace -= incomingValues; solve(right, trace)`; open selection) matches a = P V⁻¹(T D⁻¹z − t_in).
- The native examples choose the coefficient sup norm, and Mathlib's default norm on `Fin n → ℂ` is also the sup norm. The 0.625 max-row-sum and the Cbound entry sum (dual ℓ¹ row norms) are therefore the correct operator norms for that choice. `mass` is defined separately, so it is not confused with the Pi norm.

**3. Full current**
- `flux_change_bound` uses the exact decomposition from `flux_add`: q(a+e) − q(a) = q(e) + Re(⟨a,Je⟩ + ⟨e,Ja⟩). Both interference terms and the quadratic term appear, giving the constant β(2·mass a·mass e + mass e²).
- `current_change_bound` and `flux_and_current_change` keep H separate, with γ·mass(a+e)². The split is correct: [q_J(a+e) − q_J(a)] + [q_{J+H}(a+e) − q_J(a+e)]. H is a supplied bound; nothing estimates it.
- There is no Hermitian or positivity premise.
- Empty index types work. A negative β is only possible when n is empty, where both sides are 0. That explains why `flux_change_bound` needs no `hβ`, and why `hβ` is present exactly where monotonicity in δ is used (`flux_change_budget`, `flux_distance_budget`).
- The composition chain residual → amplitude → fixed-J flux → fraction is `normalized_residual_bound`. It feeds `currentBudget` into `fraction_error_budget` as α, with the numerator reference flux at the exact solution.

**4. Fractions (`Fraction.lean`)**
- `fraction_difference` is exact.
- `fraction_error_bound` uses `d ≤ |j|`, `d ≤ |jh|` and `0 < d`, giving `d² ≤ |jh||j|`. The bound |Δn|/d + |n||Δj|/d² is correct, and signed denominators are allowed.
- `denominator_nonzero` requires the strict condition δ < |j|. The zero-margin witness (j=1, jh=0, δ=1) shows the non-strict version fails.
- Lean's totalized division is never relied upon, because every estimate requires a positive margin.
- The critical scalar witness `critical_inverse_margin` supports the claim that κε ≥ 1 is only an excluded domain for the estimate, not a universal singularity.

**5. Native fidelity**
- The selected AST blocks match the source: 10, 10 and 2 statements. The observation block ends at `incoming_flux`, so the `require(incoming_flux>0)` guard (line 394) and `outgoingFluxRatio` (line 400) were not executed, as disclosed.
- `native_full_current` recomputes Re(aᴴJa) per column. `diagonal_only_current` separates from it by 0.46.
- The residual sign mutation changes exactly one Sub node.
- The documents correctly separate the unscaled-M/sup-norm bound examples from the K and C_z = C_c D⁻¹ scaling identities. They do not present those examples as a certificate for any saved physical operator.

**6. Verification records**
- The root has 51 `#print axioms` lines, and each shows only `[propext, Classical.choice, Quot.sound]`.
- There are 18 control pairs plus 6 extras, giving 24 positives. The cross-term and quadratic-error positives use identical source, so there are 23 distinct statements.
- `verify.adjudicate` accepts a rejection only with exit 1, exactly one error attributed to `contract_control` containing exactly one `⊢ False`, no warning and no bad-instrument match. The recorded mutant outputs fit that pattern.
- 8 + 18 + 24 + 1 = 51 records.
- Native counts are 13 identities/examples, 14 wrong-formula controls and 9 bound examples. The native instrument is bound by hash (`e0858317…`). Run6 marks the native record `revalidated_recorded_native_evidence`, and the run1 ERROR is retained as a failure.

**7. Scope and exclusions**
Parameter conventions, quantifiers and exclusions match the sources. No overclaim was found in COVERAGE, FIDELITY or VERIFICATION.

## Blocking findings

None.

## Optional / nonblocking observations

1. **Several controls are pure arithmetic.** The `inverse_factor`, `small_residual_large_error` and `critical_inverse_margin` pairs, and the `zero_margin` mutant, are bare numeral statements like `‖(4:ℂ) − 0‖ = 4/1`. They never reference `residual` or the canonical definitions, so their link to the contract depends on the companion witnesses in `Controls.lean` (`inverse_factor_witness`, `critical_inverse_margin`, `zero_margin_witness`).
   - The `zero_margin` mutant (`|0−1| = 0`) does not encode the mutated claim "δ ≤ |j| ⇒ jh ≠ 0". The meaningful evidence there is the positive witness.
   - This is disclosed ("canonical witnesses/arithmetic"). If controls are revisited, a small fix is to negate the mutated general statement against the existing scalar witness. Examples: `¬ ∀ δ ≤ |j|, jh ≠ 0` at (1, 0, 1), or a κ-free residual bound at `quarter`.
2. **The native current is block-diagonal by end and restricted to open outgoing channels.** Line 390 fills only same-end entries, so cross-end entries are structurally zero. COVERAGE:26 ("complete supplied outgoing matrix J") would be clearer with that stated. The Lean result is unaffected because J is arbitrary.
3. **No non-Hermitian or indefinite J is exercised.** Both `fullJ` and the native open-channel block are Hermitian positive-definite. The theorems carry no such premise, so this is only a fixture gap.
4. **The realized native denominator does not depend on the solution.** It is −Re J_in,in, so η = 0 on the residual path. The denominator-error term matters only for supplied current changes. The native `fixed_denominator_bound` examples just divide by a positive constant and do not exercise that term.
5. **Some compositions are left for the user to assemble.** `normalized_residual_bound` uses only the unperturbed A. The perturbed path stops at the amplitude (`observed_perturbed_residual_bound`), and the H term is not composed into the fraction. Neither is claimed, and `flux_distance_budget` makes the assembly straightforward.
6. **`denominator_cases` proves the two cases cover everything but not that they are disjoint.** Disjointness is trivial but is described in COVERAGE:71 without a separate proof.

## Limits of this review

- I read and reasoned about the sources only. I did not compile Lean, run Python, recompute any hash, or inspect the actual logs or objects under `_scratch`.
- Build, axiom and control outcomes are taken from the recorded JSON. Their consistency with the instruments is checked, but not independently reproduced.
- The external Mathlib and transitive cache, the 605 preserved files and 243 old objects, earlier failed runs and the guard receipts are outside the packet. Their validation is the author's.
- I did not use `S11_lean_sensitivity_review_prompt.txt` or any author PASS label as evidence.
- Native physics was read as provenance only. I did not verify the orientation convention of the native `current` matrices against the s_e factors in `SHARED_PHYSICS.md:462–465`, or whether all open outputs belong to the intended H channel. Those are premises of the physics layer, not of S1–S4.
- The packet cannot establish:
  - that any physical K or M has a given inverse bound, or that the saved system is nonsingular;
  - floating-point or interval correctness;
  - discretization, quadrature or tail error;
  - channel completeness, conservation, continuum convergence or parent-theory accuracy.
- This is one non-author review. Per L4, a second independent review is still needed before fidelity review is complete.
