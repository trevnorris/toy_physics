# Review of the P1–P4 retained-order bookkeeping packet

**Verdict: CLEAR** for the bounded P1–P4 contract. I found no blocking findings.

I did not reuse the author's PASS label or anyone else's review. I read the policy, the coverage, fidelity and verification documents, and all seven Lean sources. I also read both check instruments, the recorded contract, native, validation and resource reports, `verify.py`'s adjudication and audit code, the actual native helpers, and shared-physics §3c–§3d.

## Obligation-by-obligation assessment

**1. Rectangle, path and non-identification (P1)**
- **Independent inputs:** `rectangle` takes independent real `eta` and `sigma`.
- **Delta:** `delta` = `eta·a10 + (sigma·a01 + eta·sigma·a11)`. This matches §3c exactly: ΔA = A_η + A_σ + A_ησ, where A_η is the zero-jet contrast term. It also matches the native split: `zeroJetContrast` is grade (1,0) and `firstJet` is every grade with a σ-degree.
- **Path:** `rectangle_path` at `eta=t`, `sigma=r·t` gives the path coefficients `a00`, `a10 + r·a01` and `r·a11`, for every real `r`, including zero and negative values.
- **Conjugation:** real scalars enter only as `((x:ℝ):ℂ)`, and `Complex.conj_ofReal` is what closes `quadratic_exact`. The casts therefore preserve conjugation.
- **Non-identification:** `path_does_not_identify_rectangle` shows the path collapses the pair (−r·v, v). That pair is nonzero exactly when `v ≠ 0`, which requires a nonempty index type. The coverage document states this caveat honestly.
- **Width ratio:** identifying `r` with W̄₀/L_W is an application premise, as the fidelity record states.

**2. Full quadratic contraction (P2)**
- **All 27 terms:** I enumerated every ordered triple (i, j, k) with i+j+k = d by hand. `q0`–`q6` have 1, 3, 6, 7, 6, 3 and 1 terms, each `pair B_j a_i a_k` with the left argument conjugated. They are all present and correct.
- **Match to the source:** `q2` equals §3c's J_H^(2) term for term. It includes both cross terms, the B₁/B₂ current variations and the a₂ terms.
- **Coefficients are pinned:** `quadratic_exact` holds for every real `t`. A real polynomial identity forces each coefficient to be unique, so the definitions cannot silently redistribute terms between degrees.
- **Remainder:** `retained_with_remainder` gives the exact degrees 3–6 of the supplied degree-two polynomials. It is not a bound on the parent theory, and the documents say so. §3c's "…" in B_H and J_T,in is correctly handled as a restriction to supplied data.
- **Generality:** there are no Hermitian or positivity premises, and empty `Fintype` index types are allowed. The empty-space positive control checks this.
- **Hand check:** the scalar example gives (1, 8, 31, 72, 107, 96, 45), which is correct.

**3. Quotient (P3)**
- **Field algebra:** `quotient_equations`, `quotient_unique` and `quotient_residual` are generic over any field and each assumes `j0 ≠ 0`. I verified the residual's t³ and t⁴ terms by hand.
- **Native match:** the native `quotient` loop computes the same recurrence on complex diagonal entries.
- **Zero-leading cases:**
  - `zero_leading_obstruction` covers `j0 = 0`, `n0 ≠ 0`.
  - `zero_model_two_solutions` exhibits distinct constants for the all-zero leading equation. General singular quotients are explicitly excluded.
- **Epsilon:** `epsilon_cancels` needs only `eps ≠ 0` and relies on Lean's totalized division. `scaled_denominator_nonzero` is a separate statement. The packet does not present totalized division as a physical denominator check.

**4. Observables (P4)**
- **Total minus baseline:** `subtracted_total` keeps both cross terms. `subtracted_eq_induced_iff` is exact.
- **Baseline-free case:** `induced_coefficient` and `induced_low_coefficients` give q0 = q1 = 0 and q2 = B0[a1,a1], independent of a2, B1 and B2.
- **Witnesses:** the omitted-parent-second-order witness gives q2 = 6 = 3 + 3, which is correct. The difference 31 − 5 = 26 versus 4 matches the native total-minus-baseline versus induced fixture.
- **Scope:** no physical baseline or value is asserted.

**5. Native identification**
- **Unrestricted product:** `current.multiply` is an unrestricted bidegree convolution; I confirmed this in the source. `quadratic` computes adjoint · metric · amplitude, which matches Lean's `star x ⬝ᵥ (J *ᵥ y)` convention.
- **Path map:** `lambda_series` uses `float(ratio)**b`, so the ratio is real.
- **Quotient:** it takes each incident column's diagonal of the full contraction, keeps complex values and has no zero guard.
- **Invalid-domain fixture:** it is isolated and only asserts a nonfinite value.
- **Counts:** the native check file contains 30 identity checks, 11 wrong-formula rejections and 1 domain witness, matching the report.

**6. Mutation controls and provenance**
- **Mutant logs:** all 16 recorded mutant logs contain exactly one `error: unsolved goals ⊢ False` inside `contract_control`, and nothing else.
- **Positive logs:** every positive log is empty (the empty-file hash).
- **Adjudication:** `verify.adjudicate` requires exactly one error, one `⊢ False` and no warning. It rejects instrument failures (unknown identifiers, sorry, memory problems and similar).
- **Counts:** 7 builds + 16 rejections + 20 positives + 1 native revalidation = 44 records. The positives are 18 distinct statements because three q2 pairs reuse one positive. The root file contains 33 `#print axioms` lines, and the log shows only propext, Classical.choice and Quot.sound.
- **What the controls are:** they are numeric instances decided through canonical witnesses. They are not independent rederivations or source-replacement mutants, and the packet says so.
- **Native lineage:** the validation file records the native run1 pass inside a guard run that failed for another reason (guard exit code 1). Run7 revalidated its hashes rather than rerunning it, and the packet is transparent about this.

**7. Scope and provenance**
- **Internal consistency:** the source, instrument and report hashes are consistent across MANIFEST, validation, the native report and the contract report.
- **Dependencies:** there are six direct Mathlib import pairs and fifteen package pins. The transitive cache is declared as a pinned baseline, not a fresh rebuild.
- **Exclusions:** they are complete and consistent with §3c–§3d, which rule out strong-edge extrapolation, conversion values, conservation, channel completeness and a scattering solve.

## Nonblocking observations
1. **`path_ratio` control:** it evaluates `rectangle` at (η, σ) = (1, 2), not the path polynomial `amplitude … t` with a wrong `r²·a11`. So in Lean it is really a rectangle-evaluation control. Sensitivity to a wrong ratio rests on `rectangle_path` being a proved identity plus the native `wrong_mixed_ratio_power` check. *Optional fix:* state the control on `amplitude a00 (a10 + r•a01) (r•a11) 1` with `r = 2`.
2. **Noninjectivity witness:** `path_does_not_identify_rectangle` shows the collapse but has no nonzero instance. *Optional fix:* a scalar witness (v = 1, r = 2) where the rectangle is nonzero off the path, e.g. at η = 1, σ = 0.
3. **Order of operations:** Lean applies the path and then contracts; the native code contracts in bidegree and then applies the path. They agree because evaluation at a real ratio is a ring homomorphism that commutes with conjugation. That agreement is tested on one two-channel fixture but is not a Lean theorem. This is within the declared translation boundary.
4. **Upstream current truncation:** the native current series come from `gram`, which uses the engine's cutoff `J.multiply`. So the supplied B coefficients are themselves truncated to rectangle grades. Lean correctly treats them as supplied. The phrase "unrestricted convolution" applies only to the downstream `quadratic`.
5. **Minor control and presentation points:**
   - `zero_model_two_solutions` is stated over ℝ only.
   - Its paired control checks only the two products, not `0 ≠ 1`; the unpaired `zero_nonunique_positive` covers that.
   - No theorem combines `induced_low_coefficients` with the quotient to give the induced fraction's c2 = B0[a1,a1]/j0 from §3d. It is not claimed and would be a one-line corollary.
6. **Native quotient inputs:** it uses only the diagonal of the incident series and the complex (not real-part) numerator diagonal. This is documented; reading the result as a physical real fraction needs reality premises.

## Limits of this review
- **Nothing executed:** I had no shell. I did not compile the Lean, run any instrument, or independently recompute any SHA-256 hash. The pass/fail results and hash equalities come from the recorded reports and are only checked for internal consistency.
- **Hand verification:** I checked the proofs and statements by reading them, and the coefficient enumerations, witness arithmetic and control values by hand.
- **Outside the packet:** the transitive Mathlib cache, the historical proofs, the original logs and objects behind the hash records, and the guard implementation (`scripts/s11c_guarded_run.py`).
- **Translation boundary:** native correspondence rests on small binary-exact fixtures and is outside Lean's kernel. It says nothing about arbitrary floating-point behaviour or the production pipeline's actual operands.
- **What the packet cannot establish:** that the supplied polynomials are the physical ones, that r = W̄₀/L_W, reality or conservation of the currents, a nonzero physical denominator, or any parent-theory Taylor or strong-edge claim.

This is one of the two required independent reviews. No files were edited and nothing was committed.
