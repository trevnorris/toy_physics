**Verdict: CLEAR** for the bounded D3 contract J1–J4, with no blocking findings. This is one of the two independent non-author reviews the contract requires; J1–J4 stays open until the second one clears.

**Packet reviewed:** `MANIFEST.json`, revision "S11 D3 quadratic invariants J1–J4 fidelity contract v1" (author Codex, created 2026-09-16T15:01:18Z), aggregate SHA-256 `015d0249e7f5cf911323cdd3c787f5d9fb69a45104515fc958be5e53877e9e6c`. I took this identifier as supplied and did not recompute any hash.

## What I actually did
- **Read:** the manifest, `FORMALIZATION_POLICY.md`, the S11 `README.md`, and `D3_COVERAGE.md`, `D3_FIDELITY.md`, `D3_VERIFICATION.txt`.
- **Read in full:** all six `S11D3Invariants` modules, the audit root, the generator and both check scripts.
- **Read the records:** `_source_checks.json` and `_validation.json` in full; in `_contract_checks.json`, the axiom-audit output, both source mutations and a sample of the statement controls (all 29 records' names and outcomes via search).
- **Native sources:** `compute_q9` and its helpers in the SymPy audit, the Q9 and `P_D` sections of `S11_SHARED_PHYSICS.md`, the Wolfram `buildInvariantCensus` source, and step Move 2.
- **By hand:** I checked the polynomial coefficients, generator orientation, span equality and a sample of certificate steps (details below).
- **Not done:** I did not compile Lean, run either Python instrument, run SymPy or Wolfram, or change any file. Build, axiom and mutation results are therefore author-reported. They are consistent with the sources I read.

## Scope and assumptions checked
- **Pointwise densities only.** The contract classifies quadratic densities before any EL or divergence quotient.
- **Nothing extra assumed.** There is no positivity, symmetric-G or transverse restriction, no field equation, and D=3 is supplied rather than derived. I found no hidden hypothesis in any theorem statement.
- **No admissions.** A search of `lean/s11` for `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by`, `extern`, `unsafe` and `opaque` found nothing.
- **Standard axioms only.** The recorded audit lists 40 roots, matching the 40 `#print axioms` lines, each using only `propext`, `Classical.choice` and `Quot.sound`.

## Findings

**1. The domain really is every real quadratic form on all real 3×3 matrices.**
- `Quad` is mathlib's `QuadraticForm ℝ (Matrix (Fin 3) (Fin 3) ℝ)` (`Quadratic.lean:14`), not a list of candidate polynomials.
- `quadratic_representation` (`Quadratic.lean:109`) rewrites Q G as `Q.associated G G` and expands G over the nine matrix units (`frame_expansion`, proved entrywise).
- Its coefficients are b(eᵢ,eᵢ) on squares and b(eᵢ,eⱼ)+b(eⱼ,eᵢ) on cross terms, which is the correct expansion of b(G,G).
- The 45 monomials run over row-major index pairs (i ≤ j) in lexicographic order. This matches Q9's pinned `MONOMIAL_ORDERING`, which the source check also asserts.
- So the associated-bilinear representation covers the whole space, and it assumes nothing about which invariants exist.

**2. The rotation constraints come from admissible members of the full SO quantifier.**
- `SOInvariant` quantifies over every R with RᵀR=1 and det R=1, and every G (`Rotation.lean:12`).
- Each of the 42 constraints `e0`–`e41` (`Constraints.lean`) applies `hQ` to `rotationXY 0 1`, `rotationYZ 0 1` or `rotationXY (3/5) (4/5)`.
- Each use carries a proof that the matrix is proper (`rotationXY_proper` / `rotationYZ_proper`, general in a²+b²=1, instantiated by `norm_num`).
- Every image matrix is checked entrywise through `conjugate_apply`.
- **Spot checks I did by hand:**
  - Images of `e10`, `e32` and `e41`, plus the resulting linear equations (e41's coefficient −544/625 = 81/625 − 1 is correct).
  - The combinations `h1`, `h17` and `h24`.
- The generator (`S11_lean_d3_generate.py`) only chooses the certificate. Lean has to re-derive every equation and `linear_combination`, then close `invariant_polynomial` with `ring`.
- **Why three rotations can suffice** (my own argument, not in the packet): the two quarter-turns generate the cube's rotation group. The 3-4-5 rotation has infinite order, so invariance under it extends by continuity to every rotation about z; conjugating by the quarter-turns gives every rotation about x and y, and those generate SO(3). So a three-dimensional answer is expected, not a sampling accident.

**3. Full-O sufficiency is proved separately and does not use the finite tests.**
- `invariantForm_O` (`Rotation.lean:108`) holds for every R with RᵀR=1, using three identities: `trace_conjugate` (cyclic trace), `conjugate_mul` and `conjugate_transpose`.
- `SO_iff_O` then goes classification → `invariantForm_O` in one direction and restriction in the other.
- `Orthogonal` requires only RᵀR=1. For square matrices this is the full group O(3), including det = −1.

**4. The normalization is exact, and the coefficients are unique.**
- **Monomial definitions:**
  - `traceSquare` = x₀²+x₄²+x₈² + 2(x₀x₄+x₀x₈+x₄x₈).
  - `traceOfSquare` puts the cross coefficient 2 on the pairs (1,3), (2,6), (5,7), i.e. G₀₁G₁₀, G₀₂G₂₀, G₁₂G₂₁.
  - `frobeniusSquare` is the sum of all nine squares.
- **Each is proved equal to its trace formula:** `traceSquare_apply`, `traceOfSquare_apply`, `frobeniusSquare_apply`. So tr(G²) and (tr G)² are kept distinct.
- **The certificate's parameterization is consistent.** With `v = ![c₄/2, c₁₁/2, c₉]`, the x₀² coefficient is a+b+c (`h0`), x₀x₈ is 2a (`h8`), and x₂x₆ and x₅x₇ are 2b (`h21`, `h37`).
- **Uniqueness:** `invariantForm_injective` evaluates at E₀₁, E₀₀ and E₀₀+E₁₁ (giving v₂, then v₀+v₁+v₂, then 4v₀+2v₁+2v₂), which is a sound argument. `SO_unique` builds on it.

**5. The submodules are actual submodules, and the reflection-minus part is zero.**
- `soSpace`, `oSpace` and `oddSpace` are genuine `Submodule ℝ Quad` objects.
- `so_eq_o` is equality of submodules, and `odd_eq_bot` identifies the odd submodule directly rather than as a difference of counts.
- `odd_classification`: O-invariance under R₀ gives Q = Q∘R₀; oddness gives Q∘R₀ = −Q; together Q = 0.
- The dimensions come from `soEquiv : (Fin 3 → ℝ) ≃ₗ soSpace`, giving `census = (3,3,0)`.
- R₀ = diag(−1,1,1) is proved orthogonal with det = −1, and it matches the native R₀.

**6. The count and omission controls have admissible positive examples.**
- The positive examples all appear in the recorded run:
  - `nonzero_invariant_exists` (a nonzero invariant exists).
  - `zero_invariant` (the zero form is admissible).
  - `invariantForm_O ![-2,3,-5]` (negative coefficients).
  - Uniqueness on that same example.
- The omission controls prove that no form lies in the span of the other two. Combined with dimension 3, dropping any one form leaves a proper subspace.
- `single_entry_not_invariant` shows that the real rotation hypothesis is not vacuous.
- **The two source mutations fail for the intended reason:** "unsolved goals" in the named `_apply` theorem, showing the 3-vs-2 and −2-vs-2 coefficient mismatch.
- **The eight statement mutations each reduce to `⊢ False`.**

**7. The native link compares full spans, not just counts.**
- `S11_lean_d3_source_check.py` does four things:
  - Runs the unchanged native `compute_q9(3)`, pulled from the source with `ast`.
  - Asserts the monomial ordering.
  - For V1 and V2, checks that stacking the native basis with the three forms keeps rank 3 and that the RREFs are equal.
  - Checks V6 is 0×45, the V6 operator is the identity, `P_D` is zero and the V4 residual is zero.
- **I checked the recorded native V1/V2 RREF polynomials by hand:**
  - Row 1 is tr(G²).
  - Row 2 is ((tr G)² − tr(G²))/2.
  - Row 3 is |G|²_F − tr(G²).
  - Their span is exactly {tr², tr(G²), |G|²_F}.
- **The injected `QG_ALL` is harmless.** The native `declared_symbol` returns the same `sp.Symbol(name, real=True)` and additionally registers it in the native symbol table.

**8. The generator orientation and its control are correct.**
- **SymPy:** the native generator X has +1 at (a,b) and −1 at (b,a). δG = XG − GX is the derivative of R G Rᵀ at R = e^{tX}.
  - `action_rows[k]` is the image of monomial k, so the coefficient condition is Aᵀc = 0. The `.T` at line 517 is therefore correct.
  - The reflection rows form a diagonal matrix, so their orientation does not matter.
- **Wolfram (source inspection only):** it has the same generator sign, the same δ, and `Transpose[actionRows]`.
- **The old-orientation control:** it keeps counts 3/3/0 but fails the span comparison. This is plausible, since the untransposed system solves for the dual-weighted forms (off-diagonal weights differ by a factor of 2). I did not rerun it.

**9. The records match the stated boundaries.** Wolfram is labeled source inspection only, and the Python link is described as tested translation, not kernel-certified. No production driver, export, EL or S11c work is claimed. The mutation adjudication requires a real mathematical diagnostic and rejects instrument failures.

## Blockers
None.

## Optional improvements (none affects the bounded theorem, coverage or fidelity)
1. **Coefficient link is checked only by hand.** Nothing mechanically compares the Lean monomial coefficients with the Python `forms` rows; I checked the index correspondence manually. A one-line assertion of the expected 3×45 coefficient matrix would make that link explicit.
2. **The finite certificate itself is never mutated.** An isolated mutation dropping `e41`, or changing `c 4/2` to `c 4` in `invariant_polynomial`, would test it directly. Its soundness does not depend on such a control.
3. **"Independent" is the generator's claim.** The "42 independent constraints" wording in `D3_FIDELITY.md` comes from the generator's rank assertion; Lean does not prove independence and does not need to.
4. **Constraints build is fragile.** It took 283.8 s against a 300 s timeout. In run 3 it was reused by hash, so its compile evidence comes from the earlier run (instrument `f2c526…`), whose report is not in this packet. This reuse is disclosed and guarded by hashes.

## Limits of this verification
- I compiled nothing, ran no Python or SymPy, recomputed no hashes, and executed no Wolfram code.
- I did not independently verify the rank or solvability of the whole 42-equation certificate, only a sample. Its correctness rests on Lean accepting `Constraints.lean`, which is author-reported: PASS, empty output, with warnings treated as errors.
- The native span equality and both in-memory controls are taken from the recorded JSON plus my reading of the instrument and native source. My hand check of the recorded RREF rows agrees with them.
