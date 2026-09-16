**Verdict: CLEAR**

Reviewed packet identifier (as supplied, not independently recomputed): **S11 D3 quadratic invariants J1–J4 fidelity contract v1**, `aggregate_sha256` `015d0249e7f5cf911323cdd3c787f5d9fb69a45104515fc958be5e53877e9e6c` (`MANIFEST.json`).

This is a statement-fidelity review of J1–J4 only. It is not a J1–J4 completion claim. No required fix was found. No file was changed.

---

## Scope and assumptions checked

The contract in `lean/s11/D3_COVERAGE.md` and `lean/s11/D3_FIDELITY.md` is a classification of **pointwise homogeneous quadratic densities** on **all real 3×3 matrices**, under `G ↦ R G Rᵀ`, before any Euler–Lagrange or total-divergence quotient.

Checked against the Lean object, native Q9, and shared physics/step specification:

| Claimed assumption | What is actually formalized |
|---|---|
| Domain is every real quadratic form on all real 3×3 matrices | `Quad := QuadraticForm ℝ Mat` with `Mat := Matrix (Fin 3) (Fin 3) ℝ` (`Quadratic.lean`) |
| `G_ij = ∂_i u_j`, row-major | Lean `coordinates`; native `g_1,…,g_9`; shared physics §1 and Q9 monomial pinning |
| Full SO / full O, not a finite rotation sample as the theorem | `SOInvariant` / `OInvariant` quantify over every `Proper` / `Orthogonal` matrix (`Rotation.lean`) |
| No positivity, symmetry of `G`, transverse ansatz, field equations, or derived `D` | None of these appear as hypotheses; `Fin 3` is supplied |
| Basis is `a (tr G)² + b tr(G²) + c tr(G Gᵀ)` with those polynomials distinct | `traceSquare_apply`, `traceOfSquare_apply`, `frobeniusSquare_apply`, `invariantForm_apply` |
| Census is dimensions of actual submodules, not CAS counts | `soSpace`, `oSpace`, `oddSpace` (`Census.lean`) |
| Native link is span equality of Q9 V1/V2/V6, not dimensions alone | `_measurements/S11_lean_d3_source_check.py` |
| Wolfram is source inspection only | `Transpose[actionRows]` at `mathematica/S11_stray_longitudinal_mathematica_audit.wl:910` |

The intended native object is shared-physics **Q9** (`directives/S11_SHARED_PHYSICS.md` §Q9): `Quad(D)` of quadratic forms in the `D²` entries of `G`, action `G ↦ R G Rᵀ`, V1 = SO invariants, V2 = O invariants, V6 = `(−1)`-eigenspace of one explicit reflection on V1, `P_D` the sum of the V6 basis (empty ⇒ zero form). Step record `steps/S11_stray_longitudinal.md` Move 2 gives the D=3 counts `N_SO = 3`, `N_O = 3`, with the extras being reflection-odd and absent at D=3.

---

## J1 — object, representation, and full-group classification

**Domain is the full space.** `quadratic_representation` does not start from a candidate list. It takes an arbitrary `Q : Quad`, polarizes through `Q.associated`, expands `G` in the nine matrix units (`frame_expansion`), and writes `Q G` as a linear combination of the 45 monomials `v_i v_j` for `0 ≤ i ≤ j ≤ 8`. That is the associated-bilinear coordinate presentation of every real quadratic form on `Mat`. Over `ℝ`, `QuadraticForm` is the space of homogeneous quadratic maps; constants and linear terms are excluded by definition, which matches `D3_FIDELITY.md`.

**Necessary SO constraints come from admissible members of the full quantifier.** `invariant_polynomial` applies `hQ : SOInvariant Q` to three matrices already proved `Proper`:

- `rotationXY 0 1` and `rotationYZ 0 1` (`rotationXY_proper`, `rotationYZ_proper`)
- `rotationXY (3/5) (4/5)`, the non-permutation rational rotation

Each retained equation is `Q(R G Rᵀ) = Q(G)` at a matrix unit or pair sum, after an explicit `conjugate_apply` image check. `Constraints.lean` contains `e0`–`e41` (42 equations) and reconstructs the 42 non-free coefficients `h0`–`h44` except the free slots `c 4`, `c 11`, `c 9`. Those slots are the monomials `G00 G11`, `G01 G10`, and `G01²`, which extract `a = c4/2`, `b = c11/2`, `c = c9` from `a (tr G)² + b tr(G²) + c tr(GGᵀ)`.

Sampled conjugations match the algebra: `rotationXY 0 1` sends `E_{00}` to `E_{11}` (`e0`); `rotationXY (3/5) (4/5)` sends `E_{00}` to the rank-one matrix with entries `9/25, 12/25, 16/25` (`e41`). The rational rotation is used in the reconstructions of the diagonal-square coefficients (`h0`, `h30`, `h44`), which is where a generic angle is required. Finite tests are used only as **necessary** specializations of `∀ R, Proper R → …`.

**Full O sufficiency is independent of those tests.** `invariantForm_O` is a trace-identity argument on the whole orthogonal group:

- `trace_conjugate`: `tr(R G Rᵀ) = tr G` from `Rᵀ R = I`
- `conjugate_mul`: `R G Rᵀ · R H Rᵀ = R (G H) Rᵀ`
- `conjugate_transpose`: `(R G Rᵀ)ᵀ = R Gᵀ Rᵀ`

so `tr((RGRᵀ)²) = tr(G²)` and `tr((RGRᵀ)(RGRᵀ)ᵀ) = tr(G Gᵀ)` for every orthogonal `R`, including `det = -1`. `SO_classification` then says every SO-invariant form equals `invariantForm v`; `SO_iff_O` / `O_classification` promote that to O-invariance. Finite rotations do not stand in for this sufficiency proof.

**Normalization and unique coefficients are the displayed polynomials.**

- `traceSquare G = (tr G)²`, with the factor `2` on each paired diagonal product
- `traceOfSquare G = tr(G²)`, the trace of the matrix product, with factor `2` on `G01 G10`, `G02 G20`, `G12 G21`
- `frobeniusSquare G = tr(G Gᵀ)`, the sum of all nine squared entries

`invariantForm_injective` separates the three coefficients on `E_{01}`, `E_{00}`, and `E_{00}+E_{11}` (`Q = c`, `a+b+c`, `4a+2b+2c`). `SO_unique` is unique existence of `v`. The mutation `2 • monomial 0 4 ↦ 3 • monomial 0 4` leaves `G00 G11 * 3` versus `* 2`; `2 • monomial 1 3 ↦ -2 • monomial 1 3` leaves a sign error on `G01 G10`. Those are the exact cross coefficients, not sector-name aliases.

---

## J2 — actual submodules, equality, and the zero odd space

`soSpace`, `oSpace`, and `oddSpace` are the sets of forms satisfying the predicates, closed as `Submodule ℝ Quad`.

- `so_eq_o` is `SO_iff_O` on the carriers, so the **submodules** coincide, not merely their dimensions.
- `soEquiv : (Fin 3 → ℝ) ≃ₗ[ℝ] soSpace` from `invariantForm` / `SO_classification` / injectivity gives `so_dimension = 3` and then `o_dimension = 3`.
- `oddSpace` is `{Q | SOInvariant Q ∧ ReflectionOdd Q}` for the same `R₀ = diag(-1,1,1)` used by native Q9. `odd_classification` uses O-invariance of an SO-invariant form plus `Q(R₀ G R₀ᵀ) = -Q(G)` to force `Q = 0`. `odd_eq_bot` and `odd_dimension = 0` identify the odd submodule itself. This is not `dim SO − dim O` from CAS counts.

Zero is included (`zero_mem'`, `zero_invariant`). Arbitrary real coefficients, including negatives, are included (`invariantForm ![-2,3,-5]` in the unique-coefficient and O-invariance positives). `reflection_det` records `det R₀ = -1`, so `R₀` is an improper element; one reflection plus SO-invariance already implies oddness under every improper orthogonal map.

---

## J3 — native Q9 is a span comparison

`_measurements/S11_lean_d3_source_check.py` extracts `compute_q9`, `q9_vector`, `q9_row_to_poly`, `matrix_from_rows`, and `derivative_placeholders` from `scripts/S11_stray_longitudinal_sympy_audit.py` and runs `compute_q9(3)` only. It does not run production drivers or regenerate exports.

Independently, from the recorded bases and polynomials in `_measurements/S11_lean_d3_source_checks.json`:

- Native V1/V2 polynomials are  
  `tr(G²)`,  
  `((tr G)² − tr(G²))/2`,  
  `‖G‖_F² − tr(G²)`.
- These are an invertible change of basis from Lean’s three forms, so the **spans coincide**. Matching dimensions would not have been enough; the recorded `same_rref` / `stacked_rank` checks are actual span tests.
- The 45-vector rows match those polynomials in the pinned row-major `combinations_with_replacement` monomial order used by both Lean and Q9.
- V6 is the empty `0×45` matrix; `PD_POLY = 0`; `V6_OPERATOR` is recorded as the 3×3 identity, i.e. every V1 form is even under `R₀`.

Controls in the instrument:

- `{ (tr G)², tr(G²), G00² }` has rank 3 but is not the invariant span.
- Restoring the old generator orientation `lie_equations.extend(action_rows)` (dropping `.T`) is required to keep counts `3/3/0` while failing the span test.

Native V1 is computed from `so(3)` generators with the **transpose** of the monomial action; Wolfram has the matching `Transpose[actionRows]` at line 910. Shared physics asks for invariance under every group element; the engines use the Lie-algebra route. Lean proves the group statement. J3 only identifies the resulting subspaces. That is the specified compact link, and it is documented as tested translation, not a kernel-certified CAS bridge.

---

## J4 local controls

`Controls.lean` supplies admissible positives and the corresponding omissions:

- `nonzero_invariant_exists`, `zero_invariant`, `negative_coefficients_positive`, `unique_coefficients_positive`
- `trace_not_omittable`, `traceOfSquare_not_omittable`, `frobenius_not_omittable`
- `nonzero_odd_impossible`, `single_entry_not_invariant` (G00² fails an admissible 90° rotation)

Author-recorded contract run (`_measurements/S11_lean_d3_contract_checks.json`, `S11_lean_d3_validation.json`, `D3_VERIFICATION.txt`): 40 axiom audits on standard axioms only; ten mathematical rejections (counts 2/4/1, three omissions, odd extra, single-entry invariance, and the two coefficient identities) with `⊢ False` or unsolved identities in the named declarations; twelve positives. Syntax/import/timeout failures are excluded by the instrument. I did not rerun that suite.

---

## Blockers versus optional improvements

**Blockers / required fixes:** none for this bounded theorem, coverage, or fidelity.

**Optional (not required to clear J1–J4):**

- The `so_dimension_positive` / `o_dimension_positive` / `odd_dimension_positive` wrappers only restated the proved equalities. Nonvacuity is already carried by the existence, zero, negative-coefficient, uniqueness, and omission theorems.
- Uniqueness of the 45-coefficient polynomial representative is true on `ℝ^9` but not named as its own lemma; uniqueness of `v` is what the contract asks for.
- Native Q9 V1 is Lie-algebraic; Lean’s theorem is the group statement. The packet already treats the native link as span identification.

Out of scope, as instructed: D4/D5, dynamics, EL/divergence classes, interface physics, S11c, production/comparator/export, systematic CAS-bridge work.

---

## Limits of this independent verification

- **Read** `MANIFEST.json`, `lean/FORMALIZATION_POLICY.md`, `D3_COVERAGE.md`, `D3_FIDELITY.md`, `D3_VERIFICATION.txt`, all six `S11D3Invariants` modules, the audit root, generator, source/contract instruments and JSON records, native SymPy Q9, Wolfram `buildInvariantCensus`, shared physics Q9, and the step record.
- **Did not** compile Lean, rerun `_measurements/S11_lean_d3_source_check.py` or `_measurements/S11_lean_d3_contract_check.py`, or recompute SHA-256 / the aggregate digest. Hash agreement cited above is cross-file consistency of the supplied records, not a new hash.
- Independently checked: domain and polar representation; SO vs O quantifiers; trace-identity O sufficiency; unique `(a,b,c)` normalization; submodule census `3/3/0` with `oddSpace = ⊥`; sampled constraint images and reconstructed targets against `a (tr G)² + b tr(G²) + c tr(GGᵀ)`; native V1/V2 polynomials as an invertible re-basis of those three forms; generator `.T` / Wolfram `Transpose[actionRows]`; empty V6 / zero `P_D`; presence of same-count wrong-span and wrong-orientation controls in the instrument source.
- Lean `linear_combination` certificates and the author-reported kernel/axiom/mutation PASS remain **author-recorded**. The mathematics of the statements those certificates are supposed to prove is consistent.

**CLEAR.**
