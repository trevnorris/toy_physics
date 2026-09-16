**Verdict: CLEAR**

Reviewed packet (as supplied, hashes not recomputed): `MANIFEST.json` revision **S11 D4 quadratic invariants D4.1–D4.4 fidelity contract v1**, `aggregate_sha256` `a91cfbf327a6d180e5f8d3e9e3dd9dc0037c426b71be41c0a55ec45d7a719762`, created 2026-09-16T21:53:57.178448+00:00.

This is a statement-fidelity review of D4.1–D4.4 only. I did not compile Lean, rerun the native instrument or generator, or recompute hashes.

---

## Scope and assumptions checked

The claimed object is every real quadratic density on every real 4×4 matrix, under conjugation `G ↦ R G Rᵀ`, with `G_ij = ∂_i u_j` and independent row-major entries. No symmetric-gradient restriction, positivity, EL or divergence quotient, field equation, regularity, boundary condition, transverse ansatz, or derived physical dimension is used.

Checked against `lean/FORMALIZATION_POLICY.md`, `lean/s11/D4_COVERAGE.md`, `D4_FIDELITY.md`, the thirteen `S11D4Invariants` modules plus audit root, `_measurements/S11_lean_d4_generate.py`, the compact native check, `directives/S11_SHARED_PHYSICS.md` §7/Q9, `steps/S11_stray_longitudinal.md` Move 2, and the native Q9 helpers.

`SOInvariant` / `OInvariant` quantify over every real `R` with `RᵀR = I`, plus `det R = 1` for SO. Coefficients `(a,b,c,d)` range over all reals, including zero and negatives. The next increment (D4 odd term as divergence / zero bulk variation) is correctly excluded.

---

## D4.1 — object, 136-coefficient representation, necessity, sufficiency

**Associated-bilinear representation covers every real quadratic form on all 16 entries.**  
`Quadratic.lean` sets `Mat := Matrix (Fin 4) (Fin 4) ℝ` and `Quad := QuadraticForm ℝ Mat`. `frame_expansion` writes every matrix in the 16 matrix units. `quadratic_representation` polarizes `Q.associated` and encodes the 136 monomials `v_i v_j` for `i ≤ j` (`C(16+1,2) = 136`). Off-diagonal slots are `B(e_i,e_j)+B(e_j,e_i)`, so the expansion of `B(∑ v_k e_k, ∑ v_ℓ e_ℓ)` is complete over `ℝ` (2 invertible). This is the space of homogeneous quadratic polynomials on 16 independent entries, not a candidate list.

**All 132 necessary equations come from admissible members of the full SO quantifier.**  
Generator rotations, each proved `Proper` in `Rotation.lean`:

| Test | Lean | `a²+b²=1` |
|---|---|---|
| XY quarter-turn | `rotationXY 0 1` | `rotationXY_proper` |
| YZ quarter-turn | `rotationYZ 0 1` | `rotationYZ_proper` |
| ZW quarter-turn | `rotationZW 0 1` | `rotationZW_proper` |
| XY 3-4-5 | `rotationXY (3/5) (4/5)` | `rotationXY_proper` |

Each block goal is `have raw := hQ (rotation…) (…_proper (by norm_num)) (decode …)`: an instance of `∀ R, Proper R → ∀ G, Q(RGRᵀ)=QG`. I spot-checked the first Block0 equation: XY 90° sends `E_{00}` to `E_{11}`, so `c 0 = c 70`, matching `(-1)*c 0 + (1)*c 70 = 0`. Block3 ends with the 3-4-5 rotation; the last displayed equation has denominators 625 = 5⁴.

**Reconstruction uses those 132 equations to prove completeness.**  
Four blocks of 33 conjuncts unpack consecutively as `e0`–`e131` in `invariant_polynomial` (`Constraints.lean`). The four free slots are exactly the pairs that isolate the four forms: `(0,5)` → `a = c5/2` (cross term of `(tr G)²`), `(1,4)` → `b = c19/2` (`tr(G²)`), `(1,1)` → `c = c16` (Frobenius), `(1,11)` → `d = c26` (`G01 G23` in `P`). The other 132 coefficients are rewritten from the `e_i` by `linear_combination`, then `ring` matches `a(tr G)² + b tr(G²) + c tr(GGᵀ) + d P`. `SO_classification` packages this as `SOInvariant Q ↔ ∃ v, Q = invariantForm v`.

**Block split introduces no extra hypotheses and drops no equation.**  
Each `necessaryBlock*` has the same `(hQ : SOInvariant Q)` and the same `hc` representation. `ConstraintPolynomial.polynomial_vec` is `rfl`. The generator comment and module docs state the split is compilation-only. Reconstruction consumes all four conjunctions.

**Full-group sufficiency is independent of the finite tests.**  
`invariantForm_SO` uses `trace_conjugate`, `conjugate_mul`, `conjugate_transpose`, and `orientation_conjugate` for every proper `R`. `orientation_conjugate` is stated for **every** real 4×4 `R`, with no invertibility or orthogonality assumption:

```
orientation (conjugate R G) = R.det * orientation G
```

That is `P(RGRᵀ) = det(R) P(G)`. The finite rotations are necessity certificates; they do not stand in for this identity.

Finite tests plus reconstruction already force the four-parameter family, and that family is SO-invariant, so the two legs close.

---

## D4.2 — normalization, O / odd split, dimensions, even/odd decomposition

**Displayed forms match the named polynomials, with unique coefficients.**

| Form | Theorem | Normalization checked |
|---|---|---|
| `(tr G)²` | `traceSquare_apply` | diagonal squares plus **2** on each cross term |
| `tr(G²)` | `traceOfSquare_apply` | diagonal squares plus **2** on `G_ij G_ji` |
| `tr(G Gᵀ)` | `frobeniusSquare_apply` | sum of all 16 entry squares |
| `P` | `orientationForm_apply` | twelve monomials of the displayed Pfaffian, **+1 on `G01 G23`** |

I expanded `P = (G01-G10)(G23-G32) - (G02-G20)(G13-G31) + (G03-G30)(G12-G21)` against `orientationForm`; the twelve terms match. `invariantForm_injective` / `SO_unique` recover `(a,b,c,d)` from four matrices (`E_{01}`, `E_{00}`, `diag(1,1,0,0)`, and `G01=G23=1`). Those probes separate Frobenius, `a+b+c`, `2a+b`, and `2c+d`.

**O-invariants are exactly `d=0`.** `invariantForm_O` / `O_classification`: one reflection on the `P`-detecting matrix forces `d=0`; the three trace forms are invariant under every orthogonal `R`.

**Reflection-minus forms inside SO are exactly `a=b=c=0`.** `invariantForm_odd` / `odd_classification`. `ReflectionOdd` is the minus space of **this** reflection `diag(-1,1,1,1)` on SO-invariants, as `D4_FIDELITY.md` states. For SO-invariants that is the span of `P`, which is odd under every orientation-reversing conjugation via `det`.

**Dimensions 4/3/1 are `finrank` of the actual submodules** `soSpace`, `oSpace`, `oddSpace` of `Quad`, via linear equivalences `soMap` / `oMap` / `oddMap` (`Census.lean`), not counts of a supplied list. This matches the step record’s corrected census `N_SO(D=4)=4`, `N_O=3`, one reflection-odd extra.

**Even/odd span and disjointness, including zero.**  
`even_odd_disjoint`: O-invariance plus reflection-odd ⇒ `Q(G)=0`. Submodule disjointness is intersection `{0}`. `even_odd_span`: every SO form splits as `(a,b,c,0)+(0,0,0,d)`. `zero_invariant` puts the zero form in both; that is the overlap. Together the splitting is unique.

---

## D4.3 — native Q9 identification (exact spans, not counts)

Shared physics defines `P_D` as the **sum of emitted `V6_BASIS`**, not an epsilon formula, and forbids rescaling (`directives/S11_SHARED_PHYSICS.md` §7). Native `compute_q9` does that: `PD_POLY` is the sum of V6 rows in the emitted monomial order. Compact check (`S11_lean_d4_source_check.py` / `.json`):

- V1 / V2 / V6 vs spans of the four forms / first three / `P`: same rank, stacked rank, **and identical RREF**. Native V1 polynomials are a different basis of the same 4-space; V6’s polynomial is **string-identical** to displayed `P`.
- `PD_equals_contract_P`: `expand(PD_POLY − P) = 0`. Coefficient of `G01 G23` in the 1-dimensional RREF is +1 (pair index 26).
- Fully summed `ε_ijkl G_ij G_kl` with `ε_0123=+1` is asserted `2P`. I checked the probe `G01=G23=1` (`P=1`): permutations `(0,1,2,3)` and `(2,3,0,1)` contribute +1 each, total 2. The factor two is recorded, not divided out.
- Reflection: native `R0 = diag(-1,1,1,1)` matches Lean; `V6_OPERATOR` in the native V1 RREF chart reconstructs the reflected V1 rows.
- Same-count/wrong-span: replace `P` by `G00²`, rank stays 4, stacked rank exceeds 4.
- Historical generator orientation: in-memory drop of `.T` on `matrix_from_rows(action_rows, …)` keeps counts `[4,3,1]` and **fails** span. Native source retains  
  `lie_equations.extend(matrix_from_rows(action_rows, len(monomials)).T.tolist())`.
- Wolfram: source inspection of `Transpose[actionRows]` in `mathematica/S11_stray_longitudinal_mathematica_audit.wl` (around the Q9 census); no engine run is claimed.

The Python link is AST-extracted helper execution (`compute_q9`, `q9_vector`, `q9_row_to_poly`, `matrix_from_rows`, `derivative_placeholders`) with reconstructed `QG_ALL = g_1,…,g_25`. Native `QG_ALL` uses the same names; `G_ij = ∂_i u_j` is `coordinate_substitution`. This is tested translation, not a kernel-certified CAS bridge, and the records say so.

Native V1 is the **Lie-algebra** kernel of `so(4)` on `Quad(4)`, while Lean quantifies over finite `SO(4)`. For this compact connected linear action those spaces coincide; the compact check then compares the computed native spans to the Lean forms. That agreement is the identification, not a second Lean proof of the Lie-algebra route.

---

## D4.4 — axioms, mutations, positives, full-group sensitivity

Audit root prints 49 load-bearing declarations. The recorded axiom set is only `propext`, `Classical.choice`, `Quot.sound`. The checker also forbids `axiom` / `sorry` / `admit` in the D4 modules.

Fourteen mathematical rejections and sixteen positives are recorded. Timeouts raise and are tagged `accepted_as_mutation: False`.

| Control | What it tests | Diagnostic |
|---|---|---|
| `trace_normalization` | `2 • monomial 0 5` → `3 • …` | unsolved identity: cross coefficient 3 vs 2 in `traceSquare_apply` |
| `trace_of_square_sign` | `2 • monomial 1 4` → `-2 • …` | unsolved identity: `−2 G01 G10` vs `+2` |
| dimensions 3/4/0 | census | `⊢ False` |
| omit each of four forms | completeness | `⊢ False` |
| odd claimed O-invariant; even claimed odd | reflection split | `⊢ False` |
| `full_group_sensitivity` | `G00²` claimed SO-invariant | `⊢ False`; XY 90° sends `E00` to `E11` |
| `orientation_normalization` | `P=2` vs `P=1` | `⊢ False` |
| `determinant_sign` | reflected `P = +1` vs `−1` | `⊢ False`; uses `orientation_conjugate` then `reflection_det` |

Positives include the twelve paired true statements, nonempty odd space, zero form, `invariantForm ![-2,3,-5,-7]`, and unique coefficients. These show admissibility and coefficient signs; they do not replace the classification.

---

## Blockers vs optional

**Required fix for this bounded theorem / coverage / fidelity: none.**

Optional, not blocking:

- A conceptual Pfaffian proof of `orientation_conjugate` instead of the `det_four` + `ring` expansion.
- Native V1 is infinitesimal; the finite-group match is the compact span check, already labeled as Python/SymPy evidence.

Out of scope, as required: D4 divergence / zero bulk variation, D5, EL, dynamics, interfaces, S11c, production/comparator/export, systematic CAS bridging.

---

## Limits of this independent verification

I read the listed Lean modules, generator, instruments, native Q9 helpers, shared-physics Q9/`P_D` rule, and terminal JSON. I checked the 136-dimensional encoding, P monomial expansion, injectivity probes, first necessity equation, free-parameter map, 4×33 block unpacking, sufficiency/orientation identity, O/odd/census/even-odd lemmas, mutation source and `⊢ False` / unsolved-goal diagnostics, native RREF/span/reflection/`P_D`/`2P` logic, generator transpose, and Wolfram `Transpose[actionRows]`.

I did **not** compile Lean, so I did not kernel-check `orientation_conjugate`’s large `ring` goal, the 132 `linear_combination` proofs, or `det_four` against the full Leibniz expansion (I counted 24 terms and checked `I` and the coordinate reflection). I did **not** execute `compute_q9(4)`, the generator `--check`, or hash aggregation. Author-reported run4/run5 builds and the 44-check PASS are treated as recorded verification, not as compilation I performed.

On that reading, D4.1–D4.4 identify the intended Q9 D4 quadratic densities, prove the full-group classification with unique coefficients and census 4/3/1, and connect to native V1/V2/V6 by exact spans and a checked `P_D = P` (with `ε(G,G)=2P` recorded separately). **CLEAR.**
