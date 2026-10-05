**CLEAR**

Packet: `S11 D5 quadratic invariant D5.1–D5.4 fidelity contract v1` (MANIFEST as supplied; aggregate SHA-256 `ff8aa6e60d1403128615212ccbf601fc0c038ef82d0f64690144d9a9457c081c` not recomputed).

This review is statement fidelity of the declared D5.1–D5.4 increment only. No required fix to those claims was found. Compilation, hash recomputation, unpackaged objects, and native/Lean re-execution are author evidence and were not independently performed.

## Checked assumptions and conventions

- Object: `Quad := QuadraticForm ℝ Mat` on **all** real 5×5 matrices (`CoordinateAlgebra.lean`). No symmetric-gradient, transverse, positivity, field-equation, or nonzero restriction.
- Coordinates: 25 **row-major** entries `G 0 0, G 0 1, …, G 4 4`. Native `QG_ALL` is `g_1,…,g_25` in the same order; `coordinate_substitution` sets `g[i,j] = ∂_{x_i} u_j`, i.e. `G_ij = ∂_i u_j`.
- Action: `conjugate R G := R * G * R.transpose` (`Rotation.lean`). Native Q9 uses the matching infinitesimal `XG − GX` and finite `R G Rᵀ`.
- Claimed conclusion: unique `Q(G) = a (tr G)² + b tr(G²) + c tr(G Gᵀ)` for arbitrary real `a,b,c`, including zero and negative. Pointwise **density** classification before any Euler–Lagrange or total-divergence quotient. Dimension five is supplied; no kinetic action, unit conversion, first-variation, or bulk-response theorem.

## D5.1 — complete density space

`quadratic_representation` (`Quadratic.lean`) expands an **arbitrary** `QuadraticForm` via its associated bilinear form on the 25 frame matrices. `polynomial` / `bilinearCoefficients` use all **325** unordered pairs `i ≤ j`. Squares keep `B(e_i,e_i)`; mixed terms keep `B(e_i,e_j)+B(e_j,e_i)` once. That is a basis of homogeneous quadratics, not a three-parameter ansatz.

The bilinear identity is split into 24 mixed rows (`bilinear_cross0`–`23`) and 25 tails (`bilinear_tail24` down to public `bilinear_tail0`), then assembled with Mathlib `List.foldl1_eq_foldr1` (`BilinearExpansion.lean`, `Quadratic.lean`).

Necessary equations come from four adjacent quarter-turns and `rotationXY (3/5) (4/5)` (`S11_lean_d5_generate.py`, `ConstraintBlock16.lean`). Lean derives each selected equation from quantified `OInvariant` at an explicit matrix (`ConstraintBlock*.lean`), not from the generator’s rank assertion. Seventeen reconstruction modules target the 322 non-free coefficients with exact rational `linear_combination`s. Free slots are **6, 29, 25**:

- `c 6` = `G00 G11` → `a = c6/2`
- `c 29` = `G01 G10` → `b = c29/2`
- `c 25` = `G01²` → `c = c25`

Sampled targets match the three-form coefficients (`c0 = a+b+c`, `c12 = 2a`, `c172 = c`, `c180 = 2b`, `c323 = 0`, `c324 = a+b+c`). `Constraints.invariant_polynomial` substitutes all reconstructed coefficients and `ring`s to the three traces. `Forms.invariantForm_O` proves full-group sufficiency from `trace_conjugate`, `conjugate_mul`, `conjugate_transpose`; the five rotations are not the definition of O-invariance.

`SO_classification` / `O_classification` / `invariantForm_injective` / `SO_unique` (`Classification.lean`) give existence and uniqueness on the whole stated domain. Injectivity uses three probes (off-diagonal, one diagonal, two diagonals) that separate `c`, `a+b+c`, and `4a+2b+2c`.

## D5.2 — parity and census

`conjugate_neg` and `orthogonal_det` plus `Matrix.det_neg` in dimension five give `SO_iff_O`: equality of **invariance notions**, not the three-dimensional classification by itself. `ReflectionOdd` is oddness under `diag(-1,1,1,1,1)` with `reflection_det : det = -1`. `odd_classification`: SO invariance plus this oddness iff `Q = 0`. Zero remains admissible (`zero_invariant`); a nonzero odd invariant is impossible (`nonzero_odd_impossible`).

`so_eq_o`, `odd_eq_bot`, `soEquiv : (Fin 3 → ℝ) ≃ₗ[ℝ] soSpace`, and `census = (3,3,0)` (`Census.lean`) are dimensions of **density subspaces**, not bulk-response counts. The three forms are separately non-omittable (`Controls.lean`).

## D5.3 — compact native identification

`_measurements/S11_lean_d5_source_check.py` AST-extracts `compute_q9`, `derivative_placeholders`, `matrix_from_rows`, `q9_row_to_poly`, `q9_vector` only; it does not import the production module or run its driver. It compares full V1/V2/V6 coefficient spans to `(tr G)²`, `tr(G²)`, `tr(G Gᵀ)` in the native unordered monomial order.

The recorded native polynomials are a different spanning set of the **same** three-space (`tr(G²)`, `((tr G)² − tr(G²))/2`, `‖G‖_F² − tr(G²)`). JSON reports 3/3/0, `V6_OPERATOR = I`, `P_D = 0`, both controls detected: `G00²` same-count/wrong-span, and in-memory drop of `.T` retaining 3/3/0 but failing span. Engine files are hash-compared unchanged. Wolfram `Transpose[actionRows]` (`mathematica/S11_stray_longitudinal_mathematica_audit.wl` ~line 910) is source inspection only. This translation is **outside** Lean’s kernel. No V5, native EL, production export, or comparator claim.

## D5.4 — controls

Thirteen paired false statements match the coverage table (SO/O/odd counts 4/2/1, nonzero odd existence, each omitted basis direction, `G00²` invariance, normalizations 2/1/2 vs 4/2/1, reflection det `+1`, `soSpace ≠ oSpace`). Instrument requires one `contract_control` diagnostic, one unsolved `False`, no warning/import/resource failure; true counterparts compile. Four extra positives: zero, nonzero existence, arbitrary coefficients, unique coefficients (17 positives). `MUTATIONS` is empty: paired statements only. `ne_self_iff_false` is only on the SO/O **control** proof, as documented for run11. These are sensitivity controls, not exhaustive mistranslation detection.

Author evidence (not re-run here): 45 canonical objects, 51 axiom audits restricted to `propext` / `Classical.choice` / `Quot.sound` or empty, 30 fresh controls, 75 records. D5 modules import only `S11D5Invariants.*` and Mathlib; no older local Lean proof.

## Optional stronger results (not required)

General dimension, EL/divergence quotient, D5 bulk responses, kernel-certified CAS, canonical source-replacement mutations, V5/native EL, production/S11c, interface/spectrum/scattering/pole work.

## Verification limits

- No Lean build, axiom print, or control execution in this review.
- No SHA-256 recomputation of packet, objects, or unpackaged Mathlib/olean/historical files.
- Generator rank / `gauss_jordan_solve` were not re-executed; Lean’s quantified equations and `linear_combination`s were inspected, including sampled coefficient targets. Every `ConstraintBlock0`–`16` is imported by some reconstruction module; `e321` (rational plane rotation) is used. All 322 `e_i` occurrences in `linear_combination` lines were not exhaustively listed by hand.
- Native `compute_q9(5)` was not re-run; span identification used the instrument source plus the polynomials recorded in `S11_lean_d5_source_checks.json`.
- Mathlib pins and compiled objects are recorded, not unpacked here.

No required D5.1–D5.4 statement mismatch was found. Independent compilation/hash clearance is not claimed.
