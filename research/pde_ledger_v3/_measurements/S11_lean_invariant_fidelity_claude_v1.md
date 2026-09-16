# Verdict: CLEAR for S11 D2 invariant contract I1–I4

**Manifest aggregate hash:** `f95d1667e54b4c8139881feaa86267ebe60eb4e9d85f0f4f4bd2dff9dc2b2747`, as declared in `MANIFEST.json:6` (revision "S11 D2 invariants I1–I4 fidelity contract v1").

I found no unresolved blocking finding for I1–I4.

**What I could not check:** I had no shell in this session. That leaves two gaps:
- **Hashes:** I could not recompute the per-file SHA-256 values or the aggregate hash.
- **Builds:** I did not rebuild the Lean files. The build and axiom results below come from the recorded transcripts.

Otherwise, I compared the file set with the manifest (27 files plus `MANIFEST.json`). I cross-checked every hash that the terminal records cite against the manifest, and all of them agree. I also checked the load-bearing mathematics by hand.

## Findings by check

**I1: the object, the group action and the orientation**
- **All forms, not an assumed basis.** `Quad := QuadraticForm ℝ Mat` (`Quadratic.lean:12`). `quadratic_representation` derives the 10 coefficients from `Q.associated`, and `polynomialForm_injective`/`_surjective` make the presentation exhaustive and unique.
- **Coordinates.** t=G11+G22, s=G12−G21, x=G11−G22, y=G12+G21, with the correct inverse (`decode`, both round-trip theorems). No normalization factors are involved, so no span is affected.
- **Row-major convention.** `G 0 1`=g₂ and `G 1 0`=g₃. This matches Q9 in `S11_SHARED_PHYSICS.md:791`, SymPy `QG_ALL` (line 104) with `qg[i,j]=variables[i*n+j]`, and Wolfram's `Partition`.
- **Quantifiers.** `SOInvariant`/`OInvariant` range over every R with RᵀR=1 (and det R=1 for SO) (`Rotation.lean:15-16`).
  - `proper_rotation` shows every proper R has the form `rotation a b`.
  - `orthogonal_invariance_iff` handles det=−1 by composing with the reflection.
  - No sampled group is used.
- **Generator sign and orientation.**
  - By hand, the native X=[[0,1],[−1,0]] with δG=XG−GX gives δt=δs=0, δx=2y, δy=−2x.
  - That is the negative of the derivative of `rotation(cos θ, sin θ)` at θ=0. `INVARIANT_FIDELITY.md:42-45` states this correctly, and the sign does not change the zero constraint.
  - `coefficient_action` is the correct general identity Σᵢcᵢ(Am)ᵢ = Σⱼ(Aᵀc)ⱼmⱼ.

**I2: complete spaces and dimensions**
- **Necessary constraints only from selected rotations.** `invariant_coefficients` uses the rotations (0,1) and (3/5,4/5) only to derive necessary conditions. I checked the arithmetic by hand:
  - the quarter turn kills the coefficients of tx, ty, sx and sy;
  - the two rational evaluations give c₇−c₉ = −(7/24)c₈ and c₇−c₉ = (24/7)c₈, so c₈=0 and c₇=c₉.
- **Sufficiency for every rotation.** `invariantForm_SO` proves invariance for every proper R.
- **Dimensions and independence.** `SO_classification`, `invariantForm_injective` and the three `LinearEquiv`s give finrank 4/3/1.
- **Odd space restricted to SO.** `oddSpace` is `SOInvariant ∧ ReflectionOdd` (`Census.lean:33`), which matches Q9's "(−1)-eigenspace within V1".
- **Unique even/odd decomposition.** `even_odd_disjoint` together with `even_odd_span` gives it, including the zero form.
- **No quotient claims.** Equality is equality of polynomials, and no Euler–Lagrange or divergence quotient is claimed anywhere.

**I3: native identification**
- **Minimal SymPy fix.** The patch transposes each generator's block before stacking (`sympy_audit.py:511-517`). This is minimal and correct:
  - Transposing the stacked non-square matrix would be wrong; transposing each block is right.
  - The reflection block is diagonal for the diagonal `r0`, so leaving it untransposed is harmless.
  - V6 correctly builds the operator with coordinates of F(bᵢ) as columns.
  - No expected basis or count is inserted into the engine.
- **Wolfram already correct.** It applies `Transpose[actionRows]` inside the per-generator `Table` (`.wl:899-911`), and its `quadraticCoordinates` uses the same unnormalized cross-term convention as SymPy.
- **Native spans checked by hand.** The recorded V1 rows equal the contract span:
  - g₁²+2g₂g₃+g₄² = (t²+x²+y²−s²)/2
  - g₁g₄−g₂g₃ = det
  - (g₂−g₃)² = s²
  - (g₁+g₄)(g₂−g₃) = ts

  O = {first, det, s²} and odd = {ts}. The source check compares exact RREF matrices, which is a complete span comparison, not a count match.
- **Counterexample.** The pre-repair emitted form g₁²+g₂g₃+g₄² (probe JSON) equals Lean's `wrongNativeForm`: coefficients (½, −¼, ½, ¼) on (t², s², x², y²).
  - For R=[[3/5,−4/5],[4/5,3/5]] and G=diag(1,0), RGRᵀ=[[9,12],[12,16]]/25, so p goes from 1 to (81+144+256)/625 = 481/625.
  - `wrongNativeForm_not_SO` derives the contradiction ½=¼ from c₇=c₉.
- **Scope of the native check.** The span check covers D2 V1/V2/V6 through five extracted helpers only. The Wolfram link is source inspection only, and the documents say both of these accurately.

**I4: controls, mutations and hashes**
- **Mutations rejected for mathematical reasons.**
  - The orientation mutation fails with the unsolved goal `c i * A i j * m j = c i * m j * A j i`.
  - The odd-pairing mutation leaves the t·s term unmatched.
  - The six negated-claim controls (wrong counts 3/4/0, odd pairing omitted or made O-invariant, wrong native form accepted, rotation(1,1) accepted as proper) each fail with `⊢ False`.
  - The acceptance logic excludes unknown identifier/module, timeout and memory failures.
- **Positive controls.** All eight passed, including the admissible 481/625 witness.
- **Axiom audit.** All 46 audited theorems appear, which I cross-checked against the declarations, and each uses only the standard axioms.
- **Hash chain.** The hashes of the Lean sources, lake files, instruments and reports agree between `contract_checks.json`, `source_checks.json`, `validation.json`, `INVARIANT_VERIFICATION.txt` and `MANIFEST.json`. The probe's instrument hash matches its manifest entry. Its source hash `352fc5…` is the pre-repair file, which is not in the packet, so I could not check it.

**Boundary and H1–H4 separation**
- I found no actual inconsistency:
  - `COVERAGE.md:50` excludes invariant completeness from H1–H4.
  - The new Lake target is disclosed, and the H1–H4 hashes are not relabeled as current.
  - The README and the fidelity record disclaim D3–D5, V5 EL, PD production, comparator/export and S11c.
- **Observation (not a clearance):** at D2 the pre- and post-repair V6 rows are identical (span{ts}), so the D2 PD density does not change. V1, V2, V4 and V5 emissions do change. The effect at D3–D5 is unmeasured.

## Blockers
None.

## Optional suggestions (not closure conditions)
1. **`S11_lean_invariant_source_check.py:74-76`:** the −144/625 residual is computed from a hardcoded p, not taken from the in-memory mutant's V1 output. Yet `orientation_control` records it next to the mutant counts. Either assert that p lies in the RREF span of `wrong['V1_BASIS']`, or label the residual as independent.
2. **Orientation mutant:** it also fails on an unused-`simp`-argument linter error. The acceptance criterion correctly keys on the unsolved goal, but dropping `Matrix.transpose_apply` from the mutant (or disabling that linter there) would isolate the mathematical failure.
3. **`odd_pairing_coefficient` mutation:** it is caught at the `invariantForm_apply` definition/lemma mismatch. The substantive "omitting the odd pairing" claim is carried by `odd_pairing_is_required`; saying so in the verification text would help.
4. **Wolfram link:** it is a substring check (`'Transpose[actionRows]' in wl`). Anchoring it to the per-generator `Table` context would make it harder to satisfy by accident.
5. **Status wording:** `INVARIANT_NEXT.md:94` still says packet transfer permission is pending. Update it when the review dispositions are recorded.
