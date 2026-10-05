# Independent fidelity review: S11 D4 quadratic invariants D4.1–D4.4

## Verdict: **CLEAR**

I found no blocking problem with the mathematics, the statement fidelity, coverage or the native identification. The only obligation this review does not close is the second independent review that the policy requires. My clearance covers the mathematics and statement fidelity; the successful Lean builds are reported by the author and I did not reproduce them (see "Limits" below).

**Packet reviewed (as supplied, hashes not recomputed):**
- `MANIFEST.json`, revision "S11 D4 quadratic invariants D4.1–D4.4 fidelity contract v1".
- Aggregate SHA256 `a91cfbf327a6d180e5f8d3e9e3dd9dc0037c426b71be41c0a55ec45d7a719762`.
- The hashes recorded in the validation, contract and source records match the corresponding `MANIFEST.json` entries as text. This includes both instruments, both reports, the generator and all 14 Lean sources.

## Scope and assumptions checked

- **Object:** `Mat = Matrix (Fin 4) (Fin 4) ℝ` with 16 independent entries, and `Quad = QuadraticForm ℝ Mat`.
  - No symmetry, positivity, field-equation, regularity, boundary or transverse assumption appears anywhere.
  - Constants and linear terms are excluded because `Quad` contains only quadratic forms, and the documents say so.
- **Groups:**
  - `Orthogonal R := Rᵀ R = 1` and `Proper := Orthogonal ∧ det R = 1` (`Rotation.lean:11-12`). For square real matrices this is exactly O(4) and SO(4).
  - The action is `conjugate R G = R * G * Rᵀ`.
  - These match the spec's Q9 definitions of V1 and V2 (`S11_SHARED_PHYSICS.md:798-799`).
- **Reflection:** `diag(-1,1,1,1)`, the same matrix as the native `r0` (`sympy_audit.py:525`).
- **Coordinates:** row-major, pairs `(i ≤ j)` in lexicographic order, the same as the spec's pinned `MONOMIAL_ORDERING` (`SHARED_PHYSICS.md:791-793`). `G_ij = ∂_i u_j` matches spec line 35.

## Findings

### D4.1: every quadratic form is covered; necessity and sufficiency

1. **Representation covers every quadratic form.**
   - `quadratic_representation` (`Quadratic.lean:214`) takes an arbitrary `Q`, uses `Q G = Q.associated G G` and expands `G` in the 16-matrix frame (`frame_expansion`).
   - Coefficients are `b(Fᵢ,Fᵢ)` on the diagonal and `b(Fᵢ,Fⱼ)+b(Fⱼ,Fᵢ)` off it, which is correct. All 136 upper-triangular monomials of `polynomial` are present.
   - The only hypothesis is `Q`; no candidate family is assumed.
2. **All 132 necessary equations come from admissible members of the full SO quantifier.**
   - Every one of the 132 `have raw := hQ (…)` lines uses one of four rotations: `rotationXY 0 1`, `rotationYZ 0 1`, `rotationZW 0 1` or `rotationXY (3/5) (4/5)`.
   - Each is paired with its own `rotation*_proper (by norm_num)` proof, which discharges a²+b²=1. A mismatched lemma would not typecheck.
   - `rotation*_proper` proves both orthogonality and `det = 1` (via `det_four`, itself proved from `Matrix.det_succ_row_zero`).
   - I checked by hand, against `coordinates_rotationXY`:
     - Block 0, conjuncts 1, 5, 31 and 32.
     - The final Block 3 equation (the 3/5, 4/5 rotation at `E₀₀`; for example the c₀ coefficient is (81−625)/625 = −544/625).
3. **The reconstruction uses the equations, and the block split loses nothing.**
   - `invariant_polynomial` (`Constraints.lean:13`) destructures exactly 33+33+33+33 = 132 named equations (`e0…e131`).
   - It proves 132 coefficient identities `h_k`, one for every index except the free ones, 5, 16, 19 and 26.
   - It substitutes all of them and closes with `ring`.
   - Each `h_k` is a `linear_combination` certificate, so the kernel re-checks the rational identity; the generator's rank assertion is never trusted. By hand, I confirmed `h1`, `h31`, `h35` and `h58`.
   - The blocks are plain conjunctions with the same hypotheses `(hQ, c, hc)`. Splitting added no assumption and dropped no equation.
4. **Free coefficients are identified correctly.** a = c₅/2 (the G₀₀G₁₁ cross term of (tr G)²), b = c₁₉/2 (G₀₁G₁₀ in tr G²), c = c₁₆ (G₀₁² in tr GGᵀ) and d = c₂₆ (G₀₁G₂₃ in P). The value `h0: c0 = a+b+c` is consistent.
5. **Full-group sufficiency is proved separately from the finite tests.**
   - `orientation_conjugate` (`Orientation.lean:15`) proves `orientation (R G Rᵀ) = det R · orientation G` for every real R, with no orthogonality or invertibility hypothesis. This is the identity Pf(R A Rᵀ) = det R · Pf(A) for A = G − Gᵀ, and P is exactly Pf(G−Gᵀ) with orientation (0,1,2,3).
   - `invariantForm_SO` combines it with `trace_conjugate`, `conjugate_mul` and `conjugate_transpose`, all under the hypothesis `Orthogonal R`, and `det R = 1`.
   - `SO_classification` is an iff; `SO_unique` gives unique coefficients via `invariantForm_injective`, which evaluates at four explicit matrices.

### D4.2: reflection split and dimensions

6. **Normalizations of the four forms.**
   - `traceSquare_apply`, `traceOfSquare_apply`, `frobeniusSquare_apply` and `orientationForm_apply` tie the monomial definitions to (tr G)², tr(G²), tr(GGᵀ) and P.
   - I checked `orientationForm`'s 12 signed monomials against P by hand.
   - At G = E₀₁+E₂₃ the values are (0, 0, 2, 1), so P has unit normalization.
7. **O-invariant forms are exactly d = 0.**
   - `invariantForm_O` shows this is necessary using the reflection at E₀₁+E₂₃, and sufficient for every orthogonal R.
   - `O_classification` restricts from `SO_classification`; O-invariance implies SO-invariance via `hR.1`.
8. **Reflection-odd forms inside SO are exactly a = b = c = 0.**
   - `invariantForm_odd` and `odd_classification` prove this.
   - `oddSpace` is defined as `SOInvariant ∧ ReflectionOdd`, so it is the minus eigenspace inside the SO-invariant forms, as the fidelity record says.
9. **Dimensions 4/3/1 belong to the real submodules.**
   - `soSpace`, `oSpace` and `oddSpace` are defined from the invariance predicates themselves, not from a candidate list.
   - `so/o/odd_dimension` follow from `LinearEquiv.ofBijective`, where surjectivity comes from the classification theorems.
   - The counts agree with representation theory: M₄ = 1 ⊕ Sym₀(9) ⊕ Λ⁺(3) ⊕ Λ⁻(3) are four pairwise non-isomorphic irreducibles of real type, and the reflection swaps Λ⁺ and Λ⁻.
10. **Even/odd decomposition, including the zero form.**
    - `even_odd_disjoint` shows the intersection is {0}: Q = −Q pointwise.
    - `even_odd_span` shows `oSpace ⊔ oddSpace = soSpace`.
    - `zero_invariant` covers the zero form, and the documents state the overlap at zero explicitly.

### D4.3: native identification

11. **Full-span comparison, not just counts.**
    - `S11_lean_d4_source_check.py` asserts equal ranks, stacked rank = count and identical RREF for V1 (4), V2 (3) and V6 (1), all over 136 columns. The recorded polynomials agree with my hand expansion of P under g_{4i+j+1} = G_ij.
    - The instrument injects `QG_ALL` itself as `sp.Symbol('g_i', real=True)`. The native `declared_symbol` builds exactly that symbol (`sympy_audit.py:72-73,104`), so the injection is faithful.
12. **Native generator orientation.**
    - `delta = A·G − G·A` is the infinitesimal form of G ↦ RGRᵀ.
    - The rows of `action_rows` are images of monomials, so invariance of Σ c_p m_p is Mᵀc = 0. The transpose at line 517 is therefore correct.
    - The wrong-orientation control, run in memory, keeps the counts 4/3/1 but fails the span comparison. This is exactly the counts-versus-object distinction the policy asks for.
    - The same-count wrong span (P replaced by G₀₀²) is also rejected.
13. **Reflection operator.** `V6_OPERATOR.T * V1_BASIS == reflection_rows` is an exact 136-column identity. Because V1 is in RREF, this confirms the operator in the actual native basis.
14. **The two P conventions are distinguished correctly.**
    - `PD_POLY == P` exactly: the RREF pivot is (1,11), i.e. G₀₁G₂₃, with coefficient +1.
    - `Σ ε_{ijkl} G_ij G_kl == 2P` with `LeviCivita(0,1,2,3) = +1`. I confirmed this independently: only antisymmetric parts contribute, ¼·8·Pf = 2Pf.
    - P_D is defined by the spec's §7 rule (the sum of the emitted V6 basis, never rescaled), not by the name of the invariant.
    - The step prose (`steps/S11_stray_longitudinal.md:55`) names the D=4 extra as ε_{ijkl}G_{ij}G_{kl}. That is an inferred name, and the checked factor of 2 is recorded rather than normalized away. That is the correct handling.
15. **Evidence boundaries are stated honestly.**
    - The Wolfram link is only a source-inspection assertion that `Transpose[actionRows]` is present. I confirmed the matching generator block at `.wl:899-911`.
    - The Python/SymPy link is described as tested translation, not a kernel-certified bridge.

### D4.4: controls

16. **Mutations fail for the intended reasons.**
    - The two source mutations (trace-square cross coefficient 2→3; trace-of-square coefficient 2→−2) fail with `unsolved goals` inside the named declarations. The recorded goals show the false polynomial identities.
    - All 12 statement mutants reduce to `⊢ False`.
    - The instrument rejects results caused by an unknown identifier, a parse error, heartbeat or memory limits, or a timeout, and records timeouts separately as `accepted_as_mutation: False`.
    - The positive controls cover a nonzero odd form, the zero form, negative coefficients and unique coefficients.
    - `single_entry_not_invariant` is a real full-group sensitivity test: a quarter turn sends E₀₀ to E₁₁.
17. **Audit.** 49 declarations use only `propext`, `Classical.choice` and `Quot.sound`. Builds use `-DwarningAsError=true`, and a source scan rules out `sorry`, `admit` and `axiom`.

## Non-blocking observations (optional)

- **A. Packet-name mismatch.** The documents call the fixed manifest `_measurements/S11_lean_d4_review_packet.json`, but this packet supplies `MANIFEST.json`. Recording the aggregate hash under one name when dispositions are written would avoid ambiguity.
- **B. The odd character is only proved for one reflection.** `ReflectionOdd` uses a single reflection. The stronger statement, Q(RGRᵀ) = det R · Q(G) for every orthogonal R on `oddSpace`, follows at once from `orientation_conjugate` but is not stated as a theorem. The contract does not require it.
- **C. Shallow statement controls.** The omission and dimension mutants reuse proved lemmas, so they check that the stated claims are sensitive. They do not perturb the reconstruction certificate itself; the kernel's check of `linear_combination` covers that. A mutated certificate coefficient in `Constraints.lean` would be an optional extra demonstration.
- **D. Imprecise step prose (outside the Lean claim).** `steps/S11_stray_longitudinal.md:56-58` attributes the D=4 extra to "cross-pairings between isomorphic summands". At D=4, Λ⁺ and Λ⁻ are non-isomorphic SO(4) irreducibles; the extra invariant is the difference of their two separate norms (∝ P), not a cross-pairing. The corrected count of 4 is right, and no formal statement depends on this wording.

## Limits of my verification

- **Not done:**
  - I did not compile any Lean or run `#print axioms`.
  - I did not rerun the contract or source instruments or the generator.
  - I did not recompute any file or aggregate hash, and did not execute SymPy or Wolfram.
- **Author-reported, not reproduced:** build success, including run5's reuse of the run4 builds, the axiom output, the mutation outcomes and the native span results. I reviewed their recorded diagnostics and the code that produces them.
- **What I did check:**
  - I read all 13 modules, the audit root, both instruments, their JSON records, the validation record, the relevant generator sections, and the native Q9 source and spec sections.
  - I hand-checked:
    - the 12 signed monomials of `orientationForm`;
    - the free-coefficient map;
    - five necessary equations and four reconstruction certificates;
    - the native P polynomial;
    - the ε = 2P factor;
    - the Pfaffian transformation law;
    - the representation-theoretic 4/3/1 counts.
  - I did not check the other 127 equations or 128 certificates individually. Their correctness rests on Lean's kernel check, as reported.
- **Out of scope, as instructed:** divergence and zero bulk variation of P, EL/V5, D5, dynamics, S11c, and export/comparator work.
