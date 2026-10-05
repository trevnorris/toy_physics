## Verdict: **CLEAR**

**Packet manifest aggregate SHA-256:** `f95d1667e54b4c8139881feaa86267ebe60eb4e9d85f0f4f4bd2dff9dc2b2747`  
(SHA-256 of UTF-8 `json.dumps(files, sort_keys=True)` as recorded in `MANIFEST.json`.)

This is a non-author fidelity review of the supplied I1–I4 packet only. There is no unresolved blocking finding on the D2 invariant-space contract. Broader formalization is not a closure condition.

---

### What was reviewed

The object, Lean modules, native Q9 sources, mutation/source instruments, and terminal records were checked against `lean/FORMALIZATION_POLICY.md` and Q9 in `directives/S11_SHARED_PHYSICS.md`. Per-file hashes in `MANIFEST.json` agree with `S11_lean_invariant_validation.json`, `INVARIANT_VERIFICATION.txt`, and the two terminal JSON reports for every overlapping path.

---

### I1 — Object and coefficient action

The ambient object is Mathlib `QuadraticForm ℝ Mat` on all real 2×2 matrices, not a postulated monomial list. `quadratic_representation` reconstructs the ten coefficients from the associated bilinear form on the `t,s,x,y` frame; `polynomialForm_surjective` / `polynomial_injective` make that presentation exhaustive and unique.

Coordinates and inverse are the recorded unnormalized map

- `t=G11+G22`, `s=G12-G21`, `x=G11-G22`, `y=G12+G21`
- `G11=(t+x)/2`, `G12=(s+y)/2`, `G21=(y-s)/2`, `G22=(t-x)/2`

with `coordinates_decode` and `decode_coordinates`. Indexing is row-major `G_ij`, matching Q9 `g_((i−1)D+j)=G_ij` and the native `variables[i*n+j]` layout. The physical substitution `G_ij=∂_i u_j` is identification only; entries are independent at the census. Equality is equality of quadratic polynomials on all matrices. There is no Euler–Lagrange or total-divergence quotient in the theorems or in the native span check.

Conjugation is `R * G * Rᵀ`. `SOInvariant` / `OInvariant` quantify over every proper / orthogonal matrix. `proper_rotation` classifies every `Proper` matrix as `rotation a b` with `a²+b²=1`; `orthogonal_invariance_iff` reduces O to SO plus the specified reflection `diag(-1,1)`. No sampled subgroup replaces either quantifier.

`coefficient_action` is the general identity that monomial-image **rows** act on coefficient **columns** through the transpose. Native `compute_q9` now transposes **each generator block before stacking**; the patch does not insert expected dimensions or a candidate basis. Wolfram already has `Transpose[actionRows]` at `buildInvariantCensus` (line 910) and is unchanged. The native generator `[[0,1],[-1,0]]` is the opposite of `d/dθ rotation(cos θ, sin θ)` at the identity; the zero constraint is the same.

---

### I2 — Complete spaces

Selected rotations `(a,b)=(0,1)` and `(3/5,4/5)` appear only in `invariant_coefficients` as **necessary** constraints. Sufficiency is `invariantForm_SO`, which uses `coordinates_rotation` and `spin_two_norm` for **every** proper rotation. The classified SO space is exactly

`A t² + B s² + C(x²+y²) + D t s`,

i.e. the I2 basis `(tr G)²`, `(G12−G21)²`, `(G11−G22)²+(G12+G21)²`, `(tr G)(G12−G21)`, proved complete rather than assumed.

Linear equivalences give dimensions **4 / 3 / 1**. `oddSpace` is `{Q | SOInvariant Q ∧ ReflectionOdd Q}`: reflection-odd is taken inside the SO space, matching Q9 V6 as the `(−1)`-eigenspace within V1. `even_odd_disjoint` and `even_odd_span` give the unique even/odd decomposition, including the zero form. O invariance is `D=0`; reflection-oddness in this space is `A=B=C=0`.

---

### I3 — Compact native identification

The historical discrepancy is the concrete polynomial `G11² + G12 G21 + G22²` at checkpoint source SHA-256 `352fc502…`, with values **1** and **481/625** on `diag(1,0)` and `R=[[3/5,-4/5],[4/5,3/5]]`. `wrongNativeForm` is that polynomial; `wrongNativeForm_witness` proves those values; the source-check residual `-144/625` is the same identity.

The repaired D2 check compares **full row spaces** of native V1/V2/V6 to `{t², s², x²+y², ts}` by exact rational RREF and stacked rank, not by counts. An in-memory undo of the transpose keeps native counts **4/3** and **fails** the span test. That is the required distinction: matching dimension is not treated as fidelity evidence.

Limits are stated and kept: D2 V1/V2/V6 spans only; no Lean proof of the Python interpreter; Wolfram is source inspection of orientation, not a fresh execution; no D3–D5 census, V5/EL, PD-package production, comparator/export, or S11c clearance.

The general transpose identity and the production-site `compute_q9` change apply in every dimension that uses this routine. The packet measures and proves only the D2 census. That boundary is consistent: it is named, not used to infer later work.

---

### I4 — Controls, hashes, H1–H4 separation

Nine mathematical rejections and eight positive controls are recorded. Inspected failures are mathematical:

| Control | Observed rejection |
|---|---|
| Wrong coefficient orientation | Unsolved goal `c i * A i j * m j = c i * m j * A j i` in `coefficient_action` |
| Odd pairing omitted from `invariantForm` | Missing `t s` term in `invariantForm_apply` |
| SO/O/odd dimensions 4/3/1 flipped | `⊢ False` |
| Odd pairing in the even span / as O-invariant | `⊢ False` |
| `wrongNativeForm` as SO-invariant | `⊢ False` |
| `rotation 1 1` as proper | `⊢ False` |

Passing controls include dimensions 4/3/1, the corresponding exclusions, `¬ Proper (rotation 1 1)`, and the admissible 1 vs 481/625 witness. Tool/syntax/resource failures are excluded by the instrument. A secondary unused-`Matrix.transpose_apply` linter error on the orientation mutant is recorded and is not the acceptance criterion.

All 46 audited theorems depend only on `propext`, `Classical.choice`, and `Quot.sound`. The audit list is the full theorem set of the five modules.

H1–H4 remains a separate closed contract that excluded Q9. This packet does not relabel those historical hashes as current after the Q9 repair and the `S11Invariants` Lake target.

---

### Blockers

None.

---

### Optional suggestions (not closure conditions)

- The orientation mutant also trips `linter.unusedSimpArgs` under `warningAsError`. The algebraic mismatch is already the recorded criterion; isolating the linter would only make the log easier to read.
- SymPy leaves the `diag(-1,1)` reflection-difference block untransposed; Wolfram transposes `reflectionRows - I`. For this R0 the block is diagonal in the monomial basis in every D, so D2 V2 is unaffected, as claimed.
- `S11_lean_q9_orientation_probe.json` still carries status `CONFIRMED_SOURCE_FIDELITY_DISCREPANCY` against the pre-repair source. The fidelity record already says not to rerun it as a success test on the repaired file.

This review does not clear the second independent leg, and it does not declare I1–I4 complete as a project milestone. On the supplied packet, the I1–I4 statements, quantifiers, native spans, counterexample, and controls are faithful.
