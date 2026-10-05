**Verdict: CLEAR**

Reviewed packet identifier (as supplied, not recomputed): **S11 D3 bulk variation K1–K4 fidelity contract v1**, `aggregate_sha256` `bb5fe9bbfbc09f3a184bb6f256daba74f9ee7d6baa85cbb67b6e4eb2ddb83b54`.

This is a statement-fidelity review of the classified constant-coefficient D3 family on \(\mathbb{R}^{3+1}\), with smooth compact spacetime tests. It is not a compilation rerun, native re-execution, or hash recomputation.

---

## Definitions, assumptions, and quantifiers checked

| Object | What was checked |
|---|---|
| Family | \(Q(G)=a(\operatorname{tr} G)^2+b\operatorname{tr}(G^2)+c\operatorname{tr}(GG^T)\), \(G_{ij}=\partial_i u_j\), \(L=-Q/2\), arbitrary real \(a,b,c\); no inertia, positivity, or nonzero-coefficient hypotheses |
| Link to J1–J4 | `lagrangian` identified with `S11D3Invariants.invariantForm v` on `spatialGradient`; unique coefficients via `SO_unique` |
| Domain | Smooth fields on \(\mathbb{R}^{3+1}\) (Lean coordinates \((t,x_1,x_2,x_3)\)); tests are smooth and compactly supported in spacetime; background total action need not be finite |
| Variation | Momenta by differentiating \(L\) in jet slots; finite `relativeAction`; first derivative at \(s=0\); IBP; stationarity against all compact tests \(\Leftrightarrow\) pointwise local EL |
| Equivalence | Operator equality on **every** smooth field and **every** point, classified by \((c,a+b)\); nullness is identically vanishing first variation on every smooth background |
| Modal map | Density-derived operator for every real \(k\), including \(k=0\); identity with `S11Homogeneous.modalOperator 0 c (a+b+c)`, not a new spectrum claim |
| Native link | Actual Q9 V1 rref basis change, every V5 basis vector, arbitrary coefficients, opposite native EL sign |

Scope exclusions in the coverage contract were respected: no D4/D5, no general null-Lagrangian classification, no variable coefficients, no new spectrum/stability/interface, no S11c/production/export, no kernel-certified CAS parser.

---

## Findings

### K1 — density, momenta, finite variation, local PDE

`S11D3Bulk.lagrangian` is \(-Q/2\) with the three classified traces (`Action.lean`). `density_identity` matches `invariantForm_apply`; `all_invariant_densities` reuses `SO_unique`. Spatial rows are `J i.succ`; time-jet momentum is identically zero.

Differentiating the quadratic expansion gives
\[
P_{ri}=-(a\,\delta_{ri}\operatorname{div}u+b\,\partial_i u_r+c\,\partial_r u_i),\qquad P_{0i}=0
\]
(`momentum_eq`). Lean EL is \(-\operatorname{div} P\). Expanding and commuting mixed derivatives produces, on every smooth field,
\[
\mathrm{EL}=(a+b)\nabla(\operatorname{div}u)+c\Delta u
\]
(`eulerLagrange_eq`). Commutation is `partial_commute`, from mathlib `isSymmSndFDerivAt` on `ContDiff ℝ ∞`, not an extra field hypothesis. The three components use explicit swaps (`hcomm 2 1 1`, `3 1 2`, and cyclic).

`relativeAction` integrates only the density change under a compact variation (`Variation.lean`). `relativeAction_hasDerivAt` differentiates the proved quadratic expansion. `integrated_variation_by_parts` uses S10 `coord_integration_by_parts`. `actionStationary_iff_eulerLagrange`: compact tests give a.e. vanishing; smoothness plus full-support Lebesgue measure upgrade to a pointwise PDE. This is the correct finite contract: backgrounds need not have finite total action.

### K2 — current, bulk equivalence, nullness, dimensions

`boundaryCurrent` is \(F_i=\sum_j(u_i\partial_j u_j-u_j\partial_j u_i)\). `boundary_identity` proves \(\operatorname{div} F=(\operatorname{div}u)^2-\operatorname{tr}(G^2)\); the \(u\cdot\nabla^2 u\) terms cancel by the same smoothness commutation. `null_density_is_divergence` then gives \(L[a,-a,0]=-(a/2)\operatorname{div} F\). That is a divergence identity, not a pointwise zero density, and it does not discard boundary physics.

`BulkEquivalent` is equality of EL on all smooth fields and points. Necessity uses two smooth plane waves (`bulkEquivalent_iff`):
- transverse \(k=e_1\), \(a=e_2\) isolates \(c\);
- longitudinal \(k=a=e_1\) isolates \(a+b+c\), hence \(a+b\) once \(c\) matches.

Sufficiency is the closed local formula, not generic-\(k\) sampling. `responseMap v=!(c,a+b)` is identified with this operator equality (`bulkEquivalent_response`, `nullSpace_identification`), so \(\dim\ker=1\) and \(\dim\mathrm{range}=2\) are dimensions of the variational classification, not of an unrelated coefficient map. `response_surjective` shows every \((c,a+b)\) is attained.

`VariationallyNull` (stationarity on every smooth background) is equivalent to \(c=0\) and \(a+b=0\) (`variationallyNull_iff`), including the zero density and negative coefficients. Parameterization is \(\operatorname{span}\{(1,-1,0)\}\) (`null_parameterization`, `null_nonzero_example`). Outside that locus, `exists_nonzero_firstVariation` gives a smooth background and compact test with nonzero first variation; the admissible control instantiates this at `![-2,3,-5]`.

### K3 — modal operator and homogeneous identity

The density-derived modal operator is
\[
M a_{\mathrm{amp}}=-c|k|^2 a_{\mathrm{amp}}-(a+b)(k\cdot a_{\mathrm{amp}})k
\]
(`modalOperator`), for every real \(k\) including zero (`modal_zero_wavevector`). It is the Fourier symbol of the local EL (`eulerLagrange_planeWave`, `momentum_contraction`). Longitudinal stiffness is \(a+b+c\); transverse stiffness is \(c\). `homogeneous_operator` is the algebraic identity with `S11Homogeneous.modalOperator 0 c (a+b+c)`; \(\omega\) is dummy. This is not a new root or stability claim.

Passing numerical controls: \(v=(1,2,3)\), \(k=e_1\) gives longitudinal \(M k=-6k\) and transverse \(M e_2=-3e_2\).

### K4 — native Q9 V5 / EL correspondence

Native `compute_q9` takes the rref of the SO(3) Lie nullspace; a basis index is not assumed to name a trace form. The recorded change of basis to \((\operatorname{tr}^2,\operatorname{tr}(G^2),\mathrm{Frob})\) is
\[
\begin{pmatrix}0&1&0\\1/2&-1/2&0\\0&-1&1\end{pmatrix},
\]
invertible (\(\det=-1/2\)). Solving for general \((a,b,c)\) independently yields native weights \((a+b+c,\,2a,\,c)\), matching the recorded `general_native_basis_weights`.

Native `euler_lagrange_from_placeholders` (and Wolfram `eulerOperatorForGradientDensity`) returns \(+\operatorname{div}P\), the negative of Lean \(-\operatorname{div}P\). Hence native \(\mathrm{EL}(L)=-\)Lean \(\mathrm{EL}(L)\) and \(\mathrm{V5}(Q)=2\) Lean \(\mathrm{EL}(L)\) for \(L=-Q/2\). Displayed V5 vectors match that: V5 of \(\operatorname{tr}(G^2)\) is \(2\nabla(\operatorname{div}u)\); V5 of \(\tfrac12(\operatorname{tr}^2-\operatorname{tr}(G^2))\) is mixed-derivative commutators (null after smoothness); V5 of \(\mathrm{Frob}-\operatorname{tr}(G^2)\) is \(2\Delta-2\nabla(\operatorname{div}u)\) after commuting. Native mixed derivatives are reduced with `.doit()`; Lean proves the same commutation. Native coordinates are \((x_1,x_2,x_3,t)\) versus Lean \((t,x_1,x_2,x_3)\). Wolfram was source-inspected only (`Transpose[actionRows]` in the V1 generator). Translation is tested, not kernel-certified.

Native sign/factor mutations break V5/EL identities while leaving the homogeneous \(\mu=c\), \(B=a+b+c\) algebra intact; target-side null/non-null and wrong-\(B\) controls are present.

### Mutations and recorded evidence

Twelve intended mathematical rejections and eleven positives are recorded. Inspected diagnostics:

- `density_sign` / `density_half`: unsolved `density_identity` with \(+1/2\) or \(-1\) in place of \(-1/2\).
- `bulk_coefficient_map` (\(a+b\mapsto a-b\)): unsolved `homogeneous_operator` (`v 1 - v 0 = -v 1 + -v 0 ∨ …`); at \(v=(0,1,0)\), \(k=a=e_1\) this is \(+e_1\) versus \(-e_1\).
- `boundary_current_sign`: unsolved `boundary_identity` displaying a false expanded identity, plus an unused `partial_sub` simp warning promoted by `-DwarningAsError`. Rejection is justified by the false identity (and the smooth counterexample \(u=(x_1,0,0)\): mutated divergence \(2\), claimed density \(0\)), not by the linter warning.
- Statement mutants such as `frobenius_not_null_mutant` reduce to `⊢ False`; matching positives pass (`VariationallyNull 0`, `![1,-1,0]`, `¬ VariationallyNull ![0,0,1]`, `BulkEquivalent ![1,0,0] ![0,1,0]`, dimensions \(2/1\), \(M(0)=0\)).

No inspected rejection is a syntax/import/timeout/resource failure. Axiom audit text lists 49 roots depending only on `propext`, `Classical.choice`, and `Quot.sound`. Source/object/pin correspondence is recorded in the validation JSON; those hashes were read, not recomputed.

---

## Limits of this independent verification

Performed: read the five `S11D3Bulk` modules and audit root, J-link (`invariantForm`, `SO_unique`), S10 jet/IBP/plane-wave machinery, homogeneous modal operator, native Q9/V5/EL helpers, Wolfram EL/V1 generator, coverage/fidelity/verification texts, contract and source-check records, and independent expansion of EL, current, modal map, basis inverse, and sign/factor relations.

Not performed: Lean compilation, Python/Wolfram execution, SHA-256 recomputation, or a fresh review of unrelated H/I/E/J proofs beyond the objects this increment actually uses.

---

## Blockers versus optional improvements

**Blockers:** none.

**Optional improvements (not clearance conditions):**
- `exists_nonzero_firstVariation` is classical (`push Not` from \(\neg\forall\)) rather than an explicit compact bump. The mathematics is sound: non-null coefficients make a plane-wave EL nonzero, and compact tests then detect it. A constructed bump would be a clearer witness.
- The boundary-current mutation log mixes a genuine false identity with an unused-simp warning-as-error; the adjudication already treats them separately.

K1–K4, as stated, match the proved objects, quantifiers, and native sign/basis map.
