**CLEAR**

Packet: `S11 full D4 bulk D4C.1–D4C.4 fidelity contract v1` (`MANIFEST.json` aggregate SHA256 `061d0d96f0c3b62439875877aa4a107b53b6870a206b308ee8a04056bd8e8a11`, as supplied; not recomputed). Author: Codex. Reviewer: Grok. Increment reviewed: D4C.1–D4C.4 only.

No required fidelity fixes. Optional stronger results are correctly not claimed.

---

## Checked assumptions and conventions

The supplied family is the already classified SO(4) quadratic density

\[
Q=a(\operatorname{tr} G)^2+b\operatorname{tr}(G^2)+c\operatorname{tr}(GG^T)+\beta P,\qquad
P=F_{01}F_{23}-F_{02}F_{13}+F_{03}F_{12},\quad F=G-G^T,
\]

with action \(L=-Q/2\). Coefficients \(v=(a,b,c,\beta)\) are arbitrary constant reals, including zero and negative values; they are parameters, not varied fields. The only field varied is \(u:\mathbb{R}^{4+1}\to\mathbb{R}^4\).

Lean conventions match the contract:

| Object | Identification |
|---|---|
| Spacetime | `Point 4 = Fin 5 → ℝ`; index `0` is time |
| Spatial gradient | `G_{ij}=\partial_i u_j` via `J i.succ j` (`Even.spatialGradient`, `S11D4Odd.spatialGradient`) |
| Native coordinates | `u_functions(4)` lists `(x1,x2,x3,x4,t)`; native `x_(i+1)` is Lean spatial row `i.succ` |
| Forms | `divergence = tr G`, `transposePair = tr(G^2)`, `gradientSquare = tr(GG^T)`, `orientation = P` |
| Time momentum | identically zero; no inertia; no dimensional rescaling |

`Even.lagrangian` uses only \(a,b,c\). The full density adds the unchanged reviewed odd term `S11D4Odd.lagrangian (v 3)`.

---

## D4C.1 — density, momenta, actual first variation

**Object.** `density_identity` is \(L(v,J)=-\frac12\,\texttt{invariantForm}\,v\,(G(J))\). `all_invariant_densities` reuses `SO_unique`, so this is the entire classified SO(4) family with unique coefficients, not a sampled span.

**Momenta.** `momentum` is the actual derivative of \(L\) in the `basisJet` direction. `momentum_eq` splits even plus odd. `momentum_time` is the zero time row. Hand expansion of the even formula recovers

\[
\pi=-a(\operatorname{tr} G)I-b G^T-c G-\frac{\beta}{2}\partial P/\partial G,
\]

with odd `dualCurl` equal to \(\partial P/\partial G\).

**Variation.** `relativeAction` is the integral of the pointwise density change on a smooth background and a smooth compactly supported test. Integrability is proved for the compact linear and quadratic remainder terms only; an infinite absolute background action is not assumed. The quadratic expansion gives a real `HasDerivAt`. Integration by parts (`coord_integration_by_parts`) and the fundamental lemma (`ae_eq_zero_of_integral_contDiff_smul_eq_zero`, then continuity plus Lebesgue) give

`actionStationary_iff_eulerLagrange`: stationarity against every `TestField` iff the pointwise local EL vanishes.

Lean EL is \(-\operatorname{div}\pi\). That is the operator appearing in the first variation.

---

## D4C.2 — bulk operator and equivalence

On every smooth four-component field, `eulerLagrange_eq` (via `eulerLagrange_even` and the reviewed odd cancellation `S11D4Odd.eulerLagrange_zero`) is

\[
\mathrm{EL}=(a+b)\nabla(\operatorname{div} u)+c\Delta u,
\]

with spatial `gradDiv`/`laplacian` and actual commuting coordinate derivatives (`partial_commute`). No coefficient sign is excluded.

`BulkEquivalent` quantifies over all smooth fields and all spacetime points. Sufficiency uses the EL identity. Necessity uses two smooth waves as discriminating witnesses, not as the domain:

- transverse `planeWave 0 e₁ e₂` isolates \(c\);
- longitudinal `planeWave 0 e₁ e₁` isolates \(a+b+c\).

Controls `transverse_response` (\(-7\)) and `longitudinal_response` (\(-12\)) for \(v=(2,3,7,11)\) match that split. No generic-only sampling and no positivity hypothesis are used.

---

## D4C.3 — response, null family, currents

`responseMap` is the actual linear map \(\mathbb{R}^4\to\mathbb{R}^2\), \(v\mapsto(v_2,v_0+v_1)=(c,a+b)\). It is surjective; image and kernel both have dimension 2 (`bulk_response_dimension`, `null_dimension`).

`VariationallyNull` is zero actual first variation for every smooth background and compact test. This is equivalent to \(c=0\) and \(a+b=0\), and to the two-parameter family \((t,-t,0,\beta)\) (`variationallyNull_iff`, `null_parameterization`, `nullSpace_identification`). This classifies the already classified quadratic family; it is not a theorem about every null Lagrangian.

Currents:

- even \(J_i=\sum_j(u_i\partial_j u_j-u_j\partial_j u_i)\), \(\operatorname{div} J=(\operatorname{div} u)^2-\operatorname{tr}(G^2)\);
- odd \(K_i=\frac12\sum_j u_j M_{ij}\) with \(M=\texttt{dualCurl}\), \(\operatorname{div} K=P\).

`null_density_is_divergence` is the corresponding weighted sum. Pointwise null density and momenta remain nonzero (`even_null_density_nonzero`, odd witness reuse). The affine field `affineWitness` has \(\operatorname{div} J=2\), so boundary effects are retained. The compact homogeneous identification is \(\rho=0\), \(\mu=c\), \(B=a+b+c\), including \(k=0\) (`homogeneous_operator`, `modal_zero_wavevector`). No spectrum, stability, or boundary-value theorem is claimed.

---

## D4C.4 — controls and native translation

The sixteen paired controls are concrete statement mutations, not canonical source replacements. Sampled recorded diagnostics (e.g. `odd_momentum_sign_mutant`, `odd_momentum_factor_mutant`, `even_current_factor_mutant`) are a single `contract_control` error of the form `unsolved goals / ⊢ False`, with no warning or instrument failure. The twenty positives are eighteen distinct sources; even/odd momentum positives are repeated for separate sign and factor pairs. Extra positives cover the mixed null family, \(k=0\), negative coefficients, and existence of a nonzero actual first variation off the null family.

Native instrument (AST-selected Q9/coordinate/EL helpers only; no production import or driver):

- V1 rows in the \((\operatorname{tr}^2,\operatorname{tr}(G^2),\operatorname{tr}(GG^T),P)\) basis are \((0,1,0,0)\), \(\bigl(\tfrac12,-\tfrac12,0,0\bigr)\), \((0,-1,1,0)\), \((0,0,0,1)\);
- determinant \(-\frac12\); general native weights \((a+b+c,2a,c,\beta)\);
- native helper EL is \(+\operatorname{div}(\pi)\), opposite Lean’s \(-\operatorname{div}(\pi)\);
- native V5 of \(Q\) is twice Lean EL of \(L=-Q/2\); \(P_D=P\); zero time row;
- two actual mutations of the original `euler_lagrange_from_placeholders` append line; seven target checks (three wrong-formula rejections, four positives).

This translation is outside the Lean kernel. Wolfram strings are source inspection only.

---

## Optional, not required

- An explicit compactly supported first-variation witness; `exists_nonzero_firstVariation` is classical existence from the exhaustive iff, which is what the contract asked for.
- A general null-Lagrangian theorem, variable coefficients, D5, interface theory, spectrum/scattering/pole results, or a kernel-certified CAS bridge.

---

## Verification limits

Inspected: `MANIFEST.json`, `FORMALIZATION_POLICY.md`, `D4_BULK_COVERAGE.md`, `D4_BULK_FIDELITY.md`, `D4_BULK_VERIFICATION.txt`, all six `S11D4Bulk` modules and audit root, the classified density/odd/S10/homogeneous import closure used by the statements, native helper source, both instruments, and recorded native/formal/validation reports.

Did **not** run Lean, the native instrument, hash recomputation of the packet, Mathlib rebuild, or the mutation suite. Compilation, axiom lists, object hashes, package pins, and historical-file hashes remain author evidence. Historical hashes are preservation evidence, not a claim that every historical report is in this packet. The separate S11c calculation session was not reproduced.
