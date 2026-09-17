**Verdict: CLEAR**

**Reviewed packet identifier (as supplied, not recomputed):**  
`S11 D4 odd divergence/variation D4B.1–D4B.4 fidelity contract v1`  
`aggregate_sha256 = 1c2fb3840e56f82e506adbddb13463a6fbc645193294a3173531026f0ae838d5`

This is a statement-fidelity review of D4B.1–D4B.4 only. No required corrections. Optional improvements are listed at the end and are not blockers.

---

## Assumptions and conventions checked

The supplied odd term is the isolated XFORM_EXTRA summand. Shared specification: `L = T − W` and `W_XFORM_EXTRA` contains `(β/2) P_D` (`directives/S11_SHARED_PHYSICS.md` §2, §7, Q1, Q9). The isolated odd density is therefore `L = −(β/2) P_D`, with `β` an arbitrary real constant (no sign or nonzero premise).

Checked conventions:

| Object | Convention | Status |
|---|---|---|
| Coordinates | Lean `Point 4`: `(t,x1,x2,x3,x4)` with time = index 0. Native `u_j(x1,x2,x3,x4,t)` | Match by named spatial `xs[i]` and named `t`, not by argument position |
| `G_ij` | `∂_i u_j`; spatial row = derivative index | `spatialGradient J i j = J i.succ j`; native `g[i,j] = ∂_{x_{i+1}} u_{j+1}` |
| `F_ij` | `G_ij − G_ji` | `antisym J i j = J i.succ j − J j.succ i` |
| `P` | `F01 F23 − F02 F13 + F03 F12` | `orientation` / `orientationJet_eq` |
| `L` | `−β P / 2` | `lagrangian` |
| `M` | complementary matrix in the coverage contract | `dualCurl` matches entrywise |
| Momenta | `∂L/∂G_ij = −(β/2) M_ij`, zero time row | `momentum` by actual `deriv`; `momentum_eq` |
| Lean EL | `−∑_j ∂_j(∂L/∂J_ji)` | opposite overall sign from native Q1 / `euler_lagrange_from_placeholders` |
| Domain | all smooth `u` on `R^{4+1}`, all smooth compact `h`, all real `β` | no G-symmetry, transverse ansatz, time-independence, or `β≠0` |

Native `P_D` is the unrescaled sum of the computed `V6_BASIS` rows (`compute_q9`: `pd_poly = sum of V6 rows`). The compact check asserts `V6_DIM = 1` and `actualP − P = 0`. The recorded expanded `P` is the reviewed polynomial, not the fully summed epsilon form `2P`.

---

## D4B.1 — supplied density and actual variation

`all_odd_densities` is a direct invocation of the existing exhaustive theorem `S11D4Invariants.odd_classification`: SO-invariant and reflection-odd iff `Q = invariantForm ![0,0,0,d]`. That is the D4 quadratic classification on spatial `G`, reused rather than re-proved.

`density_identity` is the normalization link to that same family:

```
lagrangian β J = −(1/2) · invariantForm ![0,0,0,β] (spatialGradient J)
```

Because `invariantForm ![0,0,0,β] = β · orientation`, this is exactly `L = −(β/2) P`. A factor-of-two error is not absorbed.

Momentum is not a separately written matrix. It is

```
momentum β J j i := deriv (s ↦ lagrangian β (J + s • basisJet j i)) 0
```

and `momentum_eq` identifies it with `Fin.cases 0 (fun r => −(β/2) * dualCurl J r i)`. Time row (`j = 0`) is zero; spatial entries are `−(β/2) M_ri`.

`linearDensity` is the contraction of those momenta against a variation jet. `linearDensity_eq` proves it equals `variationDensity`, which is the actual directional derivative of `L` (`lagrangian_variation` / `lagrangian_increment`). So the first-variation density is the differentiated action, not an independent operator.

`eulerLagrange` is built from those same momenta, with S10’s first-variation sign `−div π`.

Witnesses prevent a zero-versus-zero reading:

- `nonzero_density`: `orientationJet witnessJet = 1` on `G01 = G23 = 1`
- `nonzero_momentum`: `momentum 2 witnessJet 1 1 = −1`, i.e. `∂L/∂G_01 = −(β/2) F23` at `β = 2`
- `witness_lagrangian (−2) = 1`, so negative `β` is live
- `current_normalization`: `K_0 = 1` for `v = (0,2,0,0)` on that jet, fixing the `1/2`

Lean index `momentum … 1 1` is derivative `x1` and field `u_1`, which is native `momenta[0,1]`. The native witness is the same `−1`.

No inertia, even densities, positivity, or modal restriction enters these definitions.

---

## D4B.2 — explicit divergence

On every `SmoothField u`:

- `dualCurl_contraction`: `∑_{ij} G_ij M_ij = 2P`
- `dualCurl_divergence_zero`: `∑_i ∂_{x_i} M_ij = 0` for each `j`
- `currentJet`: `K_i = (1/2) ∑_j u_j M_ij`
- `boundary_identity`: `div K = P`
- `lagrangian_is_divergence`: `L = −(β/2) div K`

The factor `1/2` is required by the degree-two identity `∑ G M = 2P`. It is not arbitrary.

Mixed partials are commuted only under `ContDiff ℝ ∞`, via `partial_commute` (`isSymmSndFDerivAt`) inside `dualCurl_divergence_zero`. That is the right regularity.

This is a pointwise bulk identity. Compact support is used later for the integrated first variation. Boundary contributions are not claimed to vanish for non-compact fields.

---

## D4B.3 — exhaustive odd-family bulk conclusion

`relativeAction β u h s` is the actual integral

```
∫ [ L(jet(u + s • h)) − L(jet u) ]
```

Integrability uses only the compactly supported increment (`relative_density_integrable`), via `fieldJet_add_smul` and the quadratic expansion. The background density need not be integrable.

The derivative of this integral at `s = 0` is obtained from the proved expansion (`relativeAction_hasDerivAt`), not from a separately posited EL pairing. Integration by parts is S10’s `coord_integration_by_parts` (compact support, including time). Then

```
deriv (relativeAction β u h) 0 = ∫ ∑_i eulerLagrange β u x i * h x i
```

`eulerLagrange_zero` holds for every smooth `u` and every real `β`. Therefore `firstVariation_zero` and `every_background_stationary` cover all smooth backgrounds, all smooth compact variations, and all constant `β`, including `0` and negatives. Time dependence is unrestricted: time momenta are present in the sums and are proved to be zero, not omitted.

`P = 1` and `momentum = −1` keep this from being read as “the density is identically zero.” Stationarity is bulk first variation of a divergence, not pointwise vanishing of `P`.

---

## Compact native link

The instrument extracts the original helpers (`compute_q9`, `q9_vector`, `q9_row_to_poly`, `matrix_from_rows`, `derivative_placeholders`, `u_functions`, `coordinate_substitution`, `to_coordinate`, `euler_lagrange_from_placeholders`, `q9_v5`) and does not run the production driver.

Checked identities, from the instrument source and the recorded objects:

- actual `P_D = P`, no hand-chosen substitute
- derivative momenta `= −(β/2) M`, and not `+(β/2) M` or `−β M`
- time momenta zero
- `∑ G M = 2P`, `div M = 0`, `div K = P`
- native EL of `L = −(β/2) P` is zero after `expand/doit`
- native V5 on the actual odd combination in the computed V1 basis is zero; recorded weights are `(0,0,0,1)`
- adding `G00²` yields a non-null operator whose first component is `−∂_{x1}^2 u_0` (native sign)
- mutating the EL helper in memory to append `0` makes even that non-null density return zero

Native Q1 / the helper return `+div π`. Lean EL is `−div π`. Both bulk expressions vanish, so sign/scale are taken from nonzero momenta and from the `G00²` control, not from `0 = 0`.

Wolfram evidence is the instrument’s source-string check of `Transpose[actionRows]`. No Wolfram execution is claimed, and none was performed here.

This is exact tested SymPy translation, not a kernel-certified CAS bridge.

---

## Mutations and positives

Four source mutations fail on named mathematical identities:

1. **Density normalization** (`L = −β P` instead of `−β P/2`). `density_identity` remains with unsolved goal `−β P = −(β/2) P`. At the witness with `β = 2` that is `−2` versus `−1`.
2. **Complementary sign** (`M_10 = +F23`). `momentum_eq` remains with `F23 = −F23 ∨ β = 0`. At the witness, actual `∂L/∂G_10 = +1` versus predicted `−1`.
3. **Wrong divergence index** (every `∂_i` replaced by `∂_{x1}`). `dualCurl_divergence_zero` remains an unsolved second-derivative identity. Independently, for `u = (0,0,0,x1 x2)` the `j = 0` component is `1` rather than `0`. Later rewrite/simp/`function expected` failures are downstream and were not treated as evidence.
4. **Current factor** (`1` instead of `1/2`). `boundary_identity` remains with `1 · Σ = (1/2) · Σ`. Independently, for `u = (0,x1,0,x3)` one has `P = 1` and doubled current divergence `2`. The later `change` failure in `current_normalization` is secondary.

Eight paired wrong statements reduce to `⊢ False` in `contract_control` (`P = 0`, momentum `+1`, current `2`, contraction `1`, first variation `1`, EL component `1`, `L(−2) = −1`, time momentum `1`). Their eight true partners pass. Extra positives: compact zero test field, stationarity at `β = −7` and `β = 0`, and smoothness of the affine field `(0,x1,0,x3)`.

These are mathematical rejections, not import/timeout/warning failures. Canonical universal quantifiers are retained; witnesses do not replace the proofs.

The recorded axiom audit lists the 41 named declarations and only `propext`, `Classical.choice`, and `Quot.sound`. The four new modules contain no `sorry`, `admit`, or custom `axiom`.

---

## Required corrections vs optional improvements

**Required corrections:** none.

**Optional improvements (not needed for D4B.1–D4B.4):**

- The finite relative action is `s² ∫ L(h)` after the first variation vanishes. For compact `h`, that remainder is also the integral of a compactly supported divergence and should vanish. The contract asks for first variation, which is proved; proving `relativeAction ≡ 0` would be a strengthening, not a missing claim.
- The affine positive control only proves `SmoothField`. A field-level identity `orientationJet (fieldJet affine) = 1` would sit beside the jet witness.
- The native wrong-index control uses generic non-vanishing rather than the explicit field `u = (0,0,0,x1 x2)` used in the Lean counterexample.

---

## Limits of this verification

Independently inspected: definitions, quantifiers, proof scripts, S10 derivative/IBP reuse, native helper source, compact instrument logic, recorded mutation diagnostics, and the algebraic witnesses/normalization identities.

Not done in this review:

- no Lean build
- no native instrument execution
- no hash recomputation of the aggregate or of individual files
- no Wolfram execution
- no re-proof of the historical D4 classification at `d6119f75` beyond the new uses (`odd_classification`, `orientation`, `invariantForm`)
- author-reported run2 PASS, 14 reused imports, 41 axiom audits, and object-hash correspondence are taken as recorded evidence after checking that the recorded diagnostics match the intended false identities

Out of scope, as requested: D4 even-family bulk census, variable `β`, general null-Lagrangian classification, D5, full XFORM_EXTRA spectrum, stability/interfaces, S11c, production/comparator/export, and a systematic CAS bridge.

**CLEAR** for D4B.1–D4B.4 statement fidelity. Closure of D4B.4 still requires the second independent non-author review.
