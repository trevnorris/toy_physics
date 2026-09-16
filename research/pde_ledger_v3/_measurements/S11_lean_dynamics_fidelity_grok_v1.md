# D2 odd-invariant dynamics E1–E4 fidelity review

**Packet aggregate SHA256:** `16a1dd5784c344f032716d86974b52a58c79110464e2cb1353b915c17d2581c0`  
(manifest method: SHA256 of UTF-8 `json.dumps(files, sort_keys=True)`)

**Verdict: CLEAR**

This is one independent non-author review of the fixed packet only. I did not edit files, rebuild Lean, rerun the instruments, recompute file or aggregate hashes, or execute native production/S11c jobs. Object-hash correspondence is taken from the recorded verification documents, cross-checked against the manifest and internal reports.

No unresolved blocker remains for the bounded E1–E4 contract.

---

## Scope read

Reviewed under `lean/FORMALIZATION_POLICY.md` against `lean/s11/DYNAMICS_COVERAGE.md`, `DYNAMICS_FIDELITY.md`, and `DYNAMICS_VERIFICATION.txt`.

Lean objects: `S11OddDynamics.lean` and `Action` / `Variation` / `Mixing`; S10 calculus `Action`, `PlaneWave`, `Analytic`, `Variation`; S11Invariants `Quadratic`, `Rotation`, `Classification`, `Census`, `Controls`.

Instruments and records: both check scripts, both JSON reports, `S11_lean_dynamics_validation.json`.

Native/specification path actually used for this contract: D2 Q9 V6 `P_D`, `package_build` XFORM_EXTRA−MAIN, V5 via V1 RREF coordinates, EL and the two modal routes in `scripts/S11_stray_longitudinal_sympy_audit.py`, plus Wolfram source inspection of the four recorded anchors. The rest of the production driver is provenance only.

---

## E1 — Density and object identity

`S11Invariants.oddPairing` is `invariantForm ![0,0,0,1]`, i.e. `coordinates₀ · coordinates₁ = (G₁₁+G₂₂)(G₁₂−G₂₁)` with row-major `Gᵢⱼ`. That is the already classified reflection-odd SO(2) pairing.

Lean spacetime is `(t,x₁,x₂)`: `Jet 2 = Fin 3 → Vec 2`, `spatialGradient J j i = J j.succ i`. Native functions are `uⱼ(x₁,x₂,t)` with `g[i,j] = ∂_{x_{i+1}} u_{j+1}`. The two conventions match after that index shift.

Phase: Lean `waveCovector ω k = Fin.cases (-ω) k`, so `φ = k·x − ω t`, same as the specification and Wolfram `phase`.

Beta is an unrestricted real constant, including zero; no sign or nonvanishing premise.

`lagrangian β J = -β/2 · (div J) · (skew J)` and `density_identity` identifies this with `-β/2 · oddPairing(spatialGradient J)`.

Native `compute_q9(2)` yields
`P_D = G₁₁G₁₂ − G₁₁G₂₁ + G₁₂G₂₂ − G₂₁G₂₂ = (G₁₁+G₂₂)(G₁₂−G₂₁)`.
`package_build` puts `Term(β/2, "P_D", pd_density)` in XFORM_EXTRA stiffness, and `L = T − W`, so

**XFORM_EXTRA − MAIN = −β P_D / 2**

Both sign and half-factor are the actual package difference, not a 1-dimensional normalization choice. The source instrument compares against `compute_q9` output and the live package difference; it does not insert a target polynomial into the census engine.

---

## E2 — First variation, momenta, relative action, non-nullness

`momentum` is `deriv (s ↦ L(J + s • basisJet j i)) 0`, an actual derivative of the supplied density. `momentum_eq` gives

|  | u₁ | u₂ |
|---|---|---|
| t | 0 | 0 |
| x₁ | −β c / 2 | −β d / 2 |
| x₂ | β d / 2 | −β c / 2 |

with `d = ∂₁u₁+∂₂u₂`, `c = ∂₁u₂−∂₂u₁`. This is `∂L/∂Gᵢⱼ`.

Lean EL is `Eᵢ = −∑ⱼ ∂ⱼ πⱼᵢ`. `eulerLagrange_eq` is
`E = (β/2) (∂₁c − ∂₂d, ∂₁d + ∂₂c)`, obtained from those momenta without commuting mixed derivatives of `u`.

`relativeAction` is `∫ (L(u+s h) − L(u))`, the pointwise density change. Smooth backgrounds need not have finite total action. The quadratic expansion `lagrangian_change`, compact-support integrability, `HasDerivAt` of the integrated polynomial, and S10 `coord_integration_by_parts` give

`deriv(relativeAction) 0 = ∫ E·h`.

`actionStationary_iff_eulerLagrange` upgrades vanishing against all compact tests to pointwise `E=0` by the compact-support fundamental lemma, continuity, and full-support Lebesgue measure.

`witnessField = planeWave 0 ![1,0] ![1,0] = (cos x₁, 0)` is smooth, and `E(0) = (0, −β/2)`. For every `β ≠ 0` there is an admissible compact test with nonzero first variation (`exists_nonzero_firstVariation`). `variationallyNull_iff` is **all** smooth backgrounds and compact tests iff `β=0`.

That is not a claim that every background is nonstationary for `β≠0`, and not a general null-Lagrangian or finite-domain boundary classification. The statements match the contract’s bound.

---

## E3 — Modal operator and exhaustive mixing

`modalAction β k a = L(modeJet 0 k a) = −β/2 (k·a)(r·a)` with `r = turn k = (−k₂, k₁)`.

`modalOperator M a = −β/2 ((r·a) k + (k·a) r)` is identified in two ways, for arbitrary real `ω,k,a,x`:

- `modal_variation`: amplitude derivative of `modalAction`
- `eulerLagrange_planeWave`: `E(planeWave) = cos(φ) • M a`

`polarization_decomposition` covers every amplitude when `k≠0`. Directional actions:

- `M k = −β/2 |k|² r`
- `M r = −β/2 |k|² k`

Cross pairings are both `−β/2 |k|⁴`. `Mixes` is both of those nonzero. `mixing_iff`: Mixes iff `β≠0` and `k≠0`.

`mixing_cases` is exhaustive and disjoint:

1. `β=0` (any k, including the zero intersection): `M=0`
2. `β≠0` and `k=0`: `M=0` (no nonzero polarization frame)
3. `β≠0` and `k≠0`: both cross-sector elements nonzero

This is the odd increment only. No claim about complete XFORM_EXTRA roots, positivity, or kernels.

---

## E4 — Compact native identification and controls

Selected-helper environment: AST extraction of the named constructors only, with injected symbols. No production driver, registry, emitters, or S11c.

Eight exact identities, as recorded and as checked against the helpers:

| Native | Lean |
|---|---|
| Q9 `P_D` | `oddPairing` |
| XFORM_EXTRA − MAIN | `−β P/2` |
| native EL (`+div π`, spec Q1) | `−` Lean EL |
| V5 odd combination via actual V1 RREF coordinates of `P_D` | native EL `= −(β/2)` that combination |
| route A (cosine stripped) | `−M` |
| route B (Hessian of averaged density) | `M/2` |
| period average | `modalAction/2` |
| stripped factor | `cos(phase)` |

V5 coefficients were `[0,1,0,0]`, reconstructed against live `PD_POLY`, not a hardcoded engine target.

In-memory native mutations of `Term(β/2, "P_D", …)` to the opposite sign and to a missing half-factor fail the intended identities.

Wolfram is source inspection only. The four anchors are present:

- `"XFORM_EXTRA", {{muR/2, "curl"}, {bComp/2, "div"}, {beta/2, "pd"}}`
- `lagrangianJet = Total[kineticTermsJet] - Total[stiffnessTermsJet];`
- `eulerOperatorForGradientDensity[…] := Module[` (also `+div`)
- `Integrate[Expand[planeLagrangian], {phaseVariable, 0, 2 Pi}]/(2 Pi)`

No Wolfram execution, export, or comparator clearance is claimed. Native `period_average` is `sin² → 1/2` substitution; that is the tested translation, not a kernel-certified CAS bridge.

### Ten rejected controls — mathematical failures

**Source mutations (primary polynomial/component mismatch; later tactic failures secondary):**

1. **`density_sign`** (`−β/2` → `+β/2` in `lagrangian`): `density_identity` unsolved with opposite coefficients of `β J₁₀ J₁₁`, `β J₁₀ J₂₀`, `β J₂₁ J₁₁`, `β J₂₁ J₂₀` (left `+1/2` vs right `−1/2`). The later `change` failure in `modalAction_eq` is secondary.
2. **`density_factor`** (drop the `1/2`): `density_identity` unsolved with full coefficients on the left versus half on the right. Same secondary `change` failure.
3. **`local_EL_sign`** (drop the minus on `∑ ∂π`): `eulerLagrange_eq` unsolved on both components, opposite coefficients of `∂₁c`, `∂₂d`, `∂₁d`, `∂₂c`. The later rewrite miss in `eulerLagrange_planeWave` is secondary.

**Seven false statements, each reducing to `⊢ False` in `contract_control` (exit 1, no unknown identifier / timeout / environment failure):**

4. `VariationallyNull 1`
5. `¬ VariationallyNull 0`
6. `¬ Mixes 1 ![1,0]`
7. `Mixes 0 ![1,0]`
8. `Mixes 1 0`
9. `¬ Mixes (-2) ![3,4]`
10. `eulerLagrange 1 witnessField 0 1 = 1/2`

### Eight admissible passing controls

The corresponding true statements, plus `exists_nonzero_firstVariation` at `β=1` on `witnessField`. Negative nonzero beta is included. Zero intersections are included. The witness need not have finite background action.

Axiom audit: 47 `#print axioms` lines, only `propext`, `Classical.choice`, `Quot.sound`. No `sorry`/`admit`/custom physics axioms in the new modules.

Four S10 `.olean` hashes changed on the fresh sequential rebuild with unchanged sources; S11Invariants objects did not. The E1–E4 record binds the new objects. Historical H1–H4 / I1–I4 hashes are not rewritten.

---

## Blockers

None for this bounded contract.

## Optional suggestions (not conditions of closure)

- An explicit compact bump witnessing the first variation would be more constructive than the classical existence argument from the fundamental lemma; the contract does not require it.
- Native `period_average` is a `sin²→1/2` rewrite rather than an integral; that is already the stated translation limit.
- V6 RREF happens to give coefficient 1 on `(tr G)(skew G)`; the live identity checks this, so no extra formalization is needed.

---

## Limits of this review

I did not independently recompute SHA256 hashes, rebuild Lean, or rerun either instrument. I did not execute Wolfram or production CAS. Compiled `.olean` files are not in the packet; their hashes are those recorded in `S11_lean_dynamics_validation.json` / `S11_lean_dynamics_contract_checks.json`. This review does not close E4’s second independent leg, and it does not claim H1–H4, I1–I4, D3–D5, full XFORM_EXTRA spectrum, variable beta, or S11c.

**Final verdict: CLEAR**
