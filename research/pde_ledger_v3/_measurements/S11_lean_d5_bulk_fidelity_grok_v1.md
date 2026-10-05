# D5B.1–D5B.4 independent statement-fidelity review

**Verdict: CLEAR**

Reviewer: Grok (this session). Author: Codex. Packet: `MANIFEST.json` (no `_measurements/S11_lean_d5_bulk_review_packet.json`). Scope: D5B.1–D5B.4 only.

## Evidence used

**Source-reading (this review):**
`lean/FORMALIZATION_POLICY.md`, `lean/AGENTS.md`, `lean/s11/D5_BULK_COVERAGE.md`, `D5_BULK_FIDELITY.md`, `D5_BULK_VERIFICATION.txt`; all five bulk modules plus root (`S11D5Bulk.lean`, `Action.lean`, `Variation.lean`, `Bulk.lean`, `Census.lean`, `Controls.lean`); reused D5 classification (`Forms.lean`, `Classification.lean`, `CoordinateAlgebra.lean`); S10 jet/variation (`Action.lean`, `PlaneWave.lean`, `Analytic.lean`, `Variation.lean`); `S11D4Odd/Calculus.lean`; `S11Homogeneous/Action.lean`; `S10Controls/Action.lean`; native helpers in `scripts/S11_stray_longitudinal_sympy_audit.py` (`compute_q9`, `q9_v5`, `euler_lagrange_from_placeholders`, `derivative_placeholders`, `coordinate_substitution`); `_measurements/S11_lean_d5_bulk_source_check.py` and `S11_lean_d5_bulk_contract_check.py`; recorded native/formal reports and author validation JSON; `lean/verify.py` adjudication; Wolfram file only for the two declared string anchors.

**Author-supplied execution, not independently reproduced:** recorded run4 (`D5_BULK_VERIFICATION.txt`, `S11_lean_d5_bulk_contract_checks.json`), native check report, 60 object hashes, 54 axiom lists revalidated from the run3 root log, portable D5 seed report. Compiled oleans and Mathlib sources are not in the packet. Scripts were not executed.

**Not described as a fresh build:** the 42 seeded `S11D5Invariants` objects copied from the completed portable D5 replay; the 54 `#print axioms` lists; guarded reuse of the other twelve import objects. Run4’s 32 controls are the fresh control executions in the author record.

## Domain, dimension, and normalization

Three coefficients are not five field components. `Coeff := Vec 3` is independent of `Vec 5`. Spacetime is `Point 5 = Fin 6 → ℝ` (time row `0`, spatial row `i.succ` = native `x_(i+1)`). Jets are `Jet 5 = Fin 6 → Vec 5`. The spatial gradient is a full `Matrix (Fin 5) (Fin 5) ℝ` with 25 independent entries.

Row-derivative convention holds:

```14:20:lean/s11/S11D5Bulk/Action.lean
def spatialGradient (J : Jet 5) : S11D5Invariants.Mat := fun j i => J j.succ i
...
def lagrangian (v : Coeff) (J : Jet 5) : ℝ :=
  -(1 / 2 : ℝ) * (v 0 * divergence J ^ 2 + v 1 * transposePair J + v 2 * gradientSquare J)
```

So `G_{ij} = ∂_i u_j`. `density_identity` identifies this with `L = -Q/2` for `Q = invariantForm v G`, and `invariantForm_apply` is `a (tr G)^2 + b tr(G^2) + c tr(G G^T)`. `all_invariant_densities` reuses the completed unique SO(5) classification (`SO_unique`); O(5) is the same family by that existing theorem. Coefficients are arbitrary reals, including zero and negative; nothing in the definitions restricts sign.

Momentum is the actual jet derivative, not a supplied tensor:

```57:74:lean/s11/S11D5Bulk/Action.lean
def momentum (v : Coeff) (J : Jet 5) (j : Fin 6) (i : Fin 5) : ℝ :=
  deriv (fun s : ℝ => lagrangian v (J + s • basisJet j i)) 0
...
theorem momentum_time (v : Coeff) (J : Jet 5) (i : Fin 5) : momentum v J 0 i = 0
```

`momentum_eq` gives `p_{ri} = -(a δ_{ri} div u + b ∂_i u_r + c ∂_r u_i)` on spatial rows, matching `p_{ij} = -(a δ_{ij} div + b ∂_j u_i + c ∂_i u_j)`. There is no kinetic term: time momentum vanishes, the Laplacian sums only spatial second derivatives, and `modalOperator` has no `ω²` inertia.

## D5B.1 — actual finite relative action

`relativeAction` integrates the pointwise density change, not a difference of two background integrals. `relativeAction_expansion` is the quadratic expansion; `relativeAction_hasDerivAt` differentiates that expansion; `integrated_variation_by_parts` uses S10 `coord_integration_by_parts` on compact test fields; `actionStationary_iff_eulerLagrange` upgrades compact testing to pointwise vanishing of a smooth EL by the Mathlib fundamental-lemma lemma. Domain is smooth backgrounds and smooth compactly supported variations (`SmoothField`, `TestField`). Mixed derivatives come from existing `S11D4Odd.partial_commute`.

Lean EL is the standard `-div p`:

```84:85:lean/s11/S11D5Bulk/Action.lean
def eulerLagrange (v : Coeff) (u : Point 5 → Vec 5) (x : Point 5) : Vec 5 :=
  fun i => -∑ j : Fin 6, coordDeriv j (fun y => momentum v (fieldJet u y) j i) x
```

## D5B.2 — operator, equivalence, census

On every smooth field,

```28:30:lean/s11/S11D5Bulk/Bulk.lean
theorem eulerLagrange_eq {u : Point 5 → Vec 5} (hu : SmoothField u)
    (v : Coeff) (x : Point 5) :
    eulerLagrange v u x = (v 0 + v 1) • gradDiv u x + v 2 • laplacian u x
```

so `EL = (a+b) ∇(div u) + c Δu`. `BulkEquivalent` quantifies over all smooth fields and points. Sufficiency is this identity; necessity uses one transverse and one longitudinal plane wave only. Universal equality holds iff `c` and `a+b` agree (`bulkEquivalent_iff`).

`responseMap v = ![v 2, v 0+v 1]` is surjective onto `Vec 2` (image dimension 2). Kernel dimension 1, identified with variational nullness. The null family is exactly `{ ![t,-t,0] | t : ℝ }` (`null_parameterization`), including `t = 0`. `VariationallyNull` is stationarity of the actual first variation on every smooth background, not a generic-wave witness. Non-nullness yields a nonzero first variation (`exists_nonzero_firstVariation`). `null_positive` retains unrestricted `t`.

## D5B.3 — current, witnesses, homogeneous map

```63:88:lean/s11/S11D5Bulk/Bulk.lean
def boundaryCurrent (u : Point 5 → Vec 5) (x : Point 5) : Vec 5 :=
  fun i => ∑ j : Fin 5,
    (u x i * fieldJet u x j.succ j - u x j * fieldJet u x j.succ i)
...
theorem null_density_is_divergence ...
    lagrangian ![a,-a,0] (fieldJet u x) =
      -a / 2 * (∑ i : Fin 5, coordDeriv i.succ (fun y => boundaryCurrent u y i) x)
```

`J_i = ∑_j (u_i ∂_j u_j − u_j ∂_j u_i)`, `div J = (div u)^2 − tr(G^2)`, and the null density is `-t div J / 2`. Witnesses: density `-1` on `evenJet` for `(1,-1,0)`; current divergence `2` on the matching affine field; momentum `-2` for `(1,0,0)` on that jet. Null bulk variation is not pointwise-zero density and is not absence of boundary effects; the documents do not claim otherwise.

Homogeneous identification is the modal operator, not a spectrum/stability theorem: `rho = 0`, `mu = c`, `B = a+b+c` (`homogeneous_operator`). `modal_zero_wavevector` includes `k = 0`. Negative coefficients appear in the extra positive `BulkEquivalent ![-2,1,-7] ![5,-6,-7]`.

## D5B.4 — native Q9/V5 identification

The instrument AST-selects the original helpers and does not import or run the production driver. `compute_q9(5)` produces a 3-dimensional V1 basis in the 325-monomial quadratic space. Change-of-basis rows to `(tr², tr(G²), tr(GG^T))` are `(0,1,0)`, `(1/2,-1/2,0)`, `(0,-1,1)`, determinant `-1/2`; native weights of `(a,b,c)` are `(a+b+c, 2a, c)`.

Native `euler_lagrange_from_placeholders` appends `time_coord + space_coord`, i.e. `+div(momentum)`, opposite Lean’s `-div p`. Consequently native EL of `L = -Q/2` equals `-` Lean EL, and native V5 of `Q` equals twice Lean EL of `L`. That is checked on the full three-dimensional native span and on each basis element. Native G is `G_{ij} = ∂_i u_j`. Time momenta of `L` vanish. Native sign and factor mutations of the original helper break those identities; a wrong `b`-momentum transpose, wrong `B` map, and current-factor error are rejected. Fifth-direction modal entries match `-7` (transverse) and `-12` (longitudinal). Wolfram is string-anchor inspection only. This translation is outside Lean’s kernel.

## Mutation controls

Fourteen paired controls are actual true/false mathematical statements, not syntax or resource failures. `lean/verify.py` `adjudicate` requires exit 1, exactly one `contract_control` error, exactly one `⊢ False`, and no warnings/`BAD_INSTRUMENT`. Author run4 lists 14 REJECTED and 18 PASS (17 distinct positives; momentum repeated for sign and factor).

Meaningful pairs include: response/kernel dimensions 2/1 versus 3/2 (coefficient dimension is not substituted for spatial dimension 5); `(1,1,0)` is not null (sign of the `b` slot); `(0,0,1)` is not null (Laplacian coefficient kept); null density `-1` versus `0`; momentum `-2` versus `+2` and `-4`; current divergence `2` versus `1`; transverse `-7` versus `-5`; longitudinal `-12` versus `-7`; bulk equivalence of distinct triples with the same `(c,a+b)` versus its negation; time momentum `0` versus `1` on the fifth field component; fifth-spatial-direction transverse/longitudinal `-7`/`-12` versus `0`. Extra positives: full null line, zero wavevector, negative coefficients, existence of nonzero first variation for a non-null triple.

## Coverage versus exclusions

Variable coefficients, interfaces, general-dimensional or general null-Lagrangian classification, spectrum/scattering/poles, S11c, production reruns, and a systematic CAS bridge are not claimed. D5 classification and S10/D4 calculus are unchanged reused dependencies.

## Evidence limits

This clearance is a source-and-record fidelity assessment. I did not rebuild Lean, rerun the 32 controls, rerun the native checker, or re-print axioms. Standard-axiom status of the 54 audited theorems, isolation of the sixty objects, and the precise `⊢ False` shape of each mutant are accepted as the author’s recorded run3/run4 evidence, with the 42 classification oleans treated as validated reuse rather than a new D5 classification execution.

No necessary correction was found. Optional extensions (a more constructive nonzero-variation witness; a less tautological native current-factor polynomial) are not required fixes.
