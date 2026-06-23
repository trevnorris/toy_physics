T0_FROZEN(8fa41ac51e88)

# pathA_24 T0 Freeze Artifact: Minimal Polar-OP "Little-Arrows" Parent-Action Extension

artifact_status: append-only immutable freeze
directive: `software/stage1_solver/directives/pathA_24_brane_existence_defect_structure.md`
stage: `T0 anti-circularity preregistration gate`
date: 2026-06-23
scope: reports-only; no T1 solve; no wall, spectrum, light, charge, or throat computation

This artifact freezes the target-blind rules for the minimal polar order-parameter extension. It is a contract for later rungs, not a result about whether the extension succeeds.

## 1. Verdict

`T0_FROZEN(8fa41ac51e88)`

Classification: `INDEPENDENTLY_MOTIVATED_NEW_PARENT_ACTION`, conditional on the polar-orientation postulate. The current scalar GNLS parent action does not derive this order parameter.

Reason the artifact does not emit `FAIL_NO_MINIMAL_POLAR_LAGRANGIAN`: a local minimal polar-vector Landau/Frank action can be written using only the current GNLS medium variables plus the independently motivated substructure fact that the medium's constituents carry a polar orientation.

Reason the artifact does not emit `POLAR_OP_AD_HOC`: the postulate is not justified here by wall/light/charge payoffs. Its non-payoff motivation is the 3He-style precedent that a single superfluid medium can carry orientational order as an internal material degree of freedom. The framework-circularity caveat remains active: selecting a polar vector rather than a scalar/tensor/director is not derivable from the current scalar action.

Sources: `directives/pathA_24_brane_existence_defect_structure.md` ("The mechanism under test", "Honest classification", and T0); `docs/conceptual_foundation.md` secs. 1-2; `research/pde/paper/pde.tex` sec. "Exact parent action"; `decisions/13_emergent_constants_derivation.md` sec. 8; `decisions/14_value_provenance_and_calibration_map.md` sec. 1.

## 2. Frozen Minimal Polar-OP Lagrangian

### 2.1 Existing medium facts inherited without modification

The existing GNLS medium remains the carrier:

- Bulk coordinates: `X^i=(x,y,z,w)`, `i=1..4`.
- Medium field: `psi=sqrt(rho) exp(i theta)`, with `rho=|psi|^2`.
- Existing GNLS matter action: `pde.tex`, eq. `parent-Lpsi`.
- Existing equation of state: `P(rho)=K rho^5`, `U(rho)=K rho^5/4`.
- Existing sound speed: `c_s^2(rho)=5 K rho^4/m`; at the background, `c_s0^2=5 K rho0^4/m`.
- Existing flow velocity where `rho>0`: `v_i=(hbar/m) partial_i theta - (q_star/m) A_i`.
- Existing length scale: `a`, conditionally identified with `hbar/(m c_s0)` in the current provenance map, otherwise a branch collective scale already present in the program.

No canonical document is changed by this report.

### 2.2 New field content

```freeze-action
Field:
  P^i(X,t), i=1..4, a dimensionless polar vector carried by the GNLS medium.

Polarity:
  P and -P are distinct microscopic orientations. There is no director quotient P ~ -P.

Carrier rule:
  The polar field has no standalone vacuum action. Every dynamical OP term is weighted by rho and uses the GNLS material flow v_i. When rho -> 0 the OP energy and inertia vanish with the medium.

Material derivative:
  D_t^v P^i := partial_t P^i + v^j partial_j P^i,
  v_i := (hbar/m) partial_i theta - (q_star/m) A_i.

Local sound-speed function:
  c_s^2(rho) := 5 K rho^4 / m.

Frozen T0 polar-OP Lagrangian density:
  L_pol =
      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)
    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)
    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2.

Frozen extended parent-action grammar:
  S_T0 = S_GNLS_existing + int dt d^4X L_pol.

Frozen baseline branch:
  O(4)-isotropic soft-spin polar vector with one-constant Frank stiffness,
  single-well magnitude potential |P|=1, no explicit easy axis, no brane-localized term,
  no gauge charge assigned to P except through v_i in the material derivative.
```

### 2.3 Units and dimensions

Base dimensions: length `L`, time `T`, mass `M`.

| quantity | dimension | note |
|---|---:|---|
| `d^4X dt L_density` | action | `S` has units `M L^2 T^-1` |
| `L_density` | `M L^-2 T^-2` | energy per four-volume |
| `rho` | `L^-4` | 4D number density |
| `m` | `M` | GNLS constituent mass |
| `K` | `M L^18 T^-2` | from `U=K rho^5/4` |
| `c_s` | `L T^-1` | `c_s^2=5K rho^4/m` |
| `a` | `L` | existing medium length scale |
| `P^i` | `1` | dimensionless orientation vector |
| `D_t^v P` | `T^-1` | material orientation rate |
| `partial_j P` | `L^-1` | orientational gradient |
| `m rho a^2` | `M L^-2` | OP inertia density coefficient |
| `m rho c_s^2 a^2` | `M T^-2` | one-constant Frank coefficient |
| `m rho c_s^2` | `M L^-2 T^-2` | self-potential depth density |

The three OP coefficient scales are fixed as functions of existing medium quantities. No independent OP speed, OP mass density, wall tension, light speed, charge coupling, throat scale, or Coulomb normalization is introduced.

### 2.4 Term-by-term target-blind justification

| term | included form | target-blind justification | source |
|---|---|---|---|
| Polar vector field | `P^i in R^4`, `P != -P` as states | Implements the single new substructure fact: constituents carry head-tail polar orientation. | directive "mechanism under test"; `docs/conceptual_foundation.md` sec. 2 |
| Density carrier | every term weighted by `rho` | The arrows are a property of the medium's constituents, not a second empty-space field. | `docs/conceptual_foundation.md` secs. 1-2; T0 single-medium test |
| Material derivative | `D_t^v = partial_t + v^j partial_j` | A carried internal orientation is advected by the GNLS medium flow already present in the parent action. | `pde.tex` exact current/continuity and gauge-invariant velocity |
| Kinetic term | `(1/2)m rho a^2 |D_t^v P|^2` | Minimal local quadratic inertia for an internal orientation with moment scale set by existing constituent mass, density, and length `a`. | `pde.tex` matter sector; `decisions/13` sec. 8 F1/F3; `decisions/14` sec. 1b |
| Frank stiffness | `(1/2)m rho c_s^2 a^2 |nabla_4 P|^2` | Minimal one-constant O(4)-isotropic penalty for neighboring constituents pointing differently; scale set by existing cohesion/sound energy. | `docs/conceptual_foundation.md` sec. 1; `pde.tex` EOS/sound speed |
| Self-potential | `(1/4)m rho c_s^2 (|P|^2-1)^2` | Minimal local polynomial selecting nonzero orientation magnitude while preserving O(4) and `P -> -P` degeneracy. | polar substructure fact plus current scalar single-well density action |
| Flow/rho coupling | only through `rho`, `c_s(rho)`, and `v_i` | Uses existing medium variables; avoids a decoupled director with standalone constants. | T0 operational single-medium test; `pde.tex` GNLS velocity and EOS |

## 3. Allowed vs Explicitly Excluded Invariants

### 3.1 Included/allowed in the frozen minimal action

| invariant | status | one-line reason | source citation |
|---|---|---|---|
| `rho |D_t^v P|^2` | included | Lowest-order local kinetic scalar for an orientation carried by the medium flow. | `pde.tex` current/continuity/velocity; conceptual foundation sec. 2 |
| `rho c_s^2 a^2 (partial_j P^i)(partial_j P^i)` | included | One-constant O(4)-isotropic Frank stiffness is the minimal gradient disagreement cost. | conceptual foundation sec. 1; parent EOS/sound-speed law |
| `rho c_s^2 (P^2-1)^2` | included | Minimal magnitude-selecting self-potential that does not choose an axis. | T0 target-blind minimality; parent scalar `U(rho)` remains single-well |
| coefficient dependence on `rho`, `K`, `m`, `a` via `c_s(rho)` | included | Keeps the OP material scales tied to existing medium quantities. | T0 operational single-medium test; `decisions/14` sec. 1 |
| global O(4) spatial rotations of `P^i` and `X^i` | included symmetry | The current bulk arena is isotropic before any solved wall branch selects a direction. | `pde.tex` bulk-coordinate conventions and flat metric |
| global `P -> -P` energy symmetry without identifying states | included symmetry | With no external polar source, opposite arrows are degenerate but still distinct states. | polar `+w != -w` requirement in directive and conceptual foundation |

### 3.2 Explicitly excluded from the frozen minimal action

Exclusions are on intrinsic symmetry, locality, provenance, or minimality grounds only. No exclusion below is made because it helps or hurts a later wall, light, charge, Coulomb, `e`, `alpha`, or `c_gamma=c_s` payoff.

| invariant or structure | excluded form | one-line intrinsic reason | source citation |
|---|---|---|---|
| Explicit easy-axis anisotropy | `-gamma_w (P dot w_hat)^2`, `(P dot w_hat)^4`, or any fixed `w_hat` tensor | The current scalar GNLS action supplies no fixed internal axis; adding one lowers O(4) by hand. | `pde.tex` exact parent action; directive T0 framework-circularity caveat |
| Linear polar source | `b_i P^i` | Requires an external polar vector absent from the current parent action. | parent action field content in `pde.tex` |
| Headless director quotient | `P ~ -P` | Contradicts the stipulated polar head-tail substructure. | directive "mechanism under test"; conceptual foundation sec. 2 |
| Brane-localized OP stiffness | `delta(w) |nabla_3 P|^2`, `Z(w)|nabla P|^2`, or support on a prescribed wall | Presupposes the wall/brane before T1 derives or fails it. | directive T0 before T1; `pde.tex` projection-vs-reduction firewall |
| Coupling to imposed confinement geometry | `V_conf f(P,w_hat)` or geometric profiles selecting `P` | Uses the existing imposed brane machinery to choose the OP structure. | `pde.tex` `V_conf` role; conceptual foundation sec. 2 obstacle |
| Direct gauge charge for `P` | `A_i P^i`, `F_ij P^i P^j`, or covariant phase assigned to `P` | The new datum is neutral material orientation; the existing gauge coupling already acts through `psi` and `v_i`. | `pde.tex` matter/gauge source bookkeeping and vorticity-gauge identity |
| Chiral/Lifshitz terms | one-derivative parity-odd orientation terms | Need an intrinsic handedness or epsilon-contracted structure not supplied by the polar-arrow postulate alone. | parent action has no chiral substrate datum |
| Multi-constant Frank splitting | independent splay/twist/bend constants | Symmetry allows richer elasticity, but one-constant Frank is the lowest-parameter isotropic local form. | T0 minimality/scoring axes |
| Higher-derivative gradients | `(nabla^2 P)^2`, quartic gradient terms | Excluded at first rung by locality and lowest-derivative effective-action order. | T0 minimality/locality |
| Independent displacement/elastic metric field | `u^a`, Cosserat frame, or extra solid network added beside GNLS variables | Adds another field system beyond the one polar OP carried by the medium. | T0 single-medium test; conceptual foundation one-medium rule |
| Payoff-tuned coefficient ratios | choosing coefficients to force a wall tension, photon speed, charge normalization, or Coulomb scale | Not an intrinsic invariant; all such choices violate the forbidden-information rule. | directive T0 forbidden-information rule |

## 4. Discrete Modelling-Choice List, Priors, and Branch Budget

Branch policy: the baseline branch below is the pre-committed first branch. Later rungs may run additional preregistered branches only according to the target-blind budget below. Any post-T0 structural rescue that affects a payoff requires a fresh directive/T0 or the whole run is labeled `AD_HOC_RESCUE`.

| choice | allowed values | prior ranking and target-blind reason | branch budget |
|---|---|---|---|
| OP target/magnitude | A. soft-spin polar vector `P in R^4` with single well `|P|=1`; B. hard-spin unit vector `n in S^3`; C. two-component/Ising-like construction; D. tensor/headless director | 1: A, because it is the simplest local polynomial Landau action with all T0-requested terms and no fixed axis. 2: B, because it is the hard-constraint limit of A with fewer radial dynamics but an infinite constraint. 3: C, because it introduces extra labels/fields. Prune D because it violates polarity if headless. | Run A as baseline. Run B only as a target-blind hard-spin limit/sensitivity if A's radial mode is identified as the only obstruction or if the T1 equations are ill-conditioned by amplitude collapse. Prune C unless a new T0 freezes its added labels. Prune D in this directive. |
| Alignment axis `w` | A. spontaneously selected by boundary/domain data within O(4)-isotropic action; B. explicit background axis `w_hat`; C. axis tied to pre-existing `V_conf`/brane profile | 1: A, because no axis exists in the current scalar parent action. 2: B, nonminimal because it adds a fixed tensor. 3: C, worst because it imports imposed brane structure. | Run A. Prune B/C in the minimal run. A future explicit-axis run must be separately frozen and cannot rescue this run's baseline verdict. |
| Self-potential shape | A. single radial well `(P^2-1)^2`; B. hard constraint `P^2=1`; C. component double-well in `P dot w_hat`; D. multi-well engineered potential | 1: A, because it is O(4)-invariant and polynomial. 2: B, limit of A. 3: C, requires explicit axis. 4: D, nonminimal. | Run A baseline. Run B only as the hard-spin limit noted above. Prune C/D in this directive unless a fresh T0 supplies independent microscopic provenance. |
| Frank stiffness | A. one-constant O(4) stiffness; B. multiple independent Frank constants; C. chiral/odd-gradient terms | 1: A, by isotropy and minimal parameter count. 2: B, symmetry-allowed but less minimal. 3: C, requires intrinsic chirality. | Run A. Prune B/C for this directive unless the one-constant form is mathematically ill-posed on intrinsic grounds, not because of a payoff. |
| Flow coupling | A. advective material derivative using GNLS `v_i`; B. free director time derivative `partial_t P` only; C. added co-rotational/vorticity-response coefficient; D. direct gauge charge for `P` | 1: A, because the arrows are carried by the medium. 2: C, plausible richer microstructure but needs extra closure. 3: B, less faithful to "carried by" but useful only as a negative-control in later discriminators. 4: D, not motivated by neutral orientation. | Run A baseline. Do not use B/D as physical branches. C is pruned unless a future T0 derives a material-frame coupling from substructure rather than payoff. |
| Coefficient scale closure | A. fixed medium-scale functions exactly as in the frozen block; B. same functions times independent dimensionless constants; C. arbitrary independent OP scales | 1: A, because it is the strictest single-medium closure. 2: B, honest sensitivity but adds OP-specific constants. 3: C, second-medium drift. | Run A. B/C are pruned in the main run. If later freeing coefficients is necessary for numerics, record the independent count; two or more independent OP constants trigger the single-medium drift flag. |
| Density dependence | A. local `rho`, `c_s(rho)` weighting; B. background-only `rho0`, `c_s0`; C. independent fixed OP background | 1: A, because the OP is carried locally by the medium. 2: B, less local but still medium-related. 3: C, independent. | Run A. B only as an explicitly labeled approximation if a later numerical implementation freezes background coefficients before solving; C pruned. |

Target-blind pruning criterion recorded now: prune any branch that adds fields, fixed tensors, nonlocality, independent OP constants, or brane-localized structure not forced by the polar-orientation postulate plus the existing medium action. Do not prune or promote a branch because of wall stability, wall tension, light-mode count, speed matching, charge signs, Coulomb behavior, `e`, `alpha`, or throat balance.

## 5. Provenance Class Ledger

| item | provenance class | honest note |
|---|---|---|
| Polar vector OP carried by medium | `INDEPENDENTLY_MOTIVATED_NEW_PARENT_ACTION` | New parent-action DOF, not derivable from current scalar GNLS. Independently motivated by 3He-style orientational order in a single superfluid, not by T1-T5 payoffs. |
| Polarity `P != -P` as states | `INDEPENDENTLY_MOTIVATED_NEW_PARENT_ACTION` | Directly part of the little-arrows substructure fact. |
| Soft-spin O(4) vector with `|P|=1` single radial well | `POSTULATED` | Minimal local Landau realization of the polar-vector postulate; not derived from GNLS. |
| No explicit easy axis in baseline | `POSTULATED` | Target-blind minimal symmetry choice; records the possibility that `w` may fail to emerge/stabilize later. |
| One-constant Frank stiffness | `POSTULATED` | Minimal isotropic gradient closure; richer Frank elasticity is not derived. |
| Advective material derivative using `v_i` | `INDEPENDENTLY_MOTIVATED_NEW_PARENT_ACTION` | Implements "carried by the one medium" through existing GNLS flow variables. |
| Medium-scale coefficient closure | `POSTULATED` | Strong single-medium closure tying OP scales to `rho`, `c_s`, `m`, and `a`; dimensionless equality to one is a freeze choice, not a derivation. |
| Hard-spin limit branch | `POSTULATED` | Target-blind limit of the soft-spin action. |
| Explicit easy-axis branch | `POSTULATED` only if separately frozen from independent microscopic anisotropy; otherwise `AD_HOC` | Pruned in this run. Adding it after seeing a wall failure is `AD_HOC_RESCUE`. |
| Two-component/Ising-like construction | `POSTULATED` or `AD_HOC` depending on future provenance | Pruned here because it adds labels beyond one polar vector. |
| Headless director/tensor alternative | not admitted in this directive | Contradicts the polar `+ != -` postulate and would be a different framework. |

Framework-circularity caveat: this is not a fully target-blind survey over every possible OP type. It freezes the polar-vector framework already selected by the directive. The non-payoff provenance is the 3He-style fact that a one-medium superfluid can carry orientational order; without that precedent, the polar-vector postulate would downgrade to `AD_HOC`.

## 6. Operational Single-Medium Test

Test criterion frozen by the directive: if two or more OP parameters are independent of the medium's existing parameters, the OP is functionally a second field and heads toward `FAIL_SECOND_MEDIUM_DRIFT` at later rungs.

### 6.1 Frozen parameter ledger

| OP action parameter or scale | frozen value | class | independent count | note |
|---|---|---|---:|---|
| OP inertia coefficient `I_P(rho)` | `m rho a^2` | medium-related | 0 | Built from GNLS mass density and existing length scale. |
| Frank stiffness `K_P(rho)` | `m rho c_s^2(rho) a^2` | medium-related | 0 | Uses existing sound-energy scale and length scale. |
| Self-potential depth `V_P(rho)` | `m rho c_s^2(rho)` | medium-related | 0 | Uses existing local medium energy density. |
| Magnitude minimum `P0^2` | `1` | convention/fixed normalization | 0 | Orientation magnitude is dimensionless; normalization can be absorbed into `P`. |
| Material-advection coefficient | `1` in `D_t^v` | kinematic/medium-derived | 0 | Fixed by using the GNLS flow to carry the OP. |
| Explicit easy-axis strength `gamma_w` | `0` | fixed by baseline symmetry/minimality | 0 | No fixed axis in the frozen minimal action. |
| Direct gauge charge of `P` | `0` | fixed by neutral orientation postulate | 0 | Gauge enters only through medium velocity `v_i`. |
| Brane-localized OP stiffness | `0` | fixed by no-preimposed-wall rule | 0 | No `delta(w)`, `Z(w)`, or imposed support profile. |
| Independent OP speed | none | absent | 0 | The OP's natural scale is tied to `c_s(rho)`; no `c_P` is introduced. |
| Independent OP length | none | absent | 0 | Uses existing `a`; no new correlation length is introduced. |

### 6.2 Single-medium verdict for the frozen action

Independent OP parameters: `0`.

Single-medium status: `SINGLE_MEDIUM_CLOSURE_FROZEN`, with an honesty caveat. The closure is medium-related rather than micro-derived. It follows the 3He benchmark structurally by tying orientation stiffness to medium quantities, but unlike real 3He no deeper substrate parameters such as a derived gap `Delta` or Fermi velocity `v_F` have been derived inside the current toy model. If later rungs free the inertia, stiffness, potential depth, anisotropy, or OP speed as fitted constants, the independent count must be recomputed. If the count becomes `>=2`, record `FAIL_SECOND_MEDIUM_DRIFT` pressure immediately; do not hide it under "carried by the medium" language.

## 7. Forbidden-Information Compliance Statement

I affirm that the T0 selections and exclusions above used only:

- the existing GNLS medium action, EOS, velocity, density, and length-scale provenance;
- the single new substructure fact that constituents carry a polar orientation;
- target-blind criteria: locality, O(4) symmetry, field/DOF count, minimality, and single-medium fidelity.

The following downstream payoffs were not used to admit, exclude, tune, rank, or normalize any term:

- wall existence, wall lifetime, wall tension, surface tension value, or wall-core texture;
- shear mode count, photon signature, `c_gamma`, `c_gamma=c_s`, or any `lambda_gamma` target;
- bulk shear leakage, Magnus preservation, magnetism, or sector separation;
- two charge signs, charge universality, Coulomb form, `e`, `alpha`, throat radius, or throat stability;
- dark-energy/cosmology flow behavior;
- any T1-T5 outcome label.

Temptations deliberately rejected:

- adding an explicit easy-axis anisotropy to make `+w/-w` disconnected;
- forcing a hard Ising double-well in `P dot w_hat`;
- choosing coefficient ratios to make any later wave speed equal to `c_s`;
- inserting brane-localized stiffness, `Z(w)` weighting, or a flat-core ansatz;
- coupling `P` directly to gauge fields to pre-arrange charge/Coulomb behavior;
- adding a separate displacement/Cosserat field to make a MacCullagh interpretation easier.

Binding rule for later rungs: a post-T0 structural change that affects a payoff requires either a fresh directive with a new T0 freeze or the label `AD_HOC_RESCUE` for the whole run.

## 8. Content Stamp / Hash

Hash target: the exact bytes inside the fenced block labelled `freeze-action` in section 2.2, excluding the fence lines.

Verification command from project root:

```sh
awk '/^```freeze-action$/ {on=1; next} /^```$/ && on {exit} on {print}' software/stage1_solver/reports/pathA_24_T0_freeze.md | sha256sum
```

Frozen-action SHA-256:

```text
8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064
```

Short artifact hash used in the verdict: `8fa41ac51e88`.

Append-only rule: future notes may be appended below this line only if explicitly labeled as post-freeze commentary. The frozen action block and the T0 verdict/hash lines are not to be edited during T1-T5.

---

## POST-FREEZE COMMENTARY (orchestrator audit trail — does NOT modify the frozen block or hash)

**Tri-review of T0 (2026-06-23):** orchestrator review (target-blindness, completeness vs T0 spec, dimensional homogeneity — all
three terms reduce to the action density `M L⁻² T⁻²`) + frozen-action **hash re-verified** (recomputed
`8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064`, matches §1) + a **clean adversarial agent** review. Adversarial
verdict: **`T0-SOUND-WITH-NOTES`** — no circularity defeat; the freeze forecloses no downstream FAIL and *actively protects* several
(T1c stability stays reachable — no easy axis baked in; the `A_i P^i` exclusion keeps `FAIL_DIRECTION_NO_FIELD_COUPLING` reachable;
nothing pre-bakes `BULK_SHEAR_FREE`). Three scoping notes (applied here as commentary; the frozen action is unchanged):

1. **Bulk OP mode speed `=c_s` is NOT a frozen-in `λγ=1` (cleared).** The Frank/inertia ratio gives the *bulk* spin-wave speed
   `√(K_P/I_P)=c_s` identically — but that is a legitimate "one medium → one energy/speed scale, used minimally" consequence (no
   new `c_P` is introduced, §6.1), **not** a tuning toward the forbidden `c_γ=c_s`. `c_γ` is the **wall-localized** MacCullagh shear
   speed `√(μ_brane/ρ_brane)`, a **T2 output** computed from profile-weighted integrals of the texture — a different quantity from
   the bulk spin-wave speed. So `λγ` is **not** pre-decided; T2 can still emit `c_γ≠c_s` or `FAIL_WRONG_MODE_DISPERSION`.
2. **§6.1 framing (honesty).** "Independent OP parameters: 0" is partly a consequence of the **strictest-closure convention**
   (branch A: dimensionless coefficient ratios fixed to 1) — read the §6.1 rows as `convention-fixed (branch A)`, consistent with
   the (adequate) §6.2 caveat. Single-medium legitimacy stands; only the headline was generous. Branch B (independent dimensionless
   constants) remains the recorded sensitivity; ≥2 freed ⇒ `FAIL_SECOND_MEDIUM_DRIFT`.
3. **The advective derivative `D_t^v P = ∂_t + v^j∂_j` is the MINIMAL frame; the co-rotational/Jaumann (vorticity) term `Ω·P` is
   pruned-to-a-future-T0 (load-bearing scoping for T2).** For a *vector* OP carried by a swirling flow, the co-rotational term is
   the standard choice — and it is mechanically *adjacent to the Frank→MacCullagh coupling T2 tests* (rotation coupling to the
   medium's swirl/displacement). Pruning it makes the baseline T2 a **conservative** test (harder to pass — the able-to-fail
   direction). **Therefore a baseline T2 FAIL must be scoped as `FAIL under the minimal advective frame`, NOT "the carried-vector
   light concept is dead"** — the co-rotational branch is the principled next move (a separate, independently-motivated T0 freeze;
   adding it post-failure without that is `AD_HOC_RESCUE`).

**T0 status: COMPLETE and reviewed.** None of the notes block T1. Awaiting user gate before T1.
