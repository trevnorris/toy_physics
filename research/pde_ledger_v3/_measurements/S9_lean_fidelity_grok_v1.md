# S9 bounded Lean fidelity review

## 1. Packet and coverage

**Revision:** `S9 bounded fidelity contract v1`  
**Aggregate SHA-256:** `6903a1fa2cedf302561deab7dd04c6e88f523125b372a416cd6448ed4325d0eb`  
**Author (declared):** Codex. Hashes were not recomputed.

Read `lean/FORMALIZATION_POLICY.md`, `lean/s9/COVERAGE.md`, `lean/s9/FIDELITY.md`, then the S9/S10 Lean sources, original S9 constructors, compact source/contract instruments and their filed JSON, `CLOSURE_VERIFICATION.txt`, and the S9 step record as far as C1–C4 require. Compilation was treated as deduction evidence only. No other reviewer’s findings were used. Lean, CAS, and hash tools were not executed.

**Verdict for this bounded Lean scope: CLEAR.**

Non-review obligations C1–C3 and the verification half of C4 are satisfied at the declared level. C4’s two-leg review is this process; it is not a reason to add theorems.

---

## 2. C1 — Action, census, specialization, source link

The supplied density is the original MAIN constructor.

- SymPy `construct_curl_action(rho_br * identity3, mu_R)` with default `stiffness_sign=-1` (`scripts/S9_light_requires_shear_sympy_audit.py:59-65,111`).
- Wolfram `mainLagrangian = rhoBr velocityVector.velocityVector/2 - muR curlVector.curlVector/2` (`mathematica/S9_light_requires_shear_mathematica_audit.wl:33-43`).
- Lean `lagrangian` / `jetCurl` (`lean/s9/S9Pilot/Action.lean:21-26`).

Parameter map `rho_br`/`rhoBr` → `rho`, `mu_R`/`muR` → `mu`; coordinates `(t,x,y,z)`; fields `(u1,u2,u3)`; `J j i = ∂_j u_i`. Curl components match:  
`(∂_y u_3 - ∂_z u_2, ∂_z u_1 - ∂_x u_3, ∂_x u_2 - ∂_y u_1)`.  
Normalization is `ρ/2 |u_t|² − μ/2 |curl u|²`, not an unspecified multiple.

Integrated first variation is `relativeAction` / `actionStationary_iff_eulerLagrange` (`Variation.lean:56-58,160-162`). Plane-wave reduction uses `u = a cos(k·x − ω t)` with `waveCovector = (-ω, k)` (`PlaneWave.lean:22-24`). The modal operator is

\[
M(ω,k)a = ρ ω² a − μ\bigl(|k|² a − k(k·a)\bigr),
\]

i.e. `M = ρ ω² I − μ(|k|² I − kkᵀ)` (`Action.lean:65-66`). Euler–Lagrange is `∂L/∂u − ∑_j ∂_j(∂L/∂J_j)` with `∂L/∂u = 0` (`PlaneWave.lean:81-83`). That sign is the same as both CAS residual constructors.

S10-to-S9 specialization is exact, not count-matching: `lagrangian_three`, `modalOperator_three`, `eulerLagrange_three`, `relativeAction_three`, `actionStationary_three`, and the mode-space equalities in `lean/s10/S10Pilot/Specialization.lean`. S10 MAIN stiffness reduces to `|curl|²` at D=3 by `stiffness_three`. Reuse of S10 controls is valid only through those identities; S10’s reviews did not inspect the original S9 constructors, which FIDELITY.md states.

**Census, full kernels, not determinant multiplicity.** On `ρ>0`, `μ>0`, `k≠0`:

| Case | Full amplitude kernel | Dimension | Anchor |
|---|---|---|---|
| `ω = 0` | `longitudinalSpace k = span{k}` | 1 | `zero_frequency_iff`, `longitudinalSpace_finrank` |
| `ω ≠ 0` and `ω² = (μ/ρ)\|k\|²` | `transverseSpace k = ker(k·)` | 2 | `stationary_on_cone_iff`, `transverseSpace_finrank` |
| `ω ≠ 0` and off-cone | `{0}` | 0 | corollary of `propagating_mode_iff` / `propagating_variational_mode_iff` for `a ≠ 0` |

These are subspaces of all real amplitudes, not witness vectors. `coneValue_pos` makes the cone strictly positive, so `ω=0` is disjoint from the cone. Both signs of nonzero frequency share `ω²`. The certificate exhibits a positive root; the iff covers any nonzero `ω`. Algebraic multiplicity `ρ s (ρ s − μ|k|²)²` agrees with dimensions 1+2 here; Lean’s objects are `finrank` of those kernels, and FIDELITY.md does not substitute the determinant for that.

**Source instrument limitation is honest.** `_measurements/S9_lean_source_check.py` execs selected original SymPy assigns/functions and compares the MAIN action in independent jet symbols and both native modal matrices to the handwritten Lean-shaped references, with exact zero residuals in the filed JSON. It also rejects the reversed shear sign in both action and operator; the recorded nonzero residuals are `μ|curl|²` and `2μ(|k|² I − kkᵀ)`, the expected FORM change. Wolfram is pinned by five unique native strings (coordinates, velocity, curl, MAIN Lagrangian, `Exp[I(k·x−ω t)]`); there is no fresh Wolfram run and no comparator/PIT claim. Counts are not the sole identification: the 3×3 operators and the jet action are compared.

---

## 3. C2 — Scalar phase velocity (only new claim)

The constitutive law is supplied: `v = (ħ/m) ∇θ` (`docs/model_map.md` flow form; S9 record’s weaker in-repo P2 form). Lean defines

```lean
velocity ... := fun i => (hbar / mass) * coordDeriv i.succ theta x
```

so the derivative is spatial (`Fin 3` → spacetime axes 1,2,3), not `∂_t`. The phase is the explicit real family

`θ = θ₀ + ε[A cos φ + B sin φ]`, `φ = k·x − ω t`, `θ₀` constant.

Coordinate calculus: `∂_j(A cos φ + B sin φ) = (−A sin φ + B cos φ) q_j` (`phasePerturbation_coordDeriv`). Spatial `q` is `k`. First variation in `ε` at 0 is an actual `deriv`, not a supplied amplitude (`linearVelocity`). The quadratures are

\[
\text{cosAmplitude} = \tfrac{\hbar}{m} B\, k, \qquad
\text{sinAmplitude} = -\tfrac{\hbar}{m} A\, k
\]

(`Madelung.lean:57-63`). Independent expansion of `v = (ħ/m) ε(−A sin φ + B cos φ) k` gives the same signs.

Both amplitudes lie in `span{k}` (`cosAmplitude_mem`, `sinAmplitude_mem`). For `ħ ≠ 0`, `m ≠ 0`, the cosine family is the whole longitudinal span (`cosAmplitude_range`); the sine family is the same span with a sign absorbed in `A`. For `k ≠ 0`, `span{k} ∩ transverseSpace k = {0}` (`longitudinal_transverse_eq_zero`). At `k = 0` the linear velocity is identically zero (`zero_wavevector`), and the nonzero-`k` dimension theorems are not applied. Nonvacuity: `concrete_longitudinal_velocity` at `(ħ,m,A,B,k,x)=(1,1,0,1,(0,0,1),0)` gives `(0,0,1) ≠ 0`. The `ħ=0` control shows the range is empty without a nonzero prefactor.

**Interpretation is not overstated.** The phase is `Point → ℝ`, hence single-valued; a global smooth real phase is the vortex-free chart. Uniform nonzero density is not a Lean hypothesis and is not manufactured (`COVERAGE.md:34-35`). Module header, FIDELITY.md:103-108, RESULT.md:150-153, and README.md:18-20 all refuse GNLS dynamics, Bogoliubov dispersion, one dynamical branch, Fourier/PDE completeness, a general no-spin-1/photon theorem, confinement, and microscopic emergence. Calling this the “narrow scalar-phase part of P2” with the broader P2 sentence excluded matches the S9 record (P2 remains a supplied premise; this is the kinematic `v=(ħ/m)∇θ` fragment). No sentence in the contract/proof modules was found that claims more than that.

---

## 4. C3 — Controls (records, not PASS labels)

Canonical S9 hashes in `S9_lean_contract_checks.json` match `MANIFEST.json` / `CLOSURE_VERIFICATION.txt` before and after. Instrument SHA-256 matches the pinned checker. Initial weak tests are retained in `S9_lean_contract_initial_attempt.json`; the later expansion-sign and concrete-count mutants are test-design changes, not proof repairs.

Inspected rejections are mathematical:

| Control | What failed | Why it is the intended reason |
|---|---|---|
| `action_shear_sign` | `lagrangian_variation` unsolved polynomial identity | `+μ/2\|curl\|²` vs claimed variation that still uses the minus bilinear form |
| `relative_action_expansion_sign` | `relativeAction_expansion`: `I = −I` | Wrong claimed expansion of the unchanged integral |
| `phase_wrong_coordinate` | `velocity_phasePerturbation`: rewrite fails; goal has `waveCovector … 0` (`−ω`) vs `k i` | Time derivative is not the spatial Madelung gradient |
| `phase_sine_sign` | `linearVelocity_eq`: `+ (ħ/m)A k sin` vs `− (ħ/m)A k sin` | Opposite sine quadrature |
| `transverse_count_mutant` / `longitudinal_count_mutant` | `2=3` / `1=2` → `False` | Wrong full-subspace dimensions at `k=(0,0,1)` |
| `zero_wavevector_domain_mutant` | `(0,0,1) ∈ span{0}` → `False` | Static census needs `k≠0` |
| `longitudinal_not_transverse_mutant` | unit cosine amplitude ∈ transverse → `False` | Nonvacuous exclusion |
| `nonzero_prefactor_required_for_range_mutant` | `ħ=0` amplitude `≠ 0` → `False` | Range needs nonzero prefactor |

Five matching positive fixtures compile. `phase_wrong_coordinate` also lints unused `i` under `warningAsError`; the counted diagnostic is still the rewrite/identity failure in `velocity_phasePerturbation`, not a syntax/import/timeout. Domain nonemptiness: `coneValue_pos`, concrete transverse/longitudinal action witnesses, and the nonzero Madelung example. Off-cone “kernel is `{0}`” is already the contrapositive of `propagating_mode_iff`; a separate mutant is unnecessary for clearance.

Source-check shear-sign rejection is an independent native-constructor FORM control. S10 mutation JSON in the packet is not treated as an S9 constructor check; only `Specialization.lean` licenses S10 reuse.

---

## 5. C4 and non-review completion

Filed sequential S9 rebuild (`-DwarningAsError=true`), 32 axiom audits (only `propext`, `Classical.choice`, `Quot.sound`), and no `sorry`/`admit`/custom physics axioms in canonical S9 sources. That is deduction hygiene, not fidelity. Compact source connection is documented and checked at the stated SymPy-execute / Wolfram-pin level. Policy L1 is respected: no systematic per-output CAS bridge. Broader S9 ledger, export/comparator, S11, and full GNLS remain excluded, as required.

---

## 6. Findings

No mathematical or statement-fidelity blockers.

Optional editorial (not required to clear): a named `off_cone_kernel` lemma wrapping the existing iff; a `sinAmplitude_range` twin of `cosAmplitude_range`; Wolfram route-B source anchors if a later pin is wanted. None of these changes what is proved or may be claimed.

---

## 7. Remaining limitations (in-scope, already declared)

- Wolfram is source-inspected, not re-executed; the handwritten SymPy/Lean references are a translation boundary, not a kernel-certified Python/Wolfram interpreter.
- Plane-wave ansatz only; no general PDE/Fourier completeness.
- Madelung law, uniform-density/vortex-free reading, and curl-only action remain supplied.
- This is not the full P2 / spin-1 / Bogoliubov / confinement claim, and not S9 ledger completion.
- S10 form/coefficient controls apply to S9 only through the D=3 identities in `Specialization.lean`.

**CLEAR** for the bounded S9 Lean contract at this packet revision.