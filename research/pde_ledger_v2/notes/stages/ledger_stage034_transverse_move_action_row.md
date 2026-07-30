# ledger_stage034_transverse_move_action_row

## Status

**Part V — Magnetism. V-1 (build-order 034; the FIRST stage of the 6-stage Part-V split,
user decision 2026-07-22).** The self-contained action-row cluster of the moving-throat (magnetism)
sector. Where Part IV built the *static* ±w throat as charge (030 substrate → 031 puncture-deflection
mechanism → 032 sign R1 → 033 the native-`P` departure), stage 034 asks the complementary *moving*-throat
question at the **action level only**: importing the pathA_36 transverse-vector row into G0 (which contains
**no** `u_T`) and supplying the one finite-profile moving-defect coupling, is the result a clean, stable
`(G0+δ)` amendment that damages no pre-existing G0 row? — and records the honest **YES**:

- **`TRANSVERSE_MOVE_ACTION_ROW`** — the moving-throat amendment is a well-formed transverse-vector action row

  ```text
  S_{T+move} = ∫ dt d³x [ ½ρ_br|u̇_T|²  −  ½μ_R|∇×u_T|²  +  q_T Σ_i s_i η_a(x−X_i) V_i·u_T ],
               ∇·u_T = 0,   c_γ² = μ_R/ρ_br,   q_T = λ_T τ_d,
  ```

  with `ρ_br, μ_R > 0` giving a **positive-definite transverse Hessian** (two propagating transverse
  polarizations, dispersion `ω² = c_γ²k²`, no ghost, no tachyon) and **no pre-existing G0-v0 row changed** —
  the amendment activates ONLY the previously-absent transverse DOF and the finite-profile moving coupling;
  scalar, drain, return (`F_flux`), wall, geon, and every other declared-zero G0 entry is byte-for-byte
  untouched. Verdict token **`TRANSVERSE_MOVE_ACTION_ROW`** with **`internal_inconsistency = none`** (both
  engines, both exit 0).

**⭐ Scope class: EARNED (action row; IMPORTS pathA_36 — NOT a from-scratch derivation).** This is an
EARNED action-structure result *within* the postulated G0 closure: it earns that the imported transverse row
plus one new coupling compose into a consistent, stable, non-damaging amendment. It is **NOT** a knob,
**NOT** a reduction/codimension edge, and it does **not** shrink the irreducible count. It also carries no
electric branch/sign — the moving row inherits the electric sector's open `A_E`, and `q_T`'s value is the
sim-deferred nonlinear-throat R1.

**⚠ CRITICAL TYPING (import vs earn).** The `u_T` kinetic/gradient row `½ρ_br|u̇_T|² − ½μ_R|∇×u_T|²` and the
light speed `c_γ² = μ_R/ρ_br` are **IMPORTED from Part III / stage003 / pathA_36** (the brane transverse-vector
row). Stage 034 re-instantiates that row's structure for its own dual-engine check but does **NOT** claim it as
a 034 earn and does **NOT** re-count `{u_T, μ_R, c_γ, ρ_br}` in the irreducible set (de-dup deferred to Part
VII). The magnetism-NEW content of 034 is **exactly** `{the moving coupling, q_T, τ_d}`.

**⭐ This is the FIRST stage of Part V.** After 034, magnetism is only *begun*: the native source law, the
parity census, both far-field routes, the boost-consistency comparison, the sealed landing, and the `b_T`
T-even departure are all **downstream** (035–039; see §8). 034 is the **action-row gate only** — it writes the
coupling term but does not derive the source from continuity, does not run the parity census, and does not
compute any far-field force.

## Purpose

Record, as an EARNED action-structure result, that magnetism enters the ledger as a clean `(G0+δ)`
transverse-vector amendment. The decisive object is the quadratic transverse action density: importing the
pathA_36 kinetic/gradient row into a G0 that has no `u_T`, supplying the one finite-profile moving-defect
coupling `q_T Σ_i s_i η_a V_i·u_T`, and certifying at runtime that (i) the imported row is faithfully
instantiated with the transverse dispersion `ω² = c_γ²k²`, (ii) the coupling is the correct finite-profile
current·field term with `q_T = λ_T τ_d`, (iii) the transverse Hessian is positive-definite (two polarizations,
no ghost/tachyon), (iv) no pre-existing G0 row changes (`internal_inconsistency = none`, `F_flux` explicitly
untouched), (v) the row is well-formed (variational, local, linear in `u_T`), and (vi) the full dimensional
identity closes with each action density term `[EL⁻³]` and `[S_{T+move}] = [E·T]`. Both engines run the same
six action-row teeth plus the build-global target-blindness / dual-engine / units-restored guards, the
import-vs-earn anti-overclaim guard, the computed-verdict re-derivation, and the source-to-stage manifest —
and each independently reaches `TRANSVERSE_MOVE_ACTION_ROW`.

Consumes **nothing numeric** from 030/031/032/033 — it is an independent action-row construction that CITES
the pathA_36 transverse row as provenance (see §10). It is the moving-throat twin of the static charge action
of Part IV.

## 1. The imported pathA_36 transverse row and the transverse dispersion (`ACTION_KINETIC`)

G0-v0 contains **no** transverse shear DOF `u_T`. The kinetic + gradient row supplied to the amendment,

```text
L_kin+grad = ½ρ_br|u̇_T|²  −  ½μ_R|∇×u_T|²,        ∇·u_T = 0,
```

is **imported verbatim from stage003 / pathA_36** (the brane transverse-vector row). The tooth checks that the
import is faithfully instantiated with the correct coefficients (`ρ_br/2` kinetic, `−μ_R/2` gradient), NOT that
034 earns it. For a divergence-free field the free Euler–Lagrange equation is

```text
ρ_br ü_T = −μ_R ∇×(∇×u_T) = μ_R ∇²u_T        (since ∇×∇×u_T = ∇(∇·u_T) − ∇²u_T = −∇²u_T for ∇·u_T=0),
```

so a plane wave `u_T = ε e^{i(k·x − ωt)}` with `k·ε = 0` (the two `∇·u_T = 0` transverse polarizations) obeys

```text
ρ_br ω² = μ_R k²   ⇒   ω² = (μ_R/ρ_br) k² = c_γ² k²,        c_γ² = μ_R/ρ_br.
```

The transverse dispersion `ω² = c_γ²k²` follows for both polarizations, at the brane light speed `c_γ`
imported from pathA_36. (Perturbing a coefficient — `ρ_br/2 → ρ_br/3`, flipping the gradient sign, or breaking
`c_γ² = μ_R/ρ_br` — makes the dispersion/coefficient assert fire at its own tooth.)

## 2. The magnetism-NEW moving-defect coupling (`ACTION_COUPLING`)

The one term stage 034 **supplies** is the finite-profile moving-defect coupling

```text
L_coupling = q_T Σ_i s_i η_a(x − X_i) V_i · u_T,        q_T = λ_T τ_d,
```

where `s_i = ±1` is the throat's `±w` orientation (the charge-sector sign), `η_a` is the normalized finite
mouth profile (`∫d³x η_a = 1`, `[η_a] = L⁻³`), `V_i` is the throat velocity, and `q_T = λ_T τ_d` is the
moving-throat charge built from the coupling strength `λ_T` and the active drain's time-arrow `τ_d` (the full
time-reverse maps a drain to a source). This is a proper current·field source: linear in `u_T`, with the
per-defect transverse current `J_{T,i} = q_T s_i η_a V_i`. It is **target-blind** — it carries only `q_T`, no
electric `A_E`, no `q/g` knob, and no downstream sign/landing token. (Ablating any factor — drop `V_i`, drop
`s_i`, drop `η_a`, replace `q_T` by a wrong product, or drop `τ_d` so `q_T ≠ λ_T τ_d` — fires the coupling-form
assert.)

## 3. Transverse Hessian stability — PD, two polarizations, no ghost/tachyon (`ACTION_STABILITY`)

Stability is a **computed positivity object**, not a prose claim. Per transverse component the kinetic Hessian
is `∂²L/∂u̇_T² = ρ_br I` and the Fourier-space gradient (potential) Hessian is `μ_R k² P^T` (with `P^T` the
transverse projector onto `k·ε = 0`):

```text
kinetic  Hessian PD  ⇔  ρ_br > 0    (no ghost),
gradient Hessian PSD ⇔  μ_R  > 0    (no tachyon: ω² = c_γ²k² ≥ 0).
```

With `ρ_br, μ_R > 0` both are positive-(semi)definite on the two-dimensional transverse subspace, so the row
carries **exactly two propagating transverse polarizations**, both stable — no ghost (kinetic PD) and no
tachyon (gradient PSD, `ω² ≥ 0`). (Setting `ρ_br < 0` or `μ_R < 0` makes the eigenvalue/dispersion object
actually change sign — non-PD Hessian — and the stability assert fires; it is not a prose flip.)

## 4. No G0 damage — the amendment activates only the absent transverse DOF (`G0_DAMAGE`)

The amendment is computed as a diff/overlap object over the parsed G0 row set. It touches **only** the
previously-absent transverse DOF `u_T` and the finite-profile moving coupling; every other declared-zero G0
entry — scalar, drain, **return `F_flux`**, wall `r_B`, geon, and all remaining rows — is byte-for-byte
unchanged. The diff over the pre-existing rows is empty, so `internal_inconsistency = none`: magnetism coexists
with the gravity/charge substrate without contradiction. **⚠ The active `F_flux` ledger is NOT silently
absorbed into this conservative row** — it is named as untouched and explicitly deferred (its `O(V₁V₂)`
contribution and full-force integrability are R1, handled at V-5 / Part VII). (Mutating a declared-zero G0
entry, or letting the moving coupling leak into a pre-existing row — including a silent `F_flux` alteration —
makes the diff recompute a nonzero damage and the `internal_inconsistency = none` assert fires.)

## 5. The ledger-ready `(G0+δ)` row (`LEDGER_READY_ROW`)

The imported kinetic/gradient row and the new coupling compose into a single consistent action density that is
a well-formed `(G0+δ)` amendment: **variational** (it is an action density varied to give the transverse
Euler equation), **local** (finite-profile `η_a`, no nonlocal kernel), and with the coupling **linear in
`u_T`** (a proper current·field source, not a nonlinear self-interaction). The full row is

```text
S_{T+move} = ∫ dt d³x [ ½ρ_br|u̇_T|²  −  ½μ_R|∇×u_T|²  +  q_T Σ_i s_i η_a(x−X_i) V_i·u_T ],   ∇·u_T = 0.
```

(Breaking well-formedness — making the coupling nonlinear in `u_T`, non-local, or non-variational — fires the
tooth.)

## 6. The dimensional identity — the eight dims (`FIELD_IDENTITY_UNITS`)

Units RESTORED (natural-units-independent, able-to-fail). The preregistered field is the transverse shear
**displacement** `u_T` (`[u_T] = L`); the curl candidate `b_T = ∇×u_T` is **dimensionless** (its
parity/departure is stage 039, NOT here — 034 records only its dimension). The eight dimensions:

```text
[u_T]     = L
[u̇_T]     = L T⁻¹
[∇×u_T]   = 1              (dimensionless)
[ρ_br]    = M L⁻³
[μ_R]     = M L⁻¹ T⁻²      (= E L⁻³)
[q_T]     = M T⁻¹
[η_a]     = L⁻³
[b_T]     = 1              (= [∇×u_T])
```

Each action density term is `[E L⁻³]`:

```text
½ρ_br|u̇_T|²        : [M L⁻³][L T⁻¹]²          = M L⁻¹ T⁻² = E L⁻³
½μ_R|∇×u_T|²        : [M L⁻¹ T⁻²][1]            = M L⁻¹ T⁻² = E L⁻³
q_T s η_a V·u_T     : [M T⁻¹][1][L⁻³][L T⁻¹][L] = M L⁻¹ T⁻² = E L⁻³
```

and `∫ dt d³x` over an `[E L⁻³]` density gives `[T · L³ · E L⁻³] = [E · T]`, so `[S_{T+move}] = [E·T]`
(action). (Corrupting one dimension — e.g. `[q_T] = MT⁻¹ → MT⁻²`, or `[b_T] = 1 → L⁻¹` — makes the
dimensional-homogeneity object inhomogeneous and the units tooth fires; free-carrier-independent.)

## 7. Provenance-cite discipline: import vs earn (`GUARD_IMPORT_VS_EARN`)

An anti-overclaim accounting object protects the don't-re-count discipline at runtime:

```text
earned_new_terms  ==  { moving_coupling }                     (the magnetism-NEW content of 034),
imported_terms    ⊇  { u_T_kinetic, u_T_gradient, c_γ }       (cited from stage003 / pathA_36).
```

The imported pathA_36 row (`u_T` kinetic, gradient, `c_γ`, `μ_R`, `ρ_br`; `c_γ² = μ_R/ρ_br`) is **TAGGED as
provenance**, not counted as a 034 earn. The magnetism-NEW content is exactly `{moving coupling, q_T, τ_d}`.
(Mutating the accounting to wrongly claim the imported kinetic row as a 034 earn fires the guard.)

## 8. Scope — what 034 does NOT do (deferred to 035–039)

034 is the **action-row gate only**. Named explicitly as OUT OF SCOPE:

- **The native source law** `J_T = q_T s η V` from defect-continuity (the `CONVECTION_LIKE_CONDITIONAL`
  current identity) **plus the full parity census** (`s`, `V`, `τ_d`, `q_T`, `J_T`, `u_T`, `b_T` under
  `R_w`/`P_w`/rotations/time-reversal) is **stage 035 (V-2)**. 034 writes the coupling term but does NOT derive
  the source from continuity and does NOT compute the parity census.
- **Route A (Maxwell–Darwin kernel)**, **Route B (direct shear, blind to A)**, and the **boost-consistency
  comparison** (`r_BA`, `δ_BA`, `r_cone`, `ΔU`) are **stages 036/037 (V-3/V-4)**.
- **The SEALED §4 first-match landing** — `R1_REQUIRED(electric_bc_selection)` + the three co-blockers
  (`R1_REQUIRED(direct_moving_throat)` / `(magnitude)` / `(consistency)`), the 1152-cell truth table, and the
  active-`F_flux` caveat — is **stage 038 (V-5)**.
- **The `b_T = ∇×u_T` time-reversal-EVEN departure** (`B_TIME_REVERSAL_EVEN`, the "magnetism is not exact
  Maxwell" characterization: `b_T` is T-even where a physical Maxwell `B` is T-odd) is **stage 039 (V-6)**.
  034 records only `b_T`'s dimension, NOT its parity or the Maxwell comparison.
- **The imported pathA_36/stage003 row** (`u_T` kinetic, `c_γ`, `μ_R`, `ρ_br`) is **CITED provenance** — NOT
  re-derived, NOT re-counted in the irreducible parameter set (de-dup deferred to Part VII).
- **`q_T`'s value / the throat normalization** is the sim-deferred nonlinear-throat R1. 034 introduces
  `q_T = λ_T τ_d` as a FREE-UNREDUCED / R1 parameter but does NOT resolve it (the shared Part-VII throat
  solve); this R1 is kept SEPARATE from `R1_REQUIRED(electric_bc_selection)` (the electric `A_E` branch/sign).

## 9. Source-to-stage predicate manifest

Completeness certificate: **no silently-dropped source claim.** Every source tooth in scope from the source
build (`magnetism_moving_throat_check.{py,wl}`, the 35-tooth `TOOTH_ORDER`, commit `53cf049f`) lands as
**PRESERVED** (folded as-is), **REPLACED_BY_STRONGER** (a stronger reconstruction tooth), or **SCOPED_OUT**
(deferred to its own Part-V stage, each named LITERALLY — NO wildcards). The partition is disjoint +
exhaustive, computed at runtime in **both** engines from the same canonical `(id, disposition, owner)` triples:

```text
partition = { PRESERVED: 1, REPLACED_BY_STRONGER: 8, SCOPED_OUT: 26 }   (35 total),
manifest digest (SHA-256) = 4343bd60cd974f653a0a8ac2eeced6c7aca15b1831c81be20bd358f449c454af.
```

- **PRESERVED (1):** `TARGET_BLINDNESS` (folded as the build-global source-side-dependencies guard).
- **REPLACED_BY_STRONGER (8):** the six action-row teeth `FIELD_IDENTITY_UNITS` (→ dimension object),
  `ACTION_KINETIC` (→ Hessian/dispersion), `ACTION_COUPLING` (→ differentiated source), `ACTION_STABILITY`
  (→ PD Hessians), `G0_DAMAGE` (→ parsed-G0 diff), `LEDGER_READY_ROW` (→ local/variational row), plus
  `UNITS_RESTORED` (→ whole-density firewall) and `DUAL_ENGINE_TERMS` (→ canonical term inventory).
- **SCOPED_OUT (26)** — the source teeth of V-2…V-5, each named verbatim (same partition both engines):
  - **V-2 (035):** `SOURCE_TRANSLATION_CONTINUITY`, `SOURCE_NOT_IMPORTED`, `SOURCE_BASIS`, `PARITY_RW`,
    `PARITY_PW`, `PARITY_ROTATION`, `PARITY_TIME_REVERSAL` (7)
  - **V-3 (036):** `BOOST_PROJECTOR`, `BOOST_GENERAL_VELOCITIES`, `BOOST_NEXT_ORDER` (3)
  - **V-4 (037):** `DIRECT_SOURCE`, `DIRECT_PROJECTOR`, `DIRECT_EXCHANGE_SIGN`, `DIRECT_FALLOFF`,
    `DIRECT_VELOCITY_ORDER`, `ROUTE_INDEPENDENCE`, `BOOST_COMMON_VELOCITY`, `COMPARE_COMPUTED`, `DELTA_RATIO`,
    `CONE_RATIO`, `QMAG_R1` (11)
  - **V-5 (038):** `TRUTH_TOTALITY`, `TRUTH_PRECEDENCE`, `LANDING_OWNERSHIP`, `ACTIVE_FLUX_CAVEAT`,
    `HOOK_LORENTZ` (5)

  **V-6 (039) owns NO additional source-build tooth** — it reuses V-2's parity results (`PARITY_TIME_REVERSAL`
  + `PARITY_ROTATION`) and authors its own NEW `B_TIME_REVERSAL_EVEN` departure assertions, cited-not-owned.

Both engines assert the identifier set (35 unique), the three-way disposition set, the exact counts, and the
committed digest at runtime and agree (`SOURCE_TO_STAGE_MANIFEST`; dropping one row fires it).

## 10. Consumes / cites

- **Consumes NOTHING numeric from 030/031/032/033.** 034 does not touch the puncture-deflection response matrix
  `m`, `m_gg`, the four `A_X`, the BC-ensemble data, or the native-`P` constraint signatures — it is an
  independent action-row construction.
- **Cites (provenance, NOT consumed, NOT re-counted):** the transverse-vector row `½ρ_br|u̇_T|² − ½μ_R|∇×u_T|²`,
  the field `u_T` (`∇·u_T = 0`, `[u_T] = L`), and the light speed `c_γ² = μ_R/ρ_br` — all from **Part III /
  stage003 / pathA_36** (`ledger_stage003_transverse_photons_stray_longitudinal`). De-dup of
  `{u_T, μ_R, c_γ, ρ_br}` against Part III is deferred to Part VII; they are NOT re-added to the irreducible
  count here. Also cites the charge-sector `±w` orientation `s_i` (Part IV) as the coupling's sign carrier.

## Verification

- **Dual-engine, both exit 0, genuinely independent routes.**
  `scripts/ledger_stage034_transverse_move_action_row_sympy_audit.py` — **SymPy 12 teeth**.
  `mathematica/ledger_stage034_transverse_move_action_row_mathematica_audit.wl` — **Mathematica 12 teeth**.
  Standalone, print-only, assert-zero (`raise SystemExit(1)` / `Exit[1]`), no argparse harness, no JSON/YAML
  payload, **zero file-I/O between engines**. Each engine reaches `TRANSVERSE_MOVE_ACTION_ROW` (with
  `internal_inconsistency = none`) on its own and prints its own tokens — cross-engine agreement is that they
  independently produce the same tokens, not a compare pass.
- **The `.wl` is a genuinely INDEPENDENT route** — materially distinct decomposition, not a line-port of the
  `.py`. The **SymPy** route checks transverse stability via the explicit kinetic/gradient Hessians and a
  generalized-eigenvalue dispersion; the **Mathematica** route is independent — a `NullSpace` transverse
  Fourier basis, `CoefficientArrays` quadratic forms, `PositiveDefiniteMatrixQ`, and native dispersion solves;
  each engine differentiates/decomposes its OWN action density. **No dual-engine disagreement.**
- **Per-tooth ablation** (env switch `LEDGER_STAGE034_MUTATION`): all **12 mutations FIRED_AT_OWN_ASSERT per
  engine** — the six action-row teeth, the build-global `TARGET_BLINDNESS` / `DUAL_ENGINE_TERMS` /
  `UNITS_RESTORED`, the `GUARD_IMPORT_VS_EARN` anti-overclaim guard, the `VERDICT_REDERIVATION` tooth, and the
  `SOURCE_TO_STAGE_MANIFEST` tooth. The **verdict tooth is non-tautological**: it mutates a **COMPUTED** object
  (`rho_sign = −1` → a non-PD transverse Hessian) so the pipeline **RE-DERIVES** the named token `ROW_UNSTABLE`
  (never a literal flip of `TRANSVERSE_MOVE_ACTION_ROW` — the 030 `X≡X` lesson). `GUARD_IMPORT_VS_EARN` fires
  when the imported pathA_36 kinetic row is (wrongly) claimed as a 034 earn.
- **Tri-review outcome (falsification-first — recorded transparently, not hidden).**
  - **FIDELITY:** found **1 NIT** — a dead `X≡X` clause in the `.wl` `dimensionObject` that returned
    hardcoded dimensions and then asserted them equal to themselves. **REMEDIATED:** the `.wl` now computes the
    density/action dimensions via **logarithmic derivatives of its own independently-rescaled density/action**,
    so the `.wl` `FIELD_IDENTITY_UNITS` mutation now genuinely fires → **CLEAN**.
  - **ADVERSARIAL: CLEAN with 2 DOCUMENTED non-blocking coverage notes** (verified-safe, NOT remediated):
    1. Inside `ACTION_KINETIC` the sub-conjuncts `transverse_checks == (0,0)` and `polarization_count == 2`
       are **structurally true-by-construction** (the orthogonal transverse basis / `NullSpace` construction) —
       they cannot fail individually, but they **ride a genuinely able-to-fail coefficient/dispersion assert**,
       so the tooth is not vacuous.
    2. The compound teeth carry **one mutation covering their sub-checks** (drop-`s` / drop-`V` / drop-`η_a`,
       the profile-integral, the gradient-sign) — within the directive's permissive "e.g." per-tooth-ablation
       scope.

  Arbiter re-runs of both engines reproduce the tokens and the manifest digest; the tri-review leaves the
  action row EARNED and `internal_inconsistency = none`.

## Downstream consumers

- **Part V continues** after 034: 035 (native source law + parity census) → 036/037 (Route A / Route B +
  boost-consistency comparison) → 038 (the sealed `R1_REQUIRED(electric_bc_selection)` landing) → 039 (the
  `b_T` T-even departure). 034 is the action-row foundation those stages consume by citation.
- **Parameter register:** edge **R67** (`q_T = λ_T τ_d`, moving-throat charge, `[MT⁻¹]`, FREE-UNREDUCED / R1 —
  the sim-deferred nonlinear-throat solve, kept SEPARATE from `R1_REQUIRED(electric_bc_selection)`) + edge
  **R68** (`τ_d`, active-drain time-arrow, structural/postulated, `T-odd`, not a reducible knob) + the
  **cited-not-counted** provenance edges (`u_T`, `μ_R`, `c_γ`, `ρ_br` from stage003/pathA_36; `c_γ² = μ_R/ρ_br`)
  with an explicit anti-double-count note, and `internal_inconsistency = none` as a **structural fact** (no
  `[L,T,M]`, no knob).
- **`docs/model_map.md` §3.5:** the Q-CURRENT action-row bullet — the moving throat enters as a clean stable
  transverse-vector `(G0+δ)` amendment; sign + magnitude R1, `b_T` T-even departure downstream.
- **Part VII:** the imported `{u_T, μ_R, c_γ, ρ_br}` edges enter the de-dup (counted ONCE, with Part III); the
  `q_T` / `τ_d` obligations enter the honest R1 debt alongside the electric `A_E`.

## Provenance

- **Physics source:** `software/em_charge_attribute/magnetism_moving_throat_result.md` (Deliverable #1 — the
  ledger-ready `(G0+δ)` row — plus Appendix A units) + `software/em_charge_attribute/magnetism_moving_throat_check.{py,wl}`
  (the ONE dual-engine source build, 35 teeth, commit `53cf049f`). Stage 034 extracts the **action-row cluster**
  into a focused standalone dual-engine pair; the source `argparse`/`--compare` harness, the
  `--ablate-tooth`/`--out-dir` payload plumbing, and all JSON/log file writing are STRIPPED (print-only /
  zero-file-I/O / independent-tokens contract). The stage `.wl` is a materially distinct Wolfram route, not a
  line-port of the `.py`.
- **Consumes:** nothing numeric from 030/031/032/033 (independent action-row construction).
- **Cites (provenance):** the pathA_36 / stage003 transverse-vector row (`u_T` kinetic/gradient, `c_γ`, `μ_R`,
  `ρ_br`; `c_γ² = μ_R/ρ_br`) — NOT re-derived, NOT re-counted; the Part IV `±w` orientation `s_i`.
- **Governing:** `notes/ledger_v2_blueprint.md` §5 (reshape spec) + §6 (per-tooth ablation);
  `notes/part5_magnetism_atomic_split.md` (V-1 = the action row; the tooth-allocation table);
  `docs/model_map.md` §3.5 (Q-CURRENT action-row bullet). Reshape directive + review trail:
  `research/pde_ledger_v2/_scratch/stage034_reshape_directive.md`. ⛔ **Not retained** — it lived in gitignored `_scratch/` and no copy survives; this line records that a directive existed, it is not an auditable citation.
