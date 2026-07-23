# ledger_stage036_route_a_maxwell_darwin_reference

## Status

**Part V — Magnetism. V-3 (build-order 036; the THIRD stage of the 6-stage Part-V split,
user decision 2026-07-22).** The self-contained Route-A **Maxwell–Darwin reference-kernel** cluster of
the moving-throat (magnetism) sector. Where stage 034 (V-1) built the *action row* and stage 035 (V-2)
derived the *native source law*, stage 036 (V-3) asks the *far-field* question along the **boost route**:
what interaction does the moving ±w throat produce when you **boost the electric interaction**
`U_E = s₁s₂A_E/4πR` to `O(v²/c_γ²)`? — and records the honest answer:

- **`MAXWELL_DARWIN_REFERENCE`** — inverse-Fourier-transforming the Coulomb-gauge transverse projector
  reconstructs the **Maxwell–Darwin transverse Darwin kernel** `I_ij(R) = (δ_ij + n_in_j)/8πR`, and
  Lorentz-completing the electric anchor gives the general independent-velocity interaction `U_A` to
  `O(v²/c_γ²)`. This is the **Maxwell-consistent REFERENCE** the direct moving-throat route (037) is
  tested against.

  ```text
  F⁻¹[ P^T_ij / k² ] = δ_ij/4πR − F⁻¹[ k_ik_j/k⁴ ] = (δ_ij + n_in_j)/8πR ≡ I_ij(R),
  U_A = (s₁s₂A_E/4πR)·[ 1 − (D_V + A_V)/(2c_γ²) ] + O(v⁴/c_γ⁴),
  D_V = V₁·V₂,   A_V = (V₁·n)(V₂·n),   n = (X₂−X₁)/R,   R = |X₂−X₁|.
  ```

  Verdict token **`MAXWELL_DARWIN_REFERENCE`** (both engines, both exit 0), emitted at
  **`tier_A_conditional`**.

**⭐ Scope class: EARNED (reference kernel), tier-A conditional.** This is an EARNED result — the
inverse-FT reconstruction is an exact target-blind structural derivation of the Maxwell-consistent
transverse Darwin kernel, and the Lorentz-completion is exact through `O(v²/c_γ²)`. It is **conditional**
because Route A is a **boost of the electric interaction**: the kernel and anchor **INHERIT the electric
`A_E` `R1_REQUIRED(bc_selection)`** (sign + normalization unselected; stage 032, register R63). The
kernel's STRUCTURE is exact; its overall coefficient rides the open electric `A_E`. Hence the scope is
EARNED-but-conditional (`tier_A_conditional`), **NOT** GAPS-CLOSED. It is **NOT** a knob, **NOT** a
reduction/codimension edge, and it does **not** shrink the irreducible count — the reference kernel
introduces NO new free knob (it rides `A_E` + `c_γ`).

**⚠ CRITICAL SCOPE (036 = the reference kernel ONLY, NOT the boost-consistency result).** Stage 036
produces the Maxwell-consistent **reference** the direct route is tested against — it does **NOT** compute
the direct moving-throat shear (Route B), does **NOT** run the structural COMPARISON, and does **NOT**
claim "magnetism structurally IS the boost of electricity." That boost-consistency / emergent-Lorentz
result requires the second independent route (Route B) and the comparison, and is **stage 037 (V-4)**.
Route A's verdict is the reference-kernel identity ONLY. A reader must never read 036 as claiming emergent
Lorentz.

**⚠ CRITICAL TYPING (import vs earn).** The **electric interaction** `U_E = s₁s₂A_E/4πR` (with its full
far-field coefficient `A_E`, which carries the committed `R1_REQUIRED(bc_selection)` — sign + normalization
unselected; stage 032, register R63) and the **light-cone speed** `c_γ` (`c_γ² = μ_R/ρ_br`; Part III /
stage003, re-cited via stage 034) are **CITED provenance**. Stage 036 **re-instantiates** these imported
symbols for its own dual-engine check but does **NOT** re-derive the electric sector, does **NOT** re-count
`A_E` / `c_γ` in the irreducible set, and does **NOT** re-derive the native source law (035) — **Route A
never touches `J_T`; it boosts the electric interaction.** The magnetism-NEW content of 036 is **exactly**
`{the transverse Darwin reference kernel, the Lorentz-completed O(v²/c_γ²) anchor}`.

## Purpose

Record, as an EARNED (tier-A conditional) reference-kernel result, that boosting the electric interaction
reconstructs the Maxwell–Darwin transverse Darwin kernel and the Lorentz-completed velocity anchor. The
decisive object is the **inverse-FT reconstruction** of the Coulomb-gauge transverse projector: the
transverse projector `P^T_ij(k) = δ_ij − k_ik_j/k²` divided by `k²` inverse-Fourier-transforms — via the
radial seed `F⁻¹[k⁻⁴] = −R/8π` twice-differentiated — into the transverse Darwin kernel
`(δ_ij + n_in_j)/8πR`; contracting that kernel with the independent throat velocities reconstructs the
Lorentz-completed anchor `U_A` through `O(v²/c_γ²)`, the next order `O(v⁴/c_γ⁴)` left named-but-uncomputed.
The stage certifies at runtime that (i) the reconstructed kernel equals `(δ_ij + n_in_j)/8πR`
**component-by-component** (`BOOST_PROJECTOR`), (ii) the Lorentz-completed anchor equals
`EXPECTED_A2 = −s₁s₂A_E(D_V+A_V)/(8πc_γ²R)` and reproduces the full anchor / force / radial expressions
(`BOOST_GENERAL_VELOCITIES`), (iii) the computed velocity order is exactly 2 with `O(v⁴/c_γ⁴)` named but
NOT computed (`BOOST_NEXT_ORDER`), and (iv) the whole-stage dimensional firewall closes on the real
expressions (`UNITS_RESTORED`). Both engines run these three Q-BOOST teeth plus the build-global
target-blindness / dual-engine / units-restored guards, the computed-verdict re-derivation, and the
source-to-stage manifest — and each independently reaches `MAXWELL_DARWIN_REFERENCE` at
`tier_A_conditional`.

Consumes **nothing numeric** from 030/031/032/033/034/035 — it is an independent boost-of-electric
reconstruction that CITES `A_E` (Part IV) and `c_γ` (Part III / stage003) as provenance (see §9). It is
the far-field twin of the electric interaction, boosted.

## 1. The transverse Darwin kernel — inverse-FT of the transverse projector (`BOOST_PROJECTOR`)

Route A boosts the electric interaction `U_E = s₁s₂A_E/4πR`. The Coulomb-gauge transverse projector in
momentum space is `P^T_ij(k) = δ_ij − k_ik_j/k²`. Inverse-Fourier-transforming `P^T_ij/k²` gives the
transverse Darwin kernel:

```text
F⁻¹[ P^T_ij / k² ] = F⁻¹[ δ_ij/k² − k_ik_j/k⁴ ]
                   = δ_ij/4πR − F⁻¹[ k_ik_j/k⁴ ].
```

The first term is the Coulomb kernel `F⁻¹[1/k²] = 1/4πR`. The second is obtained from the radial seed

```text
F⁻¹[ k⁻⁴ ] = −R/8π,
```

using `k_ik_j ↦ −∂_i∂_j` in position space, i.e. `F⁻¹[k_ik_j/k⁴] = −∂_i∂_j(−R/8π) = ∂_i∂_j(R/8π)`. With
`∂_iR = n_i` and `∂_i∂_jR = (δ_ij − n_in_j)/R`,

```text
F⁻¹[ k_ik_j/k⁴ ] = (δ_ij − n_in_j)/8πR.
```

Subtracting from `δ_ij/4πR`:

```text
F⁻¹[ P^T_ij/k² ] = δ_ij/4πR − (δ_ij − n_in_j)/8πR = (δ_ij + n_in_j)/8πR ≡ I_ij(R),
```

where `n = (X₂−X₁)/R`, `R = |X₂−X₁|`. This is the **Maxwell-consistent (Coulomb-gauge / Darwin) reference
kernel.** The reconstruction is certified as a **component-wise zero-difference object**: for every `i,j`,
the reconstructed `transverse[i,j]` minus `(δ_ij + n_in_j)/8πR` is symbolically zero (with
`n_z² = 1 − n_x² − n_y²` substituted, matching the build's unit-sphere reduction). ⚠ This is the primary
computed object the verdict re-derives from — the reference kernel's **coefficient `1/8π` and its
transverse `(δ_ij + n_in_j)` structure** are both load-bearing. (Doubling the reconstructed kernel
`I → 2I` — the `BOOST_PROJECTOR` mutation — or breaking transversality or corrupting the coefficient
`1/8π → 1/4π` makes the component-wise zero-difference recompute NONZERO and the tooth fires at its own
assert.)

## 2. The Lorentz-completed velocity anchor (`BOOST_GENERAL_VELOCITIES`)

Contracting the transverse Darwin kernel with the independent throat velocities `(V₁, V₂)` reconstructs
the `O(v²/c_γ²)` interaction. With `D_V = V₁·V₂` and `A_V = (V₁·n)(V₂·n)`, the contraction
`V₁ᵀ I V₂ = (D_V + A_V)/8πR`, so the `O(v²/c_γ²)` interaction piece is

```text
U_{A,2} = −(s₁s₂A_E/c_γ²)·(V₁ᵀ I V₂) = −s₁s₂A_E(D_V + A_V)/(8πc_γ²R)     (the build's EXPECTED_A2),
```

and the full Lorentz-completed anchor, force, and radial force are

```text
U_A     = (s₁s₂A_E/4πR)·[ 1 − (D_V + A_V)/(2c_γ²) ] + O(v⁴/c_γ⁴),
F_{A,2} = (s₁s₂A_E/8πc_γ²R²)·[ (V₂·n)V₁ + (V₁·n)V₂ − (D_V + 3A_V) n ],
F_{A,2,r} = − s₁s₂A_E(D_V + A_V)/(8πc_γ²R²).
```

The tooth is certified as a **symbolic zero-difference** `route_a_u − EXPECTED_A2 = 0`; the full anchor,
force, and radial expressions are reproduced and dimension-checked. (Dropping the `A_V` term
`route_a_u.subs(A_V, 0)` — the `BOOST_GENERAL_VELOCITIES` mutation — makes the zero-difference recompute
nonzero and the tooth fires.)

## 3. The order structure — honest truncation (`BOOST_NEXT_ORDER`)

The `O(1)` and `O(v²/c_γ²)` orders are **explicit** (the two terms of `U_A` above); the next term
`O(v⁴/c_γ⁴)` is **named but NOT computed** — the anchor's velocity order is exactly 2, the next order left
uncomputed as an honest R1-adjacent remainder (owned downstream at 038's higher-order hook, NOT a new
knob). The tooth is certified as a **runtime order object** `velocity_order == 2`. (Setting the claimed
order to 4 `boost_order = 4` — the `BOOST_NEXT_ORDER` mutation — makes the `velocity_order == 2` object
False and the tooth fires; this proves the tooth pins the honest truncation — 036 does NOT compute the
`O(v⁴)` term.)

## 4. Why `MAXWELL_DARWIN_REFERENCE` — a reference kernel, tier-A conditional

The kernel is **fully reconstructed in tensor form** — the transverse Darwin `(δ_ij + n_in_j)/8πR` — and
the Lorentz-completion is exact through `O(v²/c_γ²)`. But the **overall coefficient** is conditional:
Route A boosts the electric interaction `U_E = s₁s₂A_E/4πR`, whose `A_E` carries the committed
`R1_REQUIRED(bc_selection)` (sign + normalization unselected; stage 032, R63). The kernel and anchor
**INHERIT that R1**. Hence the verdict is `MAXWELL_DARWIN_REFERENCE` emitted at **`tier_A_conditional`** —
"Maxwell–Darwin reference kernel" (the structure is the transverse Darwin projector), "tier-A conditional"
(its coefficient rides the unselected electric `A_E`).

**⚠ The open electric `A_E` R1 is a scope tag, not a certification failure.** The unresolved electric
`A_E` sign/magnitude does NOT change the reference-kernel token — production is `MAXWELL_DARWIN_REFERENCE`
at `tier_A_conditional` REGARDLESS of the open `A_E`. An unresolved anchor coefficient is NOT a broken
kernel, so it must NOT re-derive to the `MAXWELL_DARWIN_REFERENCE_UNCERTIFIED` alternate (see Verification). The
tier-A conditionality is a scope tag on the production verdict, not a failure mode.

## 5. The dimensional identity — the whole-stage firewall (`UNITS_RESTORED`)

Units RESTORED (natural-units-independent, able-to-fail, free-carrier-independent). The dimensions of
every real expression in scope:

```text
[I_ij]  = L⁻¹                (the transverse Darwin kernel)
[A_E]   = E·L = M L³ T⁻²      (the electric far-field coefficient)
[c_γ²]  = L² T⁻²
[V]     = L T⁻¹
[D_V]   = [A_V] = L² T⁻²
[R]     = L
[U_A]   = E = M L² T⁻²
[F_{A,2}] = E·L⁻¹ = M L T⁻²
[s]     = 1
```

The anchor `U_A = (s₁s₂A_E/4πR)[1 − (D_V+A_V)/2c_γ²]` closes: the leading `[A_E]/[R] = ML³T⁻²/L = ML²T⁻² = E`,
and the `O(v²/c_γ²)` correction `[(D_V+A_V)/c_γ²] = (L²T⁻²)/(L²T⁻²) = 1` is dimensionless bookkeeping, so
`[U_A] = E`. ✓ The force `[F_{A,2}] = [A_E]/(c_γ²R²)·[V²] = (ML³T⁻²)(L²T⁻²)/(L²T⁻²·L²) = MLT⁻² = EL⁻¹`. ✓
Each velocity order is dimensionless bookkeeping. (Corrupting one dimension — e.g. `[I_ij] = L⁻¹ → 1`, or
`[U_A] = E → EL⁻¹` via a bad radius — makes the dimensional-homogeneity object inhomogeneous and the units
tooth fires; free-carrier-independent.)

## 6. Cite-not-count discipline: what 036 consumes vs earns

Stage 036 **imports** and re-instantiates, but does NOT re-derive or re-count:

```text
CITED (provenance, NOT re-counted; de-dup deferred to Part VII):
  A_E     (Part IV / stage 032, register R63 — the selected electric ensemble's full far-field
           coefficient; carries R1_REQUIRED(bc_selection), sign + normalization unselected)
  c_γ     (Part III / stage003 — the light-cone speed, c_γ² = μ_R/ρ_br; re-cited via stage 034)

EARNED-NEW (the magnetism-NEW content of 036):
  { the Maxwell–Darwin transverse Darwin reference kernel  I_ij = (δ_ij + n_in_j)/8πR,
    the Lorentz-completed O(v²/c_γ²) anchor  U_A }
```

The reference kernel introduces **NO new irreducible knob** — it rides the already-registered electric
`A_E` R1 (R63) and the shared `c_γ`. `A_E`, `c_γ` are re-instantiated for the dual-engine check but are
cited, not re-counted (they are NOT re-added to the irreducible set; de-dup against Part IV / Part III is
deferred to Part VII). The `O(v⁴/c_γ⁴)` next order is an honest uncomputed remainder — NOT a new knob.

## 7. Scope — what 036 does NOT do (deferred to 037–039)

036 is the **Route-A reference-kernel gate only**. Named explicitly as OUT OF SCOPE:

- **The native moving-throat source law** `J_{T,i} = q_T s_i η_a V_i` + the parity census (035, V-2 —
  DONE). Route A does **NOT** use `J_T`; it boosts the electric interaction. 036 does NOT re-derive the
  source or the census.
- **The transverse-vector action row** `S_{T+move}` (034, V-1 — DONE, `109070da`). 036 does NOT re-open
  the action, the Hessian, or `G0_DAMAGE`.
- **⚠ `BOOST_COMMON_VELOCITY`** — the single common side-by-side boost cross-check
  `F_r/F_{E,r} = 1 − v²/(2c_γ²) + O(v⁴)` — is **stage 037 (V-4), NOT 036.** It reuses BOTH routes'
  equal-velocity limit and is allocated to the comparison stage; although it is adjacent to the 036 teeth
  in the source enum, the ratified allocation places it in 037 — it is scoped OUT explicitly so a reader
  never mistakes it for a Route-A-only tooth.
- **Route B / the direct moving-throat shear** (`U_B = −s₁s₂q_T²(D_V+A_V)/8πμ_R R`, computed BLIND to A),
  the structural **COMPARISON** (`ROUTE_INDEPENDENCE`, `COMPARE_COMPUTED`), and the **ratios**
  (`r_BA = q_T²/(ρ_br A_E)`, `δ_BA = r_BA − 1`, `r_cone = c_E²/c_γ²`, `ΔU`, `QMAG_R1`) are **stage 037
  (V-4)** — 036 emits NO direct route and NO comparison. ⚠ 036 produces ONLY the reference the direct
  route is tested against; the boost-consistency / emergent-Lorentz result (the two routes matching) is
  037's claim, NOT 036's.
- **The SEALED §4 first-match landing** — the 1152-cell truth table →
  `R1_REQUIRED(electric_bc_selection)` + the three co-blockers (`R1_REQUIRED(direct_moving_throat)` /
  `(magnitude)` / `(consistency)`) and the active-`F_flux` caveat — is **stage 038 (V-5)**. 036 emits NO
  landing; its verdict is the reference-kernel identity only.
- **The `b_T = ∇×u_T` time-reversal-EVEN DEPARTURE** (`B_TIME_REVERSAL_EVEN`, "not exact Maxwell") is
  **stage 039 (V-6)** — 036 does not touch `b_T`.
- **The imported `A_E`, `c_γ`** are **CITED provenance** from Part IV / stage003 — NOT re-derived, NOT
  re-counted in the irreducible parameter set (de-dup deferred to Part VII). The `A_E` sign + magnitude
  ride `R1_REQUIRED(bc_selection)` (stage 032, R63); 036 does NOT resolve it — the reference kernel is
  exact in STRUCTURE but tier-A conditional on the open `A_E` (this is exactly why the scope is
  EARNED-but-conditional, not GAPS-CLOSED).
- **The `O(v⁴/c_γ⁴)` next order** is named but NOT computed here (an honest R1-adjacent remainder owned
  downstream at 038's higher-order hook).

## 8. Source-to-stage predicate manifest

Completeness certificate: **no silently-dropped source claim.** Every source tooth in scope from the
source build (`magnetism_moving_throat_check.{py,wl}`, the 35-tooth `TOOTH_ORDER`, commit `53cf049f`)
lands as **PRESERVED** (folded as-is), **REPLACED_BY_STRONGER** (a stronger whole-stage reconstruction
owned here), or **SCOPED_OUT** (deferred to its own Part-V stage, each named LITERALLY — NO wildcards).
The partition is disjoint + exhaustive, computed at runtime in **both** engines from the same canonical
`(id, disposition, owner)` triples:

```text
partition = { PRESERVED: 2, REPLACED_BY_STRONGER: 4, SCOPED_OUT: 29 }   (35 total),
manifest digest (SHA-256) = f8b1569834c1c5cfd404fee4fba7bc49d3d7263f332e551270ad499fd3585cac.
```

- **PRESERVED (2):** `BOOST_NEXT_ORDER` (the honest-truncation order tooth, folded as-is) + the
  build-global `TARGET_BLINDNESS` (folded as the reference-side blindness guard — no downstream
  comparison/ratio/landing/sign token enters the reference computation).
- **REPLACED_BY_STRONGER (4):** `BOOST_PROJECTOR` and `BOOST_GENERAL_VELOCITIES` (→ the stronger inverse-FT
  reconstruction + Lorentz-anchor teeth owned here), `UNITS_RESTORED` (→ whole-stage dimensional firewall
  on the real expressions), and `DUAL_ENGINE_TERMS` (→ canonical term inventory across both engines).
- **SCOPED_OUT (29)** — the source teeth of V-1/V-2 (done upstream) and V-4/V-5, each named verbatim (same
  partition both engines):
  - **V-1 (034, owned UPSTREAM — DONE):** `FIELD_IDENTITY_UNITS`, `ACTION_KINETIC`, `ACTION_COUPLING`,
    `ACTION_STABILITY`, `G0_DAMAGE`, `LEDGER_READY_ROW` (6)
  - **V-2 (035, owned UPSTREAM — DONE):** `SOURCE_TRANSLATION_CONTINUITY`, `SOURCE_NOT_IMPORTED`,
    `SOURCE_BASIS`, `PARITY_RW`, `PARITY_PW`, `PARITY_ROTATION`, `PARITY_TIME_REVERSAL` (7)
  - **V-4 (037):** `DIRECT_SOURCE`, `DIRECT_PROJECTOR`, `DIRECT_EXCHANGE_SIGN`, `DIRECT_FALLOFF`,
    `DIRECT_VELOCITY_ORDER`, `ROUTE_INDEPENDENCE`, `BOOST_COMMON_VELOCITY`, `COMPARE_COMPUTED`,
    `DELTA_RATIO`, `CONE_RATIO`, `QMAG_R1` (11)
  - **V-5 (038):** `TRUTH_TOTALITY`, `TRUTH_PRECEDENCE`, `LANDING_OWNERSHIP`, `ACTIVE_FLUX_CAVEAT`,
    `HOOK_LORENTZ` (5)

  **V-6 (039) owns NO additional source-build tooth** — it REUSES 035's parity results and authors its own
  NEW `b_T` T-even departure assertions (`B_TIME_REVERSAL_EVEN`), cited-not-owned. ⚠
  `BOOST_COMMON_VELOCITY` is scoped-OUT to 037 even though it is adjacent to the 036 teeth in the source
  enum. Accounting: **3 in-scope Q-BOOST + 3 build-global (= the 6 owned teeth: 2 PRESERVED + 4
  REPLACED_BY_STRONGER) + 29 SCOPED_OUT = 35** ✓.

Both engines assert the identifier set (35 unique), the three-way disposition set, the exact counts, and
the committed digest at runtime and agree (`SOURCE_TO_STAGE_MANIFEST`; the mutation drops the final
partition row, which changes the identifier set / count / digest and fires the tooth; dropping any one of
the 35 rows fires it).

## 9. Consumes / cites

- **Consumes NOTHING numeric from 030/031/032/033/034/035.** 036 does not touch the puncture-deflection
  response matrix, the four `A_X` BC-ensemble data, the native source `J_T`, or the parity census — it is
  an independent boost-of-electric reconstruction.
- **Cites (provenance, NOT consumed numerically, NOT re-counted):** the electric interaction
  `U_E = s₁s₂A_E/4πR` and its radial force `F_{E,r} = s₁s₂A_E/4πR²` — the far-field coefficient `A_E`
  carries `R1_REQUIRED(bc_selection)` (sign + normalization unselected) — from **Part IV (charge)**
  (stage 032, register R63); and the light-cone speed `c_γ` (`c_γ² = μ_R/ρ_br`) from **Part III /
  stage003** (re-cited via stage 034). De-dup of `{A_E, c_γ}` against Part IV / Part III is deferred to
  Part VII; they are NOT re-added to the irreducible count here. Also cites the charge-sector `±w`
  orientation `s_i` (Part IV) as the sign carrier `s₁s₂`.

## Verification

- **Dual-engine, both exit 0, genuinely independent routes.**
  `scripts/ledger_stage036_route_a_maxwell_darwin_reference_sympy_audit.py` — **SymPy 8 teeth**.
  `mathematica/ledger_stage036_route_a_maxwell_darwin_reference_mathematica_audit.wl` — **Mathematica 8
  teeth** (`TOOTH_ORDER`: `BOOST_PROJECTOR`, `BOOST_GENERAL_VELOCITIES`, `BOOST_NEXT_ORDER`,
  `TARGET_BLINDNESS`, `DUAL_ENGINE_TERMS`, `UNITS_RESTORED`, `VERDICT_REDERIVATION`,
  `SOURCE_TO_STAGE_MANIFEST`). Standalone, print-only, assert-zero (`raise SystemExit(1)` / `Exit[1]`), no
  argparse harness, no JSON/YAML payload, **zero file-I/O between engines**. Each engine reaches
  `MAXWELL_DARWIN_REFERENCE` at `tier_A_conditional` on its own and prints its own tokens — cross-engine
  agreement is that they independently produce the same tokens, not a compare pass. **No dual-engine
  disagreement.**
- **The `.wl` is a genuinely INDEPENDENT route** — materially distinct construction, not a line-port of
  the `.py` **and not the source `.wl`'s route either** (the source `.wl` `deriveRouteA` used a
  seeded-component-Hessian off a hand-supplied `−R/8π` seed, the SAME method as the stage `.py`; the stage
  `.wl` uses NEITHER). The **SymPy** route is a **Poisson-seed Hessian reconstruction**: it derives the
  `k⁻⁴` radial seed coefficient from its Poisson equation `∇²(seed·r) = −1/4πR` (solving for the
  coefficient → `SEED_K4 = −R/8π`), then takes the Cartesian Hessian `−∂_i∂_j SEED_K4` component-by-
  component to reconstruct the transverse projector. The **Mathematica** route is **materially distinct**
  — it evaluates the analytically-continued **Riesz transform**
  `F⁻¹[k⁻²ᵅ] = Γ(3/2−α)R^{2α−3}/(4^α π^{3/2}Γ(α))` at `α=2` (giving `−R/8π` by analytic continuation, NOT
  a hand-supplied seed) and `α=1` (the Coulomb scalar), builds the projector via **genuine
  Schwinger/Gaussian proper-time momentum integrals**, and organizes the velocity anchor by an
  **independent parallel/perpendicular velocity decomposition** and geometric force route. Each engine
  derives its OWN seed, its OWN kernel, and its OWN anchor by a different construction — materially
  distinct from BOTH the stage `.py` AND the source `check.wl` (NO seeded Hessian in the `.wl`).
- **Manifest digest** (SHA-256) `f8b1569834c1c5cfd404fee4fba7bc49d3d7263f332e551270ad499fd3585cac`;
  disposition partition `{PRESERVED: 2, REPLACED_BY_STRONGER: 4, SCOPED_OUT: 29}` (35 total), i.e. the 3
  in-scope Q-BOOST teeth + 3 build-global guards (6 owned) + 29 scoped-out = 35, computed at runtime in
  both engines and agreed.
- **Per-tooth ablation** (env switch `LEDGER_STAGE036_MUTATION`): all **8 teeth's mutations
  FIRED_AT_OWN_ASSERT per engine** (16 mutation runs across the two engines) — the three Q-BOOST teeth
  (`BOOST_PROJECTOR` doubles the reconstructed kernel `I → 2I`; `BOOST_GENERAL_VELOCITIES` drops the
  `A_V` term; `BOOST_NEXT_ORDER` sets the claimed order to 4), the build-global `TARGET_BLINDNESS` /
  `DUAL_ENGINE_TERMS` / `UNITS_RESTORED`, the `VERDICT_REDERIVATION` tooth, and the
  `SOURCE_TO_STAGE_MANIFEST` tooth.
- **The verdict tooth is non-tautological** — it mutates a **COMPUTED** object, not the final token. Under
  the `VERDICT_REDERIVATION` mutation it **DOUBLES the COMPUTED kernel** (`verdict_kernel = 2·KERNEL_AT_R`),
  recomputes the projector residuals, and asserts the verdict RE-DERIVES — from the computed projector
  zero-difference, the anchor zero-difference `route_a_u2 − EXPECTED_A2`, and the velocity-order object
  (never a duplicated literal) — to the authored stage-local token `MAXWELL_DARWIN_REFERENCE_UNCERTIFIED` (a
  stage-local "could-not-certify" alternate mirroring stage 035's `PARITY_INCONSISTENT` precedent, NOT a
  literal flip and NOT a physical non-Maxwell departure claim). The precedence it re-derives:

  | Computed state | Re-derived token |
  |---|---|
  | transverse Darwin kernel `(δ_ij+n_in_j)/8πR` matches (BOOST_PROJECTOR zero) **and** anchor matches `EXPECTED_A2` (BOOST_GENERAL_VELOCITIES zero) **and** velocity order = 2 (BOOST_NEXT_ORDER) | **`MAXWELL_DARWIN_REFERENCE`** (production; EARNED reference kernel, `tier_A_conditional` on the electric `A_E` R1) |
  | reconstructed kernel is not the transverse Darwin form (wrong `1/8π`, broken transversality), **or** the anchor mismatches `EXPECTED_A2` (dropped `A_V`), **or** the truncation order is misrepresented (`velocity_order = 4`) | **`MAXWELL_DARWIN_REFERENCE_UNCERTIFIED`** (the "could-not-certify" partition) |

  Three **negative controls** test the else-branch inside the verdict tooth: `projector_negative` (the
  doubled kernel), `anchor_negative` (the `A_V`-dropped anchor, only `D_V` retained), and `order_negative`
  (velocity order False) each re-derive to `MAXWELL_DARWIN_REFERENCE_UNCERTIFIED`. ⚠ The **open electric
  `A_E` R1 does NOT change the token** — production is `MAXWELL_DARWIN_REFERENCE` at `tier_A_conditional`
  REGARDLESS of the open `A_E`; an unresolved coefficient is a scope tag, not a certification failure (so
  it must NOT re-derive to `MAXWELL_DARWIN_REFERENCE_UNCERTIFIED`). Flipping the
  `MAXWELL_DARWIN_REFERENCE` literal tests nothing (the 030 `X ≡ X` lesson). `MAXWELL_DARWIN_REFERENCE`
  and `MAXWELL_DARWIN_REFERENCE_UNCERTIFIED` are the bookend-adjudicated ledger-stage verdict tokens
  (ratified in `part5_magnetism_atomic_split.md`; the source build has NO Route-A verdict enum, so none
  was overlooked); same token/partition in both engines.
- **Tri-review outcome (falsification-first — recorded transparently, not hidden).**
  - **FIDELITY and ADVERSARIAL CONVERGED on 2 NITs — both can't-fail conjuncts** (neither could ever fail,
    so neither weakened a real check): (1) a `NEXT_UNCOMPUTED_ORDER == 4` sub-clause — a literal `4 == 4`
    — ANDed into `BOOST_NEXT_ORDER`; (2) the `dimension_firewall` inventory `len == len` guard (`13 == 13`
    / `29 == 29`) riding `UNITS_RESTORED`. **Both REMEDIATED** — removed; `BOOST_NEXT_ORDER` still fires
    via the `(0,2)` / `claimed_order` / `ε⁴` order objects, and `UNITS_RESTORED` still fires via the
    dimension residuals; the **manifest digest is unchanged**
    (`f8b1569834c1…3585cac`) → **CLEAN**.
  - **One DOCUMENTED non-blocking note** (verified-safe, NOT remediated): `BOOST_GENERAL_VELOCITIES`
    bundles the force-derivation / component / radial residuals under one mutation (a compound tooth,
    within the directive's permissive per-tooth-ablation "e.g." scope) — the same class as stage 035's
    compound-tooth note.

  Arbiter re-runs of both engines reproduce the tokens and the manifest digest
  `f8b1569834c1…3585cac`; the tri-review leaves the reference kernel and Lorentz anchor EARNED and the
  verdict `MAXWELL_DARWIN_REFERENCE` at `tier_A_conditional`.

## Downstream consumers

- **Part V continues** after 036: 037 (Route B — the direct moving-throat shear, computed BLIND to A — +
  the structural boost-consistency COMPARISON + the ratios `r_BA`/`δ_BA`/`r_cone`/`ΔU`, which tests the
  direct route against 036's reference kernel) → 038 (the sealed `R1_REQUIRED(electric_bc_selection)`
  landing) → 039 (the `b_T` T-even departure). 036 is the Maxwell-consistent reference those stages
  consume by citation; the boost-consistency result is 037's, NOT 036's.
- **Parameter register:** edge **R70** (the Maxwell–Darwin transverse Darwin kernel
  `I_ij = (δ_ij + n_in_j)/8πR`, `[I_ij] = L⁻¹`, + the Lorentz-completed `O(v²/c_γ²)` anchor `U_A`;
  **DERIVED** (boost-of-electric — inverse-FT of the Coulomb-gauge transverse projector `P^T_ij/k²`);
  **explicitly tier-A conditional** — the kernel's STRUCTURE is exact but its coefficient rides the
  electric `A_E` (R63, `R1_REQUIRED(bc_selection)`) and `c_γ` (stage003) ⇒ the verdict
  `MAXWELL_DARWIN_REFERENCE` at `tier_A_conditional`; introduces **NO new irreducible knob** — do NOT
  re-add `A_E` (R63) or `c_γ` (stage003), cite them; the `O(v⁴/c_γ⁴)` next order is an honest uncomputed
  remainder, NOT a new knob).
- **`docs/model_map.md` §3.5:** the Q-BOOST / Route-A reference bullet — boosting the electric interaction
  reconstructs the Maxwell–Darwin transverse Darwin kernel `(δ_ij + n_in_j)/8πR` and the Lorentz-completed
  `O(v²/c_γ²)` anchor; verdict `MAXWELL_DARWIN_REFERENCE` at `tier_A_conditional`; the reference the direct
  route (037) is tested against.
- **Part VII:** the cited `{A_E, c_γ}` edges enter the de-dup (counted ONCE, with Part IV / Part III); the
  `A_E` obligation (`R1_REQUIRED(bc_selection)`) enters the honest R1 debt that the reference kernel
  inherits.

## Provenance

- **Physics source:** `software/em_charge_attribute/magnetism_moving_throat_result.md` (the "Q-BOOST
  (Route A)" section — the transverse Darwin kernel `I_ij(R) = (δ_ij+n_in_j)/8πR`, the Lorentz-completed
  conditional anchor `U_A` / force `F_{A,2}` / radial `F_{A,2,r}`, the explicit `O(1)`+`O(v²/c_γ²)` orders
  with the next uncomputed `O(v⁴/c_γ⁴)` — plus Appendix A "Independent derivations and dimensions" — the
  inverse-FT reconstruction `F⁻¹[P^T_ij/k²]` and the restored-units table) +
  `software/em_charge_attribute/magnetism_moving_throat_check.{py,wl}` (the Route-A block `derive_route_a` /
  `deriveRouteA`, commit `53cf049f`). Stage 036 extracts the **Route-A reference-kernel cluster** into a
  focused standalone dual-engine pair; the source `argparse`/`--compare` harness, the payload plumbing, the
  Route-B / comparison / ratios / truth-table machinery (out of 036's scope), and all JSON/log file
  writing are STRIPPED (print-only / zero-file-I/O / independent-tokens contract). The stage `.wl` is a
  materially distinct Wolfram route (Riesz analytic continuation + Schwinger/Gaussian momentum integrals +
  parallel/perpendicular velocity decomposition), not a line-port of the `.py` and not the source `.wl`'s
  seeded-Hessian route.
- **Consumes:** nothing numeric from 030/031/032/033/034/035 (independent boost-of-electric reconstruction).
- **Cites (provenance):** the electric interaction `U_E = s₁s₂A_E/4πR` and `A_E` (with
  `R1_REQUIRED(bc_selection)`) from Part IV / stage 032 (register R63); `c_γ` (`c_γ² = μ_R/ρ_br`) from
  Part III / stage003 (re-cited via stage 034) — NOT re-derived, NOT re-counted; the Part IV `±w`
  orientation `s_i` as the sign carrier.
- **Governing:** `notes/ledger_v2_blueprint.md` §5 (reshape spec) + §6 (per-tooth ablation);
  `notes/part5_magnetism_atomic_split.md` (V-3 = Route A the Maxwell–Darwin reference; the tooth-allocation
  table); `docs/model_map.md` §3.5 (Q-BOOST / Route-A reference bullet). Reshape directive + review trail:
  `research/pde_ledger_v2/_scratch/stage036_reshape_directive.md`.
