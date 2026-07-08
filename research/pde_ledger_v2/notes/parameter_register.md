# Parameter Register — knobs, dimensions, and reductions (living; feeds Part VII)

> **Purpose.** Track every parameter the ledger introduces, its dimension, where it enters, and — the load-bearing
> part — **how the parameters relate**, so the *irreducible* free-parameter count (not the nominal symbol count) can be
> read off at Part VII. Predictive surplus = held-out predictions **minus irreducible knobs**; this register is the
> honest denominator. Started at stage 005 (2026-07-07) and grown as a lightweight per-stage finalize step, because every
> stage already emits its knobs, dimensions, and carried residuals — this just aggregates them *relationally*.
>
> **Why relational, not a flat list (user, 2026-07-07):** some knobs are *manifestations* of other knobs — free only
> because a specific reduction is deferred. The final count is the **dimension of the parameter variety after the derived
> relations are imposed**, not the number of symbols. We already have two worked examples (pathA_40, pathA_41) and the
> tool to measure it (codimension / Krull-dimension count — pathA_40's `Δr=2`).

---

## Method & discipline (falsification-first, applied to the register itself)

A relation only *reduces* the count when it is **DERIVED** or **CODIM-PROVEN** — never when it is **IMPOSED**. An
imposed relation (a calibration or postulate that makes two knobs look like one) is a hidden tuning, not a reduction.

- A **named-but-unexecuted** reduction route is **reduction debt**: the knob stays counted as free, flagged with the
  exact derivation that would discharge it. Debt is honest, not a reduction.
- If a named route comes back **CLOSED-NEG** (proven not to exist), that *hardens* the knob as irreducible — real
  information, a gain not a loss.
- **Be symmetric.** Reductions can also reveal **hidden multiplicity** — one apparent knob that is secretly two
  (pathA_40's `Δr=2` did exactly this). A register that only ever finds good news is broken.

### Provenance classes
| Class | Meaning | Counts toward irreducible? |
|---|---|---|
| `ACTION` | primitive action/state constant — a genuine free input | **yes** |
| `DERIVED` | a function of other parameters (proven) | no (collapses into its inputs) |
| `CONV` | convention / units / pin choice — never free | no |
| `FREE-UNREDUCED` | free now, but a **named** reduction route is pending (reduction debt) | **yes, flagged** |
| `CALIB` / `CALIB-ANCHOR` | calibrated to a benchmark / anchor (magnitude tuned) | **yes** (a tuned knob) |
| `CANDIDATE` | labeled candidate, not yet load-bearing | tracked, not counted |
| `GAP` | a *missing* derivation (not a knob) — a deferred obligation | tracked, not counted |

### Relation-status flags (edges)
`DERIVED` (proven) · `CODIM-PROVEN` (independence/dependence established by a dimension count) · `IMPOSED` (a
calibration/postulate — a tuning, does **not** reduce) · `PENDING` (named route, unexecuted = debt) · `CLOSED-NEG`
(route proven absent → hardens the knob).

### Dimensions
Dimensions are exact `{L,T,M}` triples **only where a built stage has dual-engine-verified them**. Params from
not-yet-reshaped sectors carry `pending` and are filled in when that stage is built — the register never asserts a
dimension it has not checked. This is the register doubling as the running cross-stage dimensional ledger (the
whole-system dim check at Part VII then confirms rather than discovers).

---

## Master parameter table

Built stages: 001 (pathA_21c primitives), 002 (pathA_21c force), 003 (pathA_36 light), 004 (pathA_19 dim foundation),
005 (pathA_20 sound speed), 006 (χ_B two-phase ontology — fresh-authored), 007 (pathA_35 G0 shear-surface freeze —
the formal home of `T0_SHEAR_FROZEN(d9520d3819c3)` + `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`; the 11 computed from the
enumerated table in both engines), 008 (pathA_28 monopole/dipole return constraint-spec — **adds NO new knobs**:
the source moments `M0/D1/Q2` are constructions not parameters (`Q2` = FREE ANCHOR export), the derivative-vertex
coefficient `η_ℓ` is a recorded `branch_assumption` (tracked, not counted, no dim pinned), and `cancellation_possible`
is a scope FLAG not a quantity; the stage's obligation is edge R23). Seeded knit params (Part VI, not yet reshaped):
pathA_40 cone-lock, pathA_41 NG5.

| Param | `[L,T,M]` dim | Enters | Class | Depends on / relation | Reduction route + status |
|---|---|---|---|---|---|
| `ħ` | `M L² T⁻¹` | I-1 (004) action | `ACTION` | — | provenance is a `GAP` (`HBAR_FREE_SUBSTRATE_RELATION_MISSING`) |
| `m_GNLS` | `M` | I-1 (004) action | `ACTION` | — | `m_defect` emergence = `GAP` (`INFLOW_MASS_SOURCE_MISSING`) |
| `K` (EOS) | `M L¹⁸ T⁻²` | I-1 (004), I-2 (005) | `ACTION` | — | EOS **exponent 5** is `IMPOSED` (`EOS_CLOSURE_IMPOSED`) |
| `ρ0` | `L⁻⁴` | I-1 (004) state datum | `ACTION` (state) | — | chosen asymptotic 4D-bulk number density |
| `c_s0` | `L T⁻¹` | I-1/I-2 | `DERIVED` | `= √(5K ρ0⁴/m_GNLS)` | DERIVED from `{K, ρ0, m_GNLS}` |
| `ξ_h` | `L` | I-1 (004) | `DERIVED` | `= √2 ħ/(m_GNLS c_s0)` | DERIVED (core balance) |
| `h0` | `M L² T⁻²` | I-1 (004) | `DERIVED` | `= (m_GNLS c_s0²)/4` | DERIVED (`EOS_FROM_GNLS_FACTOR`) |
| `a` (pin) | `L` | I-1 (004) | `CONV` | `= ħ/(m_GNLS c_s0)` | pin choice — never free (`A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`) |
| `λγ = c_γ/c_s` | `1` | I-2 (005) | `DERIVED` | `= c_γ/c_s0` | **not independent** — a ratio of `c_γ` and `c_s0` (see edge R3) |
| `c_γ` (gauge cone) | `L T⁻¹` | I-2 (005), III (003), VI | `FREE-UNREDUCED` | brane: `c_γ²=μ_R/ρ_br`; bulk: `c_γ²=C_B/C_E` | Route A `PENDING`; cone-lock `IMPOSED` `λγ=1` |
| `c_E` (throat Green speed) | `L T⁻¹` | VI (040/041) | `FREE-UNREDUCED` | throat dynamic Green speed — **distinct from the Maxwell `C_E`** | cone-lock `IMPOSED` `c_E=c_γ` (R8); Route-A registered-deferred |
| `C_E, C_B` (gauge metric) | `[C_E]=M⁻¹L⁻⁴T²`, `[C_B]=M⁻¹L⁻²` (ratio `L²T⁻²`) | I-2 (005) | `FREE-UNREDUCED` (bulk) | `c_bulk²=C_B/C_E` | brane-zero-mode reduction to `c_γ` `PENDING` |
| `μ_R` (brane shear) | `M L⁻¹ T⁻²` (stages 003, 007) | **I-4 (007) freeze home**, III (003), VI | `FREE-UNREDUCED` (brane) | `c_γ²=μ_R/ρ_br`; one of the G0 "11" (erratum: count STANDS, no overcount) | Route A `PENDING` (R10, `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT`); ≠ `μ_R⁽⁴⁾` (R17 `PENDING`) |
| `ρ_br` (brane density) | `M L⁻³` (stages 003, 007) | **I-4 (007) freeze home**, III (003), VI | `FREE-UNREDUCED` (brane) | `c_γ²=μ_R/ρ_br`; one of the G0 "11" (erratum: NOT pathA_25's `varrho_br[ρ]` — that object is `OUT_OF_ACTIVE_NG5`) | Route A `PENDING` (R10) |
| `λ_Pu` (parity-repaired P–u coupling) | `M L⁻¹ T⁻²` (stage 007) | I-4 (007) freeze | `ACTION` | `L_Pu = −λ_Pu ϖ_a Ω_uᵃ`, `ϖ_a=(ŵ×P_∥)_a` — the parity-EVEN operator; re-admits the ε-contracted/chiral class T0 excluded + REQUIRES the conceded `ŵ` (structural postulate 5) | one of the G0 "11"; no reduction route named |
| `Ω_w` (bare u_w gap scale) | `T⁻¹` (stage 007) | I-4 (007) freeze | `ACTION` | `L_uw` gap term `−½ρ_br Ω_w² u_w²` | one of the G0 "11"; no reduction route named |
| `g_ℓ(w)` + width `ℓ_g` | `[g_ℓ]=L⁻¹`, `[ℓ_g]=L` (stage 007) | I-4 (007) freeze profile | `ACTION` (function: fixed Gaussian shape, ONE width knob) | `g_ℓ=exp(−(w/ℓ_g)²)/(√π ℓ_g)`, `∫g_ℓ dw=1` derived; admitted on locality/minimality grounds ONLY (target-blind G0.2) | one of the G0 "11". **Superseded as the *material-state* wall by stage006's χ_B (R21); REMAINS the constitutive-freeze profile** |
| force-magnitude norm | `1` (dimensionless coeff) | II (002) | `CALIB` | matched to inter-defect force strength | form earned; magnitude CALIBRATED; `G=GENUINE_BLOCKED` |
| `Q_E` (charge mag.) | pending (pathA_38) | IV (038) | `CALIB-ANCHOR` | — | pending pathA_38 reshape |
| `ρ_B0` | `M L⁻³` (stage 003) | III (003), VI (041) | `FREE-UNREDUCED` | on-brane compression; pathA_41 marks active irreducible | NG5 route (i) 4D→3D compression `PENDING` |
| `χ_c` | `M L⁻⁵ T²` (stage 003) | III (003), VI (041) | `FREE-UNREDUCED` | on-brane compression; pathA_41 marks active irreducible | NG5 route (i) `PENDING` |
| `B` (brane compression modulus) | `M L⁻¹ T⁻²` (stage 003) | III (003) | `FREE-UNREDUCED` (brane) | longitudinal `(∇·u)²` stiffness | (part of the light-sector departure) |
| `J` (throat winding/current) | `L² T⁻¹` (stage 003) | III (003) | `ACTION` | enters `C_J=−Jρ_B0` | — |
| `K_θ` / `κ_phase` | `M L T⁻²` (stage 003) | III (003) | `ACTION` (labeled postulate) | phase-gradient stiffness (sign is the θ-as-φ crux) | carried postulate |
| `C_J = −J ρ_B0` | `M L⁻¹ T⁻¹` (stage 003) | III (003) | `DERIVED` | `= −J ρ_B0` (sign-sensitive IBP, R15) | DERIVED |
| `B_eff = ρ_B0²/χ_c` | `M L⁻¹ T⁻²` (stage 003) | III (003) | `DERIVED` | `= ρ_B0²/χ_c` (R16) | DERIVED |
| `C_hu` | pending (pathA_41) | VI (041) | `FREE-UNREDUCED` | embedding-mixing | NG5 route (ii) embedding-overlap `PENDING` |
| `α_J` (mass-bridge coeff) | `1` | I-2 (005) | `CANDIDATE` | `m_defect = α_J ħJ/c_γ²` | labeled candidate; not load-bearing |
| `m_defect` | `M` | — | `GAP` | candidate `ħJ/c_γ² = M` (dimensional-only) | pathA_21 (deferred) |
| flux law `J_crit` | — | I-2 (005) | `GAP` | — | `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA` (deferred) |
| `χ_B` (order field) | `1` (∈[0,1]) | I-3 (006) action | `ACTION` (field) | independent scalar OP, **NOT** `\|P_∥\|²` (pin P7 — the rung-W T1-escape) | route (c) composite `χ_B=\|P_∥\|²` = named future gate, high-risk (neighbors falsified: v3 arrows wall, pathA_25, pathA_35); wall is route (a) postulated |
| `a_B` (double-well) | `M L⁻² T⁻²` | I-3 (006) action | `ACTION` | `f_B=a_B χ_B²(1−χ_B)²`, minima {0,1}, n-independent (P3/P4) | POSTULATED — parent `U(ρ)` is single-well (the wall cannot come from it) |
| `κ_B` (interface gradient) | `M T⁻²` | I-3 (006) action | `ACTION` | kink width `δ=√(κ_B/2a_B)`, `σ_wall=√(2a_B κ_B)/6` DERIVED from `{a_B,κ_B}` | — |
| `α_aniso` (P-orientation) | `M L⁻² T⁻²` | I-3 (006), dims-only | `ACTION` | `POSTULATED_ANISOTROPY` (rung_W); P dynamics = I-4 | ingredient chosen for the outcome — no independent mechanical provenance |
| `Γ_B` (conversion law) | `T⁻¹` | I-3 (006) balance | `ACTION` (law/field) | `Γ_B=Γ_return−Γ_drain`; global-return-constrained `R_0=−M_0, R_1=−D_1` (postulates, not locally asserted) | return/drain closure = deferred (pathA_28/29 lineage) |
| χ_B-gating structure | — (structural) | I-3 (006) action | `ACTION` (structural choice) | `χ_B f_shear` multiplicative gate (shear only where ordered) | structural 6th member of DRIFT(6) |
| `μ_R⁽⁴⁾` (4D shear-stiffness density) | `M L⁻² T⁻²` | I-3 (006) shear gate | `FREE-UNREDUCED` | brane `μ_R = ∫χ_B μ_R⁽⁴⁾ dw` (dim-consistent only) | projection edge R17 `PENDING` |
| `W_slab` (brane slab width) | pending (`L` expected; no audit verifies it — the stage checks only the projection `[W]=L⁻¹`) | I-3 (006) carried limitation | `FREE-UNREDUCED` | bulk-on-both-sides ⇒ kink–antikink slab; double-well selects NO width (kink admission ≠ slab stability) | old `L/a` self-selection item — "requires dynamics", sim-deferred (R19 `PENDING`) |
| `M_χ` (order mobility), `J_χ` (order transport) | `L²T M⁻¹` / `L⁻³T⁻¹` (both dual-engine-verified: P8 adjunct row; balance-row divergence) | I-3 (006) adjunct rows | `GAP`/deferred | dynamics adjunct (P8), `J_χ=0` default (P9) | tracked, not counted |
| `κ_4` (Lifshitz stabilizer) | pending (`M L³T⁻²` expected; used algebraically only in the leg-C probe, no dim check scripted) | I-3 (006) probe symbol | `GAP`/deferred | leg-C Lifshitz probe only (`k*²=−κ_phase/2κ_4`) | tracked, not counted |
| `M_n` (repair mobility) | `L⁻⁴ T M⁻¹` (dual-engine-verified, P12) | I-3 (006) throat-ontology driver | `GAP`/deferred | `J_repair ~ −M_n ∇μ` (admittance/outlet driver, not a force law) | tracked, not counted |
| `f_throat`, `f_mix` (placeholders) | dim REQUIREMENT `M L⁻² T⁻²` verified on the placeholder symbols; content undeclared | I-3 (006) action placeholders | `GAP` (`DEFERRED_PLACEHOLDER`) | mouth/wall + mixed-gauge couplings, content deferred | not knobs until given content |

*Note — stage 001 introduces NO free knobs:* `Ω_2=4π`, `Ω_3=2π²`, `⟨n_i n_j⟩=δ_ij/d` are DERIVED geometric constants.

### The G0 structural-postulate block (I-4, stage 007) — the 6 structural members of the "11"

Structural `ACTION` choices (no dimensions; counted in the freeze's `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` alongside
the 4 constants `{ρ_br, μ_R, λ_Pu, Ω_w}` + the 1 function `g_ℓ(w; ℓ_g)`; enumeration computed in both engines):
1. imposed `ŵ` axis + `w=0` surface (conceded-wall);
2. `uᵃ` = same-medium surface collective, tangentially free-slip (`u̇ᵃ ≠ vᵃ`);
3. T0 `Pⁱ` reused as the Cosserat micro-rotation reservoir (0 new DOF);
4. baseline `Pⁱ` spin-wave status = `massless` (alternates `gapped`/`slaved-rigid` named-inactive);
5. the `ŵ`-dependent parity-EVEN P–u operator (re-admits the ε-contracted/chiral class excluded by T0; requires the
   conceded `ŵ` — a structural-postulate cost, not a free choice);
6. no C5 `φ` analog / no longitudinal constraint (the flat-brane DOF=8 carries `φ`=0 — the fact the C5 crux attacks).

New-field content at G0 (`uᵃ` 3 + `u_w` 1 = 4 DOF) is SEPARATE from the 11-input drift count. The T0 couple-stress
coefficients are KEPT (0 new). Anti-absorption (2026-07-04 erratum): `{ρ_B0, χ_c, C_hu}` are the Part-VI (pathA_41)
cross-sector drift — NEVER absorbed into the 11 (guard asserted in both engines).

---

## Relations & reductions ledger (the edges)

| # | Relation | Type | Source | Effect |
|---|---|---|---|---|
| R1 | `c_s0 = √(5K ρ0⁴/m_GNLS)` | `DERIVED` | I-1/I-2 (004/005) | collapses `c_s0` into `{K, ρ0, m_GNLS}` |
| R2 | `ξ_h = √2 ħ/(m c_s0)`, `a = ħ/(m c_s0)`, `h0 = (m c_s0²)/4` | `DERIVED`/`CONV` | I-1 (004) | collapse `ξ_h, a, h0` into primitives |
| R3 | `λγ = c_γ/c_s0` | `DERIVED` | I-2 (005) | `λγ` is **not** an independent knob — the free thing is `c_γ` |
| R4 | `c_γ² = μ_R/ρ_br` (brane cone) | `DERIVED` (within pathA_36) | III (003) | brane light-cone set by brane moduli |
| R5 | `c_bulk² = C_B/C_E` (bulk gauge cone) | `DERIVED` (within pathA_20b) | I-2 (005) | bulk gauge cone set by the Maxwell metric ratio |
| R6 | brane `c_γ` ← bulk `c_bulk` (zero-mode reduction) | `PENDING` | I-2/III | debt: would relate the brane & bulk gauge speeds (`BRANE_ZERO_MODE_REDUCTION_UNDERIVED`) |
| R7 | `λγ = 1` i.e. `μ_R/ρ_br = 5K ρ0⁴/m` | `IMPOSED` (cone-lock calibration) | VI (040) | a **calibration**, not derived — does not reduce; independence proven by R9 |
| R8 | `c_E = c_γ` (the **throat Green speed** `c_E`, NOT the Maxwell coeff `C_E`) | `IMPOSED` (cone-lock calibration) | VI (040) | a second independent calibration |
| R9 | R7 and R8 are **independent** calibrations, `Δr=2` | `CODIM-PROVEN` | VI (040) | Krull/`RegionDimension` count: real-locus dim 10→8. **Hidden multiplicity** — two locks, not one |
| R10 | Route A: derive `λγ` from a nonlinear-throat `μ_R`-as-bulk-defect integral | `PENDING` | VI (040) | would discharge R7's debt (`λγ` derived, not calibrated); needs the deferred nonlinear throat |
| R11 | Route B: derive `λγ` via `h≠u_T` / thin-plate over-import | `CLOSED-NEG` | VI (040) | proven **not** a route → `λγ` reducible only via Route A (or stays calibrated) |
| R12 | NG5 route (i): compression-sector 4D→3D reduction | `PENDING` | VI (041) | would collapse `{ρ_B0, χ_c}` into bulk/brane primitives |
| R13 | NG5 route (ii): embedding-overlap reduction | `PENDING` | VI (041) | would collapse `C_hu` |
| R14 | Location-closure: every param ∈ {4D bulk, 3D brane surface, throat seam}; **no fourth arena** | `CODIM-PROVEN` (computed, able-to-fail) | VI (041) | the NG5 "drift" is **un-reduced brane-surface params**, NOT a second substance |
| R15 | `C_J = −J ρ_B0` | `DERIVED` (sign-sensitive IBP) | III (003) | collapses `C_J` into `{J, ρ_B0}` (stage 003's fidelity upgrade over pathA_36) |
| R16 | `B_eff = ρ_B0²/χ_c` | `DERIVED` | III (003) | collapses `B_eff` into `{ρ_B0, χ_c}` (the stray-longitudinal stiffness) |
| R17 | brane `μ_R = ∫χ_B μ_R⁽⁴⁾ dw` (shear-gate projection) | `PENDING` | I-3 (006) | debt: would relate the brane shear modulus to the 4D gate density (dim-consistency asserted only) |
| R18 | θ-as-Maxwell-φ `WITH_PROVENANCE` route | `CLOSED-NEG` | I-3 (006) | the Maxwell locus `K_θ=+J²ρ_B0²/ρ_br, B_eff=0, m_θ²=0` is reachable `BY_TUNING` only; the sole provenance sign-flip (`κ_phase<0`) is a Lifshitz instability = the falsified pathA_25 wall. Hardens the no-go; the θ-branch symbols are NOT admitted as live knobs |
| R19 | slab-width selection (`W_slab` from dynamics) | `PENDING` | I-3 (006) | debt: kink–antikink slab width not selected by the double-well; = the old `L/a` self-selection item (sim-deferred Gate-6 territory) |
| R20 | `δ, σ_wall` from `{a_B, κ_B}` | `DERIVED` | I-3 (006) | kink width + surface tension collapse into the well/gradient constants (single-kink admission only — see R19) |
| R21 | stage006 χ_B wall supersedes fixed-shape `g_ℓ(w)` as the *material-state* closure; the G0 freeze REMAINS the light-sector *constitutive* freeze | structural (scope split — NOT a reduction; both facts asserted in stage 007) | I-3 (006) / I-4 (007) | `ℓ_g` stays a counted knob of the constitutive freeze; stage003 consumes the frozen `L_Mac` as-is; neither artifact retro-invalidates the other |
| R22 | `[μ_R] ≠ [μ_R⁽⁴⁾]` (M L⁻¹T⁻² vs M L⁻²T⁻²) + R17 dim-consistency `[μ_R⁽⁴⁾]·L = [μ_R]` | structural (dimensional firewall, dual-engine asserted + able-to-fail — NOT a reduction; R17 remains `PENDING`) | I-4 (007) | the notational firewall: two DISTINCT symbols related only by the `PENDING` R17 projection — prevents a silent conflation that would fake a reduction |
| R23 | return-cancellation targets `R0(ω)=−M0(ω)`, `R1_i(ω)=−D1_i(ω)` (the moments any brane↔bulk return must cancel) | `PENDING` (constraint-spec obligation — bookkeeping `x−x` identity, NOT a derived cancellation; whether an admissible return delivers it = stages 009/010 + the Gate-6 `Z0_ret/Z1_ret` selector) | II-B1 (008) | the pde-ledger open-item-#9 target made precise; `cancellation_possible` stays a literal scope flag until track-3 decides |

---

## Provisional rollup (honest — a nominal-with-classes tally, NOT the final count)

The final irreducible number is a **Part VII codimension count** of the assembled constraint variety (the `Δr=2`
technique scaled up). Do **not** assert an irreducible number before then. Current provisional picture:

- **Genuinely-free medium primitives:** `{ħ, m_GNLS, K, ρ0}` (4), plus the EOS **exponent 5** as an `IMPOSED` structural
  choice.
- **Already collapsed (proven not free):** `c_s0, ξ_h, h0, a, λγ` + the stage-003 derived `C_J=−Jρ_B0` (R15),
  `B_eff=ρ_B0²/χ_c` (R16) — `λγ` in particular is *not* an independent knob, it is `c_γ/c_s0` (R3). This is the pattern
  the user predicted: apparent knobs that are manifestations of others.
- **Free-unreduced with NAMED routes (reduction debt):** the gauge normalization `c_γ` / `(μ_R, ρ_br)` / `(C_B, C_E)`
  (throat Green speed `c_E` is a *distinct* speed, R8) via Route A (R10) + brane-zero-mode (R6); the brane-longitudinal
  constants `{B, ρ_B0, χ_c}` (light-sector departure, dims verified stage 003); the NG5 trio `{ρ_B0, χ_c, C_hu}` via
  routes (i)/(ii) (R12/R13). **Every one has a named pending derivation** — none is asserted irreducible.
- **Calibrated (tuned):** the force-magnitude normalization (II); `Q_E` (IV). `G` is `GENUINE_BLOCKED`.
- **CLOSED-NEG so far:** R11 (one λγ route closed) — hardens the picture but does not yet make λγ fully irreducible
  (Route A remains).
- **Hidden multiplicity found:** R9 (`Δr=2`) — the two cone locks are independent, a sobering result the register is
  obligated to record.
- **GAPs (deferred obligations, not knobs):** `m_defect`, `ħ` provenance, flux `J_crit`, the χ_B dynamics adjunct
  (`M_χ, J_χ`).
- **The χ_B package (I-3, stage 006) — DRIFT(6), counted honestly:** `{χ_B; a_B; κ_B; α_aniso; Γ_B; gating structure}`
  are new `ACTION` inputs — **the wall is a postulated field, not a derived structure** (route (a); route (c)
  `χ_B=|P_∥|²` is a named, high-risk future gate, its neighbors already falsified). `δ`/`σ_wall` collapse into
  `{a_B, κ_B}` (R20). Two NEW `FREE-UNREDUCED` entries with named routes: `μ_R⁽⁴⁾` (projection R17) and the slab width
  `W_slab` (R19 — kink admission ≠ slab stability). `ρ_B0, χ_c` appear in BOTH the stage-006 dead θ-branch and the
  pathA_41 Part-VI drift trio — counted ONCE (Part VI); the θ-branch is dead (`THETA_BRANCH_DEAD_NOT_ADMITTED`) and
  R18 hardens it `CLOSED-NEG`.

- **Part II begun (stage 008, pathA_28): zero new knobs.** The constraint-spec introduces obligations, not
  parameters — R23 (the return-cancellation targets) is `PENDING` debt for stages 009/010 + Gate-6, `Q2` is exported
  as a FREE ANCHOR (never derived ℓ=2), `η_ℓ` is a recorded branch assumption (tracked, not counted). The frozen
  slice (`G=c=c_s=1`, `K_eos=1/500`, `(a*,L*)`) is benchmark calibration, cited not registered.
- **The G0 freeze package (I-4, stage 007) — the "11", counted honestly (computed in-engine):** 4 constants
  {`ρ_br`, `μ_R` (rows above, re-homed to I-4; Route-A R10 `PENDING`), `λ_Pu`, `Ω_w` (new `ACTION`, no routes named)}
  + 1 function `g_ℓ(w; ℓ_g)` (new `ACTION`; R21 scope split — superseded as material wall, retained as constitutive
  profile) + the 6-postulate structural block. The 2026-07-04 erratum STANDS (no `ρ_br` overcount;
  `NO_OVERCOUNT_ROUTE_A_PENDING`); the `{ρ_B0, χ_c, C_hu}` trio stays Part-VI (guarded in-engine). R22 hardens the
  `μ_R`/`μ_R⁽⁴⁾` distinction so R17's debt cannot be silently faked as discharged. Part I is now COMPLETE: every
  Part-I knob, edge, and debt above is dual-engine-verified where a stage exists.

**Reading:** the free-parameter load is real but heavily **provisional** — most of it is reduction debt with named routes
(all currently `PENDING` on the deferred nonlinear throat), not irreducible freedom. The honest question for Part VII is
how many of R6/R10/R12/R13 actually *land* (discharging debt) vs. come back `CLOSED-NEG` (hardening knobs) — and the
codimension count adjudicates it.

---

## Update protocol (per stage, lightweight)

At each stage's finalize step: (1) add any new parameter rows with dimension + class; (2) add any new relation edges with
the `DERIVED/IMPOSED/PENDING/CLOSED-NEG/CODIM-PROVEN` flag; (3) update the rollup; (4) if a stage *discharges* a pending
route (debt → DERIVED) or *closes* one (→ CLOSED-NEG), record it — those are the results that move the irreducible count.
This doc is the seed for Part VII's calibration map + the codimension count; it is not a substitute for the per-stage
dual-engine dim/able-to-fail checks (those stay in the stage scripts).

**MANDATORY verification (blueprint §5 step 7):** after every update, run a Codex read-only verification of this doc
against the stage's scripts/report (dims match, provenance classes honest, and — the highest-risk check — NO edge
mislabeled: an `IMPOSED` calibration dressed as `DERIVED`/`CODIM-PROVEN` would falsely shrink the irreducible count) and
fold the findings. The register is orchestrator-authored prose, so Codex reviews and the orchestrator folds. Template:
`_scratch/parameter_register_verify_prompt.md`.
