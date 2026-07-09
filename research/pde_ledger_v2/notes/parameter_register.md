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
is a scope FLAG not a quantity; the stage's obligation is edge R23), 009 (pathA_29 Check A flat-slab return
residual — adds the POSTULATED slab-family parameters `{d, ε0, ε1}` (geometry + DC channel transmissions, NOT
medium constants) + edge R24 (Z accounting `DERIVED`); the residual prediction is parameterized by `ε_ℓ` with the
perfect-return limits computed reversible), 010 (pathA_29 Check B slab localization p=2 + NOGO — the COMPLETING stage
of the pathA_29 joint verdict; adds **NO new counted knobs**: `k_warp` is a control-construction symbol used ONLY in
the anti-localizing warp NOGO control (tracked, not counted), and the `1/r²` localization survival on the slab family
discharges NOTHING — R19/`W_slab` and R23 stay `PENDING`; the stage's contribution is edge R25), 011 (pathA_30 II-G1a
frozen-reduction certificate — adds **NO new counted knobs**: `L0` (throat depth) = POSTULATED `ACTION`-geometry
(tracked, not counted, like stage009's `d`), `ℓ_c` (confinement length) INERT here (`δV_conf=0`; tracked, not counted),
`R_mouth` cancels out of the pinch-off root (a construction scale), and `ξ=ħ/(m c_s)` is the R2-family healing length
in the source's no-√2 convention (`DERIVED`, not new); the stage's contributions are edges R26 (frozen-reduction
validity record) + R27 (`ξ≠ℓ_c` firewall)), 012 (pathA_30 II-G1b DtN pole ladder + Robin falsifier — the COMPLETING
stage of the joint `DN_UNITTEST_BC_DEPENDENT`; adds **NO new counted knobs**: `α` (Robin cap admittance, `[α]=L⁻¹`) is
a control-construction symbol used ONLY to build the Robin counterfactual falsifier (tracked, not counted, like
stage010's `k_warp`); the frozen operator/domain/speed are CONSUMED from stage 011 (not re-derived); the stage's
contribution is edge R28 (the D/N boundary is IMPOSED → the `BC_DEPENDENT` landing, an `IMPOSED`/`PENDING` obligation,
NOT a reduction)), **013 (pathA_31 II-G2a breathing harmonic profiles + M/K projection — ⭐ the FIRST Part-II CALIBRATED
stage: adds 3 counted CALIB knobs `{μ_η, T_w, β}`** (the frozen-throat-wall constitutive packet: wall inertia `μ_η`,
wall tension `T_w`, wall inverse-length scale `β`), with `K_η = T_w β²` a DERIVED manifestation (edge R29, NOT counted)
and `r_AL` the dimensionless collective BC ratio (control-construction, tracked-not-counted). ⚠ **`β` is counted, NOT
dressed as geometry** — the source is explicit that geometry alone does not derive the wall constants
(`beta_from_R0: "geometry alone does not derive it"`; `calibration_inputs` names the "K_eta/Tw beta scale"), and
`β·L0 = 37/20` is just `β(calibrated)·L0(=37/20, the stage-011 L/a geometry)`, NOT an independent branch-pin of `β`; the
EARNED content is the STRUCTURE (harmonic-lift profiles, `M_AB`/`K_AB` by real ∫dw operator projection —
`forbidden_fit_flags` computed False via free-symbol-name ancestry), the VALUES are calibration → `BREATHING_CALIBRATED`.
Edge R30 = the named-PENDING nonlinear-throat reduction that would derive the wall response `{μ_η, T_w, β}` from the
medium (reduction debt). `c_S` NOT consumed (matter-sector deferred, `kξ≪1`)), **014 (pathA_31 II-G2b truncation
consistency — adds NO new counted knobs**: the numeric truncation controls `{FLOOR=0.9 (EPS_TRUNC=0.1), N_FINAL=16,
N_CONVERGENCE, BETA_L0_SWEEP}` are method/tolerance parameters (tracked, not counted — convergence/solver knobs, not
physics DOFs; 014 is the FIRST float-bearing stage); the frozen-wall constitutive packet `{μ_η, T_w, β, K_η}` is CONSUMED
from stage 013 via dual-site citation-integrity (NOT re-counted — no double-count); `c_S` NOT consumed (matter-sector
deferred, `kξ≪1`); the stage certifies that 013's 2-mode `{α_a, α_L}` truncation captures the low spectrum of the
combined-basis generalized eigenproblem across a computed window `β_L0∈[0.1,3.0]` ⟺ `K_η/T_w ≲ 2.6`, contributing
structural edge R31 (the truncation validity window — discharges NOTHING)). Seeded knit params
(Part VI, not yet reshaped): pathA_40 cone-lock, pathA_41 NG5.

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
| `d` (slab spacing) | `L` (stage 009: `[τ]=T` via `τ=2d/c_S` dual-engine) | II-B2 (009) postulated slab | `ACTION` (geometry of the executable flat-slab family — NOT a medium constant) | `τ = 2d/c_S` DERIVED from the Helmholtz basis | the slab family is POSTULATED; the real return geometry = the deferred nonlinear brane↔bulk closure |
| `ε0, ε1` (DC channel transmissions) | `1` (stage 009, dual-engine) | II-B2 (009) return channels | `FREE-UNREDUCED` (per-channel; independent) | residual prediction `A_res ∝ ε_ℓ/(1+ε_ℓ)`; `Z=−M0·ε0/(1+ε0)` (R24); strict `ε_ℓ→0⁺` = perfect return (prediction vanishes, orders lift to 2/4, Z→0 — COMPUTED reversibility) | named route: the track-3/Gate-6 nonlinear return would derive the transmissions from the medium (the same wall that discharges R23) — `PENDING` |
| `k_warp` (anti-localizing warp scale) | `L⁻¹` (stage 010: `k_warp·w` dimensionless, dual-engine) | II-B3 (010) NOGO control | `CANDIDATE` (control-construction symbol — the delocalizing half-line warp `μ(w)=exp(2·k_warp·w)` used ONLY to prove the classifier able-to-fail) | — | tracked, **not counted** (like the leg-C probe symbols); no medium provenance — it constructs the falsifiable counterexample, not the physics |
| `L0` (throat depth) | `L` (stage 011: domain `[0,L0]`, cap `R0(L0)=0`, dual-engine) | II-G1a (011) frozen-reduction domain | `ACTION` (straight-reference throat geometry — NOT a medium constant, like stage009's `d`) | cap endpoint SOLVED from the pinch-off `R0(s)=0` on a POSTULATED monotone taper (`R_mouth` cancels) | tracked, **not counted**; the throat depth's dynamical selection is deferred Gate-6/`W_slab` territory (R19-adjacent) |
| `ℓ_c` (confinement length) | `L` (stage 011: in `V_wall(Σ/ℓ_c)`; dual-engine `ξ≠ℓ_c` firewall) | II-G1a (011) | `CANDIDATE` (INERT here — `δV_conf=0` in the frozen `η=0` test) | distinct from the healing length `ξ` (edge R27 firewall) | tracked, **not counted**; a live confinement scale only when the wall is un-frozen (Phase-2 variable-coefficient physics) |
| `α` (Robin cap admittance) | `L⁻¹` (stage 012: `α·c_S` matches `[ω]=T⁻¹`; dual-engine) | II-G1b (012) Robin counterfactual | `CANDIDATE` (control-construction symbol — the Robin cap `∂ₛψ−αψ=0` used ONLY to build the falsifier that proves the D/N determination is not hardcodable: α→0 recovers D/N, α→∞ recovers D/D, numeric α distinct) | — | tracked, **not counted** (like `k_warp` at stage 010); no medium provenance — it constructs the falsifiable counterfactual, not the physics |
| `μ_η` (wall inertia) | `M L⁻¹` (stage 013, dual-engine: `[M_AB]=M` via `4π∫μ_η α² dw`) | II-G2a (013) breathing wall packet | `CALIB` (⭐ FIRST Part-II counted knob — the frozen-throat-wall inertia; a tuned constitutive input) | independent calibration input | reduction debt R30 (a nonlinear-throat solve would derive the wall response from the medium) — `PENDING` |
| `T_w` (wall tension) | `M L T⁻²` (stage 013, dual-engine: enters `K_AB` via `4π∫T_w α'² dw`) | II-G2a (013) breathing wall packet | `CALIB` (⭐ FIRST Part-II counted knob — the frozen-throat-wall tension; a tuned constitutive input) | independent calibration input | reduction debt R30 — `PENDING` |
| `β` (wall inverse-length scale) | `L⁻¹` (stage 013, dual-engine) | II-G2a (013) breathing wall packet | `CALIB` (counted — the wall constitutive scale, source's "K_eta/Tw beta scale"; ⚠ NOT geometry-derived: `beta_from_R0`="geometry alone does not derive it", and `β·L0=37/20` is just `β·L0(=L/a geom)`, not an independent branch-pin of `β`) | `β=√(K_η/T_w)` — one of the two independent DOF among `{T_w, K_η, β}` (with `T_w`); `K_η` is the derived third (R29) | reduction debt R30 — `PENDING` |
| `K_η` (wall stiffness) | `M L⁻¹ T⁻²` (stage 013, dual-engine) | II-G2a (013) breathing operator/K-integrands | `DERIVED` | `= T_w β²` (R29 — manifestation of `{T_w, β}`; value `calibration_tied_to_beta_squared_Tw`) | collapses into `{T_w, β}` (R29 `DERIVED`) — **not** independently counted |
| `r_AL` (collective BC ratio) | `1` (stage 013, dual-engine: `ZERO_DIM`) | II-G2a (013) `α_L` cap normalization | `CANDIDATE` (control-construction — the dimensionless `α_L(L0)=r_AL` collective length-mode ratio) | — | tracked, **not counted** (like `k_warp`/`α`); it parameterizes the collective mode, not the physics |

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
| R24 | `Z = −M0·ε0/(1+ε0)` (signed local source accounting: channel sum `Z_throat + Z_return`, sign certificate computed) | `DERIVED` (ACCOUNTING within the postulated slab family — collapses `Z` into `{M0, ε0}`; does NOT derive the `Z<0` drain-admissibility PREMISE, which stays `Z_is_premise=true`) | II-B2 (009) | stage009 sharpens open-item #9 (bounded residual `p_res=1/3`, per-channel `ε_ℓ`-contingent) — does NOT close it |
| R25 | `1/r²` localization survives the finite slab (both DC-sink completions → normalizable `m=0` zero mode → 3D-radial dsolve → `p=2`) | structural (EARNED-WITHIN-FAMILY verdict, dual-engine + able-to-fail via the anti-localizing warp → `p=3` → `RETURN_NOGO`; **NOT a reduction** — discharges NO knob) | II-B3 (010) | completes the pathA_29 joint `RETURN_RESIDUAL_PREDICTION`; explicitly does **not** discharge R19 (`W_slab` slab-width selection) or R23 (return-cancellation targets) — localization on the FAMILY ≠ the family selected by dynamics |
| R26 | frozen-reduction validity record: `L_s ψ = ψ''+(ω/c_S)²ψ` (const-coeff Helmholtz) is EXACT only on the window `{ρ0'/ρ0=0, √γ0 const, δV_conf=0, ∇Q=0, kξ≪1}`; the Bogoliubov `k⁴` term is a DEFERRED fourth-derivative intrusion `−(ħ²/4m²c_S²)ψ''''` | structural (validity certificate, dual-engine + able-to-fail — a surviving intruding term → `FAIL_OPERATOR_INTRUSION`; **NOT a reduction** — discharges NO knob; the speed is banked R1) | II-G1a (011) | records the conditions under which the frozen-throat longitudinal reduction to a Helmholtz resonator holds; the deferred variable-coefficient / BdG terms are Phase-2 physics (not dropped unconditionally) |
| R27 | `ξ = ħ/(m c_s)` (healing length) ≠ `ℓ_c` (confinement length) — distinct `[L]` symbols, never substituted | structural (dimensional/semantic firewall, dual-engine asserted + able-to-fail; analogous to R22's `μ_R≠μ_R⁽⁴⁾`; **NOT a reduction**) | II-G1a (011) | prevents a silent conflation of the two lengths that would fake a relation; `ℓ_c` INERT here (`δV_conf=0`), `ξ` is R2-family (source's no-√2 convention) |
| R28 | the D/N mouth/cap boundary pair (Dirichlet at `s=0`, Neumann at `s=L0`) is IMPOSED, not derived (`bc_provenance=imposed`, `bc_derivation_emitted=False`) → the joint verdict lands at `DN_UNITTEST_BC_DEPENDENT`, not `DN_UNITTEST_PASS` | `IMPOSED`/`PENDING` (banked calibration input, dual-engine + able-to-fail — forcing the flag True reveals the deferred `DN_UNITTEST_PASS`; **NOT a reduction** — discharges NO knob) | II-G1b (012) | the DtN ladder `Z00=−(ω/c_S)tan(L0ω/c_S)` + half-shifted resonances are a DERIVED consequence of the imposed pair (which collapses no knob); the deferred discharge is an explicit mouth/cap `V_wall` gradient derivation that would earn `DN_UNITTEST_PASS` (analogous to R23's constraint-spec obligation) — Gate-6/Phase-2 territory |
| R29 | `K_η = T_w β²` (wall stiffness = tension × inverse-length-scale²) | `DERIVED` | II-G2a (013) | collapses `K_η` into `{T_w, β}` (a manifestation, `calibration_tied_to_beta_squared_Tw`); so among `{μ_η, T_w, K_η, β}` the counted CALIB set is `{μ_η, T_w, β}` (3), NOT 4 — `K_η` is the derived third of `{T_w, K_η, β}` |
| R30 | nonlinear-throat reduction of the frozen-wall response `{μ_η, T_w, β}` from the medium | `PENDING` | II-G2a (013) | debt: would derive the calibrated wall constitutive packet (inertia/tension/scale) from the deferred nonlinear throat interior — the same wall as R10/R23's discharge; until then `{μ_η, T_w, β}` are CALIB (tuned in the frozen `G=c=c_s=1` packet `T_w=K_η=1`) |
| R31 | truncation validity window: 013's 2-mode `{α_a, α_L}` collective closure is a clean truncation of the combined-basis generalized eigenproblem `K v=ω²M v` (over `{α_a,α_L,g_1..g_N}`) ONLY for order-unity wall stiffness `K_η/T_w ≲ 2.6` (computed sweep edge — the overlap floor `o_k≥0.9` passes for `β_L0∈[0.1,3.0]`, FAILS for `β_L0≥5`; `o_1=0.99311,o_2=0.98776` at the physical `β_L0=37/20`, N-converged over `N=4/8/12/16`) | structural (validity certificate, dual-engine + able-to-fail — a gamed threshold / an everywhere-pass window is caught, per the pathA_31 v1 rejection scar; ⚠ the modal overlap does NOT guard profile-correctness — `constant_one` PASSES it, so this certifies truncation-consistency NOT profile-correctness (that is 015's HF + 013's residual); **NOT a reduction** — discharges NO knob) | II-G2b (014) | bounds where the calibrated breathing closure (013's `M_AB`/`K_AB`) is trustworthy: sharp walls (`K_η/T_w` large) need >2 modes; the phonon-limit / BdG `k⁴` stays deferred (`kξ≪1`); the 014 numeric controls `{FLOOR, N_FINAL, N_CONVERGENCE, BETA_L0_SWEEP}` are method/tolerance (tracked, not counted) |

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
- **Stage 009 (pathA_29 Check A): the slab-family package `{d, ε0, ε1}`, counted honestly.** `d` = postulated
  executable-family geometry (`ACTION`, not a medium constant); `ε0, ε1` = per-channel DC transmissions,
  `FREE-UNREDUCED` with the named Gate-6/track-3 route (the nonlinear return would derive them — the same wall as
  R23's discharge). `Z` collapses into `{M0, ε0}` via R24 (accounting; the `Z<0` premise stays a premise). The
  falsifiable residual prediction (`p_res=1/3`) is `ε_ℓ`-parameterized with COMPUTED perfect-return reversibility —
  the prediction is contingent, not baked in.
- **Stage 010 (pathA_29 Check B): zero new counted knobs — the COMPLETING stage of the pathA_29 fold.** The `1/r²`
  localization survives the slab (both DC-sink completions → normalizable `m=0` zero mode → `p=2`, edge R25),
  genuinely able-to-fail (the anti-localizing warp → `p=3` → `RETURN_NOGO`). The only new symbol, `k_warp`
  (`[k_warp]=L⁻¹`), is a control-construction `CANDIDATE` (tracked, not counted — it builds the falsifiable
  counterexample, not the physics). R25 is **structural, not a reduction**: the localization discharges no knob and
  explicitly leaves R19 (`W_slab` slab-width selection) and R23 (return-cancellation targets) `PENDING` — localization
  on the postulated FAMILY is not the family selected by dynamics.
- **Stage 011 (pathA_30 Check II-G1a): zero new counted knobs — the frozen-reduction certificate.** The frozen wall
  (`η=0`) collapses the longitudinal operator to const-coeff Helmholtz `L_s = ψ''+(ω/c_S)²ψ` (edge R26), with the speed
  BANKED (R1 at `ρ*`, not re-derived) — so the stage adds **no reduction**. New tracked-not-counted symbols: `L0`
  (throat depth, `ACTION`-geometry like `d`), `ℓ_c` (confinement length, INERT here — `δV_conf=0`); `R_mouth` cancels
  out of the pinch-off root; `ξ=ħ/(m c_s)` is R2-family (`DERIVED`). The two structural edges — R26 (the validity
  window: the reduction is exact only under `{ρ0'/ρ0=0, √γ0 const, δV_conf=0, ∇Q=0, kξ≪1}`, the BdG `k⁴` term a
  DEFERRED 4th-order intrusion) and R27 (the `ξ≠ℓ_c` firewall, analogous to R22) — **discharge nothing**; both are
  able-to-fail (a surviving intruding term → `FAIL_OPERATOR_INTRUSION`; a conflation → firewall fires). The D/N
  boundary determination + the `BC_DEPENDENT` landing are stage 012's (banked calibration, `bc_derivation_emitted=False`).
- **Stage 012 (pathA_30 Check II-G1b): zero new counted knobs — the COMPLETING stage of the pathA_30 fold.** The
  `dsolve` of stage 011's cited `L_s` → D/N BVP → the DtN `Z00=−(ω/c_S)tan(L0ω/c_S)`, the half-shifted resonance ladder
  `ω_n=πc_S(n+½)/L0`, `R_rt=1`, and the Robin counterfactual. `α` (Robin cap admittance, `[α]=L⁻¹`) is a
  control-construction `CANDIDATE` (tracked, not counted — like `k_warp`; it builds the falsifier, not the physics);
  the operator/domain/speed are CONSUMED from stage 011. R28 is **imposed, not a reduction**: the D/N boundary pair is
  a banked calibration input (`bc_provenance=imposed`, `bc_derivation_emitted=False`) → the joint verdict lands at
  `DN_UNITTEST_BC_DEPENDENT`, discharging NO knob; its deferred discharge (an explicit mouth/cap `V_wall` derivation
  earning `DN_UNITTEST_PASS`) is Gate-6/Phase-2 territory. **Completes the joint `DN_UNITTEST_BC_DEPENDENT` = (011:
  REDUCTION_CERTIFIED) ∧ (012: DtN ladder + BC_DEPENDENT landing).**
- **⭐ Stage 013 (pathA_31 II-G2a): the FIRST Part-II CALIBRATED stage — adds 3 counted CALIB knobs `{μ_η, T_w, β}`.**
  The ℓ=0 breathing reduces to a 2-mode collective closure: DERIVED harmonic-lift profiles `α_a, α_L` (proven by the
  `𝓛₀[α]=0` residual), `M_AB`/`K_AB` by real ∫dw operator projection (`forbidden_fit_flags` computed False via
  free-symbol-**name** ancestry, `operator_projection_not_static_Hessian` — NOT the legacy static Hessian), and the
  dynamical-EOM LHS. **The STRUCTURE is EARNED; the VALUES are CALIBRATED** → `BREATHING_CALIBRATED` (the 013 component
  of the joint 3-stage verdict; 014 = truncation, 015 = legacy-Hessian + HF force). The counted set is the frozen-wall
  constitutive packet `{μ_η` (inertia, `M L⁻¹`), `T_w` (tension, `M L T⁻²`), `β` (inverse-length scale, `L⁻¹`)}`, with
  `K_η = T_w β²` a DERIVED manifestation (R29 — so 3 counted, not 4). ⚠ **`β` is counted, NOT dressed as geometry** —
  the source is explicit that geometry alone does not derive the wall constants, and `β·L0=37/20` is
  `β(calibrated)·L0(=stage-011 L/a geometry)`, not an independent branch-pin of `β` (the conservative, non-shrinking
  read; the source-map's earlier "{μ_η,T_w}=2, β geometry" guess did NOT survive the provenance check). `r_AL`
  (dimensionless collective BC ratio) is a control-construction `CANDIDATE` (tracked, not counted, like `k_warp`/`α`).
  R30 is the named-PENDING nonlinear-throat reduction that would derive `{μ_η, T_w, β}` from the medium (reduction debt).
  `c_S` NOT consumed (matter-sector deferred, `kξ≪1`). **These are the first tuned knobs since Part I's `{ħ, m, K, ρ0}`
  + the brane freeze package — the gravity sector's algebra was knob-free through 008–012; the calibration enters with
  the breathing-mode wall response.**
- **Stage 014 (pathA_31 II-G2b): zero new counted knobs — the truncation-consistency component (2/3) of
  `BREATHING_CALIBRATED`; the FIRST float-bearing stage.** The combined-basis generalized eigenproblem `K v=ω²M v` over
  `{α_a,α_L,g_1..g_N}` (013's two collective modes + N sine lane modes) certifies that 013's 2-mode closure captures the
  two lowest generalized modes to the predeclared overlap floor `o_k≥0.9` (`o_1=0.99311, o_2=0.98776` at the physical
  `β_L0=37/20`), across a COMPUTED validity window `β_L0∈[0.1,3.0]` (⟺ order-unity wall stiffness `K_η/T_w≲2.6`; genuine
  FAIL rows at `β_L0≥5`) and N-converged over `N=4/8/12/16` — edge R31. The numeric truncation controls `{FLOOR=0.9
  (EPS_TRUNC=0.1), N_FINAL=16, N_CONVERGENCE, BETA_L0_SWEEP}` are method/tolerance parameters (tracked, not counted); the
  wall packet `{μ_η, T_w, β, K_η}` is CONSUMED from stage 013 (dual-site citation-integrity, no double-count); `c_S` NOT
  consumed. ⚠ **Honest caveat carried (not softened, per falsification-first):** the modal overlap does NOT guard
  profile-correctness — the `constant_one` (wrong) profile PASSES the overlap (`o_1=1.0, o_2=0.974`); the profile guard is
  013's `𝓛₀[α]=0` residual + 015's Hellmann–Feynman mismatch. The pathA_31 v1 REJECTION scar (a *gamed* truncation
  threshold) is ABSENT — the floor + window are genuinely computed and able-to-fail (adversarial per-tooth ablation
  confirmed; 2 vacuous able-to-fail legs were found + remediated to genuine coupling).
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
