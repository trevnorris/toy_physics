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

### ⭐ Two kinds of "free": declared universal constants vs reduction debt (user, 2026-07-09)
The irreducible set is NOT monolithic — split it into two first-class lists, because they score differently:
- **(a) DECLARED universal constants** — irreducible *by design*, and honest (every theory bottoms out in free constants:
  GR has `G`(+`Λ`), the SM ~19–26). Candidates = the medium primitives `{ħ, m_GNLS, K, ρ0}`. ⭐ **The scoring question for
  these is NOT "why isn't it derived" — it is "does the analog need FEWER fundamental constants than the target GR+EM need?"**
  If ~4 medium constants reproduce GR+EM structure across four sectors, that is the headline win, not a liability. ⚠ *Declaring*
  a constant fundamental is a MODELING DECISION (a user call), distinct from bookkeeping: e.g. do we accept `ħ` as fundamental,
  or keep chasing `HBAR_FREE_SUBSTRATE_RELATION_MISSING`? (Even the current primitives carry open GAPs — `ħ` provenance, `m`
  emergence — so "which do we declare fundamental vs keep reducing" is itself an open checkpoint decision.)
- **(b) REDUCTION-DEBT tunings** — `CALIB`/`FREE-UNREDUCED` knobs that are stand-ins for math/sim not yet done (the frozen-wall/
  throat packet `{μ_η,T_w,β,Vp0/ℓ_c,T_Ω,β₂}` via R30/R33/R36 — all siblings of the ONE deferred throat-interior solve; the NG5
  drift `{ρ_B0,χ_c,C_hu}` via R12/R13). **This is where the "too many knobs" risk actually lives.** ⭐ The sim is `not-yet` not
  `never` (user, 2026-07-09): the ledger is being built *sim-READY* precisely so that when solver tractability + model capability
  catch up, these convert CALIB→DERIVED in one stroke — they are placeholders with a known address, not dead ends.
- **⭐ The sharpest health metric is NOT "how many knobs" — it is "how many irreducible knobs have NO named reduction route"**
  (pure irreducible freedom that eats predictive surplus and cannot improve), vs debt-with-a-route (honest, pending). Track the
  route-less subset explicitly (currently `{ρ_B0,χ_c,C_hu}` has no *registered* reduction — the real liability — vs the throat
  packet which has R30).
- **⭐ MIDWAY KNOB AUDIT — scheduled checkpoint (user, 2026-07-09): run AFTER the entire Part-II gravity sector closes** (NOT now,
  NOT deferred to Part VII). A dry-run of the Part-VII codimension technique (pathA_40 `Δr=2`, scaled) over Parts I–II + a tally
  of held-out predictions vs the *irreducible* (route-less) count + the (a)/(b) split above — turning "are we accumulating too many
  knobs" into a number early enough to course-correct or hit a clean no-go. (Also: at that checkpoint, list every value that
  genuinely CANNOT be derived — the declared universal constants — separately from the derivable-with-more-math/sim debt.)

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
| `c_s0` | `L T⁻¹` | I-1/I-2; **LIVE units carrier at II-G4a (018)** | `DERIVED` | `= √(5K ρ0⁴/m_GNLS)` | DERIVED from `{K, ρ0, m_GNLS}` (R1); ⚠ at 018 it is the units-restoring symbol in the fingerprint physical coeffs (`z=aω/c_s`) — cited PROVENANCE, NOT a consumed value (the earned rationals + χ_Q are `c_s`-free); distinct from the frozen-wall `c_S` (011–017) |
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
| `ε0, ε1` (DC channel transmissions) | `1` (stage 009, dual-engine) | II-B2 (009) return channels; **re-expressed as `Z0_ret/Z1_ret` at II-G5b (023)** | `FREE-UNREDUCED` (per-channel; independent) | residual prediction `A_res ∝ ε_ℓ/(1+ε_ℓ)`; `Z=−M0·ε0/(1+ε0)` (R24); strict `ε_ℓ→0⁺` = perfect return (prediction vanishes, orders lift to 2/4, Z→0 — COMPUTED reversibility). ⚠ **stage023's `Z0_ret/Z1_ret` are the COORDINATE ALIASES of these SAME 2 dofs** (`ε_ℓ=Z_ℓ,ret/K_ℓ` invertible once `K_ℓ` fixed) — NOT 2 new free dofs; R42 SHARPENS this debt to "a 2-dim return-admittance nullspace survives the linear return law" (no new dofs) | named route: the track-3/Gate-6 nonlinear return would derive the transmissions from the medium (the same wall that discharges R23; R42 = the sharpened form) — `PENDING` |
| `K0c` (ℓ=0 collective stiffness) + the ℓ=1 stiffness sector `{K_eta, T_Omega}` (appearing as `K1 = K_eta+2·T_Omega`) | `M T⁻²` (stage 023, dual-engine, pathA_34 convention; scripts' `(L,M,T)`-tuple `(0,1,−2)` — ⚠ the register header order is `[L,T,M]`, in which it reads `(0,−2,1)`) | II-G5b (023) ℓ=0/1 named-constraint stiffnesses (the `named_constraints` rows; the `GENERATOR_DOFS` treat `K0c`, `K_eta`, `T_Omega` as separate directions, `K_eta`/`T_Omega` entering the physics only via `K1=K_eta+2·T_Omega`) | `FREE-UNREDUCED` **PENDING scalar-reduction — COUNTED as reduction debt** (flagged; per the register's rule pending debt stays counted until DERIVED. ⚠ these are **NOT** the registered 013/017 knobs [dims conflict, below] and **NOT** new `CALIB` knobs — CALIB set stays 6; they are the ℓ=0/1 effective-stiffness reduction-debt) | pathA_34-convention EFFECTIVE ℓ=0/1 stiffnesses; ⚠ **NOT identified with the raw 013/017 densities** — the dims `M T⁻²` do NOT match registered stage013 `K_η=T_wβ²` (`M L⁻¹T⁻²`) or stage017 `T_Ω` (`M L⁻³T⁻²`); stage017 records `K_η=T_wβ²` NON-transferable across the volume-vs-line convention (the stage016 lesson) | reduction debt (**R42-family**): an explicit scalar-reduction (profile+measure) of the ℓ=0 Gate-2 collective (§8.2, stage013) + the ℓ=1 §9.4 harmonic (stage017) sectors to the calibrated wall packet would DISCHARGE this (collapse them into `{μ_η,T_w,β,T_Ω}`) — `PENDING` (no equation yet; do NOT assert DERIVED; counted-as-debt meanwhile) |
| `k_warp` (anti-localizing warp scale) | `L⁻¹` (stage 010: `k_warp·w` dimensionless, dual-engine) | II-B3 (010) NOGO control | `CANDIDATE` (control-construction symbol — the delocalizing half-line warp `μ(w)=exp(2·k_warp·w)` used ONLY to prove the classifier able-to-fail) | — | tracked, **not counted** (like the leg-C probe symbols); no medium provenance — it constructs the falsifiable counterexample, not the physics |
| `L0` (throat depth) | `L` (stage 011: domain `[0,L0]`, cap `R0(L0)=0`, dual-engine) | II-G1a (011) frozen-reduction domain | `ACTION` (straight-reference throat geometry — NOT a medium constant, like stage009's `d`) | cap endpoint SOLVED from the pinch-off `R0(s)=0` on a POSTULATED monotone taper (`R_mouth` cancels) | tracked, **not counted**; the throat depth's dynamical selection is deferred Gate-6/`W_slab` territory (R19-adjacent) |
| `ℓ_c` (confinement length) | `L` (stage 011: in `V_wall(Σ/ℓ_c)`; dual-engine `ξ≠ℓ_c` firewall) | II-G1a (011); **LIVE at II-G2c (015)** | `CANDIDATE` at 011 (INERT — `δV_conf=0` in the frozen `η=0` test); **manifestation of the `Vp0/ℓ_c` force-scale at 015** (`δV_conf≠0`) | distinct from the healing length `ξ` (edge R27 firewall); at 015 appears ONLY via the ratio `Vp0/ℓ_c` | tracked, **not counted separately** — becomes live at 015 but only through the counted force-scale `Vp0/ℓ_c` (below), not as an independent length |
| `Vp0/ℓ_c` (breathing driving-force scale) | `M L⁴ T⁻²` (Codex-verified 2026-07-08 against the 013 dim block: the EOM RHS force-density `ρ*·Vp0/ℓ_c` has `[F_A^(HF)/w]=M T⁻²` per the `[K_AB Q]=[F_A]` balance, and with the frozen density `[ρ*]=L⁻⁴` ⇒ `[Vp0/ℓ_c]=M L⁴ T⁻²`, `[Vp0]=M L⁵ T⁻²`) | II-G2c (015) HF confinement drive | `CALIB` (⭐ the SECOND Part-II counted knob after 013 — the wall confinement potential slope that DRIVES the breathing; the RHS 013 deferred) | `Vp0` (wall confinement potential radial slope) + the now-live `ℓ_c` are its manifestations (they appear ONLY as the ratio `Vp0/ℓ_c`); `ρ*` is CONSUMED (frozen density, stage005/011) | reduction debt **R33** (a nonlinear-throat solve — the same deferred interior as R30's wall response — would derive the confinement drive) — `PENDING` |
| `α` (Robin cap admittance) | `L⁻¹` (stage 012: `α·c_S` matches `[ω]=T⁻¹`; dual-engine) | II-G1b (012) Robin counterfactual | `CANDIDATE` (control-construction symbol — the Robin cap `∂ₛψ−αψ=0` used ONLY to build the falsifier that proves the D/N determination is not hardcodable: α→0 recovers D/N, α→∞ recovers D/D, numeric α distinct) | — | tracked, **not counted** (like `k_warp` at stage 010); no medium provenance — it constructs the falsifiable counterfactual, not the physics |
| `μ_η` (wall inertia) | `M L⁻¹` (stage 013, dual-engine: `[M_AB]=M` via `4π∫μ_η α² dw`) | II-G2a (013) breathing wall packet | `CALIB` (⭐ FIRST Part-II counted knob — the frozen-throat-wall inertia; a tuned constitutive input) | independent calibration input | reduction debt R30 (a nonlinear-throat solve would derive the wall response from the medium) — `PENDING` |
| `T_w` (wall tension) | `M L T⁻²` (stage 013, dual-engine: enters `K_AB` via `4π∫T_w α'² dw`) | II-G2a (013) breathing wall packet | `CALIB` (⭐ FIRST Part-II counted knob — the frozen-throat-wall tension; a tuned constitutive input) | independent calibration input | reduction debt R30 — `PENDING` |
| `β` (wall inverse-length scale) | `L⁻¹` (stage 013, dual-engine) | II-G2a (013) breathing wall packet | `CALIB` (counted — the wall constitutive scale, source's "K_eta/Tw beta scale"; ⚠ NOT geometry-derived: `beta_from_R0`="geometry alone does not derive it", and `β·L0=37/20` is just `β·L0(=L/a geom)`, not an independent branch-pin of `β`) | `β=√(K_η/T_w)` — one of the two independent DOF among `{T_w, K_η, β}` (with `T_w`); `K_η` is the derived third (R29) | reduction debt R30 — `PENDING` |
| `K_η` (wall stiffness) | `M L⁻¹ T⁻²` (stage 013, dual-engine) | II-G2a (013) breathing operator/K-integrands | `DERIVED` | `= T_w β²` (R29 — manifestation of `{T_w, β}`; value `calibration_tied_to_beta_squared_Tw`) | collapses into `{T_w, β}` (R29 `DERIVED`) — **not** independently counted |
| `r_AL` (collective BC ratio) | `1` (stage 013, dual-engine: `ZERO_DIM`) | II-G2a (013) `α_L` cap normalization | `CANDIDATE` (control-construction — the dimensionless `α_L(L0)=r_AL` collective length-mode ratio) | — | tracked, **not counted** (like `k_warp`/`α`); it parameterizes the collective mode, not the physics |
| `{κ, χ, σ_a, σ_L}` (legacy `E_geom` Hessian constants) | energy-Hessian coefficients (stage 015: `H_legacy={aa:χ²κ+σ_a, aL:−χκ, LL:κ+σ_L}`, dual-engine) | II-G2c (015) legacy static-energy pattern basis | **NOT counted afresh** — the legacy-Hessian PATTERN basis the operator projection RECOVERS (a re-parameterization mapping to the calibrated `{μ_η,T_w,β}` closure) | ABSENT from 013's `M_AB`/`K_AB` (013 free-symbol firewall); the structure gate matches the FORM (rank/sign/zero-pattern), not the entries (`full_matrix_fit=False`) | tracked, **not counted** — counting them would double-count the already-calibrated `(a,L)` closure (edge R32); they are the comparison basis of the structure recovery, not independent physics DOFs |
| `T_Ω` (ℓ=2 angular-stiffness density) | `M L⁻³ T⁻²` (stage 017, dual-engine, pathA_32 VOLUME-density convention: enters `K₂=∫[T_w β₂'²+(K_η+6T_Ω)β₂²] a²dwdΩ`; the sourced-`[T_Ω]` corrupt probe fires `FAIL_DIMENSIONAL`) | II-G3b (017) ℓ=2 grouped-lane stiffness (first-appeared 016, count DEFERRED to here) | `CALIB` (⭐ new Part-II counted knob — the wall's ANGULAR/rotational stiffness, a genuinely new physical DOF absent from the ℓ=0 breathing sector; **calibrated NOT derived** — the source verdict ladder makes `ISOTROPY_PASS` require `T_Ω` derived from `R0`, and it is not → `ISOTROPY_CALIBRATED`) | independent calibration input (the `K₂` angular term `6·T̃_Ω`; the `λ_m=6` factor is 016's earned covariance) | reduction debt **R36** (a Gate-1 `R0`-support-equation derivation of `T_Ω`/`β₂` would earn `ISOTROPY_PASS`) — `PENDING` |
| `β₂(w)` (ℓ=2 radial profile) | `1` dimensionless (`β₂'=L⁻¹`) (stage 017, dual-engine, pathA_32 convention) | II-G3b (017) ℓ=2 radial deformation shape (first-appeared 016, count DEFERRED to here) | `CALIB` (⭐ new Part-II counted knob — the FROZEN ℓ=2 deformation profile; **calibrated NOT derived** from the Gate-1 `R0` support equation, the primary CALIBRATED-not-PASS reason, report :5) | independent calibration input; enters the response ONLY via the radial scalars `M̃/K̃/T̃_Ω` (∫·β₂²) | reduction debt **R36** — `PENDING` |
| `M̃, K̃, T̃_Ω` (ℓ=2 radial scalars) | `M`, `M T⁻²`, `M T⁻²` (stage 017, dual-engine; report :54–55) | II-G3b (017) grouped-lane `M₂`/`K₂` | `DERIVED` | radial reductions **= ∫ density·β₂² dV** (edge **R35**): `M̃=∫μ_η β₂²`, `K̃=∫[T_w β₂'²+K_η β₂²]`, `T̃_Ω=∫T_Ω β₂²` — manifestations of `{μ_η, T_w/K_η, T_Ω, β₂}` | collapse into the densities × `β₂` (R35 `DERIVED`) — **not** independently counted (the ℓ=2 analogue of R29's manifestation pattern) |
| `{B̃0,B̃2,B̃4}` support + `{Z̃0,Z̃2,Z̃4}` mixed/Maxwell scalars | D-lane support scalars (stage 017, dual-engine: `D0=K̃+6T̃_Ω−(B̃0+Z̃0)`, `D2=−(M̃+B̃2+Z̃2)`, `D4=−(B̃4+Z̃4)`; dims by order `M T⁻²`/`M`/`M T²`) | II-G3b (017) ℓ=2 port-kernel D-lane support/Maxwell content | `CALIB` (frozen support/Maxwell couplings of the ℓ=2 wall mode, FIRST load-bearing here) — ⚠ **tracked, irreducibility a Part-VII adjudication:** the isotropy VERDICT is value-INDEPENDENT of them (raw-D=0 holds symbolically for any B̃/Z̃ — they are structural port-kernel inputs, NOT tuned-for-isotropy), and the `Z̃` "Maxwell" couplings are pinned DOWNSTREAM (the ℓ=2 port kernel is consumed by 018–024 + the EM knit/Gate-4 `54/5` normalization) | the exported ℓ=2 port-kernel support content (→ 018–021 + 022/023 + 024) | downstream-pinning route (Gate-4/5 + EM knit) — `PENDING`; NOT counted in the clean `{T_Ω, β₂}` headline pending that adjudication |

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
| R32 | legacy-Hessian structure recovery: the operator-projected `M_AB`/`K_AB` (013) reproduce the STRUCTURAL SIGNATURE of the legacy static energy Hessian `H_legacy=∂²E_geom/∂Q²` (M pos-def by exact-identity Sylvester certificates; K symmetric; `K_aL<0`; rank + zero-pattern match) | structural (recovery certificate, dual-engine + able-to-fail — the `M_aa→−M_aa` / `K_aL→−K_aL` probes flip the gate; the positivity is discharged by exact-identity monotonicity certificates, NOT `is_positive` and NOT a `X≡X` form-equality; **NOT a reduction** — discharges NO knob) | II-G2c (015) | certifies the `(a,L)` closure is RECOVERED (has the form a physicist would write from a static energy functional), not re-postulated; the legacy `{κ,χ,σ_a,σ_L}` are the pattern basis (re-parameterization of the calibrated `{μ_η,T_w,β}`, NOT counted afresh — else double-count the closure); the match is structural (`full_matrix_fit=False`), not entrywise |
| R33 | nonlinear-throat reduction of the breathing driving-force scale `Vp0/ℓ_c` (the wall confinement potential slope) from the medium | `PENDING` | II-G2c (015) | debt: the same deferred nonlinear-throat interior as R30 (which would derive the wall response `{μ_η,T_w,β}`) would derive the confinement drive `Vp0/ℓ_c`; until then it is `CALIB` (the EOM RHS `F_A^(HF)` amplitude, tuned) — a sibling of R30's wall-response debt |
| R34 | ℓ=2 SO(3) covariance provenance: the five real ℓ=2 harmonics form one orthonormal SO(3)-irrep (`Gram=I₅`) with a single COMPUTED `−Δ_S²` eigenvalue `λ_m=ℓ(ℓ+1)=6` for every m (Laplace–Beltrami + Rayleigh + eigenfunction residual); the K₂ angular stiffness `K₂=K̃+λ_m·T̃_Ω` uses that computed eigenvalue | structural (covariance theorem, dual-engine + able-to-fail — the `forced_eigenvalue_probe` rejects a typed coefficient → `FAIL_NOT_COVARIANT`; **NOT a reduction** — earned angular structure, discharges NO knob) | II-G3a (016) | the earned angular content the grouped-lane isotropy (017) rides; (backfilled into the edges table 2026-07-09 — was recorded in the stage-016 rollup prose only) |
| R35 | ℓ=2 radial-scalar reduction: the grouped-lane scalars **= ∫ density·β₂² dV** — `M̃=∫μ_η β₂²`, `K̃=∫[T_w β₂'²+K_η β₂²]`, `T̃_Ω=∫T_Ω β₂²` (pathA_32 VOLUME-density convention on `a²dwdΩ`) | `DERIVED` (radial reductions — collapse `M̃/K̃/T̃_Ω` into the densities `{μ_η, T_w/K_η, T_Ω}` × the frozen profile `β₂`; the ℓ=2 analogue of R29's `K_η=T_wβ²` manifestation pattern; dual-engine dim-verified `[M̃]=M`, `[K̃]=[T̃_Ω]=M T⁻²`) | II-G3b (017) | so the counted ℓ=2 CALIB inputs are the underlying `{T_Ω, β₂}` (+ 013's `{μ_η,T_w}`, `K_η`=R29), NOT the reduced scalars — avoids double-counting the port kernel |
| R36 | Gate-1 `R0`-support-equation derivation of the ℓ=2 frozen calibration `{β₂(w), T_Ω}` from the straight-reference throat `R0(w)` | `PENDING` (the `ISOTROPY_PASS` target — the source verdict ladder requires `β₂` derived from `R0` AND `T_Ω` derived, not free-calibrated; until then `{β₂, T_Ω}` are `CALIB` → `ISOTROPY_CALIBRATED`, report :5; **NOT yet a reduction** — discharges NO knob until earned) | II-G3b (017) | the debt that would flip `ISOTROPY_CALIBRATED → ISOTROPY_PASS`; a sibling of R30/R33's deferred-interior debts but at the ℓ=2 support-equation level (the frozen wall profile + radial/support scalars, not the nonlinear throat) |
| R37 | outgoing-DtN ℓ=2 Hankel-fingerprint provenance: `Ŷ₂ᵒᵘᵗ=−3/Λ₂ᵒᵘᵗ`, `Λ₂ᵒᵘᵗ=z·h₂⁽¹⁾′/h₂⁽¹⁾`, `z=aω/c_s`, `h₂⁽¹⁾=j₂+i·y₂`, series-expands to the DERIVED fingerprint `u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴`, `v₅=a⁵/27c_s⁵` (the `27` = `1/v₅ᶻ` COMPUTED) + the sign `χ_Q=+1` outgoing / `−1` incoming (from `j₂±i·y₂`) + the standing-branch contrast (`Λ_stand(0)=+2`, no radiating term) | structural (fingerprint provenance, dual-engine + able-to-fail — a mutated Hankel derivation flips a coefficient, an incoming branch flips χ_Q, a corrupt `a`/`c_s` power fires the units-dim leg; **NOT a reduction** — earned exterior-wave structure, discharges NO knob; like R34 for 016's covariance theorem) | II-G4a (018) | the earned exterior radiative signature the quadrupole normalization (019–021) rides; ⚠ `c_s` is a **units carrier** (R1-`DERIVED`, cited PROVENANCE), NOT a consumed value — the earned rationals + `χ_Q` are `c_s`-free; the port scalars `N_n/D_n` are 019's deferred branch data; the `χ_Q` reconciliation with pathA_22b (`≈0.712`, same name / different computation) is a tracked **Part-VII** item (blueprint §8), NOT merged here |
| R38 | squared-denominator prefactor-algebra provenance: `P(ω)=D₀·N(ω)/D^cons(ω)²`, `N=N₀+N₂ω²+N₄ω⁴`, `D^cons=D₀+D₂ω²+D₄ω⁴`, series-expands to `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`, `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` (SERIES-EXTRACTED off the actual series; the `−2D₂N₀` factor-of-2 is the squared-denominator signature — plain `N/D` gives only `−D₂N₀`, provably wrong → `FAIL_PREFACTOR_ALGEBRA`) | structural (prefactor-algebra provenance, dual-engine + able-to-fail — corrupting the object to plain `N/D` fires the N/D self-check, a mutated `D_n`/`N_n` power flips a series-extracted coefficient, the swapped-in correct-object positive-control flips the discriminator to True; **NOT a reduction** — earned port-agnostic abstract algebra, discharges NO knob; like R37 for 018's fingerprint / R34 for 016's covariance) | II-G4b (019) | the prefactor structure the quadrupole normalization (020/021) rides; ⚠ the abstract `D₀,D₂,D₄` are 017's exported D-lanes + the `N₀,N₂,N₄` are `build_port_moments`' deferred Gate-6 port N-moments — carried as generic symbols, **PROVENANCE only** (the algebra is PORT-AGNOSTIC, holds for any nonzero `D0..N4`; there is NO checkable cross-stage relation, NO dual-site — the same provenance-only landing as 018, `build_port_moments` is emitted-but-never-checked); the numerical `(D_n,N_n)` remain Gate-6 sim-deferred; 019 is **units-FREE** (no `c_s`/`a`/`G`/`μ̂₀`, no dim leg — the `[P₀^phys]=1` closure is 021's) |
| R39 | `54/5=2·27/5` provenance-partition + the CALIBRATED-verdict provenance: the assembled ℓ=2 quadrupole magnitude `54Gc_s⁵/(5a⁵c⁵)` decomposes as `54/5 = 2·27/5` (a SymPy-VERIFIED rational identity BOUND to `target_rhs`/`v5_slot`, NOT a typed string nor a bare-literal tautology), where the `27` is `derived_in_gate` (018's fingerprint `1/v₅ᶻ`, cited) and the `2/5`+`G` are `external_bridge_input` (the GR Burke–Thorne `2G/5c⁵`); a 4-way PROVENANCE partition (`classify_provenance` tag-dominance `deferred>external>derived>convention`) classifies the assembled `54/5` as `external_bridge_input` (external DOMINATES the mixed earned+external tags), so the source-faithful verdict (`if g_class==derived AND mag_class==derived → QUAD_PASS; else → QUAD_CALIBRATED`) lands `QUAD_CALIBRATED` not `QUAD_PASS` | structural (provenance-partition + CALIBRATED-verdict provenance, dual-engine + able-to-fail — `3f_partition_mislabel` fires BOTH directions [G external→derived AND 27 derived→external], the dominance classifier is PROVEN by a truth-table + key-class + tag-mutation tooth [`partition_ok` alone is trivially True], the CALIBRATED verdict is proven to READ provenance by a QUAD_PASS positive control AND a REQUIRED MIXED control [one external/one derived → CALIBRATED, which an inverted rule would mislabel]; **NOT a reduction** — a classification landing / the honest CALIBRATED verdict, discharges NO knob; like R37/R38 for 018/019's fingerprint/prefactor) | II-G4c (020) | the leg that LANDS the CALIBRATED headline; ⚠ **the verdict is PROVENANCE-driven, NOT G→λG-invariance-driven** — a SEPARATE g-invariance diagnostic exposes the invariance-only TRAP (`54/5` is G-invariant yet calibrated, so an invariance-only test would MISLABEL it as earned; `invariance_only_trap_catches_54_over_5=True`). ZERO new counted knobs: `G=GENUINE_BLOCKED`/`external_bridge_input` (the thing the PDE does NOT deliver, already registered), `c`=the GR/PN vacuum light speed = the light cone `c_γ` in its GR-units-bridge role (`P₀ ∝ c_s⁵/c⁵ = 1/λγ⁵`; cited benchmark, NOT a fresh knob), the `2/5`=GR bridge, the `27`=018's derived fingerprint. 020 is UNITS-BEARING (`{c_s,a,c,G}` live in the ALGEBRA) but does NO dimensional-homogeneity gate (`[P₀^phys]=1` is 021's) — enforced by a free-symbol allowlist (no `μ̂₀` symbol; the legacy `gamma_quad_eff` μ̂₀ string DROPPED) + a structural no-dimensional-dependency cut (the 020-local verdict's signature/bytecode carries no dim parameter/dependency; the `.wl` checks its definition arity/text); consumption is PROVENANCE (018's χ_Q=+1/`27` enter the self-contained equivalence bridge, 019's `P0=N0/D0` is `Gamma5`-definitional-only, 017's D-lanes — NO dual-site) |
| R40 | μ̂₀-free `[P₀^phys]=1` dimensional closure: `[P₀^phys]=(c_s/a)²·(N₀/D₀)` is dimensionless from the SOURCED port dims `[N₀]=L⁻¹M` (density/continuity-port numerator) and `[D₀]=L⁻¹T⁻²M` (the carried reduced static conservative denominator `D₀=K−B₀−Z₀`) — `[P₀_raw]=[N₀/D₀]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`, `dimensional_ok=True`; the μ̂₀ carrier is NON-verdict (`participates_in_verdict=False`, `μ̂₀∉raw_symbol_dims`, the verdict read-set excludes `μ̂₀`/`μ_dim`/`homogeneity_pass`). The natural-units trap (the handoff `P₀=N₀/D₀` silently drops `(c_s/a)²`) is CAUGHT (`3d` `FAIL_DIMENSIONAL`); the corrupt-`[N₀]`/`[D₀]`/`[c_s]` rows FIRE `FAIL_DIMENSIONAL` (the dims that ENTER `P₀^phys`) while corrupt-`[G]` does NOT (`G∉free_symbols(P₀^phys)` — a dependency-scope diagnostic) | structural (dimensional-closure provenance, dual-engine + able-to-fail — dropping the frequency norm fires `3d`, corrupting a sourced port dim fires `3d′`, wiring the μ̂₀ back-solve into the verdict fires the read-set/demotion tooth; the v1 back-solve is a TAUTOLOGY [re-solving `[μ̂₀]` keeps `homogeneity_pass=True` under EVERY corruption → fires on nothing], correctly demoted — the μ̂₀-free gate's `[N₀]/[D₀]/[c_s]` FAIL rows are what reject it; **NOT a reduction** — a dimensional-consistency closure / the natural-units-trap catch, discharges NO knob; like R37/R38/R39 for 018/019/020's fingerprint/prefactor/provenance) | II-G4d (021) | the COMPLETING leg — **021 COMPLETES the joint `QUAD_CALIBRATED`** (018∧019∧020∧021), which STAYS CALIBRATED not PASS (020's provenance + `G=GENUINE_BLOCKED`). ZERO new counted knobs: `μ̂₀` is a free-carrier NON-verdict diagnostic (NOT a counted knob); the SOURCED port dims `[N₀]/[D₀]` are dimensional PROVENANCE (pathA_43 density-port numerator + the carried conservative denominator, Cluster-C 024/027 + pathA_34 not yet built), NOT new knobs; `c`/`G` already registered. 021 IS units-bearing AND does the `[·]` dim-homogeneity gate (the OTHER half of the 020/021 operation-level cut — 020 did algebra+provenance, 021 does dimensions). Consumption is PROVENANCE (the sourced port dims genuinely ENTER the gate so the corrupt-`[N₀]` tooth is genuine; 018's `u₂/u₄/v₅` a local frozen fixture for the `Yhat` dimensionless check; 019's `P0=N0/D0` enters `P0_raw`; 020's `target_rhs` the μ̂₀ diagnostic's rhs only — NO dual-site); the `.wl` KEEPS its native dim block (like 018/019, unlike 020) |
| R41 | cross-ℓ outgoing-DtN `−(ℓ+1)/Λ_ℓ` fingerprint provenance + the Gate-4 quadrupole non-regression: for ℓ=0,1,2 the outgoing DtN response `Ŷ_ℓᵒᵘᵗ=−(ℓ+1)/Λ_ℓ`, `Λ_ℓ=z·h_ℓ⁽¹⁾′/h_ℓ⁽¹⁾`, `z=aω/c_s`, `h_ℓ⁽¹⁾=j_ℓ+i·y_ℓ`, series-expands to the DERIVED radiative fingerprint `{ℓ=0: 1, ℓ=1: 1/2, ℓ=2: 1/27}` at leading orders `2ℓ+1={ω¹,ω³,ω⁵}`, with the static-DtN slot `Λ_ℓ(z→0)=−(ℓ+1)` DERIVED from the Hankel log-derivative (`lam_series.coeff(z,0)`, NOT the hand-set numerator) + the radiative order verified by SCANNING the imaginary series for its first nonzero power (all lower imaginary coeffs vanish); the incoming `h_ℓ⁽²⁾=j_ℓ−i·y_ℓ` flips ONLY the radiative sign; the ℓ=2 leg reproduces the completed pathA_33 quadrupole fingerprint `{u₂=1/9, u₄=4/81, v₅=1/27}` (the Gate-4 non-regression) | structural (fingerprint + non-regression provenance, dual-engine + able-to-fail — a pole-order corruption `h_mut=z·h` shifts `Λ_static: −(ℓ+1)→−ℓ` and fires (NOT the INERT outgoing→incoming branch swap, which leaves `Λ_static` unchanged — a Codex confirm-pass reasoning catch), a `+i·z` corruption fires the SCAN (the source's preselected-`2ℓ+1`-nonzero check does NOT), per-slot u₂/u₄/v₅ derivation-copy mutants each fire their own assert, the incoming branch flips the sign, `3e_break_gate4`→`FAIL_QUAD_REGRESSION`; **NOT a reduction** — earned exterior-wave structure + a consistency non-regression, discharges NO knob; like R37 for 018's ℓ=2-only fingerprint) | II-G5a (022) | the EARNED-FIRST leg of the pathA_34 2-way split — 022 lands the joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` as an EARNED-first PARTIAL (`LOCAL_AUDIT_VERDICT=CROSS_L_FINGERPRINT_OK`, exit 0; the `JOINT_LANDING_LABEL` a printed 1/2 partial). ZERO new counted knobs: `c_s` is a **units carrier** (R1-`DERIVED`, cited PROVENANCE), NOT a consumed value (the earned rationals are `c_s`-free); `a` is the `CONV` pin; ⚠ the ℓ=0/1 stiffnesses `K0c`/`K_eta+2·T_Omega` and the return admittances `Z0_ret`/`Z1_ret` are **023's** (the transfers/nullspace), NOT counted here. `static=1` + `χ_Q=1` are de-counted printed diagnostics (subsumed by `Λ_static` / by `v₅`: `27v₅−1≡27(v₅−1/27)`). Consumption: stage018's `{u₂,u₄,v₅}` is the one CHECKABLE derive-vs-typed non-regression (its typed side is 018's independently-earned literal — NO subsumed X≡X); 019/020/021 + the completed joint + 008's amplitudes + 009/010's bulk mode = PROVENANCE (NO dual-site). 022 is z-space only (NO units-restored dim leg — the `z=aω/c_s` realization + the residual dim gate are 023's). The `.wl` is RE-AUTHORED independent via built-in `SphericalHankelH1`/`H2`+`SeriesCoefficient` (the source `.wl`'s transliterated `branchData` discarded — the mirror-policy screen; UNLIKE 018/019/021's keep-native). ⚠ the FAIL headline (the native nullspace underdetermination + the `Z0_ret/Z1_ret` Gate-6 selector) is 023's — 022 owns NO nullspace and does NOT back-solve `ε_eff`/`Z` (the `FAIL_TAUTOLOGICAL` firewall preserved for 023) |
| R42 | the cross-ℓ native nullspace underdetermination departure + the sharpened Gate-6 return-selector obligation: over 11 genuine generator dofs `[OmegaU,OmegaW,Rmix,gU,gW,D0,K0c,K_eta,T_Omega,Z0_ret,Z1_ret]` the collected Gate-5 named constraints `{P0_raw, K0c, K_eta+2·T_Omega}` have constraint-Jacobian rank **3** (a genuine `sp.Matrix(rows).rank()` / constructive `NullSpace` — NO zero-padding) → native nullspace dim **8** → return-augmented rank **5** → **return-moving nullity 2**: the return admittances `{Z0_ret, Z1_ret}` survive every collected constraint (`Z_is_premise=True`, pathA_29) → `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`. The counterfactual selector `{Z0_ret=K0c, Z1_ret=K_eta+2·T_Omega}` collapses the return-moving nullity 2→0 (native nullity 8→6, rank 3→5) → `CROSS_L_RESIDUAL_PREDICTION` — a rank-collapse WITNESS that the gate is able-to-fail, NOT a proven selector | structural + `PENDING` obligation (dual-engine + able-to-fail — hardcoding the rank / zero-padding the Jacobian / faking the `.wl` `NullSpace` basis all make the audit FAIL; the selector flip is driven by the COMPUTED return-moving nullity; the `A0/A1` forward-build consumes 022's `{1,1/2}` checked vs the INDEPENDENT pathA_29 form; the `emit_epsilon_magnitude_as_derived` mutation fires `FAIL_TAUTOLOGICAL`; **NOT a reduction, NOT new free dofs** — it SHARPENS the existing `ε0/ε1` R24-family debt from "the transmissions are free" to "the linear return law leaves exactly a 2-dim return-admittance nullspace") | II-G5b (023) | the COMPLETING, FAIL-delivering leg — **DELIVERS the joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** (2/2; 022 landed the EARNED-first 1/2) and **COMPLETES the pathA_34 fold (022∧023) + closes the Gate-1–5 gravity ladder**. **Honest landing: the linear theory CANNOT pin the ℓ=0/1 return; Gate-6's nonlinear closure (sim-deferred) must supply 2 independent return equations (or an equivalent return law).** ZERO new counted **CALIB** knobs (CALIB set stays 6): `Z0_ret/Z1_ret` = coordinate ALIASES of the existing `ε0/ε1` FREE-UNREDUCED debt (row above; ZERO new free dofs, no double-count). ⚠ `K0c` + the ℓ=1 stiffness sector `{K_eta, T_Omega}` (via `K1`) = pathA_34-convention effective ℓ=0/1 stiffnesses, `FREE-UNREDUCED` **PENDING scalar-reduction — COUNTED as reduction debt** (row above; NOT `DERIVED`, NOT new `CALIB` — dims `M T⁻²`≠013/017, the stage016 convention trap; the reduction to the wall packet is PENDING). `q_free/eta_null/gain0/gain1` = control-construction symbols tracked-not-counted. Consumption: stage022's `{1,1/2}` (CHECKABLE derive-vs-typed → `A0/A1`) + the pathA_29 residual form (009/010) + `Z_is_premise` + 008's `M0/D1/R0/R1` + 017's `P0_raw` + 013/017 (context for `K0c/K1`) + `c_s` (R1)/`a` (CONV) as PROVENANCE. The `.wl` is RE-AUTHORED independent via a constructive `NullSpace` route (`Length[NullSpace[Jbase]]=8`, `MatrixRank[basis·Gᵀ]=2`, `Greturn·NullSpace[Jselector]=0`) — a materially different algorithm from the `.py`'s `augRank−rank0` (mirror-policy screen; like 020/022, UNLIKE 018/019/021's keep-native) |

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
- **⭐ Stage 015 (pathA_31 II-G2c): ONE new counted CALIB knob `Vp0/ℓ_c` — the SECOND calibration-adding Part-II stage;
  COMPLETES `BREATHING_CALIBRATED`.** The legacy-Hessian structure recovery + the Hellmann–Feynman force (the RHS 013
  deferred). The counted knob is the **breathing driving-force scale `Vp0/ℓ_c`** (the wall confinement potential slope
  driving the breathing; `Vp0` + the now-live `ℓ_c` are its manifestations — they appear ONLY as the ratio; `ρ*` is
  CONSUMED, frozen density stage005/011) — reduction debt R33 (sibling of R30's wall-response debt; the same deferred
  nonlinear-throat interior would derive it). The legacy `E_geom` constants `{κ,χ,σ_a,σ_L}` are the legacy-Hessian PATTERN
  basis (a re-parameterization of the calibrated `{μ_η,T_w,β}` closure), **NOT counted afresh** (would double-count the
  closure) — edge R32 records the structural recovery (the operator-projected `M_AB`/`K_AB` reproduce `H_legacy`'s
  structural signature; `full_matrix_fit=False`, a form match not an entry fit). ⭐ The pathA_31 v1 REJECTION's OTHER locus
  (the HF `x−x` two-route tautology) is ABSENT — the two routes are genuinely-different constructions (distributed
  projection vs Hellmann–Feynman parametric derivative) with `unsimplified_routes_identical` COMPUTED False (adversarial
  per-tooth ablation confirmed genuine; 1 vacuous HF-guard leg + 3 literal/scaffolding flags were found + remediated
  make-genuine/de-count). `c_S` NOT consumed. **⇒ Part-II counted CALIB set is now `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015)
  = 4** (with `K_η=T_wβ²` R29-derived; the legacy constants + `Vp0`/`ℓ_c` individually + `r_AL` tracked-not-counted).
- **⭐ Stage 016 (pathA_32 II-G3a): ZERO new counted knobs — the EARNED ℓ=2 SO(3) covariance theorem (the EARNED-FIRST
  leg of a 2-way split; component 1/2 of the joint `ISOTROPY_CALIBRATED`).** The five real ℓ=2 harmonics form one
  orthonormal SO(3)-irrep (`Gram=I₅`) with a single COMPUTED `−Δ_S²` eigenvalue `λ_m=ℓ(ℓ+1)=6` for every m (genuine
  Laplace–Beltrami + Rayleigh quotient + eigenfunction residual), and the K₂ **angular** stiffness `K₂=K̃+λ_m·T̃_Ω` uses
  that COMPUTED eigenvalue (the `forced_eigenvalue_probe` rejects a typed coefficient → `FAIL_NOT_COVARIANT`). The
  covariance theorem is **pure angular math** — it introduces **NO calibration** (like 011/012/014, a structural-edge
  stage). New tracked-not-counted symbols FIRST-APPEARING here but with **counting DEFERRED to 017's calibration partition**:
  `T_Ω`/`T̃_Ω` (the angular-stiffness density used in K₂ + the dim check) and `β₂(w)` (the frozen ℓ=2 radial profile) — 016
  uses them but does NOT count them (⚠ do not pre-empt 017's partition). ⚠ **The consumption is PROVENANCE + self-contained
  dimensional integrity, NOT a cross-stage relation:** pathA_32 uses a VOLUME-density / dimensionless-`β₂` convention on
  the measure `a²·dw·dΩ` (`[μ_η]=M L⁻³`, `[T_w]=M L⁻¹ T⁻²`, `[K_η]=M L⁻³ T⁻²`, `[T_Ω]=M L⁻³ T⁻²`) DIFFERENT from stage
  013's line-density / `β=L⁻¹` convention on `4π∫dw` (where `K_η=T_w β²` holds) — related by an `∫a²dΩ`≈L² bridge, NOT
  equal (a Grok compute-catch). So `μ_η`/`T_w` are cited as the same PHYSICAL wall constants counted at 013 (provenance,
  not a checkable relation), and the genuine able-to-fail integrity is pathA_32's OWN `[M₂]=M`/`[K₂]=M T⁻²` consistency
  with a one-dimension corruption of ANY sourced density (`[μ_η]/[T_w]/[K_η]/[T_Ω]`) firing `FAIL_DIMENSIONAL`. `c_S`
  NOT consumed (angular structure speed-free). New structural edge **R34** (the ℓ=2 SO(3) covariance provenance: five
  real harmonics = one SO(3)-irrep with `λ_m=6`; the K₂ angular stiffness uses the computed eigenvalue — discharges
  NOTHING, earned angular structure not a reduction). ⭐ The pathA_32 v1 trip-ups are ABSENT: the eigenvalues are COMPUTED
  (residual + `tautology_hash_collision` probe), the aggregate probe battery is intact (neuter-one flips), the sourced-`T_Ω`
  dim-mutation fires `FAIL_DIMENSIONAL`; the source's vacuous `k_coeff_equal` `λ−λ≡0` self-compare is DE-COUNTED (the K₂-
  coefficient computed-ness rides on a residual-on-the-assembled-K₂-coefficient + live `build_K2(lambdas)` binding + the
  bare forced probe). Adversarial per-tooth ablation confirmed genuine; 1 vacuous stamped-literal (`participates_in_verdict`)
  was found + remediated to a computed verdict-propagation check (fresh-agent `REVERIFY_CLEAN`). **⇒ Part-II counted CALIB
  set UNCHANGED at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) = 4** (016 adds zero; `T_Ω`/`β₂` deferred to 017).
- **⭐ Stage 017 (pathA_32 II-G3b): the THIRD calibration-adding Part-II stage — adds 2 new counted CALIB knobs `{T_Ω, β₂}`
  (+ the ℓ=2 port-kernel support/Maxwell scalars `{B̃, Z̃}` tracked/downstream-pinned); COMPLETES the joint `ISOTROPY_CALIBRATED`.**
  The grouped-P2 lane isotropy: the grouped `{20,21,22}` ℓ=2 lanes (assembled from 016's CONSUMED covariance `λ_m=6`, `K₂=K̃+λ_m·T̃_Ω`,
  Gram-diagonal `c_self=1`, via a GENUINE cross-stage dual-site — same pathA_32 convention, so a one-site `λ`/K₂-form corruption is
  CHECKABLE + fires, the coordinated escape closed by a single-`Y20` `(−Δ_S²)Y=λY` echo) respond ISOTROPICALLY: the P₂-projection
  **raw-D defects = 0 (PRIMARY)** and the normalized-u response defects = 0 (CROSS-CHECK; a pure-prefactor anisotropy MOVES raw-D but
  leaves normalized-u zero — raw-D is decisive). The response is genuinely able-to-fail (the six-probe battery {pure_prefactor,
  sector_selective, m_dependent, degenerate_beta_zero, singular_denominator, static_drop_inertia}; neuter-one flips the aggregate). The
  joint `verdict_from_gates` lands **`ISOTROPY_CALIBRATED`** (COMPLETE = 016 ∧ 017). ⭐ **The 2 new counted CALIB knobs:** `T_Ω`
  (the ℓ=2 ANGULAR/rotational stiffness density — a genuinely new physical DOF absent from the ℓ=0 breathing sector) + `β₂(w)` (the FROZEN
  ℓ=2 radial profile) — both **calibrated NOT derived** from the Gate-1 `R0` support equation (edge R36, the `ISOTROPY_PASS` target; the
  primary CALIBRATED-not-PASS reason, report :5). The radial scalars `M̃/K̃/T̃_Ω` are DERIVED manifestations `= ∫density·β₂²` (edge R35,
  the ℓ=2 analogue of R29 — NOT counted afresh). `μ_η/T_w` cited PROVENANCE (counted at 013; `K_η=T_wβ²` R29, ⚠ non-transferable across
  the volume-vs-line convention). The ℓ=2 port-kernel support/Maxwell scalars `{B̃0/2/4, Z̃0/2/4}` are frozen `CALIB` inputs FIRST
  load-bearing here but **tracked, NOT in the clean `{T_Ω,β₂}` headline** — the isotropy VERDICT is value-INDEPENDENT of them (structural
  port-kernel inputs, not tuned-for-isotropy) and the `Z̃` Maxwell couplings pin DOWNSTREAM (018–024 + the EM knit/Gate-4 `54/5`); their
  irreducibility is a Part-VII adjudication. `c_S` NOT consumed. New edges: **R35** (ℓ=2 radial-scalar reduction, `DERIVED`) + **R36**
  (the `R0`-support-equation debt, `PENDING` — the `ISOTROPY_PASS` target); **R34** backfilled into the edges table (016's covariance
  provenance, was rollup-prose only). ⭐ Tri-review CLEAN (`FIDELITY_CLEAN` + `ADVERSARIAL_ISSUES` → the `.wl` dual-site K2-form tooth was
  subsumed [reconstructed K₂ from `lambdaByChannel` instead of reading the assembled lane] → Codex rewired it to read the assembled K₂,
  mirroring SymPy → fresh-agent `REVERIFY_CLEAN` coupling meta-test). Dual-engine SymPy 118 / Mathematica 127, CWD-independent. **⇒ Part-II
  counted CALIB set is now `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6** (+ the ℓ=2 port-kernel `{B̃, Z̃}` support
  scalars tracked/downstream-pinned; `M̃/K̃/T̃_Ω` R35-manifestations; `K_η` R29-derived). ⭐ **pathA_32 fold COMPLETE (stages 016+017);
  the joint `ISOTROPY_CALIBRATED` EARNED+CALIBRATED.**
- **⭐ Stage 018 (pathA_33 II-G4a): ZERO new counted knobs — the outgoing ℓ=2 DtN Hankel-fingerprint slice (the EARNED-FIRST leg of a
  4-way split; component 1/4 of the joint `QUAD_CALIBRATED`).** The outgoing ℓ=2 DtN response `Ŷ₂ᵒᵘᵗ=−3/Λ₂ᵒᵘᵗ`, `Λ₂ᵒᵘᵗ=z·h₂⁽¹⁾′/h₂⁽¹⁾`,
  `z=aω/c_s`, series-expands to the DERIVED fingerprint `u₂=a²/9c_s², u₄=4a⁴/81c_s⁴, v₅=a⁵/27c_s⁵` (the `27`=`1/v₅ᶻ` COMPUTED — the
  `derived_in_gate` `27` that stage 020's `54/5=2·27/5` partition rides), with the sign `χ_Q=+1` outgoing / `−1` incoming COMPUTED from
  `j₂±i·y₂` (a typed `χ_Q` would be a tautology) and the standing `j₂` contrast (`Λ_stand(0)=+2`, no radiating term). ⚠ **`c_s` (the density
  sound speed) FIRST becomes a LIVE object in the Part-II gravity radiative sector here** (a physical shift from 013–017's `kξ≪1`
  matter-mode deferral) — but only as the **units-restoring carrier** in the physical coefficients (via `z=aω/c_s`): the earned rationals +
  `χ_Q` are `c_s`-FREE, so 018 does NOT consume the VALUE of `c_s`; it is `c_s0` from R1 (`c_s²=5Kρ⁴/m`, stage005), cited PROVENANCE, NOT a
  new knob (and distinct from the frozen-wall `c_S` of 011–017). New structural edge **R37** (the outgoing-DtN ℓ=2 Hankel-fingerprint
  provenance — discharges NOTHING, like R34/R26). ⚠ **The consumption is PROVENANCE ONLY (no dual-site)** — 018 is SELF-CONTAINED exterior
  spherical-Hankel algebra (the fingerprint is built from explicit `j₂`/`y₂`; the literal port-kernel `N_n/D_n` consumption is stage 019's),
  so — UNLIKE 017's genuine cross-stage dual-site — there is no checkable cross-stage relation to guard (a guard on an unused object would be
  a vacuous tooth). The `χ_Q` reconciliation (`+1` vs pathA_22b's `≈0.712`) is flagged as a Part-VII item, NOT merged. **⇒ Part-II counted
  CALIB set UNCHANGED at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6** (018 adds zero; the prefactor algebra = 019,
  the `54/5` magnitude + `G=GENUINE_BLOCKED` = 020, the μ̂₀-free dim closure = 021).
- **⭐ Stage 019 (pathA_33 II-G4b): ZERO new counted knobs — the squared-denominator prefactor algebra (the SECOND leg of the
  4-way split; component 2/4 of the joint `QUAD_CALIBRATED`).** The squared-denominator prefactor `P(ω)=D₀·N(ω)/D^cons(ω)²`
  (`N=N₀+N₂ω²+N₄ω⁴`, `D^cons=D₀+D₂ω²+D₄ω⁴`) series-expands to `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`,
  `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` (SERIES-EXTRACTED off the actual series, NOT typed — the derive-then-check firewall),
  with the N/D self-check the sharp discriminator: plain `N/D` gives `P₂=(D₀N₂−D₂N₀)/D₀²`, missing the factor of 2 (the
  computed gap `D₂N₀/D₀²≠0`) → probe `3g` fires `FAIL_PREFACTOR_ALGEBRA`; the squared-denominator object is provably the
  right one. ⚠ **This is where 018's deferred literal `N_n/D_n` port-kernel consumption lands — at the PROVENANCE level, NOT
  a value-consumption dual-site:** the abstract `D₀,D₂,D₄` are 017's exported D-lanes and the `N₀,N₂,N₄` are
  `build_port_moments`' concrete port N-moments (deferred Gate-6 branch data, emitted-but-never-checked); the algebra is
  PORT-AGNOSTIC (holds for any nonzero `D0..N4`), so — like 018 — there is NO checkable cross-stage relation to guard (a
  guard on the unused/deferred moments would be a vacuous tooth). 019 is **units-FREE** (no `c_s`/`a`/`G`/`μ̂₀`, no
  dimensional leg — the `[P₀^phys]=1` closure is 021's; the units-free algebra is enforced by a runtime free-symbol guard
  asserting the earned expressions' symbols == `{ω, D₀,D₂,D₄,N₀,N₂,N₄}`). New structural edge **R38** (the
  squared-denominator prefactor-algebra provenance — discharges NOTHING, like R37/R34). **⇒ Part-II counted CALIB set
  UNCHANGED at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6** (019 adds zero; the `54/5` magnitude +
  `G=GENUINE_BLOCKED` = 020, the μ̂₀-free dim closure = 021).
- **⭐ Stage 020 (pathA_33 II-G4c): ZERO new counted knobs — the `54/5=2·27/5` provenance partition + the CALIBRATED verdict
  label (the THIRD leg of the 4-way split; component 3/4 of the joint `QUAD_CALIBRATED`; the leg that LANDS the CALIBRATED
  headline).** The assembled ℓ=2 quadrupole magnitude `54Gc_s⁵/(5a⁵c⁵)` decomposes as `54/5 = 2·27/5` — a SymPy-VERIFIED
  rational identity BOUND to `target_rhs`/`v5_slot` (⚠ NOT a typed string, NOT a bare-literal `Rational(54,5)−2·27/5`
  tautology — a Grok/Codex genuineness catch) — where the `27` is `derived_in_gate` (018's fingerprint `1/v₅ᶻ`, cited, NOT
  re-derived) and the `2/5`+`G` are `external_bridge_input` (the GR Burke–Thorne `2G/5c⁵`). A 4-way PROVENANCE partition
  (`classify_provenance` tag-dominance `deferred>external>derived>convention`) classifies the assembled `54/5` as
  `external_bridge_input` (external DOMINATES the mixed earned+external tags), so the source-faithful verdict (`both derived →
  QUAD_PASS; else → QUAD_CALIBRATED`) lands **`QUAD_CALIBRATED` not `QUAD_PASS`** — edge **R39**. ⚠ **The verdict is
  PROVENANCE-driven, NOT G→λG-invariance-driven** — a SEPARATE g-invariance diagnostic exposes the invariance-only TRAP
  (`54/5` is G-invariant yet calibrated; an invariance-only test would MISLABEL it as earned). **ZERO new counted knobs:**
  `G` = `GENUINE_BLOCKED`/`external_bridge_input` (the thing the PDE does NOT deliver — already registered, the II-002
  force-magnitude row; 020 CLASSIFIES it, does not introduce it); `c` (the GR/PN vacuum light speed in `2G/5c⁵` and
  `54Gc_s⁵/5a⁵c⁵`) = the light cone `c_γ` in its **GR-units-bridge role** (`P₀ ∝ c_s⁵/c⁵ = 1/λγ⁵`; a benchmark units carrier
  cited from the PN corpus, NOT a fresh counted knob — the same "cited PROVENANCE, not a consumed value" logic as 018's `c_s`);
  the `2/5` = GR Burke–Thorne bridge; the `27` = 018's derived fingerprint (cited). ⚠ **020 IS UNITS-BEARING** (`{c_s,a,c,G}`
  live in the ALGEBRA — a reversal from 019's units-free slice) **but does NO dimensional-homogeneity gate** (`[P₀^phys]=1`
  is 021's) — the 020/021 cut is by OPERATION (algebra+provenance = 020; dimension = 021), enforced by a free-symbol
  allowlist (no `μ̂₀` symbol; the legacy `gamma_quad_eff` μ̂₀ STRING DROPPED) + a structural no-dimensional-dependency cut
  (the 020-local verdict's signature/bytecode carries no dim parameter/dependency). **Consumption is PROVENANCE** — 018's χ_Q=+1/`27` genuinely enter 020's SELF-CONTAINED equivalence
  bridge (a corrupt cited value breaks `forward`), 019's `P0=N0/D0` is `Gamma5`-definitional-only (absent from the `ok`
  residual), 017's D-lanes — NO theatrical cross-stage dual-site. ⭐ Tri-review CLEAN (`FIDELITY_CLEAN` + `ADVERSARIAL_ISSUES`
  → 3 LOW-severity vacuous/subsumed teeth [a near-tautological `dimensional_ok`-independence `f(x)==f(x)` in both engines; a
  subsumed `P0=N0/D0`-disjoint firewall; a subsumed tag-mutation `≠EXTERNAL` weaker than its `==DERIVED` sibling] → Codex
  remediated all 3 make-genuine, no de-counts [the structural signature/bytecode dim-cut; the `{N0,D0}⊆Gamma5 ∧ ∉residual`
  before/after firewall run-before-form; the independently-classified baseline-external "from" side] → fresh-agent
  `REVERIFY_CLEAN` coupling meta-test). Dual-engine SymPy 74 / Mathematica 82, CWD-independent; the `.wl` is a GENUINELY
  AUTHORED independent route (the source `.wl` had ZERO 020 content — native `Exponent`/`Together`/`Cancel`/`FullSimplify` +
  `Rational`/`Simplify` + a rank-based `MaximalBy` `Association` classifier + `/. G->lambdaG G`, NOT a `.py` mirror). ⚠ Grok's
  compute-verify caught the verdict-rule INVERSION at the directive stage (the default must be `else → QUAD_CALIBRATED`, not
  the inverted `PASS-unless-both-external`; a shippable rig fixed pre-build via a REQUIRED MIXED control). **⇒ Part-II counted
  CALIB set UNCHANGED at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6** (020 adds zero; the μ̂₀-free
  `[P₀^phys]=1` dim closure = 021 the COMPLETING leg).
- **⭐ Stage 021 (pathA_33 II-G4d): ZERO new counted knobs — the μ̂₀-free `[P₀^phys]=1` dimensional closure (the FOURTH,
  COMPLETING leg of the 4-way split; component 4/4 of the joint `QUAD_CALIBRATED`; the leg that LANDS the joint COMPLETE).**
  The assembled ℓ=2 quadrupole prefactor `P₀^phys=(c_s/a)²·(N₀/D₀)` is dimensionless from the SOURCED port dims `[N₀]=L⁻¹M`
  (density/continuity-port numerator) and `[D₀]=L⁻¹T⁻²M` (the carried reduced static conservative denominator `D₀=K−B₀−Z₀`):
  `[P₀_raw]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`, `dimensional_ok=True` — a **μ̂₀-FREE** gate (μ̂₀∉`raw_symbol_dims`; the
  verdict read-set excludes `μ̂₀`/`μ_dim`/`homogeneity_pass`; `participates_in_verdict=False`). The natural-units trap (the
  handoff `P₀=N₀/D₀` silently drops `(c_s/a)²`) is CAUGHT (`3d` `FAIL_DIMENSIONAL`); the corrupt-`[N₀]`/`[D₀]`/`[c_s]` rows
  FIRE while corrupt-`[G]` does NOT (`G∉free_symbols(P₀^phys)`, a dependency-scope diagnostic) — edge **R40**. ⚠ **The v1
  back-solve is a TAUTOLOGY** (re-solving `[μ̂₀]` keeps `homogeneity_pass=True` under EVERY corruption → fires on nothing),
  correctly DEMOTED to a non-verdict diagnostic — the μ̂₀-free gate's `[N₀]/[D₀]/[c_s]` FAIL rows are what reject it (a Codex
  design-review catch: the originally-proposed corrupt-`[G]` anti-v1 framing was backwards). **ZERO new counted knobs:** `μ̂₀`
  is a free-carrier NON-verdict diagnostic (NOT counted); the SOURCED port dims `[N₀]/[D₀]` are dimensional PROVENANCE
  (pathA_43 density-port numerator + the carried conservative denominator; Cluster-C 024/027 + pathA_34 not yet built), NOT new
  knobs; `c`/`G` already registered. ⚠ **021 IS units-bearing AND does the `[·]` dim-homogeneity gate** — the OTHER half of
  the 020/021 operation-level cut (020 did algebra+provenance, 021 does dimensions), the ONLY pathA_33 leg with a dimensional
  gate. **Consumption is PROVENANCE** — the sourced port dims genuinely ENTER the gate (so the corrupt-`[N₀]` tooth is genuine,
  per emitted-vs-checked); 018's `u₂/u₄/v₅` a local frozen fixture for the `Yhat` dimensionless check; 019's `P0=N0/D0` enters
  `P0_raw`; 020's `target_rhs` the μ̂₀ diagnostic's rhs only — NO dual-site. ⭐ Tri-review CLEAN (`FIDELITY_CLEAN` +
  `ADVERSARIAL_ISSUES` → the two KEY teeth [anti-v1 read-set/wired-mutant, corrupt-`[G]` scope control] GENUINE; 5 LOW-severity
  stamped/subsumed teeth → Codex remediated 2 make-genuine [`rerun_gate_logic` derived from the actual re-run,
  `participates_in_verdict` derived from the computed read-set] + 3 de-count [the QUAD-landing tautology + two subsumed
  summaries, retained as labeled prints] → fresh-agent `REVERIFY_CLEAN` coupling meta-test). Dual-engine SymPy 42 / Mathematica
  50, CWD-independent; the `.wl` KEEPS its native dim block (like 018/019, unlike 020's fresh authoring). Directive review =
  Codex→Grok→Codex bookend (Codex 6 BLOCKING folded — incl. the backwards corrupt-`[G]` discriminator, the corrupt-`[N₀]` dim
  `(1,−1,0)` not `−[D₀]`, the self-containment fixture, the `[D₀]` provenance, the `homogeneityPass`-in-guard, the `Yhat`
  structured-catch → Grok compute-verify `DIRECTIVE_CLEAN` → final Codex confirm). **⇒ pathA_33 fold COMPLETE (018∧019∧020∧021);
  Part-II counted CALIB set UNCHANGED at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6** (021 adds zero;
  the joint `QUAD_CALIBRATED` COMPLETE stays CALIBRATED not PASS, `G=GENUINE_BLOCKED`).
- **⭐ Stage 022 (pathA_34 II-G5a): ZERO new counted knobs — the cross-ℓ `−(ℓ+1)/Λ_ℓ` fingerprints + the Gate-4 non-regression
  (the EARNED-FIRST leg of the pathA_34 2-way split; component 1/2 of the joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`).** The
  cross-ℓ generalization of 018's ℓ=2-only fingerprint: for ℓ=0,1,2 the outgoing DtN `Ŷ_ℓᵒᵘᵗ=−(ℓ+1)/Λ_ℓ` series-expands to the
  DERIVED radiative fingerprint `{1, 1/2, 1/27}` at orders `{ω¹, ω³, ω⁵}`, with `Λ_static=−(ℓ+1)` DERIVED from the Hankel
  log-derivative + the radiative order SCANNED (first nonzero imaginary power, all lower vanish); the ℓ=2 leg reproduces the
  completed pathA_33 quadrupole `{u₂=1/9, u₄=4/81, v₅=1/27}` (the Gate-4 non-regression) — edge **R41**. **022 lands the joint
  `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` as an EARNED-first PARTIAL** (`LOCAL_AUDIT_VERDICT=CROSS_L_FINGERPRINT_OK`, script exit
  0; the `JOINT_LANDING_LABEL` a printed 1/2 partial). ⚠ **The FAIL headline is 023's** — 022 owns NO nullspace and does NOT
  back-solve `ε_eff`/`Z` (the `FAIL_TAUTOLOGICAL` firewall preserved; the 022-local read-set = `{cross_l_fingerprints,
  ell2_non_regression}` provably excludes any nullspace/`base_verdict`). **ZERO new counted knobs:** `c_s` is a units carrier
  (R1-`DERIVED`, cited PROVENANCE; the earned rationals are `c_s`-free), `a` the `CONV` pin; ⚠ the ℓ=0/1 stiffnesses
  `K0c`/`K_eta+2·T_Omega` and the return admittances `Z0_ret`/`Z1_ret` are **023's** (the transfers/nullspace), NOT counted
  here; `static=1`+`χ_Q=1` are de-counted printed diagnostics (subsumed by `Λ_static` / by `v₅`). **Consumption:** stage018's
  `{u₂,u₄,v₅}` is the one CHECKABLE derive-vs-typed non-regression (its typed side = 018's independently-earned literal, NO
  subsumed X≡X); 019/020/021 + the completed joint + 008's amplitudes + 009/010's bulk mode = PROVENANCE (NO dual-site). 022 is
  z-space only (NO units-restored dim leg — 023's). ⚠ **The `.wl` is RE-AUTHORED independent** via built-in
  `SphericalHankelH1`/`H2` + `SeriesCoefficient` (the source `.wl`'s transliterated `branchData` discarded per the mirror-policy
  screen; UNLIKE 018/019/021's keep-native). Directive review = Codex→Grok→Codex bookend (Codex 5 BLOCKING folded — the `.wl`
  re-authoring, the de-rigged `Λ_static`+scanned-order, the ℓ=2 double-count + χ_Q subsumption + per-slot mutants, the
  Earned/Deferred subset framing + explicit LOCAL/JOINT output, the checkable-consumption = stage018 fingerprint; then a Codex
  confirm-pass caught a directive-level reasoning error — the `Λ_static` ablation "outgoing→incoming" is INERT, so the mutant is
  a pole-order `h_mut=z·h`, and `static=1` de-counted → Grok compute-verify `DIRECTIVE_CLEAN` → final Codex confirm). Tri-review
  CLEAN (`FIDELITY_CLEAN` — independent re-derivation of the cross-ℓ fingerprints, the pole-order-vs-inert `Λ_static` mutant, the
  scan-vs-preselect order, the χ_Q subsumption, the per-slot isolation, all COMPUTED + no 019/023 leak; + `ADVERSARIAL_ISSUES` —
  both confirm-pass catches HELD live-proven, all core physics teeth genuine, 2 LOW-severity redundancy teeth [F1 `3e`
  `rerun_gate_logic` constant-length, F2 subsumed read-set-exclusion] → Codex remediated F1 make-genuine [compare the two traced
  verdicts] + F2 honest de-count [the exact-equality tooth retains the guard] + 2 `.wl` nits [symbol typo, mutation-aware parity]
  → fresh-agent `REVERIFY_CLEAN` coupling meta-test, no regression). Dual-engine SymPy 56 / Mathematica 65 → **55 / 64** (net −1
  per engine from the honest F2 de-count), exit 0, CWD-independent. **⇒ Part-II counted CALIB set UNCHANGED at `{μ_η, T_w, β}`
  (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6** (022 adds zero).
- **⭐ Stage 023 (pathA_34 II-G5b): ZERO new counted CALIB knobs (CALIB set stays 6); the ℓ=0/1 effective stiffnesses `K0c/K1` are
  new FREE-UNREDUCED reduction-debt (counted, flagged PENDING); `Z0_ret/Z1_ret` add ZERO new free dofs (ε0/ε1 aliases) — the native
  nullspace underdetermination departure that DELIVERS the joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (2/2, COMPLETING; the
  FAIL-delivering leg of the pathA_34 2-way split).** Over 11 genuine generator dofs the collected Gate-5 named constraints `{P0_raw, K0c, K_eta+2·T_Omega}`
  have Jacobian rank **3** (a genuine `sp.Matrix(rows).rank()` / constructive `NullSpace` — NO zero-padding, the v1 REJECTION
  locus repaired) → native nullspace dim **8** → return-augmented rank **5** → **return-moving nullity 2**: the return admittances
  `{Z0_ret, Z1_ret}` survive every constraint (`Z_is_premise=True`) → the EARNED FAIL; the counterfactual selector `{Z0_ret=K0c,
  Z1_ret=K_eta+2·T_Omega}` collapses the return-moving nullity 2→0 (native nullity 8→6) → `CROSS_L_RESIDUAL_PREDICTION` (the
  able-to-fail witness — the default verdict is the predictive token, so the FAIL is EARNED not baked). The `A0/A1` residuals
  FORWARD-consume 022's `{1,1/2}`, checked vs the INDEPENDENT pathA_29 form (`ε_ℓ=Z_ℓ/K_ℓ` FORWARD, never back-solved; the
  `FAIL_TAUTOLOGICAL` firewall + the `emit_epsilon_magnitude_as_derived` tooth forbid the ε_eff/Z back-solve) — edge **R42**.
  **ZERO new counted CALIB knobs (CALIB set stays 6); `Z0_ret/Z1_ret` add ZERO new free dofs (ε0/ε1 aliases); `K0c/K1` add
  COUNTED reduction-debt dofs (flagged PENDING):** `Z0_ret/Z1_ret` = COORDINATE ALIASES of the existing `ε0/ε1` FREE-UNREDUCED
  debt (no double-count, no new dof); ⚠ `K0c` + the ℓ=1 stiffness sector `{K_eta, T_Omega}` (via `K1=K_eta+2·T_Omega`) =
  pathA_34-convention effective ℓ=0/1 stiffnesses, **`FREE-UNREDUCED` `PENDING` scalar-reduction — COUNTED as reduction debt**
  (per the register rule pending debt stays counted until DERIVED; NOT `DERIVED`, NOT new `CALIB` — dims `M T⁻²`≠013 `M L⁻¹T⁻²`/017
  `M L⁻³T⁻²`, the stage016 volume-vs-line convention trap; do NOT identify with the raw densities; the reduction to the wall
  packet is PENDING); `q_free/eta_null/gain0/gain1` = control-construction, tracked-not-counted.
  **Consumption:** stage022's `{1,1/2}` = the CHECKABLE derive-vs-typed feeding `A0/A1` (corrupt `v1`→`A1_form` fires); the
  pathA_29 residual form (009/010) = the second checkable relation; `Z_is_premise` + 008's `M0/D1/R0/R1` + 017's `P0_raw` +
  013/017 (context for `K0c/K1`) + `c_s`(R1)/`a`(CONV) = PROVENANCE. ⚠ The `.wl` is RE-AUTHORED independent via a constructive
  `NullSpace` route (`Length[NullSpace[Jbase]]=8`, `MatrixRank[basis·Gᵀ]=2`, `Greturn·NullSpace[Jselector]=0`), materially
  different from the `.py`'s `augRank−rank0`. Directive review = Codex→Grok→Codex bookend (⭐ the v1 REJECTION locus — the
  SHARPEST in Part II: Codex design-review **7 BLOCKING** folded [the `.wl` decisive RE-AUTHOR; the "selector collapses the
  RETURN-MOVING nullity not the nullity" math fix — native nullity 8→6 not 0 — + corrected isolated ablation teeth; the
  counterfactual-witness-not-proven-Gate-6-selector relabel; the provenance-cut fix (022's `{1,1/2}`=`cited_earned_input` +
  `assert_not_derive` rewired to the 023-derived forward T0/T1 map + `gate4_prefactor` tag dropped); the strengthened firewall
  (`emit_epsilon` tooth + de-counted `rerun_gate_logic` + fixed `able_to_fail_bad`); the **K0c/K1 PENDING-not-DERIVED** register
  correction; the **Z0_ret/Z1_ret aliases-not-new-dofs** register correction] → Codex confirm-pass 5 BLOCKING + 2 nits →
  final-confirm 2 BLOCKING + 1 nit (all consistency-sweep gaps) → **Grok-4.5 compute-verify `DIRECTIVE_CLEAN`** [independent
  SymPy confirmed rank 3/nullity 8/return-moving 2, the selector flip with native nullity 8→6, `A1−expected_A1=0` iff `v1=1/2`,
  the dims, the `K0c/K1` dim-conflict + `Z_ret` alias conventions; validated the `.wl` constructive route incl. `Greturn·Nsel=0`
  as a genuine identity; 1 honest-scope note folded — raw nullity 8 includes `K0c/K1` self-constraint bookkeeping, verdict rides
  return-moving 2] → closing Codex confirm). Tri-review CLEAN (`FIDELITY_CLEAN` — independent SymPy re-derivation of the rank
  audit, selector flip, `A0/A1`+`v1=1/2` consumption, dims, all faithful + genuine rank not zero-padded + no ε back-solve + a
  materially-different constructive `.wl`; + `ADVERSARIAL_CLEAN` — per-tooth ablation matrix, 15 mutations both engines,
  hardcoding-the-rank/zero-padding/faking-the-`.wl`-basis all FAIL, the 4 KEY anti-rig properties CONFIRMED, 4 non-blocking
  de-count nits → Codex remediated 2 make-genuine [witness-preservation recomputes each Jacobian-row dot product from the stored
  witness; neutralized-mutation uses a cache-distinct inert context + independence check] + 2 de-count [provenance-documentation
  → labeled PROVENANCE prints; the T/ε identity → a SELF-CONSISTENCY check] → fresh-agent `REVERIFY_CLEAN` coupling meta-test, no
  regression). Dual-engine SymPy 116 / Mathematica 123 → **111 / 117** (net −5/−6 per engine from the honest de-counts), exit 0,
  CWD-independent at repo root AND `/tmp`, zero file I/O. **⇒ Part-II counted CALIB set UNCHANGED = 6** (023 adds zero); **⭐
  COMPLETES the pathA_34 fold (022∧023) + closes the Gate-1–5 gravity ladder** — the only remaining gravity item (Gate 6, the
  `{Z0_ret, Z1_ret}` nonlinear return-selector) is sim-deferred (a Part-VII open-item, R42).
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
