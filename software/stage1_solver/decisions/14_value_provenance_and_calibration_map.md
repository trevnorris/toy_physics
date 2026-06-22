# Decision 14 — Value provenance + calibration map ("every value is derived or calibration-necessity")

**Date:** 2026-06-21
**Status:** DRAFT for Codex review. Built after the user pivot to calibrate-predict (decision-13 §0 (5n)) and the user requirement:
*every value must be a named physical quantity that is either DERIVED, or fixed by CALIBRATION NECESSITY (we find what it must be for
the equations to close, freeze it, and derive it later if we can) — nothing is a free "knob."* This is the `n=5, β=3` discipline
applied to the whole Path-A constant set.
**Inputs:** two provenance sweeps (the `ANSATZ_LEDGER` classification + the pathA_19–22b reports; and the observable→constant
dependency map) + the gates 0–2 results. Companion to decision-13.

---

## 0. The rule + the classes

Every load-bearing value falls in exactly one class:

- **DERIVED** — computed from more-fundamental inputs + the field equations. (No freedom; it is what the math says.)
- **INPUT (calibration-necessity, genuine)** — an irreducible model-DEFINING choice. The real free parameters. Fixed either by
  matching a known physical input or by what closure REQUIRES, then frozen. These are the `n=5 / β=3` analogs: legitimately
  calibrated, *not* fudge factors. A small finite set — the model's true degrees of freedom.
- **BRANCH-DETERMINED (computable, NOT genuinely free)** — fully determined by the solved field configuration + the inputs; the real
  PDE *selects* it. It only *looks* like a knob because we deferred the computation. Sub-split: **computable-NOW** (from the existing
  finite-core background) vs **blocked-on-⟨named piece⟩**. **These must be COMPUTED, not calibrated** — treating them as free inflates
  the parameter count and destroys falsifiability (the user's core point).
- **BENCHMARK (external target)** — a published/physical number we calibrate *against* (not a model value). E.g. Newtonian `G`,
  the GR quadrupole `54/5`, CODATA `m_p/m_e`.

This maps onto the existing `ANSATZ_LEDGER` scheme (`research/pde_ledger/redteam_adversarial/ANSATZ_LEDGER.md`): its **(a)** =
INPUT/convention, **(b) branch-determinable** = BRANCH-DETERMINED, **(c-lit/c-rec)** = BENCHMARK. (Caveat: the ledger catalogs the
*moving-throat PDE* program (stages 022–253); the ξ factors `P0, χ_Q, g_mhat, g_G, λγ, Z_int` are the *newer Path-A* program, whose
provenance lives in the `pathA_19–22b` reports. Concept reused; sources cited per row.)

---

## 1. The value ledger

The verdict identity (pathA_22a, dimensional-skeleton-verified):
`R_norm = 0 ⟺ P0·χ_Q·g_mhat²·λγ⁵/g_G = 54/5`, with `G=(a·c_s²/m_GNLS)·g_G`, `m̂0=(c_s/(a²√m_GNLS))·g_mhat`, `c_γ=λγ·c_s`.

### 1a. Genuine INPUTS (the model's true free parameters — calibration-necessity)
| value | physical meaning | how fixed | derive-later? | provenance |
|---|---|---|---|---|
| `n=5` | quintic EoS exponent `P(ρ)=Kρ⁵` | frozen parent-action choice (defines the medium); the `β=3`-style "what it must be" call | the EoS is the model definition — not owed | ANSATZ_LEDGER §4/§5#3 (moved (b)→(a)) |
| `K` | EoS coefficient | free medium coefficient (sets the pressure scale) | calibration input | ANSATZ_LEDGER §4; pathA_19 |
| `m_GNLS` | constituent boson mass | explicit action parameter | input | pathA_19:14,44 |
| `ρ0` | asymptotic density | branch/state datum | input | pathA_19:46 |
| `μ0` | Maxwell vacuum-permeability analog | EM-sector coupling | input (EM sector) | pathA_22b:124 |
| `q_*` | charge coupling | EM-sector coupling (`q_eff=q_*/√Z_int`) | input (EM sector) | pathA_20b:17 |
| `ħ` | quantum of action | **fixed unit/external action constant (underived), NOT a calibrated model knob** | excluded from the *calibration* count but is a real external constant | pathA_19:43; pathA_20:53-57,70 |
| unit conventions, `a`-normalization | conventions | gauge/convention | not owed | ANSATZ_LEDGER §4 (a) |

→ **The model has ~6 genuine *calibration* inputs** (`n, K, m_GNLS, ρ0, μ0, q_*`), plus `ħ` (a fixed external action constant, not a
calibrated knob) and conventions. Everything below is DERIVED or BRANCH-DETERMINED from these — i.e. NOT independent.

### 1b. DERIVED (computed from inputs + equations)
| value | physical meaning | derived from | provenance |
|---|---|---|---|
| `c_s = √(5Kρ0⁴/m_GNLS)` | sound (phonon) speed | the quintic EoS (`c_s²=(1/m)dP/dρ`) | pathA_20:9-17 |
| `a = ħ/(m_GNLS·c_s0)` | healing length / mouth radius | **conditional-DERIVED** on the 4-pin/core-scale identification; otherwise a branch collective mouth-radius geometry (`A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`) | pathA_19:32,36-39,49,76 |
| `C_B/C_E = 1` ⇒ `c_bulk = 1` | EM cone is isotropic in bulk units | explicit `F_MN F^MN` expansion on the flat metric `η=diag(-1,+1,+1,+1,+1)` (Gate 2) | pathA_22b:118-124 |

### 1c. BRANCH-DETERMINED (computable, NOT free — must be computed, not calibrated)
| value | physical meaning | computable-now / blocked-on | path to derivation | provenance |
|---|---|---|---|---|
| `P0 = (c_s/a)²·N0/D0` | static quadrupole response of the throat (raw `N0/D0` has dim `T²`; `(c_s/a)²` makes it dimensionless) | **already extracted** (target-blind) from the existing finite-core background — the one branch-determined value in hand | (done; do NOT re-solve a deep/empty throat) | pathA_22a:24-25,41-43; pathA_22b:14,20 |
| `χ_Q` (=S_port) | the defect's outgoing-radiation coupling ("antenna efficiency": near-field→outgoing-wave) | **underived outgoing-DtN residual** — `S_port=1` is a current convention/DEFAULT, NOT confirmed from the existing background; feasibility of an outgoing-DtN solve on the *frozen finite-core* branch is itself OPEN (`pde.tex:2845-2849` actual branch unresolved). Tunability channel until derived | requires an actual outgoing-compact-DtN derivation; SCOPE feasibility on the frozen branch first | ANSATZ_LEDGER §1#4 (b); pathA_22a:7-11,70-71,85-88; pathA_21:90 |
| `λγ = c_γ/c_s` | speed of light ÷ speed of sound of the medium | **blocked on `β_bulk_to_brane`** (the bulk-metric→acoustic-speed identification); `C_B/C_E=1` already removed one dial | determine if the EM↔matter coupling that sets `c_γ/c_s` is specified-but-unreduced (→compute) or a genuine model gap (→ then it is an INPUT, calibrate it) | pathA_22b:130-136; pathA_20b:68 |
| `g_G` (↔ Newton `G`) | strength of the gravity-analog coupling (how defect stress sources long-range attraction) | **blocked on `W_eff`** + the mass bridge (`N_∞,3, Θ_Q, α_J`) | the `W_eff` transverse-profile derivation (Gate-0b said no shortcut — neither cancellation route holds) | pathA_21:44-48,104 |
| `g_mhat` (↔ `m̂0`) | the brane source-map scale (dimensionless residual after `m̂0` carries units) | **blocked on `W_eff`** / action-level kernel | same `W_eff` derivation (Gate 4) | pathA_22b:7-10,56 |
| `Z_int` | Maxwell localization integral (→ coupling renorm `μ0_eff=μ0/Z_int`) | **blocked** (placeholder `Z(w)` with undocumented floor=0.8; ideal integral divergent). **Off the ξ critical path** | export a genuinely decaying `Z(w)` or get floor provenance + physical w-extent | pathA_22b:66-78 |
| `L/a (=37/20)`, `c0=3/4`, `ε_r=1/20`, `γ₀`, `R0`, `J` | throat/brane geometry + profile moments | blocked on the stationary-profile / moving-throat solve | the branch solve selects them | ANSATZ_LEDGER §1#1-3,§5; pathA_21:56,65 |
| `Θ_Q`, `α_J`, `W_eff/a`, `J-selector`, branch-kernel choices | the named sub-factors inside the `g_G`/`g_mhat` mass-bridge + source-map chain | **blocked** (same `W_eff`/profile derivation as `g_G`/`g_mhat`; these ARE its constituents) | fold into / derived alongside `g_G`,`g_mhat` (Gate 4) | pathA_22a:75-81,88-89; pathA_21:63-71,97-105 |

### 1d. BENCHMARKS (external targets we calibrate against — not model values)
`Newtonian G`; `GR quadrupole 54/5` (= `2G/5c⁵` chain, Peters 1964); `CODATA m_p/m_e=1836.15…`; `g−2`, `5PN`, ringdown freqs (held-out).
Provenance: ANSATZ_LEDGER §(c-lit); prereg §H.

---

## 2. The genuine free-parameter count (the falsifiability measure)

After classification, the model's **genuine free inputs are ~6** (`n, K, m_GNLS, ρ0, μ0, q_*`). Everything in the ξ combination is
either DERIVED (`c_s, a, P0, C_B/C_E`) or BRANCH-DETERMINED (`χ_Q, λγ, g_G, g_mhat`). **None of the ξ "knobs" is a genuine free
parameter** — each is a deferred computation. This is the whole point: the apparent ">3 tunable knobs" of pathA_22a's
`TUNABILITY_CHANNEL_PRESENT` is really "3 branch-determined quantities we haven't computed," not 3 independent dials. The honest move
is to COMPUTE them (or, where genuinely a model gap like `β`, promote to a calibrated INPUT and *say so*), never to leave them as
silent free knobs.

**The under-determination that gates calibrate-predict (counting argument).** The ξ verdict has FOUR branch-determined unknowns
(`χ_Q`, `λγ`, `g_G`, `g_mhat`) but the gravity anchors give only TWO equations: Newtonian `G` (pins `g_G`, given derived `c_s, a, m`)
and the GR quadrupole (pins the single product `χ_Q·g_mhat²·λγ⁵`, given `P0, g_G`). So the gravity anchors leave the system
**under-determined by ~2**: the quadrupole match is currently ABSORBED, not predicted, and held-outs cannot be predicted until the
excess is removed. To close it we must **derive ≥2 of {`χ_Q`, `λγ`, `g_mhat`/`g_G`}** OR **add independent anchor(s)** (the EM-sector
anchor of §3 is one). The cheapest derivable pair is `χ_Q` (outgoing-DtN, feasibility-scope first) + `λγ` (settle `β`): if both close,
the two gravity anchors then pin `g_G` and `g_mhat` and the system is well-posed; if either is a genuine gap, we need the extra anchor.
This counting is the live gate on running an honest calibrate-predict and decides the next derivation targets.

---

## 3. Observable → constant dependency map (does calibrate-predict have teeth?)

Rows = observables; cells = depends (power/combo) / no / `?`. Each row flagged WORKED-OUT vs ASPIRATIONAL.

| observable | c_s | λγ | g_G(↔G) | χ_Q | g_mhat | P0/D0 bundle | μ0/q_* | Z_int | status |
|---|---|---|---|---|---|---|---|---|---|
| Newtonian G (anchor) | c_s² | no | =def | no | no | no | no | no | ASPIRATIONAL (`NEWTON_G_FORM_NOT_DERIVED`) |
| GR quadrupole 54/5 (anchor) | c_s⁵ | λγ⁵ | /g_G | χ_Q | g_mhat² | **P0=N0/D0** | no | no | PARTLY WORKED (extraction map + numeric P0; χ_Q,λγ,g_mhat²/g_G open) |
| inter-defect force | via EoS | no | ↔G | no | no | far-field | residualized | no | **WORKED (partial):** bilinear structure + reduced-3D/bulk power + LEADING far-field matter-stress attractive sign. FULL sign still residual (quantum, `V_conf`, Maxwell profile pieces not evaluated); normalization = explicit calib knob |
| g−2 | c_s (Ω_Q,σ_Q) | ? | no (only via P0) | **χ_Q** | no | **Ξ_1=P_1/P_0** | no | no | ASPIRATIONAL (reduces to Ξ_1+χ_Q; branch tangent open) |
| 5PN | c_s⁵/a⁵ | λγ³ tail | via 54/5 | χ_Q | g_mhat² | **P0,P2,P4; D0=K−B0−Z0** | no | no | ASPIRATIONAL (target surfaces; branch open) |
| ringdown / QNM | c_s | ? | no | no | no | **B_n,Z_n moments; D0** | no | no | ASPIRATIONAL (least developed) |
| moving-throat | c_s⁵/a⁵ | λγ³ tail | via 54/5 | χ_Q | g_mhat² | **full bundle** | **μ0,q_*** | **Z_int** | PARTLY WORKED (carries EM sector explicitly) |

**Synthesis:**
1. **The sharing is real (calibrate-predict has content):** the `P0/D0` modal bundle, `c_s`, `χ_Q`, and `λγ` recur across the anchor,
   5PN, g−2, moving-throat, ringdown — read from the SAME background. Held-outs are NOT per-observable free knobs.
2. **The EM sector needs its OWN anchor:** `μ0, q_*, Z_int` appear only in EM-carrying observables (moving-throat, and any EM-side
   prediction), NEVER in the gravity anchors (`G`, `54/5`). They cannot be fixed by gravity calibration → they need an EM-side
   calibration anchor (e.g. the charge/fine-structure sector), or independent derivation. Cleanest "needs its own anchor" gap.
3. **Honest status:** WORKED = inter-defect force (bilinear structure + power + *leading* far-field sign; full sign + normalization
   still open), GR-quadrupole extraction map + numeric `P0`, `c_s`. ASPIRATIONAL = Newtonian `G` itself (not derived), g−2 / 5PN /
   ringdown / multi-defect surpluses (target surfaces, refutation-firewalled as not-yet-closed). With `m̂0²·S_port` pinned=1 the anchor
   is a ~9-order MISS and pathA_22a returned `TUNABILITY_CHANNEL_PRESENT`.

---

## 4. Implications → the calibrate-predict plan

1. **Compute the branch-determined-NOW quantities instead of calibrating them.** `P0` is done. `χ_Q` is the next: it is computable
   from the existing finite-core background via the linear outgoing-DtN solve (old Gate 3) — do it, for the right reason (a
   branch-determined quantity, not a free constant). This shrinks the apparent knob set honestly.
2. **For the blocked branch-determined quantities** (`g_G, g_mhat` need `W_eff`; `λγ` needs `β`): either invest in the computation, or
   — only where it is a *genuine model gap* — promote to a calibrated INPUT and declare the dependency (the `n=5` treatment: freeze by
   necessity, derive later). Decide `β`'s status (computable vs genuine gap) FIRST.
3. **Add an EM-side anchor.** Gravity anchors cannot pin `μ0, q_*, Z_int`; the program needs one EM/charge-sector calibration anchor or
   these stay underdetermined.
4. **Minimal calibration set:** the genuine INPUTS (`n, K, m_GNLS, ρ0` for the gravity sector; `μ0, q_*` for EM) calibrated against
   `{Newtonian G, GR quadrupole, one EM anchor}`; everything else COMPUTED; held-out surplus (g−2, 5PN, ringdown, multi-defect) =
   the falsifiable predictions. Predictive content lives in the SHARED `P0/D0`-bundle + `c_s` + `χ_Q` + `λγ`.

**Bottom line for the user's question ("what are we actually calibrating?"):** ~6 genuine physical inputs, plus — only where a real
model gap exists (currently just `β_bulk_to_brane`, pending a check) — a small number of frozen-by-necessity values. The rest of the
"knobs" are computations we deferred. (Next-step recommendation revised in §5 after the GLM review: settle `β` first, not `χ_Q`.)

---

## 5. GLM tertiary review — corrections + reframing (2026-06-21)

A self-contained GLM consult on the §2 under-determination/closure argument (`_scratch/glm_consult_underdetermination.md`). Net:
one real flaw fixed, one GLM recommendation refuted by our own Gate 0b, and a sharper reframing. Digest:

- **ADOPT — `g_G`, `g_mhat` are NOT independent (GLM "Trap 1", the key catch).** Both are outputs of the SAME `W_eff` computation, and
  the verdict contains them only as the ratio `g_mhat²/g_G`. So they are ONE computation source, not two free knobs. **Never calibrate
  `g_G` and `g_mhat` independently** — that absorbs non-existent freedom and makes the quadrupole match spurious (it would always
  "succeed"). Compute the ratio / `W_eff`, or calibrate only `g_G` (from `G`) and let the ratio determine `g_mhat`.
- **REJECT (GLM didn't know) — the "ratio is cheap via `Z(w)` cancellation".** GLM asserted the `Z(w)` kernel cancels in
  `g_mhat²/g_G` "by standard KK adjoint structure" ⇒ a cheap ratio ⇒ the quadrupole becomes a free prediction. This is EXACTLY what
  **Gate 0b already tested and found `DOES_NOT_CANCEL (NOT_ESTABLISHED)`** — our specific sources establish neither the shared-factored
  scalar nor the proportional-kernel route. Reconciliation: Gate 0b is absence-of-proof, so GLM's KK intuition is a useful HINT that a
  dedicated `W_eff`/KK derivation (Gate 4) might yet prove cancellation — but it is NOT a free lunch; the gravity ratio stays blocked on
  real `W_eff` work.
- **ALREADY RESOLVED (GLM "Trap 2") — m̂ dimensionality.** Gate 0a = `MHAT_DIMENSIONFUL_CONFIRMED`; `g_mhat` is a real
  branch-determined residual, not a pure field value to be eliminated.
- **ADOPT as guards — GLM Traps 3 & 4.** (3) the `χ_Q` outgoing-DtN solve must be able to return `χ_Q≠1` (else the "derivation" is
  vacuous — assuming the `=1` default); make it a discriminating, parameter-free solve. (4) classify `β` explicitly as CONVENTION vs
  PHYSICAL PREDICTION before claiming `λγ` is "derived" — if it is a unit convention, `λγ` is assumed, not predicted.
- **REFRAMING (the user's original insight, vindicated).** Because every verdict factor (`P0, χ_Q, g_G, g_mhat`) is BRANCH-DETERMINED
  = computable, the apparent "tunability/under-determination" is largely the artifact of leaving them UNCOMPUTED and treating them as
  knobs. Computed, the dimensionless quadrupole `54/5` is closer to a **(near-)parameter-free PREDICTION than a calibration anchor** —
  genuine calibration is for the DIMENSIONFUL scales (actual `G`, masses → fix `n,K,m,ρ0,μ0,q_*`). Caveats: `λγ/β` status is open, and
  the background's own branch parameters (`R0,J`) may carry residual freedom that propagates into the factors. So the honest task is
  COMPUTE-don't-calibrate the branch-determined factors; the "under-determination by ~2" is really "2 uncomputed factors," and the
  gravity pair bottlenecks on the one `W_eff` computation (Gate 0b: un-shortcut-able from current sources).
- **REVISED NEXT STEP.** Settle `β` first (cheapest; pure analysis; classify convention-vs-prediction per Trap 4) → then `χ_Q` via the
  outgoing-DtN (non-vacuous per Trap 3). `g_G`/`g_mhat` remain the hard `W_eff` blocker (Gate 4), with GLM's KK-cancellation intuition
  as the concrete thing to TEST there. Do NOT calibrate the gravity pair independently in the interim (Trap 1).
