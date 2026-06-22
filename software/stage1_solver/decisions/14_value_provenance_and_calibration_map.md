# Decision 14 — Value provenance + calibration map ("every value is derived or calibration-necessity")

**Date:** 2026-06-21
**Status:** LIVE (Codex+GLM reviewed; updated through Gate 4, 2026-06-22 — see §6). Built after the user pivot to calibrate-predict (decision-13 §0 (5n)) and the user requirement:
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
**UPDATE (β resolved, see §5):** `λγ = c_γ/c_s` is now confirmed a **7th genuine input** (a *dimensionless* speed-ratio parameter) —
the declared parent action does not tie the photon cone to the sound cone (`BETA_GENUINE_GAP`). So the genuine free-parameter set is
`{n, K, m_GNLS, ρ0, μ0, q_*, λγ}`.
**Precision (which speed is free):** `c_s=√(5Kρ0⁴/m)` is DERIVED from `{K,ρ0,m}` — NOT a separate free parameter. Only `c_γ` is
unfixed by the action. So `λγ=c_γ/c_s` is ONE degree of freedom (all the freedom lives in `c_γ`; `c_s` is already pinned); whether
written as the dimensionful `c_γ` or the dimensionless `λγ`, it is the SAME single input — `c_s` does NOT add a second knob.
**UPDATE (Gate 4 DONE, 2026-06-22 — the gravity source-map is now a PROVEN genuine gap):** the gravity ratio `g_mhat²/g_G` is
`GENUINE_BLOCKED` (not derivable target-blind from the current action — see §1c + §6). So the source-map scale `g_mhat` (↔ `m̂0`)
joins `λγ` as an **8th genuine gap**: `g_G` is calibrated on Newtonian `G` (a benchmark anchor), and `g_mhat` is then fixed by
**consuming the GR-quadrupole `54/5` as a calibration ANCHOR** (so `54/5` is ABSORBED, not predicted). Genuine free-parameter /
gap set is now `{n, K, m_GNLS, ρ0, μ0, q_*, λγ, g_mhat}` — with the caveat that `g_mhat` is fixed by spending the `54/5` anchor,
not freely floating. Filling this gap (→ recover `54/5` as a prediction) requires ADDING a source-map/mass-bridge postulate to
the action — a modeling choice, exactly parallel to `λγ`/β, NOT a computation.

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
| `χ_Q` (=S_port) | the defect's outgoing-radiation coupling ("antenna efficiency": near-field→outgoing-wave) | **COMPUTED ≈ 0.712** (Gate 3 DONE, 2026-06-22): exact spherical-Hankel `l=2` outgoing-DtN around the frozen finite-core branch; `χ_Q = C5_defect / C5_branch-geometry-reference`; converged in even-`nw` (jackknife-stable, `nr`-independent, flat tail verified to `nw=320`). Budget = **±0.0008 numerical** + a **separate ~8.5% one-sided Z_w-reference definitional systematic** (branch-geometry = physical zero vs flat-`Z_w`). The FORM (`ω⁵`, `l=2`, `1/27` Hankel fingerprint) is derived EXACTLY (target-blind, dual-engine). NOT `S_port=1` — that default was a convention; the value is now computed. | DONE (do NOT re-solve; precision floored at ~8.5% by the reference *choice*, not by numerics) | pathA_22b Gate-3 report §; ANSATZ_LEDGER §1#4 (b); pathA_22a:7-11,70-71,85-88 |
| `λγ = c_γ/c_s` | speed of light ÷ speed of sound of the medium | **blocked on `β_bulk_to_brane`** (the bulk-metric→acoustic-speed identification); `C_B/C_E=1` already removed one dial | determine if the EM↔matter coupling that sets `c_γ/c_s` is specified-but-unreduced (→compute) or a genuine model gap (→ then it is an INPUT, calibrate it) | pathA_22b:130-136; pathA_20b:68 |
| `g_G` (↔ Newton `G`) | strength of the gravity-analog coupling (how defect stress sources long-range attraction) | **GATE 4 DONE → the ratio `g_mhat²/g_G` is `GENUINE_BLOCKED` (not derivable target-blind from the current action).** `K_stress=χ_N·ρ_∞4` derives cleanly, but `α_J` (mass bridge) sits only in this stress lane and does NOT cancel. `g_G` itself is fixed by CALIBRATING to Newtonian `G` (a benchmark anchor) — legitimate. | NOT derivable from the current action; would need a derived mass bridge (`α_J`). Calibrate `g_G` on `G`. | pathA_21:44-48,104; pathA_22b Gate-4 § overall-verdict |
| `g_mhat` (↔ `m̂0`) | the brane source-map scale (dimensionless residual after `m̂0` carries units) | **GATE 4 DONE → `BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE`** (dual-reviewed: fidelity `FAITHFUL` + adversarial `GENUINE_BLOCKED`). `m̂0` appears in the canonical action ONLY as a target-facing normalization (all 22 `m̂` sites sit in `m̂²·(branch)=54Gc_s⁵/…`); there is NO target-blind transverse source-map kernel `K_source` to derive. So `g_mhat` is **a genuine model gap → effectively a calibration INPUT** (absorbed by the quadrupole given the EM anchor), NOT a deferred computation. | reclassified INPUT/gap; derive-later requires ADDING a source-map postulate to the action (a modeling choice, cf. `λγ`/β) | pathA_22b:7-10,56; pathA_22b Gate-4 § H-findings/residual-ledger |
| `Z_int` | Maxwell localization integral (→ coupling renorm `μ0_eff=μ0/Z_int`) | **blocked** (placeholder `Z(w)` with undocumented floor=0.8; ideal integral divergent). **Off the ξ critical path** | export a genuinely decaying `Z(w)` or get floor provenance + physical w-extent | pathA_22b:66-78 |
| `L/a (=37/20)`, `c0=3/4`, `ε_r=1/20`, `γ₀`, `R0`, `J` | throat/brane geometry + profile moments | blocked on the stationary-profile / moving-throat solve | the branch solve selects them | ANSATZ_LEDGER §1#1-3,§5; pathA_21:56,65 |
| `Θ_Q`, `α_J`, `W_eff/a`, `J-selector`, branch-kernel choices | the named sub-factors inside the `g_G`/`g_mhat` mass-bridge + source-map chain | **GATE 4 confirmed BLOCKED** — `α_J` (mass bridge) and the source-map kernel are NOT derivable from the current action (no action/boundary/Hamiltonian J→m map; `m̂0` only target-facing). These are the constituents of the now-proven-blocked gravity ratio. | not derivable from the current action; would need a mass-bridge + source-map postulate | pathA_22a:75-81,88-89; pathA_21:63-71,97-105; pathA_22b Gate-4 § |

**Known numerical issue (logged 2026-06-22, separate follow-up).** The shared Maxwell assembler's centered W-lane first-derivative operator (`patha_b2b_maxwell.py`) is EXACTLY rank-deficient for ODD `nw` — a checkerboard/odd-even null mode that admits a spurious solution contaminating any quantity extracted on odd-`nw` grids. Gate 3 sidestepped it by using only EVEN `nw` (the full-rank, well-posed discretization), which is correct (excluding a singular discretization, not hiding a mode). **Other gates that may use odd `nw` should be audited**, and the shared stencil eventually fixed (non-decoupling/staggered lane operator or a clean null-projection). Not fixed in Gate 3 (Gate-3-local even-`nw` scope only).

### 1d. BENCHMARKS (external targets we calibrate against — not model values)
`Newtonian G`; `GR quadrupole 54/5` (= `2G/5c⁵` chain, Peters 1964); `CODATA m_p/m_e=1836.15…`; `g−2`, `5PN`, ringdown freqs (held-out).
Provenance: ANSATZ_LEDGER §(c-lit); prereg §H.

---

## 2. The genuine free-parameter count (the falsifiability measure)

After classification, the model's genuine free inputs/gaps are now `{n, K, m_GNLS, ρ0, μ0, q_*, λγ, g_mhat}`. Of the four ξ
factors that were once "branch-determined / compute-don't-calibrate" (`χ_Q, λγ, g_G, g_mhat`), the resolution is now SETTLED:
- `P0`, `χ_Q ≈ 0.712`: **COMPUTED** (extraction + Gate 3). These were genuinely deferred computations, now done.
- `λγ`: **genuine model gap** (`BETA_GENUINE_GAP`) → calibration INPUT (the action doesn't fix the photon/sound cone ratio).
- `g_G`, `g_mhat`: **Gate 4 (2026-06-22) PROVED the ratio `g_mhat²/g_G` is NOT derivable** target-blind from the current action
  (`GENUINE_BLOCKED`, dual-reviewed) → `g_G` calibrated on Newtonian `G`; `g_mhat` fixed by consuming the `54/5` anchor.
So the honest picture: the two factors that *could* be computed (`P0`, `χ_Q`) WERE; the other two are genuine model gaps (`λγ` by
the β-analysis, the gravity source-map by Gate 4). The "compute-don't-calibrate" program ran to completion — what remains
calibrated is calibrated because the physics genuinely doesn't fix it, not because we deferred the work.

**The under-determination that gates calibrate-predict (counting argument).** The ξ verdict has FOUR branch-determined unknowns
(`χ_Q`, `λγ`, `g_G`, `g_mhat`) but the gravity anchors give only TWO equations: Newtonian `G` (pins `g_G`, given derived `c_s, a, m`)
and the GR quadrupole (pins the single product `χ_Q·g_mhat²·λγ⁵`, given `P0, g_G`). So the gravity anchors leave the system
**under-determined by ~2**: the quadrupole match is currently ABSORBED, not predicted, and held-outs cannot be predicted until the
excess is removed. To close it we must **derive ≥2 of {`χ_Q`, `λγ`, `g_mhat`/`g_G`}** OR **add independent anchor(s)** (the EM-sector
anchor of §3 is one). The cheapest derivable pair is `χ_Q` (outgoing-DtN, feasibility-scope first) + `λγ` (settle `β`): if both close,
the two gravity anchors then pin `g_G` and `g_mhat` and the system is well-posed; if either is a genuine gap, we need the extra anchor.
This counting is the live gate on running an honest calibrate-predict and decides the next derivation targets. **Update (2026-06-22): the count is now RESOLVED.** Of the two derivations §2 said we needed, only `χ_Q` landed (Gate 3, ≈0.712); `λγ` is a genuine gap (`BETA_GENUINE_GAP`) and the gravity ratio `g_mhat²/g_G` is now PROVEN not-derivable (Gate 4, `GENUINE_BLOCKED`). So we did NOT get the ≥2 derivations that would have closed the system by derivation alone → **we must add an independent anchor.** Concretely: `g_G` is calibrated on Newtonian `G`; that leaves the quadrupole relating `g_mhat²` and `λγ⁵` — under-determined by ONE → **the EM-sector anchor (which pins `λγ`) is now LOAD-BEARING**; it is what closes the count. With it the system is well-posed but the `54/5` quadrupole is ABSORBED (it fixes `g_mhat`), not predicted. The falsifiable predictive content therefore lives ENTIRELY in the held-out surplus (g−2, 5PN, ringdown, multi-defect), which share the now-derived `χ_Q` + `P0/D0` bundle + `c_s`. Next derivation targets: the EM-sector anchor, then the held-out surplus.

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

1. **Compute the branch-determined-NOW quantities instead of calibrating them.** ✅ DONE: `P0` extracted; `χ_Q ≈ 0.712` computed
   (Gate 3). These were genuine deferred computations and they landed.
2. **The blocked quantities turned out to be genuine model gaps, not deferred computations.** ✅ RESOLVED: `λγ` = `BETA_GENUINE_GAP`
   (the action doesn't fix `c_γ/c_s`); the gravity ratio `g_mhat²/g_G` = Gate-4 `GENUINE_BLOCKED` (no target-blind source-map kernel;
   `α_J` underived + doesn't cancel). Per this section's own rule, they are now declared calibration INPUTs/anchors — not silent
   knobs. (`β` was settled first per the §5 revision; the gravity ratio was then proven blocked.)
3. **Add an EM-side anchor — now LOAD-BEARING, not optional.** With the gravity ratio blocked, `g_G` is calibrated on `G` and the
   quadrupole alone leaves `{g_mhat, λγ}` under-determined by one. The EM-sector anchor (which pins `λγ`) is what closes the count.
   It still also anchors `μ0, q_*, Z_int` (which the gravity anchors never touch). This is the next concrete derivation/calibration target.
4. **Minimal calibration set (updated):** genuine INPUTS `{n, K, m_GNLS, ρ0}` (gravity sector) + `{μ0, q_*}` (EM) + the gaps
   `{λγ, g_mhat}`, calibrated against `{Newtonian G, GR quadrupole 54/5, one EM anchor}` — note the quadrupole `54/5` is now a
   CONSUMED ANCHOR (it fixes `g_mhat`), not a prediction. Everything computable (`c_s, a, P0, χ_Q, C_B/C_E`) is COMPUTED. The
   falsifiable predictions are the held-out surplus (g−2, 5PN, ringdown, multi-defect), riding the SHARED `P0/D0`-bundle + `c_s` + `χ_Q`.

**Bottom line for the user's question ("what are we actually calibrating?") — updated 2026-06-22:** ~6 genuine physical inputs
(`n, K, m_GNLS, ρ0, μ0, q_*`) PLUS two proven model gaps that must be calibrated because the physics genuinely doesn't fix them
(`λγ` = the photon/sound cone ratio; `g_mhat` = the gravity source-map scale, Gate 4). The "compute-don't-calibrate" program ran
to completion: everything that COULD be computed was (`P0, χ_Q, c_s, a, C_B/C_E`); the residue is calibrated by necessity, not by
deferral. The cost: the GR-quadrupole `54/5` is an absorbed anchor, not a prediction — the predictive payoff is the held-out surplus.

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

### β RESOLVED (2026-06-21) — `BETA_GENUINE_GAP` (adversarial-reviewed `VERDICT_SOUND`)
The β-status analysis + a clean adversarial review concluded: the declared parent action treats EM as an INDEPENDENT fundamental gauge
field on the flat bulk metric (`pde.tex:257-262, 316-321, 357-416`), minimally coupled to the GNLS medium; nothing identifies the
bulk-metric null speed with the acoustic speed `c_s` (`pde.tex:541-577`). So `c_γ/c_s` is a genuine *physical* (invariant) ratio the
action leaves OPEN → **`λγ` is a calibration INPUT, not a prediction and not a removable convention.** It is NOT something to "compute"
— there is no derivation to do; the physics genuinely doesn't fix it. The ONLY counter-construction (where `c_γ=c_s` is forced) is the
**superseded** legacy emergent-acoustic-photon paper `em_fields.tex:160,172,475-482` (the photon as a collective mode of the medium);
Path-A canonically rejects it as non-parent (`pathA_20:35`, `pathA_20b:68`). So filling this gap would require ADDING a physical
postulate (e.g. declaring EM emergent → `c_γ=c_s`), which is a *modeling choice*, not a computation. Recommendation accepted: harden
the report addendum + `pde.tex` so this (and the m̂ dimensional inconsistency) don't re-trap a future reader.

---

## 6. GATE 4 CLOSEOUT (2026-06-22) — the gravity ratio `g_mhat²/g_G` is `GENUINE_BLOCKED`

**What was attempted.** The action-level derivation of the verdict ratio `g_mhat²/g_G` (the only gravity combination the verdict
sees — GLM Trap 1), via the cancellation test Gate-0b deferred. Directive: `directives/pathA_22b_gate4_ratio_or_blocked.md` (design-
reviewed `SOUND-WITH-FIXES` → 12 fixes → confirm-pass `SOUND-AS-IS`). A user-mediated GLM viability consult pre-judged it likely
blocked (driven by `α_J`) and recommended a CONFIRMATORY run; the directive was reframed (v3) to put the `α_J`-doom hypothesis
front-and-center and demote the `W_eff` cancellation to necessary-but-not-sufficient.

**Result (dual-reviewed: fidelity `FAITHFUL` + adversarial `GENUINE_BLOCKED`).** Verdict `BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE`:
- `K_stress(w) = χ_N(w)·ρ_∞4(w)` (the `g_G` stress lane) DERIVES cleanly, target-blind, with `N_∞,3 = ∫K_stress dw` (provenance
  pde.tex:396-416, 496-565; pathA_21b/21c). Dimensions check (`L⁻⁴ → L⁻³`).
- `K_source(w)` (the `g_mhat` source-map lane) does **NOT exist** as an independent target-blind transverse kernel. The adversarial
  agent independently swept ALL 22 `m̂` occurrences in canonical `pde.tex` — every one is target-facing (sits inside
  `m̂²·(branch data) = 54Gc_s⁵/…`); the paper itself (pde.tex:1879-1886) labels the residual `m̂_rad` an unsolved *dynamical*
  normalization. So there is no source-map kernel to derive independently of the `54/5` target.
- `α_J` (mass bridge) sits ONLY in the `g_G` stress lane (via `G_cond ∝ 1/(α_J1·α_J2)`); it has no independent route into `g_mhat`,
  so it does NOT cancel in the ratio (corroborated by pathA_21's prior `MASS_BRIDGE_FORM_NOT_DERIVED`, not circularly manufactured).
- The cancellation routes (a)/(b) were correctly `NOT_RUN_UNDEFINED_K_SOURCE` (you cannot test cancellation of a nonexistent kernel);
  negative + mutated-kernel controls confirmed the comparator is genuine (dual-engine, 10/10 `.wl`).

**Why this is the honest, valuable outcome (not a failure).** It converts "we left `g_G`/`g_mhat` as knobs" into a *derived,
adversarially-verified* statement: **the toy model as written does NOT derive its own gravitational coupling — `m̂0` and `α_J` are
imported as target-matched normalizations.** This closes Gate-0b's "did you actually try the action-level derivation?" definitively.

**Consequence for the verdict (`P0·χ_Q·g_mhat²·λγ⁵/g_G = 54/5`).** `g_G` is calibrated on Newtonian `G`; `P0`, `χ_Q` are derived;
that leaves the quadrupole relating `g_mhat²` and `λγ⁵` — under-determined by one → the **EM-sector anchor (pins `λγ`) is now
load-bearing**. With it the system closes but `54/5` is ABSORBED (it fixes `g_mhat`), not predicted. **Falsifiable content = the
held-out surplus** (g−2, 5PN, ringdown, multi-defect), riding the shared derived `χ_Q` + `P0/D0` bundle + `c_s`.

**The one path that would recover `54/5` as a prediction** (NOT a computation, a modeling choice): ADD a principled source-map
provenance + mass-bridge (`α_J`) postulate to the parent action — exactly parallel to the `λγ`/β gap. Out of scope unless the user
elects it.

**Artifacts (uncommitted at closeout):** `src/stage1_solver/patha22b_gate4.py`, `tests/test_patha22b_gate4.py`,
`tools/pathA_22b_gate4_crosscheck.wl`, report § appended to `reports/pathA_22b_minimal_combination_xi.md:322-393`.
