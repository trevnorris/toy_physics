# Part IV — Charge: RATIFIED atomic-stage split (per-Part gate SATISFIED — building)

> **⭐ Authored under the "ledger = surviving solution only" standing rule (user, 2026-07-22).** Part IV builds the
> **surviving** charge sector = the static ±w throat (puncture-deflection). It does **NOT** fold the superseded
> predecessors (leftover-scalar `u_L`-clamp, defect/antidefect wall-healing, the old `pathA_38` throat-body Coulomb
> anchor) — those are in the failures-paper backlog (`ledger_exclusions_failures_paper_backlog.md`, Exclusion 2). See
> memory `[[feedback-ledger-surviving-solution-only]]`.
>
> **Status:** ✅ **RATIFIED — building** (2026-07-22, overnight). The per-Part gate is satisfied by the user's conditional
> authorization ("if you can do the bookend cleanly and don't run into any issues, you have my permission to build it").
> **Bookend CLEAN:** Codex design-review `NEEDS_FIXES(5)+nit` → folded → Codex confirm `STILL_NEEDS_FIXES(3)` → folded →
> Codex confirm2 `CLEAN` → Grok-4.5 compute-verify `CLEAN` (all 6 SymPy checks + source-script runs exit 0) → Codex final
> `GATE_READY`. **Locked choices** (my recommendations, bookend-validated): **3 stages** (030/031/032); build ids fresh
> `030/031/032`; **IV-1 scope = fold ONLY the charge-new electric-scalar predicates**, bulk/wall/drain predicates cited as
> Part-I/II provenance (de-dup deferred to Part VII). ⚠ Concrete-source confirmation: `g0_closure_card_v0_checks.py` is a
> monolithic all-7-gate checker → stage 030 folds the electric-scalar subset only. **Build progress tracked in the
> Progress log at the bottom of this file.**
>
> Governing: `notes/ledger_v2_blueprint.md` §2 (granularity), §3 Part IV row, §5 (reshape spec — **LIGHT tier**),
> §6 (verification protocol); `docs/model_map.md` §3.4 + §5.1; resume doc `RESUME_ROADMAP.md` §3 (EM re-scope).
> Source builds: `software/em_charge_attribute/{puncture_deflection_electric_sign_result.md, g0_closure_card_v0.md,
> reports/native_p_constraint_gate.md}`.

---

## The surviving charge story, in one breath

**Charge = the static ±w throat.** A particle is a defect that punctures the brane into `±w`; the puncture bends the
brane, and the bend is the geometric field `ξ_w = ℓ h` (`[ξ_w]=L`). That embedding scalar `h` is the one gapless zero
mode of a localized parent field `H`; a mouth functional sources it orientation-oddly (`Q_χ[r_Σ,s]=s`), so the **charge
sign = the ±w orientation** — Z₂-topological, NOT additive. Two throats interact through `h` with a far-field
`F = s₁s₂ A/4πR²` (the `h`-mediated `1/R²`). The **falloff structure and `s₁s₂` form are target-blind EARNED** (the
magnitude normalization is a separate **R1-deferred** open item — an unbounded core factor `c_a`,`c_ξ` at tier-A), but
the **sign is not invariant across the admissible mouth boundary-condition ensembles** (V repel / J attract / M indefinite),
and the committed bare model does not select one → the honest terminal landing is **`R1_REQUIRED(bc_selection)`**
(neither earned NOR calibrated; needs the sim-deferred nonlinear throat/core solve). Charge coexists with gravity on the
same brane (`internal_inconsistency = none`, computed). And exact U(1)/Maxwell Gauss is proven **non-native**
(`NATIVE_P_NO_EMERGENT_GAUSS`) — a first-class characterized departure: **EM is not exact Maxwell.**

---

## RECOMMENDED split = **3 stages** (dependency order: closure → mechanism+force+landing → departure)

Each stage folds one of the three already-verified `em_charge_attribute/` script pairs (LIGHT reshape — they are
already assert-heavy print/`raise`-style engines; folding = wrap into the 6-artifact stage form + confirm the `.wl` is a
genuinely independent route).

| Build id | Stage | Source build + scripts | Headline verdict token | Scope | Reshape / notes |
|---|---|---|---|---|---|
| **030** | **IV-1** Electric scalar block + localized-`H` closure | `g0_closure_card_v0.md` §§2.3–2.4 + `g0_closure_card_v0_checks.{py,wl}` | class-(1) static PASS bundle (`PASS_LOCALIZED_H_NORM`, `PASS_REDUCED_H_INERTIA`, `PASS_STABILITY`, `PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS`, `PASS_DIMENSIONAL_HOMOGENEITY`) | **EARNED** (static prerequisites of the charge substrate): the reduced `(u_L,h)` action `S_Lh`; the parent `H` Pöschl–Teller operator `O_⊥=A†A≥0` with **one** gapless localized zero mode `f₀=1/[ℓcosh²(w/ℓ)]`, `N₀=8/(3ℓ)`, `h=P₀H`; `M_h=N₀M₄`, `K_h=N₀K₄=M_h c_E²`; the positive coupled kernel `D=B_eff K_h−C_hu²>0` (`[D]=M²T⁻⁴`) with **dimensionless witness** `D_*=7/4` and `c_±²=(3±√2)/2>0`; dim firewall. **POSTULATED** (labeled): the whole G0-v0 closure is a DRAFT-v0 postulated closure. | LIGHT. ⚠ **Fold only the charge-NEW electric-scalar sector**; the shared bulk/wall/phase pieces (`ρ,θ,r_B`) DUPLICATE Part I → cite as provenance, **de-dup DEFERRED to Part VII** (reconcile, not additive-merge). IV-1 ends at the G0 kernel stability + witness `D_*=7/4`; the response matrix `m`/`m_gg`/`z_g`/`det m` are IV-2 (puncture-deflection) objects. |
| **031** | **IV-2** Puncture-deflection: field identity → `1/R²` force → `R1_REQUIRED(bc_selection)` | `puncture_deflection_electric_sign_result.md` (all Q-blocks + App A–E) + `puncture_deflection_electric_sign_check.{py,wl}` + `puncture_deflection_electric_sign_independent_recompute.py` | **`R1_REQUIRED(bc_selection)`** | **EARNED (target-blind structure):** field identity `ξ_w=ℓh`; live coupling reduction `−J_m Q_χ H → −g_χh Q_χ h` (since `f₀(0)=1/ℓ`); the orientation-odd mouth source `η_i(k_m h − g_χh s_i)`, `Q_χ[r_Σ,s_i]=s_i` ⇒ **charge sign = ±w Z₂ orientation**; exterior `h(r)=h_A a/r` (positive holding curvature `4π(D/B_eff)a`) ⇒ the `h`-mediated `1/R` potential / `1/R²` force with the `s₁s₂` product form. The completed-square response matrix `m` (from the puncture-deflection build) gives the pair coupling `m_gg=B_eff z_g²/D` with the Robin escape factor `z_g>0` and `det m=z_g²/D>0`. **R1 (sign):** the four admissible BC classes give **decided conditional** coefficients `A_V=m_ggφ²/S_gg²>0` (repel), `A_J=−m_gg(j+g)²<0` (attract), `A_M=m_gg(q²−g²)` (indefinite), `A_MIXED=m_gg[(1−2λ)q²−2λqg−g²]` (spans); the committed bare model (G0 + U2-unresolved `𝔅` + S_hold scoped to `r_B−½`) does not select the class ⇒ `outcome_not_invariant` ⇒ terminal **`R1_REQUIRED(bc_selection)`** — the sign is **neither EARNED nor CALIBRATED**. `internal_inconsistency = none`. Co-blocker: `R1_REQUIRED(magnitude)` (`c_a`,`c_ξ` core normalizations unbounded at tier-A). | LIGHT. The whole self-contained puncture-deflection gate (one script pair). Carry the sealed **23040-cell** exhaustive §4 landing table as the able-to-fail backbone (first-matching-predicate ladder → the production tuple lands `R1_REQUIRED(bc_selection)`). |
| **032** | **IV-3** `NATIVE_P_NO_EMERGENT_GAUSS` (the characterized departure) | `reports/native_p_constraint_gate.md` + `native_p_gate_sympy.py` + `native_p_gate_dual.wl` (+ `native_p_gate_compare.py`) | **`NATIVE_P_NO_EMERGENT_GAUSS`** | **CHARACTERIZED-DEPARTURE (first-class):** the genuine Dirac–Bergmann constraint-class gate — at quadratic order native-`P` theories A and C both have **FC=0** (every constraint second-class; PB rank 8 / 12, no first-class Gauss chain); the tuned rank-drop strata that carry FC directions have a **zero Gauss descendant** (`DESCENDANT_ZERO`, hardening guard PASS). Exact U(1)/Maxwell Gauss is **non-native** → "EM is NOT exact Maxwell." ⚠ **Scope (carry honestly, does not weaken the verdict):** the regular open kinetic stratum is *fully symbolic* (FC=0 for all retained couplings); the non-common tuned rank-drop locus is *argued + fixed-seed scanned* (6 representative points per theory, 12 total across A+C), NOT an exhaustive symbolic stratification — any hypothetical missed measure-zero Gauss stratum would be a TUNED/inverse-design artifact, so the physical no-go (native `P` does not *generically* host emergent EM) is decisive independently. | LIGHT. Six able-to-fail controls anchor it (Maxwell→`FIRST_CLASS_GAUSS`, gauged-hard-unit→`MIXED`, nonconserved-current→`INCONSISTENT_PRESERVATION`, Coulomb-gauge→`SECOND_CLASS_NO_LOCAL_GAUGE`, global-U(1)→`GLOBAL_CHARGE_NO_LOCAL_GAUSS`, bare-σ→`SECOND_CLASS_RADIAL_NO_GAUSS`) + per-tooth ablation (`FIRED_AT_OWN_ASSERT`). |

**Net RECOMMENDED Part IV = 3 stages (030–032).** Build ids continue the global build-order sequence (Part III's
surviving stage is `003`; the couple-stress no-go that once held 030/031 was DROPPED to the failures paper, so those
numbers are free — but to avoid any collision/confusion with the retired 030/031, **I propose new build ids 030/031/032
are used fresh here** and note the retired numbers were never committed as stages). *(Build-id assignment is a
mechanical detail confirmable at the gate; the Part-ordered linear renumber is deferred to the pre-swap manifest pass.)*

---

## ALTERNATIVE split = **4 stages** (the gate decision to make)

Split IV-2 into its EARNED-structure half and its R1-landing half — the ledger's established idiom for joint gates
(Part II split `pathA_30`→011 earned-reduction + 012 the BC-dependent landing; `pathA_29`→009 residual + 010
localization). Costs one extra reshape: Codex splits the single `puncture_deflection_electric_sign_check.{py,wl}` into
two standalone stage-script pairs (precedent exists — 011/012 came from one `pathA_30` source).

| Build id | Stage | Headline | Scope class |
|---|---|---|---|
| 030 | IV-1 Electric scalar block + localized-`H` closure | class-(1) static PASS bundle | EARNED / POSTULATED |
| 031 | IV-2 Field identity + `χ↔h` mouth source + exterior `1/R²` | `THROAT_H_SOURCE_1_OVER_R2` (proposed token) | **EARNED (target-blind structure)** |
| 032 | IV-3 Far-field force ensembles → the honest landing | `R1_REQUIRED(bc_selection)` | **R1 (sign)** |
| 033 | IV-4 `NATIVE_P_NO_EMERGENT_GAUSS` | departure | CHARACTERIZED-DEPARTURE |

**Trade-off.** 4-stage cleanly separates the target-blind EARNED far-field structure from the R1(sign) landing (they
are different scope classes; separately citable — goal 3) and matches the Part-II joint-gate idiom. 3-stage is simpler,
maps 1:1 to the three verified script pairs (least reshape surface = least rig risk), and honors "a self-contained gate
stays one stage" + the user's standing steer toward keeping the ledger clean/minimal. **I recommend 3-stage**; flagging
4-stage so the choice is explicit at the gate.

---

## What Part IV EXCLUDES (surviving-solution compliance)

Per the standing rule, the following are **NOT folded** (already in `ledger_exclusions_failures_paper_backlog.md`
Exclusion 2 — carried as one-line "superseded predecessor" pointers only):

- **leftover-scalar `u_L`-clamp** (`NO_NATIVE_CLAMP`) — `leftover_scalar_electric_sign_result.md`;
- **defect/antidefect wall-healing** — `docs/em_phaseC_force_decomposition.md`;
- **old `pathA_38` throat-body Coulomb anchor** (`THROAT_ELECTRIC_LOCALIZED_COULOMB`) —
  `software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md`.

⚠ **The parameter register's `Q_E` row currently reads "pending (pathA_38)"** — Part IV must **re-home** it to the
surviving puncture-deflection build (as `R1_REQUIRED(magnitude)` reduction-debt, NOT a clean CALIB anchor), and note
`pathA_38` → backlog. This is a register-hygiene edit the IV-2 register update carries.

---

## Parameters Part IV adds / re-homes (register preview — the per-stage update fills these)

Charge-sector knobs from the G0 card + puncture-deflection build (dims from `g0_closure_card_v0.md` tables, `[L,T,M]`
order). Full dual-engine-verified dims land in the per-stage `parameter_register.md` update (blueprint §5 step 7).

| Param | `[L,T,M]` dim | Enters | Class (proposed) | Note |
|---|---|---|---|---|
| `H` (parent electric scalar) | `L⁻¹` | IV-1 (030) | `ACTION` (postulated PT operator) | localized; `O_⊥=A†A`, one gapless zero mode |
| `h = P₀H` (reduced scalar) | `1` | IV-1/IV-2 | `DERIVED` | zero-mode projection `h=N₀⁻¹∫2f₀H` |
| `M₄`, `K₄` (parent couplings) | `M`, `M L² T⁻²` | IV-1 | `ACTION` (postulated primitives) | `K₄=M₄c_E²`; the reduced-normalization witness `M_h*=1` is a `CONV` choice ON `M₄`, not a reduction |
| `M_h = N₀M₄` (scalar inertia) | `M L⁻¹` | IV-1 | `DERIVED` (projection relation `=N₀M₄`) | `K_h=M_h c_E²`; the `M_h*=1` value is the `CONV` normalization on `M₄` (row above) |
| `K_h = N₀K₄` (scalar stiffness) | `M L T⁻²` | IV-1 | `DERIVED` | `= M_h c_E²` |
| `c_E` (electric speed) | `L T⁻¹` | IV-1 (first entry) | `FREE-UNREDUCED` | **already seeded** (VI/040-041); **no committed cone lock** (`c_s/c_E=2` chosen to avoid one) — Part VI re-adjudicates |
| `ℓ` (embedding/healing length) | `L` | IV-1/IV-2 | `CALIBRATED`/`CONV` | `ℓ/a=ε_r=1/20` frozen handoff scale |
| `ξ_w = ℓh` (embedding displacement) | `L` | IV-2 (031) | `DERIVED` | **the field identity** (the puncture's brane-bend) |
| `A_eff = ρ_br + C_J²/κ_phase` | `M L⁻³` | IV-1 | `DERIVED` structural relation (Schur elimination of `θ_B`); the witness partition `ρ_br/A_eff=3/4`, `(C_J²/κ_phase)/A_eff=1/4` is `[POSTULATE]` (G0 tag) | eliminated-phase inertia — the *relation* is derived, the *numerical split* is postulated (do NOT read as a reduction) |
| `B_eff = ρ_B0²/χ_c` | `M L⁻¹ T⁻²` | IV-1 (consumed) | `DERIVED` | **already registered** (stage003 R16) — consumed, not re-counted |
| `C_hu` (scalar mixing) | `M T⁻²` | IV-1 (030) — **first dual-engine dim verification** | `FREE-UNREDUCED`/`PENDING` | ⚠ **NOT consumed from 041** (041 not built); 030 is C_hu's first built dim check; reduction route = Part VI/pathA_41; witness `C_hu*=1/2` is an IMPOSED edge |
| `K_m` (mouth Robin impedance) | `M L⁴ T⁻²` (`E·L²`, G0 `(4,-2,1)`) | IV-2 (031) | `ACTION` | reduced `k_m=K_m/ℓ²`, `[k_m]=E=M L² T⁻²` |
| `J_m` (mouth source) | `M L³ T⁻²` (`E·L`, G0 `(3,-2,1)`) | IV-2 (031) | `ACTION` | reduced `g_χh=J_m/ℓ`, `[g_χh]=E=M L² T⁻²` (a **decided** committed coupling, NOT R1) |
| `N_χ` (orientation-projection norm) | `L⁻¹` | IV-2 (031) | `DERIVED`/`CONV` | `Q_χ[r_Σ,s]=s` (the Z₂ orientation → sign) |
| `D = B_eff K_h − C_hu²` (physical kernel det-block) | `M² T⁻⁴` | IV-1 (witness) → IV-2 (response) | `DERIVED` | dimensionless **witness** `D_*=7/4` (`c_±²=(3±√2)/2`); the positive-definiteness is IV-1, the response objects below are IV-2 |
| `z_g` (Robin escape factor) | `1` | IV-2 (031) | `DERIVED` | `z_g=1−k_m⟨η,L_h⁻¹η⟩>0`; G0 witness `z_g=1` |
| `m_gg = B_eff z_g²/D` (pair coupling) | `M⁻¹ L⁻¹ T²` | IV-2 (031) | `DERIVED` | the two-body response amplitude |
| `det m = z_g²/D` | `M⁻² T⁴` | IV-2 (031) | `DERIVED` | `>0` (positive response) |
| `A_V/A_M/A_J/A_MIXED` (force coefficients) | `[A]=EL=M L³ T⁻²` (result:178) | IV-2 (031) | `DERIVED` (**decided** conditional formulas) | computed given a class: `A_V=m_ggφ²/S_gg²`, `A_M=m_gg(q²−g²)`, `A_J=−m_gg(j+g)²`, `A_MIXED=m_gg[(1−2λ)q²−2λqg−g²]` — NOT R1; the R1 is only WHICH class holds |
| BC-class selection + ensemble data `{φ, q, j, λ}` | — | IV-2 (031) | **R1-deferred (`bc_selection`)** | the *selection* among V/M/J/MIXED is the R1 (`outcome_not_invariant`); `{φ,q,j,λ}` parameterize the unselected classes — tracked, not counted |
| `s_i = ±1` (Z₂ orientation charge) | `1` | IV-2 (031) | `CONV`/structural (topological, **not a tuned knob**) | charge sign; ±w orientation, not additive |
| `Q_E` / charge-magnitude | pending | IV-2 (031) re-home | **R1-deferred (`magnitude`)** | **re-home from pathA_38** → puncture-deflection; `c_a`,`c_ξ` unbounded at tier-A |

New reduction-debt edges the register gains (proposed): the electric **sign** `bc_selection` and **magnitude**
normalization both discharge under the **shared sim-deferred nonlinear throat/core solve** (the same interior solve that
is gravity's `{μ_R,ρ_br}` R10 and magnetism's `q_T`) — record as `PENDING` route-debt, a sibling of R10/R30/R33.

---

## Cross-Part dependencies (so the fold stays honest)

- **Part I (medium):** IV-1's shared bulk/wall/phase pieces DUPLICATE Part I's `χ_B`/`ρ`/`θ`/`r_B` action. Do NOT
  additive-merge — cite Part I as the home, fold only the charge-NEW electric-scalar sector, flag the de-dup as a
  Part VII obligation (model_map §2 "NOT a sum", §3.4, §6).
- **Part III (light):** independent — charge is the STATIC throat; light's `u_T`/`c_γ` transverse sector is not used
  here. (Part V magnetism = the MOVING throat, which DOES reuse `u_T`/`c_γ`.)
- **Part VI (knit):** `c_E` enters here first but its cone relation to `c_γ`/`c_s` is unresolved — `pathA_40` re-adjudicates
  (no committed lock). NG5 `{ρ_B0,χ_c,C_hu}` (pathA_41) touches `B_eff`/`C_hu` consumed here.
- **Part VII (integration):** the G0-vs-Part-I de-dup; the unified action's non-variational drain/return; the shared
  throat-solve reduction that would convert the electric sign/magnitude R1 → resolved; the `F_e/F_g` hierarchy capstone
  (becomes held-out only once BOTH couplings come from one throat action — presently FIT/not-tested).

---

## Per-stage process (unchanged — blueprint §5/§6; LIGHT reshape)

Codex xhigh design-review of the per-stage reshape directive → fold → **Codex confirm-pass** (the first fold is clean) →
Grok-4.5 headless design-review → assess/fold → **final Codex confirm-pass** (the Grok fold is clean) →
**per-stage pre-exec user gate** → Codex builds the two standalone engines (`--sandbox danger-full-access`, background,
xhigh, no shell `timeout`) → dual-engine both exit 0 via independent routes → orchestrator arbiter re-run → full
tri-review on fresh agents (arbiter + fidelity + adversarial-scoped-to-reshape-integrity, per-tooth ablation) →
remediate → **update `parameter_register.md` → Codex-verify the register update → fold (MANDATORY every stage)** →
self-contained note + TeX card + registration → commit + docs/memory sync. Orchestrator authors notes/cards/registration;
**Codex codes** (incl. re-authoring each `.wl` as a genuinely independent route — not a transliteration of the `.py`).

## Gate decisions — RESOLVED (2026-07-22, via the user's conditional build authorization + clean bookend)
1. **Stage count = 3** (recommended; 4-stage alternative not taken). The R1-landing and earned-structure ride together in
   IV-2 as one self-contained puncture-deflection gate.
2. **Build ids = fresh `030/031/032`** (the retired couple-stress 030/031 were never committed as stages).
3. **IV-1 scope = charge-new electric-scalar predicates only**; bulk/wall/drain cited as Part-I/II provenance; G0-vs-Part-I
   de-dup deferred to Part VII. (Concretely: `g0_closure_card_v0_checks.py` is a monolithic all-7-gate checker; 030 folds the
   scalar-block subset.)

⚠ **Reversibility note for the user:** these are my recommendations executed under the overnight build grant. If on review
you'd prefer the 4-stage split or a different IV-1 scope, the stages are independent commits on `ledger-v2-rebuild` and can
be re-split without disturbing gravity/light.

## Progress log
- **2026-07-22 — split plan RATIFIED, bookend CLEAN** (Codex 5+3→CLEAN, Grok CLEAN, Codex GATE_READY). Building 030→031→032.
- **030 (IV-1)** — ✅ **DONE** (verdict `ELECTRIC_SCALAR_CLOSURE_STATIC`; dual-engine SymPy 16 / Mathematica 16, both
  exit 0). Directive bookend CLEAN (Codex 7→1→CLEAN, Grok CLEAN 9 SymPy checks, Codex BUILD_READY). Build → arbiter
  re-run 16/16 → tri-review: **FIDELITY_CLEAN** + arbiter + **ADVERSARIAL_CONCERNS(2)** (two `c_E`-relation value-checks
  were X≡X-tautological in production) → **remediated** (→ `PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY` +
  `PASS_REDUCED_SPEED_PRESERVATION`, both now fire on real error classes) → fresh-agent **REVERIFY_CLEAN**. 16 per-tooth
  ablations fire at own assert (both engines); `.wl` genuine independent route. Scope earned: localized-`H` spectral
  (unique ker(A)+`4/ℓ²` gap), `N₀=8/(3ℓ)`, `h=P₀H`, `M_h=N₀M₄`/`K₄=M₄c_E²`/`K_h=M_h c_E²`, split kernel
  `D*=7/4`+`c_±*²=(3±√2)/2`, `Q_s(0,0)=0`, Hessian symmetry, scalar-block dims. Bulk/wall/drain = out-of-scope
  UNDISCHARGED cross-sector G0 (Part VII de-dup). Deliverables: note + card + Part-IV appendix + register rows/edges
  (R50–R53) + coverage (canonical 30, audited 29) + PDF (82pp). Register Codex-verify: **CLEAN**. **Committed `122afe36`.**
- **⚠ Binding requirement carried to 031:** stage 031 must own a FULL able-to-fail `PASS_BARE_MOUTH_SOURCE` reconstruction
  (∫η=1, f₀(0) from the actual profile, J_m→g_χh, `Q_χ=s`, nonzero projected source) — the puncture engine currently
  assigns f₀(0)=1/ℓ directly, so 031's reshape must add the genuine reconstruction (else the predicate falls through).
- **031 (IV-2)** — ◐ **directive authored + reviewed + folded (v2), ready to resume the bookend.** Directive =
  `_scratch/stage031_reshape_directive.md`. Codex design-review returned `NEEDS_FIXES(5_blocking)` — ALL FOLDED: (1) mouth
  reconstruction now requires the odd-kernel `o_ℓ` + reflection proof `I₋=−I₊`/`Q_±=±1`/`N_χ≠0` (so `Q_χ=s` is PROVED not
  `N_χ/N_χ`) + the nonzero BARE forcing `−g_χh s` (the affine term isn't everywhere nonzero); (2) a mandatory
  source-to-stage predicate manifest (every source tooth preserved/replaced/scoped-out); (3) the landing tooth must mutate
  an UPSTREAM fact + recompute + assert a different landing (not flip the computed summary) + pin the digest literal
  `7627417a…`; (4) MIXED three-regime admissibility; (5) ensemble-data dims (`[φ]=[λ]=1`, `[q]=[j]=[g]=E`, `[S_gg]=E⁻¹`).
  ▶ **STILL OWED (next session):** the rest of the bookend — Codex confirm → Grok compute-verify → final Codex confirm —
  THEN Codex build (dual-engine + the NEW mouth reconstruction) → tri-review → finalize → commit. ⚠ Watch for tautology
  traps (the `Q_χ=s`/`N_χ/N_χ` and landing-ownership risks are the 030-lesson analogues here).
- **032 (IV-3)** — ⬜ pending (`NATIVE_P_NO_EMERGENT_GAUSS`, LIGHT reshape of `native_p_gate_*`; after 031).
