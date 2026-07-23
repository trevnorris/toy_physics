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
> `GATE_READY`. **Locked choices** (updated 2026-07-22 by user decision): **4 stages** (030/031/032/033) — 031 splits its
> EARNED mechanism/far-field structure from the `R1_REQUIRED(bc_selection)` sign landing (032); the native-P departure
> becomes 033; build ids fresh `030/031/032/033`; **IV-1 scope = fold ONLY the charge-new electric-scalar predicates**,
> bulk/wall/drain predicates cited as Part-I/II provenance (de-dup deferred to Part VII). ⚠ Concrete-source confirmation: `g0_closure_card_v0_checks.py` is a
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
the **sign is not invariant across the admissible mouth boundary-condition ensembles** (on a nondegenerate response: V repel / J attract / M indefinite),
and the committed bare model does not select one → the honest terminal landing is **`R1_REQUIRED(bc_selection)`**
(neither earned NOR calibrated; needs the sim-deferred nonlinear throat/core solve). Charge coexists with gravity on the
same brane (`internal_inconsistency = none`, computed). And exact U(1)/Maxwell Gauss is proven **non-native**
(`NATIVE_P_NO_EMERGENT_GAUSS`) — a first-class characterized departure: **EM is not exact Maxwell.**

---

## ~~RECOMMENDED split = 3 stages~~ — **SUPERSEDED (user chose the 4-stage split, 2026-07-22)**

> ⚠ **The 3-stage split below is NOT the plan.** The user reversed the overnight 3-stage recommendation in favour of the
> **ADOPTED 4-stage split** (next section): 031 splits along the EARNED/R1 scope boundary — the target-blind EARNED
> far-field structure (031) separately citable from the `R1_REQUIRED(bc_selection)` sign landing (032), matching gravity's
> `pathA_29→009/010` and `pathA_34→022/023` idiom + the blueprint §2 granularity rule ("a gate that bundles separable
> sub-derivations decomposes into several stages… each step citable"). The 3-stage table is kept below only as the source
> of the 030/whole-gate scope detail that is re-partitioned into the adopted 031/032. **Read the ADOPTED split.**

Each stage folds one of the three already-verified `em_charge_attribute/` script pairs (LIGHT reshape — they are
already assert-heavy print/`raise`-style engines; folding = wrap into the 6-artifact stage form + confirm the `.wl` is a
genuinely independent route).

| Build id | Stage | Source build + scripts | Headline verdict token | Scope | Reshape / notes |
|---|---|---|---|---|---|
| **030** | **IV-1** Electric scalar block + localized-`H` closure | `g0_closure_card_v0.md` §§2.3–2.4 + `g0_closure_card_v0_checks.{py,wl}` | class-(1) static PASS bundle (`PASS_LOCALIZED_H_NORM`, `PASS_REDUCED_H_INERTIA`, `PASS_STABILITY`, `PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS`, `PASS_DIMENSIONAL_HOMOGENEITY`) | **EARNED** (static prerequisites of the charge substrate): the reduced `(u_L,h)` action `S_Lh`; the parent `H` Pöschl–Teller operator `O_⊥=A†A≥0` with **one** gapless localized zero mode `f₀=1/[ℓcosh²(w/ℓ)]`, `N₀=8/(3ℓ)`, `h=P₀H`; `M_h=N₀M₄`, `K_h=N₀K₄=M_h c_E²`; the positive coupled kernel `D=B_eff K_h−C_hu²>0` (`[D]=M²T⁻⁴`) with **dimensionless witness** `D_*=7/4` and `c_±²=(3±√2)/2>0`; dim firewall. **POSTULATED** (labeled): the whole G0-v0 closure is a DRAFT-v0 postulated closure. | LIGHT. ⚠ **Fold only the charge-NEW electric-scalar sector**; the shared bulk/wall/phase pieces (`ρ,θ,r_B`) DUPLICATE Part I → cite as provenance, **de-dup DEFERRED to Part VII** (reconcile, not additive-merge). IV-1 ends at the G0 kernel stability + witness `D_*=7/4`; the response matrix `m`/`m_gg`/`z_g`/`det m` are IV-2 (puncture-deflection) objects. |
| **031** | **IV-2** Puncture-deflection: field identity → `1/R²` force → `R1_REQUIRED(bc_selection)` | `puncture_deflection_electric_sign_result.md` (all Q-blocks + App A–E) + `puncture_deflection_electric_sign_check.{py,wl}` + `puncture_deflection_electric_sign_independent_recompute.py` | **`R1_REQUIRED(bc_selection)`** | **EARNED (target-blind structure):** field identity `ξ_w=ℓh`; live coupling reduction `−J_m Q_χ H → −g_χh Q_χ h` (since `f₀(0)=1/ℓ`); the orientation-odd mouth source `η_i(k_m h − g_χh s_i)`, `Q_χ[r_Σ,s_i]=s_i` ⇒ **charge sign = ±w Z₂ orientation**; exterior `h(r)=h_A a/r` (positive holding curvature `4π(D/B_eff)a`) ⇒ the `h`-mediated `1/R` potential / `1/R²` force with the `s₁s₂` product form. The completed-square response matrix `m` (from the puncture-deflection build) gives the pair coupling `m_gg=B_eff z_g²/D` with the Robin escape factor `z_g>0` and `det m=z_g²/D>0`. **R1 (sign):** the four admissible BC classes give **decided conditional** coefficients (strict signs hold on a nondegenerate response `m_gg>0`/`S_gg>0`; weak in general) `A_V=m_ggφ²/S_gg²≥0` (repel; `>0` at the witness), `A_J=−m_gg(j+g)²≤0` (attract; `<0` at the witness), `A_M=m_gg(q²−g²)` (indefinite), `A_MIXED=m_gg[(1−2λ)q²−2λqg−g²]` (spans); the committed bare model (G0 + U2-unresolved `𝔅` + S_hold scoped to `r_B−½`) does not select the class ⇒ `outcome_not_invariant` ⇒ terminal **`R1_REQUIRED(bc_selection)`** — the sign is **neither EARNED nor CALIBRATED**. `internal_inconsistency = none`. Co-blocker: `R1_REQUIRED(magnitude)` (`c_a`,`c_ξ` core normalizations unbounded at tier-A). | LIGHT. The whole self-contained puncture-deflection gate (one script pair). Carry the sealed **23040-cell** exhaustive §4 landing table as the able-to-fail backbone (first-matching-predicate ladder → the production tuple lands `R1_REQUIRED(bc_selection)`). |
| **032** | **IV-3** `NATIVE_P_NO_EMERGENT_GAUSS` (the characterized departure) | `reports/native_p_constraint_gate.md` + `native_p_gate_sympy.py` + `native_p_gate_dual.wl` (+ `native_p_gate_compare.py`) | **`NATIVE_P_NO_EMERGENT_GAUSS`** | **CHARACTERIZED-DEPARTURE (first-class):** the genuine Dirac–Bergmann constraint-class gate — at quadratic order native-`P` theories A and C both have **FC=0** (every constraint second-class; PB rank 8 / 12, no first-class Gauss chain); the tuned rank-drop strata that carry FC directions have a **zero Gauss descendant** (`DESCENDANT_ZERO`, hardening guard PASS). Exact U(1)/Maxwell Gauss is **non-native** → "EM is NOT exact Maxwell." ⚠ **Scope (carry honestly, does not weaken the verdict):** the regular open kinetic stratum is *fully symbolic* (FC=0 for all retained couplings); the non-common tuned rank-drop locus is *argued + fixed-seed scanned* (6 representative points per theory, 12 total across A+C), NOT an exhaustive symbolic stratification — any hypothetical missed measure-zero Gauss stratum would be a TUNED/inverse-design artifact, so the physical no-go (native `P` does not *generically* host emergent EM) is decisive independently. | LIGHT. Six able-to-fail controls anchor it (Maxwell→`FIRST_CLASS_GAUSS`, gauged-hard-unit→`MIXED`, nonconserved-current→`INCONSISTENT_PRESERVATION`, Coulomb-gauge→`SECOND_CLASS_NO_LOCAL_GAUGE`, global-U(1)→`GLOBAL_CHARGE_NO_LOCAL_GAUSS`, bare-σ→`SECOND_CLASS_RADIAL_NO_GAUSS`) + per-tooth ablation (`FIRED_AT_OWN_ASSERT`). |

**Net RECOMMENDED Part IV = 3 stages (030–032).** Build ids continue the global build-order sequence (Part III's
surviving stage is `003`; the couple-stress no-go that once held 030/031 was DROPPED to the failures paper, so those
numbers are free — but to avoid any collision/confusion with the retired 030/031, **I propose new build ids 030/031/032
are used fresh here** and note the retired numbers were never committed as stages). *(Build-id assignment is a
mechanical detail confirmable at the gate; the Part-ordered linear renumber is deferred to the pre-swap manifest pass.)*

---

## ✅ ADOPTED split = **4 stages** (user decision 2026-07-22)

Split IV-2 into its EARNED-structure half and its R1-landing half — the ledger's established idiom for joint gates
(Part II split `pathA_30`→011 earned-reduction + 012 the BC-dependent landing; `pathA_29`→009 residual + 010
localization; `pathA_34`→022 EARNED-first + 023 the FAIL). Costs one extra reshape: Codex splits the single
`puncture_deflection_electric_sign_check.{py,wl}` into two standalone stage-script pairs (precedent exists — 011/012 came
from one `pathA_30` source; 016/017 from one `pathA_32`).

| Build id | Stage | Headline | Scope class |
|---|---|---|---|
| **030** ✅ committed | IV-1 Electric scalar block + localized-`H` closure | `ELECTRIC_SCALAR_CLOSURE_STATIC` | EARNED / POSTULATED |
| **031** | IV-2 Field identity + `χ↔h` mouth source + exterior `1/R²` + response matrix `m` | `THROAT_H_SOURCE_1_OVER_R2` (proposed token) | **EARNED (target-blind structure)** |
| **032** | IV-3 Far-field force ensembles → the honest landing | `R1_REQUIRED(bc_selection)` | **R1 (sign)** |
| **033** | IV-4 `NATIVE_P_NO_EMERGENT_GAUSS` | `NATIVE_P_NO_EMERGENT_GAUSS` | CHARACTERIZED-DEPARTURE |

**Why (user rationale, 2026-07-22).** 4-stage cleanly separates the target-blind EARNED far-field structure from the
R1(sign) landing — they are **different scope classes**, and keeping them **separately citable** for extracted papers is
the whole point of the granularity rule (goal 3). Fusing an EARNED result with its R1/FAIL landing is precisely what
gravity NEVER did (009/010, 016/017, 022/023). The overnight 3-stage recommendation over-weighted "least reshape surface"
and mis-read "a self-contained gate stays one stage" — a self-contained gate (one script pair) still decomposes when its
*derivation* has separable steps (pathA_30/32/43 each split despite being one source). The extra reshape cost is cheap
here because 031/032 are not yet built. **Boundary:** 031 owns the mechanism (field identity `ξ_w=ℓh`, the full
bare-mouth reconstruction `Q_χ=s`, exterior `h(r)=h_A a/r` → the `1/R²` `s₁s₂` structure, the completed-square response
matrix `m`/`m_gg`/`z_g`/`det m`); 032 consumes 031's `m`/`m_gg`/`S_gg` + the neutral shell `A=m_gg C` and owns the four BC ensembles + the sealed
23040-cell landing → `R1_REQUIRED(bc_selection)`.

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
`pathA_38` → backlog. This is a register-hygiene edit the **IV-3/032** register update carries (the `Q_E`/magnitude re-home
is 032's under the 4-stage split, alongside the `A_X`/ensemble/`internal_inconsistency` landing rows).

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
| `z_g` (Robin escape factor) | `1` | IV-2 (031) | formula `DERIVED`; strict `z_g>0` **POSTULATED** (Robin-admissibility witness `z_g=1`) unless independently proved | `z_g=1−k_m⟨η,L_h⁻¹η⟩`; the source only sets a positive symbol (Codex finding 3) |
| `m_gg = B_eff z_g²/D` (pair coupling) | `M⁻¹ L⁻¹ T²` | IV-2 (031) | `DERIVED` | the two-body response amplitude |
| `det m = z_g²/D` | `M⁻² T⁴` | IV-2 (031) | formula `DERIVED`; **`≥0` earned** (via `z_g²`, `D>0`); strict `>0` only at the `z_g≠0` witness (POSTULATED) | positive response at the admissible Robin witness |
| the class numerator `C` (`A=m_gg·C`) + `A_V/A_M/A_J/A_MIXED` | `[C]=E²`; `[A]=EL=M L³ T⁻²` (result:178) | **IV-3 (032)** | `DERIVED` (**decided** conditional formulas) | ⭐ 4-stage split: 031 emits the NEUTRAL shell `A=m_gg·C` (EARNED); **032 selects `C`** per class `A_V=m_ggφ²/S_gg²`, `A_M=m_gg(q²−g²)`, `A_J=−m_gg(j+g)²`, `A_MIXED=…` — NOT R1; the **sign** R1 is only WHICH class holds (strict signs CONDITIONAL on a nondegenerate response); magnitude + MIXED-param/variant are SEPARATE downstream R1 debts |
| `S_gg` (self-response) | `E⁻¹=M⁻¹L⁻²T²` | IV-2 (031) | `DERIVED` **only with explicit definition + tooth**, else `POSTULATED` | enters 032's `A_V`; the source only introduces it as a positive symbol — 031 must DEFINE it or label it POSTULATED (Codex finding 3) |
| BC-class selection + ensemble data `{φ, q, j, λ}` | `[φ]=[λ]=1`, `[q]=[j]=E` | **IV-3 (032)** | **R1-deferred (`bc_selection`)** | the *selection* among V/M/J/MIXED is the R1 (`outcome_not_invariant`); `{φ,q,j,λ}` parameterize the unselected classes — tracked, not counted |
| `internal_inconsistency = none` (Q-AMEND) | — | **IV-3 (032)** | computed | charge coexists with gravity; owned + verified at 032 (from Q-AMEND consistency), NOT 031 |
| `s_i = ±1` (Z₂ orientation charge) | `1` | IV-2 (031) | `CONV`/structural (topological, **not a tuned knob**) | charge sign; ±w orientation, not additive |
| `Q_E` / charge-magnitude | pending | **IV-3 (032)** re-home | **R1-deferred (`magnitude`)** | **re-home from pathA_38** → puncture-deflection; `c_a`,`c_ξ` unbounded at tier-A |

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

## Gate decisions — RESOLVED
1. **Stage count = 4** ✅ **(user decision 2026-07-22, reversing the overnight 3-stage recommendation).** 031 splits along
   the EARNED/R1 scope-class boundary — the target-blind EARNED mechanism/far-field structure (031) is kept **separately
   citable** from the `R1_REQUIRED(bc_selection)` sign landing (032); the native-P departure becomes 033. Matches gravity's
   `pathA_29→009/010`, `pathA_32→016/017`, `pathA_34→022/023` idiom + the blueprint §2 granularity rule. *(The overnight
   3-stage call — kept in the SUPERSEDED section above — over-weighted "least reshape surface" and mis-read "a self-contained
   gate stays one stage"; a self-contained one-script-pair gate still decomposes when its derivation has separable steps.)*
2. **Build ids = fresh `030/031/032/033`** (the retired couple-stress 030/031 were never committed as stages).
3. **IV-1 scope = charge-new electric-scalar predicates only**; bulk/wall/drain cited as Part-I/II provenance; G0-vs-Part-I
   de-dup deferred to Part VII. (Concretely: `g0_closure_card_v0_checks.py` is a monolithic all-7-gate checker; 030 folds the
   scalar-block subset.) — **unchanged by the 4-stage decision** (030 is committed as-is; the split is on IV-2, not IV-1).

## Progress log
- **2026-07-22 — split plan RATIFIED (3-stage), bookend CLEAN** (Codex 5+3→CLEAN, Grok CLEAN, Codex GATE_READY).
- **2026-07-22 — ⭐ RE-GRANULARIZED to 4 stages (user decision).** On the morning inventory review, the user flagged that
  030/031 bite off more than gravity's stages (which stayed small to cite specific parts in future papers). Adopted the
  4-stage split: 030 (committed, unchanged) → 031 mechanism/EARNED → 032 ensembles+R1 landing → 033 native-P departure.
  ▶ Next = re-carve the 031 directive into a MECHANISM directive (031) + a LANDING directive (032); rerun the bookend on
  each; the native-P directive is authored fresh as 033 when we reach it.
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
- **031 (IV-2, MECHANISM/EARNED — 4-stage)** — ✅ **DONE + committed `bba61c82`** (verdict `THROAT_H_SOURCE_1_OVER_R2`,
  EARNED within the postulated G0 closure). Owns: field identity `ξ_w=ℓh`; the FULL bare-mouth reconstruction (odd-kernel
  `o_ℓ` + reflection `I₋=−I₊` + the **independent generic `I₊>0` REFLECTION-DOMINANCE certification** so `Q_χ=s` is PROVED
  not `N_χ/N_χ`; the live bare forcing); exterior `h(r)=h_A a/r` + the exported `NONZERO_HA_REQUIRES_CORE_HOLDER`; the FULL
  response matrix `m` (heterogeneous dims, star `m_uu=K_h/D=4/7`); the target-blind neutral far-field shell `A=m_gg C`.
  **Directive bookend:** Codex `FIXES_NEEDED(6)+4nits` → confirm chain 6→2→1→CLEAN → Grok compute-verify CLEAN → `BUILD_READY`.
  **⭐ Build-1 HONESTLY HIT the `I₊≠0` STOP-gate** (the source leaves `r_{Σ,+}` an abstract `[POSTULATE]`; the source build's
  `Q_χ=s` rested on an uncertified step); STOP triple-confirmed. **User chose OPTION A** (complete the G0 postulate): a
  bounded, PEAK-NORMALIZED `r_{Σ,+}∈(0,1]` one-sided-sleeve representative (frozen `L_s`; NO new knob class) + certify `I₊>0`
  generically by **REFLECTION DOMINANCE** (`D(w)≥D(−w)` for `w>0` ⇒ `I₊=∫_{w>0}o_ℓ[D(w)−D(−w)]>0`, structural not tuned;
  even-deformation/`L_s=0`/`−w`-sleeve = load-bearing teeth; `STAGE031_STOP: I_plus_not_generic` anti-rescue fallback). Option
  A **Grok-B compute-verified CLEAN** (certificate structural; the boundedness FAIL fixed via peak-normalization). Rebuilt
  green: dual-engine SymPy 50 / Mathematica 50, both exit 0; **tri-review CLEAN both legs** (50/50 mutations fire, `I₊>0`
  parameter-generic, `N_χ/N_χ` defeated), ZERO remediation. Deliverables: note (with the full **36-tooth source-to-stage
  manifest — no silently-dropped source claim**) + card + register rows/edges **R54–R61** + coverage (canonical 31, audited
  30) + PDF (84pp); deliverables fidelity-verify `7→CLEAN` (caught a `κ`-overclaim: the exterior stiffness is a GENERIC `κ`,
  NOT `D/B_eff`; `z_g`/`z_b` opaque not DERIVED). ⚠ The `Q_χ=s`/`N_χ/N_χ` trap was defeated by the certified `I₊>0`.
- **032 (IV-3, ensembles + R1 LANDING — 4-stage)** — ✅ **DONE + committed `11253849`** (verdict
  **`R1_REQUIRED(bc_selection)`** — the electric SIGN is NEITHER earned NOR calibrated). Owns: the four DECIDED conditional
  coefficients `A_V=m_ggφ²/S_gg²` (repel) / `A_J=−m_gg(j+g)²` (attract) / `A_M=m_gg(q²−g²)` (indefinite) / `A_MIXED` (spans),
  consuming 031's `m_gg`/`S_gg` + the neutral shell by CITATION; **strict signs CONDITIONAL on the nondegenerate witness**
  (weak in general, all vanish at `z_g=0`); the MIXED three-regime admissibility; `internal_inconsistency=none` (Q-AMEND);
  the sealed **23040-cell** §4 landing ladder (digest `7627417a…`, first-match `bc_selection`, before variant/tier/magnitude)
  + genuine `LANDING_OWNERSHIP` (mutates an UPSTREAM Q-BC fact + re-adjudicates to a named different landing — passes the 030
  X≡X test); co-blocker `R1_REQUIRED(magnitude)`; the `Q_E`/magnitude re-home (pathA_38 → puncture-deflection). **Bookend:**
  Codex `SOUND` (pass-1, zero fixes; reran both 031 engines 50/50) → Grok compute-verify `CLEAN` (SymPy-verified the four
  `A_X`, MIXED regimes, ladder order, `23040=2⁹·3²·5`) → Codex `BUILD_READY`. **Build:** dual-engine SymPy 57 / Mathematica
  57 (independent Solve/Reduce + native exhaustive table), both exit 0, 58 mutations all fire; arbiter re-run CLEAN;
  **tri-review CLEAN both legs** (fidelity + adversarial 57/57, `LANDING_OWNERSHIP` genuine). **⭐ Both legs flagged the same
  weak spot** — the mandatory source-to-stage manifest tooth pinned coverage/ids/dispositions but NOT the per-disposition
  counts or the manifest digest (a count-preserving mis-disposition would slip) → **manifest-hardening remediation** (pin
  counts 22/12/10 + manifest digest `e2cfd11b…` over the sorted `(id,disposition,owner)` triples, both engines agree; add the
  count-preserving `MANIFEST_MISDISPOSITION` ablation) → arbiter + **fresh-agent re-verify CLEAN** (58/58, both engines agree
  on both digests). Deliverables authored (delegated to a clean-context agent, then Codex fidelity-verify `5→CLEAN` — caught
  prosified owner tokens, a stale `C`/`Q_E` classification, an edge-range typo, "zero remediation" wording; all folded).
  Register edges **R62–R65**; coverage canonical 31→32; PDF 86pp. ⚠ Infra: the remediation codex build was launched
  **detached via `setsid`** after two `run_in_background` reaps — watch the real launcher PID, not the `setsid` wrapper.
- **033 (IV-4, native-P departure — 4-stage)** — ✅ **DONE + committed `391da099`** (verdict
  **`NATIVE_P_NO_EMERGENT_GAUSS`** — the first-class CHARACTERIZED-DEPARTURE "EM is NOT exact Maxwell"). Exact U(1)/Maxwell
  Gauss is non-native to the P-field medium: at quadratic order native families A `{g_tA,g_sA,g_dA,g_bA}` (Hessian 8/2, 8
  constraints, PB rank 8/FC 0/SC 8) and C `{g_tC,g_sC,g_dC}` (7/4, 12 constraints, rank 12/FC 0/SC 12) are ALL second-class;
  det `(ρ_u−g_t²)³`, sole degeneracy `g_t²=ρ_u`. **The crux:** the tuned common-null points DO develop FC=2, but those
  directions have a **ZERO Gauss descendant** (`DESCENDANT_ZERO`) — so the no-go follows from `additional_G_exists=False` (no
  `∝k` Gauss candidate), NOT from an absence of FC directions. Six able-to-fail controls (maxwell→`FIRST_CLASS_GAUSS`,
  gauged-hard-unit→`MIXED`, bare-σ→`SECOND_CLASS_RADIAL_NO_GAUSS`, nonconserved→`INCONSISTENT_PRESERVATION`,
  Coulomb-gauge→`SECOND_CLASS_NO_LOCAL_GAUGE`, global-U(1)→`GLOBAL_CHARGE_NO_LOCAL_GAUSS`) prove the search is capable.
  Honest scope: symbolic OPEN stratum decisive (FC=0 all retained couplings) + tuned locus ARGUED+SCANNED (6-point boundary +
  12-total sweep, seeds A=260713/C=260715), NOT exhaustive (missed stratum ⇒ TUNED/inverse-design only ⇒ generic no-go
  stays). **Directive bookend:** Codex `NEEDS_FIXES(2B+2N)` → agent-folded → confirm `NEEDS_FIXES(1)` → fixed → Grok
  compute-verify CLEAN → `BUILD_READY`. **⭐ The source `.wl` was a MIRROR** (hardcoded the tuning the `.py` solves) — caught
  by the directive review; the stage `.wl` was **RE-AUTHORED** as a genuinely independent route (native
  `Solve`/`Reduce[ForAll]` common-null + its own six native control ablations). **Build:** dual-engine SymPy 33 / Mathematica
  33, both exit 0; verdict for A and C; the verdict tooth re-derives to `FIRST_CLASS_GENERIC_EM_CANDIDATE` under mutation
  (non-tautological, BOTH engines). **Tri-review CLEAN + adversarial CONCERNS(2)** — both non-blocking robustness smells with
  NO coverage gap (`HONEST_TUNED_SCOPE` framing-tooth is redundant with the computed-data `BOUNDARY_SCAN`/`RANDOMIZED_SWEEP`/
  `SIGNATURE` teeth; the SymPy dual-path is adversarial-verified-safe) → DOCUMENTED, not remediated (a principled split from
  the 032 manifest, which was a sole-guard gap). Source-to-stage manifest 51 rows (33 PRESERVED / 15 REPLACED_BY_STRONGER / 3
  SCOPED_OUT), SHA256 `6b191e77…` (both engines agree — already stronger than 032: pins order+counts+digest). Deliverables
  authored (delegated to a clean-context agent, then Codex fidelity-verify **3→CLEAN** — caught `½λ̇²`→`½ḃ²`, the
  "no first-class direction"→"no Gauss candidate" crux, and a manifest-partition prose error). Register edge **R66** = the
  DEPARTURE (discharges NO knob, NOT a reduction/codim — does not shrink the irreducible count) + a departure-support
  structural row + pathA_36 cross-ref; coverage canonical 32→33, Part IV row `030-033`; PDF 89pp.
- **✅ PART IV (CHARGE) COMPLETE** (030 `122afe36` / 031 `bba61c82` / 032 `11253849` / 033 `391da099`). The surviving charge
  sector: charge = the static ±w throat; `Q_χ=s` (topological, reflection-dominance `I₊>0`); `h`-mediated `1/R²` structure
  EARNED target-blind; the SIGN is `R1_REQUIRED(bc_selection)` (neither earned nor calibrated); the magnitude is
  `R1_REQUIRED(magnitude)`; exact Maxwell Gauss proven non-native (`NATIVE_P_NO_EMERGENT_GAUSS`, first-class departure). ▶ NEXT
  front = **Part V (magnetism)** (the moving throat; per `RESUME_ROADMAP.md` §3/§7).
