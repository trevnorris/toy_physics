# Part V — Magnetism: atomic-stage split (DRAFT — per-Part gate PENDING)

> **⚠ STATUS: DRAFT proposal, NOT ratified.** Authored 2026-07-22 from the RESUME_ROADMAP §3 re-scope + the memories, right
> after Part IV (charge) completed. **The per-Part user gate + a FRESH read of the magnetism build (`magnetism_moving_throat_result.md`
> + `magnetism_moving_throat_check.{py,wl}`, commit `53cf049f`) will finalize the exact stage boundaries, verdict tokens, and
> register rows** — treat the stage table below as a starting proposal, in the same way `part4_charge_atomic_split.md` was
> drafted then ratified. Do NOT build any Part V stage until this plan is ratified at the gate.
>
> **⭐ Authored under the "ledger = surviving solution only" standing rule.** Part V builds the **surviving** magnetism sector
> = the MOVING ±w throat (the electric twin in motion). It does **NOT** fold the superseded `pathA_39` `j∝sV` scope (barred;
> → failures-paper backlog). The natively-derived source `J_T=q_T sηV` REPLACES it.
>
> Governing: `notes/ledger_v2_blueprint.md` §2 (granularity), §3 Part V row, §5 (reshape spec — LIGHT tier: the build is one
> assert-heavy script pair, like puncture-deflection), §6 (verification protocol); `docs/model_map.md` §3.5 + §4 (departure
> ledger); `RESUME_ROADMAP.md` §3 (the EM re-scope) + §5 (cross-cutting facts). Source build:
> `software/em_charge_attribute/{magnetism_moving_throat_result.md, magnetism_moving_throat_check.py, ...check.wl,
> directive_magnetism_moving_throat.md}` (commit `53cf049f`; truth-table digest `983556935e…`, 1152 cells; 5-leg + GLM verified).

---

## The surviving magnetism story, in one breath

**Magnetism = the MOVING ±w throat** (⟷ gravitomagnetism = the moving drain). A moving throat drags the brane's transverse
shear `u_T` — the SAME field that carries light (`c_γ²=μ_R/ρ_br`, Part III / stage003) — via a natively-derived source
`J_T = q_T s η V` (from signed-throat/defect continuity, **NOT** the barred `j∝sV`). The far field is computed **two blind
ways**: **Route A** boosts the electric interaction (the Maxwell–Darwin reference kernel `(δ_ij+n_in_j)/8πR` to `O(v²/c_γ²)`),
and **Route B** computes the direct moving-throat shear blind to A. **They match in tensor structure, falloff, and velocity
order** — so *the model's magnetism structurally IS the boost of its electricity* (the emergent-Lorentz / boost-consistency
test PASSES at the structural level). Only the coefficient **ratio `r_BA = q_T²/(ρ_br·A_E)`** is open — it rides the
sim-deferred throat normalization `q_T` and the R1 electric coefficient `A_E` — so the magnetic **sign and magnitude INHERIT
the electric R1**: terminal **`R1_REQUIRED(electric_bc_selection)`**. Magnetism coexists with the committed sectors
(`internal_inconsistency = NONE`). And the candidate field **`b_T = ∇×u_T` is time-reversal EVEN** (a real Maxwell `B` is
T-**odd**) — a concrete first-class **characterized departure** ("not exact Maxwell"); `b_T` is correctly axial, and magnetism
requires the throat's active-drain time-arrow `τ_d`.

---

## Proposed atomic-stage split (~6 stages, build ids 034–039 — the gate may MERGE the setup pair)

Following the gravity/charge idiom: EARNED pieces kept **separately citable**, the EARNED structural result split from the
R1-landing, and the characterized departure kept **first-class and separate** (as Part III's stray-longitudinal and Part IV's
`NATIVE_P_NO_EMERGENT_GAUSS` were).

| Build id | Stage | Headline (proposed verdict token) | Scope class |
|---|---|---|---|
| **034** | **V-1** Transverse-vector action row `S_{T+move}` | `TRANSVERSE_MOVE_ACTION_ROW` | **EARNED** (action row) / imports |
| **035** | **V-2** Native moving-throat source `J_T=q_T sηV` (defect-continuity, not `j∝sV`) | `NATIVE_J_T_FROM_CONTINUITY` | **EARNED** (source law, target-blind) |
| **036** | **V-3** Route A — the Maxwell–Darwin reference kernel (boost-electric) | `MAXWELL_DARWIN_REFERENCE` | **EARNED** (reference kernel) |
| **037** | **V-4** Route B (direct shear, BLIND) + the boost-consistency COMPARISON | `BOOST_STRUCTURAL_RELATION_HOLDS` | **EARNED** (structural relation, target-blind) |
| **038** | **V-5** The landing → the magnetic sign/magnitude inherit the electric R1 | `R1_REQUIRED(electric_bc_selection)` | **R1-landing** |
| **039** | **V-6** The `b_T=∇×u_T` time-EVEN departure ("not exact Maxwell") | `B_T_TIME_EVEN_DEPARTURE` | **CHARACTERIZED-DEPARTURE** (first-class) |

**Per-stage scope sketch (to be confirmed against the build at the gate):**
- **034 (V-1):** the G0+δ transverse-vector action row `S_{T+move}=∫[½ρ_br|u̇_T|²−½μ_R|∇×u_T|²+q_T Σ sᵢ η(x−Xᵢ)Vᵢ·u_T]`,
  `∇·u_T=0`, `c_γ²=μ_R/ρ_br`, `q_T=λ_T τ_d`; transverse Hessian stable; **`internal_inconsistency=NONE`** (coexists — no
  pre-existing G0/committed row changed). ⚠ **Fold only the magnetism-NEW moving-defect coupling**; the `u_T` kinetic +
  `c_γ` are Part III / stage003 — CITE as provenance (de-dup deferred to Part VII), exactly as IV-1 cited the shared bulk.
- **035 (V-2):** the natively-derived source `J_T=q_T sηV` from signed-throat / defect continuity (`q_T` the throat charge,
  `s` the ±w orientation from the charge sector, `η` the mouth measure, `V` the throat velocity). Target-blind EARNED; it
  **REPLACES** the barred `j∝sV` (superseded, excluded).
- **036 (V-3):** Route A — boost the electric interaction to get the Maxwell–Darwin reference kernel `(δ_ij+n_in_j)/8πR`,
  general `(V₁,V₂)` to `O(v²/c_γ²)`. The Maxwell-consistent *reference* the direct route is tested against.
- **037 (V-4):** Route B — the direct moving-throat shear, computed **BLIND** to A (carries `{q_T,μ_R}`, no `A_E`) —
  independently reproduces Route A's **tensor structure + falloff + velocity order**. The COMPARISON = the boost-consistency /
  emergent-Lorentz test ⇒ **the structural relation HOLDS**. EARNED target-blind. (Only the coefficient ratio `r_BA` is left
  for 038.)
- **038 (V-5):** the landing. `r_BA=q_T²/(ρ_br·A_E)` is the sole open coefficient — it rides the sim-deferred throat `q_T`
  (R1) and the electric `A_E` (R1). ⇒ the magnetic sign/magnitude are **not independently earned or calibrated**; they
  inherit the electric R1 ⇒ terminal **`R1_REQUIRED(electric_bc_selection)`**. Co-blockers (name, do not resolve): `r_cone`
  (`c_E=c_γ?` — no committed cone lock), higher `v/c` orders, the active-flux / full-force integrability.
- **039 (V-6):** the **`b_T=∇×u_T` time-reversal-EVEN** departure — a concrete not-exact-Maxwell prediction (a real Maxwell
  `B` is T-odd); `b_T` correctly **axial**; magnetism requires the throat's **active-drain time-arrow `τ_d`**.
  CHARACTERIZED-DEPARTURE, first-class (never softened), recorded like Part III's stray-longitudinal + Part IV's `NATIVE_P`.

**⚠ Gate decision — granularity:** 034+035 (action row + source) MAY merge into one "medium+source setup" stage if the build
treats them as one block (compare: IV-1 was a single closure stage); and 036+037 could in principle merge, but keeping Route
A (reference) separately citable from Route B+comparison (the actual boost test) matches the "separately citable" goal. My
recommendation: keep all six (the boost-consistency test is the sector's whole point — each route + the comparison + the
landing + the departure are distinct, citable results). Confirm at the gate after the fresh build read.

---

## What Part V EXCLUDES (surviving-solution compliance)

- **`pathA_39` `j∝sV` scope** — the barred imported current; superseded by the native `J_T=q_T sηV`. → failures-paper backlog
  (do NOT fold). The blueprint's original Part V (`pathA_39`) scope is retired.
- The "hermetic v3" over-engineered magnetism directive that was tried + `git reset`-reverted (history only).

---

## Parameters Part V adds / re-homes (register preview — the per-stage updates fill these; continue edges after R66)

| Param | `[L,T,M]` dim (to confirm) | Enters | Class (proposed) | Note |
|---|---|---|---|---|
| `u_T` (transverse shear vector) | (Part III) | V-1 (034) — consumed | consumed (stage003) | the light-sector field; `∇·u_T=0`; NOT re-counted |
| `μ_R`, `c_γ` | (Part III) | V-1 consumed | consumed | `c_γ²=μ_R/ρ_br`; shared with light; NOT re-counted |
| `q_T = λ_T τ_d` (throat charge) | tbd | V-1/V-2 | `FREE-UNREDUCED` / R1 (throat) | rides the sim-deferred throat solve (sibling of the electric `A_E`, R10/R30) |
| `τ_d` (active-drain time-arrow) | tbd | V-1/V-6 | structural / postulated | the time-arrow magnetism requires; feeds the `b_T` T-even departure |
| `J_T = q_T sηV` (native source) | tbd | V-2 (035) | `DERIVED` (defect-continuity) | REPLACES `j∝sV`; `s`=±w orientation (charge sector), `η`=mouth measure |
| Maxwell–Darwin kernel `(δ_ij+n_in_j)/8πR` | structural | V-3 (036) | `DERIVED` (boost of electric) | the reference kernel to `O(v²/c_γ²)` |
| `r_BA = q_T²/(ρ_br·A_E)` (coefficient ratio) | dimensionless | V-4/V-5 | **R1-deferred** (`electric_bc_selection`) | the sole open coefficient; rides R1 `q_T` + R1 electric `A_E` |
| `r_cone` (`c_E=c_γ?`) | dimensionless | V-5 (038) | **R1/OPEN** | no committed cone lock (Part VI `pathA_40` re-adjudicates) |
| `b_T = ∇×u_T` (candidate mag. field) | tbd | V-6 (039) | **DEPARTURE** (T-even; discharges no knob) | first-class characterized departure, NOT a reduction — must not shrink the irreducible count |

New reduction-debt edges: the magnetic sign/magnitude (`r_BA`) discharge under the **same shared sim-deferred throat solve**
as the electric R1 (`R1_REQUIRED(electric_bc_selection)` — sibling of the charge R1s / R10/R30). Record as `PENDING`.

---

## Cross-Part dependencies

- **Part III (light):** V-1 REUSES `u_T`/`c_γ`/`μ_R` (stage003) — cite, do NOT re-derive/re-count; de-dup at Part VII.
- **Part IV (charge):** the magnetic sign INHERITS the electric R1 (`R1_REQUIRED(bc_selection)`, stage 032) — the same
  sim-deferred throat solve resolves both; `s`=±w orientation is the charge sector's topological sign. The `b_T` departure is
  the magnetic twin of `NATIVE_P_NO_EMERGENT_GAUSS` (stage 033) — both first-class "not exact Maxwell".
- **Part VI (knit):** `r_cone` / `c_E=c_γ?` is `pathA_40`'s re-adjudication (no committed lock). `q_T` throat normalization
  ties to the shared throat solve.
- **Part VII (integration):** the shared throat solve that would convert both the electric AND magnetic R1 → resolved; the
  `F_e/F_g` hierarchy capstone; the G0-vs-Part-I de-dup.

---

## Per-stage process (unchanged — blueprint §5/§6; LIGHT reshape)

Per stage: Codex→Grok→Codex directive bookend (fold via agents) → per-stage pre-exec user gate → Codex build
(`--sandbox danger-full-access`, detached `setsid` + poll — the run_in_background reap workaround) → dual-engine both exit 0,
independent routes → orchestrator arbiter re-run → fresh-agent tri-review (fidelity + adversarial per-tooth ablation) →
remediate (only real coverage gaps; document verified-safe smells) → register update + Codex-verify → self-contained note +
TeX card + registration → deliverables fidelity-verify (agent-authored + Codex-verified) → commit + tracker/memory sync
(separate commit, real hash). Orchestrator authors notes/cards/registration; Codex codes; agents fold review feedback.

## Open questions for the per-Part gate
1. **Stage count 5 vs 6** — merge 034+035 (setup) or keep separate? (Recommend keep — citability.)
2. **Build-id start = 034** (Part IV ended at 033) — confirm.
3. **Verdict tokens** — the proposed tokens are drafts; confirm against the build's actual `magnetism_moving_throat_result.md`
   Q-blocks/teeth at the fresh read.
4. **The `b_T` T-even departure** — confirm it's a standalone stage (V-6) vs folded into the V-5 landing (recommend standalone,
   matching the first-class-departure precedent).
