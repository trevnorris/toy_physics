# Part V — Magnetism: atomic-stage split (✅ RATIFIED — 6-stage split, 2026-07-22)

> **✅ STATUS: RATIFIED at the per-Part gate (2026-07-22) — all SIX stages 034–039 as tabled.** First drafted 2026-07-22 from
> the RESUME_ROADMAP §3 re-scope + memories; then **REFINED against a FRESH read of the magnetism build**
> (`magnetism_moving_throat_result.md` + a tooth-inventory audit of `magnetism_moving_throat_check.{py,wl}`, commit `53cf049f`);
> then **ratified by the user** (kept all six — the three granularity calls resolved to KEEP SEPARATE: action row / source /
> Route A / Route B+comparison / landing / departure). Verdict tokens, stage boundaries, register dims, and the **35-tooth →
> stage allocation** are build-faithful.
> **Next:** ✅ **034 (V-1) DONE + committed `109070da`** (see the Progress log at the bottom) → per-stage pre-execution gate on
> **035 (V-2)** → Codex→Grok→Codex directive → build. Per-stage builds proceed under the autonomy mode set at ratification.
>
> **Fresh-read corrections folded in:** (1) V-2 token is the build's own `CONVECTION_LIKE_CONDITIONAL` (native source earned in
> *tensor form*, conditional because `q_T` is unfixed) — not the earn-gloss `NATIVE_J_T_FROM_CONTINUITY`; (2) V-6 token is the
> build's exact `B_TIME_REVERSAL_EVEN`, and the T-even *parity is DECIDED* while the sector *landing* is R1 (V-5≠V-6, clean);
> (3) the comparison **ratios reuse BOTH routes** → they land in V-4 (computed/earned), the SEALED §4 first-match adjudication is
> V-5 (the landing); (4) V-5 names **four** emitted R1 co-blockers; (5) register dims pinned from the build's Appendix A;
> (6) build = **35 teeth** (`TOOTH_ORDER`, py/wl identical), allocated per-stage below.
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

| Build id | Stage | Headline verdict token (build-faithful) | Scope class |
|---|---|---|---|
| **034** | **V-1** Transverse-vector action row `S_{T+move}` | `TRANSVERSE_MOVE_ACTION_ROW` (`internal_inconsistency=none`) | **EARNED** (action row) / imports |
| **035** | **V-2** Native moving-throat source `J_T=q_T sηV` (defect-continuity, not `j∝sV`) + parity census | `CONVECTION_LIKE_CONDITIONAL` | **EARNED** (source law, target-blind) |
| **036** | **V-3** Route A — the Maxwell–Darwin reference kernel (boost-electric) | `MAXWELL_DARWIN_REFERENCE` | **EARNED** (reference kernel; inherits `A_E` R1) |
| **037** | **V-4** Route B (direct shear, BLIND) + structural COMPARISON + ratios (`r_BA`, `δ_BA`, `r_cone`) | `BOOST_STRUCTURAL_RELATION_HOLDS` (+ co-blocker `R1_REQUIRED(direct_moving_throat)`) | **EARNED** (structural relation, target-blind) |
| **038** | **V-5** The SEALED §4 first-match landing (1152 cells) → magnetic sign/magnitude inherit the electric R1 | `R1_REQUIRED(electric_bc_selection)` (+3 co-blockers) | **R1-landing** |
| **039** | **V-6** The `b_T=∇×u_T` time-EVEN departure ("not exact Maxwell") | `B_TIME_REVERSAL_EVEN` | **CHARACTERIZED-DEPARTURE** (first-class) |

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
- **037 (V-4):** Route B — the direct moving-throat shear `U_B=−s₁s₂q_T²(D_V+A_V)/8πμ_R R`, computed **BLIND** to A
  (`foreign_payload=None`, built FIRST in the source code to guarantee blindness; carries `{q_T,μ_R}`, no `A_E`) —
  independently reproduces Route A's **tensor structure + falloff + velocity order**. The COMPARISON = the boost-consistency /
  emergent-Lorentz test ⇒ **the structural relation HOLDS**. The ratio **expressions** `r_BA=q_T²/(ρ_br·A_E)`, `δ_BA=r_BA−1`,
  `r_cone=c_E²/c_γ²` and `ΔU` are **computed here** (DECIDED expressions — the comparison reuses BOTH routes) — but their
  *values* are unresolved (`q_T`, `A_E` open), tagged `QMAG_R1`. Route B's own magnitude co-blocker
  **`R1_REQUIRED(direct_moving_throat)`** originates in this stage. EARNED target-blind.
  ⚠ *Presentational note for the note/card:* the ledger presents Route A (V-3, reference) before Route B (V-4, test), but the
  build computes B **before** A for blindness — state this explicitly so a reader never reads the ordering as "B saw A."
- **038 (V-5):** the SEALED §4 landing — a **total first-match precedence** adjudication over the **1152-cell** truth table
  (`4×3×4×3×2×2×2`), cross-checked by an independent oracle, with a **runtime-computed** cross-engine digest (`983556935e…`
  is the production value, NOT a literal in the code). The production tuple's first match is terminal
  **`R1_REQUIRED(electric_bc_selection)`**; the complete blocker collector then emits **all four** honest R1 rows:
  `R1_REQUIRED(electric_bc_selection)` (primary) · `R1_REQUIRED(direct_moving_throat)` · `R1_REQUIRED(magnitude)` ·
  `R1_REQUIRED(consistency)`. ⇒ the magnetic sign/magnitude are **not independently earned or calibrated**; they inherit the
  electric R1. Named-not-resolved co-blockers folded into these four: `r_cone` (`c_E=c_γ?` — no committed cone lock, Part VI
  `pathA_40` re-adjudicates), higher `v/c` orders (`O(v⁴/c_γ⁴)`), the active-flux `F_flux` (`O(V₁V₂)`) + full-force
  integrability (`ACTIVE_FLUX_CAVEAT`), and emergent Lorentz = `UNDETERMINED` (`HOOK_LORENTZ`).
- **039 (V-6):** the **`b_T=∇×u_T` time-reversal-EVEN** departure — a concrete not-exact-Maxwell prediction (a real Maxwell
  `B` is T-odd); `b_T` correctly **axial**; magnetism requires the throat's **active-drain time-arrow `τ_d`**.
  CHARACTERIZED-DEPARTURE, first-class (never softened), recorded like Part III's stray-longitudinal + Part IV's `NATIVE_P`.

**⚠ Gate decisions — granularity (fresh read done; three calls for the user):**
1. **034+035 merge?** (action row + source). They are two labeled blocks in the build (Provenance/action 137–275; Q-CURRENT
   278–340) and two distinct EARNED claims ("clean G0+δ row" vs "native source law + parity"). Recommend **keep separate**
   (citability), but they MAY merge into one "medium+source setup" stage (compare IV-1, a single closure stage).
2. **036+037 merge?** The comparison/ratios reuse **both** routes, so V-4 already re-imports Route A's kernel regardless.
   Recommend **keep Route A (V-3, reference) separate from Route B+comparison (V-4, the test)** — the two *independent* routes
   matching IS the boost-consistency result, and making the independence separately-citable is the epistemic point.
3. **039 departure — standalone vs folded?** The `b_T` T-even parity is DECIDED in V-2's census; V-6 is a *characterization*
   that cites it (like 033 cited earlier charge structure). Recommend **standalone V-6** — first-class-departure precedent
   (033 `NATIVE_P`, stage003 stray-longitudinal) + RESUME §5 "record departures first-class, no softening." Alternatives:
   fold into V-5 (matches the build's "Hooks" placement) or into V-2 (where the parity is decided).

**My overall recommendation: keep all six** — each of {action row, native source, Route A reference, Route B+comparison,
SEALED landing, T-even departure} is a distinct, separately-citable result, and the two-independent-routes boost test is the
sector's whole point.

---

## Tooth allocation (build-faithful — from the fresh `check.{py,wl}` audit, 35 teeth in `TOOTH_ORDER`)

The ONE source build has **35 teeth** (identical `TOOTH_ORDER` in `.py` and `toothOrder` in `.wl`). Each ledger stage's
focused dual-engine script re-derives its own cluster (the source build is the provenance, re-authored per the mirror policy).
`TARGET_BLINDNESS`, `DUAL_ENGINE_TERMS`, and a `UNITS_RESTORED`-equivalent dimensional check are **build-global** — every
stage re-instantiates them (matches the standing per-stage dimensional-firewall rule).

| Stage | Teeth (source-build names) | # |
|---|---|---|
| **034 V-1** action row | `ACTION_KINETIC`, `ACTION_COUPLING`, `ACTION_STABILITY`, `G0_DAMAGE`, `LEDGER_READY_ROW`, `FIELD_IDENTITY_UNITS` | 6 |
| **035 V-2** source + parity census | `SOURCE_TRANSLATION_CONTINUITY`, `SOURCE_NOT_IMPORTED`, `SOURCE_BASIS`, `PARITY_RW`, `PARITY_PW`, `PARITY_ROTATION`, `PARITY_TIME_REVERSAL` | 7 |
| **036 V-3** Route A reference | `BOOST_PROJECTOR`, `BOOST_GENERAL_VELOCITIES`, `BOOST_NEXT_ORDER` | 3 |
| **037 V-4** Route B (blind) + comparison + ratios | `DIRECT_SOURCE`, `DIRECT_PROJECTOR`, `DIRECT_EXCHANGE_SIGN`, `DIRECT_FALLOFF`, `DIRECT_VELOCITY_ORDER`, `ROUTE_INDEPENDENCE`, `BOOST_COMMON_VELOCITY`, `COMPARE_COMPUTED`, `DELTA_RATIO`, `CONE_RATIO`, `QMAG_R1` | 11 |
| **038 V-5** SEALED §4 landing | `TRUTH_TOTALITY`, `TRUTH_PRECEDENCE`, `LANDING_OWNERSHIP`, `ACTIVE_FLUX_CAVEAT`, `HOOK_LORENTZ` | 5 |
| **039 V-6** `b_T` T-even departure | *(cites V-2's `PARITY_TIME_REVERSAL` + `PARITY_ROTATION`; authors its own `b_T`-axial + T-even vs Maxwell-T-odd asserts)* | census-reuse |
| **build-global** (each stage) | `TARGET_BLINDNESS`, `DUAL_ENGINE_TERMS`, `UNITS_RESTORED` | 3 |

**Count check:** 6+7+3+11+5 + 3 global = **35** ✓ (V-6 reuses V-2's parity teeth, adds authored departure asserts).

**Cross-tooth dependencies flagged by the audit (inform, don't block — ledger stages cite prior structure):** the 4 `PARITY`
teeth are one census (kept whole in V-2); `COMPARE_COMPUTED`/`DELTA_RATIO`/`CONE_RATIO`/`BOOST_COMMON_VELOCITY` reuse **both**
routes (→ all land in V-4); `LANDING_OWNERSHIP` reuses `LIVE_FACTS` spanning source+comparison+Q-MAG (→ V-5 cites V-2/V-4);
the `b_T` departure reuses `PARITY_TIME_REVERSAL`+`PARITY_ROTATION` (→ V-6 cites V-2).

---

## What Part V EXCLUDES (surviving-solution compliance)

- **`pathA_39` `j∝sV` scope** — the barred imported current; superseded by the native `J_T=q_T sηV`. → failures-paper backlog
  (do NOT fold). The blueprint's original Part V (`pathA_39`) scope is retired.
- The "hermetic v3" over-engineered magnetism directive that was tried + `git reset`-reverted (history only).

---

## Parameters Part V adds / re-homes (register preview — the per-stage updates fill these; continue edges after R66)

Dims **pinned from the build's Appendix A** units table (no longer "to confirm"):

| Param | dim | Enters | Class (proposed) | Note |
|---|---|---|---|---|
| `u_T` (transverse shear vector) | `L` | V-1 (034) — consumed | consumed (stage003) | the light-sector field; `∇·u_T=0`; NOT re-counted |
| `μ_R`, `c_γ` | `μ_R: ML⁻¹T⁻²` (`=EL⁻³`) | V-1 consumed | consumed | `c_γ²=μ_R/ρ_br`; shared with light; NOT re-counted |
| `q_T = λ_T τ_d` (throat charge) | `MT⁻¹` | V-1/V-2 | `FREE-UNREDUCED` / R1 (throat) | rides the sim-deferred throat solve (sibling of the electric `A_E`, R10/R30); parity even/even, **T-odd** |
| `τ_d` (active-drain time-arrow) | structural (T-odd arrow) | V-1/V-6 | structural / postulated | the time-arrow magnetism requires; the full time-reverse maps drain→source; feeds the `b_T` T-even departure |
| `J_T = q_T sηV` (native source) | `ML⁻²T⁻²` | V-2 (035) | `DERIVED` (defect-continuity), value cond. on `q_T` | REPLACES `j∝sV`; `s`=±w orientation (charge sector), `η`=mouth measure `[L⁻³]`; landing `CONVECTION_LIKE_CONDITIONAL` |
| Maxwell–Darwin kernel `(δ_ij+n_in_j)/8πR` | `L⁻¹` | V-3 (036) | `DERIVED` (boost of electric) | the reference kernel to `O(v²/c_γ²)`; `A_E: EL` |
| `r_BA = q_T²/(ρ_br·A_E)` (coefficient ratio) | `1` | V-4 (037) computed / V-5 landing | **R1-deferred** (`electric_bc_selection`) | the sole open coefficient; rides R1 `q_T` + R1 electric `A_E`; `δ_BA=r_BA−1` |
| `r_cone = c_E²/c_γ²` (`c_E=c_γ?`) | `1` | V-4 computed / V-5 (038) | **R1/OPEN** | no committed cone lock (Part VI `pathA_40` re-adjudicates) |
| `b_T = ∇×u_T` (candidate mag. field) | `1` (dimensionless) | V-6 (039) | **DEPARTURE** (T-even, axial; discharges no knob) | first-class characterized departure, NOT a reduction — must not shrink the irreducible count |

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

## Decisions for the per-Part gate (fresh read DONE — these are the live calls)
1. **Stage count 5 vs 6** — the three granularity calls (034+035 merge? / 036+037 merge? / 039 standalone?) are itemized in
   the "Gate decisions — granularity" block above. **Recommend all six.**
2. **Build-id start = 034** (Part IV ended at 033) — confirm. Verdict tokens are now build-faithful (V-2 `CONVECTION_LIKE_CONDITIONAL`,
   V-6 `B_TIME_REVERSAL_EVEN`) — confirm.
3. **Tooth allocation** — confirmed against the 35-tooth `TOOTH_ORDER`; the per-stage clusters + build-global set are tabled
   above. (Build-time detail; the gate just ratifies the boundaries.)
4. **Then:** ratify → per-stage pre-execution gate on stage 034 → Codex→Grok→Codex directive → build. Await the user's
   autonomy mode for the per-stage builds.

---

## Progress log
- **2026-07-22 — split plan RATIFIED (6 stages 034–039)** at the per-Part gate; all six kept separate (the three granularity
  calls resolved to KEEP SEPARATE); verdict tokens/stage boundaries/register dims/35-tooth allocation build-faithful.
- **034 (V-1, transverse-vector action row — EARNED action row)** — ✅ **DONE + committed `109070da`** (verdict
  **`TRANSVERSE_MOVE_ACTION_ROW`**, `internal_inconsistency=none`). Records, as an EARNED action-structure result within the
  postulated G0 closure, that magnetism = the MOVING ±w throat enters the ledger as a clean `(G0+δ)` transverse-vector
  amendment: the row `S_{T+move}=∫[½ρ_br|u̇_T|²−½μ_R|∇×u_T|²+q_T Σ sᵢ η_a Vᵢ·u_T]`, `∇·u_T=0`, `c_γ²=μ_R/ρ_br`,
  `q_T=λ_T τ_d`; **positive-definite transverse Hessian** (two polarizations, `ω²=c_γ²k²`, no ghost/tachyon) and **no
  pre-existing G0 row changed** (`F_flux` untouched). **Scope class EARNED (action row) / IMPORTS:** the `u_T` kinetic/gradient
  row + `c_γ²=μ_R/ρ_br` are IMPORTED from Part III / stage003 / pathA_36 — **cited, NOT re-counted** (de-dup deferred to Part
  VII); the magnetism-NEW content of 034 is EXACTLY `{the moving coupling, q_T, τ_d}`. **Action-row gate ONLY:** the sign +
  magnitude are deferred to the R1 landing (038), the `b_T` parity/departure to 039; source-from-continuity + parity census
  (035) and both far-field routes (036–037) are downstream. **Verification:** dual-engine SymPy 12 / Mathematica 12, both
  exit 0 (SymPy explicit-Hessian route; the `.wl` a genuinely INDEPENDENT NullSpace-Fourier + `CoefficientArrays` +
  `PositiveDefiniteMatrixQ` route); all 12 per-tooth mutations FIRED_AT_OWN_ASSERT in each engine; the verdict tooth is
  non-tautological (re-derives to `ROW_UNSTABLE` under a mutated Hessian, the 030 X≡X lesson). Directive Codex→Grok→Codex
  bookend clean; build + arbiter re-run; **tri-review (falsification-first): FIDELITY 1 NIT** (a `.wl` X≡X `dimensionObject`
  self-compare) **REMEDIATED → clean**, **ADVERSARIAL clean + 2 documented non-blocking notes**; Codex deliverables-fidelity-verify
  **SOUND**. Deliverables note + card (Part-V appendix input) + register edges **R67** (`q_T=λ_T τ_d`, FREE-UNREDUCED / R1,
  split from the electric `bc_selection` — the sim-deferred throat normalization) + **R68** (`τ_d`, the active-drain time-arrow,
  structural / postulated). ▶ **NEXT = stage 035 (V-2)** — the native source law `J_T=q_T sηV` from defect-continuity + the
  parity census, verdict `CONVECTION_LIKE_CONDITIONAL`.
