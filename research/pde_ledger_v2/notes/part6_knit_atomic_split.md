# Part VI — The knit: atomic-stage split (DRAFT — per-Part gate + fresh read PENDING)

> **⚠ STATUS: DRAFT proposal, NOT ratified.** Authored 2026-07-23 right after Part V (magnetism) COMPLETE, from the
> `ledger_v2_blueprint.md` §156 Part-VI row + `RESUME_ROADMAP.md` §3.6/§4/§5 + `docs/model_map.md` §3.6 + the **fresh Part-V
> results** (the cone question `r_cone` is now OPEN from stages 037/038). **The per-Part user gate + a FRESH read of the three
> source reports** (`software/stage1_solver/reports/pathA_40_cone_lock.md`, `…pathA_41_ng5_second_medium_drift.md`,
> `…pathA_42_charge_coupled_scalar.md`) **+ their scripts** will finalize the exact stage boundaries, verdict tokens, register
> rows, and the reshape/re-authoring cost. Treat the table below as a starting proposal (exactly as `part5_magnetism_atomic_split.md`
> was drafted-then-refined-then-ratified). Do NOT build any Part VI stage until this plan is ratified at the gate.
>
> **⭐ KEY subtlety — Part VI is a RE-ADJUDICATION, not just a reshape.** Parts IV/V re-scoped onto NEW builds
> (puncture-deflection, magnetism-moving-throat). Part VI carries the OLDER `pathA_40/41/42` reports, whose STATUS has SHIFTED
> given the newer G0 card + the just-completed Part-V build. In particular the **cone-lock (pathA_40) is now UNSUPPORTED**: the
> G0 card + the magnetism build do NOT establish `λ_γ=1` or `c_E=c_γ` — there is **no committed cone lock** (`c_s/c_E=2` was
> *chosen* to avoid forcing one), and Part V's stages 037/038 left **`c_E=c_γ?` explicitly OPEN via the ratio `r_cone=c_E²/c_γ²`**.
> So the fresh read must RE-ADJUDICATE pathA_40 honestly (the locks are calibration CHOICES that are NOT currently established),
> not merely transcribe its old `CONE_LOCK_CALIBRATED` token. **⭐ Authored under the "ledger = surviving solution only" rule** —
> any retired knit approaches → failures-paper backlog, not the ledger.
>
> Governing: `notes/ledger_v2_blueprint.md` §156 (Part VI row) + §5 (reshape spec) + §6 (verification); `RESUME_ROADMAP.md`
> §3.6 (the three gates) + §4 (build order) + §5 (cross-cutting cones/throat-solve facts); `docs/model_map.md` §3.6 + §4
> (departure ledger); `notes/MATHEMATICA_MIRROR_POLICY.md`.

---

## What the knit is (in one breath)

Part VI **VERIFIES "one medium, all emergent" instead of asserting it** — it asks the three questions that could still break the
single-medium claim now that all four far-field sectors exist on one shared field set:
1. **Shared cones?** (pathA_40) — are the light/gravity/electric cones locked, or are the locks unearned calibration choices?
2. **One irreducible parameter set with no smuggled second substance?** (pathA_41) — the decisive reducibility test.
3. **Does the propagating scalar stay consistent?** (pathA_42) — the charge-coupled `h`-scalar's break-risk.

This is where the knit's **real falsification power** lives (RESUME §3.6): the cone-lock is near-vacuous (it's just naming
calibration inputs), but **NG5 (pathA_41) is the decisive reducibility test** — it's the one that could expose a smuggled
second substance. None of the three is expected to cleanly "all work"; the honest outputs are a **no committed cone lock**, a
**characterized incomplete-reduction departure** (irreducible `{ρ_B0, χ_c, C_hu}`, but `no_fourth_arena=True`), and a
**sim-gated scalar departure**.

---

## Proposed atomic-stage split (~3–5 stages, build ids 040+ — the gate may split the cone-lock or add a knit-synthesis capstone)

| Build id | Stage | Headline verdict token (from blueprint §156 — RE-ADJUDICATE at fresh read) | Scope class |
|---|---|---|---|
| **040** | **VI-1** cone-lock re-adjudication (`λ_γ=1`, `c_E=c_γ`) | `CONE_LOCK_CALIBRATED` → **RE-ADJUDICATE to the honest "NO committed lock"** (`r_cone` open, `Δr=2`) | **CALIBRATED / UNCOMMITTED** (both locks are inputs, not earned) |
| **041** | **VI-2** NG5 second-medium reducibility | `SECOND_MEDIUM_DRIFT(active_irreducible={ρ_B0,χ_c,C_hu})` | **CHARACTERIZED-DEPARTURE** (incomplete 4D→3D reduction; `no_fourth_arena=True`) |
| **042** | **VI-3** charge-coupled scalar consistency | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | **CHARACTERIZED-DEPARTURE + SIM-DEFERRED** (magnitude) |

**Per-stage scope sketch (to be confirmed against the reports at the gate):**
- **040 (VI-1 — cone-lock):** the two calibrated locks — **(A)** `λ_γ=1 ⟺ μ_R/ρ_br = 5Kρ⁴/m` (light cone = gravity phonon
  cone) and **(B)** `c_E=c_γ ⟺ c_E²ρ_br = μ_R` (electric cone = light cone). RE-ADJUDICATE: neither is derived on the earned
  ledger — both are calibration INPUTS with codimension `Δr=2`; ⚠ **the G0 card + magnetism build do NOT establish them → NO
  committed cone lock**; `c_s/c_E=2` was chosen to avoid forcing one; **Part V's `r_cone=c_E²/c_γ²` (stages 037/038) left
  `c_E=c_γ?` OPEN**. The honest verdict is that the cone locks are *available calibration choices*, NOT earned facts —
  near-vacuous as a falsifier but load-bearing for Part VII's codimension count. `λ_γ=c_γ/c_s` is a DERIVED ratio, NOT a free
  knob (the free content is `c_γ`, ultimately `{μ_R, ρ_br}`).
- **041 (VI-2 — NG5, the decisive one):** the one-medium claim doesn't fully close — `{ρ_B0, χ_c, C_hu}` remain **irreducibly
  independent** (unreduced brane-surface/embedding parameters), while `{ρ_br, μ_R, c_E}` are reducible-in-principle via the
  shared throat solve. Crucially **NOT a separate substance** — `no_fourth_arena=True` (the irreducible trio lives on the ONE
  medium's brane/embedding, not a smuggled second arena). ⚠ **Anti-absorption guard:** the NG5 trio `{ρ_B0,χ_c,C_hu}` is a
  DISTINCT drift concept from the `stage007` freeze drift (operative `POST_D16_DRIFT(7)`) — the two must never be conflated /
  the trio must never be silently absorbed into the freeze count (RESUME §3.6 note). This stage carries the knit's real
  falsification content.
- **042 (VI-3 — charge-coupled scalar):** the propagating `h`-scalar doesn't cleanly break the model, but its break-risk
  magnitude is SIM-GATED: `h_EP` earned-safe on the decoupled floor; `radiation` / `universality` / `u_L_EP` / `preferred-frame`
  all `SIM_GATED`. The sharp trade-off: **`c_E→∞`** → non-radiating but preferred-frame; **`c_E=c_γ`** → Lorentz-invariant but
  a real radiating extra scalar. Records the mapped departure + the sim-gate; does NOT resolve it (that's the deferred sim).

**⚠ Gate decisions — granularity:**
1. **Cone-lock split?** 040 could split into two stages (`λ_γ=1` and `c_E=c_γ` separately) since they lock different cone
   pairs with different downstream stakes. My recommendation: **keep as one re-adjudication stage** (they are one "cone-lock"
   question), unless the fresh read shows the two locks need independent teeth/verdicts.
2. **Knit-synthesis capstone?** Consider a small **VI-4** that COLLECTS the three gate outputs into the honest one-medium
   verdict ("one medium, all emergent — with: no committed cone lock, an irreducible-but-no-fourth-arena embedding trio, and a
   sim-gated scalar departure") feeding Part VII. Recommendation: **defer to Part VII** (this synthesis IS Part VII's job) —
   but flag at the gate.
3. **Reshape/re-authoring cost:** the fresh read must confirm whether pathA_40/41/42 have existing scripts and whether their
   `.wl` are mirrors needing independent re-authoring (per `MATHEMATICA_MIRROR_POLICY` + the RESHAPE COST) — these are OLDER
   `stage1_solver` artifacts, likely `argparse --compare`/payload-mirror style like pathA_38/39 were.

---

## How Part V (magnetism) feeds Part VI (carry these into the fresh read)

- **The cone question is now OPEN from Part V.** Stages 037/038 computed `r_cone=c_E²/c_γ²` and left `c_E=c_γ?` unresolved
  (`HOOK_LORENTZ` UNDETERMINED needs `r_cone=1` AND `δ_BA=0`). So pathA_40's cone-lock re-adjudication must state that Part V
  did NOT close the cone — it left it as an open ratio. This is the concrete evidence the locks are uncommitted.
- **`c_γ` is shared by light + magnetism** (`c_γ²=μ_R/ρ_br`); `c_E` is the electric cone with no committed lock. The knit's
  cone question is exactly whether `c_E` locks to `c_γ`.
- **The shared THROAT SOLVE is the reduction lever.** The nonlinear throat solve is a shared R1 across gravity
  (`{μ_R,ρ_br}`/R10/R30), electric (`bc_selection`/sign, R63), and magnetism (`q_T`, R67) — Part V made the magnetism `q_T`
  R1 and the doubly-R1 landing (038) explicit. NG5's `{ρ_br,μ_R,c_E}` are "reducible-in-principle via the throat solve" — the
  SAME solve. Part VI states the reduction debt; Part VII's throat solve would discharge it.
- **The Part V departures** (`B_TIME_REVERSAL_EVEN`, R73) join the model_map §4 departure ledger alongside the light
  stray-longitudinal and the charge `NATIVE_P_NO_EMERGENT_GAUSS` — Part VI's `SECOND_MEDIUM_DRIFT` + `SCALAR_DEPARTURE_MAPPED`
  are two more first-class characterized departures.

---

## What Part VI EXCLUDES (surviving-solution compliance)

- Any retired knit approach / superseded cone-lock framing → failures-paper backlog (`notes/ledger_exclusions_failures_paper_backlog.md`),
  NOT the ledger. Show only the surviving re-adjudicated result.
- The actual THROAT SOLVE that would discharge the reduction debt is **SIM-DEFERRED** (Part VII records it as the central
  reduction debt) — Part VI names it, does not perform it.
- The sim that would resolve pathA_42's `SIM_GATED` scalar break-risk — named, not performed.

---

## Parameters Part VI adds / re-homes (register preview — continue edges after R73)

| Param / fact | Enters | Class (proposed) | Note |
|---|---|---|---|
| cone locks `λ_γ=1`, `c_E=c_γ` | VI-1 (040) | **CALIBRATED / UNCOMMITTED** (`Δr=2`) | NOT earned; the G0/magnetism builds don't establish them; `r_cone` (R71, Part V) is the open handle |
| NG5 irreducible trio `{ρ_B0, χ_c, C_hu}` | VI-2 (041) | **IRREDUCIBLE (unreduced embedding)**; `no_fourth_arena=True` | the decisive one-medium residue; anti-absorption-guarded vs the stage007 freeze drift; feeds Part VII codimension |
| `{ρ_br, μ_R, c_E}` reducibility | VI-2 (041) | reducible-in-principle (throat solve) | the SAME shared throat solve as gravity/electric/magnetism R1s |
| charge-coupled scalar map + sim-gate | VI-3 (042) | **DEPARTURE + SIM-DEFERRED** (magnitude) | `h_EP` earned-safe; radiation/universality/`u_L_EP`/preferred-frame SIM_GATED; the `c_E→∞` vs `c_E=c_γ` trade-off |

These are mostly **CALIBRATED / DEPARTURE / IRREDUCIBLE** rows — Part VI adds little NEW knob content; its job is to
CHARACTERIZE the knit honestly and hand Part VII a clean irreducible-count + reduction-debt statement. (Discharges no new
reduction edge; the NG5 trio ADDS to the irreducible count, it does not shrink it.)

---

## Cross-Part dependencies

- **From Part V:** `r_cone` (R71), the shared `c_γ`, the `q_T`/`A_E` R1s, the doubly-R1 magnetism landing — all feed VI-1/VI-2.
- **From Parts I–III:** `c_s²=5Kρ⁴/m` (gravity), `c_γ²=μ_R/ρ_br` (light), the `stage007` freeze drift (`POST_D16_DRIFT(7)` —
  kept DISTINCT from the NG5 trio), the medium field set.
- **To Part VII (integration):** the cone-lock adjudication (or lack thereof), the irreducible `{ρ_B0,χ_c,C_hu}` trio, and the
  scalar consistency all feed Part VII's **unified equation set + calibration map + the irreducible codimension count** (via
  the `pathA_40` `Δr=2` technique from the midway knob audit) + the **shared throat solve** as the central reduction debt. The
  2 parked knob-audit decisions (RESUME §2 — C1/C2 counting convention + R35 label) resolve **before/at Part VII**, informed
  by VI-2's irreducible count.

---

## Per-stage process (unchanged — blueprint §5/§6)

Per stage: Codex→Grok→Codex directive bookend (fold via agents) → per-stage pre-execution user gate → Codex build
(`--sandbox danger-full-access`, detached `setsid` + poll) → dual-engine both exit 0, independent `.wl` → orchestrator arbiter
re-run → fresh-agent tri-review (fidelity + adversarial) → remediate (real coverage gaps only; document verified-safe smells)
→ register update + Codex-verify → self-contained note + TeX card + registration → deliverables fidelity-verify → commit +
tracker/memory sync (separate commit, real hash). **Carry the Part-V process wins:** the ⛔ NO-CAN'T-FAIL-CONJUNCTS build-prompt
ban (killed the 037-039 NITs); the emphatic anti-give-up rule on every bookend/verify gauntlet (sub-agents are one-shot, must
busy-poll — see `feedback_subagent_marathon_infra`); Grok may be usage-limited (substitute a fresh independent compute agent).

## Open questions for the per-Part gate
1. **Stage count 3 vs 4–5** — cone-lock split (recommend keep as 1)? knit-synthesis capstone (recommend defer to VII)? Confirm
   after the fresh read.
2. **Build-id start = 040** (Part V ended at 039) — confirm.
3. **Verdict tokens** — RE-ADJUDICATE against the current post-Part-V state (esp. pathA_40 → "no committed cone lock", NOT the
   old `CONE_LOCK_CALIBRATED` as-if-settled); confirm the NG5 + scalar tokens against the actual reports at the fresh read.
4. **Reshape cost** — do pathA_40/41/42 have existing dual-engine scripts, and are the `.wl` mirrors needing re-authoring?
   (Fresh-read task; likely yes, like pathA_38/39.)
