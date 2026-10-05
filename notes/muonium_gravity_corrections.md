# Gravity-sector corrections + the muonium falsification track — what we need

**What this is:** the distilled, actionable content we need to carry forward from the muonium/gravity handoff —
the corrections to make to the *old* gravity/lepton program, and the forward-warning we plug into the ledger so
the gravity sector doesn't corner itself. **Status:** research-track + corrections, **not** a claimed prediction.
**Provenance:** the full handoff (`notes/muonium_gravity_research_track_handoff.md`) was written by Codex after
reading the ETH/PSI muonium-beam article; this doc is the orchestrator's distillation + review (2026-09-15).

**The trigger:** ETH Zurich / PSI have a cold, directional muonium (`Mu = μ⁺e⁻`) beam (LEMING); an
atom-interferometry free-fall measurement of a *second-generation* lepton is ~2–3 years out. That opens a real
falsification window for the passive gravitational response of the muon sector.

---

## 1. The methodology corrections to adopt (these are just our own discipline, in the gravity sector)

1. **⛔ No circularity.** Do **not** use assumed universality (κ=1) to infer throat geometry and then "derive"
   universality from that geometry. Solve the electron/muon branches from **nongravitational** constraints, then
   compute the gravitational response as an **output**. (This is M3 — prior art is an oracle, never a premise —
   and the target-blind / withhold-the-acceptance-criterion rule, applied to gravity.)
2. **Three-mass ledger — keep them separate.** Inertial `m_i`, **passive** `m_p` (response: how it falls — what
   LEMING measures), **active** `m_a` (source: the field it produces). Define each independently; report
   `κ^pass = m_p/m_i` and `κ^act = m_a/m_i` separately. Equality `m_i=m_p=m_a` must **emerge as a Ward-identity
   theorem** if it holds — ⛔ never assumed. (This is "keep every varying quantity live; a required freeze is the
   finding.")
3. **κ_ρ=1 is TARGET-MATCHED, not derived.** In the 1PN worldline ansatz, κ_ρ=1 is fixed by *matching the
   reduced action to Newton* — it proves the chosen effective branch has ordinary response **by construction**,
   ⛔ not that a first-principles muon branch must reproduce it. Re-label it as matched everywhere it's asserted;
   leave the species-specific κ_{ρ,s} **open** and compute it per branch.
4. **Blinding / no tuning.** Freeze the prediction (central value + interval) **before** the measurement; common
   parent-medium constants for e and μ; ⛔ no muon-specific continuous knob introduced to alter gravity; a
   post-experiment branch choice is not a prediction. (falsification-is-the-goal.)

## 2. The reset — withdrawn vs surviving

**⛔ WITHDRAWN as physical results (do-not-import):** every absolute/relative lepton throat radius, diameter,
depth, volume, aspect ratio — including `L/a≈1.85` as real geometry, any lepton radius inferred from rest mass,
"heavier ⇒ narrower/deeper", `a_j∝(2j+1)^-1` for real families, reverse-engineered deep-needle heavy-lepton
geometries, and `m_G~κ_m ρ₀ π a² L` used *quantitatively*. Any gravity conclusion resting only on these falls
with them.

**Survives:** geometry-independent identities (field identities, projection/conservation laws). **Conditional
(rerun after a real branch exists):** the `F(a,ρ)=A/a+B/a²+Ca³` functional + its virial `E_w+2E_f=3E_PV` (valid
*for that functional*, not validated as the real branch); the 11:2:5 partition; the `−57/64` breathing slope;
the absolute mass formula; D/N family scaling. **Far-field 0PN–4PN match:** survives as *"a chosen universal
effective branch reproduces the dynamics within its declared closure"* — ⛔ NOT as proof that every microscopic
species reduces to that branch. The missing task is **defect-to-worldline matching per species**.

This reset is **consistent with what the program already recognized** — the lepton tower was flagged falsified
and the throat deep-solve is off the critical path. It formalizes + extends that, it doesn't contradict it.

## 3. Grounding finding — the reset does NOT threaten the current S11 build

A targeted search (2026-09-15, **corrected** — an earlier pass wrongly said "not in `pde_ledger_v3`") shows the
withdrawn **numerical** results (`L/a`, `1.85`, `11:2:5`, `−57/64`, `2j+1`, deep-needle, `m_G~ρa²L`) are **NOT
imported into the S11 BUILD** (`research/pde_ledger_v3/scripts/*` + `*_exports.py` compute none of them) — but they
DO appear in the ledger tree as **history/register entries** (`research/pde_ledger_v3/DEFECT_REGISTER.md`,
`SESSION_REASONING.md`) and in the S16 plug, alongside the old `notes/`+`docs/` track. ⇒ **the reset is bounded
and non-blocking for the in-flight S11c-d BUILD** (nothing computed imports them), but the quarantine must cover
those ledger register/history docs too. The junction where these results would re-enter the ledger is the gravity-sector
worldtube/interior-dependence step — **S16 in `V3_STEP_PLAN.md`** — which is where the forward-warning plug is
placed (§5 below). *(Caveat: this is a targeted spot-check, not the full Phase-A dependency audit — that audit is
what makes it definitive.)*

## 4. Near-term actions (bounded, independent of the S11 hold)

1. **Dependency audit + quarantine of the OLD track.** Grep `notes/`+`docs/` **and the ledger register/history
   (`research/pde_ledger_v3/DEFECT_REGISTER.md`, `SESSION_REASONING.md`)** for `a`, `L`, `L/a`, radius,
   diameter, volume, `m_G`, `κ_ρ`, `11:2:5`, `2j+1`, `−57/64`; classify each hit
   (exact / ansatz-dependent / fitted / target-matched / reverse-engineered / obsolete); add explicit
   `do-not-import` warnings so invalid geometry isn't silently reused. Ground it in commands→`_measurements`,
   ⛔ not prose (rule 2/E1).
2. **⭐ The one load-bearing check:** verify the far-field PN ladder match is **genuinely geometry-independent**
   (doesn't secretly use withdrawn lepton geometry). If it does, the reset is bigger than the handoff claims.
   The PN work is in `notes/` (`moving_throat_*`, `4d_*pn_*`, `pathA_*`). Highest priority.
3. Re-label κ_ρ=1 as *target-matched* wherever the PN docs assert it.
4. **⭐ Structural consistency audit BEFORE any muon branch solve.** At **one fixed approximation order**, answer:
   *can a defect source gravity anomalously (`m_a ≠ m_i`), fall normally (`m_p = m_i`), and conserve total
   momentum?* This is the concrete crux of the notes tension: the 1PN reduction assumes **negligible worldtube
   boundary flux** (`notes/summaries/4d_1pn_full_summary.md:453` → normal fall) while the bridge retains an
   explicit momentum-exchange term `S_{J_i}` (`notes/summaries/4d_1pn_bridge_summary.md:701,710` → the exchange
   that could enable an active anomaly). Normal-fall and anomalous-attraction currently come from **different
   regimes**; they can coexist only if — *at the same order* — the required exchange is shown to **preserve**
   normal fall (then exhibit the compensating momentum flow) or shown **not** to (then the active anomaly must
   surface elsewhere). Do this **before** a full branch solve, ⛔ with no preference for whichever answer makes
   LEMING more useful.

## 5. The ledger plug (where we revisit this)

A forward-warning is plugged into **`V3_STEP_PLAN.md` at S16 (the worldtube reduction)** — the load-bearing
step that already flags (a) gravity depends on the interior beyond leading order, (b) the result is
**response-side (passive)** not source-side (active), and (c) it assumes a **calibrated `−Gm/r`** and supplied
mass. That is *exactly* where a species anomaly (`κ_μ≠1`) would live and where universality is currently
assumed. The plug says: before extending S16 to the interior-/species-dependent response, revisit this doc —
don't assume κ=1 to fix geometry, use the three-mass ledger, don't import the withdrawn throat geometry.

## 6. Dependency ordering + honest expectation

A genuine, blinded `κ_Mu` prediction is gated on machinery that **doesn't exist yet**: the moving-throat PDE +
throat interior (off critical path) **plus** the S11 medium chain (in flight). So this is a **deferred parallel
branch** — begin definitions/audit now, numerical prediction later. Realistically it may **not** be ready by
LEMING's 2–3 yr horizon, in which case the honest result is the handoff's **Outcome 5**: the model reproduces
the GR-like far field under a declared closure but has **not** derived equivalence from its particle ontology —
i.e. no muonium prediction. That's a legitimate outcome, not a failure to hide.

**Outcomes to keep in view:** (1) universality emerges from the real branches → strengthens the model
(consistency test); (2) small parameter-free muon deviation frozen before the measurement → genuine novel
prediction; (3) deviation already excluded by spectroscopy/kinematics/astrophysics → falsified pre-LEMING;
(4) no acceptable muon branch with common constants → the lepton construction fails (valuable falsification,
⛔ don't rescue with a species knob); (5) can't compute without importing κ=1 → no prediction (honest).

## 7. What counts as a failure vs a prediction — the passive/active guardrail

A muon throat whose geometry / `κ_μ` **surprises us is not a failure** — the geometry and `κ_μ` are *outputs* of
the branch solve, not targets, and we do **not** lock the throat radius for second-generation matter. But keep
the STATUS of a mismatch precise (⛔ "not a failure" does not auto-promote to "a win"):

- **Mismatch visible in passive (`m_p/m_i ≠ 1`)** → LEMING can measure it → a genuine **falsifiable prediction**.
  Best case.
- **Active-only mismatch (`m_a ≠ m_i`, `m_p = m_i`)** → not measurable via the muon's own fall → **not a failure,
  but not (yet) a usable prediction** — unfalsifiable *unless* the §4.4 audit shows the enabling exchange has
  other observable consequences (recoil, other species). ⛔ Active-only ≠ automatically unfalsifiable; that needs
  the calculation.
- **Mismatch that violates momentum conservation / internal consistency** → a **red flag** (error or genuine
  failure), *unless* the model **derives** the reservoir flux that absorbs it. ⛔ A leakage reservoir *permits*
  exchange; it does **not** *guarantee* the required directional flux — that must be shown, not assumed.

**The freedom is asymmetric by generation.** The **electron / first-generation** side is NOT free — the model
must reproduce ordinary-matter universality (`κ_e ≈ 1` to ~10⁻¹⁵); if it can't, *that* is a failure. The **muon**
`κ` is the output allowed to differ. **Genuine failure modes** (so the freedom has edges): can't reproduce the
muon's nongravitational observables (rest mass, charge, moment, decay) with the *same* parent constants as the
electron; a muon-specific tuning knob introduced to alter gravity; reverse-engineering the geometry from the mass.

**The discipline is symmetric:** ⛔ don't privilege `κ=1`, and ⛔ don't privilege `κ≠1` because it's more
interesting. Solve it, report what comes out, and let the mismatch's status (measurable / unfalsifiable /
inconsistent) be decided by the §4.4 audit — never by which answer we'd prefer.

## Related
Full source: [`muonium_gravity_research_track_handoff.md`](muonium_gravity_research_track_handoff.md). Ledger
plug: `research/pde_ledger_v3/V3_STEP_PLAN.md` §S16. PN/throat work: `notes/moving_throat_*`, `notes/4d_*pn_*`,
`notes/pathA_*`. Memory: the lepton-tower-falsified / throat-off-critical-path state; the dark-energy postulate
(banked in V3_STEP_PLAN) is where a *controlled* leak becomes a feature. Governing: M3 (oracle not premise),
target-blind calibrate-predict, falsification-is-the-goal.
