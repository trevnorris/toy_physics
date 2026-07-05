# The consistency knit — handoff / post-`/compact` resume note (2026-07-04)

> **⏸ UPDATED 2026-07-04c — sub-pieces (a) + (b) are DONE; this note below is the ORIGINAL scoping (kept as history).** Current state +
> the authoritative resume pointer = `STATUS.md` ▶ RESUME HERE and `decisions/13` §0. In brief: **(a) cone-lock = `CONE_LOCK_CALIBRATED`
> (`pathA_40`, `9bb58a97`, tri-reviewed EARNED). (b) NG5 = `SECOND_MEDIUM_DRIFT{ρ_B0,χ_c,C_hu}` (`pathA_41`, `453e342f`) — interpreted as
> ONE candidate medium whose 4D→3D reduction is INCOMPLETE (NOT two substances; no fourth arena); a build-rig was adversarial-caught + a
> Claude↔Codex framing re-exam corrected the routing (`ρ_br/μ_R`→pending Route-A; GLM "overcount"=red herring, pathA_35 erratum corrected).
> ✅ re-tri-review of the remediated `pathA_41` = DONE + CLEAN (2026-07-05: fidelity `FIDELITY_CLEAN` + adversarial `ADVERSARIAL_CLEAN`
> able-to-PASS + hardening fold `HARDENING_CLEAN` — caught a reproducibility bug + made location-closure able-to-fail); pathA_41 BANKED EARNED.** ▶ **NEXT = task #110** (the charge-coupled extra scalar:
> does it clearly BREAK the model or is it naturally hidden?), THEN **(c)** the `pde_ledger` assembly. Detail = [[project-brane-existence-defect-structure]]
> UPDATE 2026-07-04c + `_scratch/pathA_41_framing_codex.md`.

**Context:** all 4 force-sectors are now EARNED as independent specs (gravity `pathA_29`+PN, light `pathA_36`, charge
`pathA_38`, magnetism `pathA_39` complete). This note is the running start for the NEXT front — the **consistency knit**:
do the four sectors live in **ONE self-consistent parameter set**, and can they be assembled into the central
calibrated PDE? This is prep only; the knit is executed post-`/compact` via the usual gauntlet (Claude↔Codex setup →
directive → Codex design-review → GLM → dual-engine → tri-review). **Front door = `STATUS.md` ▶ RESUME HERE.**

---

## §1. What the knit IS — three sub-pieces
1. **The speed cone-lock (`λγ = c_γ/c_s = 1` AND `c_E = c_γ`) — the crux.** Each sector fixed a speed:
   - `c_γ² = μ_R/ρ_br` (light; pathA_36).
   - `c_s` = the speed gravitational changes propagate (`∝ ρ²`; the gravity/drain sector — NOT a ripple of the same kind as
     `c_γ`; see [[project-model-mechanics-corrections]] — three speeds: `c_s` gravity-change, `v_r` field strength (not a
     speed), `c_γ` light).
   - `c_E` = the charge/EM dynamic-Green speed (pathA_38, from `exp(iRω/c_E)/(4πR)`).
   **The question:** are these forced to coincide by the shared medium parameters (`λγ=1`, `c_E=c_γ` DERIVED), or is the
   cone-lock a **calibration input / a genuine gap**? History flags `λγ` as `BETA_GENUINE_GAP` (memory
   [[project-pathA-build]]); `c_E=c_γ` is the retired `η_T = C_hu²/(ρ_u K_h) = 1` condition (pathA_39 Stage 4 §4). Full Maxwell
   (a Lorentz-force form for magnetism) needs BOTH `c_E=c_γ` AND `c_γ=c_s`. **The pathA_39 Stage 2/4 results are
   preferred-frame UNLESS `c_E=c_γ` — so the knit is what would (or would not) restore Lorentz invariance.**
2. **The NG cross-consistency gauntlet (NG1–NG5).** The cross-sector consistency program specified in the pathA_25 directive
   (`directives/pathA_25_gnls_polar_smectic_consistency_gates.md` — NG1–NG5 + a complete parameter ledger feeding NG5). Re-pose
   it on the now-earned four-sector parameter set: are the sectors mutually consistent, or does a cross-sector relation break?
3. **Assemble the four-sector chain into the central `pde_ledger`.** The end goal (memory [[project-calibrated-pde-goal]]): a
   fully-CALIBRATED, *simulation-ready* PDE delivering GR + EM. The four earned sector specs + the cone-lock + the NG ledger
   feed `research/pde_ledger/`. Completes the SPEC; the full nonlinear sim stays deferred
   ([[project-simulation-deferred-complete-pde-strategy]]) — a no-go is still possible.

## §2. Imported (per-sector; do NOT re-derive — read the earned reports)
- **Gravity:** `c_s (∝ρ²)`; `G` calibrated (`g_G`, `GENUINE_BLOCKED` — the PDE gives the FORM not Newton's `G`); PN ladder
  1PN→4PN GR-matched; the `54/5` quadrupole (`pathA_29`/`33`); `χ_Q` (two values, un-reconciled — `pathA_22b` `≈0.712` vs
  `pathA_33` `=1`, different contexts — see the STATUS open-reconciliation note).
- **Light:** `c_γ² = μ_R/ρ_br`; `B_eff = ρ_B0²/χ_c > 0` (`pathA_36`).
- **Charge:** `c_E`, `M_h > 0`, `q_h = 2Q_E tanh(b/ℓ)/b`, `Q_E` calibrated (`pathA_38`).
- **Magnetism:** the EM field multiplet = transverse-vector + charge-coupled-scalar (`h`-branon); preferred-frame unless
  `c_E=c_γ`; `C_hu` stability bound; `a_L` (`pathA_39`).
- **The `pathA_22b` verdict count** (reference framing): `P0 · χ_Q · g_mhat² · λγ⁵ / g_G = 54/5` — `λγ` enters to the FIFTH
  power and was pinned by the EM anchor. The cone-lock either DERIVES `λγ` or confirms it stays an anchor.

## §3. First post-`/compact` action — SCOPE the knit with Codex (before any directive)
The decisive setup questions (route through Claude↔Codex, then GLM, like the pathA_39 §2 / Stage-3 setups):
1. Is `λγ = c_γ/c_s = 1` a **derivable constraint** (forced by the shared `ρ`, `μ_R`, `χ_c`, … once all sectors share one
   medium state) or a **calibration input** (a genuine gap that must be anchored, zero predictive surplus)?
2. Is `c_E = c_γ` derivable or tuned? (It restores the Lorentz/Maxwell form the pathA_39 preferred-frame results lack.)
3. What is the **able-to-fail** structure? A knit that CANNOT be satisfied in one parameter set is a **first-class
   falsification** (the four sectors don't share one medium) — the gate must be able to emit that NO-GO, not just confirm
   consistency ([[feedback-falsification-is-the-goal]], [[feedback-negative-verdict-short-circuit]]). Per
   [[feedback-calibrate-predict-methodology]]: judge by predictive surplus; a calibration that merely ABSORBS a target
   (`λγ` anchored, zero held-out surplus) must be reported plainly.
Then: build the knit directive → Codex design-review (xhigh) → GLM tertiary → dual-engine (where a second engine can verify)
→ tri-review (arbiter + fidelity + adversarial-with-ablation on fresh agents) → user gate.

## §4. Honest expectation
The knit either **CLOSES** — the four sectors live in one calibrated parameter set → "one self-consistent brane+bulk
model," the program's stated win (completes the SPEC; sim deferred) — or it hits a **NO-GO** (a cone-lock / NG relation
can't hold) → a first-class falsification. Both are first-class results; do not rescue, do not oversell. Given `λγ` is a
flagged genuine gap, a likely landing is **cone-lock = CALIBRATED (anchored, not derived)** with the surplus living in the
held-out predictions (`g−2`, 5PN, ringdown, multi-defect) riding the shared `χ_Q`/`P0` bundle — report the earned-vs-calibrated
split plainly, exactly as the sector gates did.

## §5. References (read first, post-compact)
- `STATUS.md` ▶ RESUME HERE + `software/stage1_solver/decisions/13` §0 (the `pathA_22b` verdict equation + the two-`χ_Q`
  open reconciliation live here).
- `docs/conceptual_foundation.md` §3/§4 ⭐ v8 (the four-sector chain).
- `docs/development_plan.md` (the master "all sectors → simulation-ready" plan; the knit + pde_ledger assembly are its tail).
- `directives/pathA_25_gnls_polar_smectic_consistency_gates.md` (the NG1–NG5 gauntlet spec).
- `research/pde_ledger/` (the central calibrated-PDE ledger the chain assembles into).
- Sector reports: `reports/pathA_{29,36,38,39}_*`.
- Memories: [[project-pathA-build]], [[project-calibrated-pde-goal]], [[feedback-gravity-is-flow-not-ripple]],
  [[feedback-calibrate-predict-methodology]], [[project-simulation-deferred-complete-pde-strategy]],
  [[project-brane-existence-defect-structure]].
