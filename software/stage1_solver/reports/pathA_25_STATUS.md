# pathA_25 — rolling verdict ledger (GNLS polar-smectic consistency gates)

Directive: `software/stage1_solver/directives/pathA_25_gnls_polar_smectic_consistency_gates.md` (v4, review-gauntlet SOUND).
This is the live "you are here" for pathA_25. Updated at each gate.

## ▶ RESUME HERE (2026-06-24)

> **⚠ FRONTIER PIVOTED (2026-06-24): pathA_25 gates are PAUSED.** B4 = `FAIL_NOT_CODIM1` (final; density-smectic NO-GO). The active
> frontier is now the **throat-soliton / 4D-light synthesis** — READ FIRST: repo `STATUS.md` ⭐⭐ LATEST STATE +
> `docs/conceptual_foundation.md` §4/§5/§6.1 + memory `project-light-is-4d-throat-hypothesis`.
> **UPDATE (2026-06-24): step (1) DONE — pathA_26 Derrick + open-stability check = `THROAT_DRAIN_DESTABILIZED`, interpreted NOT-a-kill**
> (conservative throat-soliton EXISTS generically; instability only at unphysical large drain = black-hole regime). See
> `reports/pathA_26_derrick.md` (+ its Interpretation §) and STATUS.md ⭐⭐. **Plan was TWO INDEPENDENT tracks: (2a)** drain-sector
> derivation; **(2b)** the R/C cubic verdict.
> **UPDATE (2026-06-24): (2b) DONE — `RC_DENSITY_SMECTIC_LIGHT_NOGO` (tri-reviewed, dual-engine, earned).** The density-smectic brane
> is a NO-GO across ALL admitted drivers: `S_R_kernel`→`FAIL_NOT_CODIM1` (C.5; bilinear+isotropic, GNLS cubic survives),
> `S_L_plus_Cdiv`→`FAIL_ADMISSION` (C.2; Goldstone integrate-out gives O(k⁰) directional EOS-stiffness softening for any λ_Cdiv≠0),
> `S_L_plus_Cpin`→`FAIL_LIGHT_STARVED` (C.6; χ_Cpin<0 DOES open a genuine codim-1 density-smectic window — the cubic IS escapable —
> but the same coupling pins P along the layer normal, `P_∥=0` at the static Ω_u=0 minimizer, starving in-plane light). So density CAN
> form a brane but only by killing light → strong earned evidence for the pivot: the brane must be a **light-confining shear surface**,
> not a density modulation. The Family-C `k→0` tautology is RESOLVED (§B.0 derived kernel). See `reports/pathA_25_gateRC_cubic.md`
> (+`_results.yaml`), `tools/pathA_25_gateRC_cubic_{sympy.py,.wl}`, prompt `_scratch/pathA_25_gateRC_execute_prompt.md` (v4).
> **(2a) drain-sector derivation is the ONLY remaining throat-soliton track — PAUSED (user, 2026-06-24); resume there next.**
> The pathA_25 gate material below is preserved as history; with (2b) done, **no live pathA_25 thread remains** (the density route is closed).

**State:** directive `pathA_25` v4 review-gauntlet SOUND + committed (`6cdaa821`). **G0 DONE** = `SECOND_MEDIUM_DRIFT_AT_FREEZE`
(5 inputs; tri-reviewed; committed `77fd0e72`, hash `f00ee99d465e`). **Gate B4 execution prompt v4 = THROUGH THE FULL REVIEW GAUNTLET**
(`_scratch/pathA_25_gateB4_execute_prompt.md`): Codex green over 3 rounds (review r1 + confirm r1 + confirm r2 `SOUND`) → single GLM
pass `SOUND_WITH_CONCERNS` (recovered via `opencode export` after the emission hung; corrections #1–5 folded — load-bearing = no
3D-BCC anchoring, 4D `S^3` weak-crystallization must be derived; #6 skipped, GLM withdrew it) → Codex re-validation (confirm r3)
`SOUND`. **▶ B4 EXECUTION IN PROGRESS** (Codex workspace-write xhigh, bg task `bh2sq3zt6`, log `_scratch/pathA_25_gateB4_exec.log`).

**B4 RESULT (2026-06-24): `FAIL_NOT_CODIM1` (preliminary, adversarial-confirmed).** Baseline `B0_Lifshitz` makes a finite-`k`
instability (Λ>1) but the nonlinear onset is a rank-2 equilateral-triad multi-`k` state, not a codim-1 lamella — the GNLS cubic
invariant (`U'''=15Kρ0²>0` + `c_s²(ρ)` dependence) gives the triad a cubic energy lowering the single-`k` lamella lacks. Adversarial
review `ADVERSARIAL_SOUND` (independently reproduced; the `Γ_T=0` rescue is a genuine **codim-1** surface in the 6-D control space, NOT
an open window → NOT `SMECTIC_CONDITIONAL`). Fidelity review `FIDELITY_ISSUES` (backbone genuinely computed, but the verdict string +
some formulas were hardcoded/pasted → being remediated). **▶ Remediation in progress** (Codex `b6z8er5qa`, danger-full-access): fix the
`.wl` idiom bug + RUN it (real dual-engine) + make both scripts COMPUTE the verdict (not assert it) + correct the report's false
"expired license" line. **Verdict will NOT change** (adversarial-confirmed); remediation makes the scripts earn it.

**B4 VERIFIED FINAL (2026-06-24):** dual-engine PASS (arbiter-reran both; `engine_agreement: PASS`), verdict COMPUTED, adversarial
SOUND, fidelity remediated, report corrected. `FAIL_NOT_CODIM1` stands. Not committed (user commits when asked).

**▶ HISTORY (these two threads are now resolved/superseded — see RESUME HERE at top):**
- **(A) R/C cubic cycle — ✅ DONE (2026-06-24): `RC_DENSITY_SMECTIC_LIGHT_NOGO`.** The pre-registered prediction VERIFIED almost
  exactly: **R = no** (bilinear/isotropic → `FAIL_NOT_CODIM1`); **C = the candidate** — and the C-branch CRUX came TRUE: the SAME χ<0
  `E_Cpin` that opens a genuine codim-1 density-smectic window pins `P` out-of-plane (`P_∥=0`) and starves in-plane light →
  `FAIL_LIGHT_STARVED`, the deep NO-GO between brane-exists & light-exists. (Cdiv added a third fail mode: O(k⁰) Goldstone EOS-softening
  → `FAIL_ADMISSION`.) Density-smectic brane route CLOSED.
- **(B) ⭐ 4D-light reframe** ([[project-light-is-4d-throat-hypothesis]], user intuition 2026-06-24): light is NOT 3D/in-brane — it's a
  4D shear field (3+1D face on brane via gapped `w`-mode → 2 vac polarizations; 4D at throats → interaction/leak; trapped = mass =
  geon). Retrofits C5 + curvature-leak + geon-mass + dissolves the A-crux P-tension. Reframes Gate L (test 4D field w/ gapped `w`-mode,
  not 3D shear). Able-to-fail on polarization count (4+1 massless = 3 pols; brane must gap exactly 1 → 2).
- Orchestrator lean: **A-then-B** (R/C is cheap + already set up; settles the density-smectic route; its result informs B — if the
  density-smectic can't form, that's evidence the brane is defined by light-confinement, not density). Do NOT run Gate L/S/B/Q/T on a
  non-existent codim-1 layer (T1 lesson) until the brane question resolves.

**Process reminders:** `codex exec --sandbox {read-only reviews | workspace-write execution} -m gpt-5.5 -c
model_reasoning_effort=xhigh`, backgrounded, verify `reasoning effort: xhigh`, never wrap the codex session in shell `timeout`.
**Baseline B4 = Family-L ONLY** (Family-C `k→0` tautology has since been FIXED in Gate R/C §B.0 — see Carry-forward "RESOLVED").
Reports→`reports/pathA_25_*`, scripts→`tools/pathA_25_*`, scratch→`_scratch/`. Commit only when the user asks; stage explicit paths.
GLM-direct = `opencode -m cloudflare-workers-ai/@cf/zai-org/glm-5.2 run '...'`; on the known output-emission hang, recover via
`opencode export <sessionID>` (full reasoning persists) — do NOT restart from scratch (re-hangs).
**Post-compact sync DONE (2026-06-26):** `decisions/13` §0 top block + repo `STATUS.md` now reframe this EM/brane-existence material as
*history* (density route CLOSED, drain-sector PAUSED); the live front is the moving-throat PDE push. This ledger stays as the EM-track record.

## Gate ledger

| Gate | Verdict | Provenance | Tri-review | Notes |
|---|---|---|---|---|
| **G0** structure freeze | `SECOND_MEDIUM_DRIFT_AT_FREEZE` (hash `f00ee99d465e…`) | — | ✅ arbiter PASS · `FIDELITY_CLEAN` · `ADVERSARIAL_SOUND` | 5 independent new inputs (`c_L1,c_L2` driver + `μ_br,J_Pu,κ_Pu` light). Honest, expected NG5 firing. User-gated: proceed. |
| **B4** smectic ground state | **`FAIL_NOT_CODIM1`** (FINAL) | dual-engine ✅ (SymPy + `.wl` both exit 0, `engine_agreement: PASS`, 31 shared exprs; arbiter-reran) | `ADVERSARIAL_SOUND` (genuine, NOT CONDITIONAL — Γ_T=0 is codim-1) · `FIDELITY` remediated (verdict now COMPUTED: lamella min −1.67e-5 vs triad min −0.090 → triad wins) | Baseline `B0_Lifshitz`. Triad beats lamella via GNLS cubic `U'''=15Kρ0²>0`. Candidate **NO-GO** for the density-smectic brane. |
| **R/C** sensitivity branches | **`RC_DENSITY_SMECTIC_LIGHT_NOGO`** (FINAL, earned) | dual-engine ✅ (both exit 0, `engine_agreement: PASS`, 72 shared exprs; arbiter-reran ×2) | `ADVERSARIAL_CONCERNS`→remediated · `FIDELITY_ISSUES`(C.6 `P_∥` hardcoded)→remediated→`REMEDIATION_CLEAN` | R→`FAIL_NOT_CODIM1`(C.5); Cdiv→`FAIL_ADMISSION`(C.2, O(k⁰) Goldstone softening); Cpin→`FAIL_LIGHT_STARVED`(C.6: opens a real codim-1 window but pins P⊥-plane, `P_∥=0` at static Ω_u=0). **Density-smectic brane route CLOSED across all admitted drivers** → brane must be a light-confining shear surface. Family-C `k→0` tautology RESOLVED (§B.0). |
| **L** light (CRUX) | _not started_ | — | — | On B4's profile. C5 left able-to-fail (no φ-analog). |
| **S** magnetism | _not started_ | — | — | Bulk reservoir on the B4 background. |
| **B** gravity bundle | _not started_ | — | — | Anisotropic-`c_s` / extra-mode / connectivity on B4. |
| **Q** charge | _not started_ | — | — | Two ±w puncture directions on B4. |
| **K** cone-lock | _not started_ | — | — | `c_γ` vs `c_s`; expected CALIBRATION_GAP. |
| **T** throat | _not started_ | — | — | Minimal T-compat required for any structure-consistent roll-up. |

## NG gauntlet state (§13)
- **NG5 (single-medium ↔ everything):** PRESSURE RAISED at G0 (5 ≥ 2 independent inputs). Recorded; not fatal per §14. The 5 inputs
  are the calibration budget — judged against held-out surplus if the structure survives the gates.
- NG1–NG4: not yet evaluable (need B4/L outputs).

## Carry-forward fixes / caveats
- **RESOLVED (2026-06-24 Gate R/C): Family-C `k→0` tautology fixed** — G0 checks now derive the Cdiv integrated-out Goldstone response and Cpin direct `O(k²)` density Hessian outside the freeze-action hash block.
- G0.6 prose ("agreed on all restored-unit reductions") is narrower than literal: kept-GNLS + T0 `L_pol` terms are hash-guarded, not
  re-dim-checked (they were checked in their own freezes).
