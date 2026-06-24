# pathA_25 — rolling verdict ledger (GNLS polar-smectic consistency gates)

Directive: `software/stage1_solver/directives/pathA_25_gnls_polar_smectic_consistency_gates.md` (v4, review-gauntlet SOUND).
This is the live "you are here" for pathA_25. Updated at each gate.

## ▶ RESUME HERE (post-/compact, 2026-06-23)

**State:** directive `pathA_25` v4 review-gauntlet SOUND + committed (`6cdaa821`). **G0 DONE** = `SECOND_MEDIUM_DRIFT_AT_FREEZE`
(5 independent inputs; tri-reviewed genuine; committed `77fd0e72`, hash `f00ee99d465e`). User-gated: **proceed to Gate B4.**

**NEXT ACTION — Gate B4 (the T1-replacement; first physics test):**
1. Draft `_scratch/pathA_25_gateB4_execute_prompt.md`, bound to directive §8 (Codex designs the numerical route; prompt states
   requirements + able-to-fail labels + controls NC-B4a/b/c/d). Question: does the **baseline `B0_Lifshitz`** (Family-L) driver make a
   **stable, codim-1, emergent-axis smectic ground state** (vs uniform / collapse / phase-separation / lattice)? Use the **full coupled
   (ρ, θ, Pⁱ) spectrum**; **derive** the layer profile (do NOT hand-draw it — the T1 lesson); verify rotationally-invariant de Gennes
   `B,K > 0` from the derived profile; show single-`k` preferred over multi-`k`.
2. **B4 is a CRUX gate → design-review the execution prompt with BOTH Codex (xhigh, read-only) AND GLM (user-mediated)** before the run
   (§18). Fold fixes.
3. Execute B4 (Codex workspace-write, xhigh) → arbiter re-run both engines + transliteration-fidelity audit + clean adversarial agent
   → user gate.

**Process reminders:** `codex exec --sandbox {read-only for reviews | workspace-write for execution} -m gpt-5.5 -c
model_reasoning_effort=xhigh`, run backgrounded, verify log shows `reasoning effort: xhigh`, never wrap the codex session in shell
`timeout`. Freeze-fidelity guard at every gate: recompute T0 `8fa41ac51e88` + G0 `f00ee99d465e`. **Baseline B4 = Family-L ONLY** (no
R/C sensitivity branches; the Family-C `k→0` derivation in the G0 scripts is tautological and must be fixed before any C-branch).
Reports→`reports/pathA_25_*`, scripts→`tools/pathA_25_*`, scratch→`_scratch/`. Commit only when the user asks; stage explicit paths.
**Post-compact sync TODO (low priority):** `decisions/13` §0 + repo `STATUS.md` active-frontier still describe the pre-G0 plan; this
ledger + memory are the fresh source of truth.

## Gate ledger

| Gate | Verdict | Provenance | Tri-review | Notes |
|---|---|---|---|---|
| **G0** structure freeze | `SECOND_MEDIUM_DRIFT_AT_FREEZE` (hash `f00ee99d465e…`) | — | ✅ arbiter PASS · `FIDELITY_CLEAN` · `ADVERSARIAL_SOUND` | 5 independent new inputs (`c_L1,c_L2` driver + `μ_br,J_Pu,κ_Pu` light). Honest, expected NG5 firing. User-gated: proceed. |
| **B4** smectic ground state | _pending_ | — | — | Baseline `B0_Lifshitz` (Family L). The T1-replacement / first physics test. |
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
- **Family-C `k→0` derivation is tautological** in the G0 dim-check scripts — fix before any Family-C sensitivity branch
  (`S_L_plus_Cdiv` / `S_L_plus_Cpin`) is run. Baseline Family-L is unaffected.
- G0.6 prose ("agreed on all restored-unit reductions") is narrower than literal: kept-GNLS + T0 `L_pol` terms are hash-guarded, not
  re-dim-checked (they were checked in their own freezes).
