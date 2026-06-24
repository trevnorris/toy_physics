# pathA_25 — rolling verdict ledger (GNLS polar-smectic consistency gates)

Directive: `software/stage1_solver/directives/pathA_25_gnls_polar_smectic_consistency_gates.md` (v4, review-gauntlet SOUND).
This is the live "you are here" for pathA_25. Updated at each gate.

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
