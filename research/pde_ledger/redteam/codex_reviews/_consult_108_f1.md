# Consult — unit 108 F1 (paper_misalignment) resolution

**Date:** 2026-05-29
**Mode:** Claude (orchestrator) + evidence agent + Codex read-only consult (`codex_session_id: 019e7279-d223-7263-95f9-b2a09371a96b`)
**Question:** Are stage-108 card Checks #2/#3 (Robin / standalone mixed-pole no-go / compensated Robin–mixed) genuinely verified at sibling stages 109–113, so stage 108 needs NO script change (direction a) — or is one absent, forcing stage 108 to carry it (direction b)?

## Evidence (both engines, real failing assertions, quantities derived not hardcoded)

- **Check 1 — Robin `chi_Q^R = 3/(3 - rho_R)`** → stage **110**
  - `scripts/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.py:29` `assert sp.simplify(chi_R - 3/(3 - rho)) == 0` (`chi_R = c5/(1/27)` from Robin DtN series)
  - `mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl:52` `expectZero["chi_Q^R - 3/(3 - rho)", chiR - 3/(3 - rho)]`
- **Check 2 — mixed-pole no-go (`kappa_W=-1/9` then `sigma_W=0`)** → stage **111**
  - `scripts/...stage111_mixed_sidechannel_pole_sympy_audit.py:32-33` `assert kappa_match + 1/9 == 0`; `assert sigma_match == 0` (`sigma_match` solved AFTER substituting `kappa=kappa_match`)
  - `mathematica/...stage111...wl:71-72` `expectZero["kappa_match + 1/9", ...]`; `expectZero["sigma_match", ...]`
- **Check 3 — compensated `chi_Q^hyb = (1-9 sigma_W gamma_W)/(1-sigma_W)`, preservation iff `gamma_W=1/9`, pairs `(sigma_W,0)` & `(4 sigma_W,1/3)`** → stage **112**
  - sympy `...stage112...py:38,48,49,51-53`; wl `...stage112...wl:50-53,75,76` plus an independent route (wl:64-66) deriving `gamma_W=1/9` from `a_0/3 + 9 a_5 = 0`

(Stage 113 has no audit script — status-only consolidation card; stage 109 only sets up/cross-references the downstream checks.)

## Verdict — CONCUR with direction (a)

Codex (read-only, opened the files itself) and the Claude evidence agent both concur: all three checks have a genuine dual-engine home at 110/111/112; none is absent from the 109–113 range. **Stage 108 needs no script change.** F1 is a paper-card over-scoping — the stage-108 card's Checks #2/#3 should be cross-referenced to stages 110/111/112 in the manual paper pass. Logged to `PAPER_CLEANUP_TRACKER`. The red-team makes no stage-108 script edit for F1.

## Side-finding (NOT stage 108 — recorded for stage 112's own fix loop)

Codex precision caveat on **stage 112**: the "preservation iff `gamma_W=1/9`" holds for the nontrivial branch only. From `chi_B = (1-9 sigma_W gamma_W)/(1-sigma_W)`, `chi_B=1` ⟹ `sigma_W(1-9 gamma_W)=0`, so `sigma_W=0` preserves for any `gamma_W`. The script does not state `sigma_W != 0`. This is a precision caveat (not a missing check). Fold the `sigma_W != 0` qualifier / case split into stage 112's own upcoming FINDINGS fix loop; do not scope-creep into batch 1.
