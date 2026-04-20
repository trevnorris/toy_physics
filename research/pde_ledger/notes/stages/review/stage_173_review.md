# Review: Stage 173 — Direct defect vs. dressing split

**Status:** Verified (SymPy PASS, hardened 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage173_direct_defect_vs_dressing_split.md`
- **Script:** `scripts/moving_throat_pde_stage173_direct_defect_vs_dressing_split_sympy_audit.py`

## Agent Review

### Agent: GPT-5 Codex — 2026-04-21
**Verdict:** PASS (HARDENED)

**Notes Derivation Review:**

The stage theorem is unchanged: the direct transfer-shape defect separates from
the dressing residual, and the direct observables remain support-blind.

**Script Review:**

The old support-blindness checks differentiated formulas that already omitted
`zeta`. The script now rebuilds `R_target` and `T^2` from the support-loaded
coherent baseline and propagates support-blindness into `N_*` from a carried
support-blind observable `R_tr`. The stage is still SymPy-only; no
Mathematica mirror exists yet.

**Issues Found:** None.

**Questions:** None.
