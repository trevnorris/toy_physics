---
unit_id: 152
batch: IV.6
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-28T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 152

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py:117-120` now contains the four anchored `expect_close` calls the directive specified, inserted between the lambda_eff prints and the pre-existing covariance consistency check:

```
expect_close("delta Pi_act scale", float(deltaPi), 0.907084414842908, tol=1e-6)
expect_close("delta Tm_act scale", float(deltaT), 0.271653979462338, tol=1e-6)
expect_close("lambda_eff^(Pi) scale", float(lam_Pi), 0.380487632771110, tol=1e-6)
expect_close("lambda_eff^(T) scale", float(lam_T), 0.378939241176339, tol=1e-6)
```

**Assessment:**
Edit matches the directive exactly (constants, ordering, tolerances, `float(...)` casts). The pre-existing `delta g from sigma1-linearized covariance check` is intact at lines 123-125. No collateral edits. Assertions are non-tautological: each binds an independently-computed mpmath quantity to a fixed numeric literal from the notes file. A regression in `Pi_star`, `Sigma_m_star`, `gprime_star`, `AT`, `BT`, or the lambda baselines would now flip the SymPy exit code.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `delta Pi_act scale = 0.9070844148429078   (target 0.907084414842908, diff 2.22e-16)`
- `lambda_eff^(Pi) scale = 0.38048763277111053   (target 0.38048763277111, diff 5.55e-16)`
- `delta g from sigma1-linearized covariance check ... diff 1.11e-16` (existing assertion preserved)

**Mathematica:** exit=0. Notable lines:
- `PASS: delta Pi_act scale`, `PASS: delta Tm_act scale`, `PASS: lambda_eff^(Pi) scale`, `PASS: lambda_eff^(T) scale`
- `Stage 152 Mathematica audit passed.`

**Output freshness:** SymPy script mtime 2026-05-28 10:01, transcript mtime 2026-05-28 11:30 (fresh). Mathematica transcript mtime 2026-05-28 11:31, .wl mtime 2026-05-11 (untouched, expected — directive only touched the .py).

## Material-change assessment

`material_change`: false. The edit only adds assertions; it does not alter any computed quantity, constant, or downstream result. No downstream units are affected.

## Side observations (non-blocking)

None.

## Verdict justification

The single F1 finding was applied verbatim per directive: four anchored `expect_close` assertions added at the specified location, each binding an independently-computed deliverable (`deltaPi`, `deltaT`, `lam_Pi`, `lam_T`) to its notes-quoted literal with `tol=1e-6`. The SymPy exec log confirms all four PASS with sub-1e-15 diffs, the pre-existing covariance consistency assertion is intact, and the Mathematica mirror continues to pass its five anchored checks. Both engines exit 0 with fresh outputs. The new assertions are non-tautological because the LHS values are computed from the canonical inputs while the RHS targets are fixed numeric literals from the notes — exactly the regression guard the auditor required.
