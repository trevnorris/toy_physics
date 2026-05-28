---
unit_id: 156
batch: IV.6
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 156

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py:127-134`, four new numeric guards were inserted directly after the preserved `g_can`/`R_can` asserts (lines 123–126): `Sigma0_can` vs 4.651033550168867 at 1e-10; `S_can` vs 0.6703621156734617 at 1e-11; `Pi_can` vs 3.871564377479002 at 1e-10; `T_hat_can` vs 1.4467083664567613 at 1e-10. Surrounding prints and `Conclusion` block untouched. The Mathematica `.wl` was not modified (correctly — directive forbade it).

**Assessment:**
Edit matches the directive verbatim (string-for-string with the "After" block). The four new tolerances mirror the Mathematica `expectApprox` tolerances one decade looser as specified. Non-tautological: each compares an independently computed quantity (`bisect`-driven `Sigma0_can`, the moment integrals `S_can`/`Pi_can`, and the closure `T_hat_can=sqrt(9 Sigma0/20)`) against the notes-anchored numeric target — they would fail under a discretization drift (wrong `N`, weights, kernel) where the original two guards would still pass. No collateral edits.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `Sigma0_can = 4.651033550168867`, `S_can = 0.6703621156734617`, `Pi_can = 3.871564377479002`, `T_hat_can = 1.4467083664567613` — exact agreement with the asserted targets; conclusion block printed cleanly (asserts silent on success, as expected).

**Mathematica:** exit=0. Notable lines: `PASS: monotone increase on scan grid`, `PASS: Sigma0_can numeric check`, `PASS: g_can`, `PASS: S_can`, `PASS: R_can`, `PASS: Pi_can`, `PASS: T_hat_can numeric check`. Unchanged from prior run (script not modified).

**Output freshness:** SymPy script mtime 1779984136 < output mtime 1779989422; Mathematica script mtime 1779945093 < output mtime 1779989509. Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. Only assertion guards were added; no derived numeric output, kernel, grid, or closure was altered. Downstream units that consume the six canonical numbers see byte-identical values.

## Side observations (non-blocking)

The shared banner "STAGE 139 — RENORMALIZED CANONICAL BRANCH" carryover in the Mathematica `.wl` (flagged in the auditor's self-test notes) was not addressed and was explicitly out of scope for this directive. Cosmetic only.

## Verdict justification

Codex (via direct orchestrator edit) applied the directive verbatim: four non-tautological numeric guards added at the specified line range, no collateral changes, no Mathematica edits. Both exec logs exit 0, the SymPy printed values agree with the new assertion targets to within tolerance, and Mathematica PASS lines are intact. F1 is fully resolved.
