---
unit_id: 155
batch: IV.6
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 155

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:94-97` now contains the four directive-specified `assert` lines anchoring `g_fp`, `S_fp`, `R_fp`, `Pi_fp` to the paper-stated constants (0.693352419668063, 0.6216013167514007, 0.2827139049082381, 1.4885734438300713) at tolerance 1e-12, inserted between the existing `print("Pi_fp=...)` (line 92) and the `dg = sp.Float(...)` block (line 99).

**Assessment:**
Edits match the directive verbatim. Constants taken from Mathematica `expectApprox` calls; non-tautological because they compare same-run moments to externally-fixed paper values. SymPy log exits 0 with the four prints unchanged from prior transcript, confirming the new asserts pass.

### F2 — paper_misalignment (banner)

**Classification:** resolved

**What changed:**
`scripts/.../sympy_audit.py:80` and `mathematica/.../mathematica_audit.wl:26` both now read `"STAGE 155 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"`. Fresh `.txt` transcripts (line 3 in each) now print the corrected `STAGE 155` banner.

**Assessment:**
Mechanical string swap matches directive. Docstring "Stage 154" cross-reference correctly left untouched per directive instruction. No collateral edits.

## Exec log assessment

**SymPy:** exit=0. Banner line 3: `STAGE 155 — FROZEN-TRACTION CO-EVOLVING FIXED POINT`. Moments printed at lines 7-10 unchanged (`g_fp = 0.693352419668063`, `R_fp = 0.2827139049082381`). Transport-law residual line 12 still well under 1e-10.

**Mathematica:** exit=0. Banner line 3 corrected. All five `expectApprox` PASS lines visible (lines 14, 16, 18, 20, 22) with diffs ≤ 2.22e-16.

**Output freshness:** SymPy `.txt` mtime 1779989421 > script mtime 1779984130; Mathematica `.txt` mtime 1779989505 > `.wl` mtime 1779945093. Both regenerated post-fix.

## Material-change assessment

`material_change`: false. The edits add assertions and correct banner strings; no derived numeric result changed. Downstream units unaffected by this stage.

## Side observations (non-blocking)

The docstring "Stage 154" vs. notes §3 "Stage 239" attribution of the transport law remains unresolved per directive; flagged by auditor as user-side resolution, not a verification blocker.

## Verdict justification

Both findings resolved with mechanical, directive-faithful edits. New SymPy asserts anchor the load-bearing moment quartet to paper-stated constants (non-tautological, independent of the transport-law identity check) and pass at ≤1e-12. Banner strings corrected in both engines and reflected in regenerated transcripts. Both exec logs exit 0 with no regressions.
