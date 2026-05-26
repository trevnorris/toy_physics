---
unit_id: 020
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 020

## Per-finding outcomes

### F1 — paper_misalignment

**Classification:** resolved

**What changed:**
Per user resolution Q5=b (TRIM), the Y20 / Gaunt lane-ratio scaffolding has been removed from both engines. The diff at `redteam/exec_logs/stage_020_diff.patch` shows:

- `scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py`: removed `from sympy.physics.wigner import gaunt` (formerly line 11), the `real_y20_square_ratio(m)` helper (formerly lines 20-26), and the three `assert_zero('Y20 overlap lane ...')` calls plus the three `lam2{0,1,2}` bindings (formerly lines 45-50). The wall-slope solve, even-gate Jacobian, dKSigma/dMSigma closed forms, the three compensated-D deficits, and the Xi1 residual all remain at lines 35-53.
- `mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl`: removed the top-of-file `ClearAll[GauntIntegral]` / `SetAttributes[GauntIntegral, Listable]` (formerly lines 6-7), the `gauntWeight` / `GauntIntegral` definitions (formerly lines 21-30 inside the Module), and the entire angular block (formerly lines 34-68, covering `overlapBase`, `lambda0`/`lambda1`/`lambda2`, `crossOne`/`crossTwo`, `angularResiduals`, and the five M1 guards including the `M1 OK` print). The Module's variable list was correspondingly pruned of `gauntWeight, overlapBase, lambda0, lambda1, lambda2, crossOne, crossTwo, angularResiduals`. The algebraic core (M2-M5: Jacobian determinant, unique-solution branch, dKSigma/dMSigma closed forms, three deficit checks, Xi1 residual) remains intact at the current lines 35-100.

**Assessment:**
The trim is exactly the direction (b) the original directive listed under "Possible directions". Every script-side artifact named in F1's `paper_missing_script_claim` (sympy `:45-50`, wl `:34-67`) is now absent, while every paper-anchored algebraic-core check (eq. stage020-wall-slopes through eq. stage020-residual-xi) is preserved. The destinations the user cited (stages 010 and 017) both already contain Gaunt machinery in their sympy scripts, so the Y20 lane content is not lost system-wide. No collateral edits visible in the diff: the only deletions are the named symbols/helpers/imports plus their declarations in the Module variable list, and no algebraic-core line was touched. Both engines still produce `STATUS: PASS`. F1 routes to `resolved` (paper card no longer overclaims relative to what the script asserts).

### F2 — tautological_check

**Classification:** resolved

**What changed:**
F2 was contingent on F1 per the original report and the directive's "Resolution gate" ("If user chose direction (b), the entire Y20 block is removed and this finding is moot"). The F1 trim deleted both `real_y20_square_ratio` (which constructed `gaunt(2,2,2,0,0,0)/base` for m=0) and the `assert_zero('Y20 overlap lane 20', lam20 - 1)` site on the SymPy side, and both the `overlapBase`/`lambda0` definitions and the `If[!TrueQ[FullSimplify[lambda0 - 1] === 0], ...]` guard on the Mathematica side. The tautological self-ratio simply no longer exists in either engine.

**Assessment:**
The defect surface is gone. There is no remaining assertion of the form `x - x == 0` masquerading as a content check. The F2 closure is mechanical — it requires no separate fix because the offending lines were deleted as part of F1.

## Exec log assessment

**SymPy:** exit=0. From `redteam/exec_logs/stage_020_sympy.log` (mtime 2026-05-21, but the script and saved output were re-touched 2026-05-25 21:42/21:49 — the saved output is the authoritative post-trim record):
- `dKSigma_from_even_gates = B01 + 27*B41 + Z01 + 27*Z41`
- `dMSigma_from_even_gates = -B21 + 3*B41 - Z21 + 3*Z41`
- `Check_D21_plus_D01_over_9 = 0` and `Check_D41_minus_2D21_over_3_minus_D01_over_27 = 0`
- `STATUS: PASS`

No `Y20 overlap lane` lines appear in the saved output, consistent with their removal.

**Mathematica:** exit=0 (inferred from `STATUS: PASS` line in saved output; no `stage_020_mathematica.log` was captured). Notable lines from `mathematica/output/.../stage020...txt`:
- `M2 even-gate Jacobian determinant residual = 0` / `M2 OK`
- `M3 solution count residual = 0`, `M3 dKSigma residual = 0`, `M3 dMSigma residual = 0` / `M3 OK`
- `M4 D01 residual = 0`, `M4 D21 residual = 0`, `M4 D41 residual = 0` / `M4 OK`
- `M5 Xi1 residual = 0` / `M5 OK`
- `STATUS: PASS`

The `M1` block (Y20 lane ratios) is correctly absent.

**Output freshness:** Confirmed.
- sympy script mtime: 2026-05-25 21:42; sympy output mtime: 2026-05-25 21:49 (output newer)
- mathematica script mtime: 2026-05-25 21:42; mathematica output mtime: 2026-05-25 21:50 (output newer)

## Material-change assessment

`material_change`: false.

The trim is purely subtractive — it removes script-only scaffolding (Y20 lane ratios and a tautological self-ratio) that no other unit consumes as a derived input. The paper-anchored algebraic-core results (wall-slope closed forms, `1/27` determinant, compensated D-deficits, Xi1 residual) are bit-identical in both engines pre- and post-trim. The user has independently confirmed that the Y20 lane verification exists in stages 010 and 017, so no claim is dropped at the project level; only the misplaced verification site is removed.

## Side observations (non-blocking)

- No `redteam/exec_logs/stage_020_mathematica.log` file is present; the Mathematica exit code is inferred from the `STATUS: PASS` terminator in `mathematica/output/.../stage020...txt`. This was flagged in the prior verification file and remains a minor orchestrator-bookkeeping gap; the freshness and content of the saved output are direct evidence of a successful run.
- The `redteam/exec_logs/stage_020_sympy.log` header date (2026-05-21) predates the current script/output mtimes (2026-05-25). The authoritative post-trim record is the saved `scripts/output/.../stage020...txt`, which is dated 2026-05-25 and contains no Y20 lane lines.
- The directive file at `redteam/directives/stage_020.md` has no `## Applied:` block from Codex because F1 was a paper_misalignment routed to user resolution and F2 was contingent on F1; the actual edits were applied via the user's Q5=b TRIM follow-up rather than the original codex apply step. The diff is complete and self-consistent regardless.

## Verdict justification

The F1 trim is mechanically clean: every script-side artifact named in the original finding is gone from both engines, the paper-anchored algebraic core is untouched, and both scripts still emit `STATUS: PASS` with all residuals identically zero. F2 was contingent on F1 and is auto-closed by the trim — the tautological `lambda_0 = 1` site no longer exists. No regressions visible in the diff, no collateral edits, no material change to downstream-visible content.
