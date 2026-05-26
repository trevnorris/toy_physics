---
unit_id: 042
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 042

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py:92` — `F_general = sp.simplify(Z_expected * S_expected / (1 - xi))` replaced by `F_general = sp.simplify(Z_overlap * S_overlap / (1 - xi))`.
- `mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl:70` — `fGeneral = FullSimplify[zExpected sExpected/(1 - xi), Assumptions -> $Assumptions]` replaced by `fGeneral = FullSimplify[zOverlap sOverlap/(1 - xi), Assumptions -> $Assumptions]`.

The diff at `redteam/exec_logs/stage_042_diff.patch` shows exactly these two single-line edits plus the directive's `## Applied: F1` block and front-matter update; nothing else was touched. The downstream lines (`F_expected` definition, prints, factorization, `expect_zero` / `expectZero` call) are untouched in both engines, matching the directive's "do not change lines 93-101 / 71-82" instruction.

**Assessment:**

The edit is correct and matches the directive verbatim. The new construction routes the F_general residual through the eigenvector-constructed overlap forms `Z_overlap` / `zOverlap` and `S_overlap` / `sOverlap` (defined earlier from `ratio_expected` / `ratioExpected`, sympy line 78, mathematica lines 62-63). The assertion `F_general - F_expected == 0` is no longer algebraically guaranteed by substitution — it only passes after the simplifier chains the `Z_overlap -> Z_expected` and `S_overlap -> S_expected` reductions (A3/A4 and B2a/B2b) through the product. If either of those reductions failed, the new `F_general - F_expected` residual would not simplify to zero, so the check is now genuinely non-tautological.

Both exec logs confirm the residual still simplifies to zero post-fix: sympy log line 79 prints `F_general - expected = 0` and mathematica log lines 32-33 print `F_general - expected = 0 / PASS: F_general - expected`. The printed factored form of `F_general` in the sympy log (lines 63-78) shows the full unsimplified rank-2 denominator `(xi-1)*(delta^2 - 2*delta*m*q^2 + ... + xi^2)`, confirming the simplifier really did chain from the constructed overlaps rather than starting from an already-collapsed `Z_expected*S_expected/dQR^2` form. The check now has real verification content.

No collateral edits in the diff. Codex did not run scripts (it stated so explicitly, and the directive forbade it); the orchestrator captured the exec logs afterwards.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Line 29-30: `row1 - expected = 0` / `row2 - expected = 0` (A1, A2 pass).
- Line 61-62: `Z_overlap - expected = 0` / `S_overlap - expected = 0` (A3, A4 pass — these are now load-bearing for the F1 fix).
- Line 79: `F_general - expected = 0` (A5 passes under the new non-tautological construction).
- Lines 92, 113, 118, 127, 134-136: tracking collapse, source-tied specialization, flat-U recovery, H_n / H_F, linear expansions — all `= 0` (A6-A12 unaffected).

**Mathematica:** exit=0. Notable lines:
- Lines 28-33: `Z_overlap - expected = 0 / PASS`, `S_overlap - expected = 0 / PASS`, `F_general - expected = 0 / PASS` (B2a, B2b, B2c pass).
- Lines 39-40: `tracking collapse = 0 / PASS` (B3).
- Lines 46-47, 52-53, 60-67: source-tied, flat-U, H_n / H_F, linear expansions — all PASS (B4-B7).
- Line 69: `Stage 042 Mathematica audit passed.`

**Output freshness:** the exec logs at `redteam/exec_logs/stage_042_{sympy,mathematica}.log` are dated `2026-05-26T01:50:20` / `01:51:04`, i.e. after Codex's `applied_at` of `2026-05-26T07:39:05Z` (which is 01:39 in mountain time, the log timestamps are local). The exec logs are therefore post-fix. The saved `scripts/output/...txt` and `mathematica/output/...txt` transcripts have older mtimes than the edited scripts (output mtimes ~1778524959, script mtimes ~1779781120). This is expected — the exec-runner under `redteam/` writes to `redteam/exec_logs/` rather than overwriting the canonical `scripts/output/` and `mathematica/output/` transcripts. The substance-bearing post-fix output is preserved in the exec logs.

## Material-change assessment

`material_change`: false.

The edit changes only the construction path of the `F_general` intermediate, not any printed numerical result or downstream-consumed quantity. `F_general` equals `F_expected` both before and after the fix (they are provably equal via A3/A4); only the algebraic route to the equality is different. No saved closed-form expression, no constant, no theorem-ledger entry changes. Downstream units that depend on Stage-25's outputs (the F_(q,r,t) closed form, the source-tied F_src, H_n^(src), H_F^(src)) see the same values as before. No upstream-stale flag needed downstream from this change alone.

## Side observations (non-blocking)

- The canonical `scripts/output/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.txt` and `mathematica/output/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.txt` transcripts have not been refreshed since Codex's edit (their mtimes are older than the edited scripts). The post-fix exec logs in `redteam/exec_logs/` confirm passing runs, but a future re-pull of the saved-output transcripts will show the pre-fix `F_general` construction in the un-displayed source. Not a verification-blocker since the assertion-level results are identical (`F_general - expected = 0` either way) and the exec logs are authoritative for this round; flagging in case the orchestrator wants to refresh the canonical outputs.
- The auditor's "Self-test notes — Soft note" about `eps > 0` declared positive while the constructive branch has `R_U < 1` (eps < 0) is not part of any formal finding and was not addressed by this fix. It does not affect the F1 verification.

## Verdict justification

The single low-severity finding (F1, tautological_check) is fully resolved: the directive's verbatim one-line change has been applied to both the SymPy and Mathematica scripts with no deviations or collateral edits, both engines re-run cleanly to exit 0, and the `F_general - expected = 0` assertion now has genuine verification content because it chains through the independently-asserted `Z_overlap -> Z_expected` and `S_overlap -> S_expected` simplifications. No regressions in the diff. No new findings introduced (per verifier scope). Verdict: `verified`.
