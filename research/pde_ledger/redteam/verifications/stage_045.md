---
unit_id: 045
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T02:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: true
---

# Verification — unit 045 (batch III.1 v2)

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:36-42` — replaced the `Coefficient[Coefficient[couplingDensity, var1], var2]` extraction chain with `D[D[couplingDensity, var1], var2]` second cross-derivatives, plus an explanatory comment block. The four `gWext / gRext / gBext / gSext` definitions now use partial derivatives instead of `Coefficient` chains.
- The branch-equation independence (`Series[branchNumRaw, {rPhi, rU, 0}]` at wl line 99 vs SymPy `branch_eq.subs(R_phi, R_U)` at py line 143) is preserved per directive (leave-alone clause).
- The F_tr-check independence is established by F3: the WL script derives `fTrFromStage044` via `FullSimplify[fContStage044 /. {rPhi -> rU, lambda0 -> lambda0DN}, Assumptions -> $Assumptions]`, while the SymPy script uses `F_cont_stage044.subs([(R_phi, R_U), (lam0, lam0_dn)])`. Both routes import the Stage-044 `F_cont` residual but use different simplification machinery.

**Assessment:**
The directive's required `D[D[...]]` replacement at wl:36-42 is applied exactly as prescribed; the diff at `redteam/exec_logs/stage_045_diff.patch` lines 14-35 shows the `cWeta / cWU / cPhiEta / cPhiU` extraction chain removed and the four partial-derivative definitions inserted. The reference amplitudes (`gW`, `gR`, `gB`, `gS`) and the four `expectZero["g_X extracted - reference", ...]` assertions are preserved unchanged. The new `D[D[...]]` route is genuinely distinct from the SymPy `coeff(...).coeff(...)` chain — they would diverge if the coupling density were nonlinear, which is the appropriate independence property. The wl exec log confirms all four `PASS: g_X extracted - reference` assertions post-fix. No collateral edits.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:114-122` — removed the `channels = [...]` list, the `M_tr_channel_sum = sp.simplify(sum(...))` construction, the `print("M_tr_channel_sum = ...")` line, and the `expect_zero("M_tr - channel_sum", ...)` assertion. The `M_mix`, `M_supp`, `M_tr` symbolic definitions and their print lines remain. An explanatory comment block notes the substantive prefactor verification lives in Stages 022/026.
- `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:77-85` — removed the `channels = {...}` list, the `mTrChannelSum = FullSimplify[Total[...]]` construction, the `Print["M_tr_channel_sum = ..."]` line, and the `expectZero["M_tr - channel_sum", ...]` assertion. `mMix`, `mSupp`, `mTr` definitions and prints remain, plus the explanatory comment.

**Assessment:**
The edit matches the directive's required change verbatim. The exec logs confirm both transcripts no longer contain a `M_tr - channel_sum = 0` line; only the `M_mix`, `M_supp`, `M_tr` values are printed (sympy log lines 35-37; mathematica log lines 29-31). No tautological assertion has been replaced or smuggled in. No collateral edits.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:171-200` — replaced the hand-written generic-`lam0` `F_track` plus its lam0=2/9 specialization check with an import of the Stage-044 `D_cont` and `F_cont` residual expressions, restated inline with a comment citing `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90,140-146`. The script then:
  1. Defines `F_track_stage044 = sp.simplify(F_cont_stage044.subs(R_phi, R_U))` and asserts `expect_zero("Stage-044 tracking F collapse", F_track_stage044 - F_track_expected)` against the closed-form generic-`lam0` tracking expression.
  2. Defines `F_tr_from_stage044 = sp.simplify(F_cont_stage044.subs([(R_phi, R_U), (lam0, lam0_dn)]))` and asserts `expect_zero("F_tr collapse from Stage-044 residual", F_tr_from_stage044 - F_tr_expected)` against the notes' closed-form `F_tr` at lam0=2/9.
  3. The `coherent normalization residual` print now uses `R_target - F_tr_from_stage044` (upstream-derived) rather than the previous `R_target - F_tr_expected` (hand-written).
- `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:119-150` — mirror import of `dContStage044` / `fContStage044` with the same two assertions (`Stage-044 tracking F collapse`, `F_tr collapse from Stage-044 residual`), routed through `FullSimplify` with `Assumptions -> $Assumptions`. Comment cites `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:81-89,128-138`.
- Cross-checked: the imported `D_cont` / `F_cont` expressions in Stage-045 match Stage-044 SymPy lines 83-91 verbatim (D_cont: line 83-86, F_cont: line 87-91), and Stage-044 Mathematica lines 81-89 (dCont) / 85-90 (fCont) verbatim. The citation line ranges in the comments are slightly off by one line at the start (82 vs 83 for SymPy) but the cited content is the correct Stage-044 object.

**Assessment:**
The user-resolved Q2 (option a) was applied as prescribed: import the Stage-044 audit's `F_cont` residual expression, substitute the tracking branch and D/N value, and compare against the notes' closed-form `F_tr`. The new assertion `F_tr collapse from Stage-044 residual = 0` PASSes in both engines (sympy log line 49; mathematica log line 44). The check is non-tautological — if Stage-044's `F_cont` formula were wrong, or if the notes' `F_tr` closed form were wrong, the assertion would fail. The subsidiary `Stage-044 tracking F collapse` PASSes too (sympy log line 48; mathematica log line 42), providing an additional anchor at generic `lam0`. The directive's prescribed plan was followed without deviation. The "coherent normalization residual" print line correctly reads `R_target - F_tr_from_stage044` now (sympy log line 50; mathematica log line 45).

### F4 — paper_misalignment

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3` — docstring `Stage 28 SymPy audit.` → `Stage 045 SymPy audit.`
- `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:31` — banner `STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT` → `STAGE 045 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT`
- `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:202` (shifted from the original line 191 by lines inserted from F3) — print `\nAll Stage-28 symbolic checks passed.` → `\nAll Stage-045 symbolic checks passed.`
- `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26` — banner `STAGE 028 — COHERENT LOCAL TRACKING` → `STAGE 045 — COHERENT LOCAL TRACKING`
- The wl footer at line 153 `Print["Stage 045 Mathematica audit passed."]` remains unchanged per directive.
- The Applied block also notes `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:232` was relabeled to "Stage 045"; verifier scope is scripts-only, but I incidentally observed line 232 reads `## 7. Best current theorem statement after Stage 045` consistent with the Applied claim (flagged in side observations).

**Assessment:**
The four required script-side label edits are applied with no collateral changes. The exec logs confirm both transcripts now banner with "STAGE 045" (sympy log line 8; mathematica log line 8) and the SymPy footer reads `All Stage-045 symbolic checks passed.` (sympy log line 52). The wl footer continues to read `Stage 045 Mathematica audit passed.` (mathematica log line 47). No math content changed.

## Exec log assessment

**SymPy:** exit=0. Notable lines (from `redteam/exec_logs/stage_045_sympy.log`):
- Line 8: `STAGE 045 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT` — new banner active.
- Lines 14-17: four `g_X extracted - reference = 0` (F1 wl-side preserved on py-side as ext-vs-ref checks).
- Line 18: `g_B g_R - g_W g_S = 0`.
- Line 21: `rho_0 - sigma_0 = 0`.
- Lines 35-37: `M_mix`, `M_supp`, `M_tr` printed; no `M_tr - channel_sum` line (F2 confirmed).
- Line 44: `tracking quadratic collapse = 0`.
- Line 47: `G_tr D/N specialization = 0`.
- Line 48: `Stage-044 tracking F collapse = 0` (F3 subsidiary anchor PASS).
- Line 49: `F_tr collapse from Stage-044 residual = 0` (F3 primary PASS).
- Line 52: `All Stage-045 symbolic checks passed.` (F4 footer active).

No `Traceback` or `AssertionError`.

**Mathematica:** exit=0. Notable lines (from `redteam/exec_logs/stage_045_mathematica.log`):
- Line 8: `STAGE 045 — COHERENT LOCAL TRACKING` — new banner active.
- Lines 10-17: four `PASS: g_X extracted - reference` (F1 D[D[...]] route PASS).
- Line 19: `PASS: g_B g_R - g_W g_S`.
- Line 23: `PASS: rho_0 - sigma_0`.
- Lines 29-31: `M_mix`, `M_supp`, `M_tr` printed; no `M_tr - channel_sum` line (F2 confirmed).
- Line 35: `PASS: tracking quadratic collapse`.
- Line 40: `PASS: G_tr D/N specialization`.
- Line 42: `PASS: Stage-044 tracking F collapse` (F3 subsidiary).
- Line 44: `PASS: F_tr collapse from Stage-044 residual` (F3 primary).
- Line 47: `Stage 045 Mathematica audit passed.`.

No `$Failed` or `FAIL:` line.

**Output freshness:**
- `redteam/exec_logs/stage_045_sympy.log` mtime 1779781909 vs `scripts/.../sympy_audit.py` mtime 1779781431 — log newer, fresh post-fix.
- `redteam/exec_logs/stage_045_mathematica.log` mtime 1779781918 vs `mathematica/.../mathematica_audit.wl` mtime 1779781430 — log newer, fresh post-fix.
- The persistent transcripts at `scripts/output/...txt` (mtime 1779475829) and `mathematica/output/...txt` (mtime 1779475839) are stale (pre-fix). Not blocking — the verifier reads the redteam exec_logs per spec — but flagged in side observations.

## Material-change assessment

`material_change`: true.

F3 is the reason: the F_tr verification is now rooted in Stage-044's `F_cont` residual rather than a hand-written `F_track`. The numerical/symbolic residual values consumed downstream are unchanged (both engines still report the same `coherent normalization residual = (R_target ... ) - F_tr`), but the structural dependency now flows from Stage-044's `F_cont` into Stage-045's `F_tr` claim. Any downstream unit that quotes the Stage-045 F_tr deliverable now indirectly depends on Stage-044's `F_cont` being correct — which is consistent (Stage-044 is upstream of Stage-045 and has its own audit). Per the standard policy the orchestrator will mark units > 045 as `upstream_stale: true`; no narrower re-audit recommendation beyond that.

F1, F2, F4 do not change derived results — F1 is an algebraic-route change in the WL engine, F2 removes a tautological assertion, F4 is purely cosmetic. The material-change flag is driven solely by F3's structural import.

## Side observations (non-blocking)

1. The `## Applied: F4` block lists `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:232` among `files_changed`. The verifier scope is scripts-only; I incidentally observed that line 232 of that notes file reads `## 7. Best current theorem statement after Stage 045`, consistent with the Applied claim. The notes-side edit was authorized per user-approved Q3 (b). Flagged for orchestrator awareness, not blocking.

2. The persistent output transcripts at `scripts/output/...txt` and `mathematica/output/...txt` have not been regenerated post-fix (mtimes stale relative to scripts). The exec-log transcripts are fresh and demonstrate the fixes pass, so verification is not blocked, but if a downstream auditor reads from the `output/` transcripts they will see pre-fix banners/assertions. Recommend the orchestrator refresh these via `$RT exec-sympy 045` / `$RT exec-mathematica 045` before promoting Stage-045 to its `verified` checkpoint state (or confirm the exec_logs are the canonical post-fix transcript record).

3. The F3 SymPy citation comment reads `:82-90,140-146` but the actual `D_cont` block at Stage-044 begins at line 83 (off by one). Minor citation imprecision; the cited content is correct.

4. The wl `Series[branchNumRaw, {rPhi, rU, 0}]` independence at lines 97-100 was preserved per directive's leave-alone clause; this verifier did not re-examine its mathematical content (out of scope for verification).

## Verdict justification

All four findings have `## Applied` blocks in the directive; all four are confirmed by reading the post-fix script files and the captured diff at `redteam/exec_logs/stage_045_diff.patch`. The exec logs show both engines exit 0 with all expected new assertions PASS — in particular `F_tr collapse from Stage-044 residual` PASS in both engines, plus the subsidiary `Stage-044 tracking F collapse` PASS, plus the F1 `D[D[...]]` route in Mathematica producing PASS on all four extracted-vs-reference checks, plus the F2 removal of the tautological `M_tr - channel_sum` assertion confirmed in both transcripts, plus the F4 relabeled banners and footers in place. The F3 fix is genuinely non-tautological — it imports Stage-044's already-audited `F_cont` expression and compares the tracking-branch + D/N specialization against the notes' closed-form `F_tr`. No regressions are visible in the diff. Material change is true because the F_tr derivation now structurally depends on Stage-044, flagging units > 045 as upstream_stale.
