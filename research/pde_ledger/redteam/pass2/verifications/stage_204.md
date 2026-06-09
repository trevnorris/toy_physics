---
unit_id: 204
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 204

The sole finding in the original report is non-script (F1 = stale committed SymPy `.txt` carrying a pre-renumber "STAGE 187" banner, captured 2026-05-11 from a revision older than the Jun-3 `.py`). The directive directed no source edit — only a re-run/overwrite of the saved output — and the orchestrator's independent re-run regenerated it. Verification confirms (A) the output refresh landed clean and the SymPy `.txt` now reads STAGE 204, (B) the audit disposition (INDEPENDENT `.wl`, no tautological/under-exercised assertions) still holds on the refreshed artifacts, and (C) `material_change: false`.

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 187" banner)

**Classification:** resolved

**What changed:**
No source change. The directive's F1 has no `## Applied:` block because it asked for "No source edit — re-run the SymPy script and overwrite the saved output," and there is no diff patch on disk (`stage_204_diff.patch` absent) — consistent with an output-refresh-only remedy. The orchestrator's independent re-run regenerated both committed transcripts. `scripts/output/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.txt` now reads `STAGE 204 — EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT PREDICTOR` at line 3 and `STAGE 204 SYMPY AUDIT PASSED` at line 221 (the stale "STAGE 187" banner at old lines 11/229 is gone). Its mtime is Jun 9 16:51, newer than the `.py` (Jun 3 15:59). The Mathematica output (already fresh per the report) was likewise refreshed (mtime Jun 9 16:51): `STAGE 204 -- EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT PREDICTOR` (line 3), `STAGE 204 MATHEMATICA AUDIT PASSED` (line 96).

**Assessment:**
Correct and complete. The refresh is exactly the prescribed remedy. The regenerated SymPy `.txt` carries the canonical stage number, every residual check prints `= 0` (34 such lines; the only non-`= 0` `=` lines are the expected primitive-table exponent tuples `e_lambda/e_c/e_gamma/e_U/e_W`, which are non-zero exponent vectors by construction, not failures), and the transcript ends "STAGE 204 SYMPY AUDIT PASSED". No FAIL/`!= 0`/Traceback anywhere. Output refresh is clean.

## Disposition re-confirmation (post-refresh)

- **Independent `.wl` (load-bearing exponents + invariance):** the auditor's INDEPENDENT verdict holds on the refreshed artifacts. The `.wl` recovers each dependent exponent by log-derivative — `sigmaDeltaRecovered = logRate[deltaGraph[tau]]` = `D[Log[deltaGraph[tau]], tau]` (wl 99–101) — plus a constancy check `D[sigmaXRecovered, tau] == 0`, whereas the `.py` posits the closed form and proves it by equating two closed forms (`delta_graph_tau_direct - delta_graph_tau_exp`, py 100/109). Distinct extraction operations for the same exponent (same for T/Keta/mu). Monomial invariance is reconstructed by inverting the timing law `deltaFromTiming[kU,tU] := Pi^2 tU/(L^2 kU)` (wl 131) applied to the lifted ray vector, rather than re-substituting the direct delta. The acknowledged shared sub-pieces (the `M_*` Stage-192 literal in M5, the predictor closed form in M7) are posit-vs-posit but secondary; their load-bearing inputs (`sigmaKeta/sigmaMu/sigmaT`) are independently recovered in M2/M3. No overturn.
- **No tautological/hardcoded/under-exercised assertions:** confirmed. Every exponent equality compares an independently-built quantity against the posited formula, so a wrong sign/factor would surface. The predictor self-test is non-tautological — both engines emit the same series (eps^2 coeff `-1/(2L0)`, eps^3 coeff `+2/(3L0)`), matching the hand-expansion of `tau_log - tau_aff`.
- **Reconciliation 16/16 MATCH, 0 misaligned:** the report's value table reconciles every derived exponent formula, the invariance identities, the primitive table, and the predictor formulas against notes/appendix; all 16 deliverables MATCH. No numeric figure-of-merit constants are emitted (only structural `= 0` residuals), so the banner relabel changes no reconciled value.

## Exec log assessment

**SymPy:** exit=0 (exec log `stage_204_sympy.log:6` `# exit_code: 0`). Notable lines: banner `STAGE 204 — EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT PREDICTOR` (log L8); terminal `STAGE 204 SYMPY AUDIT PASSED` (log L226). No Traceback/Error.

**Mathematica:** exit=0 (exec log `stage_204_mathematica.log:103` `# exit_code: 0`). Notable lines: banner `STAGE 204 -- EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT PREDICTOR` (log L8); terminal `STAGE 204 MATHEMATICA AUDIT PASSED` (log L101). No FAIL/`$Failed`.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime Jun 9 16:51, newer than the `.py` (Jun 3 15:59) and the `.wl`. The SymPy banner is now canonical "STAGE 204" (the F1 stale-187 banner is cleared). No diff patch exists, confirming no source file was touched — only the saved transcripts were regenerated.

## Material-change assessment

`material_change`: false. No source code changed (no diff patch); the only edits were the regenerated committed `.txt` transcripts (banner relabel STAGE 187 → STAGE 204 + transcript refresh). No derived result changed — the scripts are fully symbolic and emit only `= 0` structural residuals plus the (unchanged) exponent/predictor formulas — so no downstream unit is affected.

## Side observations (non-blocking)

None beyond the one reported finding. The grep for non-`= 0` `=` lines in the refreshed SymPy `.txt` surfaces only the five primitive-direction exponent rows (`e_lambda … e_W`), which are intentionally non-zero exponent tuples reported as table deliverables — not assertion failures. I concur with the auditor's INDEPENDENT verdict and add nothing.

## Verdict justification

`verified`. The lone finding is non-script and resolved: F1's stale "STAGE 187" SymPy banner is cleared by the orchestrator's independent re-run — the refreshed `.txt` now reads "STAGE 204" with every structural check `= 0` and ends "STAGE 204 SYMPY AUDIT PASSED", and the already-fresh Mathematica output was likewise refreshed and ends "STAGE 204 MATHEMATICA AUDIT PASSED". Both exec logs report exit 0. The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent on the load-bearing exponents (log-derivative recovery) and invariance (timing-law inversion), no assertion is tautological, and reconciliation is 16/16 MATCH. No source change, no diff, `material_change: false`.
