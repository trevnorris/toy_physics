---
unit_id: 209
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

# Verification — unit 209

The single finding in the original report is non-script (F1 = stale committed SymPy `.txt`, banner reading "STAGE 192", refreshed by the orchestrator's independent re-run). No Codex source edit was directed and none was needed — there is no `directives/stage_209.md` and no `exec_logs/stage_209_diff.patch` (consistent with a non-script disposition). Verification confirms (A) the output refresh landed clean with the SymPy banner now "STAGE 209", (B) the audit disposition (INDEPENDENT `.wl`, no tautology) still holds on the refreshed artifacts, and (C) `material_change: false`.

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 192" banner)

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent re-run regenerated the committed SymPy output. `scripts/output/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.txt` now reads `STAGE 209 — EXACT PAIRWISE RATIO OPTIMIZER AND MIXED-RAY WINNER THEOREM` at line 3 and `STAGE 209 SYMPY AUDIT COMPLETED SUCCESSFULLY` at line 149 (the stale "STAGE 192" banner is gone). Its mtime is 2026-06-09 16:51:54, newer than the `.py` (2026-06-03 15:59:11). The Mathematica output was already fresh and was likewise refreshed (mtime 2026-06-09 16:51:54), banner `STAGE 209 -- PAIRWISE RATIO OPTIMIZER AND MIXED-RAY WINNER THEOREM` (line 3) / `STAGE 209 MATHEMATICA AUDIT PASSED` (line 70).

**Assessment:**
Correct and complete — the orchestrator re-run is the prescribed remedy. Every SymPy result line is `= 0`: `explicit algebraic tau form = 0` / `discriminant numerator reduction = 0` (txt 37–38); `tau = 2H0 / Phi = 0` / `Phi derivative law = 0` (63–64); `quartic degree minus 4 = 0` / `quartic factorization identity = 0` / `N - (J + 2(kj-kir)S) = 0` (88–90); diagonal-neutral (103–104) and pair-symmetry (117–118) both `= 0`. The Mathematica output prints an explicit `PASS:` for every check (lines 15, 21, 23, 29, 35, 37, 39, 46, 48, 50, 56, 58, 60, 62, 68, 70, 72) with no `FAIL` anywhere (grep FAIL on both `.txt` returns nothing). Output refresh is clean.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed artifacts. The `.wl` extracts A,B,C from the assembled polynomial via `CoefficientList` (M2, residuals `{0, 0, 0}`, mathematica txt line 20) rather than positing the literals as the `.py` does; eliminates the radical via `Resultant[…, z]` to build the quartic (M5), the single strongest discriminator; and solves out N by scaled differentiation and the gradient ray by `Solve` (M4, M6). All are derive-vs-posit routes absent from the `.py`.
- **Engine cross-check on the load-bearing quartic:** confirmed. The Mathematica M5 resultant (mathematica log line 44), collected at r⁴, gives `8·H0·(2·H0·v² − 2·ki·kj·v + ki²·w)`, and at r⁰ gives `8·H0·(2·H0·v² − 2·ki·kj·v + kj²·u)`. The SymPy §VI "coefficients highest → constant" emits exactly `8·H₀·(2·H₀·v² + kᵢ²·w − 2·kᵢ·k_j·v)` (sympy txt 139–140) and `8·H₀·(2·H₀·v² − 2·kᵢ·k_j·v + k_j²·u)` (txt 150–151). The two independent routes (difference-of-squares vs resultant) land on the identical quartic.
- **Stationarity not vacuous:** confirmed. Both engines assert the derivative is non-trivial before each reduction — `M4 Phi derivative is not identically zero before reductions = True` (mathematica txt 34) and `M6 raw tau derivative is not identically zero before diagonal-neutral reduction = True` (txt 55) — so the `= 0` stationarity checks are real can-fail statements.
- **0 reconciliation misalignments:** confirmed. The report's 13-value deliverable table reconciles MATCH against the appendix/notes; every reconciled symbolic form (tau, A/B/C, disc numerator, Phi, N, derivative law, quartic Q, gradient ray, pair-symmetry) is present and identical in the refreshed outputs.

## Exec log assessment

**SymPy:** exit=0 (`# exit_code: 0`, sympy log line 156). Notable lines: `explicit algebraic tau form = 0` (8); `Phi derivative law = 0` (64); `quartic factorization identity = 0` / `N - (J + 2(kj-kir)S) = 0` (89–90); `gradient-optimal stationarity on diagonal-neutral branch = 0` (104); `equal-mix stationarity on pair-symmetric branch = 0` (118). All equalities hold to zero.

**Mathematica:** exit=0 (`# exit_code: 0`, mathematica log line 77). Notable lines: `PASS: M2 discriminant coefficient extraction residuals` (21); `PASS: M4 scaled derivative numerator minus manifest N` (37); `PASS: M5 resultant quartic degree is 4` (46); `PASS: M5 plus radical factor equals derivative numerator` (50); `PASS: M6 gradient ratio recovered from slope stationarity` (58). Every check prints PASS; no FAIL.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime 2026-06-09 16:51:54, newer than the `.py` (2026-06-03 15:59:11) and `.wl` (2026-06-02 10:47:39). SymPy banner is now canonical "STAGE 209" (the F1 stale-192 banner is cleared).

## Material-change assessment

`material_change`: false. No source code changed; the only edits were the regenerated committed `.txt` outputs (SymPy banner relabel 192→209 + transcript refresh). No derived result changed, so no downstream unit is affected.

## Side observations (non-blocking)

None beyond the single reported finding. The absent finite-set-count / bracket / promotion / winner deliverables are combinatorial/comparison wrappers (StatusNumerical), correctly not flagged by the auditor as missing CAS identities; I concur and add nothing.

## Verdict justification

`verified`. The lone finding F1 is a non-script stale SymPy banner ("STAGE 192") cleared by the orchestrator re-run — the refreshed `.txt` now carries the canonical "STAGE 209" banner (mtime newer than the `.py`) with every SymPy residual `= 0`, and the Mathematica output prints PASS on all 17 checks with no FAIL. No directive or diff exists because no Codex edit was needed, which matches the report's disposition. The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (A,B,C by `CoefficientList`, quartic Q by `Resultant`, N by scaled differentiation, gradient ray by `Solve` — none mirroring the `.py` posit route), the engine cross-check on the quartic's r⁴/r⁰ coefficients agrees, and reconciliation is 13/13 MATCH. No regressions; `material_change: false`.
