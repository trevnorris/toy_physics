---
unit_id: 210
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T16:55:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 210

Both findings in the original report are non-script (F1 = paper-side card-text lag, "Mathematica audit: none yet", USER-DEFERRED to P4-51; F2 = stale committed SymPy `.txt` carrying the P4-52 "STAGE 193" banner, refreshed by the orchestrator's independent re-run). No Codex source edit was directed and none was needed — the directive contains NO applied script edits (`applied: false`, `## Resolve before fix_loop` holds F1 for the user, F2 is an informational re-run). There is no `stage_210_diff.patch` (it does not exist), consistent with zero source edits. Verification confirms (A) the output refresh landed clean (SymPy now STAGE 210), (B) the audit disposition still holds on the refreshed artifacts, (C) F1 is correctly classed and routed to user-deferral, and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (card-text lag: "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to P4-51)

**What changed:**
Nothing in scripts (paper.tex is off-limits to the red-team). The directive holds `stage_210.tex:11` under `## Resolve before fix_loop` for the user to pick a direction; expected resolution (a) is to update the Verification line to cite `mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl`. Codex applied nothing for F1, which is the prescribed behavior. The card-text lag is USER-DEFERRED to P4-51.

**Assessment:**
Correctly classed. The `.wl` exists, runs, and passes — the refreshed Mathematica committed output ends `STAGE 210 MATHEMATICA AUDIT PASSED` (line 118) with `PASS:` on every M1-M9 check and no FAIL — so the card's "Mathematica audit: none yet" understates coverage but contradicts no result. This is a paper-prose coverage-statement staleness with zero math impact, legitimately deferred and outside the scripts-only scope. No script edit was possible or appropriate.

### F2 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 193" banner)

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent re-run regenerated the committed SymPy output. `scripts/output/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.txt` now reads `STAGE 210 — THREE-COORDINATE MIXED-SIMPLEX AND THE CANONICAL TRIPLE-SCREEN AUDIT` at line 3 and `STAGE 210 SYMPY AUDIT COMPLETED SUCCESSFULLY` at line 165; a grep for `STAGE 193` returns 0 hits (the P4-52 stale banner is fully cleared). Its mtime is 2026-06-09T16:51, newer than the `.py` (2026-06-03T15:59). The Mathematica output was already fresh and was likewise refreshed (mtime 2026-06-09T16:51), banner `STAGE 210 -- THREE-COORDINATE MIXED-SIMPLEX MATHEMATICA AUDIT`.

**Assessment:**
Correct and complete. The re-run is the prescribed remedy (no `.py` edit). Every SymPy result line is an equality `= 0`: §I lines 56-61; §II lines 88-94; §III lines 101-106; §IV lines 115-118; §V lines 147-152; §VI lines 165-167. Zero `FAIL` in the file. Output refresh is clean and the check substance is unchanged (residuals were already `= 0` pre-refresh; only the banner label and freshness gap were stale).

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed artifacts. The load-bearing gradient-optimal vector is derived (not posited) via a Lagrangian `Solve` and confirmed unique — the refreshed Mathematica output prints `M3 positive Lagrange branch = {ki/Sqrt[ki^2 + kj^2 + kk^2], ...}` (line 34), `M3 Lagrange stationarity residuals = {0, 0, 0, 0}` (line 35), and `M3 stationary branch minus gradient vector = {0, 0, 0}` (line 37), whereas the `.py` posits `avec_grad = sp.Matrix([ki/Kgrad, ...])` and verifies it. The barycenter maximizer uniqueness is `Reduce`-proved: `M6 wSigma == 2 unit-simplex condition = ai == 1/Sqrt[3] && aj == 1/Sqrt[3] && ak == 1/Sqrt[3]` with `M6 equal-mix uniqueness condition = True` (lines 77-78) — strictly stronger than the `.py`'s value-only check. The A..F coefficients are extracted by `Series`/`CoefficientList` (`M7 CoefficientList in {r,s} = {{ki^2 - 2*H0*uii, ...}}`, line 84) vs the `.py`'s whole-expression `simplify`-to-zero. Independent route confirmed — not a port.
- **Non-tautological assertions:** confirmed. Each M1-M9 row defines an object by one route and checks it against an independently-stated closed form (or solves/reduces and compares); every one prints `PASS` (lines 15, 17, 19, 24-29, 36-42, 48-52, 58-62, 68-79, 86-96, 102-108, 114-120) with no FAIL.
- **0 reconciliation misalignments:** confirmed. The report's value-reconciliation table reconciles 18/18 deliverable values MATCH against the notes (a_grad, max slope, interior ratios, w_Σ=(Σa)²−1, w_Σ≤2, barycenter 1/√3, the curvature law, the τ root map, the six A..F coefficients, the τ ratio-coordinate form, and the edge reductions); all are symbolic (no pinned floats), and each appears identically in the refreshed outputs. The only paper-side discrepancy is the non-numeric F1 coverage-statement lag.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `edge ij slope reduction = 0` (L59); `gradient-optimal normalization = 0` (L88); `w_Sigma - ((sum a)^2 - ||a||^2) = 0` (L101); `w_Sigma(equal-mix barycenter) - 2 = 0` (L105); `discriminant numerator reduction = 0` (L148); closing `STAGE 210 SYMPY AUDIT COMPLETED SUCCESSFULLY` (L170). All equalities hold to zero; exit_code 0.

**Mathematica:** exit=0. Notable lines: `M3 Lagrange stationarity residuals = {0, 0, 0, 0}` + `PASS` (L35-36); `M3 stationary branch minus gradient vector = {0, 0, 0}` + `PASS` (L37-38); `M6 equal-mix uniqueness condition = True` + `PASS` (L78-79); `M7 coefficient A residual = 0` + `PASS` (L85-86); closing `STAGE 210 MATHEMATICA AUDIT PASSED` (L123). Every M1-M9 check prints PASS; no FAIL; exit_code 0.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime 2026-06-09T16:51, newer than the `.py` (2026-06-03T15:59) and the `.wl` (2026-06-02T10:15). The SymPy banner is now canonical `STAGE 210` (the F2 stale-193 banner is fully cleared — 0 grep hits) and the Mathematica banner is `STAGE 210`.

## Material-change assessment

`material_change`: false. No source code changed; the only edits were regenerated committed `.txt` outputs (banner relabel 193→210 + transcript refresh). No derived result changed, so no downstream unit (> 210) is affected by this verification.

## Side observations (non-blocking)

The note that the diff patch `stage_210_diff.patch` is absent is expected and correct here — the directive applied no source edits — and not a defect. The `.wl` mtime (2026-06-02) predates the `.py` (2026-06-03), but both committed outputs were re-run on 2026-06-09, so the freshness criterion (output newer than its own script) is satisfied for both engines. Nothing further to flag.

## Verdict justification

`verified`. Both findings are non-script and resolved: F2's stale SymPy banner is cleared by the orchestrator re-run — the refreshed `.txt` carries the canonical `STAGE 210` banner top and bottom, zero `STAGE 193` residue, every SymPy check `= 0`, and the Mathematica output prints PASS on all M1-M9 with no FAIL; F1 is a paper-side card-text lag ("Mathematica audit: none yet") correctly USER-DEFERRED to P4-51. The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (Lagrangian `Solve` for a_grad, `Reduce`-proved barycenter uniqueness, `Series`/`CoefficientList` for A..F — all routes absent from the `.py`), every assertion is non-tautological, and reconciliation is 18/18 MATCH with 0 numeric misalignments. No Codex edit was directed, no diff exists, no regressions; `material_change: false`. Both exec logs exit 0.
