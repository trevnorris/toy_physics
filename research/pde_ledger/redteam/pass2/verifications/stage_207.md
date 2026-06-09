---
unit_id: 207
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

# Verification — unit 207

The sole finding in the original report is non-script (F1 = stale committed SymPy `.txt` carrying the pre-renumber `STAGE 190` banner, refreshed by the orchestrator's independent re-run). No Codex source edit was directed and none was needed; the directive carries no `## Applied:`/`## Blocked:` block. Verification confirms (A) the output refresh landed clean (SymPy banner is now `STAGE 207`), (B) the audit disposition holds on the refreshed artifacts (INDEPENDENT `.wl`, no tautology), and (C) `material_change: false` (only the captured `.txt` transcript changed, no derived result).

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 190" banner)

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent SymPy re-run regenerated the committed transcript. `scripts/output/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.txt` now reads `STAGE 207 — PRIMITIVE-RAY HESSIAN ENVELOPES AND CERTIFIED RAY TABLE` at line 3 and `STAGE 207 SYMPY AUDIT PASSED` at line 234 — the stale `STAGE 190` banner is gone (grep for `STAGE 190` in the refreshed `.txt` returns nothing). Its mtime is 2026-06-09T16:51:54, newer than the `.py` (2026-06-03T15:59:11). `git status` shows only this one `.txt` modified; `git diff HEAD` on the `.py` and `.wl` is empty, confirming no script-math edit. The Mathematica `.txt` was already fresh and was likewise refreshed (mtime 2026-06-09T16:51:54), banner `STAGE 207` (line 3), `STAGE 207 MATHEMATICA AUDIT PASSED` (line 128), no `FAIL`.

**Assessment:**
Correct and complete. The re-run is the prescribed remedy. The refreshed SymPy exec log (`exec_logs/stage_207_sympy.log`) carries the `STAGE 207` banner (L8), every residual `= 0` across §I–§VII (diagonal reductions L16-33, mixed-ray form L40-41, orientation law L56, all five drift rows = 0 L70-158, lower/upper comparison quadratics L175-176, turning lower/upper quadratics L201-202), and `STAGE 207 SYMPY AUDIT PASSED` / `# exit_code: 0` (L239,241). No banner/label artifact remains. Output refresh is clean.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed artifacts. The discriminator is the bracket root maps (M4/M5). The `.wl` *derives* the monotone root by `Solve`+`Reduce[…,Reals]`+`SelectFirst` branch-selection — the refreshed Mathematica log prints the explicit positive-root set with the `Sign[kappaLo]`/`Sign[kappaHi]` curvature-sign Piecewise (L90-91) and then matches the closed form (`M4 tau_lo/tau_hi selected root - closed form = 0`, L92-99) — whereas the `.py` *posits* `2*H0/(k+sqrt(k**2-2*c*H0))` and only residual-checks the comparison quadratic. Same Solve-derive-vs-posit contrast for the turning bracket (M5 L104-111: `Solve`+positive-root selection then `tau - Sqrt[-2H0/kappa] = 0` vs the SymPy posited `√(-2H0/κ)`). Method differs, not merely the simplifier → INDEPENDENT.
- **Non-tautological load-bearing checks:** confirmed. The bracket checks substitute the candidate back into the *independently stated* comparison/turning quadratic (not into its own definition), so they can fail; the elementary diagonal/mixed/drift checks would fail under mis-indexing or a mis-typed row. M2 exercises the cross-coefficient over all 10 axis pairs (log L38-77, each `= 0` + PASS).
- **0 reconciliation misalignments:** confirmed. Report reconciles 7/7 deliverable values MATCH (all symbolic closed forms / table rows, no pinned floats); every value is present and identical in the refreshed outputs (monotone root `2H0/(k+√(k²-2cH0))` py-out L163-174; turning root `√2·√H0/√a_turn` py-out L189-200; orientation `−Gam²/|Gam|`→0 py-out L48-56; five drift rows = 0 py-out L70-158).

## Exec log assessment

**SymPy:** exit=0. Notable lines: banner `STAGE 207 — PRIMITIVE-RAY HESSIAN ENVELOPES AND CERTIFIED RAY TABLE` (L8); `+e_lambda diagonal reduction = 0` … `-e_W diagonal reduction = 0` (L16-33); `canonical orientation law = 0` (L56); `lower comparison quadratic = 0` / `upper comparison quadratic = 0` (L175-176); `turning lower quadratic = 0` / `turning upper quadratic = 0` (L201-202); `STAGE 207 SYMPY AUDIT PASSED` / `# exit_code: 0` (L239,241). All residuals zero.

**Mathematica:** exit=0. Notable lines: every check prints an explicit `PASS:` — M1 diagonal (L15-33), M2 ×10 pairs quadratic-form + cross-coefficient (L39-77), M3 orientation (L83), M4 monotone root-map incl. positive-root-set Piecewise + selected-root + quadratic residual (L89-99), M5 turning root-map (L105-111), M6 drift table (L118-130); `STAGE 207 MATHEMATICA AUDIT PASSED` / `# exit_code: 0` (L133,135). No `FAIL` anywhere.

**Output freshness:** confirmed. Both committed `.txt` outputs carry mtime 2026-06-09T16:51:54, newer than the `.py` (2026-06-03T15:59:11) and `.wl` (2026-06-02T10:06:44). SymPy banner is now canonical `STAGE 207` (the F1 stale-190 banner is cleared).

## Material-change assessment

`material_change`: false. No source code changed (`git diff HEAD` on `.py`/`.wl` is empty); the only edit is the regenerated committed SymPy `.txt` transcript (banner relabel `190`→`207` plus transcript refresh, identical residuals). No derived result changed, so no downstream unit is affected.

## Side observations (non-blocking)

None beyond the single reported finding. The elimination/winner theorems are correctly carried as inequality-logic corollaries of certified-bracket admissibility (no residual identity to assert), so their absence from the scripts is not a gap; I concur and add nothing.

## Verdict justification

`verified`. The lone finding is non-script and resolved: F1's stale SymPy banner is cleared by the orchestrator re-run — the refreshed `.txt` carries the canonical `STAGE 207` banner with every SymPy residual `= 0` and `EXIT_CODE: 0`, and the Mathematica output prints PASS on every check with no FAIL. `git status` shows only the one SymPy `.txt` modified, with no `.py`/`.wl` source edit. The audit disposition holds on the refreshed artifacts: the retrofit `.wl` is genuinely INDEPENDENT (Solve+Reduce+SelectFirst branch-derived root maps vs the `.py`'s posited closed forms), the load-bearing bracket checks are non-tautological, and reconciliation is 7/7 MATCH with 0 misalignments. No regressions; `material_change: false`.
