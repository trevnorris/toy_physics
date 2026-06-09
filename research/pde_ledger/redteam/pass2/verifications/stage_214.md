---
unit_id: 214
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

# Verification — unit 214

The sole finding in the original report is non-script (F1 = stale committed SymPy `.txt`
carrying the pre-renumber `STAGE 197` banner, refreshed by the orchestrator's independent
re-run). No Codex source edit was directed and none was needed — the directive explicitly
forbids touching `.py`/`.wl` logic and the diff patch is empty. Verification confirms
(A) the output refresh landed clean (SymPy banner now `STAGE 214`), (B) the audit
disposition still holds on the refreshed artifacts AND the notes-typo corrections still
hold (150 = 5·5·6, 54, no stray 218/162/230), and (C) `material_change: false`.

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 197" banner)

**Classification:** resolved

**What changed:**
No source change (diff patch `stage_214_diff.patch` is empty). The orchestrator's
independent re-run regenerated the committed SymPy output. `scripts/output/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`
now reads `STAGE 214 — FULL INTERIOR FOUR-COORDINATE SIMPLEX OPTIMIZER` at line 3 and
`STAGE 214 SYMPY AUDIT COMPLETED SUCCESSFULLY` at line 1209 (the stale `STAGE 197` banner
at the former line 3 / line 1209 is gone). Its mtime is Jun 9 16:51, newer than the `.py`
(Jun 3 15:59). The Mathematica output was likewise refreshed (mtime Jun 9 16:51) and its
banner already read STAGE 214.

**Assessment:**
Correct and complete. The re-run is the prescribed remedy and required no source edit.
The directive's verification criteria are all met on the refreshed `.txt`: line 3 reads
`STAGE 214 — …`, the final line reads `STAGE 214 SYMPY AUDIT COMPLETED SUCCESSFULLY`,
`lifted Bezout bound = 54` (line 122), `projected one-chart Bezout bound = 150` (line 1160),
and both interior order theorems are `verified … on 924 ordered integer samples`
(lines 1205-1206). No collateral edit — content other than the banner is unchanged.

## Disposition re-confirmation (post-refresh)

- **Output refresh clean (A):** confirmed. SymPy exec log (`stage_214_sympy.log`) shows
  `STAGE 214` banner (line 8), `STAGE 214 SYMPY AUDIT COMPLETED SUCCESSFULLY` (line 1214),
  `exit_code: 0` (line 1216), `lifted Bezout bound = 54` (line 127), `projected one-chart
  Bezout bound = 150` (line 1165), and both 924-sample theorems (lines 1210-1211). The
  committed `.txt` matches.
- **Notes-typo corrections still hold (B):** confirmed. The projected one-chart bound is
  `5·5·6 = 150` in both engines (sympy out line 1160; `M5 projected product minus 5*5*6 = 0`
  / `PASS` in `.wl` out lines 82-83), the lifted four-coordinate bound is `54` (not 162),
  and a stray-token scan of the refreshed committed `.txt` finds no `218`, `162`, or `230`
  (grep clean). The degree ledger `{3,3,3,2}` (`.wl` out line 31) and eliminant degrees
  `{5,5,5}`/`{6,6,6}` (`.wl` out lines 63-64) reproduce the corrected bounds.
- **`.wl` genuinely independent:** confirmed on the refreshed artifacts. The load-bearing
  eliminants are derived by `Resultant[..., y]` lift-variable elimination (`M4 C_rs/…/S_t
  resultant minus definition = 0` → all PASS, `.wl` out lines 51-62) and checked equal to
  the `.py`'s posited closed forms (`Crs = Ms*Lr - Mr*Ls`, `Sr = Lr^2 - 4*Mr^2*Delta`).
  That derive-vs-posit split is the transliteration discriminator; the `.wl` is not a port.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `STAGE 214 — FULL INTERIOR FOUR-COORDINATE SIMPLEX
OPTIMIZER` (L8); `exact dPhi/dr compiler = 0` (L77); `lifted Bezout bound = 54` (L127);
`projected one-chart Bezout bound = 150` (L1165); `verified interior winner theorem on 924
ordered integer samples` (L1210) and `… non-improvement theorem on 924 …` (L1211);
`STAGE 214 SYMPY AUDIT COMPLETED SUCCESSFULLY` (L1214); `# exit_code: 0` (L1216). No
AssertionError, no FAIL.

**Mathematica:** exit=0. Notable lines: `M3 lifted Bezout product minus 3*3*3*2 = 0` →
`PASS` (out L45-46); `M4 S_r resultant minus definition = 0` → `PASS` (out L57-58); `M5
projected product minus 5*5*6 = 0` → `PASS` (out L82-83); `All Stage 214 Mathematica audit
checks passed.` (out L110); `# exit_code: 0` (out L111). Every `M*` check prints PASS; no
FAIL. The re-run took 72s, well under the 600s cap.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime Jun 9 16:51, newer than
the `.py` (Jun 3 15:59) and `.wl` (Jun 2). The SymPy banner is now canonical `STAGE 214`
(the F1 stale-197 banner is cleared).

## Material-change assessment

`material_change`: false. No source code changed (empty diff); the only edits were
regenerated committed `.txt` outputs (banner relabel `197`→`214` plus a fresh transcript).
No derived result changed (54, 150, all `=0` residuals, the 924 sweep counts are identical),
so no downstream unit is affected.

## Side observations (non-blocking)

- Deliverable (7) (interior winner / non-improvement order theorems, SymPy §VI, 924 samples)
  is single-engine — the `.wl` has no M-section counterpart. The auditor already noted and
  accepted this as a trivial real-interval-ordering statement; I concur and add nothing.
- The card's `\stagefield{Verification}` reads "Mathematica audit: none yet" while a passing
  `.wl` exists (card-text lag). This is paper-side prose, outside the scripts-only scope, and
  is deferred to the batch paper-cleanup step (P4-51); non-blocking, no Codex action.

## Verdict justification

`verified`. The single finding is non-script and resolved: F1's stale `STAGE 197` SymPy
banner is cleared by the orchestrator re-run — the refreshed committed `.txt` carries the
canonical `STAGE 214` banner (line 3) and completion line (line 1209) with bound 54, bound
150, and both 924-sample theorems intact, while the Mathematica output prints PASS on every
M1-M7 check (exit 0, 72s). The audit disposition holds on the refreshed artifacts: the
notes-typo corrections still hold (150 = 5·5·6, 54, no stray 218/162/230), and the `.wl`
remains genuinely independent (eliminants via `Resultant[…,y]` elimination vs the `.py`'s
posited closed forms). The diff patch is empty (output-only refresh), so there are no
regressions; `material_change: false`.
