---
unit_id: 208
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

# Verification — unit 208

The original report raised a single finding, F1 = `stale_output` (committed SymPy `.txt` carried the pre-renumber `STAGE 191` banner). It is non-script and explicitly required no Codex edit ("None mandatory ... regenerates on re-run"); the auditor wrote no directive ("no Codex-applied script edit is warranted, so no directive file is written"). Consistent with that, no `directives/stage_208.md` and no `exec_logs/stage_208_diff.patch` exist — correct, not a missing artifact. Verification therefore confirms (A) the orchestrator's output refresh landed clean (SymPy banner now `STAGE 208`), (B) the audit disposition (including the INDEPENDENT `.wl` verdict) holds on the refreshed artifacts, and (C) `material_change: false`.

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber `STAGE 191` banner)

**Classification:** resolved

**What changed:**
No source change (none directed, none needed). The orchestrator's independent re-run regenerated the committed outputs. `scripts/output/moving_throat_pde_stage208_..._sympy_audit.txt` now reads `STAGE 208 — PAIRWISE MIXED RAYS AND OFF-DIAGONAL HESSIAN SYNERGY` at line 3 and `STAGE 208 SYMPY AUDIT COMPLETED SUCCESSFULLY` at line 198; the stale `STAGE 191` banner is gone (grep for `191` across both committed `.txt` outputs returns nothing). Its mtime is 2026-06-09T16:51:54, newer than the `.py` (2026-06-03T15:59). The Mathematica `.txt` (already correct at audit time) was likewise refreshed to the same mtime, banner `STAGE 208 -- MATHEMATICA PAIRWISE MIXED-RAY AUDIT`.

**Assessment:**
Correct and complete. The refresh is the prescribed remedy for a P4-52 stale-banner artifact; the captured math content was already faithful and is unchanged. Every SymPy result line is an equality `= 0` (20 such lines, zero failures): mixed slope law (L38 of log), derivative law (L45), gradient-optimal ratio stationarity / slope / gain (L54-56), curvature decomposition + neutrality (L72-73), cross-weight derivative/stationarity/max (L80-82), the canonical-ray laws (L112-113, L134-135), envelope weighted forms (L152-153), the certified-bracket quadratic residuals (L186-187), and the gradient cross-weight formula (L200). The committed Mathematica output prints `PASS:` on every check (27 PASS) with no `FAIL`. Output refresh is clean.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed Mathematica output. The two load-bearing objects are extracted by a real method difference, not a re-simplification of the `.py`'s algebra:
  - **r_grad solved, not posited.** The output shows `Solve stationarity roots = {{r -> kj/ki}}` (log L23), `Reduce stationarity region = ...` (L24) reduced on-domain to `kj == ki*r` (L25), with the uniqueness gate `M3 unique stationary root count = 0` / PASS (L26-27) and `M3 stationarity region equals r == kj/ki = True` / PASS (L28-29). The SymPy side instead *posits* `r_grad = kj/ki` and checks stationarity at it — a count/uniqueness check is absent there. Different operation, same object.
  - **Cross-weight coefficient-extracted, not posited.** The output reports `coefficient weights {w_i, w_x, w_j} = {(1+r^2)^(-1), (2*r)/(1+r^2), r^2/(1+r^2)}` (L39) recovered via the `Series`/`Coefficient` projection on `mixedCurvature*(1+r^2)`, with `M4 coefficient-recovered cross weight = 0` / PASS (L44-45). The SymPy side posits `w_x = 2r/(1+r²)` literally. Different operation.
  The shared Rayleigh quotient / cone monomials are the legitimately-shared physical premise; the gate-defining objects of both main theorems are independently derived, so the `.wl` clears the bar.
- **0 reconciliation misalignments:** the report reconciled 16/16 symbolic deliverables MATCH; all are present and identical in the refreshed outputs (e.g. slope law SymPy txt L33-34 ↔ wl txt L16; `r_grad=kj/ki` SymPy L47-49 ↔ wl `Solve` root L23; weights SymPy L66-71 ↔ wl L39; certified bracket SymPy L170-185 ↔ wl L90-91). No pinned numeric constants (the only literals are structural), so no constant-provenance drift is possible.

## Exec log assessment

**SymPy:** exit=0. Notable lines: "mixed slope law = 0" (L38); "mixed slope derivative law = 0" (L45); "gradient-optimal ratio stationarity = 0" (L54); "mixed curvature decomposition = 0" / "diagonal neutrality when h_ij = 0 = 0" (L72-73); "cross-weight stationarity at r=1 = 0" (L81); "quadratic root relation for tau_lo/tau_hi = 0" (L186-187); "STAGE 208 SYMPY AUDIT COMPLETED SUCCESSFULLY" (L203). All residuals exactly 0.

**Mathematica:** exit=0. Notable lines: "PASS: M3 unique stationary root count" (L27); "PASS: M3 stationarity region equals r == kj/ki" (L29); "PASS: M4 coefficient-recovered cross weight" (L45); "PASS: M5 cross-weight value at r=1" (L51); "PASS: M9 lower/upper bracket closure quadratic" (L93-95). Every check (27) prints PASS; no FAIL.

**Output freshness:** confirmed. Both committed `.txt` outputs carry mtime 2026-06-09T16:51:54, newer than the `.py` (2026-06-03T15:59) and `.wl` (2026-06-02T10:04). The SymPy banner is now canonical `STAGE 208`; no `191` survives in either committed output.

## Material-change assessment

`material_change`: false. No source code changed; the only edits were the regenerated committed `.txt` outputs (SymPy banner relabel 191→208 plus transcript refresh). No derived result changed, so no downstream unit (> 208) is materially affected by this verification.

## Side observations (non-blocking)

The stage card's "Mathematica audit: none yet" text-lag (the `.wl` exists and passes) is the known card-text-lag class deferred to P4-51 with the rest of batch VI.1; the audit did not raise it as a separate finding and it is outside the scripts-only scope. Non-blocking; nothing to add.

## Verdict justification

`verified`. The lone finding F1 is non-script and resolved by the orchestrator's re-run: the refreshed SymPy `.txt` now carries the canonical `STAGE 208` banner (stale `191` eliminated everywhere), with all 20 SymPy residuals `= 0` and all 27 Mathematica checks `PASS` / no `FAIL`, both exits 0, and both outputs re-dated newer than their scripts. No directive or diff existed because none was warranted, matching the auditor's explicit disposition. The audit verdict holds on the refreshed artifacts — the `.wl` is genuinely INDEPENDENT (r_grad `Solve`/`Reduce`+uniqueness vs posit; cross-weight `Series`/`Coefficient` extraction vs posit), reconciliation is 16/16 MATCH, and no constants are pinned. No regressions; `material_change: false`.
