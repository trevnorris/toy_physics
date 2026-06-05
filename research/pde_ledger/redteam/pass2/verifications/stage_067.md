---
unit_id: 067
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 067

This was a refresh-only stage. The auditor's sole finding (F1) was a low-severity
`stale_output`: both committed `.txt` transcripts carried stale stage-number banners
(`STAGE 50` / `STAGE 050`) that disagreed with the canonical `STAGE 067` banners in
the current scripts. The auditor wrote no script edit and no directive; Codex was not
invoked. The fix is purely the orchestrator's independent re-run regenerating the two
committed outputs with canonical banners. The directive file
(`directives/stage_067.md`) is absent — confirmed; this is expected and consistent
with the "no source edit" scope. Verification proceeds from the report, the git diff
of working tree vs HEAD, and the captured exec logs.

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
Only the two `output/*.txt` files changed; the source `.py` and `.wl` are byte-for-byte
unchanged vs HEAD (`git status --porcelain` empty for both source files; `git diff --stat
HEAD` empty for both). The output diffs touch exactly one line each, the L3 banner:
- `scripts/output/moving_throat_pde_stage067_sech_gaussian_sympy_audit.txt:3`:
  `STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK` (en-dash, "50")
  → `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` (canonical hyphen, "067").
- `mathematica/output/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.txt:3`:
  `STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` ("050")
  → `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` ("067").
No other line in either output changed.

**Assessment:**
Correct and complete. The refreshed banners now match the canonical banners emitted by
the unchanged source scripts (`banner("STAGE 067 …")` at sympy L53 / `banner["STAGE 067 …"]`
at mathematica L38). All numeric/symbolic results are preserved — confirmed by diffing the
committed outputs against the exec-log bodies (excluding headers): both are byte-identical
("committed == exec log body" for sympy and mathematica). The headline deliverables are
intact: `C_res^2 = 0.994418836451529348…`, `P_res = 1.00561248776057621695…`,
`1 - C_res^2 = 0.00558116354847065129…`. Zero `FAIL` lines in either output; 13 `PASS`
lines in the mathematica output and the sympy ledger emits its full results. The fix is
the canonical pass-2 script/output-band staleness remediation with no math change — exactly
as the report's "Required change" prescribed. No source-side collateral.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` (L8 — canonical banner)
- `N_(sigma sigma) integral - 2 w_f = 0`, `N_(phi phi) integral - w_g sqrt(pi/2) = 0`
- `C_res^2 = 0.994418836451529348706428351608877628170873348983716948813464`
- `P_res   = 1.00561248776057621695172301479763550405448504648609605997534`

**Mathematica:** exit=0. Notable lines:
- `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` (L8 — canonical banner)
- `PASS: C_res^2 numeric check`, `PASS: P_res numeric check` (cross-engine NIntegrate vs mpmath)
- `PASS: constructive-branch increase up to r_*`, `PASS: constructive-branch decrease after r_*`
- `Stage 067 Mathematica audit passed.` (final ledger label canonical)

**Output freshness:** confirmed. Committed output mtimes (`1780689512` for both `.txt`)
post-date the corresponding script mtimes (sympy `.py` `1779822672`, mathematica `.wl`
`1779822689`). The committed `.txt` bodies are byte-identical to the captured exec-log
bodies, so the saved outputs are the product of the orchestrator's post-fix re-run.

## Material-change assessment

`material_change`: false.

No source script changed; no derived value changed. The only delta is a cosmetic banner
label in two transcript files. No downstream unit can depend on this; nothing to flag.

## Side observations (non-blocking)

- The mathematica output's final line was already canonical (`Stage 067 … passed`) in the
  stale transcript, while its top banner was stale — corroborating the report's note that
  the banner was edited 50/050→067 after the `.txt` was last committed. The refresh now
  makes the top banner consistent with the rest of the transcript. Non-blocking.
- The numeric values reconcile fully with the auditor's Value-Reconciliation table
  (8 values, 0 misaligned). Not a finding.

## Verdict justification

The sole finding (F1, `stale_output`) is resolved exactly as prescribed: the orchestrator's
re-run refreshed both committed transcripts so their L3 banners now read the canonical
`STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK`, output mtimes post-date the scripts,
and every numeric/symbolic result and PASS line is unchanged (verified byte-identical to the
exec-log bodies). Source `.py`/`.wl` are unchanged vs HEAD; the diff is confined to the two
banner lines; no directive was needed and none exists. Both engines exit 0 with zero FAILs.
No regressions, no collateral, no math change. Verdict: `verified`, `material_change: false`.
