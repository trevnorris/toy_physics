---
unit_id: 085
batch: III.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 085

## Per-finding outcomes

No findings were raised by the auditor (verdict `clean`, `findings_count: 0`, no directive emitted). There are no per-finding entries to evaluate. The only edit reaching the working tree was an orchestrator-direct banner relabel.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Pi_tr/C_mix - alpha_req/alpha_mix = 0`
- `zeta_req loading form - rho_alpha form = 0`
- `unblocked limit = 0`

**Mathematica:** exit=0. Notable lines:
- `PASS: Pi_tr/C_mix - alpha_req/alpha_mix`
- `PASS: zeta_req loading form - rho_alpha form`
- `PASS: unblocked limit`
- `Stage 085 Mathematica audit passed.`

All 8 assertions (A1–A8 / B1–B8 from the auditor's inventory) clear in both engines with residual = 0. The banners in both saved outputs read `STAGE 085 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS`, confirming the re-run happened after the relabel.

**Output freshness:** confirmed.
- `scripts/output/...sympy_audit.txt` mtime 2026-05-27 10:24 > script mtime 2026-05-27 10:15.
- `mathematica/output/...mathematica_audit.txt` mtime 2026-05-27 10:26 > script mtime 2026-05-27 10:15.
Both `.txt` outputs are post-fix.

## Material-change assessment

`material_change`: false.

The only edits in `git diff HEAD` for stage 085 are two single-line banner string changes: `"STAGE 68"` → `"STAGE 085"` in the `.py`, `"STAGE 068"` → `"STAGE 085"` in the `.wl`. No symbols, assumptions, expressions, assertions, or values were altered. Downstream units cannot be affected by a printed banner label.

## Side observations (non-blocking)

None. The auditor's self-test notes already cover tautology/positivity/branch-cut/missing-deliverable/constant attacks.

## Verdict justification

Auditor returned `clean` with 0 findings, so there is nothing substantive to verify against. The orchestrator-direct banner relabel is confirmed as the sole substantive change in the diff (cosmetic), both engines re-ran cleanly with the new banner stamped in the output headers, and the `.txt` mtimes are strictly newer than the script mtimes. All 8 assertions pass in each engine with residual zero. Verdict: verified.
