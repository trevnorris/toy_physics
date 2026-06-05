---
unit_id: 072
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T14:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 072

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
No source edit. The orchestrator's independent re-run regenerated both committed transcripts. The committed `.txt` files now carry the canonical banners:
- `scripts/output/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.txt:3` reads `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES` and line 31 reads `STAGE 072 THEOREM LEDGER`.
- `mathematica/output/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.txt:3` reads `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES`.

The HEAD→working diff of the two committed transcripts contains exactly three `+/-` content lines, all banner lines: `STAGE 055`→`STAGE 072` (mma L3), `STAGE 55`→`STAGE 072` (sympy L3), `STAGE 55 THEOREM LEDGER`→`STAGE 072 THEOREM LEDGER` (sympy L31). A grep for any non-banner `+/-` line returns none — every residual/PASS/closed-form line is byte-identical to HEAD.

**Assessment:**
Correct and complete. The directive required a transcript refresh only, with no source edit. Verified:
- Source `.py` and `.wl` are UNCHANGED vs HEAD: `git diff --stat HEAD` for both is empty, and the captured `stage_072_diff.patch` is 0 bytes. Source mtimes remain `2026-05-26 13:08`; both `.txt` mtimes are `2026-06-05 13:58`, i.e. newer than their scripts (freshness restored).
- The directive's verification criteria are met line-for-line (sympy L3, sympy L31 ledger header, mma L3 banners all canonical).
- No collateral edit: the only `+/-` content lines in the refreshed transcripts are the three banner lines; all eight residual asserts and all PASS lines are unchanged. No tautology was introduced (no assertion was touched at all — source is identical to HEAD).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES` (L8), `STAGE 072 THEOREM LEDGER` (L36) — canonical banners.
- `Delta0 shell leading-order matches full Delta0 = 0`, `DeltaInf shell leading-order matches full DeltaInf = 0`, `shell fail asymptotic = 0`, `shell suff asymptotic = 0` (L20-23) — all four shell residuals zero.
- `compression fail asymptotic = 0`, `compression suff asymptotic = 0` (L32-33) — comp residuals zero.
- `# exit_code: 0` (L48).

**Mathematica:** exit=0. Notable lines:
- `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES` (L8) — canonical banner.
- `PASS: Delta0 shell leading-order matches full delta0`, `PASS: DeltaInf shell leading-order matches full deltaInf`, `PASS: shell fail asymptotic`, `PASS: shell suff asymptotic` (L29,31,35,37).
- `PASS: compression fail asymptotic`, `PASS: compression suff asymptotic` (L55,57); `Stage 072 Mathematica audit passed.` (L59); `# exit_code: 0` (L60).
- The `Limit::alimv` warnings (L23,25,43) and `DeltaInf shell leading-order ratio = 2/Sqrt[5] + (5 + 2*Sqrt[5])^(-1)` (L27, exactly 1) are benign and pre-existing — the auditor already disposed of them; the subsequent `expectZero[ratio − 1]` PASS (L30-31) confirms the ratio resolves to 1.

**Output freshness:** confirmed. Both committed `.txt` mtimes are `2026-06-05 13:58:32`, newer than the `2026-05-26 13:08:42` script mtimes. Exec logs are dated `2026-06-05T13:53/13:55`. The committed transcript content matches the exec-log content (banner-only delta vs HEAD).

## Material-change assessment

`material_change`: false.

This was a transcript-banner refresh with zero source edits. No derived result, closed form, threshold surface, or asymptotic constant changed. Downstream units cannot be affected. No `Upsilon_fail/suff`, `V0_fail/suff^2`, or asymptotic value moved.

## Side observations (non-blocking)

None. The two pre-existing benign items (the `Limit::alimv` warnings and the unsimplified-but-exactly-1 `2/Sqrt[5] + 1/(5+2 Sqrt[5])` ratio) were already noted and correctly disposed of by the auditor; they are unchanged by this refresh and are not blocking.

## Verdict justification

The sole finding (F1, stale_output, low severity) is fully resolved: the orchestrator's re-run overwrote both committed transcripts so their banners now read the canonical `STAGE 072` (sympy L3 + ledger L31; mathematica L3), while the source `.py`/`.wl` remain byte-identical to HEAD (empty diff, 0-byte patch, source mtimes unchanged). The refreshed `.txt` mtimes post-date the scripts, restoring freshness. Both engines exit 0 with all eight asymptotic residuals zero and all PASS lines present and unchanged. No regression, no collateral edit, no tautology introduced. Verdict `verified`, `material_change: false`.
