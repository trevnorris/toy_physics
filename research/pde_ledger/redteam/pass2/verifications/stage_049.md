---
unit_id: 049
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T18:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 049

## Per-finding outcomes

### F1 — stale_output (unambiguous self-labels; number-only)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py:3` — docstring first line `Stage 32 SymPy audit.` → `Stage 49 SymPy audit.` (confirmed via file read and diff hunk @@ -1,6 +1,6 @@).
- `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:93` — `Print["Stage 32 Mathematica audit passed."]` → `Print["Stage 49 Mathematica audit passed."]` (confirmed via file read and diff hunk @@ -90,6 +90,6 @@).
- Both committed `.txt` outputs auto-refreshed by the orchestrator's independent re-run.

**Assessment:**
The change matches the directive's required change exactly. It is strictly label-only / strip-the-number: no assertion, helper, formula, assumption, or banner string was touched. The captured diff (`stage_049_diff.patch`) contains only the two single-line label edits — no collateral edits. The audit also noted the same single `Stage 32` self-label in each file (it cited the py docstring as line 2 vs. the directive's line 3; both point at the same docstring text line, which is now corrected). There is exactly one `Stage 32` self-label per file and both are fixed. The already-correct `STAGE 49` banners (py:51 / wl:36) were correctly left unpadded.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 49 — EXPLICIT D/N OVERLAP EXTRACTION OF THE PHYSICAL SUPPORT RATIO` (canonical banner)
- `microscopic coherent-support law = 0`, `exact twin-lane support ratio = 0`, `lowest twin half-wave = 0` (all residuals zero)
- No `STAGE 32` self-label anywhere in the log.

**Mathematica:** exit=0. Notable lines:
- `STAGE 49 — EXPLICIT D/N OVERLAP EXTRACTION OF THE PHYSICAL SUPPORT RATIO` (canonical banner)
- `PASS: microscopic coherent-support law`, `PASS: exact twin-lane support ratio`, `PASS: lowest twin half-wave` (all seven checks PASS)
- Final line `Stage 49 Mathematica audit passed.` (no longer `Stage 32`).

**Output freshness:** Both refreshed `.txt` outputs now print `STAGE 49` at output line 3 and the mathematica `.txt` final line reads `Stage 49 Mathematica audit passed.` Their content matches the post-fix exec logs verbatim (exec log dates 2026-06-05T12:08 / 12:13, post-edit). No stale `STAGE 32` banner remains in either committed transcript. Every result value (chi_n, k_n, I_n, I_n/I_0, zeta_phys, twin stiffness, x, zeta_twin, zeta_0) is identical to the audit-recorded values — only the banner/final-line label changed.

## Material-change assessment

`material_change`: false. The edit changes only human-readable stage-number self-labels in a docstring and a `Print` statement, plus the auto-refreshed output transcripts. No derived result, assertion, or carry-forward formula changed; downstream units cannot depend on a label string. No units > 049 need to be marked stale on account of this change.

## Side observations (non-blocking)

- The auditor's structural-parallelism note on the `.wl`/`.py` (borderline transliteration) is acknowledged but was correctly not raised as a finding — the load-bearing overlap integral is independently evaluated by each CAS and both agree. Nothing to act on.

## Verdict justification

The sole finding (F1, low-severity stale_output / numbering self-label) was applied exactly as the directive required: two number-only label corrections (`Stage 32` → `Stage 49`) in the py docstring and the wl final `Print`, with the committed `.txt` outputs refreshed to the canonical `STAGE 49` banner. The captured diff confirms no collateral or math edits; both engines exit 0 with all checks PASS; no `STAGE 32` artifact survives in source or outputs; all result values are unchanged. Verdict: `verified`, `material_change: false`.
