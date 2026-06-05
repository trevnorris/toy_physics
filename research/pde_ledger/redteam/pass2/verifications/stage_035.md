---
unit_id: 035
batch: II.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-04T23:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 035

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
Per `stage_035_diff.patch`, the only file touched is the SymPy script
`scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py`,
with exactly two one-token edits:
- `.py:3` (module docstring): `Moving-throat PDE Stage 18 SymPy audit.` → `Moving-throat PDE Stage 35 SymPy audit.`
- `.py:112` (final print): `print("All Stage 18 checks passed.")` → `print("All Stage 35 checks passed.")`

The directive's `## Applied: F1` block matches: `files_changed` = the sympy script, summary = relabel Stage 18→Stage 35 in docstring + final print, `deviation: none`.

**Assessment:**
The edit is exactly what the directive asked for and nothing more. The diff is purely
label-string changes: each hunk changes only the literal number `18`→`35` inside a comment
docstring and a `print(...)` string argument. No equation, value, variable, `expect_zero`
target, assertion, limit, series, or domain assumption is touched. The two hunks are the
total extent of the change (no collateral edits anywhere in the repo per the patch). This
is a cosmetic, non-load-bearing finding and it is fully closed. The output-staleness itself
(the original substance of F1) resolves on re-run — confirmed under "Output freshness" below.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 35.1 — EXACT D/N DIMENSIONLESS SHAPE FUNCTION` (canonical 2-digit banner, no residual "STAGE 18")
- `F - closed D/N form = 0`, `dF/dxi - manifestly positive form = 0`, `F(0,delta) - 1 = 0`, `softening constant - closed form = 0`, `alpha_req - closed D/N form = 0`, `alpha_crit - closed form = 0`, both near-onset series `= 0` — every residual still zero.
- `All Stage 35 checks passed.` (closing line now canonical) followed by `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines:
- `STAGE 035 — DIMENSIONLESS NORMALIZATION LOCUS` (already canonical; .wl was untouched, correctly).
- All 11 `PASS:` lines present (`PASS: F - closed D/N form`, `PASS: R_target - Pi^2 A NQ/(8 beta0)`, `PASS: dF/dxi - manifestly positive form`, `PASS: F(0,delta) - 1`, `PASS: softening constant - closed form`, `PASS: near-softening asymptotic coefficient`, `PASS: alpha_req - closed D/N form`, `PASS: alpha_crit - closed form`, `PASS: near-onset F series through O(xi^2)`, `PASS: alpha_req near-onset series through O(xi^2)`), no FAIL.
- `Stage 035 Mathematica audit passed.` then `# exit_code: 0`.

No FAIL anywhere in either engine; both exit 0.

**Output freshness:** confirmed. The committed transcripts were regenerated post-fix:
- sympy `.txt` mtime 2026-06-04 23:13:28 > sympy `.py` mtime 2026-06-04 23:09:29 (and the diff-patch banner-only change shows in it).
- mathematica `.txt` mtime 2026-06-04 23:13:28 > mathematica `.wl` mtime 2026-06-03 15:59:11.
Both `.txt` files now carry the canonical banners; the original stale-output condition is gone.

## Material-change assessment

`material_change`: false.

The only edit is two cosmetic label strings (a docstring comment and a `print` argument).
No derived symbolic/numeric result changed — every residual that read `0` before still
reads `0`, and the closed forms for `F`, `dF/dxi`, `C(delta)`, `alpha_req`, `alpha_crit`,
and both near-onset series are byte-identical in content. Nothing downstream can depend on
a SymPy banner string, so no downstream unit is affected.

## Side observations (non-blocking)

None. The auditor's value-reconciliation augmentation (9 deliverables, 0 misaligned) remains
consistent with the refreshed transcripts; the post-fix outputs reproduce every reconciled value.

## Verdict justification

The single low-severity `stale_output` finding (F1) is fully resolved. The diff is exactly
the directive's prescribed two label-only edits (`Stage 18`→`Stage 35` at `.py:3` and `.py:112`)
with zero math/value change. Both refreshed transcripts read canonical (`STAGE 35.x` / `STAGE 035`),
close with `All Stage 35 checks passed.` / `Stage 035 Mathematica audit passed.`, retain every
`= 0` residual and every `PASS:` line, show no FAIL, and exit 0. Output mtimes are newer than their
scripts, so staleness is cured. No regressions, no new findings, `material_change: false`. Verdict: verified.
