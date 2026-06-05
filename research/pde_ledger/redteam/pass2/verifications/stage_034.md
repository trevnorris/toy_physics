---
unit_id: 034
batch: II.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-04T23:20:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 034

## Per-finding outcomes

### F1 — stale_output (residual self-labels + stale transcripts)

**Classification:** resolved

**What changed:**
A single file edited, `scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`, two label-only string hunks (per `stage_034_diff.patch`):
- `.py:3` — module docstring `Moving-throat PDE Stage 17 SymPy audit.` → `Moving-throat PDE Stage 34 SymPy audit.`
- `.py:94` — `print("All Stage 17 checks passed.")` → `print("All Stage 34 checks passed.")`

Confirmed by reading the current file: line 3 reads `Stage 34`, line 94 reads `All Stage 34 checks passed.`. The `.wl` was correctly not touched (it carried no residual self-label; its staleness resolves via the orchestrator re-run only).

**Assessment:**
The edit matches the directive exactly. The diff is purely two label-string substitutions `17`→`34` with ZERO change to any equation, value, variable, assertion, or `expect_zero` target. No collateral edits anywhere in the diff. Not assertion-related, so non-tautology concerns do not apply. The fix correctly closes the only finding (cosmetic label staleness); no math was at stake.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 34.1 — EXACT SOFTENING-DEPTH SECULAR REDUCTION` (canonical banner; was `STAGE 17.1` in stale transcript)
- `alpha(lambda=A-x) - alpha(x) = 0` / `s(lambda=A-x) - s(x) = 0` / `N(lambda=A-x) - N(x) = 0` (residuals unchanged)
- `dalpha/dx - manifestly positive form = 0` / `s_x * d alpha / dx - 1 = 0` / `lambda-form vs x-form support loading = 0`
- closing `All Stage 34 checks passed.` (canonical), no FAIL.

**Mathematica:** exit=0. Notable lines:
- `STAGE 034 — SOFTENING-DEPTH NORMAL FORM` (canonical banner)
- `PASS: alpha(lambda=A-x) - alpha(x)`, `PASS: s(...)`, `PASS: N(...)`, `PASS: dalpha/dx - manifestly positive form`, `PASS: s_x * d alpha / dx - 1`, `PASS: lambda-form vs x-form support loading` — all six PASS, every `= 0` residual intact.
- closing `Stage 034 Mathematica audit passed.`

**Output freshness:** confirmed regenerated post-fix. Output `.txt` mtimes are 2026-06-04 23:13:28 (both engines), newer than `.py` (23:07:49) and `.wl` (2026-06-03 15:59:11). The refreshed SymPy transcript now reads canonical `STAGE 34.x` banners and `All Stage 34 checks passed.` (the precise stale-token failure mode is gone).

## Material-change assessment

`material_change`: false. The edit is label-only (docstring + closing print string). No derived result, value, or assertion changed; no downstream unit depends on a banner/print string. Every emitted symbolic form and every `= 0`/`PASS` residual is byte-for-content identical to the pre-fix math.

## Side observations (non-blocking)

None. The diff is minimal and exactly matches the directive's named file:line edits.

## Verdict justification

The single low-severity `stale_output` finding is fully resolved: Codex applied exactly the two named label-only edits with no collateral changes (diff confirms `17`→`34` strings only), both engines re-ran clean at exit 0 with all twelve assertions passing and no FAIL, the refreshed transcripts read canonical `STAGE 34.x`/`STAGE 034` and `All Stage 34 checks passed.`, and output mtimes postdate the scripts. No regressions, no math touched. Verdict: verified, material_change false.
