---
unit_id: 059
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 059

## Per-finding outcomes

### F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Classification:** resolved

**What changed:**
Per the captured diff (`exec_logs/stage_059_diff.patch`), exactly two lines in
`scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py` changed:
- py:3 docstring `Moving-Throat PDE — Stage 42 SymPy audit.` → `Moving-Throat PDE — Stage 59 SymPy audit.`
- py:124 `print("\nStage 42 audit passed.")` → `print("\nStage 59 audit passed.")`

Both edits are strip-the-number identical to HEAD: only the integer label changed, 2-digit
format preserved, no surrounding text altered. Confirmed against the current file (py:3 now
reads `Stage 59`).

**Assessment:**
The change matches the directive's "required change" exactly. No collateral edits: the diff
touches only the two named lines (and the regenerated `.txt` output, which is expected).
The DEFERRED cross-refs were correctly LEFT UNTOUCHED — verified in the live file:
- py:6 still reads `... on the Stage-41 branch interval,`
- py:9 still reads `... using the Stage-39 Omega_Pe series.`
- py:75 still reads `# Weak-coupling expansion using the exact Stage-39 Omega_Pe series.`
The already-correct banner py:45 (`STAGE 59`) was not padded, per directive. No assertion,
symbol, or numeric expression changed — the diff hunks are confined to the docstring header
and the closing print string.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 59 — OPERATOR-SELECTED RESIDUAL BOUNDS` (banner now matches script)
- `Omega^2 linear coefficient = 0` (residual against `(4-pi)/pi` slope holds)
- `Xi_fail*DeltaInf saturates at Pe_star diff = 0` / `Xi_suff*Delta0 ... = 7.24...E-71` (< 1e-40 tol)
- `Stage 59 audit passed.` (closing label now corrected)

**Mathematica:** exit=0. Notable lines:
- `STAGE 059 — OPERATOR BRANCH RESIDUAL BOUNDS`
- `PASS: Xi_suff - Xi_fail (ordered)`, `PASS: Xi_fail*DeltaInf ...`, `PASS: Xi_suff*Delta0 ...`, `PASS: Omega^2 linear coefficient`
- `Omega_Pe^2 linear coefficient = -1 + 4/Pi` (independent limit-of-derivative method, agrees)
- `Stage 059 Mathematica audit passed.`
The `Limit::alimv` warning is benign (assumption on the limit variable ignored), as noted by the auditor.

**Output freshness:** Both refreshed `.txt` transcripts carry post-fix run dates of
2026-06-05T12:12 (sympy) / 12:15 (mathematica), newer than the script edits, and now show the
corrected `STAGE 59`/`Stage 59 audit passed.` labels. Outputs are fresh.

## Material-change assessment

`material_change`: false. The only edits are number-only label/string corrections (docstring
and closing print). No derived value, symbol, assertion, tolerance, or numeric result changed;
all eight assertions and both engine transcripts are unchanged in substance. No downstream unit
can depend on a comment/print label, so nothing > 059 is affected.

## Side observations (non-blocking)

The Mathematica `.txt` saturation lines print `0``49.09...` (Mathematica's backtick precision
suffix on a machine-zero `0`); this is benign formatting and matches the `PASS` lines. Not a finding.

## Verdict justification

The single low-severity `stale_output` finding (F1) is fully resolved: the two unambiguous
self-labels were corrected number-only with format preserved, exactly the three deferred
cross-references (`Stage-41` at py:6, `Stage-39` at py:9 and py:75) were left untouched, the
banner was not padded, and no math/assertion/value was altered. Both engines re-run to exit 0
with all checks passing and outputs refreshed to match the corrected labels. Verdict: verified;
material_change false.
