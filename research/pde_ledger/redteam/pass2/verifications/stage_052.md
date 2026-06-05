---
unit_id: 052
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T18:10:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 052

## Per-finding outcomes

### F1 — stale_output (unambiguous filename self-label; number-only)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py:4` —
docstring filename label `moving_throat_pde_stage35_nontwin_asymmetry_threshold_sympy_audit.py`
→ `moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`. The captured
diff (`exec_logs/stage_052_diff.patch`) is exactly this single-line replacement and
nothing else; the on-disk file confirms line 4 now reads the `stage052` (3-digit) stem.
The committed SymPy `.txt` was refreshed (now prints `STAGE 52` banners at out lines 8/82
instead of the stale `STAGE 35`).

**Assessment:**
Correct and complete. The directive asked for exactly the `stage35`→`stage052` filename
self-label fix on py:4, in canonical 3-digit form, plus an output refresh — that is what
was applied. It is a strip-the-number, label-only edit: prose line 6 "Stage 52" and the
`STAGE 52` banners (py:42/124) were correctly left untouched per the directive's "do NOT
pad the already-correct banners" instruction. No assertion, symbol, or numeric expression
is touched. No collateral edits in the diff (only the .py line + the regenerated .txt).
This is identical-modulo-the-number to HEAD.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 52 — EXACT BRANCH-PRODUCT REGIME CLASSIFICATION` (banner now current)
- `dzeta_req/dPi - expected = 0`, `Delta_zeta - expected = 0`, `softening fraction - expected = 0`
- `solve(zeta_phys = zeta_req) for Omega0^2 - expected = 0` / `for Kphi0 - expected = 0`
All residual-zero identities hold; no math line differs from the pre-fix transcript content.

**Mathematica:** exit=0. Notable lines:
- `STAGE 052 — NON-TWIN ASYMMETRY THRESHOLD` (already canonical)
- Every check prints `PASS:` (boundary anchors, dZdPi independent path, Delta_zeta,
  both threshold solves, softFrac independent path, twin diagnostic, softening fraction)
- `Stage 052 Mathematica audit passed.`
The `.wl` was not edited (it was already canonical); it re-ran clean as the independent
cross-engine confirmation.

**Output freshness:** confirmed. Both `.txt` mtimes are 2026-06-05 12:22:11, newer than the
`.py` (11:52) and the `.wl` (2026-06-03). Outputs are post-fix.

## Material-change assessment

`material_change`: false. The only edit is a docstring filename string; no derived result,
assertion, symbol declaration, or numeric/symbolic value changed. No downstream unit can
depend on a comment string. Nothing downstream of 052 is affected.

## Side observations (non-blocking)

- The auditor's F1 also noted near-trivial `zeta_twin - 1` (KW/KW-1) checks and SymPy's
  `eps` declared only `positive` (omitting `eps<1`); the auditor judged both harmless
  (pure rational identities) and prescribed no change. Out of scope for this verify; no
  action needed.

## Verdict justification

The single low-severity stale_output finding was applied exactly as directed: a one-line,
number-only docstring self-label fix (`stage35`→`stage052`, 3-digit) plus an output refresh,
with no math, assertion, symbol, or value altered and no collateral edits in the diff. Both
engines re-ran to exit 0 with all PASS/residual-zero checks intact, and the `.txt` outputs
are confirmed fresh. material_change is false. Verdict: verified.
