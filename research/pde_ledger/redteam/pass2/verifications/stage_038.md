---
unit_id: 038
batch: III.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 038

## Per-finding outcomes

### F1 — stale_output (self-label) [label-only]

**Classification:** resolved

**What changed:**
The captured source diff (`stage_038_diff.patch`) touches exactly one file —
`scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py` —
and exactly two lines in the docstring:
- line 3: `moving_throat_pde_stage21_dimensionless_continuum_placement_sympy_audit.py` → `…stage038…`
- line 5: `Stage 21 SymPy audit:` → `Stage 38 SymPy audit:`

The `Stage-20 continuum formulas` cross-ref (docstring line 6, shown as diff context) is
LEFT untouched, as the directive required (it is a cross-ref to upstream stage 037 / pre-renumber
"20", owned by the deferred dedicated pass). The `.wl` was not touched (already canonical).

**Assessment:**
Strictly label-only. Only the stage-number self-tokens (`21`→`038`/`38`) changed; surrounding
docstring text is byte-identical. No equation, value, variable name, or assertion was altered —
confirmed by the diff being confined to the two docstring lines, with zero changes anywhere in the
seed forms, substitution map, `expect_zero` calls, or the derivative/sign blocks. The CROSS-ref to
037 was correctly preserved. The directive's `## Applied: F1` ("deviation: none") matches the diff
exactly. Edit is correct and complete.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 38 — DIMENSIONLESS CONTINUUM PLACEMENT AUDIT` (line 8) — canonical banner now present.
- All four dimensionless-map residuals `= 0` (lines 23–26), product relation both forms `= 0`
  (lines 35–36), nine derivative factorizations `= 0` (lines 42–50), nine sign-coefficient checks
  `= 0` (lines 51–59).
- `STAGE 38 AUDIT COMPLETE` (line 66); `# exit_code: 0` (line 70).

**Mathematica:** exit=0. Notable lines:
- `STAGE 038 — DIMENSIONLESS CONTINUUM PLACEMENT` (line 8) — canonical banner; the prior
  mixed-epoch state (`STAGE 021` header vs. `Stage 038` footer) is resolved — header and footer
  now both read 038.
- Every check carries an explicit `PASS:` (lines 24–82): all four maps, both product forms, nine
  factorizations, nine sign coefficients.
- `Stage 038 Mathematica audit passed.` (line 89); `# exit_code: 0` (line 90).

**Output freshness:** Both refreshed transcripts are dated `2026-06-05T08:09:34/40-06:00` (log
headers), newer than the source edit, and carry the canonical banners. Every prior `= 0` residual
and `PASS` from the report's assertion inventory (A1–A8 / B1–B8) is present and unchanged. Outputs
are fresh and faithful to the current scripts.

## Material-change assessment

`material_change`: false.

The only edit is two docstring label tokens in the SymPy script; the `.wl` was untouched and both
transcripts changed only in their banner/footer strings. No derived result, constant, or assertion
changed, so no downstream unit (> 038) can depend on anything that moved.

## Side observations (non-blocking)

None. (The report's A8 "partial" sign-coefficient note and the borderline second-engine
structural parallelism were already weighed and dismissed by the auditor; they are not findings and
are out of scope for this label-only verification.)

## Verdict justification

The single low-severity `stale_output` finding is fully resolved. The source diff is strictly
label-only — only the stale self-number tokens `21`→`038`/`38` in the SymPy docstring changed, the
`Stage-20` cross-ref to upstream 037 was correctly preserved, and the `.wl` was left alone. Both
refreshed transcripts show canonical banners (`STAGE 38` / `STAGE 038`), retain every `= 0`
residual and `PASS`, and both engines exit 0. No regression, no math change. Verdict: verified.
