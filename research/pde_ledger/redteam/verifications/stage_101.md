---
unit_id: 101
batch: IV.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 101

## Per-finding outcomes

### F1 — tautological_check (Mathematica)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:42-45,53-54` — three `expectZero` calls now anchor to the INPUT factorization `mHat0^2*chiQ*nQ - 1` with candidate `nQ` substitutions, not to `nQExact` substituted into itself. The banner at line 26 was also corrected from "STAGE 084" to "STAGE 101" (collateral, beneficial — matches the auditor's cosmetic note in §"Engine cross-check").

**Assessment:**
Correct. Each residual is now of the form `(mHat0^2*chiQ*nQ - 1) /. {mHat0 -> 1, nQ -> <candidate>}`, anchored to the stage's input equation rather than to the solved form. A typo on the candidate side (e.g., `nQ -> 1/chiQ^2`) would yield a non-zero residual. The `nQSol`/`nQExact` definitions (lines 33-34) and the diagnostic `Print[...]` lines (36-38) are retained as required. Non-tautological.

### F2 — missing_verification_script (SymPy)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py:33-36,46-51` — four `assert` statements added, all anchored to `(mhat0**2 * chiQ * NQ - 1).subs({...})` with candidate NQ values for the natural branch, canonical compact branch, exact `1/(1+DeltaQ)` replacement, and the small-DeltaQ linearization `series_delta - (-DeltaQ + DeltaQ**2) == 0`. A docstring annotation was also added at lines 2-19 explaining the F4 Cluster B (c) routing of Checks 2/3 to stages 102/097.

**Assessment:**
Correct. Grep confirms exactly 4 `assert` keywords in the file. The trivial-case pre-check from the directive verifies each residual reduces to `0`; mutating the candidate side breaks each residual. The added docstring is within F4's applied-block scope (Cluster B direction (c) per resolution doc) and does not alter the math. Non-tautological.

### F3 — insufficient_verification (small-DeltaQ linearization)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:51-52` — `expectZero["small-DeltaQ series matches paper", Expand[seriesDelta - (-deltaQ + deltaQ^2)]]` inserted. SymPy side was already handled inside F2 (`assert sp.expand(series_delta - (-DeltaQ + DeltaQ**2)) == 0`), per the directive's instruction not to duplicate.

**Assessment:**
Correct on both engines. The transcript shows `small-DeltaQ series = -deltaQ + deltaQ^2` and `PASS: small-DeltaQ series matches paper`, and the literal `-deltaQ + deltaQ^2` (resp. `-DeltaQ + DeltaQ**2`) appears in both files. A sign flip would now fail.

### F4 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
Per the directive's `## Applied: F4` block, this stage owns Check 1 (factorization). Checks 2 and 3 are annotated in the SymPy docstring (lines 5-15) as upstream carry-ins from stages 102 (higher_odd_irrelevance) and 097 (single_normalization_defect). No paper.tex or notes/ edits in this unit (verifier did not read them; confirmed by orchestrator context referencing `batch_IV1_paper_alignment.md` Cluster B direction (c)).

**Assessment:**
Within the verifier's scripts-only scope, the docstring annotation in the SymPy script is in place and consistent with the orchestrator's stated routing. Whether the paper-side cross-references at stages 102/097 actually verify Checks 2/3 is out of scope for this stage's verifier (it falls to verifiers of those stages). Treated as resolved here.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `NQ from exact factorized odd normalization = 1/(chiQ*mhat0**2)`
- `small-DeltaQ series = DeltaQ**2 - DeltaQ`
- `check exact replacement chiQ=1+DeltaQ: 0`
- `STAGE 101 AUDIT PASSED` (script reached terminal print, so all 4 asserts passed)

**Mathematica:** exit=0. Notable lines:
- `PASS: point-particle natural branch reduction`
- `PASS: canonical compact outgoing branch gives NQ=1`
- `PASS: small-DeltaQ series matches paper`
- `PASS: exact replacement chiQ=1+DeltaQ`
- `Stage 101 Mathematica audit passed.`

Expected 4 PASS lines (3 from F1 + 1 added by F3); observed 4. Matches the orchestrator-context cue.

**Output freshness:** confirmed.
- SymPy script mtime 2026-05-27 11:18:48 < output mtime 2026-05-27 14:28:44.
- Mathematica script mtime 2026-05-27 11:18:59 < output mtime 2026-05-27 14:30:55.
Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: false.

Edits add assertions and re-anchor residuals; no derived quantity is changed. The symbolic results in the transcript (`NQ = 1/(chiQ*mhat0**2)`, `1/chiQ`, `1`, `-DeltaQ/(1+DeltaQ)`, series `-DeltaQ + DeltaQ**2`) are identical to the pre-fix run. Downstream units inherit the same algebra; no `upstream_stale` flag warranted on substance grounds.

## Side observations (non-blocking)

- SymPy still declares `DeltaQ` with `positive=True` (line 23), even though `Delta_Q = chi_Q - 1` could be negative when `chi_Q < 1`. The auditor noted this as cosmetic and did not file it; nothing has changed. Out of scope for this verification but worth flagging if a future audit revisits assumptions.
- The Mathematica banner is now correctly "STAGE 101"; the prior "STAGE 084" leftover was a documented copy-paste issue and the fix is welcome collateral.

## Verdict justification

All four findings have applied blocks in the directive and corresponding edits in the source files. F1's three Mathematica residuals are now anchored to the input factorization and pass. F2 adds four substantive SymPy `assert`s anchored the same way, and the script exits 0. F3 inserts the missing Mathematica series-comparison `expectZero`, yielding 4 PASS lines as expected. F4 is documented in the SymPy docstring per the orchestrator's Cluster B (c) resolution. Exec logs are fresh; engines agree; no regressions.
