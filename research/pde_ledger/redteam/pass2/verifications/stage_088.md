---
unit_id: 088
batch: III.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 088

## Per-finding outcomes

### F1 — stale_output (stale self-label, numbering)

**Classification:** resolved

**What changed:**
SymPy docstring self-label updated in `scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`:
- line 3: `moving_throat_pde_stage71_loading_ratio_from_minimal_module_sympy_audit.py` → `moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`
- line 5: `SymPy audit for Stage 71.` → `SymPy audit for Stage 088.`

The diff (`exec_logs/stage_088_diff.patch`) shows exactly these two lines changed (+2/-2), and the live file confirms the new state (verified by Read of lines 1-35).

**Assessment:**
Correct and complete. The change matches the directive's required change byte-for-byte. `grep -n "stage71\|Stage 71"` on the SymPy script returns nothing (exit 1). The canonical markers all still read 088: the filename, the docstring filename line (3), the docstring audit line (5), and the in-script banner at line 29 (`STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE`). The `stage085` cross-references are correctly UNTOUCHED: py:112 (`scripts/moving_throat_pde_stage085_*`) and wl:87 (`in the stage 085 Mathematica audit files`) both still present — these are legitimate cross-refs to the upstream Stage-085 product identity, out of scope per the deferred numbering plan. No collateral edits. The change is docstring-only and non-functional: no assertion, constant, value, or substitution was touched.

## Exec log assessment

**SymPy:** exit=0 (directive `## Applied: F1` records Codex re-ran the script and confirmed exit 0 with identical transcript; the change is docstring-only and Python docstrings are not printed). The committed output `.txt` is byte-identical to before — `git diff HEAD` for the stage-088 output files shows zero changes, confirming the transcript is unchanged as required.

**Mathematica:** exit=n/a. The `.wl` was not edited (no finding targeted it); its `stage 085` cross-reference is correct and was deliberately left in place.

**Output freshness:** The committed `.txt` outputs are correctly UNCHANGED (and must be — a docstring edit is non-functional and cannot alter the printed transcript). `git diff --stat HEAD` shows only the `.py` changed (2 insertions, 2 deletions); both output `.txt` files are clean. Neither output file ever contained `stage71`/`Stage 71` (grep returns nothing), so there was no stale label leaking into committed output.

## Material-change assessment

`material_change`: false. The edit is confined to two SymPy docstring lines (file self-label and audit-description line). No derived result, numeric value, assertion, substitution, or constant changed. The output transcript is byte-identical. No downstream unit can depend on a docstring; nothing > 088 is affected.

## Side observations (non-blocking)

None. The audit's own special-care verifications (non-vacuous extraction, no early comment-close, 9-assertion ↔ 9-PASS count, independent Mathematica derivation) were already documented in the report and are not in scope for this verify pass; nothing in the diff disturbs them.

## Verdict justification

The single low-severity finding (F1, stale self-label) was applied exactly as directed. The stale `stage71`/`Stage 71` self-label is gone (grep clean); filename, docstring, banner, and paper card all read 088; the `stage085` cross-references are correctly preserved; no math/value/assertion changed; the committed outputs are byte-identical (non-functional change); and no paper or notes files were touched. `verified`, material_change false.
