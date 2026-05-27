---
unit_id: 078
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 078 (v2)

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`:

- Line 47 now reads
  ```
  thetaSuffSym = (-(45 Cosh[111 Sqrt[5]/5] + 27 Sqrt[5] Sinh[111 Sqrt[5]/5]) / (2500 - 2500 Cosh[111 Sqrt[5]/5]));
  ```
  replacing the prior `thetaFailSym * (4.21495341569977*^-2 / 3.62605617972939*^-4)` decimal-ratio bootstrap. The Stage-75 symbolic closed form is now used directly.
- Lines 44-46 carry a new provenance comment naming Stage-75 sympy output line 21 as the source.
- Lines 48-50 replace the previous comment that overpromised a "verify their ratio chi:J matches the Stage-77 ratio" check (which did not exist) with the corrected statement that the chi/J values are adopted from Stage-77 at extended precision and their independent re-derivation belongs to the Stage-077 audit.

**Assessment:**
The edit matches the directive's required change. The formula `-(45 cosh(α) + 27 √5 sinh(α)) / (2500 - 2500 cosh(α))` with α = 111√5/5 is the exact Stage-75 Theta_suff/Pe_req closed form cited in the auditor's F1 self-test. The orchestrator's one-line syntax fix (wrapping the entire RHS in outer parens) is cosmetic — Mathematica precedence on a single division was already unambiguous, so the wrap does not change semantics. The downstream consumers (thetaSuffCoeff = N[thetaSuffSym, 30] on line 56; the four expectApprox targets on lines 75-82) continue to consume the new symbolic expression. The mathematica output file confirms `Pe_suff^(chi) numeric check diff = 0``27.71` (zero at ~28 digits) and `Pe_suff^(J) numeric check diff = 0``28.36`, both PASS — i.e., the new symbolic re-derivation reproduces the prior coefficient value to the working precision. The new assertion is non-tautological: if the closed form had a sign or transcendental-function error, the diff would not collapse to zero. Engine independence on Theta_suff is now genuine, no longer a literal-decimal echo of SymPy.

### F2 — paper_misalignment

**Classification:** resolved

**What changed:**
Three banner/docstring relabels, applied by the orchestrator directly per the `## Approved by user (2026-05-27)` block (direction (b): script-side only; notes/ preserved as historical record):

- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:3` docstring: `Stage 61 SymPy audit.` -> `Stage 078 SymPy audit.`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:23` banner: `STAGE 61 — FAMILY-1 BRANCH VERDICT` -> `STAGE 078 — FAMILY-1 BRANCH VERDICT`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:32` banner: `STAGE 061 — FAMILY-1 BRANCH VERDICT` -> `STAGE 078 — FAMILY-1 BRANCH VERDICT`

**Assessment:**
All three relabels are present in the current file state and visible in the diff. Both saved output files (sympy line 3 and mathematica line 3) now display `STAGE 078 — FAMILY-1 BRANCH VERDICT`, confirming the runtime banners reflect the rename. The notes/ file is correctly left untouched per the user's resolution. No mathematical content was modified.

## Exec log assessment

**SymPy:** exit=n/a. `redteam/exec_logs/stage_078_sympy.log` is not present in the exec_logs directory (only `stage_078_diff.patch` exists). The orchestrator's directive note states "Both engines exit 0", and the canonical output `scripts/output/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.txt` (mtime 02:15) was regenerated after the .py file mtime (02:05). Notable lines from the output:

```
STAGE 078 — FAMILY-1 BRANCH VERDICT
Pe_suff^(chi) / lambda_mu^2 = 96.528524726438575954
Pe_suff_J < Pe_suff_chi  : True
Pe_fail_J < Pe_fail_chi  : True
Pe_suff_chi < Pe_fail_J  : True
```

**Mathematica:** exit=n/a. `redteam/exec_logs/stage_078_mathematica.log` is not present. The canonical output `mathematica/output/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.txt` (mtime 02:24) was regenerated after the .wl file mtime (02:20). Notable lines:

```
STAGE 078 — FAMILY-1 BRANCH VERDICT
Pe_suff^(chi) numeric check diff = 0``27.71431433552255
PASS: Pe_suff^(chi) numeric check
Pe_suff^(J) numeric check diff = 0``28.35642450198806
PASS: Pe_suff^(J) numeric check
PASS: Pe_suff^(chi) < Pe_fail^(J) (window overlap)
Stage 078 Mathematica audit passed.
```

All four `expectApprox` deltas are zero at ~26-28 digits of precision (effectively zero at working precision), and all five `expectTrue` inequalities PASS. The `peSuffChiTarget = N[thetaChiCoeffNum / thetaSuffSym, 30]` is now computed from the new symbolic closed form, so the PASS is non-tautological — it would not collapse to zero if the new symbolic form did not reproduce the Stage-75 Theta_suff coefficient.

**Output freshness:** confirmed. Sympy: .py at 02:05, .txt at 02:15 (output newer than source). Mathematica: .wl at 02:20, .txt at 02:24 (output newer than source). Both outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

No downstream-visible numeric or symbolic result changes. The new `thetaSuffSym` evaluates to the same `0.0421495341569977...` coefficient as the prior decimal-ratio bootstrap (to ~28 digits), so all four Pe ratios, all five inequality checks, and the saved output values are byte-identical to the pre-fix state at any reasonable working precision. The banner relabels are purely cosmetic. The four Pe ratios are paper-side verdict numerics, not inputs to any later stage's derivation, so no downstream unit's chain of inference is touched.

## Side observations (non-blocking)

- The mathematica and sympy exec log files were not captured into `redteam/exec_logs/` for this stage (only the diff patch is present). The output files in `scripts/output/` and `mathematica/output/` are however fresh and consistent with the orchestrator's "both engines exit 0" note. Future stages may benefit from explicit exec-log capture for the verifier's audit trail, but this does not block verification because the saved .txt outputs carry all the assertion-level PASS lines and the closing `Stage 078 Mathematica audit passed.` banner.
- The orchestrator's one-line syntax fix (outer parens around the entire RHS on line 47) is a benign wrapping — Mathematica precedence on a flat division was already unambiguous, and the wrap is purely a line-continuation hygiene measure. Worth recording only because it deviates from the verbatim multi-line directive text; no semantic effect.
- This file overwrites a stale `stage_078.md` verification from a prior (v1) audit cycle that scored four findings; v2 only has two findings (F1 mathematica_transliteration, F2 paper_misalignment), and both are resolved here.

## Verdict justification

Both F1 and F2 are fully resolved. F1's `thetaSuffSym` is now a direct symbolic re-derivation from the Stage-75 closed form `-(45 cosh(α) + 27 √5 sinh(α)) / (2500 - 2500 cosh(α))` with α = 111√5/5 (no SymPy-output decimal on the RHS), and the overpromising comment was rewritten to honestly describe the chi/J adoption. F2's three banner relabels are in place across both `.py` and `.wl`. Both saved engine outputs are fresh, banners read `STAGE 078 — FAMILY-1 BRANCH VERDICT`, and all assertions (4 `expectApprox` + 5 `expectTrue` in Mathematica; 3 inequality checks in SymPy) PASS. The new Mathematica check `Pe_suff^(chi) numeric check diff = 0``27.71` would not collapse to zero if the new symbolic form did not match the Stage-75 coefficient, so the assertion is non-tautological. No regressions visible in the diff. No new findings warranted.

stage 078 (v2): verified
