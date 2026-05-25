---
unit_id: 078
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 078

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
No direct edit at SymPy lines 41-44 or Mathematica lines 61-62; the directive itself stipulated that F1 is satisfied transitively by F4's new branch-verdict assertions and that the original ordering checks must be preserved. The original assertions are still present at `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:51-54` and `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:80-81`.

**Assessment:**
The directive's "Applied: F1" block matches what I see in the diff (no change to lines 41-44 of the .py or lines 61-62 of the .wl). The supplemental F4 assertions (verified separately below) do break the tautology: `Pe_suff_J < Pe_suff_chi`, `Pe_fail_J < Pe_fail_chi`, and `Pe_suff_chi < Pe_fail_J` each depend on the actual numeric magnitudes of both `Theta_chi` and `Theta_J` (cancellation no longer applies because the four `Pe_*` values appear in non-cancelling combinations across the inequalities). A corruption of either Theta would now be caught.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
Codex added the prescribed provenance comments at `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:26-29` and `:32-37`, naming the Stage-77 and Stage-75 upstream output files and lines, and recording the symbolic closed form for `Theta_fail/Pe_req`. Numeric values are unchanged.

**Assessment:**
The comments verbatim match the directive's required text. No collateral edits. The Mathematica half of F2 (replacing literal `expectApprox` targets with independently-computed values) was deferred to F3 per the directive and is verified there.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
At `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:37-53`, Codex replaced the four `SetPrecision[...]` literal-decimal coefficient assignments with: (i) an explicit symbolic closed form `thetaFailSym` in `Sinh[111 Sqrt[5]/5]`/`Cosh[111 Sqrt[5]/5]`, (ii) a scaled `thetaSuffSym = thetaFailSym * (4.21495341569977e-2 / 3.62605617972939e-4)`, (iii) high-precision (40-digit) numeric loads of the chi/J Theta values via `ToExpression["...`40"]`. At `:70-79`, the four `expectApprox` calls now target `peSuffChiTarget`/`peFailChiTarget`/`peSuffJTarget`/`peFailJTarget` computed from `thetaFailSym`/`thetaSuffSym` and the chi/J numerics, not the SymPy decimals.

**Assessment:**
The `Sinh`/`Cosh` closed form for `thetaFailSym` matches the form recorded in the SymPy provenance comment (and in the directive). All four `expectApprox` targets are now computed, not typed in. The .wl output file confirms all four checks PASS (`diff = 0.` for two, `diff = 0\`\`25.6...` and `0\`\`26.3...` for the other two — both are arithmetic-equivalent to zero at working precision).

Codex applied one explicit deviation: dropped the `100` factor from `thetaSuffSym` (the directive had `thetaSuffSym = 100 thetaFailSym * (4.21495e-2 / 3.626e-4)`, Codex wrote `thetaSuffSym = thetaFailSym * (4.21495e-2 / 3.626e-4)`). The deviation is correct and necessary: the SymPy script uses `Theta_suff_coeff = 4.21495e-2` (not 4.21495); the `100` factor in the directive would have made `thetaSuffSym ≈ 4.21` (an Upsilon-scale value) and broken `Pe_suff^(chi) numeric check` since `peSuffChiTarget = thetaChi/thetaSuffSym` would then be off by 100×. The SymPy provenance comment `Theta_suff = Upsilon_suff / 100` is internally consistent with `4.21495e-2` and confirms Codex's interpretation. The deviation is documented in the directive's "Applied: F3" block.

One residual weakness (non-blocking, noted as side observation): `thetaSuffSym` is derived as `thetaFailSym × (sympy_theta_suff/sympy_theta_fail)`, so it is not fully independent of the SymPy decimals — the *ratio* is imported. However, the absolute scale of `thetaFailSym` is independent (symbolic closed form), so a sign or transcendental-function error in the closed form would still be caught. This matches the directive's structure and is acceptable.

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:56-77`, Codex appended three branch-verdict assertions plus three `print` lines, verbatim matching the directive. At `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:83-86`, Codex appended three `expectTrue` calls before the existing `Print[""] / Print["...passed."] / Exit[0]` block.

**Assessment:**
The three SymPy assertions (`Pe_suff_J < Pe_suff_chi`, `Pe_fail_J < Pe_fail_chi`, `Pe_suff_chi < Pe_fail_J`) are non-tautological. The first two require `Theta_J < Theta_chi`; the third requires the actual numeric ordering of `Theta_chi/Theta_suff_coeff` vs `Theta_J/Theta_fail_coeff`. All three would fail if `Theta_chi` or `Theta_J` were corrupted in sign or magnitude. The Mathematica mirror is exact. Output files confirm all three pass in both engines:
- SymPy: `Pe_suff_J < Pe_suff_chi  : True`, `Pe_fail_J < Pe_fail_chi  : True`, `Pe_suff_chi < Pe_fail_J  : True`
- Mathematica: `PASS: Pe_suff^(J) < Pe_suff^(chi)`, `PASS: Pe_fail^(J) < Pe_fail^(chi)`, `PASS: Pe_suff^(chi) < Pe_fail^(J) (window overlap)`

## Exec log assessment

**SymPy:** exit=n/a. `redteam/exec_logs/stage_078_sympy.log` is missing. The canonical output `scripts/output/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.txt` is present, fresh (mtime 23:31 vs script mtime 23:28), and contains all the expected new lines, implying a successful run. Notable lines:

```
Pe_suff^(chi) / lambda_mu^2 = 96.528524726438575954
Pe_fail^(chi) / lambda_mu^2 = 11220.544162625905301
Pe_suff^(J)   / lambda_mu^2 = 22.006222633075413597
Pe_fail^(J)   / lambda_mu^2 = 2558.0189234920526360
Pe_suff_J < Pe_suff_chi  : True
Pe_fail_J < Pe_fail_chi  : True
Pe_suff_chi < Pe_fail_J  : True
```

**Mathematica:** exit=n/a. `redteam/exec_logs/stage_078_mathematica.log` is missing. The canonical output `mathematica/output/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.txt` is present, fresh (mtime 23:31 vs script mtime 23:30), and ends with `Stage 078 Mathematica audit passed.` All ten PASS lines are present (four `expectApprox`, two original `expectTrue`, three new F4 `expectTrue`, plus the closing banner). Notable lines:

```
Pe_suff^(chi) numeric check diff = 0.
PASS: Pe_suff^(chi) numeric check
Pe_fail^(chi) numeric check diff = 0``25.648956084907972
PASS: Pe_fail^(chi) numeric check
PASS: Pe_suff^(J) < Pe_suff^(chi)
PASS: Pe_fail^(J) < Pe_fail^(chi)
PASS: Pe_suff^(chi) < Pe_fail^(J) (window overlap)
```

**Output freshness:** confirmed. SymPy output mtime 23:31 > script mtime 23:28. Mathematica output mtime 23:31 > script mtime 23:30. Both `.txt` outputs were re-generated after Codex's edits.

## Material-change assessment

`material_change`: false.

No printed symbolic content or numeric value changed downstream-visibly. The four `Pe_*^*` numeric values (`96.528...`, `11220.544...`, `22.006...`, `2558.019...`) are byte-identical to the pre-edit outputs. New assertions and provenance comments do not alter any downstream-consumable derivation; they only add internal validation. The F3 deviation (no `100` factor) reflects an arithmetic correction to the directive itself, not a change to the operational coefficient values.

## Side observations (non-blocking)

- F3's `thetaSuffSym` independence is partial: `thetaSuffSym = thetaFailSym × (sympy_suff/sympy_fail)` imports the ratio of SymPy decimals. The directive itself prescribed this structure (`thetaSuffSym = 100 thetaFailSym * (sympy_suff/sympy_fail)`), so Codex did not introduce the weakness. A stronger independence would derive `thetaSuffSym` from its own Stage-75 closed form. Flagging for future hardening; not a blocking finding.
- The `100` factor removal in F3 is a substantive math correction to the directive that the auditor's directive arithmetic missed. The corresponding `Applied: F3 / deviation` block accurately records this. Worth keeping in mind when the orchestrator reviews the diff trail.
- Exec logs for stage 078 are missing from `redteam/exec_logs/`. Canonical `.txt` outputs are present and fresh, which is sufficient to confirm passing runs, but the orchestrator may wish to verify why log capture didn't run for this stage.

## Verdict justification

All four findings are `resolved`. Codex applied each finding as directed (F1 by transitive coverage from F4 per the directive's explicit instruction; F2 via the prescribed provenance comments; F3 by replacing literal decimals with a symbolic closed form and computed `expectApprox` targets; F4 by appending the three branch-verdict assertions in both engines). The one deviation Codex took (dropping the `100` factor in F3's `thetaSuffSym`) is mathematically correct and aligns the independent target with the operational `Theta_suff` scale; without it the new `expectApprox` would have failed by two orders of magnitude. Canonical refreshed outputs confirm all assertions PASS in both engines, no downstream numeric values changed, and the new assertions are substantive (non-tautological).

stage 078: verified
