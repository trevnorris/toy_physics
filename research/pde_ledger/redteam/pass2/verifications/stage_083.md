---
unit_id: 083
batch: III.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 083

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy (`scripts/...sympy_audit.py:64-75`): the two tautological `expect_zero` "defining-equation residual" blocks (old py:71-81, of the form `denom·(numer/denom) − numer ≡ 0`) are deleted along with the false "non-tautological" comment. Replaced by two `expect_close` numeric anchors comparing the closed forms `Delta0_F1`/`DeltaInf_F1` against external literals `1.73302079021525e-4` (atol 1e-16) and `2.01447565540522e-2` (atol 1e-15), under an honest comment "Numeric anchors for the Family-1 closed forms reported in the stage notes."
- Mathematica (`mathematica/...audit.wl`): the "Independent BVP derivation" comment block (old wl:62-71) is replaced by an honest "Redundant structural checks ... the independent numeric anchors later in the script pin the reported values" (wl:62-64). The byte-identical `A_F1 independent vs closed-form` tautology (old `aF1Indep = (kappaF1+(Pi/2)^2)/...`, since `(Pi/2)^2 ≡ Pi^2/4`) and its `expectApprox` are fully deleted (wl:86-91 now only checks the genuine `yRoot` residual). The `Omega(Pe)` comment (wl:127-131) is corrected to state plainly it is "tautological for the typed definition ... not as independent verification." The two `Delta` defining-equation `expectZero` residuals remain but are now honestly relabeled as redundant structural checks; the genuine numeric anchor battery (wl:154-164) is retained.

**Assessment:**
Correct and complete. The directive's preferred minimal path was taken on the `.wl` (correct misleading comments + delete the byte-identical `A_F1` tautology + keep the real numeric anchors), which exactly matches the offered "EITHER delete/correct comments ... OR perform genuine DSolve" choice. Adversarially: the new SymPy anchors compare the closed form against FIXED EXTERNAL literals (matching notes:81/83 and the `.wl`'s own pins at wl:154-155), not self-derived from `Delta0_F1` — a `cosh↔sinh` paste error in the closed form would shift its numeric value and blow the diff past the 1e-16/1e-15 tolerance, flipping PASS→FAIL. The literal precision (~15 sig figs) versus the tolerance is consistent: exec log shows actual diffs of 1.6e-19 and 4.0e-17, both comfortably inside tolerance and far below what any real perturbation would produce. SymPy now has a genuine anchor on `Delta_0`/`Delta_inf` where previously it had none. No tautology survives unlabeled.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
A reusable `expect_close` helper was added (`scripts/...sympy_audit.py:32-36`). Nine new numeric assertions were added for the stage's headline deliverables: four `Pe` windows (py:118-121) against literals `96.5285247264385`/`11220.5441626259`/`22.0062226330754`/`2558.01892349205`, and five `zeta` window/ceiling values (py:156-160) against `2.46622291347846`/`2.46752913273870`/`2.44257571477179`/`2.46752736855058`/`2.46752922945601`. Tolerances copy the `.wl` battery (1e-10..1e-7 for Pe, 1e-12 for zeta).

**Assessment:**
Correct and complete. The SymPy engine now asserts the stage's actual Output (previously print-only). Adversarially non-tautological: each `Pe`/`zeta` value is COMPUTED from the closed forms and operator coefficients (`Pe_minus_chi = Xi_chi_coeff * Delta0_F1`, `zeta_F1 = A_F1*Omega^2`) and asserted against a FIXED external literal — perturbing `Theta_chi_coeff` or the closed forms (as the directive's self-test predicts) would shift the computed value outside tolerance and fail. `grep -c expect_close` = 11 (2 Delta + 4 Pe + 5 zeta), meeting the directive's ≥11 bar. Exec log shows all nine new window checks PASS with diffs in the 1e-15..1e-11 range. No deliverable value moved — printed outputs are byte-identical to the report's reconciliation table.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Delta_0(F1) numeric anchor diff = 1.6008...E-19` (anchor passes, well inside 1e-16)
- `Pe_+^(chi) numeric check diff = 1.0389...E-11` (inside 1e-7)
- `zeta_max^(F1) numeric check diff = 2.1796...E-15` (inside 1e-12)
- All deliverable values unchanged, e.g. `zeta_max^(F1)= 2.46752922945601223332958450157`.

**Mathematica:** exit=0. 17 PASS lines, including `PASS: Delta_0(F1) numeric check`, the four `Pe` and five `zeta` checks, and `PASS: zeta_F1 monotone increasing`. No `aF1Indep`/"A_F1 independent" line remains (the tautology is gone). The two `Delta` defining-equation residuals still report `= 0` but are now honestly labeled redundant.

**Output freshness:** confirmed. `.py` and `.wl` mtimes = 1780695509; committed `sympy.txt` = 1780696024, `mathematica.txt` = 1780696044 — both outputs regenerated AFTER the scripts. Committed SymPy `.txt` contains 11 anchor/check lines; committed `.wl` `.txt` has 17 PASS lines and 0 occurrences of the deleted `A_F1 independent`/`Independent BVP` strings. Committed outputs match the exec logs.

## Material-change assessment

`material_change`: false. The edits add assertions and correct/delete misleading comments only; no formula, coefficient, or output value changed. Every deliverable (`Delta_0`, `Delta_inf`, four `Pe`, five `zeta`, ceiling, `Pi/C_mix`) is byte-identical to the pre-fix output and to the auditor's reconciliation table (21/21 match). No downstream unit's inputs are perturbed.

## Side observations (non-blocking)

- The two `.wl` `Delta` defining-equation `expectZero` residuals (wl:70-84) are still structurally tautological but are now HONESTLY labeled as redundant structural checks, and real numeric anchors (wl:154-155) cover the same constants — this matches the directive's accepted minimal path, so it is not a defect. The `Omega` residual is similarly honest now.
- The `Theta_chi_coeff`/`Theta_J_coeff` operator selectors and the `136900` prefactor remain unanchored upstream (flagged by the existing in-script `TODO(upstream-anchor)` comment, pre-existing, not part of either finding). Not blocking.

## Verdict justification

Both findings are fully resolved. F1: every tautology is either deleted (SymPy residuals, `.wl` `A_F1` byte-identical check) or honestly relabeled (`.wl` `Delta`/`Omega` residuals, with real numeric anchors carrying the actual coverage), and the misleading "non-tautological" / "Independent BVP derivation" comments are corrected. SymPy gains genuine external-literal anchors on `Delta_0`/`Delta_inf` that would catch a `cosh↔sinh` paste error. F2: nine new `expect_close` window asserts make the SymPy engine actually verify the stage's headline Output, all against fixed external literals and all passing with tiny residuals. Both engines exit 0, outputs are fresh and match the logs, and no deliverable value moved — material_change is false.
