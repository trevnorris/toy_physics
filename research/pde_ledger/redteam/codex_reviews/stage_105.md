---
stage: 105
reviewer: codex
verdict: FINDINGS
findings_count: 2
files_reviewed: [stage_105.md, stage_105.tex, moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py, moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt, moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl, moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt]
---
# Codex review — stage 105

## What the edit was supposed to do
F1 required the Mathematica audit to stop being a line-by-line port of the SymPy audit and derive \(\chi_Q=1\) through a distinct algebraic path, with different variable names and a Reduce-based uniqueness check. It also required the deformed branch to be derived from the operator identity rather than by directly assigning \(Y=-3/\Lambda\). F2 identified a stage-label paper/script mismatch; the current files use Stage 105 labels, matching the file path and paper label.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:37 | canonical \(\sigma_Q=4a^5/(27c_s^5)\) | yes | yes, DtN coefficient normalization | PASS |
| 2 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:45-47 | retarded branch \(\omega^2,\omega^4,\omega^5\) coefficients | yes | yes, outgoing DtN fingerprint | PASS |
| 3 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:49-55 | solve outgoing match for \(\chi_Q=1\) | yes | yes, central paper claim | PASS |
| 4 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:64-66 | deformed normalized branch coefficients | yes | yes, deformation coefficient relation | PASS |
| 5 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:42 | canonical \(\sigma_Q=4a^5/(27c_s^5)\) | yes | yes | PASS |
| 6 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:49-52 | unfactored ratio equals partial-fraction form | only for numerator algebra; shared denominator errors pass | supporting check only | PASS |
| 7 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:62-64 | retarded branch coefficient checks | yes | yes | PASS |
| 8 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:69-72 | Reduce-based \(\chi_Q=1\) check | not cleanly; malformed extraction still simplifies to PASS | intended to, but does not prove uniqueness as written | FINDING |
| 9 | moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:89-92 | polynomial-inverse deformed coefficients | yes | yes, but directive naming constraint is violated | FINDING |

## Findings
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:69`
- **What's wrong:** `chiReduce = Reduce[cRet5/I == aThroat^5/(27*cSound^5), chiQ, Reals];` does not include the declared positive assumptions. The saved output confirms this: Mathematica emits `ReplaceAll::argt` and prints a disjunction with degenerate `aThroat == 0` branches before reporting `chi_Q - 1 = 0`. That is not the directive's requested unique real solution proof; it is a malformed `ToRules`/`ReplaceAll` expression later simplified under `$Assumptions`.
- **What it should be:** Reduce the assumed problem directly and assert equivalence to `chiQ == 1`, e.g. `Reduce[aThroat > 0 && cSound > 0 && cRet5/I == aThroat^5/(27*cSound^5), chiQ, Reals]` followed by an explicit check that the result is `chiQ == 1`.

### R2 — other
- **Where:** `moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:80,87`
- **What's wrong:** The directive explicitly required the substrings `lamDef` and `yDef` to disappear from the Mathematica file. They remain as prefixes in `lamDeformed = -3 + z^2/3 + z^4/9 + I*xiQ*z^5/9;` and `yDeformedSeries = Expand[yAnsatz /. ySolved];`. This is not a mathematical failure of the coefficient solve, but it is not faithful to the applied directive and would fail the stated verifier substring check.
- **What it should be:** Rename those symbols to names without the forbidden substrings, such as `lambdaPoly` and `seriesFromInverse`.

## Bottom line
The SymPy checks and most Mathematica coefficient checks are mathematically anchored to the Stage 105 claim and are corroborated by the saved outputs. The blocking issue is the Mathematica \(\chi_Q\) proof: the output itself shows the Reduce result is not a clean unique-solution derivation, and the script reaches PASS only after simplifying a malformed replacement expression. Combined with the directive-forbidden `lamDef`/`yDef` substrings still present in the `.wl`, this stage cannot receive a trustworthy PASS.
