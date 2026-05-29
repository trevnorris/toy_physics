---
stage: 115
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage115_core_balance_compensation_sympy_audit.py, moving_throat_pde_stage115_core_balance_compensation_sympy_audit.txt, moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl, moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.txt, stage_115.md, stage_115.tex, stage_appendix_part04.tex]
---
# Codex review — stage 115

## What the edit was supposed to do
The directive required the Mathematica audit to stop being only a line-by-line transliteration of the SymPy route. It asked for an independent parent-overlap reparametrization using \(\mathfrak r=\lambda/\sqrt{K_sK_q}\), \(\mathfrak g=g_q\sqrt{K_s}/(g_s\sqrt{K_q})\), and \(1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2\). The load-bearing goal was to derive \(\sigma_c=\sigma_*\) from that parent-family branch, while preserving the existing balance and outlet fingerprint checks.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl:53 | parent-family residual equals balance equation times factor | yes; wrong factor, wrong \(\mathfrak g\), or wrong \(\sigma_c\) denominator leaves a nonzero residual | yes; anchors parent family to the card’s coupling balance | PASS |
| 2 | moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl:69 | explicit \(\mathfrak g_-\) is one Solve root of parent family | yes; wrong half-discriminant or wrong family equation fails | yes; checks the appendix/root form used to realize the branch | PASS |
| 3 | moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl:80 | substitute parent-derived \(g_q\) and verify \(\sigma_c=\sigma_*\) | yes; wrong branch translation, target, or core formula fails | yes; this is the compensated-branch conclusion quoted by the stage | PASS |

## Bottom line
PASS. I checked the new Mathematica assertions against the directive, the stage card’s coupling-balance claim, and the saved transcript. I tried to break the checks by inspecting whether the first residual was an \(X-X\) construction, whether Solve ordering made the root check vacuous, whether the final \(\sigma_c=\sigma_*\) substitution reused `gQBranch`, and whether the parent-ratio domain introduced an unjustified branch assumption. The checks are not tautological: the equivalence factor at line 56 is independently necessary, the root check compares Solve output to the closed-form parent root, and the final substitution uses `gQFromFrakMinus`, not the original `gQBranch`; the saved Mathematica output shows all three new checks at `0` with `PASS`, and the existing SymPy/Mathematica transcripts still corroborate the original stage checks.
