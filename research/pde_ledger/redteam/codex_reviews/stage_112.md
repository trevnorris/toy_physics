---
stage: 112
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [stage_112.md, stage_112.tex, moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py, moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.txt, moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl, moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.txt, moving_throat_pde_stage112_hybrid_robin_mixed_compensation.md, moving_throat_pde_stage092_dynamic_geometry_obstruction.md]
---
# Codex review — stage 112

## What the edit was supposed to do
The directive held the label mismatch as a paper-misalignment item pending user resolution, and required a real Mathematica independence improvement for the former transliteration finding. The requested mathematical fix was to add a Stage-92 linearized cross-check on the nontrivial compensated branch, using the branch-selection data \((b,a_0,a_5)=(0,3\sigma_W,-\sigma_W\gamma_W)\) and the preservation condition \(a_0/3+9a_5=0\) to rederive \(\gamma_W=1/9\). The existing chi-based branch checks were supposed to remain as a complementary route, not the only derivation.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:62 | `a0Def - 3*sigma == 0` | Yes; wrong branch or constant deformation normalization would fail. | Yes; checks the Stage-92 \(a_0=3\sigma_W\) data used for the paper's branch. | PASS |
| 2 | moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:63 | `a5Def + sigma*gamma == 0` | Yes; wrong odd sign or subtraction of the canonical \(1/9\) would fail. | Yes; checks the Stage-92 \(a_5=-\sigma_W\gamma_W\) data. | PASS |
| 3 | moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:64-66 | solve `a0Def/3 + 9*a5Def == 0` and assert `gammaFromLinear == 1/9` | Only generically; it silently divides out the \(\sigma=0\) factor. | Not over the declared assumptions; the paper's iff needs \(\sigma_W\ne0\) or a degenerate exception. | FINDING |

## Findings
### R1 — symbol_assumption_error
- **Where:** `moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:64-66`
- **What's wrong:** `gammaFromLinear = FullSimplify[gamma /. First[Solve[a0Def/3 + 9*a5Def == 0, gamma, Reals]], Assumptions -> $Assumptions];` followed by `expectZero["independent: gamma_W from a_0/3 + 9 a_5 = 0", gammaFromLinear - 1/9];` asserts a generic solution as an iff result. From the immediately preceding verified data, the condition is \(\sigma(1-9\gamma)=0\), not \(1-9\gamma=0\), under the script's assumptions. The script assumes `sigma != 1`, but not `sigma != 0`; at \(\sigma=0\), the side-channel amplitude vanishes, the preservation condition is identically true, and \(\gamma\) is unconstrained. The saved output corroborates that Mathematica printed `gamma_W from linearized preservation = 1/9` and `PASS`, but that is only the generic branch of `Solve`, not the full declared real domain.
- **What it should be:** Either add an explicit nontrivial-loading assumption such as `sigma != 0` before deriving the iff, or split the check: verify that \(\sigma=0\) is the degenerate no-side-channel case where \(\gamma\) is not determined, and verify \(\gamma=1/9\) only under \(\sigma\ne0\).

## Bottom line
The new Mathematica block is mostly a real independent check: I tried breaking the \(a_0\) and \(a_5\) assertions by changing the branch, sign, and canonical subtraction in the algebra, and those would not survive. The load-bearing gamma assertion, however, proves only the generic \(\sigma_W\ne0\) case while the script's declared assumptions still allow \(\sigma_W=0\), where the paper's “iff \(\gamma_W=1/9\)” statement is false without a degenerate-case qualifier.
