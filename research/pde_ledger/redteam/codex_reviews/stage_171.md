---
stage: 171
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [stage_171.md, stage_171.tex, moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py, moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.txt, moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl, moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.txt]
---
# Codex review — stage 171

## What the edit was supposed to do
F1 was supposed to repair the Mathematica second-engine check for the factor-heavy \(Z\) and \(N\) obstruction bundles. The directive's first proposed repair was itself tautological, so the orchestrator reworked it: keep the paper's collected closed-form bundle targets, compare them against Mathematica's `D[]`-derived total differentials, and add independent `Series`-coefficient routes for the same bundles. F2 corrected traceability labels from Stage 154 / Stage 170 to Stage 171; it did not add a mathematical assertion.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:123 | `D[]`-derived `zCombExact` equals collected \(dZ_2+dZ_0/9\) target | Yes; wrong collected coefficients leave uncancelled `dQ`, `dS`, `dG`, or `dDelta` residuals | Yes, outgoing-transfer bundle in the \(K_A\) decomposition | PASS |
| 2 | moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:130 | `Series` route for \(z_2+z_0/9\) equals same collected Z target | Yes; this is not `zCombExact - zCombExact` and would fail on a wrong target | Yes, same \(K_A\) outgoing-transfer bundle through an independent route | PASS |
| 3 | moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:138 | `D[]`-derived `nCombExact` equals collected \(dN_0+p_0dZ_0\) target | Yes; wrong `dQ`, `dP`, or `dDelta` slopes leave residuals | Yes, outgoing-transfer bundle in the \(G_A\) decomposition | PASS |
| 4 | moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:143 | `Series` route for \(n_0+p_0z_0\) equals same collected N target | Yes; it extracts a first-order coefficient from perturbed expressions, not from the summed differential | Yes, same \(G_A\) outgoing-transfer bundle through an independent route | PASS |

## Bottom line
PASS. I tried to break the edited checks by changing the \(Z\) bundle's \(1/(9\delta)\) `dQ` contribution, the \(-q/(9\delta^2)\) `dDelta` contribution, the \(N\) bundle's \(2p/\delta^2\) `dP` slope, and the \(p_0q/\delta^2\) `dDelta` term; each attack would leave an independent residual coefficient, so these are load-bearing rather than tautological. The saved Mathematica output shows all four edited bundle checks as `= 0` with `PASS`, and both saved outputs show the Stage 171 banner. I read the directive, the stage card, both scripts, and both saved transcripts; the nonzero denominator assumptions are justified by the rational formulas and do not smuggle in the result.
