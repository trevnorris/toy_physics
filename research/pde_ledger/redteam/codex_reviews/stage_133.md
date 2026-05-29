---
stage: 133
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.py, moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.txt, moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl, moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.txt, stage_133.md, stage_133.tex]
---
# Codex review — stage 133

## What the edit was supposed to do
The directive required replacing the Mathematica hand-ansatz construction, which duplicated SymPy's `cCoeff`/`aCoeff` algebra, with an independent `DSolveValue` derivation of the scalar D/N mouth-layer solution. The downstream checks were supposed to remain: ODE residual, both boundary conditions, the mouth derivative kernel \(u'(0)/G\), and the static-shell limit. The purpose was to make Mathematica an independent derivation of the kernel used in the paper's two-channel fixed-point law.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:54 | DSolve-derived `u` satisfies scalar ODE | Yes; the residual is recomputed from the returned solution, not formed as `X-X`. | Yes, supports the D/N boundary-layer response behind \(\mathcal S(\Pi,\kappa)\). | PASS |
| 2 | moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:55 | `u(0) = 0` | Yes; it substitutes the independently solved profile into the Dirichlet condition. | Yes, checks the D side of the D/N layer. | PASS |
| 3 | moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:56 | `u'(1) = 0` | Yes; it differentiates the solved profile and evaluates the Neumann condition. | Yes, checks the N side of the D/N layer. | PASS |
| 4 | moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:57 | `u'(0)/gSrc` equals target kernel | Yes; `sKernel` comes from `DSolveValue`, while `sTarget` is a separate closed-form formula. | Yes, this is the scalar kernel in the quoted two-channel fixed-point law. | PASS |
| 5 | moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:63 | static-shell limit `S(Pi,0)=1` | Yes; the series limit is compared against the constant `1`. | Yes, checks the retained normalized static limit of the kernel. | PASS |

## Bottom line
PASS. I tried to break the fix by looking for reused hand coefficients, shared construction of both sides, and a copied SymPy ansatz in the Mathematica path: `cCoeff`/`aCoeff` are absent, `u` is obtained by `DSolveValue`, and the kernel target is only used after the independent solve. I also manually re-derived \(u'(0)/G\) from the D/N ODE and got the same \(\mathcal S(\Pi,\kappa)\), including the static limit \(1\). The saved Mathematica output shows `0` plus `PASS` for all five checks, and the SymPy output corroborates the same scalar kernel identity.
