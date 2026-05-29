---
stage: 163
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py, moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.txt, moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl, moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.txt, stage_163.md, stage_163.tex]
---
# Codex review — stage 163

## What the edit was supposed to do
The directive required only the Mathematica audit to stop being a pure transliteration of the SymPy script. It had to add an implicit-function derivation of \(g_-'(r)\) from \(F(g,r)=0\), and a Series-based derivation of the microscopic \(\delta r\), \(\delta g\), and \(\delta_\perp\) formula from the parent ratios. The existing checks were to remain in place, and the saved Mathematica output had to show the new PASS lines plus the pre-existing PASS lines.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl:40 | `gPrime` equals implicit-function derivative from `fComp` | yes: a sign/error in `gPrime` or in the partial-derivative ratio leaves a nonzero residual | yes: verifies the tangent direction used to define the off-family normal coordinate | PASS |
| 2 | moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl:76 | Series-derived `delta r` equals hand logarithmic form | yes: changing any parent-ratio exponent/sign makes the residual nonzero | yes: verifies the parent-action perturbation entering microscopic \(\delta_\perp\) | PASS |
| 3 | moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl:77 | Series-derived `delta g` equals hand logarithmic form | yes: changing any \(g_q,g_s,K_s,K_q\) exponent/sign makes the residual nonzero | yes: verifies the parent-action perturbation entering microscopic \(\delta_\perp\) | PASS |
| 4 | moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl:79 | Series-route microscopic `delta_perp` equals explicit formula | yes: altering the normal projection or any coefficient in the explicit formula gives a nonzero residual | yes: directly checks the card’s explicit microscopic off-family drift scalar | PASS |

## Bottom line
PASS. I tried to break the added checks by treating them as possible \(X-X\) identities, sign-copy checks, or checks of adjacent quantities only. The implicit derivative check uses partial derivatives of `fComp`, not `D[gMinus,r]` on both sides; the Series route derives `deltaRSubst` and `deltaGSubst` from parent ratios before comparing to the hand forms; and the final Series-route \(\delta_\perp\) check compares that independently derived perturbation against the explicit microscopic formula. The Mathematica saved output corroborates all four new checks with `= 0` and `PASS`, and the pre-existing SymPy and Mathematica outputs still show the stage’s transport and tangent/normal checks passing.
