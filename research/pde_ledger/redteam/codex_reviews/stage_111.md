---
stage: 111
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.py, moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.txt, moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl, moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.txt, stage_111.md, stage_111.tex]
---
# Codex review — stage 111

## What the edit was supposed to do
The directive identified the Mathematica audit as a line-by-line transliteration of the SymPy audit. It required adding an independent Mathematica re-derivation of `chi_Q^mix` that bypasses `lambdaMix`'s `l5` extraction and instead computes the pole's `z^5` contribution directly from the geometric series. The existing four checks were to remain intact, with a fifth `expectZero["chi_Q^mix routes agree", chiMix - chiMixAlt]` added and shown passing.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl:75 | `chi_Q^mix routes agree` | Yes. It compares `chiMix` from `l0,l5` extracted out of `lambdaMix` against `chiMixAlt` from a separate pole-series coefficient route. | Yes. It independently corroborates the mixed-pole odd-normalization expression used by the stage's compensated-branch check. | PASS |

## Bottom line
PASS. I tried to break the added check as a tautology, but `chiMixAlt` is not assigned from `chiMix`, `l5`, or `l0`; it recomputes the pole's `z^5` coefficient from `sigma/(1 - kappa*z^2 - I*gamma*z^5)` and uses the explicit `L0 = -3 - sigma` route. A sign error in the pole subtraction, a wrong `z^5` coefficient extraction, or dropping the `/I` normalization would make line 75 fail. The assumptions only require real symbols and `sigma != -3`, so they do not smuggle in the `sigma = 0` absence result. The saved Mathematica transcript shows the independent route as `(3 - 27*gamma*sigma)/(3 + sigma)` and the new residual `chi_Q^mix routes agree = 0` followed by `PASS`; the SymPy transcript and the existing Mathematica checks still corroborate the even-branch and odd-normalization identities.
