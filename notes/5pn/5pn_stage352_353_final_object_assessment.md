# Stage 352–353 — Final Numerical Object Assessment

## Goal

Decide whether the current file stack already contains the **last missing numerical object**: the PDE-selected orbit-lock point, rather than only its exact theorem surface.

## What is already numerically located

The explicit Family-1 support/source branch **is** numerically located on the refreshed exact `\Lambda_{EM}` geometry branch. The current reduced operator family gives

- `\Lambda_\ell = L/\ell \approx 36.9497315424`
- `\kappa \approx 2457.5087899001`
- `\zeta_{\max} \approx 2.4675297457`
- on the `\chi`-weighted wall-depth extraction,
  - `Pe_* \approx 11155.7265863`
  - `\zeta_{\rm phys}(Pe_*) \approx 2.4675296479`
  - `\rho_{\alpha,\max} \approx 3.4675296479`

So the explicit support/source side is not only algebraically under control; it is numerically located and safely above the canonical isotropic demand.

## What the current files still do **not** numerically locate

The coherent-branch finish packet is still only present **symbolically** in the current scripts/notes. The exact actual-branch quantities are already compiled as

- `R_tr = (\chi_0 + \delta_U + 1) / ((\chi_0+1)(\delta_U+1))`
- `\epsilon = \epsilon_W (9\delta_U+11)/(11(\delta_U+1))`
- `R_target = \Lambda (1-\epsilon_\eta)(1-\epsilon)^2 / (Z_W(1+\chi_0)^2)`
- `N_Q = 1/\chi_Q` on the natural source-map branch

but these remain formulas in the branch variables, not a numerically selected branch point.

The decisive script verdict from the current stack is still:

> *“No numerical PDE-selected point is present yet in the notes/scripts.”*

So the final unlocated object is specifically the **orbit-lock / coherent placement point** — equivalently the realized values of `(R_tr, R_target, \epsilon_\eta)` or the realized drift verdict `(\delta\ln R_tr, \delta\ln R_target, \delta\ln \epsilon_\eta)` on the completed moving-throat branch.

## Bottom line

We have already numerically located the explicit support/source branch, but **we have not numerically located the final PDE-selected orbit-lock point** from the present file stack alone.

So the remaining gap is not more reduced algebra. It is the actual numerical or analytic solve of the completed moving-throat branch strongly enough to read off the coherent placement variables.

## Deliverables in this turn

- `5pn_stage352_final_object_status.py`
- `5pn_stage352_final_object_status_output.txt`
- `5pn_stage352_353_final_object_assessment.md`
