---
stage: 117
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [stage_117.md, audit_directive_117.md, stage_117.tex, moving_throat_pde_stage117_outlet_core_status_sympy_audit.py, moving_throat_pde_stage117_outlet_core_status_sympy_audit.txt, moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl, moving_throat_pde_stage117_outlet_core_status_mathematica_audit.txt]
---
# Codex review — stage 117

## What the edit was supposed to do
The directive required applying only F4: replace the capstone classification table's hardcoded boolean columns with booleans computed from the earlier residuals and branch solves. The compensated survivor specifically had to be wired to the section-5 residual `delta_core - delta_core_expected`, while keeping the final survivor-set assertion. F1-F3 were explicitly blocked: the Mathematica transliteration and the `kappa_c`/`gamma_c` normalization tautologies were not to be mechanically fixed without user-supplied upstream derivations.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:184 | SymPy capstone booleans wired to section results | Yes; bad roots/residuals change the booleans | Partly: reduced outlet classification, not independent D/N realization | PASS |
| 2 | moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:216 | SymPy nontrivial survivor set | Yes; changed booleans alter `nontrivial_survivors` | Yes, for the reduced-model survivor classification | PASS |
| 3 | moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:136 | Mathematica capstone booleans mirror SymPy | Only as the same reduced algebra; not independently | No independent second-engine corroboration | R1 |
| 4 | moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:166 | Mathematica nontrivial survivor set | Only as the same reduced algebra; not independently | No independent second-engine corroboration | R1 |

## Findings
### R1 — transliteration
- **Where:** `moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:136`
- **What's wrong:** `evenOkScale = TrueQ[Sort[beta /. betaSolutions] === {-1, 1}];` begins a capstone block that is still a line-by-line port of the SymPy block, including `nontrivialCompensated = TrueQ[Normal[Series[deltaCore - deltaCoreExpected, {z, 0, 5}]] === 0]`. This was already called out as F1 and blocked in the directive. The saved Mathematica output corroborates that the copied check passes, but it is not an independent derivation; a wrong `deltaCoreExpected` target or carried-forward `kappa0/gamma0` substitution copied into both engines would pass in both.
- **What it should be:** Either derive the Mathematica side through an independent path from the Schur complement / D/N BVP / bare mixed-channel premises, or explicitly mark the `.wl` file as duplicate corroboration rather than a second-engine audit.

## Bottom line
The SymPy F4 wiring is no longer the old hardcoded survivor table: the booleans are tied to prior solves and the section-5 residual, and the saved output shows the expected five-row table. The stage still cannot receive a clean PASS because the Mathematica “second engine” remains the known transliteration from F1, so its passing output does not independently protect the paper-side claim from a copied wrong reduced target.
