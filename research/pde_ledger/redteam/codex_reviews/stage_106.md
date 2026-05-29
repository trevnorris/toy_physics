---
stage: 106
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py, moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.txt, moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl, moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.txt, stage_106.md, stage_106.tex]
---
# Codex review — stage 106

## What the edit was supposed to do
The directive required replacing tautological branch-identity checks with checks against the independent hardcoded target literals. It also required re-authoring the Mathematica audit along a structurally independent one-pole series path, including the omega-series evidence and canonical closure. Finally, it required adding a first-order `Delta_Q` sensitivity check. The unresolved F1 paper-misalignment warned that the stage card still asks for higher-odd and outgoing DtN fingerprint checks unless the paper is explicitly changed to delegate them upstream.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:64 | `K0_target*K4_target - 4*K2_target^2` | Yes | Yes, target-literal consistency | PASS |
| 2 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:68 | `Gamma5_target - 9*sqrt(K2_target^5/K0_target^3)` | Yes | Yes, target-literal consistency | PASS |
| 3 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:88 | `Delta_Q` first-order slope | Yes | Yes, reduced failure sensitivity | PASS |
| 4 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:93 | `Delta_Q` zeroth coefficient | Yes | Yes, canonical limit | PASS |
| 5 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:41 | `sigma_Q^can` normalization | Yes | Partial, one-pole setup only | PASS |
| 6 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:63 | omega^5 coefficient form | Yes | Partial; does not derive `chi_Q=1` or DtN fingerprint | FINDING |
| 7 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:74 | target-literal identities | Yes | Yes, target-literal consistency | PASS |
| 8 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:87 | `N_Q=1`, canonical branch | Yes | Yes, product closure | PASS |
| 9 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:104 | `Delta_Q` sensitivity checks | Yes | Yes, reduced failure sensitivity | PASS |
| 10 | moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:47 | higher-odd “verification” by printed coefficients | No, print-only | No pass/fail assertion | FINDING |

## Findings
### R1 — paper_misalignment
- **Where:** `moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:5`
- **What's wrong:** The script states that paper Checks (ii) and (iii) are “verified upstream,” but `stage_106.tex:21-25` still requires this stage to check higher odd terms and the outgoing `l=2` DtN fingerprint. The SymPy output contains no `z` expansion, no `Lambda_2^{out}`, no `Yhat_2^{out}`, and no higher-odd assertion.
- **What it should be:** Either the paper card must be changed to explicitly delegate those checks, or this stage must add real assertions for the Hankel/DtN expansion and the higher-odd onset.

### R2 — insufficient_verification
- **Where:** `moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:47`
- **What's wrong:** The comment says the omega series verifies that higher odd terms begin beyond omega^5, but lines 52-60 only print coefficients. If `omega7Coeff` were wrong, absent, or zero, the script would still pass.
- **What it should be:** Add assertions, e.g. coefficient checks for zero odd terms below omega^5 and an expected nonzero omega^7 coefficient.

### R3 — insufficient_verification
- **Where:** `moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:63`
- **What's wrong:** The Mathematica audit checks only that the one-pole omega^5 coefficient equals `I chiQ sigmaQcan/4`. It never matches this coefficient to the canonical `Gamma5`/`K0` relation or solves/derives `chi_Q = 1`; the result block instead says `chi_Q = 1` is a carry-in.
- **What it should be:** Assert the coefficient-to-`Gamma5` normalization relation required by F3/M1, then derive the canonical `chi_Q=1` closure rather than importing it.

## Bottom line
The replacement target-literal and `Delta_Q` assertions are not tautological and are corroborated by the saved outputs, but the stage still does not verify two paper-card checks. The most important problem is that the scripts substituted comments and printouts for the required DtN fingerprint and higher-odd assertions, while the paper card still assigns those checks to stage 106.
