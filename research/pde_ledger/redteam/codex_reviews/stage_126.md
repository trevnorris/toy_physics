---
stage: 126
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [audit_directive_126.md, verify_directive_126.md, stage_126.tex, moving_throat_pde_stage126_positive_source_families_sympy_audit.py, moving_throat_pde_stage126_positive_source_families_sympy_audit.txt, moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl, moving_throat_pde_stage126_positive_source_families_mathematica_audit.txt]
---
# Codex review — stage 126

## What the edit was supposed to do
The directive held F1 as a notes-side typo resolution, then required script fixes for F2-F4. F2 was supposed to add nonnegativity checks for the positive mouth-source family on `z in [0,L]`, `xi in [0,1]`. F3 was supposed to make the SymPy interval check raise on failure, and F4 corrected the banner from Stage 109 to Stage 126.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage126_positive_source_families_sympy_audit.py:51 | `min sigma_match on [0,L] = 0` | Yes for wrong endpoint/sign changes | Yes for the self-matched endpoint source | PASS |
| 2 | moving_throat_pde_stage126_positive_source_families_sympy_audit.py:55-66 | `xi=0`, `xi=1`, and corner positivity checks | Only for endpoint failures; an interior negative perturbation can pass | Only partially; does not check full `sigma_xi >= 0` on the box | FINDING |
| 3 | moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:64-75 | boundary-value positivity checks | Only for endpoint failures; an interior negative perturbation can pass | Only partially; does not check full `sigmaXi >= 0` on the box | FINDING |
| 4 | moving_throat_pde_stage126_positive_source_families_sympy_audit.py:80-83 | raise if `2/pi < g_- < pi/4` is false | Yes | Yes; brackets the lower branch between family endpoint biases | PASS |
| 5 | moving_throat_pde_stage126_positive_source_families_sympy_audit.py:13 / moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:26 | Stage 126 banner | N/A label fix | Traceability only | PASS |

## Findings (only if verdict = FINDINGS)
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage126_positive_source_families_sympy_audit.py:55-66`; `moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:64-75`
- **What's wrong:** The added positivity checks only test `xi=0`, `xi=1`, and one corner. For example, SymPy checks `min_xi0 = ... sigma_xi.subs(xi, 0)`, `val_xi1 = ... sigma_xi.subs(xi, 1)`, and `sigma_corner = sigma_xi.subs([(z, L), (xi, 0)])`. A bad family with an added term like `-A*xi*(1-xi)*sin(2*pi*z/L)/L` can preserve those endpoint/corner checks and normalization while becoming negative in the interior of `[0,L] x [0,1]`. The saved outputs corroborate the endpoint checks, but not the full paper claim that the convex family is positive.
- **What it should be:** Assert the global claim directly, e.g. prove/check `sigma_xi(z) >= 0` for `0 <= z <= L`, `0 <= xi <= 1`, `L > 0`, or add a non-tautological structural check that `sigma_xi` is affine in `xi` plus endpoint minima. In Mathematica, a `Resolve[ForAll[...]]`/successful `Minimize`-style check would match the directive’s intent.

## Bottom line
The interval hardening and banner fix are fine, and the saved transcripts show the new checks executing. The problem is F2: the paper’s load-bearing claim is positivity of the whole convex source family, but the edited scripts only verify boundary slices and a corner. That is not enough for an adversarial audit, because a sign-changing interior source can pass the added assertions while violating the stage’s positive-family claim.
