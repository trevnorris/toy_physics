---
stage: 108
reviewer: codex
verdict: FINDINGS
findings_count: 2
files_reviewed: [moving_throat_pde_stage108_robustness_classes_sympy_audit.py, moving_throat_pde_stage108_robustness_classes_sympy_audit.txt, moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl, moving_throat_pde_stage108_robustness_classes_mathematica_audit.txt, stage_108.md, stage_108.tex]
---
# Codex review — stage 108

## What the edit was supposed to do
The directive identified that the scripts only checked the beta=1 additive preservation locus, while the notes claimed the general beta-dependent locus `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27`; direction (a) was to extend both scripts with a scale+argument+additive class and re-solve the even-match compensation. It also required fixing a Mathematica precedence bug so `chi_arg(beta=1) - 1` substitutes beta=1 before subtracting 1. Finally, both script banners had to be updated from stage 91/091 to stage 108.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage108_robustness_classes_sympy_audit.py:25 | SymPy banner says stage 108 | yes | metadata, not math | PASS |
| 2 | moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:26 | Mathematica banner says stage 108 | yes | metadata, not math | PASS |
| 3 | moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:55 | parenthesized `chi_arg(beta=1) - 1` | yes | pure argument beta=1 odd normalization | PASS |
| 4 | moving_throat_pde_stage108_robustness_classes_sympy_audit.py:84-90 | unique beta-parameterized even-match solution | yes | compensated branch even coefficients | PASS |
| 5 | moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:95-96 | unique beta-parameterized even-match solution | yes | compensated branch even coefficients | PASS, but see R2 |
| 6 | moving_throat_pde_stage108_robustness_classes_sympy_audit.py:98-101 | solved locus equals `S(1-beta^5)/9 - Sigma0/27` | yes | general preservation locus | PASS |
| 7 | moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:104-107 | solved locus equals `S(1-beta^5)/9 - sigma0/27` | yes | general preservation locus | PASS, but see R2 |
| 8 | moving_throat_pde_stage108_robustness_classes_sympy_audit.py:102 | substitute solved locus back into same `chi_gen == 1` | no | only self-consistency of `Solve` | FINDING R1 |
| 9 | moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:108 | substitute solved locus back into same `chiGen == 1` | no | only self-consistency of `Solve` | FINDING R1 |
| 10 | moving_throat_pde_stage108_robustness_classes_sympy_audit.py:103-106 | general locus reduces to beta=1 Class C | yes | special-case reduction | PASS |
| 11 | moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:109-112 | general locus reduces to beta=1 Class C | yes | special-case reduction | PASS, but see R2 |

## Findings
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage108_robustness_classes_sympy_audit.py:102`; `moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:108`
- **What's wrong:** `expect_zero('general preservation locus check', chi_gen.subs(Sigma5, chi_pres_gen) - 1)` and `expectZero["general preservation locus check", (chiGen /. sigma5 -> sigma5PresGen) - 1];` substitute a value obtained by solving that same equation. If the model for `chi_gen`/`chiGen` were wrong but internally self-consistent, this check would still pass.
- **What it should be:** Substitute the hard-coded paper locus directly, and assert the iff structure, e.g. check `chi_gen.subs(Sigma5, S*(1-beta**5)/9 - Sigma0/27) - 1 == 0` and factor `chi_gen - 1` to show it is proportional to `Sigma5 - (S*(1-beta**5)/9 - Sigma0/27)` under the nonzero denominator branch.

### R2 — transliteration
- **Where:** `moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:83-112`
- **What's wrong:** The added Mathematica Class D block is a line-by-line translation of the SymPy block: build `lambdaGen`, extract `l0g/l2g/l4g/l5g`, define `m2g/m4g`, solve for even compensation, define `chiGen`, solve for `sigma5`, and compare to the same hard-coded target. It corroborates that a second CAS simplifies the same copied algebra, but it is not an independent derivation; a copied modeling error would survive in both engines.
- **What it should be:** Use Mathematica to check the preservation statement independently, preferably with `Reduce`/`FullSimplify` on `chiGen == 1 <=> sigma5 == sNorm*(1-beta^5)/9 - sigma0/27` under the declared assumptions, or derive `chiGen` from a direct series coefficient route that is not just the SymPy intermediate chain rewritten.

## Bottom line
The corrected Mathematica precedence check is real and the saved output shows `chi_arg(beta=1) - 1 = 0` with `PASS`; the banner fixes are also corroborated in both transcripts. The hard-coded general locus comparison would fail under the old beta=1-only model, a wrong beta power, or a wrong odd normalization, so that assertion is meaningful. The failure is that the subsequent “general preservation locus check” is solved-locus self-substitution, and the Mathematica copy does not provide independent protection against a shared modeling mistake.
