---
stage: 144
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py, moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.txt, moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl, moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.txt, stage_144.md, stage_144.tex]
---
# Codex review — stage 144

## What the edit was supposed to do
F1 was supposed to correct copy-paste banner labels from stage 127 to stage 144 in both engines and transcripts. F2 was supposed to add `Sigma0` reporting plus numerical target assertions for the compensated branch inequalities, `Pi_*`, `That(Pi_*)`, `Sigma0(Pi_*)`, `Pi_match`, and `That(Pi_match)`. F3 and F4 were explicitly unresolved in the directive; F4 already identifies the Mathematica audit as a line-by-line SymPy transliteration.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:16,70 | Stage banner relabelled to 144 | Yes, transcript would still show 127 | Directive label fix, not mathematical | PASS |
| 2 | moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:40-43 | Report `Sigma0(Pi_*)` and `Sigma0(Pi_match)` from the current closure | Yes as output, but print-only | Supports numerical ledger reporting | PASS |
| 3 | moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:53-56 | Assert `g_+^F1 > 1` and `2/pi < g_-^F1 < 1` | Yes | Exercises branch exclusion/bracketing | PASS |
| 4 | moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py:58-67 | Assert numerical targets for `Pi_*`, `That_*`, `Sigma0_*`, `Pi_match`, `That_match` | Yes | Exercises the paper card's numerical fixed-point record | PASS |
| 5 | moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:27,73 | Stage banner relabelled to 144 | Yes, transcript would still show 127 | Directive label fix, not mathematical | PASS |
| 6 | moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:59-62 | Report `Sigma0(Pi_*)` and `Sigma0(Pi_match)` | Yes as output, but print-only | Supports numerical ledger reporting | PASS |
| 7 | moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:64-71 | Assert same branch and target checks as SymPy | Locally yes, but not independently | Same numerical claim, copied derivation | FINDING |

## Findings
### R1 — transliteration
- **Where:** `moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl:64-71`
- **What's wrong:** `If[!(Abs[N[piStar - 1.50882951349316\`30, 30]] < tol), fail["Pi_* drift", N[piStar, 30]], pass["Pi_* matches notes target"]];` and the surrounding added assertions are the Mathematica analogue of the SymPy block, using the same definitions, same root equations, same seeds, and same target constants. The saved Mathematica output shows the added `PASS:` lines, but it is not an independent second-engine derivation; a shared wrong algebraic target copied into both engines would pass in both. This is exactly the unresolved F4 risk, and the directive contains no applied/rework note accepting that policy risk.
- **What it should be:** Either record an explicit accepted resolution that this stage permits transliteration, or replace one engine with an independent derivation/check of the branch law or closure before counting it as second-engine corroboration.

## Bottom line
The SymPy-side F2 checks are non-tautological and the saved output corroborates the new numerical assertions. The banners are also corrected in both transcripts. The blocking issue is that the Mathematica additions are still a line-by-line mirror of the SymPy computation, so the second engine does not provide independent evidence for the paper-side numerical branch claim.
