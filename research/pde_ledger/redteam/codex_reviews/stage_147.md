---
stage: 147
reviewer: codex
verdict: FINDINGS
findings_count: 5
files_reviewed: [stage_147.md, stage_147.tex, moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py, moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.txt, moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl, moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.txt]
---
# Codex review — stage 147

## What the edit was supposed to do
The directive was to replace print-only or tautological stage-147 checks with assertions: numerical anchors for \(A_T\), \(B_T\), and \(|A_T|/B_T\), a chain-rule consistency check for \(A_T\), and a centered-kernel structure check. For Mathematica it also required replacing the old `R_q(g_minus)-1/4` tautology and fixing the stage banner. F3 was supposed to add centering / moment support for the rigidity-kernel claim, but the inserted code only resubstitutes cached \(g_*\) and \(S_*\) formulas.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:56 | `A_T` vs paper literal | Yes, against hard-coded paper value | Yes, coefficient quote | PASS |
| 2 | moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:59 | `B_T` vs paper literal | Yes, against hard-coded paper value | Yes, coefficient quote | PASS |
| 3 | moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:62 | `|A_T|/B_T` vs literal | Yes, against hard-coded ratio | Adjacent/redundant to coefficients | PASS |
| 4 | moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:73 | `A_T` closed form vs chain route | Only if the copied algebra diverges | Weak; self-consistency, not independent derivation | FINDING |
| 5 | moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:97 | `Wcenter` constant offset | No meaningful physics failure; built from same `Wcenter` definition | Weak; does not test centering source moments | FINDING |
| 6 | moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:109 | `g_*` resubstitution drift | No; repeats line 25 expression | No moment-orthogonality check | FINDING |
| 7 | moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:111 | `S_*` resubstitution drift | No; repeats line 26 expression | No moment-orthogonality check | FINDING |
| 8 | moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:67 | `A_T` vs paper literal | Yes, outside tolerance | Yes, coefficient quote | PASS |
| 9 | moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:69 | `B_T` vs paper literal | Yes, outside tolerance | Yes, coefficient quote | PASS |
| 10 | moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:71 | `|A_T|/B_T` vs literal | Yes, outside tolerance | Adjacent/redundant to coefficients | PASS |
| 11 | moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:78 | `A_T` closed form vs chain route | Only if the copied algebra diverges | Weak; self-consistency, not independent derivation | FINDING |
| 12 | moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:90 | `Wcenter` centered form | No meaningful physics failure; also samples one `x` | Weak; does not test centering source moments | FINDING |
| 13 | moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:96 | `g_*` resubstitution drift | No; repeats line 41 expression | No moment-orthogonality check | FINDING |
| 14 | moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:98 | `S_*` resubstitution drift | No; repeats line 42 expression | No moment-orthogonality check | FINDING |

## Findings (only if verdict = FINDINGS)
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:69`
- **What's wrong:** `dSigma_dPi_at_star = 1/(1 - S_star/4) + Pi_star * Sp_star / (4*(1-S_star/4)**2)` is then used to rebuild the same algebra used in `AT` at lines 33-38. I tried to break it by assuming the derivative formula itself has a missing factor or sign; the check would pass if both the closed form and this “chain” route share that same copied mistake.
- **What it should be:** Differentiate an independently defined `T(Pi) = sqrt(9*(Pi/(1-Sformula(Pi)/4))/20)` with respect to `Pi`, then compare `-dT/dPi / dgPi/dPi` at `Pi_*` to `A_T`.

### R2 — tautological_check
- **Where:** `moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:84`
- **What's wrong:** `Wcenter = sp.simplify(AT*(c-gminus) + BT*(Kq-Sformula.subs(Pi, Pi_star)))` defines the kernel in exactly the form later asserted at lines 93-98. I tried the failure mode “the centered kernel derivation is wrong but the final expression was typed in the desired form”; this assertion still passes.
- **What it should be:** Derive the kernel from an independent variation/projection formula, or directly check the stated centering moments for `c-g_*` and `K_q-S_*` under the declared source-moment inner product.

### R3 — insufficient_verification
- **Where:** `moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:88`
- **What's wrong:** `wCenterConst = FullSimplify[(wCenter - (aT*c + bT*kq)) /. x -> 1/2]` checks one sample point after `wCenter` was already defined as `aT*(c - gMinus) + bT*(kq - sStar)`. A nonconstant erroneous residue vanishing at `x = 1/2` would pass.
- **What it should be:** Check the full symbolic residual in `x`, and source `wCenter` from an independent derivation rather than the asserted target form.

### R4 — tautological_check
- **Where:** `moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:107`
- **What's wrong:** `g_star_resub = sp.N(gPi.subs(Pi, Pi_star), 40)` and `S_star_resub = sp.N(Sformula.subs(Pi, Pi_star), 40)` repeat the exact definitions of `g_star` and `S_star` from lines 25-26. This can catch mutation, not mathematical centering or moment orthogonality.
- **What it should be:** Compute the moments from the actual inner-product/integral definitions and compare those independent moment values to `gFormula(Pi_*)` and `Sformula(Pi_*)`.

### R5 — transliteration
- **Where:** `moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:63`
- **What's wrong:** The Mathematica assertion block is a direct port of the SymPy block: same literals, same hand-written chain rule, same target-form kernel definition, same resubstitution drift checks. The saved `.txt` files corroborate that both engines print `PASS`, but a wrong target copied into both engines would still pass in both.
- **What it should be:** Use the second engine for an independent derivation path, especially for the chain rule and centering/moment checks.

## Bottom line
The numerical coefficient anchors are real and the saved outputs show the expected `PASS` lines, including the fixed Mathematica banner. The important failure is that the new kernel and moment checks do not independently test the stage’s centered first-order rigidity-kernel claim: they mostly reassert definitions already typed into the scripts. The stage therefore remains insufficiently verified despite passing transcripts.
