---
stage: 143
reviewer: codex
verdict: FINDINGS
findings_count: 2
files_reviewed: [moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py, moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.txt, moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl, moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.txt, stage_143.md, stage_143.tex]
---
# Codex review — stage 143

## What the edit was supposed to do
The directive required the scripts to gate, not merely print, the endpoint limits, limiting constants, traction asymptotic, and positivity pieces used to prove \(\mathfrak g_\Pi<1\). It also required Mathematica to replace a hardcoded \(S_\infty=1\) and algebraic substitutions with actual limits of the dynamical \(s_Q,r_Q,\sigma_0,\widehat T\). Finally, it required an independent Mathematica `Reduce` proof of numerator positivity to break the prior SymPy/Mathematica transliteration.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:50 | `pi**2 - 2*pi > 0` | yes | yes, decomposition lemma | PASS |
| 2 | moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:51 | `pi**2/2 - 4 > 0` | yes | yes, decomposition lemma | PASS |
| 3 | moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:53-54 | exponential remainder cubic coefficient | not for global positivity | no, only local Taylor data | FINDING R1 |
| 4 | moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:61-62 | endpoint limits \(g_0, g_\infty\) | yes | yes | PASS |
| 5 | moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:77-79 | \(R_\infty,S_\infty,\widehat T/\sqrt\Pi\) constants | yes | yes | PASS |
| 6 | moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:60-65 | `Reduce[num > 0]` | yes | yes, independent \(g_\Pi<1\) proof | PASS |
| 7 | moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:71-72 | constant positivity pieces | yes | yes, decomposition lemma | PASS |
| 8 | moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:73-74 | exponential remainder cubic coefficient | not for global positivity | no, only local Taylor data | FINDING R2 |
| 9 | moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:82-83 | endpoint limits \(g_0, g_\infty\) | yes | yes | PASS |
| 10 | moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:91-102 | dynamic limits and limiting constants | yes | yes | PASS |

## Findings
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:53-54`
- **What's wrong:** `expect_equal("exp remainder leading term is Pi**3/6", exp_rem_series.coeff(Pi, 3), sp.Rational(1, 6))` checks only one Taylor coefficient. I tried replacing the intended positive remainder mentally with an expression like `Pi**3/6 - Pi**4`: this check still passes, but the expression is negative for large positive `Pi`. It does not assert the directive's required \(e^\Pi-1-\Pi-\Pi^2/2>0\) for all finite \(\Pi>0\).
- **What it should be:** A direct proof/check of the full inequality for `Pi > 0`, for example via a positivity assertion on the whole remainder or a derivative/Taylor-remainder argument that gates global positivity, not just the cubic coefficient.

### R2 — insufficient_verification
- **Where:** `moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:73-74`
- **What's wrong:** `expectEqual["exp remainder leading term is piM^3/6", Coefficient[expRemSeries, piM, 3], 1/6];` has the same defect. It corroborates the saved output's `PASS`, but only for the cubic coefficient, not for \(e^{\pi_M}-1-\pi_M-\pi_M^2/2>0\) on `piM > 0`.
- **What it should be:** A direct `expectPositive`/`Reduce`-style proof of `Exp[piM] - 1 - piM - piM^2/2 > 0` under `piM > 0`.

## Bottom line
The Mathematica hardcoded-limit rework is materially better, and the new `Reduce[num > 0, piM, Reals]` check independently supports the paper's core \(g_\Pi<1\) claim; the saved transcripts show those checks passing. The remaining blocker is that both engines still claim to gate the exponential-remainder positivity while only checking its leading Taylor coefficient. That is not a global positivity proof and can pass for a wrong finite-\(\Pi\) remainder.
