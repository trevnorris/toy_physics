---
stage: 118
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [stage_118.md, stage_118.tex, moving_throat_pde_stage118_parent_core_sympy_audit.py, moving_throat_pde_stage118_parent_core_sympy_audit.txt, moving_throat_pde_stage118_parent_core_mathematica_audit.wl, moving_throat_pde_stage118_parent_core_mathematica_audit.txt]
---
# Codex review — stage 118

## What the edit was supposed to do
The directive identified a sign mismatch for the mixed coupling and, under resolution direction (a), required the uniform-core lambda convention to be \(\lambda=-q_*v_{w0}I_{sq}\). It also required added coverage for three paper-side deliverables: \(K_q\), \(g_s\), and the general bilinear lambda form \(\lambda=-q_*v_{w0}J_sI_q\). The stage card says the audit must support overlap formulas for \(K_s,K_q,\lambda,g_s,g_q\) from the parent-action throat-core ansatz.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage118_parent_core_sympy_audit.py:97 | `K_q closed form` | No; `K_q` is defined as the same closed form at line 82 | No; it does not connect \(K_q\) to the D/N mode stiffness integral | FINDING |
| 2 | moving_throat_pde_stage118_parent_core_sympy_audit.py:100 | `g_s closed form` | Yes; it depends on the shell moment/J_s reduction | Yes | PASS |
| 3 | moving_throat_pde_stage118_parent_core_sympy_audit.py:102 | `lambda uniform closure` | Yes; it checks the minus sign and closed-form overlap product | Yes | PASS |
| 4 | moving_throat_pde_stage118_parent_core_sympy_audit.py:104 | `lambda from bilinear` | No; it follows from the immediately preceding definitions | No; it does not use the extracted `sq_coeff` | FINDING |
| 5 | moving_throat_pde_stage118_parent_core_mathematica_audit.wl:105 | `K_q closed form` | No; `kQ` is defined as the same closed form at line 90 | No; it does not connect \(K_q\) to `chiGrad` | FINDING |
| 6 | moving_throat_pde_stage118_parent_core_mathematica_audit.wl:108 | `g_s closed form` | Yes; it depends on the shell moment/jS reduction | Yes | PASS |
| 7 | moving_throat_pde_stage118_parent_core_mathematica_audit.wl:110 | `lambda uniform closure` | Yes; it checks the minus sign and closed-form overlap product | Yes | PASS |
| 8 | moving_throat_pde_stage118_parent_core_mathematica_audit.wl:111 | `lambda from bilinear` | No; it follows from the immediately preceding definitions | No; it does not use the extracted `sqCoeff` | FINDING |

## Findings
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage118_parent_core_sympy_audit.py:97`; `moving_throat_pde_stage118_parent_core_mathematica_audit.wl:105`
- **What's wrong:** `expect_zero("K_q closed form", K_q - (Zq/mu0) * sp.pi**2 * c_s**2 / (4*L_W**2))` and `expectZero["K_q closed form", kQ - (zQ/mu0)*Pi^2*cSound^2/(4*lW^2)];` compare \(K_q\) to the exact expression used to define it. I tried to break the check by making the D/N stiffness-to-\(K_q\) connection wrong; the check would still pass because it never uses `chi_grad` / `chiGrad`.
- **What it should be:** Derive or check \(K_q\) through the independently computed stiffness moment, e.g. `K_q - (Zq/mu0)*c_s**2*chi_grad` and the Mathematica analogue.

### R2 — tautological_check
- **Where:** `moving_throat_pde_stage118_parent_core_sympy_audit.py:104`; `moving_throat_pde_stage118_parent_core_mathematica_audit.wl:111`
- **What's wrong:** `expect_zero("lambda from bilinear", lam_uniform + qstar * v0 * J_s * I_q)` and `expectZero["lambda from bilinear", lamUniform + qStar*v0*jS*iQ];` are true by construction because `I_sq_uniform = J_s * I_q` and `lam_uniform = -qstar * v0 * I_sq_uniform`. I tried to break the actual bilinear extraction by ignoring or changing `sq_coeff`; this assertion would not notice.
- **What it should be:** Connect lambda to the extracted coefficient, avoiding division assumptions, e.g. `expect_zero("lambda from bilinear", lam_uniform*varrho_s*A_q - sq_coeff*J_s*I_q)` and the Mathematica analogue.

### R3 — transliteration
- **Where:** `moving_throat_pde_stage118_parent_core_mathematica_audit.wl:90-111`
- **What's wrong:** The Mathematica audit is a line-by-line mirror of the SymPy definitions and targets for the edited section. It independently confirms that Mathematica simplifies the same self-referential residuals to zero, but it is not an independent derivation that could catch copied wrong targets.
- **What it should be:** The second engine should derive the same deliverables through independent intermediates, especially `chiGrad` for \(K_q\) and `sqCoeff` for lambda.

## Bottom line
The saved outputs do show `0` / `PASS` for the edited checks, but two of the new checks are self-fulfilling. The most important defect is that the added `lambda from bilinear` check does not actually test the bilinear coefficient extracted in section IV; it only restates the lambda definition from section V. That leaves the paper-side bilinear deliverable insufficiently verified despite passing transcripts.
