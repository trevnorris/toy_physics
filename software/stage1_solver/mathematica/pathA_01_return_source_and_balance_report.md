# Path-A 01 Return Source and Balance Verification

Genuine checks pass: PASS
Verifier including construction-restatements/guards pass: PASS

## Genuine Checks

| Check | Result | Detail |
| --- | --- | --- |
| reciprocity | PASS | Taylor series gives L_cross=k1 eta drho and L_self=-(1/2) k2 rho0 eta^2; equation-kernel Hessian is symmetric with kernel -k1. |
| gauge_invariance | PASS | Three live identities pass: E_w, C_a, and the hydrodynamic covariant current are gauge invariant. |
| dimensional_consistency | PASS | Matter return terms and all wall-PDE left-hand terms share the same source dimension vector. |
| reduction_consistency | PASS | The local eps^2 coefficient of promoted S_Sigma reproduces mu_eta, T_w, T_Omega, and the unintegrated quadratic terms. |
| target_blind_artifacts | PASS | Note, report markdown, Mathematica verifier, SymPy verifier, and SymPy diagnostics contain no forbidden target-leakage tokens. |

## Construction-Restatements and Guards

These are excluded from the genuine derivation gates and carry `_not_a_physics_gate` suffixes.

| Check | Result | Detail |
| --- | --- | --- |
| no_numerator_knob_corollary_not_a_physics_gate | PASS | Corollary of reciprocity: return kernels are adjoints of the forward kernels, so they add no new free numerator magnitude. The P/N0 symbol-placement scan is only a construction restatement. |
| honest_open_items_disclosure_guard_not_a_physics_gate | PASS | Disclosure-presence guard: the Posited or Open section still names constitutive, gauge, boundary, and no-solve limits. |
| ibp_bookkeeping_identity_not_a_physics_gate | PASS | Bookkeeping identity: symbolic IBP derives dA=partial_w(T_w,Sigma,R R0') and the compact K_eta algebra agrees after dropping the named boundary density. |
| declared_maxwell_no_direct_eta_variation_not_a_physics_gate | PASS | Declared fixed-domain Z(w),H(w) Maxwell action has zero direct eta variation by construction. |
| allowed_gauge_data_placeholder_not_a_physics_gate | PASS | Construction placeholder: the allowed reduced gauge data list contains only J_cov_eta, E_w, and C_a names, not raw potentials or chi. |

Non-physics guards/restatements pass: PASS

Diagnostics JSON: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/pathA_01_return_source_and_balance/pathA_01_return_source_and_balance_diagnostics.json`