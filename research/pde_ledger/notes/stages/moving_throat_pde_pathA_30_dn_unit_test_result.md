# pathA_30 frozen-wall D/N unit test result

Verdict: `DN_UNITTEST_BC_DEPENDENT`.

This records the Phase-1 scaffold §9/§11.1 check for the straight frozen throat. The reduced operator is the phonon-limit scalar Helmholtz channel on `s in [0,L0]`, with outward-mouth DtN `-omega*tan(L0*omega/cS)/cS` and pole ladder `omega_n = pi*cS*(n + 1/2)/L0, n=0,1,2,...`. The Robin counterfactual is live and destroys the half-shift in the Dirichlet-cap limit.

BC provenance remains `imposed` because the explicit `V_wall` mouth/cap gradient derivation was not emitted; therefore this is the honest `BC_DEPENDENT` rung, not a full `PASS`.

Dimensional retrofit: the SymPy and Mathematica engines walk the sourced EOS speed law `5*K*rho_star^4/m`, the DtN tangent argument `omega*L0/c_s`, and the `Z00` prefactor. The sourced `[K]` is fixed from the 4D EOS chain `[P]=(-2,1,-2), [rho]=(-4,0,0), [K]=[P]-5[rho]=(18,1,-2)`, not back-solved from `c_s`; corrupting `[K]` flips the recomputed verdict to `DN_UNITTEST_FAIL_DIMENSIONAL`, while the unmutated run remains `DN_UNITTEST_BC_DEPENDENT`.

Primary artifacts:
- `software/stage1_solver/reports/pathA_30_dn_unit_test.md`
- `software/stage1_solver/reports/pathA_30_results.yaml`
- `software/stage1_solver/tools/pathA_30_dn_unit_test_sympy.py`
- `software/stage1_solver/tools/pathA_30_dn_unit_test.wl`
