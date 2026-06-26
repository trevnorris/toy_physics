DN_UNITTEST_BC_DEPENDENT

## Reduced Operator

Frozen reduction background: `R=R0, eta=0; matter perturbation remains dynamic as exp(-I*omega*t)`.

Emitted reduced operator:

`L_s psi_hat = (cS**2*Derivative(psi_hat(s), (s, 2)) + omega**2*psi_hat(s))/cS**2`

Equivalent ODE artifact:

`psi''(s) + (omega/cS)^2 psi(s) = 0`

Operator check: `operator_is_helmholtz=True`; speed trace: `c_s^2=5*K*rho_star^4/m`; domain check: `[0,L0]`.

## Reduction Certificate

- Background: `straight finite throat, eta=0, s in [0,L0], R0(L0)=0`, `rho0(s)=rho_star`, `A_M0=0`.
- Projection measure: `sqrt(gamma0)(s)=A_perp0 on the straight reference throat`, derivative `0`, so no first-derivative measure term.
- `rho0` fate: `uniform straight reference, O(rho0'/rho0)=0`.
- `c_s(s)` fate: `c_s is the bulk sound speed sqrt(5*K*rho_star^4/m), not wall/healing renormalized`.
- `V_conf(Sigma0)` fate: `absorbed into the stationary reference chemical potential / zero off the wall`; `delta_V_conf=0 in the frozen eta=0 unit test`.
- Quantum potential: background `grad Q=0`; perturbation `delta Q` gives BdG term `hbar**2*omega**4/(4*cS**4*m**2)` with ratio `hbar**2*omega**2/(4*cS**4*m**2)` and is deferred only under `k*xi << 1 with xi=hbar/(m*c_s); equivalently omega << m*c_s^2/hbar`.

## DtN Derivation

- General ODE solution from `dsolve`: `C1*sin(omega*s/cS) + C2*cos(omega*s/cS)`
- D/N coefficient matrix: `[['0', '1'], ['omega*cos(L0*omega/cS)/cS', '-omega*sin(L0*omega/cS)/cS']]`
- D/N determinant: `-omega*cos(L0*omega/cS)/cS`
- BC-applied solution: `psiM*(sin(omega*s/cS)*tan(L0*omega/cS) + cos(omega*s/cS))`
- DtN before final simplification: `-omega*tan(L0*omega/cS)/cS`
- Derived outward-mouth DtN: `-omega*tan(L0*omega/cS)/cS`

## Pole Ladder

Pole equation: `Eq(cos(L0*omega/cS), 0)`.

Solved ladder: `omega_n = pi*cS*(n + 1/2)/L0, n=0,1,2,...` with `halfshift=True`.

## Static Limit

`-L0*omega**2/cS**2 - L0**3*omega**4/(3*cS**4) + O(omega**6)`

This is the small-omega expansion of the dynamic DtN; no separate static solve is used.

## Round Trip

`r_D=-1`, `r_N=+1`, `R_rt=-exp(2*I*L0*omega/cS)` and after substituting the D/N ladder `R_rt=1`, `phi0=0 mod 2*pi`.

## Robin Counterfactual

- Robin coefficient matrix: `[['0', '1'], ['-(alpha*cS*sin(L0*omega/cS) - omega*cos(L0*omega/cS))/cS', '-(alpha*cS*cos(L0*omega/cS) + omega*sin(L0*omega/cS))/cS']]`
- Robin determinant: `(alpha*cS*sin(L0*omega/cS) - omega*cos(L0*omega/cS))/cS`
- Robin denominator core: `-(alpha*cS*sin(L0*omega/cS) - omega*cos(L0*omega/cS))/cS`
- Robin DtN: `omega*(alpha*cS*cos(L0*omega/cS) + omega*sin(L0*omega/cS))/(cS*(alpha*cS*sin(L0*omega/cS) - omega*cos(L0*omega/cS)))`
- Numeric alpha: `2/L0` gives `-omega*(L0*omega*sin(L0*omega/cS) + 2*cS*cos(L0*omega/cS))/(cS*(L0*omega*cos(L0*omega/cS) - 2*cS*sin(L0*omega/cS)))`
- Guard: `{'robin_determinant_emitted': True, 'recovers_DN_at_alpha0': True, 'recovers_DD_at_alpha_inf': True, 'halfshift_destroyed_for_DD': True, 'numeric_alpha_distinct': True, 'dtn_mismatch': True}`

## BC Provenance

`bc_provenance=imposed`, `bc_derivation_emitted=False`.

The explicit `V_wall` mouth/cap gradient derivation is not emitted, so the verdict remains `DN_UNITTEST_BC_DEPENDENT` rather than `DN_UNITTEST_PASS`.

## Engine Agreement

`engine_agreement=True` via Mathematica `FullSimplify[a-b==0]` checks. Details: `{'status': 'pass', 'digest_matches': True, 'mathematica_yaml': 'software/stage1_solver/_scratch/pathA_30_mathematica_results.yaml', 'details': {'dtn': True, 'pole_denominator': True, 'robin_dtn': True, 'robin_denominator': True, 'static_series': True, 'dd_limit': True, 'alpha0_limit': True}, 'mathematica_route': 'transfer_matrix_resolvent'}`.

## Dimensional Check

`dim_check=pass`.
