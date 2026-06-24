RC_DENSITY_SMECTIC_LIGHT_NOGO

## Gate R/C Verdict

Conditional on the frozen postulated Family R and Family C density drivers, no branch produces an admissible light-viable codim-1 density smectic. Family R remains in the B4 cubic no-go class; Family Cdiv fails the G0 admission premise through a derived long-wavelength Goldstone response; Family Cpin creates an open density-smectic window for the negative sign, but the computed static C.6 minimizer has zero in-plane polar order.

Per-branch verdicts:
- `S_R_kernel`: `RC_S_R_kernel_FAIL_NOT_CODIM1` at C.5 - R is bilinear and isotropic; the kept GNLS cubic survives at the R roton off a codimension-one surface.
- `S_L_plus_Cdiv`: `RC_S_L_plus_Cdiv_FAIL_ADMISSION` at B.0/C.2 - integrating out transverse P Goldstones shifts the k->0 EOS density stiffness directionally for any nonzero lambda_Cdiv.
- `S_L_plus_Cpin`: `RC_S_L_plus_Cpin_FAIL_LIGHT_STARVED` at C.6 - negative chi_Cpin has an open density-smectic anisotropy window, but the computed static Omega_u=0 minimizer gives P_parallel=0.

## Engine Status

- SymPy: `timeout 600 python3 software/stage1_solver/tools/pathA_25_gateRC_cubic_sympy.py` exited `0` and wrote this report/YAML.
- Mathematica: `timeout 600 math -script software/stage1_solver/tools/pathA_25_gateRC_cubic.wl exited 0 and asserted PASS`.
- Engine agreement: `PASS` on 72 shared symbolic expressions.

## Freeze Fidelity

- T0 hash: `8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064`
- G0 hash: `f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290`
- Exact T0 freeze-action bytes embedded in G0: `true`.
- The corrected Family-C check changes only verification logic outside the hash-bearing freeze-action block.

## Corrected G0 Family-C k->0 Check

The P-sector static block used here is three transverse Goldstones with stiffness `5*K*a**2*rho0**5 k^2` plus one magnitude mode with gap `10*K*rho0**5` and the same gradient stiffness.

`E_Cdiv` gives mixed Hessian entries `-k*lam*sin(theta)` and `-k*lam*cos(theta)`. Integrating out the transverse Goldstone plus magnitude mode gives

`Delta A_Cdiv = -lam**2*(a**2*k**2*sin(theta)**2 + a**2*k**2*cos(theta)**2 + 2*sin(theta)**2)/(5*K*a**2*rho0**5*(a**2*k**2 + 2))`.

The computed low-`k` limit is `(5*K*a*rho0**4 - lam*sin(theta))*(5*K*a*rho0**4 + lam*sin(theta))/(5*K*a**2*rho0**5)`, so the directional shift from the EOS stiffness is `-lam**2*sin(theta)**2/(5*K*a**2*rho0**5)` and the liminf shift is `-lam**2/(5*K*a**2*rho0**5)`. This is an admission failure for any nonzero `lambda_Cdiv`; if `lambda_Cdiv^2 >= 25*K**2*a**2*rho0**8`, it is also a phase-separation failure.

`E_Cpin` has no rho-P quadratic mixed block (`0`) and contributes the direct density Hessian `chi*k**2*cos(theta)**2`. Its low-`k` limit is `5*K*rho0**3`, so it preserves the EOS stiffness and reaches the morphology test. The angle-dependent roton is `k_*^2=(4*c_L1*m*rho0 - 4*chi*m*rho0*cos(theta)**2 - hbar**2)/(8*c_L2*m*rho0)` with threshold `c_L1>(8*sqrt(5)*sqrt(K)*sqrt(c_L2)*m*rho0**(5/2) + 4*chi*m*rho0*cos(theta)**2 + hbar**2)/(4*m*rho0)`; the soft direction is parallel to `P0` for negative `chi_Cpin` and perpendicular for positive `chi_Cpin`.

## Family R

`A_Rho(k) = (4*AR*k**4*m*rho0 - 8*AR*k**2*kR**2*m*rho0 + 20*K*kR**4*m*rho0**4*exp(k**2/kR**2) + hbar**2*k**2*kR**4*exp(k**2/kR**2))*exp(-k**2/kR**2)/(4*kR**4*m*rho0)` and `omega_rho^2=k**2*(4*AR*k**4*m*rho0 - 8*AR*k**2*kR**2*m*rho0 + 20*K*kR**4*m*rho0**4*exp(k**2/kR**2) + hbar**2*k**2*kR**4*exp(k**2/kR**2))*exp(-k**2/kR**2)/(4*kR**4*m**2)`.
The roton is determined by `-x*(4*AR*m*rho0*x**4 - 16*AR*m*rho0*x**2 + 8*AR*m*rho0 - hbar**2*kR**2*exp(x**2))*exp(-x**2)/(2*m*rho0) = 0` with threshold `-(20*K*m*rho0**4 + hbar**2*kR**2*x**2)*exp(x**2)/(4*m*rho0*x**2*(x**2 - 2))` over the negative annulus `0 < x < sqrt(2)`.
The R quadratic term has computed cubic contribution `0`. The cubic invariant is `-3*(-40*K*m*rho0**4 + hbar**2*k**2)/(32*m*rho0**2)`, zero only on `k_Rstar^2 = 40*K*m*rho0**4/hbar**2`. Therefore R is `FAIL_NOT_CODIM1` in every finite-k region off that codimension-one surface.

## Family Cpin Morphology

Predeclared competitors: single-k lamella, rank-2 equilateral triad, rank-3 orthogonal regular star, rank-4 orthogonal regular star. Each regular competitor uses one common amplitude and re-optimizes `P0` over the declared high-symmetry orientations.

For `chi_Cpin < 0`, the computed open density-smectic certificate is `eps>0, beta>0, gamma>0, zeta > 2*(135*beta*eps + 4*gamma**2)/(135*beta)`. In that region the lamella minimum is `-eps**2/(6*beta)` and the declared multi-k competitors have nonnegative minima, so C.5 passes as a density-only statement.
Omitted BCC/FCC-like cubic-active competitors are flagged as `unresolved_by_scope_but_large_zeta_gap_covered`. The computed gap lower bound is `3/2`, giving finite large-`zeta` cover threshold `(4*avg2O*eps*u4O + c3O**2*gamma**2)/(3*u4O)` for any finite omitted competitor of the recorded generic form.
The synthetic anisotropic-shell positive control uses `{'eps': '1/100', 'beta': '1', 'gamma': '1', 'zeta': '3', 'lamella_min': '-1/60000', 'rank2_triad_min': '0', 'rank3_min': '0', 'rank4_min': '0', 'NC_RC_Cescape_pass': True}` and returns `NC_RC_Cescape_pass=true`.

For `chi_Cpin > 0`, the triad can lie in `P0`-perp and keeps the baseline cubic energy `A**2*(45*A**2*beta - 8*A*gamma - 12*eps)/8`; the sample triad minimum is `-(10 + sqrt(130))**2*(2*sqrt(130) + 29)/27000000` below the lamella sample `-1/60000`, so that sign is `FAIL_NOT_CODIM1`.

## C.6 Light-Substrate Test

The quotient energy is `(curlU**2*muBr + kappaPu*mismatch**2)/2` with Hessian `[['muBr', '0'], ['0', 'kappaPu']]` and positive-definite conditions `['mu_br > 0', 'kappa_Pu > 0']`.
The orientation energy minimized for `chi_Cpin<0` is `(KPSigma*pParallel**2*qParallel**2 + OmegaU**2*kappaPu + 4*OmegaU**2*muBr - 2*OmegaU*kappaPu*pParallel + chiAbs*gRhoSqSigma*pParallel**2 + kappaPu*pParallel**2)/2` with stationarity `KPSigma*pParallel*qParallel**2 - OmegaU*kappaPu + chiAbs*gRhoSqSigma*pParallel + kappaPu*pParallel=0`, so `P_parallel(Omega_u)=OmegaU*kappaPu/(KPSigma*qParallel**2 + chiAbs*gRhoSqSigma + kappaPu)`.
The static brane-formation ground state has no imposed rotation drive (`Omega_u=0`), hence `P_parallel=0` and `P_parallel^2=0`. C.6 returns `FAIL_LIGHT_STARVED`.
The driven control sets nonzero `Omega_u` and obtains `P_parallel=1/2`, making the `RC_S_L_plus_Cpin_CODIM1_CONDITIONAL` fork reachable without changing the static verdict.

## Controls

- `NC-RC-regress`: executed, pass=`True`; `F_tri_neg(zeta=0)` restores the baseline triad energy.
- `NC-RC-pos`: executed, pass=`True`; the Z2 shell minimizer returns a single-k lamella when the cubic is absent.
- `NC-RC-Cnull`: executed inside `NC-RC-regress`, pass=`True`.
- `NC-RC-Rtriad`: executed, pass=`True`; `<eta_T^3>/A^3=3/2` and `<eta_T |grad eta_T|^2>/(A^3 k^2)=3/4`.
- `NC-RC-Cescape`: executed, pass=`True`; the predeclared anisotropic benchmark returns a single-k lamella.
- `NC-RC-C6-fork`: executed; static pass=`True`, driven conditional pass=`True`.

## Units

MLT checks include `[A_rho]=(1, 6, -2)`, `[Gamma]=(1, 10, -2)`, `[lambda_Cdiv]=(1, 3, -2)`, `[chi_Cpin]=(1, 8, -2)`, `[kappa_Pu]=(1, -1, -2)`, `[mu_br]=(1, -1, -2)`, `[chi_Cpin <|grad rho|^2>_Sigma]=(1, -1, -2)`, and `[omega^2]=(0, 0, -2)`.

## Short-Circuit Discipline

`S_L_plus_Cdiv` stops at admission/C.2. `S_R_kernel` stops at C.5. `S_L_plus_Cpin` reaches C.6 only in the negative-sign density-smectic region and stops at `RC_S_L_plus_Cpin_FAIL_LIGHT_STARVED`; no de Gennes layer constants or downstream gates are reported.
