FAIL_NOT_CODIM1

## Gate B4 Smectic Verdict

Baseline `B0_Lifshitz` does produce a finite-`k` density softening over an open region, but the generic nonlinear weak-crystallization onset is not a single-`k` codim-1 lamellar stack. A rank-2 equilateral-triad modulation on the same 4D `|k|=k_*` shell has a nonzero cubic invariant and lowers the fixed-mean energy before any single-`k` lamella can be the ground state. The failure is therefore `FAIL_NOT_CODIM1`, conditional on the frozen postulated driver family `B0_Lifshitz`.

This verdict is now produced by both engines from computed branch conditions: finite-`k` softening, convex/no-phase-separation, boundedness, nonnegative minimized light sector, and a fixed-mean lamella-vs-triad Landau minimization. It is not a hardcoded line-item conclusion.

No final layer profile or de Gennes `B,K` extraction is reported: directive F requires a hard stop once C.5 shows a multi-`k` state is preferred. The single-`k` profile below is used only as an excluded variational branch.

## Engine Status

- SymPy: `timeout 600 python3 software/stage1_solver/tools/pathA_25_gateB4_smectic_sympy.py` exited `0` and wrote `software/stage1_solver/reports/pathA_25_gateB4_results.yaml`.
- Mathematica: `timeout 600 math -script software/stage1_solver/tools/pathA_25_gateB4_smectic.wl` exited `0` and wrote `software/stage1_solver/_scratch/pathA_25_gateB4_smectic_mathematica.json`.
- Engine agreement: Mathematica imports the SymPy JSON machine dump and asserts exact agreement on 31 shared expressions, including `A_rho(k)`, `omega_rho^2`, `k_*^2`, `c_L1^crit`, threshold `k_*^2`, triad averages, `Gamma_T`, the `Gamma_T=0` surface, the B4d benchmark minimum, SH `M=1..4` control energies, boundedness, light-sector zeros, and the computed lamella-vs-triad minima.

## Freeze Fidelity

Both scripts recompute the freeze blocks before physics:

- T0 hash: `8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064`
- G0 hash: `f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290`
- The exact T0 `freeze-action` bytes are asserted embedded inside the G0 block.
- Frozen baseline lines for `B0_Lifshitz`, `lambda_Cdiv=0`, `chi_Cpin=0`, `L_Mac`, `L_Pu`, and no prescribed layer profile/normal are asserted present.

B4 background conditions used throughout: `A_mu=0` up to pure gauge, `V_conf=0`, no prescribed support, no fixed layer normal, no boundary potential, fixed conserved mean `bar rho=rho0`, and `rho(X)>=0`.

## Control Predeclarations

NC-B4b positive method control is the 4D, mean-zero, `Z2`-symmetric conserved Swift-Hohenberg/Brazovskii shell:

`F = int d^4X [ (r/2) phi^2 + (1/2)((Delta+k0^2)phi)^2 + (u/4)phi^4 ]`, with shell coefficient `r<0`, `u>0`.

For `M=1..4` independent shell directions, both scripts construct the shell Landau energy from `<phi^2>` and `<phi^4>`, minimize in the amplitude, and obtain:

- `M=1`: `-r^2/(6u)`
- `M=2`: `-r^2/(9u)`
- `M=3`: `-r^2/(10u)`
- `M=4`: `-2r^2/(21u)`

The method therefore returns a stable single-`k` stripe in the declared 4D control window.

NC-B4d multi-`k` benchmark is the 4D shell Landau model with a rank-2 equilateral-triad cubic invariant and predeclared coefficients `tau=1/100`, `gamma=1`, `u=1`. The benchmark has `F_lam,min=0` at positive `tau`, while the independently minimized triad gives `F_tri,min=-3.614683057877558e-05`, so the method returns `FAIL_NOT_CODIM1`.

Other controls:

- NC-B4a: driver off gives `omega^2=(rho0 k^2/m)(5K rho0^3 + hbar^2 k^2/(4m rho0))>0`, hence `FAIL_NO_MODULATION`.
- NC-B4c: a pure attractive negative-compressibility kernel softens first at `k=0`, hence `FAIL_PHASE_SEPARATION`.

## Units

The scripts restore MLT dimensions:

- `[omega^2]=T^-2`
- `[c_L1]=M L^8 T^-2`
- `[c_L2]=M L^10 T^-2`
- density Hessian kernel `[A_rho(k)]=M L^6 T^-2`
- if reached, de Gennes `[B]=M L^-2 T^-2` and `[K_bend]=M T^-2`

The T0 length `a` drops out of the density stability region because the `rho-P` quadratic block vanishes.

## Coupled Linear Spectrum

With `rho=rho0+eta`, fixed mean, and `P0=(1,0,0,0)`, the quadratic density kernel is

`A_rho(k)=5K rho0^3 + (hbar^2/(4m rho0)-c_L1) k^2 + c_L2 k^4`.

The density-phase quadratic Lagrangian has canonical entry `hbar eta dot(theta)` and energy Hessian block

`diag(A_rho(k), rho0 hbar^2 k^2/m)` for `(eta,theta)`.

Thus

`omega_rho^2(k) = (rho0 k^2/m) A_rho(k)`.

The roton minimum is at

`k_*^2 = (c_L1 - hbar^2/(4m rho0))/(2 c_L2) = c_L1/(2c_L2) - hbar^2/(8m rho0 c_L2)`,

when `c_L1 > hbar^2/(4m rho0)`. The finite-`k` softening threshold is

`c_L1^crit = hbar^2/(4m rho0) + 2 sqrt(5K rho0^3 c_L2)`,

with `k_*^2` at threshold equal to `sqrt(5K rho0^3/c_L2)`.

The P-sector block is independent at quadratic order. The scripts compute the stiffnesses, divide by the computed `p_inertia=m rho0 a^2`, and get:

- Three transverse `O(4)->O(3)` Goldstones: `omega_P,T^2 = c_s^2 k^2 = 5K rho0^4 k^2/m`.
- One longitudinal/magnitude mode: `omega_P,L^2 = c_s^2 k^2 + 2 c_s^2/a^2`.
- All `rho-P` and `theta-P` Hessian blocks vanish by explicit mixed-derivative checks. The `rho c_s^2(rho)` dependence multiplies `(partial P)^2`, `(dot P)^2`, or `(P^2-1)^2`, so its first mixed terms are cubic or higher about constant unit `P`.

## Light Sector C.0

About the uniform state, `Sigma_n[rho0]` is empty and the light package has no quadratic contribution. This empty-support statement is labeled in the YAML as `analytic_argument_not_computed`.

On any formed stripe support, the static light energy minimized over the layer fields is

`E_light = sum_n int d^3sigma [ (1/2) mu_br |curl u|^2 + (1/2) kappa_Pu |delta P_parallel - Omega_u|^2 + inherited nonnegative T0 surface terms ] >= 0`

for `mu_br>=0` and `kappa_Pu>=0`. The minimum is attained by `u=0` and `delta P_parallel=Omega_u=0`, giving identically zero integrand on every `Sigma_n[rho]`. Therefore variation of the support `Sigma_n[rho]` multiplies zero and cannot shift or destabilize the density ground state. `J_Pu` is inertial and does not enter the static energy. If G0 allowed negative `mu_br` or `kappa_Pu`, the light package would have its own collapse independent of B4; that is not used to rescue the baseline.

## Finite-k vs k=0

At `k=0`, `A_rho(0)=5K rho0^3>0`; long-wavelength compressibility does not go negative. Same-mean macroscopic phase separation is also excluded for the baseline density sector because the scripts compute `d^2(K rho^5/4)/d rho^2 = 5K rho^3 >= 0` on `rho>=0`, so Jensen's inequality makes the uniform local term minimal at fixed mean.

The energy is bounded below at fixed mean: minimizing the Fourier symbol gives

`min_q [(1/2)c_L2 q^2 - (1/2)c_L1 q] = -c_L1^2/(8c_L2)` at `q=k^2=c_L1/(2c_L2)`,

while the positive `U(rho)~rho^5` controls the resulting `L2` term and the quantum-pressure term is nonnegative.

## Nonlinear Branches and C.5 Morphology

The finite-`k` instability lives on the full `S^3` direction shell by O(4) invariance, so the selected axes are spontaneous, not box-pinned.

For a single-`k` lamella,

`eta_L = A cos(k n) + O(A^2)`,

with zero mean and chemical potential fixing the conserved density. This branch has no cubic invariant:

`<eta_L^3>=0`, `<eta_L |grad eta_L|^2>=0`.

For a rank-2 equilateral triad on the same `S^3` shell, choose `k1+k2+k3=0`, `|ki|=k_*`, and

`eta_T = A [cos(k1.X)+cos(k2.X)+cos(k3.X)] + O(A^2)`.

The checked averages are

- `<eta_T^2>/A^2 = 3/2`
- `<eta_T^3>/A^3 = 3/2`
- `<eta_T |grad eta_T|^2>/(A^3 k_*^2) = 3/4`

The baseline cubic energy coefficient is therefore

`Gamma_T = 15K rho0^2/4 - 3 hbar^2 k_*^2/(32m rho0^2)`.

Equivalently, `Gamma_T=0` only on the codimension-one surface

`k_*^2 = 40 K m rho0^4/hbar^2`.

For generic baseline parameters `Gamma_T != 0`, the triad phase/sign is chosen so `Gamma_T A^3<0`. The lamella has no cubic lowering. Hence an independently minimized rank-2 modulation beats every single-`k` lamella in an open neighborhood of the ordering transition. This is already sufficient to exclude a codim-1 ground state; 3- and 4-direction competitors may lower the energy further, but cannot restore the lamellar verdict.

The computed normalized baseline comparison used for the verdict sets `epsilon=1/100`, positive quartic stabilizer `beta=1`, and a generic off-surface sample with `Gamma_T=891/256`. Independent amplitude minimization gives

- single-`k` lamella: `F_L,min = -1/60000 = -1.6666666666666667e-05`
- rank-2 triad: `F_T,min = -0.09020779478363558`

The triad minimum is strictly lower. At `epsilon=0`, the computed triad minimum is already negative while the lamella minimum is zero, giving the open-neighborhood certificate.

## Stability Region

After C.0/C.1 reductions, the density-region control is

`Lambda = (c_L1 - hbar^2/(4m rho0))/(2 sqrt(5K rho0^3 c_L2))`.

Finite-`k` softening requires `Lambda>1`. A pass-class codim-1 smectic would also require an open region where the single-`k` branch is lower than all multi-`k` states. The baseline fails that requirement generically because the computed triad minimum is below the computed lamella minimum whenever `Gamma_T != 0` in the checked open neighborhood. `Gamma_T=0` remains only the codimension-one tuning above, not an open stability region, and cannot support `SMECTIC_GROUND_STATE_STABLE` or `SMECTIC_CONDITIONAL`.

## Short-circuit

C.1 and C.2 pass the finite-`k` and no-phase-separation checks. C.3 gives bounded density energy and consistent constant `P`, but C.5 fails: a rank-2 multi-`k` state is preferred. Per directive F, the gate stops with `FAIL_NOT_CODIM1`; no stable codim-1 layer profile, `B`, `K`, or rigid-tilt de Gennes extraction is produced.
