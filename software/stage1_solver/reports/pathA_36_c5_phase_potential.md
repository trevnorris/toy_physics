# pathA_36 C5 Phase-Potential Derivation

artifact_status: computed dual-engine result  
directive: `software/stage1_solver/directives/pathA_36_c5_phase_potential.md`  
frozen input: `T0_SHEAR_FROZEN(d9520d3819c3)`  
machine result: `software/stage1_solver/reports/pathA_36_c5_phase_potential_results.yaml`

## Verdict

`FAIL_CAUCHY_STRAY_LONGITUDINAL`

The main single-medium branch, with `delta_rho_B = -rho_B0 div u` and finite conjugate-density stiffness, has one physical
longitudinal propagating DOF per finite `k`.  The engines derive a second-class constraint pair, not a Maxwell first-class Gauss
chain.

The physical obstruction is twofold:

- finite density stiffness contributes `B_eff = rho_B0^2/chi_c`, a Cauchy longitudinal modulus;
- the provenance-fixed order-parameter phase sign is conventional, `K_theta = -kappa_phase < 0`, while the Maxwell square requires
  the electric sign `K_theta = C_J^2/rho_br > 0`.

The tuned Maxwell fixture does work algebraically, but only with
`K_theta = J^2 rho_B0^2/rho_br`, `B_eff = 0`, and `m_theta^2 = 0`.  That is `BY_TUNING`, not `WITH_PROVENANCE`.

## Primitive Start

Both engines start from the directive's primitive quadratic Lagrangian:

`L = 1/2 rho_br dot(u)^2 - 1/2 mu_R (curl u)^2 - 1/2 B (div u)^2 + J dot(theta) delta_rho_B + 1/2 K_theta (grad theta)^2`

No pre-completed `1/2 epsilon (...)^2` is used as input.  In the longitudinal finite-`k` convention used by the engines,

`L_L = 1/2 rho_br dot(u_L)^2 - C_J k u_L dot(theta) + 1/2 K_theta k^2 theta^2 - 1/2 B_eff k^2 u_L^2`

with the derived sign

`C_J = -J rho_B0`.

## Dirac-Bergmann Result

For the first-order longitudinal system, the canonical momenta are

`p_u = rho_br dot(u_L)`,  
`pi_theta = J k rho_B0 u_L`.

Thus the primary constraint is

`Phi_1 = pi_theta - J k rho_B0 u_L`.

Constraint preservation gives the secondary constraint

`Phi_2 = -k (J p_u rho_B0 + k kappa_phase rho_br theta)/rho_br`

on the provenance-fixed conventional-sign branch.  The constraint bracket is

`{Phi_1, Phi_2} = k^2 (J^2 rho_B0^2 + kappa_phase rho_br)/rho_br`,

so the constraints are second-class.  With two configuration variables and two second-class constraints,

`N_phys = (4 - 2)/2 = 1`.

For finite stiffness, `B_eff = rho_B0^2/chi_c`, the derived pole is

`omega^2 = k^2 kappa_phase rho_B0^2 / (chi_c (J^2 rho_B0^2 + kappa_phase rho_br))`,

with positive residue and a bounded reduced Hamiltonian.  It is a real stray longitudinal mode, not a ghost and not a gauge-removal.
The independent initial-data count is 2 functions per finite `k`.

If the density stiffness is absent and `B_eff=0`, the same wrong-sign `K_theta` branch has one second-class physical zero mode:
`FAIL_C5_LONGITUDINAL_ZERO_MODE`.

## Tuned Locus

On the algebraic Maxwell locus,

`K_theta = C_J^2/rho_br = J^2 rho_B0^2/rho_br`,  
`B_eff = 0`,  
`m_theta^2 = 0`,

the engines derive `{Phi_1, Phi_2}=0`, two first-class constraints, no longitudinal pole, and 0 physical longitudinal DOF.  The local
generator is

`G[chi] = (rho_br/C_J) (chi Phi_2 - dot(chi) Phi_1)`,

giving `delta u_L = k chi`, `delta theta = -(rho_br/C_J) dot(chi)`.  No inverse `k`, inverse `omega`, or on-shell condition is used.
The longitudinal Hamiltonian is bounded on the Gauss constraint surface.

This is not `WITH_PROVENANCE`: the frozen definitions do not force the electric sign `K_theta>0`, do not force
`K_theta = J^2 rho_B0^2/rho_br`, and do not remove the finite `rho_B0^2/chi_c` Cauchy term.

## Branches

| branch | result | derived count |
|---|---:|---:|
| (b) slaved, finite stiffness, conventional `K_theta` | `FAIL_CAUCHY_STRAY_LONGITUDINAL` | 1 longitudinal DOF, 2 initial-data functions |
| (b) slaved, curl-only/no stiffness, conventional `K_theta` | `FAIL_C5_LONGITUDINAL_ZERO_MODE` | 1 longitudinal zero-mode DOF |
| (b) slaved, tuned Maxwell locus | `C5_RESOLVED_MAXWELL_BY_TUNING` | 0 longitudinal DOF |
| (a) independent `delta_rho_B` with continuity, integrated out | `FAIL_CAUCHY_STRAY_LONGITUDINAL` | same effective slaved finite-stiffness sector |
| (a) continuity removed ablation | `FAIL_EXTRA_SCALAR_DOF` | extra propagating `theta` scalar |

Branch (a) proof: the finite-frequency continuity equation gives
`omega (delta_rho_B + rho_B0 k u_L)=0`, hence `delta_rho_B = -rho_B0 k u_L` in the fixed-number sector.  Substitution gives the same
Josephson cross-term and the same `B_eff` increment `rho_B0^2/chi_c`.  If continuity is removed, Gaussian integration instead gives
`1/2 chi_c J^2 dot(theta)^2`, which is the extra-scalar/second-medium ablation.

## Transverse Sector

The transverse block remains

`L_T = 1/2 rho_br dot(u_T)^2 - 1/2 mu_R k^2 u_T^2`

for each of two transverse polarizations.  Both engines derive 2 massless transverse photons with

`omega^2 = (mu_R/rho_br) k^2`, `c_gamma^2 = mu_R/rho_br`.

The `epsilon != rho_br` control closes the longitudinal square with `epsilon = 2 rho_br`, but fails as
`FAIL_TRANSVERSE_DISTURBED` because the transverse speed shifts to `mu_R/(2 rho_br)`.

## Controls

All directive controls fired:

| control | status |
|---|---|
| 1 no-theta | `FAIL_C5_LONGITUDINAL_ZERO_MODE` |
| 2 Cauchy bulk, no theta | `FAIL_CAUCHY_STRAY_LONGITUDINAL` |
| 3 positive `K_theta` detuned, `B=0` | `FAIL_C5_LONGITUDINAL_ZERO_MODE` |
| 3 `K_theta <= 0` | `FAIL_C5_LONGITUDINAL_ZERO_MODE` |
| 3 positive `K_theta` detuned with `B != 0` | `FAIL_CAUCHY_STRAY_LONGITUDINAL` |
| 3 positive `K_theta` with negative residue | `FAIL_GHOST_OR_NEGATIVE_NORM` |
| 3 `B != 0` on the square locus | `FAIL_SECOND_CLASS_NOT_MAXWELL` |
| 3 `epsilon != rho_br` | `FAIL_TRANSVERSE_DISTURBED` |
| 4 decoupled theta, slaved/no theta kinetic | `FAIL_C5_LONGITUDINAL_ZERO_MODE` |
| 4 decoupled theta, independent density without continuity | `FAIL_EXTRA_SCALAR_DOF` |
| 5 transverse | `PASS_TRANSVERSE_UNDISTURBED` |
| 6 provenance ablation | fixed coefficients fail; free locus gives `C5_RESOLVED_MAXWELL_BY_TUNING` |
| 7 compressibility absent vs included | `FAIL_C5_LONGITUDINAL_ZERO_MODE -> FAIL_CAUCHY_STRAY_LONGITUDINAL` |
| 8 theta mass | `FAIL_SECOND_CLASS_NOT_MAXWELL` |

## Dimensional Firewall

Both engines restored units inline for 18 expressions: brane inertia, MacCullagh curl energy, Cauchy bulk term, Josephson density
term, slaved density, signed phase gradient, density compressibility, theta mass, `C_J`, IBP cross-term, all three electric-square
pieces, the Maxwell locus `C_J^2 = rho_br K_theta`, `c_gamma^2`, the independent-branch theta kinetic, the slaved compressibility
increment, and transverse `omega^2`.

Three able-to-fail ablations fired:

- dropping `rho_B0` from the Josephson cross-term;
- omitting the gradient from the phase stiffness;
- multiplying by `chi_c` instead of dividing the density stiffness by `chi_c`.

## Drift

`DRIFT(6)` for branch (b): new field `theta`; new constants `B`, `J`, `rho_B0`, `K_theta`, `chi_c`.

`DRIFT(7)` for branch (a): adds the independent field `delta_rho_B`.

No `AXIS_RE_ADMITTED` or `U_W_COLLISION` tag fired.

## Dual-Engine Agreement

`ENGINE_AGREE`.

Compared headline fields include the main constraint classification, main DOF count, initial-data count, pole count, curl-only
subcase, tuned-locus first-class count, branch (a) integrated verdict, transverse count/speed, all controls, and dimensional ablation
count.  The SymPy and Mathematica engines each assemble these quantities from their own primitive Lagrangian and Poisson-bracket
calculation.

## Run Log

Commands run from `/var/projects/toy_physics`:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_36_c5_sympy.py
exit 0

timeout 600 math -script software/stage1_solver/tools/pathA_36_c5.wl
exit 0

timeout 600 python3 software/stage1_solver/tools/pathA_36_c5_sympy.py --compare
exit 0
```

Scripts:

- `software/stage1_solver/tools/pathA_36_c5_sympy.py`
- `software/stage1_solver/tools/pathA_36_c5.wl`
