# pathA_35 Gate L: Light on the Shear-Surface Brane

artifact_status: computed gate result
directive: `software/stage1_solver/directives/pathA_35_shear_surface_brane_gates.md`
stage: `Gate L`
date: 2026-06-26
frozen input: `T0_SHEAR_FROZEN(d9520d3819c3)` from `pathA_35_G0_freeze.md`

## Verdict

`FAIL_COUPLE_STRESS_NOGO`

Provenance: `CONDITIONAL_ON(both)` because the tested structure is conditional on the imposed axis and the postulated MacCullagh package.

The remediation changed the derivation provenance and corrected the slaved-rigid coefficient, but it did not change the overall verdict.

No `FREE_LIGHT_OK_CONDITIONAL` pass was reached on any frozen configuration. The slaved-rigid branch resolves the live-`P` mode-count plus
closure pair, but the frozen package still has no C5 `phi` analog; a `phi=u_w` rescue collides with the required `u_w` gap, while an
independent `phi` is a fresh G0, not a Gate-L pass.

The verdict is now aggregated from the computed sub-hurdle statuses. An all-sub-hurdles-pass fixture is included in the dual-engine payload
and aggregates to `FREE_LIGHT_OK_CONDITIONAL`.

Machine-readable result: `software/stage1_solver/reports/pathA_35_gateL_light_results.yaml`.

## Computational Setup

Both engines use the flat-brane Fourier basis with `k` along one in-plane axis and build the quadratic action from the frozen `S_G0` terms:
`L_Mac`, `L_Pu`, the brane kinetic term, `L_uw`, and the T0 `P` kinetic/Frank terms projected through `g_ell(w)`. The surface-projected
inherited T0 arrow coefficients are derived as

- `J_P = m rho a^2 ell_g`
- `Gamma_P = m rho c_s^2 a^2 ell_g`

with `int dw g_ell(w)=1` and mode weight `int dw ell_g g_ell(w)=ell_g`. These are not new constants; they are the T0 inertia/stiffness
integrated across the frozen one-width confinement scale. They were dimension-checked inline. The MacCullagh feed-forward speed recorded
for later gates is

`c_gamma^2 = mu_R/rho_br`.

The transverse L(b) Hessian is now taken as the second variation of the derived quadratic energy, not typed into the audit:

`[[mu_R k^2, lambda_Pu k], [lambda_Pu k, Gamma_P k^2]]`.

For the slaved-rigid branch, the script performs the frozen substitution
`P_parallel = w_hat x Omega_u`, hence `varpi = w_hat x P_parallel = -Omega_u`. Retaining the bilinear coupling gives the derived

`mu_eff = mu_R - 2 lambda_Pu`

and transverse dispersion

`omega^2 = ((mu_R - 2 lambda_Pu) k^2 + Gamma_P k^4)/rho_br`.

Thus exact nondispersive equality is false at finite `k`; the branch is clean only as a low-`k` EFT mode with positive effective modulus and

`k^2 << (mu_R - 2 lambda_Pu)/Gamma_P`.

## Configuration A: G0 Baseline

Frozen fields: massless `P^i`, parity-repaired `P-u` coupling, no `phi`, gapped `u_w`.

| sub-hurdle | result | computed driver |
|---|---|---|
| L(a-i) traction-not-torque | PASS, `ARROWS_SUPPLY_TRACTION` | `P-u` virtual work gives rank-2 cut traction; Frank-only rank is 0 |
| L(a-ii) hidden modes | FAIL: `FAIL_HIDDEN_PROPAGATING_MODE` | determinant count: 5 gapless positive branches at `lambda_Pu=0`, i.e. 3 extra `P` spin waves; nonzero `lambda_Pu` gives 2 low-`k` negative `omega^2` branches |
| L(a-iii) C5 | FAIL: `FAIL_C5_LONGITUDINAL_ZERO_MODE` | no-`phi` branch leaves the kinetic longitudinal zero mode |
| L(b) bounded/closure | FAIL: `FAIL_NOT_BOUNDED_BELOW` | transverse `u/P` Hamiltonian minor is `k^2 (Gamma_P mu_R k^2 - lambda_Pu^2)` |
| L(c) leak | PASS | flat direct traction is 0; flat one-`k` indirect vorticity source is 0 |
| L(d) `u_w` gap | PASS | coupled flat scalar spectrum has `u_w` descendant `Omega_w^2` |

For L(b), the Hamiltonian block for one transverse pair is

`[[mu_R k^2, lambda_Pu k], [lambda_Pu k, Gamma_P k^2]]`.

For `0 < k^2 < lambda_Pu^2/(Gamma_P mu_R)`, the principal minor is negative. This is an energy-matrix failure, not a dispersion-only
claim. The `lambda_Pu -> 0` sanity check restores the minor to `Gamma_P mu_R k^4`. The gapped-`P` leg also fires from the derived response
`P = -lambda_Pu Omega_u/(M_gap_P + Gamma_P k^2)`: the low-`k` couple-stress divergence goes to 0 and the closure residual remains
`2 mu_R Omega_u`.

Configuration A verdict: `FAIL_COUPLE_STRESS_NOGO`.

## Configuration B: Slaved-Rigid Escape

Frozen branch: `P_parallel = w_hat x (nabla_parallel x u)`, no independent `P` spin modes; all other G0 choices inherited.

| sub-hurdle | result | computed driver |
|---|---|---|
| L(a-i) traction-not-torque | PASS, `POSTULATED_SURFACE_ELASTICITY` | independent arrow modes are gone; clean `k^2` traction is the postulated surface modulus |
| L(a-ii) hidden modes | PASS with low-`k` dispersion caveat | no extra propagating `P` modes; `k^4` correction remains |
| L(a-iii) C5 | FAIL: `FAIL_C5_LONGITUDINAL_ZERO_MODE` | frozen no-`phi` decision is inherited |
| L(b) bounded/closure | PASS conditional on `mu_R - 2 lambda_Pu > 0` | slaving gives closure residual 0 and energy coefficient `k^2(mu_R - 2 lambda_Pu + Gamma_P k^2)` |
| L(c) leak | PASS | same flat direct/indirect leak result as A |
| L(d) `u_w` gap | PASS | `Omega_w>0` leaves no massless `u_w` descendant |

Configuration B verdict: `FAIL_C5_LONGITUDINAL_ZERO_MODE`.

Branch B does escape the L(a-ii)/L(b) live-arrow horn: the determinant count gives 2 gapless positive branches, 1 gapped `u_w` branch, and
0 extra gapless `P` branches; closure is algebraic. The costs are (i) the `k^4` dispersion scale above, (ii) the leading traction
provenance becomes `POSTULATED_SURFACE_ELASTICITY`, and (iii) the frozen no-`phi` C5 failure still kills Gate L.

## Section 2.6 Attribution

The four-way no-go materialized in the following computed sense:

- Live massless `P` gives the needed reservoir candidate, but it produces 3 hidden propagating spin waves and the gyroscopic `P-u`
  Hamiltonian is unbounded at low `k`.
- Gapping `P` removes the live modes but leaves a nonzero low-`k` closure residual, `2 mu_R Omega_u`.
- Slaving `P` removes the hidden modes and closes angular momentum, but only with derived
  `mu_eff = mu_R - 2 lambda_Pu`, a `k^4` correction, and traction provenance downgraded to the postulated surface modulus.
- The frozen package has no C5 `phi`; an independent variational `phi` fixture removes the longitudinal mode in the control, but it was not
  frozen in G0. The `phi=u_w` route collides with L(d), because the scalar-potential descendant would have to be massless while `u_w` must
  remain gapped.

Therefore the overall label is `FAIL_COUPLE_STRESS_NOGO`, not a collapse onto L(b) alone.

## Controls

All required controls were real computations and fired:

- Frank-only reference: rank-0 `u` traction, `FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION`.
- Cauchy reference: 3 propagating elastic modes, `FAIL_CAUCHY_STRAY_LONGITUDINAL`.
- No-`phi` C5 branch: constrained physical zero mode, `FAIL_C5_LONGITUDINAL_ZERO_MODE`.
- Independent variational `phi` fixture: removes the longitudinal mode, but is fresh-G0 only.
- Raw `div u=0` projector: fails for no variational provenance.
- Omitted couple-stress reservoir: nonzero closure residual.
- Large-gap `P` leg: nonzero low-`k` closure residual.
- Bent interface: nonzero direct `T_na`; the actual flat wave `u_y=A cos(kx-omega t), u_w=0` gives `v_w=partial_t u_w=0` and hence
  `T_na=0`.
- Frank-only indirect channel: zero induced `P` source; the coupled channel uses
  `P=-lambda_Pu partial_x u_y/(Gamma_P k^2-J_P omega^2)` and the T0 current `J_P(partial_t P)partial_a P`; the flat one-`k` curl is 0,
  while the nonplanar able-to-fail control gives nonzero advective curl.
- Ungapped `u_w`: trips `FAIL_BENDING_MASSLESS_FIFTH_FORCE`.
- Dimensional ablations: dropped `m` from `T_wa` and corrupted `P-u` gradient structure both fired.

## Dimensional Firewall

Result: `PASS`, `ENGINE_AGREE`.

The engines checked 42 expressions inline, including the projected T0 coefficients, the traction operators, the full projected interface
traction `T_na` with `T_wa=m rho v_w v_a`, the principal symbols, `c_gamma^2`, the slaved `k^4` dispersion, Hamiltonian minors, closure
residuals, indirect advective currents, and the coupled `u_w` scalar spectrum.

Able-to-fail dimensional ablations:

- `rho v_w v_a` computed as `L^-2 T^-2` instead of `M L^-2 T^-2`.
- `lambda_Pu P u` computed as `M T^-2` instead of `M L^-1 T^-2`.

## Run Log

Commands run from `/var/projects/toy_physics`:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_35_gateL_sympy.py < /dev/null
timeout 600 math -script software/stage1_solver/tools/pathA_35_gateL.wl < /dev/null
timeout 600 python3 software/stage1_solver/tools/pathA_35_gateL_sympy.py --compare < /dev/null
```

All exited `0`. The comparison wrote `ENGINE_AGREE` over the derived Hessian, slaved reduction, mode counts, leak quantities, sub-hurdle
statuses, and aggregated verdict.

Scripts:

- `software/stage1_solver/tools/pathA_35_gateL_sympy.py`
- `software/stage1_solver/tools/pathA_35_gateL.wl`
