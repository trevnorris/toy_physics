# MT15 M1b BdG + Wall Derivation Report

Generated from `software/stage1_solver/mathematica/mt15_02_bdg_wall_derivation.wls`.

This is a pre-freeze, target-blind machinery run.  The BdG matter sector and wall weak-form sector are derived.  The mixed-Maxwell port sector is explicitly posited, not derived, because the canonical 1-D mixed eigenproblem and Path-A `S_eta^(A)` kernel remain open.

**MACHINERY, NOT PHYSICS:** The exported `varpi_alpha`, `c_alpha`, and `B0/B2/B4` are eigenmodes/couplings of the operator linearized around a NON-EQUILIBRIUM target-blind background (stationary residual norm `243.39250922131095`). They are correctly-derived MACHINERY, NOT physical frequencies/couplings; physical values require the self-consistent background at M1c.

## Commands

```bash
timeout 600 wolframscript -script software/stage1_solver/mathematica/mt15_02_bdg_wall_derivation.wls
timeout 600 python research/pde_audit/scripts/stage_v2_22c_end_to_end_smoke_pipeline.py --solver-output software/stage1_solver/mathematica/runs/mt15_02_bdg_wall_derivation/mt15_02_v2_22b_bdg_wall_packet.json --tol 1e-9 --out-report software/stage1_solver/mathematica/runs/mt15_02_bdg_wall_derivation/mt15_02_v2_22c_pipeline_report.json --out-profile-manifest software/stage1_solver/mathematica/runs/mt15_02_bdg_wall_derivation/mt15_02_v2_22a_profile_manifest.json --out-v21-manifest software/stage1_solver/mathematica/runs/mt15_02_bdg_wall_derivation/mt15_02_v2_21_manifest.json --out-observable-packet software/stage1_solver/mathematica/runs/mt15_02_bdg_wall_derivation/mt15_02_v2_22c_observable_packet.json --out-tolerance-budget software/stage1_solver/mathematica/runs/mt15_02_bdg_wall_derivation/mt15_02_v2_22c_tolerance_budget.json
```

## Background

Geometry reuses the M1a open throat:

`R0(w)=a+(R_exit-a)(3x^2-2x^3)`, with `L=2.0`, `a=1.0`, `R_exit=1.65`.

The stationary matter profile is target-blind engineering smoke, not a full 2-D coupled Newton solution:

`rho0=0.035+0.18 Exp[-(r/(0.72 R0(w)))^2] (1+0.06 Cos[Pi w/L])^2`, `psi0=Sqrt[rho0]`.

The smoke stationary residual norm is `243.39250922131095`.  M1b uses this only as a well-defined background for the operator eigensolve; M1c tightens background self-consistency.

## Canonical Operator Map

Compact source quote, `notes/moving_throat_pde_program_compact.md` L1406-L1428: `i hbar partial_t [delta psi, delta psi*] = L_BdG[...] + C_A + C_eta`, with wall source `delta V_conf = -V_wall'(Sigma0/lc) eta/lc`.

EOS source, compact L578-L582: `P=K rho^5`, `U=(K/4)rho^5`, `h=dU/drho=(5K/4)rho^4`.

Term-by-term implementation:

- Kinetic: `-(hbar^2/(2m))*(laplacian-l(l+1)/r^2)`, with `l=2`, `l(l+1)=6`; matches `p2_tangent.py` L240-L245.
- Confinement: `V_radial*(r/R0)^4 + 0.5*V_axial*((w-L/2)/L)^2`; matches the implemented confinement entering `coupled_branch.py` L263-L270.
- Quintic enthalpy: `h=(5K/4)rho0^4`; matches compact L578-L582 and `physics.py` L46-L47.
- Linearization: `delta(h psi)=h delta_psi + (rho0 dh/drho)(delta_psi+delta_psi*)`, with `dh/drho=5K rho0^3`.
- Chemical potential: `-mu` on the single-particle diagonal; matches `coupled_branch.py` L263-L270.
- Gauge coupling: `q A0` retained; spatial covariant derivative terms are represented in the operator ledger but vanish because this smoke background uses `Ar=Aw=0`.  Source: compact L569-L573 and `coupled_branch.py` L253-L258.
- Wall drive: `delta V_conf=-4 V_radial r^4/R0^5 eta`; source compact L1080-L1085, L1424-L1428, and `p2_tangent.py` L165-L183.

## Derived BdG Spectrum

**MACHINERY, NOT PHYSICS:** The exported `varpi_alpha`, `c_alpha`, and `B0/B2/B4` below are from the correctly assembled operator around the NON-EQUILIBRIUM target-blind background, not physical frequencies/couplings.

Finite BdG matrix: `360 x 360` from `12` radial cells and `15` axial points.  Positive real modes found: `180`; exported modes: `3`.

| mode | varpi | relative eigen-residual | c_alpha | adapter I_eta_phi | lambda_B convention |
|---|---:|---:|---:|---:|---:|
| 1 | 6.4326845472035865 | 8.600892158865581e-14 | 0.18690592450983112 | -0.4698931123942154 | 0.3977626391617057 |
| 2 | 7.695104597592589 | 9.198683525361953e-14 | 0.07599805468245249 | -0.10224563334254178 | 0.7432890011825245 |
| 3 | 10.752603138823837 | 5.711111504737992e-14 | 0.05731756462299828 | 0.10053474922243924 | 0.5701268970809256 |

All exported `varpi_alpha > 0`.  The residuals are direct checks of `L_BdG v_alpha = hbar varpi_alpha v_alpha`, so these modes come from the eigensolve, not from hand-entered frequencies.

`lambda_B` is only the packet normalization-convention factor in `c_alpha = lambda_B * I_eta_phi`.

`c_alpha` here is a defensible realization of the `delta V_conf` wall<->BdG coupling: symplectically-normalized mode density-response projected through the drive onto the wall shape. It is derived as an engineering realization; exact equivalence to the canonical Schur-complement reduced coupling is pending M1c.

Derived BdG moments:

- `B0 = 0.0009701851412361938`
- `B2 = 2.229517295898049e-05`
- `B4 = 5.229950077800426e-07`
- Stieltjes identity: `B0*B4 - B2^2 = 1.0327248218049946e-11 >= 0` by Cauchy-Schwarz for any positive `(c,varpi)`; this is a structural sanity identity that catches a coding bug only, NOT a physics validation gate.

## Derived Wall Coefficients

Basis: `{w/L, (w/L)^2, (w/L)^3}`.  Constitutive placeholders are target-blind:

- `mu_eta = 1+0.08 Cos[Pi w/L]^2`
- `T_w = 1.1+0.10 Sin[Pi w/L]`
- `T_Omega = 0.35+0.03 Cos[2 Pi w/L]`
- `K_eta = 0.45+0.08*256*x^4*(1-x)^4`
- `V_l = K_eta + l(l+1)T_Omega`, `l=2`
- `Y_L(0)=0`

Weak-form matrices:

```text
M =
[[0.6974081534333318, 0.5261122301499979, 0.4229177220612085],
 [0.5261122301499979, 0.4229177220612085, 0.35377392156969134],
 [0.4229177220612085, 0.35377392156969134, 0.3041438720748873]]

K =
[[2.317854347926562, 1.8946128582829593, 1.6358182203438723],
 [1.8946128582829593, 1.8286215370424828, 1.7484519833049588],
 [1.6358182203438723, 1.7484519833049588, 1.7924323068387131]]
```

`M` eigenvalues: `1.3848558225405911`, `0.03916772355384407`, `0.00044620147499199335`; therefore `M>0`.

The packet scalar wall coordinate is the lowest generalized wall shape normalized to the V2 profile weight:

- `K = 1.4830544424409784`
- `M = 0.46742489897033024`
- profile weight norm `= 0.9999999999999999`

## Posited Mixed Ports

Status: `posited_not_derived`.

These values are sane arbitrary placeholders, not target-fit and not reverse-engineered:

- `Omega_U = 3.25`
- `Omega_W = 4.35`
- effective `R = 0.08`
- effective `g_U = 0.18`
- effective `g_W = 0.13`
- `I_eta_u = 0.9999999999999999`
- `I_eta_w = 0.9846958899840685`
- `I_u_w = 0.9846958899840685`
- packet factors: `lambda_U = 0.18000000000000002`, `lambda_W = 0.1320204555764961`, `lambda_R = 0.08124335727784375`
- `Delta_eff = 199.86250624999994 > 1e-12`

The packet metadata includes `mixed_ports_status: "posited_not_derived"`.

## V2 Packet And Residual

Packet:

`software/stage1_solver/mathematica/runs/mt15_02_bdg_wall_derivation/mt15_02_v2_22b_bdg_wall_packet.json`

V2 validation:

- `validation_pass = True`
- `mechanical_pipeline_pass = True`
- `open_gate_pass = True`
- `stability_gate_pass = True`
- `target_packet_pass = False`
- `branch_target_realization_claimed = False`

Loud residual label:

**BdG + wall DERIVED; mixed ports POSITED.  `R_norm`, `R_pole`, `R_P2`, and `R_P4` are target-blind PARTIALLY-POSITED diagnostics, NOT clean falsification tests and NOT physics results. They are partially/largely determined by the POSITED mixed ports (`g_U=0.18`, `g_W=0.13`, `R=0.08`, `Omega_U=3.25`, `Omega_W=4.35`), so they are NOT clean tests of the derived BdG/wall sectors.**

V2 smoke residuals:

- `R_norm = -0.11056372954240365`
- `R_pole = -0.6564301341214106`
- `R_P2 = 2.4149712110438795e-05`
- `R_P4 = 1.2307316834719538e-05`
- `R_tail = 14.625`

The proximity `|R_norm| ~= 0.1106` to the GR target constant `54G c_s^5/(5a^5c^5)=0.110592` is a calibration artifact of the near-zero posited `N0` (`N0 ~= 4.8e-5`, making the `mhat0^2*S_port*N0/D0` transfer term about `2.8e-5`), NOT validation and NOT tuning.

The downstream calibrated direct-coefficient control still passes, proving the V2 pipeline remains mechanically discriminating.
