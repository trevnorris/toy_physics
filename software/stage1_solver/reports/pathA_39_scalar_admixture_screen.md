# pathA_39 Stage 0+1 Scalar-Admixture Screen

Computed headline: `FAIL_OBSERVABLE_SCALAR_ADMIXTURE`.

Robust floor: `FAIL_EXTRA_H_BRANON` is import-forced by `M_h>0` and pathA_38 `q_h != 0`; its `h`-branon charge residue is `4*QE**2*tanh(b/ell)**2/(Mh*b**2)` and it is independent of `a_L`.  The `FAIL_OBSERVABLE_SCALAR_ADMIXTURE` headline additionally includes the density-admixture residue proportional to `q_L^2`, whose `OBSERVABLE` upgrade rests on the sim-deferred `a_L != 0` branch.

The decisive Stage 0+1 branch uses the pathA_36 density block with imported `B_eff=rho_B0^2/chi_c`, the pathA_38 `h` block with `Kh=Mh*cE^2`, and the Stage-1 projected `q_L=sCharge*aL*Nu`.  The scalar block is not tuned to Maxwell.

## Scalar Block

`Phi=(u_L,h)` and `x=omega^2/k^2`.

```text
D_x = [['rhoBr*x - rhoB0**2/chiC', '-Chu'], ['-Chu', '-Mh*cE**2 + Mh*x']]
G_scalar_x = [['chiC*(-Mh*cE**2 + Mh*x)/(-Chu**2*chiC - Mh*cE**2*chiC*rhoBr*x + Mh*cE**2*rhoB0**2 + Mh*chiC*rhoBr*x**2 - Mh*rhoB0**2*x)', 'Chu*chiC/(-Chu**2*chiC - Mh*cE**2*chiC*rhoBr*x + Mh*cE**2*rhoB0**2 + Mh*chiC*rhoBr*x**2 - Mh*rhoB0**2*x)'], ['Chu*chiC/(-Chu**2*chiC - Mh*cE**2*chiC*rhoBr*x + Mh*cE**2*rhoB0**2 + Mh*chiC*rhoBr*x**2 - Mh*rhoB0**2*x)', 'chiC*(rhoBr*x - rhoB0**2/chiC)/(-Chu**2*chiC - Mh*cE**2*chiC*rhoBr*x + Mh*cE**2*rhoB0**2 + Mh*chiC*rhoBr*x**2 - Mh*rhoB0**2*x)']]
A_qq = (4*Chu*QE*b*chiC*qL*tanh(b/ell) - Mh*b**2*cE**2*chiC*qL**2 + Mh*b**2*chiC*qL**2*x + 4*QE**2*chiC*rhoBr*x*tanh(b/ell)**2 - 4*QE**2*rhoB0**2*tanh(b/ell)**2)/(b**2*(-Chu**2*chiC - Mh*cE**2*chiC*rhoBr*x + Mh*cE**2*rhoB0**2 + Mh*chiC*rhoBr*x**2 - Mh*rhoB0**2*x))
```

## Decisive Residues (C_hu -> 0 Decoupled Limit)

| pole | speed^2 | residue to charge | residue to mass |
|---|---:|---:|---:|
| `density_cs` | `rhoB0**2/(chiC*rhoBr)` | `Nu**2*aL**2*sCharge**2/rhoBr` | `qM**2/rhoBr` |
| `h_branon` | `cE**2` | `4*QE**2*tanh(b/ell)**2/(Mh*b**2)` | `0` |

This table is the `C_hu -> 0` decoupled limit.  The full mixed-eigenpole residues are recorded under `poles_generic_mixing` (`generic_poles` in the engines), and the `root_minus`/`root_plus` charge residues are included in the dual-engine comparison.  The density residue is nonzero for the parameterized Stage-1 `q_L`, while the `h` pole supplies the import-forced `FAIL_EXTRA_H_BRANON` floor.

## Controls

| control | status | verdict/result |
|---|---:|---|
| `inject_qL_epsilon` | `FIRED` | `FAIL_OBSERVABLE_SCALAR_ADMIXTURE` |
| `closed_steady_current_wire_limit` | `FIRED` | `WIRE_LIMIT_NO_LONGITUDINAL_RESPONSE` |
| `B_eff_positive_vs_zero` | `FIRED` | `FAIL_OBSERVABLE_SCALAR_ADMIXTURE -> FAIL_EXTRA_H_BRANON` |
| `Mh_positive_qh_nonzero` | `FIRED` | `FAIL_EXTRA_H_BRANON` |
| `cE_import_match` | `FIRED` | `IMPORT_MATCH` |
| `mixing_on_with_derived_qL_zero` | `FIRED` | `FAIL_OBSERVABLE_SCALAR_ADMIXTURE` |

## Stage 1 Projection

`q_A^T = Nu*aT*sCharge` and `q_L = Nu*aL*sCharge` from the declared ansatz.  The coefficients `aT`, `aL`, bulk compactness, anisotropy, and operator motion perturbation remain conditional/sim-deferred.

## Provenance

Imported:
- `B_eff`: rho_B0^2/chi_c
- `c_gamma_squared`: mu_R/rho_br
- `q_h`: 2*QE*tanh(b/ell)/b
- `c_E`: cE from pathA_38 dynamic Green exp(I*R*omega/cE)/(4*pi*R)
- `f_h`: [['1/(ell*cosh(w/ell)**2)'], ['1/(ell*cosh(w/ell)**2)'], ['0']]
- `f_u`: pathA_36 shear profile f_u(w); normalization kept as Nu because the report does not give a closed profile
- `M_h`: positive imported property from pathA_38; coefficient kept as symbolic Mh>0

Declared / parameterized / deferred:
- `imported_exact`: B_eff=rho_B0^2/chi_c; c_gamma^2=mu_R/rho_br; q_h=2*QE*tanh(b/ell)/b; c_E=cE
- `imported_positivity_parameterized`: M_h>0 as symbolic Mh; f_u normalization Nu
- `declared_scan_parameters`: q_L from Stage-1 ansatz coefficient aL; C_hu scalar mixing
- `sim_deferred_values`: aT; aL; bulk compactness; operator motion perturbation deltaO_V

## Dual Engine

`ENGINE_AGREE` over `41` compared quantities, including `generic_root_minus_charge_residue` and `generic_root_plus_charge_residue`.

Run commands:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_39_scalar_admixture_screen_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_39_scalar_admixture_screen.wl
timeout 600 python3 software/stage1_solver/tools/pathA_39_scalar_admixture_screen_sympy.py --compare
```
