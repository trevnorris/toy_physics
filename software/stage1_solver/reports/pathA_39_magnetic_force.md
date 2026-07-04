# pathA_39 Stage 2 Magnetic Force

Computed headline: `MAGNETIC_FORCE_DERIVED` with dual-engine `ENGINE_AGREE`.

The Stage-2 interaction is now obtained by integrating out the transverse and longitudinal fields.  The `R` powers are measured from the resulting expressions, not supplied to the classifier.  The result is conditional on the declared Stage-1 moving-source amplitudes `q_A^T=Nu*aT*s` and `q_L=Nu*aL*s`; the values of `aT` and `aL` remain sim-deferred.

## Derivation

In the static projector eigenbasis, the quadratic operators and their inverses are:

```text
O_T = k2*muR,  G_T = 1/(k2*muR)
O_L = k2*rhoB0**2/chiC,  G_L = chiC/(k2*rhoB0**2)
```

The exchange integral used by both engines is:

```text
U_12=-int d^3k/(2*pi)^3 j1(-k).G(k).j2(k)
F^-1[1/k^2] = 1/(4*pi*R)
F^-1[1/k^4] seed = -R/(8*pi)
F^-1[k_i k_j/k^4] contracted with V1,V2 = -(Aprod - Ddot)/(8*pi*R)
```

## Kernel

With `R=X2-X1`, `n=R/R`, `D=V1.V2`, and `A=(V1.n)(V2.n)`:

```text
K_T compact = (D + A)/(8*pi*R)
K_L compact = (D - A)/(8*pi*R)
U_12 compact = -s1*s2*Nu^2/(8*pi*R)*[aT^2*(D+A)/mu_R + aL^2*(D-A)/B_eff]
F_12 compact = -grad_R U_12, derived symbolically from the computed U_12

K_T = (Aprod + Ddot)/(8*pi*R)
K_L = -(Aprod - Ddot)/(8*pi*R)
U_T = -Nu**2*aT**2*s1*s2*(Aprod + Ddot)/(8*pi*R*muR)
U_L = Nu**2*aL**2*chiC*s1*s2*(Aprod - Ddot)/(8*pi*R*rhoB0**2)
U_12 = Nu**2*s1*s2*(Aprod*aL**2*chiC*muR - Aprod*aT**2*rhoB0**2 - Ddot*aL**2*chiC*muR - Ddot*aT**2*rhoB0**2)/(8*pi*R*muR*rhoB0**2)
F_12 = [['Nu**2*s1*s2*(-aL**2*chiC*muR*(-3*Aprod*nx + Ddot*nx + arad*v2x + brad*v1x) + aT**2*rhoB0**2*(-3*Aprod*nx - Ddot*nx + arad*v2x + brad*v1x))/(8*pi*R**2*muR*rhoB0**2)'], ['Nu**2*s1*s2*(-aL**2*chiC*muR*(-3*Aprod*ny + Ddot*ny + arad*v2y + brad*v1y) + aT**2*rhoB0**2*(-3*Aprod*ny - Ddot*ny + arad*v2y + brad*v1y))/(8*pi*R**2*muR*rhoB0**2)'], ['Nu**2*s1*s2*(-aL**2*chiC*muR*(-3*Aprod*nz + Ddot*nz + arad*v2z + brad*v1z) + aT**2*rhoB0**2*(-3*Aprod*nz - Ddot*nz + arad*v2z + brad*v1z))/(8*pi*R**2*muR*rhoB0**2)']]
```

Reinserting the same derived Green tensor into the static EOM gives the consistency residual `U_eom - U_integrate = 0`.

Measured falloff: potential `R^-1`, point force `R^-2`.

## Sign Diagnostic

For side-by-side parallel throats (`V_i.n=0`), the scalar-to-transverse magnitude ratio is `aL**2*chiC*muR/(aT**2*rhoB0**2)`.  It is not a cancellation ratio: both stable channels have the same attractive sign for like sources.

Transverse like-current coefficient: `-Nu**2*aT**2/(8*pi*R**2*muR)`.
Longitudinal like-current coefficient: `-Nu**2*aL**2*chiC/(8*pi*R**2*rhoB0**2)`.
Total like-current coefficient: `-Nu**2*(aL**2*chiC*muR + aT**2*rhoB0**2)/(8*pi*R**2*muR*rhoB0**2)`.

Landing: `NO_CANCELLATION_BOTH_CHANNELS_ATTRACTIVE`.  Outward radial force means repulsion; inward means attraction.

| case | like charge, like currents | like charge, opposite currents | opposite charge, like currents | opposite charge, opposite currents |
|---|---|---|---|---|
| additive T+L | `ATTRACT` | `REPEL` | `REPEL` | `ATTRACT` |

The transverse `u_T` contribution gives the magnetic like-current attraction.  The longitudinal `u_L` contribution is an unavoidable attractive scalar-current admixture under the same stable-sign assumptions.

## Lorentz-Form Diagnostic

`PREFERRED_FRAME_UNLESS_cE_EQUALS_cGAMMA`.  The residual is `(cE**2*rhoBr - muR)/rhoBr`, so a Lorentz-force-form rewrite requires the additional condition `c_E=c_gamma`.  This is reported as a diagnostic, not a Stage-2 fail.

## Controls

| control | status | verdict/result |
|---|---:|---|
| `muR_propagator_perturbation` | `FIRED` | `RE_DERIVED_RESPONSE_CHANGED` |
| `source_projection_scale_perturbation` | `FIRED` | `RE_DERIVED_RESPONSE_CHANGED` |
| `projection_tensor_perturbation` | `FIRED` | `RE_DERIVED_RESPONSE_CHANGED` |
| `propagator_functional_perturbation_kminus4` | `FIRED` | `RE_DERIVED_FALLOFF_CHANGED` |
| `target_readback_fixture` | `FIRED` | `FAIL_TARGET_READBACK` |
| `wrong_falloff_fixture` | `FIRED` | `FAIL_WRONG_FALLOFF` |
| `noncompact_source_fixture` | `FIRED` | `FAIL_WRONG_FALLOFF` |
| `V_equals_zero` | `FIRED` | `NO_MOVING_SOURCE` |
| `neutral_plus_minus_composite` | `FIRED` | `ZERO_MONOPOLE_CURRENT_SOURCE` |
| `charge_flip_s_to_minus_s` | `FIRED` | `SOURCE_SIGN_FLIPS` |
| `ghost_wrong_sign_transverse` | `FIRED` | `MU_R_TO_MINUS_MU_R_REDERIVED_FLIPS_ATTRACTION_TO_REPULSION` |

## Dimensional Firewall

Units-restored checks passed for `8` expressions.  The able-to-fail ablations fired for omitted velocity, omitted compact brane-source normalization, mass-current/charge-current mixing, and using `c_s` where the imported transverse `c_gamma` is required.

## Provenance

Imported:
- `B_eff`: rho_B0^2/chi_c
- `c_gamma_squared`: mu_R/rho_br
- `pathA_36_transverse_sector`: 1/2 rho_br dot(u_T)^2 - 1/2 mu_R (curl u_T)^2
- `q_h`: 2*QE*tanh(b/ell)/b
- `c_E`: cE from pathA_38 dynamic Green exp(I*R*omega/cE)/(4*pi*R)
- `f_h`: [['1/(ell*cosh(w/ell)**2)'], ['1/(ell*cosh(w/ell)**2)'], ['0']]
- `f_u`: pathA_36 shear profile f_u(w); normalization kept as Nu
- `M_h`: positive imported property from pathA_38; not used to tune the Stage-2 vector kernel
- `stage0_plus_1`: pathA_39 scalar-admixture screen with q_A^T=Nu*aT*sCharge and q_L=Nu*aL*sCharge

Declared / parameterized / deferred:
- `imported_exact`: B_eff=rho_B0^2/chi_c; c_gamma^2=mu_R/rho_br; u_T MacCullagh shear sector from pathA_36; q_h=2*QE*tanh(b/ell)/b; c_E=cE; Stage 0+1 q_A^T and q_L projection form
- `declared_stage2`: two compact brane-localized moving throats; j_T=q_A^T V with q_A^T=Nu*aT*s; j_L=q_L V with q_L=Nu*aL*s; static O(V1 V2) exchange limit in the medium rest frame
- `sim_deferred_values`: aT; aL; bulk compactness and finite-mouth form factors; anisotropy tensor; operator motion perturbation deltaO_V; absolute hierarchy magnitude

## Dual Engine

`ENGINE_AGREE` over `73` compared quantities, including the inverted propagators, inverse-Fourier kernels, Green-tensor consistency residual, force components, measured falloff powers, corrected sign codes, dimensional ablation count, and target-readback controls.

Run commands:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_39_magnetic_force_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_39_magnetic_force.wl
timeout 600 python3 software/stage1_solver/tools/pathA_39_magnetic_force_sympy.py --compare
```
