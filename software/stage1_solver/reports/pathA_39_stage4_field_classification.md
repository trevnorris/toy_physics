# pathA_39 Stage 4 Field-Coupling Classification

Computed headline: `primary=FIELD_SCALAR_VECTOR_DEPARTURE` with flags `{'scalar_sector_stable': True, 'density_charge_coupled': True, 'operator_parity_contamination': True}` and dual-engine `ENGINE_AGREE`.

This is a Stage-4 field-content classification from the assembled inverse propagator, not a relabel of Stages 0-3.  The healthy real-branch departure is conditional on the computed scalar stability bound `B_eff K_h - C_hu^2 > 0`.

## Assembled Q

Coordinates: `Phi=(u_T1,u_T2,u_L,h)`.

```text
Q(omega,k) = [['-k**2*muR + omega**2*rhoBr', '0', '0', '0'], ['0', '-k**2*muR + omega**2*rhoBr', '0', '0'], ['0', '0', 'omega**2*rhoBr - k**2*rhoB0**2/chiC', '-Chu*k**2'], ['0', '0', '-Chu*k**2', '-Mh*cE**2*k**2 + Mh*omega**2']]
J_charge = [['Nu*aT*sCharge'], ['Nu*aTp*sCharge'], ['Nu*aL*sCharge'], ['2*QE*tanh(b/ell)/b']]
J_q record = {'q_A_T_components': ['Nu*aT*sCharge', 'Nu*aTp*sCharge'], 'q_L': 'Nu*aL*sCharge', 'q_h': '2*QE*tanh(b/ell)/b', 'q_M': 'qM', 'mass_source_is_not_counted_as_charge_residue': True}
R=J^T Q^-1 J = Nu**2*sCharge**2*(aT**2 + aTp**2)/(-k**2*muR + omega**2*rhoBr) + chiC*(4*Chu*Nu*QE*aL*sCharge*tanh(b/ell)/b + Nu**2*aL**2*sCharge**2*(-Mh*cE**2 + Mh*x) + 4*QE**2*(rhoBr*x - rhoB0**2/chiC)*tanh(b/ell)**2/b**2)/(k**2*(-Chu**2*chiC - Mh*cE**2*chiC*rhoBr*x + Mh*cE**2*rhoB0**2 + Mh*chiC*rhoBr*x**2 - Mh*rhoB0**2*x))
```

The mass channel `q_M` is recorded as a separate source to `u_L`; it is not counted as charge residue.

## DOF And Constraints

| branch | kinetic rank | longitudinal class | first class | second class | transverse | u_L | h | total physical DOF |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| `real_provenance_fixed` | `4` | `SECOND_CLASS_PAIR` | `0` | `2` | `2` | `1` | `1` | `4` |
| `maxwell_counterfactual` | `2` | `FIRST_CLASS_MAXWELL_CHAIN` | `2` | `0` | `2` | `0` | `0` | `2` |
| `clean_coexistence` | `4` | `SECOND_CLASS_PAIR` | `0` | `2` | `2` | `1` | `1` | `4` |
| `aL_to_0` | `4` | `SECOND_CLASS_PAIR` | `0` | `2` | `2` | `1` | `1` | `4` |
| `large_C_hu` | `4` | `SECOND_CLASS_PAIR` | `0` | `2` | `2` | `1` | `1` | `4` |

The real branch uses the imported pathA_36 `SECOND_CLASS_PAIR` consequence for `u_L`.  The Maxwell counterfactual removes the `h` block and imposes the imported `FIRST_CLASS_MAXWELL_CHAIN`, giving the required DOF drop from 4 to 2.

## Scalar Stability

`det(stiffness)/k^4 = -Chu**2 + Mh*cE**2*rhoB0**2/chiC`.

Finite-`k` stiffness eigenvalues: `k**2*(Mh*cE**2 - sqrt(4*Chu**2 + (-Mh*cE**2 + rhoB0**2/chiC)**2) + rhoB0**2/chiC)/2` and `k**2*(Mh*cE**2 + sqrt(4*Chu**2 + (-Mh*cE**2 + rhoB0**2/chiC)**2) + rhoB0**2/chiC)/2`.

Real branch: `scalar_sector_stable=True` from `computed det(stiffness)/k^4=-Chu**2 + Mh*cE**2*rhoB0**2/chiC is positive; with positive diagonal stiffness entries this gives both finite-k stiffness eigenvalues > 0`.

Large-`C_hu` control: `det(stiffness)/k^4 = -3*Mh*cE**2*rhoB0**2/chiC`, `scalar_sector_stable=False`, class `FIELD_SCALAR_SECTOR_UNSTABLE`.

In the `C_hu -> 0` limit the scalar roots reduce to the density speed `B_eff/rho_br` and the h speed `c_E^2`; Stage 4 keeps the mixed eigenroots instead of importing the decoupled residues.

## Poles And Residues

| branch | pole | speed^2 | charge residue | mass residue |
|---|---|---:|---:|---:|
| `real_provenance_fixed` | `u_T1` | `muR/rhoBr` | `Nu**2*aT**2*sCharge**2/rhoBr` | `` |
| `real_provenance_fixed` | `u_T2` | `muR/rhoBr` | `Nu**2*aTp**2*sCharge**2/rhoBr` | `` |
| `real_provenance_fixed` | `scalar_root_minus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `chiC*(4*Chu*Nu*QE*aL*sCharge*tanh(b/ell)/b + Nu**2*aL**2*sCharge**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr)) + 4*QE**2*(-rhoB0**2/chiC + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh))*tanh(b/ell)**2/b**2)/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |
| `real_provenance_fixed` | `scalar_root_plus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `chiC*(4*Chu*Nu*QE*aL*sCharge*tanh(b/ell)/b + Nu**2*aL**2*sCharge**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr)) + 4*QE**2*(-rhoB0**2/chiC + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh))*tanh(b/ell)**2/b**2)/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |
| `maxwell_counterfactual` | `u_T1` | `muR/rhoBr` | `Nu**2*aT**2*sCharge**2/rhoBr` | `` |
| `maxwell_counterfactual` | `u_T2` | `muR/rhoBr` | `Nu**2*aTp**2*sCharge**2/rhoBr` | `` |
| `clean_coexistence` | `u_T1` | `muR/rhoBr` | `Nu**2*aT**2*sCharge**2/rhoBr` | `` |
| `clean_coexistence` | `u_T2` | `muR/rhoBr` | `Nu**2*aTp**2*sCharge**2/rhoBr` | `` |
| `clean_coexistence` | `scalar_root_minus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `0` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |
| `clean_coexistence` | `scalar_root_plus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `0` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |
| `aL_to_0` | `u_T1` | `muR/rhoBr` | `Nu**2*aT**2*sCharge**2/rhoBr` | `` |
| `aL_to_0` | `u_T2` | `muR/rhoBr` | `Nu**2*aTp**2*sCharge**2/rhoBr` | `` |
| `aL_to_0` | `scalar_root_minus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `4*QE**2*chiC*(-rhoB0**2/chiC + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh))*tanh(b/ell)**2/(b**2*(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))))` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |
| `aL_to_0` | `scalar_root_plus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `4*QE**2*chiC*(-rhoB0**2/chiC + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh))*tanh(b/ell)**2/(b**2*(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))))` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(4*Chu**2*Mh*rhoBr + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |
| `large_C_hu` | `u_T1` | `muR/rhoBr` | `Nu**2*aT**2*sCharge**2/rhoBr` | `` |
| `large_C_hu` | `u_T2` | `muR/rhoBr` | `Nu**2*aTp**2*sCharge**2/rhoBr` | `` |
| `large_C_hu` | `scalar_root_minus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `chiC*(8*sqrt(Mh)*Nu*QE*aL*cE*rhoB0*sCharge*tanh(b/ell)/(b*sqrt(chiC)) + Nu**2*aL**2*sCharge**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr)) + 4*QE**2*(-rhoB0**2/chiC + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh))*tanh(b/ell)**2/b**2)/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC - sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |
| `large_C_hu` | `scalar_root_plus` | `(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh*rhoBr)` | `chiC*(8*sqrt(Mh)*Nu*QE*aL*cE*rhoB0*sCharge*tanh(b/ell)/(b*sqrt(chiC)) + Nu**2*aL**2*sCharge**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr)) + 4*QE**2*(-rhoB0**2/chiC + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*Mh))*tanh(b/ell)**2/b**2)/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` | `chiC*qM**2*(-Mh*cE**2 + (Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2))/(2*rhoBr))/(-Mh*cE**2*chiC*rhoBr - Mh*rhoB0**2 + chiC*(Mh*cE**2*rhoBr + Mh*rhoB0**2/chiC + sqrt(16*Mh**2*cE**2*rhoB0**2*rhoBr/chiC + (Mh*cE**2*rhoBr - Mh*rhoB0**2/chiC)**2)))` |

The clean branch keeps `q_A^T` transverse charge but sets `q_L=q_h=0`; the two scalar charge residues compute to zero.  The `a_L->0` branch sets `q_L=0` but keeps `q_h!=0`, so the computed h-only pole charge residue keeps the scalar-vector departure while the computed density-only pole charge residue leaves `density_charge_coupled=false`.

## Controls

| control | status | derived class/result | recorded artifacts |
|---|---:|---|---|
| `real_provenance_fixed` | `FIRED` | `FIELD_SCALAR_VECTOR_DEPARTURE` | `Q=[['-k**2*muR + omega**2*rhoBr', '0', '0', '0'], ['0', '-k**2*muR + omega**2*rhoBr', '0', '0'], ['0', '0', 'omega**2*rhoBr - k**2*rhoB0**2/chiC', '-Chu*k**2'], ['0', '0', '-Chu*k**2', '-Mh*cE**2*k**2 + Mh*omega**2']]`, `J=[['Nu*aT*sCharge'], ['Nu*aTp*sCharge'], ['Nu*aL*sCharge'], ['2*QE*tanh(b/ell)/b']]`, `dof=4`, `class=FIELD_SCALAR_VECTOR_DEPARTURE` |
| `maxwell_counterfactual` | `FIRED` | `FIELD_EXACT_MAXWELL_STRUCTURE` | `Q=[['-k**2*muR + omega**2*rhoBr', '0', '0'], ['0', '-k**2*muR + omega**2*rhoBr', '0'], ['0', '0', '0']]`, `J=[['Nu*aT*sCharge'], ['Nu*aTp*sCharge'], ['0']]`, `dof=2`, `class=FIELD_EXACT_MAXWELL_STRUCTURE` |
| `clean_coexistence` | `FIRED` | `FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY` | `Q=[['-k**2*muR + omega**2*rhoBr', '0', '0', '0'], ['0', '-k**2*muR + omega**2*rhoBr', '0', '0'], ['0', '0', 'omega**2*rhoBr - k**2*rhoB0**2/chiC', '-Chu*k**2'], ['0', '0', '-Chu*k**2', '-Mh*cE**2*k**2 + Mh*omega**2']]`, `J=[['Nu*aT*sCharge'], ['Nu*aTp*sCharge'], ['0'], ['0']]`, `dof=4`, `class=FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY` |
| `aL_to_0` | `FIRED` | `FIELD_SCALAR_VECTOR_DEPARTURE` | `Q=[['-k**2*muR + omega**2*rhoBr', '0', '0', '0'], ['0', '-k**2*muR + omega**2*rhoBr', '0', '0'], ['0', '0', 'omega**2*rhoBr - k**2*rhoB0**2/chiC', '-Chu*k**2'], ['0', '0', '-Chu*k**2', '-Mh*cE**2*k**2 + Mh*omega**2']]`, `J=[['Nu*aT*sCharge'], ['Nu*aTp*sCharge'], ['0'], ['2*QE*tanh(b/ell)/b']]`, `dof=4`, `class=FIELD_SCALAR_VECTOR_DEPARTURE` |
| `large_C_hu` | `FIRED` | `FIELD_SCALAR_SECTOR_UNSTABLE; restored bound -> FIELD_SCALAR_VECTOR_DEPARTURE` | `Q=[['-k**2*muR + omega**2*rhoBr', '0', '0', '0'], ['0', '-k**2*muR + omega**2*rhoBr', '0', '0'], ['0', '0', 'omega**2*rhoBr - k**2*rhoB0**2/chiC', '-2*sqrt(Mh)*cE*k**2*rhoB0/sqrt(chiC)'], ['0', '0', '-2*sqrt(Mh)*cE*k**2*rhoB0/sqrt(chiC)', '-Mh*cE**2*k**2 + Mh*omega**2']]`, `J=[['Nu*aT*sCharge'], ['Nu*aTp*sCharge'], ['Nu*aL*sCharge'], ['2*QE*tanh(b/ell)/b']]`, `dof=4`, `class=FIELD_SCALAR_SECTOR_UNSTABLE` |
| `import_fidelity` | `FIRED` | `{'B_eff_corrupt': 'FIELD_CLASSIFICATION_UNDERDETERMINED', 'K_h_corrupt': 'FIELD_CLASSIFICATION_UNDERDETERMINED', 'c_E_corrupt': 'FIELD_CLASSIFICATION_UNDERDETERMINED'}` | B_eff/K_h/c_E corruptions each run through the extractor; branch summaries below |
| `dof_count_discriminator` | `FIRED` | `real 4 -> Maxwell 2` | real and Maxwell branch DOF records |

Import-fidelity corruption branch summaries:

| branch | imported value corrupted | assembled scalar Q[2,2] | assembled scalar block | class |
|---|---|---:|---:|---|
| `import_fidelity_B_eff_corrupt` | `B_eff` | `-k**2*(deltaB + rhoB0**2/chiC) + omega**2*rhoBr` | `[['-deltaB + rhoBr*x - rhoB0**2/chiC', '-Chu'], ['-Chu', '-Mh*cE**2 + Mh*x']]` | `FIELD_CLASSIFICATION_UNDERDETERMINED` |
| `import_fidelity_K_h_corrupt` | `K_h` | `omega**2*rhoBr - k**2*rhoB0**2/chiC` | `[['rhoBr*x - rhoB0**2/chiC', '-Chu'], ['-Chu', '-Mh*cE**2 + Mh*x - deltaKh']]` | `FIELD_CLASSIFICATION_UNDERDETERMINED` |
| `import_fidelity_c_E_corrupt` | `c_E` | `omega**2*rhoBr - k**2*rhoB0**2/chiC` | `[['rhoBr*x - rhoB0**2/chiC', '-Chu'], ['-Chu', '-Mh*cBad**2 + Mh*x']]` | `FIELD_CLASSIFICATION_UNDERDETERMINED` |

## Provenance

Imported blocks and values:
- `B_eff`: rho_B0^2/chi_c from pathA_36 finite-compressibility branch
- `c_gamma_squared`: mu_R/rho_br from pathA_36 transverse sector
- `longitudinal_real_constraint_class`: SECOND_CLASS_PAIR with 1 physical longitudinal DOF
- `longitudinal_maxwell_constraint_class`: FIRST_CLASS_MAXWELL_CHAIN with 0 physical longitudinal DOF and 2 first-class constraints
- `q_A_T`: q_A^T components parameterized as Nu*aT*sCharge and Nu*aTp*sCharge
- `q_L`: Nu*aL*sCharge from Stage 0+1 source projection; a_L sim-deferred
- `q_h`: 2*QE*tanh(b/ell)/b
- `q_M`: mass source to u_L, recorded separately from charge residues
- `M_h`: positive symbolic Mh from pathA_38 branon zero-mode normalization
- `c_E`: cE from pathA_38 dynamic Green exp(I*R*omega/cE)/(4*pi*R)
- `K_h`: M_h*c_E^2
- `C_hu`: sim-deferred Stage 0+1 scalar mixing coefficient Chu
- `magnetic_sign`: Stage 2 source of truth: transverse and longitudinal like-current channels attract
- `operator_parity_contamination`: true, magnitude sim-deferred, from Stage 3

Declared / sim-deferred:
- `declared_stage4`: one assembled Q(omega,k) over (u_T1,u_T2,u_L,h); charge source vector decomposed as q_A^T, q_L, q_h with q_M tracked as mass-channel source; real branch uses the imported SECOND_CLASS_PAIR longitudinal class
- `conditional_or_sim_deferred`: C_hu stability bound C_hu^2 < B_eff K_h; a_L and density charge coupling magnitude; operator-parity contamination magnitude; c_E=c_gamma and lambda_gamma knit

## Dimensional Firewall

Passed `14` assembled-action/source checks.  The able-to-fail ablations fired for missing `k^2` in `C_hu`, using charge density for `q_L`, and counting `q_M` as a charge source.

## Dual Engine

`ENGINE_AGREE` over `77` independently derived and compared feature quantities: assembled `Q`, `J`, kinetic rank, imported Dirac constraint consequences, DOF per sector, scalar stability sign/eigenvalue codes, pole speeds, `J^T Q^-1 J`, residues, residue-derived feature flags, control classes, and import-fidelity guards.

Run commands:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_39_stage4_field_classification_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_39_stage4_field_classification.wl
timeout 600 python3 software/stage1_solver/tools/pathA_39_stage4_field_classification_sympy.py --compare
```

## Remediation

REMEDIATED: fixes 1-5, verdict/numbers unchanged.  The large-`C_hu` control carries the stable fixture label but is classified unstable from the computed determinant sign; Mathematica now derives the discrete integers before comparison; both engines use the same classifier predicates; scalar/density/h flags are residue-derived; expected-landing checks are non-fatal records.

## Honest Scope

This closes the magnetism sector as a spec-level field-content skeleton.  It does not perform the nonlinear throat solve, the `lambda_gamma` knit, the `c_E=c_gamma` Lorentz-cone test, or the simulation-deferred operator-contamination magnitude.  It does not un-earn Stages 0-3.
