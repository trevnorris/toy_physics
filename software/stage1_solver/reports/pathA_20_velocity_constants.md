# Path-A 20 velocity constants summary

## Verdicts

- `c_gamma/c_s`: `C_GAMMA_RATIO_UNDERDETERMINED`. The carried ratio is `lambda_gamma=c_gamma/c_s`; `tail=(c/c_s)^3=lambda_gamma^3`.
- `flux_law_verdict`: `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`. No accepted `J_crit` law is produced in this step.
- `hbar` provenance: `HBAR_PROVENANCE_UNDETERMINED`. The `h/2pi` split is bookkeeping-useful for complete windings, not a provenance split.

## S1 sound speed

- Formula: `c_s^2(rho)=(1/m_GNLS)*dP/drho=5*K*rho^4/m_GNLS`.
- Dimension: `L T^-1` machine-checked in SymPy and Mathematica.
- State dependence: c_s(rho) proportional to rho^2; d ln c_s / d ln rho = 2.
- Profile status: rho, v_b, and c_s are one stationary profile through continuity plus quantum-Bernoulli; c_s=1 denotes the asymptotic c_s0 pin.
- Provenance: derived relative to imposed EOS P=K*rho^5, not from an hbar-free microscopic EOS derivation in this step.

Source anchors: `part01_parent_geometry.tex:194-203`; `pde.tex:344-352`.

## S2 velocity structure

| Velocity | Role | Status |
|---|---|---|
| `v_b` | background condensate flow velocity, v_b=(hbar/m_GNLS) grad(theta) in the ungauged lane; profile variable | `[v_b]=L T^-1` checked |
| `c_s` | bulk density/phonon sound speed; profile c_s(x) with asymptotic c_s0 | `[c_s]=L T^-1` checked |
| `c_gamma` | photon/gauge-wave speed from the gauge-wave kinetic operator; brane light-cone speed | `[c_gamma]=L T^-1` checked; ratio to `c_s` underdetermined |

`c=c_gamma` result: c=c_gamma from the massless wave-sector ceiling: omega^2=c_gamma^2 k^2 gives group velocity c_gamma; a trapped transverse mode has omega^2=c_gamma^2(k_parallel^2+k_perp^2), so d omega/dk is bounded by c_gamma and approaches it at high drive.

Bound-mode clock: The trapped-mode rest oscillation is omega0=c_gamma*k_perp from the wave boundary condition. A boosted wave-operator solution has phase exp[-i*omega0*gamma*(t-v*x/c_gamma^2)], so along the packet center x=v*t the internal clock advances at omega0/gamma. No E=m_defect*c_gamma^2 or Compton premise is used.

Constants vs profiles: Constants/input labels in this step: K, m_GNLS, hbar, conserved/no-leakage J label, and asymptotic rho0,c_s0. Profiles: rho(x), v_b(x), c_s(x), and possibly c_gamma(x) until the gauge ratio is fixed.

Mass bridge recorded only as candidate form: `m_defect=alpha_J*hbar*J/c_gamma^2, or alpha_J*h*J_nu/c_gamma^2 for a cycle-count rate; candidate conversion only`. This does not collapse `M` and does not derive `alpha_J`.

The localized Maxwell sources expose `Z(w)`, projection, and measured-vs-flux closure data (`part01_parent_geometry.tex:225-389`; `pde.tex:355-565`). The legacy EM acoustic reuse (`em_fields.tex:149-184`, `482-499`, `692-705`) is a prior, not the kinetic-operator proof required here.

## S2b flux law

Verdict: `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`.

- Conservation statement: In a steady no-leakage region, J=int rho v_b.dSigma is surface-independent, but local v_b can accelerate as rho and area vary.
- No-net-accretion status: `NO_NET_ACCRETION_BC_UNDERIVED`.
- Missing branch data: solved `R0`, `psi0`, `A_M0`, confinement/wall data, quantum-pressure contribution, leakage/topology BC, support/gauge spectra, and overlap/DtN data.
- Consequence: pathA_21 must consume the verdict, not an unconditional choked flux.
- Environment dependence: any solved flux law inherits `P0=K*rho0^5` and `c_s0^2=5K*rho0^4/m_GNLS`; this was dimension-checked but not converted into an accepted choked law.

A conditional ideal Euler-nozzle algebra check was recorded but not accepted as the branch law: with upstream rest flow and no `Q`/`V_conf` variation, `c_s,* / c_s0=3^(-1/2)`, `rho_* / rho0=3^(-1/4)`, and `Jcrit/(rho0 c_s0 A_*)=3^(-3/4)`.

Source anchors: `pde.tex:2515-2566`, `2847-2879`; `brane_bulk_ontology.tex:1998-2042`.

## S3 hbar and h/2pi

Verdict: `HBAR_PROVENANCE_UNDETERMINED`.

- Current action role: hbar is an explicit action/PDE coefficient in the present parent theory.
- Anti-tautology gate: hbar=m_GNLS*c_s0*a is a pin rearrangement unless a is independently fixed by an hbar-free substrate/action relation.
- `h/2pi` assessment: Partially meaningful only in winding/rate bookkeeping: h is the natural action per complete phase winding, while hbar remains the angular phase/PDE coefficient. It does not split charge and mass provenance by itself; J cycle-vs-angular status and the 2*pi placement are deferred.
- `n^2` assessment: A vortex kinetic-energy ladder proportional to kappa^2 would scale like n^2, so higher windings are energetically disfavored; this is recorded as a conditional winding-sector prediction, not a mass spectrum derivation.

Role catalog dimensions checked: circulation `kappa=int v_b dl=h*n/m_GNLS`, phase momentum `p=hbar grad(theta)`, quantum pressure `Q`, and candidate bridge `hbar J/c_gamma^2`.

Source anchors: `part01_parent_geometry.tex:174-219`, `270-289`; `pde.tex:429-469`; `brane_bulk_ontology.tex:668-671`, `1169-1180`.

## Residuals

- `EOS_CLOSURE_IMPOSED`: CARRIED_FORWARD. c_s is derived only relative to the imposed stiff-polytropic EOS; deriving the EOS from a deeper substrate remains outside pathA_20. Source: part01_parent_geometry.tex:194-203; pde.tex:344-352 state P=K*rho^5 and c_s^2=(1/m_GNLS)dP/drho.
- `C_GAMMA_RATIO_UNDERDETERMINED`: BLOCKS_NUMERIC_C_GAMMA_OVER_C_S. Carry lambda_gamma=c_gamma/c_s and tail factor (c/c_s)^3=lambda_gamma^3 into pathA_21/pathA_22 until the gauge/density kinetic operator and localization profile fix it. Source: part01_parent_geometry.tex:225-389 and pde.tex:355-565 give a localized Maxwell action, projection law, Z(w) renormalization, and measured-vs-flux closure issue; em_fields.tex:149-184, 482-499, 692-705 only give the legacy weak-field acoustic reuse.
- `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`: FLUX_LAW_VERDICT. No accepted choked-flux law or nontransonic alternate law is available here; pathA_21 must consume this verdict rather than an unconditional Jcrit. Source: pde.tex:2515-2566 requires a branch data set {R0, psi0, A_M0, wall data, spectra, overlaps, a, c_s, mhat}; pde.tex:2847-2879 states the actual stationary throat profile and DtN data remain branch-dependent.
- `NO_NET_ACCRETION_BC_UNDERIVED`: CARRIED_FORWARD. J_in=J_out requires throat-bottom topology/BC input; Gauss flux conservation is local no-leakage conservation, not a no-net-accretion proof. Source: part01_parent_geometry.tex:298-330; pde.tex:511-539; brane_bulk_ontology.tex:1998-2042.
- `HBAR_FREE_SUBSTRATE_RELATION_MISSING`: BLOCKS_HBAR_EMERGENT. S3 verdict is HBAR_PROVENANCE_UNDETERMINED; hbar remains explicit in the PDE and in the candidate mass bridge. Source: part01_parent_geometry.tex:174-219 and pde.tex:318-408 keep hbar as an action coefficient; pathA_19 only gives the pin relation a=hbar/(m_GNLS*c_s0).
- `H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED`: CARRIED_FORWARD. Use h for complete winding/cycle-count relations and hbar for angular/PDE coefficients; pathA_21 must decide whether J is cycle-rate or angular-rate and where the 2*pi sits in alpha_J. Source: pde.tex:429-469 gives psi=sqrt(rho) exp(i theta) and v_i=(hbar/m_GNLS)partial_i theta; brane_bulk_ontology.tex:668-671 and 1169-1180 treat circulation as quantized/integer-labeled.

## Algebraic harness summary

- Dimensional checks: 21 consistent, 0 inconsistent, 21 total.
- Algebraic checks: 5 consistent, 0 inconsistent, 5 total.
- Acceptance status: `PASS_WITH_NAMED_RESIDUALS`.

