THROAT_DRAIN_DESTABILIZED

**pathA_26 Derrick + Open-System Stability Report**

Computed phase labels: Phase A = `A_CANDIDATE_MIN`, Phase B = `B_RESCUABLE`, Phase C = `C_GENERICALLY_UNSTABLE`.
The top-line follows the directive decision table from those computed labels and the approximation gate.

**Functional And Sources**

The fluid GNLS/Hamiltonian and n=5 EOS are anchored at `notes/inner_throat/inner_throat_hard_mode.md:331`, `notes/inner_throat/inner_throat_hard_mode.md:365`, and `notes/inner_throat/inner_throat_hard_mode.md:1544`.
The scalar wave surrogate and wave Hamiltonian are anchored at `notes/inner_throat/inner_throat_hard_mode.md:459`, `notes/inner_throat/inner_throat_hard_mode.md:479`, and `notes/inner_throat/inner_throat_hard_mode.md:1540`.
The fixed smooth confinement gates/regulators are from `notes/inner_throat/inner_throat_freeze_sheet.md:54`, `notes/inner_throat/inner_throat_freeze_sheet.md:62`, `notes/inner_throat/inner_throat_freeze_sheet.md:70`, `notes/inner_throat/inner_throat_freeze_sheet.md:89`, `notes/inner_throat/inner_throat_freeze_sheet.md:101`, and `notes/inner_throat/inner_throat_freeze_sheet.md:113`.
The 4D tube measures are from `notes/inner_throat/inner_throat_freeze_sheet.md:196`, `notes/inner_throat/inner_throat_freeze_sheet.md:202`, `notes/inner_throat/inner_throat_freeze_sheet.md:206`, and `notes/inner_throat/inner_throat_freeze_sheet.md:213`.
Open-system flux/back-pressure diagnostics are from `notes/inner_throat/4d_next_steps.md:1184`, `notes/inner_throat/4d_next_steps.md:1210`, `notes/inner_throat/4d_next_steps.md:1217`, `notes/inner_throat/4d_next_steps.md:1224`, the fixed-density reservoir template from `notes/inner_throat/4d_next_steps.md:3145`, and the back-pressure observable from `notes/inner_throat/4d_next_steps.md:3246` and `notes/inner_throat/4d_next_steps.md:3262`.
The internal 4D self-flow and D/N ladder are from `notes/lepton_mass_notes.md:188`, `notes/lepton_mass_notes.md:201`, `notes/lepton_mass_notes.md:266`, and `notes/lepton_mass_notes.md:303`; the three distinct flux objects are separated at `notes/lepton_mass_notes.md:1043` and `notes/lepton_mass_notes.md:1047`.

**Object And Closures**

The tested object is `E*(a,L)=Omega_fluid,excess + I*omega(a,L) + P_vac*V + sigma*A`.
The fixed-mu fluid closure uses `mu=U'(rho0)` and the sharp-wall depletion contributes `+P(rho0)*V`; in the execution slice `P(rho0)=0.002`, `P_vac=0.001`, and `sigma=0.003`.
The fixed-action wave branch is `omega=sqrt(c_w^2*(pi^2/a^2+(pi/2)^2/L^2)+mu_in^2)`, with `chi_perp=pi`, `chi_DN=pi/2`, `mu_in=0.2`, and the binding threshold `mu_out=2.5`.

**Forms**

| term | form | a,L exponents/form | sign | source |
| --- | --- | --- | --- | --- |
| fluid_grand_potential_depletion | `P(rho0)*(4*pi/3)*a^3*L` | `[3, 1]` | + | notes/inner_throat/inner_throat_hard_mode.md:1544; notes/inner_throat/inner_throat_freeze_sheet.md:202 |
| vacuum_pressure_geometry | `P_vac*(4*pi/3)*a^3*L` | `[3, 1]` | + | notes/inner_throat/inner_throat_freeze_sheet.md:191; notes/inner_throat/inner_throat_freeze_sheet.md:202 |
| brane_tension_side | `sigma*4*pi*a^2*L` | `[2, 1]` | + | notes/inner_throat/inner_throat_freeze_sheet.md:191; notes/inner_throat/inner_throat_freeze_sheet.md:206 |
| brane_tension_caps | `sigma*(8*pi/3)*a^3` | `[3, 0]` | + | notes/inner_throat/inner_throat_freeze_sheet.md:191; notes/inner_throat/inner_throat_freeze_sheet.md:209 |
| fixed_action_wave | `I*sqrt(c_w^2*(pi^2/a^2+(pi/2)^2/L^2)+mu_in^2)` | `algebraic; asymptotic pieces (-1,0) and (0,-1), not a monomial` | + | notes/inner_throat/inner_throat_hard_mode.md:479; notes/lepton_mass_notes.md:303 |
| phaseB_internal_4D_self_flow | `Phi^2/(8*pi^2*rho*a^2)` | `[-2, 0]` | + | notes/lepton_mass_notes.md:188; notes/lepton_mass_notes.md:201; notes/lepton_mass_notes.md:1047 |
| phaseB_soft_slab | `sigma_ret*(L/d_slab)^2` | `[0, 2]` | + | docs/conceptual_foundation.md:321; docs/conceptual_foundation.md:325 |
| phaseB_hard_slab | `constraint L<=d_slab` | `inequality constraint` | constraint | docs/conceptual_foundation.md:325; docs/conceptual_foundation.md:326 |
| optional_Willmore_side | `lambda_W*(16*pi/9)*L` | `[0, 1]` | + | notes/inner_throat/inner_throat_hard_mode.md:586; notes/inner_throat/4d_next_steps.md:353 |
| phaseC_conductance | `G_c=a^3/L` | `[3, -1]` | + | notes/inner_throat/4d_next_steps.md:1210; notes/inner_throat/inner_throat_freeze_sheet.md:200 |
| phaseC_work_kernel | `G_work=a^2*L/ell0^5` | `[2, 1]` | s in {+,-} enters F_nc | notes/inner_throat/4d_next_steps.md:1224; notes/inner_throat/inner_throat_freeze_sheet.md:206 |

The wave term is intentionally not monomialized.  The scalar surrogate, Maxwell `F^2`, and transverse brane-shear all share the same two-derivative dispersion structure for this scaling test; polarization count changes degeneracy, not the `a,L` powers.  The exception is a future vector/shear calculation with extra boundary-localized gauge constraints, which this necessary-condition test does not include.

**Phase A**

Sharp-wall stationary point: `a*=1.892410005412`, `L*=1.81206770991`, `L/a=0.957544984823`.
Analytic envelope Hessian eigenvalues: `[0.311311733728, 1.616343493965]`.
Virial residuals: `a=0.0`, `L=0.0`.

Coarse fixed-regulator envelope Hessian eigenvalues: `[0.960114275551, 1.942629378481]`.
The numeric grid is a cheap constrained smooth-gate depletion check, not a high-resolution PDE solve; it uses grid `[48, 48]` and relaxed depletion `eta=0.95`.
Fluid-only contrast: `FLUID_ONLY_COLLAPSE_NO_INTERIOR_STATIONARY` because both positive-coupling derivatives are positive for `a,L>0`.

**Phase B**

Phase B label: `B_RESCUABLE`.
The certified open positive region is `{'B_flow': '(0, 0.0001)', 'sigma_ret': '(0, 0.0001)', 'lambda_W': '(0, 0.0001)', 'd_slab': 2.5}`.
The hard slab check gives `interior_not_boundary` for `L <= 2.5`.
The optional Willmore side term uses the averaged side mean curvature `H=2/(3a)`, giving `(16*pi/9)*lambda_W*L` plus fixed-regulator edge corrections, so it is not scale-invariant under independent `(a,L)` variation.

**Phase C**

The slaved drain closure is `Phi_w=kappa_c*(a^3/L)*mu_drive`, `Delta_h=zeta*Phi_w`.
The single work-kernel is `G_work=a^2*L/ell0^5`, so after slaving `F_nc=s*(kappa_c*mu_drive)^2/rho0*(2*a^7/L, a^8/L^2)`.
`K_open=-dF_nc/dq`, `M=diag(rho0*V,rho0*V/4)`, and passive damping is `diag(c_a,c_L)`.

For `s=+1`, the nonzero-drain sample fixed point is `[1.897036263838, 1.824510457739]` and the Jacobian spectrum is `[{'re': -0.017712259, 'im': 0.1814321003}, {'re': -0.017712259, 'im': -0.1814321003}, {'re': -0.0302034748, 'im': 0.1438780504}, {'re': -0.0302034748, 'im': -0.1438780504}]`.
For `s=-1`, the nonzero-drain sample fixed point is `[1.887848509896, 1.800012627417]` and the Jacobian spectrum is `[{'re': -0.0176380603, 'im': 0.188353188}, {'re': -0.0176380603, 'im': -0.188353188}, {'re': -0.0316423622, 'im': 0.1465903422}, {'re': -0.0316423622, 'im': -0.1465903422}]`.
The stable near-zero-drain thresholds are `gcrit_plus=0.0013108765` and `gcrit_minus=0.0060561151`, where `g=(kappa_c*mu_drive)^2/rho0`.

At the box-required center gain `g_open=88.82874001` evaluated at `q*=[1.892410005412, 1.81206770991]`, the Phase-C failure is static DIVERGENCE, not flutter: for `s=+1`, `det(H+K_open)=-66344569.698162` and `lambda_min(sym(H+K_open))=-32837.374848`; for `s=-1`, `det(H+K_open)=-66344104.227761` and `lambda_min(sym(H+K_open))=-6228.919444`.
Because the stiffness has a negative determinant and a negative symmetric-part eigenvalue at the required box gain, passive damping `C>=0` cannot remove the negative-stiffness direction; the damping value only sets timescales, so the damping axis need not be swept to `c_max`.
Robust stable region exists: `False`.
Instability proof: Stable fixed points exist only in a tiny near-zero-drain corner. The most favorable in-box 30%-axis-radius ball center already has kappa_c=mu_drive=3.07, hence g_open=88.8287, far above max(gcrit)=0.00605612; therefore no qualifying ball can be stable.
The conclusion is robust to the magnitude of `M`: `M` changes timescales in the Jacobian, while the no-robust-ball proof uses the steady point and open stiffness thresholds.

**Approximation Gate**

Gate checks: `{'scale_separation': True, 'binding': True, 'candidate_interior_to_numeric_grid': True, 'hessian_sign_agreement': True}`.
Scale separation min ratio: `36.2413541982`.
Binding margin fraction: `0.2466213649`.

**Lambda-Ray Regression**

The reduced-ledger regression reproduces `F=A/a+B/a^2+C*a^3` with exponents `[-1, -2, 3]` and virial `E_w + 2*E_f = 3*E_PV`.
This is only the reduced ledger identity from `notes/lepton_mass_notes.md:50`, `notes/lepton_mass_notes.md:60`, and `notes/lepton_mass_notes.md:72`; it is not a claim that the full 4D continuum has the same Lambda-ray volume scaling.

**Engine Agreement**

Engine agreement status: `AGREE`.
Agreement details: `{'status': 'AGREE', 'mathematica_json': '/var/projects/toy_physics/software/stage1_solver/_scratch/pathA_26_derrick_mathematica.json', 'mismatches': []}`.

**Falsifier**

Named falsifier: `FALSIFIER_C_ROBUST_DRAIN_BALL`.
Input that would flip the top-line: A measured non-conservative drain closure with a 30%-axis-radius LHP-stable gain ball for nonzero drain would change Phase C to C_STABILIZABLE.  For this kernel that is equivalent to reducing the effective open-force normalization below eta_open=4.55037e-06, or otherwise changing the measured work-kernel/Jacobian so the robust ball fits inside the O(1) box.

**Scope**

This is a necessary-condition test over the two collective coordinates `(a,L)` and the named `j=0` wave branch.  A positive conservative result is not an existence proof.  The scalar wave stands in only for the two-derivative scaling of the trapped support sector.  The grid is coarse and regulator-fixed.  The Phase-C drain law is parameterized, not derived from a PDE.  Stability of the `(a,L)` family does not exclude off-family shape or field negative modes.

---

## Interpretation (orchestrator + user, 2026-06-24 — NOT a computed result; read before citing the top-line)

**`THROAT_DRAIN_DESTABILIZED` must NOT be read as "the gravitating particle is killed."** Two distinct findings:

1. **Conservative existence is the solid, real result.** Phase A = `A_CANDIDATE_MIN` is *generic* (an independent scan found an interior
   positive-definite minimum in 75/75 coupling combinations spanning ~2 decades each). The throat-soliton energy balance — trapped-wave
   outward pressure (`ω∝1/a`: squeezing the throat raises the frequency, hence the restoring push) vs brane-tension + vacuum-cost inward —
   **closes**. The throat-soliton is NOT Derrick-forbidden.

2. **The Phase-C instability is at UNPHYSICALLY LARGE drain.** The destabilization is a static divergence (`det(H+K_open)<0`) only ABOVE a
   critical drain `gcrit≈0.006`; for weaker drain the throat is stable (`gcrit>0`). The `DESTABILIZED` top-line came from a directive
   choice (demand robustness over an O(1) drain box, `g` up to ~89) that does NOT reflect that **gravity is the weakest force**. Prior repo
   work shows the drain/gravity sector is the runt (`E_w:E_f:E_PV=11:2:5`, feed ≈11% of rest energy; `g_G` calibrated-not-derived, bare
   `P0~10⁻⁹`; analog horizon only at the transonic limit `v_b→c_s`; trans-brane current `J_w=0` on the exact trapped mode). So a real
   particle's drain is tiny → it sits deep in the **stable corner**; the instability plausibly marks the **black-hole/transonic regime**,
   where the drain *should* dominate (a consistency feature, not a bug).

**The precise gap (never done):** nobody mapped a particle's physical drain strength onto this dimensionless `g`-axis to PROVE
`g_phys ≪ gcrit`. The drain law here is parameterized, not derived; `g≈89` is a synthetic worst-case, not a particle. Closing this
rigorously = the drain-sector derivation (the named falsifier `FALSIFIER_C_ROBUST_DRAIN_BALL`). **Deeper, separate open question:** with
`J_w=0` on the exact mode, whether the "gravity = drain" mechanism is even *nonzero / right-sized* is itself unbuilt (needs the dynamic
AC→DC rectification specced-not-built in `notes/inner_throat/4d_next_steps.md`). Net: conservative throat-soliton ✓ exists; stable at the
physical (tiny) drain ✓ (structurally — pin pending); the real frontier is the **drain sector itself**.
