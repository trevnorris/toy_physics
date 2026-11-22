# Correcting the force on a slowly moving draining mouth

**The earlier push calculation overestimated the resistance by a factor of 5,978 in the largest finite box at the lowest tested frequency.** Its pressure integral used a one-sided field gradient at the sharp junction where the flare opens into the receiving bulk. Measuring the actuator force from carrier momentum conservation brings push, steady tow and cycle-averaged work into agreement on the same rigid geometry. Each box kept its original boundary conditions. The comparison reused all 22 archived harmonic fields. [Original readings](mouth_rigid_drag_table.md) · [Corrected diagnostic plot](figures/mouth_finite_wall_force.png)

### Force, work and steady tow

These are model-unit resistance coefficients in units of $10^{-6}$ at angular frequency 0.000225, mouth radius 1 and nominal intake 0.25. Positive resistance means the actuator supplies power. For the solver’s $e^{+i\omega t}$ convention, the plotted push observable $q=\omega\operatorname{Im}(1/\chi_a)$ is the negative of the tabulated resistance $\operatorname{Re}(F/V)$. The last two columns measure the same signed work coefficient, $2\overline P/|V|^2$, on two interior windows.

| Box radius / depth | Old wall resistance | Corrected push | Brane + bulk work | Inner window | Wider window |
|---|---:|---:|---:|---:|---:|
| 9 / 10 | +542.11157 | −6.36253 | −6.36253 | −6.42743 | −6.42597 |
| 13.5 / 15 | +547.95068 | −0.57728 | −0.57728 | −0.58308 | −0.58160 |
| 18 / 20 | +548.61829 | +0.0917743 | +0.0917743 | +0.0895538 | +0.0910407 |

In the largest box the corrected harmonic resistance is $+9.17743\times10^{-8}$ and the steady tow gives $+9.17740\times10^{-8}$. The two remote measurements of harmonic momentum and energy use the same boundary field gradients. The two interior work measurements check the answer on different spatial supports. Their spread matters when subtracting nearly equal port fluxes. In the largest box, the corrected readings at the two sampled dynamical frequencies lie below their respective 1% force thresholds at separations 10 and 20. Applying those thresholds to an open exterior requires an exterior convergence test. [All 22 corrected rows, controls and source hashes](../../../research/pde_ledger_v3/_measurements/mouth_finite_wall_force_audit.json)

<details>
<summary>Signed push observable at all six frequencies on mesh 34</summary>

The entries below give $q$ in units of $10^{-6}$. The measurement record also includes the four independent resolution-control rows.

| Frequency | 9 / 10 | 13.5 / 15 | 18 / 20 |
|---:|---:|---:|---:|
| 0.000225 | +6.362529 | +0.577283 | −0.091774 |
| 0.000637 | +6.362542 | +0.577285 | −0.091776 |
| 0.002 | +6.362674 | +0.577303 | −0.091797 |
| 0.005 | +6.363445 | +0.577411 | −0.091917 |
| 0.01 | +6.366202 | +0.577798 | −0.092349 |
| 0.02 | +6.377240 | +0.579343 | −0.094094 |

</details>

### Where the error entered

The finite-element equations enforce zero relative normal current through the wall in their weak form. The old force integral treated the one-sided gradient as satisfying that condition pointwise. The flare exit makes a reentrant corner in the fluid cross-section. A virtual displacement equal to one within 0.4 mouth radii of the exit and tapered to zero by 0.7 accounts for 99.7% of the discrepancy on mesh 34 in the largest box. The error concentrates at that exit under quadrature and resolution checks. It is already present in the steady tow, which rules out an error requiring oscillatory damping. [Localization instrument and results](../../../research/pde_ledger_v3/_measurements/mouth_wall_force_localization.json)

The corrected solver evaluates the change in laboratory carrier momentum plus the outward momentum flux across both translating reservoir boundaries. It includes advection relative to the moving cut and the incident pressure in the separate stream experiment. As a check on the correction at angular frequency 0.01 in the 9 / 10 box, the revised native solve gives $a/F=0.08156169-0.0000042350i$ in the push and $a/U_s=0.02036563-0.0002925221i$ in the stream. An oscillating stream still produces a lag at zero intake. [Corrected force assembly](../tools/mouth_carrier_free_center.py) · [Force regression tests](../tests/test_mouth_rigid_force.py)

### Does drag vanish?

Across adjacent box sizes the measured brane work magnitudes give decay exponents 2.994 and 2.999. The receiving bulk gives 3.787 and 3.907. These trends are consistent with inverse-cube and inverse-fourth-power finite-box corrections. The net changes sign between radii 13.5 and 18. Three sizes do not uniquely determine the powers, crossover radius or infinite-domain intercept. The present result removes the apparent positive wall-resistance plateau as evidence against coherent drag cancellation. The remaining test is to hold this force definition fixed and vary the exterior matching radius at fixed frequency before taking the slow limit. Equilibrium and universal free fall of a material mouth remain separate questions. [Box-scaling calculation](../tools/mouth_drag_box_scaling.py) · [0PN roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

**Independent local review.** The momentum and energy reviewer reported “Scoped clear for the finite-box force correction, table and figure.” The source reviewer reported “Scoped clear for the corrected finite-box rigid-carrier force definition, its two regression tests, and the audit report.”
