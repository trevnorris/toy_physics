# Step 1. Does a slowly moving mouth lose energy to drag?

**Archived table.** The later finite-box audit identified a pressure-readout error at the flare exit and replaced the push resistance in all 22 rows. The signed reservoir and interior work measurements below still apply. [Corrected forces, stream and full diagnostic plot](mouth_finite_wall_force_audit.md)

**Result.** The wall force inferred from a weak push has a low-frequency plateau at all three reservoir sizes. The cycle-averaged energy carried through the brane and bulk accounts for a much smaller force, with opposite signs in the smallest and largest boxes. The wall plateau therefore fails the independent work check. The steady, unbounded-mouth drag remains unresolved. This measures a coherent carrier and a prescribed rigid mouth. The supersolid brane, its shear waves and the self-supported flare are physical requirements for the full model. [Corrected roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

### Push in resting reservoirs

All runs use the same radius-one cubic mouth, nominal stationary intake 0.25, and fixed far phases on the translating finite chart. Both reservoirs retain their steady throughflow. The push imposes zero incident sideways stream in the laboratory. The main series keeps the original mesh-34 nodes while moving the outer boundaries. Mesh labels 26 and 48 are independent resolution controls at the original box size. The observable is the response to external force on the entrance, $\chi_a=a/F$. Peak amplitudes and the $e^{+i\omega t}$ convention give

$$
q=\omega\,\operatorname{Im}(1/\chi_a)=-\operatorname{Re}(F/V),\qquad
\gamma_E=2(\overline P_{\rm brane}+\overline P_{\rm bulk})/|V|^2 .
$$

If the wall force and independently measured energy transfer balance on the same harmonic drive, $-q$ equals $\gamma_E$. Positive $\gamma_E$ means net excess power outward. A negative value means the reservoirs together supply power in that preparation. The sign by itself does not identify microscopic heat.

**Every force coefficient below is in units of $10^{-6}$.** “Brane” and “bulk” denote the two signed outward power coefficients. Their sum is $\gamma_E$. The last two columns give the corresponding total from two spatial integrals around the mouth. Box dimensions are radial extent / depth in mouth-radius units. [Frequency and box-size plot](figures/mouth_rigid_drag_table.png) · [All complex inputs, guard outcomes and computed rows](../../../research/pde_ledger_v3/_measurements/mouth_rigid_drag_table.json)

| Mesh | Box | $\omega$ | Push $q$ | Brane | Bulk | $\gamma_E$ | Inner window | Wider window |
|---:|:---:|---:|---:|---:|---:|---:|---:|---:|
| 26 | 9 / 10 | 0.002 | -501.75407 | +9.66456 | -15.99114 | -6.32657 | -6.43965 | -6.43551 |
| 26 | 9 / 10 | 0.01 | -501.75139 | +9.67188 | -16.00197 | -6.33009 | -6.44321 | -6.43908 |
| 34 | 9 / 10 | 0.000225 | -542.11157 | +9.69821 | -16.06074 | -6.36253 | -6.42743 | -6.42597 |
| 34 | 9 / 10 | 0.000637 | -542.11156 | +9.69824 | -16.06078 | -6.36254 | -6.42744 | -6.42599 |
| 34 | 9 / 10 | 0.002 | -542.11148 | +9.69851 | -16.06118 | -6.36267 | -6.42757 | -6.42612 |
| 34 | 9 / 10 | 0.005 | -542.11100 | +9.70011 | -16.06356 | -6.36345 | -6.42835 | -6.42690 |
| 34 | 9 / 10 | 0.01 | -542.10931 | +9.70583 | -16.07203 | -6.36620 | -6.43113 | -6.42967 |
| 34 | 9 / 10 | 0.02 | -542.10253 | +9.72876 | -16.10600 | -6.37724 | -6.44225 | -6.44079 |
| 34 | 13.5 / 15 | 0.000225 | -547.95068 | +2.88080 | -3.45808 | -0.57728 | -0.58308 | -0.58160 |
| 34 | 13.5 / 15 | 0.000637 | -547.95068 | +2.88082 | -3.45810 | -0.57728 | -0.58308 | -0.58160 |
| 34 | 13.5 / 15 | 0.002 | -547.95071 | +2.88100 | -3.45831 | -0.57730 | -0.58310 | -0.58162 |
| 34 | 13.5 / 15 | 0.005 | -547.95089 | +2.88209 | -3.45950 | -0.57741 | -0.58321 | -0.58173 |
| 34 | 13.5 / 15 | 0.01 | -547.95154 | +2.88598 | -3.46378 | -0.57780 | -0.58360 | -0.58212 |
| 34 | 13.5 / 15 | 0.02 | -547.95413 | +2.90161 | -3.48095 | -0.57934 | -0.58516 | -0.58368 |
| 34 | 18 / 20 | 0.000225 | -548.61829 | +1.21576 | -1.12399 | +0.09177 | +0.08955 | +0.09104 |
| 34 | 18 / 20 | 0.000637 | -548.61829 | +1.21577 | -1.12400 | +0.09178 | +0.08956 | +0.09104 |
| 34 | 18 / 20 | 0.002 | -548.61836 | +1.21591 | -1.12412 | +0.09180 | +0.08958 | +0.09106 |
| 34 | 18 / 20 | 0.005 | -548.61877 | +1.21674 | -1.12482 | +0.09192 | +0.08970 | +0.09118 |
| 34 | 18 / 20 | 0.01 | -548.62022 | +1.21968 | -1.12733 | +0.09235 | +0.09013 | +0.09161 |
| 34 | 18 / 20 | 0.02 | -548.62605 | +1.23152 | -1.13743 | +0.09409 | +0.09187 | +0.09336 |
| 48 | 9 / 10 | 0.002 | -521.25020 | +9.72219 | -16.10867 | -6.38648 | -6.41846 | -6.41810 |
| 48 | 9 / 10 | 0.01 | -521.24826 | +9.72952 | -16.11954 | -6.39002 | -6.42201 | -6.42165 |

### What the work check resolves

At $\omega=0.01$ in the 9 / 10 box on mesh 34 the direct wall reaction implies mean actuator power $+0.000271055\,|V|^2$. The two interior energy measurements give $-0.00000321556\,|V|^2$ and $-0.00000321484\,|V|^2$. Changing mesh from 26 to 48 in that box preserves the sign of the net reservoir power. Extending the box at fixed inner resolution reduces that power and then changes its sign. In the 18 / 20 box the independently computed **steady** force coefficient from the full momentum balance is $+9.17740\times10^{-8}$, close to the $+9.17743\times10^{-8}$ reservoir coefficient at $\omega=0.000225$. This establishes a slow limit in that finite preparation. A bound on open-domain drag requires a transparent or independently verified exterior and an independently balanced mouth force.

For the roadmap’s 1% force-tolerance pilot, the reference source held at rest gives allowable resistances of $7.80906\times10^{-5}$ at 10 mouth radii and $2.76092\times10^{-5}$ at 20. Both largest-box reservoir readings are smaller than the respective threshold. The wall-versus-reservoir discrepancy is $5.48527\times10^{-4}$ at both sampled frequencies. These data cannot certify the pilot range. [Sensitivity inputs and arithmetic](../../../research/pde_ledger_v3/_measurements/mouth_rigid_drag_checks.json)

The energy current for the carrier Hamiltonian is $\mathbf J_E=-\operatorname{Re}[(\partial_t\psi)^*\nabla\psi]$ in the laboratory. A cut moving with the frame carries relative flux $\mathbf J_E-e\mathbf V_{\rm cut}$. Both reservoir boundaries move in the experiment. Interior and outer measurements track flux independently of differentiating the field on the sharp wall. Cycle-mean storage vanishes on the bounded periodic state. The unknown second-order field enters this energy current through the time derivative of a periodic correction multiplied by the stationary gradient, whose cycle mean vanishes. The separately solved mean particle throughput does change. Restoring a common chemical energy per particle shifts the two individual reservoir readings by opposite amounts. Their sum is unchanged to the measured particle-balance precision. The signed and half-strength powers follow perturbation theory and were also checked on reconstructed fields. Four fully nonlinear trajectories would provide a separate amplitude test. [Energy calculation](../tools/mouth_rigid_push_energy.py) · [Mean throughput and amplitude controls](mouth_push_mean_throughput.md)

### Stream past a mouth

In a separate drive the far brane stream oscillates and the far bulk rests. At mesh 34 and $\omega=0.01$, $a/U_s=0.02040828-0.00020124i$. Dividing that response by $a/F$ on the same preparation reconstructs the total stream-force coefficient; its quadrature divided by frequency was −0.35737 at 0.01 and −0.35821 at 0.02. At zero intake the mouth still has a stream lag while its push resistance is zero on the tested mesh. An unsteady incident-pressure force is consistent with these data. The full Froude–Krylov contribution would require integrating the undisturbed incident pressure on a specified surface. The stream’s 1% lag supplies no estimate of steady drag. [Separate push and stream measurements](mouth_drag_response_audit.md)

**Independent check.** A local reviewer recovered the signed susceptibility from the applied forces and complex velocities in all 22 saved runs. An independent momentum calculation compared the direct wall, two interior windows and both reservoirs, including the largest-box steady tow. The review verdict was “Scoped clear for the 22 signed harmonic rows, both nested-window power estimates, and the independent matched-bath DC momentum comparison.” [Independent steady and harmonic readings](../../../research/pde_ledger_v3/_measurements/mouth_rigid_drag_checks.json)

**Next decision.** Test a momentum-conserving force on the moving mouth with a coherent exterior that transmits outgoing waves. Verify the exterior by moving the matching surface at fixed frequency, then repeat the slow limit and steady tow. This distinguishes wave radiation or physical momentum export from work by a distant, moving reservoir. The equilibrium test in Step 2 remains a separate prerequisite for interpreting any measured impedance as that of a physical mouth.
