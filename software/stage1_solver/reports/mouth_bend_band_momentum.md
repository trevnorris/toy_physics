# Sideways momentum at the entrance and middle of a draining bend

**The draining fluid carries more sideways momentum at mid-depth. The wall material’s own momentum remains unmeasured.** At mass intake 0.0373999726942 and peak displacement \(10^{-5}\), we integrated the component of the fluid’s laboratory-frame momentum in phase with the peak bend velocity \(\dot b\). The entrance covers depths \(0\) to \(H/4\), the middle \(3H/8\) to \(5H/8\), on the same draining arm. The figures below give momentum divided by peak bend speed in model units, on the finer grid. [Carrier calculation](../tools/mouth_bend_carrier_band_current.py) · [Guarded outputs on both grids](../../../research/pde_ledger_v3/_measurements/mouth_bend_band_momentum.json)

| Angular frequency | Entrance, full section | Middle, full section | Entrance, outermost \(R/8\) | Middle, outermost \(R/8\) |
|---:|---:|---:|---:|---:|
| 0.01 | 1.07474 | 1.92013 | 0.38397 | 0.94918 |
| 0.02 | 1.08758 | 1.93391 | 0.38710 | 0.95269 |

At 0.02 the coarser grid gave 1.10178 and 1.94839 across the full sections. Reducing the sampled thickness from \(R/8\) to \(R/32\) on the finer grid reduced the two near-wall readings to 0.10761 and 0.32307. The thinnest shell spans two radial cells on that grid. These volumes contain flowing material. The solver prescribes the wall motion, and its material action takes the local wall inertia as an input. There is no separate wall momentum to read from this preparation. The spectral calculation already includes the fluid reaction. [Reservoir field equations](../tools/mouth_open_reservoirs.py#L158) · [Wall action and effective inertia](../../../research/pde/paper/pde.tex#L775)

**Would two wall readings decide the peak count?** We tested two positive distributions of normally moving wall material at the same total radial mass 59.3830741481. One follows the supporting wave’s local pressure. The other puts some material just beyond the measured middle band, where the wall still moves appreciably, and matches both band integrals. These rows are conditional predictions. [Wall calculation](../tools/mouth_bend_band_momentum.py)

| Assumed wall mass profile | Entrance \(P_x/\dot b\), one arm | Middle \(P_x/\dot b\), one arm | Bend inertia, both halves | Peaks |
|---|---:|---:|---:|---:|
| Proportional to wave pressure | 0.626832 | 2.497242 | 7.802142 | 1 |
| Redistributed at the same band readings | 0.626832 | 2.497242 | 13.771688 | 2 |

Those two readings and the radial normalization allow limiting bend inertias from 4.37578 to 14.57709. The fine-grid transition on the tested branch is **13.00479**, at sidewall stiffness 0.424660352295 and the previously measured radiation loss. Both spectral grids give the peak counts in the table for the equal-area readout with 1:3 widths, a common shift and zero extra dephasing. The pressure row integrates the bend shape through each pressure cell. The earlier midpoint approximation gave 7.802236. An independent reviewer verified the bounds and constructed other matching profiles of mass 12 and 14 with one and two peaks. [Independent calculation](../tools/mouth_bend_band_momentum_check.py) · [Spectrum and its assumptions](mouth_bend_coefficients.md)

**The deciding measurement.** Set \(u=w/H\) and \(f(u)=\sin^2(\pi u)\). Under the assumed local normal inertia law, the conjugate bend mass is

\[
M_{\rm wall}=\sum_{\rm halves}\int_0^1 f(u)\,\frac{dP_{x,\rm wall}(u)}{\dot b}.
\]

Here \(dP_{x,\rm wall}\) is the sideways projection of the **normal wall-material momentum** in a depth slice, in phase with \(\dot b\). One full-depth integral suffices if the wall halves respond identically, in which case its transition value is **6.50240** per half. A material equation or an independent wall-sensitive measurement must supply that response. An independent total holding force can also determine the inertia and stiffness after subtracting the modeled field reactions. This fixed-end test leaves the freely translating mouth at 0PN open. Working result at base `49dca23e`, independently reviewed at the scopes recorded in the manifest. [Held-force alternative](mouth_bend_holding_force.md)
