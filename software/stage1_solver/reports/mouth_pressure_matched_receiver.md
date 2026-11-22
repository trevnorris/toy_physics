# Can the receiving bulk lose shear while it stays dense?

The free-fall calculation needs the sideways force delivered to the ordered material as constituents enter the bulk. A fixed crystal carried force and recorded no constituent passages in 99 finite runs. We therefore tested whether the same four-dimensional material can supply an ordered surface and a dense, rearranging receiver at a common external pressure. [Fixed-layer measurement](mouth_dense_ordered_interface.md) · [Free-fall roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

**The reference receiver.** At thermal energy 1.024 bond depths, one saved sample lost 84.31% of its original neighbors. Four independent directions selected before the hold still showed strong periodic density variations at the end. Neighbor replacement alone cannot identify shear loss. A separate held-strain experiment used the three archived preparations, two new noise histories each, unit constituent mass, Langevin damping rate 1, and zero imposed external pressure. All six undriven samples stayed connected. Over the last quarter of a 48.88171-native-unit hold, mean density was 1.68214 and mean pressure was −0.07043. The intervals below account for the three parent preparations. [Structure](../tools/mouth_receiver_structure.py) · [Pressure and shear measurements](../tools/mouth_receiver_pressure_readout.py)

| Native observable at 1.024 | Mean | Approximate 95% interval |
|---|---:|---:|
| Late/prompt shear, strain ±0.02 | 0.41591 | −0.49001 to 0.58899 |
| Late/prompt shear, strain ±0.01 | 0.43490 | −1.29958 to 0.75803 |
| Compression compliance, pressure steps ±0.25 | 0.023518 | 0.019155 to 0.027881 |

**Tagged material.** Across the six undriven paths the mean final root-mean-square displacement was 0.53021 bond lengths after removing box dilation and center-of-mass motion. Between 0 and 51 of the 128 sites moved at least one bond length, depending on the preparation and thermal path. On average 23.2070% of the original neighbor pairs were absent at the end. These displacements describe motion within a homogeneous sample. A measurement of material crossing an ordered interface requires that interface in the experiment. [Tagged-site reconstruction](../tools/mouth_receiver_tracer.py)

**Finding two preparations.** At 1.28 all sites stayed connected at the final four readings, but mean densities fell to 0.99980–1.10947. At 2.048 they fell to 0.10842–0.16611. Across eight temperatures, none met the prespecified density, connectivity, structural-order and motion criteria in all three histories. The criteria included a final-quarter density of at least 1.2 and a final tagged displacement of at least one bond length RMS. At **1.088**, one stayed connected at density 1.24236, lost 89.30% of its original neighbor pairs, and had weak periodic intensity at both endpoints (strongest peak 0.08134 → 0.09019). The other two had strong final peaks of 0.71516 and 0.76797. This isolated result motivated an independently prepared comparison at the same temperature. [Screen and thresholds](../tools/mouth_receiver_material_screen.py)

**Physical units.** The earlier two-sheet assignment of site mass, pitch and binding energy would put a 48.88171-unit mouth transit at **840.17332 native units**. The shorter holds cover 5.818% of that transit under this assignment. A neutral tracer measurement timed against the passing material could check the mapping. The thermostat supplies heat in these tests. The actual conversion energy and ambient bulk pressure remain open. [Conditional clock and measured sample energies](../tools/mouth_receiver_energy_clock.py)

**Same-temperature comparison.** All seven preparations received the same zero drive, opposite full and half strains (±0.02 and ±0.01 in the brane–fourth-direction plane), and opposite pressure steps (±0.25). The table gives the late signed shear response divided by the prompt signed response. The density-wave peak measures periodic order. Each row is one thermal history.

| Preparation at 1.088 | Final density-wave peak, undriven | Late/prompt shear, full | Late/prompt shear, half | Compression compliance |
|---|---:|---:|---:|---:|
| Sample 87, cooled A | 0.08873 | -6.85% | -26.95% | +0.12269 |
| Sample 87, cooled B | 0.09006 | +4.26% | -1.54% | +0.07593 |
| Sample 87, cooled C | 0.08499 | +3.46% | -11.36% | +0.13632 |
| Independent sample 41, cooled | 0.83559 | +78.77% | +81.42% | +0.03409 |
| Sample 87, ordered continuation A | 0.71661 | +46.50% | +39.20% | -0.00256 |
| Sample 87, disordered continuation B | 0.07181 | +5.44% | -25.98% | +0.10548 |
| Sample 87, ordered continuation C | 0.73768 | +60.27% | +56.42% | +0.02840 |

All 21 branches of the three cooled sample-87 preparations remained dense and weakly ordered at the measured endpoints. Across these three preparations the full-strain late/prompt ratio was **0.00115** (95% interval −0.14128 to +0.17669), the half-strain ratio was **−0.13687** (−0.42352 to +0.22550), and compression compliance was **0.11165** (+0.03296 to +0.19034). Only the full-strain interval fit within the prespecified ±0.25 screening margin. Negative single-history estimates, including the small negative compliance on ordered continuation A, remain unresolved individually.

Independent sample 41 already had a strong density-wave peak of 0.74513 before cooling, versus 0.05867–0.07515 for the three disordered precursors. Its large retained shear showed why bath temperature alone could not select the receiver. A subsequent test prepared that same independent geometry in the disordered state. [Independent follow-up](mouth_independent_receiver.md) · [Signed comparison](../tools/mouth_receiver_preparation_readout.py) · [Measurements and independent checks](../../../research/pde_ledger_v3/_measurements/mouth_pressure_matched_receiver.json)

Working baseline `a91a633f`. Formal repository review pending.
