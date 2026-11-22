# Can a warm 4D bulk lose shear while retaining compression?

**At the highest temperature the bulk rearranged while remaining connected and resisting compression. Its sideways restoring force over one transit is still unresolved.** Sample 41 had already shown that internal motion can lower the measured sideways force while leaving the material in the same local energy minimum. This experiment asks whether warming the less ordered sample 87 can instead produce sustained neighbor exchange in dense material within one downstream transit. [Previous warm and cold measurements](mouth_bulk_isotropic_4d.md)

We kept the same 128 constituents, four periodic spatial dimensions, direction-independent Morse potential and zero external pressure used to select the cold sample. Each of three independent thermal histories starts from those coordinates. Each history warms in steps of 0.008, 0.064, 0.256, 0.512 and 1.024 in units of the bond well depth, spending 48.8817124312 model time units preparing at each step. At each temperature matched copies branch off for a separate hold of the same duration. The preparation itself supplies the starting point for the next temperature. The noise is common to all signed drives within one history.

- **Shear.** Apply opposite deformations along the easiest cold shear direction. Hold the deformed shape while allowing the volume to respond at zero pressure. Repeat at half amplitude and compare the remaining stress in the final quarter of the hold.
- **Compression.** Apply small opposite pressure steps about zero while holding shape fixed. Measure the signed change in log volume. A resolved positive compliance gives a finite-time bulk modulus. Compare it with the thermal pressure scale of a dilute gas at the measured density.
- **Material exchange.** Track density, connected cluster size, persistent and replaced neighbors, displacements after subtracting cell expansion and translation, and potential energy per constituent. The sweep measures finite-time material states at zero applied pressure. Expansion or dispersal would have to be distinguished from loss of shear in dense material.

All 75 planned trajectories passed their guards. The cold force probe gave shear stiffness 23.70723 and compression compliance 0.007764814. Independent energy calculations gave 23.70716 and 0.007764813. The warm test used the same cold-selected shear direction and a volume controller with response time 0.5. The table gives mean signed slopes over the last quarter of the hold, averaged over the three histories. The last column reports the fraction of starting neighbors absent at the end of the undriven hold. [Experiment](../tools/mouth_bulk_temperature_steps.py) · [Readout](../tools/mouth_bulk_temperature_readout.py) · [Measurements and independent check](../../../research/pde_ledger_v3/_measurements/mouth_bulk_temperature_steps.json)

| Thermal energy / bond depth | Late shear at full strain | Late shear at half strain | Original neighbors absent |
|---:|---:|---:|---:|
| 0.008 | 24.17 | 24.11 | 0% |
| 0.064 | 40.71 | 40.58 | 0–0.42% |
| 0.256 | 62.24 | 62.29 | 0–4.63% |
| 0.512 | 62.60 | 643.80 | 0.96–21.49% |
| 1.024 | 17.51 | 23.93 | 40.14–84.31% |

At 0.008 all three undriven histories retained their original neighbors at every saved frame. At 0.256 individual late restoring slopes were positive, with sizable differences between histories. At 0.512 one history changed from 64.49 to 1826.30 on halving the shear strain. Its measured slope depends strongly on the probe amplitude. The three-history compression compliance interval at this temperature also crosses zero.

At 1.024 all 128 sites belonged to one cluster at every saved undriven frame. Final-quarter densities were 1.611–1.707 in model units. Each history compressed under the finite pressure difference. Inverting the mean compliance gives an apparent bulk modulus of 38.34, with an approximate three-history 95% interval of 20.45–306.85. The ideal thermal pressure scale at the measured mean density was 1.71. The mean internal pressure differences between the two driven branches ranged from 0.320 to 0.628 for an imposed difference of 0.500. The three signed shear slopes were +65.07, −49.80 and +37.25 at full strain and +94.97, −26.31 and +3.14 at half strain. The full-strain mean’s approximate 95% interval runs from −131.35 to +166.36. The independent audit reproduced the signed reductions and checked pressure and shear directly from all 15 saved undriven endpoints to within 1.1 × 10⁻⁸.

**Consequence for 0PN.** Heating can produce extensive neighbor exchange in a connected sample that supports compression. The measured shear still permits both substantial restoring stress and relaxation near zero over one transit. The condition for advancing to the mouth test remains open. A further neutral-material test would need a more precise upper bound on late shear across independent histories, with agreement at both probe amplitudes, before testing time step, pressure control and additional shear directions. The sweep imposed the temperature. The energy available during actual conversion and the free-mouth acceleration remain open.

Working checkpoint `a91a633f`. Independent numerical review complete within the scope above. Formal repository review is pending.
