# Does a companion produce a different mouth stress from a matched stream at 0PN?

**Same-reservoir comparison.** The target has radius 1 and prescribed inward number rate $J=0.25$. The companion has the same preparation. Reuse the guarded pair, target-alone and companion-alone fields at separations 8 and 16 in the cylindrical reservoir of extent 48, with 4 lateral cells per radius and 16 depth intervals. For each separation solve two additional target-alone controls with the density fixed at the same distant boundary and a uniform incident carrier phase gradient. One matches the source-only axial speed interpolated at the geometric center. The other matches the source-only speed weighted over the opening by the isolated target intake. The uniform bare stream is a plane wave; its speed was checked with the force readout's spatial derivative. All four new target runs verified the field equation and the intake. [Solver](../tools/companion_uniform_stress.py) · [Inputs, guarded outcomes and output](../../../docs/_measurements/companion_uniform_stress.json)

**Measured difference.** Positive is toward the companion. The excess below is the pair minus its matched stream, with the same isolated target baseline in both. The combined column includes carried sideways momentum and mixed density stress through the full opening. Its ppm unit is relative to $JU$ for the indicated incident-speed convention. The stress column is the mixed-stress contribution alone, in solver force units.

| Separation | Incident speed match | Combined excess (ppm) | Stress excess |
| ---: | --- | ---: | ---: |
| 8 | Center | 147.192371 | $-5.98593\times10^{-12}$ |
| 16 | Center | 38.5382443 | $-1.07370\times10^{-13}$ |
| 8 | Intake-weighted | 7.92816932 | $-6.37585\times10^{-12}$ |
| 16 | Intake-weighted | 0.984187992 | $-1.33024\times10^{-13}$ |

For the center match, doubling separation multiplies the paired flux by $0.243959596$, the combined excess by $0.0638810012$, and the stress excess by $0.0179370108$. Reconstructing the combined flux from the interior layers gives center-matched excesses of 147.930385 and 38.6947192 ppm. These are two distances on one mesh. They establish no continuum exponent or physical force bound. [Comparison outputs](../../../docs/_measurements/companion_uniform_stress.json)

**Ideal-fluid check.** On a sphere in steady, constant-density, inviscid potential flow, superpose a radial target and an undisturbed exterior point sink. At matched center speed $U$, pressure and carriage each give the same target–incident cross term as for a uniform stream: their *inward* momentum contributions are $-\rho Q_tU/3$ and $+4\rho Q_tU/3$. The companion-only Bernoulli pressure varies across the sphere and its integrated axial difference starts at $d^{-5}$; its companion-only carriage cancels it in the closed-surface total. The analytic sphere and the finite-slab aperture are distinct cuts. [Mathematica derivation and SymPy cross-term check](../../../docs/_measurements/companion_uniform_stress.json)

**Physical status.** The matched flow accounts for the leading observed carrier flux in this test. An independent pull through the brane would require a companion-dependent **total** traction at matched material slip. The carrier archive supplies neither that slip nor the reaction holding the geometric mouth still. Additional material stresses, outgoing momentum and the free acceleration remain open. Working on base fd667ea9, uncommitted pending review. [Rim test and missing readings](one_mouth_rim_slip_comparison.md)
