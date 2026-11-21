# Does the neutral patch reproduce the brane’s pressure law?

**One scale matches the resting compression modulus. The same 24-site patch stiffens too rapidly to match the carrier’s pressure changes at other densities.** This tests a shared-material hypothesis away from a mouth, at zero intake. We used the relaxed central-bond preparation from the attachment test. The carrier’s homogeneous energy gives the target pressure increment at conserved constituent number: \(p(n)-p(1)=nW'(n)-W(n)=(n^5-1)/5\), in solver units. [Carrier energy](mouth_carrier_0pn_methods.md) · [Neutral preparation](neutral_patch_sites.md)

**Calibration.** We assigned one pitch of depth per layer, one pitch to the unresolved direction and the reference triangular area per site, with pitch 0.01 mouth radii. At rest density one, the cell represents 2.07846 × 10⁻⁷ carrier number in the same numerical four-volume. This association is a hypothesis. We strained the three brane directions equally and kept the fourth-direction thickness fixed. In the primary case the layers were separated along the fourth direction, so their gap stayed fixed. In the alternative orientation their gap changed with the brane strain. Fifteen edge sites followed the imposed strain, while nine interior sites could rearrange in their own planes. We also tested moving every site with the strain. The conjugate pressure is the mean stress for this loading path. All comparisons subtract the rest load, including any gap preload.

Only a **1% density increase** set the multiplicative scale. For the relaxed fixed-gap cell, the fit assigned one patch energy unit a value of 4.46280 × 10⁻¹⁰ in carrier energy units, conditional on the chosen cell volume. Independent energy differences verified the boundary forces.

| Change in density | Carrier pressure change | Relaxed fixed-gap prediction | Relative discrepancy in magnitude |
|---|---:|---:|---:|
| −2% | −0.01921584 | −0.01839252 | −4.28% |
| +0.5% | +0.00505025 | +0.00501500 | −0.70% |
| +1%, calibration | +0.01020201 | +0.01020201 | fitted |
| +2% | +0.02081616 | +0.02110262 | +1.38% |
| +5% | +0.05525631 | +0.05820886 | +5.34% |

Let \(\eta=\log n\). The ratio \(p''(0)/p'(0)\), with derivatives taken in \(\eta\), is independent of an overall force scale. At log step 0.0005 it was 7.87112 with the relaxed fixed gap and 8.47240 with the relaxed brane-direction gap. Both affine controls also exceeded the carrier value of 5. This rules out exact matching with one multiplier in all four tested finite-cell preparations. [Calculation](../tools/neutral_patch_carrier_compression.py) · [Measurements and independent review](../../../research/pde_ledger_v3/_measurements/neutral_patch_carrier_compression.json)

**Meaning for 0PN.** At intake 0.25 the archived deep wall from depth 1.175 to 2 sampled densities 0.995856–0.999050 on mesh 34. At nine evenly spaced densities in that interval the 1% calibration underestimated the magnitude of the homogeneous pressure change by 1.54–1.99% for the relaxed fixed-gap cell. A separate calibration to the infinitesimal resting modulus reduced that to 0.14–0.60%. Mesh 26 gave essentially the same interval. These percentages pertain to homogeneous compression. The total wall reaction also depends on shear, motion, coverage and where the other layer transfers momentum. The force per actual wall area and the physical free acceleration remain open. [Operating-density check](../tools/neutral_patch_compression_range.py)

**Decision.** Before changing the bond law, test a neutral sample with periodic boundaries in the represented in-plane directions, preserving the bond law and testing size dependence. This would distinguish a finite-edge contribution from a persistent error in the compression curve. Even agreement would require an independently justified amount of material per wall area before a balanced mouth push could predict free fall. Working base `1dc45d09`. An independent reviewer reproduced the compression results with a separate energy calculation and checked the operating-density comparison. The measurement record contains both scoped verdicts. Formal review and commit are pending.
