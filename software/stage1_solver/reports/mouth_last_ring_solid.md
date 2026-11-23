# Bonded crystal at the downstream rim of the radius-four mouth

With the archived fields and displacement footprint held fixed, the bonded crystal 97.7569% down the rejected wall has **positive radial–axial stiffness and an outward residual force** when integrated across the registered material triangles. A native four-point-per-axis element rule reports an inward radial force of −0.00186587 on the same material. This sign difference matters for deciding which correction to make before releasing the wall.

The test holds the wall at downstream depth 0.79475445 and exit radius 3.50085165. The original last free material node moves by a signed radial or axial distance, tapering linearly to zero over its six surrounding triangles. All refinements preserve that displacement pattern and the saved strain field. Splitting the integration where the reference crystal changes triangle resolves gradients that the original integration sites miss. The independent fluid-grid comparisons hold the crystal and its registration fixed.

| Saved-field integration | Radial force | Axial force | Lower stiffness | Higher stiffness |
|---|---:|---:|---:|---:|
| Original rule | −0.001866 | +0.000106 | 551.75 | 2167.22 |
| 32 uniform subdivisions per edge | +0.195968 | +0.010884 | 558.81 | 2180.96 |
| 64 uniform subdivisions per edge | +0.190588 | +0.010644 | 557.03 | 2179.84 |
| 32 depth subdivisions aligned to material and order boundaries | +0.191377 | +0.010604 | 555.38 | 2181.97 |
| 64 aligned depth subdivisions | +0.191375 | +0.010680 | 555.05 | 2180.55 |

Stiffnesses are the two eigenvalues of the finite-step energy curvature under radial and axial motion, including their cross term. The interpolated order has slope jumps, so the exact infinitesimal tangent can differ slightly. The original row uses signed steps of 0.0001, the remaining rows 0.00005. At 64 aligned subdivisions, the mixed stiffness is +629.22. Doubling the step gives eigenvalues 555.00 and 2180.13. At the main ring the minimum stretch is 0.323664. For radial steps of ±0.0001 it becomes 0.323677 and 0.323651 respectively. For axial steps it becomes 0.323809 and 0.323515. One ring farther into the solid the minimum on the unperturbed footprint is 0.301901. All sampled volumes stay positive. The main footprint reaches order 0.255 where it enters the conversion front. The two adjacent hats also have positive two-direction stiffness, with radial forces +0.302782 upstream and −0.227638 inward of the target.

At the target ring the original radial elastic stress and spatial-order terms are +0.818957 and −0.820823. The integration aligned to both material charts gives +1.042197 and −0.850822 for the same fields. Most of the original radial sign error comes from integrating the elastic stress.

Independently rematching the carrier and terminal order on the same saved wall changes the solid force below. The integral uses 64 aligned depth subdivisions and signed steps of 0.0001 for all three readouts.

| Carrier and conversion mesh | Radial force | Axial force | Lower stiffness | Higher stiffness | Minimum order |
|---|---:|---:|---:|---:|---:|
| Saved native | +0.191375 | +0.010680 | 555.00 | 2180.13 | 0.254748 |
| Globally doubled | +0.115917 | +0.147799 | 553.08 | 2183.64 | 0.262719 |
| Further refined at rim | +0.140119 | +0.181955 | 554.48 | 2182.58 | 0.263584 |

The shrinking-footprint control samples progressively less material. Its radial force falls from +0.35148 at half the original width to +0.04387 at one sixteenth. The finite fixed-footprint force above is the appropriate comparison across resolutions. Its scope is the saved solid and its terminal order loading. Establishing a stationary sharp corner additionally requires releasing interior displacements and resolving the total closed-surface fluid force. [Calculation](../tools/mouth_last_ring_solid.py) · [Measurements](../../../research/pde_ledger_v3/_measurements/mouth_last_ring_solid.json)
