# Closed axial momentum balance around the radius-four rim

**Enclosing the whole rim explains most of the sensitivity to measurement thickness. It leaves the absolute axial force unresolved.** We measured draining minus resting momentum on the same wound, angularly relaxed solid. Each volume encloses both sides of the full disk and the surrounding receiver. Its two end caps and outer cylindrical side lie outside the mobile conversion layer. Increasing its depth also encloses more of the held, sloping duct wall. [Commands and measured force components](../../../research/pde_ledger_v3/_measurements/mouth_closed_rim_work.json)

| Carrier and front grid, 33 × 32 solid | Axial shell | Outward flux | Sum of both sealed-face tractions | Residual |
|---|---|---:|---:|---:|
| 13 × 10 × 14 | Near | −0.510153 | +0.502440 | −0.007713 |
| 13 × 10 × 14 | Far | −0.463688 | +0.457604 | −0.006084 |
| 26 × 20 × 28 | Near | −0.499924 | +0.447171 | −0.052753 |
| 26 × 20 × 28 | Far | −0.454026 | +0.401387 | −0.052639 |

The integration weights equal one across the outlet. The near axial taper spans 0.172–0.224 below it and 0.158–0.277 above. The far taper spans 0.224–0.276 below and 0.277–0.486 above. The table uses the radial taper from 4.300 to 4.839 at disk height. Moving that taper out to 4.839–5.593 gave refined residuals of −0.052837 and −0.052731. The same transfer to the finer **solid** mesh at fixed carrier resolution had changed the residuals by less than 0.000056. [Independent sharp-surface reconstruction](../../../research/pde_ledger_v3/_measurements/mouth_closed_rim_review.json)

The doubled **carrier** grid subdivided exactly the same physical integration regions. The solid was held while the fields rematched. Current changed from 11.61684 to 11.53320. All four refined weak mass balances closed within 1.21 × 10⁻⁹. The thickness dependence at the inner radial taper dropped from 0.001629 to 0.000114, while the common force level shifted by 0.045–0.047. One-sided gradients used to recover wall traction retain a spurious normal-current trace on the weakly impermeable wall.

A common translation convected the crystal and all order pins. Its signed energy derivative, including both sealed-wall vortex terms, verified the native weak flux within 4 × 10⁻⁷. The drive-induced crystal stress cancels on the outer measurement surface. Assigning the incident or receiving axial slope to the shared order pin gives different Eulerian constraint works, −0.01732 to −0.07894 on the nearer refined shell. Those are alternative assignments. Adding either to the common translation would change the imposed variation. The residual cannot establish a physical force discontinuity or supply a free-mouth force. The next discriminator is a conservative wall-traction reconstruction tied to a specified mechanical motion of the order–crystal junction.
