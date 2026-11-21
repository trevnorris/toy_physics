# Centered and off-center attachments at the same shear order

**At half shear order, the off-center attachment exerted 25% less sideways force in the central binding well. A second contact law with identical force curves along each axis gave equal forces.** Both measurements use the upstream site of the same held 0PN mouth, radius 1 and intake 0.25, at half the reference conversion drive. “Off center” describes the alignment of two local material features within their binding region. [Earlier single-attachment test](mouth_surviving_attachment.md)

**Matched displacements.** At each depth, two copies started at zero sideways extension. One feature was centered and the other held one eighth of a material pitch off center in the perpendicular direction. We gave them the same sideways displacement of \(10^{-5}\) pitch, then reversed and halved the step on fresh copies. The table gives restoring force magnitudes in units of \(10^{-6}\) model force per unit contact area. Both attachments stayed connected.

| Local shear order | Centered, either law | Off center, central well | Off center, independent directions |
|---|---:|---:|---:|
| 0.75 | 2.5000 | 2.2222 | 2.5000 |
| 0.50 | 2.5000 | 1.8750 | 2.5000 |
| 0.30 | 2.5000 | 0.7639 | 2.5000 |

In the central well the binding strength depends on distance from the binding center in the plane of the two displacements. The off-center attachment starts closer to the edge. Its small-step force divided by the centered force approaches \(1-(2\eta/\chi)^2\), where \(\eta\) is the measured offset in pitches. The comparison law applies the same one-dimensional binding force independently in each direction. Reversing the step reversed the sideways force in both laws.

**What happens when the feature can move.** We followed two more copies through conversion, starting from the same equilibrium with an elastic anchor whose unloaded position was one eighth of a pitch off center. One kept the unit-stiffness anchor and the other was released. At order 0.30 their actual offsets were 0.10584 and 0.02156 pitch. Their instantaneous sideways forces in the central well were 50.21% and 97.93% of the centered reading. The drag and anchoring are trial material parameters. [Calculation](../tools/mouth_matched_attachment.py)

**An independent material input.** Holding the off-center feature at half order took slightly less force along its original offset during either sign of sideways displacement. The decrease was \(5.00\times10^{-11}\) model force per unit area in the central well and zero in the independent-directions control. Integrating this change along the offset reproduced the difference in sideways holding work. An independent material law or force reading with both displacements present would test this coupling. The present brane fields supply the local order and transit, while the attachment force and its contribution to 0PN motion remain conditional.

Working base `49dca23e`. Two independent local reviewers checked the recorded scopes. Formal review and commit are pending. [Readings, commands and reviews](../../../research/pde_ledger_v3/_measurements/mouth_matched_attachment.json)
