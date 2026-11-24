# Radius-four checkpoint 10: a fold at the edge of ordered material

The 24 × 24 unpinned relaxation is paused at checkpoint 10. The held entrance radius is 4, depth is 0.81914036 and exit radius is 3.30637727. The inward current is 10.70857269, of which 66.7562% at the incident disk cut is dense and subsonic including swirl. The fully melted shaft extends from about physical depth 0.2846 to the exit. The test uses one quantum of circulation, phase drop 0.4510084282857874 and the original fluid and order energies. [Measurements and reproducing command](../../../research/pde_ledger_v3/_measurements/mouth_r4_checkpoint10_boundary.json) · [Checkpoint diagnostic](../tools/mouth_r4_checkpoint10_boundary.py)

## What the shape is trying to do

Each value is the negative centered difference of the energy with respect to one solver coordinate, at fixed material offsets and transported fields. Positive and negative give the preferred direction of an isolated displacement. The exit coordinate is four times the change in log exit radius. The width modes change the profile at fixed endpoints.

| Generalized force | Checkpoint 8 | Checkpoint 9 | Checkpoint 10 |
|---|---:|---:|---:|
| Log depth | +0.162504 | +0.151886 | +0.150565 |
| Scaled log exit radius | −0.656703 | −0.642728 | −0.641236 |
| Middle width | +0.609296 | +0.695823 | +0.705499 |
| Entrance-to-exit width asymmetry | −0.590523 | −0.659160 | −0.666834 |
| Quartic middle-width mode | +0.442165 | +0.463347 | +0.465603 |
| Polar–equatorial distortion | +0.009859 | +0.009103 | +0.009007 |

All six forces retain their signs. The inward force per unit exit radius weakens from −0.790568 to −0.775757. At checkpoint 9 the drive to widen the middle becomes the largest component. The increasingly negative asymmetry component favors widening the entrance-side half and narrowing the exit-side half. The exit narrows over these three states. They establish neither runaway closure nor a self-supporting equilibrium.

## Why the permitted remap cannot repair the fold

The backtracking log identifies spatial order node 51603 at radius 3.63538, depth 0.181112 and polar angle 0.121295 radians. The production check encloses eight folded material prisms on checkpoint 10. Only the enclosure for prism **(depth 5, radial 17, polar 0), triangle 0** includes this node. A tighter full-prism test using the actual spatial basis of 51603 proves that its positive support is disjoint from every nonpositive determinant on the saved state. After the full shape trial the node had order 0.0000364030 and a one-sided lower-bound derivative of −7.20 × 10⁻⁹. That trial failed the conservative refreezing check. Its full field was not archived.

The implicated prism has a real fold at solid vertex **(5, 18, 1)**, at radius 3.550163, depth 0.202342 (24.7017% of mouth depth) and polar angle 0.261776 radians. There χ = 0, the signed physical volume ratio is −0.00168337 and the singular stretch magnitudes are 0.002417, 0.809405, 0.809540 and 1.062819. These characterize the inverted label map in melted material. At the other angular end of the same edge, vertex **(5, 18, 0)** has χ = 0.000115411, positive volume ratio 0.251629 and minimum stretch 0.368459. Continuity gives positive order inside the prism near that vertex. Every one of its six nodes therefore belongs to an active element. Fixing those six nodes fixes the entire prism map, including its negative determinant. This precludes the proposed zero-energy repair under the specified constraints.

**Decision.** Retain checkpoint 10. No label change, production equivalence claim or restart is authorized by this result. For continued evolution of the existing action, test the actual basis support and order derivative on a rejected trial. If every label must have positive volume, the mixed prism needs a different local discretization. A change in condensation energy would answer a separate physical question about equilibrium.

Independent checks reproduced the signed force history, the folded prism, the active-vertex incidence and the disjoint support of node 51603 on checkpoint 10. [Earlier melted-label controls](mouth_r4_melted_collar.md)
