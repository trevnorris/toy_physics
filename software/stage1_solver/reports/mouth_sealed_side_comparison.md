# What changes when the lateral wall of a draining mouth is sealed?

**Sealing the lateral wall reverses the carrier’s stiffness on the radius-1 test flare.** We kept the geometry, coherent material on both sides, and two reservoir phases identical and changed the location where current can cross. This directly tests whether the 95.62% lateral discharge in the earlier calculation matters for shape support. [Previous permeable-wall calculation](mouth_exchange_carrier.md)

| At the same reservoir phase difference | Permeable lateral wall | Sealed lateral wall |
|---|---:|---:|
| Total intake | 1.89673719 | 0.80912217 |
| Through the lateral wall | 1.81363280 | 0 by construction |
| Through the terminal melting disk | 0.08310440 | 0.80912216 |
| Carrier energy slope for a widening neck | +0.766960 | +0.222293 |
| Carrier stiffness for the same change | +0.193441 | −0.500150 |

The two solved fields are the steady density and phase at a far phase difference of 0.4510084282857874. The intake adjusts in response. Positive displacement widens the entrance and narrows the flare at fixed length and volume. The generalized force is minus the energy slope. Both branches favor narrowing the entrance at the original shape. On the sealed branch, the carrier also lowers the second derivative of the energy, which increases the stiffness needed from other sources. The change of boundary condition includes pressure and density-gradient stresses as well as flow. It does not isolate a pure Bernoulli force. The original shape was a turning surface of constant mean curvature.

**Working boundary for the free-profile tests.** The intact slab face and the lateral duct both have zero normal carrier current. The bulk outside the duct remains present and can exert pressure. On the permeable terminal disk the carrier has the same conservative Route A phase delay as before. The amplitude is matched at the disk, with free amplitude traces on the sides of the sealed wall. Equal surface tensions on the wall and on the terminal order interface are an explicit additional assumption for the relaxation. Their work and the wave energy are separate from the carrier-only table.

At the finer mesh (64) the energy derivative from signed steps and the independently integrated two-sided stress agreed in slope to \(4.7\times10^{-6}\). The stiffness from stress was −0.500140. Doubling the step gave −0.500119 from energy. At mesh 40 the corresponding result was −0.500107. This comparison tests the consequences of a seal. The formation of a seal and the absolute physical conductance of the terminal front are open material questions. [Sealed and permeable field results](../../../research/pde_ledger_v3/_measurements/mouth_sealed_side_comparison.json) · [Boundary selection in the roadmap](../../../docs/zero_pn_free_fall_roadmap.md)
