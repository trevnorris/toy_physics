# Radius-four relaxation with exact order support

At checkpoint 10 the 24 × 24 crystal has inverted material labels inside the melted shaft. The existing action permits these wherever crystalline order is exactly zero. The new guard asks whether there exists any point in a material prism with both positive interpolated order and nonpositive signed volume. It tests the actual support of every nonnegative spatial basis function, including knots, prism faces and the polar axis. Rational subdivision certifies separated regions. An exact polynomial feasibility test resolves any remaining intersection, including a tangency at zero volume. The energy, constitutive law and physical boundary conditions are the same as in the [checkpoint-10 calculation](mouth_r4_checkpoint10_boundary.md).

On the original checkpoint, old and new evaluations agreed bit for bit on total energy (92185.53614224629), all 10,074 internal crystal forces, all 12 signed shape-probe energies, all six resulting shape forces, every field derivative and all eight current readings or current differences. Setting order to 0.001 at folded crystal corner (5, 18, 1) triggered rejection in the production virtual-work evaluator. [Measurement record](../../../research/pde_ledger_v3/_measurements/mouth_r4_exact_support.json)

## Reconstructed full shape step

The saved optimizer history reproduced the checkpoint-10 full-step rejection verbatim. On the rematched field, node 51603 had order 0.000036403035309499604 and a signed lower-bound energy derivative of −7.1967443814220655 × 10⁻⁹. An independent check certified disjoint support on all eight folded prisms, both for this node alone and for the complete order field. The rejection came from including a nearby node in a bounding box.

| Held-mouth reading | Checkpoint 10 | Full shape proposal |
|---|---:|---:|
| Depth | 0.81914036 | 0.82008883 |
| Exit radius | 3.30637727 | 3.27757267 |
| Inward feed current | 10.70857269 | 10.57367995 |
| Dense, subsonic share of inward disk current, including swirl | 66.7562% | 66.1106% |
| Depth from which the full disk has zero order through the exit | 0.2846 | 0.2846 |

The full shape proposal uses the checkpoint's internal crystal offsets and rematched carrier and order fields. Its feed and bulk currents agree within 1.36 × 10⁻¹³. The production force target is 0.0004. The entrance stays held at radius 4 with one quantum of circulation, phase drop 0.4510084282857874 and fluid potential \\(U(n)=(n^5-5n+4)/20\\).
