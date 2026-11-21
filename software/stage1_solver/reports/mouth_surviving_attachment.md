# Does one surviving attachment soften during conversion?

**The same upstream attachment lost 74.20% of its entry stiffness before breakup in the order-dependent run. In the constant-stiffness run it kept its full stiffness until breakup.** The force readings distinguish the two preparations. The carrier and order fields do not independently select a contact law. Both tests use the held radius-one 0PN mouth at intake 0.25, the same wall and ordered sites, and the same tag one eighth of a pitch off center. [Original junction and work histories](mouth_tagged_contact.md)

**Force measurement.** At each stage we froze the ordered endpoint, conversion and registration. We moved only the upstream wall anchor by \(\pm10^{-5}\) model lengths on separate copies, repeated at half strength, and subtracted the unperturbed force on the ordered endpoint. The table gives the increase in force per unit area for the positive full step, in units of \(10^{-6}\). Order is evaluated at the surviving attachment.

| Stage | Local order | Order-dependent attachment | Constant attachment |
|---|---:|---:|---:|
| Entry | 0.96916 | 4.84579 | 5.00000 |
| Junction center reaches half order | 0.62379 | 3.11893 | 5.00000 |
| Immediately before and after the neighboring attachment breaks | 0.47866 | 2.39331 | 5.00000 |
| Immediately before this attachment breaks | 0.25000 | 1.25000 | 5.00000 |
| Immediately after it breaks | 0.25000 | 0 | 0 |

Halving the displacement halved the force change. A second reading used only the older whole-junction forces and known relaxed offsets. Interpolating the affine force curve to the neighbor’s unloaded position isolates the upstream force even while both attachments are engaged. This check gives the same slopes without reading either stored spring stiffness. Both instruments read the same simulated material. [Calculation](../tools/mouth_surviving_attachment.py) · [Guarded readings and review scopes](../../../research/pde_ledger_v3/_measurements/mouth_surviving_attachment.json)

**Can the patch stiffness discriminate?** The earlier trial wall law gives stiffness per area proportional to order. A narrowing capture window also gives a surviving fraction proportional to order for uniformly registered sites. With equivalent parallel loads this constrains the mean stiffness of survivors. An individual off-center site can still soften. [Wall law](mouth_bend_static_shear.md)

For a concrete test, extend the earlier force that keeps two material features aligned to an isotropic well. On a copy held at registry \(\eta\) and zero sideways extension, its sideways stiffness is \(-F_{\rm registry}/\eta=S[1-(2\eta/\chi)^2]\) while attached. At \(S=0.25\) the off-center tag’s slope falls from 0.23336 at entry to 0.18180 at its neighbor’s loss in the reference history. It tends to zero at the edge. Yet the average among survivors stays \(2S/3\), giving patch stiffness proportional to \(\chi\). This prediction assumes a common metric and force law in the two directions and fixed registration. [Underlying registry force](../tools/mouth_contact_crossing_force.py)

The discriminating substrate measurement is a sideways push on a centered and an off-center surviving site at the same local order, with registration measured or held. For the off-center site, compare its measured sideways slope with the independently measured force holding it in registry divided by its offset, using the sign convention above. The work available for 0PN motion still depends on that material response.

Working base `49dca23e`. Two independent local reviewers checked the recorded scopes. Formal review and commit are pending.
