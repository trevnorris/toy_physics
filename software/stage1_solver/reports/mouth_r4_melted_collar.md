# Can the draining radius-four mouth survive loss of its crystalline collar?

The coarse survival test resumes from saved update 6 with the entrance held at radius four, one conserved winding, distant phase drop 0.4510084282857874 and the existing continuous order energy. The carrier uses the previous two-sided geometry with a sealed lateral wall. In this model an open, current-carrying slit can persist even where the crystal loses order. Mechanical equilibrium requires an independent balance of shape and material forces.

## What can be removed without changing the action

Where the spatial order field is exactly zero on a neighborhood, the crystal contributes zero energy and zero material force there. Its **one-sided response to restoring order** can still be finite. On the eight zero-order nodes surrounding the old limiting corner, the original solver predicts an energy slope of 0.2076556451473667 for increasing their order together. Direct perturbation gave 0.20765564515135004.

The test requires positive orientation wherever order is positive, including the ordered portions of mixed material cells. It uses no order threshold. Every oriented quadrature site contributes its original crystal load to the order field. In a melted region where the label map has folded, re-forming a crystal would require an oriented reference. Every such reference has nonnegative strain energy. We therefore solve the unrestricted order problem with that lower bound on the unknown cost. An accepted state must stay exactly melted on every affected order node, and the lower-bound first variation must be nonnegative there for every free node. A failure of that certificate leaves reordering unresolved.

## Independent fixed-field controls

| Test on the saved radius-four preparation | Result |
| --- | --- |
| Compare the old and new evaluators on the accepted state | Total energy, crystal energy, every material force, and every crystal-to-order load agreed exactly. The independently assembled full order gradient agreed within 8.89 × 10⁻¹⁶. |
| Move the limiting corner just past zero determinant | The signed chart determinant reached −3.3998961 × 10⁻⁹ between integration points. The whole-element check admitted the fold at zero order and detected it after order rose to 0.001. |
| Move one label farther inside the disordered shaft | 6,153 integration points crossed zero volume. Energy and material forces were unchanged. The order load changed by up to 0.02301817. All 128 affected nodes remained at zero order with a minimum lower-bound reaction of +0.00161120. |
| Reconstruct current from the original carrier field | Inward entrance, outward bulk and weak disk currents were 10.923297040080078, 10.923297040080131 and 10.923297040080053. |

The saved states establish that the shaft had lost order before its material coordinates nearly folded. They leave the direction of front propagation unresolved. A survival verdict requires a trajectory of order, exit radius, depth, intake and dense subsonic current share together with shape, material and field residuals. Any conclusion about physical collapse also depends on the prescribed carrier wall topology. [Original stop and spatial profile](mouth_r4_fold_gate.md) · [Control measurements](../../../research/pde_ledger_v3/_measurements/mouth_melted_collar_controls.json)

**Independent review**

> Scoped clear for continuing the held 24 × 24 radius-four relaxation with orientation required throughout ordered material and a tested one-sided lower bound for reordering in melted cells. Physical survival, continuum convergence and free release require the endpoint evidence.

**Checkpoint views.** The equatorial cutaway preserves the physical proportions. Both field sections enlarge depth by 2.7× and use fixed scales across the three saved updates. The exit narrows from 3.35189 at updates 6 and 7 to 3.32269 at update 8. Across all three, the entire shaft has zero crystalline order beyond physical depth 0.28460. The dark center of the density map follows the vortex core. [Three-frame sequence](figures/r4_melted_collar/checkpoint_sequence.gif) · [Checkpoint 8 at full resolution](figures/r4_melted_collar/checkpoint_08.png) · [Render inputs and measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_survival_views.json)
