# Force on the first bend while the same material slips, shears and crosses the gates

**After one measured material transit, the signed waist-to-middle wall force was −3.38% of its early value in the reference strip.** The twelve sections that supplied this clock also transmitted the sideways force. The result depends on continuing the neutral material's positive shear stiffness into the compressed strip. A separate energy test found that the existing microscopic bond law favors sideways displacement there even with the wall attachments held fixed. [Guarded measurements, tag histories and local reviews](../../../research/pde_ledger_v3/_measurements/mouth_joint_strip_hold.json)

We held the first bend of the radius-one cubic mouth at intake 0.25. Longitudinal Morse forces and carrier drag set the site positions and arrival times in an open strip. Nine contact families per twelve-site section sampled the passing material and loaded the same exterior brane stress field. We connected corresponding sites along the flow using the shear stiffness measured on the neutral sheet, scaled by local shear order. At first order a sideways bend leaves the longitudinal motion unchanged in this model. Wall coverage, smooth longitudinal sliding and the correspondence between neighbors are material assumptions. [Joint equations and force readout](../tools/mouth_joint_strip_hold.py)

| Time since the bend was displaced | Model time | Signed band force / early force |
|---|---:|---:|
| Early reference | 0.0625 | +100% |
| Mean material transit from waist to middle | 6.92037 | −4.31% |
| Mean material transit from waist to depth 2 | 13.32770 | −3.38% |
| Earlier full carrier transit | 28.43103 | −3.24% |

The refined early force was +4.34174 × 10⁻⁶ for a displacement of 10⁻⁵. The cut stress independently reproduced the integrated wall reaction. The undisplaced copy gave zero, and reversing and halving the displacement reversed and halved the force. At the middle gate, the mean angular displacement of tagged site 12 relative to the average of its neighbors eight sections upstream and downstream was +5.86 × 10⁻⁸. The measurement file includes each tag's three crossing times and the wall force at those crossings. The table reads at elapsed transit durations. This strip extends to depth 2, beyond the middle outlet of the previous row test. [Earlier passage experiment](mouth_connected_passage.md)

**What sets the late force.** On the reference spatial grid, removing inter-section shear changed the full-material-transit fraction from −3.38% to −4.80%. Extending only the inlet from 8 to 24 bond lengths gave +0.71%. Extending only the outlet gave −4.10%. Each fraction uses that preparation's own early force and measured transit.

**Can the microscopic bonds support the centered strip?** We reconstructed staggered triangular rows on a measured transport snapshot and displaced small groups sideways with both the wall tips and exterior held fixed. All nine trial displacements lowered the total energy to second order on each of two carrier meshes. For the eight-spacing disturbance at the middle, the bond curvature was −44,384.57 and the wall contacts added +255.44, leaving −44,129.13. The positive reference-shear law gave +1,336.74 for the same trial including the wall. The centered microscopic configuration therefore has a direction in which it can lower its energy. The joint force prediction describes the reference-shear model. [Energy calculation](../tools/mouth_joint_strip_stability.py)

**Next 0PN test.** At mid-depth, keep the measured feed, axial compression, hoop size and wall contacts. Release the lateral site positions from opposite weak perturbations. Check whether the driven material reaches a reproducible rearranged state and measure its force under small opposite shear steps. That would supply a shear response for repeating the held-bend test. The force on a freely accelerating mouth remains open.

Working base `49dca23e`. Independent local review covered the conditional force, tag clocks and constrained energy test. Formal project review and commit are pending.
