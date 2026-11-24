# Where the unpinned radius-four crystal first folds

The 24 × 24 relaxation stopped at saved update 6. One further accepted update raised the largest absolute shape force from 0.74525349 to 0.77343822. The minimum chart determinant stayed at $3.80182\times10^{-16}$ at the same vertex, and the material solver again rejected a step that would fold a crystal element. At the preserved state the mouth still had depth 0.81770231, exit radius 3.35189208 and inward port current 10.92329704. These are measurements on an unbalanced shape. The entrance radius, one quantum of circulation, phase drop 0.4510084282857874 and carrier potential $U(n)=(n^5-5n+4)/20$ were fixed throughout.

### The folding corner

The limiting vertex lies at 24.6607% of the physical depth in the first polar cell, 15.00° from the winding axis. Its spherical radius is 3.55416 and its perpendicular distance to that axis is 0.91981. The determinant comes from one side of the lower triangular prism with depth/radial/angular indices (5, 17, 0). [Material measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_fold_constitutive.json)

| Quantity at saved update 6 | Measured value |
| --- | ---: |
| Crystalline order and its spatial gradient | Both exactly zero in the containing order-field cell |
| Physical volume ratio $J$ | $2.51\times10^{-16}$ |
| Principal stretches, smallest first | $3.60\times10^{-16}$, 0.81032, 0.81046, 1.06249 |
| Coherent carrier density $n$ | 0.18460748 |
| Crystal energy density $W$, per undeformed volume | 35.39064 |
| Actual crystal energy density $\chi W$ | 0 |

The weakest stretch runs almost along the depth axis (absolute axial component 0.99834). At the center of this prism $J=0.12346$. All points of the parent prism satisfy $0\le\chi\le0.000139858$, and one vertex has $\chi=0.00011624036$. The upper bound encloses the entire mapped cell. Direct reconstruction with exact sums of the saved material offsets gives a positive chart determinant of $4.70\times10^{-16}$. Rounding the landed node positions first gives $3.80\times10^{-16}$. The physical $J$ in the table uses the former convention and the production reference volume.

### Order from entrance to exit

Across the full shaft, crystalline order is at most 0.1 from 7.4932% depth onward, at most 0.01 from 18.5081%, and identically zero from 34.8054% through the exit. At the requested 97.765% section it is zero over the full disk, including its edge. It is also zero around the coarse mesh’s last free material ring, which landed near 97.1716%. The fluid continues to drain through this section. [Section plot](figures/mouth_r4_fold_section.png) · [Order profiles and bounds](../../../research/pde_ledger_v3/_measurements/mouth_r4_fold_order.json)

The reconstructed seed for the held, undriven preparation had zero of 30,085 shaft nodes at or below 0.1 order. By its first surviving spatial snapshot, update 60, 29,413 were at or below 0.1 and 25,562 were exactly zero. The later coarse relaxation began with $\chi=0.00130405$ at the eventual folding vertex. It reached zero by coarse update 1, before the determinant approached machine precision. The available snapshots establish widespread order loss before the fold. Earlier spatial snapshots would be needed to distinguish a front arriving from the exit from local nucleation. Optimization updates do not measure propagation in physical time.

**Decision.** With $\chi=0$, the local $-\chi\ln J$ resistance drops out even while the solver requires positive orientation of every crystal element. The run reached that kinematic limit inside material that had already lost its crystalline stiffness. This leaves continuum collapse and the fate of a mobile exit open. The calculation stopped before any change to the mesh or material law. Radius four remains held, and the blind radius-three prediction remains available.

**Independent review**

> Scoped clear for stopping the coarse radius-four descent at accepted update 6, for loss of the crystal’s compressive penalty in a zero-order neighborhood of the nearly singular shaft corner, and for prior loss of order at the held exit ring. The equilibrium exit radius, physical collapse and free fall remain unresolved.
