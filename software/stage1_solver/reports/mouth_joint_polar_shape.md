# Can the radius-four mouth carry its own load?

**The mouth stayed open and draining through the measured deformation. Its crystal still has a direction of lower energy, and the axial load at the exit rim remains unresolved. Free translation stays gated.**

One quantum of circulation threads the ordered mouth and continues into the superfluid bulk. We released depth, outlet radius, three interior width coordinates and a polar–equatorial distortion alongside the bonded core and collar. The conversion layer and both reservoirs respond. The angular wall displacement vanishes at both ends. The coarse solve took three energy-lowering steps. We transferred the final shape to the finer solid mesh and re-equilibrated its meridional material displacements and terminal fields. [Shape measurements](../../../research/pde_ledger_v3/_measurements/mouth_joint_polar_shape.json)

| At the same displaced shape, in solver units | 24 × 24 solid | 33 × 32 solid |
|---|---:|---:|
| Depth / minimum radius, reached at the circular exit | 0.795904 / 3.501226 | 0.795904 / 3.501226 |
| Largest absolute force among six shape coordinates | 1.21839 | 1.75746 |
| Largest free meridional material-node force | 1.12 × 10⁻⁵ | 4.71 × 10⁻⁶ |
| Principal stretches / minimum volume ratio sampled | 0.32920–1.70343 / 0.11846 | 0.30430–1.73581 / 0.12705 |
| Inward current / percentage retained against same-lattice unwound control | 11.61687 / 84.2427% | 11.61684 / 84.2426% |
| Incoming current both dense and subsonic, including swirl | 70.4773% | 70.4769% |

The coarse force on the width-asymmetry coordinate went from −19.00221 to +0.10421, and the polar force from +0.0840594 to −0.00202946. At depth 0.270608 the polar and equatorial radii were 3.837413 and 3.837151. Moving the crystal between pole and equator while holding the wall fixed lowered energy on both meshes. That independent coordinate had forces −0.00302046 and −0.00301984. This establishes an unbalanced internal load at the last verified shape. [Internal material test](mouth_polar_material_work.md)

The receiving-order residual reached 9.89 × 10⁻⁹ at this shape. The solid and front agreed on exchanged work within 1.5 × 10⁻¹³. The independent fine run from the original shape stopped on a material residual of 0.00003406 before accepting a displacement. After the tabulated state the coarse trial reached 0.00068307 and the force-controlled fine trial reached 0.00081288. Both exceed the material tolerance of 0.00003. These stops leave collapse and stability undecided. [Solver and independent order checks](../../../research/pde_ledger_v3/_measurements/mouth_joint_polar_shape.json)

**On the refined solid, the integrated axial stress in the draining state minus the resting state was +0.85088 on the incident side and +0.37667 on the receiving side.** Including side flux and crystalline work gave +0.41439 and +0.41374 at axial half-thickness 0.12. At half-thickness 0.03429 the same full disk gave +0.52913 and +0.41784. The interior agrees much more closely. The rim shares current with material outside the disk and carries reactions from the prescribed order and lattice. A physical force discontinuity has yet to be established. The next test should let the interior crystal move angularly at fixed wall, then measure force and constraint work on a surface enclosing all sides of the rim. Radius 3 remains blind. [Disk measurements and independent review](mouth_disk_pillbox.md)
