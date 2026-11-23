# An open, draining radius-four mouth with its entrance held

The mouth stayed open when its depth, exit radius, three depth-profile coordinates, polar deformation and internal crystal displacements relaxed together. Both solid meshes passed the predeclared force tolerances. The result applies to a stationary, bonded lattice with an ordered interior and one circulation quantum. The entrance radius was held at four. Motion of the center relative to the brane remains an open question. [Equations and preparation](mouth_whole_r4_relaxation.md) · [Endpoint measurements](../../../research/pde_ledger_v3/_measurements/mouth_whole_r4_endpoints.json)

| Measurement in solver units | Coarse solid, 24 × 24 | Fine solid, 33 × 32 |
|---|---:|---:|
| Depth | 0.792701 | 0.792279 |
| Narrowest radius, at the exit | 3.501793 | 3.502069 |
| Polar / equatorial radius at 34% depth | 3.837574 / 3.837222 | 3.837752 / 3.837399 |
| Largest shape force, full / half probe | 0.00037549 / 0.00037521 | 0.00014200 / 0.00014219 |
| Largest crystal-node / ring force component | 0.00000614 / 0.00002501 | 0.00002414 / 0.00002546 |
| Sampled principal stretches | 0.34554–1.58371 | 0.31987–1.60743 |
| Minimum sampled volume ratio | 0.14708 | 0.12709 |
| Inward current / fraction of unwound intake | 11.53793 / 84.227% | 11.53918 / 84.230% |
| Dense, subsonic share in the feeding-side terminal cells, including swirl | 72.119% | 72.125% |

**The earlier solid stall cleared.** On the fine mesh the last free ring, at 97.765% of the depth, had radial and axial forces of −0.00001073 and −0.00000435. Every accepted map retained orientation. At the sample stations 0.0722 inside the disk edge, the minimum stretch on the fine mesh went from 0.32594 before descent to 0.32577 afterward. The coarse endpoint gave 0.40817 there. The compressed zone persisted during descent and its strength still depends on solid resolution. The two calculations establish first-order stationarity in the allowed coordinates. A test of stability to general perturbations would also need the coupled stiffness, including symmetry-breaking directions.

**The current and force checks answer different questions.** The resting wound control had zero through-current. The fine driven state gave the same variational feed and bulk currents to within $7.1\times10^{-12}$. Independently integrating the physical current from the two sides of the disk gave 11.55207 and 10.94991, a difference of 5.22% of intake. On those two faces, dense, subsonic material carried 71.36% and 74.74% of the positive current. Both solids used the same carrier mesh. The fine near and far axial readings were +0.54323445 and +0.48981577. Including work in the intervening material brought them within $3.4\times10^{-11}$. This checks integration against the same displacement. Including the slit faces directly left axial sums of −0.05760 and −0.05614. The enclosed order, crystal and fixture reactions have yet to be identified independently. These sums alone cannot assign a force to a freely moving mouth.

**Release decision.** Preserve both endpoints for held-state measurements. Keep radius three blind. Resolve the physical face flux and the residual enclosed force on the held fine shape before attempting capture, push inertia or translation. Widening the mobile order band on the held fine shape restored a radial force of −1.74663 on the last free ring. Physical matching of the conversion layer now determines whether this shape can support itself. Review velocity damping, vortex anchoring and material mobility before designing a release. [Order band and entrance fixture](mouth_whole_r4_supports.md) [Protocol](../../../docs/zero_pn_free_fall_roadmap.md#after-both-radius-four-relaxations-finish)

Independent endpoint review, covering both guarded exits, 69 pinned source and input files and 31 accepted snapshots at solver revision `fa5bed3d` · “Scoped clear for first-order stationarity of the two saved, fixed-entrance, axisymmetric radius-four preparations and for positive weak drainage. Local edge strain convergence and closed physical momentum remain unresolved.”
