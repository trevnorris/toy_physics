# Can the solved elements determine the axial force on the radius-four exit lip?

**The two sides of the lip still fail the independent wall-force test.** We used the saved draining and resting fields on the same fixed wall, with an ordered interior and one circulation quantum. We recovered axial forces on the duct and receiving sides of the last 0.03429 of depth. Two independent calculations used the element fields before comparing their predictions with the prior signed wall-displacement measurements.

| Rim-refined mesh, axial force | Duct side | Receiving side |
|---|---:|---:|
| Direct signed work, displacement extended 4 cells into material | −0.04765548 | +0.02361719 |
| Direct signed work, same wall motion extended 8 cells | −0.04623509 | +0.02339762 |
| Force assembled on wall nodes alone | −0.04631421 | +0.02799741 |
| Earlier pointwise wall stress, integrated over the same motion | −0.00453854 | +0.06160983 |

Forces use draining minus resting at the same winding. The table uses four-point integration per axis and reports the reaction opposing the energy derivative. Including the forces on displaced **interior** nodes reproduces the signed work to within 6.84 × 10⁻⁹ across both meshes, integration orders, and three wall strips. But the same terminal wall motion gives different work when carried into four or eight cells. A single surface traction must miss at least one reading by approximately 0.00071020 on the duct side and 0.00010978 on the receiving side. Halving the signed step changes either by less than 0.000000010. [Nodal forces and signed controls](../../../research/pde_ledger_v3/_measurements/mouth_terminal_element_recovery.json)

**What accounts for the difference.** Integrating axial momentum separately within each element isolates the coherent carrier as the largest contributor. For the eight-cell motion on the refined mesh, its interior contribution relative to the directly sampled wall stress is −0.04602 (duct) and −0.03713 (receiving). On the duct side, +0.03469 from stress jumps between elements partly offsets −0.08071 from stress variation inside elements. The corresponding receiving-side terms are −0.06009 and +0.02296. Other energy terms account for the remainder. This identifies where the work enters the calculation. How much belongs to the physical support of the conversion front remains open. [Element momentum budget](../../../research/pde_ledger_v3/_measurements/mouth_terminal_element_sources.json)

| Closed axial residual with independently measured exterior, near / far | Doubled grid | Rim-refined grid |
|---|---:|---:|
| Original integrated stress | −0.05342362 / −0.05330769 | −0.07774070 / −0.07783445 |
| Assembled wall nodes, shared rim counted once | −0.02303795 / −0.02344811 | −0.04042705 / −0.04078059 |
| Conservative face flux, extrapolating the interior source to the boundary | −0.05395794 / −0.05301948 | −0.06880558 / −0.06848992 |

All three conservative flux reconstructions balance every cell to within 4 × 10⁻¹⁶. All miss both independent four-cell wall forces on the finer grid. Replacing the exterior *as well* with the projected flux recovers the imposed cell balance by construction. Keeping the exterior measurement independent exposes the grid drift. [Conservative flux measurements](../../../research/pde_ledger_v3/_measurements/mouth_rim_conservative_faces.json)

**Decision.** Keep the radius-four lip fixed. The 0.176506 minimum stretch about 0.0722 inside the edge comes from denser sampling of the earlier locally relaxed patch. This force diagnostic does not establish a corner singularity. The next discriminator is a fixed-wall test with the same signed wall shift and the same physical crystal and receiving-order constraints. Let the fields re-equilibrate, measure work on those constraints separately, and check whether two ways of moving the interior mesh give the same wall reaction. Free release requires a common wall force that also converges across grids. Radius three stays blind. [Compression record](mouth_terminal_force_balance.md) · [Independent review](mouth_rim_traction_recovery_review.md)
