# Does the puncture close, and where does the occupied wave go?

**Closure remains unresolved at both entrance radii.** Both model mouths still carry current through a finite terminal opening in the sampled states. The coarse radius-2 calculation gives a separate, positive result about the wave. Its saved oscillation is trapped on a six-node shear-bearing island while the mouth still drains. [Release and signed work](mouth_exit_release.md)

### Wave at the saved states

We rebuilt the mass and stiffness matrices at each saved order field and applied them to its saved wave vector. The kinetic energy on each connected component comes from the mass quadratic form. Radiated power comes from the outgoing traction acting on the vector at the port. No new eigenmode was selected. The table reports power per unit stored energy. The fine radius-2 entry is an intermediate state of the fixed-geometry order relaxation. [Wave operator](../tools/mouth_strained_wave.py) · [Outgoing port](../tools/mouth_exchange_shear.py) · [Measured component masses and physical outgoing power](../../../research/pde_ledger_v3/_measurements/mouth_wave_trapping.json)

| Entrance, mesh | Exit radius | Shear component occupied by the saved mode | Cavity kinetic share | Radiated power / stored energy |
|---|---:|---|---:|---:|
| 1, 48 | 0.46202 | Connected to port | 0.87243 | 0.022085 |
| 1, 72 | 0.44978 | Connected to port | 0.87995 | 0.020363 |
| 2, 48 | 1.01986 | Six nodes, isolated from port | 1 | 0 |
| 2, 72, iteration 2425 | 1.03697 | Connected to port | 0.91768 | 0.022626 |

On the coarse radius-2 grid the other component has 1,934 nodes and contains the port. Its kinetic share is zero, and the carrier through-current is 1.263651993023878. Iteration 1150 recorded one connected wave region at frequency 2.913391924010481. Iteration 1175 recorded the detached mode at 1.8125577812311795. The stored diagnostics bracket the change. To determine how the prepared wave action divides at the split would require time-resolved fields or resolved component projections across that event. The single-mode relaxation supplies no such partition. [Component selection](../tools/mouth_live_waist_active_wave.py) · [Order and wave measurements](../../../research/pde_ledger_v3/_measurements/mouth_exit_release.json)

### What would decide closure

The chart fixes the entrance, pins the terminal order at one half and represents the exit by $b(1)=1.5R\exp(x_{\rm exit})>0$. At an interior pinch or a closed terminal disk, a material section with positive initial radius would map to the axis. With finite positive axial stretch the bonded-core Jacobian $J=(r_{\rm cut}/r_{\rm cut,0})^3 L/L_0$ vanishes there. Both the stored-energy density and the incremental wave assembly use the inverse of this map. The carrier mesh also retains a distinct cap linking the two reservoirs. Replacing the log radius by a radius coordinate leaves these obstructions in place. This identifies a boundary of the chart. The energy cost of physically reaching it is still open. [Shape and aperture](../tools/mouth_sealed_shape_relax.py) · [Order boundary](../tools/mouth_live_waist_order.py) · [Core map](../tools/mouth_poly_bonded_core.py) · [Material attachment](mouth_exit_bond_audit.md)

To establish a sealed downstream exit with the entrance held, follow the minimum section to zero under free material and shape forces, continue through contact with a physical rule for the joining lattice and carrier phase, and verify that steady excess through-current over an intact sheet at the same bias disappears. Showing the entire puncture heals to a flat sheet would additionally require releasing the entrance constraint. To identify wave trapping across that event, one must also account for the prepared wave energy remaining in each supported region, transferred to the material and radiated through the port. A reproducible positive aperture with restoring work against feasible variations would instead establish local arrest. No sampled state meets either test.
