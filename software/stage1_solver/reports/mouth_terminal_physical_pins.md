# Radius four: which constraint stays fixed when the exit wall is displaced?

The test uses the same ordered radius-four mouth, one circulation quantum and the same two saved carrier meshes. The signed axial displacement rises from zero at depth 0.76046874 to its maximum at 0.77761160 and returns to zero at the disk at 0.79475445. Its continuation into the material spans either four or eight cells of the original doubled grid on **both** meshes.

**Bonded duct.** The duct-wall crystal and the wall share a position. The proposed axial step has a normal component of 0.20475660 per unit at the peak. A step of 0.000020 sweeps 0.00001105071 of four-dimensional volume into or out of the crystal domain. Holding those physical crystal sites stationary while shifting the bank would move the wall through its own material. A control subtracts the eight-cell displacement from the four-cell displacement. The wall displacement then vanishes everywhere. This tests whether moving only the interior mesh produces an apparent force with the crystal, wall and receiving constraints all fixed. It yields a difference between the two extensions at first order. Individual duct wall reactions would require specifying how the actuator moves the bonded material.

**Receiving bank.** The saved equations allow a separate virtual shift of this carrier-field boundary. At each moved node we resample the prescribed receiving order at its physical position. The crystal stays fixed on the duct side. Resampling changes the energy compared with carrying the old nodal values along with the mesh. This contribution is measured separately. A stationary external support has zero mechanical displacement work. Its field representation can still contribute to the calculated variation.

| Four-cell minus eight-cell reaction, draining minus resting | Doubled grid | Rim-refined grid |
|---|---:|---:|
| Bonded duct, fixed-wall interior motion | +0.00062939 | +0.00017420 |
| Receiving bank, physical order sampled at moved nodes | −0.00011981 | −0.00002156 |
| Receiving bank, saved order evaluated directly in laboratory coordinates¹ | −0.00016597 | −0.00007831 |

These are **predictions from the saved stationary fields**. The receiving values are the mean of the positive- and negative-step derivatives. The receiving-order interpolation has different one-sided slopes. On the finer grid the four-cell predictions for negative and positive steps are +0.02545792 and +0.02572784, respectively, compared with +0.02547948 and +0.02574939 for eight cells. The refined duct prediction consists of −0.00001714 from the carrier at fixed nodal values, +0.00018039 from sampling the field at fixed crystal sites, and +0.00001095 from the imposed order. [Separated contributions](../../../research/pde_ledger_v3/_measurements/mouth_rim_pinned_constraint_work.json)

¹ The independent check evaluates the saved piecewise-linear receiving-order field directly at the moved integration points in laboratory coordinates. Keeping that field fixed off the nodes changes the apparent rate of convergence while preserving the sign of the predicted gap. [Independent review and controls](mouth_rim_physical_pin_review.md)

The signed test lets the free density, phase and conversion order settle for each displacement and compares both step sizes with these predictions. The terminal wall probe vanishes at the common disk rim, where the enclosing-surface probe has unit weight. A terminal-strip result alone cannot update the previously measured whole-rim drift or establish equilibrium for free release. [Previous axial measurements](mouth_terminal_element_recovery.md)
