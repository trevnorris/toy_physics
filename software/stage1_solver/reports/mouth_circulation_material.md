# What happens to order and strain around a static vortex?

**One quantum of circulation leaves the previously melted radius-1 site energetically free to reorder.** The material prefers a small local increase of crystalline order on both held preparations tested there. Reordering an entire surrounding patch costs energy. These two observations leave whole-mouth healing open.

The homogeneous coherent carrier with asymptotic density one has radial amplitude \(f(s)\), where \(s\) is distance perpendicular to the vortex sheet. The existing carrier law gives
\[
 f''+s^{-1}f'-N^2 f/s^2=(f^8-1)f/2,\qquad n=f^2 .
\]
Regularity at the core and the far-bath density determine the profile. For one quantum the half-density radius is 1.5479513 and the homogeneous long-wave speed \(c=n^2\) equals \(|N|/s\) at 2.0499349. For two quanta these distances are 3.1811918 and 3.9142614. Changing the outer radius from 32 to 64 and tightening collocation tolerance preserves the unit-winding distances within \(2\times10^{-9}\). The integrated radial and azimuthal stresses satisfy cylindrical force balance within \(2.2\times10^{-11}\) in the finer unit case. In this neutral control a flat equilibrium order profile \(\chi(w)\) coexists with the same vortex on each side. The azimuthal current is tangent to order contours and the canonical order momentum vanishes. The core therefore depletes coherent density without requiring the crystal to melt. [Native vortex and material probes](../tools/mouth_circulation_material.py)

**Order on the draining mouth.** At fixed density and meridional phase, the direct swirl energy is \(nN^2/(2s^2)\). Its derivative with respect to axisymmetric order is zero. The wound density and meridional flow can also change the reciprocal carrier–order stress. To measure that indirect effect we integrated it over the original 145-by-73 order grid using the reconstructed wound and unwound fields. Both followed the same crystal registration. Full variation of the native carrier energy with respect to the tested order direction at fixed reconstructed fields agreed with the analytical derivative. [Order-work calculation](../tools/mouth_winding_order_readout.py)

| Entrance and order probe | Starting order at the probe center | Energy slope at zero winding | Energy slope at one quantum |
|---|---:|---:|---:|
| 1, previously melted site | 0 | −0.0001483874 | −0.0001485177 |
| 1, surrounding patch | 0 | +0.7793832 | +0.7793786 |
| 2, same material-coordinate site | 0.9790313 | −0.0002518583 | −0.0002519652 |
| 2, surrounding patch | 0.9790313 | +0.0658974 | +0.0658984 |

Slopes are energy per positive order-probe amplitude. The reference site is at depth fraction 106/144 and radius fraction 43/72. The continuous patch raises order only within its compact footprint, leaving fully ordered nodes fixed. Two carrier meshes reproduce all signs for one quantum. With two quanta the coarser radius-1 calculation still gives −0.0001486970 at the melted site and +0.7793707 for the wider patch. A single material node occupies less than one carrier element; integrating its footprint avoids falsely assigning it zero carrier work from a missed cell center. These are held-shape derivatives on reconstructed fields. The total energy and an order profile after free relaxation are separate measurements. [Original site selection](mouth_exit_release.md)

**Collar.** The fixed prestrained registration has minimum sampled volume Jacobians of 0.0037405 and 0.0090115 across the two collars at entrances 1 and 2. The outer collars reach principal stretches from 0.04705 to 5.02472 and from 0.01573 to 9.76640. They remain orientation preserving in the finite-strain energy. The held geometry has zero incremental displacement by construction, and the law supplies no independently measured damage threshold. No balanced vortex-loaded state was available for a test of its released collar strain.
