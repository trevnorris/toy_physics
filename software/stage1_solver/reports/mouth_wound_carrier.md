# A coherent vortex continued through both sides of a held exit

A single winding supplies outward work at both reference radii on the straight, centered vortex preparation. The same winding strongly depletes the carrier at the outlet. At unchanged reservoir phase difference, the normal flow alone exceeds the local long-wavelength sound speed throughout the sampled disk. The comparison uses the existing sharp conversion boundary, whose quiet planar calibration cannot be transferred to these states.

**Geometry and conserved quantity.** Let \(r^2=x^2+y^2+z^2\) in the three brane directions, \(s=\sqrt{x^2+y^2}=r\sin\vartheta\), and \(w\) denote depth. We put \(\Theta=N\phi+\psi(r,\vartheta,w)\) on both sides of the exit. In solver units, \(\oint d\Theta=2\pi N\), \(v_\phi=N/s\), and a linked vortex spans the \((z,w)\) plane. Both far reservoirs have this same straight vortex; their residual phases differ by 0.4510084282857874. The curved side has independent reflecting carrier traces on the inside and outside. The disk passes mass through the existing receiving half-layer of phase-delay length 0.04/12. The angular momentum flux through the disk equals \(N\) times the mass flux. Vortex preparation, stability against offset or splitting, and a localized return path are open. [Carrier solver](../tools/mouth_wound_carrier.py)

**Finite-core calculation.** Write \(\sqrt n=f_N(s)v\), where the positive straight-vortex solution satisfies
\[
 f_N''+f_N'/s-N^2 f_N/s^2=(f_N^8-1)f_N/2,
 \qquad f_N\sim a s^{|N|},\qquad f_N(\infty)=1 .
\]
Subtracting the energy of this same uniform vortex from all six spatial regions removes a common, box-extensive energy. Integration by parts gives amplitude-gradient energy \(f_N^2|\nabla v|^2/2\), potential \(f_N^{10}(v^{10}-5v^2+4)/20\), and boundary energy \(f_N\partial_n f_N(v^2-1)/2\). On the two sealed side faces the normals are opposite and the amplitudes may differ. The residual phase minimizes the native positive carrier energy including order entrainment in the core and the matching cap current. The three-dimensional reduced mesh retains both polar and meridional gradients. At every signed exit displacement we vary \(v,\psi\) while holding the preparation of the crystal fixed. [Density, speed and force reduction](../tools/mouth_wound_speed_readout.py)

| Entrance radius | Closing load at zero winding | Additional outward force, one quantum | Additional outward force, two quanta | Intake, zero → one quantum | Smallest sampled total Mach, one quantum |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.344251 | +0.520650 | +0.278595 | 0.289965 → 0.00304997 | 708.65 |
| 2 | 0.385372 | +0.770849 | +0.612592 | 1.331302 → 0.140538 | 14.727 |

The table compares signed 0.6% exit variations, 13 meridional cells per patch direction and 10 polar cells per hemisphere. Subtracting the zero-winding carrier reaction on the same mesh isolates the force to add to the archived total load. Enlarging the radial and depth boundaries from (28,22) to (38,31) gives one-quantum increments +0.517927 and +0.762542. The smallest nonzero quantum clears both loads but overbalances both held shapes. The nonmonotone one- and two-quantum readings preclude a quadratic interpolation to a balancing circulation.

At the equatorial rim sample the actual outlet radii are 0.44978 and 1.03697. The calculated total speeds are 2.42658 and 1.22859, while local homogeneous sound speeds \(c=n^2\) are 0.00342423 and 0.0834244. The minimum sampled absolute normal Mach numbers across the corresponding disks are 278.87 and 5.1201. The hydrodynamic Mach ratio becomes singular as density vanishes on the vortex axis; no conclusion relies on its maximum there. These checks establish failure of the normal calibration in this preparation. They do not establish a dynamical vortex instability or a physical steady conversion state.
