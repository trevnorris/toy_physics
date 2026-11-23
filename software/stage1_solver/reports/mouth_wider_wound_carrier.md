# Where does a widened exit carry its inward current?

A straight vortex of one quantum extends from the feeding brane into the receiving bulk. Both reservoirs retain the original meridional phase difference. The two-sided carrier adjusts its density and meridional phase at each held exit. For a reference meridian of depth \(L\), the path is \(r(w;b)=r_{\mathrm{ref}}(w)+(b-b_{\mathrm{ref}})w/L\) inside the throat. The entrance radius, depth, other meridional coefficients and material labels are held fixed. The exit coordinate \(b\) therefore tests the response of the entire registered collar. The disk has intrinsic three-volume \(4\pi b^3/3\). [Solver and signed work](../tools/mouth_wider_wound_carrier.py)

At 13 meridional by 10 polar cells per subdomain, the denominator of each flux fraction below is the integrated positive inward current through the exit. The normal Mach number is \(|j_w|/n^3\) since local long-wave sound speed is \(n^2\). Total Mach includes the radial, polar and azimuthal velocities. The density cutoff for the last column is one half of the far-field carrier density. It characterizes carrier density independently of crystalline order. Force here is the carrier contribution by itself, from the signed derivative of two-sided excess energy for exit shifts of ±0.6%. The total includes order and crystalline elastic work.

| Entrance | Exit \(b\) | Carrier force | Intake | Minimum normal Mach | Flux below normal sound speed | Flux below total sound speed | Flux dense and below normal sound speed |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1.040 | +0.317275 | 0.141966 | 5.0342 | 0% | 0% | 0% |
| 2 | 1.042 | +0.316909 | 0.142914 | 4.9783 | 0% | 0% | 0% |
| 2 | 1.35 | +0.232948 | 0.330539 | 1.0551 | 0% | 0% | 0% |
| 2 | 1.65 | +0.099838 | 0.581375 | 0.30275 | 27.93% | 0% | 4.30% |
| 2 | 2.00 | −0.091985 | 0.928636 | 0.08428 | 50.99% | 0% | 17.56% |
| 2 | 2.45 | −0.305288 | 1.409449 | 0.01538 | 69.40% | 8.25% | 36.94% |
| 2 | 3.50 | −0.391878 | 2.423349 | 0.0000695 | 88.18% | 33.54% | 64.22% |
| 1 | 0.455 | +0.263155 | 0.003219 | 265.80 | 0% | 0% | 0% |
| 1 | 1.00 | +0.278874 | 0.065978 | 3.7039 | 0% | 0% | 0% |
| 1 | 1.35 | +0.213787 | 0.135695 | 0.38958 | 23.08% | 0% | 0% |
| 1 | 1.65 | +0.130911 | 0.188797 | 0.07552 | 50.52% | 0% | 3.10% |

The same-exit, same-bias unwound intakes are 1.337372 at radius 2 and exit 1.040, 5.004661 at exit 3.50, and 0.295930 at radius 1 and exit 0.455. Thus widening the radius-2 exit to 3.50 restores 48.42% of the matched unwound intake. Even at exit 1.040, 91.85% of radius-2 inward flow is farther than half the outlet radius from the vortex axis. In this geometrically outer annulus every sampled point still crosses faster than the local normal sound speed. A 16 by 12 carrier refinement with ±0.3% exit steps gives a carrier force of +0.326529 there and minimum normal Mach of 4.9711. The measured flux at exit 1.040 falls from 0.141966 to 0.141074 under that refinement. At radius 1 and exit 0.455, the same refinement changes the carrier force from +0.263155 to +0.267077, the minimum normal Mach from 265.80 to 256.79, and the intake from 0.0032190 to 0.0030891.

The local disk phase jump in this sweep uses the inherited sharp receiving half-layer. These data measure where the boundary condition carries mass and whether its calibration transfers. Determining the transmitted stress requires the separate resolved-layer calculation. Stability of a displaced vortex remains open.
