# What resists thinning around one active drain at 0PN?

**Result.** At fixed inward rate the brane carrier pushes outward less strongly near the drain than in the undisturbed bath. In the slowly varying far-field calculation, squeezing the slab weakens that outward push further. The net force that restores the two faces to their original separation remains undetermined.

**Measured on the saved slab.** Take an isolated radius-one drain, inward number rate $Q=0.25$, and distant density and prescribed width both equal to one. Write the face positions as $w_\pm=Z\pm W/2$. If $P_\pm$ denote the carrier's outward normal momentum fluxes above the bath reference, virtual work gives $(P_+-P_-)\delta Z+(P_++P_-)\delta W/2$. On the intact reflecting faces the carrier tensor gives $P_\pm=(n^5-1)/5-\partial_w^2 n/4$. I evaluated the normal density derivative directly and independently recovered the normal field curvature from the steady carrier equation. Both reconstructions used the same archived fields. The table gives sphere averages on the largest grid domain ($L=48$), in solver stress units. [Stress tensor](../../../docs/drain_driven_initial_acceleration.md#held-two-drain-attraction-measured-in-the-surrounding-brane) · [Readings on three grids](../../../docs/_measurements/active_drain_width.txt)

| Distance in mouth radii | Carrier load conjugate to width, positive outward | Same load from the field equation |
|---:|---:|---:|
| 2 | $-4.01340058316\times10^{-5}$ | $-4.01285958603\times10^{-5}$ |
| 8 | $-5.10997670147\times10^{-8}$ | $-5.10997663328\times10^{-8}$ |
| 16 | $-3.05760382526\times10^{-9}$ | $-3.05760383636\times10^{-9}$ |

**Response to a prescribed squeeze.** In a broad, locally flat exterior annulus assume steady radial subsonic flow, uniform through the depth, slow variation relative to the slab depth and density healing length, and fixed distant chemical potential. At leading order omit density-gradient stress and the quantum potential. With $I=Q/(4\pi r^2)$, inward speed $v=I/(Wn)$, $P=n^5/5$, and $(n^4-1)/4+v^2/2=0$, differentiation gives

$$
\left.\frac{\partial P}{\partial W}\right|_{I,\mu_\infty}
=\frac{n^5v^2}{W(n^4-v^2)}>0 .
$$

At eight radii, narrowing the annulus by 1% lowered pressure by $9.8096406052\times10^{-10}$ relative to the unperturbed annulus. The central slopes from the $\pm2\%$ and $\pm1\%$ width pairs were $9.67047729153\times10^{-8}$ and $9.66467529412\times10^{-8}$. Both engines gave the differential slope $9.66274245529\times10^{-8}$ to the quoted precision. These are local imposed-width calculations. The archived fixed-width fields supplied the baseline loads in the table. [Calculations and assumptions](../../../docs/_measurements/active_drain_width.txt)

**Why the physical width is still open.** In the uniform quadratic material energy, sealing a column against exchange and lateral dilation gives a width stiffness larger than or equal to the freely density-relaxed value by $B(1/W_0-C/B)^2$ for $B>0$. The active-drain equations select neither ensemble. Even adding $\lambda(W-W_0)^2/2$ would leave every pinned-width baseline load unchanged while changing its incremental stiffness by $\lambda$. This illustrates the missing datum. Bulk loading, exchange and material stress must be included on the same faces. [Width equation and limits of the historical breathing calculation](../../../research/pde_ledger_v3/steps/S11bB_interface_assembly.md) · [Live normal response and bulk loading](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md)

**Implication for gravity.** The measured carrier load extends into the exterior. If the *net* width load retains a tail $-A/r^4$, the diagnostic equation $(K-D\nabla^2)\delta W=f$ with constant $K,D>0$ gives a forced tail $-A/(Kr^4)$. For $K=0$, the exterior contains a particular term $A/(2Dr^2)$ and may contain a $1/r$ term. Its coefficient depends on the complete signed width source, including the inner boundary. The exterior carrier reading alone cannot establish inverse-square attraction. [Width diagnostic](drain_density_thickness.md)

**Decisive measurement.** On the same active preparation, apply equal and opposite laboratory forces to the two faces of a broad exterior patch. Compare the changes in $W=w_+-w_-$ and $Z=(w_++w_-)/2$ under reversed and half-amplitude loads, including density and the bulk reaction. Repeat at increasing patch sizes to resolve the long-wavelength response. The fixed-width solver cannot perform this measurement. It needs independently located material faces and a normal momentum law or measured net traction versus separation at that intake. Working result at fd667ea9, pending review.

## What an applied pressure patch would measure

**One-drain test status.** The saved isolated intake is $0.25$. Its archive contains the carrier, current and fixed grid. The solver sets the depth before solving, and the adopted shear-holding faces have no identification with a chosen density contour. The live bulk normal-load response and the brane's normal constitutive response are both open. The earlier single-interface tension and uniform breathing calculations do not supply those responses at this intake. Consequently there is no physical reading of thickness under an applied load in this preparation. [Saved input and solver evidence](../../../docs/_measurements/width_pressure_patch.txt) · [Definition of the faces](../../../research/pde_ledger_v3/steps/S12_order_source.md) · [Width requirement R-S7-01](../../../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#r-s7-01--the-slab-width-flat-direction-under-a-non-uniform-width)

**Calibration.** Consider two initially flat, equivalent faces in a laterally unbounded, homogeneous slab. Apply outward pressure $p$ on each face over the same 3D ball of radius $L$, with zero load outside and unperturbed width at infinity. A squeeze has $p<0$. Each face receives a force of magnitude $4\pi L^3|p|/3$ in opposite laboratory directions. Let $b=\delta W$ and suppose the net static linear response has constant coefficients $K\ge0$, $D>0$, with $(K-D\nabla^2)b=p$ inside the ball and zero on the right outside. Matching across its edge gives the width change at its centre,

$$
\frac{b(0)}{p}=
\begin{cases}
L^2/(2D), & K=0,\\[2pt]
[1-(1+s)e^{-s}]/K,\quad s=L\sqrt{K/D}, & K>0.
\end{cases}
$$

For two independent interfaces with equal tension $\sigma$ and only area energy, expanding the areas of $Z\pm W/2$ gives $D=\sigma/2$. Doubling the patch radius at fixed pressure then quadruples the central displacement. Positive local stiffness gives the same small-patch trend, followed by saturation at $p/K$ for patches much broader than $\sqrt{D/K}$. A single patch size cannot distinguish the two responses. The single-wall tension has no established value as the active slab's width stiffness. [Independent algebra checks](../../../docs/_measurements/width_pressure_patch.txt) · [Scope of the wall tension and live bulk law, §§8.2 and 8.5](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md)

**Handoff.** On a freely responding active slab, centre both patches at the same exterior location and keep them clear of the opening and distant boundaries. At each size, use the same pressure magnitude on the two faces and both pressure signs, then repeat at half strength. Record the change in their separation and the midpoint separately, with the same intake and far reservoir preparation. Interpret the size law only on scales where the homogeneous approximation applies. Measure any drain-generated width profile with the test pressure removed, since the patch itself produces a $1/r$ exterior response in the tension-only calibration. The physical squeeze test stops at the missing material input. Working result at fd667ea9, pending review.
