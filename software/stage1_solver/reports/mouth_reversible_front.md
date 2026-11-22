# Can ordered brane convert to coherent bulk reversibly?

**A Hamiltonian material law supports a stationary conversion layer with current and zero net steady power. Its planar interface can translate freely in a homogeneous medium.** This gives a candidate for the energetic part of drainage. The calculation resolves the fields on a neutral patch with coherent material on both sides. The same law applies to a patch at either calibration mouth, but a plane has no mouth radius, geometric rim or trapped wave. [Guarded measurements](../../../research/pde_ledger_v3/_measurements/mouth_reversible_front.json) · [Patch calculation](../tools/mouth_reversible_front.py)

We retained the coherent carrier energy \(U(n)=(n^5-5n+4)/20\), writing \(f=\sqrt n\), and the interface energy \(W(\chi)=\chi^2(1-\chi)^2\). The tested extension gives the crystalline order its own reversible momentum density \(\pi\). The trial Hamiltonian density is

\[
\mathcal H=\tfrac n2|\nabla\theta-\pi\nabla\chi/n|^2+
\frac{\pi^2}{2I}+\tfrac12|\nabla f|^2+U(n)
+\tfrac\kappa2|\nabla\chi|^2+W(\chi).
\]

The canonical time part of the action is \(-n\dot\theta+\pi\dot\chi\). Variation gives \(\mathbf j=n\mathbf u=n\nabla\theta-\pi\nabla\chi\), \(\pi=I(\partial_t+\mathbf u\cdot\nabla)\chi\), and \(\partial_t\pi+\nabla\cdot(\pi\mathbf u)=\kappa\nabla^2\chi-W'\). Thus order changes as material moves across the sheet and exerts a reciprocal force on the carrier. Positive \(I=1\) and \(\kappa=2\) are common patch parameters in interface-width units. An order-amplitude oscillation on a uniform patch would fix \(I\) from \(I\omega^2=W''+\kappa k^2\) in this model. The ordered and bulk phases retain their respective shear properties; shear strain is zero in this longitudinal test.

In the stationary plane the phases at the two density-one reservoirs differ by 0.4510084282857874, with both stationary in time. Both density and order respond to the unknown current \(j\). The two pure order values are set far from the interface. Reflection about the midpoint selects a centered representative of the sheet for the field solve; all perturbations of the sheet are allowed in the subsequent stiffness calculation. The order stress at the outer ends is \(5.06\times10^{-14}\) in the 32-width sample and \(5.70\times10^{-21}\) in the 48-width sample. The separated phases have equal energy density. In this stationary state the Noether energy current is \(\mu j\), where \(\dot\theta=-\mu\), and the reservoir contributions cancel. There is no dissipative order term in the trial Hamiltonian.

| Preparation | Current per area | Integral of carrier velocity | Extra phase from order |
|---|---:|---:|---:|
| 32 widths, positive drive | +0.014019626 | +0.448671214 | +0.002337215 |
| 32 widths, reversed drive | −0.014019626 | −0.448671214 | −0.002337215 |
| 32 widths, zero drive | 0 | 0 | 0 |
| 32 widths, zero order inertia limit | +0.014092657 | +0.451008428 | 0 |
| 48 widths, positive drive | +0.009363090 | +0.449447732 | +0.001560697 |

The last two columns sum to the imposed phase difference. Normal momentum flux across the positive-drive patch varies by less than \(1.6\times10^{-13}\) in the conserved-stress calculation. The local carrier-pressure, quantum-stress and order-stress sum agrees within \(5.9\times10^{-11}\) on the 32-width sheet. A phase delay through the sheet after subtracting \(\int u\,dx\), together with an independent loss measurement, can test the proposed order inertia.

For the longitudinal stability test, minimize the positive flow and momentum energy at fixed end phases. The phase energy becomes \((\Delta\theta)^2/(2C)\), with \(C=\int(n^{-1}+I\chi'^2/n^2)dx\). The second variation includes the nonlocal positive term \(j^2(\delta C)^2/C\). On 160, 320 and 640 linear elements the lowest static eigenvalue is 0.00095109, 0.00023803 and 0.00005952. Its mode agrees with translation of the order sheet with overlap exceeding 0.99999998. The next eigenvalue approaches 1.50012 on the finest mesh. This is a fixed-port energy stiffness, in a joint density-amplitude and order norm. The translational freedom of the infinite homogeneous problem follows from its symmetry; the finite ends weakly constrain translation. No negative direction was resolved in these planar static tests. A curved mouth also requires the radius and depth of the conversion surface, the trapped shear-wave occupation, two-sided stresses and the forces where the surface joins the slab and exits into the bulk. A whole-mouth intake or translational inertia cannot be calculated by multiplying a patch reading by a chosen area. Radius 3 was reserved.

The reviewer concluded “Scoped clear for the traction-balanced planar Hamiltonian front, its reversible current and unresolved translation.” An independent full-interval minimization allowed the center to move and recovered current 0.01401963053 with 512 cells. Starting 1.75 widths away retained that separation with an energy difference below \(8\times10^{-13}\) in the tested domain.
