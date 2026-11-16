# Holding one draining mouth with a companion present

**Working result, November 12, 2026.** Two held, inward-draining mouths produce an inward momentum contribution directed toward the companion in the steady three-dimensional exterior calculation. The change in the force needed to hold the physical mouth still remains undetermined. It depends on the momentum passing through the sleeve and slab faces and on the change in the other stresses around the core.

## Alone and with the companion

Place the target at the origin and its companion at $d\hat x$. Compare stationary states with the same target, ambient material, and drain and return laws. For the exterior calculation, assume constant density $\rho$, inviscid potential flow, and prescribed volume intake rates $Q_1,Q_2>0$. Subscripts 1 and 2 denote target and companion. Rest at infinity and inward volume flux give the velocity potential $Q_i/(4\pi r_i)$ of each point intake, with $\mathbf v=\nabla\sum_i Q_i/(4\pi r_i)$. Integrate the difference of the **full** steady currents $\rho\mathbf v\mathbf v+p\mathbf 1$, with $p=p_\infty-\rho|\mathbf v|^2/2$, on a sphere of radius $b<d$ enclosing only the target. The field of the companion varies across the sphere. [Exact integration](../software/stage1_solver/tools/held_drain_pair.py) · [Independent surface quadrature](../software/stage1_solver/tools/held_drain_pair.wl)

The target alone gives zero axial influx. With the companion, the change in inward current of axial momentum is

$$
\Delta I_{\rm matter,x}=\frac{\rho Q_1Q_2}{4\pi d^2}.
$$

Advection supplies $4/3$ of this value and the pressure change contributes $-1/3$. The integral is independent of $b$ in this exterior. For $d=10$ and $20$ in the same length units, the computed values of $4\pi\Delta I_{\rm matter,x}/(\rho Q_1Q_2)$ were $1/100$ and $1/400$, each checked at $b=1$ and $3$. This recovers the earlier far-field matter term by integrating the spatially varying companion field. A physical use of the point-intake exterior requires a control sphere outside the core in a region where that exterior applies. The native source-to-flux map, full coefficient and other stress contributions remain open. [Run evidence](_measurements/held_drain_pair.txt) · [Earlier matter-stress calculation](../software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md) · [Native slow-limit status](../research/pde_ledger_v3/steps/S14a_2_slow_limit.md)

## What the holder would measure

Take a volume around the target, closed across the sleeve and the remaining slab faces. Let $I_{\rm side,x}$ be the **total** axial momentum entering through its lateral boundary. Let $L_x$ be the outward axial momentum current through the other cuts, including carried momentum and stress. Partition each surface once and include any independent volume forces as $F_{\rm other,x}$, signed toward the companion. In two stationary states the difference in the support force is

$$
\Delta H_x=\Delta L_x-\Delta I_{\rm side,x}-\Delta F_{\rm other,x},
\qquad \Delta=\text{pair}-\text{isolated}.
$$

A negative $\Delta H_x$ means the holder must pull farther away from the companion to keep the mouth still. Matching the computed matter term to $I_{\rm side,x}$ requires the native profiles and the remaining stresses. The signed change in $L_x$ requires the carried momentum and traction on the actual cuts. An isolated scalar intake rate supplies neither vector integral. The same accounting must include any orientation-dependent stress when interpreting a total load as drain attraction. [Native momentum inventory](../research/pde_ledger_v3/directives/S15_decision_list.md) · [Face exchange requirement](../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md)

The face law narrows the next calculation. For instance, consider the upper graph face $w=h$ on a fixed spherical annulus outside the core. If its mechanical outward momentum current is $p\hat n$, its axial projection is $-p\,\partial_x h\,d^3x$. Here $p$ is a signed normal stress amplitude. Write the companion-induced dipole changes in height and amplitude as $\epsilon h_1(r)\cos\theta$ and $\epsilon p_1(r)\cos\theta$ about radial profiles $h_0,p_0$, with the polar axis toward the companion. To first order, the integrated **upper-face stress contribution per radial interval** changes by

$$
\left.\frac{\partial}{\partial\epsilon}\frac{dL^{\rm stress}_{+,x}}{dr}\right|_0
=-\frac{4\pi r^2}{3}\left[p_1 h_0'+p_0\left(h_1'+\frac{2h_1}{r}\right)\right].
$$

Thus even shear-free bulk pressure can matter through the face slope. The lower graph has the opposite outward orientation and its own response. This projection assumes the normal-loading form of the isolated steady model for the comparison. Computing the native holding force requires the companion-induced axial change in **total** current through both faces and the sleeve, including material transport and any core or edge reactions. The available laws leave the load amplitude, surface response and transfer data open. The [Newtonian gravity requirements](../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#what-two-inward-drains-must-reproduce-in-the-newtonian-limit) record that missing momentum exchange for the future model. [Steady face-loading premises](../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md) · [Source excerpts and measured projection](_measurements/held_drain_pair.txt)

Repository base `d544ef73`. Local working calculation pending independent review.
