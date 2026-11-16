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

## Do unequal mouths exert equal and opposite 0PN forces?

Put mouths 1 and 2 at rest a distance $d$ apart, with $+x$ pointing from 1 to 2. Each has its own radius $a_i$ and the same physical radius-to-intake law $D(a_i)$. For a reading on mouth $i$, subtract its reading alone at the same position. For a whole-region reading $Y$, use $\delta Y=Y_{12}-Y_1-Y_2+Y_0$, where $0$ is the common bath without either mouth. This counts the background once and keeps all signs in the same frame.

**Momentum crossing between the two mouths.** In the steady constant-density, potential-flow exterior above, set $\mathbf v=\mathbf v_1+\mathbf v_2$. Subtracting the two isolated advective and pressure currents from the paired one gives

$$
 \delta\Pi_{xx}=\rho\bigl(2v_{1x}v_{2x}-\mathbf v_1\cdot\mathbf v_2\bigr).
$$

Put a plane between the mouths at distances $\ell$ and $b$ from the two centers, with $\ell+b=d$, outside both openings. For cylindrical radius $s$ on that plane, the integral of $\delta\Pi_{xx}\,2\pi s\,ds$ from zero to infinity is $-\rho Q_1Q_2/(4\pi d^2)$ for any such plane. At fixed $d$ the interaction stress at a distant closing surface of radius $R$ is $O(R^{-4})$, so its integrated flux vanishes as $R\to\infty$. The steady exterior stores no changing momentum. The currents **into** the left and right regions are therefore $+F_*$ and $-F_*$, with $F_*=\rho Q_1Q_2/(4\pi d^2)>0$. Each points toward the companion. Rotational symmetry about the line of centers cancels the transverse components on integration. The integral at the midpoint was also evaluated directly; the result agrees with the earlier sphere integral around a single drain. Under steady flux matching $Q_i=D(a_i)/\rho$, so $F_*=D(a_1)D(a_2)/(4\pi\rho d^2)$ for any positive $D$. [Exact plane integrals](_measurements/unequal_mouth_pair_balance.txt) · [Single-drain control surface](_measurements/held_drain_pair.txt)

**Forces on the physical openings.** Enclose both mouths, their assigned bound material, and the intervening brane. Denote the total outward $x$-momentum flux on its external cuts by $\Phi_x$, including carriage and traction at the far boundary and on any open slab or bulk cuts. A shared internal cut cancels. Let $B_x$ be any independently applied volume force other than the holders, with face exchange counted only once. In stationary held preparations the intrinsic drive on each mouth is the negative of the paired-minus-alone holding force. For a free pair take the corresponding mouth momentum rates, and let $P_{\rm rest,x}$ be the momentum in the region excluding the two translating bodies. The joint momentum balance gives, with support forces removed from the reported mouth forces,

$$
 F_{1x}+F_{2x}
 =\delta B_x-\delta\Phi_x-\delta\dot P_{\rm rest,x}.
$$

Storage vanishes for genuinely stationary held fields. For a free pair the storage rate and the cut fluxes must be evaluated in the released states. Equal and opposite forces at inverse-square order require the right side to be $o(d^{-2})$. For a held pair this can be checked directly by taking $d^2(\Delta H_{1x}+\Delta H_{2x})$ at large separation. Its limit must vanish while $-d^2\Delta H_{1x}$ tends to a positive value. These are total mechanical support forces. The native balance retains open stress and face-exchange terms. Its enclosing surfaces alone do not fix this coefficient. [Moving-region balance](../research/pde_ledger_v3/steps/S16_moving_region.md#L38) · [Physical momentum partition](drain_driven_initial_acceleration.md#momentum-arriving-and-momentum-leaving)

In the flux-matched, single-monopole limit, write the two zero-acceleration drives as $e_1F_*$ and $-e_2F_*$ at leading order. The unknown intakes factor out, so equal and opposite attraction requires $e_1=e_2>0$. As a conservation check, suppose mouth 1 retains the entire incident sideways momentum and mouth 2 passes its entire share on, with no storage or other applied force. In units of $F_*$ the two drives are $+1$ and $0$, and the exported $x$-momentum rate is $-1$. Total momentum balances even though the two drives differ. This example assigns the transfers for illustration. The remaining physical requirement is the **combined signed outward flux plus storage rate, less independent applied forces, at order $d^{-2}$**. Its measurement or derivation for unequal mouths tests reciprocity without fitting $D(a)$ or specifying an inertial mass. [Signed export check](_measurements/unequal_mouth_pair_balance.txt) · [Independent inertia requirement](../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#how-drain-rate-and-retained-momentum-must-scale-for-free-fall)

## What can be read for the two held unequal mouths?

The saved stationary fields cover a radius-1 mouth on the left at $x=-d$ and a radius-1.5 mouth on the right at $x=0$. Positive $x$ points right. For each opening I integrated the inward current times the local $x$ velocity and the mixed stress across its actual footprint. I subtracted the isolated reading for that same opening at the same position, on the same grid, with the same boundary law and reservoir. The left mouth used a prescribed rate of $0.25$. The right used a speed profile dilated from the unit-mouth calibration. This comparison reuses those fields and assumes no new physical intake law. [Paired preparation](../software/stage1_solver/tools/lateral_paired_acceleration.py#L115) · [Isolated off-center preparation](../software/stage1_solver/tools/lateral_target_rate.py#L126) · [Isolated centered preparation](../software/stage1_solver/tools/lateral_radius_momentum.py#L103)

| Separation | Extra load on left | Extra load on right | Signed sum | Sum as % of mean load magnitude |
| ---: | ---: | ---: | ---: | ---: |
| 4 | $+0.00105122206413$ | $-0.00105086764686$ | $+3.54417266901\times10^{-7}$ | $+0.0337204702$% |
| 8 | $+0.000261909823373$ | $-0.000261178004586$ | $+7.31818787465\times10^{-7}$ | $+0.279807232$% |

These are **brane-side** loads in solver force units at reservoir extent 48, with four lateral cells per unit length and 16 intervals through the slab. The signed sum uses the same laboratory axis for both mouths. The denominator for the percentage is $(|L_1|+|L_2|)/2$ using the extra loads. At separation 8 the off-center mouth alone had load $+2.36198917227\times10^{-7}$. A second reconstruction of the full momentum current using normal field differences on the same solutions gave fractional sums of $+0.000990501193$% and $+0.245252843$% at separations 4 and 8. At separation 8 and the same spacing, moving the reservoir from extent 24 to 48 reduced the primary fractional sum from $+2.46463237$% to $+0.279807232$%. On both apertures the additional sideways load points toward the companion. The dependence on reservoir position and stress readout leaves the asymptotic coefficient of their sum undetermined. [Full readout including each control and both stress evaluations](_measurements/unequal_mouth_holding_readout.txt)

**Physical holding forces.** Neither packet records an external support reaction. Both centers are solver inputs. For a stationary mouth, let $E_i$ be the paired-minus-alone rate of $x$ momentum passed on from the receiving side of the cut, including traction and any other applied-force correction with the sign in the balance above. Then $\Delta H_i=-L_i+E_i$. Full retention would give outward holding forces equal to the negatives of the two tabulated loads. If both mouths instead pass on their entire load, both additional holding forces vanish. The same stored feeding-side fields admit either accounting until the mouth law fixes $E_i$. Determining the two physical reactions requires the actual sideways constraint force on each opening or an independently specified net momentum-transfer law for each radius. A measurement of just the sum cannot distinguish equal nonzero attraction from zero force on both mouths. [One-body force balance](em_u1_body_definition.md#L56) · [Force readout and boundary terms](drain_driven_initial_acceleration.md#what-would-turn-the-measured-held-force-into-a-release-acceleration)
