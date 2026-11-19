# A pressure-balanced flow guide that gives the electric sign

**One open guide can both emit the bend and receive a repulsive force.** Suppose the mouth receives the net force from a widening stream through a force-only joint. The inlet follows the local sheet normal, the outlet retains an independent bulk direction, and the upstream steering reaction belongs to a separate body. Equal and opposite momentum leaves with the fluid or reaches a separate support. Assume intersecting thrust axes or an independently balanced torque. These conditions describe a possible boundary law. The physical draining throat has yet to supply them. [Mouth attachment](wall_load_to_bend.md) · [Existing interface conditions](../g0_closure_card_v0.md#8-the-only-branch-slot-𝔅)

**Momentum and energy.** Write \(m=A_o/A_i>1\) for the widening ratio, \(v_i=D/(\rho A_i)\), \(v_o=v_i/m\) and \(p_o=p_{\rm amb}\). Bernoulli gives \(p_i-p_{\rm amb}=\rho(v_o^2-v_i^2)/2\). Each thrust includes flow and pressure, \(T=Dv+(p-p_{\rm amb})A\), so \(\Delta T=T_i-T_o=\rho A_i(v_i-v_o)^2/2>0\). The same normal recoil does work on the height \(Z=\ell\langle h\rangle\). A translating mouth also sees the Cartesian component of the inclined inlet thrust. In the sharp limit the combined virtual work gives \(L=\ell\Delta T\) and \(B=-\ell T_o\), for the definitions in the [force report](drain_bend_response.md). Like mouths repel. The pressure head \(\rho v_o^2/2\) above outlet ambient must be supplied upstream. At rest the steady flow exports energy; under translation the difference in pressure-inclusive energy flux equals the stream force times the guide velocity.

**Ownership changes the result.** Join an upstream elbow to the mouth, turning a bulk-aligned stream into the required inlet direction at the same speed and pressure. Its isolated vertical load is zero, so the emitted bend is unchanged. Its lateral reaction restores \(B=+L\). Both assemblies have identical guide intake, pressures and isolated bend. Also, if both cuts turn together, their net force is purely normal and the sharp guide’s tangent response is zero. This is a momentum reaction on whichever part turns the stream; the hole itself carries no imposed spin.

**Finite-width test on one brane.** Retain the earlier annular restoring energy \(k_m\int\eta h^2/2\), replacing its prescribed linear source by the stream recoil. The annulus samples both height and inclination. Write \(z\) for monopole transmission and \(u=\langle\partial_x h_{\rm dipole}\rangle/g\) for its response to a unit incident gradient. Scalar and contact forces together give \(L=\ell z\Delta T\). The additional inlet force gives \(B=L-\ell T_i u\); adding the normal recoil again would duplicate its work.

| Radius | \(z\) | \(u\) | \(B/L\) at \(m=2\) |
|---|---:|---:|---:|
| 1 | 0.919127347040 | 0.985741822475 | −4.362378922 |
| 2 | 0.957190292880 | 0.992653901125 | −4.185248474 |

For any fixed \(m>1\), \(B/L=1-(u/z)(m^2+1)/(m-1)^2\). The ratios are negative and unequal at these radii. Choosing \(D(2)/D(1)=2.77162005175\) to equalize emission for geometrically similar inlets still leaves a **1.04232256434** ratio between the two ordered mixed-pair force magnitudes at \(m=2\).

For a radius-dependent widening, common emission \(L_0>0\) and common repulsion ratio \(\beta=-B/L>0\) require
\[
\frac{m^2+1}{(m-1)^2}=(1+\beta)\frac{z}{u},\qquad
D^2=\frac{2\rho A_iL_0}{\ell z(1-1/m)^2}.
\]
These conditions also suffice within the guide assumptions when the first right-hand side exceeds one. They can be checked using independently obtained outlet areas and flows. Setting those quantities from the desired force would reproduce it by construction.

**Can the known wall pull select the flow?** If a single straight-wall pull \(P=C+4\pi\sigma R^2\) is transmitted in series to the brane and balances \(\Delta T\), it requires \(D^2=2\rho A_iP/(1-1/m)^2\) and upstream pressure head \(P/[A_i(m-1)^2]\). Doubling \(R\) at fixed \(m\) changes this head by a factor between \(1/8\) and \(1/2\). The emitted strength is then \(\ell zP\), larger at the larger radius even for curvature alone under this annular law. This balance specifies the compatible head and rate. Their physical selection and the end receiving the reaction remain open. Working result on base e0a94594. [Wall pull](order_rim_radius.md) · [Guarded calculations](../../../docs/_measurements/drain_bend_response.json)
