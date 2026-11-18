# Two-drain exterior energy at 0PN

**Deliverable:** compare the ideal pair's interaction energy with the attractive load's work, including energy transport and pressure work at both moving cuts.

**Result within the ideal benchmark.** As the drains approach, their exterior interaction kinetic energy **increases** by the same amount as the mechanical work of the attractive load. Both inner regions together must make available twice that work as an interaction increment. The calculation identifies a required energy transfer; it supplies no physical conversion law or free acceleration.

**Preparation and subtraction.** Keep each physical radius $a_i$ and its associated symbolic intake $D_i=D(a_i)>0$ fixed. The two constant-density potential exteriors have signed rates $q_i=-D_i$ in three brane dimensions and $q_i=+D_i$ in four bulk dimensions. Prescribe weak motion along their separation axis, with $d=X_2-X_1>0$ and $\dot d=V_2-V_1$. Retain terms linear in these prescribed velocities; no response dipoles or acceleration law are selected. Use the same exterior domain, excluding two disjoint spheres, for the pair, both isolated fields and the zero-flow background. Thus the interaction energy is

$$
K_\times=\rho\int_{\rm common\ exterior}\mathbf v_1\cdot\mathbf v_2\,d^nx
=\frac{q_1q_2}{\rho\Omega_n(n-2)d^{n-2}},\qquad n=3,4.
$$

Here $\Omega_n$ is the unit-sphere area. Green's identity gives this expression independently of the enclosing radii while the spheres remain disjoint. It is a common-domain cross energy, not a sum of unknown interior energies. Subtracting isolated energies on different domains would additionally require the self-field energies in the omitted regions.

**Both channels and mechanical work.** Summing their interaction energies gives

$$
K_\times=\frac{D_1D_2}{4\pi\rho_3d}
+\frac{D_1D_2}{4\pi^2\rho_4d^2}.
$$

The benchmark's attractive load magnitude is $A=D_1D_2[1/(4\pi\rho_3d^2)+1/(2\pi^2\rho_4d^3)]$. Assigning each whole load to its associated mouth, the mechanical work rate is

$$
P_{\rm mech}=A(V_1-V_2)=\dot K_\times.
$$

Both are positive during approach. This kinetic cross energy therefore does not supply the mechanical gain by decreasing. The bulk term belongs to the ideal four-dimensional spreading geometry; it is not a post-Newtonian correction.

**Where the energy crosses.** For each channel, the moving potentials obey $\partial_t\phi_i=-\mathbf V_i\cdot\mathbf v_i$. Bernoulli pressure gives the cross part of the lab energy current as $\mathbf J_\times=-\rho[(\partial_t\phi_1)\mathbf v_2+(\partial_t\phi_2)\mathbf v_1]$. Subtracting the transported energy of each moving cut gives the outward energy rate from inner region $i$ into the exterior,

$$
J_i=q_i\mathbf V_j\cdot\mathbf v_j(X_i),\qquad
\mathbf F_i=-q_i\mathbf v_j(X_i),\qquad j\ne i.
$$

The surface integrals give $J_1+J_2=\dot K_\times$. Adding mechanical work, the net energy required from each inner region is $J_i+\mathbf F_i\cdot\mathbf V_i=P_{\rm mech}$, with both channels included. Consequently,

$$
\Delta P_{\rm supply}-\Delta P_{\rm other,out}
-\frac{d\Delta E_{\rm inner}}{dt}=2P_{\rm mech}.
$$

This is the increment relative to the matched isolated controls, not twice the total drain power. It could use stored energy or redirect existing conversion power from another output. Those channels must be identified and counted once. No physical supplier is inferred. [Energy obligations](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#6-energy-input-and-forcepower-pairing)

**Checks and limits.** SymPy's Green construction and Mathematica's Gauss-volume construction agree. Direct surface quadratures retain the full companion field at two cut-radius pairs, unequal rates, both flow signs and two motion allocations. Omitting unsteady pressure or moving-cut transport produces a nonzero energy defect. [Scripts, literal readings and guarded outcomes](../../../docs/_measurements/two_drain_energy.json). The [physical mouth response](mouth_conversion_drag.md), conversion power and inertia remain open. Base `fd667ea9291b4430a96afe3fd932926232f5507d`; working, uncommitted pending independent review.
