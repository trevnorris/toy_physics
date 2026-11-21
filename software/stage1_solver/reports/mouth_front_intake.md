# Can the conversion layer select a resting mouth’s intake?

**Working result.** A weak local head drop selected inward mass flow through a stationary, attached layer in two conditional 0PN models. Pressure supplied the dissipated energy. The existing brane equations leave the layer’s attachment, relaxation rate and connection to the receiving reservoir open. [Measurements and independent reviews](../../../research/pde_ledger_v3/_measurements/mouth_front_intake.json)

### Pressure-controlled layer

Let \(\chi\) fall from one in ordered brane to zero toward the bulk. We tested the earlier relaxation hypothesis \(D_t\chi=-M\mu_\chi\), with material force density \(\mu_\chi\nabla\chi\), constant mass density \(\rho\) and no other conversion drive. Its equal-depth order energy alone gives zero steady crossing for a free planar front between homogeneous phases. As one possible attachment we added a conservative spring localized across the layer and tied its preferred profile to the stationary rim. This allows the front to shift under load while material passes through. Pressure here is the incompressible mechanical load. The parent carrier would need its full specific head on the same two sides. [Order energy and mobility](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md) · [Earlier steady balance](../../../docs/_measurements/rim_curvature_head.txt)

For area \(A\), normal speed \(u\), intake \(D=\rho A u\) and local head drop \(H=\Delta p/\rho\), the pressure-controlled solution obeys

\[
H=\frac{u}{\rho M}\int(\partial_z\chi)^2 dz,\qquad P_{\rm loss}=D H .
\]

The stationary anchor took the normal load and did zero work. Together with the equilibrium tension, density and area, two decay rates of weak order disturbances in constrained resting material predicted a zero-head conductance \(dD/dH=14.4\) in trial units. The boundary solver took pressure as input. At \(H=0.002\) and \(0.001\), its respective \(D/H\) values were **14.400669** and **14.400167**.

At the positive baseline \(H=0.002\), independent pressure steps gave \(dD/dH=14.402008\). If the [ideal exterior](mouth_stream_power.md) transmits its entire head change \(U^2/2-UV\) to this layer, the predicted four-reversal coefficient is \(-14.402008\). New local-head solves gave \(-14.402116\) at full strength and \(-14.402033\) at half strength. A fixed local head gives zero for the same comparison. These controls test the intake law. The transfer of head from actual reservoirs to the mouth still needs independent data.

### An observable at the rim

In the second model the conversion surface spans a three-ball of radius \(R\) and attaches only around its edge, where its depth displacement vanishes. Here the projected area is \(A=4\pi R^3/3\). For a thin layer with small slopes, uniform material properties and zero sideways flow, a stationary depth profile \(h\) satisfies \(u=-\nu\nabla^2 h\), with \(\nu=M\kappa\) and equilibrium tension \(\sigma=\kappa\int(\partial_z\chi_0)^2 dz\). The curvature relaxation is the same mathematical mechanism studied for a material order boundary. [Allen and Cahn](https://doi.org/10.1016/0001-6160(79)90196-2)

For an axisymmetric profile, the edge slope gives the total mass rate

\[
D=-4\pi\rho\nu R^2 h'(R).
\]

A fixed-edge decay time \(\tau_1=R^2/(\pi^2\nu)\) determines \(\nu\) independently if material motion during decay is suppressed or accounted for. Uniform pressure also permits a center-depth estimate. For the archived intake shape proportional to \((1-r^2/R^2)^3\), that shortcut returns **\(187/64\)** times the true rate in the continuum. Coarse and fine saved profiles gave factors of 2.8753 and 2.8959. The saved carrier current was prescribed. These shape predictions are available for comparison with a future independent depth reading.

**Gravity test.** Measure the rim slope, relaxation time and actual inlet-to-outlet head at the same mouth. A uniform normal load on a bowed fixed-edge surface integrates to zero sideways force. A small prescribed load gradient \(\partial_x\Delta p=g\) over the uniformly loaded starting cap gives \(F_x=(2/5)g A h(0)\) to first order in \(g\). Whether a companion creates this gradient, and how its reaction reaches a free center, remain open. Working base 49dca23e. The pressure solves and profile calculation passed. The independent edge review is analytic. Its Wolfram checker hit the guard’s task limit. Formal review is pending.
