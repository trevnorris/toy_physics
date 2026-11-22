# Can the same material select a steadily draining depth and flare at two radii?

**The stationary power balance rules out a positive-intake equilibrium at either calibration radius for the passive conversion law and the Step 4 reservoirs.** The obstruction survives arbitrary changes to depth and flare. Radius 3 remains reserved. The mechanism that seals intact faces remains parked. [Two-radius energy audit](../../../research/pde_ledger_v3/_measurements/mouth_free_shape_power.json) · [Previous radius comparison](mouth_radius_pivot.md)

### What was allowed to change

We chose radii 1 and 2 for calibration and radius 3 for a blind prediction of depth, profile, intake and inertia. Both calibration mouths see the same two baths of unit carrier density, a fixed spatial phase difference of 0.4510084282857874, sealed intact faces and an open melting exit. The reservoir radius is 28, roof height 22 and slab thickness 1. The reference solutions have positive intakes of 0.249617766 and 1.371519070 respectively. These are measured on held shapes. The test for a free shape below permits the intake to change by **any** amount as the shape relaxes. Its only throughput requirement is a strictly positive steady current from ordered feed to disordered bulk. The zero-current state cannot generate the drain's exterior inverse-square flow.

### The obstruction before any shape fit

The material law being tested is finite-mobility passive conversion, \(D_t\chi=-M_\chi\mu_{\rm eff}\), where the order fraction goes from one in the feed to zero in the bulk. On a stationary shape its energy loss is

\[
P_\chi=\int (\mathbf u\cdot\nabla\chi)^2/M_\chi\,d^4x .
\]

Choose the smooth order marker \(H(\chi)=3\chi^2-2\chi^3\). Let \(n\) be number density and \(D>0\) the through-current in units where constituent mass is one. Mass conservation, the two pure-phase end values of \(H\) and sealed remaining boundaries give \(D=-\int n H'(\chi)\mathbf u\cdot\nabla\chi\,d^4x\). By Cauchy–Schwarz,

\[
 P_\chi\,\mathcal I\ \ge D^2,\qquad
 \mathcal I=\int M_\chi n^2[H'(\chi)]^2\,d^4x .
\]

This bound applies to finite, regular conversion layers, with finite positive mobility where material converts. The integral is taken over the mobile material; locked faces have no conversion. At the two *reference* throughputs the required products are at least 0.0623090 and 1.88106 in solver units. For any other positive intake, spreading conversion across a larger finite flare can reduce the energy required per constituent, but it leaves a strictly positive loss. On a separate saved radius-1 order profile at intake 0.25 the lower bound is 0.0124923 power, against 0.0128076 directly integrated. That one-mouth cost does not fix the cost at radius 2. [Direct passive-loss integral](mouth_coupled_conversion.md)

The coherent port energy current is \(-\operatorname{Re}(\dot\psi^*\nabla\psi)\). The prescribed phases are constant in time, so the native port power is exactly zero. A common temporal phase frequency gives equal energy per constituent at the two ends and zero net power at matched current. The archived port currents agree within \(1.9\times10^{-10}\) at both training radii. The two pure order phases have equal energy in the neutral material law. A stored wave can exert pressure on a stationary surface, although sustained power drawn from its finite energy would change the prepared occupation. Energy conservation therefore excludes a **steady, passive, draining** shape under these conditions. This is a necessary test before solving for a stationary depth. Neither the equilibrium profile nor its free acceleration was measured. [Carrier energy current](mouth_stress_energy.md) · [Reviewable calculation](../tools/mouth_free_shape_power.py)

A positive control distinguishes this obstruction from an inability of the solver to carry current through a curved front. In a separate pressure-driven annulus the finer run supplied 0.00462897 energy per constituent. Its two readings, reservoir power and conversion loss, were both \(8.09684\times10^{-5}\) per inner area. That annulus has different end conditions and its own normalization. A net energy drop can fund conversion. Establishing a stable full flare still requires two-sided traction, wave occupation and force balance including both ends. [Powered curved interface and release](mouth_coupled_conversion.md)

### What the held mouths say about the required geometry

Sixteen signed shape changes at the two calibration radii test the near-lockstep scaling in the slender reference mouths. At fixed radius, a local 1% increase in axial depth lowers intake per unit carrier inertia by about 2.04% at radius 1 and 2.08% at radius 2. Widening the middle of the flare by 1% at fixed depth and fixed ends lowers the ratio by about 1.33% and 1.43%. In the latter test intake increases, but inertia increases faster. These measured slopes support the proposed depth sign near the tested shapes. They leave open the combined changes in width, depth, occupation and mass of a material-selected shape. [Signed steps, intake, inertia and geometric work](mouth_selected_shape_sensitivity.md)

### What would decide between a family and one elementary radius?

A radius held during a calculation can carry an unmeasured radial reaction. To identify a *free* radius, release that reaction as well as those at the axial ends. Let \(F(y;\lambda)=0\) denote the full stationary balance including depth, radius and profile, with \(\lambda\) specifying the physical preparation. A differentiable family at identical \(\lambda\) obeys \(F_y\,dy/dR=0\). After removing overall translations and the arbitrary carrier phase, a tangent that changes the radius is a neutral deformation. A continuously adjustable family with positive restoring stiffness in the radius direction therefore needs a varying physical label such as conserved trapped-wave action. Alternatively there can be multiple isolated stable roots at the same preparation. Finding neither a powered branch nor a free radial root here gives no evidence for a unique or quantized elementary radius, or for clusters of elementary mouths.

**Next independent input.** Determine how much specific energy an incoming constituent can release as it loses order into the actual receiving phase, at the intended common pressure and material preparation. Include any prepared strain energy in the ordered feed and any continuing boundary work. The present equal-well model supplies zero. A measured positive energy budget, or an independently derived reversible conversion law, would reopen the free-shape calculation. Use the same law and independently specified wave preparation at radii 1 and 2, test local and endpoint forces and stability, then freeze the radius-3 predictions before solving there. Capture tracking intake to 0.21% in a held rigid family satisfies a useful reciprocity check at that accuracy. Establishing action–reaction for unequal released mouths still requires their forces and the momentum left in the medium. [Reservoir and material controls](mouth_face_conversion_lock.md) · [Free-fall roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

The independent reviewer concluded “Scoped clear for the passive-conversion energy exclusion at fixed-phase reservoirs, the conditional depth and flare sensitivities at radii 1 and 2, and preserving the radius-3 prediction.” [Guarded audit](../../../research/pde_ledger_v3/_measurements/mouth_free_shape_power.json)
