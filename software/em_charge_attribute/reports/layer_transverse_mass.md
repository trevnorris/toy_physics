# Which mass sets the long-wavelength reflection phase?

**A single mass per unit area sets the leading phase for a lossless layer with a shear-traction-free bulk edge.** Measured from the entrance to the layer, the apparent reflection depth tends to \(\Sigma_\perp/\rho_0\). Here \(\Sigma_\perp=\int_a^b\rho_\perp(q)\,dq\) is the mass per unit cross-sectional three-area that participates in the transverse motion, and \(\rho_0\) is the transverse inertia per four-volume on the fully ordered side. The current model leaves the physical value of \(\Sigma_\perp\) open. [Earlier reflection test](conversion_layer_reflection.md) · [Live material inertia](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L129)

**Derivation.** Hold radius and intake fixed in the same planar, gapless shear benchmark. Neglect through-flow effects on propagation and loss or exchange during a wave cycle. The rim is \(q=0\), the varying layer occupies \(a<q<b\), and the uniform incident medium has speed \(c_0\). Write \(T=\mu(q)U'\) for the shear stress amplitude and normalize \(U(a)=1\). The wave equation gives \(T'=-\omega^2\rho_\perp U\). Since \(T(b)=0\), \(T(a)=\omega^2\int_a^b\rho_\perp U\,dq\). The zero-frequency displacement is uniform. A sufficient regularity condition, allowing stiffness to vanish at the edge, is finite \(\int_a^b M(q)/\mu(q)\,dq\), where \(M(q)=\int_q^b\rho_\perp(s)\,ds\). The integral equation for \(U\) then gives \(U=1+O(\omega^2)\) uniformly.

With \(k=\omega/c_0\), \(\mu_0=\rho_0c_0^2\) and time convention \(e^{-i\omega t}\), the continuous rim reflection phase gives
\[
\frac{\arg R}{2k}
=a+\frac{\Sigma_\perp}{\rho_0}
+k^2\left[\frac{c_0^2}{\rho_0}\int_a^b\frac{M(q)^2}{\mu(q)}\,dq
-\frac{1}{3}\left(\frac{\Sigma_\perp}{\rho_0}\right)^3\right]+O(k^4).
\]
Both terms in the bracket follow by expanding the entrance stress and matching to incident and reflected waves. A separate mass per area attached at the bulk edge adds its mass to the leading readout. Application to a draining mouth requires checking its end load and that wavelengths long compared with the layer coexist with an elastic response on that timescale.

**Check at equal mass.** Let \(x=(b-q)/h\), \(h=b-a\), and \(\mu=\mu_0 x\). Two positive interior inertia profiles \(\rho_\perp/\rho_0=x\) and \(x^2(3-2x)\) both give \(\Sigma_\perp=\rho_0h/2\). Their leading reflection depths are exactly \(h/2\). Their first corrections are \(k^2h^3/48\) and \(3k^2h^3/224\), respectively. SymPy expanded the displacement and integrated the distributed momentum. Mathematica independently expanded the ratio of stress to displacement and checked both profiles at three nonzero wavelengths using a regular series. Its initial differential-equation solver failed near the free edge; its extrapolated values were discarded. [Commands and outputs](../../../docs/_measurements/layer_transverse_mass.json)

**Independent calibration.** The carrier kinetic term uses the total number density \(n\) and carrier velocity. The order balance counts \(\chi_B n\). Neither identifies the velocity of the transversely displaced material throughout a converting layer. The homogeneous brane kinetic term and the rule for momentum carried out with converted material leave the stored transverse momentum response open. [Order and kinetic terms](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md#L93) · [Homogeneous kinetic term](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L376) · [Carried momentum](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L217)

For the same inward drain and fixed cuts, compare opposite weak sideways drives. Independently track the nearly uniform transverse material velocity and the change in stored transverse momentum of that patch. Their odd-response ratio \((P_{\perp,+}-P_{\perp,-})/[A(v_+-v_-)]\) gives \(\Sigma_\perp\) in this regime. If momentum is inferred from impulse, include all applied and boundary tractions and subtract signed net advective export over the same time interval. Match or remove any end-support inertia. The existing carrier records supply no such paired reading. [Available one-mouth fields](rim_bend_readouts.md)

Working result on base e0a94594, uncommitted pending independent review.
