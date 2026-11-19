# What axial momentum does the receiving bulk carry away from one rounded end?

**A steady spreading receiver carries a finite axial momentum rate near the opening. In the open-bath benchmark that rate goes to zero at large distance.** The local value changes with the discharge profile even at the same mouth radius \(a\), mass intake \(D(a)\), and far-bulk density. Including the neutral carrier's density-gradient stress leaves a second independent boundary datum. The physical conversion law and force delivered to the slab remain open. [Previous common-displacement test](common_face_work.md) · [Guarded calculations and inputs](../../../docs/_measurements/rounded_end_bulk_flux.json)

**Receiving-side calculation.** Use a four-dimensional, steady, irrotational, incompressible receiver of density \(\rho_4\), at rest at infinity. Enclose the rounded end by an inner hemispherical cut of radius \(a\); the plane outside its equator reflects the flow. This is an ideal exterior geometry. Assume all the intake reaches this receiver, with no further material sources outside the cut. The cut prescribes a discharge profile for each test; it does not impose a pressure or identify the actual cap as permeable. For an upper throat put \(\mu=\cos\theta\) toward \(+w\). Choose the regular outward speeds \(D(a)[1+\epsilon(4\mu^2-1)]/(\pi^2\rho_4 a^3)\) at the inner cut. Both \(\epsilon=0\) and \(1/4\) are positive everywhere and have exactly the same total throughput.

On the outer hemisphere of radius \(R\ge a\), solve Laplace's equation and integrate carriage plus Bernoulli pressure relative to the resting bath. For a mirrored throat the sign \(s\) reverses in fixed lab coordinates. The axial momentum rate is

\[
\Phi_w(R)=s\,\frac{2D(a)^2}{\pi^3\rho_4 R^3}
\left[\frac13-\frac{2\epsilon}{5}\frac{a^2}{R^2}
+\frac{\epsilon^2}{7}\frac{a^4}{R^4}\right].
\]

At \(R=a\) the shaped profile gives **\(407/560\)** of the uniform profile's rate. At any \(R\), the difference \(\Phi_w(a)-\Phi_w(R)\) is exactly the signed pressure transfer to the intervening plane. The distant hemisphere exports zero axial momentum although the mass throughput stays \(D(a)\). SymPy integrated the curved surface. Mathematica checked the bulk equation and recovered the same result from the exterior plane. For a self-similar discharge at fixed \(R/a\) and \(\rho_4\), a nonzero radius-independent **finite-cut** rate would require \(D(a)\propto a^{3/2}\). No such law was imposed or derived. [SymPy](../rounded_end_bulk_flux.py) · [Mathematica](../rounded_end_bulk_flux.wl)

**Density and stress.** For the stationary neutral carrier, take weak radial outflow and a density perturbation \(\eta=n-n_0\) of second order in the speed. Let \(m\) be constituent mass and \(k=\hbar^2/(4m)\). Continuity, the stationary phase equation, and the full carrier stress give the upper-hemisphere result through that order,
\[
\Phi_w(R)=\frac{2D(a)^2}{3\pi^3 m n_0R^3}+4\pi kR^2\eta'(R).
\]
Two exterior solutions can differ by \(\delta\eta=A K_1(R/\xi)/R\), where \(K_1\) is the decaying modified Bessel function and \(\xi^2=k/[n_0 h'(n_0)]>0\) for enthalpy per constituent \(h\). Both approach the same far density at the same throughput; \(A\) is of the perturbative order. Their finite-cut fluxes differ by \(-4\pi k A R K_2(R/\xi)/\xi\). Even radial symmetry leaves the inner density or its normal derivative to be selected by matching. [Carrier tensor and momentum equation](../../stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md)

**Boundary relation to carry forward.** On the stationary physical conversion layer, use the same normal \(\hat n\) from lumen to bulk on both sides. The axial momentum delivered by the carrier is \(\mathbf e_w\cdot(\Pi_{\rm lumen}-\Pi_{\rm bulk})\hat n\). The missing law must connect the local material transfer and density traces to that signed traction and specify how much reaches the slab after layer storage, stress transported along the layer to its rim, and other support reactions are accounted for. Each channel enters once. The permeable E3 condition sets the carrier traction jump to zero. The impermeable E2 condition excludes normal carrier flow through the cap and allows a normal reaction. Neither identifies the physical active conversion support. Until a matched law or independent load measurement is available, the calculation stops here. [Contact alternatives](../g0_closure_card_v0.md#L431) · [Physical matching balance](../../../docs/drain_driven_initial_acceleration.md#matching-a-conversion-layer-to-the-sleeve-and-receiving-bulk-at-0pn)

Working result on base e0a94594, pending independent review.
