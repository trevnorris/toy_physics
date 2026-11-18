# Can the exterior bulk close the 0PN holding-force measurement?

The stored one-mouth runs give the carrier momentum entering the opening under opposite sideways flows. The signed slope is 0.2499999853 at speed magnitude 0.02 and 0.249999985299 at 0.01 in solver units, using the boundary-law stress. Negative drives are reflected partners. The receiver packets used for lateral stress have inward number rates 0 and 0.125 and an acceleration drive. The feeding packets have rate 0.25 and steady incident flow. Combining them would mix preparations. [Signed feed readings](../../../docs/_measurements/held_crossflow_budget.txt) · [Receiver protocols](../../../docs/_measurements/open_sleeve_lateral_balance.txt)

**One exterior condition can eliminate the bulk export.** Consider a full four-dimensional bulk exterior governed by the steady neutral carrier equation, with spherical baseline outflow, uniform chemical potential, and no external forcing or distributed exchange in the outer annulus. Let $n_0(R)>0$ be its number density, $s$ its outward speed, $S=2\pi^2R^3$, and $\dot M=Smn_0s$ the conserved mass throughput. For a weak drive along $x$, use $\mu=\hat r\cdot\hat x$. Write the first-order fractional density dipole as $b(R)\mu$ and the velocity dipole as $u(R)\mu\hat r+t(R)(\hat x-\mu\hat r)$. With $L=\log(n_0/n_{\rm ref})$ for constant $n_{\rm ref}>0$, primes denoting radial derivatives and $k=\hbar^2/(4m)$, the outward momentum rate in the drive direction is

$$
\delta\Phi_{{\rm bulk},x}=
\frac{\dot M}{4}(u+3t+s b)
+\frac{S k n_0}{4}(L'b'-bL'').
$$

The ordinary pressure and density-gradient stress each contain the curvature $b''$ of the sideways density change. It cancels in their sum. SymPy projected the logarithmic-density stress on the sphere. Mathematica differentiated the wave amplitude in four Cartesian directions. Both retained the baseline density gradients and obtained the same rate. [Matter stress and EOS](pathA_21c_force_from_noether_stress_tensor.md) · [Both calculations](../../../docs/_measurements/mouth_bulk_far_flux.json)

For a finite positive far density, bounded $b$, and $u,t\to0$, the first term vanishes at infinity. If also $R^3(L'b'-bL'')\to0$, the remote bulk carries **zero net first-order sideways momentum**, even with a local response near the mouth. On a matched steady, source-free annulus the integrated rate is the same at each radius. A finite box with a fixed bath value does not establish the derivative condition. A partial bulk surface also requires the fluxes on its remaining boundaries.

**Physical status.** The proposed force test still lacks a matched brane-material and bulk solution or a direct holder reaction. On such a solution, integrate the actual flux over two enclosing surfaces, include all other cuts and separately applied forces, and verify the same signed balance at both speeds. No throat interior is needed for that measurement. The scalar exterior result cannot by itself identify the receiver of the physical converting mouth. Working on base fd667ea9, uncommitted pending review. [Physical bulk and support inputs](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L149) · [Towing comparison](mouth_tow_and_companion.md)
