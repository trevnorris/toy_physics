# Does the first change of flare shape bring it into equilibrium?

**The flowing carrier amplifies one small reshaping of the flare. Positive wall tension and the occupied ground wave both resist that reshaping. None of the five sampled profiles balance the radial load at all eight measured depths.** The result gives a way to distinguish support laws on an isolated mouth, before a gravity or unequal-radius measurement.

### What was allowed to move

The radius-one mouth drained at prescribed intake 0.25. We held the entrance and exit geometry fixed and changed the radius in five of eight measured interior bands, opening some while narrowing others. The signed coordinate $q$ gives a maximum local change of $|q|$ mouth radii. We selected the pattern from the archived cubic shape so that at the reference mesh the first-order work from constant tension, affine backing pressure and one trapped ground mode cancels. This isolates a residual carrier force. The deformed states include all of these load terms. We re-solved the coherent carrier and the trapped ground mode at $q=0,\pm0.0001,\pm0.0002$ on two resolutions. [Test and geometry](../tools/mouth_free_radial_direction.py) · [Guarded measurements](../../../research/pde_ledger_v3/_measurements/mouth_free_radial_direction.json)

At the undeformed shape the carrier exerted $+3.84715\times10^{-5}$ units of force along this pattern. A stress integral and the virtual work of the fluid agreed. A separate difference of the fixed-intake energy included the work of the reservoirs. On the finer carrier mesh the half-strength pair gave the following slope of net force with deformation.

\[
\frac{dF_q}{dq}\simeq
+0.10097 -1910.04\,\sigma -0.05315\,J
-1.95169\,P_0-2.96220\,g .
\]

Here $\sigma$ is uniform mechanical surface tension, $J=E/\omega$ is trapped-wave action, and $P_0+gw$ is the assumed backing pressure. Positive slope drives a further displacement in the same direction. The wave coefficient was −0.04702 on the coarser wave grid, and an independent finer wave calculation gave −0.05601. Its restoring sign persisted. The near-zero wave work at the reference is a property of the chosen discrete projection, with measurable leakage on refinement. A negative total slope indicates a restoring contribution around the trial shape. Stationarity also requires the force itself to vanish in every allowed direction. [Independent derivatives and work checks](../../../research/pde_ledger_v3/_measurements/mouth_free_radial_band_check.json)

### What eight separate force measurements ruled out

At each shape we let tension and action take any nonnegative values and independently varied the offset and slope of the backing pressure. The following are the smallest achievable values of the largest pressure imbalance among the eight bands, in model pressure units, from carrier volume work on mesh 34.

| Largest signed deformation $q$ | Minimum worst pressure imbalance |
|---:|---:|
| −0.0002 | $9.11924\times10^{-5}$ |
| −0.0001 | $9.47516\times10^{-5}$ |
| 0 | $9.88188\times10^{-5}$ |
| +0.0001 | $9.32101\times10^{-5}$ |
| +0.0002 | $9.30164\times10^{-5}$ |

Each bound has an independent signed combination of band measurements that no permitted choice of coefficients can cancel. For comparison, the largest wall-stress versus volume-work difference in any band across these five fine-mesh shapes was below $5.5\times10^{-9}$. The independent coarser mesh gave the same exclusion. The historical smaller bound on the cubic allowed negative tension and negative wave action. The present bounds apply to the five sampled profiles and this class of support laws. [Pressure calculation](../tools/mouth_free_radial_band_check.py) · [Exact-decimal checks and source data](../../../research/pde_ledger_v3/_measurements/mouth_free_radial_band_check.json) · [Earlier cubic result](mouth_free_equilibrium_inputs.md)

### Decision for 0PN

Prioritize Step 2. The finite-box drag calibration and outgoing-wave standards are available for a supported geometry. A physical equilibrium still requires an independently prepared load-bearing wall, its mechanical tension, the trapped-wave occupation and the reactions where the wall joins the rest of the material. The next discriminator is the excess normal force on a weakly curved patch of that same wall, measured relative to a flat patch at the same material state and wave preparation. Opposite curvatures and a second curvature scale would check whether one surface tension predicts the response. The known conversion-layer energy can supply the tension only if that layer carries the wall load. With calibrated loads, let several depths move together and reserve other depth bands to test the resulting shape. The radius-4 shortfall becomes a physical universality test when one such equilibrium family predicts mouths of different radii. [What is known about the wall location](mouth_surface_registration.md) · [Five-stage roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

**Independent review.** The shape reviewer reported “Scoped clear for the fixed-intake five-shape test, conditional restoring slopes and eight-band exclusion at the sampled shapes.” A second reviewer reported “Scoped clear for the eight-band equilibrium bounds.”
