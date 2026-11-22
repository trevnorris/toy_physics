# Can a pressure gradient support the existing flare?

**An arbitrary constant pressure gradient with depth leaves a measurable imbalance on the specified cubic wall.** At radius 1 and intake 0.25 we compared eight radial displacements in the archived, steadily draining preparation. Allowing the uniform wall tension, ground-wave action, backing pressure and its gradient all to vary leaves a minimum worst discrepancy of \(2.07892\times10^{-5}\) in model pressure on mesh 48. The two wall ends can exert forces, but all eight displacements vanish in neighborhoods of both ends. Endpoint work therefore cannot remove this discrepancy on the prescribed profile.

### What the free wall would have to satisfy

Write \(R(w)\) for the radius of a spherical section in the three brane coordinates. In four spatial dimensions an isotropic membrane has area \(4\pi\int R^2\sqrt{1+R'^2}\,dw\). With fixed axial endpoints a radial variation \(\eta=\delta R\) does total effective virtual work

\[
\delta W_{\rm resist}-\delta W_{\rm load}
=4\pi\int R^2(\sigma K-q)\eta\,dw+
4\pi\left[\frac{\sigma R^2R'}{\sqrt{1+R'^2}}\,\eta\right]_{\rm ends},\qquad
K=\frac{2}{R\sqrt{1+R'^2}}-\frac{R''}{(1+R'^2)^{3/2}} .
\]

Here \(q=p_{\rm carrier}+Jp_1-p_{\rm back}\) is the net outward load in the tested elastic law. \(p_1\) is the load per unit trapped-wave action. The stationary carrier pressure already includes its flow contribution. The formula assumes constant tension, no bending elasticity and radial virtual displacements. Joining to other wall pieces or a load-bearing edge adds their endpoint work. Fixing an endpoint radius sets \(\eta=0\) there and leaves its holding reaction to be accounted for. On a radially mobile junction the net endpoint force must balance. A freely moving axial endpoint would require axial work as well.

### Eight radial measurements

For each bump \(b_i(w)\), use the measured carrier pressure \(p_i\), curvature \(K_i\), and unit-wave pressure \(p_{1i}\), all weighted by \(R^2 b_i\). A backing \(P_0+gw\) has exactly the mean \(P_0+g\bar w_i\) with that same weight. We minimized the largest absolute error in \(p_i=P_0+g\bar w_i+\sigma K_i-Jp_{1i}\) over four unrestricted real coefficients. A signed combination of bands eliminates all four columns, has unit sum of absolute weights, and gives the same error bound as the fitted residual. The optimum on mesh 48 even uses negative tension and negative action. Requiring physically nonnegative values cannot improve the fit. [Calculation](../tools/mouth_free_equilibrium_inputs.py) · [Guarded data, coefficients and signed certificates](../../../research/pde_ledger_v3/_measurements/mouth_free_equilibrium_inputs.json)

| Carrier mesh | Minimum worst error from wall stress | From volume work |
|---:|---:|---:|
| 26 | \(2.05750\times10^{-5}\) | \(2.05323\times10^{-5}\) |
| 34 | \(2.07922\times10^{-5}\) | \(2.07906\times10^{-5}\) |
| 48 | \(2.07892\times10^{-5}\) | \(2.07916\times10^{-5}\) |

These finite-discretization results leave non-affine backing, a different equilibrium profile and independently supported extra stresses open. A chosen pressure profile could support any specified smooth curve through \(p_{\rm back}=p_{\rm carrier}+Jp_1-\sigma K\). This inverse relation explains why fitting a pressure at each depth cannot determine a physical mouth.

### Independent inputs for a release

The neutral order equation supplies the planar interfacial energy per three-dimensional area, \(\sqrt{2a_B\kappa_B}/6\). Using it as the membrane’s mechanical tension requires locating a load-bearing wall at that interface or measuring the force connecting them. The absolute bound-mode action is also free in the archived wave equation. A prescribed formation pulse can set an occupation conditional on that pulse. The steady-drain tests supplied no independently selected occupation at the reference intake. [Order energy](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md) · [Material surface](mouth_surface_registration.md) · [Formation](mouth_formation_capture.md) · [Wave power from the drain](mouth_wave_drain_power.md)

Prepare the two end junctions independently of gravity and measure their allowed radial motion and reactions. On the same preparation read the bulk-side normal traction by depth and the mode energy and frequency, or establish a stable gain–loss balance selecting its action. With a calibrated wall law, re-solve the carrier and wave on each trial free shape at the chosen supply. Test the resulting balance on unused bands and its recovery from opposite weak displacements. A direct opening-force measurement on a reproducibly self-supported mouth is another route. At present the material preparation supplies neither that mouth nor the loads needed to predict it. This blocks a physical flexible-mouth release and leaves the rigid carrier drag test available. [Equilibrium work on the reference shape](mouth_occupied_radial_balance.md) · [Development plan](../../../docs/zero_pn_free_fall_roadmap.md)

The independent five-band calculation and endpoint-work check returned “Scoped clear for the affine-backing exclusion on the prescribed cubic shape, its four-dimensional membrane variation and the conditional physical-release gate.” [Independent readings](../../../research/pde_ledger_v3/_measurements/mouth_rigid_drag_checks.json)
