# Can the mouth drain while an intact face stays quiet?

**Uniform mobility gave nearly the same conversion per area at a supported flat face and a curved drain. A curvature threshold can seal a flat face, but it does not reliably select a mouth.** We used the same coherent carrier, equal-energy order wells and passive conversion law in both preparations. The curved control is a four-dimensional annulus with a free conversion surface and a specified shear mirror. The flat control has a stationary restraint whose force is measured. [Prior energy and release tests](mouth_coupled_conversion.md)

### Conversion rate and energy

Positive current is toward the bulk, per unit three-dimensional front area. All entries below are stationary pressure-driven controls. The wave-occupied curved reference has trial action 0.001.

| Matched at the reservoirs | Flat current | Curved current | Flat / curved |
|---|---:|---:|---:|
| Static pressures, equal path span | 0.0134672 | 0.0138028 | 0.975686 |
| Specific energy drop, equal span | 0.0138058 | 0.0138028 | 1.000219 |
| Static pressures and saved spatial phase difference | 0.0135538 | 0.0138196 | 0.980765 |
| Both endpoint densities, phases and specific energies | 0.0138243 | 0.0138196 | 1.000338 |

The last two rows use +0.450142783 radians. Matching the phases required solving for the flat path length. In the last row the flat length was 32.4876 interface widths compared with 32 for the curved path. The flat ends had measurable quantum stresses. In that fully matched pair, the reservoirs supplied **0.00462897 energy per constituent** to each preparation, all accounted for by conversion loss. The homogeneous order wells supplied zero net energy. The earlier **0.05123** was a cost on the saved radius-one throat in its own units. A metastable brane could in principle release energy as it converts. The written wells give no such energy difference, and a uniform release would also be available on intact faces. [Matched pressures and head](../tools/mouth_matched_flat_rate.py) · [Simultaneous phase and energy control](../tools/mouth_matched_reservoir_phase.py)

The flat restraint carried a force of −0.00449917 per area in the first comparison. Without support, a free planar interface with equal wells and vanishing end gradients requires zero steady material crossing in its own frame. At unequal bath pressures that leaves an unbalanced normal load. Physical steady leakage from a freely moving slab and any change to the exterior inverse-square law remain unmeasured. A gravity test would need the change in face flux caused by the drain after subtracting the undisturbed background.

### Curvature threshold

We tried a kinetic law with $b=\delta_B\sqrt{k_1^2+k_2^2+k_3^2}$ on each regular order surface. Mobility rises smoothly from zero at $b_c$ to its uniform value at $2b_c$. The thresholds are trial material parameters. The lock seals an intact plane by construction. The stationary locked face still requires normal support under unequal pressures.

At $b_c=0.006$, the fully responding spherical control drained at **0.00294686** with the wave and **0.00293468** without it. At the same threshold, **all sampled points** of the saved mouth’s half-order surface were locked on three meshes using two curvature estimates. At thresholds 0.002 and 0.004, the finest mesh retained about 73% and 0.76% of that surface’s uniform conductance under the same local driving force. These are fixed-geometry conductances. A smoothly curved, unpunctured spherical patch also activates the law. The threshold has no way to identify a puncture of the same curvature. [Coupled curvature test](../tools/mouth_curvature_lock.py) · [Saved mouth and intact controls](../tools/mouth_order_bend_gate.py)

### One connected carrier

A separate weak-drive test put the opening region and a flat annulus between the same two uniform phase boundaries. At rest with positive carrier density $n_0(w)$, number conservation gives $\nabla\cdot(n_0\nabla\theta_1)=0$. With reflecting outer sides and continuous phase, the unique solution has $j_w=\Delta\theta/\int dw/n_0$ throughout the cross-section. The saved phase difference across 32 length units gave **0.0140670** current per projected area even on the flat annulus. This is current through a fixed plane in a neutral reference cell. A locked material surface there must move with the carrier. Remaining stationary would require zero normal current. This first-order incompatibility does not assume a shape instability or model the nonlinear throat. [Joined-cell calculation](../tools/mouth_locked_face_phase.py)

**Next discriminator.** Compare the normal mass current versus a weak applied phase difference in an intact bent patch and a punctured patch with matched local curvature. Measure the local phase gradient, surface velocity and normal support alongside conversion. This tests whether the puncture changes coherent transport itself or the driving phase that reaches the face. A change in conversion mobility alone leaves the neutral carrier equation unchanged. [Guarded measurements and independent checks](../../../research/pde_ledger_v3/_measurements/mouth_face_conversion_lock.json)

**Independent review.** “Scoped clear for the supported-flat versus free-spherical controls under matched reservoir loads and phases, and for stationary annular drainage under the trial curvature lock.” A separate reviewer found “Scoped clear for the joined neutral phase-cell incompatibility.”
