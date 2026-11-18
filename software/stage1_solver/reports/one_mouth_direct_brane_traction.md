# Direct tangential loading of the brane at 0PN

A traction-controlled calculation can bypass the mouth’s contact coefficient. Prescribe the change in **net tangential force per unit area on the brane** and read the change in material velocity at the same point. I ran this test on the homogeneous exterior of the radius-one opening. The physical response of the inward-draining brane remains open. [Rotational traction at a boundary](pathA_23_stage3_noleak_closure.md#L64) · [Draining-state material inputs](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L123)

**Test and readings.** The reference brane starts at rest and has unit density and transverse wave speed. The load is the projection of a sideways push onto the tangent plane around the spherical rim. Its signed equatorial value is \(q(t)=q_* (t/T)^3\) for \(0\le t\le T=0.001\). The geometric opening stays fixed. The solver imposes traction on the shear field and obtains the equatorial velocity from the material momentum equation using the radial gradient on the brane side. The loads below are in solver traction units and the readings in solver velocity units. [Direct-drive solver](../tools/brane_direct_traction.py) · [Run manifest](../../../docs/_measurements/brane_direct_traction.txt)

| Net equatorial traction at \(T\) | Equatorial material velocity at \(T\), 1536 cells |
|---:|---:|
| \(+2\times10^{-5}\) | \(+1.99999580748\times10^{-5}\) |
| \(-2\times10^{-5}\) | \(-1.99999580748\times10^{-5}\) |
| \(+1\times10^{-5}\) | \(+9.9999790374\times10^{-6}\) |
| \(-1\times10^{-5}\) | \(-9.9999790374\times10^{-6}\) |
| \(0\) | \(0\) |

For this outgoing shear mode, SymPy and Mathematica give the transfer from equatorial traction to equatorial velocity. With rim radius \(a\), brane density \(\rho\), wave speed \(c_T\), \(\tau=a/c_T\), and Laplace variable \(s\), it is

\[
\frac{\widehat v}{\widehat q}=\frac{1}{\rho c_T}
\left(1+\frac{1}{s\tau(1+s\tau)}\right).
\]

The maximum relative velocity errors against Mathematica at the three sampled times fell from \(0.0005441065230921567\) to \(0.0001367379358602001\) to \(0.00003427336119408571\) on successive grid doublings. The zero-load readings were exactly zero. The nearly equal numerical traction and velocity at the end of this short ramp reflect the chosen material units and \(T/\tau=0.001\). [Independent continuum calculation](../tools/brane_direct_traction_check.wl) · [Guarded outputs](../../../docs/_measurements/brane_direct_traction.txt)

**Requirement for the draining test.** Keep the radius, inward intake, feed state and mouth support fixed across the signed loads and the zero-load run. Record the perturbation of the brane-material velocity at the equator relative to that baseline. Calibrate the net traction at the same point, including any change in contact reaction. For example, on the trial linear-slip branch with the mouth held fixed, \(\delta q_{\rm contact}=-\gamma\,\delta v\). Using the actuator load alone would include this contact in the inferred impedance. The archived carrier evolves density and phase. A draining-state material evolution that supplies stress and inertia about that background, or an independent measurement on that preparation, is still needed. Direct loading would then determine the brane response over the probed time or frequency range. The earlier imposed-mouth accelerations would still need a material observation to fit contact strength. [Carrier archive](../../../docs/_measurements/mouth_material_readout.txt) · [Rim contact law and calibration](one_mouth_exterior_shear.md#what-the-interface-equations-determine-at-one-rim-point)

### Can the same-intake subtraction give the draining response?

The saved radius-one carrier preparation has inward number rate \(0.25\) and zero sideways reservoir speed. Its field is the carrier amplitude and phase. A brane mass-density profile, material velocity and stress history are absent. The map from the carrier’s normal number current to the brane material loss \(j_n\) is also open. The uniform shear action supplies the reference test above. The existing nonuniform operators take supplied background profiles and do not supply a momentum or stress law for this draining preparation. [Packet and calculation](../../../docs/_measurements/draining_rim_material_balance.txt) · [Current and mass conventions](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L205) · [Scope of the material equations](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L377)

Even granting the same **pointwise** \(j_n\) and brane density in both preparations, \(\nabla\cdot(\rho_{\rm br}\mathbf V)=-j_n\) only fixes the divergence of the mass current. We constructed a continuously differentiable mass-current perturbation supported in a band around the rim. Its divergence and radial current at the rim are both zero. Its equatorial sideways current has arbitrary amplitude \(b\), giving a velocity change \(b/\rho_{\rm br}(a)\) there. SymPy and Mathematica verified this in spherical and Cartesian coordinates. This checks the freedom in continuity alone. The prescribed traction would select the motion through a material momentum equation. Matching intake also leaves a change \(j_n\,\delta V_x\) in the sideways momentum carried by exchanged material wherever loss occurs, even with \(\delta j_n=0\). That is a transport term to account for in the balance. [Exact checks](../../../docs/_measurements/draining_rim_material_balance.txt) · [Carried momentum convention](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L220)

For the forced-minus-unforced velocity at fixed intake, the missing input is the **linear transverse momentum and stress response about this draining state**. A field implementation needs the material background and the changes in momentum storage, transport, shear stress and any exchange partners under weak loading. For the rim observable alone, an independently calibrated relation between net traction history and velocity change on this preparation would suffice. Prescribing net traction already bypasses the contact coefficient. The physical draining solver and its paired rim-velocity measurement remain open. [Live momentum and stress requirements](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L165)

Working result pending review. Repository base `fd667ea9`. No new commit.
