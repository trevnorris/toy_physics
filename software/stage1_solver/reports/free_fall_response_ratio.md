# Does the unknown receiving-side force scale cancel in 0PN free fall?

**Momentum conservation allows the same steady pull and two different slow accelerations.** This is a test of what follows from the current equations. The two cases are alternative laws for the unresolved sideways response of one specified drain. The radius, intake, exterior and centered radial preparation are common. The calculation determines no physical coefficient for the archived mouth. [Open mouth inputs](../../em_charge_attribute/u1_body_dynamics_inputs.yaml#L48)

**Held mouth.** Write \(U\) for the weak incident speed toward the companion and \(CU\) for the steady pull on a fixed mouth. There is an explicit positive example. In the ideal constant-density exterior, superpose a uniform stream on a spherical sink of radius \(R\) and inward speed \(w\). Integrating the carried momentum and Bernoulli pressure over the sphere gives \(C=4\pi R^2\rho w\), the material mass intake. Take the single-holder benchmark with the same other loads in both cases. The physical drain can have a different \(C\). [Exterior assumptions](held_force_angular_falloff.md) · [Moving-mouth and remote-flow conditions](mouth_tow_and_companion.md)

**Moving mouth.** Let \(q(z,t)\) be the incremental sideways displacement of a bound material or wave response over depth \(L\). Couple \(q(0,t)\) to the geometric mouth displacement \(X(t)\). Choose a free far end for this lateral stress, and zero direct force from the incident stream on this sector. Material still drains. Its quadratic energy is \(\tfrac12\int_0^L\lambda\mu(\dot q^2+c^2q_z^2)\,dz\). Each positive \(\lambda\) gives the same signal speed \(c\), motion for prescribed \(X(t)\), and zero incremental energy and stress in the centered, held state. Only the normalization of this *translation perturbation* varies. Whether the background and the perturbation share a normalization in a microscopic throat is open. [Projected field momentum and traction](wall_wave_cut.md)

The force the rim exerts on this field is \(-\lambda\mu c^2q_z(0)=\frac{d}{dt}\int_0^L\lambda\mu\dot q\,dz\). The field exerts the opposite force on the rim. For a slowly varying drive proportional to \(e^{st}\), \(sL/c\ll1\), the field follows the mouth and its force per acceleration approaches \(\lambda m\), with \(m=\mu L\). Write \(M\) for the common slow inertia of the other moving fields. With the steady pull as the only incremental drive and the common drag absent, the two free responses are

\[
\left.\frac{a}{U}\right|_{\lambda=1}=\frac{C}{M+m},
\qquad
\left.\frac{a}{U}\right|_{\lambda=2}=\frac{C}{M+2m}.
\]

Both conserve momentum. For positive \(C\) and \(m\) the second accelerates less. If the sector also adds steady pull \(\lambda B U\), the ratio \((C+\lambda B)/(M+\lambda m)\) is scale independent exactly when \(BM=Cm\) in this family. The passive example has \(B=0\). These are slow-response coefficients. An abrupt release can excite internal transients.

**0PN requirement.** The physical mouth must select the sideways pull and co-moving inertia together, or directly supply their ratio. Target-independent free fall in the same source flow requires a common positive ratio across physical mouth radii. An effective law measured at the mouth could supply this information. The need for a full interior solution remains open. [Radius and source requirements](../../../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#how-drain-rate-and-retained-momentum-must-scale-for-free-fall)

SymPy used the continuous field equation. Mathematica independently integrated the exterior flux and varied the displacement of a two-element field, including its free endpoint. Both obtained the same slow gains and zero momentum-balance residual. Working result on base cda9a8eb, pending independent review. [Scripts, commands and literal results](../../../docs/_measurements/free_fall_response_ratio.json)
