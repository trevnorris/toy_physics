# How much sideways momentum reaches a draining mouth at 0PN?

On the spherical surface in the relaxed viscous exterior, a held mouth with strong normal contact takes in sideways momentum at rate $\dot M U/3$ to leading order. Here $U$ is the incident brane-material speed and $\dot M=\rho Q_t$ is the positive material mass intake. Material can slide along the surface even when its sideways-induced *normal* velocity there vanishes. SymPy integrated the momentum and pressure/viscous stress on arbitrary spheres. Mathematica reconstructed the velocity from a stream function and obtained the stress change from the momentum equation in the volume between spheres. [Both guarded outputs](../../../docs/_measurements/mouth_conversion_momentum.json)

**Surface calculation.** Keep the radius $a$ and radial drain fixed. Write $V$ for the mouth's lab speed, $W=U-V$, and $\mu$ for the cosine to the incident axis. At the surface, the relative outward radial velocity is $-Q_t/(4\pi a^2)+b\mu$, where $b$ is the response to the weak sideways drive. The calculation uses the previous incompressible exterior with constant viscosity $\eta$ and zero incremental tangential traction. Take the linear response at finite intake, with $|b|\ll Q_t/(4\pi a^2)$, and retain terms through first order in $\epsilon=\dot M/(4\pi\eta a)$. [Exterior solution and contact law](rim_material_entrainment.md)

Integrating inward mass flux times lab axial velocity over the whole sphere gives

$$
J_{\rm in}=\dot M\left(V+\frac{W}{3}+b\right).
$$

Changing the material velocity contributes $\dot M(V+W/3+2b/3)$. Redistributing the intake over the existing radial flow contributes $\dot M b/3$, even though total intake stays fixed. As a conservation check, the sum of mechanical stress and inward carriage in the mouth frame is $4\pi\eta a(W-b)+\dot M(W+b)/2$ at every exterior radius through this order. This is the supply to everything enclosed. [Angular integrals and two-surface check](../../../docs/_measurements/mouth_conversion_momentum.json)

**What reaches the opening itself.** At each conversion event the model sends material into the bulk with its local in-plane brane velocity. If conversion occurs at the sampled sphere with the same local mass flux, all of $J_{\rm in}$ is carried onward. The carriage then gives zero additional impulse on the opening at conversion. The mechanical rim traction still acts. The physical conversion surface and the spherical integration surface have no established pointwise identification. Matching total intake alone cannot establish their momentum match. [Local conversion rule](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#L220)

**Open measurement.** Measure the isolated mouth's steady drive force per speed $C_{\rm drift}$ in resting brane. With a held companion, compare the pull $-\Delta H$ to $C_{\rm drift}U_s$, using the companion-only incident **material** speed $U_s$ at the target. A resolved residual calls for checking the in-plane traction and the flux-weighted material velocity at the actual mouth. The archived separation-8 and separation-16 runs save carrier fields and currents. They supply neither an independent material velocity nor an actuator reaction. Free acceleration remains open. Working on base fd667ea9, uncommitted pending review. [Archive inventory](../../../docs/_measurements/mouth_conversion_momentum.json) · [Matched carrier comparisons](companion_uniform_stress.md)

A subsequent comparison found that an independently flowing receiving bulk changes this prediction. Towing tests the response to both incident streams together. A held run with brane flow alone separates their force coefficients. [Towing and companion comparison](mouth_tow_and_companion.md)
