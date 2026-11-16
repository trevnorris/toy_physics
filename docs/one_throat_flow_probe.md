# Holding force and acceleration of one draining intake

**Result, November 12, 2026.** For a spherical intake that turns the incoming fluid fully into the bulk, reversing a weak sideways flow reversed the required holding force. Letting the outlet carry the incident sideways momentum through gave zero holding force at both speeds. The exterior flow, intake rate and body mass were identical. A separate small force probe gave the same incremental acceleration in both cases. These are results for the specified three-dimensional benchmark. The native throat’s outlet stress and translational response remain open. [Run output and source excerpts](_measurements/one_throat_flow_probe.txt)

## Measurements

The radius is $a$, the volume intake rate is $Q>0$, and the exterior density is $\rho$. Units are $a$ for length, suction speed $s=Q/(4\pi a^2)$ for velocity, and displaced mass $m_d=4\pi\rho a^3/3$ for mass. Force is in $m_d s^2/a$ and acceleration in $s^2/a$. I chose body mass $M=m_d$ for this test. Each row starts with the mouth at rest and a settled flow. The columns give the support force, then the initial free acceleration with and without an applied probe force of $+0.015$. The drain stays on throughout. [Exact results](_measurements/one_throat_flow_probe.txt)

| Outlet’s sideways velocity at rest | Incident flow | Holding force | Acceleration without probe | Acceleration with probe |
| --- | ---: | ---: | ---: | ---: |
| Zero | +0.01 | −0.03 | +0.02 | +0.03 |
| Zero | 0 | 0 | 0 | +0.01 |
| Zero | −0.01 | +0.03 | −0.02 | −0.01 |
| Same as incident flow | +0.01 | 0 | 0 | +0.01 |
| Same as incident flow | 0 | 0 | 0 | +0.01 |
| Same as incident flow | −0.01 | 0 | 0 | +0.01 |

## What sets the two responses

The exterior is incompressible, inviscid potential flow. The spherical surface admits fluid at uniform inward normal speed $s$ relative to the moving mouth, with tangential slip. Both outlet cases discharge mass at $\rho Q$ with zero additional lateral outlet traction. If the mouth moves at $V$, the prescribed mean outlet velocity along the flow axis is $V+k(U-V)$, with $k=0$ or $1$. The model gives the intake a constant mechanical mass and neglects internal tilt, internal momentum storage relative to the body, and radiation. These are benchmark choices.

Solving for the exterior velocity and integrating the convective and pressure currents gives inward force $\rho Q U-(m_d/2)A$ at the instant of release. The outlet carries away $k\rho Q U$ then. The support force $H$ for steady holding and acceleration $A$ with the support removed therefore satisfy

$$
H=-(1-k)\rho Q U,\qquad
(M+m_d/2)A=F_{\rm probe}+(1-k)\rho Q U .
$$

In these units the effective mass inferred from the probe is $0.015/0.01=1.5$. One unit was supplied as the body’s mass. The other half came from the pressure needed to accelerate the exterior fluid. The result was computed with two symbolic engines. All 24 entries agreed, including the negative-probe control. The no-drain limit of the pressure term also gives the classical spherical added mass. [SymPy calculation](../software/stage1_solver/tools/one_throat_flow_probe.py) · [Wolfram calculation](../software/stage1_solver/tools/one_throat_flow_probe.wl) · [MIT sphere derivation](https://ocw.mit.edu/courses/2-016-hydrodynamics-13-012-fall-2005/7cfc1ce5fb298d4f54beaf36b9f00109_add_mass_sphere.pdf)

## What the drain model still owes

The first quantity to derive for the native mouth is the change in **total sideways momentum per unit time leaving through the sleeve and slab faces** when the imposed flow reverses. That includes carried momentum and traction. The native order-conversion rate alone supplies neither. For the probe, the existing one-body calculation supplies part of the translational inertia. Core and surface contributions, coupling to tilt, and the acceleration-dependent exchange remain unresolved. The assumed $M$ and computed spherical added mass above cannot supply those values. [Release balance](drain_driven_initial_acceleration.md) · [One-body results, phases B1 and B2](../software/em_charge_attribute/reports/u1_body_dynamics.md) · [Face and inertia requirements](../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md)

Repository base `d544ef73`. Local working result, pending independent review.
