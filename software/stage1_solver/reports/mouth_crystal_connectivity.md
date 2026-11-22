# Can crystalline connectivity select the open exit?

**Decision.** The coherent crystalline-field test admits a conducting defect but does not select it at the mouth. Park microscopic sealing. For conditional Steps 3 and 4, assume **sealed intact slab faces and an open downstream melting exit**. The origin of this boundary condition, physical support of the flare and the energy source for sustained conversion remain open. [Previous carrier barrier and smooth-exit test](mouth_bent_puncture_transport.md)

We compared an intact bent face with a smooth converting patch at the same shear order and curvature. In addition to the coherent carrier phase, we introduced a complex crystalline density-wave amplitude. Its phase specifies the positions of crystalline layers. We let its amplitude relax, so a dislocation could have a smooth core. We used a local tube of radius eight interface widths on a three-sphere of radius 100 widths, translationally invariant along the defect. The trial normal phase stiffness was
\[
K_\perp/K_0=1-4\chi(1-\chi)f^2 .
\]
Here $f$ measures crystalline registry in addition to the common ordered fraction $\chi$. The positive crystal energy penalizes gradients and favors $f=1$. Neither a mouth-specific energy nor a defect source was inserted. The material coefficients and the local relation to phase stiffness are hypotheses. [Field calculation](../tools/mouth_crystal_connectivity.py)

| Crystalline preparation | Result in the trial law |
|---|---|
| Intact face, registration free to choose | Uniform crystal minimizes the energy. Normal conductance is zero. |
| Smooth melting exit, registration free to choose | The same uniform crystal minimizes the energy. Normal conductance is zero. |
| One imposed positive or negative winding at the exit | A resolved dislocation opens a finite path, costing 7.741 model energy units per unit defect length. |
| The same winding on intact material | The same conductance and energy cost are obtained. |

For the imposed defect, independent current and phase variational bounds give **55.47% to 67.55%** of the transparent patch conductance. An independent calculation resolving lateral and normal currents gives **63.18%**. [Independent field calculation](../tools/mouth_crystal_connectivity_check.py) The defect phase concerns the crystal, while the transported current belongs to the coherent carrier. A spherical opening does not require a crystalline dislocation around it. The percentage characterizes this prepared patch and gives no independently selected mouth intake.

The measurement establishes what such a dislocation *could* do under the proposed coupling. A theory that predicts where dislocations form, their persistence during melting and their relation to carrier stiffness would be needed to adopt it as the sealing law. The gravity calculations now use the declared boundary assumption. [Measurements and scoped review](../../../research/pde_ledger_v3/_measurements/mouth_crystal_connectivity.json) · [Free-fall roadmap](../../../docs/zero_pn_free_fall_roadmap.md)
