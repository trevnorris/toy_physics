# One rounded mouth can turn the current into the bulk

At 0PN a held-open, isolated mouth admits a stationary flow with one phase from the feeding slab into the receiving bulk. The current crosses inward at every sampled radius from the center to just inside the lip. At a tangent attachment it follows the lip toward the bulk, so the crossing at the edge itself is zero. The flat-feed continuation in the [previous test](mouth_responding_rim.md) had given outward flow near this edge.

## What the solve selects

Assume the longitudinal material remains coherent as transverse shear order is lost. A single neutral carrier obeys

\[
-\tfrac12\Delta_4\psi+\tfrac14(|\psi|^8-1)\psi=0,
\qquad \mathbf j=\operatorname{Im}(\psi^*\nabla\psi).
\]

The radius-one opening connects a laterally fed slab to a bulk reservoir around a rounded reflecting lip. Both distant baths have unit density. Adjusting their phase difference to 0.41696727 radians produced intake 0.25. The local density, phase slope and direction of discharge were free throughout the opening. The rounded wall was held fixed and its attachment angle chosen for tangency at this intake with the trial coefficient \(K=0.4\). The conversion surface \(w=h(r)\) was then reconstructed from the current using \(-K\Delta_3h=q\), without feeding a surface force back into the carrier.

Here \(q=j_w-h'j_r\) measures inward crossing per unit *projected* area. On the finest mesh it was 0.186798 at the center and 0.000918654 at \(0.99R\). All 2,000 sampled interior radii had positive crossing. A separate integration of radial current on the same field gave inward transfer in all 51 shells. At the lip \((j_r,j_w)=(-0.0450431,+0.00224735)\). Tangency sets \(q(R)=0\) in the continuum. The numerical residual there fell from \(8.74\times10^{-5}\) to \(7.07\times10^{-6}\) with refinement and larger baths. The integrated crossing was 0.25000158, and the two bath readings agreed with 0.25 to within \(9\times10^{-12}\). [Flow figure](figures/mouth_curved_phase.png) · [Fields, commands and checks](../../../research/pde_ledger_v3/_measurements/mouth_curved_phase.json)

An [independent analytic example](../tools/mouth_confocal_phase.py) gives strictly inward crossing throughout the interior and exact grazing in a different rounded geometry. Halving or doubling the receiver’s density stiffness relative to the feed retained the inward sign on a fixed tangent cap. Three finite meshes had a single phase throughout every cell and a positive fixed-bath energy Hessian. The fixed-wall state is energetically stable under the specified baths. Its relaxation time and the physical origin of coherence at conversion remain open. [Independent stability check](mouth_curved_phase_stability.md)

## What the conversion layer must supply

The homogeneous smooth carrier gives equal pressure on the two sides of the same surface. If the earlier trial tension law \(p_f-p_b=-\sigma\Delta_3h\) also applies, it requires an integrated pressure difference \((\sigma/K)Q=5/24\) at \(\sigma=1/3\). A smooth pressure gradient through the opening supplies zero jump. The next isolated-mouth test is whether loss of order supplies the missing normal stress while preserving the through-flow. This is an input for a future two-mouth force balance. The present calculation measures no acceleration. [Interface balance](../tools/mouth_curved_phase_balance.py)

Working result on base `49dca23e`. Independent local checks completed. Formal project review remains pending.
