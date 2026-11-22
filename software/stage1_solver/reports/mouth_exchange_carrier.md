# Carrier work on the radius-one neck and flare

**The carrier supplied positive stiffness and a first-order force toward a narrower entrance on the tested permeable flare.** At fixed spatial phase difference 0.45100843, the density-relaxed calculation gave energy slope **+0.76696** and stiffness **+0.19344** for the archived volume- and length-preserving change. The generalized force exerted by the carrier is minus the slope. A local pressure decrease at a constriction alone does not give the sign of the integrated stiffness. [Measurements and source hashes](../../../research/pde_ledger_v3/_measurements/mouth_exchange_carrier.json)

### Matching the same material on both sides

In the stationary reversible Hamiltonian, number flux is conserved across a thin transition, while the phase gradient exceeds the physical velocity by the entrained order momentum. Integrating the excess gradient across a planar kink gives a conservative phase jump. With unit constituent mass and ambient density, the trial extension to a curved front is

\[
 [\theta]=\lambda_0 j_n/n_\Sigma^2,\quad
 E_\Sigma=\int n_\Sigma^2[\theta]^2/(2\lambda_0)\,dA,\quad
 \lambda_0=I_\chi/(6\delta).
\]

Here the exit-side phase minus the core phase defines the jump; positive current leaves the ordered side. We assumed that the density varies slowly through the transition and shares a common amplitude trace. The kink shape is held registered during the deformation. Setting the dimensionless ratio \(I_\chi/(n_0\delta^2)=1\) as on the planar patch and choosing \(\delta/R=0.04\) gives \(\lambda_0/R=0.006667\). The physical width and order inertia in mouth units remain uncalibrated. In pure material the energy density is \(|\nabla\sqrt n|^2/2+(n^5-5n+4)/20+n|\nabla\theta|^2/2\), with no dependence on which side of the registered surface the material occupies. Solving for both fields includes pressure, flow, and density-gradient work on both sides. The interface stores energy reversibly; its steady dissipated power is zero. The neutral surface tension is counted in the shape calculation. [Hamiltonian material patch](mouth_reversible_front.md) · [Finite-drive solver](../tools/mouth_exchange_carrier_density.py)

The ordered unit-thickness slab opens into the same radius-1 to 1.5 turning flare as the wave. The intact slab face is sealed. Coherent receiving material fills the region immediately outside the entire side wall and beyond its terminal disk. The phase is fixed on the remote feeding rim, receiving rim and roof, with unit density at those ports. The receiving rim starts above the sealed plane. Radial distance to the side port is 28 and roof height is 22. The sharp interface passes current through its side and disk. At zero deformation the new through-current is 1.89674, with 1.81363 crossing the side and 0.0831044 crossing the disk. These are outputs on a prescribed meridian at fixed phase. The older prescribed intake and the intake of a free material mouth are different preparations.

### Work and controls

Two separate readouts used the signed stationary-energy difference and integrated variational stress. In the weak-flow limit the bulk derivative contains \((|\nabla\theta|^2\mathbf1/2-\nabla\theta\otimes\nabla\theta):\nabla\mathbf V\), with the corresponding area derivative on the permeable sheet. The density-relaxed version includes the potential and density-gradient stress and the measured local density in the phase terms. The prescribed shape speed reaches both the aperture and the downstream junction. Its continuation through the two regions moves integration coordinates at fixed bath position. This keeps the end work in the force readout and avoids integrating a one-sided gradient at a corner. On-shell variation removes the implicit field derivatives. [Two-sided assembly](../tools/mouth_exchange_carrier.py)

| Fixed-port comparison | Energy slope | Energy stiffness |
|---|---:|---:|
| Full drive, mesh 64, signed 0.005 steps | +0.766956 | +0.193441 |
| Full drive, mesh 64, signed 0.010 steps | +0.766943 | +0.193453 |
| Full drive, mesh 40, signed 0.005 steps | +0.767749 | +0.192992 |
| Half drive, mesh 64, signed 0.005 steps | +0.191871 | +0.047767 |
| Leading quadratic response evaluated at full drive, mesh 80 | +0.767479 | +0.190395 |

Integrated full-drive stress at the reference gives +0.766960. Differencing the stresses gives stiffness +0.193445. Under the documented coordinate variation, the feeding slab contributes +0.12805, the ordered conduit −0.00742, receiving material beside it +0.07870, material above the exit −0.00310 and reversible sheet storage −0.00279 to that stiffness. The partition depends on the continuation of the virtual displacement through the fluids. Only their sum is the work conjugate to the tested change. At leading order, changing to a fully transparent internal mark gives total stiffness +0.18881 and taking \(\lambda_0=0.1/6\) gives +0.19150. Moving the ports inward to 14 and 11 gives +0.23180. These sensitivities put the finite-reservoir assumption in the numerical result.

A second, independent volume-preserving disturbance holds both end radii fixed and expands one portion of the wall while narrowing another. Its full-drive energy slope is +0.00064127. At mesh 80 the quadratic carrier gives +0.00061768 there; the corresponding transparent-mark null test leaves a mesh residual of +0.00003206. The exact transparent continuum with fixed aperture has no material mark in its energy. The quadratic primary slope persists at +0.77013 in that same transparent limit, where only changing the open aperture affects the continuum carrier energy. [Held-end readout](../tools/mouth_exchange_carrier_interior.py)

The numbers use the carrier energy and mouth-radius units above. They can be included with \(J\,\omega'\), \(J\,\omega''\) and the measured \(-2.70831\,\sigma\) only in a common normalization and material preparation. Actual interface registration, mouth-unit tension, wave action and transfer of wave stress to the two-sided medium remain physical inputs. [Reference shape](mouth_conversion_routes.md)

The first mesh-24 quadratic run emitted its signed fields, then failed the unextrapolated full-step force comparison. It is excluded from the table. The finer independent runs checked both signed step sizes and retained the original force tolerance on the step-extrapolated readout.
