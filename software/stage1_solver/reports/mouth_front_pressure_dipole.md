# Does the companion’s pressure push a bowed mouth sideways at 0PN?

**The held drain sees higher feeding-side pressure toward its companion.** On the finest saved grid, doubling the separation reduced its pressure moment to 0.25015 of its earlier value. How much of that load reaches the rim depends on the pressure from the receiving bulk. A measurement of how the conversion surface tilts around the rim could recover the difference without resolving the throat interior. [Recorded runs and scoped reviews](../../../research/pde_ledger_v3/_measurements/mouth_front_pressure_dipole.json)

### What the feeding field supplies

Both mouths have radius \(R=1\) and prescribed inward number rate 0.25. With \(\xi\) positive toward the companion, take the pressure change from the same target alone, \(\delta p_f=\delta(n^5/5)\), on its flat feeding face. Integrate \(M_f=\int_{B_R}\xi\,\delta p_f\,d^3x\) over the round opening. An independent integration using exact intersections of grid cells with the ball agreed within 0.56% on the finer mesh. [Field calculation](../tools/mouth_front_pressure_dipole.py) · [Independent integration](../tools/mouth_front_pressure_cells.py)

| Reservoir extent, cells per radius | \(M_f\) at \(8R\) | \(M_f\) at \(16R\) | Ratio |
|---|---:|---:|---:|
| 48, 4 | 4.88898 × 10⁻⁶ | 1.19729 × 10⁻⁶ | 0.24490 |
| 96, 4 | 4.90310 × 10⁻⁶ | 1.22601 × 10⁻⁶ | 0.25005 |
| 96, 6 | 5.01709 × 10⁻⁶ | 1.25504 × 10⁻⁶ | 0.25015 |

These two distances test finite-range falloff. The prescribed normal speed \(j_w/n\) is reflection-even by construction. Density and inward current can redistribute despite that restriction. Including carriage and density-gradient stress gives a **different** observable. At \(8R\) on the fine grid the first moment of the full normal carrier momentum flux was \(1.63279\times10^{-5}\) from the field equation and \(1.57698\times10^{-5}\) from the saved depth layers. Substituting either for an isotropic pressure jump would require an interface law. [Carrier stress convention](../../em_charge_attribute/reports/rim_normal_return.md)

### What the rim could measure

Assume a quasistatic, thin conversion surface with constant tension \(\sigma\), small slopes, a radial initial depth \(h_0(r)\), and a fixed spherical edge as its sole support. Let \(P=p_f-p_b=-\sigma\nabla^2 h\) be the actual normal pressure jump there. The sideways pressure force is \(-\int P\partial_\xi h\,d^3x\). Integrating by parts, retaining the change in shape, gives the linear force on the rim

\[
\delta F_\xi=\frac{s_0}{R}(M_f-M_b),\qquad s_0=-h_0'(R),\qquad M_f-M_b=\sigma A s_1,\qquad A=4\pi R^3/3 .
\]

Here \(s_1\) is the dipole of the change in inward rim slope, \((3/4\pi)\int \cos\theta\,\delta(-\partial_r h)\,d\Omega\). The moment identity follows from Green’s identity with the harmonic coordinate \(\xi\). Looking only at the two points nearest and farthest from the companion also picks up higher odd angular patterns. [Nonuniform load and shape check](../tools/mouth_front_pressure_check.py) · [Angular readout check](../tools/mouth_front_rim_pressure_readout.py)

**Discriminating test.** If the measured flat-face pressure change reaches the bowed surface and the bulk contributes zero pressure moment, \(\sigma s_1\) should be \(1.198\times10^{-6}\) at \(8R\) and \(2.996\times10^{-7}\) at \(16R\) on the fine grid. Equal feeding and receiving moments give zero. In the first case an independently measured slope \(s_0=0.1\), for example, would give a sideways pressure force about 0.65% of the flat opening’s carried-plus-stress sideways momentum rate. The latter is a comparison scale. The current archives contain neither the responding shape nor a matched receiving pressure. Measure the full angular slope and tension, or both local pressures on the same front. A free-fall claim additionally needs the total force transmitted to the mouth and its independently determined inertia. Working result on `49dca23e`, formal review pending. [Receiving inlet comparison](mouth_inlet_trace.md)
