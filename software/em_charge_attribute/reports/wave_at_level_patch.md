# Does the supporting wave prevent the conversion layer from leveling off?

**Reflection and wave pressure allow a level patch.** In a lossless shear model, a standing wave survives where the face becomes parallel to the resting brane and its projected normal vanishes. The wave continues to load the inner rim. For the tested flare the wave energy at fixed action decreases as the signed count is lost. These results concern one scalar polarization on prescribed geometry. The pressure balance of a freely draining mouth is still open. [Geometric unwinding path](mouth_count_unwinding.md) · [Guarded calculations](../../../docs/_measurements/wave_at_level_patch.json)

**Exact local test.** Follow a constant-width material strip through the level point. Along the north meridian of the flare at tilt \(\beta=5/2\), the height gradient vanishes halfway in arclength between \(r=2R\) and \(3R\). Let \(s\) increase along that strip from a reflecting node at 0 to the shear-free end at \(L\). Give the transverse inertia and stiffness the same linear order factor \(\chi=1-s/L\), with constant wave speed \(c\). The intrinsic equation \((\chi U')'+(\omega/c)^2\chi U=0\) admits

\[
 U(s)=J_0[j(1-s/L)],\qquad \omega=cj/L,\qquad j=2.404825557696
\]

where \(j\) is the first positive zero of \(J_0\). With unit end amplitude, the level point has \(U(L/2)=0.669929738985\) and \(\chi=1/2\). Its cycle-averaged longitudinal wave stress \(\chi(\rho_0\omega^2 U^2+\mu_0 U'^2)/4\) is positive there. The free edge has zero linear shear traction, and the inner node takes a mean outward force \(E/L\) for stored energy \(E\). The local stress across a cut in the strip and the force on the inner rim are separate readings. [Order-dependent shear and tested inertia law](conversion_layer_reflection.md) · [Bessel zeros](https://dlmf.nist.gov/10.21)

**Three-dimensional collar test.** Use the previous flare on \(2\le r\le3\), with \(R=c=1\) and slope parameter \(\epsilon=0.2\). Join it to a fixed node at \(r=1\) by multiplying the tilt by the smooth step \(G(r-1)\) on \(1<r<2\), where \(G(x)=10x^3-15x^4+6x^5\). Set \(\chi=1\) on this join and \(\chi=G(3-r)\) in conversion. Solve the scalar wave equation on the induced metric with stiffness and inertia proportional to \(\chi\), using the natural free edge at \(r=3\). The table gives the lowest axisymmetric mode and the integrated mean pressure on the inner node per unit wave energy. The positive intake from the geometric test is a matched background. Flow advection, losses and coupling to other modes are omitted. [Python spline solve](../wave_at_level_patch.py) · [Independent Mathematica solve](../wave_at_level_patch.wl)

| Tilt \(\beta\) | Upper-mouth layer count | Frequency \(\omega\) | Rim force per energy |
|---|---:|---:|---:|
| 0 | 1 | 0.650163 | 1.063 |
| \(5/2\) | \(1/2\) | 0.569861 | 0.897 |
| 4 | 0 | 0.498680 | 0.776 |

**Which way does the wave push?** At fixed wave action \(J=E/\omega\), the generalized force increasing tilt is \(-J\,\partial_\beta\omega\). A centered difference with half-step \(0.01\) at \(\beta=5/2\) gave \(+0.0481767478\,J\) with splines and \(+0.0481767402\,J\) with an independent polynomial basis. The respective rim loads per energy were \(0.8969553790\) and \(0.8971655580\). The sign of the tilt force includes the chosen join and the whole mode. To decide whether the physical throat resists this deformation requires the opposing work from wall tension, bending, the rim attachment and material flowing into the bulk along the same path. A positive energy barrier could stabilize a preparation. A universal electric charge would also require a mechanism that preserves the signed unit and couples it to the electric field. [Wall and rim terms](order_rim_radius.md) · [Electric-source test](mirrored_mouth_probe.md)

Working calculation on base e0a94594, uncommitted pending independent review.
