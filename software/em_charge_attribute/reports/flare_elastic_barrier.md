# Can wall elasticity protect the signed mouth count?

**A finite barrier appears on the tested deformation for some preparations.** Most importantly, the bending cost comes from the region joining the flare to the fixed mouth. The converting part of the flare actually reduces its bending energy as the count unwinds. A freer attachment could change the result. [Previous unwinding geometry](mouth_count_unwinding.md) · [Energy curve for one preparation](flare_elastic_barrier.png)

Use the same radius \(R\), tilt \(\beta\), inner join and order-weighted wave as in the standing-wave test. The upper-mouth count is one through \(\beta=2\), half at \(2.5\), and zero at \(3\). For constant effective wall tension \(\sigma>0\) and bending modulus \(\kappa_b\ge0\), write the change from the untilted collar as

\[
 \frac{\Delta E}{\sigma R^3}=\Delta A+k\,\Delta B+j\,\Delta\hat\omega,\qquad
 k=\frac{\kappa_b}{\sigma R^2},\quad j=\frac{Jc}{\sigma R^4}.
\]

\(A\) is the wall area divided by \(R^3\), \(B=R^{-1}\int H^2dA_{\rm phys}\) uses the average of all three principal curvatures, and \(\hat\omega=\omega R/c\). \(J\) is held fixed during a tilt. The wall acts over the \(R\le r\le3R\) collar. The model supplies no matched effective moduli or receiver attachment for this draining preparation. Work beyond the cut, material-flow work and dissipation remain open. [Wall functional](order_rim_radius.md) · [Wave ensemble and geometry](wave_at_level_patch.md)

| Stage | Tension coefficient \(\Delta A\) | Bending coefficient \(\Delta B\) | Wave coefficient \(\Delta\hat\omega\) |
|---|---:|---:|---:|
| Count starts to fall | 9.342621 | 4.503866 | −0.056657 |
| Half the count remains | 14.224334 | 5.840227 | −0.080302 |
| Count reaches zero | 19.913330 | 6.980132 | −0.104489 |

At complete loss the bending change splits into \(+7.417089\) in the join and \(-0.436957\) in conversion. Along this path the zero-count endpoint costs energy for \(j<190.58+66.80k\) (rounded). A higher \(j\) can still leave an intervening barrier. At the halfway point the *local* restoring force dominates for \(j<220.12+51.39k\). Python integrated the graph Hessian and Mathematica took the divergence of the normalized gradient. Their wall coefficients agree to the displayed precision. [Graph area variation](https://web.stanford.edu/~ochodosh/GenericLondon2023.pdf) · [Run manifests and literal readings](../../../docs/_measurements/flare_elastic_barrier.json)

**An example with a lower-energy endpoint.** As a conditional way to select the action, require this isolated three-term energy to be stationary under uniform dilation at zero tilt, scaling the entire collar and conversion width together. This gives \(j=(3A_0+kB_0)/\hat\omega_0\). For the illustrative stiffness \(k=5\), \(j=582.37187\) on the fine mesh. The energy rises by \(0.418529\,\sigma R^3\) near \(\beta=1.10115\), before any count is lost. It ends \(6.03733\,\sigma R^3\) below the initial energy. At \(k=4\) the same dilation prescription gave no positive energy change at the sampled tilts through complete loss. Uniform dilation is a property of this isolated trial family. The physical open mouth can exchange work with its surroundings and need not obey that prescription. [Peak refinement and figure source](../plot_flare_elastic_barrier.py)

**Implication for charge.** The example gives energetic resistance along one path. The minimum barrier over allowed deformations, its size at different physical radii, and the coupling to the electric field remain unknown. The next discriminator is to relax the shape of the join at the same neck, intake and material parameters and check how much barrier survives. [Electric probe](mirrored_mouth_probe.md)

Working calculation on base e0a94594, uncommitted pending independent review.
