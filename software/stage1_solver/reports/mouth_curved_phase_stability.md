# Fixed-bath stability of the rounded-mouth carrier

At transported rate $Q=0.25$, all three represented fields admit a single continuous phase. Prescribing the complex carrier values $1$ and $e^{i\delta}$ at the distant baths, with $|\delta|<\pi$, selects a unique energy minimizer on the bounded continuum domain with the wall fixed.

## The stored fields

The lower bounds on $u=\operatorname{Re}(e^{-i\delta/2}\psi)$ hold over **every point** of every Q2 cell and are rounded downward. The two eigenvalues use the consistent mass matrix and the solver's four-point tensor quadrature, with zero perturbations on both baths. $L$ is the scalar operator below and $H$ is the real and imaginary energy Hessian.

| Mesh $n$ | Bath radius | Q2 cells | Lower bound on $u$ | Lowest eigenvalue of $L$ | Lowest eigenvalue of $H/2$ |
|---:|---:|---:|---:|---:|---:|
| 14 | 6 | 1036 | 0.9784777884093507 | 0.1386395339957829 | 0.13870391267854518 |
| 26 | 6 | 3562 | 0.9784775199315556 | 0.1386389093127812 | 0.1387032870384152 |
| 34 | 9 | 6018 | 0.9783459254334127 | 0.06146748781778277 | 0.06148555389221552 |

The independent check transforms the stored binary real and imaginary nodal values into exact tensor Bernstein coefficients. Rational alternating-series bounds enclose the phase rotation. The nonnegative Bernstein basis sums to one on the reference square, so the smallest enclosed coefficient bounds the field on the entire mapped cell. The eigenvalues are numerical estimates. Shifting below the sampled minimum potential targets the algebraically lowest modes. The largest relative eigenpair residual was $2.9277563222005586\times10^{-11}$. [Check script](../tools/mouth_curved_phase_stability.py) · [Carrier and mesh](../tools/mouth_curved_phase.py)

## The fixed-bath continuum energy

Write $n=|\psi|^2$, $W(n)=n^5/20-n/4$, $E[\psi]=\int (|\nabla\psi|^2/2+W(n))$, and $L=-\Delta/2+W'(n)$. A minimizer exists by the direct method on the bounded connected domain, with fixed traces, coercive growth in $H^1\cap L^{10}$, and almost-everywhere convergence along a subsequence. Projection onto the unit disk decreases the energy wherever the modulus exceeds one. Reflecting $u$ to $|u|$ preserves the energy and bath traces, so a minimizer can be chosen with $u\ge0$ and $W'(n)\le0$. In its weak equation $Lu=0$, testing with $(\cos(\delta/2)-u)_+$ gives $u\ge\cos(\delta/2)>0$.

For each bounded real perturbation $a$ vanishing at the baths, substitution of $a^2/u$ into that same equation gives

$$
\langle a,La\rangle=\tfrac12\int u^2|\nabla(a/u)|^2 .
$$

The reflecting walls give zero boundary term. Truncation extends positivity to finite-energy perturbations. Indeed the mixed-boundary Poincaré inequality gives a lower bound proportional to $\int a^2$, since $\cos(\delta/2)\le u\le1$ and $a/u$ vanishes at the baths. For a complex perturbation $\eta$, put $m=|\psi+\eta|^2$. Stationarity yields the exact change

$$
E[\psi+\eta]-E[\psi]
=\langle\operatorname{Re}\eta,L\operatorname{Re}\eta\rangle
+\langle\operatorname{Im}\eta,L\operatorname{Im}\eta\rangle
+\int [W(m)-W(n)-W'(n)(m-n)] .
$$

Convexity of $W$ on nonnegative densities makes the last integrand nonnegative. This proves uniqueness of the fixed-bath minimizer and strict energetic stability. The discrete eigenvalue check tests the actual Q2 equations independently of the continuum test function $a^2/u$, which need not lie in Q2. The calculation leaves relaxation, phase-slip dynamics, moving-wall stability, and the infinite-reservoir spectral gap open.
