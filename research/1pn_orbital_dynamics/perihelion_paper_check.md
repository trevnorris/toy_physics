# A paper check of the perihelion shift

Working note, shelved for the brane-feedback stage. The question is whether the older paper's proposed change in inertia would rotate an eccentric orbit, and what an independent calculation of the brane would have to supply.

Take a static spherical source and a test body. Let $\mu$ be the strength of the Newtonian field, $c_s$ the speed in the proposed correction, and $r$ the distance to the source. The paper assumes the following Lagrangian per unit reference mass. The dimensionless $\beta$ is left free throughout the orbit calculation. [Kinetic term in the original paper](paper/1pn_orbital_dynamics.tex#L91)

$$
L=\tfrac12\left(1+\frac{\beta\mu}{r c_s^2}\right)v^2+\frac{\mu}{r}.
$$

Write $\mathbf P=\partial L/\partial\mathbf v$ for the canonical momentum per unit reference mass. To first order in $1/c_s^2$, the Hamiltonian correction is $H_1=-\beta\mu |\mathbf P|^2/(2r c_s^2)$. It is already of first order, so evaluate it on the unperturbed ellipse of semimajor axis $a$ and eccentricity $0<e<1$. There $|\mathbf P|^2=\mu(2/r-1/a)$. A time average using eccentric anomaly $E$, $r=a(1-e\cos E)$, and $dt/T=(1-e\cos E)dE/(2\pi)$ gives

$$
\left\langle r^{-1}\right\rangle=\frac1a,\qquad
\left\langle r^{-2}\right\rangle=\frac1{a^2\sqrt{1-e^2}},\qquad
\langle H_1\rangle=-\frac{\beta\mu^2}{2c_s^2a^2}
  \left(\frac2{\sqrt{1-e^2}}-1\right).
$$

The angular-momentum action $J=\sqrt{\mu a(1-e^2)}$ is conjugate to the direction of perihelion $\varpi$. Differentiate the averaged energy with respect to $J$ at fixed $\Lambda=\sqrt{\mu a}$, and multiply by the Kepler period $T=2\pi\sqrt{a^3/\mu}$. The resulting rotation between successive closest approaches is

$$
\Delta\varpi
=T\left.\frac{\partial\langle H_1\rangle}{\partial J}\right|_{\Lambda}
=\frac{2\pi\beta\mu}{c_s^2 a(1-e^2)}
\quad\text{radians per orbit, to first post-Newtonian order.}
$$

Zero $\beta$ gives a closed Kepler ellipse. Positive $\beta$ advances perihelion and negative $\beta$ makes it regress. The assumed correction conserves energy and angular momentum. This orbit average reproduces the project's radial-equation derivation. [Original orbit calculation](mathematica/derive_precession_general.wl)

At $\mu=GM_\odot$ and $c_s=c$, the result for Mercury is $(\beta/3)\times 42.98$ arcseconds per century. The factor 42.98 is the leading solar relativistic contribution. The old paper obtained $\beta=3$ by matching the Schwarzschild result. Its conditional density and added-mass estimates sum to $1+1/2$, and it assigned the remaining $3/2$ by that match. A prediction from moving drains still requires a brane calculation that yields the correction and its coefficient independently of the orbit data. The linear 0PN response at release supplies neither. [S17 coefficient record](../pde_ledger_v3/steps/S17_post_newtonian_ladder.md#L23) · [Clifford Will, Eq. 2 and Mercury value](https://arxiv.org/pdf/1802.05304)
