# Independent source calculation

The source review preceded receipt of the proposed calculation. The aperture is a three-ball in a flat face. Use a common normal from feed to receiver. With a continuous differentiable single-carrier phase, equal constituent mass, and pointwise continuous normal number current, carriage cancels between the two faces. The feed gives the incoming load $mJ(U)U$. Define the signed jump $b=(\partial_N\log n)^+-(\partial_N\log n)^-$ and $k=\hbar^2/(4m)$. The mixed-stress drive in the positive incident-flow direction is $k\int n\,\partial_x b\,d^3x$. A total support reading has the opposite sign to the total drive.

Let $n=n_0+U\mu p+o(U)$, $b=b_0+U\mu t+o(U)$, and $\mu=x/r$. On a sphere the mean of $\mu^2$ is $1/3$, and $\partial_x(\mu t)=\mu^2t'+(1-\mu^2)t/r$. Integration gives

$$
 F'_{\rm mixed}={4\pi k\over3}\left([r^2n_0t]_0^a+
 \int_0^a r^2(p b_0'-n_0't)\,dr\right).
$$

The origin term vanishes by regularity. The upper endpoint can survive at a finite radius. A smooth taper to zero for the dipolar jump removes it. Indicators, singular rim stresses, and moving-domain terms require an explicit edge balance. This radial endpoint differs from a sleeve's axial end traction. The latter can contribute separately to the mechanical support balance.

For the stipulated local conservative translated-density energy $e(r,n)$, fixed coefficients, one trace, and no additional arguments, two-sided variation gives $k b=e_n$. Then $k t=e_{nn}p$ and $k b_0'=e_{nr}+e_{nn}n_0'$. With zero radial endpoint the fixed-field translation gives $F'=(4\pi/3)\int r^2p e_{nr}dr$, where the radial partial holds density fixed. A force or independent coupling already represented by this stress jump counts once.

At each radius two nearby spherically symmetric preparations of the same material law can measure a local secant in $n$ when both two-sided jumps are observed and the two local densities differ. A differentiable law yields its tangent as these trace densities converge. At a fixed finite separation the derivative at either preparation generally remains unknown. Other uncontrolled changing arguments of an actual surface law can enter the measured secant. Neither the weak density dipole under signed sideways flow nor independent rim/receiving loads follow from rest data alone. The archived boundary laws prescribe the feeding face and have no independent receiving jump or support reading.

The script `independent.py` integrates Cartesian monomials on the three-ball and checks variation of an explicit energy at fixed carrier field. In units $k=1$ with slopes divided by $4\pi$, a nontapered test gives $24/35=2/3+2/105$ (endpoint plus interior). Two tapered conservative laws have identical slope jumps at both rest densities but give drive slopes differing by $16/17325$. Each equals its own direct energy translation. Literal stdout is in `witness-run/stdout`; `witness-run/limit-validation.json` records 1 GiB and zero swap. No solver run was used.
