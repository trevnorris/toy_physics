# Which mouth response sets the 0PN carrier force?

Consider a flat, stationary, radius-$a$ opening in a neutral carrier. Assume the same nonzero field trace, constituent mass and pointwise normal current on both sides. Positive $x$ is toward the companion. Across this cut the carried sideways momentum passes through. The carrier force on the interface comes from the jump in mixed stress. With $N$ pointing from the feeding side into the opening, put $J=(\partial_N\ln n)^+-(\partial_N\ln n)^-$ and $K=\hbar^2/(4m)$. Then $F_x=K\int n\,\partial_x J\,d^3x$. A stationary external hold balances this force and any independent physical rim or other forces. [Reviewed transmission identity](../../../docs/held_drain_pair.md#what-passes-on-if-the-carrier-phase-is-continuous-at-the-mouth)

**What to measure.** For a spherical reference and weak incident velocity $U$, expand $n=n_0+U n_1\mu+\cdots$ and $J=J_0+U J_1\mu+\cdots$, with $\mu=x/r$. Only the dipole enters the linear integral. Angular integration and a radial integration by parts give

$$
F'_x(0)=\frac{4\pi K}{3}\left\{\int_0^a r^2\big[n_1 J_0'-n_0'J_1\big]\,dr +a^2 n_0(a)J_1(a)\right\}.
$$

The edge term vanishes for a smooth jump supported inside the opening. Given the rest matching $J_0$ and the measured density response, the one new carrier datum is the indicated weighted projection of $J_1$, including its edge value when applicable. A measurement needs the normal density slopes immediately on both sides at zero sideways flow and at opposite weak speeds, repeated at half speed. The archived slab ends at the feeding face. [Brane solver](../../../software/stage1_solver/tools/mouth_feed_weak_flow.py#L68)

**Can rest data fix this response?** Add a trial real surface energy per unit three-area $\delta e=\kappa h(\mathbf x-X\hat x)(n-n_0(\mathbf x-X\hat x))^2/2$, with $h$ a mouth weight vanishing smoothly at its edge. For any nonnegative $\kappa$, its energy, first density derivative and translation force vanish on the reference trace. It preserves the local matched rest slopes and normal current. For a given weak density trace the added force slope is $-\kappa\int h\,(\partial_U n)_0\,\partial_x n_0\,d^3x$. In the explicit spherical example $h=(1-r^2/a^2)^3$ inside, $n_0=1+b(r/a)^2$, $n_1=sr/a$, two calculations give $-128\pi a^2 b\kappa s/3465$. The comparison allows different local receiving slopes under drive. Global receiving states remain undetermined.

On the archived radius-one fields with imposed inward rate 0.25, the added odd force per speed per unit trial $\kappa$ was −1.69608024596 × 10⁻⁵ at $|U|=0.01$ and −1.69606371567 × 10⁻⁵ at 0.02. The negative-speed traces use the benchmark reflection symmetry. Independent signed translations of the surface energy reproduced both values. These are finite-speed sensitivities of an illustrative law. The physical surface response and any separate rim reaction remain open. [Model input](../../../software/em_charge_attribute/u1_body_dynamics_inputs.yaml#L56)

The [calculation](../../../software/stage1_solver/tools/mouth_stress_response.py) used archived fields from baseline 4b49fe61. [Literal results](force.stdout) · [Calibration feasibility](calibration.md)
