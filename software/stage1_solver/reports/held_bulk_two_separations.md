# Measuring the receiving-bulk dipole at two held separations

**Result.** Stationary angular velocity fields recover the separation slope of the receiving-bulk dipole in the 0PN potential exterior. For two equal mouths with uniform outward normal velocity imposed over each entire spherical opening, the extrapolated slope at \(d=4a\) is \(+0.00395231085946\,v_0a^3\). For two transparent radial sources with the same radii, intakes and far conditions, the slope vanishes. The angular discharge condition selects the value. The physical mouth's condition is open. [Earlier energy test](receiving_bulk_energy_response.md) · [Scripts and literal readings](../../../docs/_measurements/held_bulk_two_separations.json)

**Preparation.** Both centers are at rest in every solve. Set their common radius to \(a\), keep their common intake \(q=D(a)>0\) fixed, and release that mass into an unbounded, constant-density 4D receiver with fluid at rest at infinity. Units are \(a\) for length and \(v_0=q/(2\pi^2\rho a^3)\) for velocity. The transparent-source control superposes the two radial potentials. It permits the normal velocity to vary over either enclosing sphere. The second calculation imposes total normal velocity \(v_0\) everywhere on each mouth. The needed angular correction is solved for both mouths together. Both calculations have the same integrated discharge. Uniform local discharge is an added assumption.

**Readout.** On a sphere centered on the target let \(\mu\) be the direction cosine toward the companion, \(v_N\) the outward velocity, and \(M(R)=4\langle\mu v_N\rangle_{S^3}\) the full angular average with that weight. In a source-free harmonic annulus, the axial dipole potential is \((G r+B/r^3)\mu\). Higher angular modes integrate to zero and \(\langle\mu^2\rangle=1/4\). Two measuring radii in the same field therefore give

\[
M(R)=G-\frac{3B}{R^4},\qquad
B(d)=\frac{M(R_o)-M(R_i)}{3(R_i^{-4}-R_o^{-4})}.
\]

I used \(R_i=1.2a\) and \(R_o=1.6a\). The subtraction cancels the companion's regular dipole \(G\), including its own boundary response. A reading at just the near and far poles can mix higher odd angular modes into the result. This construction assumes a full spherical annulus outside all conversion regions, with harmonic flow there.

**Measured uniform-discharge branch.** \(B\) below is in \(v_0 a^4\) and the secant in \(v_0a^3\). Negative \(B\) means the induced dipole adds outward normal velocity on the companion-facing side relative to the regular incident field. The induced component grew stronger as separation decreased.

| Held separations \(d/a\) | \(B\) at smaller separation | \(B\) at larger separation | Secant |
|---|---:|---:|---:|
| 3.900 and 4.100 | −0.00564929490686 | −0.00485712295247 | +0.00396085977194 |
| 3.950 and 4.050 | −0.00543588220197 | −0.00504043761011 | +0.00395444591868 |
| 3.975 and 4.025 | −0.00533322461667 | −0.00513558238546 | +0.00395284462426 |

Two separations give a secant. Halving the separation interval twice tested the smooth midpoint limit. Extrapolating the last two secants gave the result above. Independently differentiating the Mathematica boundary equations at \(4a\) gave approximately \(+0.00395231105199\,v_0a^3\). On the finest angular solve the largest error in either normalized total intake was \(8.882\times10^{-16}\), and the largest departure from the imposed normal speed was \(5.329\times10^{-14}\) in velocity units. The transparent-source control retained a separation-dependent angular flow signal on each cut but yielded \(B=0\) to roundoff after the two-radius subtraction. [Guarded checks](../../../docs/_measurements/held_bulk_two_separations.json)

**Physical requirement.** The saved paired fields at separations 4 and 8 contain brane-side data. The available finite receiver uses one mouth under an acceleration drive. A physical determination needs matched receiving-side angular velocities at nearby held separations with the same physical mouth preparation, or an independently derived mouth condition that selects those fields. Static differences determine a quasistatic response. If the receiver has memory, a moving-history response is also needed for the energy test. The earlier energy example assumes \(B(d_*)=0\). The uniform-discharge branch has nonzero \(B(4a)\), so substituting this slope alone into that example would omit other terms. [Paired archive](../../../docs/held_drain_pair.md#what-can-be-read-for-the-two-held-unequal-mouths) · [Receiver preparation](../../../docs/_measurements/open_sleeve_lateral_balance.txt)

Working calculation at base cda9a8eb, pending independent review.
