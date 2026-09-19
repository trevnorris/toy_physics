# Local finite-pencil contour diagnostic

The complete 16-point contour passed in **517.42 seconds (8m37s)**, with four
clean single-thread workers and empty stderr. All 16 matrices have full rank
645. Every point includes all 80 rows, 160 terms, 35 sources and four incoming
columns; 2,845,568 new momentum nodes were evaluated. The native rules and
physical input remain unchanged. Worker peak RSS was at most 669,808 KiB.

On the circle centred at 1-0.01i with radius 0.02, both sampled 8/16 determinant
windings are zero. The maximum 16-node phase increment is **2.70022 radians**;
this still motivates one midpoint refinement. Zero sampled winding is not a
certified empty spectrum, as the saved degree16 alias control demonstrates.

Full inverse moment norms (orders 0–3) decrease from
8.73872e-2, 5.75245e-3, 3.21324e-4, 1.48556e-5 at 8 points to
7.02241e-6, 3.72828e-7, 1.83022e-8, 9.72474e-10 at 16 points.
The actual open-response moment norms at 16 points are at most 9.23198e-14.
Moments use the fixed recorded coefficient frames and are not assumed to be
projectors or simple-pole residues. Small source-to-observation moments alone
do not establish absence of an inverse pole.

The maximum sampled matrix condition is 24,087.7 and minimum singular value
8.97099e-5 in the fixed seed frames. Ten complete end clusters close numerically
within 4.25579e-13, while all frequency-dependent forcing/observation and phase
maps remain in construction. Complex-frequency amplitudes are analytic
coordinates, not physical gain/loss measurements.

Saved-operand acceptance checks every source/input/worker/artifact hash, all
16 source and 80-row censuses, actual measures and direct actions, independent
solves and both inverse identities, end maps and full moment arrays. Complete
original packets and logs remain preserved. No quadrature or solve was repeated
for acceptance. The next selected comparison adds only the 16 angular midpoints
for a 32-point contour. No physical pole set, principal-part classification,
global completeness or certified empty spectrum has been computed.
