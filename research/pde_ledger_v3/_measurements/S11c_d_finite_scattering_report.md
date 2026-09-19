# S11c-d finite scattering pilot

Implementation assembles all five fields and80 nonlocal rows from accepted
source factors, including the actual source derivatives in all35 generic source
amplitudes. It forms all four incident columns with the original physical current
matrices and full outward bases. Exact local operators and full symbolic source
operands are retained; the numerical adapter changes no upstream physics.

Focused boundary checks find five independent outward directions at each end:
two propagating transverse directions and three evanescent thickness directions.
Both trace-basis condition numbers are below5.94. Modal value/derivative map
residuals are below4.47e-16. Independent polynomial checks through derivative
order3 differ by at most3.56e-15. These checks validate the finite boundary
instrument; they do not prove a transparent nonlocal boundary condition.

Implementation is committed at `efc61982`. The325-unknown pilot has launched with a15-minute budget, one native numerical thread
and2 GiB ceiling. Save complete operator matrices before solving all incoming
columns together. Compare matrix actions with separate native-factor contractions
and actual measure mutations. Inspect residuals, rank, conditioning and amplitudes
before choosing more resolution. No physical scattering or continuum-grade
expansion is accepted yet. The regulator, source/momentum cutoffs, modest
quadrature and approximate modal boundary remain explicit numerical choices.

The completed adaptive sequence is published/annex-verified at4ba95bcf/a6bfb7a6.
Future work follows exploratoryAcceptanceV1 and the focused completion plan;
no broad tail/regulator refinement campaign is queued.
