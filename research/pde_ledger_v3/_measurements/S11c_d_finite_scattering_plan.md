# First finite S11c-d scattering construction

Under `exploratoryAcceptanceV1`, the next result is a complete finite numerical
response on the approved LAB_HELD/RHO4_CONSTANT example. The independent
quadrature sequence is accepted at `a6bfb7a6`; do not repeat it.

Use a polynomial collocation representation for all five fields on a stated
finite interval. Assemble the four accepted local derivative matrices and all
80 source-factorized nonlocal rows, retaining their actual input field, equation,
source/profile factors, momentum limits and measures. Form numerical matrices
by applying the computed operator to the chosen basis; no fitted potential or
typed mode-conversion coefficient substitutes for these rows. Keep inherited
equation/field units and the bound eta/sigma instance explicit.

At each end select the computed physical-sheet outward-decaying candidates
and current-defined outgoing propagating modes. Use their full computed bases,
including all degenerate directions. If they span the five field traces,
compute a modal value-to-derivative boundary map and its reconstruction residual;
derive the inhomogeneous boundary data by inserting every incoming mode. This
is an approximate finite modal radiation boundary condition, not a proved
transparent boundary for the full nonlocal operator or its cut continuum.
Do not assume that an incomplete/singular trace basis is invertible. Retain
every excluded candidate and the selection reason.

The initial pilot asks whether this full finite problem can be assembled and
solved affordably and whether its residuals/conditioning permit further use.
Choose one modest collocation and momentum rule, with explicit source/profile
rules and regulator. Measure construction time and memory before launching any
larger batch. Include the finite source cutoff in the model approximation;
outgoing exterior and branch-continuum omissions are boundary/domain uncertainty,
not certified zero tails. Save matrices before the solve, then solve all incident
columns together. Retain raw residuals, scaled residuals, rank/singular values,
modal amplitudes and all physical current matrices. No automatic unitary-flux
condition is imposed on this potentially non-Hermitian retained operator.

Initial pilot:65 Chebyshev coefficients per field (325 unknowns), four incoming
columns, interval/source bound48, momentum bound4, regulator0.2, outer8 and
concentration-aware inner2/2 momentum rules, source/profile128. This deliberately
modest momentum/source/profile rule is instrument/cost evidence, not the accepted
fine quadrature. Use one native thread,2 GiB address-space ceiling and a900-second
wall budget; save each layout and partial matrix. No total runtime is extrapolated
from earlier integral-only jobs. Inspect this pilot before any production choice.

Check polynomial differentiation and source-basis evaluation independently,
actual channel reconstruction and input directions, and the finite matrix action
against a separate contraction of its computed term operands. Preserve all
nonlocal terms, with a measure mutation capable of changing their action. These
are focused implementation checks, not another exhaustive quadrature campaign.

The pilot is not accepted physical scattering or the required continuum-grade
expansion. Inspect its result before selecting one practical resolution change
and targeted boundary/regulator checks. Report open-channel conversion only for
actually available channels; evanescent thickness fields are near-field response,
not an invented outgoing flux channel. Poor conditioning or a material boundary
inconsistency is recorded and resolved before physical interpretation.

Use saved source operands, new-stage authority/source hashes, durable scratch,
one native numerical thread for the initial pilot and a completion/error hook
if needed. Stop the pilot on its explicit time/memory budget and preserve partial
matrices; an incomplete pilot remains unfinished. The budget is an instrument
cost test, not permission to report missing entries as zero. Commit each usable
checkpoint and continue toward the complete numerical response under the focused
completion plan; only a distinct upstream physical repair requires user input.
