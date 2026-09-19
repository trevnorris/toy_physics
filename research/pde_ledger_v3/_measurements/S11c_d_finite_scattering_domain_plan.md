# Finite boundary and regulator sensitivity of the scattering response

The selected resolution study is validated at `c1d466e8`. Its dominant response
changes are below the practical reporting target; small reflection and current
deficits remain unresolved. Test the more consequential finite-boundary and
regulator assumptions next, preserving the approved physical inputs.

Use three sequential cases, all with momentum16/4/4, source/profile256/512,
momentum cutoff4, profile cutoff14 and every native contribution:

1. At position/source bound48 and regulator0.2, increase97 to129 coefficients
   per field. This supplies a matching-rule baseline for the domain comparison.
2. Keep129 coefficients and regulator0.2, extending the position/source interval
   from[-48,48] to[-64,64]. Derive the source cutoff substitution from the accepted
   position-domain adapter, preserving all80 rows,35 generic sources and six
   profiles. Reevaluate every contribution; hold no nonlocal term.
3. At the larger interval and unchanged rules, reduce regulator0.2 to0.1.
   Recompute the source-derived Abel width and all actual regulator occurrences.
   This measures finite-regulator sensitivity, not an Abel limit.

Boundary-anchored complex amplitudes include the propagation phase to their
reference planes. Compare them at a common origin: derive the input/output
phase maps by evaluating the actual accepted modal plane waves at each boundary.
Test the maps against direct modal-vector evaluation and projection, including
a wrong-sign control. Normalize the common-grid fields to the same incident
phase. Retain full current matrices and verify the rephased quadratic forms.
Do not exponentially propagate evanescent amplitudes back to the origin: retain
them at their actual boundaries as domain-dependent diagnostics.

Generate the variable-domain constructor from the accepted finite constructor
with a reverse whole-function AST join: only the interval argument and the two
actual regulator bindings change. The quadrature, local/nonlocal coefficients,
boundary maps and solve remain identical. Derive cutoff rebinding and common-grid
inspection from their accepted routines with similarly explicit substitutions.
An unchanged-parameter comparison must reproduce saved matrices/outputs from
their actual operands; a small bounded prefix suffices for quadrature wiring.

Keep the prior1% target for resolved observables, amplitude reporting resolution
1e-4 and normalized-current resolution1e-6. Report actual complex-amplitude,
current and field changes, raw/scaled equations, boundary residuals, rank and
conditioning. Small deviations from unit survival are unresolved unless their
numerical spread supports the claimed precision; never force current conservation.

The final97-coefficient case previously took39.4 seconds of construction.
Increasing to129 and changing the finite domain/regulator is expected to take
several minutes for this set; this is an estimate. Bound total child computation
to900 seconds, one numerical thread and2 GiB per child. Save each complete case
and all unique batch partials. Do not restart completed cases if a later case
fails. No day-long sweep or rigorous infinite-limit certification is queued.

After this set, assess the actual observable spread before any further change.
Stability here supports a scoped finite numerical response; it does not prove
transparent nonlocal boundaries, resolve a tiny reflected/lost signal, supply
the required continuum grade expansion, or complete the pole/control/export work.
