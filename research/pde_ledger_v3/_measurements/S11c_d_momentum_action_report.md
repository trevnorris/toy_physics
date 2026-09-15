# S11c-d finite-action momentum quadrature

The accepted finite source quadrature is published at 36dce2e2 and its annex
verification is committed at daf404a7. This stage consumes those bound source
amplitudes and frequencies, all 80 factorized native integrals, all six nested
profile operands and all three source-derived Abel transfer pairs.

The implementation streams the remaining momentum quadrature in bounded
batches. It evaluates the source and profile integrals directly at requested
momenta and preserves native limit order. Both saved whole-action grids are
checked first, covering both Gaussian fields, three positions, 300 components
and every original nonlocal term. Three concentration-aware momentum grids
then record changes at fixed source/profile orders, finite bounds and regulator.

Every bound row and completed layout/grid is saved. Partial numerical sums are
saved every 64 batches. Actual operand metadata is exercised before expensive
integration; complete emission replay and source/packet hashes follow it.
The 32 MiB phase and cache budgets are recorded separately from process RSS.

Focused checks pass with exit zero and empty stderr: all 70 bound source
evaluations agree with their literal operands within 2.30e-14 scaled error; all
six profile evaluations match the native evaluator. Three Abel panel tests,
including a peak near the finite boundary, agree with the saved exact primitive
within 1.20e-14. Omitted panels and changed weights produce nonzero responses.
The complete accepted engine AST joins after removing only the new helper.

No full-action momentum result has been accepted yet. Finite refinement alone
establishes no infinite-domain tails/interchange, Abel weak limit, two-ended
matching, scattering or pole solve. Approved inputs and the retained solver/
export contract remain in force.
