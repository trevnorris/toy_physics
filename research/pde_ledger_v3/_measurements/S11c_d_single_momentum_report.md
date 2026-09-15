# S11c-d one-momentum refinement

The complete finite momentum result is published at 20fe9381 and annex-verified
at 621d81db. Its saved term decomposition identifies the one-momentum layout as
the source of the remaining 3.64e-5 / 1.05e-5 changes. The corresponding two-
and three-momentum contributions changed by at most 7.72e-9 and 2.23e-12.
The decomposition reproduces full-action differences to 2.90e-16.

The prepared calculation refines all 40 one-momentum rows on both approved
fields. Outer orders 64/96/144/216 and source orders 128/192/256 distinguish
momentum from source quadrature error. Independent adaptive Gauss-Kronrod outer
integration uses the same source operands. Complete actions explicitly retain
the accepted finest two/three-momentum terms; their convergence is not presumed.

Only the source Gauss node/weight rules are cached, with exact comparisons and
read-only arrays. Native numerical thread pools are limited to one. Frozen
sources and every completed numerical record remain in repository scratch.
Focused checks pass with exit zero and empty stderr: 128 point evaluations
reconstruct all 40 rows at three positions for both fields. Point-sum, accepted
layout and complete-action residuals are exactly zero. Both source-rule caches
match the original nodes/weights exactly; actual weight mutations respond.
No new refinement or adaptive result has been accepted yet. Finite-grid tests establish
no infinite-domain tail/interchange, Abel weak limit, scattering or pole result.
Approved inputs and the retained solver/export contract remain unchanged.
