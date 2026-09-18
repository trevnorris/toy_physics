# S11c-d acceptance for the spare-time toy model

Status: user approved on 2026-09-18 (`exploratoryAcceptanceV1`). The user asked
to refocus future work on the proposed practical numerical standard and to let
existing runs finish. This is a change in completion criteria and priorities,
not a physical repair or permission to change the model.

## Authority and scope

Apply this document to future S11c-d planning, implementation and reporting.
It supersedes requirements to establish rigorous full-operator tail, Abel-limit
and outgoing-inverse certificates **before attempting numerical scattering**.
That includes the prerequisite language in numerical-action plan step 5,
quadrature-limits plan step 5 and the matching plan's convergence interpretation.
The retained builder contract remains verbatim historical authority, interpreted
with this acceptance addendum and `S11c_d_NONLINEAR_POLE_CONTRACT.md`.

The complete program's physical outputs, reduced equations, full coupling
vertex, approved profiles and parameters, units, independent grades, source
identities, actual boundary/channel construction and corrected nonlinear-pole
mathematics remain required. Numerical evidence does not turn a false identity
into a true one. No absent integral, mode, pole search, expansion or export may
be replaced by a fingerprint, a typed answer, zero or an empty set.

## What is sufficient for a numerical result

1. Compute the complete two-ended response on the approved development example,
   including every required open channel and evanescent matching contribution.
   Derive it from the accepted reduced local and nonlocal operators. A finite
   discretized solve is initially a scoped numerical result; retain the separate
   continuum expansion required by the original physics and export contract.
2. Check equations, boundary matching, current normalization and relevant
   controls. Use existing accepted comparisons where their operands and settings
   apply. New algorithms need focused checks that can detect sign, measure,
   indexing or omitted-term mistakes; unchanged machinery does not need another
   exhaustive instrument campaign.
3. Assess the actual reported observables with a small selected set of resolution,
   domain and regulator changes. Vary these independently before attributing an
   effect. Reuse prior evidence as context, not as a scattering error bound from
   two Gaussian action witnesses. Include a different numerical route where it
   is informative and affordable, targeting the dominant uncertainty.
4. Set the intended reporting precision before the acceptance comparisons. A
   working default is roughly two significant figures (about 1% stability) for
   well-resolved nonzero observables. This is a numerical precision goal, not
   the withheld physical falsification threshold. Record absolute tolerances
   in the declared unit frame too. Near zero, compare with the measured numerical
   uncertainty and state a resolution limit; a relative test with an arbitrary
   denominator floor cannot establish a tiny conversion, a sign or an absence.
   Small effects that carry the physical claim need tighter targeted checks.
5. Report observed numerical spread, residuals, conditioning, parameter/domain
   restrictions and any untested limit. Do not call a refinement difference a
   rigorous upper error bound. A finite regulator result stays labelled with
   that regulator until numerical limit evidence supports a stronger statement.
   An unstable feature stays unresolved; continue independent useful work rather
   than silently lowering its tolerance or declaring all of S11c-d complete.

Pole work remains a separate targeted bounded search. Preserve sheet, decay,
width, channel-closure, multiplicity, complete Laurent and response-map checks
where they support the claim. A resolved numerical candidate is not a general
meromorphic-operator theorem. Zero residue does not exclude a higher-order pole;
unsuccessful searches do not establish an empty spectrum. Thresholds and
coalescences require explicit treatment if encountered or claimed, not a global
survey before obtaining an ordinary-parameter scattering result.

## What becomes optional follow-up

Uniform estimates on trial/test unit balls, complete Banach-space operator-error
budgets, rigorous outgoing-inverse bounds, formal infinite-tail/Abel proofs,
global exceptional-locus coverage and general analytic Fredholm/Riesz premises
are optional certification work. Their absence limits the strength and domain
of a claim; it is no longer a blanket gate on a numerical calculation. Lean
handoffs remain useful guides and correctness checks. A theorem or certified
claim still requires its actual hypotheses; no such claim is manufactured here.

All originally required computed outputs, controls, continuum bookkeeping and
the final four-case/export integration remain program work. Do not silently
delete them while deferring the additional certification layer.

## Effort and stopping rules

Prioritize the first complete scattering result, then its important observables.
Stop refining an intermediate contribution once its uncertainty is demonstrably
irrelevant at the chosen reporting precision; a small action change alone does
not establish that after a possibly ill-conditioned solve. Record the numerical
sensitivity or an appropriate estimate before making that decision.

For each new multi-hour job, first record the physical question it resolves,
the observable and precision at stake, reusable evidence, cheaper alternatives,
an initial measured cost range, resource limits and a stopping decision. Use
one selected comparison batch and inspect it before extending the sweep. If
uncertainty remains material, target its source or change the numerical method;
do not automatically repeat day-long grids in successive directions. A planning
estimate is not a hard deadline or authority to discard unfinished work. Routine
choices remain autonomous under the user's standing authorization.

## Transition from the active run

Let the current wider-box triple adaptive run finish under its original source
pins, settings and guards. Preserve and validate its saved results; publish any
accepted transcript through DataLad/git-annex and verify hashes. A failure is
retained as evidence and handled from saved operands, not reclassified as a pass.
Do not stop/restart its workers or change their tolerances for this policy.

Afterward follow `../_measurements/S11c_d_focused_completion_plan.md` rather than
automatically launching more quadrature, box or regulator sweeps. The existing
watcher reads the execution/report checkpoints; its loaded notification text
and original message hash remain intact. New stages must pin this addendum and
their own current plan alongside unchanged consumed physical sources. Preserve
old frozen plans and packets; never substitute a changed source hash to rebase
them. Keep durable scratch, silent local completion/error hooks, at most four
authorized single-thread workers, stepwise commits and annex publication.
