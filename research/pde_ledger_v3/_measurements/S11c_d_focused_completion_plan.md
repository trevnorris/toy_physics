# S11c-d focused route to numerical physical results

This plan implements the user-approved 2026-09-18
[acceptance change](../directives/S11c_d_EXPLORATORY_ACCEPTANCE.md).
It supplements the source-pinned matching/numerical-action plans without editing
their frozen files. First complete one approved LAB_HELD / RHO4_CONSTANT example;
then perform required bookkeeping, controls, other cases and export integration.

## 0. Finish the active independent quadrature comparison

The wider-box triple adaptive production launched at `b5597574` continues
unchanged. Accept its actual result only after the original saved-operand,
worker, residual, metadata and source checks. Publish, annex-verify and commit
accepted evidence. If it fails, preserve partials and diagnose the issue without
an automatic full rerun. A physical disagreement still triggers the user's
repair checkpoint. Instrument fixes should reuse completed work.

This ends the planned fixed-box instrument/refinement sequence. Further such
work needs a material connection to the scattering observable, not merely a
smaller raw integral residual. No new day-long sweep is queued by this plan.

## 1. Build the first complete finite scattering solve

Use accepted both-end modes, full subspaces and currents, reduced local orders
0–3, all 80 nonlocal operands, nested profiles and approved inputs. The two
Gaussian action probes validate parts of the instrument; they do not span the
scattering solution or certify a new discretization. Construct actual boundary
matching and all required channel insertions from those operands.

The next implementation checkpoint is a small finite-domain discretization
that assembles the complete operator and boundary data and solves for both
incident ends. Record unknowns, basis/mesh, physical and numerical parameters,
coefficient/field/limit census, cost and memory, equation and boundary residuals,
and a numerical conditioning/sensitivity diagnostic. Include every open
reflected/transmitted channel and the required evanescent modes. Keep full
physical current matrices; do not assume conservation or invent absent channels.

Prototype on a small problem first. Reuse cached transforms, accepted groups or
other saved operands only when their actual arguments/settings apply; source
values sampled at old frequencies are not a uniform interpolation certificate.
An approximate transform, projection or alternative quadrature must be tested
at the arguments it consumes. Retain full symbolic coefficient dependence.

Deliver a scoped first numerical response with conversion and transverse
survival diagnostics. Explicitly distinguish a finite-contrast truncated-pencil
response from the required retained-order continuum expansion. Lack of a uniform
outgoing-inverse theorem no longer blocks this finite numerical construction.

## 2. Check the observables, with a bounded initial comparison set

Select one numerical resolution change, one domain change and one regulator
change around the baseline, reusing accepted evidence where applicable. These
are an initial diagnostic set, not a proof that any three comparisons suffice.
Use results to choose an additional targeted comparison only if needed; avoid a
Cartesian sweep of every order, cutoff and regulator. Include an affordable
independent route or physical control that tests the most consequential risk.

Before production, record the numerical precision target, units, near-zero
reporting rule, measured pilot cost and the claim each comparison could change.
Use the acceptance addendum's approximate 1% default only for resolved nonzero
observables. A small conversion, a cancellation or poor conditioning may require
more accuracy. Report numerical spread as empirical evidence, with every finite
cutoff/regulator retained. If sensitivity stays significant, identify its source
and reassess method/cost before another multi-hour batch.

The milestone is a numerically supported continuum response on the declared
instance after the required eta/sigma re-expansion, with enough precision for
its stated physical claims. Uniform norms and certified infinite limits are
optional follow-up. Unresolved numerical limits are explicit scope restrictions.

## 3. Compute the remaining physical outputs

Finish amplitude/flux baseline, interference and quadratic bookkeeping, weak
coefficients and required section 5 controls, using the computed response.
Carry out a separate bounded, profile-dependent frequency pole search with
targeted candidate refinements under nonlinearPoleV2. Track full principal parts
and actual source/observation coupling; do not infer capture rates, pole absence
from zero residue, or global spectral coverage from generic samples. Difficult
regions may be reported unresolved without blocking unrelated outputs, but an
unexecuted mandatory construction remains outstanding.

## 4. Integrate and finish

Integrate the numerical path and transparent symbolic continuum/weak objects
with the existing engine and own-row export contract. Run the remaining cases
with the validated method and practical per-observable tolerances. Preserve
case-specific channel/threshold/domain failures. Run required export closure,
semantic and output checks; publish via DataLad/git-annex and commit. Do not
repeat every historical development grid for every case by default.

## Deferred certification and operating rules

The rigorous T1–T4 operator norm/inverse application and global analytic spectral
certification are optional. Keep their hypotheses visible when a particular
formula or claim relies on them. New physical bugs still require a repair
discussion; the acceptance change does not relax equation correctness.

Each new multi-hour job gets a short question/cost/precision/stopping record.
Use the existing one-supervisor, up-to-four-worker silent completion workflow.
Save operands and recover validation/emission from them without physics reruns.
Commit each substantive step; do not request fresh permission for already
authorized routine work. No elapsed-time promise is made for unbuilt solves.
