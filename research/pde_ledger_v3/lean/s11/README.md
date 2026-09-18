# S11 Lean contracts

The homogeneous contract below is the first bounded S11 contract under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Read
[COVERAGE.md](COVERAGE.md) for H1–H4 and [FIDELITY.md](FIDELITY.md) for the precise
claims, native action/operator correspondence and limits. Local verification
passed; see [VERIFICATION.txt](VERIFICATION.txt). Both independent fidelity
reviews returned CLEAR, and **H1–H4 are complete**. The
[closure record](FIDELITY_REVIEW.md) retains the reviewed revision and findings.

The supplied homogeneous action adds compression to the existing curl stiffness.
The proofs reuse S10's calculus and geometry, classify full amplitude kernels
including frequency coincidence, recover B=0 and handle k=0 separately. The
bulk threshold theorem is kinematic; interfaces, bound states, leakage and
nonuniform scattering remain separate work.

| Source | Responsibility |
|---|---|
| [Action.lean](S11Homogeneous/Action.lean) | Actual density, integrated variation, local PDE, phase-average normalization and modal operator. |
| [Spectrum.lean](S11Homogeneous/Spectrum.lean) | Exhaustive full kernels, D3 counts, positive frequencies, coefficient coincidence and boundary limits. |
| [Threshold.lean](S11Homogeneous/Threshold.lean) | Phase-matching identity and below/on/above-grazing sign classification. |
| [S11Homogeneous.lean](S11Homogeneous.lean) | Selected load-bearing axiom audits. |

From the parent `lean/` directory, the resource-conscious verification command is

```sh
python3 ../_measurements/S11_lean_contract_check.py
```

It builds only the three S11 modules and the audit root, sequentially, then runs
isolated controls. Logs live in `s11/_scratch/verification/`; the durable result
is `../_measurements/S11_lean_contract_checks.json`. `--reuse-build` requires
identical local transitive sources, dependency pins, commands and output-module
hashes. The 4096 MiB setting limits Lean's allocator, not OS resident memory.

The compact source check is
`python3 ../_measurements/S11_lean_source_check.py`. It executes the selected
native D3 MAIN constructor and two small modal routes, not the production audit.
No new installation or dependency update is needed.


## D2 invariant-space contract I1–I4

The separately authorized [invariant contract](INVARIANT_NEXT.md) proves the
complete quadratic invariant spaces under SO(2) and O(2) conjugation, with
dimensions 4 and 3 and a one-dimensional reflection-odd complement. Local
verification passes and both independent fidelity reviews returned CLEAR.
**I1–I4 are complete**; see the
[review and closure record](INVARIANT_FIDELITY_REVIEW.md).
Read [INVARIANT_FIDELITY.md](INVARIANT_FIDELITY.md) for the object conventions,
minimal native Q9 transpose correction and its exact D2 span check, and
[INVARIANT_VERIFICATION.txt](INVARIANT_VERIFICATION.txt) for the evidence.

From `lean/`, its separate sequential verification command is

```sh
python3 ../_measurements/S11_lean_invariant_contract_check.py
```

The compact source command is
`python3 ../_measurements/S11_lean_invariant_source_check.py`. These use separate
reports from H1–H4. The original orientation-probe script is historical and
expects the pre-repair native source. The new contract does not certify D3–D5,
EL/total-divergence classes, production exports, or S11c calculations.

## D2 odd-invariant dynamics E1–E4

The next bounded increment is [DYNAMICS_COVERAGE.md](DYNAMICS_COVERAGE.md), with
its precise object and normalization map in
[DYNAMICS_FIDELITY.md](DYNAMICS_FIDELITY.md). **E1–E4 is complete: local verification passed and both independent reviews
returned CLEAR.** See [DYNAMICS_VERIFICATION.txt](DYNAMICS_VERIFICATION.txt) and
[DYNAMICS_FIDELITY_REVIEW.md](DYNAMICS_FIDELITY_REVIEW.md). It addresses the actual bulk first variation of the D2
odd term and its exhaustive longitudinal/transverse mixing criterion.

`S11OddDynamics/Action.lean` identifies `-beta P/2`, its derivative-defined
momenta, local PDE and modal operator. `Variation.lean` derives the finite
relative-action variation using existing S10 calculus. `Mixing.lean` identifies
the zero loci and supplies a smooth background with nonzero compact-test first
variation whenever beta is nonzero. `S11OddDynamics.lean` audits the new theorems.

Separate commands, run sequentially from `lean/`, are

```sh
python3 ../_measurements/S11_lean_dynamics_source_check.py
python3 ../_measurements/S11_lean_dynamics_contract_check.py
```

The latter binds the nine unchanged local dependency sources to their compiled
objects, then checks the new modules, axioms and mutations with one Lean worker.
It supports hash-checked `--reuse-build`. These checks do not regenerate pinned
S11c exports or certify the complete XFORM_EXTRA spectrum.

## D3 invariant completeness J1–J4

The bounded contract is [D3_COVERAGE.md](D3_COVERAGE.md), with the exact object
and native span map in [D3_FIDELITY.md](D3_FIDELITY.md). **J1–J4 is complete:
local verification passed and both independent fidelity reviews returned
CLEAR.** See [D3_VERIFICATION.txt](D3_VERIFICATION.txt) and
[D3_FIDELITY_REVIEW.md](D3_FIDELITY_REVIEW.md).

`S11D3Invariants` classifies quadratic forms on all real 3×3 gradients under
the full SO/O conjugation action. Its proved census is 3/3/0, with the unique
trace-square, trace-of-square and Frobenius-square parameterization. This is a
density classification before any EL or total-divergence quotient.

The compact native check and one-worker formal checks run sequentially:

```sh
python3 ../_measurements/S11_lean_d3_source_check.py
python3 ../_measurements/S11_lean_d3_contract_check.py
```

The formal check supports hash-checked `--reuse-build`. It also checks that the
finite coordinate certificate matches its generator. Neither command reruns a
production audit or modifies S11c exports.

## D3 bulk variation K1–K4

The authorized increment [D3_BULK_COVERAGE.md](D3_BULK_COVERAGE.md) connects
the completed three-density classification to bulk equations. **K1–K4 is
complete:** local verification passes and both independent reviewers returned
CLEAR. See [D3_BULK_FIDELITY_REVIEW.md](D3_BULK_FIDELITY_REVIEW.md). The object and sign map
are recorded in [D3_BULK_FIDELITY.md](D3_BULK_FIDELITY.md), with current evidence
in [D3_BULK_VERIFICATION.txt](D3_BULK_VERIFICATION.txt).

`S11D3Bulk` reuses S10 calculus and the D3 classification. Its checked results
are the actual first variation, the one-dimensional null family, the two
independent bulk responses and their identification with homogeneous stiffness.
It includes an explicit divergence current and compact native Q9 V5 checks.
It does not add new spectral or interface calculations.

This is a stronger verification of existing mathematics and native CAS results,
not a new physical discovery. For
`Q = a (tr G)² + b tr(G²) + c tr(G Gᵀ)` and `L = -Q/2`, Lean proves that
the bulk operator is `(a+b) grad(div u) + c Delta u` for every smooth field.
Its entire null family is `c = 0, a+b = 0`: three independent densities give
exactly two independent bulk responses. The explicit divergence current
explains the missing response; it does not establish that boundary effects
vanish. The added assurance is completeness, checked conventions, meaningful
negative controls and two independent fidelity reviews.

## D4 invariant completeness

The bounded [D4_COVERAGE.md](D4_COVERAGE.md) contract is complete, with local
verification passed and both independent fidelity reviews CLEAR. Lean proves
the full SO(4)/O(4)/reflection-odd quadratic density classification, unique
coefficients, census 4/3/1 and complete even/odd split. Compact native checks
match the complete spaces and establish `P_D=P`, with the fully summed
epsilon contraction equal to `2P`. All fourteen mathematical mutations were
rejected and sixteen positives passed. See [D4_FIDELITY.md](D4_FIDELITY.md)
and [D4_VERIFICATION.txt](D4_VERIFICATION.txt); review provenance and optional
finding dispositions are in [D4_FIDELITY_REVIEW.md](D4_FIDELITY_REVIEW.md).
The divergence/zero bulk effect of the D4 odd term is covered by the separate
D4B increment below.

## D4 odd-density bulk variation

The bounded [D4_ODD_COVERAGE.md](D4_ODD_COVERAGE.md) increment is complete:
local verification passes and Claude and Grok independently returned CLEAR.
It proves the explicit divergence current and zero bulk first variation of
the constant-coefficient odd term, using the completed D4 classification.
The current retains its factor of one half; nonzero density and
momentum controls distinguish a bulk cancellation from pointwise vanishing.
Four modules and the audit root pass, with 41 standard-axiom audits, twelve
mathematical rejections, twelve positives and compact native correspondence.
See [D4_ODD_FIDELITY.md](D4_ODD_FIDELITY.md) and
[D4_ODD_VERIFICATION.txt](D4_ODD_VERIFICATION.txt), with review provenance and
finding dispositions in [D4_ODD_FIDELITY_REVIEW.md](D4_ODD_FIDELITY_REVIEW.md).
This is a statement about constant coefficients and compact variations;
possible boundary effects remain.

## Tail, Abel and inverse-stability estimates

The authorized [T1–T4 contract](ANALYTIC_ERROR_COVERAGE.md) addresses the first
question from the S11c calculation session. It proves conditional integral
tail/Abel and bounded-inverse perturbation estimates; operator-specific error
constants and an outgoing inverse margin remain application obligations.
The fresh recorded suite passed four modules plus the audit root, 41 standard-axiom
audits, twelve mathematical rejections and sixteen positives. The compact native
source check also passed. Claude and Grok independently cleared the bounded contract. Sixteen positive
executions represent fifteen distinct statements. Review findings and limits
are recorded in [ANALYTIC_ERROR_FIDELITY_REVIEW.md](ANALYTIC_ERROR_FIDELITY_REVIEW.md).
See [ANALYTIC_ERROR_FIDELITY.md](ANALYTIC_ERROR_FIDELITY.md)
for the compact native normalization link and the explicit limits of the claim.

## Variable coefficients and interface terms

The authorized [VC1–VC4 contract](VARIABLE_COEFFICIENT_COVERAGE.md) addresses
question 2 using the reviewed D3 family and D4 odd density. It proves the
coefficient-gradient corrections, weighted currents and an actual flat-interface
normal-slice integration identity. Local verification passes: 31 standard-axiom
audits, twelve mathematical rejections, sixteen positives and compact native
checks. Both independent fidelity reviews are CLEAR; the bounded contract is
complete. See [review and closure](VARIABLE_COEFFICIENT_FIDELITY_REVIEW.md) and
[the fidelity boundary](VARIABLE_COEFFICIENT_FIDELITY.md); no full S11c operator
or multidimensional transmission result is claimed.

## Nonlinear-pencil finite core NP1–NP4

VC1–VC4 is complete at `9865050f`. The user authorized proceeding to the
calculation session's third question after that checkpoint. The
[assessment](POLE_ASSESSMENT.md) reconciles the already adopted
`nonlinearPoleV2` correction with its analytic and physical application
obligations. [POLE_COVERAGE.md](POLE_COVERAGE.md) bounds this increment to typed
modal algebra, finite Laurent/circle integrals and discriminating actual small
pencils. [POLE_FIDELITY.md](POLE_FIDELITY.md) records the precise proof boundary.

Local verification passes: seven guarded canonical builds, 55 standard-axiom
audits, 17 mathematical rejections and 21 positives. The compact read-only
identification passes 19 checks and four controls. Claude and Grok independently
returned CLEAR; the bounded contract is complete. See
[review dispositions](POLE_FIDELITY_REVIEW.md),
[POLE_VERIFICATION.txt](POLE_VERIFICATION.txt), and the
[calculation-session handoff](POLE_HANDOFF.md). Run the separate
`../_measurements/S11_lean_pole_contract_check.py` from `lean/` for the one-worker
suite; `--reuse-build` requires its recorded source/command/object guards.
This increment does not rerun the old native diagnostics or S11c production,
prove analytic Fredholm theory, or certify a physical pole/scattering solve.
Its README and lakefile additions do not update the historical VC packet.

## Full D4 bulk equivalence (complete)

`S11D4Bulk` implements the bounded [D4C.1–D4C.4 contract](D4_BULK_COVERAGE.md):
actual variation of the full four-coefficient D4 density, exhaustive bulk
operator equivalence, the two-dimensional response and null spaces, and both
boundary currents. It reuses the completed D4 odd proof. Recorded verification
passes 67 standard-axiom audits, sixteen mathematical rejections and twenty
positive executions, plus compact native correspondence. Claude and Grok both
cleared the fixed packet; author closure validation confirms unchanged proofs,
objects and protected historical inputs. See
[D4_BULK_FIDELITY_REVIEW.md](D4_BULK_FIDELITY_REVIEW.md) and
[D4_BULK_FIDELITY.md](D4_BULK_FIDELITY.md). This adds no D5 or S11c calculation
scope and does not assert absence of boundary effects.

## D5 invariant completeness (complete)

The authorized [D5.1–D5.4 contract](D5_COVERAGE.md) targets the complete
SO(5)/O(5)/reflection-odd quadratic density classification on all real 5×5
gradients, with census 3/3/0. It uses the odd-dimensional conjugation identity,
an exhaustive finite quadratic representation and checked necessary rotation
constraints, with separate full-group sufficiency. See
[D5_FIDELITY.md](D5_FIDELITY.md) for conventions and the compact native span
check. Local verification passes: 45 guarded canonical objects, 51 standard-axiom
audits, thirteen mathematical rejections and seventeen positives. Native complete
spans/reflection and live-source/preservation checks pass. Claude and Grok both
returned CLEAR with no required fix; the bounded contract is complete. See
[D5_FIDELITY_REVIEW.md](D5_FIDELITY_REVIEW.md) and
[D5_VERIFICATION.txt](D5_VERIFICATION.txt).
This increment does not include D5 bulk variation or S11c calculations.
