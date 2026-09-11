# S11c-d response to the S10 Lean stratum-coverage concern

2026-09-11. Read-only inspection of the frozen end-resolvent engine and the
inventories committed at `718e5ced`. The current regular-case regeneration was
already running when the user supplied this concern. This note does not change
that computation, establish exceptional-domain coverage, or authorize an
upstream repair. The separate Lean trial is untouched.

The concern identifies real unfinished coverage, but the inspected root and
nullspace routines do not reduce each root set or eigenspace to one witness.
The current results concern complete finite root sets at specified parameter
bindings, followed by numerical mode checks. They do not cover the entire
parameter variety. Sampling each component or stratum would improve coverage;
one sample alone would still not certify an identity throughout that stratum.
Using a representative for a rank or related property requires establishing
that the property stays constant there, or identifying and handling its
exceptional subloci. The user's subsequent note from the S10 session agrees
with these boundaries and calls for preserving them without repeating the
per-root and full-basis checks already present. That session's edits to
`directives/S10_SHARED_PHYSICS.md` and
`scripts/S10_anisotropic_strata_comparator.py` are outside this construction
and are left untouched.

## 1. Native end roots: each isolated root and both normal lifts

[EndSpectrumCoverage.isolate](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2905)
uses square-free factors, enumerates their roots, constructs separate rational
disks with exact sufficient isolation bounds, checks pairwise disjointness,
and compares the multiplicity-weighted root count with the polynomial degree.
Square-free factorization is not an irreducible parameter-space decomposition;
the degree-count/disjoint-disk certificate is what covers this finite root set.

[The downstream constructor](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2962)
loops over every `(root, disk)` and both signs of the normal square root. It
evaluates the original rational physical matrix separately for each candidate,
then emits its singular values, left/right spaces, residuals, projector and
sheet record. The new [resolvent constructor](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3708)
likewise visits every native record. These are per-root checks at the bound
inputs, not a single evaluation of an aggregate root set.

## 2. Sheets and cut banks: explicit paths, incomplete region coverage

[BulkContinuationAudit](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3173)
tests positive/negative-frequency propagating and evanescent reference points,
upper/lower rays, local path ordering, winding loops, and separate banks at
three offsets from selected targets. Momentum-bank targets use unresolved
native candidates when available; otherwise they are computed branch-ray
probes. Both physical matrices are evaluated on each defined bank pair.

This is not an exhaustive decomposition of the joint complex frequency and
momentum domain into regions or path classes. A successful bank/path check
does not establish coverage on the rest of a cut or across every singularity.
All 24 branch-intersection controls have unique `(packet, label)` records with
their own vertices and unresolved status; they are not merely an aggregate
count. Their unresolved status records a failed transport domain, not a solved
threshold operator. The 48 unresolved native sheet labels also remain separate.

## 3. Degenerate spaces: full matrix bases, numerical rank limits

[FullSectorModes.solve_sample](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2479)
and the native constructor take every SVD column assigned to the nullspace.
For a two-dimensional space the right and left residuals are matrix residuals
on both columns. Frequency normalization and derivative pairings operate on
the full two-by-two systems. Oblique projectors use the full right/left bases
after checking that their overlap has rank equal to the computed nullity.
The rectangular jets also evolve matrix invariant pairs for the whole cluster.
The new Laurent pairing uses the same full-subspace construction.

This excludes the specific single-vector failure mode. It is not an exact
basis-completeness certificate independent of numerical rank: the SVD
threshold and singular values are emitted, and the projector has an
idempotence residual, but no separate exact projector-rank certificate is
computed. At `718e5ced`, the native dataset has 336 one-dimensional and 96
two-dimensional spaces across 432 candidates, all with recorded
algebraic-minus-geometric multiplicity zero. The older 528-candidate jet
dataset is a different census from these native candidates.

## 4. Exceptional strata: detected at chosen inputs, not enumerated

The native constructor computes exact gcds with the denominator norm, normal
threshold and radical-branch factor after parameter binding. It records
algebraic/geometric multiplicity differences and singular overlaps. The new
resolvent records singular derivative pairings, thresholds and unresolved
local contour separation when encountered. Its circles deliberately avoid
the computed branch/denominator loci.

None of those operations solves the exceptional parameter loci themselves.
The four case labels and PIT samples are not a census of threshold, root
coalescence, denominator intersections, modal-gap closures or Jordan strata.
Defective generalized modes and their resolvents remain unconstructed.

A concrete guard limitation also needs attention before introducing such
inputs: `REGULAR_ALGEBRAIC_MODE_COVERAGE` checks finite-root coverage,
denominator/normal-lift regularity and positive numerical nullity, but does
not explicitly require algebraic/geometric multiplicity agreement. Its name
must not be used as a defective-mode coverage certificate. This inspection
found no nonzero multiplicity difference in the committed native candidates;
the limitation is therefore not evidence that one of those candidate records
was silently omitted.

## 5. Continuum and bound poles: separate future constructions

[Shared physics section 3b](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:479)
requires conditional, profile-dependent bound-pole existence and distinguishes
it from continuum conversion. The retained solver/export contract already
requires bounded frequency searches, explicit sheets and isolating contours;
an unsuccessful search remains unresolved. The new local normal-momentum
residues do not implement that frequency-pole search or assert bound states.

Before any claim of full sheet/spectrum coverage or a section 3b result, the
remaining work must explicitly derive exceptional-locus conditions from the
operator, distinguish structural multiplicity from new coalescence, retain
their intersections and admissibility domains, and construct targeted
evaluations or unresolved records on those loci. It must also account for
the tested sheet regions/path classes and supply dedicated profile-frequency
pole searches and closure tests. Generic witnesses cannot discharge these
requirements. This is a dependency to plan at the next stop, not an expansion
silently inserted into the frozen regular-case run.

The [inspection evidence](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_stratum_coverage_evidence.json)
pins the engine and committed inventories, records the native nullity and
multiplicity census, and lists all 24 individual unresolved path controls.
No new physical premise, upstream numerical discrepancy, or completed
exceptional-domain certificate was established by this inspection.
