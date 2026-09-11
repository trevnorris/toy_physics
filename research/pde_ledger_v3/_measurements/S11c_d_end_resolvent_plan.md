# S11c-d end-resolvent construction plan

2026-09-11. Begin after user-requested checkpoint `718e5ced` (five transcripts
saved through DataLad/git-annex; twelve other files through Git). The Lean trial
belongs to a separate user session and is outside this build.

Authorities: the S11c-d program brief, shared physics sections 2/3a/3b/6, the
build directive, and the retained solver/export contract. Continue build, run,
report, stop. No review legs, comparator, Wolfram engine or downstream stage.

## Objective and dependency

Construct the source-to-field inverse of each computed constant-end physical
pencil on its continued bulk-radical curve. The existing bank records evaluate
the operator; this step computes the resolvent, regular normal-momentum pole
residues and cut-bank resolvent jumps needed for subsequent Green-function and
variable-profile matching work. These normal-momentum poles at fixed real
frequency are distinct from the profile-dependent frequency poles of section
3b. Do not export them as bound states, Riesz data or scattering channels.

## Execution

1. Develop on the pinned LAB_HELD/RHO4_CONSTANT symbols and native spectrum from
   the completed joint-sheet run. Bind the actual reduced five-field operator
   and derive its total normal-momentum derivative on the radical relation.
   Preserve import wiring, Fourier reduction and every existing emitted result.
2. Evaluate the inverse on both existing banks. Emit each inverse, operator
   provenance, left/right inversion residuals, coefficient-frame conditioning
   and the inverse-jump identity residual with both operands. Restore inverse
   dimensions from the source equation and field dimensions.
3. For each native algebraic candidate, compute the derivative pairing of its
   actual left/right nullspaces and solve the Laurent cancellation equations
   for the regular pole residue. Emit singular pairings and exceptional domains
   explicitly. Use closed normal-momentum contours on the candidate's local
   radical lift to compute residues independently, with two radii and two node
   counts. Record branch/other-root/denominator clearance, matrix conditioning,
   pole-count diagnostics, Laurent moments and refinement residuals. Numerical
   clearance is not an exact contour certificate. Retain all original sheet
   labels; a local algebraic lift does not choose a global physical sheet.
4. At cut-bank targets coinciding with native normal-momentum poles, retain the
   raw inverse jump and explicitly compute the appropriate local pole
   subtraction on each bank. Keep regular-part refinement and any unresolved
   subtraction domains visible. Finite-offset samples alone do not establish
   a limiting continuum spectral measure.
5. Integrate this construction into reference/left/right, physical-input and
   algebraic-PIT packets for all four cases. Run the full engine once after
   focused checks. Inventory new objects and preservation of the prior
   continuation, root and inverse-Fourier evidence. Atomically publish completed
   outputs under scripts/out without writing through annex pointers. Report
   source hashes, actual child resources, residuals, dimensions and limits.

Heavy matrices use fingerprints of already computed objects and whole-object
SHA digests. Scalar diagnostics and small residual matrices remain inspectable.
The source grades, background origin, declared L/T/M frame and retained-model
status accompany evaluated records. No placeholder export is permitted.

## Completion and stop boundary

This is the constant-end resolvent dependency within the spectrum/continuation
program. It does not solve the complete variable-profile resolvent, its contour
pinches, a global continuum measure, energy-current normalization or scattering.
Exceptional/defective modes remain explicit domain limitations. Broader TODOs
are removed only when their full construction is actually complete.

If the calculation contradicts an earlier operand or requires a new physical
premise, an upstream repair or a change to the agreed approach/work breakdown,
preserve the evidence and stop to explain it before making that change. The
supplied-input, shear-normalization and c2 cross-engine operand/sign debts
remain open. Preserve the eight-point solver/export contract byte-for-byte.

## Focused implementation choices

The first physical/PIT contour checks evaluated all 18 native candidates per
packet and agreed with the modal Laurent residues. The closest PIT bank
samples exposed cancellation in double-precision inverse products. Bank
operators/inverses now use 40/60-digit evaluation of the original rational
pencil; exact bound frequency and stored pole-root coordinates are retained
for pole-centered bank targets, and their differences from the earlier
rounded coordinates are emitted. Both original double and refined inverse
identity residuals remain visible. Local subtraction residues are recomputed
at the same precisions from the actual nullspaces and derivative pairing.
The contour route remains a separate double-precision/refinement computation.

Heavy matrix metadata uses explicit row dimensions plus column dimension
offsets, matrix shape, computed axis-encoding residuals, and evaluated tensor
grade/homotopy support. This reconstructs each entry's restored dimensions
without repeating full leaf-path tables. Nonmatrix metadata retains the prior
format; the new inventory checks both encodings. Source grade support and
background origin remain separate from evaluated tensor grades.

## User-supplied stratum-coverage concern during regeneration

The S10 Lean overlap map prompted a read-only inspection documented in
[the coverage note](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_stratum_coverage_note.md).
The frozen regular-case run retains its scope. Each isolated native candidate
and each computed nullspace basis column is checked, but generic witnesses do
not cover exceptional parameter loci or every sheet region. The next stop must
address explicit locus/region coverage and the regular-summary guard limitation
before extending claims to threshold, coalescent, defective or bound-pole data.
No additional physical construction is inserted into this run.

## Completed regular-case checkpoint

The four focused checks and one fresh four-case run completed on the same
frozen source. The full run took 7265.559 seconds, peak child RSS 1727696 KiB,
exit 0 and empty stderr. Inventories account for 432 residues, 1728 contours,
384 bank inverses and 144 local pole subtractions, with no metadata gaps or
native nullity differences. Earlier spectral, sheet and inverse-Fourier
records are preserved within the documented carrier-order equivalence.
All five completed transcripts are published under scripts/out, with the old
annex payload preserved. The 189492142-byte main output exceeds the brief's
original size target; metadata/output compaction remains a mechanical debt.
The construction report and run record contain the actual residuals, sizes,
provenance and limitations. All ten broad TODOs and the absent export remain.
No exceptional-domain expansion, downstream work or new commit follows this
checkpoint; the user-supplied coverage audit remains part of its stop boundary.
