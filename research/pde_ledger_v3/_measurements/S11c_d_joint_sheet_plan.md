# S11c-d joint bulk-sheet continuation plan

2026-09-11. Start from the user-requested preservation checkpoint `cc6f8b3c`.
The six prior transcripts were saved through DataLad/git-annex and the other
files through Git. This next construction follows build, run, report, stop;
its resulting changes are not included in that earlier commit.

Authorities are the program brief, the shared physics sections 2/3a/3b and the
build directive, with S11b section 1b supplying the frequency-continuation
prescription. Preserve the eight-point solver/export contract in the builder
report. No review leg, comparator, Wolfram engine or downstream stage runs.

## Objective

Extend the fixed-positive-frequency bulk-radical chart with explicit joint
frequency/normal-momentum paths. Compute continuation from the actual reduced
real-axis branch operands, including negative-frequency propagating seeds and
evanescent seeds. Frequency rays follow S11b's upper-rim prescription. General
joint paths carry their path history; a path-dependent lift is not a global
physical-sheet label. Keep sheet transport separate from decay, mode direction,
full resolvent singularities and flux normalization.

## Execution

1. Use the source-pinned LAB_HELD/RHO4_CONSTANT reference cache from the completed
   end-spectrum run. Extract the radical relation and real-axis branch operands
   from the computed reduced symbols/branch reduction record. Compute the
   positive/negative-frequency and propagating/evanescent joins to the algebraic
   pencil before extending its domain. A discrepancy with an earlier operand
   triggers the stop condition below.
2. Construct piecewise-linear path transport on the actual algebraic radical
   curve. Derive segment branch loci and the implicit differential equation
   from that relation. Record path clearance, finite resolution, refinements,
   radical residuals and independent differential-equation transport residuals.
   A branch intersection remains unresolved. Never reselect by decay at a
   complex frequency.
3. Compute both signs of frequency, upper/lower frequency rays, a local joint
   path-order comparison, branch-point intersections, and winding loops. Record
   encounters with S11b's downward frequency cuts where that fixed-real-momentum
   chart applies. Compute separate limiting cut-bank data under offset
   refinement, including the fixed-frequency candidates left unresolved on a
   branch ray. Do not replace those original labels by choosing a bank silently.
   Evaluate the actual rational physical end matrix on the two bank values;
   report its jump and denominators as operands. Bank samples alone do not
   construct the full continuum spectral measure or profile resolvent.
4. After focused checks, integrate the computation into all reference/left/right
   end records for each of the four cases. Run the full engine once, inventory
   the new records, and check preservation of the previous root census and
   inverse-Fourier evidence. Publish completed transcripts atomically under
   `scripts/out/` without writing through annex links. Record source hashes,
   commands, actual child RSS, wall time, dimensions/grades and limitations.

Heavy matrices receive computed fingerprints and SHA digests; path operands
and scalar residuals remain inspectable. Numerical coordinates are coefficients
in the declared L/T/M frame, with their restored dimensions recorded. Algebraic
PIT and the explicit physical input remain separate.

## Completion boundary

This step targets the bulk radical's path and cut-bank construction. The full
generic-sheet TODO also covers contour transport for the complete nonlocal
profile resolvent, including other singularities and possible pinches. Retain
that TODO unless the computed scope actually supports removing it. Full end
spectra additionally retain exceptional denominator/threshold domains, mixed
degeneracies and generalized modes at defective roots. Current/scattering,
pole/Riesz, bookkeeping, controls and export remain later constructions.

If the work exposes an incorrect earlier operand, requires a new physical
premise, or requires changing the agreed approach or work breakdown, preserve
the evidence and stop to explain it before making that repair. The existing
supplied-input, shear-normalization and c2 cross-engine debts remain open.

## Completed bounded checkpoint

2026-09-11. Steps 1–4 completed with four final focused checks and one fresh
four-case run. The run and all inventories exit 0 with empty stderr. There are
24 continuation packets, 504 paths (480 transported, 24 deliberate unresolved
branch intersections) and 192 bank pairs. Earlier numerical spectrum and
inverse-Fourier evidence is preserved. Five completed transcripts were
published atomically under `scripts/out/`; the previous annex payload is intact.
The [construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_report.md)
and run record carry the actual counts, residuals, resources and hashes.

No stop-condition discrepancy arose. The full profile-resolvent contour,
continuum measure and exceptional/defective spectrum remain unresolved, so
both relevant broad TODOs stay open. All ten program TODOs and the absent
own-row export are recorded in the builder report. No review, downstream run
or subsequent commit was performed.
