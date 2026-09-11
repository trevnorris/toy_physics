# S11c-d builder checkpoint: joint bulk-sheet paths and cut banks

2026-09-11. The planned bulk-radical continuation construction is implemented,
run across all four cases and published. The prior regular end-spectrum
checkpoint was first saved at user request as `cc6f8b3c`: six `.out` files through
DataLad/git-annex, the other files through Git. This subsequent build is an
unreviewed runnable checkpoint with no further commit. The complete S11c-d
engine and export remain unfinished.

[JointBulkSheetPath](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2023)
transports the actual reduced bulk radical along explicit frequency/momentum
paths, starting from the separately reduced real-axis branch operand. It
computes segment branch loci, two adaptive root transports, implicit-ODE
transport, cut encounters and dimensioned residuals.
[BulkContinuationAudit](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3153)
computes real-axis matrix joins, upper/lower rays, local path-order checks,
winding and branch-intersection controls, and refined frequency/momentum cut
banks. Both bank values are evaluated in the original rational physical matrix;
their matrices, jump and denominator operands are emitted with fingerprints,
digests, input frames, grades and restored L/T/M dimensions.
The [plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_plan.md)
and [construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_report.md)
detail the computation and its domain.

The fresh full run completed in `5779.514` seconds, peak child RSS `1725656` KiB,
exit 0 and empty stderr. Source hashes before/after match. The
[104,669,285-byte main transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
has 38,828 unique tags, one completion marker, no new metadata gaps or nonfinite
objects, and three empty dimensional-constraint records. It contains 24
continuation packets: 504 paths, of which 480 transport and 24 deliberate branch
intersections remain unresolved, plus 192 bank pairs. All 48 original native
PIT unresolved candidate labels remain explicit.

The physical real-axis matrix joins have 400 exact zero residual entries.
PIT joins evaluate both actual operands at 40/60 digits; their maximum residuals
are `2.870e-42` at `[-2,-2,1]` and `9.724e-63` at `[-4,-1,1]`. The maximum
ODE/root-transport endpoint difference is `6.397e-12` at `[0,-1,0]`; the maximum
sampled ODE radical residual is `3.386e-11` at `[0,-2,0]`. The root refinements
agree in these samples. Local path-order and double-loop return residuals are
at most `4.450e-16` and `4.441e-16`, both at `[0,-1,0]`. The four final focused
checks and three full inventories exit 0 with empty stderr.

The 24 native finite-root certificates and their 432 candidate records remain
intact. Among 8,016 native payloads including metadata, 24 raw differences are
only carrier-association/metadata ordering; their values, dimensions and grades
match. All numerical native records, 12 earlier full-sector symbols and 528
legacy candidate records are unchanged. All 96 inverse-Fourier carriers retain
zero inverse/source-image/remainder/branch residuals; projections of the 222
integral and eight row residual fingerprints remain zero. No upstream repair
was indicated by these checks.

The [run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_runs.json)
retains commands, source/input/cache provenance, resources and publication
hashes. The [joint inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_inventory.json),
[spectrum inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_spectrum_inventory.json)
and [Fourier inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_inverse_inventory.json)
retain the computed evidence. All five successful transcripts are published in
`scripts/out/` (106,367,981 bytes total); atomic replacement preserved the old
annex payload. They are working files for the next user-requested DataLad save.
The engine and both new instruments compile. Only `run` changed among existing
definitions; none was removed. The directive-named
`reduction/derived_or_declared.py` and `reduction/engine_output_checks.py` remain
absent and were not run. The solver/export contract below is unchanged.

All ten live TODOs remain. Explicit bulk-radical paths do not yet construct a
global physical sheet or the full profile-resolvent contour, including other
singularities and contour pinches. Finite-offset banks are not a completed
continuum measure. Double-precision path clearance and sampled ODE residuals
are diagnostics, not global certificates. Full spectra still require
exceptional threshold/denominator domains, mixed degeneracies and generalized
modes at defective roots. Negative-frequency branch checks do not extend the
positive-frequency channel-input API or provide scattering data.

Closed nonlocal current/flux normalization, complete two-ended scattering,
poles/Riesz/overlap, survival, flux bookkeeping, weak coefficients, Section 5
controls and own-row export remain later constructions.
`scripts/S11c_d_exports.py` remains absent because its required objects have not
been computed. Section 1 inputs remain SUPPLIED and unfalsifiable here; the
separate shear-normalization and c2 cross-engine operand/sign debts remain
open. No review leg, comparator, Wolfram engine, downstream stage or subsequent
commit ran. Stop at this build/run/report checkpoint.

## Retained user-approved solver/export contract


1. Preserve `EdgeReduction`, the positional three-parent fold, and the exact
   direct-lookup manifest. All numerical assembly must consume the computed
   reduced rows, including the full nonlocal terms and full coupling vertex.
2. Accept explicit, independent dimensionless profile functions w(xi), m(xi),
   their derivatives, asymptotic limits and tail information. A selected smooth
   step with an independently adjustable localized modulus bump is a numerical
   instance; it does not replace the interface class. Store the profile formula
   and digest in every case record. The current preflight input is recorded in `S11c_d_channel_preflight_input.json`.
3. Require a complete parameter map in a declared L/T/M unit frame, real
   continuum frequency and tangential momentum, small contrast, and
   sigma_W = eta_bg W_0/L_W for evaluations on the physical homotopy. Retain
   independent eta/sigma grades in the symbolic calculation. Test actual
   reference/end channel availability before attempting flux normalization.
   Do not manufacture an incident channel by assigning a sector label.
4. Compute both-end modes, left/right normalization, the S11b-derived current,
   and the variable-profile matching problem. Re-expand the continuum response
   to the retained rectangle; do not present a finite-contrast numerical
   solution as a higher-order continuum prediction. Retain evanescent matching
   modes, channel degeneracies and domain failures explicitly.
5. Numerical pole searches have an explicit profile, parameter map, sheet,
   bounded search region and isolating contours. Evaluate the retained operator
   without the continuum re-expansion. Record boundary/quadrature resolution,
   domain size, precision, root residuals, contour-count evidence, and changes
   under refinement. A bounded search does not establish a global pole set.
   An unsuccessful or inconclusive search is unresolved, not an empty pole set.
   Compute residues/projectors and sheet/decay/width/closure tests only for
   actually resolved candidates; emit spectral overlap, not capture probability.
6. Separate transparent symbolic expressions from evaluated numerical records.
   Symbolic operator/continuum/weak-coefficient exports remain differentiable
   SymPy expressions, compacted with algebraic equivalence checks. Numerical
   mode and pole datasets retain their input bindings, domain, convergence
   evidence, dimensions and truncated-model status. They are not stand-ins for
   a generic symbolic profile-dependent root function. This is the export
   distinction motivating the user-approved contract; the downstream consumer
   will need to bind the appropriate representation explicitly.
7. Fingerprints summarize already constructed/evaluated objects. They do not
   evaluate nonlocal integrals or replace a spectral solve. Algebraic PIT and
   physical numerical evaluation are separate records. All completed roots use
   fresh lowerCamel write-keys and the existing bind-closure/minimal-delta guards.
8. Finish the one-case path, then implement controls and bookkeeping, then run
   all four cases once and write the complete export. No review legs,
   comparator, Wolfram engine, downstream stage, or commit belongs to this lane.
