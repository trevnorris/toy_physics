# S11c-d SymPy builder checkpoint

The preceding two-frequency current proof is committed as `5b007ddb`, with
source/reports in Git and its `.out` in DataLad/git-annex. The earlier
mechanical-load repair remains `c643112a`. This subsequent reference-subspace
continuation is uncommitted; no new upstream repair is indicated.

The native `ModalCurrentSubspaces` builder computes full left/right spaces,
nonlinear frequency pairings, and physical slab-plus-bulk current matrices for
all **18 isolated root/lift candidates** in LAB_HELD/RHO4_CONSTANT REFERENCE:
14 scalar and four two-dimensional spaces. All 18 frequency pairings have
full computed rank. Its current/energy reconstruction retains the actual
row-power bridge, projection defects, interface exchange, upper boundary and
radical derivatives. It builds on the preceding 802-residual acoustic/slab
proof without altering its constructor or the import/Fourier wiring.

Eight candidates have disk-certified bulk decay. Two physical-sheet,
exactly real-normal two-dimensional spaces admit signed physical right-current
normalization, giving four normalized basis directions. The maximum
current-reconstruction residual is 8.21e-13 and the signed-current normalization
residual is 3.24e-16 (rounded upward, double-precision reference-unit norms).
The row-dual pairing and physical field current remain separate typed objects;
their canonical adjoint-field export join remains ahead. No bare group velocity,
sector label or generic spectral witness substitutes for the calculated forms.

The [subspace report](S11c_d_modal_subspace_report.md) and
[inventory](S11c_d_modal_subspace_checkpoint.json) record scope, source pins,
tolerances and residuals. The new focused transcript is **1,134,351 bytes**,
with 1,330 object/metadata pairs, 2,658 indexed source assignments, 1,298 checked
tensors and 2,794 numerical residual scalars in 27 families. Its 250 quadratic
extraction scalar residuals are zero. SHA/PIT and literal payload checks,
metadata coverage and restored dimensions passed. The final run took 128.69
seconds and 145,324 KiB peak RSS, with exit zero and empty stderr. Correcting
factor unit labels reproduced every computed tensor exactly. Publication under
`scripts/out/` used an atomic replacement; the new transcript is uncommitted.

Next are source-level mass-transfer, bulk-load and boundary-work controls,
the canonical adjoint-field current join, and extension to both ends and all
cases; see the [continuation plan](S11c_d_modal_current_plan.md). The four-case
production transcript remains the committed repair checkpoint. No new full
run or incomplete export was produced under the retained one-case contract.
All ten broad TODOs remain, including completed-stage integration, scattering,
profile-frequency bound poles, survival, bookkeeping, controls and export.
Threshold, defective, sheet, denominator and convergence limitations and
supplied upstream debts remain explicit. Constant-end momentum roots remain
separate from section 3b's conditional profile-frequency bound poles. No
S10/Lean or authority edits, review legs, comparator, Wolfram, downstream run
or push occurred. Only the user-requested initial checkpoint was committed.

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
