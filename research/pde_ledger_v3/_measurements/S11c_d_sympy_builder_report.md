# S11c-d builder checkpoint: regular algebraic end spectra

2026-09-11. Regular finite algebraic spectrum coverage at the bound inputs is
implemented and verified across all four cases. The focused checks, full run,
inventories and canonical publication completed. This is an unreviewed runnable
checkpoint extending `c5af3181`; the full S11c-d engine and export remain unfinished.

[EndSpectrumCoverage](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2747)
consumes the computed physical five-field operator and canonical sector pencil.
It constructs their coordinate pullback, clears rational denominators, computes
the radical polynomial and exceptional-locus gcds, isolates its finite roots,
and evaluates left/right nullspaces in the original physical matrix. The output
retains both normal-momentum lifts, multiplicities, sheet and classifier domain
records, equation residuals and oblique nullspace projector residuals. Heavy
objects have compact fingerprints and digests with input, grade and dimension
metadata. These are finite-input algebraic modes, not flux-normalized channels.
The [plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_plan.md)
and [construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_report.md)
detail the computation and remaining domains.

The full run took `3963.933` seconds, peak RSS `1724396` KiB, with exit 0 and
empty stderr. The [95,506,147-byte transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
has 33,436 unique tags with no metadata gaps and three empty dimensional
constraint records. All 24 native polynomial certificates account for their
roots: 432 `(k,q)` candidate records across reference/left/right, physical-input
and algebraic-PIT packets. Each degree-11 polynomial has nine distinct radical
roots and 18 candidates. All isolation and separation checks hold, with zero
degree-count and pullback residuals. The largest root-refinement difference is
`5.796e-50` in inverse-time reference units. Dimensioned field/projector residual
maxima are retained in the construction report and inventory.

The 12 physical-input packets have no unresolved sheet labels; the 12 native
PIT packets retain 48. These classify the existing fixed-frequency chart and
do not count open channels. All 12 previous full-sector symbols and all 528
legacy candidate records are unchanged, including their earlier unresolved
labels. All 96 carrier inverse-Fourier checks remain intact, with 288 zero
inverse residual entries and zero projections of the 222 integral and eight
row residual fingerprints. No additional upstream repair was needed.

The [run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_runs.json)
retains commands, source/cache linkage, hashes, publication and resource records;
the [spectrum inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_inventory.json)
and [Fourier preservation inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_inverse_inventory.json)
retain the computed checks. All six successful transcripts are in `scripts/out/`
(98,218,875 bytes total). Atomic publication preserved the previous annex payload;
the user subsequently requested a preservation checkpoint. Its six outputs
are saved through DataLad/git-annex; source, instruments, reports and inventories
use ordinary Git. This local save conveys no review clearance. The engine and two
new instruments compile. Only `run` changed among pre-existing engine definitions.
The directive-named `reduction/derived_or_declared.py` and
`reduction/engine_output_checks.py` remain absent. No review leg, comparator,
Wolfram engine, downstream stage or commit ran in this build.

All ten live TODO items remain. Full spectra still require generic sheet
continuation, cut-bank/continuum treatment, exceptional threshold/denominator
loci, mixed degeneracies and generalized modes at defective roots. The other
constructions are closed nonlocal current/flux normalization; complete two-ended
scattering; poles/Riesz/overlap; survival; flux bookkeeping; weak coefficients;
Section 5 controls; and own-row export. The spectrum constructor is a regular
algebraic checkpoint within that program. `scripts/S11c_d_exports.py` remains
absent because its required constructions are unfinished.

Section 1 inputs remain SUPPLIED and unfalsifiable here. The separate
shear-normalization and c2 cross-engine operand/sign debts remain open. The
[previous inverse-Fourier report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_report.md)
retains its historical checkpoint evidence. Transform and coordinate/input
checks do not supply the Section 5 physical profile-FORM controls.

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
