# S11c-d builder checkpoint: carrier inverse Fourier reconstruction

2026-09-10. The all-carrier inverse-Fourier construction is complete across all
four cases. The full run, output inventories and canonical publication completed.
This is an unreviewed runnable checkpoint; the complete S11c-d engine and export
remain unfinished.

[FourierCarrierReconstruction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1006)
computes independent source and reduced-image inverses for every profile and
jet carrier, including both middle-leg arguments. It derives weak kernels and
normalizations by integration, retains the regulated zero-jet tail through
inversion, and reconstructs both half-lines and the subtraction origin. The
source route for the two full rows now uses the actual imported Fourier
bindings. The forward Fourier reduction, import fold and spectral construction
are preserved. The [construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_report.md)
identifies the computation sites, scope and focused mutation evidence.

The full run took `3706.142` seconds, peak RSS `1723024` KiB, with exit 0 and
empty stderr. The [81,914,994-byte transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
contains 25,306 unique tags. All 96 carrier inverses have zero residuals (288
dimensioned entries); the 96 source-image and 24 branch residuals are also zero.
All 222 integral and eight row residual fingerprints have zero projections.
The five-slot census, metadata coverage and all three dimensional checks pass.
Deliberate source-weight, source/image-phase and Abel-removal mutations produce
the corresponding live residuals; vanishing tangential jets remain zero.
These are mathematical transform controls, not the Section 5 physical controls.

All 12 full-sector pencil symbols and all 528 root/nullity/sheet records match
checkpoint `aff093da` literally. The 528 regular jets and 64 unresolved sheet
labels are preserved. No additional upstream repair was needed for this step.
Exact commands, source/cache linkage, hashes and resource records are in the
[run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_runs.json),
with the [carrier inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_inventory.json)
and [spectral inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_spectrum_inventory.json).
All successful transcripts are under `scripts/out/`; atomic publication
preserved the previous annex payload. The engine and both new instruments
compile. The directive-named `reduction/derived_or_declared.py` and
`reduction/engine_output_checks.py` remain absent. No review leg, comparator,
Wolfram engine, downstream stage or commit ran during the build.

The user subsequently requested a local preservation checkpoint. The three
output payloads (83,389,349 bytes total) are saved through DataLad/git-annex;
the source, instruments, reports and inventories use ordinary Git. This save
does not confer review clearance or complete the engine or export.

`ALL_CARRIER_INVERSE_FOURIER_ROUNDTRIPS` is removed from the live TODO. Ten
constructions remain: full end spectra; generic sheet continuation; closed
nonlocal current and flux normalization; complete two-ended scattering;
poles/Riesz/overlap; survival; flux bookkeeping; weak coefficients; Section 5
controls; and own-row export. The next spectrum work must go beyond the current
regular reference/end mode jets and their PIT domains. `scripts/S11c_d_exports.py`
remains absent because its required roots are uncomputed.

Section 1 inputs remain SUPPLIED and unfalsifiable here. The separate
shear-normalization and c2 cross-engine operand/sign debts remain open. The
[prior mode-jet report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_rectangular_mode_jet_report.md)
and [fixed-frequency repair report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_repair_report.md)
retain their historical evidence.

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
