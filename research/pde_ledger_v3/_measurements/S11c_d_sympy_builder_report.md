# S11c-d SymPy builder checkpoint

The bulk exceptional-geometry and generalized-threshold work was committed as
`18f3236a`: 38 ordinary Git files and eight `.out` files through DataLad/git-annex.
The verified full threshold transcript remains unchanged. Its construction and
coverage boundaries are recorded in the
[threshold report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_modes_report.md).

The next [current construction plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_nonlocal_current_plan.md)
now has a runnable one-case preflight, but normalization is paused under the
user's instruction to stop at a required change of work. The reconstructed
S11b mass closure matches the reduced mass row exactly. The five mechanical-load
coefficients are the negatives of the independently reconstructed S11b face-work
load. The conservative thickness-stiffness coefficient anchors the same row
orientation in both constructions. This touches the carried face-force/closure-
fold sign debt; no sign was selected, inherited row changed, or upstream repair
launched. A bounded sign audit is proposed before validated current normalization.

[SlabEnergyBalance](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1870)
computes the actual slab energy balance, its material density-rate boundary
correction and chemical functional derivative.
[ClosedAcousticEnergy](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2049)
derives acoustic energy/current, depth integrals and both face closures. The
[focused report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_nonlocal_current_report.md)
records 25 zero scalar residuals, five nonzero mechanical comparison coefficients,
and the zero mechanical sum and row-orientation residual. It states the precise
scope: REFERENCE, LAB_HELD/RHO4_CONSTANT, with symbolic closure parameters.

The focused run took 195.725 seconds and 200,248 KiB peak child RSS. Its output
is 1,738,791 bytes with 146 unique tags, no metadata gaps or nonfinite objects,
and eleven unchanged source pins. Exit zero records completed emission, not
agreement of the mechanical comparison. All 38 previous top-level definitions
are unchanged. The new classes are exercised only by the focused instrument;
the main four-case output remains the committed threshold producer. Compilation,
whitespace, artifact and contract-retention checks passed.

All ten broad TODOs remain. Combined nonlocal current, left/right and flux
normalization, two-ended scattering, profile-frequency poles, survival,
bookkeeping, weak coefficients, controls and own-row export are unfinished.
The previous finite-root, threshold and selected-sheet coverage limitations
remain explicit. Section 1 premises, c2 operand/sign debt and the separate
shear-normalization debt remain supplied. The new current work is uncommitted.
No S10/Lean edit, review, comparator, Wolfram, downstream run, new full four-case
regeneration, export or push occurred. The retained contract below is byte-identical.

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
