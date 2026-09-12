# S11c-d SymPy builder checkpoint

The finite generalized-threshold-mode checkpoint is complete, verified and
published. It extends the previous finite-slice coverage work after Git checkpoint
`3332fe0a`; both extensions remain uncommitted. The full S11c-d program and
own-row export are still unfinished.

[NormalTaylorChains](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3440)
now computes exact full left/right chain spaces from the physical matrix germ.
[ThresholdModeAudit](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3562)
constructs polynomial normal modes, joins both radical lifts to the reduced
Fourier seed, derives local frequency unfolding, and tracks both full mode
spaces through approach points, bypasses and loops. Exact local exception-disk
certificates exclude additional enumerated exceptional frequencies in the
recorded neighborhood. The [construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_modes_report.md)
gives computation-line references, residuals and scope boundaries.

The four-case run completed in 10,051.411 seconds with 1,733,148 KiB peak child
RSS, exit zero, empty stderr and 27 unchanged source/input pins. Three focused
LAB_HELD/RHO4_CONSTANT runs use the same final source. The
[main transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
is 83,919,481 bytes, atomically published with three focused threshold outputs.
The old annex data is preserved. [Run provenance](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_modes_runs.json)
records hashes, resources, verification and development repairs.

There are 24 finite threshold points, 48 complete left/right chain spaces and
96 individual chains. Each space has kernel counts and basis ranks `[2,4,4]`,
with two chains of length 2. All 7,320 exact scalar chain/block/polynomial-mode/
radical/multiplicity residuals are zero. All 24 local disk certificates and
connections are defined. The 192 approach points and 2,448 path nodes have full
nullity 2 under the stated 53-bit SVD tolerance; all 144 normal and 336 bulk paths
are defined. Residuals remain separated by restored dimensions in the
[full census](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_full_summary.json).
Twelve zero-frequency intersections remain explicitly unresolved.

Inventories found 229,636 unique tags with no new metadata gaps or nonfinite
objects. All 124,804 previous tags are accounted for, with no unclassified
physical changes. The prior 432 modes, 432 residues, 1,728 contours, 504 joint-sheet
paths and their unresolved domains are preserved. Fourier reconstruction and
solved dimensional residuals remain zero. Lossless decoding restores all payloads
and 229,634 source-index assignments. The named generic triage helpers were not
available; the dedicated inventories actually run are recorded in the report.

All ten broad TODOs remain. Next is the reduced nonlocal S11b current and mode/
flux normalization, followed by two-ended scattering, profile-frequency poles,
survival, bookkeeping, weak coefficients, physics controls and export. Remaining
singular, mixed/defective and global sheet domains stay explicit. These local
normal-threshold results do not establish a global sheet atlas or a section 3b
profile-frequency bound pole. Section 1 and the c2 operand/sign and
shear-normalization debts remain supplied premises. No new premise, upstream
repair or change of work breakdown was required. No S10/Lean edit, review,
comparator, Wolfram, downstream run, export, commit or push was performed.
The retained solver/export contract below is byte-identical.

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
