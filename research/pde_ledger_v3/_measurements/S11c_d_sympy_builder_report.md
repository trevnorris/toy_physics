# S11c-d builder checkpoint: paused exceptional-domain construction

2026-09-11. Work resumed from `56595cf7`, preserving the S11c-d `f22cb682`
producer. The source is runnable and the frozen reference check completed.
Work then stopped under the user's instruction to pause when the breakdown
needs changing. The right-end check was interrupted. **No fresh four-case run
was started: the main transcript still represents `f22cb682`, not this source.**

The reason is a separate bulk-continuum threshold. The new end-mode resultants
find exceptional points of the discrete end roots. The completed reference
operands also give a collision of bulk normal branch points at frequency
coefficient `sqrt(5)` (`2.2360679775` in the declared time unit), where none of
the fifteen current end-mode exclusion conditions vanishes and the end-spectrum
polynomial has no radical root at zero. Bulk branch/denominator strata need
an independent family of records before intersecting them with end-mode loci.
This is a coverage finding; it does not establish an upstream operator error.
The [scope probe](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_bulk_threshold_scope_probe.json)
and [paused plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_exceptional_strata_plan.md)
record the next required change.

The [regularity criteria](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2967)
now include algebraic/geometric multiplicity agreement, ranks of both full
bases, and the normal-derivative pairing rank. A development reference run
retained all 18 candidates; each new criterion was true. Software controls
removing required record fields prevent regular coverage from being reported.
Existing root isolation, import wiring, Fourier reduction, sheet paths and
constant-end resolvent computations are retained.

[EndExceptionalSlice](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3193)
computes square-free multiplicities, end-mode exclusion polynomials, exact
real-locus gcds/intervals, and generic rank identities on a declared frequency
slice with other carriers bound. Determinant order and identities for all
minors of the required size bound nullity from both sides on the recorded
regular domain. This is a frequency-slice certificate, not a parameter-variety
or global-sheet certificate. The reference factors have degrees seven and two
with multiplicities one and two; 1+25 minor identities have zero remainders
and reconstruction residuals.

The [targeted threshold computation](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3342)
uses the original rational physical matrix and exact full bases. At the
computed end-mode threshold (`0.27386127875` in the declared time unit), both
radical roots give matrix rank three, both nullities/basis ranks two, one
coincident normal lift and normal-derivative pairing rank zero. Equation
residuals are emitted. Generalized threshold modes, exceptional-point physical
sheet membership and the zero-frequency intersection remain unresolved.

The frozen reference run used engine SHA
`ed359f16e78deccffdbce8f15cbb7e9735d60344e8fe0b38d5dd5e8ea457419e`,
completed in 363.169 seconds at 265,436 KiB peak child RSS, exit zero, empty
stderr and unchanged source hashes. Its
[reference transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_exceptional_strata_reference.out)
is 254,459 bytes with 179 unique tags, empty dimensional constraints and no
metadata gaps, unmatched objects or nonfinite payloads. Source-index decoding
restores all 176 preceding tag assignments. The
[run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_exceptional_strata_runs.json)
separates this completion from the deliberately interrupted right-end run.

The [output codec](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_output_codec.py:19)
shares exact payload text and records source-line positions compactly. Every
tag, dimension/grade payload, residual and digest survives decoding. The old
189,492,142-byte transcript encodes to 52,531,413 bytes with zero payload or
source-line assignment differences. Updated inventory readers accept both
representations; the codec's command-line expansion writes the previous
representation to a separate new file. The main annex payload is unchanged.
The eight-point solver/export contract below is byte-identical. Compilation
and whitespace checks completed. The directive-named reduction triage tools
remain absent and were not run.

All ten broad TODOs remain, including full sheet/exceptional coverage, closed
nonlocal current and flux normalization, scattering, profile-frequency
poles/Riesz/overlap, survival, bookkeeping, weak coefficients, Section 5
controls and the own-row export. `scripts/S11c_d_exports.py` remains absent.
These are constant-end results, not section 3b profile-dependent bound poles.
Section 1 and the carried c2 operand/sign debts remain supplied premises.
S10/Lean files were not edited here; new changes in those areas appeared during
this work and were left untouched. No review leg, comparator, Wolfram engine,
downstream stage or commit ran. The remaining focused checks, full regeneration
and main-output publication await the coverage-plan revision.

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
