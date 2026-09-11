# S11c-d builder checkpoint: constant-end resolvents and normal-momentum residues

2026-09-11. The user-requested prior checkpoint was committed as `718e5ced`:
five `.out` files through DataLad/git-annex and twelve other files through Git.
The subsequent constant-end resolvent construction is implemented, run across
all four cases, inventoried and published. It remains an unreviewed runnable
checkpoint with no later commit. The complete S11c-d engine and export remain
unfinished.

[EndResolventAudit](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3363)
extends the existing continuation builder. It derives the total normal-momentum
derivative of the computed physical end operator, constructs full-subspace
Laurent residues, and independently integrates inverse matrices on two radii
with 32/64 contour nodes. It evaluates both cut-bank inverses at 40/60 digits
and explicitly subtracts computed local normal-momentum poles where applicable.
Original double-precision operands, coordinate refinements, sheet labels,
residuals, grades and restored L/T/M dimensions remain visible. These are
fixed-frequency end data; profile-frequency bound poles and their Riesz data
remain separate future constructions. The
[plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_plan.md)
and [construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_report.md)
record the computation and its limits.

The single fresh full run completed in `7265.559` seconds, peak child RSS
`1727696` KiB, exit 0 and empty stderr. Source hashes before/after match the
current files. The transcript contains 24 packets, 432 regular Laurent
residues, 1,728 contours and 384 bank inverses. It records 144 local pole
subtractions, no nullity differences and no unresolved subtraction at these
bank targets. The maximum modal/Cauchy residue difference is `5.607e-12` at
`[2,2,-1]`, and the derivative-projector difference is `1.928e-11` at `[0,0,0]`.
Near-pole cancellation produces a raw double inverse-jump residual of `2.256`
at `[3,2,-1]`; the refined 60-digit calculation gives at most `4.632e-52` at
that dimension. Those are finite-input residuals, not global error bounds.

All four final focused checks and four full inventories exit 0 with empty
stderr. The main output has 109,612 unique tags, one completion marker, empty
dimensional constraints and no new metadata gaps or nonfinite objects.
All 24 native root certificates and 432 candidate records are preserved.
Native and continuation comparisons each retain 24 carrier-association/metadata
ordering differences with no semantic difference. The 504 earlier paths,
192 bank pairs, 48 unresolved sheet labels, 12 legacy full-sector symbols and
528 legacy candidate records remain intact. All 96 inverse-Fourier carriers
retain zero inverse/source-image/remainder/branch residuals; projections of
222 integral and eight row residual fingerprints remain zero. No upstream
operand mismatch was indicated.

The [run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_runs.json)
contains the source/input/cache pins, resources, inventories and publication
hashes. The [main transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
is 189,492,142 bytes; the four focused outputs bring the total to 203,816,095
bytes. This exceeds the original tens-of-MB target: the new numerical objects
occupy 24.6 MB, their metadata 50.4 MB, and the emission index is now 14.4 MB.
Heavy objects still use fingerprints/digests; metadata/output compaction remains
a mechanical size debt. All five files are in `scripts/out/`; atomic publication
preserved the old annex payload. They await the next user-requested DataLad save.
The engine and new instruments compile. Only `run` changed among existing
definitions, one class was added and none removed. The directive-named
`reduction/derived_or_declared.py` and `reduction/engine_output_checks.py` remain
absent and were not run. The solver/export contract below is unchanged.

The user's S10 Lean concern is addressed in the
[stratum-coverage note](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_stratum_coverage_note.md).
Native checks visit every isolated root and both normal lifts; nullspace
residuals use all computed basis columns under recorded numerical tolerances.
Each of the 24 unresolved branch-intersection tests has its own path record.
These checks do not enumerate all joint sheet regions or exceptional parameter
strata. A representative witness also requires constant rank/property on its
stratum or explicit handling of exceptional subloci. The regular-coverage
summary does not explicitly gate on algebraic/geometric multiplicity agreement,
although all committed native candidate differences are zero. This guard and
targeted exceptional-locus coverage need attention before extending the domain.
No checks were rerun solely for this read-only inspection.

All ten broad TODOs remain. The full profile resolvent, Fourier Green-function
reconstruction, global sheet/continuum measure, exceptional/defective modes,
closed nonlocal current/flux normalization, scattering, profile-frequency
poles/Riesz/overlap, survival, bookkeeping, weak coefficients, Section 5 controls
and own-row export remain unfinished. `scripts/S11c_d_exports.py` remains absent
because its required objects have not been computed. Section 1 inputs remain
SUPPLIED and unfalsifiable here; shear-normalization and c2 cross-engine
operand/sign debts remain open. The other session's S10 and Lean work is
untouched. No review leg, comparator, Wolfram engine, downstream stage or
subsequent commit ran. Stop at this build/run/report checkpoint.

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
