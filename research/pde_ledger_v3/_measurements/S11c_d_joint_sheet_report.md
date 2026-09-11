# S11c-d joint bulk-sheet paths and cut-bank construction

2026-09-11. The planned bulk-radical path and cut-bank construction completed
after preservation checkpoint `cc6f8b3c`. All four focused checks, the fresh
four-case regeneration, inventories and atomic publication completed. This is
an unreviewed runnable checkpoint; the complete S11c-d engine and export remain
unfinished. This checkpoint was subsequently saved as `718e5ced` using
DataLad/git-annex for its outputs. The later
[end-resolvent checkpoint](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_report.md)
regenerated the canonical main output. This report retains the historical
joint-sheet evidence; its original main transcript is linked below from the
preserved annex payload.

## Computation

[JointBulkSheetPath](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2023)
consumes the actual reduced quadratic radical relation and a reduced real-axis
branch operand. Substitution of the affine path ansatz gives each segment's
radicand polynomial. Its zeros supply dimensionless path-parameter clearance
and branch-intersection diagnostics. Adaptive root transport uses two step
fractions; a separately integrated implicit differential equation supplies an
endpoint comparison and sampled equation residual. The path record contains
the vertices, seed, cut encounters, resolution and ODE tolerances, endpoints,
argument change and residuals. A branch intersection remains unresolved.

[BulkContinuationAudit](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3153)
receives the computed reduced branch bindings and physical end matrix. It
computes the branch-to-radical scale from the quadratic coefficient, compares
the separately reduced branch operands, and checks the real-axis physical
matrix against its algebraic representation for both frequency signs in the
propagating and evanescent domains. It then constructs upper/lower frequency
rays, local frequency/momentum path-order checks, winding loops, branch-locus
intersections and separate frequency/momentum cut banks.

The two bank sides use decreasing real offsets `10^-3`, `10^-5`, `10^-7` times
the computed branch scale. The engine evaluates the original rational physical
matrix on each transported bank value and emits both matrix fingerprints,
their computed jump, denominator coefficients and radical refinement
differences. Momentum-bank targets include the existing unresolved native
candidates when present; their old labels are retained. Heavy matrices have
whole-object SHA digests and numerical tensor fingerprints. Every object has
grade and restored L/T/M metadata in its declared numerical frame.

## Focused evidence

The physical reference source/algebraic comparison has four zero 5-by-5 matrix
residuals and zero radical residuals. The final focused path run completed in
`58.902` seconds with peak RSS `169792` KiB, exit 0 and empty stderr. It contains
39 paths: 36 transported and three branch-intersection records. Its metadata
inventory has no gaps, unmatched objects or nonfinite objects; dimensional
constraints are empty.

The largest physical-reference ODE/root-transport endpoint difference is
`6.397e-12` at `[0,-1,0]`; the largest sampled ODE radical-equation residual is
`3.386e-11` at `[0,-2,0]`. The root-transport refinements agree literally in
these samples. The local path-order residual is `4.450e-16` at `[0,-1,0]`.
One frequency loop changes the root sign; a second loop starting at that
transported root returns within `4.441e-16` at `[0,-1,0]`. Each loop records a
downward-cut encounter and leaves the supplied fixed-real-momentum frequency
cut chart while its radical lift remains defined. This measured path dependence
does not supply a new sheet-selection rule.

The final PIT reference run completed in `27.849` seconds, peak RSS `169640`
KiB. The two physical end checks took `16.721` and `48.079` seconds, with peak
RSS `170032` and `170040` KiB. All exit 0 with empty stderr, unchanged source
hashes and empty dimensional constraints. Their combined inventory contains
102 paths: 96 transported and six deliberate branch intersections. There are
36 bank pairs across the three finite offsets, and no metadata gaps or
nonfinite objects. The PIT packets preserve four unresolved candidate labels.

The PIT real-axis matrices are evaluated independently at 40 and 60 decimal
digits. Their maximum scaled 60-digit comparison residual is `9.724e-63` in
the declared coefficient frame; the maximum 40-digit residual is `2.870e-42`
at `[-2,-2,1]`. This numerical route avoids large exact algebraic-number-field
simplification. The physical-input matrix comparisons remain exact.

The [run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_runs.json)
retains commands, source/cache provenance and measured child resources. The
[plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_plan.md)
defines the construction scope and stop condition.

## Four-case evidence and preservation

The full run exited 0 with empty stderr in `5779.514` seconds (96.3 minutes),
peak child RSS `1725656` KiB. Source hashes before and after match the published
engine and inputs. The
[104,669,285-byte main transcript](/var/projects/toy_physics/.git/annex/objects/Mx/pQ/MD5E-s104669285--9a91e625b06d0613d85cd5947b81fd25.out/MD5E-s104669285--9a91e625b06d0613d85cd5947b81fd25.out)
has SHA-256 `e51274711cf0fa7ea59758c3d11eadbd35dc02639b1b2f72c234f208bf4e1030`.
It has 38,828 unique tags, one completion marker, no duplicate tags, no new
metadata gaps or nonfinite objects, and three empty dimensional-constraint
records. The prior main transcript was 95,506,147 bytes; heavy-object
fingerprints keep this extension's increase to 9,163,138 bytes.

All 24 physical-input/PIT reference/left/right packets completed. Their 504
paths comprise 480 transported paths and 24 deliberate branch intersections
left unresolved. There are 192 bank pairs and 16 recorded fixed-real-momentum
frequency-cut encounters. All 48 original native PIT candidate labels remain
unresolved; the additional bank records do not replace them.

The 16 physical real-axis matrix joins contain 400 exact zero residual entries.
The 16 PIT test points emit 32 matrices at 40/60 digits, comprising 800 residual
entries. Their maxima are `2.870e-42` at `[-2,-2,1]` and `9.724e-63` at
`[-4,-1,1]`, respectively. The largest 60-digit scaled residual is also
`9.724e-63`, against the declared `1e-30` comparison tolerance. Reduced-seed and
real-axis radical residuals are zero. Across all transported paths, the maximum
ODE/root endpoint difference is `6.397e-12` at `[0,-1,0]`, and the maximum
sampled ODE radical residual is `3.386e-11` at `[0,-2,0]`. Root step-refinement
differences are zero in these samples. The maximum local path-order and
double-loop return residuals are `4.450e-16` and `4.441e-16`, both at `[0,-1,0]`.

The earlier 24 degree-11 root certificates still account for all 432 native
`(k,q)` candidates, with zero degree-count residuals and isolated, disjoint
root disks. Of 8,016 native spectrum payloads including metadata, 24 differ
textually only in ordering: twelve physical-input carrier associations and
their metadata. Key/value, metadata path, dimension and grade comparisons are
identical; the inventory retains both the raw differences and this restricted
ordering comparison. All numerical native spectrum evidence is unchanged.
All 12 earlier full-sector symbols and 24 legacy mode packets (528 candidate
records) match the committed baseline literally.

The unchanged inverse-Fourier instrument finds all 96 carrier inverses, 288
inverse residual entries, 96 source-image residuals, 288 remainder entries and
24 branch residuals intact and zero. All 666 projections of the 222 integral
residual fingerprints and 132 projections of the eight row residual
fingerprints are zero; the slot census and metadata coverage remain complete.

The [joint inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_inventory.json),
[focused inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_focused_inventory.json),
[spectrum inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_spectrum_inventory.json)
and [Fourier inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_joint_sheet_inverse_inventory.json)
retain these computations. The three full inventories exited 0 with empty
stderr in `31.785`, `152.228` and `6.853` seconds; peak child RSS was `185288`,
`206284` and `184972` KiB. The run record retains commands, source/input/cache
hashes, runtime versions and publication hashes.

All five successful transcripts are in `scripts/out/` (106,367,981 bytes total).
Atomic replacement left the previous annex payload unchanged. They were subsequently
saved through DataLad/git-annex at `718e5ced`. The engine
and both new instruments compile. Only `run` changed among pre-existing engine
definitions; no definition was removed. The directive-named
`reduction/derived_or_declared.py` and `reduction/engine_output_checks.py` remain
absent and were not run. The eight-point solver/export contract is preserved
byte-for-byte in the builder report.

## Limits and next dependency

This computation supplies explicit bulk-radical lifts and finite-offset bank
data. It does not assign a path-independent physical sheet in joint complex
frequency/momentum space. S11b's frequency-cut chart is applied only at fixed
real momentum. General joint paths retain their history, including winding;
no decay requirement is re-imposed at complex frequency. Negative-frequency
branch tests do not extend the positive-frequency channel-input API or create
negative-frequency scattering data.

The clearance calculation uses double-precision roots in a dimensionless path
parameter, not an exact complex-domain certificate. ODE residuals are sampled
checks, not global error bounds. Bank refinement is not a certified limiting
continuum spectral measure. The full nonlocal profile-resolvent contour,
other singularities and contour pinches remain uncomputed. Accordingly both
the generic-sheet and full-spectrum TODOs remain open, as do exceptional
denominator/threshold domains, mixed degeneracies and defective-root modes.

Current/flux normalization, complete scattering, poles/Riesz/overlap, survival,
bookkeeping, weak coefficients, Section 5 controls and own-row export remain
later constructions. Section 1 supplied inputs and the separate shear and c2
cross-engine operand/sign debts remain open. No upstream repair has been
required by these checks. No review leg, comparator, Wolfram engine or
downstream stage is part of this build.
