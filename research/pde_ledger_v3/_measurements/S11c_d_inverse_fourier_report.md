# S11c-d carrier inverse-Fourier construction

2026-09-10. The focused checks, fresh four-case run and output inventories
completed. The regenerated transcript and focused evidence are published under
`scripts/out/`. This extends checkpoint `aff093da` and conveys no review clearance.

## Computed construction

[FourierCarrierReconstruction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1006)
computes a Gaussian weak Fourier kernel, its mass, characteristic-function
limit and second moment. The delta ansatz weight and contour factor are
computed from those mathematical operands. Every imported profile and jet
definition is instantiated at the actual transfer or middle-leg argument.
No imported Fourier normalization is supplied as a number.

The [source inverse](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1116)
uses the actual three-dimensional definition, profile ansatz, coordinate and
momentum Jacobians, and computed weak kernels. The separate
[image inverse](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1147)
consumes the existing `hat()` result. It integrates localized/jet terms and
computes residues of the regulated Abel part in a complex transfer variable.
It retains the regulator until inversion, then records the positive half-line,
negative half-line and symmetric subtraction-origin values. Profile endpoint
limits temporarily become explicit coefficient symbols for CAS operations;
their definitions are restored in the reconstructed values.

[Source-image linearity](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1105)
separates the computed Fourier coefficient from the normal measure Jacobian.
The raw and resulting operands are emitted. The source route of
[EdgeReconstruction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1284)
now obtains hat factors from the source definitions and independently eliminates
the tangential constraints. It previously shared `hat()` with the forward route.
Computed zero source integrals are evaluated explicitly. The forward `hat()`,
action contraction, profile/field ansatz, import fold and spectral constructions
are unchanged.

These are weak inverses on the supplied smooth, tangentially homogeneous
interface class, with the additional half-line-tail premise for transformed
zero jets. They are not inverses on arbitrary three-dimensional backgrounds.
The check retains all five source slots and all actual carrier arguments.
Every new physical record has grade and restored L/T/M metadata; heavy row
comparisons retain their carrier fingerprints and SHA digests.

## Focused evidence

On `LAB_HELD / RHO4_CONSTANT`, all 24 carrier inverse records have literal
zero residuals at the three coordinate domains: 72 entries with dimension
`[0,0,0]`. All 24 source-image differences and all 72 inverse-remainder census
entries are zero. The 28 action-integral and two row residual fingerprints
have zero projections. Dimensional constraints and unresolved dimensions are
empty; no new metadata leaf is marked `ZERO_MAP`. The check took `142.683`
seconds and peak RSS `127272` KiB, with exit 0 and empty stderr.

The independent mutation run took `23.267` seconds and peak RSS `122200` KiB,
with exit 0, empty stderr and empty dimensional records. Doubling a source
coefficient exposes the live profile or normal-jet difference at all three
domains. Reversing either source or image phase exposes the two half-line
differences. Removing the Abel terms exposes the right endpoint, left endpoint
and their mean at the origin. The two tangential jets remain computed zeros.
These are mathematical transform controls, not Section 5 physical profile-FORM
or scattering/flux controls.

## Four-case evidence and publication

The full run took `3706.142` seconds, peak RSS `1723024` KiB, exit 0 and empty
stderr. Source hashes stayed fixed throughout. Its
[81,914,994-byte transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
has 25,306 unique tags. All 96 carriers have literal zero inverse residuals
(288 entries at dimension `[0,0,0]`), zero source-image differences and empty
inverse remainders. All 24 branch reconstruction residuals are zero. The 222
action-integral and eight row residual fingerprints have zero projections.
The independently read source census matches the emitted carriers, branches,
integrals and all five slots for both rows in all four cases. There are no
coverage or metadata gaps, and all three dimensional records are empty.

All 12 full-sector pencil symbol payloads and all 528 root/nullity/sheet
records match checkpoint `aff093da` literally. The 528 regular rectangular
jets remain defined (336 scalar and 192 two-dimensional nullspaces), and the
64 unresolved sheet labels remain explicit. These are the inherited PIT
records, not a completed generic end-spectrum or physical channel census.
The new checks exposed no additional upstream repair.

Commands, source snapshots, cache provenance, resource measurements and hashes
are in the [run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_runs.json).
The [carrier inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_inventory.json)
and [spectral inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_inverse_fourier_spectrum_inventory.json)
retain the measured counts and comparison scope. The successful
[focused production transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_inverse_fourier_production.out)
and [mutation transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_inverse_fourier_mutations.out)
are also at canonical output paths. Publication used same-directory temporary
files and atomic replacement; the previous annex payload's hash is unchanged.

The engine and both new instruments compile. The two directive-named checking
utilities remain absent. No review leg, comparator, Wolfram engine, downstream
stage or commit ran during the build. The user subsequently requested a local
preservation checkpoint, with all three output payloads in DataLad/git-annex
and the source, instruments, reports and inventories in ordinary Git. This
save conveys no review clearance.

`ALL_CARRIER_INVERSE_FOURIER_ROUNDTRIPS` is removed from the
live TODO. Ten constructions remain: full end spectra, generic sheets, closed
nonlocal current/flux normalization, complete scattering, poles/Riesz/overlap,
survival, flux bookkeeping, weak coefficients, Section 5 controls and own-row
export. The export remains absent because its required roots are uncomputed.
