# S11c-d constant-end inverse symbols and normal-momentum residues

2026-09-11. The user-requested prior checkpoint was saved as `718e5ced` using
DataLad/git-annex for its five transcripts and Git for the other twelve files.
The new end-resolvent construction has four completed focused checks and one
fresh four-case run. All inventories and atomic publication completed. The
canonical main output now contains the new construction; its committed
predecessor is preserved in git-annex. This build remains unreviewed and has
no later commit. The complete S11c-d engine and export remain unfinished.

## Computed construction

[EndResolventAudit](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3363)
extends the existing `BulkContinuationAudit`. It consumes
the same computed reduced five-field end operator and intercepts the actual
bank matrices and transported endpoints emitted by the parent. The earlier
operator, spectrum and continuation computations are retained. It derives the
total normal-momentum derivative from the operator and implicit radical
relation, then constructs inverse matrices and local Laurent data.

For every native candidate, the [local Laurent construction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3489)
uses the actual operator's left/right nullspaces to give
the coefficients of a Laurent ansatz. Projecting its constant equation onto
the left nullspace supplies the derivative-pairing system for the residue.
The resulting derivative projector and both Laurent equation residuals are
emitted. Branch points, other normal-root coordinates and denominator-resultant
zeros supply numerical contour clearance. The [Cauchy construction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3441)
uses two radii and 32/64-node circles to
compute inverse-matrix Cauchy residues independently, together with derivative
integrals, determinant winding, second Laurent moments and refinement residuals.
These are fixed-real-frequency normal-momentum data; they are not the
profile-dependent frequency poles or bound-state Riesz data of section 3b.

Both original cut-bank operators are inverted. The inverse jump and the
separately evaluated inverse-jump identity operand are emitted with their
residual. At a bank target coinciding with a native normal-momentum pole on the
selected local lift, its computed Laurent term is explicitly subtracted.
The raw inverse, pole contribution, regular part and their finite-offset
refinement remain separate. Frequency-bank records do not substitute an empty
frequency-pole set or perform that normal-momentum subtraction.

The closest PIT bank samples exposed cancellation in double-precision inverse
products. The final bank construction evaluates the original rational pencil
at 40/60 digits and recomputes the relevant pole residues at those precisions.
For momentum banks it retains the exact bound frequency and restores matching
pole targets from the stored 50-digit normal roots; the differences from the
older rounded coordinates and operators are emitted. The radical is evaluated
on its actual relation and matched to the earlier transported bank. The
original double-precision inverses and identity residuals remain as operands.
The Cauchy route stays a separate double-precision/refinement calculation.

Heavy matrices have whole-object SHA digests and numeric fingerprints. Their
metadata restores each entry from explicit row dimensions and column dimension
offsets, with shape, computed axis-encoding residual and evaluated tensor
grade/homotopy support. Nonmatrix objects retain the earlier leaf metadata.
Unbound source-grade support, background origin, input binding and L/T/M frame
are recorded separately. Finite-input inversion is retained-operator data;
no continuum re-expansion or flux normalization is claimed.

## Focused evidence

All four final focused runs use the same engine hash, pinned producer symbol
caches and input. They exit 0 with empty stderr and unchanged before/after
source hashes. Each emits empty dimensional constraints.

| LAB_HELD / RHO4_CONSTANT check | Wall seconds | Peak child RSS, KiB | Output bytes |
|---|---:|---:|---:|
| Physical reference | 142.474 | 183240 | 3945664 |
| Physical left end | 100.854 | 179836 | 3130955 |
| Physical right end | 143.199 | 179996 | 3134720 |
| Reference PIT | 131.881 | 185088 | 4112614 |

The combined inventory accounts for 72 native candidates and residues, 288
contours, 36 bank pairs and 72 bank inverses. It finds no duplicate tags,
metadata gaps, unmatched objects or nonfinite objects. There are no nullity
differences from the native spectrum and no unresolved local-pole subtractions.
The PIT packet retains its four original unresolved sheet labels.

The largest modal/Cauchy residue difference is `5.607e-12` at `[2,2,-1]`;
the derivative-projector comparison is `1.928e-11` at `[0,0,0]`. The PIT
pole-count residual is at most `2.015e-12`. The original double inverse-jump
identity residual reaches `2.256` at `[3,2,-1]` near a pole; the 60-digit
calculation gives at most `4.632e-52` at that dimension. Its 40-to-60-digit
inverse refinement difference is at most `4.852e-25` at `[3,2,-1]`.
Left/right inverse residual maxima in the PIT packet are `3.168e-56` at
`[-1,1,0]` and `3.104e-56` at `[0,0,0]`. These are dimensioned finite-input
residuals, not global numerical error certificates.

The initial double-precision development checks are superseded. An early
40/60-digit pass exposed mpmath's rejection of negative slice indices; using
explicit nonnegative slice starts repaired that implementation error. Only
the four successful final-source focused transcripts were published.

## Full run, preservation and publication

The fresh four-case run exited 0 with empty stderr in `7265.559` seconds
(121.1 minutes), peak child RSS `1727696` KiB. All pinned source hashes agree
before/after the run and with the current files. The main transcript has
109,612 unique tags, one completion marker, three empty dimensional-constraint
records, and no new metadata gaps, unmatched objects or nonfinite records.

The [resolvent inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_inventory.json)
accounts for 24 packets, 432 native candidates and regular Laurent residues,
1,728 contours, 192 bank pairs and 384 bank inverses. All nullities agree with
the native records. There are 144 explicitly subtracted local normal-pole
terms and no unresolved local-pole subtractions at these bank targets. The
48 original unresolved sheet labels remain unchanged. The full-run maxima
for modal/Cauchy residues, derivative projectors, precision refinement and
inverse identities equal the focused maxima above. The maximum contour
coefficient condition is `1.206e6`, and the maximum coefficient-frame inverse
residual is `9.000e-11`; those are finite-input numerical diagnostics.

The [joint inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_joint_inventory.json)
retains all 504 paths (480 transported, 24 individually unresolved branch
intersections) and 192 original bank pairs. The native and continuation
comparisons each have 24 raw carrier-association/metadata ordering differences
and zero semantic differences. All 8,016 native payloads otherwise match,
including every numerical candidate record. The
[spectrum inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_spectrum_inventory.json)
retains 24 finite-root certificates, all 432 native candidates, 12 earlier
full-sector symbols and 24 legacy packets containing 528 candidate records.
The [Fourier inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_inverse_inventory.json)
retains all 96 carriers with zero inverse/source-image/remainder/branch
residuals. All 666 projections from 222 integral residuals and 132 projections
from eight row residuals are zero. All four inventories exit 0 with empty
stderr. No upstream operand mismatch was indicated.

The [main transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
is 189,492,142 bytes, SHA-256
`b04b50c87e227498a7acf48dec0a8a46e3e2288bb09045f479c9de2523a1d617`.
This exceeds the brief's original tens-of-MB target. New resolvent object
payloads occupy 24,576,504 bytes and their metadata 50,412,584 bytes; the
existing emission-index mechanism now occupies 14,406,893 bytes. Heavy
objects retain fingerprints and digests rather than full symbolic solutions.
Output/metadata compaction remains a mechanical size debt; no format rewrite
was inserted into this completed run.

The four focused transcripts are published beside the main output as
`S11c_d_end_resolvent_{reference,pit,left,right}.out`; all five total
203,816,095 bytes. Copy, fsync, digest verification and atomic replacement
preserved the old 104,669,285-byte annex payload. The new transcripts are
working files for the next user-requested DataLad/git-annex save. The
[run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_resolvent_runs.json)
retains commands, source/input/cache provenance, resources, inventory and
publication hashes, and the focused/development distinction.

The engine and both new instruments compile. Only `run` changed among
pre-existing definitions; `EndResolventAudit` was added and none was removed.
The directive-named `reduction/derived_or_declared.py` and
`reduction/engine_output_checks.py` remain absent and were not run. The
eight-point solver/export contract remains byte-for-byte unchanged.

## Remaining boundary

The user's S10 Lean concern prompted a [read-only stratum-coverage inspection](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_stratum_coverage_note.md)
while the source was frozen. Root checks visit every isolated candidate, and
nullspace residuals use every computed basis column. The 24 branch-intersection
controls have individual records. Sheet-region enumeration and deliberate
exceptional-parameter-locus evaluation remain absent. The regular-coverage
summary also does not explicitly gate on algebraic/geometric multiplicity
agreement, although all 432 committed native differences are zero. Those
limits must be addressed before extending the scope to exceptional or
defective modes; the inspection does not alter this regular-case run.

The full variable-profile resolvent and its contour pinches, normal Fourier
Green-function reconstruction and a limiting continuum measure remain open.
A local algebraic lift does not determine a global physical sheet. Numerical
contour clearance, winding and refinement diagnostics do not replace exact
complex-domain certification. Exceptional denominator/threshold and defective
mode domains remain explicit limitations, and the native sheet labels remain
unchanged.

Closed nonlocal energy current and flux normalization, complete scattering,
profile-frequency poles/Riesz/overlap, survival, bookkeeping, weak coefficients,
Section 5 controls and the own-row export remain unfinished. All ten broad
TODOs remain; no placeholder export or empty uncomputed pole set is supplied.
Section 1 inputs remain supplied and unfalsifiable here. The shear-normalization
and c2 cross-engine operand/sign debts remain open. No upstream repair was
indicated by the preservation checks. No review leg, comparator, Wolfram engine or
downstream stage runs in this lane. The separate Lean trial is untouched.
The other session's S10 work is also untouched. Stop at this build/run/report
checkpoint before extending the scope.
