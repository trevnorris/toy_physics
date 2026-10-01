# Selected uniform near-unity result — 2026-10-01

The exact modal/acoustic match is regular on the tested transverse subspace at
both ends. The surrounding finite-depth sweep is **UNRESOLVED**, so this is a
partial result, not completion of the near-unity check. No further run was
launched. Evidence and the metadata/byte audit are in
`S11c_d_near_unity_uniform_completion.json`; the immutable run remains under
`_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01`.

This is the strict rest-bulk, LAB_HELD/RHO4_CONSTANT effective-coefficient family
at omega 3 and tangential momenta (1/5, 1/10). Equality means the selected modal
speed equals the effective bulk speed. It is not primitive R7 calibration, a
draining background, or a defect-loss calculation.

## What completed

| End | Exact matching bulk speed | Normal momentum | Selected current magnitudes |
| --- | --- | --- | --- |
| LEFT | sqrt(3/2) | +/-sqrt(595)/10 | 0.2195335965, 32.9300394777 |
| RIGHT | sqrt(150/101) | +/-sqrt(601)/10 | 0.2206377121, 33.4266133829 |

Both normal directions and both transverse polarizations passed at each end's
own match. The original five-row pencil residuals passed; the selected lift has
rank two. Current eigenvalues have the group-direction sign and remain nonzero;
the normalized current matrices are +/-identity to floating roundoff. Full
limiting pencil rank three/nullity two is a numerical diagnostic, not a complete
outgoing-basis certificate.

At these four end/direction checks, all eight named face projections on both
faces and both harmonic legs are exactly zero in the readable saved outputs.
The selected bulk-normal current, depth current and interface-power matrices
are zero. Native pencil omission moves the residual by about 0.0414; a contributing
current-entry omission moves it by about 0.0041. Direction reversal disagrees with
the independent direction, and all four native eW-velocity row controls move by
1. These are source/control checks, not a defect leakage measurement.

All four exact grazing-limit operations completed. Radiating and evanescent
approaches give identical finite selected limits, with zero saved differences
for the pencil, currents and faces. Each LEFT sign saved 322/295 unique scalar
limit evaluations on the two paths; each RIGHT sign saved 338/311, including raw
and reduced order operands. These are source-internal consistency results, not
independent validation or proof of neighborhood-wide smoothness.

## What stopped the sweep

Of 24 scheduled end evaluations, two exact matching evaluations passed and 22
remained unresolved: four successful sign checks and 44 unresolved sign checks.
Every unresolved call stopped in `domain_at_point`, before its current/face
point checks. No nongrazing native sheet-reversal probe was reached.

The predicate is `not finite(bound) or bound.is_zero is not False`. It requires
an explicit symbolic nonzero property; an undecided property is refused. The
message “nongrazing source denominator outside admitted domain” must not be
read as proof of a pole. For example, a saved rejected value is

`2250 + sqrt(266)*(-1365 + 900*I)/4 + 7500*I`.

Its imaginary part is `7500 + 225*sqrt(266)`, strictly positive, so this particular
denominator is finite and nonzero. The worker did not separately save the
`is_zero` and finiteness flags. This source/readable-output inspection does not
certify every rejected factor or promote the other points to passed status.
The next concrete dependency, if separately authorized, is exact nonzero
certification of the refused denominators using saved operands, followed only
by the unfinished point/control checks. No tolerance relaxation, blanket
acceptance of unknown flags, completed-work replay or defect sweep is implied.

## Cost, preservation and acceptance boundary

Worker time was 515.239 seconds (8.6 minutes); peak cgroup memory was 260,624,384
bytes (248.55 MiB). Swap and all memory-event counters stayed zero. Actual
enforcement verified 8 GiB native/cgroup caps in the 16 GiB pool, CPU 15, one
native thread, 32 tasks, and a 4 GiB host reserve. There was no deadline and no
restart. The exit code was 2 because the scientific status was unresolved,
not because a time or memory limit terminated the run.

All 41 input/source posthashes, 32 source snapshots, nine original input copies,
5,383 opaque journal blobs and every receipt reference passed byte/hash checks.
All 71 started operations have terminal receipts: 27 complete and 44 unresolved.
Strict scientific stderr is empty and stdout is byte-identical to checks.json.
The result tree has 32 files / 46,689,296 bytes. Inspection restored no scientific
payload. Runtime scratch stays ignored and uncommitted; no new tracked generated
file exceeds 1 MiB. Shared guard, scientific sources, Lean and the protected
builder were not edited.

Both literal independent method verdicts remain **CLEAR FOR THIS SELECTED
UNIFORM METHOD**, for packet
`bf55b921706cc1d2c3df225fe87ec8a20405df782de771f6e50f8d60506f3752`.
They do not clear the worker or its results. Actual constant-end source/profile
joins completed, but unsupplied nonuniform and direct mixed-grade composition
remain conditional dependencies. The separate upstream assessment is not
presumed delivered or clear. No loss, calibration, drain response, complete mode
census, Green/FORM/A11/A12 or defect-sweep acceptance is claimed.
