# Fixed-input outgoing prescription: saved-output validation complete

The fixed-input, separated-point outgoing-prescription candidate is built and
validated. The explicitly approved corrected validator completed at
01:38:06 UTC on 2026-09-28 (19:38:06 MDT on September 27). The
[current checkpoint](S11c_d_outgoing_prescription_checkpoint.json) records
`VALIDATED_FIXED_INPUT_SEPARATED_POINT_OUTGOING_PRESCRIPTION_CANDIDATE`.
This completes the bounded prescription stage: the saved full 5x5 symmetric
principal-value Fourier integrals and signed-current delta terms pass their
saved-output checks on the unchanged physical slice, for distinct normal
positions. The integrals remain unevaluated.

This does not accept a complete outgoing Green operator, full FORM,
two-asymptote response, A11/A12, a coincident-position distributional extension,
complex-frequency retarded equivalence or a radiating witness. The fixed
physical slice remains evanescent.

## Validated evidence

All 30 structural/source/unit/domain checks, 26 denominator certificates,
50 decaying-tail checks, both complete pole-block residual/mutation bundles
and seven PV/delta assembly checks passed. The validator used only the new
saved outputs. It did not reconstruct derivatives or tails, import a source
producer, solve modes/roots/LU, or evaluate integrals.

Both physical pole blocks have nullity two: block 17 has current direction
-1 and block 16 direction +1. Each maximum residual norm is
`5.744387398727955e-16`, below the `1e-8` bound. Every saved residual, projected
operand, norm and mutation roundtrip difference is zero. Both direction
mutations have norm `2.8284271247461903`; reversing the delta coefficient
moves it by norm `3.7717975547957288`. Exact cofactor residue identities and
complete block shapes passed. Full operands and returns are retained.

The checked joins include the source inverse/determinant/adjugate, saved
regular density and phase/Fourier measure, source units and actual reference
grades, and excluded opposite-sheet lifts. The denominator checks join actual
saved polynomial coefficients and their signs; the tail checks join saved
series and leading terms for all 25 entries at both ends of the momentum
axis. The domain explicitly requires
`Ne(s11cdNormalPosition, s11cdSourceNormalPosition)`. The zero polynomial
contact part follows the decaying tails; it is not a diagonal distributional
extension. The contour-sign identity remains a wiring check, not independent
physics evidence.

## Preservation and completion checks

The [production checkpoint](S11c_d_outgoing_prescription_unlimited_checkpoint.json)
preserves 1069 files / 6,786,185 bytes and 226 completed construction operations
(217 restored, nine new). Block 17's completed non-journal bundle was restored
byte-for-byte. Its historical validation-pending status is preserved; the
current checkpoint above supersedes that status without rewriting history.

The [corrected validation checkpoint](S11c_d_outgoing_prescription_validation_v2_checkpoint.json)
preserves 1854 files / 3,490,903 bytes and 441 completed validation operations:
217 restored returns and 224 new operations. All 217 restored return files
are byte-identical to their originals. Seventeen newly serialized argument
tuples differ in bytes; the checkpoint records both routes and does not claim
input byte identity. Every original input/return remains unchanged, and the
ordered journal bypassed the completed functions.

All 1069 candidate files, 901 original validator files and 498 constructor
input routes still match their pinned hashes. Every completed validation
operation's saved input, return and receipts was checked. Worker stdout is
byte-identical to `checks.json`; scientific and infrastructure stderr are
empty. Coordinator, supervisor, guard, native containment and effective
limits agree. All job PIDs and the service cgroup are gone. The hook was
armed before launch and finished with its completion notification queued;
any delayed duplicate event requires no new launch.

The [first validator failure and correction](S11c_d_outgoing_prescription_validation_disposition.md)
remain preserved. That validator rejected SymPy `BooleanTrue` using Python
identity `is True`, after the first tail's mathematical checks had passed.
The user explicitly approved the corrected attempt. It restored the old
false-flag return, corrected only the literal-Boolean comparison in one new
operation, confirmed the other six checks unchanged, then finished the
remaining validation. No automatic retry or scientific reconstruction ran.

## Measured cost and exact artifacts

| Completed run | Worker time | Guard time | Peak whole-job memory |
|---|---:|---:|---:|
| No-deadline saved-return constructor | 3394.435 s (56 m 34 s) | 3395.462 s | 132,718,592 bytes (126.57 MiB) |
| Approved corrected saved-output validator | 3.930 s | 4.420 s | 89,706,496 bytes (85.55 MiB) |

These are the final two runs, not a total including earlier stopped attempts.
The remaining block-16 adjugate derivative accounted for 3344.348 constructor
seconds. The constructor had 1698 resource samples; the validator had four.
Both had zero swap/memory events and at most three tasks. Actual resource
limits remained 2 GiB, zero swap, one CPU, nice 15, 32 tasks and one native
thread. The completed constructor alone used the approved indefinite-duration
exception. Validation used the unchanged shared guard with 900 seconds outer
and 840 seconds native.

The primary files below remain under
`_scratch/s11c/s11c-d-outgoing-prescription-20260927/unlimited/complete`:

| File | SHA256 |
|---|---|
| `outgoing-prescription-candidate.pickle` | `7702b16291588d2ea823f1503b63c5d1009f111d896980eb209b9164a55037a8` |
| `prescription-domain.pickle` | `2b296cda6d33fdb0089fb625cd8a7d045cb217dd4363b4b61f83f4858ad0de8e` |
| `spatial-sign-identity.pickle` | `d2a186f3e545fc6519d3b6f3f089704f851077a36f7d407116b5c7c850918882` |
| `checks.json` | `534ea04d84106a39bac46404a5ea93c4c98cfb066e915e2fbb625bd7129c5cbf` |

The current checkpoint links exact production and validation artifact hashes.
Source/guard/manifest/gate/authorization snapshots are in sibling
`unlimited-complete-source` and `validation-v2-complete-source`. Large machine
checkpoints have explicit `.gitattributes` diff exclusions. Binary scratch
artifacts remain at their original paths. All earlier failed/stopped runs and
the original failed reference validator remain intact. Lean, shared guard and
the protected builder suffix are unchanged.

## Next method boundary

The earlier literal reviews remain Claude **NEEDS REVISION** and Grok
**CLEAR FOR THIS BOUNDED STAGE**, with the bounded local disposition preserved.
Saved-output validation is not a new external method review or fresh
independent CLEAR.

The next dependency is a concrete two-asymptote response plan. It must retain
the actual unequal left/right asymptotic fields and grades, specify how
nondecaying forcing and half-line/distributional tails are handled, and keep
the full nonlocal operator and cross-interface reconstruction/cancellation
checks. Simply integrating nondecaying forcing against this reference
prescription does not complete that method. The
[constructor plan](../directives/S11c_d_FORM_constructor_plan.md) and
[execution plan](S11c_d_scattering_form_execution_plan.md) retain those
obligations and the separate FORM/export/A11/A12 requirements.
No new physics job, replay or external submission is authorized by this
completed validation or its checkpoint.
