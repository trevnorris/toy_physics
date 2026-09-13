# S11c-d both-end current source — runtime repaired, RIGHT joins nonzero

The all-material-symbolic RIGHT source construction for LAB_HELD/RHO4_CONSTANT
and its exact arithmetic checks have completed on this box. The runtime refactor preserves
the algebra: **127 cross-product residuals are exactly zero**. The source
comparison itself has **10 nonzero scalar entries**: all five entries of each
mass and mechanical face-join row, at `(epsilon,eta,sigma)=(0,1,0)` and lambda
order one. These are retained terms, not omitted higher-order remainders.
This establishes a source-consistency discrepancy; its cause has not been
assigned to the downstream reconstruction or the inherited closed rows.

Of the 33 top-level residual scalars in 21 source families, the remaining
23 are exactly zero. Both independently computed closure-row increments,
both reconstructed face rows, their literal differences, all original source
operands and denominator exclusions remain in the diagnostic packet.
The reference current mapping and both-end frequency packets are preserved.

## Computation and representation

`ClosedAcousticEnergy.cancel_background_expression` (native line 2082)
collects field amplitudes and background powers structurally, preserving
background-dependent denominator factors and all material symbols. Only then
does it cancel each material coefficient (line 2077). The source projections
remain at their original sites. The row comparisons are computed at native
lines 2228–2229. No action, profile, material parameter or response kernel was
changed, and no extra background truncation was introduced.

The [resumable instrument](S11c_d_end_current_resumable.py) saves source
operands and completed stages atomically. A resume requires matching source,
input and operand hashes. A changed emitter may seed construction packets
only under identical native/loader/input pins, and may seed exact proofs only
under an unchanged identity-function AST. Seed paths and packet hashes are
recorded explicitly. This does not relabel an interrupted computation as a
completed residual.

Exact proof arithmetic combines fractions before direct sparse polynomial
cross-multiplication. An uncombined mechanical numerator/denominator had
artificial degrees 154/146; exact fraction combination reduced them to 21/13
in 0.66 seconds. Original before/after expressions and their uncombined
exclusions remain separate from the combined proof operands. All 127 literal
cross-product residuals computed zero without numerical material binding.

Proof metadata derives zero units from the physical rows and denominator
factors. Untruncated rational expressions carry their exact numerator and
denominator grade/lambda supports plus the full fraction's restored unit;
they are not mislabeled as finite polynomial series. Heavy operands use the
native carrier/PIT fingerprint and SHA representation.

The final source/proof replay and emission completed in **345.70 seconds**,
peaking at **201436 KiB**, and wrote **9520298 bytes**. This timing excludes
previous construction and proof attempts; it is not a fresh uncached build
benchmark. The run exited **1 after completing emission**, through the guard
on the computed retained source discrepancies. Publication validation completed: **85 source objects, 1651 cancellation
objects, 368 native metadata paths and 3484 tags**. The
[RIGHT transcript](../scripts/out/S11c_d_end_current_source_right.out) is
atomically published; its SHA256 is
`7d0987a2a552da82ddaf2dd864b1881bc710b6bb90aabf892162a51d22f6149d`.
The [RIGHT inventory](S11c_d_end_current_source_right_checkpoint.json) pins
its sources and full validation census.

## Preserved LEFT baseline and historical runtime evidence

The earlier LEFT construction produced 33 top-level and four nested face
exact-zero residuals. Its 85 source objects remain preserved in
[the original LEFT output](../scripts/out/S11c_d_end_current_source_left.out)
and [inventory](S11c_d_end_current_source_left_checkpoint.json). That run took
236.99 seconds and peaked at 200532 KiB. The fresh runtime regression completed in **123.03 seconds**, at **201636 KiB**,
with all 33 top-level and four nested face source residuals and all **82** exact
arithmetic identities zero. All **85** source objects agree with the baseline:
83 are structurally identical; the complete 52-entry parameter map agrees by
symbolic key; two face-load expression differences have zero exact residual.
The initial positional map comparison was corrected without changing either
source packet. The new [LEFT runtime transcript](../scripts/out/S11c_d_end_current_source_left_runtime.out)
has **2625518 bytes**, SHA256
`d60af16865aef8437dc629f8849ba3550a0d5e5484432852c7b57670c9819234`.
Its [inventory](S11c_d_end_current_source_left_runtime_checkpoint.json) records
1066 cancellation objects, 368 native metadata paths and 2314 tags. The original
LEFT output was not overwritten.

The original RIGHT stop was a computation bottleneck, not a residual result.
Its [immutable inventory](S11c_d_end_current_source_stop.json) is unchanged.
The inventory's mechanical label for frozen line 2175 was one line off:
that line is the mass increment; the mechanical increment is line 2176.
Later attempts also located dense cancellation in the face load.

The 8 MiB stack attempt segfaulted in the final diagnostic GCD after saving
125 operations. The 256 MiB retry completed all 127 operations and saved the
acoustic objects, then segfaulted during eager expansion of an equivalence
check. Stack exhaustion was not established as a common cause. A sparse
verification attempt saved 96 zero identities before being stopped to combine
its repeated denominator factors. The next attempt computed all 127 zero
identities but stopped during rational-grade metadata emission; those results
were preserved and the metadata representation repaired. No OOM was observed.
Frozen sources, raw operands, progress/resource logs and failure traces remain
under `/tmp/s11cd-right-source-resume-20260913`; the
[runtime plan](S11c_d_right_current_runtime_plan.md) records the sequence.

## Next boundary

Publication and the LEFT regression are complete. Stop before changing a
physical source formula. Trace the two first-contrast row differences back through the
chemical driver, mass-rate normalization and face reconstruction to the
actual reduced closed rows. Establish whether repair belongs in S11c-d or an
earlier source stage before editing either. A verified RIGHT physical current
and flux normalization retain this dependency. Independent end-pencil and
frequency data remain valid within their documented domains.

Global parameter/sheet/exceptional coverage, complete scattering, section 3b
profile-dependent frequency poles, controls/bookkeeping and export remain
open. No S10/Lean or authority edit, full four-case rerun, review/comparator/
Wolfram/downstream run, commit or push belongs to this checkpoint.
