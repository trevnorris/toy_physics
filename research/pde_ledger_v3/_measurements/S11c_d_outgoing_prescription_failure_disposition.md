# Outgoing prescription: preserved stop and proposed continuation

The single guarded run **stopped without producing an outgoing kernel**.
The worker took 14.543 seconds; the supervisor took 14.784 seconds. Peak
whole-job memory was 115,191,808 bytes (about 110 MiB), with zero swap and zero
memory events. The 2-GiB / zero-swap / one-CPU / nice-15 / 32-task / one-thread
limits were verified. This was an implementation check failure, not a resource
failure. No automatic retry occurred.

The complete [failure checkpoint](S11c_d_outgoing_prescription_failure_checkpoint.json)
pins all 181 saved files / 2,264,672 bytes, all 43 completed operation returns,
the actual exception and posthashes of all 48 consumed inputs. All source
hashes were unchanged. Scientific stdout is empty, stderr contains the
798-byte traceback, and there is no checks JSON, tail result, residue block
or kernel artifact to accept. Earlier accepted reference-kernel outputs and
the original validation failure remain unchanged.

## Actual failure

The review-requested source-denominator check was implemented too narrowly:
it accepted only a nonzero constant times a power of the positive radical.
The actual saved determinant denominator reduces on the acoustic branch to

```text
D(r) = (-28512000000000 + 5760000000000*i)*r^2
     + (576000000000 + 5760000000000*i)*r
     + 288000000000,                 r > 0.
```

This is the literal expression in `production/complete/determinant-denominator.json`.
Its imaginary part is `5760000000000*r*(r+1)`, strictly positive for `r>0`.
Thus the recorded expression is not evidence of a real-axis denominator zero;
it is outside the monomial-only sufficient condition. This algebraic
observation motivates the correction; no CAS or scientific pickle restoration
was run during failure inspection.

## Concrete bounded continuation

The separate [continuation worker](S11c_d_outgoing_prescription_continue.py)
keeps the failed worker unchanged. Its new nonzero certificate accepts a
polynomial denominator when its real or imaginary component has coefficients
of one weak sign and at least one coefficient of that strict sign. On a
strictly positive radical, every monomial is positive, so such a component
cannot vanish. The worker saves the actual polynomial, coefficients, component
and sign before the guard. It stops if this sufficient proof is unavailable;
there is no root search, numerical tolerance or assumption that a denominator
is harmless. The displayed coefficients are not hard-coded into the worker.

The continuation restores the **43 complete saved operation returns in order**
and bypasses those operation functions. Their original argument/return hashes
are pinned. The first new operation is the corrected nonzero proof using the
saved denominator reduction; subsequent work is the unfinished adjugate
denominator checks, tails, mode-block residues, coverage joins and integral
assembly. Cheap reconstruction of local bindings and checks of restored
identities still occur; source-producing operations and the 43 completed
journaled calls are not rerun. Every reused return is labelled and linked to
its original operation in the new journal.

The original `production` directory, reviewed packet, failed worker, manifest,
gate and all operands/results remain intact. The new worker, manifest, proposed
gate, launcher and completion message target a fresh `continuation` directory.
The proposed gate is **PENDING_EXPLICIT_CONTINUATION_APPROVAL** and cannot launch.
Preparation uses only standard-library metadata/hash/AST operations.

## Approval boundary

The authorized first attempt has been consumed. The user's explicit no-retry
instruction therefore requires a new approval for this corrected saved-return
continuation. Proposed limits remain one 900-second / 2-GiB / zero-swap /
one-CPU / nice-15 / 32-task / one-thread guarded worker around the existing
supervisor, with the completion hook armed first. No automatic second attempt
or external review submission is included. The original literal review verdicts
remain Claude NEEDS REVISION and Grok CLEAR FOR THIS BOUNDED STAGE; no fresh
independent CLEAR is claimed for this elementary certificate correction.

Until that approval and actual successful validation, the outgoing prescription,
two-asymptote response, full Green operator, FORM and A11/A12 remain unfinished.
