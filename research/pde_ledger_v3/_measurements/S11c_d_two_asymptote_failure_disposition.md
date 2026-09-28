# First two-asymptote attempt: preserved coefficient-guard failure

The authorized attempt stopped at `block17-regular-0-0`, its first new local
regular inverse coefficient, after **24.807482 worker seconds**. It raised
`ValueError: denominator leading coefficient nonzero undecided`. No end lift
or forcing result is accepted. The failure is preserved in
`_scratch/s11c/s11c-d-two-asymptote-20260927/production`; no retry has run.
The [completion checkpoint](S11c_d_two_asymptote_failure_checkpoint.json)
records the actual logs, resource samples, operation receipts, artifact hashes
and all source posthashes.

The worker selects the first structurally nonzero denominator-series entry
and then requires its SymPy `is_zero` property to be Python `False`. That
guard failed. Its selected order, coefficient expressions and zero-property
values were local variables inside the incomplete operation and were not
saved before the guard. The exception alone does **not** establish whether
the selected coefficient is zero, nonzero, or undecidable. In particular,
this inspection does not assume the property was `None`, or authorize
weakening the guard.

The full incomplete input is preserved at
`production/complete/operations/0061-block17-regular-0-0/input.pickle`
(9,710 bytes, SHA-256
`5e6d6a95b6e3336d3ab976ab82d29b742ca68aa3b0b65100f64126ae863fc1b4`).
It has a started receipt but no completed return. All **61 completed
operations**, their complete input/return receipts, and all **279 files /
3,575,556 bytes** in the partial output are retained. All **101 source/input
hashes** and the twelve launch snapshots match their pinned bytes.

Seventeen saved 5-by-5 source/end-pencil/grade/branch-derivative summaries
report exact zero residuals. These are worker-reported results, not a new
independent saved-object validation. The native census records 25 cells,
71 local terms and 160 nonlocal terms. The completed runtime grade-inventory
join and its operands are retained. The premise summary still says
`UNVERIFIED_NONLOCAL_PLANE_JET_EXTENSION`: no regulator/profile occurrence at
zero grade was detected, but neither translation invariance nor moving
momentum derivatives through ordered integrals and the Abel limit was
established. End coefficients, end equations, carrier reconstruction,
responsive mutations and final domain assembly were **not reached**.

This was not a resource-limit stop. The shared guard enforced 900 seconds
outside the existing supervisor and 840 native seconds, 2 GiB memory, zero
swap, one CPU, nice 15, 32 tasks and one native thread. The whole guard took
25.384841 seconds; peak memory was **124,932,096 bytes** (119.145 MiB), with
zero swap and zero memory-event counts. Fourteen resource samples report at
most three tasks. Coordinator, guard, child and supervisor all record exit
1. Scientific stderr contains the 1,114-byte traceback; stdout is empty and
`checks.json` does not exist, so stdout/checks identity is inapplicable.
There is no scientific acceptance based on process status.

The original reviewed packet and raw reports remain unchanged: Claude
**NEEDS REVISION**, Grok **CLEAR FOR THIS BOUNDED INGREDIENT**. The earlier
local review disposition is preserved; no fresh independent CLEAR is
claimed and no review submission is proposed here.

## Prepared next diagnostic, awaiting separate approval

The [diagnostic worker](S11c_d_two_asymptote_entry_diagnostic.py) consumes
only the saved incomplete entry tuple. It uses the same scalar simplifier,
convolution and polynomial-series helpers as the failed worker. It computes
only the missing branch/numerator/denominator series inside that unfinished
operation, persists each complete intermediate before examining zero flags,
records the actual flag values/types and the original guard outcome, then
stops. It performs no new nonzero proof, quotient, regular coefficient,
end-field continuation, numerical probe or integral evaluation. It does not
replay any of the 61 complete operations. This is an instrumented calculation
of the incomplete entry, not merely metadata inspection or validation of an
already saved denominator series.

The [pending gate](S11c_d_two_asymptote_entry_diagnostic_pending_gate.json)
and [preparation receipt](S11c_d_two_asymptote_entry_diagnostic_preparation.json)
pin the diagnostic, all partial results and original inputs, and the unchanged
shared guard/supervisor. A separate new `entry-diagnostic` directory and the
same 900/840-second, 2-GiB/zero-swap/one-CPU/nice15/32-task/one-thread limits
are prepared. The hook targets session
`01a0e01b-ef84-7192-817f-584cda5d339b` and must arm before the job. Missing
explicit diagnostic approval prevents gate activation and launch. The
completion event expressly forbids retry/replay, so this calculation needs
new approval; no constructor correction or continuation is authorized by it.

No scientific pickle was restored during this completion inspection or
preparation. Accepted reference/outgoing artifacts, previous failures,
incident history, the protected builder suffix and shared guard are intact;
Lean work is untouched. Full Green/FORM, two-asymptote response, A11/A12,
forcing-domain pairing, diagonal extension, complex-frequency retarded
equivalence and a radiating witness remain open. The fixed slice is
evanescent.
