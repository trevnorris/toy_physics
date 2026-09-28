# Incomplete-entry diagnostic completed

The user explicitly approved the prepared diagnostic. The
[authorization](S11c_d_two_asymptote_entry_diagnostic_authorization.json) and
[activated gate](S11c_d_two_asymptote_entry_diagnostic_gate.json) permitted
only the missing intermediate calculation of `block17-regular-0-0`, with no
quotient, end-lift continuation or replay of the parent's 61 complete
operations. The hook was armed before science and queued completion to
session `01a0e01b-ef84-7192-817f-584cda5d339b`.

The diagnostic completed in **4.771488 worker seconds**. Its
[checkpoint](S11c_d_two_asymptote_entry_diagnostic_checkpoint.json) verifies
actual logs, containment, exact stdout/checks byte identity, empty scientific
stderr, every operation input/return receipt and all source hashes. The
diagnostic evidence is available; no end-lift or forcing result is accepted.

## Actual reason for the previous guard failure

The denominator constant coefficient is the exact saved `Integer(0)`.
The selected first-order coefficient has `is_zero = None`, with actual type
`builtins.NoneType`. The original nonzero guard therefore still fails. Its
full saved expression is:

```text
sqrt(36445)*(-27911338644360 - 79738064747200*I)/81
+ sqrt(555)*(-259257681200800 + 79429384856000*I)/81
```

These are the actual diagnostic JSON literals, not substituted examples.
By inspection of that expression its real part is

\[
\operatorname{Re}d_1 =
-\frac{27911338644360\sqrt{36445}
       +259257681200800\sqrt{555}}{81}<0.
\]

Both radicands and integer multipliers are positive, so the displayed
coefficient is nonzero. This identifies a concrete exact sign-certificate
route. The diagnostic did **not** implement a new nonzero certificate,
weaken the guard, calculate the quotient, or return a regular coefficient.
There was no independent scientific pickle roundtrip during completion
inspection. The negative-real-part argument here is an inspection of the
saved expression; it is not a claim that the constructor has been corrected
or validated.

The zero flags for the next two saved denominator coefficients are also
`None`; none was used to bypass the leading-coefficient guard. All four
numerator and denominator series coefficients are saved. The denominator
return completed before flag observation began.

## Preserved evidence and costs

The run root is
`_scratch/s11c/s11c-d-two-asymptote-20260927/entry-diagnostic`.
The checkpoint pins all **59 output files / 627,266 bytes**, comprising
**12 completed journal operations** with complete operands and returns.
The key routes below are relative to its `complete/` directory.

| Artifact | SHA-256 |
|---|---|
| `operations/0009-numerator-series/value.pickle` | `12e36dad89709ae2e51198e99c7b926ca2094291b990350ced357bfda3127670` |
| `operations/0010-denominator-series/value.pickle` | `2374af9c4c4a9e80b1b84b4b1d66c6d8625e4d56d895f83eeae027e6c43ffaf8` |
| `operations/0011-observe-original-denominator-guard/value.pickle` | `f988dfc3960de2ff3ffd2322b911e20460698c5b0c24b12a93fd9e6c698981fc` |

`denominator-observation.json` contains full expression text and `srepr`,
selected order, flag values and Python types. `operation-index.json` and
`artifact-index.json` supply all remaining operand/return routes and hashes.
The original incomplete input is unchanged at SHA-256
`5e6d6a95b6e3336d3ab976ab82d29b742ca68aa3b0b65100f64126ae863fc1b4`.
Its original failed operation still has no completed return.

All **392 manifest records / 386 unique pinned files** match their prehashes,
posthashes and current bytes. This includes all 279 parent partial artifacts
and all 101 original source/input records. Fourteen launch snapshots match
their original sources. None of the 61 completed parent operations was
replayed.

The unchanged shared guard around the normalization supervisor enforced
900 seconds outside / 840 native seconds, 2 GiB memory, zero swap, one CPU,
nice 15, 32 tasks and one native thread. Whole-guard time was 5.075145 seconds;
peak memory was **70,053,888 bytes** (66.809 MiB). Four samples show zero
swap, zero memory-event counts and at most three tasks. Coordinator, guard,
child and supervisor completed successfully; acceptance of this limited
diagnostic rests on the saved evidence and checks above, not exit codes.

## Next specific dependency

Before any separately approved constructor continuation, implement a
source-bound exact nonzero certificate for the actual selected coefficient.
For this entry the displayed negative real part supplies the certificate:
extract from the saved coefficient, preserve the exact reconstruction and
sign/domain operands, and require them to establish nonzero. Do not replace
the failed guard with a structural comparison or a floating-point tolerance.
If an exact certificate cannot be established, stop with the evidence saved.

Any future continuation must restore the 61 parent complete returns and all
applicable diagnostic intermediates instead of recomputing them. In
particular, the numerator and denominator series now have completed returns.
This next correction/continuation is specified, **not implemented or launched**
by the diagnostic approval. No additional physics job or external review was
submitted. A fresh scope and pinned gate would be needed for that execution.

Literal historical reviews remain Claude **NEEDS REVISION** and Grok
**CLEAR FOR THIS BOUNDED INGREDIENT**, with the earlier local disposition;
there is no fresh independent CLEAR. Accepted reference/outgoing results,
all prior failures and incident history remain intact. The shared guard,
protected builder suffix, Lean paths and `S11_lean_*` work are untouched.
No newly generated tracked data exceeds 1 MiB. Full Green/FORM,
two-asymptote response, A11/A12, forcing-domain pairing, diagonal extension,
complex-frequency retarded equivalence and a radiating witness remain open.
