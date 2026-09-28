# Reference-kernel implementation review disposition

Both fixed-packet reviews finished before this adjudication or any worker
edit. Claude's literal verdict is **CLEAR FOR THIS BOUNDED INGREDIENT**.
Grok's literal final `text` verdict is **NEEDS REVISION**, with one blocking
finding about extracting LU column entries. Grok's separate `thought` field
is not its final report. Neither verdict is rewritten as a new reviewer verdict.

The original 16-file packet, archive, prompts, reports, stderr and run receipts
remain at `_scratch/s11c/s11c-d-reference-kernel-20260926/build-review`.
Their hashes and all reviewed source hashes matched before the local edit.
Both processes exited zero; that alone supplied no clearance. Claude stderr
was empty. Grok stderr contained configuration warnings and unsuccessful Read
tool calls during file discovery. Its final report cites the supplied worker,
source excerpts, guard and conventions and identifies a specific code line.
The report's factual claims were adjudicated against local source, not accepted
from its exit code or repeated verdict text.

## LU finding and bounded disposition

Grok says `columns[j][i]` returns a 1-by-1 matrix and requests
`columns[j][i, 0]`. The actual `/usr/bin/python3` package resolves to
`/home/trevnorris/.local/lib/python3.10/site-packages/mpmath`.
Its `matrices/matrices.py`, `_matrix.__getitem__`, converts an integer index
on a one-column vector to `(key, 0)` and returns a scalar when neither index
is a slice. Thus the reported failure does not occur for this package's
column vectors. Reading this definition did not import a scientific library
or run a calculation.

Nevertheless, the exact one-line revision requested is adopted: the worker
now spells out `columns[j][i, 0]`. It is equivalent for the installed package,
preserves the column-solve orientation, and makes the scalar extraction
explicit. No other worker behavior or physical input changed. An AST check
and exact one-line comparison against the preserved reviewed worker verify
this limited delta. The actual source-bound LU comparisons remain required
inside the single guarded production job. This closes the sole blocking
finding locally; it is not a fresh independent CLEAR verdict.

Grok's optional uncertainty about the `reductionState` wrapper is resolved by
`save_reduced_action_cache` in the pinned audit source: its saved packet has
the explicit `reductionState` key containing `dict(vars(reduction))`. The
existing packet/hash joins and runtime checks still govern restoration.

## Method boundary and optional findings

Both final reports accept the limited ingredient stop. The Fourier convention
and source branch data do not by themselves supply the full coupled real-axis
singular integral. The regular inverse and spectral density may therefore be
constructed under the approved first-stage scope. They do not complete the
outgoing Green operator, FORM, or A11/A12.

Claude identifies continuation from the selected radiation branch with the
`exp(-i omega t)` convention, saved REFERENCE normal-mode/current evidence,
and behavior of the coupled pencil along a limiting-absorption path as
dependencies to inspect for a later outgoing prescription. This is recorded
as a review suggestion for that later method decision, not as a proven contour,
a compulsory new zero census, or a global stability theorem. No such
operation enters this worker. Its determinant-nonzero domain and pending
prescription status remain explicit.

The other observations do not block this stage: the outer timeout retains
operation receipts and its own outcome even if the native failure receipt
is absent; the decoder is only a producer-import restriction for hash-pinned
local packets, not a security sandbox; norm thresholds are specific to the
saved unit frame; exact physical-input/checkpoint paths are already manifest
members; branch evaluation and additional joins are optional future evidence.
The profile Abel mass is a separate operand and is not equated to the normal
Fourier mass. A small mutation response or ill-conditioned fixed probe must
stop without selecting an easier replacement or retrying.

The scope remains one guarded job, 900 seconds, 2 GiB, zero swap, one CPU,
nice 15, 32 tasks and one native thread, using the existing supervisor and
local completion hook. Scientific acceptance still requires the actual
resource records, strict stderr, source posthashes, persisted operands and
residuals. No accepted scientific result or Lean file was changed.
