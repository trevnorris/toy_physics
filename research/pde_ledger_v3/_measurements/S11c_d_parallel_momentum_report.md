# S11c-d parallel numerical execution

The four-worker execution preflight is accepted. Recovery exited zero with empty
stderr in 57.4 seconds, reusing all numerical operands. All four 65,536-node
prefixes agree exactly in values, measure mutations, mass and counters; literal
native frequency-prefix replay agrees. All six complete underresolved records,
native terms and full actions agree exactly. This validates scheduling only;
no source/profile refinement result or physical limit has yet been accepted.

The initial preflight reused a payload encoder across independent output files.
Two encoder resets repair it; the full preflight AST rejoins after their removal.
The original second stream still fails standalone decoding as required, and its
combined-stream payloads agree with both recovered independent transcripts.
Both recovered transcripts are byte-identical: 1,530,835 bytes, 5,754 tags,
2,875 write keys and 26,705 metadata paths. Every source, original artifact and
copied packet hash passes. The engine, codec, serial constructor/emitter and
numerical worker implementation remain unchanged.

One workspace difference is fully accounted for: the reused serial process
reports 12,677,760 bytes accumulated from earlier grids while that isolated
worker reports 6,386,304 bytes. Each serial peak equals the running maximum of
the worker peaks. Per-worker measured RSS was at most 122,340 KiB (119.5 MiB).
Parallel prefix elapsed time was 17.09 seconds. Serial file timestamps indicate
about 30.07 seconds, or approximately 1.76x faster; the serial monotonic timer
was not saved, so this is an approximate comparison, not an exact benchmark or
production-duration prediction.

The earlier resume test reproduces an uninterrupted 65,536-node prefix exactly
from node 16,384, including frequency counts; an altered count is rejected.
The authorized serial handoff preserved 709 partials and copied the latest
11,616,256-node accumulator byte-for-byte. Original logs and snapshots remain.
The old watcher is disarmed and that serial run must not be restarted.

Next: launch four production workers, one per field/source-profile setting,
with one native thread and a 2 GiB address-space ceiling each. Retain the native
summation order, deterministic aggregation, independent logs/checkpoints and
silent completion/error wake-up. The accepted focused transcripts remain in
repository scratch; publish the production physical transcript only after its
full validation. See the acceptance and handoff checkpoints.
