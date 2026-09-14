# S11c-d matching channels

Both physical-end channel packets are validated for the supplied LAB_HELD /
RHO4_CONSTANT instance. Each retains all 18 isolated-root/lift candidates and
22 basis directions, including every closed, wrong-sheet and unresolved record.
Each end has two current-normalized two-dimensional subspaces: two incoming and
two outgoing basis directions. Classifier operands remain explicit; no sector
label is inferred from channel availability.

The engine evaluated all four eligible mode-subspace pairs at each end using
the polarized slab current and the computed convergent bulk-depth integral.
Field-basis and flux-basis currents, the basis map and both reconstruction routes
are emitted separately. The largest basis-change residual norm is 6.33e-16,
the largest accepted source diagonal-block difference is 4.82e-16, and the
largest Hermitian residual is 3.64e-16. Cross-root current norms are 2.83e-17
at LEFT and 3.53e-17 at RIGHT; these were computed, not assigned zero. Every
full-rank/nullity and incoming/outgoing orientation-set residual is zero.

The [checkpoint](S11c_d_matching_channels_checkpoint.json) records 500 tags,
2,763 metadata paths and 248 fresh write-keys. Every computed payload,
fingerprint, unit and grade was replayed against the transcript. All source,
source-function, end-input and artifact joins passed. The 309,046-byte
[transcript](../scripts/out/S11c_d_matching_channels.out) has SHA-256
`92e7750c02cfe95efef5ad74be7f01db16e8555ab5abba65a0946252f41d0bf6`.
Construction/validation took 748.91 seconds including supervision, at 282,156
KiB peak RSS, with exit zero and empty stderr. The local watcher delivered its
completion event without model polling.

The first attempt stopped after 1.58 seconds before calculating a channel,
because the accepted end packets predated logging/pre-emission storage changes.
Its source snapshots are preserved. The repaired check requires exact agreement
of the four functions actually consumed at each end, in addition to historical
snapshot and all other source pins. All eight function joins pass. Removing the
new matching class leaves the entire engine AST identical to the accepted end
snapshots; no physical constructor or Fourier reduction changed.

This completes steps 1–2 of the [matching plan](S11c_d_variable_profile_matching_plan.md).
The next construction is the source-driven local/nonlocal operator assembly.
There is no interior solve, complete S-matrix, continuum re-expansion or
profile-frequency pole result yet. Existing exceptional-domain limits and c2
operand debt remain.

Publication is committed at `3a5d3d25`; the canonical 309,046-byte annex
payload SHA-256 is `92e7750c02cfe95efef5ad74be7f01db16e8555ab5abba65a0946252f41d0bf6`.
