# S11c-d matching-channel construction

The first variable-profile matching checkpoint is implemented in the existing
engine's `TwoEndedMatchingChannels`, committed at `0f5498ba`. The focused run
stopped after 1.58 seconds at the full-file normalization-helper pin, before any
channel calculation. Accepted end packets predate logging/pre-emission storage
changes. The instrument now requires exact AST agreement for every consumed
helper function, while still verifying full frozen snapshots and other source
pins. All eight end/function joins pass; the failed attempt is preserved under
repository `_scratch/s11c/s11c-matching-channels-20260914/complete`. No physical
constructor changed and no channel-current result is accepted yet. The fresh
`retry-02` run is active with its own completion/error watcher after compatibility
repair commit `2953fa48`.

The construction consumes the accepted LEFT/RIGHT current/adjoint packets,
preserves every isolated-root/lift basis and classifier/domain record, and
computes the full cross-mode physical current from the native polarized slab
and bulk operands. It emits field/flux bases, current matrices, coordinate maps,
source diagonal-block joins, basis-change and Hermitian residuals. Cross-root
current entries are calculated explicitly. The instrument checks every serialized
object and its units/grades, the full candidate/rank census and the orientation
sets against the accepted end records before conditional publication.

Removing the new class leaves the entire engine AST identical to both accepted
normalization snapshots. The physical producer calculations and Fourier reduction
are unchanged. The [matching plan](S11c_d_variable_profile_matching_plan.md)
continues with reduced local/nonlocal operator assembly, variable-profile boundary
matching and continuum re-expansion. No complete S-matrix or profile-frequency
pole result is supplied by this channel checkpoint. Existing exceptional-domain
limits and c2 operand debt remain.
