# S11c-d matching-channel construction

The first variable-profile matching checkpoint is implemented in the existing
engine's `TwoEndedMatchingChannels`, committed at `0f5498ba`. The focused run
is active under repository `_scratch/s11c/s11c-matching-channels-20260914/`, with
a local completion/error watcher. No channel-current result is accepted yet.

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
