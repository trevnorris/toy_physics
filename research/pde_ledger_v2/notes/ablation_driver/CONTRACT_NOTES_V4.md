# Findings and decisions for ablation-driver requirements v4

## What v4 removed from the contract

1. **The public `verify` and `rerun` operations are gone.** The CLI now has only `validate`, `run`, and
   `repair`. With `rerun` went its special exit status, row selection, generated one-row list, separate
   output package, source-package verification, and caller target/artifact snapshots. A9 still has a
   complete truth table to apply; it does not need a driver verb that grades the driver's own evidence.

2. **Replay closure is gone rather than renamed.** The contract no longer requires executable or
   interpreter digests, imported-library closure, lockfiles, dependency archives, replay capsules,
   `environment.json`, public target/artifact entry snapshots, `SHA256SUMS`, or a package manifest.
   It makes no promise that a commit can re-run a row. External toolchain, OS, kernel, and hardware
   bytes are explicitly outside the commit set.

3. **The old C-6 replay procedure is gone.** The replacement answers v4's actual question: a completed
   run leaves the accepted config/list, results, per-row captures, restore proof, completion marker, and
   a short `EVIDENCE.md` scope statement for commit. That statement names the repository at the
   containing commit, the included configuration/environment, the toolchain command names, and POSIX
   process/filesystem/signal semantics as what a reader needs beyond the row. It also says that the
   package is not self-contained replay evidence.

4. **Findings 3, 4, and 11 from `CONTRACT_NOTES.md` no longer generate infrastructure.** Finding 3's
   replay-package closure and finding 4's platform/capsule boundary were consequences of withdrawn R8.
   Finding 11 remains the ordinary limitation of any finite public fixture suite; v4 supplies no
   generality-proof requirement, so the contract adds no opaque-list service or anti-special-case
   protocol.

## What remains, and why

5. **A closed child environment and a declared stable project read set remain for R5, not replay.**
   `stable_inputs` is the reduced successor to `replay_inputs`: it covers project files/data read by
   children, plus the driver automatically, and recovery state authenticates their entry digests so
   resume can reject drift. The digest inventory is not committed. It does not capture external
   executables, libraries, or a transitive platform. Direct argv, no ambient environment inheritance,
   no network/time/unseeded-random decisions, and stable toolchain behavior are the determinism
   conditions needed for resumed and uninterrupted result bytes to be equal. This is the load-bearing
   part identified by old finding 5.

6. **Exact config/list copies, prefix banking, stable capture paths, and `COMPLETE` remain for R5.**
   They decide list/config sameness, make an interrupted prefix recognizable, prevent a banked row from
   running again, and exclude attempt IDs/timestamps from result equality. They are resume state, not
   replay state.

7. **Authenticated recovery snapshots remain only under `recovery_path` for R7.** They are required to
   repair after an uncatchable process death and are expressly excluded from the commit set. Per-target
   and per-artifact SHA-256 rows in `restore.tsv` remain because C-7 and R7 require digest proof,
   including an answer for multiple targets.

8. **Raw banked captures and their stable paths/digests remain for R4/A5 and R9.** R4 requires that
   every mutant's captures survive without overwrite; A5 requires correspondence to that mutant; and
   C-5 derives tally, failure, verdict, moved-value, and emitted-artifact fields from those bytes.
   Interrupted attempts remain separate because overwriting them would recreate old finding 9.

9. **The fixed observation grammars and C-9 legacy mapping remain.** They are required by A7 and R9,
   not by replay. Without them the legacy tallies/statuses/moved values would again be implementation
   guesses, as old finding 8 explained. The outcome vocabulary remains deliberately exit-based so it
   does not claim what a checker read or that a mutant was scientifically caught.

10. **Target/mutant digests and exact application counts remain for R2, R6, and R7.** They bind
    applicability to the target's run-entry bytes, distinguish zero from multiple matches, and connect
    mutation/restoration claims to the real configured path. They are not retained as historical
    snapshots for replay.

## Remaining satisfiability and scope findings

11. **Literal R3 is still unsatisfiable for arbitrary unconstrained commands.** A child can read or
    mutate host state the driver cannot enumerate or isolate. The contract makes R3 meaningful over
    declared target/artifact/reset/evidence/recovery state and prohibits other mutable inputs. If R3
    means hermetic noninterference over the whole host, v4 still lacks the sandbox/process model needed
    to satisfy it.

12. **Literal unconditional R5 is still unsatisfiable for nondeterministic children or a changing
    toolchain.** The contract states deterministic-child and stable-toolchain behavior as conformance
    preconditions and mechanically checks the project inputs it can name. This is not a replay promise.
    If R5 must hold despite clocks, randomness, network state, or changed interpreter/library behavior,
    no driver at this scope can provide it.

13. **R7's uncatchable class still ends at surviving recovery storage.** `repair` can recover after
    process death only when authenticated recovery state survives. Host/disk loss that removes that
    state remains outside the property.

14. **Revised R8 and A8 are satisfiable and non-vacuous.** `EVIDENCE.md` narrows the claims, named
    captures and commit-tree digests make their support inspectable, and A8 can fail on a missing file,
    a changed commit target/input, or an overclaim. It no longer asks execution to prove replay.

15. **The C-9/manifest wording conflict — ✅ RESOLVED 2026-07-29, after this note was written.** A7 is
    exact on every mapped field, not byte equality with the legacy 13-column table.
    `DIMENSION_REWRITE.md` §12b(b) now says precisely that — **"exactly on every mapped field"**, "⚠ Not
    byte equality, which was never achievable", with the mapping frozen as part of the contract — so the
    wording correction C-9 required has been made; the committed legacy table remains untouched.

No v4 requirement is vacuous under these readings. A9 remains a real schema rejection, R8 remains a
real commit-content/prose gate, and the other acceptance legs retain their able-to-fail predicates.

## Decisions not expressly requested by §C

16. The producer may not mutate a target, and a checker may not mutate a target or declared artifact;
    artifact observation occurs before checkers. This preserves attribution of producer output (old
    finding 7).

17. Interrupted, unbanked captures use a namespace separate from final banked captures and are not part
    of the required commit set. This resolves the R4/R5 collision identified by old finding 9.

18. A9 consumers apply C-4/C-5 directly; the driver does not expose a self-verifying operation. This
    keeps the three-session fixture/implementation separation intact.

19. The pre-driver stage023 evidence remains only the A7 mapped oracle and is not retroactively failed
    for lacking files that this future driver will leave. Its prose must continue to describe its own
    omissions honestly.

## Size observation

20. The v4 contract is still **601 lines** (measured 2026-07-29, after the confirm-pass corrections that
    added lines; it was 597 when this finding was written), down from 618 but still materially larger than the
    manifest's informal “~100-line” description. The replay surface is gone; the remaining size comes
    from requirements that v4 retains: three exact operations, deterministic resume and banking,
    crash repair with multi-path digest proof, list/result schemas, raw-capture correspondence,
    observation extraction, a total outcome truth table, and the exhaustive stage023 mapping. This is
    an observation, not a line-count target. A roughly 100-line implementation may still be possible
    only if it realizes all of that behavior; the contract does not weaken R1–R7 or R9 to make the
    estimate true.
