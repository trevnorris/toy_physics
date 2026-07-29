# Findings against ablation-driver requirements v3

1. **§C-9 / A7 has a live normative contradiction.** The new fixed result schema cannot be byte-equal
   to the committed legacy 13-column schema. `DIMENSION_REWRITE.md` §12b(b)'s word “exactly” must be
   corrected to “exactly under CONTRACT.md §C-9's mapped projection.” The legacy file must remain
   unchanged.

   ✅ **RESOLVED 2026-07-29 — marker only; the finding above is unedited and reads in its original
   present tense.** The manifest was corrected: §12b(b) now requires the retrofit to reproduce the
   committed table **“exactly on every mapped field”** and states **“⚠ Not byte equality, which was
   never achievable”**, with the legacy→new column mapping frozen in `CONTRACT.md` §C-9 rather than
   chosen at comparison time (`manifests/DIMENSION_REWRITE.md:1344-1349`). The contradiction is closed;
   the "must be corrected" sentence above is a record of what was demanded, not live work.

2. **R3 is unsatisfiable literally for arbitrary commands.** No driver can prove that *no* prior state
   anywhere on the host is observable by an unconstrained producer. The contract makes the useful
   property decidable for declared artifacts/reset paths and prohibits undeclared mutable inputs, but
   fixture success proves only the planted state classes. If R3 is intended literally rather than as a
   closed-read-set contract, it requires a hermetic sandbox/filesystem/process model not named in v3.

3. **R8 is not a result-schema property.** It requires raw captures, entry target/artifact bytes, config,
   include-list, dependency bytes, environment prerequisites, restore proof, and the driver—not merely
   a row. Even then, completeness of a declared transitive read set is not mechanically provable
   without tracing or hermetic execution. The contract therefore makes omissions fail `rerun`, and A8
   must be performed against the actual committed tree. A pre-commit scratch fixture cannot establish
   this property.

   ✅ **OVERTAKEN by v4 — marker only; the finding above is unedited.** R8's original replay property
   was withdrawn by user decision (`REQUIREMENTS.md:8`), and v4 names this finding's replay-package
   closure as a consequence of that withdrawal, generating no infrastructure
   (`CONTRACT_NOTES_V4.md:23-24`). The `rerun` verb it turns on no longer exists
   (`CONTRACT_NOTES_V4.md:5-8`). Revised R8/A8 are now judged satisfiable and non-vacuous
   (`CONTRACT_NOTES_V4.md:83-85`; `CONTRACT.md:461-465`).

4. **Literal “committed evidence alone” has a platform boundary.** A commit cannot contain the reader's
   kernel/hardware, and proprietary or host-only runtimes cannot be replayed merely from version text.
   CONTRACT.md requires sensitive executable state to be committed as a replay capsule and treats
   irreducible host services as explicit prerequisites. If “alone” excludes even those prerequisites,
   R8 is unsatisfiable.

   ✅ **OVERTAKEN by v4 — marker only; the finding above is unedited.** v4 names this finding's
   platform/capsule boundary as a consequence of withdrawn R8 (`CONTRACT_NOTES_V4.md:23-24`). Replay
   capsules are **gone rather than renamed**, and external toolchain, OS, kernel and hardware bytes are
   explicitly outside the commit set (`CONTRACT_NOTES_V4.md:10-14`; `CONTRACT.md` §C-6, which now states
   the commit set and names prerequisites rather than promising replay).

5. **R5 needs determinism outside the driver.** Byte-identical result rows are achievable because the
   schema excludes clocks, PIDs, attempt IDs, and absolute paths. Reaching the same outcome after resume
   still cannot be guaranteed for a producer/checker that reads time, randomness, network state, or an
   undeclared file. The closed environment/read-set rule is therefore necessary, not optional metadata.

6. **R7's uncatchable class must mean surviving-storage process death.** `repair` can recover after
   `SIGKILL` only if its authenticated recovery state survives. It cannot repair host/disk loss that
   destroys that state. If v3 intends host loss too, R7 needs separately stored recovery evidence and an
   availability requirement.

7. **v3 did not ask who may mutate declared artifacts or targets.** Without a rule, a checker could
   create the artifact later attributed to the producer, or a producer could rewrite its mutant. The
   contract therefore observes artifacts before checkers and makes target/checker artifact mutation an
   integrity failure.

   ✅ **SETTLED — marker only; the finding above is unedited.** The rule it asked for was adopted and
   v4 retained it: the producer may not mutate a target, a checker may not mutate a target or declared
   artifact, and artifact observation occurs before checkers — recorded as preserving attribution of
   producer output, "old finding 7" (`CONTRACT_NOTES_V4.md:98-100`, decision 16).

8. **v3 did not ask how extraction claims become fields.** Tallies, first failures, checker statuses,
   mismatch names, and moved values are otherwise transcriptions with no grammar. The contract adds
   fixed `exit-only`, ledger tally/result, and ledger dimension-record observation grammars. Supporting a
   new prose format requires a future contract version, not an implementation guess.

   ✅ **SETTLED — marker only; the finding above is unedited.** The fixed observation grammars and the
   C-9 legacy mapping were added and **v4 keeps them**, required by A7/R9 and not by replay — v4 cites
   "old finding 8" as the reason (`CONTRACT_NOTES_V4.md:55-58`). They live in `CONTRACT.md` §C-5
   (outcome truth table) and §C-9 (stage023 legacy mapping).

9. **v3 did not define interrupted, unbanked captures.** R4 (“never overwritten”) and R5 (exact result
   bytes) conflict if a resumed attempt reuses a partial capture path or embeds a fresh attempt ID. The
   contract fixes stable banked paths, excludes attempt IDs from results, and requires interrupted
   attempt bytes to be retained separately.

   ✅ **RESOLVED — marker only; the finding above is unedited.** v4 states it in terms: interrupted,
   unbanked captures use a namespace separate from final banked captures and are not part of the
   required commit set, which "resolves the R4/R5 collision identified by old finding 9"
   (`CONTRACT_NOTES_V4.md:102-103`, decision 17; `CONTRACT.md:358,518`).

10. **The existing stage023 evidence is deliberately not R8 evidence.** Its summary states that the raw
    captures, run log, pristine run target, and driver were not committed. The legacy table remains a
    valid A7 mapped oracle; its missing source captures cannot be repaired by normalization or prose.

    ✅ **SETTLED by v4 decision — marker only; the finding above is unedited.** The pre-driver stage023
    evidence remains only the A7 mapped oracle and is **not retroactively failed** for lacking files
    the future driver will leave; its prose must keep describing its own omissions honestly
    (`CONTRACT_NOTES_V4.md:108-110`, decision 19). With R8's replay property withdrawn there is no
    longer a standard it fails. The oracle's exactness condition is stated at
    `manifests/DIMENSION_REWRITE.md:1344-1349`.

11. **Finite visible fixtures cannot exclude an input-hash special case.** A build can dispatch on the
    frozen fixture/list bytes and fabricate their public evidence while rejecting unseen conforming
    inputs. The contract makes such a build nonconforming, but the stated acceptance set alone cannot
    distinguish it. An opaque post-build conforming list or code review is required if this threat is
    in scope; no §C wire-format decision can make a finite known suite prove generality.

    ⚖ **DECIDED, LIMITATION STANDS — marker only; the finding above is unedited.** v4 took the
    decision: **no opaque-list service and no anti-special-case protocol** is added, because v4 supplies
    no generality-proof requirement (`CONTRACT_NOTES_V4.md:23-27`, which also records that "finding 11
    remains the ordinary limitation of any finite public fixture suite"). ⚠ The consequence of that
    decision is that **the limitation remains open and unmitigated, by choice** — not closed. What
    makes it acceptable, so a future reader can re-open it on evidence rather than on unease: it is the
    ordinary limitation of *any* finite public fixture suite and is **not specific to this
    deliverable**; and the standing mitigation is only the three-session separation — fixtures are
    frozen by a session that does not build the driver, and the building session may neither author nor
    weaken them (`REQUIREMENTS.md:27`; `manifests/DIMENSION_REWRITE.md:1333`). ⛔ **That separation does
    not exclude the threat above:** the frozen inputs are visible to the building session, so a build
    can still dispatch on their bytes. Re-opening this needs an opaque post-build conforming list or
    code review, as the finding says — not a further wire-format decision.

12. **The manifest's “~100-line” description is not compatible with treating it as a cap.** Exact
    resume, crash repair, five public operations, schema verification, capture hashing, observation
    extraction, multi-target proof, and R8 replay closure are all required behavior. If “~100-line” is
    normative rather than an informal size estimate, it conflicts with v3's accepted scope and should
    be removed rather than used to weaken the contract.

    ✅ **SETTLED — marker only; the finding above is unedited.** "~100 lines" is now recorded as
    non-normative wherever it appears: v4 calls the contract's own length "an observation, not a
    line-count target" and refuses to weaken R1–R7/R9 to make the estimate true
    (`CONTRACT_NOTES_V4.md:112-122`, size observation 20); the manifest records that the survey's
    sizings "are not budgets, and the driver's has not held"
    (`manifests/DIMENSION_REWRITE.md:1434-1437`); and `STATUS.md:69-71` states the "~100 lines" sized
    the hand-rolled loops the driver replaces, not the deliverable.
