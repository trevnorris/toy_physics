# Ablation-driver conformance fixtures report

> ⛔⛔ **FULLY SUPERSEDED — USER DECISION, 2026-07-29/30.** The driver was re-scoped small and the whole
> conformance-fixture apparatus was cut: **no contract, no frozen fixtures, no three-session shape, no byte
> authority.** `fixtures_v4/` and `fixtures_v4.accepted.sha256` still exist on disk but gate nothing. Live
> spec: `REQUIREMENTS.md` (v5), live set **R1, R2, R3, R4, R6, R7**. Read this file as history only.

## ⚠ STATUS BANNER (added 2026-07-29) — read before using any part of this report

- **Written against requirements v1.** The accepted requirements are now **v4**
  (`REQUIREMENTS.md`, "Ablation driver — requirements and acceptance intent (v4)"). Everything below
  predates the v3 division-of-labour correction, the v4 scope-down that **withdrew R8's original replay
  property** — R8 still exists under that name, now stating the narrower honest-evidence rule "the
  evidence claims exactly what it supports" (`REQUIREMENTS.md:108`) — and `CONTRACT.md`.
- ⛔ **The suite it describes is SUPERSEDED and must be REDONE, not extended.** It was frozen against v1,
  covers only A1–A6/A8, and lives at `_scratch/ablation_driver/fixtures/` — **gitignored, not copied out,
  and not guaranteed to survive**. Nothing below should be treated as a frozen acceptance artifact for
  the driver that is now to be built.
- ⛔ **"A8" below is the v1 item, not v4's.** Below, A8 is *"the fixtures can fail"* — a **suite-quality**
  item. In v4, **A8 is the honest-evidence requirement**: "R8 holds — every claim in the run's committed
  output is checkable against the committed files alone, and each stated prerequisite is named rather
  than assumed." v4 also grades fixture-suite quality **separately**, explicitly not as driver
  acceptance. Do not read the table below as evidence for v4's A8.
- ✅ **What remains valuable is the "Untestable or underspecified as written" section** — its ten
  specification findings, which are about the requirements rather than about the suite. They are
  unedited and are the reason this report was kept.

## Scope

This session produced fixtures, producer/checker programs, expected outcomes, non-empty floors, and
failure exhibits only. It produced no driver, driver prototype, reference implementation, example
implementation, or acceptance substitute. No fixture was run against a driver because none exists in
this session.

The frozen suite is under `fixtures/`; `fixtures/MANIFEST.md` is the acceptance index and
`fixtures/SHA256SUMS` is its byte-level freeze. A1--A5 execute only in disposable fixture copies. A6
points directly at the real committed include-list and pins both its Git blob and parse fingerprint.

## Fixture coverage

### A1 — artifact reset is load-bearing

`a1_artifact_reset` has three ordered rows: nonzero/no-emission, zero/emission, then
nonzero/no-emission. It is run both without a starting artifact and with a planted residual. The
successful row emits a value distinct from `pristine/artifact.txt`; the failing rows report a PASS/FAIL
tally and named first failure but emit no artifact. The checker returns different verdict text and exit
status for absent, fresh, and residual bytes.

The floor requires both trials, all six applied rows, exact exit/emission patterns, non-empty captures,
absence at both failing rows, movement at the emitting row, and claim scope consistent with R9. A seeded
cache sentinel also makes R3's cache clearing observable.

Direct component exercises established the required shapes:

- failing producer: exit 23, 1 PASS/1 FAIL, first failure
  `INTENTIONAL_PRODUCER_FAILURE`, no artifact;
- successful producer: exit 0, 2 PASS/0 FAIL, exact artifact `artifact=FRESH_EMIT\n`;
- checker with no artifact: exit 12, `CHECKER|verdict=ARTIFACT_ABSENT`;
- checker with fresh artifact: exit 0, `CHECKER|verdict=FRESH_EMIT`.

### A2 — exact-once mutation floor

`a2_mutation_floor` supplies one old text occurring zero times and another occurring twice. The
producer and checker create invocation journals if they are called, so an incorrect attempt to run
either row leaves evidence.

The floor is two distinguishable application-outcome rows with observed counts 0 and 2, zero completed
mutants, and no invocation journals. Thus a no-op, silent skip, outcome collapse, or replace-all behavior
cannot pass.

### A3 — exact resume

`a3_resume_exact` supplies three deterministic rows and a producer completion journal. Its second
producer can pause before completing, giving an external conformance run a stable interrupted state
after row `one` is banked. The resumed result bytes must exactly equal an uninterrupted reference; the
completion journal must be exact text `one\ntwo\nthree\n`.

The floor requires a non-empty interrupted prefix, all three final rows in include-list order, unique
keys, exact result bytes, non-empty captures, and exactly one producer completion per key.

### A4 — restore on clean and interrupted paths

`a4_restore` pins the pristine target SHA-256 and makes its producer record the SHA-256 of the mutated
target bytes it actually observed. The clean trial completes normally. The interrupted trial pauses
after that marker and is terminated with SIGTERM.

The floor requires both mutation evidence and subsequent byte/digest equality. A target that was never
mutated cannot pass merely by already having the pristine digest.

### A5 — captures are not clobbered

`a5_captures` emits unique producer stdout, producer stderr, checker stdout, and checker stderr for each
of `red`, `green`, and `blue`. `EXPECTED.md` pins all twelve payloads byte-for-byte.

The floor is three completed rows and twelve non-empty mutant-and-stream-specific files. Empty files,
reconstructed output, duplicated payloads, combined streams, and final-row-only captures all fail.

### A6 — real include-list parses

`a6_real_include_list` names the required committed path directly; there is no synthetic copy. The
current worktree file matches Git blob `3124d4842d6f6f0a6e7d6d0865397345aa6e867d`. Independent TSV
parsing measured 51 data rows, five fields per row, and a 6,942-byte canonical JSONL stream with
SHA-256 `666e0a6ca0db22fd48bb74d374172040280dbc3e3d5c3004bda8bd6103f7627e`.

The floor requires all 51 losslessly parsed rows, not a successful open or line count. This fixture
makes no claim about applying the production mutations or about the production corpus.

### A8 — the fixtures can fail

`a8_failure_exhibits/FAILURE_EXHIBITS.md` freezes thirteen deficient cases and the discriminator each must
trip. I exercised them in disposable temporary directories or, for deliberately malformed result
states, evaluated the planted observation against the frozen exact floor. Every case was rejected:

| Case | Measured deficient observation | Frozen discriminator that failed |
|---|---|---|
| `F-A1-preexisting` | producer 23; checker 9 `RESIDUAL_ARTIFACT` | checker must be 12 `ARTIFACT_ABSENT` |
| `F-A1-between` | producer 23; checker 0 `FRESH_EMIT` | artifact must be absent after the non-emitting row |
| `F-R3-cache` | producer 31; `STALE_CACHE_PRESENT` | expected per-row producer exit/tally and cache absence |
| `F-R9-overclaim` | producer 23, no artifact, checker 12; planted claim `detected` | absent-artifact row cannot claim checker detection |
| `F-A2-zero-skip` | one observed row instead of two | exact two-row key floor |
| `F-A2-collapse` | both counts absent; both outcomes `NOT_APPLIED` | required counts 0/2 and distinguishable outcomes |
| `F-A2-replace-all` | producer 97; invocation journal present | producer must not run; zero completed mutants |
| `F-A3-rerun` | journal `one,one,two,three` | exact completion journal and unique-key floor |
| `F-A3-drop` | final keys only `two,three` | three-row floor and exact uninterrupted equality |
| `F-A4-clean` | final target SHA-256 `f1c1c197…f142c98` | pristine SHA-256 is `3b612246…61fe9c` |
| `F-A4-interrupt` | mutation marker present; final target still `f1c1c197…f142c98` | interrupted cleanup must restore pristine bytes |
| `F-A5-clobber` | four files remain, all BLUE; RED/GREEN absent | twelve-file floor and per-key byte payloads |
| `F-A6-whitespace` | 51 rows tokenized to widths 7, 9, 11, or 13 | every parsed row must have five preserved TSV fields |

No persistent throwaway driver stubs were written. The failure demonstrations were disposable probes
that invoked the fixture programs or inspected planted deficient observations, and their temporary
directories were deleted by the temporary-directory owner after execution.

## Untestable or underspecified as written

These are specification findings. They are not fixture failures.

1. **No fixed driver invocation or result schema is specified.** The suite can freeze inputs and
   semantic outcomes, but it cannot name a single executable acceptance command or mechanically locate
   every R8 field in an arbitrary result file. R8's “readable from the commit” and R9's claim-scope
   judgment therefore require an agreed external invocation/evidence mapping before they can be fully
   automatic. The fixtures do not invent one.

2. **R1 and A6 disagree about where the target file is named.** R1 says every include-list row names the
   target file. The required real A6 list has exactly five columns
   `axis,key,record,old_text,new_text` and no target-file column. Lossless parsing is decidable and
   frozen here; R1's per-row target-file assertion is not decidable for that input without an additional
   run-level contract that the requirements do not state.

3. **A6 full execution conflicts with the synthetic/self-contained constraint.** The real list can be
   parsed without making a production claim. Applying its 51 mutations would additionally require a
   production target, producer, checker, and pristine artifact, which this directive expressly excludes
   from synthetic fixtures. No parse-only driver operation is required by R1--R9. Consequently this
   suite makes A6's parser compatibility decidable, but cannot force a full mutation run against that
   list without either extending the public contract or violating fixture isolation.

4. **A3 does not define an interruption boundary.** “Interrupted after some rows” does not name a signal,
   row limit, or observable stop point. The fixture chooses SIGTERM while row `two` is paused and requires
   row `one` already be flushed. The exact-resume property is decidable under that trial. The literal
   phrase “no row runs twice” is ambiguous for an in-flight row: a row begun but not banked must be
   attempted again to finish, while R5 clearly forbids re-running already banked rows. This fixture
   measures completed producer rows and banked result keys.

5. **R7's “including a run that ended early” is unbounded.** Cleanup after catchable SIGTERM is
   executable and frozen in A4. Cleanup after SIGKILL, machine loss, or interpreter destruction cannot
   be performed by the terminated process. The requirements need to bound the early-ending class if A4
   is intended to cover more than the SIGTERM trial.

6. **R3's complete internal ordering is not externally distinguishable.** A1 proves artifact deletion
   and cache clearing occurred before the producer, and A4 proves mutation was visible at the real target
   before its producer. It does not distinguish two internal orders that lead to the same observable
   bytes—for example, artifact deletion immediately before versus immediately after cache removal. An
   exact total-order acceptance requires an external event trace or watcher contract not specified
   here.

7. **R8 leaves several grammars undefined.** “Tally,” “named assertion,” “verdict text,” and “emitted
   value” have no serialization or extraction rules. A1 supplies unambiguous synthetic instances and
   exact expected facts, but a generic conformance checker cannot know whether differently formatted
   evidence carries all required information without reading and interpreting it.

8. **R8's commit-readability property cannot be demonstrated solely in a pre-commit scratch run.** It
   becomes decidable only against the exact set of files selected for a commit. This session was
   forbidden to commit, correctly, so it freezes what evidence must survive but does not claim the
   future commit contains it.

9. **A7 is intentionally not a synthetic fixture in this directive.** It is a production retrofit
   comparison, while this directive requested fixtures for A1--A6 and A8 and prohibited production
   corpus assertions. A7 is decidable only in the later production run against the committed stage023
   tables; no A7 result is claimed here.

10. **R9 has no exhaustive claim vocabulary.** The A1 oracle decisively rejects the named overclaim
    “detected” for absent-artifact rows and requires the underlying producer/artifact/checker facts.
    Semantically equivalent overclaims using other language still require human interpretation unless a
    controlled result vocabulary is added to the acceptance contract.
