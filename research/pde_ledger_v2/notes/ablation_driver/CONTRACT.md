# Ablation driver public contract

> ⛔⛔ **SUPERSEDED — USER DECISION, 2026-07-29/30. This contract is NOT a specification and NOT a gate.**
> The driver was re-scoped to *mutate a declaration, confirm the declared assert fires, record it,
> reviewed by one fresh agent* — **no contract, no frozen fixtures, no three-session shape.** The live
> specification is `REQUIREMENTS.md` (v5), whose live requirement set is **R1, R2, R3, R4, R6, R7**.
> ⛔ Do not implement to this document, and do not cite it as binding. Nothing here obliges a builder:
> the five operations, the exit-code table, the JCS/TSV wire formats, the outcome truth table, the
> banked-prefix resume semantics, the restoration-event proof and the evidence commit set were all cut.
> ⭐ **ONE PART IS STILL WORTH READING — §C-9, the stage023 legacy→new field mapping table — as a REFERENCE,
> not an authority.** ⛔ **Nothing here is authoritative, §C-9 included.** The live spec (`REQUIREMENTS.md`
> A7) downgraded the retrofit from a blocking gate to a cross-check and explicitly **does not require**
> field-by-field agreement via this table; it calls it "the useful reference for which legacy column
> corresponds to which new observation". Read it that way: the *mapping* is informative, the *schema it maps
> into* is gone, and a disagreement is a finding to report, never something either side is adjusted into.
> ⚠ Retained as history because it is the concrete measure of what the tooling-defence tower cost:
> ~600 lines of interface specification, none of which could catch a wrong derivation.

**Contract version:** `ablation-driver-v1`

**Applies to:** ~~`REQUIREMENTS.md` v4, 2026-07-29~~ — superseded by v5, 2026-07-29/30 (see banner)

This document fixes only the driver's public invocation, accepted inputs, emitted evidence, and the
meaning of that evidence. It does not prescribe the driver's implementation or internal file layout.
The words **must**, **must not**, and **exactly** are normative.

All project paths governed below are POSIX paths relative to the project root. Absolute project paths,
`..` path components, symlinks in a target/artifact/reset path, tabs, CR, NUL, and paths that escape the
project root are invalid. An external executable named through `PATH` is not a project path and is not
copied into run evidence. All text control files are UTF-8 without a BOM and use LF line endings.
SHA-256 means the SHA-256 digest of the exact file bytes, rendered as 64 lowercase hexadecimal
characters.

JSON values use RFC 8785 JSON Canonicalization Scheme (JCS). Thus a JSON-valued TSV cell has one
byte representation. A standalone JSON file is its compact JCS value followed by exactly one LF; a
JSON-valued TSV cell is bare JCS with no surrounding whitespace. JSON `null` is not an inapplicability
marker. Indented JSON examples below are expanded only for readability.

## C-1. Operations, argv, output, and exit status

The public executable is
`research/pde_ledger_v2/scripts/ablation_driver.py`, invoked from the project root. It accepts exactly
these three forms:

```text
research/pde_ledger_v2/scripts/ablation_driver.py validate --config CONFIG.json
research/pde_ledger_v2/scripts/ablation_driver.py run --config CONFIG.json
research/pde_ledger_v2/scripts/ablation_driver.py repair --config CONFIG.json
```

There are no positional aliases, short options, optional flags, or environment-variable substitutes
for these arguments. Every successful operation writes exactly one JCS JSON value plus LF to stdout
and writes nothing to stderr. A failure writes no success JSON; its stderr diagnostic is human-readable
but non-normative.

`validate` validates the configuration and losslessly parses the whole include-list without mutating
anything or creating evidence/recovery state. On success it writes one JCS JSON value plus LF to stdout:

```json
{
  "axis_counts": {"AXIS": 1},
  "columns": ["axis", "key", "record", "old_text", "new_text"],
  "list_sha256": "64-lowercase-hex",
  "row_count": 1,
  "rows": [
    {
      "cells": ["AXIS", "key", "record", "old bytes as text", "new bytes as text"],
      "row_id": "target:AXIS:key",
      "target": "target"
    }
  ],
  "schema": "ablation-list-validation-v1"
}
```

`axis_counts` has one member for every distinct `axis`, `columns` and every `cells` array retain input
header/column order, and `rows` retain data-row order. This output is the A6 observable; a success-only
or hash-whitelist validator is nonconforming.

`run` starts a new run when the configured evidence directory does not exist. When it contains a
matching, structurally valid incomplete run, the same invocation resumes it. Any other pre-existing
evidence directory is an integrity error. A successful run emits every include row, restores the
restoration set, writes C-6's `COMPLETE` marker, and writes a JCS summary containing
`schema="ablation-run-summary-v1"`, project-relative string `evidence`, unsigned integer `rows`,
`outcome_counts` with exactly the five C-5 tokens and unsigned integer values, and
`restored=true`. Those are its only members. Row outcomes never by themselves make `run` fail.
`run` does not silently stand in for `repair`: an incomplete run with an active attempt and no terminal
restore event exits `66` until `repair` records one, even if the target bytes happen already to match.

`repair` uses the named configuration and its surviving recovery state to restore every member of the
restoration set after an uncatchable termination. It is idempotent while valid recovery state remains:
an already-restored member is verified, not treated as an error. It appends repair proof to
`restore.tsv` and writes a JCS summary having exactly `schema="ablation-repair-summary-v1"`,
project-relative string `evidence`, unsigned integer `restore_event`, and `restored=true`.

Exit statuses have exactly these meanings:

| Code | Meaning |
|---:|---|
| `0` | The requested operation completed and all checks belonging to that operation passed. |
| `64` | Invalid operation, argv, or option spelling. |
| `65` | Invalid config, include-list, result truth-table row, or observation grammar. |
| `66` | Evidence/recovery digest mismatch, changed stable input, or unsafe resume. |
| `70` | A child could not be started or an evidence/capture operation failed; the observable restoration paths were restored. |
| `75` | A catchable signal interrupted work; the observable restoration paths were restored and the in-flight row is unbanked. |
| `76` | Restoration or its verification failed. This status takes precedence over every other failure. |

If the process itself is uncatchably killed, it naturally cannot emit a contractual exit status;
`repair` is then the contractual next operation. Child timeout is recorded as exit `124`; a child
terminated by signal `N` is recorded as `128+N`. A child start failure is not a child exit and therefore
does not bank a row.

Every `run` exit after recovery state has been established first attempts to restore the complete
restoration set and append the corresponding C-7 event. A failed restoration changes the exit to `76`.

There is deliberately no public replay or evidence-verification operation. Fixtures and readers apply
C-4 through C-9 directly. The withdrawn re-runnable-from-the-commit property does not reappear as
another verb.

*Choice and tradeoff — fixed verbs, observable validation output, and disjoint exit meanings make every
acceptance leg executable; they foreclose alternate CLIs and treating mutant rejection as driver
failure.*

## C-2. Run configuration and execution semantics

`CONFIG.json` is one JCS JSON object with exactly these members; unknown or missing members are invalid:

```json
{
  "artifacts": [
    {
      "name": "sidecar",
      "observation_format": "ledger-dimension-v1",
      "path": "path/to/sidecar"
    }
  ],
  "checkers": [
    {
      "argv": ["python3", "path/to/checker.py"],
      "name": "dimensions",
      "report_format": "ledger-result-v1",
      "timeout_seconds": 600
    }
  ],
  "contract": "ablation-driver-v1",
  "environment": {"LC_ALL": "C", "PATH": "/usr/bin:/bin", "TZ": "UTC"},
  "evidence_path": "path/to/evidence",
  "include_list": "path/to/include_list.tsv",
  "producer": {
    "argv": ["python3", "path/to/producer.py"],
    "report_format": "ledger-tally-v1",
    "timeout_seconds": 600
  },
  "recovery_path": "path/to/recovery-state",
  "reset_paths": ["path/to/__pycache__"],
  "stable_inputs": ["path/to/input-or-directory"],
  "targets": [{"name": "stage", "path": "path/to/target.py"}],
  "workdir": "."
}
```

The result path is exactly `evidence_path/results.tsv`; the banked capture root is exactly
`evidence_path/captures`; restore proof is exactly `evidence_path/restore.tsv`. They have no separate
override. `recovery_path` is mutable crash-recovery state and is never a substitute for committed
evidence.

Names use `[A-Za-z0-9_.-]+`; target, checker, and artifact names are unique in their respective arrays.
There is at least one target and one checker. `artifacts`, `reset_paths`, and `stable_inputs` may be
empty. Target, artifact, reset, evidence, and recovery paths are mutually distinct and non-overlapping.
`workdir` may contain any of them. A stable-input directory may contain producer, checker, driver, or
data paths, but must not contain a target, artifact, reset, evidence, or recovery path. A declared
artifact is absent or a regular file at run entry and may only become absent or a regular file during a
row.

Every `argv` is a nonempty JSON array of nonempty strings. It is executed directly, never through a
shell: no quoting pass, globbing, variable expansion, pipelines, redirection, or command substitution
is performed. Relative argv paths resolve from `workdir`. The child environment is exactly the
configured `environment`; the driver does not inherit undeclared variables. Environment keys and
values are strings without NUL, keys are nonempty and contain no `=`, and the configuration must
provide any required `PATH`. `timeout_seconds` is an integer from 1 through 600 inclusive.

The producer runs once for an exactly-applied row. Every checker then runs once, in configuration-array
order, even when the producer exits nonzero. Producer and checker stdout and stderr are separate byte
streams. A checker may inspect any declared input available to its argv—target, stdout made available
by its own command arrangement, an artifact, or another committed input. No outcome rule assumes that
a checker reads an artifact.

Before every producer, including the first and a resumed row:

- the selected target has its run-entry bytes except for exactly that row's replacement;
- every declared artifact is absent;
- every reset path is absent; and
- no banked row is invoked again.

The selected target remains at its real configured path throughout producer/checker execution. The
producer and checkers must not change a target. Artifact state is observed immediately after the
producer and before the first checker; checkers must not change a declared artifact. A changed target
or artifact is an operational/integrity error and the row is unbanked.

At `run` entry, every target and every then-existing artifact is snapshotted byte-for-byte into
authenticated recovery state before the first mutation. Targets and entry artifacts form the
**restoration set**. These snapshots are operational recovery material, not committed replay material.
`reset_paths` are explicitly declared disposable cache state: they are removed before each producer
and are absent after restoration. A path not declared as a target, artifact, reset, evidence, or
recovery path is not authorized as mutable run state.

`report_format` is either `exit-only`, `ledger-tally-v1` for a producer, or
`ledger-result-v1` for a checker:

- `exit-only` makes no textual tally/verdict claim.
- `ledger-tally-v1` recognizes at most one UTF-8 stdout line of the exact form
  `TALLY LABEL: P pass + F fail = T checks`, where `LABEL` is nonempty, `P`, `F`, and `T` are unsigned
  decimal integers, and `T=P+F`. The first earlier stdout line beginning `FAIL` followed by one or more
  spaces is the first failed assertion. `F=0` requires no such earlier line and emits `-`; `F>0`
  requires at least one. With no matching tally, both tally and first-failure fields are inapplicable;
  malformed or multiple matching tallies are evidence errors.
- `ledger-result-v1` recognizes at most one UTF-8 stdout line beginning `RESULT|`, followed by unique
  nonempty `name=value` fields separated by `|`, including one `status` field. A field name matches
  `[A-Za-z_][A-Za-z0-9_]*`; its value is nonempty and contains no `|`; the first `=` separates name from
  value. It also recognizes each line `MISMATCH NAME:` followed by any text, retaining the nonempty
  bytes between `MISMATCH ` and the first `:` as `NAME`, in encounter order. No `RESULT|` line yields
  status `NO_RESULT_LINE`; malformed/multiple result lines are evidence errors.

An artifact `observation_format` is `none` or `ledger-dimension-v1`.
`ledger-dimension-v1` recognizes unique records of the form
`DIM|axes=AXES|name=NAME|exponents=VALUE`, where all three values are nonempty and contain no `|`. For
a row it selects the matching record whose `NAME` equals that row's `record`, compares its `VALUE` with
the run-entry artifact's record, and emits the structured moved-value observation defined in C-4.
Missing post-producer artifacts emit no moved-value observation. Missing/duplicate selected records in
an otherwise emitted artifact are evidence errors.

`stable_inputs` is retained solely for R5. It declares the complete project-tree read set of producer
and checker execution apart from the selected target, declared artifacts, copied configuration, and
copied include-list. The driver file is added automatically. A directory expands recursively to its
regular files in bytewise path order. The driver authenticates their entry path/digest set in recovery
state and requires the same set and digests before and after every attempted row and on resume.

Network access, clock-dependent decisions, unseeded randomness, and reads from mutable undeclared
inputs are nonconforming. External executables, interpreters, libraries, the OS, and hardware are not
made part of this set; R5 therefore also presumes that the external toolchain's behavior does not
change across the uninterrupted/resumed comparison. For equal target/artifact entry state, stable
inputs, configuration, include-list, environment, and toolchain behavior, a conforming producer and
checker must produce the same exits, captures, and artifacts. This is a command/configuration
precondition: the driver cannot prove determinism of an arbitrary child.

*Choice and tradeoff — direct argv, a closed environment, ordered checkers, declared reset state, and
fixed report grammars make R3, R5, and extraction decidable over declared state; they foreclose shell
command strings, ambient environment dependencies, undeclared mutable state, and arbitrary prose
parsing. The stable-input rule is not a promise that the evidence is a replay capsule.*

## C-3. Include-list wire format

The include-list is TSV and its header contains each of these required names exactly once:

```text
axis	key	record	old_text	new_text
```

Its data rows are the complete mutation set. The driver must neither scan a target for candidate sites
nor synthesize, omit, reorder, or add a row.

The five required columns may appear in any order/header positions so that the committed stage023 list is
accepted losslessly. The only optional driver-defined columns are `target` and `line`; all other unique
header names are opaque extra columns. Header names are nonempty and unique.

Every data row has exactly as many fields as the header. Quoting and backslash escaping do not exist;
tabs and line breaks therefore cannot occur inside a value. `axis` and `key` match
`[A-Za-z0-9_.-]+`. `record` and `old_text` are nonempty. `new_text` may be empty, permitting deletion.
`old_text` and `new_text` must differ. `line`, when present, is either empty or a positive unsigned
decimal line number.

If the header has no `target` column, the config must have exactly one target and every row resolves to
it. If the column exists, every row's value is a configured target name; an empty target cell is
invalid. A multi-target run therefore requires the column.

Row identity is the resolved target name, `axis`, and `key`; that triple must be unique. `record` is an
observation selector/metadata value, not identity. Extra columns are never discarded: their
header/value pairs are carried into the result row's `include_extra` object, in input header order
before JCS canonicalization.

Applicability is measured against the selected target's run-entry text:

- with an empty/absent `line`, `match_count` is the number of non-overlapping exact occurrences of
  `old_text` in the whole target;
- with `line=N`, it is the number of non-overlapping exact occurrences wholly within logical line `N`.

Every target is a regular UTF-8 file without a BOM. Logical lines are LF-delimited, numbered from one,
with an unterminated final line allowed. Mutation replaces the UTF-8 bytes of the sole counted
occurrence and no other bytes. A row is **unusable for execution** when `match_count` is not exactly
one. That is not a list-validation failure: `run` must bank the zero/multiple application outcome and
must invoke neither producer nor checker for it.

Syntactically blank lines, short/long rows, duplicate headers/identities, an unknown target, invalid
UTF-8, invalid line number, identical old/new text, or a missing required value are list errors.

*Choice and tradeoff — named required columns accept the real five-column list while exact field
preservation defeats accept-and-ignore parsing; the TSV choice forecloses embedded tabs/newlines and
defines list whitespace as data rather than formatting.*

## C-4. Result and capture wire formats

### `results.tsv`

The result file has exactly this tab-separated header and no extra columns:

```text
schema	list_sha256	row_number	row_id	target	axis	key	record	line	old_text	new_text	include_extra	match_count	applied	target_entry_sha256	target_mutant_sha256	producer_exit	producer_tally	first_failed_assertion	artifacts_emitted	moved_values	checker_exits	checker_verdicts	captures	outcome
```

It has one LF-terminated data row per banked include row, in include-list order. Empty cells are
invalid. The sole whole-cell inapplicability sentinel is the single byte `-`. It is never accepted as a
synonym for an empty array/object/string.

Field grammar is:

| Field | Exact representation |
|---|---|
| `schema` | Literal `ablation-result-v1`. |
| `list_sha256` | SHA-256 of the exact include-list bytes. |
| `row_number` | One-based unsigned decimal include-list position, without leading zeroes. |
| `row_id` | `target:axis:key`; colons are impossible in its components. |
| `target`, `axis`, `key` | Resolved token values. |
| `record`, `old_text`, `new_text` | JCS JSON strings, preserving the list values exactly. |
| `line` | Positive unsigned decimal, or `-` for an empty/absent list cell. |
| `include_extra` | JCS object of every non-required, non-`target`, non-`line` include column with every value encoded as a JSON string; `{}` when none. |
| `match_count` | Unsigned decimal integer without leading zeroes. |
| `applied` | Exactly `true` or `false`; no `yes`, `no`, `1`, or `0`. |
| `target_entry_sha256` | SHA-256 of the exact bytes in which matches were counted. |
| `target_mutant_sha256` | SHA-256 of the exactly-once mutant, or `-` when not applied. |
| `producer_exit` | Integer `0..255`, or `-` when not applied. |
| `producer_tally` | JCS object `{"fail":F,"pass":P}`, or `-` when not reported/not applied. |
| `first_failed_assertion` | JCS JSON string containing the complete recognized line, or `-`. |
| `artifacts_emitted` | JCS array of emitted artifact names in config order; `[]` is meaningful. |
| `moved_values` | JCS array of moved-value objects in artifact-config then artifact-record order; `[]` is meaningful. |
| `checker_exits` | JCS array of `{"exit":N,"name":"NAME"}` in config order. |
| `checker_verdicts` | JCS array of verdict objects in config order. |
| `captures` | JCS array of capture-reference objects in role order. |
| `outcome` | Exactly one C-5 token. |

A moved-value object is exactly:

```json
{"after":"{13, -7, 5}","artifact":"sidecar","before":"{1, 0, 0}","moved":true,"name":"sourced_dims.a"}
```

`moved` is `true` iff `before` and `after` differ byte-for-byte, otherwise `false`. Absence of an
artifact is represented by absence of its name from `artifacts_emitted` and by no moved-value object;
it is never called an unmoved value.

An `exit-only` checker verdict object is exactly
`{"mismatch_names":[],"name":"NAME","status":"EXIT_ZERO"}` or the same with `EXIT_NONZERO`.
A `ledger-result-v1` verdict object is exactly
`{"mismatch_names":["NAME",...],"name":"NAME","status":"STATUS"}`. Arrays retain encounter order and
preserve duplicates; an empty array is `[]`.

A capture reference is exactly
`{"bytes":N,"path":"PATH","role":"ROLE","sha256":"DIGEST"}`. Applied rows contain, in this order:
producer stdout/stderr; each checker's stdout/stderr in checker order; then every emitted artifact in
artifact order. Empty stdout/stderr are still retained zero-byte captures. Roles are
`producer.stdout`, `producer.stderr`, `checker.NAME.stdout`, `checker.NAME.stderr`, and
`artifact.NAME`.

The fixed banked capture paths are:

```text
captures/ROW_NUMBER/producer.stdout
captures/ROW_NUMBER/producer.stderr
captures/ROW_NUMBER/checker-NN-NAME.stdout
captures/ROW_NUMBER/checker-NN-NAME.stderr
captures/ROW_NUMBER/artifact-NAME
```

`ROW_NUMBER` is six-digit zero-padded and `NN` is the two-digit one-based checker position. A run with
more than 999999 rows or 99 checkers is invalid. Captures from an interrupted, unbanked attempt must be
retained separately and must never replace these banked captures; they are not referenced by a result
row.

The other files a completed run leaves for commit, and the claims made for them, are fixed in C-6.

*Choice and tradeoff — a closed header, one sentinel, canonical JSON collections, and raw capture
references make row claims checkable without creating a replay package; they foreclose ad hoc result
columns, blank-as-meaning, alternate boolean spellings, and committing only a summarized TSV.*

## C-5. Complete outcome truth table

First apply these row invariants. A row violating any one is a schema violation and has **no valid
outcome**:

1. Identity/list fields exactly reproduce the resolved include row and `list_sha256`.
2. `applied` is `true` iff `match_count=1`.
3. If `applied=false`, `target_mutant_sha256`, `producer_exit`, `producer_tally`,
   `first_failed_assertion`, `artifacts_emitted`, `moved_values`, `checker_exits`,
   `checker_verdicts`, and `captures` are all `-`.
4. If `applied=true`, `target_mutant_sha256` is a digest different from
   `target_entry_sha256`; `producer_exit` is an integer; `artifacts_emitted` and `moved_values` are
   arrays; `checker_exits`/`checker_verdicts` contain exactly the configured nonempty checker sequence;
   and `captures` contains exactly the required stream and emitted-artifact references with matching
   bytes/digests.
5. `producer_tally` is `-` or is derived from the producer capture under its report grammar;
   `first_failed_assertion` is `-` or is the capture's derived first failure. Every artifact and checker
   observation is likewise derivable from its own capture and configured grammar.

For any row satisfying those invariants, exactly one predicate below is true, and `outcome` must be its
token:

| `outcome` | Complete predicate over that row |
|---|---|
| `NOT_APPLIED_ZERO` | `match_count=0` and `applied=false`. |
| `NOT_APPLIED_MULTIPLE` | `match_count>=2` and `applied=false`. |
| `PRODUCER_NONZERO` | `match_count=1`, `applied=true`, and `producer_exit` is in `1..255`; checker exits do not alter this predicate. |
| `CHECKER_NONZERO` | `match_count=1`, `applied=true`, `producer_exit=0`, and at least one member of `checker_exits` has `exit` in `1..255`. |
| `ALL_ZERO` | `match_count=1`, `applied=true`, `producer_exit=0`, and every member of the nonempty `checker_exits` has `exit=0`. |

The predicates are exhaustive because `match_count` is an unsigned integer, exactly-one application
has one integer producer exit, and the checker list is nonempty with integer exits. They are disjoint
by construction. In particular:

- `PRODUCER_NONZERO` claims only a producer exit, not that the mutation was scientifically “caught”;
- `CHECKER_NONZERO` claims only checker exit behavior, not what any checker read; and
- `ALL_ZERO` claims only zero process exits, not that a mutant “survived” or is semantically acceptable.

Unknown tokens, `PRODUCER_NONZERO` with producer exit zero, `CHECKER_NONZERO` with all exits zero,
`ALL_ZERO` with any nonzero exit, or either application token with execution fields are schema
violations. The driver must never bank one, and an A9 consumer applying this contract must reject one;
no separate driver verb is required to restate the table.

*Choice and tradeoff — exit-based, mutually exclusive predicates make R9 mechanical without assuming
checker inputs; they foreclose stronger causal words such as caught, killed, survived, or detected.*

## C-6. Commit set, evidence claims, and named prerequisites

A successfully completed `run` leaves this commit-candidate set in `evidence_path`:

```text
EVIDENCE.md
config.json
include_list.tsv
results.tsv
restore.tsv
captures/...
COMPLETE
```

`config.json` and `include_list.tsv` are byte copies of the accepted inputs. `COMPLETE` is exactly
`ablation-evidence-v1` plus LF and is created only after every list row is banked and the final
successful restoration event is durable. It is R5 completion state, not a package-integrity manifest.

Interrupted-attempt captures may additionally remain under `attempts/`; they are never result-row
evidence and need not be committed. `recovery_path`, reset caches, recovery snapshots, external
executables/libraries, stable-input digest inventory, and host state are not in the commit set. There
is no `SHA256SUMS`, `environment.json`, target/artifact entry-snapshot directory, caller-snapshot
package, dependency archive, or replay capsule.

`EVIDENCE.md` is UTF-8 Markdown with this fixed content, except that the final list has one item for
each distinct producer/checker `argv[0]` in first-use order. Each item is that value encoded as a JCS
JSON string:

```text
# Ablation evidence scope

This directory records one ablation run's observations; it is not a self-contained replay package.
The claims are the fields in `results.tsv`, the bytes and digests referenced under `captures/`, and the restoration events in `restore.tsv`.
The accepted run configuration and include-list are included as `config.json` and `include_list.tsv`.
Audit also requires the repository tree at the commit containing this directory: each target and any pre-existing artifact must match the final restored state.
Reproducing a row additionally requires that repository tree, compatible external toolchain behavior, POSIX filesystem/process/signal semantics, and the child environment recorded in `config.json`; stable-input entry digests and external executable, library, OS, kernel, or hardware bytes are not captured here.
Accordingly this evidence does not claim that the commit alone can re-run a row or that a new run will repeat an outcome.

Command entry points required for reproduction:
- "python3"
```

The last heading is followed by at least one item because the producer and checker lists are nonempty;
duplicate command names appear once. A committed summary may derive claims from this set, but must not
say that an omitted capture, recovery snapshot, toolchain, or run log is committed or independently
auditable.

R8 holds for a committed run only when the whole commit-candidate set is committed and the repository
at that commit passes the target/artifact checks stated in `EVIDENCE.md`. Raw captures support the
extracted textual fields and emitted-artifact observations; C-5 limits outcome claims to the row's
recorded process exits. A later target edit, a missing referenced capture, or prose that claims more
than these files contain is an R8/A8 failure, not a reason to add replay infrastructure.

The pre-driver stage023 package remains the A7 oracle under C-9. It is not retroactively judged by this
run-output rule, and its missing captures must still be described honestly.

*Choice and tradeoff — the run leaves observations and an explicit scope statement, while naming the
repo, configuration, toolchain behavior, and POSIX semantics a reader additionally needs. This answers
v4's evidence question without promising replay from a commit.*

## C-7. Digest and restoration evidence

`restore.tsv` is in the evidence root and has exactly this header:

```text
event	operation	kind	name	entry_state	entry_sha256	pre_restore_state	pre_restore_sha256	restored_state	restored_sha256	restored
```

`event` is a one-based unsigned decimal restoration-event number. `operation` is `run`, `signal`,
`error`, or `repair`. `kind` is `target` or `artifact`; `name` is the configured name.
`entry_state`, `pre_restore_state`, and `restored_state` are `file` or `absent`. A state of `file` has
the corresponding SHA-256; `absent` has `-`. For each event, there is exactly one row for every
configured target followed by every configured artifact, in config order, sharing the event and
operation. Targets always have `entry_state=file`. Entry state/digest name the original run-entry
state. `restored` is exactly `true` iff restored state and digest equal entry state and digest,
otherwise `false`.

An operation may report restoration success only if its final event has the exact configured
target-and-artifact set and every row has `restored=true`. Proof is never collapsed into one aggregate
digest. A multi-target run therefore has one proof row per target, including configured targets not
selected by a banked row, plus its artifact rows.

The recovery state must exist and authenticate every restoration-set snapshot before any target is
mutated. Restoration copies those run-entry bytes; it must never derive them from Git and must never
invoke `git checkout`, `git restore`, or `git stash`.

*Choice and tradeoff — per-event, per-path SHA-256/state proof exposes partial multi-target or artifact
restoration and supports repair; it forecloses a single aggregate “restored” assertion and any
Git-derived rollback.*

## C-8. Row identity, banking, list sameness, and equality

A `row_id` is exactly `resolved-target:axis:key`. The include-list uniqueness rule makes it
unambiguous.

A row is **banked** only when all of the following are simultaneously true:

- its complete LF-terminated result row is present in `results.tsv` at its include-list position;
- every referenced capture exists at its final path and matches its recorded size/digest;
- the row passes all C-4 invariants and exactly one C-5 predicate; and
- all earlier include rows are banked.

Thus banked rows are a contiguous prefix. A partial result line, captures without a valid row, a row
without all captures, or an in-flight child is unbanked. Resume skips the prefix and starts with the
first unbanked include row. Captures belonging to any banked or interrupted attempt are never
overwritten; an unfinished attempt may be retained outside the fixed banked paths, but may not be
mistaken for banked evidence.

“Present” above means the LF-terminated row and referenced capture bytes are visible to an independent
reader after the writer has flushed them; data held only in process memory is not banked.

Two include-lists are the **same list** iff their exact file bytes have the same SHA-256 and are
byte-for-byte equal. Semantic TSV equivalence, newline normalization, column reordering, or whitespace
normalization is not sameness. Resume additionally requires byte-equal accepted config after removing
no fields—`evidence_path` and `recovery_path` included—matching entry snapshot digests, and the exact
C-2 stable-input path/digest set recorded at run entry.

Two result files are **equal** iff their complete bytes are equal, including header, row order, JSON
serialization, sentinels, and final LF. Result rows intentionally contain no timestamps, process IDs,
attempt numbers, absolute paths, or run-generated random identifiers. This is the equality R5 requires
between uninterrupted and resumed runs.

*Choice and tradeoff — content identity, prefix banking, and byte equality make resume corruption and
silent list drift visible; they foreclose “equivalent after normalization” resumes and nondeterministic
metadata in results.*

## C-9. Stage023 legacy mapping for A7

The A7 comparison joins the legacy `results.tsv`, the committed legacy `include_list.tsv`, and the new
result on `(axis,key)`. Both sides must have exactly 51 joined rows, exactly 22
`A1_DECLARATION` and 29 `A2_BINDING`, unique joined keys, and equal `record`. New mutation text must
equal the corresponding legacy include-list `old_text`/`new_text`. The A7 config names the sidecar
artifact exactly `stage023_sidecar` with `ledger-dimension-v1`, uses `ledger-tally-v1` for the producer,
and has exactly one checker named `stage023_comparator` with `ledger-result-v1`.

The complete mapped projection is:

| Legacy field/value | New field/value |
|---|---|
| `axis` | `axis`, exact string |
| `key` | `key`, exact string |
| `record` | decoded `record` JSON string, exact string |
| `stage_exit` | `producer_exit`, decimal integer |
| `pass_count`, `fail_count` | decoded `producer_tally.pass`, `.fail` |
| blank `first_fail` | `first_failed_assertion=-` |
| nonblank `first_fail` | decoded `first_failed_assertion` JSON string, exact string |
| `sidecar_written=yes` | `stage023_sidecar` is present in `artifacts_emitted` |
| `sidecar_written=no` | `stage023_sidecar` is absent from `artifacts_emitted` |
| `record_moved=no_sidecar` and blank `emitted_value` | no `moved_values` member for `stage023_sidecar` |
| `record_moved=yes` | exactly one matching moved-value member with `moved=true` |
| `record_moved=no` | exactly one matching moved-value member with `moved=false` |
| nonblank `emitted_value` | that member's `after`, exact string |
| `cmp_exit` | exit of `stage023_comparator` |
| `cmp_status` | that checker's decoded `checker_verdicts.status`, exact string |
| `mismatch_names=none` | that checker's `mismatch_names=[]` |
| any other current `mismatch_names` cell | a one-element `mismatch_names` array containing that exact cell |

For this import/comparison only, legacy `yes` and `no` normalize to JSON booleans `true` and `false`
under the two membership tests above; the new result grammar never accepts `yes`/`no`. A legacy blank
normalizes only where the table says so; no other blank is legal. `none` is a legacy mismatch-list
sentinel, not the new schema's `-`. The table is exhaustive: an unlisted legacy spelling or an
inconsistent combination such as `no_sidecar` plus a nonblank emitted value is a comparison error.

The new outcome is compared with the outcome derived from the normalized legacy exits:

- `stage_exit != 0` derives `PRODUCER_NONZERO`;
- `stage_exit = 0` and `cmp_exit != 0` derives `CHECKER_NONZERO`;
- `stage_exit = 0` and `cmp_exit = 0` derives `ALL_ZERO`.

Applied to the committed legacy bytes, this derivation yields 16 `PRODUCER_NONZERO`, 35
`CHECKER_NONZERO`, and zero `ALL_ZERO` rows; the same file has 35 `sidecar_written=yes` /
`record_moved=yes` rows and 16 `sidecar_written=no` / `record_moved=no_sidecar` rows.

The new row must also have `match_count=1` and `applied=true`; the legacy result has no fields that can
establish application. New-only provenance, digest, capture, list-digest, and schema fields have no
legacy counterpart and are not fabricated.

Therefore A7 is an **exact mapped-field comparison**, not byte equality and not exact reproduction of
the legacy 13-column file. ✅ **DONE, 2026-07-29 — the manifest correction this paragraph called for has
been made.** `manifests/DIMENSION_REWRITE.md` §12b(b) now states A7 as reproducing the committed
`results.tsv` — 22 `A1_DECLARATION` and 29 `A2_BINDING` rows — **“exactly on every mapped field”**, adds
“⚠ **Not byte equality, which was never achievable**”, and makes the legacy-column-to-new-field mapping
part of the driver's contract, frozen with it rather than chosen at comparison time. The committed legacy
file itself must not be edited, and was not.

*Choice and tradeoff — an explicit total mapping preserves every legacy observation without pretending
that absent provenance exists; it forecloses a byte-for-byte “exactly” claim, which the manifest no
longer makes.*
