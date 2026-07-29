# Ablation driver — requirements and acceptance intent (v4)

**Status: USER DECISION, 2026-07-29** (`manifests/DIMENSION_REWRITE.md` §12b(b)). Leads the build queue;
existing hand-rolled harnesses switch over to it; every subsequent stage reuses it. It is built before
stage027 begins, and 027 is the first conversion to use it — that is the build queue, not the conversion
order (§8).

⚠ **v4 — USER DECISION 2026-07-29: scope down, R8 is dropped.** The contract round showed the scope had
outgrown the "~100-line driver" the manifest records: a 618-line contract with five operations, replay
capsules and closed read sets. Tracing the cost, three of its twelve findings descended from **R8**
("committed evidence alone lets a reader re-run a row"), which the requirements author had added in
response to an evidence summary whose claims were true of its run directory and false of its commit.
⭐ That was a **prose** defect with a **prose** fix; converting it into an infrastructure requirement is
what pulled in transitive read sets, replay capsules and platform prerequisites. R8 is therefore replaced
by the honest-scope rule below, and the retrofit (A7) remains the real gate.

⚠ **v3 was a division-of-labour correction, not a third attempt at the same document.** v1 and v2 both
failed design review for one underlying reason: they tried to fix an **interface** — argv, wire format,
sentinels, outcome predicates, column mappings — and that is design. Design of deliverable code is
Codex's; `docs/development_pipeline.md` (Roles) says the orchestrator "never pre-designs the
computational route." Checklist item 1b pushed that way, because fixtures need something to bind to. The
resolution is to **sequence**, not to design:

1. **Codex authors the CONTRACT** — everything in §C — as a standalone deliverable, with no
   implementation.
2. **A different session freezes the fixtures** against that contract.
3. **A third session implements** the driver, and may neither author nor weaken those fixtures.

Fixtures still precede the implementation and are still independently authored, which is what 1b
protects. Prior findings: `codex_verify_persisted.log` (v1), `codex_review_v2.log` (v2),
`FIXTURES_REPORT.md` §"specification findings".

**What it is.** A per-tooth ablation runner: it takes an orchestrator-supplied list of mutations, applies
them one at a time to the real file at its real path, runs a producer and one or more checkers after
each, and banks one result row plus retained captures per mutation.

**What it is not.** It never decides *what* to mutate. Tools whose `init` enumerates the sites were
rejected partly for that: the target list belongs to the orchestrator (`development_pipeline.md`
checklist item 6).

---

## C — what the CONTRACT must settle

⛔ This is a list of questions, not answers. The contract's content is Codex's to author; what is
required here is only that each of these has **one** answer, written down, before fixtures are frozen —
because each was found, by review, to be a place where two parties would otherwise invent different
things and the fixtures would end up binding to the implementation they grade.

1. **Operations and invocation.** What operations exist (at minimum: validating a list, running one, and
   repairing a tree after an uncatchable kill), their exact argv, and their exit-code meanings.
2. **Run configuration.** How the target, producer command, checker commands, declared artifacts, working
   directory, result path and capture path are supplied; how repeated checkers are ordered and encoded;
   shell-vs-argv execution and environment.
3. **The include-list wire format.** Required and optional columns, uniqueness, where the target is named
   when a list has no target column, what "unusable" means, and what happens to columns the driver does
   not use.
4. **The result wire format.** Header and column order, whether extra columns are allowed, boolean
   spellings, the sentinel for inapplicable, and how multi-valued fields — several checkers, several
   artifacts, several moved values — are serialized.
5. **The outcome truth table.** Not a vocabulary: a **complete predicate per outcome**, over the row's own
   fields, such that exactly one outcome is derivable from any row and an unsupported claim is a schema
   violation rather than a judgement. ⚠ Do not assume every checker reads an emitted artifact; some
   inspect the target, stdout, or another input.
6. **What the committed evidence claims, and what it names as prerequisites.** Settle which files a run
   is expected to leave behind for commit, and what a reader is told they additionally need in order to
   reproduce a row (the repo at that commit, the toolchain, the run configuration). ⛔ The contract does
   **not** have to make a row re-runnable from the commit alone — that requirement was withdrawn. It has
   to make the evidence state its own limits, so no claim in it is true of the run directory and false of
   the commit.
7. **Digest evidence.** Where restore proof is emitted, its algorithm and format, and what it means when a
   run has more than one target.
8. **Row identity and equality.** What makes a row "banked", what makes two lists "the same list", and
   what equality of two result files means.
9. **The legacy mapping for A7.** The committed `notes/stage023_step_h_evidence/results.tsv` uses
   `stage_exit`, `pass_count`/`fail_count`, `first_fail`, `sidecar_written` and has no outcome column. The
   contract must state the mapping to the new schema, and how legacy blanks and `yes`/`no` spellings
   normalize. ✅ **DONE, 2026-07-29:** the contract states it as a mapped comparison (§C-9), and
   `manifests/DIMENSION_REWRITE.md` §12b(b) has been corrected to match — it now says the retrofit
   reproduces that file **"exactly on every mapped field"**, explicitly **not** byte equality, with the
   mapping frozen as part of the driver's contract.

---

## R — requirements (properties the driver must have)

**R1 — the list is an input.** The driver never scans for candidate sites and never adds a row.

**R2 — the mutation is visible at the real file's real path while the producer runs.** This workstream's
controls hash the file's own bytes and derive artifact paths from `__file__`; a mutant living anywhere
else makes those controls pass vacuously or fail uniformly.

**R3 — no state from a previous mutant is observable by the next** — no artifact it emitted, no stale
bytecode.

**R4 — captures are per mutant and are never overwritten.**

**R5 — resume re-runs no already-banked row, completes every unbanked row, and yields the same result
file as an uninterrupted run.**

**R6 — a mutation that does not apply exactly once is a recorded outcome**, distinguishing zero from
multiple matches, never banked as a mutant that ran.

**R7 — restore on every path the driver can observe** (normal end, error, catchable signal), with digest
proof; plus a separately invocable repair for the uncatchable class. ⛔ Never by `git checkout`,
`restore` or `stash` — this workstream ablates uncommitted work and those commands destroy it.

**R8 — the evidence claims exactly what it supports.** Every claim in what the run leaves behind is true
of the committed files, not merely of the run directory; anything a reader would additionally need in
order to reproduce a row is named. ⚠ This replaces the withdrawn re-runnable-from-the-commit property:
the defect it was written against was a claim that outlived its evidence, not an inability to replay.

**R9 — no row's outcome asserts more than its own fields support.**

---

## A — acceptance intent

Each item states **what must be demonstrated by execution**. The mechanics of each demonstration are
settled by §C's contract, and each must be able to fail: for every item, the fixture author states which
implementation deficiency it catches, and shows it catching it.

**A1** — R3 is load-bearing: a residual artifact that would otherwise change a mutant's outcome does not.
**A2** — R6's zero-match and multi-match floors fire, distinguishably.
**A3** — resume is exact under a named interruption boundary.
**A4** — restore holds after a clean run, after that interruption, and after an uncatchable kill via R7's
repair operation.
**A5** — every mutant's captures survive the run and correspond to that mutant.
**A6** — the real committed include-list validates, with a floor that an accept-and-ignore parser fails:
all 51 rows, the two axis counts, and the exact mutation text of named rows reproduced.
**A7** — the retrofit reproduces the committed stage023 tables under §C-9's mapping, exactly on every
mapped field. ⛔ A disagreement is a finding to report, never a difference to adjust either side into.
**A8** — R8 holds: every claim in the run's committed output is checkable against the committed files
alone, and each stated prerequisite is named rather than assumed. A claim that requires a file the run
does not leave behind fails this.
**A9** — R9 holds: a row whose fields do not support its outcome is rejected by the contract's own truth
table.

⚠ **Fixture-suite quality is graded separately and is not driver acceptance.** The requirement that every
negative fixture be shown failing against a deficient variant grades the *suite*; a driver cannot fail
it, and the fixture session correctly reported it satisfied before any driver existed. It stays a
requirement on the fixture deliverable, stated there.

---

## Evidential note — what motivates R3 and R8, and what survives

⚠ Both motivations were measured in this session's gitignored scratch tree. State them as motivation,
never as reproducible evidence:
- **R3** — in the pre-fix run of stage023's declaration axis, most checker verdicts were freshness
  failures decided by a residual artifact rather than comparisons. The committed `results.tsv` is the
  **corrected** run; the defective one is not in git.
- **R8** — several claims in this stage's hand-written evidence summary were true of the run directory
  and false of the commit, and were corrected under review. Only the corrected text is committed.

---

## Process constraints

⛔ **Acceptance tooling cannot grade itself** (`development_pipeline.md` checklist 1b) — hence the
three-session sequence above.
⛔ **No expected value or reason** in any directive or in any file a build session can read.
Scripts run under `timeout 600`; a 124 is a failure to report, never a reason to raise the cap.
