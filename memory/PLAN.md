# Research Memory Migration and Implementation Plan

Status: planning
Created: 2026-08-24
Target directory: `memory/`

## Purpose

Build a small, durable research memory that lets an agent answer:

- Where is the current treatment of a topic?
- What does a script or audit do, and what does it produce?
- Which conclusions are measured, derived, provisional, open, disputed, or
  superseded?
- Which older source was revised, and why is another source currently
  preferred?

The memory is a navigation and synthesis layer over the repository. It is not
a second formal model of the research and it is not an authority independent
of the cited source files.

The useful content in `atlas/` and `graph/` will be migrated before those
directories are deleted. Deletion is a separate, gated cutover step.

## Governing Principles

1. Repository sources remain authoritative. Memory pages summarize and point
   to them.
2. The memory is intentionally incomplete. It tracks selected canonical
   sources, not every file matching an extension.
   The confirmed normal discovery roots are `research/` and `software/`, and
   discovery is limited to Git-tracked content. Root `notes/` and `docs/` are
   excluded.
3. Every substantive statement cites a repository-relative source path and,
   when available, a stable TeX label or Markdown heading.
4. AI-generated prose never upgrades the evidentiary status of a claim.
5. Newer does not automatically mean better. Explicit revision relationships,
   applicable evidence, reviewed source adjudications, and configured lineage
   outrank dates. Review of memory wording is a separate property.
6. A detected contradiction is recorded as a candidate conflict until the
   sources explicitly resolve it or it is reviewed.
7. Old information is retained as superseded or historical when it helps
   explain the research lineage.
8. Deterministic code owns discovery, hashes, state transitions, and linting.
   AI owns summaries, topic linking, and candidate-conflict detection.
9. Synchronization state must be rebuildable from the repository and memory
   pages. It is an operational cache, not research truth.
10. Prefer one useful representation of a fact. Do not reproduce the old
    graph as claim pages, edges, canvases, dashboards, and exports.

## Planned Layout

Only this plan is created initially. The remaining structure is created after
Phase 0 decisions are recorded.

```text
memory/
├── PLAN.md
├── README.md                 # how agents and humans use the memory
├── index.md                  # small entry point; read first
├── conflicts.md              # open and recently resolved conflicts
├── topics/                   # synthesis across selected sources
├── sources/                  # replace-on-update capsule for each source
├── scripts/                  # script catalogs grouped by research domain
├── prompts/                  # versioned semantic extraction instructions
├── tools/                    # deterministic scanner/status/lint helpers
└── _meta/
    ├── schema.md             # authoritative page and status rules
    ├── config.yaml           # source allowlist, exclusions, classifications
    ├── state.json            # generated hashes and sync state
    └── atlas-migration.yaml  # disposition of legacy atlas/graph content
```

There will be no `raw/` directory. Existing repository paths are cited
directly.

There will be no general `claims/` directory in the first version. Important
claims live in source capsules and topic pages. Standalone claim records may be
introduced later only if actual use demonstrates a need for them.

## Page Model

### Source capsules

Each selected source unit gets one capsule that is regenerated when that unit
changes. A source unit may be one file or a declared bundle, which is required
for componentized papers, step records with paired engines, and software result
families. Every bundle member is hashed independently. A capsule contains:

- Source path, content hash, source kind, and relevant declared date.
- Short purpose and scope.
- Key statements with exact source anchors where possible.
- Important assumptions, limitations, and open questions.
- Measurement or script provenance stated by the source.
- Related topic IDs.
- Explicit revision or supersession relationships found in the source.

A capsule reports what its source says. It does not silently reconcile that
source with the rest of the repository.

### Topic pages

Topic pages synthesize several source capsules. They contain:

- A concise current working position.
- Statement evidence and memory-review qualification.
- The preferred sources and the explicit reason they are preferred.
- Superseded positions that remain relevant.
- Open questions and candidate conflicts.
- Links to relevant scripts and source capsules.

Topic pages are the main answer surface for future agents. The number of topic
pages should remain small enough that `index.md` is genuinely useful.

### Script catalogs

Scripts are grouped by domain rather than placed in one repository-wide table.
Only meaningful entry points, validators, builders, and audit tools are
cataloged by default. Helpers, caches, generated code, scratch probes, and
vendored material are excluded unless a selected entry point depends on them
in a way an agent needs to understand.

Each catalog entry contains:

- Repository path and language.
- What object it computes, builds, validates, or compares.
- Inputs and literal output artifacts.
- Guards or invariants it checks.
- Related source and topic pages.
- The source/commit against which the description was last refreshed.

The description must not turn a script's computed output into an unsupported
conclusion.

### Status dimensions

Avoid one overloaded `status` property. Pages and statements use separate
dimensions:

```yaml
lifecycle: current | superseded | deleted | retired
memory_review: ai_draft | reviewed
```

Evidence is statement-level because one page commonly mixes measured,
derived, provisional, open, and disputed material.

These dimensions answer different questions:

- `lifecycle` describes place in the document lineage.
- Statement-level `evidence` describes support represented by the cited source.
- `memory_review` describes only whether the memory wording itself was checked.

The schema will define when each value may be used. The AI may default to
`provisional` and `ai_draft`; it may not infer `measured` or `reviewed` merely
from confident prose.

## Source Precedence and Conflict Policy

The current working position is selected in this order:

1. An explicit `revises`, `supersedes`, correction, or retraction relationship
   in the sources.
2. A recorded reviewed adjudication comparing the same object and scope.
3. The canonical source designated in `_meta/config.yaml` for that document
   lineage.
4. Recency, but only among sources of comparable role, evidence, source
   adjudication, scope, and lineage.
5. If none of the above resolves the difference, retain both positions and
   record a candidate conflict.

AI-detected conflicts must include:

- The two source paths and exact anchors when possible.
- The statements that appear incompatible, paraphrased conservatively.
- Whether the conflict may instead be a scope, convention, regime, or notation
  difference.
- The precedence rule, if any, that resolves it.
- Review state and resolution history.

## Phase 0: Scope and Semantic Contract

Phase 0 prevents uncontrolled ingestion and establishes what the memory is
allowed to say.

### 0.1 Inventory source classes

Create a deterministic inventory of candidate sources, grouped into at least:

- Canonical papers and maintained TeX sources.
- Current audit records, result reports, and handoff documents inside the selected corpus.
- Measurement records and literal engine output.
- Working notes and exploratory material.
- Scripts and executable audit entry points.
- Prompts, directives, reviews, generated reports, scratch files, and archives.

The inventory records counts and representative paths; it does not ingest the
files.

### 0.2 Establish an explicit corpus policy

Write `_meta/config.yaml` with:

- Included roots and explicitly selected files.
- Excluded roots and path patterns.
- Source-kind classifications.
- Canonical document lineages.
- Script entry-point rules.
- Maximum file-size or chunking policy.
- Treatment of generated, archived, scratch, and review artifacts.

The confirmed candidate roots are `research/` and `software/`. Discovery uses
the selected committed Git tree and reads its blob objects; mutable index and
working-tree bytes, ignored files, and untracked material are not eligible for
normal ingest. Root `notes/` and `docs/` are excluded. Tracked process history
such as directives, prompts, and red-team material is citation-only by default.

The initial corpus should be a small representative set, approximately 20–50
sources, rather than the whole repository. It must cover at least one current
paper, one changing audit/software family, one conflict or supersession,
one status firewall, and representative Python and Wolfram Language scripts.

### 0.3 Define the evidence contract

Write the precise rules for `measured`, `derived`, `provisional`, `open`, and
`disputed`, consistent with the repository's measurement discipline in
`CLAUDE.md`.

In particular:

- A memory page cannot create evidence absent from the source.
- A measured statement must point to the command/output record that supports
  it when the source provides one.
- A script description states what is computed and emitted; interpretation is
  attributed to its step record or paper.
- Missing evidence is recorded as missing, not repaired by AI reasoning.

### 0.4 Define stable citations

Choose citation rules in this preference order:

1. Repository path plus TeX `\label{...}`.
2. Repository path plus Markdown heading.
3. Repository path plus a distinctive nearby identifier.
4. Line number only as a disposable convenience, never as the sole stable
   anchor.

### 0.5 Inventory the legacy atlas

Populate `_meta/atlas-migration.yaml` with one disposition for every valuable
legacy category and, where practical, every unique hand-maintained legacy
record:

```yaml
disposition: migrate | covered_by_source | generated_duplicate | obsolete
target: memory/topics/example.md
reason: concise explanation
```

The inventory must explicitly consider:

- Status firewalls and invalid-inference warnings.
- Open gates and unresolved research questions.
- Source crosswalks and canonical TeX anchors.
- Claim/theorem summaries.
- Physical and mathematical ontology summaries.
- Equation and derivation summaries.
- Script/query knowledge.
- Future-paper and citation-backfill queues.
- Query-validation examples.

Generated Obsidian pages, Canvas files, Bases, JSON mirrors, HTML/SVG/DOT
views, and generator metadata are expected to be classified as generated
duplicates unless they contain unique information.

### 0.6 Define a retrieval benchmark

Before implementation, write a small set of real questions the memory must
answer. It should include:

- Where is the latest treatment of a named active topic?
- What is currently open versus closed in that topic?
- Why must projection not be described as reduction?
- Which source revised an older conclusion?
- What does a representative Python script do?
- What does a representative Wolfram Language audit do and emit?
- Which sources support a selected paper-facing statement?

Record the expected source paths, not canned prose answers. This benchmark is
used during migration and before deleting the legacy system.

### Phase 0 completion gate

Phase 0 is complete only when:

- The initial corpus is explicit and bounded.
- Evidence and precedence rules are unambiguous.
- Every legacy content category has a migration policy.
- The retrieval benchmark exists.
- No source ingestion has occurred outside the approved corpus.

## Phase 1: Skeleton

Create the planned directories and initial files:

- `README.md`, `index.md`, and `conflicts.md`.
- `_meta/schema.md`, `_meta/config.yaml`, and empty generated state.
- Page templates or prompt contracts for source, topic, conflict, and script
  content.
- A short root-agent pointer to `memory/index.md`, without duplicating the
  schema in another instruction file.

Do not initialize the successful-sync baseline to `HEAD` before initial ingest
succeeds. An empty memory must remain visibly uninitialized.

## Phase 2: Deterministic Synchronization Core

Implement a small CLI with these operations:

```text
memory init
memory status
memory update [--paths ...]
memory lint
memory query <question>
```

One future agent skill may expose these operations; there should not be five
independently maintained skills.

The deterministic layer must:

- Discover only normalized configured source units.
- Enumerate normal inputs from Git-tracked paths under `research/` and
  `software/`; never recurse through the raw filesystem as the source list.
- Compare stored identities with committed target-tree blobs; inspect the
  working tree only to report dirty tracked members.
- Use Git diffs as an optimization, not the sole correctness mechanism.
- Detect committed changes, dirty tracked files, newly tracked files,
  deletions, and renames. Untracked files remain outside normal discovery.
- Track extractor/schema versions so semantic migrations can request
  reprocessing even when source content is unchanged.
- Materialize immutable committed inputs and an explicit worklist before
  invoking AI processing.
- Seal each semantic task with frozen prompts/schema, direct source inputs,
  prior IDs/history, hard output budgets, and staged dependency identities.
- Launch semantic writers in a Gitless isolated filesystem where the live
  repository is absent; mount only the task packet read-only and its declared
  output read-write, and require the isolation attestation at finalization.
- Stage and lint pages before journaled publication; update state last and
  recover or roll back interrupted publication.
- Preserve the prior successful state when a batch fails partway through.
- Never edit source research files.

Config also declares derived topic/script/conflict/index tasks, their exact
generated regions, input/reverse dependencies, direct-source coverage, order,
and output budgets. Capsules may nominate sources for synthesis, but every
load-bearing derived statement must receive committed direct-source bytes or a
deterministic excerpt in its sealed task.

`_meta/state.json` should record at least:

- Last fully processed commit.
- Per-source content hash and last processed commit.
- Source kind and lifecycle.
- Derived source capsule and related topic IDs.
- Extractor/schema version.
- Last processing result and error, if any.

The state format should be machine-owned and deterministic. AI output belongs
in Markdown pages, not in the synchronization logic.

## Phase 3: Initial Population and Atlas Migration

### 3.1 Populate source capsules

Process the approved Phase 0 corpus. Review representative capsules before
scaling the batch. A changed source replaces its capsule body from the current
source rather than repeatedly patching old prose.

### 3.2 Migrate valuable legacy knowledge

Migrate in this order:

1. Source crosswalks and source anchors needed by later pages.
2. Status firewalls and invalid-inference warnings.
3. Open gates and unresolved questions.
4. Selected high-value claim, equation, derivation, physical, and mathematical
   summaries.
5. Useful query examples and future-paper/citation queues.

Legacy material must be re-anchored to the original TeX, Markdown, measurement,
or script source. New memory pages must not rely on `atlas/` or `graph/` as
their final evidentiary citation.

Retain legacy node IDs in migration metadata where useful for auditability, but
do not reproduce graph edges or require those IDs for normal memory queries.

### 3.3 Build topic pages

Create only the topics needed by the initial corpus and retrieval benchmark.
Each topic page should synthesize source capsules, surface status boundaries,
and show unresolved conflicts without trying to encode a full dependency
graph.

### 3.4 Build script catalogs

Catalog representative entry points first. Expand by domain only when a script
is important to an active source or query. Record literal output locations and
validation targets rather than generic descriptions such as "does symbolic
math."

### 3.5 Refresh navigation

Build `index.md` from topic and source metadata, with a small curated overview
followed by deterministic lists of:

- Current major topics.
- Recently changed sources.
- Open conflicts.
- Open gates.
- Script domains.
- Memory freshness and pending changes.

## Phase 4: Incremental Update Workflow

For each update:

1. Load configuration and prior state.
2. Compute the scoped source delta from hashes and Git information.
3. Show the proposed worklist.
4. Regenerate capsules for changed or added sources.
5. Mark deleted sources and preserve relevant lineage.
6. Re-synthesize only topics linked to changed capsules.
7. Detect candidate conflicts against the affected topics.
8. Refresh script entries when selected scripts change.
9. Run lint and retrieval checks.
10. Atomically record the successful state.

`memory status` must warn a querying agent when relevant sources are pending or
the last update failed. `memory query` should read `index.md`, check freshness,
then load only the relevant topic, source, conflict, and script pages.

## Phase 5: Legacy Deletion and Cutover

Do not delete `atlas/` or `graph/` merely because the new files exist. Delete
them only after all of the following pass:

- Every legacy content category has a recorded disposition.
- Every item marked `migrate` has a valid target page.
- Migrated statements cite original repository sources rather than legacy
  generated pages.
- Status firewalls and open gates pass focused completeness review.
- The retrieval benchmark is answered from `memory/` without reading
  `atlas/` or `graph/`.
- `memory lint` reports no missing sources, broken internal links, orphaned
  derived pages, or stale state entries.
- Repository code and documentation have been checked for live dependencies on
  `atlas/` and `graph/`.
- A fresh agent can start at `memory/index.md` and locate the benchmark answers.
- The complete migration is committed so the legacy files remain recoverable
  from Git history.

Then perform the cutover as a separate reviewable change:

1. Remove live references and obsolete maintenance commands.
2. Delete `atlas/` and `graph/`.
3. Run the repository's relevant checks and `memory lint`.
4. Re-run the retrieval benchmark.
5. Commit the deletion separately from the migration.

## Phase 6: Controlled Expansion and Polish

After the cutover:

- Add sources by explicit corpus-policy changes.
- Add standalone claim pages only if repeated queries demonstrate that topic
  sections are insufficient.
- Add reminders or hooks only after the manual update workflow is reliable.
- Track query failures and improve navigation based on observed use.
- Periodically prune generated wording and topic sprawl while retaining source
  lineage.

## Validation Requirements

The lint layer should eventually check:

- Every cited repository path exists, or is explicitly marked deleted.
- Frontmatter matches the schema.
- Required page lifecycle and memory-review fields are present, and key
  statements carry evidence status.
- Every source capsule is present in state and vice versa.
- Stored content hashes match processed content.
- Topic source lists and state reverse dependencies agree.
- Internal links resolve.
- Superseded pages identify their successor when known.
- Conflicts cite at least two positions or explain why only one source remains.
- No generated page cites `atlas/` or `graph/` after cutover.
- No excluded or generated memory file is recursively ingested.

## Success Criteria

The system is successful when:

- A small source change reprocesses only its capsule, linked topics, and
  affected index/conflict sections.
- A clean repository reports no pending work, while dirty tracked sources are
  visible and untracked files remain ineligible until deliberately added.
- A newer low-authority note does not silently override an applicable source
  selected by explicit relation or reviewed adjudication.
- Explicit revisions preserve the older position and explain the transition.
- Script descriptions identify concrete inputs, outputs, and validation roles.
- Agents answer the retrieval benchmark from a small number of memory files
  with source citations.
- The memory remains understandable without its update tool and rebuildable if
  `_meta/state.json` is lost.
- `atlas/` and `graph/` can be removed without losing unique, useful research
  knowledge or day-to-day retrieval capability.

## Explicit Non-Goals for the First Version

- Reconstructing the full legacy graph.
- Ingesting every TeX, Markdown, Python, shell, or Wolfram Language file.
- Formal AST or PDF parsing.
- Automatically adjudicating scientific disagreements.
- Treating AI summaries as measurements or reviewed conclusions.
- Maintaining an Obsidian vault, reader site, Canvas, or graph visualization.
- Creating a page for every atomic claim.
- Automatically committing generated memory changes.
