# Research Memory Semantic Schema

Status: Phase 0 authoritative contract
Schema version: 2
Last updated: 2026-08-24

## 1. Authority and scope

This schema governs generated and curated content under `memory/`. Repository
sources under the selected `research/` and `software/` corpus remain
authoritative. Memory pages are navigation and synthesis aids; they cannot
create evidence, settle an open calculation, or upgrade the claim made by a
source.

The governing order is:

1. The cited repository source for scientific and software content.
2. This schema for what a memory page may say and how it must say it.
3. `_meta/config.yaml` for corpus selection, source-unit boundaries, and
   canonical lineages.
4. `_meta/state.json` for rebuildable synchronization state.

`atlas/` and `graph/` are migration inputs only. A migrated statement must cite
an original `research/` or `software/` source before the legacy system is
deleted. A legacy node ID may be retained as migration provenance, but it is
never the final evidence citation.

Meta documents (`PLAN.md`, `_meta/*.md`, prompts, and measurement records) are
not ordinary memory pages and do not require the page frontmatter below.

## 2. Terms

- **Source file:** one Git-tracked repository path admitted by
  `_meta/config.yaml`.
- **Source unit:** the unit of refresh. It is either one source file or a
  declared bundle of related files.
- **Bundle member:** a concrete repository path in a source unit. Members have
  roles and are hashed separately.
- **Source capsule:** a replace-on-refresh report of what one source unit says.
- **Topic page:** a synthesis across source capsules and direct source
  citations.
- **Script catalog:** a domain-grouped inventory of executable entry points and
  result families.
- **Candidate conflict:** two or more statements that may be incompatible but
  have not yet been adjudicated.
- **Substantive statement:** a scientific conclusion, status judgment,
  software/result interpretation, precedence decision, or claim that an
  artifact computed or validated something. Pure navigation text is not a
  substantive statement.

Repository paths are POSIX-style and relative to the repository root. Page and
statement IDs are lowercase kebab-case and stable after publication. Dates use
`YYYY-MM-DD`; a partial source date may remain `YYYY` or `YYYY-MM` and must not
be made artificially precise.

## 3. Common page frontmatter

Every source capsule, topic page, script catalog, and conflict register begins
with these fields:

```yaml
---
schema_version: 2
id: <stable-kebab-case-id>
title: <human-readable title>
type: source_capsule | topic | script_catalog | conflict_register | index
lifecycle: current | superseded | deleted | retired
memory_review: ai_draft | reviewed
sources:
  - research/or/software/path
content_owner: ai_generated | human_curated | mixed
last_updated: YYYY-MM-DD
---
```

Rules for the common fields:

- `sources` contains every direct repository source used by the page, in
  sorted, duplicate-free order. It does not contain other memory pages.
- Page lifecycle and memory review describe different things and must never be
  collapsed into a single `status` or numeric priority.
- Evidence is statement-level. A page commonly mixes measured, derived,
  provisional, open, and disputed statements, so frontmatter must not assign
  one evidentiary status to the whole page.
- `memory_review` describes the current memory wording, not whether a cited
  paper, result, adjudication, or script was reviewed.
- `last_updated` is the date the memory wording was refreshed. It is not a
  source date and cannot be used for precedence.
- `generated_from_commit` is required on `ai_generated` pages and on mixed
  pages with generated regions. It is the full committed source snapshot read
  from the transaction, never an inferred working-tree identity.
- `superseded_by` is required when `lifecycle: superseded` and a successor is
  known. `supersedes` is optional on the successor. Both contain page or
  statement IDs, not free-form titles.
- `reviewed_by`, `reviewed_at`, and `review_record` are required when
  `memory_review: reviewed`. `review_record` cites the repository or memory
  record that documents review of this wording.

`open` and `disputed` are not lower-confidence versions of `derived` or
`measured`; they mean the object is unresolved or has unresolved incompatible
treatments.

## 4. The three independent status dimensions

### 4.1 Lifecycle

- `current`: the page or statement is part of the current memory view. It does
  not mean the physics is correct or closed.
- `superseded`: a newer or preferred item explicitly replaces this item for a
  stated scope. Preserve the old text and identify the successor and basis.
- `deleted`: the source path is no longer present in the selected Git-tracked
  corpus. Preserve enough metadata and lineage to explain existing citations.
- `retired`: the source still exists, but its unit was deliberately removed
  from the configured corpus. Preserve its last served page and state as audit
  history; normal query does not rank it as current.

Do not set `superseded` merely because another source is newer, longer, or more
confidently written. Partial revision is statement-level: do not supersede an
entire source capsule when only one conclusion changed.

### 4.2 Evidence

- `measured`: the prepared evidence supplies the exact recorded invocation,
  readable literal stdout or a tracked transcript containing the relevant
  operands and residual, and a separate source that owns the exact
  interpretation. Cite the executable/invocation record, output record and
  stable tag or nearby identifier, and interpretive source. A prose
  measurement report qualifies as output evidence only when it embeds that
  literal chain. Code without its recorded invocation, an asserted zero, a
  pass label, an exit code, or a guard without the operands/residual does not
  establish a measured conclusion.
- `derived`: the cited source explicitly derives the statement analytically or
  proves it inside declared assumptions/closure. The memory cites the
  derivation or theorem anchor and states the applicable scope. AI reasoning
  performed while writing the capsule is not a source derivation.
- `provisional`: the source asserts or motivates the statement, but the memory
  cannot point to the evidence required for `measured` or the explicit
  derivation required for `derived`; it is also the default for a new AI
  synthesis when support is unclear.
- `open`: the cited source explicitly leaves the object unresolved, deferred,
  or dependent on an unbuilt calculation, unselected branch, or missing datum.
  Do not fill the gap by inference.
- `disputed`: apparently incompatible statements about the same object, scope,
  regime, and conventions remain unresolved. The conflict record must cite all
  live positions.

These values are not a universal linear ranking. A literal software measurement
is strong evidence about that run; an analytic derivation may be the applicable
evidence for a theorem. Compare evidence only for the same object and scope.
`provisional` never outranks applicable explicit derivation or measurement, and
`open` or `disputed` is never presented as settled.

For computed material, preserve the repository's evidence boundary:

- A script catalog says what objects the program constructs, prints, compares,
  and guards. It does not turn its own `PASS` or exit code into a physical
  conclusion.
- A scientific interpretation is attributed to the paper, step record, result
  report, or adjudication that makes it.
- Cross-engine agreement is evidence about the compared objects over the
  comparator's actual coverage. It is not automatically independent physical
  evidence when engines share premises or an action.
- Missing commands, literal outputs, operands, residuals, or interpretation
  records are reported as missing rather than reconstructed in prose.

### 4.3 Memory review

- `ai_draft`: generated or changed memory wording has not received the recorded
  review required to call the wording checked.
- `reviewed`: the exact memory wording and its citations have been checked, with
  `reviewed_by`, `reviewed_at`, and `review_record` recorded.

An AI may initialize or reset `memory_review` to `ai_draft`; it may not infer `reviewed`
from a source's review history or from its own rereading. Any substantive change
inside a reviewed generated region resets the affected statement and the page
roll-up to `ai_draft`. Unchanged human-curated regions may retain their own
recorded review status.

### 4.4 Statement status

Every substantive statement uses this compact form. This includes scientific
conclusions, status judgments, artifact-behavior claims, limitations, open
questions, revision/supersession claims, preferred positions, and every
conflict position or resolution. Only neutral navigation and member metadata
may appear without a statement status block.

```markdown
### <stable statement ID> — <short name>

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

<Conservative statement, including assumptions and scope.>

Sources:

- `research/example/paper.tex` — `\label{sec:exact-anchor}`
- `research/example/record.md` — heading `## Interpretation`
```

A paragraph may inherit the page default only when it is non-load-bearing
context. Do not use the inherited default to hide mixed or weaker evidence.

## 5. Stable citations

Every substantive statement cites the concrete source member that supports it.
Use the first available anchor type:

1. Repository path plus TeX `\label{...}`.
2. Repository path plus exact Markdown heading text or explicit `{#anchor}`.
3. Repository path plus a distinctive nearby identifier such as a step ID,
   emitted tag, function name, theorem name, or result key.
4. A line number only as an additional disposable convenience, never as the
   only anchor.

Canonical display examples:

```markdown
- `research/pde/paper/pde.tex` — `\label{sec:status-parent}`
- `software/stage1_solver/STAGE1_VERDICT.md` — heading `## The test`
- `research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py` — function `main`
- `software/em_charge_attribute/puncture_deflection_electric_sign_check.py` — function `section4_adjudicate`
- `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` — heading `## Concise result (body; ≤2 pages)`
```

Clickable Markdown links may accompany this display, but the visible
repository-relative path and stable anchor must remain present for linting.
When a source unit is a bundle, cite the responsible member—not merely the
bundle ID or capsule. If a statement combines several members, cite each role:
for example the engine, literal transcript, comparator, and interpretive step
record.

If the prepared packet supplies no stable anchor, use the exact sentinel
`anchor-unavailable` and keep the statement `evidence=provisional` (or `open`
only when the source explicitly says it is open). A statement using that
sentinel cannot be `measured` or `derived`.

For a deleted path, cite the last tracked path plus the stored last-seen commit
or blob hash and mark the citation deleted. Do not silently redirect it to a
different file.

## 6. Source units and bundle rules

A file unit contains exactly one selected path. A bundle groups files that must
be understood and refreshed together, such as a componentized paper, a step
record with paired engines and outputs, or a software run/result family.

Each resolved bundle member has:

```yaml
- path: research/or/software/path
  role: <configured stable role string>
  read_mode: semantic | excerpt | identity_only
  mode: <Git mode>
  object_type: blob
  blob_oid: <full Git blob object ID at the processed commit>
  blob_size: <bytes>
```

Rules:

- Every member is Git-tracked and admitted by `_meta/config.yaml`, except a
  deleted member retained for lineage.
- Every member is hashed independently. Any member change refreshes its owning
  capsule and linked topic/script regions.
- Exact paths and the schema-defined prefix/recursion/name selectors may define
  a bundle in config. Generic or Git pathspec globs are not supported. A
  capsule and state record contain the resolved member list used for refresh.
- Process history is citation-only by default. A directive, prompt, review, or
  red-team record becomes a member only when required to understand a selected
  result or dispute.
- Measurement files normally support an owning unit; they do not automatically
  receive their own capsule.
- Helpers and tests are dependencies of a selected script family unless they
  are meaningful entry points in their own right.
- A source capsule reports each member's role and keeps member-specific
  citations. It never synthesizes a bundle into a stronger claim than its
  interpretive source makes.
- Roles are stable strings copied from config, not an inferred closed
  vocabulary. Common roles include primary, paper member, governance, step
  record, paper card, engine, comparator, control, premises, measurement,
  result, adjudication, entry point, dependency, and supporting history.

## 7. Source precedence and supersession

Select a preferred current position only in this order:

1. **Explicit source relationship.** A source explicitly says it revises,
   corrects, supersedes, retracts, or overturns another statement. Cite the
   relationship and limit it to the declared scope.
2. **Reviewed evidence adjudication.** A recorded review compares the same
   object and scope and selects the treatment with the stronger applicable
   evidence. Cite the adjudication and evidence records.
3. **Configured canonical lineage.** `_meta/config.yaml` designates the
   canonical source for that document lineage. This chooses the retrieval
   default; it does not manufacture scientific support.
4. **Recency tie-break.** Use an explicit source date only between sources of
   comparable role, evidence, source adjudication, scope, and lineage.
5. **Unresolved.** Retain all positions and create a candidate conflict.

Git commit time, file modification time, page length, filename suffix, or a
numeric AI-assigned priority is not a scientific precedence rule.

A supersession record contains:

- old and new page or statement IDs;
- both source citations and exact anchors;
- `basis: explicit_source | reviewed_adjudication | configured_canonical`;
- scope of replacement;
- concise explanation;
- date and review state.

Never delete the older position merely because it is superseded. If the basis
is ambiguous or only inferred from chronology, record a candidate conflict
instead.

## 8. Candidate-conflict policy

AI conflict detection is deliberately conservative. Before opening a
candidate, check whether the difference is explained by object identity,
notation, sign convention, coordinate/gauge choice, approximation order,
regime, branch, source role, or historical versus current scope.

An AI may:

- create `conflict_state: candidate`;
- paraphrase positions conservatively with source anchors;
- list plausible scope/convention explanations;
- apply an explicit precedence rule whose basis is directly cited;
- draft a source-declared resolution while leaving memory wording
  `memory_review: ai_draft`.

An AI may not:

- select a winner because one file is newer;
- rederive the physics to adjudicate the conflict;
- mark a candidate reviewed;
- erase the losing or historical position;
- describe different scopes as contradictions without recording the scope
  test.

A conflict becomes resolved only by an explicit source resolution or a
recorded reviewed adjudication. Preserve the resolution basis and history.
Recently resolved conflicts remain visible; older ones may move to a historical
section but are not silently removed.

### 8.1 Legacy migration accounting

Every inventory item with `disposition: migrate` receives one structured record
inside its target page's declared generated region:

```markdown
#### Migration record — <legacy ID>

- `legacy_id`: `<exact inventory ID>`
- `migration_disposition`: `migrated | blocked`
- `target_statement_ids`: [`<page-namespaced statement ID>`, ...]
- `inventory_sha256`: `<frozen inventory digest>`

Sources:

- `research/or/software/path` — `<stable anchor>`
```

Every required original source/anchor supplied by the sealed migration task is
present in that record or explicitly accounted as blocked. A legacy ID merely
appearing in prose, a comment, or curator text is not completion. The record
points to current status-bearing statements; it does not use Atlas wording as
evidence.

`migration_items_accounted` means all such records validate. It is not
`cutover_ready`. Cutover additionally requires full state freshness, served
lint and hashes, no Atlas/graph evidence links or normal references, a recorded
pass for every retrieval benchmark item, and a fresh-agent retrieval check.

## 9. Ownership and update boundaries

Ownership values have operational meaning:

- `ai_generated`: the synchronizer may replace the generated body from current
  sources. Source capsules use this mode by default. Manual edits inside the
  generated body will be overwritten. This frontmatter value explicitly marks
  the entire body as generated, so whole-body pages do not also use region
  comments.
- `human_curated`: automation may lint but does not rewrite the body.
- `mixed`: generated regions are replaceable; all unmarked regions are
  human-owned and preserved byte-for-byte.

Mixed pages use exact markers:

```markdown
<!-- BEGIN GENERATED:<stable-region-id> -->
...replaceable content...
<!-- END GENERATED:<stable-region-id> -->
```

Region IDs are unique within a page. Nested regions are forbidden. A missing or
unbalanced marker is a lint error and blocks update rather than risking human
content.

Ownership by artifact:

- `_meta/config.yaml`, this schema, the benchmark, and curated migration
  decisions are human-curated policy.
- `_meta/state.json`, content hashes, resolved bundle membership, last processed
  commit, and deterministic index lists are machine-owned.
- Source-capsule summaries and candidate-conflict drafts are AI-generated.
- Topic working-position regions, script entries, conflict candidates, and
  index summaries may be generated inside mixed pages.
- Review dispositions, manual precedence adjudications, and curator notes are
  human-owned unless a review workflow explicitly records another owner.

The deterministic layer may update machine-owned frontmatter fields but cannot
change statement evidence, memory review, scientific precedence, conflict
resolution, or curator notes. The semantic extractor may conservatively set
statement evidence, initialize `memory_review: ai_draft`, and draft source
relationships; it cannot set `reviewed`.

### 9.1 Sealed semantic tasks

Every semantic task is self-contained. Its sealed packet includes frozen
schema and prompt copies; target and refresh identities; page/task metadata;
complete current and required prior member identities; semantic blobs,
deterministic excerpts, and identity-only metadata; allowed citations;
stable-ID registry; desired lifecycle; related/reverse dependencies; staged
base page for mixed ownership; delegated frontmatter fields; dependency page
hashes/freshness; and an output budget.

The writer runs in a Gitless isolated filesystem with only the packet mounted
read-only and one declared output mounted read-write. Live repository paths are
absent. The packet contains a pre-run isolation contract; a trusted launcher
creates the post-run attestation only after exit. Publication requires that
attestation and verification that only the declared output changed. Prompt
prose forbidding a live read is not a substitute for this boundary.

For a mixed page, the transaction seals the base page and digests every
unmarked byte and non-delegated frontmatter value. Finalization rejects any
change outside the declared generated region/delegated fields.

Derived topic and conflict statements require target-snapshot direct-source
content or deterministic excerpts for every load-bearing position and
resolution. Source capsules may nominate original paths, stable IDs, and
navigation, but their `ai_draft` paraphrases do not supply semantic authority.
A memory dependency records page hash, generation commit, unit digest when
applicable, policy digests, lifecycle, and freshness at target. Stale or
pending pages may be linked with warnings but cannot support a current working
position.

### 9.2 Output budgets and stable IDs

Config supplies enforceable ceilings and each task copies its applicable
budget. Initial defaults are:

- source capsule: 1,800 words and 12 substantive statements;
- topic generated region: 1,000 words and 8 substantive statements;
- script entry: 250 words; at most 10 entries and 2,500 words per catalog;
- conflict entry: 350 words; at most 5 entries and 1,800 words per register;
- index navigation region: 700 words.

When the source exceeds a budget, prioritize benchmark-relevant results,
status boundaries, open gates, and invalid-inference guards. Do not overflow or
silently create extra pages.

Statement IDs are globally unique and page-namespaced as
`<page-id>--<semantic-key>`. Refreshes receive the prior ID registry and reuse
the same ID for the same semantic object. New IDs follow the task's
deterministic allocation rule and never derive identity from mutable wording
alone.

## 10. Page templates

### 10.1 Source capsule

```markdown
---
schema_version: 2
id: source-<unit-id>
title: <source-unit title>
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
  - research/example/record.md
content_owner: ai_generated
last_updated: YYYY-MM-DD
generated_from_commit: <full Git commit SHA>
source_kind: paper | governance | step_family | result | software_project | verification
source_unit:
  id: <config unit ID>
  shape: file | bundle
  entrypoint: research/example/record.md
  unit_digest_sha256: <canonical unit digest>
  members:
    - path: research/example/record.md
      role: step_record
      read_mode: semantic
      mode: "100644"
      object_type: blob
      blob_oid: <full Git blob object ID>
      blob_size: <bytes>
extractor_version: <version>
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

## Source-unit map

## Key statements

### <statement ID> — <short name>

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

<Statement with its assumptions and limits.>

Sources:

- `<member path>` — `<stable anchor>`

## Computed evidence represented by the source

List commands, engines, literal transcript/result paths, stable output tags,
and the separate source that interprets them. Say `not supplied` where absent.

## Assumptions, exclusions, and open questions

## Revision and supersession relationships

## Related topics and scripts
```

For a monolithic paper, `shape: file` and `members` contains the paper path. For
a componentized paper, the TeX entry point and every resolved `\input` member
belong to the bundle. For a step family, distinguish the interpretive record,
paper card, independent engines, comparator, exports, and literal outputs.

### 10.2 Topic page

```markdown
---
schema_version: 2
id: topic-<topic-id>
title: <topic title>
type: topic
lifecycle: current
memory_review: ai_draft
sources:
  - research/example/current.tex
  - research/example/historical.tex
content_owner: mixed
last_updated: YYYY-MM-DD
generated_from_commit: <full Git commit SHA>
---

<!-- BEGIN GENERATED:working-position -->
## Current working position

### <statement ID> — <short name>

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

<Working position, scope, and citations.>

## Why these sources are preferred

Name the applicable precedence rule for each preference. Do not write only
"newer" or "higher priority."

## Status boundaries and invalid inferences

## Superseded or historical positions

## Open questions and candidate conflicts

## Relevant source capsules and scripts
<!-- END GENERATED:working-position -->

## Curator notes

Human-owned; generated updates preserve this section.
```

Keep topic count small. A topic is justified by repeated retrieval value, a
benchmark question, an important status firewall/open gate, or synthesis across
several selected sources—not merely because a term occurs in the corpus.

### 10.3 Script catalog

```markdown
---
schema_version: 2
id: scripts-<domain-id>
title: <domain> script catalog
type: script_catalog
lifecycle: current
memory_review: ai_draft
sources:
  - research/example/run_audit.py
  - research/example/step_record.md
content_owner: mixed
last_updated: YYYY-MM-DD
generated_from_commit: <full Git commit SHA>
domain: <domain-id>
---

<!-- BEGIN GENERATED:entries -->
## `<repository/path/to/entrypoint>`

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

- Language: Python | Wolfram Language | shell | other
- Source unit: `<unit-id>`
- Role: builder | validator | comparator | runner | reporter | visualizer
- Computes/builds: <named objects, not a prose conclusion>
- Inputs: <paths, arguments, fixtures, imported exports>
- Emits: <stdout tags and tracked output/result paths>
- Guards/controls: <what can fail and the measured scope>
- Invocation: <exact recorded command, or `not recorded`>
- Interpretation source: <paper/step/result path and anchor, or `none`>
- Related topics: <memory links>
- Refreshed against: <commit and member hashes>

Sources:

- `<script path>` — function/tag/nearby identifier
- `<output path>` — emitted tag or nearby identifier
- `<interpretive source>` — heading or TeX label
<!-- END GENERATED:entries -->

## Curator notes
```

Catalog the meaningful entry point and its family, not every helper. If a
runner coordinates paired engines, comparator, fixtures, and a result report,
one entry may describe the family and list those dependencies. Never summarize
a script as “proves” or “shows” a scientific conclusion unless an interpretive
source makes that scoped statement.

### 10.4 Conflict register and entry

```markdown
---
schema_version: 2
id: conflicts
title: Research memory conflict register
type: conflict_register
lifecycle: current
memory_review: ai_draft
sources:
  - research/example/older.md
  - research/example/newer.md
content_owner: mixed
last_updated: YYYY-MM-DD
generated_from_commit: <full Git commit SHA>
---

<!-- BEGIN GENERATED:candidates -->
## <conflict ID> — <neutral title>

Conflict state: `candidate | open | resolved | scope_difference`

### Positions

#### <page-namespaced position A ID> — <neutral name>

Status: `lifecycle=current` · `evidence=disputed` · `memory_review=ai_draft`

<Conservative position A paraphrase.>

Sources:

- `<path>` — `<anchor>`

#### <page-namespaced position B ID> — <neutral name>

Status: `lifecycle=current` · `evidence=disputed` · `memory_review=ai_draft`

<Conservative position B paraphrase.>

Sources:

- `<path>` — `<anchor>`

### Scope/convention test

Record object, regime, branch, notation, sign, gauge/coordinate, approximation,
and source-role differences that may explain the apparent incompatibility.

### Precedence check

Record each applicable schema rule and its cited basis. If none resolves it,
write `unresolved; no preferred position`.

### Disposition and history

#### <page-namespaced disposition ID> — <neutral disposition name>

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

For an unresolved entry, write `unresolved; no preferred position`. For a
resolution, record `resolution_basis: explicit_source |
reviewed_adjudication`, the preferred statement ID, date, original-source
citations, replacement scope, and why the old position remains visible.
<!-- END GENERATED:candidates -->

## Curator dispositions
```

For `candidate` and `open`, live incompatible positions use `disputed`. For
`resolved` and `scope_difference`, preserve each position's own justified
evidence and lifecycle instead of assigning an aggregate disputed status; the
separate disposition statement carries the evidence value justified by its
explicit source resolution/adjudication. If all live conflicts are resolved,
each retained entry still keeps its historical positions. A `scope_difference`
disposition explains why no actual contradiction remains and cites the direct
scope evidence.

### 10.5 Index

The index uses the common frontmatter with `type: index`, `content_owner:
mixed`, and one `GENERATED:navigation` region. Its generated region contains
neutral navigation and freshness warnings, not a second scientific synthesis.
Any unavoidable substantive description uses the statement form and direct
original-source citations. The unmarked reader guidance and curator notes are
preserved byte-for-byte.

## 11. Minimal lint invariants

Phase 1/2 linting must enforce at least:

- common required page fields and controlled values, plus statement-level
  evidence values;
- unique page IDs and statement/conflict IDs;
- `reviewed` provenance and `superseded` successor requirements;
- repository-relative, Git-tracked source paths, or an explicit deleted record;
- concrete member citations for bundle-backed statements;
- resolved bundle membership and per-member full hashes on source capsules;
- balanced, unnested generated-region markers;
- byte-for-byte preservation of unmarked mixed-page content and
  non-delegated frontmatter;
- isolation attestation, declared-output-only verification, and applicable
  output-budget compliance;
- no AI-written `reviewed` transition without a review record;
- no conflict resolution without an explicit-source or reviewed basis;
- no generated final citation to `atlas/` or `graph/` after cutover;
- no script interpretation presented as measured without command, literal
  output, and interpretation provenance;
- no `measured` or `derived` statement using `anchor-unavailable`;
- internal links, source-capsule/state agreement, and configured corpus scope.

When lint cannot decide whether an evidentiary or precedence requirement is
met, it warns and leaves the content `ai_draft`; it does not silently upgrade or
delete it.
