# Index-refresh task

Apply `00-snapshot-contract.md`. This task refreshes only the `navigation`
region of the mixed research-memory index after staged capsules, topics,
script catalogs, and conflicts are complete.

## Inputs

The transaction task supplies the staged base index plus deterministic page
metadata and links for served/staged current pages, unit freshness, recent
changes, pending units, open gates, and conflicts. Use those lists as
authoritative. Do not inspect source content or invent categories from titles.

## Navigation content

Produce a compact entry point in this order:

1. snapshot/freshness notice, including target commit and any pending or dirty
   warning supplied by the task;
2. major topic links with one sentence about retrieval scope, not a new
   scientific conclusion;
3. current source capsules or canonical lineages needed for direct lookup;
4. open gates and candidate/open conflicts with links;
5. script-catalog domains;
6. recently changed units from the deterministic transaction list.

The index routes readers; it does not adjudicate or restate detailed physics.
When a one-sentence description is substantive, give it statement status and
original mapped citations or replace it with neutral navigation text. Preserve
warnings about `ai_draft`, pending sources, and identity-only evidence. Do not
call the memory current when the transaction reports partial coverage.

Rank exact and configured retrieval defaults before broader related pages, but
do not describe retrieval order as scientific authority. Historical,
superseded, disputed, and open material remains discoverable through the
appropriate topic/conflict links.

## Output contract

Start from the staged base index. Replace only the content between the exact
`BEGIN GENERATED:navigation` markers and preserve every unmarked body byte and
both markers. Update only delegated machine-owned frontmatter fields. The page
roll-up remains `memory_review: ai_draft`; copy `generated_from_commit`
exactly.

Use only relative links explicitly supplied in the task, verify that each link
target is in the staged page set, and keep the region short enough to serve as
the first file a new agent reads. Never cite a transaction path, `atlas/`, or
`graph/`.
Obey the index budget (default ceiling: 700 words).
