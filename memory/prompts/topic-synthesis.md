# Topic-synthesis task

Apply `00-snapshot-contract.md`. This task refreshes one declared
`working-position` region in a mixed topic page from prepared source capsules
plus committed direct-source content or deterministic excerpts for every
load-bearing statement.

## Inputs

The transaction task must supply the topic ID/title, staged base page,
generated region ID, relevant capsule paths and freshness identities,
reverse-dependency reasons, canonical-lineage records, and committed direct
source content or excerpts covering every proposed substantive statement. Do
not retrieve additional sources based on terminology. A capsule may nominate
sources, IDs, and navigation only; it cannot supply semantic support for a new
paraphrase. If direct evidence for a load-bearing statement is absent, omit it
and return a missing-input report for that statement. Do not transfer a draft
capsule's paraphrase and citation.

A stale, pending, superseded, deleted, or retired memory dependency may be
linked with an explicit warning, but cannot support a current working position.
Reject a dependency whose page or policy hash differs from the sealed task
identity.

When the sealed task contains `migration_requirements`, emit one exact
schema-version-2 structured migration record for every requirement inside the
generated region. Bind its legacy ID to page-namespaced target statement IDs
and every supplied original citation/anchor. If direct support is missing, use
`migration_disposition: blocked` and report the missing input; a bare legacy ID
or unrelated citation never counts as migrated.

## Synthesis

Write a small answer surface organized as:

1. `Current working position`
2. `Why these sources are preferred`
3. `Status boundaries and invalid inferences`
4. `Superseded or historical positions`
5. `Open questions and candidate conflicts`
6. `Relevant source capsules and scripts`

Select a current retrieval default only through the precedence rules in the
shared contract. State the exact basis and scope. Configured canonical lineage
may determine which source a reader opens first; label it as navigation and do
not present it as proof. Never infer preference from recency, target commit,
step number, word count, evidence vocabulary alone, or `PASS`.

Retain status firewalls, invalid-inference warnings, and open gates when they
affect likely queries. Do not collapse projection into reduction, architecture
into result, an engine run into interpretation, or an existence result into a
uniqueness/result claim unless the sources explicitly do so.

When capsules appear incompatible, perform the scope/convention test. If the
packet cannot resolve the difference, include both scoped positions and link a
candidate conflict; do not choose a winner. Preserve explicit supersession and
the older position. Mark every generated substantive statement
`memory_review: ai_draft` with its own evidence value.

## Output contract

Start from the staged base page. Replace only the content between the exact
`BEGIN GENERATED:working-position` and matching end marker. Preserve all
unmarked body bytes and marker lines. Update only task-delegated machine-owned
frontmatter values; `sources` is sorted and contains every direct original
source used, `generated_from_commit` is copied from the transaction, and the
page roll-up remains `memory_review: ai_draft`.

Keep the generated region concise enough to scan. Prefer several precise
status statements over a narrative survey. Memory-page links belong only in
the final navigation section; substantive citations show original mapped
`research/` or `software/` paths and stable anchors.
Obey the topic budget (default ceiling: 1,000 words and eight substantive
statements).
