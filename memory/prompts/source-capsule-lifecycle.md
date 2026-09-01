# Source-capsule deletion and retirement task

Apply `00-snapshot-contract.md`. This task handles an existing source unit that
is `source_deleted` or `retired_from_corpus`. It does not summarize a new
source snapshot.

## Inputs

The sealed task must supply the prior served capsule, prior statement-ID
registry, prior successful state, complete last-seen member identities and
committed member content/excerpts needed to retain each statement, the last
source commit, the lifecycle action commit, the desired page lifecycle
(`deleted` or `retired`), deletion/retirement scope and reason, and every
affected derived-page ID. Stop if any prior statement lacks its last-seen
direct-source support in the packet.

`deleted` means the configured source disappeared from the target tree.
`retired` means the source still exists but its unit was deliberately removed
from the configured corpus. Neither action is scientific supersession and
neither selects a different claim as true.

## Output contract

Preserve the prior capsule's retrieval-useful historical statements and stable
IDs, changing their lifecycle to the task's desired lifecycle without changing
their scientific evidence value. Mark all generated wording
`memory_review: ai_draft`. Retain last-seen original citations and annotate
them with the supplied last source commit/blob identity. Do not redirect a
missing path or infer replacement content.

Set page lifecycle and lifecycle metadata exactly from the task. Record the
action, scope, reason, last source commit, lifecycle action commit, and affected
derived pages in a status-bearing statement. Normal navigation must not present
the page as current. Obey the capsule output budget.

For a partially deleted bundle, use the normal kind-specific capsule prompt
with the task's surviving committed members plus last-seen metadata/content for
the missing members. Include a status-bearing structural-deletion statement;
do not label the whole unit deleted when its configured primary source remains.
