# Conflict-detection and register-refresh task

Apply `00-snapshot-contract.md`. This task refreshes the `candidates` region of
the mixed conflict register from a bounded set of prepared topic/capsule
nominations plus committed direct-source content or deterministic excerpts for
every position, scope claim, and resolution.

## Inputs

The transaction task supplies the staged base register, existing conflict IDs
and history, affected statements, their original citations, applicable
canonical-lineage records, and any explicit source resolution or reviewed
adjudication. Each position and resolution must have target-commit direct
source bytes in the packet; a topic/capsule may nominate it but cannot lend
semantic authority. If direct support is absent, omit the proposed entry and
return a missing-input report. Reject stale or hash-mismatched memory
dependencies. Do not scan unrelated pages or sources for contradictions.

## Detection and adjudication limits

Create a candidate only when two positions appear incompatible for the same
object and potentially overlapping scope. Before doing so, compare:

- object identity and definitions;
- notation and sign conventions;
- coordinate/gauge choice;
- approximation order and regime;
- selected branch and boundary data;
- source role and historical/current scope.

For each entry retain neutral paraphrases of all positions and their original
source anchors. Record plausible scope explanations without deciding that one
is correct. `Conflict state: scope_difference` requires supplied source
evidence that the scopes differ; otherwise leave it candidate/open.

Resolve only from an explicit source correction/resolution or a supplied
recorded reviewed adjudication of the same object and scope. Record
`resolution_basis`, preferred statement ID, date, citations, replacement
scope, and why the older position remains visible. Configured canonical
lineage can choose retrieval order but cannot resolve science. Recency, a
successful command, `PASS`, stronger prose, or your own derivation cannot
resolve a conflict.

Preserve existing unresolved and recently resolved entries unless the task
explicitly carries a lifecycle action. A deleted/retired source does not erase
history; use its supplied last-seen identity and deletion status. Never mark
the wording reviewed.

## Output contract

Start from the staged base register. Replace only the content between the exact
`BEGIN GENERATED:candidates` markers, retaining existing stable conflict IDs
for the same positions. Preserve all unmarked body bytes and markers. Update
only delegated machine-owned frontmatter; the page and all generated entries
remain `memory_review: ai_draft`.

Each entry follows the schema-version-2 order: neutral heading, conflict state,
status-bearing positions, scope/convention test, precedence check, and a
separate status-bearing disposition/resolution statement. Use `disputed` only
on candidate/open positions. Resolved or scope-difference entries retain each
historical position's justified evidence and give the resolution/scope
statement its own source-justified evidence value. If no rule resolves it,
write exactly `unresolved; no preferred position`. Every substantive position
and resolution cites original mapped paths, never transaction paths, capsules,
`atlas/`, or `graph/`.
Obey the register budget (default ceilings: 350 words per entry, five entries,
and 1,800 words for the generated region).
