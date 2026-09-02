# Deterministic Scanner and State Decisions

Status: reviewed and folded once; builder-ready
Date: 2026-08-24
Updated: 2026-09-01 for git-annex storage and the PDE-ledger corpus exclusion
Review record: `memory/_measurements/scanner_decision_reviews.md`

This is the builder-facing contract for the deterministic memory tooling. It
contains the verified fold from two independent review legs.

## Committed-source boundary

1. Normal discovery is limited to committed objects under `research/` and
   `software/`. Root `notes/` and `docs/` are excluded.
   The roots `research/pde_ledger/`, `research/pde_ledger_v2/`, and
   `research/pde_ledger_v3/` are also absolute exclusions: normal discovery,
   exact membership, selectors, derived direct sources, supporting lineages,
   and Atlas migration may not restore them.
2. `atlas/` and `graph/` are visible only to the explicit migration workflow.
3. Enumerate the target commit once with NUL-safe Git plumbing equivalent to
   `git ls-tree -rz --full-tree <commit> -- research software`.
4. Record each entry's mode, object type, object ID, size, and path. Regular
   blobs are supported. Symlinks, gitlinks, and other types fail unless a unit
   explicitly declares identity-only handling.
5. Do not use Git glob pathspecs; this repository's Git does not support the
   required magic on `ls-tree`.
6. Dirty tracked paths are reported separately and never alter the committed
   baseline. Untracked paths are ineligible until committed.
7. `.gitignore` controls new/untracked eligibility, not committed-tree
   membership. Explicit unit selectors and exclusions control committed paths.

## Executable source-unit configuration

8. Config contains one normalized `units` list. Every unit has a stable ID,
   kind, unique capsule target, lifecycle, and explicit members/selectors.
9. Exact members declare `path`, `role`, `read_mode`, and whether required.
10. The only selector dialect is:
    `prefix` + `recursive` + optional `suffixes`/`basenames`/`name_prefixes`.
    Matching is segment-aware, case-sensitive, and tested locally after the
    committed tree is enumerated.
11. Required selectors resolving empty on initial configuration fail. A member
    present in prior state but absent at target becomes a deletion action.
12. One path has at most one primary owner. Other units may name it only as a
    `shared_dependency`.
13. Unit IDs and capsule targets are unique. Topic, script-catalog, conflict,
    and index pages are derived aggregations with reverse dependencies, not
    shared capsule targets.
14. Directives, prompts, review history, tests, helpers, and measurements are
    dependencies by default. Exact unit membership may override a broad
    semantic exclusion when the source is load-bearing, except for the three
    absolute PDE-ledger exclusions above.
15. Scientific unit identity is explicit. Filename prefixes never infer that
    similarly named records are separate or equivalent research steps.

## Read modes and bounded evidence

16. `semantic` material is staged in full and may be summarized.
17. `identity_only` material is fully hashed but its contents are not supplied
    to the AI. This is the default for large transcripts, exports, and binary
    artifacts whose interpretation already exists in a selected record.
18. `excerpt` material is supplied only through declared stable tags/headings
    and a bounded context policy. No arbitrary model-chosen chunking is used.
19. Config declares default per-member and per-unit semantic byte limits. A
    semantic member exceeding a limit fails unless it has an explicit smaller
    override or excerpt rule.
20. Tracked `.out` files are valid evidence dependencies. They are never
    blanket-excluded by extension. When a transcript is stored as a git-annex
    symlink, the committed tree contains only the annex pointer; configure that
    member as `identity_only` rather than resolving worktree-local annex
    content. Its owning selected record and executable engines remain the
    portable semantic evidence.

## Identity and drift

21. Each member identity includes path, role, read mode, Git mode/type, blob
    object ID, and blob size.
22. A unit digest is SHA-256 over a versioned canonical serialization of its ID,
    semantic configuration, and sorted resolved member identities.
23. State and transactions pin content digests for config, schema, extraction
    prompts, and the tool/extraction contract. A manual version string alone is
    insufficient.
24. Rename detection is informational. Correctness follows old/new configured
    membership and blob identity.

## Immutable transaction inputs

25. `update --prepare` materializes the exact target-commit blobs into an
    ignored transaction snapshot. The AI reads only that snapshot, never live
    `research/` or `software/` paths.
26. The transaction maps snapshot paths back to visible repository-relative
    citations. Generated pages cite original paths, not staging paths.
27. Identity-only members appear only in the transaction manifest. Excerpts are
    generated deterministically from committed blobs and declared anchors.
28. The extraction prompt explicitly forbids reading live candidate roots.

## Staging, publication, and recovery

29. AI-generated pages are written under the transaction staging tree. Live
    memory pages remain untouched during extraction.
30. Finalize lints and hashes the complete staged page set before publication.
31. Finalize takes a repository-local lock, then rechecks HEAD, prior state,
    config/schema/prompt/tool digests, transaction completeness, and staged page
    hashes.
32. Publication uses a durable journal and backups: write journal; replace
    target pages; write state last; mark complete; then clean backups.
33. Handled publication failure restores every old page and prior state.
34. Every CLI command first recovers an interrupted journal: roll back when the
    new state was not committed, or finish cleanup when it was.
35. Query and lint reject live pages whose recorded hash/generation is not in
    successful state. They never read staging pages.
36. State and journal writes use same-directory temporary files, flush, and
    atomic replacement.

## Partial and full baselines

37. State records a processed commit and digest per unit.
38. A filtered update may process selected configured units, but does not mark
    unrelated pending units complete.
39. The global `last_fully_processed_commit` advances only when every active
    unit matches the same target commit and current policy digests.
40. `--paths` selects owning units; it never expands beyond configured units or
    silently changes the global baseline.
41. A unit removed from config is `retired_from_corpus`; a configured source
    missing from the tree is `source_deleted`. Both remain auditable until
    dependent pages are deliberately retired. A hard corpus exclusion removes
    the semantic capsule and queryable IDs from served state; ignored sealed
    transactions may retain non-served provenance.

## CLI behavior

42. One repository-local Python CLI provides `init`, `status`, `update`,
    `lint`, and `query`.
43. It uses subprocess argument arrays, never shell strings built from config.
44. Config/frontmatter YAML uses `SafeLoader` plus a custom mapping constructor
    that rejects duplicate keys.
45. `init` creates only missing skeleton/state files.
46. `status` is read-only and reports unit actions, partial/full freshness,
    dirty tracked members, policy drift, missing members, and transaction state.
47. `update --prepare` creates a transaction; `update --finalize <transaction>`
    publishes it. Neither command commits Git changes.
48. `lint` has explicit modes: served pages validate against successful state;
    staged pages validate against their transaction target.
49. `query` retrieves relevant served pages/excerpts for an agent to answer. It
    does not call a model or claim to adjudicate science.
50. No command parses PDFs, edits candidate sources, commits, or deletes legacy
    directories.

## Lint and query contract

51. Lint validates required frontmatter, vocabulary, unique page/statement IDs,
    unique unit targets, state coverage, resolved internal links, and ownership
    markers.
52. Original-source citations are checked against the applicable committed
    blobs, including TeX labels, Markdown headings/anchors, or declared
    distinctive identifiers when present.
53. Lint verifies membership metadata and blob identities, but never infers or
    promotes scientific status.
54. Superseded statements identify a successor when known. Conflicts retain all
    positions and a resolution basis if resolved.
55. After cutover, no normal page cites `atlas/` or `graph/` as evidence.
56. Query warns on any pending unit unless relevance is provable. Newly
    introduced vocabulary must not suppress a freshness warning.
57. Exact ID/path/title matches rank before metadata/body matches. Current
    reviewed memory wording may rank higher, but historical/disputed pages
    remain discoverable.

## Ownership

58. Config, schema, migration decisions, and reviewed adjudications are curated
    policy.
59. State, transactions, resolved membership, hashes, and deterministic lists
    are machine-owned.
60. Source capsule bodies are generated. Topic/script/conflict/index pages may
    contain explicit generated regions; unmarked content is preserved.
61. Regeneration never heuristically merges prose. Missing or unbalanced region
    markers block publication.

## Implementation and acceptance

62. Put the CLI in `memory/tools/memory.py` and focused tests in
    `memory/tools/tests/`.
63. Separate Git/tree identity, selector resolution, transaction publication,
    YAML/frontmatter lint, and retrieval into testable functions.
64. Tests use isolated temporary Git repositories and do not receive frozen
    corpus result counts.
65. Demonstrate: ignored-volume independence; selector root/nested behavior;
    duplicate/overlap rejection; mode/type hashing; dirty/live isolation;
    policy drift rejection; addition/deletion/retirement; partial-baseline
    behavior; staged failure rollback and interrupted-journal recovery;
    anchor/link lint; and stale-query warnings.
