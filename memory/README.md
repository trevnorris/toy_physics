# Research memory

This directory is the repository's small navigation and synthesis layer for
selected sources under `research/` and `software/`. Repository sources remain
authoritative. Start every lookup at [`index.md`](index.md).

## Freshness model

Normal memory updates read the committed Git tree. Dirty tracked files should
be reported by status, but they are not incorporated into the committed memory
baseline. Untracked files are outside normal discovery. Consequently, memory
pages describe the last successfully processed committed snapshot, not
necessarily the current working tree.

Before relying on the memory, run these commands from the repository root:

```bash
python3 memory/tools/memory.py status
python3 memory/tools/memory.py lint --served
```

`status` reports initialization, the last successful snapshot, pending
committed changes, dirty tracked inputs, and prior update failures. Served
`lint` checks page structure, source scope and citations, generated-region
boundaries, links, and synchronization consistency. Until initial ingest
succeeds, the index remains visibly unpopulated rather than implicitly current.

## Reading boundaries

- Treat `memory_review: ai_draft` as unchecked memory wording, not as a
  scientific review decision.
- Follow citations to the concrete `research/` or `software/` source for
  authority. A memory page cannot create evidence or settle an open result.
- Generated regions may be replaced by synchronization. Preserve unmarked
  regions on mixed pages.
- Do not infer that a script's pass label proves a scientific interpretation;
  look for the separate source that interprets the recorded computation.

The authoritative rules are in [`_meta/schema.md`](_meta/schema.md), and the
selected corpus policy is in [`_meta/CORPUS.md`](_meta/CORPUS.md) and
[`_meta/config.yaml`](_meta/config.yaml).

## Legacy migration

The legacy Atlas and graph are temporary migration inputs and remain in place.
Their valuable status boundaries, open gates, and crosswalks must be
re-anchored to original `research/` or `software/` sources before cutover.
Deletion is gated on successful population, linting, and the
[`retrieval benchmark`](_meta/retrieval-benchmark.md); legacy pages are not
final evidence for new memory content.
