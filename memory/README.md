# Research memory

This is a small, human-readable guide to the toy-model research in
`research/` and `software/`. The original files remain authoritative. Start
with [`index.md`](index.md), then follow links to a topic, source capsule, or
script catalog.

The memory intentionally excludes root `notes/`, root `docs/`, and all three
PDE-ledger trees:

- `research/pde_ledger/`
- `research/pde_ledger_v2/`
- `research/pde_ledger_v3/`

## Updating it

The updater is content-based. It does not store or compare Git commit IDs.
Rebases and unrelated commits therefore have no effect.

```bash
python3 memory/update.py status
```

For each changed or new source unit:

1. Read the listed source files.
2. Edit its capsule directly under `memory/sources/`.
3. Update any affected topic or script pages reported by `status`.
4. Cite the original repository paths, not another memory page, for factual
   claims.
5. Ask Grok for a second look when the update contains meaningful scientific
   synthesis or conflict resolution.
6. Record completed work:

```bash
python3 memory/update.py mark-units <unit-id> [...]
python3 memory/update.py mark-pages <page-id> [...]
python3 memory/update.py lint
```

`mark-*` only records SHA-256 hashes of current file contents. Git is used
solely to determine which source files are tracked, so ignored and untracked
files are never pulled into the catalog accidentally.

Source units are curated reading bundles, not transitive build-dependency
graphs. For generated scientific results, regenerate the tracked report/result
first; the memory summarizes that published bundle and a few meaningful entry
points rather than hashing every helper used to build it.

## Reading it

- Prefer source pages marked current and newer explicit revisions over older
  material.
- Preserve older useful information as superseded rather than silently
  deleting it.
- Treat unresolved disagreements as conflicts; do not manufacture consensus.
- A script catalog says what code does and emits. A script printing `PASS`
  does not by itself establish a physical conclusion.
- Follow citations into `research/` or `software/` before relying on a detail.

The short conventions are in [`_meta/schema.md`](_meta/schema.md). The selected
source bundles and page dependencies are in [`catalog.json`](catalog.json).
The remaining Atlas migration checklist is in [`MIGRATION.md`](MIGRATION.md).
