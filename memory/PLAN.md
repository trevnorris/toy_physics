# Research memory plan

The goal is a useful LLM-facing wiki for this toy model, not a publication or
provenance system.

## Keep

- concise source capsules for selected papers and maintained software results;
- a small set of cross-source topic pages;
- short script catalogs;
- an explicit conflict page;
- original source-path citations;
- content hashes so unchanged material is not reread.

## Do not rebuild

- a graph database;
- atomic claim records for every statement;
- commit-bound snapshots or transaction archives;
- model attestations, isolated review packets, or publication gates;
- multiple overlapping pages for the same topic;
- catalogs of root `notes/`, root `docs/`, or any PDE ledger.

## Workflow

1. `python3 memory/update.py status` compares current tracked source contents
   with `_meta/state.json`.
2. Codex reads only the changed source units and writes the affected Markdown
   pages directly.
3. Grok performs a lightweight second review when the prose resolves conflicts
   or makes important scientific distinctions.
4. Codex corrects the pages, runs `memory/update.py lint`, marks the completed
   units/pages, and commits the normal Markdown and state changes.

There is no stored Git SHA. Git only supplies the tracked-file list. A rebase,
merge, or unrelated commit does not make memory stale; changing source bytes
does.

## Remaining work

1. Finish the seven unpopulated software/source capsules.
2. Populate the eight topic pages and four script catalogs.
3. Fold the useful Atlas status firewalls into those topic pages using
   [`MIGRATION.md`](MIGRATION.md).
4. Refresh `conflicts.md` and make `index.md` a useful entry point.
5. Ask a fresh agent a handful of real retrieval questions.
6. Delete `atlas/` and `graph/` only after every migration checklist item is
   either represented or deliberately discarded.

That is the entire system. If maintaining the helper ever becomes harder than
reading a few changed files, simplify it again.
