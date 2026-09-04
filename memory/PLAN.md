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

## Initial buildout (completed)

- [x] Populate the selected paper, software, and proposal capsules.
- [x] Populate the eight topic pages and four script catalogs.
- [x] Migrate the useful Atlas status firewalls with original-source citations.
- [x] Refresh `conflicts.md` and make `index.md` the entry point.
- [x] Pass independent retrieval checks using only the memory.
- [x] Retire `atlas/` and `graph/` after completing the migration checklist.

Ongoing work is source-driven maintenance: when `status` reports a changed
unit, reread that unit and refresh only its affected pages.

That is the entire system. If maintaining the helper ever becomes harder than
reading a few changed files, simplify it again.
