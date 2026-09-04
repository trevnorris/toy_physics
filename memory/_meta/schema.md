# Memory conventions

Repository files under `research/` and `software/` are authoritative. Memory
pages summarize and navigate them.

## Page frontmatter

Use only fields that help a reader:

```yaml
---
title: Moving-throat dynamics
type: topic          # source | topic | script_catalog | conflict_register | index
status: current      # current | superseded | conflict | draft
sources:
  - research/pde/paper/pde.tex
last_updated: 2026-09-03
superseded_by:       # optional path
---
```

Do not put Git commit IDs, Git blob IDs, transaction IDs, model attestations,
or tool-version hashes in pages. Content hashes belong only in
`_meta/state.json`.

## Source capsules

A source capsule should answer:

- What is this source trying to do?
- What are its main claims or results?
- What assumptions and limitations matter?
- What does it revise or supersede?
- What remains open?

Important statements cite an original path plus a TeX label, Markdown heading,
or nearby code identifier when practical.

## Topic pages

A topic page gives the current working picture across a few source capsules.
It names disagreements and explains why one treatment is preferred. Prefer an
explicit correction or supersession statement; use dates only when the sources
otherwise have comparable roles and scopes.

## Script catalogs

For each meaningful entry point, record its path, role, inputs, outputs, and
what it checks. Do not enumerate every helper. Do not turn a successful program
exit into a scientific conclusion.

## Conflicts

Keep both sides, cite both original sources, and mark the conflict open or
resolved. Never silently erase the older position.

## Review

AI-written pages are drafts until read. A lightweight Grok check is enough for
important synthesis in this toy-model repository; no stored review packet or
attestation is required.
