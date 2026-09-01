# Scanner Decision Review Fold

Date: 2026-08-24
Artifact reviewed: `memory/_meta/SCANNER_DECISIONS.md` (author draft)

Two independent review legs were completed before builder launch. Neither leg
edited the decision list or implementation.

## Leg 1 — corpus and scanner audit

Verified blocking findings:

- Git glob pathspec assumptions were unusable on the installed Git.
- Config was descriptive rather than executable and had overlapping v3 units.
- Blob identities did not ensure the AI read committed rather than dirty bytes.
- Atomic state did not protect live pages from partial publication.
- Partial updates could falsely advance a global commit baseline.
- Policy bodies were not digest-pinned.
- Deletion and corpus retirement were conflated.
- Large tracked evidence lacked read bounds; important `.out` files were
  omitted.
- Page review status and scientific source authority were conflated.
- Page-level evidence could promote mixed-content pages.
- YAML duplicate keys, snapshot-specific lint, unique targets, and anchor
  validation needed explicit controls.

The leg also verified that the sixteen selected monolithic paper entry points
have no component `\input` trees, nested `.git` directories are not candidate
submodules, and tracked `.out` evidence must remain eligible.

## Leg 2 — independent scanner review

Independently verified:

- Source units needed normalized IDs, roles, targets, selectors, and shared
  dependency declarations.
- Prepare must materialize immutable target-tree blobs and pin mode/type/size.
- Generated pages needed staging plus journaled publication and recovery.
- Per-unit processed commits were required for partial updates.
- Large evidence needed identity-only or deterministic excerpt treatment.
- Selector semantics, deletion/tombstone transitions, duplicate-key rejection,
  anchor lint, and conservative stale-query warnings were underspecified.

## Fold disposition

All verified blockers were accepted. The folded decision list now requires:

- one committed-tree enumeration and a local prefix/suffix selector dialect;
- normalized executable units;
- immutable source snapshots for AI extraction;
- semantic/identity/excerpt read modes and byte limits;
- config/schema/prompt/tool digests;
- staged pages, publication journal, rollback, and crash recovery;
- per-unit and global freshness separation;
- explicit source deletion versus corpus retirement;
- snapshot-aware anchor/link lint and served-page hash checks.

One recommendation was not adopted: adding root `STATUS.md` and `CLAUDE.md` as
normal memory source inputs. They remain orchestration controls outside the
user-approved semantic corpus boundary of `research/` and `software/`.
