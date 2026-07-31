# archive/ — superseded apparatus, preserved not deleted

Moved here **2026-07-30** by `git mv`, so history is intact and `git log --follow` works.

⛔ **Nothing here is live.** ⛔ Do not build on it, cite it as current, or restore it without a decision.
⭐ **It was archived because the METHOD changed**, not because it was wrong — see
`docs/derivation_walkthrough_plan.md` §0 and `research/pde_ledger_v2/walkthrough/DECISIONS.md`.

## `census/` — the DERIVED-vs-DECLARED census

A backward audit over the finished corpus: a provenance schema, a pilot on stages 016/023/043, and two
blocking review rounds. It worked and it was retired because it made the physics unfollowable for its one
reviewer — after eleven commits, no physics had been verified.

⭐ **Its findings are NOT here — they were extracted before archiving** and live in
`research/pde_ledger_v2/_scratch/NEXT_SESSION.md`:
- **zero cross-artifact citations resolve to a locus** (measured, both pilot stages, both engines);
- a **false provenance attribution** — stage016's engines claim `M̃`/`K̃`/`T̃_Ω` are
  `CONSUMED-from-011/012/013`; those stages contain none of them;
- a **live contradiction between tracked documents** over whether `K0c`/`K_eta`/`T_Omega` were ever
  identified — a tier-1-vs-tier-3 disagreement;
- a **wrong locus in four tracked files** (`:355-366` should be `:314-325`).

## `manifests/` — the 44-manifest integration-test system

4 of 44 manifests were ever extracted. Superseded as a route by the walkthrough.
Contains `composite_build.py`, both schema versions, `MANIFEST_README.md`, `EXTRACTION_PROTOCOL.md`,
`LEDGER_WIDE_PLAN.md`, the mutators, examples, reports, and the 4 extracted stage manifests.

⭐⭐ **Two files were deliberately LEFT BEHIND at `research/pde_ledger_v2/manifests/`:**
- **`DIMENSION_REWRITE.md`** — the active conversion workstream's canonical doc, cited by many files;
- **`DIM_ORDER_DECISION.md`** — despite the manifest-sounding title, its content is a retraction about
  *audit-script* dimension order (`[L,T,M]` vs `(L,M,T)`), which is live dimension-rewrite material.

⚠ Their paths are **unchanged on purpose**. Moving them would break every citation.

## ⚠ One known dangling reference

`research/pde_ledger_v2/schemas/validate_dimension_survey.py` pattern-matches locus citations of the form
`composite_build.py:N` (`:121`, `:1678`). It does **not import** the file, so nothing breaks at runtime —
but any citation it validates now points into `archive/`. ⛔ Left as-is rather than edited, because that
validator's own live/retired status is undecided.
