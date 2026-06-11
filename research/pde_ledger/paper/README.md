# Moving-Throat PDE Derivation Companion

This project now has two paper entry points built from the same paper source
hierarchy:

- `pde_ledger.tex`: the authoritative archive ledger.
- `pde_ledger_reader.tex`: the collaborator-facing reader ledger.

## Compile

From this directory:

```bash
latexmk -pdf pde_ledger.tex
latexmk -pdf pde_ledger_reader.tex
```

or, if `latexmk` is unavailable:

```bash
pdflatex pde_ledger.tex
pdflatex pde_ledger.tex
pdflatex pde_ledger_reader.tex
pdflatex pde_ledger_reader.tex
```

## Build roles

- The archive build keeps the full frontmatter firewalls, stage ledger, stage
  appendices, reproducibility map, source-file index, and fill workflow.
- The reader build keeps short orientation frontmatter, the same Parts I--VIII,
  and compact provenance/verification appendices.

## Included structure

- `pde_ledger.tex`: archive entry point.
- `pde_ledger_reader.tex`: reader entry point.
- `document_setup.tex`: shared package and macro setup.
- `main_parts.tex`: shared Part I--VIII include list.
- `macros.tex`: theorem/status/stage macros.
- `frontmatter/`: archive governance chapters plus reader-only orientation chapters.
- `parts/`: eight main body parts; this is the shared theorem-block
  synthesis layer.
- `stages/`: canonical per-stage stage cards, provenance anchors, and audit
  anchors.
- `appendices/`: archive appendices plus compact reader-only summaries.  The
  archive stage appendices carry inline theorem-block synthesis narrative and
  `\input` all 253 stage cards from `stages/`.

## Current source rule

The active paper source is two-layered.  The theorem-block synthesis layer
lives inline in `parts/` and in `appendices/stage_appendix_part*.tex`; the stage
appendices are not mere input wrappers.  The stage-card layer lives in
`stages/stage_NNN.tex`; those files are the canonical per-stage provenance and
audit anchors.  All 253 stage cards are `\input` into the archive stage
appendices.

For provenance disputes beneath either paper layer, escalate to the
PDE-ledger-relative per-stage notes and audit record, especially
`notes/stages/` and its review/per-stage notes subtrees, plus the associated
script and notebook evidence.

Generated stage count: 253.
