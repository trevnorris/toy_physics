# Moving-Throat PDE Derivation Companion

This project now has two paper entry points built from the same theorem-block
source:

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
- `parts/`: eight main body parts.
- `stages/`: generated or archival stage templates; not the active narrative source for either build.
- `appendices/`: archive appendices plus compact reader-only summaries.

## Current source rule

The canonical stage narrative currently lives in
`appendices/stage_appendix_part*.tex`. The `stages/` directory should be treated
as template inventory unless and until an explicit generation pipeline is put
back in charge of the appendix text.

Generated stage count: 253.
