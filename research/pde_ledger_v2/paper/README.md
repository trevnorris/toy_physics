# PDE Ledger Archive Skeleton

This directory is the empty LaTeX skeleton for the rebuilt ledger.

## Compile

From this directory:

```bash
latexmk -pdf pde_ledger.tex
```

or, if `latexmk` is unavailable:

```bash
pdflatex pde_ledger.tex
pdflatex pde_ledger.tex
```

## Included Structure

- `pde_ledger.tex`: archive entry point.
- `document_setup.tex`: shared package setup.
- `main_parts.tex`: placeholder for future Part includes.
- `macros.tex`: theorem/status/stage macros.
- `frontmatter/`: placeholder governance chapters.
- `parts/`: empty directory for future theorem-block chapters.
- `stages/`: empty directory for future canonical stage cards.
- `appendices/`: infrastructure appendices reset to zero-stage placeholders.

Generated stage count: pending.
