# v3 ledger — paper

⛔ **v3 is SELF-CONTAINED.** `document_setup.tex` and `macros.tex` are **copies** of v2's, not links.
Do not `\input` anything from `research/pde_ledger_v2/paper/`.

## Compile

```bash
pdflatex -interaction=nonstopmode pde_ledger_v3.tex   # twice, for the TOC
```

`latexmk` is not installed on this machine; two `pdflatex` passes are the working recipe.

## Structure — and the ordering is the METHOD, not a presentation choice

⭐⭐ **Parts run REQUIREMENTS-FIRST: each force sector first, the substrate LAST.**
Each sector states what it *needs*; the brane and bulk are defined at the **knit** (Part V) by asking
whether one medium supplies every requirement at once. ⛔ **Do not reorder to put the medium first** —
that is the direction v3 explicitly moved away from (`CHARTER.md` §0).

| | |
|---|---|
| `parts/part01_light` … `part04_charge_magnetism` | **half one** — the sectors' requirements |
| `parts/part05_the_knit` | **half one's verdict** — can one medium satisfy them all? A no-go here IS the falsification |
| `parts/part06_simulation_handoff` | **half two** — what only a simulation can settle. ⛔ Nothing requiring a sim is DONE here |
| `steps/` | one card per step, named by **v3 step ID** (⛔ not v2's `stage_NNN`) |
| `appendices/registry_provenance` | ⭐ where **S0.5** lives — it is *bookkeeping repair, not physics*, so it gets no card in a physics part |
| `appendices/source_map` | per-step input loci |

## Status

**One card exists: `steps/S9_light_requires_shear.tex`.** Everything else is a stub that fixes the
ordering before content is written into it. ⚠ Stubs say *"Not yet walked"* — ⛔ do not read a stub as an
empty result.

## ⚠ The built PDF is NOT tracked

Deliberate, and revisit when the ledger is near done: a 19-page skeleton regenerated at every step is
~300 KB of binary churn per commit with no review value. ⚠ v2 *does* track `pde_ledger.pdf`, so this is
a departure from that precedent, made knowingly.
