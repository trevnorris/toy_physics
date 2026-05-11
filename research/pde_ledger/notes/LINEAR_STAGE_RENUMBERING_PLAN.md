# Linear stage renumbering tracker

Created: 2026-05-11
Updated: 2026-05-11

## Goal

Make the PDE ledger linear again after importing the derivation-only projected
Maxwell work from `notes/em_projected`.  The canonical paper, scripts, and
verification references should list in one global order:

- `stage_001`--`stage_003`: unchanged front-end geometry, wall, and BdG stages.
- `stage_004`--`stage_020`: projected Maxwell derivation imported from
  `notes/em_projected` derivation notes through step 18.  The source directory
  intentionally has no step 06 file.
- `stage_021`: the retained reduced one-port normal form that closes the
  projected EM packet into the grouped-response bridge.
- `stage_022`--`stage_253`: the previous post-EM ledger, shifted forward by
  17 stage numbers.

No reader-facing paper path should use the old local Stage 004 substage scheme.
The computational and branch-search work from `notes/em_projected` step 19 and
later remains outside the paper derivation.

## Current numbering policy

| Current stage range | Source |
|---|---|
| `001--003` | Existing PDE-ledger stages, unchanged. |
| `004--020` | Projected Maxwell derivation notes through step 18; no step 06 file exists. |
| `021` | Retained reduced one-port normal form. |
| `022--253` | Former stages `005--236`, shifted by `+17`. |

The total canonical stage count is now `253`.

## Projected Maxwell block

| Stage | Imported source |
|---:|---|
| `004` | Projected Maxwell bundle index / README. |
| `005` | Covariant projected Maxwell law. |
| `006` | Vector projected Maxwell form. |
| `007` | Projection/reduction comparison. |
| `008` | Projected Maxwell extension. |
| `009` | Near-throat projected Maxwell packet (`notes/em_projected` step 07). |
| `010` | Push-bundle master identities (`notes/em_projected` step 08). |
| `011` | Grouped `P_2` bridge (`notes/em_projected` step 09). |
| `012` | Projected Maxwell primitive bridge (`notes/em_projected` step 10). |
| `013` | Mouth-Taylor master identities (`notes/em_projected` step 11). |
| `014` | Mouth-Taylor gate bridge (`notes/em_projected` step 12). |
| `015` | Parent throat action master (`notes/em_projected` step 13). |
| `016` | Parent throat action candidate (`notes/em_projected` step 14). |
| `017` | Parent throat action weak-axisymmetric limit (`notes/em_projected` step 15). |
| `018` | Parent throat action bundle master (`notes/em_projected` step 16). |
| `019` | Parent throat action isotropic bundle (`notes/em_projected` step 17). |
| `020` | Parent throat action weak-axisymmetric packet (`notes/em_projected` step 18). |
| `021` | Reduced one-port normal form retained from the prior EM reduction stage. |

## Required consistency surfaces

- Paper stage files: `paper/stages/stage_001.tex` through
  `paper/stages/stage_253.tex`.
- Paper flow files: part introductions, result-anchor map, reader status
  summary, stage appendices, stage ledger, source-file index, dependency and
  reproducibility maps.
- SymPy audits and transcripts: top-level `scripts/moving_throat_pde_stageNNN_*`
  files and `scripts/output/*.txt`.
- Mathematica audits and transcripts:
  `mathematica/moving_throat_pde_stageNNN_*` files and
  `mathematica/output/*.txt`.
- Numerical stress artifacts whose names or prose mention old stage numbers.
- Metadata notes: citation map, provenance index, verification coverage,
  integration progress notes, and any active review/status summaries.

## Execution status

- Stage TeX files have been promoted to `stage_004.tex`--`stage_021.tex`, and
  former stages `005--236` have been shifted to `022--253`.
- Top-level SymPy audit filenames now follow the same order as the paper.
- The reduced one-port Mathematica audit is now Stage `021`; former
  Mathematica audits after the old EM stage have been shifted forward.
- Stage appendix Part I and the consolidated stage ledger list `001--253`
  linearly.
- Result-anchor ranges and Part I--IV handoff prose have been updated for the
  new stage numbers.
- Active metadata files are being updated to describe the linear numbering;
  historical integration notes should either carry a clear supersession banner
  or be rewritten before they are used as current guidance.

## Verification completed

1. SymPy transcripts regenerated:
   `TOTAL: 241  PASS: 241  FAIL: 0  SKIPPED: 0`, including the repo-level
   master audit and `240` stage-level audits.
2. Mathematica transcripts regenerated through the one-kernel-at-a-time runner:
   `TOTAL: 165  PASS: 165  FAIL: 0`.
3. Numerical stress transcripts are current for the Python and Mathematica
   numerical harness sets.
4. `\StageFile{...}` and stale-numbering scans were run against paper-facing
   files and active metadata.
5. `research/pde_ledger/paper/pde_ledger.tex` was rebuilt in the
   post-renumbering verification pass.
