# EM projected integration progress

> Supersession note, 2026-05-11: this tracker originally recorded the
> temporary local Stage 004 substage import.  The current canonical ledger no
> longer uses that scheme.  The projected Maxwell derivation now occupies
> linear Stages `004--021`, and the old post-EM ledger starts at Stage `022`.

Created: 2026-05-11

## Working scope

Integrate the derivation-only portion of `notes/em_projected` into the PDE
ledger. Use the derivation notes through step 18 as source material; the source
directory has no step 06 file.  Retain the reduced one-port normal form as the
Stage 021 closure after the projected EM packet.  Exclude the step 19
branch-export diagnostic and step 20+ computational, runtime, frontier, and
falsification work from the paper derivation.

## Implementation checklist

- [x] Promote projected Maxwell derivation-only notes to linear
      `paper/stages/stage_004.tex` through `paper/stages/stage_020.tex`.
- [x] Retain the reduced one-port normal form as `paper/stages/stage_021.tex`.
- [x] Shift the previous post-EM ledger to `stage_022.tex` through
      `stage_253.tex`.
- [x] Add matching ordered SymPy audits under `research/pde_ledger/scripts/`.
- [x] Retain the reduced one-port Mathematica mirror as Stage `021`; shift the
      former post-EM Mathematica scripts to the new linear order.
- [x] Update Part I prose.
- [x] Update result-anchor map and dependency map.
- [x] Update Stage Appendix Part I.
- [x] Update provenance and verification coverage notes.
- [x] Update source-file index.
- [x] Run targeted new SymPy audits.
- [x] Run targeted new Mathematica audits if `math` is available.
- [x] Build the LaTeX paper.

## Notes

Projection-first Maxwell is now a normal stage sequence rather than a local
Stage 004 subsection.  Filesystem order, paper order, SymPy script order, and
transcript order all use the same global numbers.  The imported `notes/em_projected`
portion ends at Stage `020`, Stage `021` supplies the retained one-port closure,
and the step 19 branch-export diagnostic and step 20+ computational work remain
excluded from the paper derivation.

## Verification completed

- Full SymPy transcript regeneration completed after the linear renumbering:
  `TOTAL: 241  PASS: 241  FAIL: 0  SKIPPED: 0`, including the repo-level
  master audit and `240` stage-level audits.
- Mathematica transcript regeneration completed for the current Mathematica
  coverage set: `TOTAL: 165  PASS: 165  FAIL: 0`.
- File-for-file Mathematica mirrors for projected Maxwell Stages `004--020` are
  still not claimed.  The retained Mathematica mirror is the reduced one-port
  normal form at Stage `021`, plus the shifted post-EM Mathematica audits.
- Numerical stress summaries are current for the Python and Mathematica
  numerical harness sets.
- The LaTeX paper was rebuilt after transcript regeneration and stale-path
  scans in the post-renumbering verification pass.
