# EM projected integration progress

Created: 2026-05-11

## Working scope

Integrate the derivation-only portion of `notes/em_projected` into the PDE
ledger. Use steps 1--18 as source material. Exclude the step 19 branch-export
diagnostic and step 20+ computational, runtime, frontier, and falsification
work from the paper derivation.

## Implementation checklist

- [x] Create ordered Stage 004 substage TeX files.
- [x] Rework `paper/stages/stage_004.tex` as the public wrapper.
- [x] Add matching ordered SymPy audits under `research/pde_ledger/scripts/`.
- [x] Add matching Mathematica mirrors under `research/pde_ledger/mathematica/`.
- [x] Update Part I prose.
- [x] Update result-anchor map and dependency map.
- [x] Update Stage Appendix Part I.
- [x] Update provenance and verification coverage notes.
- [x] Update source-file index.
- [x] Run targeted new SymPy audits.
- [x] Run targeted new Mathematica audits if `math` is available.
- [x] Build the LaTeX paper.

## Notes

The public global Stage 004 number remains stable. Projection-first Maxwell is
inserted as ordered Stage 004 substages using `stage_004_001_*` names so local
filesystem ordering matches derivation order.  The root audit files keep the
same `stage004_00N` order and basename intent as the TeX substages; no hidden
`step_*` source subfolder is used.

## Verification completed

- `bash research/pde_ledger/scripts/run_all_audits.sh stage004_0 --force`
  passed all 18 ordered Stage 004 SymPy audits after the derivation-only
  boundary was tightened.
- File-for-file Mathematica mirrors for the imported EM-projected scripts are
  not yet claimed.  The only retained Mathematica substage mirror is the legacy
  reduced one-port adapter, and
  `bash research/pde_ledger/mathematica/run_all_audits.sh stage004_018 --force`
  passed that adapter.
- `pdflatex -interaction=nonstopmode -halt-on-error pde_ledger.tex` completed
  cleanly after the normal cross-reference reruns.
