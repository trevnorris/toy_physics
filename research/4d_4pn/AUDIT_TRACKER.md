# 4d_4pn Audit Tracker

This file tracks the `research/AUDIT_WORKFLOW.md` closeout for the `4d_4pn` paper
bundle.

## Bundle Inventory

- `paper/4d_4pn.tex`
- `scripts/*.py`:
  - `4pn_onebody_audit.py`
  - `4pn_quartic_legendre_audit.py`
  - `4pn_local_scaffold_audit.py`
  - `4pn_local_target_import_audit.py`
  - `4pn_local_hamiltonian_to_ordinary_audit.py`
  - `4pn_hamiltonian_chart_generic_frame_lift_audit.py`
  - `4pn_hamiltonian_chart_canonical_slice_audit.py`
  - `4pn_generic_frame_ordinary_translation_audit.py`
  - `4pn_generic_frame_aligned_seed_lift_audit.py`
  - `4pn_local_referee_master_sympy_audit.py`
  - `4pn_tail_bridge_audit.py`
  - `4pn_tail_hereditary_bridge_audit.py`
  - `4pn_conditional_referee_master_audit.py`
- `notes/*.md`: migrated stage notes plus the consolidated `4d_4pn_full_notes.md`
- `mathematica/*.wl`: grouped dual-CAS mirror of the referee-facing theorem chain

## Status

- [x] Existing SymPy and note archive migrated into `research/4d_4pn/`
- [x] Local scripts README added
- [x] Mathematica manifest / README / output-block updater added
- [x] All SymPy scripts rerun from bundle-local paths
- [x] Mathematica mirror created
- [x] All Mathematica scripts rerun and embedded `Output:` blocks normalized
- [x] Paper verification section updated to match actual archive layout
- [x] Paper rebuilt cleanly on two passes

## Remaining Non-Bundle Duplicates

- There is still a duplicate copy of `notes/4d_4pn_full_notes.md` at repo root-level
  `notes/4d_4pn_full_notes.md`; the paper bundle now carries the migrated version in
  `research/4d_4pn/notes/4d_4pn_full_notes.md`.

## Current Archive Layout

- `scripts/*.py`: full stagewise SymPy referee archive for the local and hereditary chain
- `mathematica/*.wl`: grouped Mathematica mirror for the one-body gate, quartic
  compiler, local scaffold/target import, Hamiltonian chain, ordinary chain, tail bridge,
  and conditional full-4PN replay
- `mathematica/wl_notes.txt`: official WL manifest used by `update_output_blocks.py`
