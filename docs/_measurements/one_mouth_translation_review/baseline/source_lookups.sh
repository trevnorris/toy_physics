#!/bin/bash
set -ex
sha256sum _scratch/s11c/one_mouth_translation_review01/r0/section.md
sed -n '149,161p' research/pde_ledger_v3/steps/S16_moving_region.md
sed -n '2207,2213p;2219,2221p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
sed -n '444,449p;456,474p;735,739p' docs/drain_driven_initial_acceleration.md
sed -n '9,31p' software/stage1_solver/src/stage1_solver/mouth_probe.py
sed -n '88,93p;100,112p' docs/held_drain_pair.md
rg -n '^Positive speeds|^per_speed|^U=|^abs_U=' docs/_measurements/held_crossflow_budget.txt
rg -n '^  "status"|^  "scope"' docs/_measurements/mouth_release_onset.json
rg -n '^held acceleration|^support-removal residual' research/pde_ledger_v3/_measurements/mouth_release_review/local_witness.stdout.txt
