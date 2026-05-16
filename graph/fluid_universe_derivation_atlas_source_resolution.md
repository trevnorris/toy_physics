# Fluid Universe Derivation Atlas - Source Resolution

This file records the cleanup from generated, bare source filenames to repository
paths. Existing maintained TeX papers are treated as canonical. Markdown-only
derivation ledgers remain linked deliberately so they mark future papers that
still need to be written.

## Canonical Paper Targets

| Legacy source | Canonical target | Source note retained |
|---|---|---|
| `4d_summary.md`, `4d.tex` | `research/4d/paper/4d.tex` | `notes/summaries/4d_summary.md` |
| `4d_em_fields_summary.md`, `4d_em_fields.tex` | `research/4d_em_fields/paper/4d_em_fields.tex` | `notes/summaries/4d_em_fields_summary.md` |
| `4d_plasma_summary.md`, `4d_plasma.tex` | `research/4d_plasma/paper/4d_plasma.tex` | `notes/summaries/4d_plasma_summary.md` |
| `4d_1pn_bridge_summary.md`, `4d_1pn_bridge.tex` | `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` | `notes/summaries/4d_1pn_bridge_summary.md` |
| `4d_1pn_full_summary.md`, `4d_1pn_full.tex` | `research/4d_1pn_full/paper/4d_1pn_full.tex` | `notes/summaries/4d_1pn_full_summary.md` |
| `4d_2pn_summary.md`, `4d_2pn.tex` | `research/4d_2pn/paper/4d_2pn.tex` | `notes/summaries/4d_2pn_summary.md` |
| `4d_2_5pn_summary.md`, `4d_2_5pn.tex` | `research/4d_2_5pn/paper/4d_2_5pn.tex` | `notes/summaries/4d_2_5pn_summary.md` |
| `4d_3pn_summary.md`, `4d_3pn.tex` | `research/4d_3pn/paper/4d_3pn.tex` | `notes/summaries/4d_3pn_summary.md` |
| `4d_4pn_summary.md`, `4d_4pn.tex`, `4d_4pn_full_notes.md` | `research/4d_4pn/paper/4d_4pn.tex` | `notes/summaries/4d_4pn_summary.md`, `research/4d_4pn/notes/4d_4pn_full_notes.md` |
| `pde_ledger.tex` | `research/pde_ledger/paper/pde_ledger.tex` | `notes/moving_throat_pde_program_compact.md` |

## Paper Components

The moving-throat stage notes have maintained TeX components in the PDE ledger:

| Legacy source | Canonical target |
|---|---|
| `moving_throat_pde_stage001_geometry_lift.md` | `research/pde_ledger/paper/stages/stage_001.tex` |
| `moving_throat_pde_stage002_breathing_reduction.md` | `research/pde_ledger/paper/stages/stage_002.tex` |
| `moving_throat_pde_stage003_bdg_coupling.md` | `research/pde_ledger/paper/stages/stage_003.tex` |
| `notes/em_projected/step_01_projected_maxwell_readme.md` through `notes/em_projected/step_18_parent_throat_action_weak_axisym_packet_notes.md` | `research/pde_ledger/paper/stages/stage_004.tex` through `research/pde_ledger/paper/stages/stage_020.tex` |
| `moving_throat_pde_stage021_reduced_one_port_normal_form.md` | `research/pde_ledger/paper/stages/stage_021.tex` |
| `moving_throat_pde_stage022_grouped_p2_normalization_bridge.md` | `research/pde_ledger/paper/stages/stage_022.tex` |
| `moving_throat_pde_stage023_full_grouped_bundle.md` | `research/pde_ledger/paper/stages/stage_023.tex` |
| `moving_throat_pde_stage024_overlap_isotropy.md` | `research/pde_ledger/paper/stages/stage_024.tex` |
| `moving_throat_pde_stage025_minimal_isotropic_normalization.md` | `research/pde_ledger/paper/stages/stage_025.tex` |

## Future Paper Backlog

These sources intentionally remain markdown because no maintained TeX paper target
exists yet.

| Current source | Future-paper note |
|---|---|
| `notes/5pn/5pn_notes_full.md` | Future 5PN paper needed. The legacy `4d_5pn_summary.md` filename was not present in this checkout. |
| `notes/moving_throat_pde_program_compact.md` | Future compact moving-throat/theorem-status paper or maintained paper section needed. |
| `notes/pde_audit_full.md` | Future PDE-audit paper or maintained ledger target needed. |
| `notes/atom_work.md` | Future atomic reduced-sector paper needed. |
| `notes/lepton_work.md` | Future lepton conditional-theorem paper needed. |
| `notes/lepton_mass_notes.md` | Future lepton mass paper needed. |
| `notes/g2/g2_full_output.md` | Future g-2/anomaly paper needed. |
| `notes/moving_throat_notes_full.md` | Future full moving-throat derivation paper or maintained ledger target needed. The legacy `moving_throat_output_full.md` filename was not present. |

## Provenance Rules

- `source_kind: paper` means the graph should link to the TeX paper.
- `source_kind: paper_component` means the graph should link to a TeX component included by a maintained paper.
- `source_kind: future_paper_note` means the markdown file is the current source and should be replaced only after a paper exists.
- `legacy_*` fields are retained as provenance for generated anchors and should not be used as canonical targets.
