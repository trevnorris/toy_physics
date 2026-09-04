---
title: Research conflicts
type: conflict_register
status: current
sources:
  - research/4d/paper/4d.tex
  - research/4d_1pn_full/paper/4d_1pn_full.tex
  - research/4d_2pn/paper/4d_2pn.tex
  - research/4d_2_5pn/paper/4d_2_5pn.tex
  - research/4d_4pn/paper/4d_4pn.tex
  - research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex
  - research/brane_bulk_ontology/paper/brane_bulk_ontology.tex
  - research/em_fields/paper/em_fields.tex
  - research/pde/paper/pde.tex
  - research/pde_audit/notes/stage_v2_26_program_status_after_audit.md
  - research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md
  - software/em_charge_attribute/puncture_deflection_electric_sign_result.md
  - software/em_charge_attribute/reports/native_p_constraint_gate.md
  - software/stage1_solver/README.md
  - software/stage1_solver/STAGE1_VERDICT.md
last_updated: 2026-09-03
---

# Conflicts

This register keeps genuine corrections and unresolved choices visible. A
difference in scope or claim status is not automatically a contradiction.

## Resolved or scoped corrections

| Issue | Current reading | Earlier reading retained for context |
|---|---|---|
| Electric-charge ontology | Use the oriented branch sign $\eta_Q$, microscopic charge $q_\star=\eta_Qe_\star$, and normalized brane coupling $q_{\rm eff}=q_\star/\sqrt{Z_{\rm int}}$. Circulation is magnetic/vortical; breathing does not change charge. | `research/em_fields/paper/em_fields.tex` (`eq:q-def`, `subsec:gauss-breathing`, `eq:q-from-volume`) and `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` (`subsec:charge_flux`) used circulation, geometry, or breathing as electric-charge motivation. Those charge definitions are superseded in scope; unrelated formulas remain usable. Current sources: `research/4d/paper/4d.tex` (`sec:fields`, `eq:qeff-canonical`) and `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` (`Electric charge variables`). |
| Isotropic wake mixing | The full isotropic projector gives a real calibration $\alpha^2=3/4$. | An earlier truncated contraction suggested an imaginary or negative solution. `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` (`subsec:nbody-derivation`) records the corrected full contraction. |
| Stage-1 outcome | The frozen effective-closure branch was run and robustly missed, with $R_{\rm norm}\approx-10.8$. This rejects that branch under its frozen assumptions, not all of Path A. | `software/stage1_solver/README.md` predates the run and says no physical branch had been solved. Use it for architecture; use `software/stage1_solver/STAGE1_VERDICT.md` for the later outcome. |

## Open choices and unresolved gates

| Status | Issue | Current reading and sources |
|---|---|---|
| Open | Electric boundary class | The sign is not unique before the mouth ensemble is selected: held value is positive, fixed source is negative, fixed monopole is indefinite, and mixed conditions span outcomes. Current landing: `R1_REQUIRED(bc_selection)`. Source: `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` (`Q-BC`, `Q-FORCE / Q-COMBINE`, and `§4 landing`). |
| Open | Origin of electromagnetism | Localized Maxwell theory is assumed in the parent action. Native $P^a$ has no generic quadratic Gauss constraint in the tested regular stratum; the added spin-ice route is an identification, not a dynamical construction of the throat carrier. Sources: `research/4d/paper/4d.tex` (`sec:fields`) and `software/em_charge_attribute/reports/native_p_constraint_gate.md` (`Decision-table result`). |
| Open | Moving-throat parent | The strict parent does not provide autonomous wall dynamics. The linear $S_\eta^{(2)}$ closure is consistent, while a nonlinear $S_\Sigma[R]$, material law, source/port law, and same-branch exporter remain unselected. Sources: `research/pde/paper/pde.tex` (`sec:status-parent`, `sec:linearized-pde`) and `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` (`What the audit says is missing`). |
| Open | Quadrupole normalization | The angular/STF source map is fixed, but the radial/axial/source-to-outgoing normalization and transport are not. The same open effective quadrupole datum controls 2.5PN reaction and the 4PN hereditary tail. Sources: `research/4d_2_5pn/paper/4d_2_5pn.tex` (`sec:conditional-theorem-gap-normalization`) and `research/4d_4pn/paper/4d_4pn.tex` (`eq:4pn-tail-exact-bridge`, `sec:tail-no-new-datum`). |
| Conditional, not contradictory | PN matches versus branch misses | The PN results are algebraic matches after declared closures; branch tests ask whether a particular frozen PDE/effective branch realizes the needed data. A miss can coexist with correct conditional algebra. Sources: `research/4d_1pn_full/paper/4d_1pn_full.tex` (`sec:taxonomy`), `research/4d_2pn/paper/4d_2pn.tex` (`sec:status-open`), and `software/stage1_solver/STAGE1_VERDICT.md` (`Interpretation — what it means (and does NOT mean)`). |
| Conditional, not contradictory | “No new 4PN datum” versus open tail transport | The canonical 4PN bridge needs no second quadrupole-normalization constant; it inherits the open 2.5PN normalization. Physical transport and branch realization are still open. Source: `research/4d_4pn/paper/4d_4pn.tex` (`sec:tail-no-new-datum`, `sec:fixed-open-still-open`). |

See [Status and reading rules](topics/status-and-reading-rules.md) before
promoting any conditional match or program output to a physical claim.
