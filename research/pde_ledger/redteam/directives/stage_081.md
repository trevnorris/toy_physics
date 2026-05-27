---
unit_id: 081
batch: III.4
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 081 (resolved orchestrator-direct)

The single F1 `paper_misalignment` (banner relabel "Stage 64" -> "Stage 081" on the SymPy script) was resolved as direction (b) on 2026-05-27 (script-side relabel only; legacy numbering in `.wl` was already correct) and applied by the orchestrator directly:

- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:3` docstring `Stage 64` -> `Stage 081`
- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:28` banner `STAGE 64` -> `STAGE 081`

No Codex invocation required. Verifier will confirm the new banner appears in transcripts and the script exits 0.

Codex: do nothing here. Do not edit `paper.tex`, `notes/`, or the scripts to "fix" this paper_misalignment unless a follow-up directive (after the user chooses a direction) explicitly authorises it.

## F1 — paper_misalignment

**Subtype:** target_mismatch (script self-identifies as the wrong stage)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_081.tex:1` quote: `\section[Stage 081]{Stage 081: Family--1 Demand Window in the Branch-Product Variable}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage081_family1_pi_thresholds.md:1` quote: `# Moving-Throat PDE — Stage 081: Final Explicit Family-1 Quadrupole-Demand Window in the Branch-Product Variable `Pi_tr``
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex:140` quote: `081 & Branch-product demand window & \StatusExactClosure{} & Thresholds in \(\Pi_{\rm tr}\) with blocking parameter.`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:38` quote: `banner["STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS"];`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:3` quote: `"""SymPy audit for Stage 64.`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:28` quote: `banner("STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS")`
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt:3` quote: `STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS` (echo of banner)

The SymPy script's filename, location, all of its math content (Stage-35 inversion, Stage-63 zeta thresholds, branch-product `Q`), and its companion Mathematica script all correspond to Stage 081. Only the docstring at line 3 and banner at line 28 carry the obsolete label "Stage 64". This is a leftover from commit `0d09ef6 fully reorder the pde ledger`, in which stages were renumbered (the pre-renumber identifier was apparently "Stage 64").

## Resolve before fix_loop

The SymPy script `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py` self-identifies as "Stage 64" in its docstring (line 3) and banner (line 28), but the paper card, notes, appendix row, and companion Mathematica script all call this Stage 081. The script's math content is correct and aligned with the Stage 081 paper card; only the labels disagree.

Which direction should be applied?

- (a) Script is stale (almost certainly the case — global renumbering committed in `0d09ef6` updated the filename but not the in-file string literals) -> change SymPy script line 3 from `"""SymPy audit for Stage 64.` to `"""SymPy audit for Stage 081.` and line 28 from `banner("STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS")` to `banner("STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS")`. Re-run `redteam exec-sympy 081` to refresh the transcript. No paper-side edit.
- (b) Paper is stale (extremely unlikely — Stage 64 already exists as a different stage post-renumber) -> would require renaming the file, the Mathematica script, and the notes file, and updating the paper card heading; flagged as inconsistent with the rest of the ledger.
- (c) Both labels refer to the same stage under different naming conventions (informational rename only) -> still requires choosing one canonical label across the ledger, and the rest of the ledger already standardises on "Stage 081".

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. Recommended direction: (a).
