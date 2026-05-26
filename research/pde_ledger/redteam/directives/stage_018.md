---
unit_id: 018
batch: I.2
created_at: 2026-05-25T18:06:25-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 018

This directive contains a single `paper_misalignment` finding. **Codex must apply nothing on this unit until the user resolves the question below.** The orchestrator will halt and surface the question; do not edit paper.tex, notes/, or the scripts to "fix" this finding speculatively.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_018.tex:26-28` quote: "Stage~018 exports the parent-wall bundle bridge \eqref{eq:stage018-msigma}--\eqref{eq:stage018-ksigma}."
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_018.tex:14-22` quote: the only body equations are `M_\Sigma=\int dw\,\mu_\eta\beta_2^2` and `K_\Sigma=\int dw\,[T_w(\beta_2')^2+(K_\eta+6T_\Omega)\beta_2^2]`.
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:58` quote: "018 & Parent throat action bundle master & \StatusExactClosure{} & Bundle-level parent-action identities used by the projected electromagnetic response."
- No `notes/stages/moving_throat_pde_stage018_*.md` file exists; the docstring-referenced `step_16_parent_throat_action_bundle_master_notes.md` is not present anywhere in the repository.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:1-2` quote: `"""Master-note audit for step_16_parent_throat_action_bundle_master_notes.md."""`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:25-88` — assertions covering one-pole numerator identity, two `KSigma` closures, compatibility cross-closure, even-gate determinant `= 1/27`, closed-form wall-stiffness slope `dKSigma = B01+Z01+27(B41+Z41)`, closed-form wall-inertia slope `dMSigma = -(B21+Z21)+3(B41+Z41)`, and residual amplitude `Xi1 = N01/N0 - 27(B41+Z41)/(KSigma-B0-Z0)`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:90-95` — the only script-side checks that directly map to the paper's bridge equations are the Gaussian-profile reductions `MSigma_example = sqrt(pi)` and `KSigma_example = 3*sqrt(pi)/2`, with `mu_eta = T_w = K_eta + 6 T_Omega = 1` and `beta = exp(-w^2/2)` hardcoded.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl:3-15` — header explicitly enumerates claim families M1–M8 (one-pole numerator, KSigma closures, compatibility, even-gate determinant, wall-slope solve, Xi1 residual, Gaussian integrals), mirroring the SymPy scope.

## Resolve before fix_loop

The script (both engines) verifies five families of identities, but the stage-018 paper card only declares the two `M_Σ` / `K_Σ` bridge integrals as `\stagefield{Output}`. The script's families 1–4 (one-pole numerator, KSigma closures + compatibility, even-gate determinant + wall-slope solve, Xi1 residual) have no anchor in either the .tex card or any existing notes file. The docstring-referenced `step_16_parent_throat_action_bundle_master_notes.md` has been deleted from the repository.

**Which is correct?**

Possible directions (the user picks one):

- **(a) The paper card is too terse and should be expanded.** Stage 018 genuinely owns all five script-side claim families. The paper card and a fresh `notes/stages/moving_throat_pde_stage018_*.md` need to be written so the script's assertions trace to declared deliverables. If this is the direction: a follow-up directive will name the .tex and notes edits and Codex will be told to update the script docstring's stale `step_16_*` reference; no script-side algebra changes.
- **(b) The script's extra families belong to a different stage.** Families 1–4 actually live in (e.g.) stage 015 ("Parent throat action master"), 016 ("Parent throat action candidate"), 017 ("Parent throat action weak-axisymmetric law"), or 021 ("Reduced Maxwell/mixed one-port normal form"). If this is the direction: a follow-up directive will tell Codex to delete A1–A17 (sympy) and M1–M7-mut (mathematica) from the stage-018 scripts and (separately, if needed) verify the deleted claims appear in the named upstream stage's scripts. The stage-018 audit then reduces to A18/A19 / M8 only, and even those should be lifted from a single Gaussian profile to a symbolic check of the bridge identity itself.
- **(c) The paper card is correct as-is and the script is over-scoped.** Stage 018's responsibility ends at exporting the two bridge integrals; the extra algebra is scaffolding from an earlier draft. If this is the direction: a follow-up directive will tell Codex to delete A1–A17 / M1–M7-mut, restructure A18/A19 / M8 to verify the abstract `M_Σ = ∫ μ_η β² dw` and `K_Σ = ∫ [T_w (β')² + (K_η + 6 T_Ω) β²] dw` identities symbolically (not just under the Gaussian collapse `μ_η = T_w = K_η + 6 T_Ω = 1`), and update the docstring.

In addition to choosing (a), (b), or (c), the user should also indicate whether the dangling `step_16_parent_throat_action_bundle_master_notes.md` reference in the SymPy docstring (line 2) reflects a missing notes file that should be restored, or simply a stale reference that should be replaced with the correct stage-018 notes path once one exists.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
