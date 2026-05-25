---
unit_id: 011
batch: I.1
created_at: 2026-05-25T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 011

This directive contains exactly one finding (`paper_misalignment / paper_missing_script_claim`). It requires user resolution before any code or paper edits are applied. Codex must not edit `paper/`, `notes/`, or the stage 011 scripts in response to this directive.

## F1 — paper_misalignment

**Subtype:** `paper_missing_script_claim`

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_011.tex:37-39` quote:
  > Stage~011 exports the `z_0`-cancellation and the projected `P_2` compatibility bridge \eqref{eq:stage011-dcompat}.
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:44` quote:
  > `P_2`-lane bridge from projected electromagnetic slots into the grouped normalization language.

**Script side (delivers more than the paper card enumerates):**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py:132-155` — transported-target K-surface and compatibility variant (paper card mentions only the fixed-target form).
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py:168-176` — constant-prefactor branch factorizations imposing `P_2 = 0`, `P_4 = 0`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py:179-196` — real-Y20 weak-axisymmetric lane lambdas `(1, 1/2, -1)` and grouped trace + `b = 3a`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py:199-204` — static Xi1 prefactor slope `Xi1^{(proj,static)} = na/N_0 + za/D_0`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py:207-211` — per-lane `u_2` slope `(D_0 z_{2a} - D_2 z_a)/D_0^2`.
- Mathematica counterparts at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl:97-153` mirror items 1, 3, 4, 5 above.

The stage 011 script's core (sympy A1-A9, A13, A15-A17 and Mathematica M1-M5) faithfully verifies what the paper card declares as Output. The items above are additional verified identities the paper card does not anchor.

## Resolve before fix_loop

The stage 011 script verifies more than the stage card declares as `Output`. The card promises only (i) the K-eliminated compatibility surface `C = N_0/P_{0,target} - 3 S^2/T`, (ii) its first variation `δC = n_0/P_{0,target} - 6 S z_2/T + 3 S^2 z_4/T^2`, and (iii) the `z_0`-cancellation property. The script additionally verifies the transported-target variant, a constant-prefactor branch factorization, the real-Y20 weak-axisymmetric lane lambdas, the static Xi1 prefactor slope, and a per-lane `u_2` slope. Which direction is correct?

Possible directions (the user picks one):
- **(a) Paper is correct as scoped — trim the script.** Remove sympy L132-155 (transported-target variant), L168-176 (constant-prefactor factorization), L179-211 (Y20 / Xi1 / u2 lane). Remove Mathematica `.wl` L97-153. Verify no downstream stage script imports or quotes these results from stage 011 before deleting. If a downstream stage does need them, move the relevant block into that downstream stage's audit script (which presumably has a paper-side anchor that declares the result).
- **(b) Script is correct — extend the paper card.** Add to `stage_011.tex` an Output block (or supplemental equations) that enumerates: the transported-target compatibility variant (`compat_transport = D_{0,target} - 3 S^2/T` with its `z_0`-independent first variation), the constant-prefactor branch identities `P_2 = (N_2 - 2 D_2 N_0/D_0)/D_0` and `P_4 = (N_4 - N_4^{*})/D_0`, the real-Y20 grouped lane signature `(λ_{20}, λ_{21}, λ_{22}) = (1, 1/2, -1)` with `xbar = x^{(0)}` and `b = 3a`, the static Xi1 slope `Xi1^{(proj,static)} = na/N_0 + za/D_0`, and the per-lane `u_2` slope `(D_0 z_{2a} - D_2 z_a)/D_0^2`. Optionally also add per-stage notes under `notes/em_projected/step_011_*.md` so future audits have a derivation trail.
- **(c) Split across stages.** Keep only the core K-eliminated compatibility surface in stage 011 (matching the current card), and move the Y20 / Xi1 / u2 / constant-prefactor work into the stage(s) that actually consume them (likely 012 "primitive bridge", 022 "grouped normalization bridge", or 023 "full grouped bundle"). This requires both script-side relocation and paper-side anchoring at the destination stage.

Once the user has chosen, a follow-up directive should specify exactly which file:line ranges to edit and which paper card(s) (if any) to amend. The orchestrator will not invoke Codex on stage 011 until that follow-up arrives.
