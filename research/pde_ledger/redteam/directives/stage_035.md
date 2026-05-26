---
unit_id: 035
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 035

This directive contains a single `paper_misalignment` finding. **Do not apply any edit.** The orchestrator is holding for user resolution. The user must choose a direction before any follow-up directive authorizes a paper- or notes-side edit. Codex must not edit `paper/`, `notes/`, or any prose document on its own initiative, and the scripts should be left as-is until the user has decided.

## F1 — paper_misalignment

**Subtype:** `target_mismatch` (also acts as `notes_contradicts_script`, since the notes carry the same wrong polynomial as the paper card)

**Paper side:**

- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_035.tex:65-73` quote (boxed `eq:app-stage035-F-derivative`):

  ```
  \frac{\partial F}{\partial\xi}
  =
  \frac{(9\delta+11\xi)^3
  \left(81\delta^3+72\delta^2+206\delta^2\xi+297\delta\xi^2+138\xi^3\right)}
  {81(1-\xi)^2(9\delta^2+18\delta\xi+11\xi^2)^3}>0
  ```

- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md:85-87` quote:

  ```
  dF/dxi
  = (9 delta + 11 xi)^3
    [ 81 delta^3 + 72 delta^2 + 206 delta^2 xi + 297 delta xi^2 + 138 xi^3 ]
    / [ 81 (1 - xi)^2 (9 delta^2 + 18 delta xi + 11 xi^2)^3 ]
  > 0
  ```

**Script side:**

- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py:66-70` quote:

  ```
  dF_target = sp.simplify(
      (9 * delta + 11 * xi) ** 3
      * (81 * delta**3 + 72 * delta**2 + 189 * delta**2 * xi + 297 * delta * xi**2 + 121 * xi**3)
      / (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)
  )
  ```

- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:65-68` quote:

  ```
  dFTarget = FullSimplify[
    (9*delta + 11*xi)^3*(81*delta^3 + 72*delta^2 + 189*delta^2*xi + 297*delta*xi^2 + 121*xi^3)/
      (81*(1 - xi)^2*(9*delta^2 + 18*delta*xi + 11*xi^2)^3),
    Assumptions -> $Assumptions
  ];
  ```

- The Mathematica saved transcript (`mathematica/output/.../moving_throat_pde_stage035_..._mathematica_audit.txt:12`) prints:

  ```
  dF/dxi = ((9*delta + 11*xi)^3*(81*delta^3 + 297*delta*xi^2 + 121*xi^3 + 9*delta^2*(8 + 21*xi)))/...
  ```

  i.e., `9*delta^2*(8 + 21*xi) = 72 delta^2 + 189 delta^2 xi`. The `expectZero` assertion `dF - dFTarget` passes (residual `= 0`, line 13), confirming the engine-derived `D[f, xi]` reduces to the script's literal polynomial.

## Resolve before fix_loop

The paper card and notes box the bracket polynomial of `partial F / partial xi` with coefficients `206 delta^2 xi` and `138 xi^3`, while both SymPy and Mathematica scripts (and independent hand derivation, verified numerically at `delta = 0` where `F = 121/[81(1-xi)]` so `dF/dxi = 121/[81(1-xi)^2]`) give `189 delta^2 xi` and `121 xi^3`. Which is correct?

Possible directions (the user picks one):

- (a) **Scripts are correct → update paper and notes.** Direction: in `paper/stages/stage_035.tex` line 71 replace `\left(81\delta^3+72\delta^2+206\delta^2\xi+297\delta\xi^2+138\xi^3\right)` with `\left(81\delta^3+72\delta^2+189\delta^2\xi+297\delta\xi^2+121\xi^3\right)`. In `notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md` line 86 replace `[ 81 delta^3 + 72 delta^2 + 206 delta^2 xi + 297 delta xi^2 + 138 xi^3 ]` with `[ 81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3 ]`. No script change. No re-run of either engine is required; the script assertions already pass with the correct polynomial. This is the direction the auditor recommends based on independent hand derivation and engine agreement.

- (b) **Paper is correct → update scripts.** Direction: in both `scripts/moving_throat_pde_stage035_..._sympy_audit.py` line 68 and `mathematica/moving_throat_pde_stage035_..._mathematica_audit.wl` line 66, replace `189 delta^2 xi` with `206 delta^2 xi` and `121 xi^3` with `138 xi^3`. Re-run sympy and Mathematica via `redteam exec-sympy 035` and `redteam exec-mathematica 035`. The auditor predicts both will FAIL with a nonzero residual under this direction, because the engine-derived `D[f, xi]` does not equal the paper's polynomial (see the `delta = 0` numerical check in the report). If they FAIL, return to option (a).

- (c) **Both are derived from a third source that contradicts both → flag for deeper review.** If the user believes there is an upstream Stage-034 derivation that produces a different bracket polynomial than either form, escalate: re-derive `partial F / partial xi` from scratch using the Stage-034 definition of `N_-(x)` and check which polynomial form (if either) reproduces the result.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. After the user's choice, a follow-up directive may be issued that explicitly authorizes either a paper-side edit (direction a) or a script-side edit (direction b).
