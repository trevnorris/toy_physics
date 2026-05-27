---
unit_id: 074
batch: III.4
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
---

# Codex directive -- unit 074 (resolved orchestrator-direct)

The single F1 `paper_misalignment` (alpha value mismatch) was resolved by the user as direction (a) on 2026-05-27 and applied by the orchestrator directly:

- `paper/stages/stage_074.tex:31`: `\frac{128}{\sqrt5}` -> `\frac{111}{\sqrt5}`
- `notes/stages/moving_throat_pde_stage074_family1_healing_lock.md:117`: `179/sqrt(5)` -> `111/sqrt(5)`
- `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:63`: `179/sqrt(5)` -> `111/sqrt(5)` (same typo, collateral fix)
- `scripts/.../stage074..._sympy_audit.py`: added `expect_zero("alpha_ref - 111/sqrt(5)", alpha_ref - 111/sqrt(5))` after the `kappa_ref` assertion.
- `mathematica/.../stage074..._mathematica_audit.wl`: added `expectZero["alpha_ref - 111/sqrt(5)", alphaRef - 111/Sqrt[5]]` after the `kappaRef` assertion.

No Codex invocation required. Verifier will confirm the new assertion appears and both scripts exit 0.

Do NOT touch paper.tex, notes/, or any prose document.
Do NOT edit the scripts to "fix" a paper_misalignment.
Do NOT run python or mathematica.

## F1 -- paper_misalignment

**Subtype:** value_mismatch

**Paper side:**

- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_074.tex:26-31` quote:
  > `\kappa=\frac95\Lambda_\ell^2=\frac{12321}{5}, \qquad \alpha=\sqrt\kappa=\frac{128}{\sqrt5}`

- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage074_family1_healing_lock.md:113-117` quote:
  > `alpha = sqrt(12321/5) = 179/sqrt(5) ~ 49.6407091.`

**Script side:**

- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:59-66` script-output quote:
  > `alpha = 111*sqrt(5)/5`,  `alpha (numeric) = 49.640709100495331260`

- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.wl:48-55` script-output quote:
  > `alpha = 111/Sqrt[5]`,  `alpha (numeric) = 49.64070910049533126028365544583433242664`

## Resolve before fix_loop

The paper card states (boxed) `alpha = sqrt(kappa) = 128/sqrt(5)`, but the same equation also boxes `kappa = 12321/5`. Since `sqrt(12321/5) = 111/sqrt(5)` (because `111^2 = 12321`), the paper's three claims `kappa = 12321/5`, `alpha = sqrt(kappa)`, and `alpha = 128/sqrt(5)` are mutually inconsistent. The notes likewise state `alpha = 179/sqrt(5)`, which is also wrong as a literal (`179^2 = 32041 != 12321`), but the notes' decimal `~49.6407091` is correct and matches the scripts' `111/sqrt(5)`. The scripts compute the arithmetically correct value `111/sqrt(5)` but do not assert it.

Question for the user: which value of `alpha` is intended?

Possible directions (the user picks one):

- (a) `alpha = 111/sqrt(5)` is correct (consistent with `kappa = 12321/5` and with both engines' output). In `paper/stages/stage_074.tex` line 31, change `\frac{128}{\sqrt5}` to `\frac{111}{\sqrt5}`. In `notes/stages/moving_throat_pde_stage074_family1_healing_lock.md` section 5, change `179/sqrt(5)` to `111/sqrt(5)`. No script change required.
- (b) `kappa` is wrong and the correct `kappa` is `(128)^2/5 = 16384/5`, with `alpha = 128/sqrt(5)`. This would require revising the entire chain `kappa = (9/5) Lambda_ell^2`, the Stage 54 coefficient `4 chi_s^2 + (4/5) Lambda_ell^2`, and the carry-forward `Lambda_ell = 37`. Highly unlikely (would invalidate Stages 56, 073, 075+), but listed for completeness.
- (c) Both literals are wrong relative to a third intended derivation -- flag for deeper review.

Recommended direction: (a). Both engines independently produce `111/sqrt(5)`; the notes' decimal also corresponds to `111/sqrt(5)`; only the two prose literals (paper `128`, notes `179`) disagree, and they disagree with each other as well as with the engines, which is a strong signature of a copy-paste typo on the paper-prose side.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## Self-test trace (auditor-side, for the user's reference)

- Verified `111^2 = (100+11)^2 = 10000 + 2200 + 121 = 12321`. Confirmed.
- Verified `128^2 = 16384 != 12321`.
- Verified `179^2 = (180-1)^2 = 32400 - 360 + 1 = 32041 != 12321`.
- Verified numeric: `111/sqrt(5) = 111/2.2360679... ~ 49.64070910...`. Matches the scripts and the notes' decimal.
- Confirmed no script-side assertion currently fails: the scripts only print `alpha`, they do not assert any literal, so the discrepancy is invisible to the current verifier and was not caught by the v1 audit.
