---
unit_id: 020
batch: I.2
created_at: 2026-05-25T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 020

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_020.tex:53-55` quote: "Stage~020 exports the weak-axisymmetric parent packet \eqref{eq:stage020-wall-slopes}--\eqref{eq:stage020-residual-xi}." (entire Output paragraph; no Y20 / Gaunt / lane-ratio anchor anywhere in the card)
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:62` quote: "020 & Parent throat action weak-axisymmetric packet & ... & Weak-axisymmetric packet exported from the parent throat-action bundle." (no Y20 mention)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:45-50` quote:
  ```
  lam20 = real_y20_square_ratio(0)
  lam21 = real_y20_square_ratio(1)
  lam22 = real_y20_square_ratio(2)
  assert_zero('Y20 overlap lane 20', lam20 - 1)
  assert_zero('Y20 overlap lane 21', lam21 - sp.Rational(1, 2))
  assert_zero('Y20 overlap lane 22', lam22 + 1)
  ```
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl:34-67` quote: Y20 lane ratio assertions (`lambda0 - 1`, `lambda1 - 1/2`, `lambda2 - (-1)`, `crossOne - 0`, `crossTwo - 0`).

## Resolve before fix_loop

The paper card and appendix row for stage 020 export only the algebraic packet `delta M_Sigma`, `delta K_Sigma`, `D_{n1}`, `K_1`, `H_even`, `Xi_1`, the wall-slope solve, and the residual `Xi_1`. They do not mention spherical-harmonic / Gaunt lane ratios. Both engine scripts, however, devote a substantial block to asserting `lambda_{20} = 1`, `lambda_{21} = 1/2`, `lambda_{22} = -1`, and `gaunt(2,2,2,0,m,m) = 0` for `m = 1, 2`. Is the Y20 lane block part of stage 020's deliverable (and the paper card needs to be amended), or is it scaffolding that belongs to a different stage (e.g., 023 "Full grouped bundle and projectors", which the appendix row at stage_appendix_part01.tex:68 says includes "isotropic branch test" — that may be where the Y20 ratios live)?

Possible directions (the user picks one):
- (a) Y20 ratios belong here → add a paragraph and equation to `paper/stages/stage_020.tex` (e.g., after the "Even compensation" paragraph) that anchors the three lane ratios and the selection-rule vanishing; no script change. Then F2 below should be applied to fix the lambda_{20} tautology.
- (b) Y20 ratios belong elsewhere → remove sympy lines 13 (the `gaunt` import), 20-26 (`real_y20_square_ratio` definition), and 45-50 (the three `assert_zero` calls). Mirror: remove `wl` lines 6-7 (the `ClearAll[GauntIntegral]` and `SetAttributes`), 21-30 (`GauntIntegral` definition), and 34-68 (the angular block including the `M1 OK` print). After removal, F2 is moot.
- (c) Y20 ratios are partially load-bearing here but in a sense the card has not articulated → user writes a short note in `notes/stages/moving_throat_pde_stage020_*.md` that nails down why the ratios appear, then chooses (a) or (b) given that anchor.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:20-26, 45, 48`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl:34-35, 49-51`

**Issue:**

In `real_y20_square_ratio(m)` (sympy lines 20-26), `base` is `gaunt(2,2,2,0,0,0)` and the return for `m = 0` is `gaunt(2,2,2,0,0,0) / base`, which is identically `1` by construction. The downstream assertion `assert_zero('Y20 overlap lane 20', lam20 - 1)` (line 48) therefore cannot fail no matter what the Gaunt code returns. Same defect in Mathematica at lines 34-35 (`overlapBase`, `lambda0`) and lines 49-51 (assertion on `lambda0 - 1`). The non-trivial assertions in the block (`lam21`, `lam22`, same-sign vanishing) DO have content; only `lam20`/`lambda0` is the tautology.

**Resolution gate:** Apply only if user has chosen direction (a) above (Y20 ratios stay in the script). If user chose direction (b), the entire Y20 block is removed and this finding is moot.

**Required change (if F1 direction (a) is chosen):**

In `sympy_audit.py`, anchor `gaunt(2,2,2,0,0,0)` to its known closed form before forming the ratios. After line 32 (the import of `gaunt`) is fine; concretely, add after line 44 (`Xi1 = sp.expand(...)`) and before line 45:

Before (line 45-50):
```
    lam20 = real_y20_square_ratio(0)
    lam21 = real_y20_square_ratio(1)
    lam22 = real_y20_square_ratio(2)
    assert_zero('Y20 overlap lane 20', lam20 - 1)
    assert_zero('Y20 overlap lane 21', lam21 - sp.Rational(1, 2))
    assert_zero('Y20 overlap lane 22', lam22 + 1)
```

After:
```
    # Anchor the Y20-self overlap to its closed form, so the lane-ratio
    # denominators are not validated against themselves.
    gaunt_base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    assert_zero('Y20 base gaunt closed form', gaunt_base - sp.sqrt(5/sp.pi)/(7*sp.sqrt(sp.pi)/sp.sqrt(1)) * sp.Integer(1))
```

The exact closed form of `gaunt(2,2,2,0,0,0)` is `sqrt(5/(4*pi)) * 2/7`. Codex: compute this with sympy and replace the closed-form expression in the assertion above with the value sympy actually returns from `sp.sqrt(sp.Rational(5,4)/sp.pi) * sp.Rational(2,7)` — do not invent an answer. If you cannot determine the exact closed form mechanically without running sympy, append `## Blocked: F2` with the question instead.

In `mathematica_audit.wl`, mirror the same change: introduce
```
gauntBase = FullSimplify[GauntIntegral[2, 2, 2, 0, 0, 0]];
If[!TrueQ[FullSimplify[gauntBase - Sqrt[5/(4 Pi)]*(2/7)] === 0],
  Print["FAIL: M1 gauntBase closed form"]; Exit[1]];
```
before line 35, and leave the existing `lambda0` check in place (it is still tautological, but the new `gauntBase` check now anchors the base value). Codex: insert the new lines but do not delete the existing `lambda0` check; the orchestrator wants both for trace continuity.

**Claim manifest:** N/A (this is not a missing-script finding).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 020` and `redteam exec-mathematica 020` and confirm:
1. New `Y20 base gaunt closed form` assertion appears at sympy line ~46.
2. New `gauntBase` assertion appears at wl around line 35.
3. Both scripts exit 0 with `STATUS: PASS`.
4. The new check could in principle fail if the closed-form value were perturbed.
