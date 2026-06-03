---
unit_id: 233
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-02T22:06:47-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 233

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py:115-126`

**Issue:** The "numerical recovery of the carried budgets" block is a circular round-trip. It defines `Pbar_num`, `robust_budget`, `nonempty_budget` as literals (lines 115–117), then computes `Pcrit_robust = Pbar_num*(1+robust_budget)` (line 119) and `recovered_robust = Pcrit_robust/Pbar_num - 1` (line 122). Algebraically `recovered_robust = Pbar_num*(1+robust_budget)/Pbar_num - 1 = robust_budget` for ANY inputs, so `assert_close(recovered_robust, robust_budget)` (line 125) cannot fail. Same for the nonempty branch (lines 120, 123, 126). The check verifies `b == b` and recovers nothing. The notes (§9, deliverable 8) require a genuine recovery of the two budgets from independently sourced data.

**Required change (in priority order — apply the first that is achievable in-repo):**

1. **Preferred — independent `P_crit`:** If an independent value of `P_crit` (the carried critical pressure) is available anywhere in this repo (search the Stage-240/Stage-241 source notes and any earlier stage script that defines a numeric `Pcrit`/`P_crit`), introduce it as a literal `Pcrit_num = <that value>` with an inline comment citing the source file. Then compute the budgets from the genuine identity `budget = Pcrit_num / Pbar_num - 1` and assert each equals the carried `robust_budget` / `nonempty_budget`. Note: a single `Pcrit_num` cannot reproduce two different budgets at one `Pbar_num`; the robust and nonempty budgets correspond to two different ceiling readings, so each branch needs its own independent critical-pressure literal (`Pcrit_robust_num`, `Pcrit_nonempty_num`) sourced from the data, NOT back-computed from the budget. Remove lines 119–120 (`Pcrit_* = Pbar_num*(1+budget)`).

2. **Fallback — non-trivial cross-relation:** If no independent `Pcrit` literal exists in-repo, do NOT keep the round-trip. Instead replace it with a check that ties the two carried budgets through a relation neither defines by construction. Acceptable: confirm that the two budgets, taken as `1+robust_budget` and `1+nonempty_budget`, satisfy the same `Pcrit/Pbar - 1` form with a SINGLE shared `Pbar_num` only if the two implied `Pcrit` values are independently carried; if they are not, this fallback is not achievable either.

3. **If neither is achievable:** append `## Blocked: F1` stating "No independent `P_crit` value found in-repo to break the round-trip; the carried budgets cannot be re-derived without an external critical-pressure source. Need user/notes to supply `P_crit`." Do NOT leave the existing tautological round-trip in place silently — if blocked, also delete the misleading print lines that claim "Recovered ... budget" (lines 134–135) or mark them clearly as "carried (not independently recovered)".

**Self-test reminder (already done by auditor, re-check before applying):** substituting your new definition of `recovered_robust` must NOT collapse to `robust_budget` for arbitrary inputs — i.e. if you perturb `Pcrit_robust_num` by a small delta the assertion must FAIL. If it still passes for any `Pcrit`, the round-trip is not broken.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 233` and confirm: (a) `Pcrit_robust`/`Pcrit_nonempty` are no longer `Pbar_num*(1+budget)`; (b) the recovery uses an independently sourced critical-pressure literal (or the block is `## Blocked`); (c) the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py`
- summary: Replaced the tautological budget round-trip with independent carried critical-pressure literals from the Stage 224 ceiling dictionary.
- deviation: none

## F2 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_mathematica_audit.wl` (NEW FILE — does not exist yet)

**Issue:** Stage 233 has only a SymPy script. Every claim is a closed-form algebraic identity plus one numeric relation, all independently verifiable in Mathematica. The dual-engine rule requires a second-engine derivation wherever Mathematica CAN verify the stage. It can here. Write a NEW independent-route Mathematica script at the target path above. The `_mathematica_audit.wl` suffix and the `mathematica/` directory are mandatory (matches the 224 existing `.wl` files).

**Independent route, not transliteration:** Use native Mathematica primitives (`FullSimplify`, `Simplify`, `Reduce`, `Refine`, `Together`/`Apart`) and a DIFFERENT decomposition than the `.py`. Do not echo the `.py`'s variable choreography line-for-line. Verify each in-file check with an explicit zero/equality test that exits nonzero on failure (e.g. ``If[FullSimplify[lhs - rhs] =!= 0, Print["FAIL Mn"]; Exit[1]]``). Strip `ConditionalExpression[0, ...]` wrappers from any `Solve`/`Reduce` output before the zero test, and for any pole/limit use `1/expr == 0` rather than `=!= Infinity`.

**Claim manifest** — the new `.wl` must independently verify each of these:

- **M1.** Compiler definitions reproduce the rigid-mouth collapse. With `Theta1 = dlnRtr`, `Xi1 = dlnNstar - Bstar*dlnRtr`, `ceta = epsEtaStar/(1-epsEtaStar)`, `R1 = -ceta*dlnEpsEta - Xi1`: confirm symbolically (not by bare substitution of 0) that the `dlnRtr`-coefficient of `Xi1` is `-Bstar` and of `Theta1` is `1`. Independent route: extract coefficients via `Coefficient[Xi1, dlnRtr]` and `Coefficient[Theta1, dlnRtr]` rather than substituting `dlnRtr -> 0`.
- **M2.** Track-locked specialization: imposing `dlnRtr = 0` gives `Theta1 = 0` and `Xi1 = dlnNstar`. Derive `Theta1 = 0` from M1's coefficient result (`Theta1` has no constant term and coefficient 1 on `dlnRtr`), not from a bare `Theta1 /. dlnRtr->0`.
- **M3.** Prefactor identity: with `P0 = N0/D0`, `P1 = (N01*D0 - N0*D01)/D0^2`, `XiLoad = N01/N0 - D01/D0`, confirm `FullSimplify[P1/P0 - XiLoad] == 0`.
- **M4.** Operator-rigidity specialization: `XiLoad /. D01 -> 0` equals `N01/N0` (confirm `FullSimplify[(XiLoad /. D01->0) - N01/N0] == 0`).
- **M5.** Transported ceiling equivalence: `gateRhs = Pcrit*mhat0^2/(DeltaNorm + Tquad) - 1` equals `Pcrit/PbarExpr - 1` with `PbarExpr = (DeltaNorm + Tquad)/mhat0^2` (confirm `FullSimplify[gateRhs - (Pcrit/PbarExpr - 1)] == 0`, with `mhat0 > 0`, `DeltaNorm >= 0`, `Tquad > 0`, `Pcrit > 0`).
- **M6.** Calibrated-branch ceiling: imposing `DeltaNorm = 0` on the $\bar P_0$ form yields `Pcrit*mhat0^2/Tquad - 1`. Derive this by starting from the `Pcrit/PbarExpr - 1` form, set `DeltaNorm -> 0`, and `FullSimplify` against `Pcrit*mhat0^2/Tquad - 1` — i.e. route through the $\bar P_0$ form, NOT by re-substituting into `gateRhs` and re-asserting the same expression (this addresses report F3).
- **M7.** $\bar P_0$-form: `Pbar = (DeltaNorm + Tquad)/mhat0^2` implies the gate RHS equals `Pcrit/Pbar - 1` (subsumed by M5; state explicitly).
- **M8.** Numerical budget relation: confirm `Pcrit/Pbar - 1` reproduces the carried budgets `0.367930328492646` and `0.737619063660757` at `Pbar = 0.002069792318062885`, USING THE SAME independent critical-pressure input(s) that F1 settles. If F1 ends `## Blocked` (no independent `Pcrit`), then M8 must instead only assert the algebraic identity `Pcrit/Pbar - 1` symbolically and state in a comment that the numeric budgets are carried-not-recovered (do not reproduce the F1 round-trip in Mathematica).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 233`, confirms the new `.wl` exists at the path above, all M1–M8 checks print no FAIL and the script exits 0, and confirms the `.wl` is not a transliteration of the `.py` (different decomposition for M1/M2/M6 at minimum).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_mathematica_audit.wl`
- summary: Added the independent Mathematica audit covering M1 through M8 with coefficient, Pbar-form, and independent critical-pressure checks.
- deviation: none

## F3 — insufficient_verification (folded into F2)

**Target:** addressed by F2's M2 and M6 independence requirements; no standalone `.py` edit required.

**Issue:** In the `.py`, line 95's calibrated-branch check re-substitutes `Delta_norm: 0` into `gate_rhs` and asserts equality with the same hand-written form (confirms `subs`, not an identity), and line 45's `assert_zero(Theta1_rm)` is a bare zero-substitution into an identity map. These are low-severity filler.

**Required change:** None in the `.py` standalone. Ensure the F2 `.wl` derives M2 and M6 by the independent routes specified (coefficient extraction for M2; route-through-$\bar P_0$ for M6) so the second engine does not merely echo the trivial substitutions.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_mathematica_audit.wl`
- summary: Covered the low-severity SymPy filler checks through the F2 Mathematica M2 coefficient extraction and M6 Pbar-form derivation.
- deviation: none

## F4 — stale upstream stage labels (low, optional, script-side comment/print only)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py:20` ("Exact Stage 188 branch-observable compiler")
- `:32` (`print("\nStage 188 observable compiler:")`)
- `:113` ("Numerical recovery of the carried Stage 224 budgets")
- `:128` (`print("\nStage 223 / Stage 224 carried numbers:")`)

**Issue:** These comments/print strings cite stale upstream stage numbers. The canonical citations per the card and notes are: the observable compiler is from **Stage 239** (not 188), the compatibility-point `Pbar` is from **Stage 240** (not 223), and the carried budgets/ceiling are from **Stage 241** (not 224). This is the known project-wide incomplete-renumber from the EM-extension realignment; the math is unaffected.

**Required change:** Update only these four comment/print strings to cite Stage 239 (compiler, lines 20, 32), Stage 240 (compatibility point, line 128), and Stage 241 (budgets/ceiling, lines 113, 128) to match the card/notes. Do NOT change any symbol, value, or assertion. Do NOT perform any other renumbering anywhere — this is a targeted four-string fix for this file only. If unsure which exact upstream number applies to a given string, prefer matching the notes §9 wording (Stage 239 / Stage 240 / Stage 241) and append `## Applied: F4` noting the mapping you used.

**Verification command:**
`redteam exec-sympy 233` exits 0 and the four strings now read 239/240/241; no assertion or numeric value changed.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py`
- summary: Updated the four targeted comment and print labels to Stage 239 for the compiler, Stage 240 for Pbar, and Stage 241 for budgets and ceiling.
- deviation: none
