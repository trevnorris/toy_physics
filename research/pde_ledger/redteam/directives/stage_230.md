---
unit_id: 230
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T17:36:17-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 230

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl`

**Issue:** Stage 230 has only a SymPy verifier. It is checkpoint:false and NOT status-only, and every step is a native Mathematica primitive (rational-function `Simplify`/`Together`, `D`, `Limit`, sign via `Reduce`/`Resolve`, `Log`, `N`). The dual-engine rule requires an independent Mathematica audit because Mathematica CAN verify this stage. Create a NEW, INDEPENDENT derivation — do NOT transliterate the `.py`.

**Anti-transliteration guard (mandatory):**
- Use Mathematica-native primitives throughout: `D[...]` for derivatives, `Limit[..., dir -> ...]` for the endpoint limits, `Reduce`/`Resolve[ForAll[...]]` or `Simplify[... < 0, assumptions]` for sign claims, `Solve` for the threshold inversions, `Log` and `N` for numerics.
- Do NOT mirror the `.py` variable choreography. In particular, do NOT reproduce the SymPy `dynamic_ceilings` helper line-for-line; instead express each ceiling directly as a piecewise/`Which` construct and let Mathematica evaluate it. Prove monotonicity by `Reduce[D[Splus,R] < 0 && R >= 0, R]` (a genuine sign proof over the half-line) rather than by reading off a constant's sign.
- Prove the sign-flip and onset thresholds by `Solve`, not by substituting a pre-defined value back in (avoid repeating F2's tautology on the Mathematica side).
- Establish the static-first "everywhere" claim by combining a `Reduce`-proved monotonicity with the `Limit[..., R -> Infinity]` infimum, then `N`-compare the infimum against the static budgets — do NOT just check the same four sample points the `.py` checks.

**Claim manifest** (each must be independently verified by the new `.wl`):

- **M1** — Exact classifier and its monotonicity. `R_ND(xi,delta) = 72 delta^2 (1-xi) / ((9 delta + 11 xi)(9 delta^2 + 18 delta xi + 11 xi^2))`. Verify (a) `R_ND(0,delta) = 8/(9 delta)`; (b) `Limit[R_ND, xi -> 1, Direction -> "FromBelow"] = 0`; (c) `D[R_ND, xi] < 0` for `0 <= xi < 1, delta > 0` (prove via `Reduce`, do not just match a hand-written derivative).
- **M2** — Affine share compiler. With share weights `w_num = R/(1+R)`, `w_den = 1/(1+R)` and the carried per-unit rigid slopes (carry-forward literals; see note below) `s_plus_num = -0.508643465308977`, `s_plus_den = -0.301516097158113`, `s_minus_num = -0.334368725711457`, `s_minus_den = +0.411024574532864`, define `S_pm(R) = (R*s_pm_num + s_pm_den)/(1+R)`. Verify `D[S_plus,R] < 0` and `D[S_minus,R] < 0` for `R >= 0` (via `Reduce`).
- **M3** — Sign-flip threshold. Solve `S_minus(R) == 0` for `R` and verify the unique nonnegative root equals `R_star = s_minus_den/(-s_minus_num)`, with `N[R_star, 30] ≈ 1.229255438463336`. Verify `S_plus(R_star) < 0`, `S_minus(R) > 0` for `0 <= R < R_star`, `S_minus(R) < 0` for `R > R_star` (via `Reduce`).
- **M4** — Onset threshold. Solve `R_ND(0,delta) == R_star`, i.e. `8/(9 delta) == R_star`, for `delta`; verify the solution equals `delta_dyn_star = 8/(9 R_star)` with `N[delta_dyn_star, 30] ≈ 0.723111617875019`.
- **M5** — Dynamic ceilings in `|eps Xi_1|`. With carried `RQ_minus = 30.199907560250075`, `RQ_plus = 36.171186483269487`, `RQ_req = 21.854566296358396`, set `ell_minus = Log[RQ_minus/RQ_req] ≈ 0.323428979934714`, `ell_plus = Log[RQ_plus/RQ_req] ≈ 0.503852964869151`. Verify the endpoints: `B_both(R=0) ≈ 1.671064893775584` and `B_nonempty(R=0) = Infinity`; `Limit B_both as R->Infinity ≈ 0.967282389363822`; `Limit B_nonempty as R->Infinity ≈ 0.990581810705233`. Also reproduce the four sample rows at `R = 0, 1, R_star, 10` (`S_plus, S_minus, B_both, B_nonempty`) matching the SymPy sample table.
- **M6** — Static-first theorem (the card's Output). With carried static budgets `B_stat_both = 0.367930328492646`, `B_stat_nonempty = 0.737619063660757`, verify `inf_{R>=0} B_both = (R->Infinity limit) > B_stat_both` and `inf B_nonempty (over finite region) = (R->Infinity limit) > B_stat_nonempty`, where "infimum at R->Infinity" is justified by the M2 monotonicity proof.

Carry-forward note (do NOT re-derive, just carry as literals — they are anchored upstream): the four per-unit slopes, the three R_Q figures, the two Xi_1 unit norms, and the two static budgets all originate in `scripts/moving_throat_pde_stage228_..._sympy_audit.py` (slopes lines 418-421, R_Q base line 368, R_Q_req line 432, Xi_1 norms lines 294/298, static budgets lines 441-442) and the classifier in `scripts/moving_throat_pde_stage229_..._sympy_audit.py`. Use the SAME numeric literals; the second engine independently re-derives the *downstream algebra* (compiler, thresholds, ceilings, inequalities), not these frozen inputs.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 230` and confirms the new `.wl` exists at the path above, the M1-M6 checks appear, and the script exits 0; `R_star`, `delta_dyn_star`, the endpoint ceilings, and the two static-first inequalities must agree with the SymPy output.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M6 with Reduce/Solve/Limit-based checks and dynamic ceiling comparisons.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.py:119`

**Issue:** Line 119 asserts `sp.simplify(onset.subs(delta, delta_dyn_star) - R_star) == 0`. Because `delta_dyn_star` is defined on line 114 as exactly `8/(9*R_star)` and `onset = 8/(9*delta)` (line 63), substituting `delta -> 8/(9*R_star)` returns `R_star` by pure construction; the residual is identically zero regardless of any physics. The check verifies only that `8/(9·x)` is its own inverse, not the onset-threshold derivation. The real numeric anchoring of `delta_dyn_star` is already at line 115.

**Required change:**
Replace the definitional round-trip with an independent solve of the onset-threshold equation. Specifically, change line 119 from:

```python
    assert sp.simplify(onset.subs(delta, delta_dyn_star) - R_star) == 0
```

to a block that solves `onset(delta) = R_star` for `delta` and confirms the solution equals `delta_dyn_star` (this exercises the inversion the notes §4.2 actually claims, rather than assuming it):

```python
    delta_solutions = sp.solve(sp.Eq(onset, R_star), delta)
    assert len(delta_solutions) == 1
    assert sp.simplify(delta_solutions[0] - delta_dyn_star) == 0
```

Here `onset` is the symbolic `8/(9*delta)` from line 63 and `delta` is the symbol from line 47. Do not change lines 63, 111, or 114. Leave line 118 (`S_-(R_*) = 0`) untouched.

**Self-test (already performed by auditor):** `solve(8/(9*delta) - R_star, delta)` returns the single root `delta = 8/(9*R_star)`, which is exactly `delta_dyn_star`; substituting the anchored `R_star ≈ 1.2292554` gives `delta ≈ 0.7231116`, a nonzero literal, so the replacement assertion is non-trivially satisfiable and matches line 115's value. No new constant is introduced.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 230` and confirms the new solve-based check appears at/near line 119 and the script exits 0 with all checks passing.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.py`
- summary: Replaced the onset-threshold substitution round trip with a solve-based inversion check for `delta`.
- deviation: none

## Authorized notes renumber (USER-AUTHORIZED 2026-06-02)

The user authorized renumbering stale VII.1 notes-prose stage labels to canonical. The audit logged a notes-side renumbering drift here (notes reference pre-renumber stage numbers; the card/appendix/scripts use canonical). Notes-only cleanup in THIS fix loop (Codex applies notes/*.md; Claude reviews):
- In `notes/stages/moving_throat_pde_stage230_..._sympy_audit.md`, renumber every stale stage-number reference to match the canonical numbering used in this stage's SymPy script comments + the paper card (self-reference → Stage 230; cited upstream stages → the numbers the .py comments use). Math/content unchanged.
- Do NOT touch scripts, paper.tex, or appendix. Acceptance: notes stage labels match the .py comments + card. Append `## Applied: notes-renumber`.

## Applied: notes-renumber

- files_changed:
  - `notes/stages/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.md`
- summary: Renumbered stale notes references from Stage 245/246/247 to canonical Stage 228/229/230 labels and updated the supporting filename.
- deviation: none
