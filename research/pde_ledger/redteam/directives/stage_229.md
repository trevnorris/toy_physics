---
unit_id: 229
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T17:34:21-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 229

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

The `paper_misalignment` finding(s) (F1) have been RESOLVED by the user (2026-06-02) — see the `## RESOLVED` block at the END of this directive and apply the authorized notes-only edit there as part of this fix loop. (Codex applies notes/*.md; Claude reviews.)

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`math -script <path>` for the new Mathematica file) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex or the appendix. The ONLY authorized prose/notes edit is the notes-only change specified in the `## RESOLVED` block at the END of this directive (user-authorized 2026-06-02); apply exactly that and make no other notes/prose edits.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md:292` quote: "189\xi^3+297\delta\xi^2+333\delta^2\xi+81\delta^3-72\delta^2"
- same file `:302` quote: "363\xi^2+594\delta\xi+333\delta^2 > 0"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.py:96` quote: "expected_P = 121 * xi**3 + 297 * delta * xi**2 + 333 * delta**2 * xi + 81 * delta**3 - 72 * delta**2"

## Resolve before fix_loop

The notes (line 292) write the crossover cubic leading term as `189\xi^3`, but the script asserts `121\xi^3` (line 96). The math points decisively to `121`: expanding the classifier-equals-one condition gives the $\xi^3$ term as $11\xi\cdot11\xi^2=121\xi^3$, the script's in-script re-derivation (`numerator` at line 95-97) forces `121`, and the notes' OWN derivative at line 302 (`363\xi^2`) is $\partial_\xi(121\xi^3)$, not $\partial_\xi(189\xi^3)=567\xi^2$ — so the notes contradict themselves. Which is correct?

Possible directions (the user picks one):
- (a) Script is correct (strongly indicated) → change notes:292 leading term `189\xi^3` to `121\xi^3`. No script change. This is a notes-only edit; per the RESOLVED block below, Codex applies it and Claude reviews.
- (b) Notes are correct → then both the script's `expected_P` (line 96) AND the notes' own derivative (line 302) are wrong; flag for deeper review, since this would imply the classifier $\mathcal R_{ND}$ form is itself mis-stated.
- (c) Defer; treat as a pure notes typo to be fixed in the notes-sync pass.

F1 is RESOLVED by the user (direction (a); see the ## RESOLVED block at the end): apply the authorized notes edit (189→121), then proceed with F2.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md`
- summary: Corrected the notes crossover-cubic leading coefficient from `189\xi^3` to `121\xi^3`.
- deviation: none

## F2 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_mathematica_audit.wl`

**Issue:** No Mathematica script exists for unit 229. The unit is `is_status_only_candidate: False` and `is_checkpoint: False`, so the dual-engine rule requires a `.wl` because Mathematica CAN independently verify every claim here (rational-function reduction, symbolic log-derivatives, one-sided limits, a cubic from clearing denominators, derivative positivity, and real-root location). Create an INDEPENDENT re-derivation — NOT a transliteration of the `.py`.

**Required change:**
Create the `.wl` at the Target path. Anti-transliteration guard — you MUST:
- Use native Mathematica primitives only: `Together`, `Cancel`, `Numerator`, `Denominator`, `Factor`, `D[Log[...], xi]`, `Limit[..., xi -> 1, Direction -> "FromBelow"]`, `CoefficientList`, `Resolve`/`Reduce`, `NSolve`/`Root`.
- Derive the cubic by `Numerator[Together[RND - 1]]` (then `Expand`/`Factor`), NOT by writing a literal `expected_P` and comparing — let Mathematica produce the polynomial, then read off / assert its coefficients via `CoefficientList[poly, xi]`.
- Build `s_-`, `N_-`, and `F` from the constants `kappa0sq = 8/Pi^2`, `kappa1sq = 16/(9 Pi^2)` directly; do not import the SymPy intermediate forms.
- Use a different check idiom than the SymPy `assert simplify(...)==0`: an `expectZero[expr_] := If[Simplify[expr] =!= 0, Print["FAIL ..."]; Exit[1]]` helper, and for poles use the reciprocal trick (`1/expr == 0`) rather than `=!= Infinity`. Strip any `ConditionalExpression[0, ...]` from `Simplify`/`Reduce` results before the zero test.
- Declare symbols real/positive via `Assuming[xi > 0 && xi < 1 && delta > 0, ...]` consistent with the notes' stable-branch domain $0\le\xi<1,\ \delta>0$.

**Claim manifest** (the new `.wl` must independently verify each):
- M1: $N_-(A\xi, A\delta) = \dfrac{8\beta_0}{\pi^2 A}\,F(\xi,\delta)$ with $F=\dfrac{(9\delta+11\xi)^4}{81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2}$, built from $\kappa_0^2=8/\pi^2,\ \kappa_1^2=16/(9\pi^2)$ and $s_-=\dfrac{[\kappa_0^2(x+\Delta K_{ax})+\kappa_1^2 x]^2}{\kappa_0^2(x+\Delta K_{ax})^2+\kappa_1^2 x^2}$.
- M2: $F = F_{\rm num}\,F_{\rm den}$ with $F_{\rm num}=\dfrac{(9\delta+11\xi)^4}{81(9\delta^2+18\delta\xi+11\xi^2)^2}$, $F_{\rm den}=\dfrac{1}{1-\xi}$.
- M3: $L_{\rm num}:=\partial_\xi\ln F_{\rm num}=\dfrac{72\delta^2}{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)}$ and $L_{\rm den}:=\partial_\xi\ln F_{\rm den}=\dfrac{1}{1-\xi}$ (both derived by `D[Log[...],xi]`, then `Simplify`).
- M4: $\mathcal R_{ND}:=L_{\rm num}/L_{\rm den}=\dfrac{72\delta^2(1-\xi)}{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)}$.
- M5: onset $\mathcal R_{ND}(0,\delta)=\dfrac{8}{9\delta}$.
- M6: near-softening $\displaystyle\lim_{\xi\to1^-}\mathcal R_{ND}=0$ (use `Direction->"FromBelow"`); and $\displaystyle\lim_{\xi\to1^-}L_{\rm num}=\dfrac{72\delta^2}{(9\delta+11)(9\delta^2+18\delta+11)}$ (finite).
- M7: crossover cubic — let $\mathcal P:=\mathrm{Numerator}[\mathrm{Together}[\mathcal R_{ND}-1]]$ (taken with the sign so the $\xi^3$ coefficient is positive); assert via `CoefficientList[P, xi]` that $\mathcal P = 121\xi^3+297\delta\xi^2+333\delta^2\xi+81\delta^3-72\delta^2$. (NOTE: this independently re-derives the coefficient that F1 disputes; the correct value is `121`. Do NOT hardcode the literal cubic as the source — derive it, then verify the coefficients match these targets.)
- M8: $\partial_\xi\mathcal P = 363\xi^2+594\delta\xi+333\delta^2$, and this is $>0$ for all $\xi\ge0,\ \delta>0$ (verify via `Resolve[ForAll[{xi,delta}, xi>=0 && delta>0, 363 xi^2 + 594 delta xi + 333 delta^2 > 0]]` returning `True`, or by confirming the discriminant in $\xi$ is negative and the leading coefficient positive).
- M9: $\mathcal P(0,\delta)=9\delta^2(9\delta-8)$, giving the threshold $\delta=8/9$ (root of $9\delta-8$).
- M10: sample crossover roots in $(0,1)$ of $\mathcal P(\xi,\delta_0)=0$ via `Root`/`NSolve`: $\delta_0=1/4\Rightarrow\xi_*\approx0.107223051105697$; $\delta_0=1/2\Rightarrow\xi_*\approx0.081847937860074$; $\delta_0=3/4\Rightarrow\xi_*\approx0.032505121082825$; verify $\mathcal R_{ND}>1$ just left and $<1$ just right of each root; and verify $\mathcal R_{ND}<1$ for all $\xi\in\{1/100,1/5,3/5,9/10\}$ at $\delta=1$ (always-denominator slice).

**Self-test (already performed by auditor):** every `D[..., xi]` in M3/M8 acts on an expression that genuinely depends on $\xi$ (none is identically zero). The $\xi=0$ substitution in M9 reduces to $81\delta^3-72\delta^2=9\delta^2(9\delta-8)$, nonzero for generic $\delta$. The hand-expansion of M7's numerator yields leading coefficient `121` (= $11\cdot11$), so the Mathematica-derived cubic must match `121`, providing the independent cross-check against the F1 notes typo.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 229` and confirms the new `.wl` exists, exits 0, contains the M1-M10 checks, and that its independently-derived cubic (M7) and classifier (M4) agree with the SymPy residuals (engine cross-check).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering the selected-branch reduction, log-slope classifier, crossover cubic, monotonicity, threshold, and sampled crossover roots.
- deviation: none

## RESOLVED — F1 paper_misalignment (USER-AUTHORIZED 2026-06-02)

Direction: correct the notes to match the verified SymPy script (authoritative; script:96; the cubic is re-derived in-script and the notes' own next-line derivative `363 ξ²` already equals ∂(121 ξ³), confirming 121). Notes-only; Codex applies, Claude reviews. Do NOT change the script, paper.tex, or appendix.
- In `notes/stages/moving_throat_pde_stage229_..._sympy_audit.md`, the crossover-cubic leading coefficient `189 ξ³` (`189\xi^3`) → `121 ξ³` (`121\xi^3`). ONLY the leading 189→121 changes; the rest of the cubic (`297δξ² + 333δ²ξ + 81δ³ − 72δ²`) and the `363 ξ²` derivative stay as-is.
- Acceptance: stale `189\xi^3` gone; `121\xi^3` present. Append `## Applied: F1`.
