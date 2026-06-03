---
unit_id: 241
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-03T08:46:38-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 241

Apply each finding below in order (F1, F2, and the F3 authorized notes edit). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F3 is a `paper_misalignment` the USER HAS RESOLVED — see `## RESOLVED — F3` below: notes typo, script canonical. Apply F3 by making ONLY the one authorized notes edit specified there (`193/369`→`125/369` on line 577). Do NOT change the script.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

The ONLY authorized prose/notes edit is the F3 notes correction in `## RESOLVED — F3`; do not touch `paper.tex` or any other prose document.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py:61` (and the `eps_sel` definition at line 58)

**Issue:** Line 58 defines `eps_sel = sp.simplify(1 - sp.Rational(3, 2) * varrho)`. Line 61 asserts `assert_zero(eps_sel - (1 - sp.Rational(3, 2) * varrho), "eps_* = 1 - 3 varrho / 2")`. The subtrahend is a verbatim copy of the `eps_sel` definition, so the residual is `X - X = 0` by construction and the assertion can never fail. The stage's actual claim (notes section 1, line 164) is that the selected-branch reduction `eps_* = 1 - 3varrho/2` is the *solution* of the Stage-240 support law `varrho = 2(1 - eps_*)/3`. The script must verify that derivation, not a self-assignment. The companion check on line 62 (`sigma_sel - (4/(3*varrho) - 2)`) is genuine — leave it alone.

**Required change:**
Keep the `eps_sel` definition at line 58 (it is the working value used downstream). Replace the line-61 assertion so the residual is built from the Stage-240 support law rather than from a literal copy of `1 - 3*varrho/2`. Acceptable forms (Codex picks one; you design the exact code):
- `assert_zero((2 * (1 - eps_star) / 3 - varrho).subs(eps_star, eps_sel), "eps_* solves Stage-240 support law varrho = 2(1-eps_*)/3")`, OR
- solve the support law for `eps_star` and assert the solution equals `eps_sel`: `assert_zero(sp.solve(2 * (1 - eps_star) / 3 - varrho, eps_star)[0] - eps_sel, ...)`.

The acceptance criterion is that the new residual passes through `varrho = 2*(1 - eps_*)/3` (notes line 164) and contains no second literal copy of `1 - 3*varrho/2`. Do not change the label text materially (a short clarifying suffix is fine).

**Verification command:**
The verifier will run `redteam exec-sympy 241` and confirm (a) the new check at/near line 61 references `2 * (1 - eps_star) / 3`, (b) there is no longer a bare `eps_sel - (1 - sp.Rational(3, 2) * varrho)` residual, and (c) the script exits 0 with `[ok]` for the selected-branch reduction.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py`
- summary: Replaced the tautological selected-branch assertion with a residual obtained by substituting `eps_sel` into the Stage-240 support law.
- deviation: none

## F2 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.wl` (NEW FILE — `.wl` lives in `mathematica/`, the `_mathematica_audit.wl` suffix is mandatory)

**Issue:** Stage 241 has only a SymPy script (card line 11: "Mathematica audit: none yet."). The unit is not a checkpoint and not status-only, so both engines are required. Every claim is verifiable with native Mathematica primitives (`Together`, `Factor`, `FullSimplify`, `D`, `Reduce`/`Solve`, exact rational arithmetic). Write a NEW independent-route Mathematica script.

**Independence requirement (do NOT transliterate the `.py`):**
The SymPy script imports the crossover conditions `eps_* = 1/(2+beta^2)` and `eps_* = beta/(1+beta+beta^2)` as literals (lines 121-132) and then checks consistency. The Mathematica route must instead DERIVE those crossovers: build the weights from `eps_*` and `beta`, form the differences `w_W - w_Lambda` and `w_Umag - w_Lambda`, and `Solve`/`Reduce` each `== 0` for `eps_*` to obtain the crossover, then substitute the support law to obtain `varrho_WLambda`, `varrho_ULambda`. This is a different decomposition (solve-for-the-zero) than the SymPy script's (assert-the-known-form), so it is a genuine second derivation, not a port. Use `Together`/`Factor`/`FullSimplify` (not a line-by-line echo of the SymPy `assert_zero` choreography).

**Claim manifest** (each must be independently verified; halt-and-fix if any fails):

- **M1** — Selected-branch reduction: solving `varrho == 2 (1 - epsStar)/3` for `epsStar` yields `epsStar == 1 - 3 varrho/2`; and with `epsStar = 1 - 3 varrho/2`, `sigma := 2 epsStar/(1 - epsStar)` simplifies to `4/(3 varrho) - 2`.
- **M2** — Twin-window inclusion: with `sigmaSel = 4/(3 varrho) - 2`, `FullSimplify[sigmaSel - (1/varrho - 2)] == 1/(3 varrho)` and `FullSimplify[(2/varrho - 2) - sigmaSel] == 2/(3 varrho)` (both positive for `0 < varrho < 2/3`).
- **M3** — Weight identities. With `N = 5(1 - epsStar)^2 + 6 epsStar^2 (1 + beta^2)`, `wLambda = epsStar^2 (1+beta^2)/N`, `wZ = (1 - 2 epsStar + (2+beta^2) epsStar^2)/N`, `wChi = 2 wZ`-candidate `= 2(1 - 2 epsStar + (2+beta^2) epsStar^2)/N`, `wW = epsStar(1-epsStar)/N`, `wUmag = beta epsStar(1-epsStar)/N`:
  - `wChi - 2 wZ == 0`;
  - `wZ - wLambda == (1 - epsStar)^2/N`;
  - `wZ - wW == (beta^2 epsStar^2 + 3(epsStar - 1/2)^2 + 1/4)/N`;
  - `wW - wUmag == epsStar(1 - epsStar)(1 - beta)/N`.
- **M4** — Crossover derivation (the independence-critical step): `Solve[wW - wLambda == 0, epsStar]` gives `epsStar -> 1/(2 + beta^2)` (drop the trivial `epsStar -> 0` root); `Solve[wUmag - wLambda == 0, epsStar]` gives `epsStar -> beta/(1 + beta + beta^2)` (drop `epsStar -> 0`). Then substituting the support law `epsStar = 1 - 3 varrho/2` and solving for `varrho` gives `varrhoWLambda = 2(1+beta^2)/(3(2+beta^2))` and `varrhoULambda = 2(1+beta^2)/(3(1+beta+beta^2))`.
- **M5** — Factorized sign laws over `D = 18 beta^2 varrho^2 - 24 beta^2 varrho + 8 beta^2 + 33 varrho^2 - 24 varrho + 8`, with `epsStar -> 1 - 3 varrho/2` substituted into the weights:
  - `FullSimplify[(wLambda - wW) - (2 - 3 varrho)(2 + beta^2)(varrhoWLambda - varrho)/D] == 0`;
  - `FullSimplify[(wLambda - wUmag) - (2 - 3 varrho)(1 + beta + beta^2)(varrhoULambda - varrho)/D] == 0`.
- **M6** — Threshold ordering: `FullSimplify[varrhoULambda - varrhoWLambda] == 2(1+beta^2)(1-beta)/(3(1+beta+beta^2)(2+beta^2))` (positive for `0<beta<1`); `FullSimplify[2/3 - varrhoULambda] == 2 beta/(3(1+beta+beta^2))` (positive).
- **M7** — Numeric windows: `varrhoWLambda /. beta -> 0` is `1/3`; `varrhoWLambda /. beta -> 2/11` is `125/369`; `varrhoULambda /. beta -> 0` is `2/3`; `varrhoULambda /. beta -> 2/11` is `250/441`. Also `D[varrhoWLambda, beta] == 4 beta/(3 (beta^2 + 2)^2)` and `D[varrhoULambda, beta] == -2(1 - beta^2)/(3(1 + beta + beta^2)^2)`.
- **M8** — Representative region orderings at `beta -> 1/10` and three sample `varrho` (one in each of `(0, varrhoWLambda)`, `(varrhoWLambda, varrhoULambda)`, `(varrhoULambda, 2/3)`): confirm the three orderings `wChi > wZ > wLambda > wW > wUmag` (Region I), `wChi > wZ > wW > wLambda > wUmag` (Region II), `wChi > wZ > wW > wUmag > wLambda` (Region III). Use exact rationals and verify each pairwise difference is positive (do not rely on machine floats).

Use exact `Exit[1]` on any failed check (e.g. a `checkZero[expr_, label_] := If[FullSimplify[expr] =!= 0, Print["FAIL ", label]; Exit[1], Print["[ok] ", label]]` pattern) so the script exits nonzero on failure. Preemptively strip `ConditionalExpression[0, ...]` from any `Solve`/`Reduce` output, and for the `beta -> 0` substitutions use the limit/direct substitution rather than a pole check.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 241` and confirms M1-M8 each emit `[ok]`, the script exits 0, and the file is NOT a transliteration of the `.py` (M4 must solve for the crossover, not import the literal).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.wl`
- summary: Added an independent Mathematica audit deriving the crossover roots with `Solve` and checking M1-M8 with exact symbolic and rational residuals.
- deviation: none

## F3 — paper_misalignment (notes_contradicts_script)

**Subtype:** notes_contradicts_script

**Paper side (notes prose):**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md:577` quote: "`\frac13 < \varrho_{W\Lambda} < \frac{193}{369}\approx 0.338753`."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py:184` quote: "`varrho_WLambda.subs(beta, beta_max) - sp.Rational(125, 369)`" (asserts `varrho_WLambda(beta=2/11) = 125/369`).

## RESOLVED — F3 (user direction 2026-06-02: correct notes to script)

Direction (a): notes numerator typo, script canonical. For `beta = 2/11`, `varrho_WLambda = 250/738 = 125/369 = 0.338753` — which equals the notes' OWN printed decimal `0.338753`; `193/369 = 0.523` is impossible here. **Authorized notes edit (Codex APPLIES; Claude reviews):**

In `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md` line 577, change `\frac{193}{369}` → `\frac{125}{369}` (leave the `\approx 0.338753` decimal). No script change. The new F2 `.wl` (M7) independently computes `varrhoWLambda /. beta -> 2/11 == 125/369` and cross-engine-corroborates.

## Applied: F3

- files_changed:
  - `notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md`
- summary: Corrected the authorized notes numerator from `193/369` to `125/369` while leaving the matching decimal unchanged.
- deviation: none
