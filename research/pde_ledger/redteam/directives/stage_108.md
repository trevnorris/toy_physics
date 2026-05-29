---
unit_id: 108
batch: IV.2
created_at: 2026-05-28T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 108

Apply each non-paper_misalignment finding below (F2, F3, F4) in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For the `paper_misalignment` finding F1, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" F1 unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_108.tex:22-25` quote: "Check pure scale and pure argument deformations separately. / Check Robin and standalone mixed-pole limits before imposing compensation. / Check that the compensated branch preserves the even coefficients as well as the odd normalization."
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:27` quote: "Stages 107--113: low-frequency outlet deformations and the compensated Robin--mixed branch."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage108_robustness_classes.md` (covers ONLY Classes A/B/C and the preservation submanifold; no Robin / mixed-pole / compensated material).

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py` (no Robin / standalone-mixed-pole / compensated-Robin–mixed checks anywhere).
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl` (same).

## Resolve before fix_loop

The stage-108 card's `Checks` items #2 and #3 require Robin, standalone mixed-pole, and compensated-Robin–mixed verification (`chi_Q^R = 3/(3-rho_R)`; `kappa_W=-1/9` then `sigma_W=0`; `chi_Q^hyb = (1-9 sigma_W gamma_W)/(1-sigma_W)` preserved iff `gamma_W=1/9`). Neither stage-108 script tests any of these, and the stage-108 notes file does not mention them. The appendix says this material spans "Stages 107--113." Question: do the Robin/mixed-pole/compensated checks belong to stage 108, or to a sibling stage (≈109–113)?

Possible directions (the user picks one):
- (a) They belong to a sibling stage → the stage-108 card's `Checks` #2/#3 are block-level over-scoping; re-scope/annotate the card (paper-side edit), no script change. Confirm which sibling stage's scripts actually verify Robin/mixed/compensated.
- (b) They belong to stage 108 → extend BOTH scripts with Robin (`chi_Q^R = 3/(3-rho_R)`), standalone mixed-pole no-go (`kappa_W=-1/9` then `sigma_W=0`), and compensated-Robin–mixed (`chi_Q^hyb`, even solutions `(rho_R,kappa_W)=(sigma_W,0)` and `(4 sigma_W,1/3)`, `gamma_W=1/9`) checks; a follow-up directive will then specify the claim manifest.
- (c) Shared block checklist deliberately repeated on every stage in 107–124 → leave both card and scripts as-is and record the convention; no edit.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:83-112`

**Issue:** The `.wl` Class D block ports the `.py`'s exact `Series`→`Coefficient`→`Solve`→assemble route (and copies the inline comments verbatim). To satisfy the second-engine independent-derivation requirement, the `.wl` must derive the same Class D result by a structurally different route and confirm it matches the existing series-route `chiGen`.

**Required change:**
After the existing Class D block computes `chiGen` (currently around line 100), insert an INDEPENDENT-route cross-check that does NOT use `Series` of the rational function. Build the raw expansion coefficients of `lambdaGen` directly as symbolic expressions in `{sNorm, beta, sigma0, sigma2, sigma4, sigma5}` and impose the even-fingerprint constraints as plain equations:

```wolfram
(* Independent route: build chiGen from raw L-coefficients without Series. *)
L0raw = -3*sNorm + sigma0;
L2raw = sNorm*beta^2/3 + sigma2;
L4raw = sNorm*beta^4/9 + sigma4;
L5raw = sNorm*beta^5/9 + sigma5;
solAlt = First[Solve[{-L2raw/L0raw == 1/9, L2raw^2/L0raw^2 - L4raw/L0raw == 4/81},
                     {sigma2, sigma4}, Reals]];
chiGenAlt = FullSimplify[(27*(-L5raw/L0raw)) /. solAlt, Assumptions -> $Assumptions];
expectZero["chiGen independent-route agreement", chiGenAlt - chiGen];
```

Place this block AFTER the line that defines `chiGen` (so `chiGen` is in scope) and BEFORE the existing `sigma5PresGen` solve at line 102, OR immediately after `chiGen` is printed at line 101 — either ordering is acceptable as long as `chiGen` is already bound and the new `expectZero` runs. Do not delete or alter the existing series-route Class D lines; this is an additive independent-route confirmation.

**Self-test (verified by auditor):**
- `L0raw,L2raw,L4raw,L5raw` are the literal coefficients of `lambdaGen = sNorm*(lambdaOut /. z->beta z) + sigma0 + sigma2 z^2 + sigma4 z^4 + I sigma5 z^5`: constant `-3 sNorm + sigma0`; `z^2`: `sNorm beta^2/3 + sigma2`; `z^4`: `sNorm beta^4/9 + sigma4`; `z^5` imaginary part: `sNorm beta^5/9 + sigma5`. These match the appendix expansion.
- Solving the two even-match equations for `sigma2,sigma4` and forming `27*(-L5raw/L0raw)` yields `chiGenAlt = 3(sNorm beta^5 + 9 sigma5)/(3 sNorm - sigma0)`, identical to the series-route `chiGen` (see Mathematica output line 36). So `chiGenAlt - chiGen` simplifies to 0 — falsifiable: if either route had a wrong coefficient the residual would be nonzero.
- No `D[...]`/derivative introduced; no `Series` in the new route. The two routes share only the physical premise (`lambdaOut` literal), not the algorithm.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 108` and confirm the new `chiGen independent-route agreement` line appears, prints `0`, `PASS`, and the script exits 0.

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:30-32`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:35-37`

**Issue:** Class A asserts `Y_scale - Y_can == 0`, but `Y_scale = (-3*S)/(S*Lambda_out)` reduces to `-3/Lambda_out = Y_can` by cancellation of `S` before the series, so the assertion cannot fail for any `Lambda_out`. Anchor it to the paper's literal canonical fingerprint (appendix eq:app-part04-Yout-dtn) so it is falsifiable.

**Required change:**

SymPy — replace lines 30-32:
```python
# Class A: pure scale
Y_scale = sp.series((-3*S)/(S*Lambda_out), z, 0, 6).removeO()
Y_can = sp.series((-3)/Lambda_out, z, 0, 6).removeO()
expect_zero('pure scale invariance', Y_scale - Y_can)
```
with:
```python
# Class A: pure scale -- anchor normalized scaled response to the literal canonical fingerprint
# (appendix eq:app-part04-Yout-dtn) so the check is falsifiable, not a bare S-cancellation.
Y_scale = sp.series((-3*S)/(S*Lambda_out), z, 0, 6).removeO()
Y_can_literal = 1 + z**2/sp.Integer(9) + 4*z**4/sp.Integer(81) + I*z**5/sp.Integer(27)
expect_zero('pure scale invariance', sp.expand(Y_scale) - Y_can_literal)
```

Mathematica — replace lines 35-37:
```wolfram
yScale = Expand[Normal[Series[(-3*sNorm)/(sNorm*lambdaOut), {z, 0, 5}]]];
yCan = Expand[Normal[Series[(-3)/lambdaOut, {z, 0, 5}]]];
expectZero["pure scale invariance", yScale - yCan];
```
with:
```wolfram
yScale = Expand[Normal[Series[(-3*sNorm)/(sNorm*lambdaOut), {z, 0, 5}]]];
(* anchor to literal canonical fingerprint (appendix eq:app-part04-Yout-dtn) *)
yCanLiteral = 1 + z^2/9 + (4*z^4)/81 + (I/27)*z^5;
expectZero["pure scale invariance", yScale - yCanLiteral];
```

**Self-test (verified by auditor):**
- `-3/Lambda_out` expanded to O(z^6) equals `1 + z^2/9 + 4 z^4/81 + i z^5/27` (hand-expansion of `1/(1 - z^2/9 - z^4/27 - i z^5/27)`), so `Y_scale - Y_can_literal` reduces to 0. Falsifiable: a wrong `Lambda_out` coefficient changes `Y_scale` and the residual becomes nonzero.
- The literal matches appendix line 321 exactly; no new constant introduced.

**Verification command:**
`redteam exec-sympy 108` and `redteam exec-mathematica 108`; the `pure scale invariance` line still prints `0`/`PASS` and the scripts exit 0.

## F4 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:68-70,102`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:78-81,108`

**Issue:** Two solve-then-resubstitute round-trips cannot fail. Replace the Class C round-trip with a paper-anchored locus check; demote the redundant Class D round-trip to a print (the Class D submanifold anchor already exists and stays).

**Required change:**

SymPy — replace lines 68-70:
```python
chi_pres = sp.solve(sp.Eq(chi_add, 1), Sigma5)[0]
print('Sigma5 preservation locus =', sp.simplify(chi_pres))
expect_zero('preservation locus check', chi_add.subs(Sigma5, chi_pres) - 1)
```
with:
```python
chi_pres = sp.solve(sp.Eq(chi_add, 1), Sigma5)[0]
print('Sigma5 preservation locus =', sp.simplify(chi_pres))
# Anchor the Class C locus value to the notes submanifold reduced to beta=1 (Sigma5 = -Sigma0/27),
# instead of re-substituting the solved value back into the same equation.
expect_zero('Sigma5 locus (Class C) = -Sigma0/27', sp.simplify(chi_pres) - (-Sigma0/sp.Integer(27)))
```

SymPy — at line 102, replace:
```python
expect_zero('general preservation locus check', chi_gen.subs(Sigma5, chi_pres_gen) - 1)
```
with:
```python
# (Round-trip check demoted to a print; the submanifold anchor above is the load-bearing test.)
print('general preservation locus check =', sp.simplify(chi_gen.subs(Sigma5, chi_pres_gen) - 1))
```

Mathematica — at lines 78-81, replace:
```wolfram
sigma5Pres = FullSimplify[sigma5 /. First[Solve[chiAdd == 1, sigma5, Reals]], Assumptions -> $Assumptions];
Print["Sigma5 preservation locus = ", fmt[sigma5Pres]];
expectZero["Sigma5 preservation locus + sigma0/27", sigma5Pres + sigma0/27];
expectZero["preservation locus check", (chiAdd /. sigma5 -> sigma5Pres) - 1];
```
with:
```wolfram
sigma5Pres = FullSimplify[sigma5 /. First[Solve[chiAdd == 1, sigma5, Reals]], Assumptions -> $Assumptions];
Print["Sigma5 preservation locus = ", fmt[sigma5Pres]];
expectZero["Sigma5 locus (Class C) = -sigma0/27", sigma5Pres - (-sigma0/27)];
```
(Note: the existing `Sigma5 preservation locus + sigma0/27` assertion at line 80 is already a valid anchor and is equivalent to the new one; if Codex prefers, it may instead keep line 80 and only DELETE the round-trip line 81. Either keeping the anchored line 80 + deleting line 81, or the replacement above, is acceptable — the goal is: no round-trip, and the locus value `-sigma0/27` remains asserted.)

Mathematica — at line 108, replace:
```wolfram
expectZero["general preservation locus check", (chiGen /. sigma5 -> sigma5PresGen) - 1];
```
with:
```wolfram
Print["general preservation locus check = ", fmt[FullSimplify[(chiGen /. sigma5 -> sigma5PresGen) - 1, Assumptions -> $Assumptions]]];
```

**Self-test (verified by auditor):**
- Class C: `solve(chi_add==1, Sigma5)` gives `Sigma5 = -Sigma0/27`, so `chi_pres - (-Sigma0/27) == 0`. Falsifiable: if `chi_add`'s closed form were wrong the locus would shift. Matches notes submanifold `S(1-beta^5)/9 - Sigma0/27` at beta=1 → `-Sigma0/27`. No new constant.
- Class D: the submanifold anchor `chi_pres_gen - (S(1-beta^5)/9 - Sigma0/27)` (sympy 98-101 / wl 104-107) is left untouched and remains the load-bearing Class D check, so demoting the round-trip removes no real coverage.
- No `diff`/`D`, no parity-sensitive integral; both edited assertions/prints depend only on already-bound symbols.

**Verification command:**
`redteam exec-sympy 108` and `redteam exec-mathematica 108`; confirm `Sigma5 locus (Class C) = -Sigma0/27` (or the equivalent `... + sigma0/27`) asserts and prints `0`/`PASS`, the demoted Class D line prints `0` (no longer an assertion), the submanifold anchor still `PASS`es, and the scripts exit 0.
