---
unit_id: 119
batch: IV.3
verifier_model: claude-opus-4-8-1m
verify_date: 2026-05-29T12:45:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 119

## Per-finding outcomes

### F1 — tautological_check (de-tautologization)

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py:53-56` — new assertion `tube-length law (rc -> rhat**2 link)` inserted after the existing `tube-length law` check: `L_sel.subs(rc, rhat**2) - sp.pi*a*sp.sqrt((1+rhat**2)/3)/2`.
- Mathematica `mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:60-63` — analogous `tube-length law (rC -> rHat^2 link)` assertion `(lSel /. rC -> rHat^2) - (Pi*a*Sqrt[(1 + rHat^2)/3])/2`. Section III `Clear`/`$Assumptions` (lines 52-53) were extended to include `rHat`.

**Assessment:**
Genuinely non-tautological. `L_sel` is the quantity *solved* from `kappa0 == (1+rc)/3` (so it carries the section's own `rc`), while the assertion substitutes the notes-mandated identification `r_c = rhat**2` and compares against an *independently written* closed form in `rhat`. The load-bearing point is that `rhat` is the SAME symbol object used in Sections I-II — verified: in the `.py`, `rhat` is declared once via `sp.symbols` (line 26) and never rebound; in the `.wl`, `rHat` is never bound to a value anywhere in the file (only `==`/`^2` usages), so it stays a free symbol shared across sections. The check therefore ties the L_W section to the dimensionless-ratio family parameter rather than restating the line-50 algebra; were the notes identification different (e.g. `r_c = rhat` or `2*rhat**2`), the residual would be non-zero. Not the X−X pitfall: the target is not defined as `L_sel` immediately prior.

Leak check (`.wl` `rHat` extension): clean. Section III sets `$Assumptions = Element[{a, lW, rC, rHat}, Reals] && a>0 && lW>0 && rC>0`, but Section IV (lines 67-70) does its own `Clear[...]` and *fully reassigns* `$Assumptions` to a fresh list (which independently and correctly declares `rHat` real). The Section III positivity block does not survive into Section IV, and `rHat` carries no stale value. No leak.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy `...sympy_audit.py:75-84` — two new assertions `T_m (+ branch) match` / `T_m (- branch) match` comparing `Tm_sol_plus`/`Tm_sol_minus` against `2*sqrt(2)*sqrt(K_s)*sqrt(Zq)/(J_s*sqrt(L_W)*c_s*sqrt(mu0)*(2*rhat ± sqrt(1+rhat**2)))`.
- Mathematica `...mathematica_audit.wl:88-99` — inline `stripCE[expr_] := expr /. ConditionalExpression[v_, _] :> v;` helper, then two assertions comparing `stripCE[tMPlus]`/`stripCE[tMMinus]` against the same factor-of-2-rearranged notes target.

**Assessment:**
The asserted quantities (`Tm_sol_plus`/`tMPlus` etc.) are obtained by `solve`ing the two-branch matching condition `ghat_expr == rhat ± sqrt(1+rhat^2)/2`, where `ghat_expr` itself descends from the substituted parent-action `Kq`/`gq` formulas — an independent chain. The comparison target is the notes §5 closed form `2 sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W)(2 rhat ± sqrt(1+rhat^2)))`, written out literally and NOT defined as the solved quantity immediately before. So this is solved-quantity-vs-independent-notes-form, not a self-comparison. The matching uses the Section II two-branch law `rhat ± sqrt(1+rhat^2)/2`, so the assertion closes the loop from parent-action substitution → branch matching → notes traction law. Confirmed via the printed solver outputs (sympy log lines 32-33; Mathematica log lines 38-39 showing the `ConditionalExpression` heads that `stripCE` removes) followed by residual `0`. Not tautological.

## Exec log assessment

**SymPy:** exit=0. 8 substantive zero-residual checks, all `= 0`:
- `tube-length law (rc -> rhat**2 link) = 0` (F1)
- `T_m (+ branch) match = 0`, `T_m (- branch) match = 0` (F2)
- plus the 5 pre-existing checks (dimensionless law, ±branch, tube-length law, ghat explicit). Closes with `All Stage 119 symbolic checks passed.`

**Mathematica:** exit=0. 8 `PASS:` lines:
- `PASS: tube-length law (rC -> rHat^2 link)` (F1)
- `PASS: T_m (+ branch) match`, `PASS: T_m (- branch) match` (F2)
- plus the 5 pre-existing. The `T_m` lines show `ConditionalExpression[...]` in the printed solutions; `stripCE` reduces them so the residual is `0`. Closes with `Stage 119 Mathematica audit passed.`

Engine independence confirmed: SymPy `expect_zero` uses `sp.simplify(sp.expand(...))`; Mathematica `expectZero` uses `FullSimplify[Together[Expand[...]], Assumptions -> $Assumptions]` and additionally must strip the `ConditionalExpression` head (the sympy solver did not emit a conditional, so it needs no analog). Different mechanics, not a transliteration.

**Output freshness:** confirmed. Script mtimes are 2026-05-27 16:35; the committed `.txt` outputs and the orchestrator's exec-log captures are 2026-05-29 12:36 — newer than the scripts, post-fix. Both committed `.txt` outputs contain the three new check lines (sympy output L20/L29/L30; Mathematica output L24-25/L35-38). Note: `stage_119_diff.patch` is 0 bytes and `git status` shows the files clean — the F1/F2 edits were already committed in `b4e02d8` ("batch IV.3 close"), so the working-tree diff is empty by design; this is not a missing-edit signal (the edits are present in the read-back script bodies and in the committed outputs).

## Material-change assessment

`material_change`: false. Both findings only ADD assertions (consistency-link and notes-form matches); no derived expression, constant, or solved result was altered. `L_sel`, `ghat_expr`, `Tm_sol_plus/minus` are computed exactly as before. Nothing downstream units could depend on changed.

## Side observations (non-blocking)

- The auditor's original report flagged the partial coverage of notes deliverables III and V (L_W upstream derivation, T_m boxed form); F1/F2 close the asserted-anchor gap. The deeper "kappa_c = 1/3 from D/N eigenvalue selection" derivation is still not exercised in-script (F1 option (a) was applied as the rc→rhat^2 family link, the auditor's preferred concrete suggestion), which is the intended scope of the low-severity finding — not a defect.
- `kappa0` remains a misleading name (it is `4 L_W^2/(pi^2 a^2)`, not the κ_c eigenvalue), as the auditor noted; renaming was out of directive scope and not required.

## Verdict justification

Both findings are `resolved`. F1's new assertion genuinely ties Section III's `rc` to the family-wide `rhat` symbol via the notes identification `r_c = rhat^2` (shared free symbol confirmed, no value binding, no `.wl` assumption leak because Section IV reassigns `$Assumptions`); it is non-tautological. F2's new assertions compare the independently-solved `T_m` branches against the literally-written notes §5 closed form (factor-of-2 rearranged), not a self-comparison, with the Mathematica side correctly stripping `ConditionalExpression` via the established `stripCE` idiom. Both engines exit 0 with 8 substantive checks each, use genuinely different simplification mechanics, and the saved outputs are fresh and contain the new lines. No regressions in the (committed) edits. Verdict: verified.
