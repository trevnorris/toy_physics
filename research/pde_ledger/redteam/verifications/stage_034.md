---
unit_id: 034
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 034

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py:83-92` — the two prior `expect_zero` calls on `gBreq_sq_over_varpi2 - (alpha_x - alpha_mix)` and `alpha_mix + gBreq_sq_over_varpi2 - alpha_x` are gone. In their place is a fresh `sp.solve` of the loading equation against `alpha_lam` (the lambda-form closed expression) producing `gBreq_lambda`, followed by a single `expect_zero("lambda-form vs x-form support loading", gBreq_lambda.subs(lam, A - x) - gBreq_sq_over_varpi2)`.
- `mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl:91-97` — the two prior `expectZero` calls on `gBReqSqOverVarpi2 - (alphaX - alphaMix)` and `alphaMix + gBReqSqOverVarpi2 - alphaX` are gone. In their place: `gBSolutionLambda = Solve[alphaMix + gBsqOverVarpi2 == alphaLambda, gBsqOverVarpi2, Reals]`, a `Length` guard via `fail`, `gBReqLambda = FullSimplify[gBsqOverVarpi2 /. First[gBSolutionLambda], …]`, and a single `expectZero["lambda-form vs x-form support loading", (gBReqLambda /. lambda -> A - x) - gBReqSqOverVarpi2]`.
- The `Print["alpha_mix = …"]` and `Print["g_B,req^2 / varpi^2 = …"]` informational prints are preserved, as the directive instructed.
- Diff is restricted to the two named line ranges; no collateral edits to docstrings, banners, variable names, `$Assumptions`, or other assertions.

**Assessment:**
Both engines now drop the two tautological residual checks (which were guaranteed zero because `gBreq` had been constructed as the CAS solution of that exact linear equation) and replace them with one cross-form check. The new check constructs `gBreq_lambda` from the closed-form `alpha_lam` and then verifies that substituting `lambda -> A - x` recovers the x-form solution. The residual does NOT collapse by construction: it would be non-zero if `alpha_lam` were defined incorrectly or if the closed-form `alpha_x` did not algebraically equal `alpha_lam.subs(lam, A-x)`. This matches the directive verbatim — the SymPy block is character-for-character what the directive prescribed, and the Mathematica block reproduces the prescribed `Solve` / `Length` guard / `FullSimplify` / `expectZero` pattern exactly.

Both saved outputs show the new line:
- SymPy output line 56: `lambda-form vs x-form support loading = 0`
- Mathematica output lines 22-23: `lambda-form vs x-form support loading = 0` followed by `PASS: lambda-form vs x-form support loading`

And the two prior tautological-check lines (`solved g_B,req^2/varpi^2 - (alpha(x) - alpha_mix) = 0`, `alpha_mix + g_B,req^2/varpi^2 - alpha(x) = 0` and the analogous Mathematica `PASS:` lines) are absent. The pre-existing A1-A5 / A8-A12 assertions (the lambda-to-x equivalence of `alpha`, `s_-`, `N_-`, derivative identity, reciprocity) still appear and still pass.

The Mathematica file's second `$Assumptions` block (lines 79-82) does not include `lambda`, but the directive explicitly flagged and accepted this — the substitution `lambda -> A - x` eliminates `lambda` before `FullSimplify` runs, and the saved output confirms the residual reduces to exactly 0.

## Exec log assessment

**SymPy:** exit=n/a (no exec log captured at `redteam/exec_logs/stage_034_sympy.log`). Substituting saved output (`scripts/output/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.txt`, mtime 2026-05-21 17:33, newer than the script mtime). Notable lines:
- `alpha(lambda=A-x) - alpha(x) = 0` (line 26)
- `s_x * d alpha / dx - 1 = 0` (line 44)
- `lambda-form vs x-form support loading = 0` (line 56)
- `All Stage 17 checks passed.` (line 57)

**Mathematica:** exit=n/a (no exec log captured at `redteam/exec_logs/stage_034_mathematica.log`). Saved output (`mathematica/output/…audit.txt`, mtime 2026-05-21 17:33, newer than the script mtime). Notable lines:
- `PASS: alpha(lambda=A-x) - alpha(x)` (line 10)
- `PASS: s_x * d alpha / dx - 1` (line 19)
- `PASS: lambda-form vs x-form support loading` (line 23)
- `Stage 034 Mathematica audit passed.` (line 25)

Since `expect_zero` raises `AssertionError` and `expectZero` exits 1 via `fail` on any non-zero residual, the fact that both files completed through the final "passed" banner implies exit 0. The orchestrator simply did not write the per-engine log files, but the saved txt outputs are post-fix (mtimes 2026-05-21 17:33:45 and 17:33:50, both newer than the script mtime 2026-05-21 17:32:47).

**Output freshness:** confirmed. Script mtime: 1779406367 (2026-05-21 17:32:47 UTC). SymPy output mtime: 1779406425 (~58 s later). Mathematica output mtime: 1779406430 (~63 s later). Outputs were regenerated after the edits landed.

## Material-change assessment

`material_change`: false.

The closed-form derivations `alpha_x`, `s_x`, `N_x`, `dalpha_dx`, and the symbolic `gBreq_sq_over_varpi2 = alpha_x - alpha_mix` are unchanged in both engines. Only the trailing pair of residual checks was swapped out for one cross-form check. No symbol bindings, no variable renames, no `$Assumptions` modifications, no docstring or banner edits. Downstream units that consume the softening-depth closed forms (or the support-loading expression `gBreq`) see identical symbolic results.

## Side observations (non-blocking)

- The new SymPy/Mathematica check is algebraically reducible to the already-asserted A1/A8 identity (`alpha_lam.subs(lam, A-x) - alpha_x == 0`): the residual `gBreq_lambda.subs(lam, A-x) - gBreq_sq_over_varpi2` expands to `(alpha_lam - alpha_mix).subs(lam, A-x) - (alpha_x - alpha_mix) = alpha_lam.subs(lam, A-x) - alpha_x`. So the new check is not adversarially independent of A1/A8 — it re-exercises the same content via a different code path. The directive's auditor proposed this exact wording; it is non-tautological in the directive's sense (the residual is not zero by construction), so I am not blocking on this. Flagging in case a later auditor wants to add a numerically-substituted point check on top.
- The CAS-solved sympy output (lines 51-55 of the saved output) prints `gBreq_sq_over_varpi2` in an expanded numerator form rather than the simpler `alpha_x - alpha_mix`. This is just a `sp.simplify` quirk; the residual at line 56 is still 0. Not a finding.

## Verdict justification

Codex applied F1 verbatim to both engines, the diff is limited to the two prescribed line ranges, the two tautological assertions are gone, and both saved outputs (regenerated post-edit) show the new `lambda-form vs x-form support loading = 0` line passing. No collateral edits, no regressions, no symbol changes that would propagate to downstream units. Verdict: verified.
