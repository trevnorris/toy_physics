---
unit_id: 031
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 031

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:71-76`: Inserted three new lines between the existing `expect_zero("generic quotient/HF identity", ...)` (line 69) and the `print("dP0_-/dalpha =")` (line 75):
  ```
  dP_direct = sp.diff(P0_sel, alpha)
  dP_physical = beta0 * (ds_expected * lam_minus + s_minus**2) / lam_minus**2
  expect_zero("dP0_-/dalpha direct identity", sp.simplify(dP_direct - dP_physical))
  ```
  and changed the subsequent pprint argument from a re-typed expression to `sp.simplify(dP_physical)`. The generic-identity assertion remains intact.

**Assessment:**
Edit matches the directive verbatim (modulo the cosmetic pprint refactor, which uses the now-defined `dP_physical` instead of re-typing the closed form — equivalent in content). The new assertion is non-tautological: it invokes the physical `P0_sel = beta0 * s_minus / lam_minus` built from the explicit eigenvalue expressions, takes `sp.diff(...)` against `alpha`, and compares to the closed form `beta0*(ds_expected*lam_minus + s_minus**2)/lam_minus**2`. A sign error in `s_minus` or in `ds_expected` would cause this check to fail because `dP_direct` is computed independently of those symbols. The saved sympy output shows `dP0_-/dalpha direct identity = 0` at line 18, between `generic quotient/HF identity = 0` (line 17) and `dP0_-/dalpha =` (line 19), exactly as the directive's verification clause requires.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:89-94`: Replaced the tautological `(A*(A+DK) - alpha*T0).subs(alpha, alpha_crit)` check with a denesting-then-simplify check on `lam_minus.subs(alpha, alpha_crit)`. Codex did NOT use `sp.radsimp` as the directive suggested; instead it introduced `root_crit = A**2 x1 + (A+DK)**2 x0` and a custom `expr.replace(...)` that matches any `Pow` with exponent 1/2 whose base simplifies to `root_crit**2`, replacing the sqrt with `root_crit`. This deviation is acknowledged in the `## Applied: F2` block ("`sp.radsimp` did not reduce ... in this environment").
- `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:61`: Replaced `expectZero["alpha_crit condition", (a*(a + dK) - alpha*t0) /. alpha -> alphaCrit]` with `expectZero["lam_-(alphaCrit)", FullSimplify[lamMinus /. alpha -> alphaCrit, Assumptions -> $Assumptions]]` exactly as specified.

**Assessment:**
Both engines now assert the substantive physical claim `lambda_-(alpha_crit) = 0` (the selected eigenvalue collapses at the threshold), not the tautological rearrangement.

On the SymPy deviation: the custom denester is mathematically sound. The radical inside `lam_minus.subs(alpha, alpha_crit)` after `sp.together` is `sqrt(numerator/T0**2)` where `numerator = T0**2 * R(alpha_crit)**2`. Codex's `replace` only fires when the sqrt's base equals `root_crit**2` modulo `sp.simplify`; if the F3 perfect-square identity (verified separately at line 100) failed, the matcher would not fire and the residual would survive `sp.simplify` non-zero. The check is therefore genuinely substantive and not made-to-pass by the rewrite — it depends on the perfect-square structure that F3 establishes. Since all parameters are positive, `sqrt(root_crit**2) = root_crit` is valid (no |.| sign issue). The output shows `lam_-(alpha_crit) = 0` (sympy line 71) and `lam_-(alphaCrit) = 0 / PASS: lam_-(alphaCrit)` (mathematica lines 31-32), satisfying the directive's verification clause.

A genuine concern was whether the SymPy implementation might mask a sign error in `T0` or `alpha_crit`: if `T0 = (A+DK)*x0 - A*x1` (the typo example from the audit report), then `alpha_crit = A(A+DK)/T0_wrong`, `R(alpha_crit)` would no longer evaluate to `root_crit/T0_wrong`, the matcher would not match `root_crit**2`, the radical would survive, and `sp.simplify` of the unreduced rational expression would not produce 0. So the check is robust against that class of error.

### F3 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:97-101`: Replaced the hand-typed 9-term polynomial `radcrit = A**4*x0**2 + ...` with `radcrit_derived = sp.expand(T0**2 * (R**2).subs(alpha, alpha_crit))`. The assertion now expands the difference `radcrit_derived - (A**2 * x1 + (A + DK) ** 2 * x0) ** 2`.
- `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:64-65`: Replaced the hand-typed `radCrit = ...` with `radCritDerived = Expand[t0^2 * (r^2 /. alpha -> alphaCrit)]` and updated the assertion to `Expand[radCritDerived - (a^2*x1 + (a + dK)^2*x0)^2]`.

**Assessment:**
Edits match the directive verbatim in both engines. The LHS is now derived in-script from the physical radicand `R**2 = (DK + alpha*(x0-x1))**2 + 4*alpha**2*x0*x1`, evaluated at `alpha = alpha_crit`, and weighted by `T0**2`. The RHS is the closed-form perfect square. The check therefore anchors the threshold radical to the actual quadratic-formula discriminant, not to a hand-supplied polynomial. A typo in `R`, `T0`, or `alpha_crit` would now propagate to `radcrit_derived` and the difference would be non-zero. Output shows `threshold radical square identity = 0` in both engines (sympy line 72, mathematica line 33). Confirmed via grep that neither file contains the tail term `dK^4*x0^2` / `DK**4*x0**2` of the old hand-typed polynomial.

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The same Mathematica edit applied for F3 (lines 64-65). Directive explicitly states F4 is addressed by the F3 Mathematica patch; no additional edits required or applied.

**Assessment:**
The most direct evidence of transliteration (the identical hand-typed 9-term `radCrit` polynomial shared verbatim between `.py` and `.wl`) is eliminated. The Mathematica file now derives the perfect-square identity from its own symbolic `r^2 /. alpha -> alphaCrit` rather than copying SymPy's polynomial. PART II already does independent work via `D[p0Sel, alpha]` (which was unchanged). Verbatim transliteration in PARTs I, III, V remains, but that is intrinsic to the algebraic problem and is out of scope per the directive's minimal-change instruction. The grep above confirms `dK^4*x0^2` is gone from the Mathematica file.

## Exec log assessment

**SymPy:** exit=0 (inferred). `exec_logs/stage_031_sympy.log` is not present; per the verifier prompt, the saved output file `scripts/output/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.txt` (mtime 2026-05-21 17:22, newer than the script mtime 17:20) is the freshness record. Notable lines:
- `ds_-/dalpha exact formula = 0` (line 5)
- `generic quotient/HF identity = 0` (line 17)
- `dP0_-/dalpha direct identity = 0` (line 18) — new F1 assertion
- `lam_-(alpha_crit) = 0` (line 71) — new F2 assertion
- `threshold radical square identity = 0` (line 72) — new F3 assertion
- `lambda_- * lambda_+ - T0*(alpha_crit-alpha) = 0` (line 93)
- `STAGE 14 AUDIT COMPLETE` banner present (line 121); since `expect_zero` raises `AssertionError` on any non-zero residual, reaching the completion banner implies all assertions passed and the script exited 0.

**Mathematica:** exit=0 (inferred). `exec_logs/stage_031_mathematica.log` is not present; the saved output `mathematica/output/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.txt` (mtime 2026-05-21 17:22, newer than the .wl mtime 17:18) is the freshness record. Notable lines:
- `PASS: ds_-/dalpha exact formula` (line 6)
- `PASS: exact monotonicity identity` (line 13) — unchanged, PART II already correct per F1's source
- `PASS: lam_-(alphaCrit)` (line 32) — new F2 assertion
- `PASS: threshold radical square identity` (line 34) — F3 derived form
- `PASS: lambda_- * lambda_+ - T0*(alpha_crit-alpha)` (line 42)
- `STAGE 14 AUDIT COMPLETE` banner present (line 48); 10 `PASS:` lines, 0 `FAIL:` lines. `Exit[0]` is in the script at line 83, and `expectZero` calls `fail[...]; Exit[1]` on any non-zero residual, so reaching the banner indicates exit 0.

**Output freshness:** Confirmed. Script mtimes are 17:18 (.wl) and 17:20 (.py); output mtimes are both 17:22, post-fix. Diff patch mtime also 17:22.

## Material-change assessment

`material_change`: false.

The fixes tighten existing checks and replace hand-typed expressions with derived ones, but they do not alter any quoted result: `alpha_crit = A(A+DK)/T0`, `lambda_+^(crit) = ((A+DK)^2 x0 + A^2 x1)/T0`, `ds_-/dalpha`, `dP0_-/dalpha`, and `P0_-` are all unchanged in the saved outputs. Downstream units that depend on `alpha_crit`, the perfect-square identity, or the prefactor derivative will see the same symbolic forms.

## Side observations (non-blocking)

- The F2 SymPy implementation is more elaborate than the directive's `sp.radsimp(...)` suggestion. Codex's `## Applied: F2` `deviation` note explains this is environmental (sp.radsimp did not reduce). The custom `Pow.replace` with a `root_crit**2` matcher is a localized denester. It is sound under positivity assumptions; the only risk is a future maintainer changing `root_crit` to a different (but equivalent) form that fails the `expr.base - root_crit**2` simplifier — at which point the matcher would silently no-op and the residual sqrt would survive. This is a maintainability concern, not a correctness concern, and not part of any finding.
- The pretty-printed `lambda_+^(crit)` in the SymPy output (lines 78-88) is not auto-denested by `sp.simplify`; the Mathematica equivalent at line 36 is fully simplified to `((a + dK)^2*x0 + a^2*x1)/(dK*x0 + a*(x0 + x1))`. The two are algebraically equal (the original audit verified this), so this is purely a cosmetic engine difference and is not part of any finding.
- PARTs I, III, V of the Mathematica file remain structurally parallel to the SymPy file. The directive scoped F4 to the F3 Mathematica edit only, so this is expected and not a defect.

## Verdict justification

All four findings are `resolved`. F1 adds the substantive direct-derivative check the audit specified; F2 replaces the tautological `alpha_crit` check with `lambda_-(alpha_crit) = 0` in both engines (the SymPy implementation deviates from the suggested `sp.radsimp` to a custom radical denester but remains substantively correct and non-tautological); F3 derives the threshold radical polynomial from the physical `R**2(alpha_crit) * T0**2` in both engines; F4 is structurally addressed by the F3 Mathematica edit. Saved outputs confirm all new assertions evaluate to 0, both scripts reach their completion banners (implying exit 0), and no hardcoded 9-term polynomial remains in either file. No regressions introduced by the diff. No findings were blocked or skipped.
