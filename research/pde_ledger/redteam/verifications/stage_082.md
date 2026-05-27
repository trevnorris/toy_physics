---
unit_id: 082
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 082 (v2 audit cycle)

This verification covers iteration 2 of the audit cycle for stage 082, addressing the three findings (F1, F2, F3) raised by the 2026-05-27 auditor report. The prior iteration's verification (4 findings: tautological_check, insufficient_verification, mathematica_transliteration, hardcoded_result) is superseded by this report.

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim, severity medium)

**Classification:** resolved

**What changed:**
A new closed-form pin block for `zeta_phys` was inserted in both engines, after the F3 derivative-factorization block and before the Family-1 strength prints:

- SymPy: `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:94-138` — symbolic `Omega_Pe_expr = pi*Pe*(2*Pe*exp(Pe) + pi) / ((4*Pe^2 + pi^2)*(exp(Pe) - 1))`, `zeta_phys_closed = Omega_Pe^2 * (kappa + pi^2/4) / (kappa + y^2)`, an `Omega_Pe -> pi/2` limit assertion via `sp.limit(..., Pe, sp.oo)`, an mpmath bisection for `y_F1` (smallest positive root of `y*tan(y) = 37` on `(1.5, 1.55)`), evaluation of `zeta_phys_F1_limit = (pi^2/4) * (kappa_F1 + pi^2/4) / (kappa_F1 + y_F1^2)` at Family-1 `(eta_F1, kappa_F1) = (37, 12321/5)`, and an `assert` comparing against `zeta_max_F1_reference = 2.467529229455835` to `< 1e-10`.
- Mathematica: `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:88-113` — parallel block using `Limit[OmegaPe, peSym -> Infinity]`, `FindRoot[ySym*Tan[ySym] - 37 == 0, {ySym, 1.527}, WorkingPrecision -> 30]`, evaluation at `kappaF1 = 12321/5`, and an `If[TrueQ[diff < 10^-10], pass[...], fail[...]]` gate against the same reference constant.

Banner also relabelled `STAGE 65` -> `STAGE 082` in both scripts (a long-standing cosmetic side observation from the prior verification — now resolved as a collateral edit).

**Assessment:**
The edit closes the F1 missing-claim cleanly. The closed form for `zeta_phys` is now instantiated as written in eq. 082-zeta-phys (`Omega_Pe^2 (kappa + pi^2/4) / (kappa + y(eta)^2)`), with the transport overlap `Omega_Pe(Pe)` symbolically defined per the paper's transport block and verified to limit to `pi/2`. The Family-1 numerical pair `(37, 12321/5)` is instantiated, the Robin-root condition `y*tan(y) = eta` is solved numerically, and the Pe->oo limit is cross-checked against stage 084's pre-computed `zeta_max^(F1)` reference. The deviation from `sp.nsolve` to `mpmath.findroot(..., (1.5, 1.55), solver='bisect')` is well-justified: SymPy's Newton iteration is unstable near `pi/2` (where `tan` diverges) and can jump to a far-away root; bracketed bisection on `(1.5, 1.55)` is correct because the smallest positive root of `y tan y = 37` lies in that bracket (`y tan y` is monotonically increasing from 0 at `y = 0` through arbitrarily large values as `y -> pi/2^-`, so a unique bracketed root exists). Mathematica's `FindRoot[..., {y, 1.527}, WorkingPrecision -> 30]` is independently stable.

The two engines independently converge to the same `y_F1` to 30 digits (sympy: `1.52948248371469964992710762240`; mathematica: `1.5294824837146996499271076224024730...`) and the same `zeta_phys_F1_limit` (sympy: `2.4675292294560122333`; mathematica: `2.46752922945601223332958450157053...`). Both report `|zeta_phys(F1, Pe->oo) - zeta_max^(F1)| ≈ 1.77e-13`, well under the `1e-10` threshold. Engine cross-agreement to ~13 digits is the substantive cross-check.

Non-tautological: the `Omega_Pe -> pi/2` limit assertion involves a non-trivial transcendental limit of a rational-in-`exp(Pe)` expression that would fail under any sign error in the numerator or denominator; the Family-1 pin compares against a value computed in a separate stage's script via a different code path.

### F2 — paper_misalignment (script_missing_paper_claim, severity low)

**Classification:** resolved

**What changed:**
Per the directive, F2 is subsumed by F1 (no separate script edit beyond F1's block).

Sub-question (i): the Family-1 numerical pair `(eta, kappa) = (37, 12321/5)` is now instantiated inside F1's block (sympy line 123 `kappa_F1 = sp.Rational(12321, 5)` together with `y_F1` from the `y tan y = 37` solve; mathematica line 104 `kappaF1 = 12321/5` and `yF1` from `FindRoot[ySym*Tan[ySym] - 37 == 0, ...]`).

Sub-question (ii): the Upsilon_w convention (`100` vs `117`) is resolved upstream by the Q2 paper-side update. `paper/stages/stage_075.tex:7` now reads `\Upsilon_w=100\Theta_w` (verified), so the script's existing `Upsilon_w = 100 Theta_w` is consistent with the canonical paper text. The carry-forward TODO comment at sympy 142-146 / mathematica 115-118 was rewritten from a provenance TODO into a citation block referencing the stage 075 update.

**Assessment:**
Both sub-questions are closed by the upstream changes (F1 for (i); Q2 paper edit for (ii)). No separate script edit was required, which matches the directive's "Required change: none beyond F1". The comment rewrite is a documentation clarification, not a new claim.

### F3 — insufficient_verification (severity low)

**Classification:** resolved

**What changed:**
- SymPy: `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:77-92` — removed the tautological `dR_quad/dzeta_phys + 1 == 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-) == 0` checks. Replaced with `sp.fraction(sp.together(sp.simplify(sp.diff(zeta_req, Pi_tr))))` to factor `dzeta_req/dPi_tr`, plus two new `expect_zero` calls asserting `numerator - C_mix*(1 - eps_blk) == 0` and `denominator - (C_mix - eps_blk*(2*C_mix - Pi_tr))**2 == 0`.
- Mathematica: `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:71-86` — parallel replacement using `Numerator[Together[...]]` / `Denominator[Together[...]]` with the same two `expectZero` factorization checks under `$Assumptions`.

**Assessment:**
The new assertions are non-tautological and exercise the strict-positivity content of notes section 4. The numerator-equals-`C_mix*(1 - eps_blk)` check would fail under any algebra error in the quotient rule (e.g., a sign flip in `eps_blk*(2*C_mix - Pi_tr)` or a mis-collected `eps_blk*C_mix` term). The denominator-equals-`(C_mix - eps_blk*(2*C_mix - Pi_tr))^2` check verifies that the denominator of `dzeta_req/dPi_tr` factors cleanly into the squared form of `zeta_req`'s own denominator — a non-trivial property. Together, the two assertions establish that `dzeta_req/dPi_tr` is sign-controlled on the physical branch (where `eps_blk ∈ [0, 1)`): the numerator is `> 0` and the denominator is a square `> 0`, hence the derivative is strictly positive.

Both engines independently confirm the factorization (sympy output lines 36-37; mathematica output lines 21-24 with PASS lines). The original two tautological checks are entirely gone from both scripts (verified by searching for `dR_quad/dzeta_phys + 1` and finding no occurrences in the post-fix scripts).

## Exec log assessment

**SymPy:** exit=0. The manifest records `last_codex_exit: 0`, and the prompt confirms both engines exit 0. The standalone `redteam/exec_logs/stage_082_sympy.log` file is absent (orchestrator did not emit a separate log file), but the canonical script output `scripts/output/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.txt` was freshly regenerated post-fix (mtime `2026-05-27 02:17`, vs. script mtime `2026-05-27 02:17`). Notable lines from the regenerated output:

- `Omega_Pe -> pi/2 as Pe -> oo = 0` (the limit-assertion residual)
- `y_F1 (root of y tan y = 37, smallest positive) = 1.52948248371469964992710762240`
- `zeta_phys(Pe->oo, kappa_F1, y_F1) = 2.4675292294560122333`
- `|zeta_phys(F1, Pe->oo) - zeta_max^(F1)| = 1.77233329100813225006336971399E-13`
- `PASS: Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10.`
- `numerator(d zeta_req/d Pi_tr) - C_mix*(1 - eps_blk) = 0`
- `denominator(d zeta_req/d Pi_tr) - (C_mix - eps_blk*(2*C_mix - Pi_tr))**2 = 0`

All `expect_zero` residuals display `= 0`; the assert for the Family-1 pin would have raised AssertionError on any failure but the printed PASS line confirms it passed.

**Mathematica:** exit=0. The manifest records `last_codex_exit: 0`. The canonical output `mathematica/output/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.txt` was freshly regenerated (mtime `2026-05-27 02:19`, vs. script mtime `2026-05-27 02:15`). Notable lines:

- `PASS: Omega_Pe -> Pi/2 as Pe -> oo` (line 27)
- `y_F1 (root of y tan y = 37, smallest positive) = 1.5294824837146996499271076224024730266867806218758702464039`30.` (line 28)
- `zeta_phys(Pe->oo, kappa_F1, y_F1) = 2.46752922945601223332958450157053542039`20.` (line 29)
- `|zeta_phys(F1, Pe->oo) - zeta_max^(F1)| = 1.772333295845015705...e-13` (line 30)
- `PASS: Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10` (line 31)
- `PASS: numerator(d zeta_req/d Pi_tr) - C_mix*(1 - eps_blk)` (line 22)
- `PASS: denominator(d zeta_req/d Pi_tr) - (C_mix - eps_blk*(2*C_mix - Pi_tr))^2` (line 24)

Absence of `FAIL:` lines (the `fail[...]` function calls `Exit[1]` on any check failure) and the script's terminating `Exit[0]` confirm a clean run.

**Output freshness:** confirmed. SymPy output mtime `02:17`, script mtime `02:17` (both updated in same batch); Mathematica output mtime `02:19`, script mtime `02:15` (output is newer). Both `.txt` files show the post-fix `STAGE 082` banner and the new F1/F3 blocks.

**Engine cross-check:** `y_F1`, `zeta_phys_F1_limit`, and the diff to the upstream reference all agree between SymPy and Mathematica to the printed precision (~13 digits for the diff; 20+ digits for the values). Engines independently confirm the F1 pin.

## Material-change assessment

`material_change`: false.

The F1 block computes the Pe->oo Family-1 limit of `zeta_phys` and cross-checks it against the *existing* upstream constant `zeta_max^(F1) ≈ 2.467529229455835` (sourced from stage 084 .wl, which itself was verified at v2). Stage 082 consumes this constant; it does not redefine it. The F3 change replaces tautological assertions with non-tautological factorization checks of the same `dzeta_req/dPi_tr` expression — no new symbolic result is introduced; the printed forms of `zeta_req`, `Q`, `Pi_suff`, `Pi_fail`, `R_quad`, `Xi_F1` are all unchanged. The banner relabel and the rewritten provenance comment are documentation-only. Downstream stages > 082 do not depend on any new value from this stage and do not need re-audit on the basis of stage 082's content alone.

## Side observations (non-blocking)

1. The directive's prose at line 121 mentioned the expected `y_F1 ~ 1.5269...`, but the actual smallest positive root of `y tan y = 37` is `~1.52948`. Both engines independently confirm `1.52948...`. The downstream cross-check against `zeta_max^(F1)` to `~1.77e-13` confirms the value is correct — the directive's prose was an approximation (or typo) but the scripts compute the correct value.
2. The directive file `redteam/directives/stage_082.md` does not contain `## Applied: Fn` blocks for F1, F2, or F3 (a minor process slip by Codex), but the user's prompt to the verifier enumerates the applied edits and the diff matches the directive. Not blocking — the substance is correct and the regenerated outputs confirm correctness.
3. The carry-forward provenance comment (sympy 142-146 / mathematica 115-118) was rewritten from a TODO to a citation block referencing the Q2 stage 075 resolution. Cosmetic clarification, no new claim.
4. SymPy uses `import mpmath` mid-file (line 115) rather than at module top. Stylistically unusual but functionally correct and well-justified by the inline comment explaining the `sp.nsolve` instability near `pi/2`. Not a verification concern.
5. The prior verification (May 25, 4 findings: tautological_check, insufficient_verification, mathematica_transliteration, hardcoded_result) is superseded by this report. The earlier "F1 (tautological_check)" was a different finding raised by an earlier auditor cycle and is unrelated to the current F1 (paper_misalignment); to avoid confusion, this file overwrites the prior content entirely.

## Verdict justification

All three findings from the 2026-05-27 audit are resolved. F1's closed-form pin instantiates the missing paper-side `zeta_phys` formula and cross-checks it against the upstream Family-1 constant from stage 084 at the Pe->oo limit, with both engines independently agreeing to ~13 digits — well within the 10^-10 tolerance. F2 is structurally subsumed by F1 (Family-1 numerical pair `(37, 12321/5)` is instantiated) and by the upstream Q2 paper-side edit (Upsilon_w convention unified at `100 Theta_w`). F3's tautological derivative checks are replaced with a non-tautological numerator/denominator factorization of `dzeta_req/dPi_tr` that exercises the strict-positivity content underlying the threshold theorems of notes section 4. Both engines exit 0 with all PASS lines, the regenerated `.txt` outputs reflect the post-fix state, and no regressions are visible in the diff. Verdict: verified.
