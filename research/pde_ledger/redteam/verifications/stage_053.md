---
unit_id: 053
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 053

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py:74-76` — the literal `expected_linear = sp.simplify(2 / pi - sp.Rational(1, 2))` was replaced with `linear_coeff = sp.simplify(series_small.coeff(alpha, 1))`, and the assertion now reads `expect_zero("linear coefficient - (4-pi)/(2pi)", linear_coeff - (4 - pi) / (2 * pi))`.
- `mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl:70` — the literal `linearCoeff = FullSimplify[2/Pi - 1/2];` was replaced with `linearCoeff = FullSimplify[Coefficient[seriesSmall, alpha, 1], Assumptions -> ell > 0];`. Lines 69 (`seriesSmall = ...`) and 78 (`expectZero["linear coefficient - (4-Pi)/(2Pi)", linearCoeff - (4 - Pi)/(2 Pi)]`) are unchanged.

**Assessment:**
Both edits match the directive's "required change" verbatim. The new SymPy `linear_coeff` is extracted from `series_small`, which itself depends on the computed `Omega_alpha = I_alpha / I_W` (i.e., on the actual integration). The Mathematica `linearCoeff` is now `Coefficient[seriesSmall, alpha, 1]`, where `seriesSmall = Series[omegaAlpha, {alpha, 0, 2}]` is built from the integrated `omegaAlpha`. The post-fix assertions are therefore non-tautological — they would fail if the integrand or the closed-form check on line 63/65 produced an incorrect Omega_alpha. The SymPy transcript prints `small-alpha series = alpha**2*(-4/pi**2 + 1/12 + 1/pi) + alpha*(-1/2 + 2/pi) + 1` and `linear coefficient = (4 - pi)/(2*pi)`; the Mathematica transcript prints `small-alpha series = 1 + alpha^2*(1/12 + (-4 + Pi)/Pi^2) + alpha*(-1/2 + 2/Pi)` and `linear coefficient = -1/2 + 2/Pi`. Both assertions evaluate to 0 / PASS. No collateral edits beyond the orchestrator's preemptive `ConditionalExpression` strip in the `expectZero` helper (noted by the orchestrator as a generic idiom; no substantive effect here since none of the residuals in this script are `ConditionalExpression`-wrapped).

## Exec log assessment

**SymPy:** exit=0 (inferred from refreshed transcript; orchestrator log line 61 confirms `refresh sympy` succeeded and stage advanced to `codex_applied`). Notable lines:
- `small-alpha series = alpha**2*(-4/pi**2 + 1/12 + 1/pi) + alpha*(-1/2 + 2/pi) + 1`
- `linear coefficient = (4 - pi)/(2*pi)`
- `linear coefficient - (4-pi)/(2pi) = 0`

**Mathematica:** exit=0 (transcript ends with `Stage 053 Mathematica audit passed.`; orchestrator log line 62 confirms `refresh mathematica` succeeded). Notable lines:
- `small-alpha series = 1 + alpha^2*(1/12 + (-4 + Pi)/Pi^2) + alpha*(-1/2 + 2/Pi)`
- `linear coefficient = -1/2 + 2/Pi`
- `PASS: linear coefficient - (4-Pi)/(2Pi)`

(The two `Limit::alimv` warnings on the `omega0` / `omegaInf` lines are unchanged from the pre-fix audit; they reflect the existing pattern of passing `Assumptions -> ell > 0` while the limit variable is `alpha` — they do not block any assertion and predate this directive.)

**Output freshness:** Confirmed. Script mtimes 2026-05-22 17:36:14; output mtimes 17:37:18 (SymPy) and 17:37:25 (Mathematica). Both outputs were regenerated after the post-Codex edits.

## Material-change assessment

`material_change`: false.

The new linear-coefficient extraction produces the same value `(4-pi)/(2*pi) = 2/pi - 1/2` that the prior literal carried, and the substantive downstream numbers (`Omega_max = pi/2`, `A_I,max = pi^2/4`, `Omega_alpha(alpha->0) = 1`, `Omega_alpha(alpha->+infty) = pi/2`, `Sigma_total = L`) are unchanged. No quoted constants or derived results are altered, so no downstream unit needs re-audit on account of this fix.

## Side observations (non-blocking)

- The orchestrator-applied `ConditionalExpression -> e` strip in the Mathematica `expectZero` helper (`.wl:22-27`) is a generic idiom fix and changes no assertion in this script. Flagged only because it appears in the diff alongside F1; it is not a finding.
- The pre-existing `Limit::alimv` warnings on the alpha-limit lines are independent of F1 and out of scope here.

## Verdict justification

The directive's single finding F1 was applied exactly as specified in both engines: the hardcoded `2/pi - 1/2` literals are gone, replaced by genuine extractions from the computed small-alpha series. Both refreshed transcripts show the new assertion passing and confirm the linear coefficient is read from the integrated `Omega_alpha`, so the check is now non-tautological. No regressions, no scope creep, and no downstream-affecting numerical changes.
