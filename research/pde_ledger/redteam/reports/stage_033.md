---
unit_id: 033
batch: II.1
auditor_model: claude-opus-4-7
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage033_microscopic_normalization_equation.md
  paper_appendix: present
---

# Audit unit 033 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_033.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage033_microscopic_normalization_equation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (only the `\input{stages/stage_033}` line on row 104; no inline summary for this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.txt`

## What the paper claims

Stage 033 writes the selected-branch normalization problem in microscopic coupling variables and packages a stability gate plus an exact necessary onset condition for the universal target. The `\stagefield{Output}` line states verbatim: "Stage~033 outputs the coupling-level definitions \eqref{eq:app-stage033-defs}, the stability gate \eqref{eq:app-stage033-stability}, the microscopic target \eqref{eq:app-stage033-micro-target}, the onset condition \eqref{eq:app-stage033-onset}, and the weak-loading expansion \eqref{eq:app-stage033-weak-loading}." Concretely the deliverables are:

1. **Microscopic definitions** (eq `app-stage033-defs` and `app-stage033-alpha-beta`): `A = K_0 - g_U^2/Omega_U^2`, `Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`, `X = Omega_U^2 g_W + g_R g_U`, with `beta_0 = X^2/Delta_0^2` and `alpha_0 = g_B^2/varpi^2 + X^2/(Omega_U^2 Delta_0)`.
2. **Stability gate** (eq `app-stage033-stability` and `app-stage033-alpha-crit`): `Delta_0 > 0`, `A > 0`, `alpha_0 < alpha_crit` with the closed form `alpha_crit = 9 pi^2 A(A + DeltaK_ax) / [8(11 A + 9 DeltaK_ax)]` after using `kappa_0^2 = 8/pi^2`, `kappa_1^2 = 16/(9 pi^2)`.
3. **Microscopic target** (eq `app-stage033-micro-target`): the selected-branch equation `X^2/Delta_0^2 * s_-(alpha_0)^2 / (kappa_0^2 lambda_-(alpha_0)) = N_Q^{target}` with `alpha_0` as above. This is an equation to be solved at physical parameter values, not an identity.
4. **Onset condition** (eq `app-stage033-onset`): `N_-(0) <= N_Q^{target}`, equivalently `X^2 <= N_Q^{target} A Delta_0^2 / kappa_0^2`. The notes additionally express this as an onset stiffness bound `K_0 >= g_U^2/Omega_U^2 + kappa_0^2 X^2 / (N_Q^{target} Delta_0^2)`.
5. **Weak-loading expansion** (eq `app-stage033-weak-loading` and `-DN`): `N_-(alpha_0) = beta_0 kappa_0^2/A + alpha_0 beta_0 kappa_0^2 (4 A kappa_1^2 + DeltaK_ax kappa_0^2)/(A^2 DeltaK_ax) + O(alpha_0^2)`, with finite-throat constants substituted to `8 beta_0/(pi^2 A) + 64 alpha_0 beta_0 (8A + 9 DeltaK_ax) / (9 pi^4 A^2 DeltaK_ax) + O(alpha_0^2)`.

`\stagefield{Checks}` adds: stability decomposes into three independent conditions; onset follows from monotonicity of `N_-`; weak-loading coefficient is positive for positive `A`, `DeltaK_ax`, `beta_0`. The notes (`§4`) also state the explicit positivity formula `dN_-/dalpha_0 = beta_0 [2 s_- (ds_-/dalpha_0) lambda_- + s_-^3] / (kappa_0^2 lambda_-^2) > 0` for every stable branch point.

## What the script claims to verify

After the prior v1 audit's fixes (directive `stage_033.md` `applied: true`, three findings applied), the scripts now verify:

1. The monotonicity identity `dN/dalpha = beta_0 (2 s_- ds lambda_- + s_-^3)/(kappa_0^2 lambda_-^2)` (Stage 16.1) — this implicitly tests the lemma `dlambda_-/dalpha_0 = -s_-`.
2. The closed-form refined threshold `alpha_crit = 9 pi^2 A(A+DeltaK) / (8(11A + 9 DeltaK))` against the generic form `A(A+DeltaK)/((A+DeltaK) kappa_0^2 + A kappa_1^2)` (Stage 16.2).
3. `N_-(0) = beta_0 kappa_0^2 / A` (Stage 16.3).
4. The weak-loading derivative at `alpha_0 = 0` against both generic and finite-throat-closed forms (Stage 16.4).
5. The microscopic K0-onset value: `sp.solve(N_-(0)|_{mic} == NQ, K0)` matches `gU^2/OmegaU^2 + kappa_0^2 Chi^2/(NQ Delta_0^2)`. Mathematica mirrors this with `Solve[n0Mic == NQ, K0]` plus a back-substitution check `N_-(0) at K0_onset == NQ` (Stage 16.5).
6. The fully substituted stability gate denominator: `together(alpha_crit_mic - alpha_0_mic)` produces a denominator that differs from `8 varpi^2 Omega_U^2 Delta_0 (11 A_mic + 9 DeltaK)` only by a parameter-free constant (which the output transcripts show as `9 pi^2` in sympy and `-9 pi^2` in mma) (Stage 16.6).
7. Numerical cross-check (Mathematica only): the Stage 16.1 monotonicity identity and the Stage 16.6 gate identity are re-evaluated at two rational-parameter rule sets at 30-digit precision.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| 1. Microscopic defs (`A`, `Delta_0`, `X`, `beta_0`, `alpha_0`) | sympy lines 81-85, mma lines 89-93 (direct definitions) | match |
| 2a. `alpha_crit` closed form `= 9 pi^2 A(A+DeltaK)/(8(11A+9DeltaK))` | sympy `expect_zero("alpha_crit - closed finite-throat form")` (line 60), mma line 65 | match |
| 2b. Conditions `Delta_0 > 0`, `A > 0` as stability requirements | `A`, `DeltaK`, `varpi`, `OmegaU`, `OmegaW`, `K0`, `NQ` declared positive in script (sympy line 40, mma line 30-31, 86-87); `Delta_0` is derived, not constrained | partial (script enforces the positivity at the level of symbol domains; paper presents them as gate conditions to be checked, not enforced; this is normal convention for an algebraic-identity audit) |
| 3. Microscopic target `X^2/Delta_0^2 s_-^2/(kappa_0^2 lambda_-) = NQ` (general `alpha_0`) | Components verified (`Nminus = beta_0 s_-^2/(kappa_0^2 lambda_-)` by construction; `beta_0 = X^2/Delta_0^2` by definition at line 84/92); the full equation `N_-(alpha_0) = NQ` is an equation to be solved at physical parameter values, not an identity | match (full identity-style assertion would not be physically meaningful; script-side structural verification plus the alpha_0=0 boundary check (deliverable 4) cover the paper's intent) |
| 4. Onset condition `N_-(0) <= NQ` ↔ `Chi^2 <= NQ A Delta_0^2/kappa_0^2` ↔ K0 stiffness bound | sympy line 87 solves `K0_onset` from `N_-(0)|_mic = NQ`, then line 97-100 checks equality with closed form; mma line 95-97 mirrors via `Solve`, line 107 back-substitutes, line 108-111 checks closed form | match (boundary case is verified algebraically; the inequality direction follows from monotonicity verified in deliverable 5's prerequisite, Stage 16.1) |
| 5a. Weak-loading expansion (generic): `beta_0 kappa_0^2/A + alpha_0 beta_0 kappa_0^2 (4A kappa_1^2 + DeltaK kappa_0^2)/(A^2 DeltaK)` | sympy line 70 `coef1_target`; line 73 `expect_zero(coef1 - coef1_target)` | match |
| 5b. Weak-loading expansion (DN-substituted): `8 beta_0/(pi^2 A) + 64 alpha_0 beta_0 (8A + 9 DeltaK)/(9 pi^4 A^2 DeltaK)` | sympy line 71 `coef1_target_closed`; line 74 `expect_zero(coef1 - coef1_target_closed)`; mma lines 76-78, 82 | match |
| Implied: monotonicity formula `dN/dalpha = beta_0 (2 s_- ds lambda_- + s_-^3)/(kappa_0^2 lambda_-^2)` (notes §4) | sympy line 51-54, mma line 54-60 | match (note: the identity is verified; strict positivity at all stable branch points is not asserted, but the form of the identity makes positivity follow from `s_-, lambda_- > 0`) |
| Implied: positivity of weak-loading coefficient (Checks bullet 3) | Not directly asserted; follows from positive `A`, `DeltaK`, `beta_0` declarations and the closed form's manifest positivity | match (not load-bearing for the symbolic audit) |

`paper_alignment: aligned` — every primary deliverable is exercised; the `partial`/structural items are appropriate to an algebraic-identity audit and not red-team failures.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `expect_zero(dN - dN_formula)` | implied monotonicity formula (notes §4) | yes (implicitly tests `dlambda_-/dalpha_0 = -s_-`) |
| A2 | sympy | 60 | `expect_zero(alpha_crit - alpha_crit_target)` | claim 2a (alpha_crit closed form) | yes (non-trivial after `kappa` substitution) |
| A3 | sympy | 65 | `expect_zero(N0 - beta0*kappa0_sq/A)` | claim 4 prerequisite (`N_-(0)`) | yes |
| A4 | sympy | 73 | `expect_zero(coef1 - coef1_target)` | claim 5a (generic weak-loading) | yes |
| A5 | sympy | 74 | `expect_zero(coef1 - coef1_target_closed)` | claim 5b (DN-substituted weak-loading) | yes |
| A6 | sympy | 97-100 | `expect_zero(K0_onset - target_closed_form)` (K0_onset from `sp.solve`) | claim 4 (onset stiffness bound) | yes |
| A7 | sympy | 120-122 | `assert den_ratio.is_number` | Stage 16.6 gate denominator structure (script-internal scaffolding around claim 2a/1) | yes (substantive guard; the ratio is parameter-free `9 pi^2`) |
| A8 | sympy | 125-128 | `expect_zero(alpha_crit_mic - alpha0_mic - gate_num_target/gate_den_claim)` | same as A7 (decorative) | **no — tautological** (`gate_num_target = gate_num_actual/den_ratio` by construction, so the residual is identically zero) |
| M1 | mma | 60 | `expectZero[dN - dNFormula]` | mirrors A1 | yes |
| M2 | mma | 65 | `expectZero[alphaCrit - alphaCritClosed]` | mirrors A2 | yes |
| M3 | mma | 69 | `expectZero[n0 - beta0*kappa0Sq/A]` | mirrors A3 | yes |
| M4 | mma | 81 | `expectZero[coef1 - coef1Target]` | mirrors A4 | yes |
| M5 | mma | 82 | `expectZero[coef1 - coef1Closed]` | mirrors A5 | yes |
| M6 | mma | 107 | `expectZero[(n0Mic /. K0 -> k0Onset) - NQ]` | claim 4 (back-substitution of solved K0_onset) | yes |
| M7 | mma | 108-111 | `expectZero[k0Onset - (gU^2/OmegaU^2 + kappa0Sq*chi^2/(NQ*delta0^2))]` | claim 4 (closed form vs Solve result; non-tautological after F2 fix) | yes |
| M8 | mma | 122-125 | `If[!NumericQ[denRatio], fail[...]]` | Stage 16.6 gate denominator structure | yes (substantive guard) |
| M9 | mma | 128-131 | `expectZero[alphaCritMic - alpha0Mic - gateNumTarget/gateDenClaim]` | decorative mirror of A8 | **no — tautological** by construction |
| M10 | mma | 147-152 | `If[Abs[monotonicityNumeric] > 10^-20, fail, pass]` over two rational rules | mirrors A1 (independent numeric derivation) | yes |
| M11 | mma | 153-158 | `If[Abs[gateNumeric] > 10^-20, fail, pass]` over two rational rules for gate identity | mirrors M9 (numeric instantiation of the same tautological identity) | **no — tautological** (`gateNumTarget/gateDenClaim ≡ alphaCritMic - alpha0Mic` symbolically; substituting rationals gives 0 by construction) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:124-128`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl:126-131,146-160`

**What's wrong:**
The Stage 16.6 final identity check `expect_zero("alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim", alpha_crit_mic - alpha0_mic - gate_num_target / gate_den_claim)` (sympy line 125-128) is tautological by construction. From lines 105-106 and 124:

```python
gate_diff = sp.cancel(sp.together(alpha_crit_mic - alpha0_mic))      # = num_actual / den_actual
gate_num_actual, gate_den_actual = sp.fraction(gate_diff)
den_ratio = sp.simplify(gate_den_actual / gate_den_claim)
gate_num_target = sp.simplify(gate_num_actual / den_ratio)
```

Substituting back:
- `gate_num_target = gate_num_actual / den_ratio = gate_num_actual * gate_den_claim / gate_den_actual`
- `gate_num_target / gate_den_claim = gate_num_actual / gate_den_actual = gate_diff = alpha_crit_mic - alpha0_mic`

So `alpha_crit_mic - alpha0_mic - gate_num_target / gate_den_claim ≡ 0` in any commutative algebra — it would pass even if `alpha_crit_mic` were perturbed (because `gate_num_target` would be re-derived to compensate). The Mathematica mirror at lines 128-131 reproduces the same construction. The substantive content is in the *guards* at sympy line 120 (`assert den_ratio.is_number`) and mma line 122-124 (`If[!NumericQ[denRatio], fail[...]]`), which DO carry meaning (they fail if the claimed denominator is structurally wrong). The labelled `expect_zero` / `expectZero` is decorative.

The Mathematica numerical cross-check for the gate (lines 153-158) inherits the same tautology: `gateNumTarget/gateDenClaim - alphaCritMic + alpha0Mic` is symbolically zero by construction, so substituting any rational rule and computing `N[..., 30]` produces zero independent of any algebra error in `alphaCritMic` or `alpha0Mic`. (By contrast, the monotonicity numerical check at lines 147-152 IS substantive because `dN - dNFormula` was only confirmed zero via `FullSimplify`; numeric substitution at concrete points is a structurally distinct check that would catch a `FullSimplify` mis-evaluation.)

**Why this matters:**
A reader of the transcripts sees "PASS: alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim" and concludes that the gate-identity factorization has been verified end-to-end. In fact, only the structural part (denominator factor matches up to a parameter-free constant) has been verified; the rest of the assertion is a self-cancellation. If a future maintainer introduces a sign error or wrong constant in `alpha_crit_mic` or `alpha0_mic`, the `den_ratio.is_number` guard catches structural denominator mismatches but the cosmetic equality below it cannot detect a numerator error. The numeric cross-check for the gate identity provides no additional safety net.

The reason this is severity `low` rather than `medium`: (a) the substantive guard `den_ratio.is_number` IS in place and IS non-tautological; (b) the prior v1 audit's F1 directive specifically prescribed this two-stage structure (guard + reconstruct-and-recheck), so removing only the cosmetic step does not undo the v1 fix; (c) downstream stages do not depend on the precise polynomial form of `gate_num_target`. The risk is documentation-level (a maintainer trusts a misleading transcript label) rather than physics-level.

**Required change:**
Two options; pick (a) for the minimum fix, (b) for a stronger fix.

(a) Re-label the decorative assertion to acknowledge its construction. Replace the label `"alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim"` with `"gate_num_target/gate_den_claim - gate_diff (tautological by reconstruction; substantive check is den_ratio above)"`. Same change in mma. Drop the numeric `gate-identity` block from `mathematica/...mathematica_audit.wl:146-160` (keep `monotonicity numeric residual` rules — those are substantive).

(b) Replace the construction with an *independent* polynomial-identity check: instead of defining `gate_num_target` as `gate_num_actual/den_ratio`, write `gate_num_target` as a symbolic polynomial-in-couplings whose coefficients are read out (using `sp.Poly(..., gens=[K0, gU, gB, gR, gW, OmegaU, OmegaW, varpi, DeltaK])`) and check that `expand(alpha_crit_mic * gate_den_claim - alpha0_mic * gate_den_claim) - gate_num_target_independent` is identically zero. This way both the numerator polynomial form and the denominator are pinned independently. Mma analogue: use `PolynomialReduce` / `CoefficientList` to build the same polynomial-coefficient identity.

For Codex auto-application, (a) is the safe minimum (cosmetic relabel + drop the redundant numeric block); (b) requires more design and would be deferred to a follow-up directive if the user wants stronger coverage.

**Verification:**
After fix (a): the new sympy and mma assertion labels reflect that `gate_num_target/gate_den_claim` is tautologically equal to the difference by construction; the gate-identity numeric block in the mma script is absent; both scripts still exit 0 and the substantive guards still appear in the transcript.

## Independent-derivation check (Mathematica)

The Mathematica script's analytic block remains a transliteration of the SymPy script (same `r`, `lambdaMinus`, `sMinus`, `nMinus`, same microscopic substitutions, same gate construction). The v1 F3 directive added a numerical cross-check block (mma lines 138-160) that addresses the transliteration concern *for the monotonicity identity* (`monotonicityNumeric` evaluates `dN - dNFormula` at two rational rules; this is structurally distinct from the analytic `FullSimplify` because numeric substitution would catch an erroneous `FullSimplify` reduction). The cross-check for the gate identity at lines 153-158, however, is still tautological by construction of `gateNumTarget` (see F1), so it does not extend the independent-derivation coverage to Stage 16.6. A new `mathematica_transliteration` finding is NOT raised because (i) the monotonicity numeric check IS independent, satisfying the second-engine policy for the most physics-loaded identity, and (ii) the Stage 16.6 weakness is captured under F1 as a tautology rather than a transliteration.

## Engine cross-check

Both engines report PASS on every assertion. Numerical outputs printed agree up to formatting:
- `N_-(0) = 8 beta_0 / (pi^2 A)` — sympy output line 38, mma output line 13.
- `alpha_crit = 9 pi^2 A (A + DeltaK) / (8 (11 A + 9 DeltaK))` — sympy output line 32, mma output line 10 prints `(9 A (A+DeltaK) Pi^2)/(88 A + 72 DeltaK)` (same value: `88A + 72 DeltaK = 8(11A+9DeltaK)`).
- `alpha0(mic) = g_B^2/varpi^2 + Chi^2/(Omega_U^2 (Omega_U^2 Omega_W^2 - 88 g_R^2/(9 pi^2)))` — sympy line 55, mma line 25.
- `K0_onset = g_U^2/Omega_U^2 + 648 pi^2 (...)^2 / (NQ (9 pi^2 Omega_U^2 Omega_W^2 - 88 g_R^2)^2)` — sympy line 57, mma line 27 (different surface form: mma writes `(7744 gR^4 gU^2 NQ + 72 OmegaU^2 (...) Pi^2 + 81 gU^2 NQ OmegaU^4 OmegaW^4 Pi^4) / (NQ (88 gR^2 OmegaU - 9 OmegaU^3 OmegaW^2 Pi^2)^2)`; both reduce to the same closed form, confirmed by mma's M7 assertion which passes).
- `denominator ratio` differs in sign (`+9 pi^2` sympy vs `-9 pi^2` mma) — both numeric/parameter-free, both pass the guard.
- Mma's `N::meprec` warnings during numerical evaluation indicate precision-limit warnings while computing the (large symbolic) residuals to 30 digits; the final reported magnitudes (`0``78.83...` etc.) are at 78-digit residuals (much smaller than the `10^-20` tolerance). No engine disagreement.

Outputs are fresh: sympy script mtime `1779406179`, sympy output mtime `1779406300` (output newer by 121 s); mma script mtime `1779406179`, mma output mtime `1779406310` (output newer by 131 s). `outputs_fresh: true`.

## Verdict justification

The paper's five `\stagefield{Output}` deliverables (microscopic definitions, stability gate including `alpha_crit` closed form, microscopic target structure, K0 onset, weak-loading expansion in both generic and DN-substituted forms) are all covered by non-tautological script checks after the v1 audit's three fixes (`applied: true`). The micro-target equation is appropriately covered piecewise (definitions + boundary + structural identity rather than as a single global identity, which would be physically meaningless). The monotonicity identity from notes §4 is verified analytically in both engines and numerically in Mathematica with a structurally distinct path (the only true independent-derivation step in the Mathematica script). The K0_onset block is now solve-based in both engines, eliminating the v1 F2 tautology, and the back-substitution check (mma M6) provides a non-trivial consistency test.

The one residual issue is the Stage 16.6 final `expect_zero`/`expectZero` (and its numerical mma echo): these are tautological by construction of `gate_num_target = gate_num_actual / den_ratio`. The v1 F1 fix introduced a substantive guard (`den_ratio.is_number` / `NumericQ[denRatio]`) that does carry the verification load, but left the cosmetic identity in place where a reader might mistake it for an independent check. This is low-severity (documentation-level, not physics-level), and I file it as F1 here primarily so that a future maintainer is not misled by the transcript labels.

Verdict: `findings` (1 low-severity tautology). No stop-cold. `paper_alignment: aligned`.

## Self-test notes

I checked the four required traps. (1) Variable independence: every `sp.diff(EXPR, alpha0)` / `D[expr, alpha0]` has `alpha0` genuinely present in `EXPR` (`s_minus` and `Nminus` both depend on `alpha0` through `R` and `lambda_minus`); no zero-by-construction derivatives. (2) Symmetry/parity: no integrals in this unit; trap N/A. (3) Trivial-case pre-check: at `alpha_0 = 0` I hand-evaluated `R = DeltaK`, `lambda_- = A`, `s_- = (sigma + delta_kappa)/2 = kappa_0^2`, so `N_-(0) = beta_0 kappa_0^2 / A` matches the script; I also hand-checked the F1 tautology by substituting `gate_num_target/gate_den_claim` and confirmed the residual is structurally zero. (4) Paper round-trip: the F1 fix proposal (option (a) cosmetic relabel + drop redundant numeric block) introduces no new paper-side claim and does not change any verified identity; the suggested option (b) would re-pin the numerator polynomial form, which is also script-internal scaffolding without paper-side counterpart. No new `paper_misalignment` risk.
