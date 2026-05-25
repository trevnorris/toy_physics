---
batch_id: III.4
label: Part III.4 — Family-1 geometry, thresholds, quadrupole
started: 2026-05-22
completed: 2026-05-25
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: false
---

# Batch III.4 — Final summary

All 12 stages reached `verified` status. Both engines (SymPy + Mathematica) exit 0 on every unit; all 40 findings closed substantively. **First batch where every stage carried at least one finding** — no clean-on-first-read stages in 073-084. No `material_change: true` flag; all derivation-route rewrites left printed symbolic content and numeric values byte-identical.

## Per-stage outcomes

| Stage | Audit findings | Codex iterations | Verifier verdict | Notes |
|---|---|---|---|---|
| 073 | 4 (taut×2, hardcoded, math_translit info) | 1 | verified | Mathematica precedence bug in `eta /. (len/ell) -> 37 - 37` fixed via parentheses (`Rule` had lower precedence than `Plus`, making the residual trivially 0); both `eta` and `Lambda_ell` reconstructions now built from symbolic `K_m` / `L/ell` and a symbolic `Lambda_ell - L/ell` identity check; F4 transliteration logged informationally per directive |
| 074 | 2 (taut, insuff_verif) | 1 | verified | `chi_lock = Lambda_ell/2` tautology replaced with the physical substitution chain `chi_def = m_psi*c_s*L/hbar` → `subs(c_s)` → `subs(L/ell)`; four-line provenance comment anchoring `4` and `4/5` Family-1 EL coefficients inserted; `alpha_numeric = 49.640709...` byte-identical pre/post |
| 075 | 3 (taut, hardcoded, math_translit) | 1 | verified | F1 round-trip now uses constructed `Upsilon_fail`/`Theta_fail` instead of `100*Theta_w == 100*Theta_w`; F2 added free-symbol `Delta_0`/`Delta_inf` identity checks (alpha_sym/eta_sym must be proved by simplify, not be self-consistent literals); F3 discharged via F2's independent-derivation leg |
| 076 | 3 (taut, hardcoded, math_translit) | 1 | verified | `U` now derived from `Integrate[P/rho^2]`; `mu_star` and `c_sw` routed through `sp.solve`/`Solve`; n=5 enthalpy identity paired with an explicit n=3 fail-check (`raise AssertionError if residual_n3 == 0`); F2 docstring typo (`(25/4)` → `25`) human-resolved as a stale comment — math, assertion, and paper all consistent on `25`; `thetaHealTarget` → `thetaHealReduced` rename breaks line-by-line `.wl`/`.py` mapping. One controlled codex deviation (`P = K*rho^n_poly` instead of literal `P = K*rho^(1+1/n)`) — verifier confirmed necessary, the directive's literal form was internally inconsistent with the `h = m cs^2/4` identity |
| 077 | 3 (insuff_verif×2, taut) | 1 | verified | Symbolic `1 - alpha_r * S(xi_*)^2 = 0` identity added in both engines (proves `xi_* = atanh(2/sqrt(alpha_r) - 1)` is the cut point rather than postulating it); tautological `xi_* numeric check` removed from `.wl`; SymPy `expect_close` helper added for the per-value Jensen-floor tolerance checks; four cross-engine constants now verified to ~1e-50 |
| 078 | 4 (taut, hardcoded, math_translit, insuff_verif) | 1 | verified | Provenance comments on `Theta_chi_coeff`/`Theta_J_coeff` literals; Mathematica replaced literal-decimal coefficients with symbolic `Sinh`/`Cosh` closed forms + high-precision `ToExpression["...`40"]` loads; three branch-verdict assertions added (`Pe_suff_J < Pe_suff_chi`, `Pe_fail_J < Pe_fail_chi`, `Pe_suff_chi < Pe_fail_J`). One controlled codex deviation (removed the directive's spurious `100` factor on `thetaSuffSym`) — verifier confirmed this was a math error in the directive itself |
| 079 | 3 (hardcoded, taut, math_translit) | 1 | verified | Bypassed literal `zeta0`/`zetaInf` decimals (copied byte-for-byte from SymPy output); now uses `Limit[aF1*omega^2, pe → 0/Infinity]` with `expectApprox(target=0, tol=10^-40)` per the directive's permitted fallback (high-precision numeric atoms can't be `=== 0`-tested); F3 added a `D[omega, pe]` symbolic-derivative slope check returning `(4-Pi)/(2*Pi)` |
| 080 | 3 (taut, hardcoded, insuff_verif) | 1 | verified | SymPy: four `zeta numeric check` lines now exercise the Stage-61 Pe constants instead of being trivially-passing `c*lambda^2 → oo`; Mathematica: four `zetaTarget*` independent-path computations replace literal SymPy-output copies; four ordering PASS lines (chi-pair, J-pair, J≤chi suff, J≤chi fail) added per the prose ledger |
| 081 | 4 (hardcoded, taut×3) | 1 | verified | `qq = piOfZeta / cMix` now derived via `Solve[zeta == zetaExpr, piTr]`; five `expectApprox` magic-number self-checks replaced with residuals against `1 + zeta_*` functional form; blocking-ceiling check rewritten as reciprocal identity `epsCeiling*zetaMaxF1 - 1`. **One orchestrator-applied hot-fix**: codex's `Solve` introduced `ConditionalExpression` wrappers that the `expectZero` helper failed on; standard 3-line `ConditionalExpression[e_, _] :> e` strip retrofitted to the helper and to `piOfZeta`/`qq` assignments (downstream `zeta → 0/1` substitutions need the bare form). Substantive PASS unchanged (residuals were 0 modulo wrapper before; 0 plain after) |
| 082 | 4 (taut, insuff_verif, math_translit, hardcoded) | 1 | verified | Two `Xi_F1` arithmetic self-checks demoted to `print` (no `PASS:` follow); two new derivative assertions (`dR_quad/dzeta_phys + 1 = 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)`) added in both engines; `.wl` `zetaReq` now `Solve`-derived from `qMap` with `ConditionalExpression` strip per project idiom; TODO(provenance) on `Lambda_ell = 37` / `lambdaEll = 37` literals |
| 083 | 4 (math_translit, hardcoded, insuff_verif, taut) | 1 | verified | `delta0Residual`, `deltaInfResidual`, `omegaResidual`, `y_F1`/`aF1Indep` cross-engine identity checks and `dzetaDpe` monotonicity sign-check added; SOURCE-ANCHOR comment blocks for `Theta_chi_coeff`/`Theta_J_coeff`/`136900` literals; legacy 11 numeric self-checks kept per directive option (b) since F1's identity residuals now carry the substantive verification. One benign deviation: `nsolve(prec=80)` bump for numeric stability |
| 084 | 3 (taut, hardcoded, insuff_verif) | 1 | verified | Two `(37^2 - 1369)`-style tautologies replaced with cross-route `(xiF1FromUpsilon /. upsilonW → 100*thetaW) - xiF1FromTheta` consistency check; two hardcoded float-difference `expectApprox` calls replaced with four `expectZero[..., If[TrueQ[...], 0, 1]]` ordering inequalities (`zetaMinusJ`/`zetaPlusJ` now in the assertion set); added `FindRoot`/`Limit[zetaPhys, Pe→Infinity]` block returning `2.467529229456...` matching `zetaMaxF1` to ~14 digits |

## Findings breakdown

- `tautological_check`: 14
- `hardcoded_result`: 12
- `mathematica_transliteration`: 7
- `insufficient_verification`: 7

Total: 40 across 12 stages (avg 3.33/stage, highest yet — III.3 was 2.25/stage). Three categories tied or close: `hardcoded_result` rose sharply (1 in III.3, 12 in III.4) because the Family-1 numerology stages (075-084) accumulate many literal constants (`alpha_r=10`, `1/20`, `4.06863235...`, `0.927552032...`, `136900`, `1369`, etc.) and the auditors held each to either a derivation or a provenance comment.

## Toolchain hardening landed mid-batch

- **`ConditionalExpression` strip retrofit on stage 081** (orchestrator-applied). Codex's directive-mandated `Solve[zeta == zetaExpr, piTr]` returned `ConditionalExpression[..., conditions]` under aggressive `$Assumptions`. The `expectZero` helper's `=== 0` test failed on the wrapped residual, and downstream `qq /. zeta -> {0, 1, ...}` substitutions hit the ConditionalExpression's boundary and returned `Undefined`. Standard 3-line strip pattern (`res = res /. ConditionalExpression[e_, _] :> e`) added to the helper, plus matching strip on `piOfZeta` and `qq` assignments. This is the same pattern documented in `codex.md` after III.2 — the directive for 081 did not preemptively include it, and codex did not preemptively add it. **Action item**: consider updating the directive template so that any stage where the directive prescribes a `Solve`/`Reduce` call adjacent to downstream substitutions also prescribes the strip. III.4 found one such gap; previous batches' directives happened to prescribe it explicitly.

## Process observations

1. **Zero `material_change: true` this batch** (vs. one each in III.2/III.3). All derivation-route rewrites — including the substantial Solve/Integrate/Limit/Series additions on 075-084 — left printed symbolic forms and numeric outputs byte-identical pre/post. The transition from "every stage is dirty" (40 findings across 12 stages) to "no material change" is the expected pattern when the audits are surfacing assertion-strength gaps rather than algebra errors.

2. **First batch with no clean stages.** III.1 had 2 clean (042, 048); III.2 had 1 (056); III.3 had 2 (061, 066); III.4 has 0. The Family-1 numerology cluster (075-084) packs many literal constants into compact stages, which the auditors flagged uniformly as `hardcoded_result`. Worth tracking whether this pattern persists into III.5+.

3. **One orchestrator-applied mid-batch hot-fix** (stage 081 `ConditionalExpression` strip). Caught by the fix_loop's post-codex sanity-exec FAIL grep; resolved by direct edit on the `.wl`, manual `math -script` re-run to confirm exit 0, then `$RT set-status 081 codex_applied` + `$RT capture-diff 081` + resume of fix_loop on 082-084. The orchestrator's intervention was a known canonical pattern from project memory (`feedback_mathematica_script_idioms`), not novel debugging.

4. **Two controlled codex deviations**, both verified necessary:
   - Stage 076 F1: codex used `P = K*rho^n_poly` instead of the directive's literal `P = K*rho^(1+1/n_poly)`. Verifier confirmed the directive's literal form is internally inconsistent with the `h = m cs^2/4` identity; codex's form makes the identity hold iff `n=5`. Acceptable.
   - Stage 078 F3: codex removed the directive's spurious `100` factor on `thetaSuffSym`. Verifier confirmed this was a math error in the directive itself — the operational SymPy `Theta_suff_coeff` is `4.21495e-2` (Upsilon/100), so the `100` factor would have broken the new `expectApprox`. Codex caught the error; acceptable.

5. **Exec_logs absent across all 12 verifications** (same as III.3). Verifiers fell back to the freshly-regenerated canonical output `.txt` transcripts. Non-blocking. Same action-item as III.3: consider making `fix_loop.sh` also copy the sanity-exec output into `redteam/exec_logs/` for the verifier's audit trail.

6. **One human-resolved docstring discrepancy** (stage 076 F2). Auditor flagged that the docstring claimed `Theta_w = (25/4) lambda_mu^2 rho_w^2` while the assertion and paper used `25 lambda_mu^2 rho_w^2`. Human review (cross-referenced `paper/appendices/stage_ledger.tex:221`, `paper/appendices/stage_appendix_part03.tex:130`, `paper/stages/stage_077.tex:24`, `paper/stages/stage_082.tex:46`) confirmed the docstring was a stale typo; the math, assertion, and paper are correct on `25`. Directive simplified to "delete the `/4` from the docstring" before fix_loop launch. No paper change, no downstream cascade.

## Prompt hardening landed this batch

None. The III.2 `codex.md` additions (`ConditionalExpression` strip, `1/pi1 == 0` infinity test) held in 11 of 12 stages. Stage 081 needed the strip retrofit, but that was a per-stage application of the documented pattern, not a new pitfall.

## Tracker updates landed in this commit

- `notes/MATHEMATICA_MIRROR_POLICY.md`: 12 new Independent-Mirror Set entries (073-084); snapshot bumped to "2026-05-25 (batch III.4 close)".
- `notes/CHECKPOINT_TRUST_AUDIT.md`: snapshot date updated; tier unchanged.
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: snapshot date bumped; no new constants surfaced.
- `notes/PAPER_CLEANUP_TRACKER.md`: new P4-37 (III.4 batch) row and change-log entry dated 2026-05-25.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: III.4 is out of the linear projected-EM core range (004-021); only a date bump + completed-checks bullet noting zero material_change.
- `notes/STAGE_VERIFICATION_COVERAGE.md`: snapshot bumped; cumulative count updated to 84 / 253; new III.4 row in the per-batch coverage table.
