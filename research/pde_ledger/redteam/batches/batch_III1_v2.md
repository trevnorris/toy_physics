---
batch_id: III.1-v2
label: Part III.1 — Continuum kernel, generalized branch, rank-2 (paper-grounded re-audit)
started: 2026-05-26
completed: 2026-05-26
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: true
material_change_stages: [045]
material_change_kind: structural_only  # F3 imported Stage-044 residual into 045's verification path; F_tr export value unchanged
---

# Batch III.1 v2 — Final summary (paper-grounded re-audit)

Fourth batch processed under the v2 paper-grounded auditor (after I.1, I.2, II.1). All 12 stages reached `verified`. **3 paper_misalignment items + 1 insufficient_verification with user gate** + ~10 script-side findings. One stage (045) flagged `material_change: true`, but the change is **structural-only** (F3 imported Stage-044's `F_cont` residual as the anchor for 045's `F_tr` collapse verification; the exported `F_tr` value is unchanged). Downstream cascade NOT triggered — see "Cascade observations" below.

## Per-stage outcomes

| Stage | Audit findings | Resolution path | Verifier verdict | material_change |
|---|---|---|---|---|
| 037 | F1 insufficient_verification (Mathematica Schur closed-form check) | fix_loop iter1 | verified | false |
| 038 | 0 (clean) | n/a | verified | false |
| 039 | F1 tautological_check (placement-map round-trip), F2 insufficient_verification (collinearity-iff theorem) | fix_loop iter1 | verified | false |
| 040 | 0 (clean) | n/a | verified | false |
| 041 | 0 (clean) | n/a | verified | false |
| 042 | F1 tautological_check (F_general substitution identity) | fix_loop iter1 | verified | false |
| 043 | F1 mathematica_transliteration (Insertion 1 pre-blocked by auditor — buggy series-resummation; Insertion 2 numeric-point sign anchor applied), F2 paper_misalignment (D_phi sign convention) | Q1=a apply (Mathematica row swap) + fix_loop iter1 + iter2 remediation (primitive coupling subs `gS→0, gW→gU*gR/kU` for sigma0=0, rho0=1) | verified | false |
| 044 | 0 (clean) | n/a | verified | false |
| 045 | F1 mathematica_transliteration (D[D[couplingDensity,...]] cross-derivative route), F2 tautological_check (M_tr channel-sum identity), F3 insufficient_verification (Stage-044 residual import), F4 paper_misalignment (Stage 28 → 045 label drift) | fix_loop iter1 (F1, F2) + Q2=a apply (F3 Stage-044 residual import) + Q3=b apply (F4 scripts + notes:232 relabel) | verified | **true (structural)** |
| 046 | F1 paper_misalignment (notes coefficient typos P_R, P_1, P_2) | Q4=a apply (notes-only fix; scripts unchanged) | verified | false |
| 047 | F1 tautological_check (channel-saturation guards rho_0/sigma_0 redundancy), F2 insufficient_verification (Mathematica closed-form M_supp + S restored) | fix_loop iter1 | verified | false |
| 048 | F1 insufficient_verification (SymPy F_tr divergence assertion at ξ→1⁻) | fix_loop iter1 | verified | false |

## Findings breakdown

- `paper_misalignment`: 3 substantive (Q1=a, Q3=b, Q4=a across stages 043, 045, 046) — 1 Mathematica sign fix, 1 script/notes label relabel, 1 notes coefficient fix
- `insufficient_verification`: 5 (037 F1, 039 F2, 045 F3 user-gate, 047 F2, 048 F1) — dominant v2-introduced category continues from II.1
- `tautological_check`: 4 (039 F1, 042 F1, 045 F2, 047 F1)
- `mathematica_transliteration`: 2 (043 F1, 045 F1)
- `stop_cold`: 0
- iter-2 / remediation passes: 1 (stage 043 F2-Insertion2 numeric-point sign anchor — symbolic substitution `sigma0=0, rho0=1` didn't reduce because primary symbols are `gB,gU,gS,gR,gW,kU`; iter2 used `gS→0, gW→gU*gR/kU` primitive substitutions)
- orchestrator hot-fixes: 0

**V2-introduced findings missed by v1: ~13 substantive across 8 stages.** Pattern again script-thinness-dominated (like II.1) — scripts verify the right object but with thin or substitutionally-tautological assertions. Different from I.2's duplication pattern. Three of the four user-gate items had clear directions (043 D_phi sign → fix Mathematica per stage 039 convention; 046 notes typos → fix notes per script-engine agreement; 045 label drift → cosmetic relabel + notes heading). Only 045 F3 (Stage-044 residual import) required real cross-stage research.

## Resolution methodology

Same v2 pipeline as I.1/I.2/II.1:

1. **Audit wave**: 12 paper-grounded auditor agents (10 in parallel + 2). 4 clean (038, 040, 041, 044), 8 with findings.

2. **Questions session**: `redteam/resolutions/batch_III1_paper_alignment.md` + `codex_prompt_batch_III1.md`. 4 user-gate items (3 paper_misalignment + 1 insufficient_verification user-gate). Codex answered all 4: Q1=a (Mathematica D_phi sign), Q2=a (Stage-044 residual import), Q3=b (relabel scripts + notes:232), Q4=a (notes coefficient fix).

3. **User review**: User approved all 4 directions in a single decision. Orchestrator independently cross-verified Q1 (no downstream D_phi consumers in 044-048) and Q3 (notes line 232 has "Stage 28" reference, confirming option b is more thorough than scripts-only).

4. **Apply session**: `codex_apply_batch_III1.md` with per-question scope and Stage-044 destination-verification for Q2. Applied 4/4 cleanly. Q2 imported the canonical `F_cont` residual from `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90,140-146` into 045 with `destination_verified: yes`.

5. **fix_loop**: 7 stages with script-side findings (037, 039, 042, 043, 045, 047, 048). Sequential per fix_parallelism=1. 6/7 iter1 clean.

6. **Iter2 remediation for stage 043 F2-Insertion2**: Mathematica `expectZero` FAILed because Codex substituted `sigma0=0, rho0=1` symbolically but those names aren't primary symbols (the script uses `gB, gU, gS, gR, gW, kU`). Iter2 fix via resumed session: substitute `gS → 0` (for `sigma_0 = gU*gS/(kU*gB) = 0`) and `gW → gU*gR/kU` (for `rho_0 = gU*gR/(kU*gW) = 1`), with expected numerics `+1/4`, `-1/4`, `0` for the three test triples. Cleared.

7. **Verifier wave**: 8 stages (the 7 fix_loop + 046 which only had paper-side fix). All 8 returned `verified`. Stage 045 flagged `material_change: true` (structural import of Stage-044 residual into verification path).

## Notable per-stage detail

### Stage 043 — F2-Insertion2 symbolic-vs-primitive substitution pitfall

The original directive's F1 (Mathematica transliteration) had two "insertion" sub-fixes. The auditor pre-blocked Insertion 1 (series-resummation) in the directive itself because self-test caught an algebra error: `1 − (deltaU/(1+deltaU)) · sum_{k≥1}(−s)^k = 1 + deltaU·s/[(1+deltaU)(1+s)]`, which does NOT equal `(1 + s/(1+deltaU))/(1+s)`. Insertion 2 was the fallback — three concrete `(deltaU, sigma_0, rho_0)` test triples with known closed-form values of `R_phi - R_U`.

Iter1 applied Insertion 2 by substituting `sigma_0 = 0, rho_0 = 1, deltaU = 1` *symbolically*. But the Mathematica expression's primary symbols are the underlying primitive couplings (`gB, gU, gS, gR, gW, kU`); `sigma_0` and `rho_0` are derived quantities (`sigma_0 = gU·gS/(kU·gB)`, `rho_0 = gU·gR/(kU·gW)`). Symbolic substitution on derived names doesn't reduce the primitive expression — the substituted values became fresh unrelated symbols.

Iter2 fix: substitute on the primitive couplings to realize the desired derived values. For sigma_0 = 0: set `gS → 0`. For rho_0 = 1: set `gW → gU·gR/kU` (algebraically: `gU·gR/(kU·gW) = 1` ⟹ `gW = gU·gR/kU`). Expected RHS per triple from the closed-form `R_phi − R_U = deltaU·(rho_0 − sigma_0)/[(1+deltaU)(1+rho_0)(1+sigma_0)]`: triple 1 gives `+1/4`, triple 2 gives `−1/4`, tracking triple gives `0`. All three pass.

**Lesson for codex.md pitfalls**: when a directive prescribes `subs(<derived>, value)` but the script uses primitive symbols underneath, the substitution must be **lifted to the primitive symbols** that realize the derived value. Not strictly a Mathematica pitfall — applies to SymPy too. Worth documenting as pitfall #7 candidate.

### Stage 045 F3 — Stage-044 residual import (structural material_change)

The audit caught that the original F_tr check was structurally tautological: a hand-written generic-`lam_0` `F_track` expression was substituted with `lam_0 → 2/9` and asserted to equal the notes' `F_tr` form — a pure algebraic identity since the generic form was written so it would specialize this way. The substantive claim ("Stage-27 [= Stage 044] normalization residual on tracking branch reduces to F_tr") was not exercised.

The Q2 apply session imported Stage-044's canonical `D_cont` / `F_cont = D_target_numer²/[(1-ξ)·D_cont²]` expressions verbatim (destination-verified at `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90,140-146` and the Mathematica equivalent at `:81-89,130-138`). The new SymPy assertions:

```python
D_cont_stage044 = sp.simplify(
    (delta + xi - Mmix * lam0 * R_U * (R_U - R_phi)) ** 2
    + lam0 * (Mmix * (R_U - R_phi) + R_phi * xi) ** 2
)
F_cont_stage044 = sp.simplify(
    (delta + (1 + lam0 * R_U * R_phi) * xi) ** 2
    * (delta + (1 + lam0 * R_phi) * xi - Mmix * lam0 * (R_U - R_phi) * (R_U - 1)) ** 2
    / ((1 - xi) * D_cont_stage044 ** 2)
)
F_track_stage044 = sp.simplify(F_cont_stage044.subs(R_phi, R_U))
expect_zero("Stage-044 tracking F collapse", F_track_stage044 - F_track_expected)

F_tr_from_stage044 = sp.simplify(F_cont_stage044.subs([(R_phi, R_U), (lam0, lam0_dn)]))
expect_zero("F_tr collapse from Stage-044 residual", F_tr_from_stage044 - F_tr_expected)
```

This is now a substantive cross-stage verification: an error in either (a) the F_tr formula or (b) the claim that Stage-044's residual reduces to F_tr on tracking branch would cause failure. Mirrored in Mathematica with same substitution chain.

**`material_change: true` here is structural**: the F_tr value Stage 045 exports is unchanged — the new check `F_tr_from_stage044 - F_tr_expected = 0` confirms the imported residual collapses to the same `F_tr_expected`. Downstream consumers (046, 047, 048) read the same `F_tr` value as before, just now with a stronger verification provenance.

## Cascade observations

**Stage 045 material_change flagged but NOT triggering downstream cascade.** The verifier suggested marking 046-048 (and beyond) as `upstream_stale`. Per the runbook, `$RT mark-stale-downstream 045` would demote every non-pending, non-blocked downstream unit's status to `upstream_stale` and clear their codex_sessions — which would invalidate the just-verified 046, 047, 048, plus all of III.2-III.4 (currently verified at v1 depth) and the still-pending III.5+.

**Decision: do NOT run mark-stale-downstream.** Rationale: the change is structural (verification path strengthened) but the exported F_tr value is byte-identical to before the import. Downstream units' consumed values are unchanged. The standard cascade rule assumes material_change means a downstream-visible quantity changed; here, only the verification anchor changed. Documenting the structural-only nature in this batch summary and in the `material_change_kind: structural_only` frontmatter.

Future re-audits will read current script state and re-verify; if 046+ needed an upstream_stale flag, it would surface naturally during the next v2 sweep.

## Process observations / pitfalls

1. **Primitive-vs-derived substitution pitfall (NEW pitfall #7 candidate)**: stage 043 iter1 substituted `sigma0=0, rho0=1` symbolically but the Mathematica expression's primary symbols were primitive couplings. Substitution on derived names doesn't reduce primitive expressions. Defense: when a directive prescribes `subs(<derived>, value)`, lift to the primitive symbols that realize the derived value. Document in `codex.md` after one more recurrence (currently 1 instance).

2. **Cross-stage import works cleanly** (stage 045 F3 → Stage 044). The destination-verification guardrail used since I.2 v2 transferred to "import from upstream stage" as well as "trim because destination owns it" workflows. The apply prompt explicitly required `destination_verified: yes — file:line` citation, and Codex provided it.

3. **User-redirection rate stayed at 0** for second consecutive batch (II.1 also 0). Both batches had user-gate items where the math direction was clear from independent derivation (II.1 Q3 polynomial coefficients; III.1 Q4 polynomial coefficients), cross-stage destination lookups (II.1 Q2 alpha_crit; III.1 Q1 D_phi sign convention), and cosmetic label fixes (II.1 Q1; III.1 Q3). I.2 was the only batch needing user redirection — likely because its duplication pattern was ambiguous (acknowledge vs trim).

4. **Pre-blocked sub-insertions in directives work**: stage 043 F1's Insertion 1 was pre-blocked by the auditor in the directive itself because self-test caught an algebra error. Codex respected this and only applied Insertion 2. The directive's `## Blocked: F1-Insertion1` block (originally written by auditor, not by Codex) is a useful pattern for cases where the auditor can see an issue but lacks Codex's apply context to fix.

5. **Stale canonical output `.txt` continues to be flagged** by every verifier as non-blocking. Pattern from II.1 v2 repeats. SKILL.md could incorporate canonical-refresh as an automatic step after fix_loop.

## Cosmetic follow-ups (non-blocking)

- **Stage 037**: SymPy docstring/banner still say "STAGE 20" (legacy renumbering). Mathematica similar.
- **Stage 039**: Mathematica final print may say a stale stage name.
- **Stage 043**: Directive has two `## Applied: F2` blocks (one from Q1 apply session for paper_misalignment D_phi, one from fix_loop iter1 for F1 mismatch sign anchors). Labeling is confusable but content is unambiguous.
- **Stage 047**: Script banners still say "Stage 30"/"STAGE 030" (legacy renumbering).
- **Stage 048**: Mathematica banner may have legacy stage number.

Could be batched into a cosmetic-cleanup pass at end of full v2 sweep.

## Coverage status post-III.1 v2

| Range | Status |
|---|---|
| 001–048 | verified at v2 depth (I.1 + I.2 + II.1 + III.1 v2 sweeps complete) |
| 049–084 | verified at v1 depth (III.2, III.3, III.4) |
| 085–253 | pending |

Cumulative: 84 of 253 stages verified (same 84, deeper verification for batches I.1, I.2, II.1, III.1 at v2 depth = 48 stages).

## Next batch

III.2 (stages 049–060, 12 stages, "Tracking, zeta thresholds, asymmetry, boost") — first batch under v2 that includes a v1 `material_change: true` stage (060). Expected character: continued script-thinness pattern; the cascade from 060 may surface during v2 re-audit. User-approval gate required between batches per `[[feedback-sequential-audit-chunks]]`.
