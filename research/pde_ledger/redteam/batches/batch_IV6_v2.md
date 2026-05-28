---
batch: IV.6
range: 151-163
total_stages: 13
verified: 13
findings_count: 19
findings_resolved: 19
findings_blocked_legitimate: 0
material_change_count: 1
clean_stages: [153, 159, 162]
status_only: [153]
dirty_stages: [151, 152, 154, 155, 156, 157, 158, 160, 161, 163]
checkpoints: []
audit_date: 2026-05-27
verify_date: 2026-05-28
status: closed
---

# Red-team batch IV.6 — Correction, coevolution, traction, off-family

## Summary

13-stage audit unit for IV.6 (`Part IV.6 — Correction, coevolution, traction, off-family`). 12 dual-engine stages and 1 status-only stage (153, `full_mouth_correction_status`). No checkpoints.

19 findings distributed across 10 dirty stages — 19 resolved, 0 blocked. 13 stages verified end-to-end; both engines exit 0; outputs fresh.

One material change: **stage 151 SymPy was rewritten from a symbolic-integration approach to mpmath numerical integration** with concrete `Pi_star`, `r1`, `r2` anchors because `sp.integrate(Sigma_star * cos(pi*x/2), (x, 0, 1))` with symbolic `Pi_star` was pathologically slow (>30 minutes CPU before manual kill). The verifier classified the rewrite as `material_change: true`; downstream impact is zero since the script verifies algebraic identities, not specific numeric outputs consumed downstream.

Twelfth consecutive zero-redirection batch.

## Per-stage findings tally

| Stage | Status | Findings | Engines | Notes |
|-------|--------|----------|---------|-------|
| 151 | dirty | 3 | SymPy + Mathematica | F1 tautological (delta_g definition vs integration), F2 insufficient (deltaPi/deltaT not asserted), F3 transliteration (Series-based Mathematica path) |
| 152 | dirty | 1 | SymPy + Mathematica | F1 insufficient (4 scale anchors absent in SymPy) |
| 153 | clean | 0 | — | Status-only; legitimate carry-forward |
| 154 | dirty | 1 | SymPy + Mathematica | F1 transliteration + banner; Series-based Mathematica path |
| 155 | dirty | 2 | SymPy + Mathematica | F1 insufficient (4 fixed-point asserts), F2 paper_misalignment/banner |
| 156 | dirty | 1 | SymPy + Mathematica | F1 insufficient (4 canonical-point numerics absent in SymPy) |
| 157 | dirty | 3 | SymPy + Mathematica | F1 transliteration (Solve-based independent delta C), F2 insufficient (trivial mul-by-zero check), F3 stale_output informational |
| 158 | dirty | 3 | SymPy + Mathematica | F1 paper_misalignment (Cluster C: forward-carry items 2-3), F2 tautological (Ms law), F3 insufficient (composed delta Mq + delta Pi laws) |
| 159 | clean | 0 | — | All assertions match paper; H1 + body citations renumbered via Cluster A/B |
| 160 | dirty | 1 | SymPy + Mathematica | F1 transliteration (direct chain-rule total-differential replaces Series + Coefficient) |
| 161 | dirty | 3 | SymPy + Mathematica | F1 tautological (eps_gamma derivation), F2 tautological (linearized slippage), F3 transliteration + banner |
| 162 | clean | 0 | — | All three load-bearing identities pass; engines agree to ~28 sig figs |
| 163 | dirty | 1 | SymPy + Mathematica | F1 transliteration (implicit-function + Series chain-rule independent paths) |

**Totals:** 19 findings, 10 dirty stages, 3 clean stages.

## Three user-gate clusters (all `(Recommended)`)

### Cluster A — Mass renumbering (mechanical, 20 edits)

- 18 banner edits at −17 offset across 9 stages × 2 engines (154 137→154, 155 138→155, 156 139→156, 158 141→158, 159 142→159, 160 143→160, 161 144→161, 162 145→162, 163 146→163).
- 2 notes H1 edits at −85 offset: 159 (Stage 244→159), 160 (Stage 245→160).

Stages 151, 152, 157 had no banner offset (151/152 use title-only convention; 157 banner was correct).

### Cluster B — Body-text forward-stage citation re-attribution (53 edits across 13 notes files)

Three offset rules surfaced (content-verified):
- **−51** for pre-renumber 188-199 → current 137-148 (IV.4/IV.5 references)
- **−85** for pre-renumber 239-248 → current 154-163 (IV.6 internal cross-references)
- **−102** for pre-renumber 221 → current 119 (IV.3 parent compensation family) and 249-250 → current 147-148 (IV.5 rigidity kernel + representative families)

Notes are now self-consistent. The −85 offset is new in IV.6 — IV.5's body citations used −102 for the 220-251 range; IV.6's pre-renumber sequence numbered 239-248 as a contiguous block that became current 154-163.

### Cluster C — Stage 158 paper-card Checks downgrade

Stage 158's paper card listed three `\stagefield{Checks}` items but only item 1 was script-side verified. Items 2 (even-preservation) and 3 (tangent motion δ⊥=0) describe the broader transport program — even-preservation is verified downstream in stage 159 (canonical-even gate in `hybrid_outlet_projection`), tangent motion δ⊥=0 in stages 162/163 (parent compensation rigidity + off-family normal coordinate).

Items 2-3 rewritten as forward-carry citations of `\ref{stage:159}` and `\ref{stage:162}` / `\ref{stage:163}`. **First forward (downstream) carry-forward in the v2 batches** — IV.4 stage 134 and IV.5 stage 144 both carried items upstream.

## Orchestrator catches (rework loop)

1. **Stage 151 SymPy symbolic-integration hang.** Directive prescribed `sp.integrate(Sigma_star * c_kernel, (x, 0, 1))` with symbolic `Pi_star`. SymPy's symbolic integrator hung (>30 min CPU, killed). Rewrote with mpmath numerical integration at `Pi_star = 1.50882951349316`, `r1 = 1.7`, `r2 = -0.9`, verifying all moment-shift identities to ~40 dps. **Material change** (verifier flagged), but downstream impact zero — script verifies algebraic identities not numeric outputs. **First v2 batch where the auditor-prescribed approach was infeasible at the engine level and the orchestrator had to redesign the verification.**

2. **Stage 154 Mathematica multivariate `Series` keeps cross-products.** Directive prescribed `Normal[Series[piExpr, {dSigma0, 0, 1}, {dR, 0, 1}, {dS, 0, 1}]]` — but this multivariate Series gives the bivariate-truncation hull, not the joint linear part, so it retained `dSigma0*dR`, `dSigma0*dS`, `dR*dS` cross-terms and the `dPi identity` assertion failed with residual `-(dS*dSigma0*rStar) - dR*(dS*(dSigma0+sigma0)+dSigma0*sStar)`. Switched to a single-epsilon parameterization `piExprEps = piExpr /. {dSigma0 -> epsLin*dSigma0, dR -> epsLin*dR, dS -> epsLin*dS}` and `Normal[Series[piExprEps, {epsLin, 0, 1}]] /. epsLin -> 1`, which correctly extracts only the joint linear part.

3. **Stage 161 directive variable-substitution typo.** Directive prescribed `expect_zero("d eps_gamma = d ln gamma0 - d ln(1+r_c)", depsg_direct.subs(dgamma0, (1+rc)*dln_gamma0/9) - (dln_gamma0 - drc/(1+rc)))`, but `depsg_direct` retained the free `gamma0_sym`. The correct expression substitutes `gamma0_sym -> (1+rc)/9` first (yielding `depsg_branch`), then substitutes `dgamma0`. Caught when SymPy first run failed with residual `drc*(rc+1 - 9*gamma0)/(rc+1)^2`. Mirrored fix in Mathematica.

## Verification

All 13 verification files written under `redteam/verifications/stage_*.md`. Final verdicts:

- `verified` (10 from finding-bearing stages): 151, 152, 154, 155, 156, 157, 158, 160, 161, 163
- `verified` (3 from clean stages): 153, 159, 162
- `needs_rework` (0)
- `blocked_unfixable` (0)

Material change: 1 (stage 151).

## Cumulative

Range 001-163 paper-aligned at v2 depth. **175/253 stages red-team verified (69.2%).** Twelfth consecutive zero-redirection batch.

Next batch (sequential-audit-chunks rule, awaits explicit user authorization): **V.1 = stages 164-175** ("Microscopic logs, drifts, bundle inversion"), 12 stages.
