---
consult: batch V.2 (first-pass audit, stages 176–187)
date: 2026-05-30
mode: codex-chat read-only ephemeral
codex_session: 019e77af-7315-76e3-8a89-751065348cf2
questions: 3
verdicts: 2 CONCUR, 1 DISPUTE (all resolved; none conceptual; none escalated to user)
prompt: redteam/tmp_prompts/consult_batch_v2.md
---

# Claude+Codex consult — red-team batch V.2

Three narrow math-coverage / how-it's-checked calls from the V.2 first-pass audit. No
paper_misalignments, no conceptual items, nothing escalated to the user (these are the kind of
mirror-policy / feasibility decisions the user delegated to Claude+Codex).

| Q | Topic | Verdict | Action taken |
|---|---|---|---|
| Q1 | 178 F1 — apply independent route vs. accept sanctioned mirror | **CONCUR — apply** | 178 directive applied as-is (no edit) |
| Q2 | 181 F1/F2 — which path closes the Mathematica `FullSimplify→0` | **CONCUR — path (b), recompute from closed forms** | Added `## Consult resolution` steer to `directives/stage_181.md` |
| Q3 | 182 F2 — genuine independent route vs. `## Blocked: F2` | **DISPUTE the Block — lighter genuine route exists** | Added `## Consult resolution` steer to `directives/stage_182.md` |

## Q1 — Stage 178 mirror disposition → CONCUR: apply the independent route

Policy default requires rewriting one-for-one mirrors to native Mathematica primitives before
close (`MATHEMATICA_MIRROR_POLICY.md:75`). 178's `.wl` is mirrored choreography (routes §4
through `pExpected`/`dExpected` then `nuDirect`, wl:73/96). The prescribed
`nuFromData = Coefficient[Series[Log[pA^2/dA^2], eps], …]` route is a **real** native
Series/Coefficient extraction from component data — NOT an X−X port — and `pA`/`dA` both carry
`eps*lam` (wl:62), so the coefficient is non-vacuous. **Do not sanction 178 as a policy mirror;
apply the additive `nuFromData` check.**

## Q2 — Stage 181 Mathematica feasibility → CONCUR: path (b), recompute from closed forms

The fragile point is real: `lamNorm`/`rTarget`/`t2Loaded` are computed under the first
assumption context (wl:32), then `$Assumptions` is reset without `zeta`/`zW`/`omegaW2`/
positivity (wl:55). Since the file already proves `t2Loaded == t2Direct` (wl:48), rebuilding
`t2DirectPert`/`rTargetPert` directly from the closed forms (with local assumptions for `s` and
drift symbols) is the lowest-risk symbolic path — the directive names it as a fallback
(stage_181.md:73,108). First-order log-derivatives of rational products → no reason to expect
non-closure. **Do not raise the 600s cap.**

## Q3 — Stage 182 F2 independence → DISPUTE `## Blocked`: a concrete lighter route exists

The full 8-log→5-Σ inverse is under-determined (directive notes this, stage_182.md:116), but
that is not a blocker. The `.wl` has the direct microscopic forms (wl:62) and five linear Σ
definitions (wl:45); Mathematica can **gauge-fix / Solve five of the eight logs** (linear
`SolveAlways`/`CoefficientArrays`), substitute into `xi1Direct`, `Collect` on the Σ symbols,
and verify each coefficient — a route with a real fail mode (no solution / leftover free logs /
coefficient mismatch) that stays in scope and avoids the forbidden upstream `zeta`
reconstruction (stage_182.md:31). Derive the split via
`sigmaChi -> (sigmaTr - (1 + chi0) sigmaDel)/(1 + deltaU)`. **Reserve `## Blocked: F2` only if
that concrete linear reduction genuinely fails.**

---

## Iteration-2 addendum — stage 185 F1 (checkpoint), codex_session 019e77e6

The verify wave flagged **181 F1** and **185 F1** as `needs_rework`. 181's delta was mechanical (no
consult). 185 F1 needed a second read because the verifier's own proposed fix
(`C_tr_star·(D/(chi0s·deltaUs)) − 1`) is itself X−X (since `C_tr_star` IS that ratio). A focused
read-only consult settled the genuinely load-bearing anchor:

**Route (b) — reconstruct the observable defects independently from the slippage drifts.** Build
`Theta_1`/`Xi_1` from `chi1_indep = chi0s·Sigma_chi`, `deltaU1_indep = deltaUs·Sigma_delta`,
`Sigma_Z`, `Sigma_eps`, `Sigma_delta` (the appendix Stage 182/183 slippage laws) — NOT via
`C_tr_star`/`A_tr_star`. Then `Theta1 − (−C_tr_star·Sigma_tr)` has residual `(C_typed − C_true)·Sigma_tr`,
nonzero (Sigma_tr ≠ 0) for a wrong typed coefficient → genuinely load-bearing. Orchestrator
cross-checked the algebra: `Theta1_indep = −(chi0s·deltaUs/D)·[(1+chi0s)Sigma_delta + (1+deltaUs)Sigma_chi]
= −C_tr_star·Sigma_tr` only when `C_tr_star` is correct (using `Sigma_tr = (1+deltaUs)Sigma_chi +
(1+chi0s)Sigma_delta`). Applied as the 185 iteration-2 delta; the two tautological "coefficient form"
anchors are deleted, F2 (det) and F3 (banner) kept as-is.
