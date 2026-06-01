---
consult: V.3
codex_session: 019e843e-9934-7ed0-88ca-8e2aec9c3181
date: 2026-06-01
mode: read-only
topics: [stage_189 Section II de-tautologization]
verdict: CONCUR (full)
---

# Claude+Codex consult — V.3

Read-only consult (`codex-chat -s read-only -C /var/projects/toy_physics --ephemeral`).
Prompt: `redteam/tmp_prompts/consult_189.md`.

## Q — stage 189 Section II "selected-branch identity"

The iter-1 F1 fix produced `Rtarget_oneport = Lambda0*(1-epseta)/T2_direct` then
`expect_zero(Rtarget_oneport*T2_direct - Lambda0*(1-epseta))`. Because R_target is back-defined,
the residual is identically 0 for ANY `T2_direct`/`Lambda0` (the `.subs(Lambda0, Lambda0_val)`
cancels too) — still tautological. Is the genuine content the direct-slope bridge
`δln T_A^2 = ε λ_A Ξ_1` (notes §2.1-2.3)? Should the selected-branch identity be demoted to a
definition?

## A — CONCUR (full)

1. **DEMOTE** the Section II product assertion to a printed DEFINITION of R_target (rank-2
   bookkeeping; notes lines 215-226). The `Lambda0` substitution does not make it physical.
2. **ADD the genuine direct-slope bridge**: perturb the coherent-branch `T2_coh`, compute
   `Series[log(T2_pert/T2_0)]` (first-order slope), compare to `eps*lambda_A*Xi1_closed` where
   `Xi1_closed` is the INDEPENDENTLY-defined first grouped defect scalar, copied from the upstream
   defect law — NOT defined from that derivative. Reference pattern: stage 181
   (`...stage181...sympy_audit.py:76` eps split + Xi1 def, then `:110/:114` derivative check).
3. (Q1) (b) IS genuine content and is NOT the same as Section I's abstract `dn - Bstar*dr = Xi1`
   (Section I only verifies the observable-packet compiler; the concrete slope test links that
   abstract scalar to the actual continuum transfer shape).
4. (Q2) No additional independent content in the two-`T_A^2`-forms equality (R_target is defined
   by it) unless R_target is imported from a separate continuum-kernel derivation.
5. (Q3) Demoting to a definition is the honest move.
6. **Mathematica mirror**: same — NO `Rtarget = Lambda0(1-epseta)/T2; Simplify[Rtarget*T2-...]`
   assertion; use a printed definition + the non-tautological first-order slope check.

## Disposition

Applied as **189 iteration-2** (both engines). The iter-1 verify verdict (`verified`) was
OVERRIDDEN on orchestrator review — the verify agent flagged the "structurally identity by
construction" smell but passed it; the orchestrator rejected per the no-rubber-stamp rule
(same class as the V.2 checkpoint-185 catch).
