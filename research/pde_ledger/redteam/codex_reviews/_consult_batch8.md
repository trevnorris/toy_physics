---
consult: batch_8
date: 2026-05-29
mode: read-only ephemeral
codex_session_id: 019e7760-3aba-7340-8d81-f3ace2add17d
stages: [175]
outcome: 4/4 CONCUR — no dispute, no conceptual escalation
---

# Batch-8 Claude+Codex consult — stage 175 (R1 only)

One read-only consult covering the single open finding on stage 175 (the LAST of the 29
FINDINGS stages). R1 = the Mathematica `Sigma_N` differential block is a line-by-line
transliteration of the SymPy block (same expression construction, both using the same
`dlog = d/deps log(.)|_{eps->0}` primitive); the V.1 disposition accepted it as a policy
mirror (F3-step3), and the re-review (`codex_reviews/stage_175.md`) rejects that waiver.

## SUMMARY

All four design questions resolved by unconditional CONCUR. Decision: **option (B) —
SUPPLEMENT** the `.wl` Sigma_N block with ONE independent `dlogSeries` (Series+Coefficient)
series-route check, comparing series-DIRECT vs SHAPE target, leaving the existing `dlog`
lines untouched as corroboration. Mandatory escape clause: if `dlogSeries` does not land
`=== 0` robustly at exec time (CAS-normalization limit, not a real algebraic discrepancy),
accept the `dlog` block as a SANCTIONED policy mirror with a written MIRROR_POLICY
exception — recorded honestly as "waived with justification," NOT "independence achieved."
No user escalation (this is a how-it's-checked / coverage decision, not a conceptual one).

## Per-question resolution

**Q1 — is the Sigma_N differential slope singly-routed? — CONCUR.**
The differential slope `2 dln(P/Delta) - dK = dln(Lambda^2/K)` is, across both engines, only
ever computed via the mirrored `dlog` route. `N0 - Lambda^2` (wl:61/py:83) is genuinely
independent but is the STATIC factorization — it does not extract the first-order slope. The
`Common-shape Sigma_N + dK` (wl:113/py:156) and `Xi_load + dK` (wl:118/py:164) checks are
downstream specializations of `sigmaNDirect`, so they inherit (do not add independence to)
the slope route. ⇒ R1 is a legitimate finding; the slope needs an independent second route.

**Q2 — replace vs supplement? — CONCUR with (B) SUPPLEMENT.**
Add ONE new `.wl` line:
`expectZero["Sigma_N - dln(Lambda^2/K) [series route]", 2*dlogSeries[exprPoverDeltaPhys] - kappa - dlogSeries[(lambda^2/k)/.subsEps]]`
and leave the existing `dlog` lines untouched. The new line becomes the load-bearing
independent differential check; the old line remains corroborating mirror coverage. Option
(A) (redefine `sigmaNDirect`/`sigmaNShape` to use `dlogSeries`) is also valid but has
unnecessary blast radius — downstream `sigmaNCommon`/`Xi_load` would inherit the new
extraction path, exposing 3 checks to the feasibility risk instead of 1.

**Q3 — feasibility + anti-guard — CONCUR.**
For a nonzero analytic argument at eps=0 (a ratio of sums/products of `Exp[s_i eps]` is
analytic), `Coefficient[Normal[Series[Log[e],{eps,0,1}]],eps]` is EXACTLY
`(D[Log[e],eps] /. eps->0) = e'[0]/e[0]`. So `dlogSeries` should land the same rational
first-order log-slope as `dlog`, modulo Mathematica simplification strength/assumptions ⇒
the DIRECT-minus-SHAPE difference should still `FullSimplify` to 0. ANTI-GUARD CONFIRMED: the
meaningful check is series-route DIRECT (`2*dlogSeries[P/Delta] - kappa`) vs the SHAPE target
(`dlogSeries[(lambda^2/k)]`), NOT `dlogSeries[e] - dlog[e]` on the SAME argument (that only
validates two extraction methods against each other — a differentiation-method tautology, not
the Sigma_N physics/factorization identity). Keep `-kappa` (= -delta_K) symbolic.

**Q4 — escape clause — CONCUR with honesty caveat.**
Do NOT force a route that makes a previously passing audit silently change meaning or become
brittle. If `dlogSeries` fails to reduce robustly AND the failure is a CAS-normalization
limitation (not a demonstrated algebraic discrepancy), the honest fallback is a written
mirror-policy exception in MATHEMATICA_MIRROR_POLICY. CAVEAT (recorded): that must be logged
as a SANCTIONED MIRROR / accepted-risk closure ("R1's structural-independence ask waived with
justification"), NOT as "independent series-route coverage achieved." The exec-backed answer
comes from the `codex-invoke` run (Codex runs Mathematica); the directive carries both the
primary `dlogSeries` prescription and this mandatory escape clause.

## Directive implications (for the rewrite)

- findings_count = 1 (R1 only). Do NOT re-prescribe F1/F2 (already resolved + PASS; re-touching
  the F1 Sigma_N lines risks reintroducing the simplify-commutes near-tautology that V.1 already
  removed).
- Scope = `.wl` ONLY (the SymPy `dlog` route is the reference engine; do not touch the `.py`).
- Primary edit: define `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps]`
  (with `$Assumptions` as needed) and ADD the single supplemental series-route `expectZero` line
  for the Sigma_N block (series-DIRECT vs SHAPE; keep `-kappa` symbolic). Leave the existing
  `dlog`-based `Sigma_N - dln(Lambda^2/K)` line in place.
- Mandatory escape clause: if the series route does not robustly reduce to `=== 0`, append
  `## Blocked: F1` (banner/observation) and the orchestrator records the sanctioned-mirror
  exception in MATHEMATICA_MIRROR_POLICY honestly (waived, not achieved).
