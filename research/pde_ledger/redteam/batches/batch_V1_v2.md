---
batch: V.1
range: 164-175
total_stages: 12
verified: 12
findings_count: 22
findings_resolved: 22
findings_blocked_legitimate: 0
material_change_count: 1
clean_stages: []
status_only: []
dirty_stages: [164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175]
checkpoints: []
audit_date: 2026-05-28
verify_date: 2026-05-28
status: closed
---

# Red-team batch V.1 — Microscopic logs, drifts, bundle inversion

## Summary

12-stage audit unit for V.1 (`Part V.1 — Microscopic logs, drifts, bundle inversion`). All 12 are dual-engine computational stages; no checkpoints, no status-only units.

22 findings across all 12 stages — 22 resolved, 0 blocked. 12 stages verified end-to-end; both engines exit 0; all outputs fresh.

This batch was unusually homogeneous: **every one of the 12 stages carried a stale `STAGE N−17` script banner** (the legacy Cluster A offset), and **9 of 12 carried a `mathematica_transliteration` finding** — both the dominant patterns the mirror policy predicts. The substantive (non-cosmetic) work concentrated in five stages: 165 (SymPy short-circuited its headline `δln L_W = δln a` via `if False else 1`), 166 (inversion slopes only printed, not asserted), 169 (load-bearing Family-1 numeric combination only printed), 170 (paper_misalignment — the weak-axisymmetric `(1,½,−1)` signature deliverable was unverified), and 175 (a literal `X−X` Σ_N tautology + an unasserted `Σρ=1` aggregate).

**One paper_misalignment (stage 170 F1)** required genuine user resolution — the **first non-cluster paper_misalignment in the v2 run**. It is *not* a Cluster C downstream carry-forward: the `(1,½,−1)` inheritance into the outlet maps δκ_W/δγ_W (with amplitudes κ1/γ1) is stage 170's *own* Sec. 5 deliverable; stages 171/173 verify the same lane *pattern* only on other quantities (microscopic obstructions / load), not on 170's outlet maps. The user chose **direction (a): add the missing check** to both engines (no paper-side edit). Resolved additively — no carried/derived result changed.

No stop-cold verdicts (no `UNFIXABLE`, no `CRITICAL_DOWNSTREAM`). **Thirteenth consecutive batch clear of stop-cold and of any downstream-invalidating change.**

## Per-stage findings tally

| Stage | Status | Findings | Engines | Notes |
|-------|--------|----------|---------|-------|
| 164 | dirty | 1 | SymPy + Mathematica | F1 transliteration → independent `Series`/`Coefficient` route on healing-locked monomials; banner −17 |
| 165 | dirty | 2 | SymPy + Mathematica | F1 insufficient (SymPy `if False else 1` short-circuit of headline δln L_W=δln a → real log-derivative assert); F2 tautological channel-closure demoted to prints |
| 166 | dirty | 2 | SymPy + Mathematica | F1 insufficient (general inversion slopes only printed → 4 asserts added both engines); F2 transliteration → independent `Inverse[Mmat]` route + banner |
| 167 | dirty | 1 | SymPy + Mathematica | F1 banner −17 mislabel |
| 168 | dirty | 1 | SymPy + Mathematica | F1 banner −17 mislabel |
| 169 | dirty | 3 | SymPy + Mathematica | F1 tautological distributive check → 3 per-coefficient numeric asserts vs paper Family-1 weights; F2 transliteration accepted as policy mirror; F3 banner −17 |
| 170 | dirty | 2 | SymPy + Mathematica | **F1 paper_misalignment** ((1,½,−1) signature unverified → user dir (a): added Sec.5 signature block both engines); F2 transliteration → `D[…,eps]`+direct-`Solve` route (removed `du2sym`/`dP0sym` idiom); F3 banner −17 |
| 171 | dirty | 2 | SymPy + Mathematica | F1 transliteration (directive route was tautological → orchestrator rework: collected-literal target + independent `Series` route); F2 banner + docstring |
| 172 | dirty | 2 | SymPy + Mathematica | F1 transliteration → implicit-differentiation route; F2 banner −17 |
| 173 | dirty | 2 | SymPy + Mathematica | F1 transliteration → `Coefficient[Series[…]]` route (removed `d41Hidden` general-solve); F2 banner −17 |
| 174 | dirty | 1 | SymPy + Mathematica | F1 transliteration → `D[…,eps]` perturbation derivation of b01/z01/n01; banner −17 |
| 175 | dirty | 3 | SymPy + Mathematica | F1 tautological `X−X` Σ_N (→ minimal resolution: keep load-bearing `Σ_N−dln(Λ²/K)`); F2 insufficient (added `Ξ_load` Σρ=1 aggregate); F3 transliteration banner + F1 mirror (step3 policy-accepted) |

**Totals:** 22 findings, 12 dirty stages, 0 clean stages, 0 status-only.

## Clusters

### Cluster A — script-banner renumbering (mechanical, −17 offset, all 12 stages)

Every stage's `.py`/`.wl` `banner(...)` carried `STAGE N−17` (164→147 … 175→158). 11 stages were fixed inside their per-finding directives; **4 residual banners not flagged as findings** (164 `.py`, 165 `.py`+`.wl`, 174 `.py`) were mass-fixed in place under the pre-authorized Cluster A "fix all" direction. No notes-H1 or paper-card-banner offsets surfaced this batch.

### Cluster B — body-text forward-stage citation re-attribution

**None this batch.** No notes prose in 164-175 cited legacy upstream stage numbers requiring re-attribution.

### Cluster C — paper-card `Checks` downgrade

**None this batch.** Stage 170's `Checks` item for the `(1,½,−1)` signature was initially a candidate, but content-disambiguation showed the deliverable is stage 170's *own* (not verified downstream), so it was resolved by adding the script check (direction a), not by a carry-forward citation.

## Mathematica mirror policy — transliteration disposition

9 of 12 stages flagged `mathematica_transliteration` (75% — the highest v2 share, reflecting the mechanically-similar linear/rational algebra throughout V.1). Independent second-engine routes added to **7**: 164 (`Series`/`Coefficient` on monomials), 166 (`Inverse[Mmat]`), 170 (`D[…,eps]`+direct `Solve`), 171 (collected-literal + `Series` route, after rework), 172 (implicit differentiation), 173 (`Coefficient[Series[…]]`), 174 (`D[…,eps]` perturbation). **2 accepted as policy mirrors** (no forced rewrite): 169 F2 (grouped-invariant closed form; F1 numeric checks restore substantive coverage) and 175 F3-step3 (`dlogSeries` route; banner + F1 fix applied). Both dispositions are consistent with the mirror-policy default that transliteration is expected, not exceptional.

## Orchestrator catches (rework loop)

1. **Stage 166 round-trip vector residual.** The F2 matrix-inverse cross-check's prescribed `expectZero["matrix round-trip", Mmat . solVec - {…}]` produces a length-4 vector; `expectZero` tests `res === 0`, which is False for a list, so the first run FAILed with residual `{0,0,0,0}` (every component already zero). Scalarized to `Total[(…)^2]`. Re-run clean (19 PASS).

2. **Stage 175 F1 near-trivial cross-check.** The applier's cross-check variant (`2 dln(P/Delta) − 2 dln Lambda` etc.) reduced on run to a simplify-commutes identity (`P/Delta` and `Lambda := simplify((P/Delta)·subs_hat)` are value-equal). Switched to the directive's blessed minimal resolution: keep only the load-bearing `Σ_N − dln(Λ²/K)` (sensitive to `kappa = δK`) and drop the two near-trivial lines in both engines.

3. **Stage 171 F1 tautological directive route.** The directive's prescribed "independent reconstruction" rebuilt `zCombFormula`/`nCombFormula` from the *same* `D[z2,…]` calls that build `zCombExact`, making the bundle checks `X−X`. The first verification caught it (`needs_rework`). Reworked to: collected closed-form literal target (load-bearing — engine's own `D[]` bundle vs the collected coefficients) **plus** an independent `Series`-coefficient second route per bundle. Re-run clean (23 PASS); re-verified `verified` by a fresh agent.

## Verification

All 12 verification files written under `redteam/verifications/stage_*.md`. Final verdicts:

- `verified` (12): 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175.
- `needs_rework` → reworked → re-`verified`: 171 (1 rework loop).
- `blocked_unfixable` (0).

Material change: 1 (stage 170 — *additive* new Sec. 5 weak-axisymmetric deliverable coverage; no carried/derived result altered, zero downstream propagation). The other additive checks (165/166/169/175) were classified `material_change: false`.

## Cumulative

Range 001-175 paper-aligned at v2 depth. **187/253 stages red-team verified (73.9%).** Thirteenth consecutive batch clear of stop-cold; first non-cluster paper_misalignment (170 F1, resolved additively).

Next batch (sequential-audit-chunks rule, awaits explicit user authorization): **V.2 = stages 176-187** ("Load shape, transfer shape, coherent slippage"), 12 stages.
