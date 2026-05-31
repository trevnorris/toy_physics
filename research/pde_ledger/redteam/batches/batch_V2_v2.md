---
batch: V.2
range: 176-187
total_stages: 12
verified: 12
findings_count: 24
findings_resolved: 24
findings_blocked_legitimate: 0
material_change_count: 0
clean_stages: [179, 180]
status_only: []
dirty_stages: [176, 177, 178, 181, 182, 183, 184, 185, 186, 187]
checkpoints: [185]
consult: redteam/codex_reviews/_consult_V2.md
audit_date: 2026-05-30
verify_date: 2026-05-30
status: closed
---

# Red-team batch V.2 — Load shape, transfer shape, coherent slippage

## Summary

12-stage audit unit for V.2 (`Part V.2 — Load shape, transfer shape, coherent slippage`),
the **first batch of the resumed first pass** after the IV.x/V.1 orchestrator-direct integrity
remediation (batches 1–8) closed. All 12 are dual-engine computational stages; **185 is a
checkpoint** (higher verify bar); no status-only units.

**Resumed under the restored contract: Codex is the fix-applier** (the last first-pass batches
III.5–V.1 had run orchestrator-direct with Codex bypassed — the very drift the remediation fixed).
Per-stage loop: clean Claude audit agent → directive → orchestrator review + Claude↔Codex consult
→ `codex-invoke` (Codex applies + runs + iterates ≤600s) → orchestrator independent exec re-run +
output refresh → clean Claude verify agent.

**2 clean** (179, 180 — confirmed by clean-verify agents), **10 dirty**. ~24 findings, **all
resolved, 0 blocked, 0 stop-cold, ZERO paper_misalignment** (contrast V.1's 1 — no user-resolution
items this batch). Both engines exit 0 on every stage; all outputs fresh. **Material change: 0** —
every fix only adds/strengthens checks; no derived value, constant, identity target, or paper number
moved, so no `upstream_stale` propagation.

The two dominant patterns mirror V.1: **`mathematica_transliteration` (7 stages)** and the same
**stale `STAGE N−17` script-banner cluster** (the legacy global-renumber offset, 176→159 … 187→170).
The substantive (non-cosmetic) work concentrated in the "passes-for-the-wrong-reason" findings —
three HIGH-severity (181 F1, 182 F1, 184 F1) and several X−X / vacuous-check de-tautologizations.

**Two verify-wave catches drove iteration-2 reworks** (the verifier's no-rubber-stamp rule earning
its keep): **181 F1** and the **checkpoint 185 F1** (see "Iteration-2 reworks" below). Everything
else verified on the first verify pass.

## Per-stage findings tally

| Stage | Status | Findings | Notes |
|-------|--------|----------|-------|
| 176 | dirty | 3 | F1 tautological (rigidity corollary substituted into the *constructed* Σ_factored → reduce the independent Σ_exact under rigidity); F2 transliteration (comment-only: documents the existing `D[Log]` vs `series` extraction divergence); F3 banner −17 |
| 177 | dirty | 3 | F1 transliteration → independent factored-`(M,I,H)` load-factor route + banner; F2/F3 insufficient — **informational, no edit** (headline grouped collapse / prefactor-slope are trivial-by-linearity; any "fix" would be tautological — confirmed sound by verify) |
| 178 | dirty | 2 | F1 transliteration → Mathematica-native `Coefficient[Series[Log[pA²/dA²]]]` ν_r route (consult Q1: apply, not sanctioned-mirror); F2 banner −17 |
| 179 | **clean** | 0 | clean-verify confirmed; cosmetic banner drift only |
| 180 | **clean** | 0 | clean-verify confirmed; cosmetic banner "163" only |
| 181 | dirty | 3+1 | **F1 insufficient HIGH** (Ξ₁ typed-in, "ζ-drops-out" untested → perturbation-derived Ξ₁ + ζ-blindness; **iter-2** fixed a vacuous ζ-free Mathematica M2); F2 tautological (R₁ self-def → R_target perturbation); F3 transliteration (independent `sSupport` grouping + ported spoiled negative-control); F4 banner −17 |
| 182 | dirty | 3 | **F1 insufficient HIGH** (support-blindness checks differentiated w.r.t. never-wired symbols → structural free-symbol-absence assertion); F2 transliteration → consult-Q3 gauge-`Solve` linear Σ-coefficient route (NOT Blocked); F3 tautological (X−X self-subtractions → notes-§3 equalities) + banner |
| 183 | dirty | 1 | F1 tautological (triple-rigidity tested only the trivial ⟸ → added the ⟹ content via `expect_nonzero` on branch prefactors) |
| 184 | dirty | 2 | **F1 tautological HIGH** (Mathematica drifts hardcoded = targets, `expr−expr==0` → `SeriesCoefficient[Log[composite/ref]]` route); F2 insufficient (dead composites become load-bearing, same root). SymPy correctly untouched (already correct) |
| 185 | dirty (ckpt) | 3 | F1 insufficient (coefficients multiplied onto a proven-zero → **iter-2**: independent slippage-law `Theta₁`/`Xi₁` reconstruction); F2 insufficient (det M_*^(τ,κη,μ)=1+χ0* unasserted → differentiated-minor det check); F3 banner −17 |
| 186 | dirty | 2 | F1 tautological (ε_η-orbit X−X → `Solve`-derived K_η scaling vs paper literal 2C−U); F2 transliteration → M_* rows from `Coefficient[Log[monomial]]` + re-solve |
| 187 | dirty | 2 | F1 insufficient (`Exp[row]`/`Log` X−X → log-ratios built from §1 monomial primitives); F2 transliteration → same, implemented natively (deletes the `Exp[row]` block) |

**Totals:** ~24 findings, 10 dirty, 2 clean, 0 status-only, 0 blocked.

## Clusters

- **Cluster A — script-banner renumbering (−17 offset):** the same legacy pattern as V.1. Present
  on 176/177/178/179/180/181/182/186 (and the checkpoint 185 at −17). Fixed inside the per-finding
  directives (or folded into F1) where dirty; the clean stages 179/180 carry residual cosmetic
  banner drift, recorded as non-blocking side observations (gate no assertion).
- **Cluster B — body-text forward-stage citation re-attribution:** none (red-team is scripts-only;
  no notes prose touched).
- **Cluster C — paper-card `Checks` downgrade:** none (zero paper_misalignment this batch).

## Mathematica mirror policy — transliteration disposition

7 stages flagged `mathematica_transliteration` (177, 178, 181-F3, 182-F2, 186-F2, 187-F2, plus
176-F2 as a documented comment-only near-mirror). **ALL received genuine independent routes — none
accepted as a sanctioned policy mirror this batch** (contrast V.1's 2). Routes added: 177 (factored
`(M,I,H)` load-factor), 178 (`Coefficient[Series[Log[pA²/dA²]]]`), 181 (`sSupport` grouping +
spoiled control), 182 (gauge-`Solve` linear Σ-coefficient extraction), 184 (`SeriesCoefficient[Log]`
drifts, de-tautologizing hardcoded equals), 186 (`Coefficient[Log[monomial]]` M_* rows + re-solve),
187 (monomial-primitive log-ratios replacing `Exp[row]`). 176-F2 is comment-only (the existing
`D[Log,eps]` vs `series` divergence is genuine; documented, not rewritten). The mirror-policy
default (transliteration is expected; rewrite unless sanctioned) held — no sanctioned-mirror
acceptances were needed.

## Claude+Codex consult (`redteam/codex_reviews/_consult_V2.md`, codex_session 019e77af)

Three narrow math-coverage calls, no paper_misalignments, nothing conceptual, none escalated:
- **Q1 (178):** CONCUR — apply the independent `nuFromData` route; do not sanction as a mirror.
- **Q2 (181):** CONCUR — recompute the perturbed quantities from closed forms (the post-`$Assumptions`-reset block drops ζ); don't raise the cap.
- **Q3 (182):** DISPUTE the `## Blocked` default — a concrete gauge-`Solve`/`CoefficientArrays` linear route exists; reserve Block only if it genuinely fails. (Codex applied it; not blocked.)

A second iteration-2 consult (session 019e77e6) settled 185 F1 (see below).

## Iteration-2 reworks (verify-wave catches)

1. **181 F1 — vacuous ζ-free Mathematica M2 (`partial` → resolved).** Codex's Q2 closed-form
   deviation over-applied: it set `t2LoadedPert = t2DirectPert` (ζ-free), so `D[xi1Loaded, zeta]`
   was identically 0 — the support-loaded ζ-cancellation was never exercised on the Mathematica
   side (SymPy's was genuine). Iter-2 rebuilt `t2LoadedPert` from the ζ-bearing `rTargetLoaded`
   perturbation **and added a `FreeQ[t2LoadedPert, zeta]` guard** that fails if ζ is lost. Re-verified
   `verified` (4/4).
2. **185 F1 — checkpoint coefficient not load-bearing (`partial` → resolved).** The re-pointed law
   still multiplied `C_tr,*`/`A_tr,*` onto `(Σ_tr − Σ_tr_compiled)`, a quantity proven zero at
   sympy:166, and the "coefficient form" anchors were byte-identical X−X. The verifier's own proposed
   normalization `coeff·(1/formula)−1` would *also* be X−X (since the coefficient IS that formula) —
   caught on orchestrator review. A focused consult (019e77e6, route (b)) settled the genuinely
   load-bearing anchor: **reconstruct `Theta₁`/`Xi₁` independently from the slippage drifts**
   (`chi1_indep=chi0s·Σ_chi`, `deltaU1_indep=deltaUs·Σ_delta`, …) so the residual becomes
   `(C_typed−C_true)·Σ_tr` (Σ_tr ≠ 0) — a wrong typed coefficient now fails. Re-verified `verified`
   on the checkpoint bar (`wrong_coeff_fails: yes`, both identities hand-checked).

## Orchestrator catches

- **No pre-Codex directive-review catches this batch** (contrast batch 7's two 148 catches) — the
  directive-review agent independently re-derived 185 F2 (`det = 1+chi0s`), 187 F1 (monomial
  log-ratios reproduce the script's rows), and 184's series route, and found no X−X "independent
  route", fabricated literal, or zero-derivative trap.
- **Post-verify catch (185 F1):** rejected the verifier's proposed normalization fix as itself X−X;
  routed to a 2nd consult for the genuine independent anchor. This is the checkpoint higher bar +
  no-rubber-stamp rule working as designed.

## Verification

All 12 verification files under `redteam/verifications/stage_*.md`. Final verdicts:
- `verified` (12): 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187.
- `needs_rework` → reworked → re-`verified`: 181 (1 loop), 185 (1 loop, checkpoint).
- `blocked_unfixable` (0).

Material change: **0** (`material_change: false` on all 12 — additive/strengthening checks only;
SymPy reference engines untouched where the defect was Mathematica-only, e.g. 184).

## Cumulative

Range 001-187 paper-aligned at v2 depth. **199/253 stages red-team verified (78.7%)** (was 187 after
V.1). **First batch of the resumed first pass**, run under the restored Codex-as-fix-applier
contract; zero stop-cold, zero paper_misalignment, zero material change.

Next batch (sequential-audit-chunks rule, awaits explicit user authorization): **V.3 = stages
188-200** (`Part V.3 — Branch observables, isotropic target, home stretch`), 13 stages (checkpoint
200). The planned full end-to-end **second pass** remains a later cross-check, only after the first
pass reaches stage 253.
