---
batch_id: III.2-v2
label: Part III.2 — Tracking, zeta thresholds, asymmetry, boost (paper-grounded re-audit)
started: 2026-05-26
completed: 2026-05-26
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: false
material_change_stages: []
material_change_kind: none  # All verifier verdicts material_change=false; Stage 050 has a paper-card extension (new boxed S_n^(max) equation) but the script's existing S_n_max formula did not change, so no export-value cascade.
orchestrator_hotfix_stages: [058]
---

# Batch III.2 v2 — Final summary (paper-grounded re-audit)

Fifth batch processed under the v2 paper-grounded auditor (after I.1, I.2, II.1, III.1). All 12 stages reached `verified` with `material_change: false`. **2 paper_misalignment items + 14 script-side findings**. One stage (058) required an orchestrator hot-fix after Codex's iter2 prescription hung sympy for 7+ hours on a symbolic dsolve+BC simplify; replaced with an equivalent kernel-integral identity check. No downstream cascade required.

## Per-stage outcomes

| Stage | Audit findings | Resolution path | Verifier verdict | material_change |
|---|---|---|---|---|
| 049 | F1 tautological_check (Mathematica overlapFormula regression from v1 F2 fix) | fix_loop iter1 | verified | false |
| 050 | F1 mathematica_transliteration, F2 paper_misalignment (paper_missing_script_claim — S_n^(max) ceiling), F3 insufficient_verification (S_0 doubling) | Q1=(a) apply (paper card extended) + fix_loop iter1 (F1, F3) | verified | false |
| 051 | 0 (clean, checkpoint) | n/a | verified | false |
| 052 | 0 (clean) | n/a | verified | false |
| 053 | 0 (clean) | n/a | verified | false |
| 054 | F1 insufficient_verification (window monotonicity certificate) | fix_loop iter1 | verified | false |
| 055 | F1 tautological_check (Mathematica xFloor), F2 symbol_assumption_error (SymPy domain comment) | fix_loop iter1 | verified | false |
| 056 | 0 (clean) | n/a | verified | false |
| 057 | F1 paper_misalignment (script_missing_paper_claim — `∂_Pe ζ > 0`), F2 tautological_check (Solve-based kappaMax/kappaReq/y_req), F3 insufficient_verification (`partial_kappa` sign sweep), F4 mathematica_transliteration (physical-operator A_K derivation) | Q2=(a) apply (Stage-057 script edit) + fix_loop iter1 (F2-F4) | verified | false |
| 058 | F1-F4 all insufficient_verification (IVT bracket existence, Green-kernel BVP, kernel monotonicity, IFT weak-coupling slope) | fix_loop iter2 + **orchestrator hot-fix** (BVP dsolve replacement) | verified | false |
| 059 | F1 insufficient_verification (`Xi_fail <= Xi_suff` ordering) | fix_loop iter2 | verified | false |
| 060 | 0 (clean, v1 material_change held) | n/a | verified | false |

## Findings breakdown

- `paper_misalignment`: 2 substantive (050 F2 user-gate Q1=a, 057 F1 user-gate Q2=a) — one `paper_missing_script_claim` (scripts had a ceiling theorem the paper card didn't advertise; chose to extend the paper) and one `script_missing_paper_claim` (paper claimed Pe-monotonicity scripts didn't verify; chose to verify locally at Stage 057).
- `insufficient_verification`: 8 (050 F3, 054 F1, 057 F3, 058 F1-F4, 059 F1) — dominant v2-introduced category continues; especially heavy on 058 (4 of 4 findings).
- `tautological_check`: 3 (049 F1, 055 F1, 057 F2)
- `mathematica_transliteration`: 2 (050 F1, 057 F4)
- `symbol_assumption_error`: 1 (055 F2 — Codex resolved with a documentary comment)
- `stop_cold`: 0
- iter-2 / remediation passes: 2 (stage 058 codex iter2 — followed by orchestrator hot-fix; stage 059 codex iter2)
- **orchestrator hot-fixes: 1 (stage 058 — heavy BVP dsolve replacement)**

**V2-introduced findings missed by v1: ~16 substantive across 7 stages.** Pattern still script-thinness-dominated (consistent with II.1 and III.1 v2). Stage 058 was the densest finding cluster (4 of 16) — that stage's v1 scripts essentially only checked the Delta closed form numerically at endpoint sweeps; v2 added IVT existence, kernel construction, kernel monotonicity, and IFT slope verification.

Both user-gate items had unambiguous directions:
- **Q1=(a) (050 F2)**: extend paper card. Cross-check confirmed no downstream stage (051-072) consumes `S_n^(max)`, so option (b) "relocate" had no natural target. Paper card now boxes `S_n^{twin}(x;ε) < S_n^{max}(ε) := 1 + (1-ε)/((2n+1)^2 - ε)` with label `eq:app-stage050-Sn-max` and updated Output line.
- **Q2=(a) (057 F1)**: add local Pe-monotonicity script-side check. Cross-check confirmed Stage 056's scripts only verify the identity `dOmega_Pe/dPe = Cov(chi_0,s)/I_W` but NOT the sign `> 0` (the sign argument is in notes prose only). So the carry-forward chain has an upstream script-side gap; option (b) "rely on Stage 056 carry-forward" was not sound. Stage 057 scripts now sweep `Pe in {1/10, 1/2, 1, 2, 5, 10}` at `(kappa, y) = (1, π/4)` with `if val <= 0: raise AssertionError`, plus a carry-forward comment citing notes 056 §4.

## Orchestrator hot-fix (stage 058) — detail

**What Codex iter2 did:** Codex applied F1-F4 with no blocked findings. For F2 (Green-kernel BVP construction), it added a full `sp.dsolve` of `-Phi''(x) + α^2 Phi(x) = Sigma_Pe(x)` in sympy, then `sp.solve` for the two integration constants against the Robin (x=0) / Neumann (x=1) BCs, then `sp.simplify(phi(1) - phi(0) - Delta)`. The Mathematica mirror used `DSolve` + `FullSimplify`.

**What happened:** the sympy run hung at 100% CPU for 7+ hours on the symbolic BVP simplify. The orchestrator killed the process and diagnosed the bottleneck (heavy symbolic machinery on a 2nd-order BVP with exponential RHS over a denominator `Pe² - α²`).

**Hot-fix applied directly by orchestrator** (not via Codex):
1. **Sympy F2 replacement**: replaced the dsolve+BC simplify block (sympy:94-110) with a numerical sweep `Delta = integral(K(x) * Sigma_Pe(x), x, 0, 1)` on 4 concrete `(α, η, Pe)` tuples. `sp.integrate` with concrete substituted values is fast. Same Green-function physical content; same identity verified.
2. **Mathematica F2 removal**: removed the `DSolve` block entirely. The equivalent identity already exists at the pre-existing line 84 (`delta independent integral matches combination form`) which compares `Integrate[kernel*sigmaPe, {x, 0, 1}]` (Green-function side) vs the `Ic`/`Is` combination form (closed-form side). No new check needed.
3. **Pe == α singularity guards**: `Delta` has a removable 0/0 at `Pe = α` (denominator factor `Pe² - α²` cancels in the numerator but `subs()` doesn't take the limit). Sympy monotonicity sweep (L143-145) and IVT sweep (L165-168) now `continue` past collisions. Mathematica emits benign `Power::infy`/`N::meprec` warnings at the same singular points but all `expectZero` assertions still pass.

**New pitfall #8 candidate (consider promoting to `codex.md`)**: heavy-machinery BVP checks via `dsolve`/`DSolve` are not worth the symbolic cost. The natural equivalent is the Green-function identity `Delta = ∫ K · Σ` which can be verified either symbolically (when sympy can integrate it) or numerically on a concrete `(α, η, Pe)` grid. The dsolve path adds nothing the integral path doesn't already cover.

Codex sessions on stage 058 cleared (codex-reset) and the `## Applied: F2` block in the directive now describes the original (replaced) dsolve approach; the divergence is recorded in a new `orchestrator_hotfix_2026_05_26` block in the directive front-matter.

## Cascade observations

No `material_change: true` stages this batch. Stage 050's paper card was extended (new boxed equation), but the script's `S_n_max` formula is byte-identical to before; no downstream stage references `S_n^(max)` so no cascade.

Stage 060 (the v1 `material_change: true` stage carried forward from earlier) returned **clean (0 findings)** under v2 audit — the v1 gain-definition `Xi_micro = Λ²L²/(Θ T_X)` propagated correctly through paper card, notes, and both engines, and into downstream Stage 061 consumer. No v2 retrofit needed for this stage; the v1 material_change has now been confirmed sound at v2 depth.

## Resolution methodology

Same v2 pipeline as I.1/I.2/II.1/III.1 with one new escalation:

1. **Auditor wave (12 fresh agents)** — paper-grounded read of `paper/stages/stage_NNN.tex` + `notes/stages/moving_throat_pde_stageNNN_*.md` + `paper/appendices/stage_appendix_part03.tex` BEFORE opening scripts. Each agent writes `redteam/reports/stage_NNN.md` and (if findings) `redteam/directives/stage_NNN.md`.
2. **Consolidation** — `redteam/resolutions/batch_III2_paper_alignment.md` collected the 2 `paper_misalignment` items with options (a/b/c/skip) for each.
3. **Codex recommendation pass** — fresh Codex session reads paper/notes/scripts/downstream and fills `## Recommendation` blocks; Q1 and Q2 both came back direction=(a) with well-cited rationales.
4. **Orchestrator destination-verification cross-check** — independent grep of (a) downstream consumers of `S_n^(max)` (0 found, confirming Q1=a) and (b) Stage 056's script-side sign assertion for `dOmega_Pe/dPe` (absent — only the identity is checked, not the sign — confirming Q2=a).
5. **User gate** — surfaced both recommendations with destination-verification evidence; user approved (a) for both.
6. **Codex apply session (paper-edit auth)** — applied Q1 paper-card edit + Q2 Stage-057 script edit. Both re-ran scripts and confirmed exit 0 with new PASS lines.
7. **fix_loop on stages with remaining script-side findings** (049, 050 F1+F3, 054, 055, 057 F2-F4, 058, 059) — sequential `$RT codex-invoke <NNN> <directive>` with sanity refresh of both engines.
8. **Orchestrator hot-fix on 058** — replaced unworkable BVP dsolve block with kernel-integral identity, ran scripts manually to confirm clean exit, captured outputs.
9. **Verifier wave (7 fresh agents)** — all returned `verified`, `material_change: false`.

## Findings tally vs. prior v2 batches

| Batch | Stages | Findings | paper_misalignment | insufficient_verification | tautological_check | other | iter-2/3 | orchestrator hot-fix |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| I.1 | 12 | 10 | 7 | 0 | 0 | 3 | 0 | 1 (stage 004) |
| I.2 | 11 | 10 | 3 | 1 | 2 | 4 | 1 | 0 |
| II.1 | 13 | 18 | 3 | 8 | 5 | 2 | 1 remediation | 0 (codex remediation only) |
| III.1 | 12 | 13 | 3 | 5 | 4 | 1 | 1 iter-2 | 0 |
| **III.2** | **12** | **16** | **2** | **8** | **3** | **3** | **2 iter-2** | **1 (stage 058)** |

Cumulative v2 findings to date: ~67 across batches I.1+I.2+II.1+III.1+III.2. The `insufficient_verification` category is now the clear v2 leader (22 of ~67).

## Verification artifacts

- Audit reports: `redteam/reports/stage_049.md` ... `stage_060.md` (all 12)
- Directives: `redteam/directives/stage_049.md`, `050.md`, `054.md`, `055.md`, `057.md`, `058.md`, `059.md` (7 with findings)
- Verifications: `redteam/verifications/stage_049.md` ... `stage_059.md` (7 with findings; clean stages do not get verification files)
- Apply log: `redteam/resolutions/batch_III2_paper_alignment.md` (Q1, Q2 both applied)
- Codex logs: `redteam/codex_logs/batch_III2_recommendation.txt`, `batch_III2_apply.txt`, `049_iter2.txt`, ... `059_iter2.txt`
- Fix-loop log: `redteam/fix_batch_III.2_v2.log` (resumed across 049-055-057-058-059 with one HALT/recovery on 058)
- Hot-fix detail: `redteam/directives/stage_058.md` frontmatter `orchestrator_hotfix_2026_05_26` block
