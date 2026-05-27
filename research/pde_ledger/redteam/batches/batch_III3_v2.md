---
batch_id: III.3
batch_label: "Microclosure, gain thresholds, equilibrium, walls"
stage_range: 061-072
stage_count: 12
audit_pass: v2 (paper-grounded)
audit_start: 2026-05-26T17:36:00Z
audit_close: 2026-05-26T19:35:00Z
verdict_summary:
  verified: 12
  blocked: 0
material_change_kind: none
orchestrator_hotfix_stages: [064]
checkpoint_stages: [069]
v1_material_change_carry_forward_clean: [068]
---

# Batch III.3 v2 close

## Per-stage table

| Stage | Status | v2 audit findings | Material change | Notes |
|---|---|---|---|---|
| 061 | verified | 0 | none | clean at v2 |
| 062 | verified | 3 (F1 paper_misalignment, F2 paper_misalignment, F3 mathematica_transliteration) | none | 2 substantive paper_misalignment items resolved (a)/(a) via user gate |
| 063 | verified | 0 | none | clean at v2 |
| 064 | verified | 2 (F1 insufficient_verification, F2 mathematica_transliteration) | none | orchestrator hot-fix (new pitfall #9 candidate: Integrate + constant factoring) |
| 065 | verified | 0 | none | clean at v2 |
| 066 | verified | 1 (F1 insufficient_verification) | none | added 4 SymPy monotonicity sign assertions |
| 067 | verified | 1 (F1 paper_misalignment / banner relabel) | none | orchestrator-applied directly (codex declined due to stale Resolve block) |
| 068 | verified | 3 (F1 tautological_check, F2 mathematica_transliteration, F3 insufficient_verification) | none | **v1 material_change confirmed clean at v2** — Solve-derived Wfail_res/Wfail_match preserved |
| 069 | verified | 0 | none | CHECKPOINT clean at v2 |
| 070 | verified | 2 (F1 mathematica_transliteration, F2 insufficient_verification) | none | codex applied edits via unusual "do NOT run/edit prose" directive instruction |
| 071 | verified | 0 | none | clean at v2 |
| 072 | verified | 1 (F1 paper_misalignment / banner relabel) | none | applied via codex-invoke; 5 string replacements |

## Findings breakdown (13 v2-introduced findings)

- **paper_misalignment**: 4
  - 062 F1: `script_missing_paper_claim` (second equality of boxed `G_micro` + Cauchy-Schwarz bound) → user-gate Q1, direction (a) extend scripts. Downstream 063 directly consumes the factored form.
  - 062 F2: `notes_contradicts_script` (σφ coupling sign convention) → user-gate Q2, direction (a) flip script sign to match notes. Stage 064 already follows the notes' minus convention.
  - 067 F1: `notes_contradicts_script` (banner STAGE 50/050 → STAGE 067 relabel) → unambiguous direction, orchestrator-applied directly.
  - 072 F1: `notes_contradicts_script` (banner STAGE 55/055 → STAGE 072 relabel) → unambiguous direction, applied via codex-invoke.
- **mathematica_transliteration**: 4 (062 F3, 064 F2, 068 F2, 070 F1)
- **insufficient_verification**: 4 (064 F1, 066 F1, 068 F3, 070 F2)
- **tautological_check**: 1 (068 F1)

`mathematica_transliteration` share: 4/12 dirty stages (33%) — down from III.3 v1's 9/12 (75%). Many v1-flagged transliterations were already remediated in v1 and held at v2.

## Resolution methodology

### Paper-alignment user-gate resolutions (Q1, Q2 — both 062)

1. Audit identified 4 `paper_misalignment` items; 2 (banner relabels) routed directly to codex/orchestrator with unambiguous direction; 2 (062 F1, F2) sent to codex-chat for recommendations with destination-verification questions.
2. Codex returned direction (a) for both with cited rationale (063 directly consumes factored form for Q1; 064 follows notes' sign convention for Q2).
3. Orchestrator independently spot-checked citations: confirmed `C_sp_sq` at scripts/063 line 36 + Cauchy saturation at line 117; confirmed `−Λ σ φ` at scripts/064 lines 63 and 146.
4. User approved both as (a).
5. Apply pass via codex-chat: edited both stage 062 scripts (added second-equality + Cauchy assertion, flipped sign, added susceptibility route for F3). Both engines exited 0.

### Stage 064 orchestrator hot-fix (new pitfall #9 candidate)

Codex iter1's F1 fix integrated `hFun*chiSigmaFun^2` and tried to assert `FullSimplify[Integrate[hFun*chiSigmaFun^2] - gPhi^2*Integrate[chiPhi^2/hFun]] == 0`. `FullSimplify` could not pull the constant `gPhi^2` outside `Integrate[...]` with unspecified symbolic `chiPhi[y]` and `hFun[y]`. Residual surface form: `-gPhi^2*Integrate[c^2/h] + Integrate[gPhi^2*c^2/h]` — algebraically zero but doesn't reduce.

Hot-fix: verify integrand equality first via `FullSimplify[hFun*chiSigmaFun^2 - gPhi^2*chiPhi^2/hFun] == 0`, then define `thetaGeneral = lambdaGeneral = gPhi^2*i1Integral` directly. Integrand equality implies integral equality, bypassing Mathematica's symbolic-Integrate limitation.

**Pitfall #9 candidate** (cross-engine but symptomatic in Mathematica): `Integrate[]` with symbolic unspecified functions does not factor constant multipliers. Verify integrands BEFORE comparing integral values, or factor manually.

### Stage 067 / 072 banner relabels

- 067: auditor wrote a `## Resolve before fix_loop` block in the directive even though the audit report explicitly said direction is unambiguous. Codex correctly halted. Orchestrator applied the 3 string replacements directly and marked directive applied.
- 072: codex-invoke applied 5 string replacements cleanly via the standard fix loop; no Resolve block in this directive.

### Stage 070 unusual directive

Auditor's directive contained "Do NOT run python or mathematica. Only edit files." instruction — non-standard. Codex applied F1 and F2 edits as prescribed but did not run scripts (or edit directive prose). Orchestrator confirmed edits via codex_logs/070_iter1.txt and exec'd scripts independently — both engines exit 0.

## v1 material_change carry-forward (stage 068)

Stage 068's v1 verification carried `material_change: true` because the v1 fix lifted `Wfail_res`/`Wfail_match` from postulated values to `Solve`-derived expressions on resonance-corrected premises. The v2 audit re-validated this carry-forward:
- SymPy section 2 still uses `sp.solve` on matched and profile Peclet balances.
- Mathematica section 2 upgraded from `Solve` to `Reduce` (more rigorous, same outputs).
- Four threshold expressions and two band widths symbolically unchanged from v1.
- The four `Pe_req / [matched-or-profile Delta]` expressions identical to v1 → no further cascade.

**Verdict: clean at v2 depth.** Same pattern as stage 060 in III.2 v2.

## Findings vs prior batches

| Category | I.1 v2 | I.2 v2 | II.1 v2 | III.1 v2 | III.2 v2 | III.3 v2 |
|---|---:|---:|---:|---:|---:|---:|
| paper_misalignment | 7 | 3 | 3 | 3 | 2 | 4 |
| insufficient_verification | 0 | 1 | 8 | 5 | 8 | 4 |
| mathematica_transliteration | 1 | 2 | 2 | 1 | 2 | 4 |
| tautological_check | 0 | 1 | 1 | 1 | 3 | 1 |
| symbol_assumption_error | 0 | 1 | 1 | 1 | 1 | 0 |
| hardcoded_result | 1 | 1 | 0 | 1 | 0 | 0 |
| stale_output | 1 | 1 | 1 | 0 | 0 | 0 |
| Other | 0 | 0 | 2 | 1 | 0 | 0 |
| **Total** | **10** | **10** | **18** | **13** | **16** | **13** |

User-redirection rate this batch: **0** (4th consecutive batch with zero redirections — codex's recommendations all held up after orchestrator destination-verification).

## Cumulative coverage as of 2026-05-26 (III.3 v2 close)

- Stages verified: 96 of 253 (across 8 batches I.1/I.2/II.1/III.1/III.2/III.3 + 2 pending v1-only batches III.4 forward).
- Cumulative v2-introduced findings: 80 (10 I.1 + 10 I.2 + 18 II.1 + 13 III.1 + 16 III.2 + 13 III.3).
- Cumulative substantive findings closed (v1 + v2): ~299.

## Notable non-blocking observations

1. **Stage 070 Mathematica typo**: comment says `"(analytic 8/15 = 0.5333)"` for `I_g`, but actual integral of `(d²/dξ² sech)²` is `14/15 ≈ 0.9333`. Print-only; both assertions use `IgNum` consistently so no false pass. Flag for future cleanup.
2. **Saved `.txt` snapshots stale**: stages 066, 067, 068, 072 have `scripts/output/*.txt` and `mathematica/output/*.txt` files older than their script mtimes. The orchestrator's `redteam/exec_logs/` captures are the authoritative post-fix transcripts. User may want to refresh canonical snapshots out-of-band.
3. **Banner drift across stages 061-072**: every stage 061-072 still carries pre-renumbering "STAGE NN" banner strings (e.g., 061 → STAGE 44, 062 → STAGE 45, etc.). Only the directly-flagged ones (067, 072) were relabeled this batch. The audit reports for clean stages noted this cosmetically. Candidate for batched cosmetic sweep later.

## Next batch suggestion

III.4 v2 (stages 073-084, "Family-1 geometry, thresholds, quadrupole"). Notable: high concentration of `hardcoded_result` findings expected (Family-1 numerology cluster 075-084 in v1). Stage 081's directive prescribed `Solve` adjacent to substitutions where `ConditionalExpression` wrapper bit in v1 — codex.md pitfall #4 should preempt that.
