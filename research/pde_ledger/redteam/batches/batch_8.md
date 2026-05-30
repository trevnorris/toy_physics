---
batch: 8
pass: IV.x/V.1 orchestrator-direct integrity remediation
range: [175]
total_stages: 1
verified: 1
findings_count: 1
findings_resolved: 1
findings_blocked_legitimate: 0
material_change_count: 0
material_change_stages: []
clean_stages: []
status_only: []
dirty_stages: [175]
checkpoints: []
status_only_candidates: []
consult: redteam/codex_reviews/_consult_batch8.md
audit_date: 2026-05-29
verify_date: 2026-05-29
status: closed
final_batch: true
remediation_complete: true
---

# Red-team remediation batch 8 (FINAL) — stage 175

## Summary

Eighth and **FINAL** batch of the IV.x/V.1 orchestrator-direct integrity remediation.
One stage whose first-pass fix had been applied orchestrator-direct (Codex bypassed) was
Codex-reconciled and re-verified:

- **175** (V.1, `wall_normalized_load_shape`) — the **LAST of the 29 FINDINGS stages**.

With this close, **ALL 29 findings stages flagged by the Codex re-review are now
remediated and `verified`.** The only remaining planned work is the full end-to-end
**second pass** (a fresh adversarial run of the whole pipeline as a cross-check).

175 is a dual-engine computational stage (SymPy + Mathematica). **NO checkpoints and NO
status-only units** in the batch-8 range — `is_checkpoint: false` and
`is_status_only_candidate: false` confirmed in `redteam/MANIFEST.yaml` (V.1 carries no
checkpoints; nearest checkpoints are 096 and 105).

1 finding — **resolved, 0 blocked**. 175 REMEDIATION-`verified` end-to-end; both engines
exit 0/0; committed outputs fresh. **Material change: 0** (`material_change: false`) — the
fix only ADDS a corroborating independent check; no derived value, constant, identity
target, or paper number moved.

## Consult

One Claude+Codex read-only consult (`redteam/codex_reviews/_consult_batch8.md`, Q1–Q4):
**4 of 4 unconditional CONCUR.** No DISPUTES, no conceptual items, nothing escalated to
the user (this is a how-it's-checked / cross-engine-coverage decision, not a conceptual
change to any published claim).

| Q | Topic | Verdict |
|---|---|---|
| Q1 | Is the Sigma_N differential *slope* singly-routed across engines? | CONCUR — yes; `N0 - Lambda^2` is independent but *static*; common-shape/`Xi_load` are downstream of `sigmaNDirect` |
| Q2 | Replace `dlog`→`dlogSeries` vs SUPPLEMENT with one new line | CONCUR — **option B SUPPLEMENT** (smaller blast radius) |
| Q3 | Feasibility of `dlogSeries` + anti-guard | CONCUR — `Series`+`Coefficient` == `D[Log]` for analytic `e`; compare series-DIRECT vs SHAPE (NOT series-vs-D[Log] on same arg); keep `-kappa` symbolic |
| Q4 | Escape clause if it won't reduce to 0 | CONCUR (+ honesty caveat) — fallback = sanctioned MIRROR_POLICY exception recorded as **"waived with justification," NOT "independence achieved"** |

The escape clause was **available but NOT triggered** — the series route landed `=== 0`,
so R1's structural independence was genuinely **ACHIEVED**.

## Per-stage findings tally

| Stage | Status | Findings | material_change | Notes |
|-------|--------|----------|-----------------|-------|
| 175 | dirty | 1 | false | F1 R1: the Mathematica `Sigma_N` differential block was a line-by-line transliteration of the SymPy block (both via the same `dlog = D[Log[.],eps]/.eps->0` primitive). SUPPLEMENTED with a Mathematica-native `dlogSeries` (Series+Coefficient) independent slope route; the existing `dlog` line and the SymPy `.py` left untouched. |

**Totals:** 1 finding, 1 dirty stage, 0 clean, 0 status-only, 0 blocked.

## Finding closed (detail)

### 175 — independent series-route for the Sigma_N differential slope

- **F1 (mathematica_transliteration / R1):** the Mathematica `Sigma_N` differential block
  (`exprPoverDeltaPhys`/`sigmaNDirect`/`sigmaNShape`, wl:95-98) was a line-by-line port of
  the SymPy block — both engines extracted the first-order log-slope via the SAME
  `dlog` primitive (`D[Log[.],eps]/.eps->0` in Mathematica, `sp.diff(sp.log(.),eps).subs`
  in SymPy) — so the differential SLOPE identity `2 dln(P/Delta) - dK = dln(Lambda^2/K)`
  was singly-routed across engines. In V.1 this had been ACCEPTED as a policy mirror (the
  old "175 F3-step3" `dlogSeries` step was waived per MATHEMATICA_MIRROR_POLICY default);
  the Codex re-review REJECTED that waiver for the differential block.

  **Fix (option B SUPPLEMENT):** added a Mathematica-native extractor
  `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps]` (after
  the existing `dlog` helper) and exactly ONE new check immediately after the existing
  `dlog`-based line:
  `expectZero["Sigma_N - dln(Lambda^2/K) [series route]", 2*dlogSeries[exprPoverDeltaPhys] - kappa - dlogSeries[(lambda^2/k)/.subsEps]]`.
  The new check compares the series-route DIRECT slope (`2*dlogSeries[P/Delta] - kappa`)
  against the SHAPE target (`dlogSeries[(lambda^2/k)]`) — NOT `dlogSeries` vs `dlog` on the
  SAME argument (that would be a differentiation-method tautology, not the physics
  identity); `-kappa` (= `-delta_K`) kept symbolic. The existing `dlog` line was LEFT
  UNTOUCHED as corroboration, and the SymPy `.py` was LEFT UNTOUCHED as the reference
  engine.

  `Series`+`Coefficient` is structurally distinct from SymPy's `sp.diff(sp.log(...))`, and
  the residual landed `=== 0`, so the mandatory escape clause (accept the sanctioned mirror
  with a written MIRROR_POLICY exception) was AVAILABLE but NOT triggered — independence was
  **ACHIEVED, not waived** (the converse of the V.1 F3-step3 disposition). Both engines
  exit 0: **SymPy 13 PASS (UNCHANGED — `.py` untouched); Mathematica 13→14 PASS** (+1 = the
  new `[series route]` line). Codex applied iter=1, deviation: none.

  **Note (non-tautology scope, recorded for honesty):** as with the kept `dlog` line, the
  load-bearing degree of freedom in the new check is the `-kappa` (= `-delta_K`) term plus
  the independent extraction-method coverage; the `2 dln(P/Delta)` vs `2 dln(Lambda)` part
  is value-equal (`lambda` is the simplified `(p/delta).subsHat`). The R1 fix adds a
  *structurally independent slope route*, not new physics — which is exactly what the
  finding asked for. A wrong SHAPE target (e.g. `lambda^2*k`) or a wrong `-kappa`
  coefficient breaks the new check.

## Orchestrator catches (directive-review, stage 175)

**None.** Unlike batch 7 (two catches on the 148 draft), the 175 directive draft was clean
on orchestrator review — the consult had fully specified the edit (helper + one
supplemental line, exact live variable names, anti-guards, escape clause), and the
agent-drafted directive encoded it faithfully. No edits were needed before `codex-invoke`.

The orchestrator's independent exec re-run (reliability gate) confirmed both engines exit
0/0 and that the SymPy `.py` and its committed output are byte-identical (untouched
reference engine), with the Mathematica output gaining exactly the two expected lines
(`Sigma_N - dln(Lambda^2/K) [series route] = 0` / `PASS:`).

## Verification

Verification file written at `redteam/verifications/stage_175.md`. Final verdict:

- `verified` (1): 175.
- `needs_rework` (0); `blocked_unfixable` (0).

Material change: 0 (`material_change: false` — additive corroborating check only; SymPy
reference engine untouched).

## Tracker sync

Six prose trackers synced for batch 8:

- `MATHEMATICA_MIRROR_POLICY` — new batch-8 entry: 175's V.1-accepted policy mirror
  (F3-step3 waiver) UPGRADED to a genuine independent `Series`+`Coefficient` route
  (achieved, not waived); 175 added to the Independent-Mirror Set.
- `CHECKPOINT_TRUST_AUDIT` + `CHECKPOINT_CONSTANT_PROVENANCE` — no checkpoints in range; no
  new constant (the `dlogSeries` fix adds no literal; `-kappa`/targets unchanged);
  cumulative checkpoint trust unchanged (105 at `strong`).
- `PAPER_CLEANUP_TRACKER` — **P4-49**: NO new paper-card items (no `.tex`/paper number
  touched); flagged ONE misfiled notes-review artifact
  (`notes/stages/review/stage_175_review.md` body points at pre-renumber
  `stage022_grouped_p2_normalization_bridge` source files) for separate orchestrator/notes
  repair — same class as the batch-7 trio.
- `EM_PROJECTED_INTEGRATION_TRACKER` — no impact (out of EM-projected range).
- `STAGE_VERIFICATION_COVERAGE` — 175 → verified (remediation); the `[series route]`
  independent check now sits alongside `N0 - Lambda^2`, common-shape, and `Xi_load`;
  Mathematica PASS 13→14; `material_change: false`.

## Cumulative — INTEGRITY REMEDIATION COMPLETE (first pass NOT done)

Range 001-175 paper-aligned at v2 depth and remediation-hardened. **All 29 FINDINGS stages
of the IV.x/V.1 orchestrator-direct integrity remediation are now remediated and
`verified`** (batches 1–8). The orchestrator-direct integrity recovery is **closed**.

**This remediation was a DETOUR within the first pass — the first full pass is NOT yet
complete.** The first full pass has reached only stage 175 (V.1); **stages 176–253 are
still PENDING / un-audited** (batches V.2 176–187, V.3 188–200, VI.1 201–218, VII.1
219–230, VII.2 231–242, VIII.1 243–253 = 78 stages — see `redteam/BATCHES.md`).

**Next (awaits explicit user authorization):** RESUME the first full pass at **V.2 =
{176–187}** — where the team was headed when the integrity issue surfaced — then run
forward through V.3 / VI.1 / VII.1 / VII.2 / VIII.1 to stage 253. The planned full
end-to-end **second pass** (a fresh adversarial re-run of the whole pipeline as a
cross-check, per the "full second pass planned" project note; intermediate retrofit
cross-checks were skipped in favor of that single comprehensive re-run) comes only AFTER
the first pass completes at stage 253 — it is NOT the immediate next step.

**Backlog flag (not blocking):** the bulleted "Current Independent-Mirror Set" in
`MATHEMATICA_MIRROR_POLICY.md` is not backfilled for the recent IV.x/V.x remediation stages
that gained independent routes (e.g. 144/147/166/169 and the 105/108/116 batch-1/2 stages);
their dispositions live in the dated prose entries. 175 was added there this batch. A future
cleanup could backfill the bulleted set for the whole remediation run, or leave the prose
entries as the system of record.
