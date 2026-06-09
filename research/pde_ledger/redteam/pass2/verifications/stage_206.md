---
unit_id: 206
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 206

All three findings are non-source-edit items: F1 = paper-side card lag (USER-DEFERRED to P4-51), F2 = a SymPy section-V cross-ref label drift "Stage 239" → "Stage 205" (DEFERRED to the content-keyed numbering pass `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`; math correct), F3 = stale committed SymPy `.txt` (banner "STAGE 189"), refreshed by the orchestrator re-run. The captured Codex diff patch (`exec_logs/stage_206_diff.patch`) is empty — Codex applied nothing, which is exactly what the directive prescribes (`applied: false`, "Apply nothing on this unit until a follow-up directive authorizes a direction"). Verification confirms (A) the SymPy output refresh landed clean (now STAGE 206), (B) the audit disposition holds on the refreshed artifacts, (C) F1→P4-51 and F2→numbering plan are correctly classed deferrals with no Codex/paper edit and no math impact, and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (card Verification field "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to P4-51)

**What changed:**
Nothing in scripts. The diff patch is empty; no edit to `paper/stages/stage_206.tex:11` or `paper/appendices/stage_appendix_part06.tex:918` (both off-limits to the scripts-only red-team and to Codex without authorization). The finding routes the stale card text to the paper-cleanup queue (P4-51) to cite the present, passing `.wl`.

**Assessment:**
Correctly deferred. A complete Mathematica `.wl` exists and passes every check (Mathematica output M1-M7, F3a, F3b all `PASS:`, no `FAIL:`), so the card understates coverage but contradicts no result. This is a prose↔artifact lag with zero math impact — outside the scripts-only scope and properly parked at P4-51. The directive's `## Resolve before fix_loop` correctly held for the user gate rather than letting Codex auto-edit paper/.

### F2 — notes_contradicts_script (SymPy section-V labels "Stage 239", should be "Stage 205")

**Classification:** resolved (correctly routed — DEFERRED to the content-keyed numbering pass)

**What changed:**
Nothing. The three section-V labels remain in the `.py` exactly as reported: line 131 (`# V. Collapse to the Stage 239 quadratic log predictor under exact curvature`), line 133 (`subbanner("V. Collapse to the Stage 239 quadratic log predictor")`), and line 145 (`"Stage 206/239 log-predictor collapse"`), plus the intermediate `print("Stage-239 tau_log2 =")` at line 142. Codex did not touch them (diff empty), consistent with the directive's hold and the routing of this label drift to `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (content-keyed, never offset-sweep).

**Assessment:**
Correctly classed as a numbering cross-ref label drift, not a math defect. The collapse identity itself is verified by both engines: SymPy `Stage 206/239 log-predictor collapse = 0` (output line 100) and the printed collapse target `-2·log(Φ₀)/(L₀ - √(L₀² - 2·L₁·log(Φ₀)))` (output lines 93-98), which is exactly the oriented quadratic log predictor the notes attribute to Stage 205. The "239" is a stale label only — pointing a reader at the wrong upstream stage — and the math is unaffected. Deferral to the dedicated content-keyed numbering pass is the established disposition for this class; no Codex edit here is correct.

### F3 — stale_output (committed SymPy `.txt` banner "STAGE 189")

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent re-run regenerated the committed SymPy output. `scripts/output/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.txt` now banners `STAGE 206 — CERTIFIED RAY RANKING AND LOCAL BRACKETING` (line 3) and `STAGE 206 SYMPY AUDIT PASSED` (line 115) — the stale "STAGE 189" banner is gone. The `.txt` mtime is 2026-06-09T16:51, newer than the `.py` (2026-06-03T15:59). Git status shows only this one `.txt` as modified, confirming the refresh is the sole edit and carries no spurious collateral.

**Assessment:**
Correct and complete; the refresh is the prescribed remedy. Every SymPy result line is `= 0` or the required boolean: §I lines 20-23, §II lines 40-43, §III lines 64, 71, §IV lines 86-88, §V line 105, §VI booleans lines 111-117 (`pairwise ordering theorem = True`, `pairwise ordering negation satisfiable = False`, sieve cases as required). The refreshed `.txt` still carries the "Stage 239" section-V labels (lines 86, 100), which is internally consistent — it faithfully reflects the current (still-deferred) `.py`; the F2 relabel will propagate to the `.txt` when the numbering pass applies it.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed artifacts. The load-bearing root map is derive-vs-posit: the `.wl` runs `Solve` for the quadratic roots and selects the physical branch by its zero-curvature limit (Mathematica output line 14 `M1 solve roots = {(k - Sqrt[-2*c*H0 + k^2])/c, (k + Sqrt[...])/c}`, line 15 `M1 selected root = (k - Sqrt[...])/c`), then confirms the posited closed form equals the Solve-selected branch (line 19 `M1 closed root - Solve-selected branch = 0`, PASS line 20). The `.py` instead posits the closed form and checks the quadratic residual (output line 20 `quadratic residual = 0`). Same independence on the turning root: `.wl` `Solve`+`SelectFirst` (output line 67 `M7 turning roots = {-(Sqrt[2]*Sqrt[H0]/Sqrt[a]), Sqrt[2]*Sqrt[H0]/Sqrt[a]}`, lines 68-71) vs `.py` posited `sqrt(2H0/a)`. Pairwise ordering uses genuinely different methods: `.wl` runs `Resolve[ForAll[...]]` QE (output line 78 `F3a pairwise certified interval ordering = True`) vs `.py`'s constructive slack-gap + `sp.ask` (output line 111). The `.wl` additionally carries strict-sign assertions absent from the `.py` (line 33 `M3 curvature derivative sign = True`, line 41 `M4 strict endpoint descent on simple-root branch = True`) — doing strictly more, not echoing. INDEPENDENT.
- **Non-vacuous checks:** confirmed. The degenerate-envelope collapse (`.py` §II output lines 42-43, `.wl` M5 output line 52 `M5 degenerate-envelope collapse = 0`) is a real substitution, not an identity-by-construction; the F3a "ordering without separation is not tautological = False" (Mathematica output line 81) confirms the QE premise can fail, so the implication is a can-fail statement; the M6 `eta^2 coefficient cancellation = 0` (line 59) confirms the leading-η width law is the genuine content.
- **0 value misalignments:** confirmed. The report reconciles 10 deliverable symbolic values MATCH against the notes; the lone discrepancy is the F2 stage-number *label* (not a value). The collapse target, root map, width law, turning roots, and ordering forms all agree across both engines (e.g. root map SymPy `2*H0/(k + sqrt(-2*H0*c + k^2))` ≡ Mathematica line 16 `(-2*H0)/(-k - Sqrt[-2*c*H0 + k^2]*Sign[k])` under `k>0`).

## Exec log assessment

**SymPy:** exit=0. Notable lines: `quadratic residual = 0` (L20); `implicit derivative identity = 0` (L22); `small-envelope width law = 0` (L64); `Stage 206/239 log-predictor collapse = 0` (L105); `pairwise ordering theorem = True` / `pairwise ordering negation satisfiable = False` (L111-112); terminal `STAGE 206 SYMPY AUDIT PASSED` (L120), `exit_code: 0` (L122). All equalities hold to zero; banner is now canonical STAGE 206.

**Mathematica:** exit=0. Notable lines: `PASS: M1 closed root - Solve-selected branch` (L20); `PASS: M3 curvature derivative sign` (L34); `PASS: M4 strict endpoint descent on simple-root branch` (L42); `PASS: F3a pairwise certified interval ordering` (L79); `PASS: F3a ordering without separation is not tautological` (L81); terminal `STAGE 206 MATHEMATICA AUDIT PASSED` (L94), `exit_code: 0` (L96). Every check prints PASS; no FAIL.

**Output freshness:** confirmed. SymPy `.txt` mtime 2026-06-09T16:51 > `.py` 2026-06-03T15:59, banner now STAGE 206 (F3 stale-189 cleared). Mathematica `.txt` mtime 2026-06-09T16:51 > `.wl` 2026-06-02T11:26, banner already/now STAGE 206. Captured Codex diff patch is empty — no source edit, as the directive prescribes.

## Material-change assessment

`material_change`: false. No source code changed (diff patch empty); the only edit was the regenerated committed SymPy `.txt` (banner relabel STAGE 189 → STAGE 206 + transcript refresh). No derived result changed — every symbolic deliverable is identical to the pre-refresh content — so no downstream unit is affected. The two deferred items (F1 paper card → P4-51; F2 "Stage 239" label → numbering pass) are documentation/label-only and likewise carry no derived-result change.

## Side observations (non-blocking)

The refreshed SymPy `.txt` still carries "Stage 239" at lines 86 and 100, matching the still-deferred `.py`. This is correct and internally consistent — F3's refresh reflects the current script, and F2's 239→205 relabel is owned by the content-keyed numbering pass, which will re-propagate to the `.txt` at that time. Not a defect; noted only so a later reader does not double-flag it.

## Verdict justification

`verified`. All three findings are non-source-edit and correctly resolved: F3's stale "STAGE 189" SymPy banner is cleared by the orchestrator re-run (refreshed `.txt` banners STAGE 206, mtime newer than `.py`, every SymPy check `= 0` and every Mathematica check PASS with no FAIL); F1 is a paper-side card lag correctly USER-DEFERRED to P4-51 with no Codex/paper edit; F2 is a numbering cross-ref label drift ("Stage 239" → "Stage 205") correctly DEFERRED to the content-keyed numbering pass, with the underlying collapse identity verified by both engines. The captured Codex diff is empty, exactly matching the directive's "apply nothing" hold. The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (Solve/SelectFirst/Resolve-ForAll QE + strict-sign assertions the `.py` lacks), the checks are non-vacuous, and reconciliation is 10/10 value MATCH with 0 value misalignments. No regressions; `material_change: false`.
