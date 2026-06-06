# Pass-2 Batch III.5 — stages 085–090 (Part III.5: quadrupole cancellation, loading ratio, verdict)

**Result: 6/6 verified, all `material_change: false`, 0 stop-cold/blocked, 0 Codex deviations
beyond one sanctioned helper-add. Two checkpoints in range (089 & 090) — BOTH cleared the
higher bar.**

## Headline
- **NO genuine `paper_misalignment` anywhere → ZERO paper/notes edits.**
- **Value reconciliation: 59 deliverable values batch-wide, 0 misaligned** (085=13, 086=13,
  087=5, 088=8, 089=14, 090=6).
- 3 clean stages: **085, 086, 090** (090 = checkpoint, clean, higher bar cleared by a
  substantive own-assertion `zeta_req = rho_alpha − 1` + branch-ordering inequalities).
- 3 findings stages: **087** (1), **088** (1), **089** (2). All `needs_user_resolution: false`
  (script-side fixes only).

## Findings & fixes

### 087 — F1 `insufficient_verification` (low; .py + .wl)
The "cross-check against upstream stage-086" asserts compared each window literal
(`rho_suff/rho_fail/rho_max`) against the SAME literal re-typed in the same file (py:58 vs
py:73; wl:57 vs wl:61) — a hollow self-comparison that cannot detect a mistyped literal (both
sides move together). Output records confirmed it (sympy diff ~1e-16, mma `diff = 0.`). The
comments overclaimed an "upstream cross-check."
**Fix (both engines):** replaced the three self-comparisons with can-fail structural relations
— `rho_suff < rho_fail`, `rho_fail < rho_max`, `0 < rho_max − rho_fail < 1e-6` (the tight
constructive-ceiling gap = 9.671731e-8) — and reworded the overclaiming comments. The `.wl`
`zeta_*` numeric checks (wl:73-75) were ALSO self-referential (target = `rho_X − 1` re-typed)
→ re-anchored to `rhoSuff − 1 / rhoFail − 1 / rhoMax − 1` so they genuinely test the
`epsBlk→0` substitution of `zetaReq`. Load-bearing `unblocked zeta_req` (py:55/wl:44) and
`d zeta_req exact formula` (wl:43) untouched; the three literal VALUES unchanged. **iter-2
followup:** the SAME overclaim survived in the SymPy top docstring (py:17-20) → reworded to
match (docstring-only, non-printing; transcript byte-identical to the post-F1 run). `.wl`
verify-confirmed STILL INDEPENDENT (not a transliteration). 8 PASS (mma).

### 088 — F1 `stale_output` / stale self-label (low; .py docstring only, 0-seat)
SymPy docstring carried the pre-renumber OWN-number self-label: py:3
`...stage71_...` and py:5 `SymPy audit for Stage 71.` (088 − 17 = 071 EM-extension drift).
**Fix:** → `stage088` / `Stage 088` (docstring-only). Banner (py:29), filename, paper card,
`.wl` banner all already canonical 088; `stage085` cross-refs (py:112, wl) correct → untouched.
Committed outputs **byte-identical** (docstrings are not printed). **Both first-pass
fragilities confirmed ALREADY FIXED (no new finding):** (a) the `omega**2 → u` substitution
now lands on the atomic `omega**2` in `Y_rho` (output shows genuine extracted residues
`c0=1/rho_alpha`, `c1=(rho-1)/rho`, not the `0` a no-op would give); (b) the upstream ref is
`stage 085` (space, no `_*`) so no embedded `*)` closes the comment early — assertion count 9 =
PASS-line count 9 (no silent partial run).

### 089 — F1 + F2 `tautological_check` (CHECKPOINT, higher bar; de-tautologized BOTH, not deleted)
- **F1 (SymPy line 77 + .wl additive):** `Q` was baked with `eps_blk = sp.Integer(0)`, so
  `Q ≡ 1 + zeta` and `expect_zero(Q − (1 + zeta))` was a pure X−X self-cancel (the script's own
  comment 87-89 flagged it). **Fix:** made `eps_blk` symbolic, introduced general
  `Q_gen = (1+(1−2 eps_blk) zeta)/(1 − eps_blk zeta)`, set `Q = simplify(Q_gen.subs(eps_blk,0))`
  (still `1+zeta` → downstream `rho_*` byte-identical), and asserted the eps→0 reduction on the
  GENERAL form: `Q_gen.subs(eps_blk,0) − (1+zeta)` (a structural transcription error in Q now
  fails). Added a parallel `.wl` reduction assert `expectZero["Q(zeta;0)=1+zeta reduction",
  q[zetaRed,0] − (1+zetaRed)]` for engine symmetry. **Sanctioned Codex deviation:** the `.wl`
  had no `expectZero` helper (only `expectTrue`/`expectApprox`) → Codex added a minimal one
  (`FullSimplify[Together[Expand[expr]]]`, pass iff `=== 0`); verify-confirmed correct.
- **F2 (both engines):** the boxed Output `Pe_req = 0` was verified by a literal checked against
  itself (`expect_zero(Pe_req)` with `Pe_req=sp.Integer(0)`; `expectApprox[peReq,0]`) = `0==0`.
  **Fix:** removed the self-check; replaced with a CAN-FAIL positivity assertion on the
  zero-bias success margin `zeta_F1(0) − zeta_min` (= A_F1 − 1/3 = 0.6667185954688619953260008,
  the quantity that FORCES Pe_req=0), then constructed/printed `Pe_req = 0` as the consequence.
  ⭐ **ORCHESTRATOR FALSE-POSITIVE/SAFETY GUARD:** the audit agent's drafted directive proposed
  a `sp.Piecewise((0,cond),(sp.nan,True)) → expect_zero` gate — FRAGILE, because in SymPy
  `abs(complex(sp.nan)) > tol` evaluates to `False`, so a failed precondition would pass
  SILENTLY. The orchestrator rewrote the directive to the explicit-raise margin form; Codex
  applied it correctly (explicit `raise`/`fail`, not the silent-nan gate).
- **Checkpoint verdict: clears the higher bar.** After the fixes, both engines carry substantive
  non-tautological assertions, agree to displayed precision, paper alignment exact (boxed
  `Pe_req = 0` Output line retained); NO remaining named tautology; NO new pinned constant
  (`A_F1`, `rho_*` thresholds derived/carried with provenance; the Mathematica side independently
  re-derives the upstream `Pe` via `FindRoot` — robust route, not the latent `nsolve`-near-`tan`
  pitfall #10). 12 PASS (mma; was 11: +1 reduction +1 margin −1 `peReq`).

## Infra / process
- **Seat policy held:** 088 = 0-seat (.py docstring); 087 + 089 = 2 `.wl`-touching Codex
  sessions (within the 2-seat cap), ran concurrently; 087 iter-2 = 0-seat. Orchestrator `exec-*`
  run SEQUENTIALLY after all Codex done (no overlap with any build).
- **Output refresh:** `$RT exec-*` writes `exec_logs/` only → orchestrator sed-refreshed every
  committed `output/*.txt`. 088 sympy + mma outputs byte-identical (non-printing change). All
  exec exit 0.
- **Arbiter grep on refreshed committed outputs: CLEAN** — no pre-renumber self-epoch label
  (068–073 band) appears as a self-banner/closing in ANY of the 6 III.5 committed outputs
  (findings + clean). Only canonical `STAGE 0{85..90}` banners + `app-stage089-Pe-zero`
  (canonical self paper-eq ref).
- **Deferred numbering CROSS-refs** (content-keyed, never offset-sweep → `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`):
  086 py:37 `Stages 63-64`→080/081; 089 comments `Stage-62`/082/075/074 (cross-refs);
  090 notes provenance-tracker `Stage 73` (out of audit scope). 088 `stage085` refs are correct
  (no defer). Codex held this scope (0 cross-ref touched).

## Trackers
6 prose trackers synced (PAPER_CLEANUP **P5-09** = ZERO paper/notes edits). Pass-1
`MANIFEST.yaml` untouched (isolation held).

## NEXT
**HALT at the III.5 boundary for the user gate.** NEXT = IV.1 (091–102) — first IV.x batch
(`$RT next-batch` → IV.1).
