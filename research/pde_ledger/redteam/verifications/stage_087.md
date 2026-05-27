---
unit_id: 087
batch: III.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 087

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim + banner target_mismatch)

**Classification:** resolved

**What changed:**
- Sympy script banner at `scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:47` now reads `"STAGE 087 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE"` (previously `"STAGE 70 — …"`).
- Mathematica banner at `mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:32` now reads `"STAGE 087 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE"` (previously `"STAGE 070 — …"`).
- Sympy docstring (lines 2–21) explicitly names this stage a checkpoint-consolidation and points readers to upstream sources `scripts/moving_throat_pde_stage081_*`, `stage082_*`, `stage085_quadrupole_demand_cancellation_*`, and `stage086_family1_loading_ratio_window_*` as the locus of the actual cancellation chain.
- Mathematica comment block at `wl:46-56` mirrors this carry-forward narrative naming the upstream stage_085 / stage_086 scripts.

**Assessment:**
Resolution per user direction (a) "status/checkpoint consolidation" is correctly applied in both engines. The banners no longer carry the stale stage-70 label, and the docstring/comment block makes the carry-forward intent explicit. No collateral edits to upstream files or to paper.tex/notes/. No new symbols introduced (consistent with the consolidation classification). The change is non-material to downstream derivations.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- Sympy: an `expect_close(name, value, target, tol)` helper was added at lines 40–44. The previous tautological residual checks `expect_zero("zeta_suff - 2.4662...")`, `…zeta_fail - 2.4675...`, `…zeta_max - 2.4675...` are removed entirely. They are replaced (lines 73–75) by three `expect_close("rho_X^(...) vs stage-086", rho_X, sp.Float("...", 30), sp.Float("1e-13", 30))` calls that cross-check the carried literals against the upstream stage-086 quoted values. The `expect_zero("unblocked zeta_req", ...)` line (line 55) is retained as a definitional consistency probe of the simplification, not a physics anchor.
- Mathematica: three new `expectApprox["rho_suff^(chi) vs stage-086", rhoSuff, 3.46622291347846, 10^-14]` cross-checks were added at lines 61–63 with the same intent (rho_X carried literals vs stage-086 reference values). The `rhoSuff/rhoFail/rhoMax` literals are now declared at extended precision (`\`20` precision marker, lines 57–59) so the cross-check has room to fail on transcription drift.

**Assessment:**
The sympy residuals are now genuine numeric differences (~1e-16 in the saved output: `rho_suff^(chi) vs stage-086 diff = 1.214e-16`, `rho_fail vs stage-086 diff = 1.070e-16`, `rho_max vs stage-086 diff = 5.367e-17`), not algebraic identities — they reflect the high-precision vs target-precision floating delta and would actually fail under any non-cosmetic transcription drift in the literals. The Mathematica `rho_X vs stage-086` checks at lines 61–63 are likewise real anchors. Resolution direction (i) is implemented per the user-resolution markdown.

### F3 — mathematica_transliteration

**Classification:** resolved (won't-fix per F1 = (a))

**What changed:** No code change. The resolution markdown (`redteam/resolutions/batch_III5_paper_alignment.md:51`) records the won't-fix decision: consolidation stages don't need a second independent derivation; both engines running the same algebra by design is acceptable because the algebra is itself a carry-forward restatement, not a fresh derivation.

**Assessment:**
Consistent with the original directive's F3 text ("If F1 direction = (a) (status/checkpoint), F3 can be closed as won't-fix"). The Mathematica script does retain one genuinely independent assertion — `expectZero["d zeta_req exact formula", dZeta - dZetaExpected]` at line 43 — which provides cross-engine coverage on the symbolic derivative shape.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `unblocked zeta_req = 0`  (definitional, OK as retained probe)
- `rho_suff^(chi) vs stage-086 diff = 1.214391841494943541093365767218853890273E-16`
- `rho_fail^(chi) vs stage-086 diff = 1.070317295670973629640094647445676925917E-16`
- `rho_max^(F1)   vs stage-086 diff = 5.366711433453002476816172462459761803208E-17`

Real residuals, well within the 1e-13 tolerance. Non-tautological.

**Mathematica:** exit=0. Notable lines:
- `PASS: d zeta_req exact formula` (independent derivative residual)
- `PASS: rho_suff^(chi) vs stage-086 diff = 0.` / `rho_fail vs stage-086 diff = 0.` / `rho_max vs stage-086 diff = 0.` (extended-precision rho literal matches machine-precision target to within tolerance)
- `PASS: zeta_suff/fail/max numeric check` retained from pre-fix file (see side observations).

**Output freshness:** SymPy script mtime 1779898699 < output mtime 1779899092. Mathematica script mtime 1779898705 < output mtime 1779899180. Both saved `.txt` outputs were regenerated after the script edits.

## Material-change assessment

`material_change`: false.

No symbolic identity, numeric threshold, or boxed result downstream of stage 087 depends on these scripts in their previous form. The edits relabel banners, document upstream provenance, and convert tautological residuals into real anchors against carried literals. Downstream stages 088–096 carry forward `rho_alpha = alpha_req/alpha_mix` as a *symbolic* input, not a residual produced here.

## Side observations (non-blocking)

The Mathematica script retains three legacy `expectApprox["zeta_suff numeric check", zetaSuff, 2.4662229134784601...]` calls at lines 73–75. Because `zetaSuff` is computed as `N[zetaReq /. {rhoAlpha -> rhoSuff, epsBlk -> 0}, 30] = rhoSuff - 1`, and the target literal is the same `rho_suff - 1` typed out to 25 digits, these calls still compare a quantity to (an extended-precision version of) itself. The resolution markdown explicitly stated F2's apply block should add "four calls comparing `zeta_suff/zeta_fail/zeta_max` and the unblocked `zeta_req(eps_blk=0)` against the notes-quoted upstream targets" — the sympy file removed the three tautological zeta_X residuals and added rho_X anchors; the Mathematica file added the rho_X anchors but did not remove the legacy zeta_X residuals. The new rho_X cross-checks at wl:61–63 already provide the same anchor coverage on the load-bearing literals, so the leftover zeta_X residuals are cosmetic noise rather than illusory coverage. Not blocking on this per the verifier's no-new-findings policy and the user's explicit statement that F2 is resolved; flagging for the post-batch sweep if cleanup is desired.

## Verdict justification

All three findings are addressed per the explicit user-resolution markdown at `redteam/resolutions/batch_III5_paper_alignment.md`. F1 (a) is mechanically visible in the relabelled banners and the new docstring/comment naming upstream stage_081/082/085/086 sources. F2 (i) is implemented in sympy (tautological residuals removed, real `expect_close` anchors added) and in Mathematica (new `rho_X vs stage-086` anchors at wl:61–63, real residuals in the exec log). F3 is correctly closed as won't-fix under F1 = (a). Both engines exit 0; the saved outputs are fresher than the scripts; no downstream-relevant numerical result was changed.
