---
batch_id: III.5
range: 85-90
stage_count: 6
label: Part III.5 — Quadrupole cancellation, loading ratio, verdict
audit_pass: v2 (first-pass under paper-grounded auditor)
audit_date: 2026-05-27
verify_date: 2026-05-27
verified: 6
clean_on_first_read: 1   # 085
dirty: 5                  # 086, 087, 088, 089, 090
findings_count: 15
findings_resolved: 15
user_redirection_rate: 0
material_change_count: 0
orchestrator_hotfixes: 2  # 088 sympy subs canonicalization + 088 mathematica comment *) bug
codex_invocations: 0
---

# Batch III.5 v2 — Red-Team Audit Close

Snapshot date: `2026-05-27`. Batch closed.

## Why this matters

Batch III.5 is the cancellation/closure batch for the Family-1 outgoing branch: it takes the quadrupole-demand product-language cancellation (stage 085), the Family-1 loading-ratio window (stage 086, 087 consolidation), the contact-plus-pole loading-ratio extraction from the minimal isotropic module (stage 088), the explicit Family-1 verdict at zero transport bias (stage 089 CHECKPOINT), and the updated-reduced-status verdict (stage 090 CHECKPOINT). Both checkpoint stages 089 and 090 are paper-grounded `Output` lines that downstream Part IV consumes.

This is the **first** III.5 audit pass (no prior v1 for this batch — auditor used the v2 paper-grounded prompt directly).

## Coverage

| Stage | SymPy | Mathematica | Checkpoint | Pre-Audit | Post-Audit |
|---|---|---|---|---|---|
| 085 | present | present | — | pending | verified, 0 findings |
| 086 | present | present | — | pending | verified, 2 findings (tautological×2) |
| 087 | present | present | — | pending | verified, 3 findings (paper_misalignment + tautological + transliteration) |
| 088 | present | present | — | pending | verified, 3 findings (tautological×2 + transliteration) |
| 089 | present | present | ✓ | pending | verified, 4 findings (paper_misalignment + transliteration + tautological + hardcoded) |
| 090 | present | present | ✓ | pending | verified, 3 findings (tautological + script_missing_paper_claim + insufficient_verification) |

## Finding breakdown

- `paper_misalignment`: 2 substantive (087 F1 consolidation status, 089 F1 chain closure) + 12 banner relabels (orchestrator-direct sweep, every script's banner stale from the global renumber)
- `tautological_check`: 7 (086 F1+F2, 087 F2, 088 F1+F2, 089 F3, 090 F1)
- `mathematica_transliteration`: 3 (087 F3 marked won't-fix per F1=(a); 088 F3 rewrote to paper-form path with `Limit + subtraction`; 089 F2 rederived peSuffChi/peFailChi via FindRoot)
- `hardcoded_result`: 1 (089 F4 — sympy literal Pe values; resolved via provenance comment per pitfall #10)
- `script_missing_paper_claim`: 1 (090 F2 — resolved by F1's new anchors)
- `insufficient_verification`: 1 (090 F3 — added Pe_req = 0 carry-forward in both engines)

15 substantive findings closed across the 5 dirty stages.

## Resolution pattern

Per the III.4 lesson ("Codex stalled mid-Q1 consultation"), III.5 bypassed Codex entirely. The 3 user-gate questions (087 F1, 087 F2 conditional, 089 F1) were consolidated into `redteam/resolutions/batch_III5_paper_alignment.md` with orchestrator-direct recommendations + destination-verification greps. The user approved all three as recommended (a/i/a). Apply was orchestrator-direct (no codex-chat session).

## Hot-fixes (2)

### Stage 088 SymPy — subs canonicalization

The directive prescribed:
```python
Y_rho_u = Y_rho.subs(omega**2/Omega_Q**2, u)
```
But `Y_rho` after `sp.simplify` becomes `(Omega_Q**2*rho_alpha - omega**2)/(rho_alpha*(Omega_Q**2 - omega**2))` — the combined ratio `omega**2/Omega_Q**2` is no longer a syntactic subexpression. Substitution returned `Y_rho` unchanged; downstream `sp.limit((1-u)*Y_rho_u, u, 1)` produced `c1_extracted = 0`, and `c0_extracted - 1/rho_alpha` was non-zero, failing the first assertion.

**Fix** (sympy lines ~55-56):
```python
Y_rho_u = sp.simplify(Y_rho.subs(omega**2, u * Omega_Q**2))
```
Substitute on `omega**2` (which IS a syntactic subexpression of the simplified form) then re-simplify. Same approach applied to `Y_Q_paper_u` extraction. After fix: `c1_extracted = (rho_alpha-1)/rho_alpha`, `c0_extracted = 1/rho_alpha` — both correct.

### Stage 088 Mathematica — `*)` substring closes comment prematurely

The orchestrator's apply included a comment:
```
(* Stage-085 product identity: Pi_tr = rho_alpha * C_mix (verified upstream
   in mathematica/moving_throat_pde_stage085_*). Substitute rho_min. *)
```
The substring `stage085_*)` contains an embedded `*)` that Mathematica reads as the end of the outer comment. The trailing fragment `. Substitute rho_min. *)` was then parsed as code, raising `Syntax::sntx: Invalid syntax in or before "in mathematica/moving_throat_pde_stage085_*). Substitute rho_min. *)"`. Mathematica's parser recovered to subsequent code and reached `Exit[0]` cleanly, so the orchestrator's `rc=0` check passed. However, the F1 assertion `expectZero["Pi_tr_from_rho - (4/3) C_mix", ...]` and the regime trichotomy at lines 90-91 were skipped: only 7 PASS lines appeared in the log instead of the expected 9.

**Verifier caught it** from the missing PASS line count — directly confirming the verifier-prompt warning that "passing exec log is necessary but not sufficient. A script that now passes because Codex made the assertion tautological still fails verification."

**Fix** (wl line 87):
```
   in the stage 085 Mathematica audit files). Substitute rho_min. *)
```
Strip the `*)` substring from the comment body. Re-run produced all 9 PASS lines and the final `Stage 088 Mathematica audit passed.`

### New pitfall #11 candidate

**Mathematica comment bodies cannot contain `*)` substrings.** The first `*)` ends the comment; subsequent text is parsed as code (which usually fails syntax and triggers parser-recovery skip-ahead). The failure is silent at the rc=0 level — the script still exits 0 because Mathematica recovers and reaches the final `Exit[0]`. **The verifier must count expected PASS lines, not just check rc=0**, to detect this class of silent partial run.

Promote to `codex.md` if recurs in IV.1+.

## Banner-relabel sweep (12 scripts)

Every III.5 script carried a stale banner from the pre-renumber numbering. Auditor-flagged for 085 (.py+.wl), 086 (.py+.wl), 087 (.py+.wl), 088 (.py+.wl), 089 (.py+.wl), 090 (.py+.wl) — orchestrator-direct sweep across all 12 in a single pass. Bonus: 089 sympy docstring at line 3 had a stale filename `moving_throat_pde_stage72_*`; fixed to `stage089` alongside.

| Old banner | New banner |
|---|---|
| STAGE 68 / 068 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS | STAGE 085 |
| STAGE 69 / 069 — FAMILY-1 PURE LOADING-RATIO WINDOW | STAGE 086 |
| STAGE 70 / 070 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE | STAGE 087 |
| STAGE 71 / 071 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE | STAGE 088 |
| STAGE 72 / 072 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH | STAGE 089 |
| STAGE 073 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION | STAGE 090 |

## Material-change assessment

Zero `material_change: true` flags. All derived numerics unchanged:
- `rho_alpha = 4/3` ✓
- `zeta_req = 1/3` ✓
- `Pi_tr = (4/3) C_mix` ✓
- `Pe_req = 0` ✓
- Family-1 window `rho_suff^(chi) = 3.46622291347846`, `rho_fail^(chi) = 3.46752913273870`, `rho_max^(F1) = 3.46752922945601` ✓
- `A_F1 ≈ 1.00005192880220` ✓

Edits restructure verification paths (independent extraction in 088, chain closure in 089, derive-from-coefficients in 090) and add cross-engine anchors. No downstream Part IV consumer sees a value change.

## User-gate resolutions

All three approved as orchestrator-recommended:
- **Q1 (087 F1) — (a) status/checkpoint consolidation**: paper card explicitly designates 087 as consolidation; cancellation chain proven upstream.
- **Q2 (087 F2) — (i) expect_close vs upstream stage 069/082**: same approach as 086 F1, strictly stronger than print demotion.
- **Q3 (089 F1) — (a) strengthen scripts**: chain closure with sp.limit/Mathematica Limit is ~6 lines per engine; honors checkpoint bar; downstream consumers in Part IV depend on Pe_req=0 being established.

Six consecutive zero-redirection batches (II.1, III.1, III.2, III.3, III.4, III.5).

## Codex availability note

Per the III.4 stall lesson, III.5 bypassed Codex entirely. The orchestrator-direct math-authority workflow held up because:
1. Audit evidence was already paper-grounded and conclusive.
2. The 3 user-gate questions each had a clear (a/b/c) option set with explicit destination-verification greps.
3. The substantive edits were specified in detail by the auditor directives.

For IV.1+ the orchestrator should re-test Codex availability (fresh sessions appearing under `~/.codex/sessions/YYYY/MM/DD/` within ~30s of invocation). If Codex is healthy, the math-authority delegation is preferred for novel derivations; if not, orchestrator-direct is the proven fallback.

## Output paths

- Reports: `redteam/reports/stage_{085..090}.md`
- Directives: `redteam/directives/stage_{086..090}.md` (085 had no directive — clean)
- Verifications: `redteam/verifications/stage_{085..090}.md`
- Resolutions: `redteam/resolutions/batch_III5_paper_alignment.md`
- Audit prompts: `redteam/tmp_prompts/audit_prompt_{085..090}.md`
- Verify prompts: `redteam/tmp_prompts/verify_prompt_{085..090}.md`
- Updated scripts: `scripts/moving_throat_pde_stage{085..090}_*_sympy_audit.py`
- Updated Mathematica: `mathematica/moving_throat_pde_stage{085..090}_*_mathematica_audit.wl`
- Refreshed outputs: `scripts/output/...` and `mathematica/output/...` (all 12 mtimes post-fix)
- Updated MANIFEST: `redteam/MANIFEST.yaml` (status: verified, iteration_count: 1, last_verify_date: 2026-05-27 for all 6)

## Cumulative coverage state

102 of 253 stages red-team verified (was 96 → +6 from III.5). Entire range 001-090 now paper-aligned at v2 depth. Next batch: IV.1 (stages 091-102, "Part IV.1 — Grouped p2 geometry, decoupling, contamination") — 12 stages, includes no checkpoints per BATCHES.md.
