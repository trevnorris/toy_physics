---
batch_id: I.2-v2
label: Part I.2 — Maxwell bridge, parent throat action, reduced one-port (paper-grounded re-audit)
started: 2026-05-25
completed: 2026-05-25
stages_total: 11
stages_verified: 11
stages_blocked: 0
material_change_any: true
---

# Batch I.2 v2 — Final summary (paper-grounded re-audit)

Second batch processed under the v2 paper-grounded auditor (first was I.1). All 11 stages reached `verified`. 4 stages flagged `material_change: true` (013, 014, 015, 018 — script-side scope reductions that removed duplicated content; downstream consumers should rebind to the destination stages that genuinely own each piece).

## Per-stage outcomes

| Stage | Audit findings | Resolution path | Verifier verdict | material_change |
|---|---|---|---|---|
| 013 | F1 paper_misalignment, F2 insufficient_verification (Xi anchor), F3 insufficient_verification (K1/H_even guards) | Q1=b trim apply + F2 fix_loop; F3 blocked-legitimate (K1/H_even trimmed) | verified | **true** |
| 014 | F1 paper_misalignment, F2 tautological_check (comp_surface round-trip) | Q2=b trim apply + F2 fix_loop | verified | **true** |
| 015 | F1 paper_misalignment (high — wall-only/Y20/grouped duplicates 017), F2 insufficient_verification (concrete IBP parity), F3 tautological_check (wall-only K1/H_even), F4 mathematica_transliteration | Q3=b trim apply + F2 fix_loop + F4 iter2 EL re-derivation; F3 blocked-legitimate | verified | **true** |
| 016 | 0 (clean) | n/a | verified | false |
| 017 | 0 (clean) | n/a | verified | false |
| 018 | F1 paper_misalignment (medium) | Q4=b trim apply | verified | **true** |
| 019 | F1 symbol_assumption_error | fix_loop | verified | false |
| 020 | F1 paper_misalignment (medium), F2 tautological_check (lambda_0 self-ratio, contingent) | Q5=b trim apply (F2 cascade-resolved) | verified | false |
| 021 | F1 paper_misalignment (script_missing_paper_claim — inverse direction) | Q6=a add apply + iter2 closed-form rework | verified | false |
| 022 | 0 (clean) | n/a | verified | false |
| 023 | F1 tautological_check (N0_target self-substitution), F2 insufficient_verification (Schur derivation) | fix_loop F1 + codex remediation F2 (option α, sign convention from paper §2) | verified | false |

## Findings breakdown

- `paper_misalignment`: 6 (Q1-Q6 across stages 013, 014, 015, 018, 020, 021)
- `tautological_check`: 3 (stage 014 F2, stage 020 F2 cascade-resolved, stage 023 F1)
- `insufficient_verification`: 4 (013 F2 Xi anchor, 013 F3 K1/H_even guards blocked-legitimate, 015 F2 IBP profile, 023 F2 Schur derivation)
- `symbol_assumption_error`: 1 (019 F1 Reals)
- `mathematica_transliteration`: 1 (015 F4 — K_eta block resolved via EL re-derivation; wall-only portion auto-resolved by F1 trim)

**V2-introduced findings missed by v1: 6 substantive** (the six paper_misalignment items). Pattern this batch was notably different from I.1: I.1 was dominated by **a** (expand paper) and **c** (acknowledgement) recommendations. I.2 was dominated by **b** (trim scripts), because the cross-stage check revealed duplication: stage 015's wall-only/Y20/grouped blocks duplicated stage 017; stage 018's one-pole/gate/Xi_1 duplicated stages 019/020; stage 013's δP_n/sieve duplicated stages 010/014; stage 020's Y20 duplicated stages 010/017. The duplication came from the EM-projected scripts being file-for-file ports of `notes/em_projected/step_NN_*` master notes, while the paper cards were later compact summaries that distributed content across multiple stages.

## Resolution methodology

Same pattern as I.1:
1. **Questions session**: `redteam/resolutions/batch_I2_paper_alignment.md` + `codex_prompt_batch_I2.md`. Codex answered all 6, skipped 0.
2. **User review with critical re-direction**: User pushed back on Codex's (c) acknowledgement recommendations for Q1 and Q2, citing "each stage builds on prior" principle. Claude cross-checked Codex's "already owned by stage X" claims by reading destination paper cards + scripts. Found: stage 017's script literally contains the same `real_y20_square_ratio`, `D01_full`, `K1_full`, `wall_only` patterns as 015; stage 019's script verifies one-pole closure; stage 020's script verifies even-gate determinant + Xi_1. Codex was right that 015/018's content was duplicated, but wrong to recommend acknowledgement (c) for 013/014 — the cleaner answer was (b) trim, with destinations stage 010 (δP_2/δP_4 added in I.1 v2 sweep) and stage 014 (sieve) actually owning the content. User-revised: Q1 and Q2 flipped from (c) to (b) before apply.
3. **Apply session**: `codex_apply_batch_I2.md` with explicit **destination-verification guardrail** (Codex must grep destination script to confirm equivalent assertion exists before deleting from source). All 6 applied (`0 revised, 0 blocked`). Each apply block records `destination_verified: yes — <file:line>`.
4. **Output refresh**: serial Mathematica + parallel sympy for 6 modified stages.
5. **fix_loop**: 5 stages (013, 014, 015, 019, 023) processed sequentially per `fix_parallelism: 1`. Iter1 each except 015 (iter2 for Mathematica F4 fix — see below).
6. **Verifier wave**: 8 agents in parallel (stages with findings). Returned 6 verified, 1 needs_rework (021), 1 blocked_unfixable (023 F2).
7. **Codex remediation pass**: combined prompt for 021 F1 iter2 (closed-form RHS) + 023 F2 (sign-convention decision from paper). Codex read the paper for 023, chose option α (`+R U W` → `−R` off-diagonal in frequency space → existing `+2*g_U*g_W*R_mix` numerator is correct), added Schur derivation matching existing code. Both clean after remediation.

## Notable per-stage detail

**Stage 015 F4 iter2**: Iter1 introduced Mathematica K_eta via Euler-Lagrange linearization using `Dt[..., Constants -> {...}]`. The `Dt` call produced an unevaluated `Dt[TwR0, w, Constants -> ...]` residual because the Constants directive didn't fully isolate `TwR0` from `w` dependency. Iter2 fix: switched to ordinary `D` with explicit temporary `twR[w]` slot variable, side-stepping the `Dt` lurking-dependency trap. The resulting derivation is genuinely independent of SymPy's `Series` path (different intermediate representations; the residual collapse to zero arises from definitional IBP identity, not from sharing intermediate variables). Worth noting as a potential pitfall #6 candidate, though the trigger (using `Dt` with `Constants` for EL derivations) is narrow.

**Stage 013 F3 blocked-legitimate**: F3 (K_1/H_even coefficient guards for the 1/9, 2/3, -1/27 literals) was correctly blocked by fix_loop because the K_1/H_even definitions no longer exist in trimmed stage 013 (they live in stage 014 after Q1's trim). The coefficient-guard concern is **displaced to stage 014's scope** — the verifier flagged this as a follow-up: stage 014 still has hardcoded coefficients `1/9`, `2/3`, `-1/27` in its K_1/H_even definitions without coefficient-specific anchor checks. Optional future fix_loop on stage 014 could add the guard pattern from the original F3 prescription, adapted to stage 014's variable scope. Not blocking this batch.

**Stage 023 F2 resolution path**: Codex correctly refused to apply the directive's `+Rmix` prescription in iter1 (detected sign conflict with existing `Q_expr`'s `+2*g_U*g_W*Rmix`). Remediation pass had Codex read `paper/stages/stage_023.tex` `eq:app-stage023-full-lagrangian` and notes file to determine the sign by physics, not by directive fiat. Result: paper's `+R U W` Lagrangian → frequency-space spring matrix off-diagonal is `-R` → existing code is right → applied Schur derivation with `-R` matrix off-diagonal, matching the existing `Q_expr` cross sign. This is a clean instance of "Codex as math authority" working as designed: directive was wrong, Codex deferred to the paper, sign was resolved by re-deriving from the Lagrangian.

## Process observations

1. **User re-direction was load-bearing this batch.** Codex's initial (c) for Q1/Q2 would have left duplicated assertions in place; user's "each step builds on prior" principle pushed toward (b) which produced a much cleaner result. Cross-checking Codex's "destination owns X" claims by actually reading the destination paper + script (not just trusting the recommendation) caught the under-recommendation. Pattern worth preserving for I.3+ batches: when Codex recommends (c) acknowledgement, ask "is the destination ALREADY both paper- and script-anchored? if yes, why not (b) trim?"

2. **Destination-verification guardrail in the apply prompt prevented orphan trims.** Adding `"GREP destination script to confirm equivalent assertion exists before deleting from source"` to the apply prompt's procedure meant Codex's `destination_verified:` blocks made the safety argument explicit. No orphan trims occurred; if a destination paper claim existed without a destination script assertion, Codex would have blocked rather than silently delete coverage. (None did, this batch.)

3. **Sequential fix_loop is slower but cleaner.** `fix_parallelism: 1` meant ~5 minutes per stage × 5 stages = ~25 minutes serial fix_loop. Parallel fix_loop would race on Mathematica single-seat + MANIFEST writes. Fine trade-off; the verifier wave is parallel so the overall batch wall-clock is still hours not days.

4. **Verifier wave catches under-applied fixes.** Stage 021's tautological RHS slipped past Codex apply ("Dcorr.subs(Gamma_port, ...) - (-I * N0 * a^5/(27 c_s^5) * omega^5)" reduces to X-X=0 since Dcorr itself contains the N0 factor). The verifier agent caught it and prescribed the closed-form RHS, which then went through a remediation pass cleanly. Worth a paragraph in `prompts/codex.md`: when adding "composed assertion" of two upstream pieces, neither side may use the bare symbol of the upstream piece — must use the closed form on at least one side.

5. **One toolchain pitfall worth tracking (not yet promoted to codex.md):** `Dt[..., Constants -> {list}]` may leave unevaluated `Dt[symbol, w, ...]` residuals if `symbol` has a lurking dependency on `w` that isn't in the Constants list. Workaround: use ordinary `D` with explicit slot variables (`twR[w]` etc.) rather than `Dt` with Constants. Narrow trigger (Euler-Lagrange in Mathematica using Dt for the chain rule); leave for now and document if it recurs.

6. **Cumulative v1+v2 finding density**: I.1 v2 added ~10 net findings vs v1 (7 paper_misalignment + 3 new script-side). I.2 v2 added ~10 net findings vs v1 (6 paper_misalignment + 4 script-side). Roughly 1 v2 paper-misalignment per ~2 stages; ~1 script-side miss per ~4 stages. Stable rate so far.

## Tracker updates landed in this commit

- `notes/MATHEMATICA_MIRROR_POLICY.md`: snapshot bumped to I.2 v2 close.
- `notes/CHECKPOINT_TRUST_AUDIT.md`: snapshot bumped; checkpoints in I.2 are 022, 023 — both verified.
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: snapshot bumped; no new constants surfaced from trims (only duplicates removed).
- `notes/PAPER_CLEANUP_TRACKER.md`: new P3-11 row for batch I.2 v2 sweep; change log entry.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: v2 paper-grounded re-audit row added (stages 013-021 are in this tracker's scope; v2 sweep significantly tightened script-vs-paper alignment).
- `notes/STAGE_VERIFICATION_COVERAGE.md`: snapshot bumped to "batch I.2 v2 close"; cumulative count unchanged (84/253 — same stages, deeper verification).

## Open follow-ups (not blocking)

- **Stage 014 coefficient guards for K_1/H_even literals 1/9, 2/3, -1/27** — displaced from stage 013 F3 by the Q1 trim. Could be added in a future stage 014-only fix_loop. Auditor's original F3 prescription is in `redteam/reports/stage_013.md`; adapt to stage 014's variable scope.
- **Stage 015 docstring** still references the nonexistent `step_13_parent_throat_action_master_notes.md`. Cosmetic; can be cleaned up in any future edit. (Same pattern for stage 018 referencing nonexistent `step_16_*_notes.md`.)
- **Pitfall #6 candidate (Mathematica Dt+Constants residual leak)** — leave undocumented for now; promote to `codex.md` if it recurs in I.3+ batches.
