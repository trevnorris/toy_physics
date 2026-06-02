---
batch: retro
range: 121-123
total_stages: 3
verified: 3
findings_count: 3
findings_resolved: 3
findings_blocked_legitimate: 0
material_change_count: 0
clean_stages: []
status_only: []
dirty_stages: [121, 122, 123]
checkpoints: []
consult: none
audit_date: 2026-06-01
verify_date: 2026-06-01
status: closed
---

# Red-team retro-sweep {121, 122, 123} — dual-engine retrofit

## Summary

3-stage RETRO-SWEEP, not a new audit unit. Stages 121, 122, 123 were already
`verified` in batch IV.3 but were **SymPy-ONLY** (no `.wl`). Under the
dual-engine rule (a Mathematica `.wl` is **REQUIRED wherever Mathematica CAN
independently verify** — the test is "**is it possible**," NOT "is it
necessary"), each was retrofitted with a **NEW independent-route Mathematica
`.wl`**. **Codex designed + wrote the `.wl`; Claude reviewed (audit + verify).**

Outcome: **3/3 verified**, **`material_change: false` on all 3**. The SymPy
`.py` reference engines were **UNCHANGED**; **paper/notes UNTOUCHED**. **0
transliterations accepted, 0 iteration-2 reworks, 0 blocked, 0 stop-cold.**
The single "finding" per stage was the missing-`.wl` dual-engine gap itself,
closed by the new independent route. **No checkpoints in range. No new
constants. No EM-projection change.**

This sweep is "verified-on-disk / uncommitted" at the time of this record.

## Headline — the dual-engine gap for already-verified non-status-only stages is now CLOSED

Per the inventory, 121/122/123 were the **only 3** already-`verified`
non-status-only stages still missing a `.wl`. With the three new independent
routes landed, that gap is now **CLOSED**. The 11 status-only single-engine
stages compute nothing a `.wl` could check and **legitimately remain
single-engine**.

**Cumulative coverage is UNCHANGED at 200/253 verified (79.1%)** — a
retro-sweep adds second engines to already-`verified` stages, it does **NOT**
add new verified stages.

## Per-stage tally

| Stage | Name | New `.wl` | PASS | Notes |
|-------|------|-----------|-----:|-------|
| 121 | geometric_r_selection | `mathematica/moving_throat_pde_stage121_geometric_r_selection_mathematica_audit.wl` | 6 | M1 `r_geom` closed-form derived by `Solve`-ing the Stage-99 length law + explicit positive-branch `Select`; M2 tube-length round-trip; M3 `r_F1` surd certified via `FullSimplify`/`RootReduce`; M4 `r_c`; M5 `Omega_W` definitional-parity (acknowledged low-information); M6 threshold exact-0. Genuine independent route. |
| 122 | mouth_source_compensation_test | `mathematica/moving_throat_pde_stage122_mouth_source_compensation_test_mathematica_audit.wl` | 9 | `g_±` DERIVED by `Solve`-ing the compensation quadratic at `r_F1` + sign-based branch selection (not hard-typed). The traction-ratio de-tautology is preserved cross-engine: the constant `C` (`cStage`) is carried as a FREE POSITIVE SYMBOL, the cancellation is COMPUTED, and the residual stays symbolic as (g_nat−1)/g_± — collapsing to 0 only when g_nat→1 is substituted last (NOT an X−X self-check). |
| 123 | parent_normalized_branch_values | `mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl` | 6 | The `Xi_v` and `Xi_T` LAWS emerge from `Reduce`-based inversions (for `v_w0` and `T_m`); numerics computed from the derived laws at upstream-anchored `r_F1`/`g_±`. CRITICAL guard satisfied: the un-squared λ NEGATIVE sign (stage-118 convention) is preserved, so `Xi_v(F1) ≈ −1.01675633282526` (negative, not the +1.0168 mis-pass); `v_w0` declared REAL (not positive); healing lock `c_s→ℏ/(2 m ℓ)` applied only in the `Xi_T` inversion. |

**Totals:** 3 stages, 3 dirty (the missing-`.wl` gap), 0 clean, 0 status-only,
0 blocked. **3 new independent-route `.wl` (121, 122, 123).** 21 PASS across the
three new scripts.

## Mathematica mirror policy — the dual-engine retrofit

All 3 new `.wl` are **GENUINE INDEPENDENT routes** (0 transliterations
accepted). All 3 are added to the Independent-Mirror Set in
`MATHEMATICA_MIRROR_POLICY.md` (inserted in stage order between `090` and
`175`). The labor split was strictly enforced: **Claude reviews** (audit +
verify); **Codex writes ALL script code** (designs and writes the new `.wl`);
the directives stated only the requirement + acceptance criteria, never script
code.

## SymPy side

**UNCHANGED.** The SymPy `.py` reference engines for 121/122/123 were not
touched — this sweep only adds the second engine. (The IV.3 SymPy-side
de-tautologization of 122's reciprocal-traction check — derived from the
stage-119 proportionality law `g = C/T_m` plus the `g_nat = 1` ansatz — remains
as it was at IV.3 close; the new `.wl` preserves that de-tautology cross-engine
via the free-positive-symbol `cStage`.)

## Paper / notes

**UNTOUCHED.** No `.tex` derived value, constant, or identity target moved.
The verify agents did observe pre-existing pre-renumber stale-numbering in the
paper/notes prose (`stage_122.tex` header "Stage 139", `stage_123.tex` header
"Stage 140"; notes reference "Stage 221/223") — these are out of red-team
SCRIPT scope and are logged to `PAPER_CLEANUP_TRACKER.md` (appended to the
existing banner-drift / stale-numbering cluster entry P4-51) for a later
paper-side pass. No edit proposed here.

## Consult

None. No iteration-2 reworks; no math directions needed Claude+Codex
resolution; nothing escalated to the user.

## Verification

The three new-`.wl` reviews are recorded with the stage verification artifacts
(`redteam/verifications/stage_121.md` … `stage_123.md`). Final verdicts:
- `verified` (3): 121, 122, 123.
- `needs_rework` (0); `blocked_unfixable` (0).

Material change: **0** (`material_change: false` on all 3 — second-engine
addition only; the SymPy reference engines and the paper/notes are unchanged,
so no derived value, constant, identity target, or paper number moved, and no
`upstream_stale` propagation).

## Cumulative

**200/253 stages red-team verified (79.1%)** — UNCHANGED by this retro-sweep
(a retro-sweep adds second engines to already-`verified` stages, it does NOT
add new verified stages). What changed: the dual-engine gap for already-
`verified` non-status-only stages is now **CLOSED** (121/122/123 were the only
3; confirmed against the inventory). The 11 status-only single-engine stages
compute nothing a `.wl` could check and legitimately remain single-engine.

The 6 prose trackers under `notes/` were synced with a dated retro-sweep entry:
`MATHEMATICA_MIRROR_POLICY` (121/122/123 added to the Independent-Mirror Set +
gap-closed note), `STAGE_VERIFICATION_COVERAGE` (count unchanged at 200/253,
gap closed), `CHECKPOINT_TRUST_AUDIT` (no checkpoints → no change),
`CHECKPOINT_CONSTANT_PROVENANCE` (no new constant → no change),
`PAPER_CLEANUP_TRACKER` (stale-numbering observation appended to P4-51),
`EM_PROJECTED_INTEGRATION_TRACKER` (out of EM-projected range → no change).
