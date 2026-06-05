---
unit_id: 069
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T14:10:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 069

Checkpoint stage — higher bar applied. The sole finding (F1, `stale_output`, low) is a
REFRESH-ONLY fix: no source edit, the two committed `.txt` transcripts re-generated so the
banner self-label reads the canonical `STAGE 069`.

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
No source edit (the directive explicitly required none, and the captured
`stage_069_diff.patch` is 0 bytes). The orchestrator's independent re-run overwrote both
committed transcripts. `git diff HEAD` on each `.txt` shows a single-line change, line 3:
- `scripts/output/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.txt:3` —
  `STAGE 052 …` → `STAGE 069 — FINAL REDUCED SUPPORT/SOURCE VERDICT`
- `mathematica/output/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.txt:3` —
  same `STAGE 052` → `STAGE 069`

Both diffs are `2 +- / 2 -` total (one line each); nothing else in either transcript changed.

**Assessment:**
Correct and complete. The two source scripts already emit the canonical banner
(`scripts/...sympy_audit.py:65` `banner("STAGE 069 …")`;
`mathematica/...mathematica_audit.wl:33` `banner["STAGE 069 …"]`), and `git diff --stat HEAD`
confirms neither source `.py`/`.wl` is modified — exactly the expected refresh-only shape.
The refreshed committed line 3 of both `.txt` files reads `STAGE 069 — FINAL REDUCED
SUPPORT/SOURCE VERDICT`. The `Stage-49`/`Stage-51` lines surviving in the FINAL LEDGER of
both transcripts are cross-references (matched-branch → upstream, resonance-family →
upstream), identical between the committed `.txt` and the fresh exec logs; per orchestrator
context these are correctly deferred to the dedicated numbering plan and are not a defect.
No collateral edit. No source change ⇒ no downstream impact.

## Exec log assessment

**SymPy:** exit=0. Banner now `STAGE 069 — FINAL REDUCED SUPPORT/SOURCE VERDICT`. Substantive
checkpoint checks pass symbolically:
- `matched fail edge from W_match(Delta_inf) = 0`, `matched success edge from W_match(Delta_0) = 0`
- `W_match decreasing in Delta_eff > 0 -> 1`; `matched window width = 0`
- `failure-side width - Pe_req(1-C_res^2)/(C_res^2 Delta_inf) = 0`, `success-side … = 0`
- `P_res - 1 - (1-C_res^2)/C_res^2 = 0`; interior-point ordering all strictly positive
  (`failure-band point - matched fail edge > 0`, `resonance fail edge - failure-band point > 0`, …)

**Mathematica:** exit=0. Banner now `STAGE 069`. Independent route (Cres2Prim primitive,
`Pres = 1/Cres2`, Solve-derived gap):
- `PASS: matched window width`, `PASS: failure-side width - Pe_req(1-C_res^2)/(C_res^2 Delta_inf)`,
  `PASS: success-side width - …`, `PASS: P_res - 1 - (1-C_res^2)/C_res^2`
- `PASS: Pres-PresGap consistency via Solve`; all five interior-point positivity checks `PASS`.

No numeric constant is pinned/re-asserted: `grep` for `0.994418…`/`1.00561…`/`0.56` in both
source files returns NONE; `C_res^2`/`P_res` are carried as the free symbols
`Pres_gap`/`Cres2` (SymPy) and `Cres2Prim` (Mathematica). Checkpoint bar met — load-bearing
identities (window width, side-band widths, `P_res-1=(1-C²)/C²`, interior ordering) are real
algebraic zeros over free positive parameters, not re-assertions.

**Output freshness:** confirmed. Exec logs dated `2026-06-05T13:5x`; committed `.txt` line 3
now matches the fresh banner; `git diff HEAD` on both `.txt` is exactly the `STAGE 052→069`
banner line and nothing else, with the assertion body identical to the exec logs (22 SymPy
result/PASS lines, 23 Mathematica PASS lines).

## Material-change assessment

`material_change`: false. The only change is the stale banner self-label in two committed
transcripts (`STAGE 052` → `STAGE 069`). No source edit, no derived result altered, no symbol
or constant introduced. No downstream unit can depend on a transcript banner string.

## Side observations (non-blocking)

- The working tree shows sibling-stage edits (061–072) from the in-progress batch III.3; only
  the two stage-069 `.txt` files and unchanged stage-069 sources are in scope here, and those
  are exactly as expected.
- `stage_069_diff.patch` is 0 bytes — correct for a refresh-only stage (the patch captures
  source edits, of which there are none); the `.txt` refresh is visible via `git diff HEAD`.

## Verdict justification

verdict = verified. The single low-severity finding (stale `STAGE 052` banner) is resolved by
the orchestrator's re-run: both committed transcripts now carry the canonical `STAGE 069`
banner with the assertion body unchanged, the source `.py`/`.wl` are untouched vs HEAD (empty
diff patch; clean `git diff --stat HEAD`), both engines exit 0, every substantive checkpoint
identity passes, and no numeric `C_res²`/`P_res` is pinned (kept symbolic). No regression and
no material change.
