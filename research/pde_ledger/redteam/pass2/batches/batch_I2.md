# Batch I.2 (stages 013-023) — pde_ledger SECOND PASS

Date: 2026-06-04

Part I.2 — Maxwell bridge, parent throat action, reduced one-port.

## Method

- v2 paper-grounded auditor PLUS the **exhaustive script→doc
  value-reconciliation augmentation** (see
  `redteam/pass2/RECONCILIATION_AUGMENTATION.md`) — every audit agent read both
  its rendered `audit_prompt_NNN.md` and the augmentation doc.
- One clean per-stage audit agent per stage (no context contamination); 11
  agents run in parallel (audit agents execute nothing → 0 Mathematica seats).
- Independent exec reliability gate (orchestrator re-runs the affected engine;
  refreshes the committed output). The exec re-run is the ARBITER.
- Codex as the sole fix-applier (orchestrator does not hand-apply directive
  fixes).
- Clean per-stage verify agent for the one finding stage.

## Result

All 11 stages reached `verified` at v2 depth + value-reconciliation
augmentation. **11/11 verified, `material_change=false` on all 11. No
stop-cold, no blocked.**

| Outcome | Stages |
| --- | --- |
| Clean (10) | 013, 014, 015, 016, 017, 018, 019, 020, 022, 023 |
| Findings → resolved → verified (1) | 021 |

## The one resolution

- **021** — single low-severity `stale_output`. The committed Mathematica
  output `.txt` was a pre-numbering-reconciliation capture (banner read
  `STAGE 004 — MAXWELL + MIXED-SECTOR REDUCTION`), and the live `.wl:195`
  retained a leftover stale label `Print["Stage 004 Mathematica audit
  passed."]`. The numbering-reconciliation Phase-1 scan keyed on the `banner[...]`
  call (fixed at `.wl:35` → `STAGE 021`) and on docstrings, so this bare `Print`
  literal slipped through. **Fix: one label-only `.wl:195` edit `Stage 004` →
  `Stage 021` + output refresh.** ZERO math / equation / assertion change; no
  result value changed; no paper or notes edit; **NOT a `paper_misalignment`**.
  Codex-applied (iter 1, exit 0), orchestrator-re-run (exit 0, banner now
  `STAGE 021`, all `PASS:`/`= 0` residuals unchanged), Claude-verified.
  `material_change: false` → no downstream invalidation.

## Value reconciliation (pass-2 augmentation)

Applied on all 11 stages; **every emitted deliverable value reconciles —
0 misaligned** across the batch. Values checked per stage: 013=6, 014=5, 015=4,
016=9, 017=6, 018=2, 019=5, 020=11, 021=11, 022=14, 023=13 (86 deliverable
values total). No MISMATCH and no MISSING-DELIVERABLE anywhere in I.2, so the
augmentation surfaced no `paper_misalignment` this batch (consistent with the
guard against false-positive flooding on terse cards).

## Dual-engine / mirror status

All 11 stages already carried both an independent SymPy `.py` and Mathematica
`.wl` from the first pass; the re-audit confirmed genuine independence (e.g.
021's `.wl` uses `LinearSolve` / `SphericalHankelH1` / `VariationalMethods`'
`EulerEquations` against the SymPy hand-built closed forms / `euler_equations`).
**No new `.wl` built, no mirror reclassification, 0 sanctioned mirrors
introduced.**

## Checkpoints

022 and 023 (the two checkpoint stages in range) re-verified clean with the
value-reconciliation augmentation; no checkpoint constant changed or was
introduced, no trust impact.

## Status

11/11 verified, all `material_change=false`; no stop-cold, no blocked. **No
paper/notes edits this batch** — the single finding (021) was a script-internal
stale label + output refresh, not a prose item. See
`notes/PAPER_CLEANUP_TRACKER.md` (P5-03). HALT at the I.2 boundary for the user
gate before I.3.

## Key lesson (reinforced)

The numbering reconciliation was declared "mechanically EXHAUSTED (001-253)",
but its deterministic scan was keyed on `banner[...]` calls / docstrings /
self-titles — a bare `Print["Stage NNN ..."]` summary literal (021 `.wl:195`)
was outside that key and survived. The pass-2 value/label re-audit is what
caught it. Cost to fix is trivial (one label-only line), but it confirms the
second pass earns its keep on residual label drift the dedicated numbering pass
could not have keyed on.
