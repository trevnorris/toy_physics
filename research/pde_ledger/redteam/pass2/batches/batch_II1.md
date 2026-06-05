# Batch II.1 (stages 024-036) — pde_ledger SECOND PASS

Date: 2026-06-04

Part II.1 — Overlap isotropy through continuum kernel.

## Method

- v2 paper-grounded auditor PLUS the **exhaustive script→doc value-reconciliation
  augmentation** (`redteam/pass2/RECONCILIATION_AUGMENTATION.md`) — every audit
  agent read both its rendered prompt and the augmentation doc.
- 13 clean per-stage audit agents in parallel (audit agents execute nothing → 0
  Mathematica seats).
- Independent exec reliability gate — orchestrator re-ran BOTH engines for every
  findings stage (sequential, ≤1 seat) and refreshed the committed transcripts.
  The exec re-run + an arbiter grep are the ground truth (the audit agents
  under-/over-called the residual-label scope; the grep resolved it).
- Codex as the sole fix-applier (incl. the 024 notes edit). 3 Codex waves of 2
  (≤2 seats; exec deferred until builds finished).
- Clean per-stage verify agent for each of the 6 findings stages.

## Result

All 13 stages reached `verified` at v2 depth + value-reconciliation augmentation.
**13/13 verified, `material_change=false` on all 13. No stop-cold, no blocked.**

| Outcome | Stages |
| --- | --- |
| Clean (7) | 025, 026, 027, 028, 029, 030, 031 |
| Findings → resolved → verified (6) | 024, 032, 033, 034, 035, 036 |

Checkpoints **024** and **036** both re-verified (024 with a resolved notes typo).

## Resolutions

- **024** — TWO findings, both resolved:
  - **F1 `paper_misalignment` (value_mismatch) — RESOLVED via Claude+Codex
    (non-conceptual, notes-only, published card UNAFFECTED).** notes:213 stated the
    unit-sphere sixth-moment prefactor as `4π/122`; the textbook value is `4π/105`
    (`1/(d(d+2)(d+4)) = 1/(3·5·7)` in d=3), which BOTH engines use (py:128 `4*pi/105`;
    the `.wl` integrates the moment directly) and which the notes' OWN downstream
    `κ_* = √5/(7√π)` requires (a `122` prefactor would rescale κ_* by 105/122). The
    published `.tex` card is silent on the prefactor → unaffected. Codex corrected
    notes:213 `122`→`105` (the single hunk); no script change. `material_change:false`.
  - **F2 `stale_output`** — committed transcripts predated the scripts; refreshed
    (the Mathematica banner was stale `STAGE 007` → now `STAGE 024`).
- **032/033/034/035/036** — each a single `stale_output` (NOT paper_misalignment).
  The committed transcripts predated the June-3 numbering commit (`e2a4780`), which
  fixed the scripts' main `banner(...)` headers to canonical `STAGE 3N.x` but left
  residual stale **self-labels** (the SymPy module docstring + the closing
  `All Stage NN checks passed.` print, which prints INTO the transcript; 032 also a
  `.wl` closing Print; 033 a `.wl` sub-stage comment). These are the exact class as
  the I.2 stage-021 `.wl:195` fix. **Label-only** corrections `Stage {15,16,17,18,19}`
  → `Stage {32,33,34,35,36}` (matching each file's own banner format) + both-engine
  re-run. ZERO math/equation/assertion/value change; `material_change=false`.
  (Audit agents for 032/033/034 under-called these "no directive needed"; the arbiter
  grep on the source caught the residual labels and the directives were authored to
  cover them. 036's directive was incomplete — augmented to add the missed
  print-summary fix.)

## Value reconciliation (pass-2 augmentation)

Applied on all 13 stages; **177 deliverable values checked batch-wide, 1 misaligned**
— the single misalignment was 024's `4π/122` notes typo (resolved above). Per stage:
024=16(1), 025=14, 026=12, 027=19, 028=17, 029=12, 030=18, 031=11, 032=14, 033=16,
034=8, 035=9, 036=11. No MISSING-DELIVERABLE anywhere.

## Dual-engine / mirror status

All 13 already carried both independent engines from pass 1; the re-audit confirmed
genuine independence (e.g. 024's `.wl` integrates the sphere moments directly via
`Integrate[n[i]…Sin θ]` vs SymPy's combinatorial pairing-sum; `LinearSolve`/quotient-rule
vs `series` routes). **No new `.wl`, no mirror reclassification, 0 sanctioned mirrors.**

## ⚠️ DISCOVERY — numbering drift persists in the II.1 script/output band (NOT exhausted)

The pass-2 re-audit surfaced that the "mechanically exhausted (001-253)" numbering
reconciliation did **not** clear the script/output band for II.1. Beyond the
unambiguous self-labels resolved above, the following stale labels REMAIN and were
deliberately **left untouched** (they are content-dependent and partly ambiguous —
fixing them is the careful, gated, per-reference numbering work, NEVER an offset-sweep,
NOT the red-team loop; per memory `numbering-drift-root-cause`):

- **Clean stages' committed OUTPUTS are themselves stale** (predate `e2a4780`), with
  stale self-banners the audit agents did not flag (math verified clean regardless):
  028 out `STAGE 011 — LOADED PROFILE SELECTION`; 030 out `STAGE 13 AUDIT COMPLETE`;
  031 out `STAGE 14 AUDIT COMPLETE`; 025 out `Stage 8 Mathematica audit passed.`;
  026/027 out cross-ref `STAGE-8`/`STAGE-8/9` + self `Stage 9`/`Stage 10` passed lines.
- **Ambiguous self-vs-cross multi-epoch refs in findings-stage sources** (left as-is):
  024 `.wl:293 FINAL STAGE-007 LEDGER` + `.py "Stage-6 … bundle/stack"` (024 carries
  self-labels from BOTH the −17 epoch [→007] and another [→6], and "Stage-6 grouped
  bundle" could instead be a cross-ref to stage 023); cross-refs 033 `.py:50 Stage-15`
  (=stage 032), 036 `.py:106/127 Stage-18` (=stage 035). Canonical 3-digit cross-refs
  to `Stage 022`/`Stage 023` are already correct and untouched.

**Recommendation (user decision at the gate):** a dedicated careful per-reference
numbering pass over the script/output band (refresh all stale committed outputs +
content-map each remaining stale self/cross label), separate from the red-team loop.
The red-team's own job — verifying the MATH — is complete for II.1 (13/13).

## Status

13/13 verified, all `material_change=false`; no stop-cold/blocked. Paper/notes edits:
the single 024 notes typo (Codex-applied, Claude-reviewed) — see PAPER_CLEANUP **P5-04**.
Pass-1 `MANIFEST.yaml` untouched (isolation held). HALT at the II.1 boundary for the
user gate (incl. the numbering-drift decision above) before III.1.
