---
unit_id: 132
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 132

## Per-finding outcomes

### F1 — paper_misalignment

**Classification:** resolved

**What changed:**

The user-resolved Cluster B (orchestrator-direct, Codex bypassed) edited the notes file
`/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md`
at two locations:

- Line 2 (H1): `# Moving-Throat PDE — Stage 234: Mouth Boundary-Layer Status After Explicit Source-Law Extraction` → `# Moving-Throat PDE — Stage 132: Mouth Boundary-Layer Status After Explicit Source-Law Extraction` (Cluster A's mass H1 renumber, +102 offset removed for 132).
- Line 6 (carry-forward attribution): `After Stages 180–182, the mouth-source side is no longer an abstract profile problem.` → `After Stages 129–131, the mouth-source side is no longer an abstract profile problem.`

Git diff against HEAD confirms both edits are exactly these two text changes; no other content in the notes file was altered. The paper card `paper/stages/stage_132.tex` is untouched (correct — its own `Verification: SymPy audit: none yet. Mathematica audit: none yet.` was already truthful).

**Assessment:**

Both edits fully address F1.

- The H1 now matches the filename slug `stage132`, the paper section header `Stage 132`, and the appendix `\input{stages/stage_132}` row. The stage-234 misanchor is gone.
- The carry-forward attribution now references stages numerically *less than* 132. Per the resolution doc, the actual upstream content lives at:
  - Stage 130 — closed form of the Family-1 mouth-bias factor `g_Π = 2Π(2Π e^Π + π) / [(4Π² + π²)(e^Π − 1)]`.
  - Stage 131 — parent micro-threshold root yielding `Π_* ≈ 1.50882951349316`.
  - Stage 129 — mouth boundary-layer setup (`σ_Π(z) = Π e^{−Πz/L} / [L(1 − e^{−Π})]`).
  The chosen citation "Stages 129–131" is a tighter range than the resolution-doc draft (which said "Stages 130–131"); 129 is included because the explicit source family `σ_Π` itself originates there, so the attribution is correctly inclusive of all three load-bearing upstream derivations for the items boxed in stage 132's notes. The chain now points strictly backward in the ledger and supports the boxed prose claim of stage 132.

No script edits were required (status-only unit; the audit explicitly stated this is a `paper_misalignment` requiring user resolution, with no script side to verify). No collateral edits made beyond the two prescribed lines.

## Exec log assessment

**SymPy:** exit=n/a. No script exists under the `stage132` slug (status-only unit); the paper card declares `SymPy audit: none yet.` truthfully. Audit log path `scripts/output/moving_throat_pde_stage132_mouth_boundary_layer_status_sympy_audit.txt` is absent by design.

**Mathematica:** exit=n/a. Same — no `.wl` under the `stage132` slug; paper card declares `Mathematica audit: none yet.` truthfully.

**Output freshness:** N/A — no scripts, hence no outputs.

## Material-change assessment

`material_change`: false.

The fix is purely textual within a status-only notes file. No numeric constant changed; the cited carry-forward stages (129–131) are the same stages already producing the numeric values (`σ_Π`, `g_Π`, `Π_* ≈ 1.50882951349316`) inside batch IV.4. No downstream unit's numeric anchor moved. The orchestrator can leave `upstream_stale` as-is for units > 132 with respect to this fix; the value `Π_* ≈ 1.50882951349316` quoted at line 30 of the notes is preserved exactly.

## Side observations (non-blocking)

- The notes file uses an en-dash in "Stages 129–131" (U+2013), matching the original "Stages 180–182" punctuation; consistent with batch-wide notes style. Non-blocking.
- The paper card's `\stagefield{Downstream use}` still cites Stages 133–145 as consumers; those downstream consumers are unchanged in numeric terms (Cluster B did not touch their carry-forward anchors), so no follow-up to those units is implied by this fix.

## Verdict justification

The single finding F1 was a `paper_misalignment` flagged by the auditor because the notes H1 and carry-forward chain pointed to numerically-later stages (234 / 180–182), breaking the upstream ordering for a status-only unit. The user-resolved Cluster B made the exact two-line edit prescribed: H1 corrected to "Stage 132", and the carry-forward attribution corrected to "Stages 129–131" — strictly less than 132 and matching the actual IV.4 upstream stages where `σ_Π`, `g_Π`, and `Π_*` are derived. Both edits are visible in the working tree against HEAD, with no collateral changes. The verdict is `verified`.
