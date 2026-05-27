---
unit_id: 120
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00-06:00
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage120_core_parameter_status.md
  paper_appendix: present
---

# Audit unit 120 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_120.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage120_core_parameter_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_120}` row at line 1274 and the `MTDC-T8` anchor block at lines 1175-1179 are relevant)
- sympy: (missing — intentionally; manifest flags `is_status_only_candidate: true`)
- mathematica: (missing — intentionally; manifest flags `is_status_only_candidate: true`)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

The stage 120 card is a status-only ledger entry under Part IV, anchor `MTDC-T8`. The `\claimstatus{}` line marks it as `\StatusExactClosure{} / \StatusOpen{}`, and the `\stagefield{Verification}` line says verbatim: "SymPy audit: none yet. Mathematica audit: none yet." The body claim is the single line in the `\begin{quote}` block: "Records the surviving finite parent-level target after overlap extraction." The `Downstream use` field is explicit that the card "is a derivation ledger entry, not an unconditional actual-branch theorem" and that the open-status tag must be carried forward when cited. The `Inputs` enumerates four upstream items: compensated outlet conditions, the two-channel shell/mixed Schur complement, the finite D/N mixed-tube geometry, and parent-action overlap data — all of which are derived/verified upstream in stages 114-119. The accompanying notes file goes further and states a sharper *structural* gate ("compensated canonical outgoing quadrupole branch exists iff `1 + r^2 = 4 (g - r)^2`, with `L_W = (pi a / 2) sqrt((1+r^2)/3)`"), but this is presented in the notes as a conditional realization condition, not a verified actual-branch theorem.

## What the script claims to verify

No SymPy or Mathematica scripts exist for this unit; the manifest's `is_status_only_candidate: true` flag explicitly permits this, and the paper card itself acknowledges no script-side verification ("SymPy audit: none yet. Mathematica audit: none yet."). The unit is a paper-side ledger entry that consolidates upstream verified results into a status statement; the carry-forward content (the gate equation `1 + r^2 = 4 (g - r)^2` from the notes) is actually exercised by neighboring scripts:
- `scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:20` defines `LW_formula = sp.pi*a/2 * sp.sqrt((1+r**2)/3)` — the L_W formula from the notes.
- `scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:55-56` calls `expect_zero("lower branch balance relation", 1 + r**2 - 4 * (gminus - r) ** 2)` and the analogous upper-branch check — directly the gate equation from the notes.

So the structural content the notes hint at is independently verified by upstream/downstream scripts; the status-only carve-out is legitimate.

## Paper ↔ script cross-check

| Paper deliverable | Coverage | Notes |
|---|---|---|
| "Records the surviving finite parent-level target after overlap extraction" | match (status-only) | Card explicitly status-only; no script required for this unit. |
| Inputs: compensated outlet conditions | match (upstream) | Carried from earlier outlet-compensation stages in the same part. |
| Inputs: two-channel shell/mixed Schur complement | match (upstream) | Carried from stage 114 (`*_concrete_core_schur_*`). |
| Inputs: finite D/N mixed-tube geometry | match (upstream) | Carried from stage 116 (`*_dn_mixed_tube_realization_*`). |
| Inputs: parent-action overlap data | match (upstream) | Carried from stage 119 (`*_parent_balance_*`). |
| Notes-side sharp gate `1+r^2 = 4(g-r)^2` | match (verified downstream) | Verified in stage 125 sympy at lines 55-56. |
| Notes-side `L_W = (pi a/2) sqrt((1+r^2)/3)` | match (verified downstream) | Verified in stage 121 sympy at line 20. |
| Verification status: "SymPy audit: none yet. Mathematica audit: none yet." | match (truthful) | No script present; card does not overclaim. |
| Status tag `\StatusOpen{}` | match | Card is honest about the conditional / status-only nature. |

`paper_alignment` set to `aligned`: every load-bearing paper-side statement either (a) is an honest status disclaimer, (b) is an input citation traceable to an upstream script-verified stage, or (c) is a structural content statement verified in a neighboring script.

## Assertion inventory

No SymPy or Mathematica scripts exist for this unit, so there are no assertions to inventory. The status-only carve-out applies. The empty inventory is consistent with both the manifest flag (`is_status_only_candidate: true`) and the paper card's own self-disclosure (`SymPy audit: none yet. Mathematica audit: none yet.`).

## Findings

None.

## Independent-derivation check (Mathematica)

N/A — no Mathematica script exists. Status-only carve-out applies; no upstream Mathematica result is asserted to be carried forward in this card (the Inputs list cites upstream stages, not constants pinned by this card alone).

## Engine cross-check

N/A — neither engine has a script for this unit.

## Verdict justification

The card is a faithful status-only ledger entry. Attacks attempted:
1. **Does the card overclaim?** No — it explicitly marks `\StatusOpen{}` and says "SymPy audit: none yet. Mathematica audit: none yet." The body line "Records the surviving finite parent-level target after overlap extraction" is descriptive, not assertive of a new theorem.
2. **Does the card contradict its notes?** No — the notes state a sharper structural gate (`1+r^2 = 4(g-r)^2`, `L_W = (pi a / 2) sqrt((1+r^2)/3)`), and the card's body is a coarser summary, but the notes themselves frame the gate as a conditional realization condition, not as a verified actual-branch theorem. The card's `\StatusOpen{}` is consistent with that conditional framing.
3. **Is the carry-forward chain legitimate?** Yes — the structural identities the notes invoke are independently exercised by `scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:20` (L_W formula) and `scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:55-56` (gate equation, both branches). The four Inputs listed are each carried from numbered upstream stages within the same part.
4. **Does the card assert a new derivation that no script supports?** No — every quantitative claim in the card or its notes is either flagged as conditional/open or is verified in a neighboring script.
5. **Numbering drift.** The card's section heading is `\section[Stage~137]{Stage~137: Core-Parameter Extraction Status}` while its label is `\label{stage:120}`, and the notes header reads "Stage 222: Core-Parameter Extraction Status". This is a pre-reorder artifact consistent across all neighboring cards in Part IV (e.g., stage 115 is "Stage 134", stage 121 is "Stage 138") and is not a content finding for this audit — it does not change what is being claimed.

Verdict: `clean`. No findings; no Codex directive required.

## Self-test notes

I checked: (1) whether the card asserts a new derivation — it does not, it is an open-status ledger entry whose Verification field is explicit about "none yet"; (2) whether the notes' sharp gate `1+r^2 = 4(g-r)^2` and `L_W` formula appear without script support — they do appear in upstream/downstream sympy scripts (stage121 line 20, stage125 lines 55-56), so the status-only carve-out is supported; (3) whether the Inputs list cites items that are verified upstream — yes, all four trace to numbered stages in the same part. Status-only carve-out applies cleanly.
