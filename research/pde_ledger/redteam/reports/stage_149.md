---
unit_id: 149
batch: IV.5
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage149_mouth_rigidity_status.md]
  paper_appendix: present
---

# Audit unit 149 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_149.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage149_mouth_rigidity_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_149}` line on L1332 references this stage; no narrative row)
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 149 is a status-summary ledger entry in the "Finite mouth-profile corrections" block of Part IV. The card explicitly carries `\claimstatus{\StatusExactClosure{} / \StatusOpen{}}` and states the headline summary in a quotation: "Mouth corrections are controlled by a finite rigidity kernel, not by branch ambiguity." The card's `\stagefield{Verification}` is "SymPy audit: none yet. Mathematica audit: none yet." — the paper card itself does NOT assert any new computed identity for this stage. The card lists three `\stagefield{Checks}` items (deformation centering, rigidity-kernel control, one-step nonlinear correction inside the reduced regime). The notes file expands the body into five enumerated bullets summarizing the first-order moment representation, the bias retuning formula $\delta\Pi$ and traction shift $\delta\widehat T_m$, and quotes upstream-derived constants $A_T\approx-4.27263956256927$, $B_T\approx0.134875005736706$, with $|A_T|/B_T\approx 31.68$. The notes are explicit that "the stage remains a perturbative rigidity summary, not yet a full nonlinear mouth theorem for arbitrary finite deformations" and reframe Stages 197–199 / 228 / 247 results as a status update — i.e., a carry-forward, not a new derivation.

## What the script claims to verify

No script exists for this unit. There are no assertions to evaluate. The unit's manifest entry has `is_status_only_candidate: True`, the paper card explicitly states "none yet" for both engines, and the notes frame the stage as a summary of earlier work rather than a new theorem. The status-only carve-out described in the prompt therefore applies.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Status summary: mouth corrections controlled by finite rigidity kernel, not branch ambiguity | (no script) | n/a — status-only |
| Carry-forward of first-order moment representation $\Sigma_\epsilon = (1-\epsilon)\Sigma_* + \epsilon\varsigma$ | (no script; upstream stages 197–199 derive the first-order machinery) | n/a — carry-forward |
| Carry-forward of $\delta\Pi$, $\delta\widehat T_m$ formulas and constants $A_T$, $B_T$ | (no script; constants originate in upstream stages, per notes) | n/a — carry-forward |
| Checks list (deformation centered, rigidity kernel controls, one-step correction inside regime) | (no script; these are status assertions about prior derivations) | n/a — status-only |

The paper card and notes do not introduce any new numeric identity that demands fresh in-stage verification. The card's `\stagefield{Verification}` line is itself "none yet," so the paper makes no claim of script verification. The unit is consistent with status-only treatment.

`paper_alignment: aligned` — the script side (empty) does not contradict the paper side (no script claimed). The paper makes no assertion that the absent scripts can fail to support.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | (none) | — | — | — | — |

No assertions to inventory.

## Findings

None.

## Independent-derivation check (Mathematica)

No `.wl` script exists. Not applicable.

## Engine cross-check

No engines present. Not applicable.

## Verdict justification

Stage 149 is a status-only ledger entry whose paper card explicitly declares "SymPy audit: none yet. Mathematica audit: none yet." The notes corroborate that the unit is a perturbative summary of upstream results (Stages 197–199, 228, 247) and not a fresh derivation. The numeric constants $A_T$, $B_T$ named in the notes are sourced upstream, not freshly derived in this stage; the notes label them with "approximately" and frame them as carry-forward values. The `is_status_only_candidate: True` manifest flag, the status carve-out in the prompt, and the explicit paper-side "none yet" together mean that missing SymPy/Mathematica scripts are not findings here. No `paper_misalignment` exists because the paper does not claim a verified identity the script fails to honor (the script does not exist; the paper does not claim one does). Attacks attempted: (1) Does the paper card or notes assert a fresh, in-stage identity not present upstream? No — the boxed expression "what small non-exponential correction does the real moving-throat mouth layer induce around $(\Pi_*, \widehat T_{m,*})$?" is framed as the *open question*, not a proved result, and the carry-forward constants are labeled as upstream-derived. (2) Are the listed `\stagefield{Checks}` (deformation centering, rigidity-kernel control, one-step correction within regime) phrased as new theorems demanding in-stage verification? They are phrased as status conditions that reference the upstream first-order machinery; the card itself flags them as ledger items, not new derivations. (3) Does the notes file disagree with the paper card on any numeric or symbolic? No — the card is terse and consistent with the notes' more detailed bullets. Verdict: `clean`.

## Self-test notes

Checked: (1) status-only carve-out applicability — confirmed by manifest flag, paper card's explicit "none yet" verification field, and notes' "perturbative rigidity summary" framing; (2) whether constants $A_T \approx -4.27263956256927$, $B_T \approx 0.134875005736706$ in notes are presented as new derivations — no, they are quoted with "approximately" and the notes mention Stages 197–199/228/247 as their context; (3) paper round-trip — no script exists, so no risk of script-side claim outrunning paper. No directive written (zero findings).
