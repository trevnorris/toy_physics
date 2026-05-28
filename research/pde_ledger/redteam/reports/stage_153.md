---
unit_id: 153
batch: IV.6
auditor_model: claude-opus-4-7-1m
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
  notes_stage_files: [moving_throat_pde_stage153_full_mouth_correction_status.md]
  paper_appendix: present
---

# Audit unit 153 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_153.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage153_full_mouth_correction_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read intro list at L25-33 and the `\input{stages/stage_153}` insertion at L1340; the appendix has no per-stage prose row for 153 beyond the input line and the block summary "Stages 146--153: first-order mouth-profile rigidity and finite mouth-only correction.")
- sympy: (missing — no `scripts/moving_throat_pde_stage153_*.py` exists)
- mathematica: (missing — no `mathematica/moving_throat_pde_stage153_*.wl` exists)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 153 is explicitly labelled `\stagefield{Verification}{SymPy audit: none yet. Mathematica audit: none yet.}` in the .tex card. The card is a status / derivation-ledger entry inside the "Finite mouth-profile corrections" block whose `\stagefield{Purpose}` is described as "a finite mouth-profile corrections ledger step" and whose `\stagefield{Inputs}` import the regular canonical mouth branch, positive deformations of the source profile, the first-order rigidity kernel, and the full mouth potential. The single quoted bottom-line result is qualitative: "Branch ambiguity is gone; finite mouth-only correction motivates full co-evolution." The checks list three qualitative properties (centred first-order deformations, rigidity kernel controlling non-exponential corrections, one-step nonlinear correction within the reduced mouth-layer regime) but no new closed-form identity. The notes supply quantitative numerical estimates (δg ≈ -0.0648, δS ≈ -0.0388, δΠ ≈ 0.9071, δT_m ≈ 0.2717; corrected canonical point Π_corr ≈ 2.4159, T_m,corr ≈ 1.1731; one-step Π_1 ≈ 2.5391, T_m,1 ≈ 1.2104), but the notes themselves explicitly describe these as a one-step consistency check, "not a proof of full nonlinear convergence" — i.e., a status snapshot, not a theorem.

## What the script claims to verify

No SymPy or Mathematica script exists for stage 153. The audit prompt confirms `is_status_only_candidate: True`, and the paper card's Verification field corroborates this with the explicit "none yet" annotation. Stage 153 carries forward results from upstream stages (the rigidity kernel, the canonical regular branch selection, the full mouth potential, the source-profile deformations enumerated in Inputs) and packages them as a "branch ambiguity is gone" ledger entry. No new identity is asserted in the card that would require its own verification engine here.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| "Branch ambiguity is gone; finite mouth-only correction motivates full co-evolution." (qualitative status statement) | — (no script; status-only) | match (carry-forward; upstream branch-selection stages own the underlying identities) |
| Check 1: first-order profile deformations centred before covariance use | — (no script) | match (upstream rigidity-kernel stage owns the centring claim) |
| Check 2: rigidity kernel controls non-exponential corrections | — (no script) | match (upstream rigidity-kernel stage owns this) |
| Check 3: one-step nonlinear correction stays within reduced mouth-layer regime | — (no script; notes describe it qualitatively as a consistency check) | match (status carry-forward) |

The numerical estimates quoted in the notes (δg, δS, δΠ, δT_m, etc.) are not asserted as paper-card deliverables in `stage_153.tex` — they appear only in the notes as quantitative context. The notes themselves frame them as a status snapshot ("not a proof of full nonlinear convergence"). No `paper_misalignment` arises because the card does not commit those numerics as verification targets.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | — | — | (no scripts exist for this unit) | n/a | n/a |

## Findings

None.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` file exists.

## Engine cross-check

Not applicable — neither engine is present.

## Verdict justification

Verdict: `clean` with `findings_count: 0`. Stage 153 is a legitimate status-only / ledger card: the .tex Verification field explicitly self-declares "SymPy audit: none yet. Mathematica audit: none yet.", and the manifest flags `is_status_only_candidate: True`. The card's bottom-line is qualitative ("branch ambiguity is gone; finite mouth-only correction motivates full co-evolution") and the Inputs list four upstream artifacts (regular canonical mouth branch, positive source-profile deformations, first-order rigidity kernel, full mouth potential) that all live in earlier stages of the same Part IV block ("Stages 146--153: first-order mouth-profile rigidity and finite mouth-only correction"). The quantitative numbers that do appear (δg, δS, δΠ, δT_m, Π_corr, T_m,corr, Π_1, T_m,1) are confined to the notes and are framed there as a one-step consistency check rather than a theorem; the paper card does not commit them as verification targets, so no `missing_*` finding is warranted under the status-only rule of the audit prompt. The notes' reference to a "Stage 249 rigidity kernel" appears to be a numbering / labelling note inside the working notes and is not load-bearing for the paper card's claim, which only references "the first-order rigidity kernel" via Inputs. I attempted to break the alignment by (a) hunting for a hardcoded quantitative claim in the .tex that no upstream stage could plausibly own, (b) checking that the Inputs really do match the deliverables enumerated in the Checks list and the qualitative Output, and (c) looking for any new identity introduced only here; all three attacks failed. Carry-forward is legitimate.

## Self-test notes

I checked: (1) the status-only carve-out — confirmed Verification field says "none yet" and notes frame numerics as a consistency-check status, not a theorem; (2) paper-vs-notes consistency — the .tex card's qualitative claim and the notes' qualitative interpretation agree; (3) Inputs vs Checks vs Output — all three are coherent and reference upstream artifacts (rigidity kernel, mouth potential, source deformations, canonical branch) rather than introducing new identities. Conclusion: no finding warranted; status-only carry-forward is legitimate.
