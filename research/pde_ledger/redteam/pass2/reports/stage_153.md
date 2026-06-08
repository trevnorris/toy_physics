---
unit_id: 153
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00-06:00
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage153_full_mouth_correction_status.md]
  paper_appendix: present
---

# Audit unit 153 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_153.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage153_full_mouth_correction_status.md` (only stage153 notes file present)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows at lines 19, 31, 1175-1179, 1340 reference this unit / its anchor MTDC-T8)
- sympy: (missing — by design; `is_status_only_candidate: true`, manifest path `null`/`exists: false`)
- mathematica: (missing — by design; manifest path `null`/`exists: false`)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 153 is explicitly a **status / derivation-ledger card** ("v58 Status"), not a fresh computation. The card's `\stagefield{Verification}` states verbatim: "SymPy audit: none yet. Mathematica audit: none yet." Its `\stagefield{Purpose}` is "a finite mouth-profile corrections ledger step" whose "audit target is the verification output quoted below," and the quoted bottom line is: "Branch ambiguity is gone; finite mouth-only correction motivates full co-evolution." The card carries no numeric results itself; the supporting numbers live in the companion notes file, which records: a definite first-order non-exponential correction around the unique regular lower compensated branch; the residual `R_*(x)=Phi_*(x)-Pi_* x` with `R_*(0)=R_*'(0)=0`, `R_*''(0)=-3 Sigma_m^* Pi_*/(1-e^{-Pi_*})<0` (tangent-matched, sublinear → source broadening); the actual first-order Family-1 moment shifts `delta g_act ≈ -0.0648069688`, `delta S_act ≈ -0.0388718369`, retuning `delta Pi_act ≈ 0.9070844148`, `delta Tm_act ≈ 0.2716539795`; the corrected canonical point `Pi_corr ≈ 2.4159139283`, `Tm_corr ≈ 1.1731380336`; and a one-step nonlinear Picard check giving `g_1 ≈ 0.6844235741`, `S_1 ≈ 0.6163331306`, retuned point `Pi_1 ≈ 2.5391484761`, `Tm_1 ≈ 1.2103694208`. All of these are presented as carried-forward consolidation of the Stage 152 result, projected against the Stage 147 rigidity kernel. The downstream-use line confirms it feeds Stages 154-163 and carries a status tag, "not an unconditional actual-branch theorem."

## What the script claims to verify

No script exists for this unit, by design. There are no assertions to evaluate. The audit therefore reduces to: (a) confirm the unit is genuinely a status/consolidation card whose results are produced and verified upstream, not a stage that should own executable scripts; and (b) reconcile every value the card/notes state against the upstream stage(s) they cite.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| Bottom-line status: "branch ambiguity gone; finite mouth-only correction motivates co-evolution" | none (status card) | n/a — consolidation, no own computation |
| Residual `R_*` tangent-match + negative curvature (broadening) | none in-stage | carried from Stage 152 / Stage 151 derivation |
| `delta g_act`, `delta S_act`, `delta Pi_act`, `delta Tm_act`, `Pi_corr`, `Tm_corr` | none in-stage | carried verbatim from Stage 152 (boxed there) |
| one-step nonlinear `g_1, S_1, Pi_1, Tm_1` | none in-stage | carried verbatim from Stage 152 §3 |

`paper_alignment` = aligned: every stated value traces to an upstream stage's stated result, and the status tag the card carries matches the consolidation role.

## Assertion inventory

No assertions — no SymPy or Mathematica script exists for this unit (status-only by design). Table omitted (empty).

## Findings

None.

The single status-only carve-out test from the audit prompt — "is a `missing_sympy`/`missing_mathematica` finding valid?" — fails to trigger here: every numeric result the card/notes report is a verbatim carry-forward of values that ARE derived and (per upstream cards) verified at Stage 152 (the "Actual Family-1 Mouth Correction and One-Step Nonlinear Check" stage) using the Stage 147 rigidity coefficients `A_T`, `B_T`. The values that would need engine support are owned upstream, so no engine is required at 153. The card itself explicitly disclaims any in-stage audit ("SymPy audit: none yet. Mathematica audit: none yet."), consistent with `is_status_only_candidate: true` and `is_checkpoint: false` in the manifest. This matches the pattern of the other notes-only/status stages (103/113/120/124).

## Independent-derivation check (Mathematica)

N/A — no `.wl` script exists (status-only by design).

## Engine cross-check

N/A — neither engine present.

## Verdict justification

Clean. I read the paper card, the single stage153 notes file, and the part-04 appendix rows, then attacked the "should this have scripts?" question by tracing every reported number to its source. The card is unambiguously a status/derivation-ledger entry (it self-declares "v58 Status" and "SymPy audit: none yet / Mathematica audit: none yet"), and every value in the companion notes is a verbatim carry-forward from Stage 152, which itself rests on the Stage 147 rigidity kernel — both of which are present and own the actual derivations. I confirmed no stale `168π²`/`168%`/`100%` artifact appears anywhere in the stage153 docs. I verified the internal arithmetic the card/notes assert: `Pi_corr = Pi_* + delta_Pi_act = 1.50882951349316 + 0.907084414842908 = 2.41591392833607`, matching the stated `Pi_corr ≈ 2.4159139283` (and Stage 152's boxed value); the truncated values in the 153 notes (`-0.0648069688`, `-0.0388718369`, `0.9070844148`, `0.2716539795`, `1.1731380336`, etc.) are consistent leading-digit truncations of the full-precision Stage 152 boxed values. No paper_misalignment, no missing-script finding, no internal inconsistency. The unit consolidates correctly.

## Self-test notes

Checked: (1) carry-forward provenance — every 153 value traces to a Stage 152 boxed result (which cites Stage 147 `A_T`, `B_T`), so no orphaned/unsupported constant; (2) internal arithmetic — `Pi_corr = Pi_* + delta_Pi_act` reproduces the stated 2.41591392833607 exactly; (3) stale-constant scan — no `168`/`168π²`/`100%` artifact present in either stage153 doc; (4) status-only legitimacy — manifest `is_status_only_candidate: true`, `is_checkpoint: false`, card self-declares no audit, so missing engines are not findings. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

No scripts emit values for this unit (status-only, no `.py`/`.wl`/output `.txt`). The reconciliation therefore checks every numeric deliverable the stage153 notes/card state against the upstream stage that derives it (Stage 152, drawing on Stage 147). Base of reconciliation is the doc text plus the upstream notes, as the augmentation permits when outputs are absent.

| value | source (stage153 notes line) | upstream/.tex/.md location | status |
|---|---|---|---|
| `delta g_act ≈ -0.0648069688` | notes:31 | stage152 notes:37 (`-0.0648069687666328`) | MATCH (truncation) |
| `delta S_act ≈ -0.0388718369` | notes:32 | stage152 notes:38 (`-0.0388718368650403`) | MATCH (truncation) |
| `delta Pi_act ≈ 0.9070844148` | notes:35 | stage152 notes:62 (`0.907084414842908`) | MATCH (truncation) |
| `delta Tm_act ≈ 0.2716539795` | notes:36 | stage152 notes:68 (`0.271653979462338`) | MATCH (truncation) |
| `Pi_corr ≈ 2.4159139283` | notes:42 | stage152 notes:78 (`2.41591392833607`) | MATCH (truncation); also `= Pi_* + delta_Pi_act` |
| `Tm_corr ≈ 1.1731380336` | notes:44 | stage152 notes:86 (`1.17313803363654`) | MATCH (truncation) |
| `g_1 ≈ 0.6844235741` | notes:54 | stage152 notes:112 (`0.684423574065325`) | MATCH (truncation) |
| `S_1 ≈ 0.6163331306` | notes:55 | stage152 notes:113 (`0.616333130570251`) | MATCH (truncation) |
| `Pi_1 ≈ 2.5391484761` | notes:60 | stage152 notes:121 (`2.53914847609768`) | MATCH (truncation) |
| `Tm_1 ≈ 1.2103694208` | notes:61 | stage152 notes:122 (`1.21036942084359`) | MATCH (truncation) |
| `R_*''(0) = -3 Sigma_m^* Pi_*/(1-e^{-Pi_*}) < 0` | notes:20-21 | symbolic; consistent with Stage 151/152 residual setup | MATCH (symbolic; sign verified < 0) |

INTERNAL (scaffolding / definitional, no finding): `Phi_*(x)=4 Sigma_m^* T_s - Sigma_m^* T_q`, `R_*(x)=Phi_*(x)-Pi_* x`, `R_*(0)=R_*'(0)=0`, `Sigma_1(x) ∝ e^{-Pi_* x - R_*(x)}` — definitional/structural carriers, not standalone deliverable numerics.

reconciliation: complete; 11 deliverable values checked, 0 misaligned
